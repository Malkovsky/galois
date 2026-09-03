// Standalone top-level Reed-Solomon benchmark for klauspost/reedsolomon.
//
// This runner intentionally avoids CGo. It emits Google Benchmark-compatible
// JSON so its verified single-threaded results can be merged with the C++
// reference suite after collection on the same pinned host.
package main

import (
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"math"
	"os"
	"regexp"
	"runtime"
	"sort"
	"time"

	"github.com/klauspost/reedsolomon"
	"golang.org/x/sys/unix"
)

const (
	moduleVersion = "v1.14.2"
	moduleCommit  = "af9e2b1b1bad1889954523347758996aafd9c805"
	shardBytes    = 1024
	defaultSeed   = int64(0x6d2b79f5)
)

type benchmarkCase struct {
	data     int
	recovery int
}

var comparisonCases = []benchmarkCase{
	{8, 248},
	{16, 240},
	{32, 224},
	{64, 192},
	{128, 128},
	{192, 64},
	{224, 32},
	{240, 16},
	{248, 8},
}

type operation int

const (
	encode operation = iota
	decodeMax
)

type benchmarkInstance struct {
	caseIndex  int
	kind       operation
	name       string
	iterations int
}

type outputDocument struct {
	Context    map[string]any `json:"context"`
	Benchmarks []outputRow    `json:"benchmarks"`
}

type outputRow struct {
	Name                   string  `json:"name"`
	FamilyIndex            int     `json:"family_index"`
	PerFamilyInstanceIndex int     `json:"per_family_instance_index"`
	RunName                string  `json:"run_name"`
	RunType                string  `json:"run_type"`
	Repetitions            int     `json:"repetitions"`
	RepetitionIndex        *int    `json:"repetition_index,omitempty"`
	Threads                int     `json:"threads"`
	Iterations             int     `json:"iterations"`
	RealTime               float64 `json:"real_time"`
	CPUTime                float64 `json:"cpu_time"`
	TimeUnit               string  `json:"time_unit"`
	BytesPerSecond         float64 `json:"bytes_per_second"`
	AggregateName          string  `json:"aggregate_name,omitempty"`
	AggregateUnit          string  `json:"aggregate_unit,omitempty"`
}

type measurement struct {
	iterations     int
	realNS         float64
	cpuNS          float64
	bytesPerSecond float64
}

type fixture struct {
	encoder      reedsolomon.Encoder
	shards       [][]byte
	originalData [][]byte
	erased       []bool
}

func nextRandom(state *uint32) uint32 {
	*state ^= *state << 13
	*state ^= *state >> 17
	*state ^= *state << 5
	return *state
}

// maxErasurePattern exactly mirrors benchmarks/rs_benchmark_cases.h. Its
// logical order is [recovery][data].
func maxErasurePattern(data, recovery int) []bool {
	erased := make([]bool, data+recovery)
	for i := 0; i < recovery; i++ {
		erased[i] = true
	}

	state := uint32(0x6d2b79f5) ^ uint32(data*257)
	for i := len(erased) - 1; i > 0; i-- {
		position := int(nextRandom(&state) % uint32(i+1))
		erased[i], erased[position] = erased[position], erased[i]
	}
	return erased
}

// physicalErasurePattern maps [recovery][data] into klauspost's
// [data][parity] order.
func physicalErasurePattern(data, recovery int) []bool {
	logical := maxErasurePattern(data, recovery)
	physical := make([]bool, data+recovery)
	for i := 0; i < data; i++ {
		physical[i] = logical[recovery+i]
	}
	for i := 0; i < recovery; i++ {
		physical[data+i] = logical[i]
	}
	return physical
}

func newFixture(testCase benchmarkCase) (*fixture, error) {
	encoder, err := reedsolomon.New(
		testCase.data,
		testCase.recovery,
		reedsolomon.WithMaxGoroutines(1),
		reedsolomon.WithInversionCache(false),
	)
	if err != nil {
		return nil, err
	}

	extensions, ok := encoder.(reedsolomon.Extensions)
	if !ok {
		return nil, errors.New("klauspost encoder does not expose aligned allocation")
	}
	shards := extensions.AllocAligned(shardBytes)
	if len(shards) != testCase.data+testCase.recovery {
		return nil, fmt.Errorf("aligned allocation returned %d shards", len(shards))
	}

	for shard := 0; shard < testCase.data; shard++ {
		for offset := range shards[shard] {
			shards[shard][offset] = byte(
				(shard*131 + offset*17 + testCase.data*3 + testCase.recovery) & 255,
			)
		}
	}
	if err := encoder.Encode(shards); err != nil {
		return nil, err
	}
	valid, err := encoder.Verify(shards)
	if err != nil || !valid {
		return nil, fmt.Errorf("encoded fixture verification failed: valid=%v err=%v", valid, err)
	}

	originalData := make([][]byte, testCase.data)
	for i := range originalData {
		originalData[i] = append([]byte(nil), shards[i]...)
	}

	return &fixture{
		encoder:      encoder,
		shards:       shards,
		originalData: originalData,
		erased:       physicalErasurePattern(testCase.data, testCase.recovery),
	}, nil
}

func (f *fixture) makeDecodeWorkSets(iterations int) [][][]byte {
	sets := make([][][]byte, iterations)
	for iteration := range sets {
		sets[iteration] = make([][]byte, len(f.shards))
	}
	f.resetDecodeWorkSets(sets)
	return sets
}

func (f *fixture) resetDecodeWorkSets(sets [][][]byte) {
	for _, work := range sets {
		for shard := range work {
			if f.erased[shard] {
				work[shard] = f.shards[shard][:0]
			} else {
				work[shard] = f.shards[shard]
			}
		}
	}
}

func (f *fixture) poisonMissingData() {
	for shard := range f.originalData {
		if !f.erased[shard] {
			continue
		}
		for offset := range f.shards[shard] {
			f.shards[shard][offset] = 0xa5
		}
	}
}

func (f *fixture) verifyData() error {
	for shard := range f.originalData {
		if len(f.shards[shard]) != len(f.originalData[shard]) {
			return fmt.Errorf("data shard %d has length %d", shard, len(f.shards[shard]))
		}
		for offset := range f.originalData[shard] {
			if f.shards[shard][offset] != f.originalData[shard][offset] {
				return fmt.Errorf("data shard %d differs at byte %d", shard, offset)
			}
		}
	}
	return nil
}

func reconstructBatch(
	workSets [][][]byte,
	count int,
	reconstruct func([][]byte) error,
) error {
	for iteration := range workSets[:count] {
		if err := reconstruct(workSets[iteration]); err != nil {
			return err
		}
	}
	return nil
}

func threadCPUTime() (time.Duration, error) {
	var value unix.Timespec
	if err := unix.ClockGettime(unix.CLOCK_THREAD_CPUTIME_ID, &value); err != nil {
		return 0, err
	}
	return time.Duration(value.Sec)*time.Second + time.Duration(value.Nsec), nil
}

func cpuPinning() string {
	var affinity unix.CPUSet
	if err := unix.SchedGetaffinity(0, &affinity); err != nil {
		return "unknown"
	}
	if affinity.Count() == 1 {
		return "pinned"
	}
	return "not_pinned"
}

func measure(instance benchmarkInstance, iterations int) (measurement, error) {
	testCase := comparisonCases[instance.caseIndex]
	fixture, err := newFixture(testCase)
	if err != nil {
		return measurement{}, err
	}

	var wallElapsed, cpuElapsed time.Duration
	var lastDecodeSet [][]byte
	if instance.kind == encode {
		cpuStart, err := threadCPUTime()
		if err != nil {
			return measurement{}, err
		}
		wallStart := time.Now()
		for iteration := 0; iteration < iterations; iteration++ {
			if err = fixture.encoder.Encode(fixture.shards); err != nil {
				return measurement{}, err
			}
		}
		wallElapsed = time.Since(wallStart)
		cpuEnd, err := threadCPUTime()
		if err != nil {
			return measurement{}, err
		}
		cpuElapsed = cpuEnd - cpuStart
	} else {
		// ReconstructData restores missing slice lengths. Run bounded batches so
		// header reset and allocation stay outside the accumulated timed regions
		// without allocating one 256-entry header array per total iteration.
		const batchSize = 512
		fixture.poisonMissingData()
		workSets := fixture.makeDecodeWorkSets(min(iterations, batchSize))
		remaining := iterations
		for remaining > 0 {
			count := min(remaining, batchSize)
			fixture.resetDecodeWorkSets(workSets[:count])
			cpuStart, err := threadCPUTime()
			if err != nil {
				return measurement{}, err
			}
			wallStart := time.Now()
			if err = reconstructBatch(
				workSets, count, fixture.encoder.ReconstructData,
			); err != nil {
				return measurement{}, err
			}
			wallElapsed += time.Since(wallStart)
			cpuEnd, err := threadCPUTime()
			if err != nil {
				return measurement{}, err
			}
			cpuElapsed += cpuEnd - cpuStart
			lastDecodeSet = workSets[count-1]
			remaining -= count
		}
	}

	if instance.kind == encode {
		valid, verifyErr := fixture.encoder.Verify(fixture.shards)
		if verifyErr != nil || !valid {
			return measurement{}, fmt.Errorf(
				"post-encode verification failed: valid=%v err=%v", valid, verifyErr,
			)
		}
	} else {
		if err := fixture.verifyData(); err != nil {
			return measurement{}, err
		}
		for shard, erased := range fixture.erased {
			if erased && shard < testCase.data && len(lastDecodeSet[shard]) != shardBytes {
				return measurement{}, fmt.Errorf("data shard %d was not reconstructed", shard)
			}
			if erased && shard >= testCase.data && len(lastDecodeSet[shard]) != 0 {
				return measurement{}, fmt.Errorf("parity shard %d was unexpectedly reconstructed", shard)
			}
		}
	}

	cpuPerIteration := float64(cpuElapsed.Nanoseconds()) / float64(iterations)
	return measurement{
		iterations:     iterations,
		realNS:         float64(wallElapsed.Nanoseconds()) / float64(iterations),
		cpuNS:          cpuPerIteration,
		bytesPerSecond: float64(testCase.data*shardBytes) / (cpuPerIteration * 1e-9),
	}, nil
}

func warm(instance benchmarkInstance, duration time.Duration) error {
	if duration <= 0 {
		return nil
	}
	const chunkIterations = 256
	var elapsed time.Duration
	for elapsed < duration {
		measurement, err := measure(instance, chunkIterations)
		if err != nil {
			return err
		}
		elapsed += time.Duration(measurement.cpuNS * chunkIterations)
	}
	return nil
}

func calibrate(instance benchmarkInstance, target time.Duration) (int, error) {
	if target <= 0 {
		return 1, nil
	}
	iterations := 1
	for {
		measurement, err := measure(instance, iterations)
		if err != nil {
			return 0, err
		}
		elapsed := measurement.cpuNS * float64(iterations)
		if elapsed >= float64(target.Nanoseconds()) {
			return iterations, nil
		}
		estimate := int(math.Ceil(float64(iterations) * float64(target.Nanoseconds()) / math.Max(elapsed, 1) * 1.15))
		if estimate < iterations*2 {
			estimate = iterations * 2
		}
		if estimate > iterations*10 {
			estimate = iterations * 10
		}
		if estimate == iterations {
			return iterations, nil
		}
		iterations = estimate
	}
}

func instanceName(kind operation, testCase benchmarkCase) string {
	operationName := "Encode"
	if kind == decodeMax {
		operationName = "DecodeMax"
	}
	return fmt.Sprintf(
		"RS/Klauspost/Native/%s/K:%d/R:%d/bytes:%d",
		operationName,
		testCase.data,
		testCase.recovery,
		shardBytes,
	)
}

func makeInstances(filter *regexp.Regexp) []benchmarkInstance {
	instances := make([]benchmarkInstance, 0, len(comparisonCases)*2)
	for _, kind := range []operation{encode, decodeMax} {
		for caseIndex, testCase := range comparisonCases {
			name := instanceName(kind, testCase)
			if filter == nil || filter.MatchString(name) {
				instances = append(instances, benchmarkInstance{
					caseIndex: caseIndex,
					kind:      kind,
					name:      name,
				})
			}
		}
	}
	return instances
}

type splitMix64 struct {
	state uint64
}

func (random *splitMix64) next() uint64 {
	random.state += 0x9e3779b97f4a7c15
	value := random.state
	value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9
	value = (value ^ (value >> 27)) * 0x94d049bb133111eb
	return value ^ (value >> 31)
}

func shuffle(values []int, seed int64) {
	random := splitMix64{state: uint64(seed)}
	for i := len(values) - 1; i > 0; i-- {
		position := int(random.next() % uint64(i+1))
		values[i], values[position] = values[position], values[i]
	}
}

func statistic(values []float64, name string) float64 {
	switch name {
	case "mean":
		var total float64
		for _, value := range values {
			total += value
		}
		return total / float64(len(values))
	case "median":
		copyValues := append([]float64(nil), values...)
		sort.Float64s(copyValues)
		middle := len(copyValues) / 2
		if len(copyValues)%2 == 1 {
			return copyValues[middle]
		}
		return (copyValues[middle-1] + copyValues[middle]) / 2
	case "stddev":
		if len(values) < 2 {
			return 0
		}
		mean := statistic(values, "mean")
		var squares float64
		for _, value := range values {
			difference := value - mean
			squares += difference * difference
		}
		return math.Sqrt(squares / float64(len(values)-1))
	case "cv":
		mean := statistic(values, "mean")
		if mean == 0 {
			return 0
		}
		return statistic(values, "stddev") / mean
	default:
		panic("unknown statistic " + name)
	}
}

func appendAggregates(document *outputDocument, instance benchmarkInstance, rows []outputRow) {
	realValues := make([]float64, len(rows))
	cpuValues := make([]float64, len(rows))
	throughputValues := make([]float64, len(rows))
	for i, row := range rows {
		realValues[i] = row.RealTime
		cpuValues[i] = row.CPUTime
		throughputValues[i] = row.BytesPerSecond
	}

	for _, aggregate := range []string{"mean", "median", "stddev", "cv"} {
		aggregateUnit := "time"
		if aggregate == "cv" {
			aggregateUnit = "percentage"
		}
		document.Benchmarks = append(document.Benchmarks, outputRow{
			Name:                   instance.name + "_" + aggregate,
			FamilyIndex:            int(instance.kind),
			PerFamilyInstanceIndex: instance.caseIndex,
			RunName:                instance.name,
			RunType:                "aggregate",
			Repetitions:            len(rows),
			Threads:                1,
			Iterations:             len(rows),
			RealTime:               statistic(realValues, aggregate),
			CPUTime:                statistic(cpuValues, aggregate),
			TimeUnit:               "ns",
			BytesPerSecond:         statistic(throughputValues, aggregate),
			AggregateName:          aggregate,
			AggregateUnit:          aggregateUnit,
		})
	}
}

func main() {
	runtime.GOMAXPROCS(1)
	runtime.LockOSThread()
	defer runtime.UnlockOSThread()

	var (
		filterText      string
		listOnly        bool
		outputPath      string
		repetitions     int
		warmup          time.Duration
		minimumTime     time.Duration
		fixedIterations int
		seed            int64
	)
	flag.StringVar(&filterText, "filter", "", "regular expression selecting benchmark names")
	flag.BoolVar(&listOnly, "list", false, "list matching benchmark names and exit")
	flag.StringVar(&outputPath, "output", "", "write Google Benchmark-compatible JSON to this file")
	flag.IntVar(&repetitions, "repetitions", 15, "number of measured repetitions per benchmark")
	flag.DurationVar(&warmup, "warmup", 100*time.Millisecond, "discarded CPU warmup per benchmark")
	flag.DurationVar(&minimumTime, "min-time", 100*time.Millisecond, "minimum CPU time per repetition")
	flag.IntVar(&fixedIterations, "iterations", 0, "fixed iterations per repetition; skips warmup and calibration")
	flag.Int64Var(&seed, "seed", defaultSeed, "random interleaving seed")
	flag.Parse()

	if repetitions <= 0 || fixedIterations < 0 {
		fmt.Fprintln(os.Stderr, "repetitions must be positive and iterations non-negative")
		os.Exit(2)
	}

	var filter *regexp.Regexp
	var err error
	if filterText != "" {
		filter, err = regexp.Compile(filterText)
		if err != nil {
			fmt.Fprintln(os.Stderr, "invalid filter:", err)
			os.Exit(2)
		}
	}
	instances := makeInstances(filter)
	if listOnly {
		for _, instance := range instances {
			fmt.Println(instance.name)
		}
		return
	}
	if len(instances) == 0 {
		fmt.Fprintln(os.Stderr, "no benchmarks matched")
		os.Exit(2)
	}

	hostName, _ := os.Hostname()
	document := outputDocument{
		Context: map[string]any{
			"date":                time.Now().Format(time.RFC3339),
			"host_name":           hostName,
			"executable":          os.Args[0],
			"num_cpus":            runtime.NumCPU(),
			"cpu_pinning":         cpuPinning(),
			"library_version":     "standalone-go-runner",
			"library_build_type":  "release",
			"json_schema_version": 1,
			"go_version":          runtime.Version(),
			"go_os":               runtime.GOOS,
			"go_arch":             runtime.GOARCH,
			"reedsolomon_module":  "github.com/klauspost/reedsolomon",
			"reedsolomon_version": moduleVersion,
			"reedsolomon_commit":  moduleCommit,
			"gomaxprocs":          1,
			"max_goroutines":      1,
			"inversion_cache":     false,
			"interleave_seed":     seed,
		},
	}

	schedule := make([]int, 0, len(instances)*repetitions)
	for repetition := 0; repetition < repetitions; repetition++ {
		for instanceIndex := range instances {
			schedule = append(schedule, instanceIndex)
		}
	}
	shuffle(schedule, seed)

	rowsByInstance := make([][]outputRow, len(instances))
	for _, instanceIndex := range schedule {
		instance := &instances[instanceIndex]
		if instance.iterations == 0 {
			if fixedIterations > 0 {
				instance.iterations = fixedIterations
			} else {
				if err := warm(*instance, warmup); err != nil {
					fmt.Fprintln(os.Stderr, instance.name+":", err)
					os.Exit(1)
				}
				instance.iterations, err = calibrate(*instance, minimumTime)
				if err != nil {
					fmt.Fprintln(os.Stderr, instance.name+":", err)
					os.Exit(1)
				}
			}
		}

		measurement, err := measure(*instance, instance.iterations)
		if err != nil {
			fmt.Fprintln(os.Stderr, instance.name+":", err)
			os.Exit(1)
		}
		repetitionIndex := len(rowsByInstance[instanceIndex])
		row := outputRow{
			Name:                   instance.name,
			FamilyIndex:            int(instance.kind),
			PerFamilyInstanceIndex: instance.caseIndex,
			RunName:                instance.name,
			RunType:                "iteration",
			Repetitions:            repetitions,
			RepetitionIndex:        &repetitionIndex,
			Threads:                1,
			Iterations:             measurement.iterations,
			RealTime:               measurement.realNS,
			CPUTime:                measurement.cpuNS,
			TimeUnit:               "ns",
			BytesPerSecond:         measurement.bytesPerSecond,
		}
		document.Benchmarks = append(document.Benchmarks, row)
		rowsByInstance[instanceIndex] = append(rowsByInstance[instanceIndex], row)
	}

	for instanceIndex, instance := range instances {
		appendAggregates(&document, instance, rowsByInstance[instanceIndex])
	}

	encoded, err := json.MarshalIndent(document, "", "  ")
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	encoded = append(encoded, '\n')
	if outputPath == "" {
		_, err = os.Stdout.Write(encoded)
	} else {
		err = os.WriteFile(outputPath, encoded, 0o644)
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}
