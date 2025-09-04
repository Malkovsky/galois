build/Release/benchmarks \
    --benchmark_format=json \
    --benchmark_out=benchmarks/benchmarks.json

python scripts/parse_benchmarks.py \
    --metric cpu_time \
    --output docs/images/benchmarks.svg \
    --title "Matrix multiplication" \
    benchmarks/benchmarks.json
