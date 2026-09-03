package main

import "testing"

func TestReconstructBatchUsesOnlyLivePartialBatch(t *testing.T) {
	workSets := make([][][]byte, 512)
	calls := 0
	err := reconstructBatch(workSets, 1, func([][]byte) error {
		calls++
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	if calls != 1 {
		t.Fatalf("got %d reconstruction calls, want 1", calls)
	}
}

func TestMaxErasurePatternMatchesSharedCases(t *testing.T) {
	expectedMissingData := []int{7, 16, 28, 49, 55, 49, 28, 15, 8}
	expectedMissingRecovery := []int{241, 224, 196, 143, 73, 15, 4, 1, 0}
	expectedPatternHash := []uint64{
		0x3deebfec14a438dd,
		0x88ba2078ea24ac03,
		0x6bd8a738779e020d,
		0x3c9666b43a92059d,
		0xfc7f499b142ff235,
		0x698e1a23b5dbc50f,
		0xd830fecf61e12bcb,
		0x83f4651310e2e295,
		0x3745ca6111157135,
	}

	for caseIndex, testCase := range comparisonCases {
		pattern := physicalErasurePattern(testCase.data, testCase.recovery)
		if len(pattern) != testCase.data+testCase.recovery {
			t.Fatalf("case %d: got pattern length %d", caseIndex, len(pattern))
		}

		missingData := 0
		missingRecovery := 0
		patternHash := uint64(14695981039346656037)
		for shard, erased := range pattern {
			value := byte(0)
			if erased {
				value = 1
			}
			patternHash ^= uint64(value)
			patternHash *= 1099511628211
			if !erased {
				continue
			}
			if shard < testCase.data {
				missingData++
			} else {
				missingRecovery++
			}
		}

		if missingData != expectedMissingData[caseIndex] ||
			missingRecovery != expectedMissingRecovery[caseIndex] {
			t.Fatalf(
				"case %d: got data/recovery losses %d/%d, want %d/%d",
				caseIndex,
				missingData,
				missingRecovery,
				expectedMissingData[caseIndex],
				expectedMissingRecovery[caseIndex],
			)
		}
		if missingData+missingRecovery != testCase.recovery {
			t.Fatalf("case %d: erasure total does not equal R", caseIndex)
		}
		if patternHash != expectedPatternHash[caseIndex] {
			t.Fatalf(
				"case %d: got erasure fingerprint %#x, want %#x",
				caseIndex,
				patternHash,
				expectedPatternHash[caseIndex],
			)
		}
	}
}
