package rhealpixorb

import (
	"reflect"
	"testing"
	// rhealpix "github.com/sixy6e/go-rhealpix"
)

func TestMergeRangesWithGap128(t *testing.T) {
	tests := []struct {
		name     string
		ranges   []Uint128Range
		maxGap   uint64
		expected []Uint128Range
	}{
		{
			name:     "Nil or Empty Slices",
			ranges:   nil,
			maxGap:   1,
			expected: nil,
		},
		{
			name: "Single Range",
			ranges: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 20},
			},
			maxGap: 1,
			expected: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 20},
			},
		},
		{
			name: "Disjoint Ranges Outside Gap",
			ranges: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 20},
				{MinHigh: 0, MinLow: 25, MaxHigh: 0, MaxLow: 30},
			},
			maxGap: 1,
			expected: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 20},
				{MinHigh: 0, MinLow: 25, MaxHigh: 0, MaxLow: 30},
			},
		},
		{
			name: "Contiguous Ranges (Gap = 1)",
			ranges: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 20},
				{MinHigh: 0, MinLow: 21, MaxHigh: 0, MaxLow: 30},
			},
			maxGap: 1,
			expected: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 30},
			},
		},
		{
			name: "Unsorted & Overlapping Ranges",
			ranges: []Uint128Range{
				{MinHigh: 0, MinLow: 50, MaxHigh: 0, MaxLow: 100},
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 30},
				{MinHigh: 0, MinLow: 25, MaxHigh: 0, MaxLow: 60},
			},
			maxGap: 1,
			expected: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 100},
			},
		},
		{
			name: "64-bit Low Word Carry Across Boundary (Low Overflow)",
			// first range ends near max uint64. Adding maxGap=5 causes Low to roll over to 3 and High to increment to 1.
			// second range starts at High=1, Low=2, which falls within the gap
			ranges: []Uint128Range{
				{MinHigh: 0, MinLow: 100, MaxHigh: 0, MaxLow: 0xFFFFFFFFFFFFFFFE},
				{MinHigh: 1, MinLow: 2, MaxHigh: 1, MaxLow: 50},
			},
			maxGap: 5,
			expected: []Uint128Range{
				{MinHigh: 0, MinLow: 100, MaxHigh: 1, MaxLow: 50},
			},
		},
		{
			name: "128-bit Sorting Across High and Low Fields",
			ranges: []Uint128Range{
				{MinHigh: 2, MinLow: 5, MaxHigh: 2, MaxLow: 10},
				{MinHigh: 1, MinLow: 500, MaxHigh: 1, MaxLow: 600},
				{MinHigh: 1, MinLow: 10, MaxHigh: 1, MaxLow: 100},
			},
			maxGap: 1,
			expected: []Uint128Range{
				{MinHigh: 1, MinLow: 10, MaxHigh: 1, MaxLow: 100},
				{MinHigh: 1, MinLow: 500, MaxHigh: 1, MaxLow: 600},
				{MinHigh: 2, MinLow: 5, MaxHigh: 2, MaxLow: 10},
			},
		},
		{
			name: "Fully Subsumed Range",
			ranges: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 100},
				{MinHigh: 0, MinLow: 25, MaxHigh: 0, MaxLow: 50},
			},
			maxGap: 0,
			expected: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 100},
			},
		},
		{
			name: "Exact Contiguous (maxGap = 0)",
			ranges: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 20},
				{MinHigh: 0, MinLow: 21, MaxHigh: 0, MaxLow: 30},
			},
			maxGap: 0,
			expected: []Uint128Range{
				{MinHigh: 0, MinLow: 10, MaxHigh: 0, MaxLow: 20},
				{MinHigh: 0, MinLow: 21, MaxHigh: 0, MaxLow: 30},
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := MergeRangesWithGap128(tt.ranges, tt.maxGap)
			if !reflect.DeepEqual(got, tt.expected) {
				t.Errorf("MergeRangesWithGap128() = %v, want %v", got, tt.expected)
			}
		})
	}
}

func TestAdd128Carry(t *testing.T) {
	// test boundary arithmetic helper directly
	high, low := add128(0, 0xFFFFFFFFFFFFFFFF, 1)
	if high != 1 || low != 0 {
		t.Errorf("add128 overflow failed: got (%d, %d), want (1, 0)", high, low)
	}

	high, low = add128(5, 0xFFFFFFFFFFFFFFFE, 3)
	if high != 6 || low != 1 {
		t.Errorf("add128 overflow with gap failed: got (%d, %d), want (6, 1)", high, low)
	}
}
