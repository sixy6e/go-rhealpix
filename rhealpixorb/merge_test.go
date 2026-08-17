package rhealpixorb

import (
	"reflect"
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
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

// ============================================================================
// RANGE MERGING SECURITY TESTS
// ============================================================================

func TestMergeRangesDisjointTrunks(t *testing.T) {
	gap := uint64(100)

	t.Run("64-bit Disjoint Trunk Isolation", func(t *testing.T) {
		// R7... vs R8... trunks at Res 12
		cellR7, _ := rhealpix.PackCellID64(4, 1, []uint8{7})
		cellR8, _ := rhealpix.PackCellID64(4, 1, []uint8{8})

		r7Min, r7Max := cellR7.SubtreeRange(12)
		r8Min, r8Max := cellR8.SubtreeRange(12)

		ranges := []Uint64Range{
			{Min: uint64(r7Min), Max: uint64(r7Max)},
			{Min: uint64(r8Min), Max: uint64(r8Max)},
		}

		merged := MergeRangesWithGap(ranges, gap)
		if len(merged) != 2 {
			t.Fatalf("expected 2 unmerged ranges across disjoint trunks, got %d", len(merged))
		}
	})

	t.Run("128-bit Disjoint Trunk Isolation", func(t *testing.T) {
		cellR7, _ := rhealpix.PackCellID128(4, 1, []uint8{7})
		cellR8, _ := rhealpix.PackCellID128(4, 1, []uint8{8})

		r7Min, r7Max := cellR7.SubtreeRange(20)
		r8Min, r8Max := cellR8.SubtreeRange(20)

		ranges := []Uint128Range{
			{MinHigh: r7Min.High, MinLow: r7Min.Low, MaxHigh: r7Max.High, MaxLow: r7Max.Low},
			{MinHigh: r8Min.High, MinLow: r8Min.Low, MaxHigh: r8Max.High, MaxLow: r8Max.Low},
		}

		merged := MergeRangesWithGap128(ranges, gap)
		if len(merged) != 2 {
			t.Fatalf("expected 2 unmerged 128-bit ranges across disjoint trunks, got %d", len(merged))
		}
	})
}

func TestMergeRangesResolutionIsolation(t *testing.T) {
	gap := uint64(100)

	t.Run("64-bit Resolution Level Isolation", func(t *testing.T) {
		cell, _ := rhealpix.PackCellID64(4, 1, []uint8{7})
		rRes10Min, rRes10Max := cell.SubtreeRange(10)
		rRes12Min, rRes12Max := cell.SubtreeRange(12)

		ranges := []Uint64Range{
			{Min: uint64(rRes10Min), Max: uint64(rRes10Max)},
			{Min: uint64(rRes12Min), Max: uint64(rRes12Max)},
		}

		merged := MergeRangesWithGap(ranges, gap)
		if len(merged) != 2 {
			t.Fatalf("expected 2 unmerged ranges across different resolution headers, got %d", len(merged))
		}
	})

	t.Run("128-bit Resolution Level Isolation", func(t *testing.T) {
		cell, _ := rhealpix.PackCellID128(4, 1, []uint8{7})
		rRes10Min, rRes10Max := cell.SubtreeRange(10)
		rRes12Min, rRes12Max := cell.SubtreeRange(12)

		ranges := []Uint128Range{
			{MinHigh: rRes10Min.High, MinLow: rRes10Min.Low, MaxHigh: rRes10Max.High, MaxLow: rRes10Max.Low},
			{MinHigh: rRes12Min.High, MinLow: rRes12Min.Low, MaxHigh: rRes12Max.High, MaxLow: rRes12Max.Low},
		}

		merged := MergeRangesWithGap128(ranges, gap)
		if len(merged) != 2 {
			t.Fatalf("expected 2 unmerged 128-bit ranges across different resolution headers, got %d", len(merged))
		}
	})
}
