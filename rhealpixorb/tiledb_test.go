package rhealpixorb_test

import (
	"testing"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
	"github.com/sixy6e/go-rhealpix/rhealpixorb"
)

// ============================================================================
// GEOJSON PIPELINE INTEGRATION TESTS
// ============================================================================

func TestBoundingBoxToTileDBRanges64_Integration(t *testing.T) {
	el := rhealpix.NewWGS84()
	targetRes := uint8(12)

	// Canberra / ACT Bounding Box [minLon, minLat, maxLon, maxLat]
	canberraBBox := orb.Bound{
		Min: orb.Point{148.7, -35.6},
		Max: orb.Point{149.4, -35.1},
	}

	ranges, err := rhealpixorb.BoundingBoxToTileDBRanges(el, canberraBBox, targetRes)
	if err != nil {
		t.Fatalf("BoundingBoxToTileDBRanges failed: %v", err)
	}

	if len(ranges) == 0 {
		t.Fatalf("expected non-zero ranges for Canberra BBox, got 0")
	}

	for i, r := range ranges {
		minCell := rhealpix.CellID64(r.Min)
		maxCell := rhealpix.CellID64(r.Max)

		// invariant: Min <= Max in 1D integer space
		if r.Min > r.Max {
			t.Errorf("Range %d has Min > Max: Min=%d (0x%X), Max=%d (0x%X)",
				i, r.Min, r.Min, r.Max, r.Max)
		}

		// invariant: Resolution Header Alignment
		if minCell.Resolution() != targetRes {
			t.Errorf("Range %d Min resolution header mismatch: got %d, want %d (Min Hex: 0x%X)",
				i, minCell.Resolution(), targetRes, r.Min)
		}
		if maxCell.Resolution() != targetRes {
			t.Errorf("Range %d Max resolution header mismatch: got %d, want %d (Max Hex: 0x%X)",
				i, maxCell.Resolution(), targetRes, r.Max)
		}

		// invariant: Single-Range Trunk Uniformity (Min and Max must share the same SUID trunk)
		minTrunk := r.Min >> 52
		maxTrunk := r.Max >> 52
		if minTrunk != maxTrunk {
			t.Errorf("Range %d spans across distinct SUID trunks: MinTrunk=0x%X (0x%X), MaxTrunk=0x%X (0x%X)",
				i, minTrunk, r.Min, maxTrunk, r.Max)
		}
	}

	// invariant: Consecutive Range Trunk Isolation
	for i := 0; i < len(ranges)-1; i++ {
		currMaxTrunk := ranges[i].Max >> 52
		nextMinTrunk := ranges[i+1].Min >> 52

		// if consecutive ranges sit in different trunks, verify they were NOT merged
		if currMaxTrunk != nextMinTrunk {
			// ensures current.Max cannot reach next.Min across different trunk boundaries
			if ranges[i+1].Min <= ranges[i].Max {
				t.Errorf("Ranges %d and %d overlap across trunk boundaries! Range %d Max=0x%X (Trunk 0x%X), Range %d Min=0x%X (Trunk 0x%X)",
					i, i+1, i, ranges[i].Max, currMaxTrunk, i+1, ranges[i+1].Min, nextMinTrunk)
			}
		}
	}
}

func TestBoundingBoxToTileDBRanges128_Integration(t *testing.T) {
	el := rhealpix.NewWGS84()
	targetRes := uint8(20)

	// Canberra / ACT Bounding Box [minLon, minLat, maxLon, maxLat]
	canberraBBox := orb.Bound{
		Min: orb.Point{148.7, -35.6},
		Max: orb.Point{149.4, -35.1},
	}

	ranges, err := rhealpixorb.BoundingBoxToTileDBRanges128(el, canberraBBox, targetRes)
	if err != nil {
		t.Fatalf("BoundingBoxToTileDBRanges128 failed: %v", err)
	}

	if len(ranges) == 0 {
		t.Fatalf("expected non-zero 128-bit ranges for Canberra BBox, got 0")
	}

	for i, r := range ranges {
		minCell := rhealpix.CellID128{High: r.MinHigh, Low: r.MinLow}
		maxCell := rhealpix.CellID128{High: r.MaxHigh, Low: r.MaxLow}

		// invariant: Min <= Max in 128-bit integer space
		minIsGreater := r.MinHigh > r.MaxHigh || (r.MinHigh == r.MaxHigh && r.MinLow > r.MaxLow)
		if minIsGreater {
			t.Errorf("Range 128 %d has Min > Max: Min=%s, Max=%s", i, minCell.Hex(), maxCell.Hex())
		}

		// invariant: Resolution Header Alignment
		if minCell.Resolution() != targetRes {
			t.Errorf("Range 128 %d Min resolution header mismatch: got %d, want %d (Min Hex: %s)",
				i, minCell.Resolution(), targetRes, minCell.Hex())
		}
		if maxCell.Resolution() != targetRes {
			t.Errorf("Range 128 %d Max resolution header mismatch: got %d, want %d (Max Hex: %s)",
				i, maxCell.Resolution(), targetRes, maxCell.Hex())
		}

		// invariant: Single-Range Trunk Uniformity
		minTrunk := r.MinHigh >> 52
		maxTrunk := r.MaxHigh >> 52
		if minTrunk != maxTrunk {
			t.Errorf("Range 128 %d spans across distinct SUID trunks: MinTrunk=0x%X, MaxTrunk=0x%X",
				i, minTrunk, maxTrunk)
		}
	}

	// invariant: Consecutive Range Trunk Isolation
	for i := 0; i < len(ranges)-1; i++ {
		currMaxTrunk := ranges[i].MaxHigh >> 52
		nextMinTrunk := ranges[i+1].MinHigh >> 52

		if currMaxTrunk != nextMinTrunk {
			nextGeCurr := ranges[i+1].MinHigh > ranges[i].MaxHigh ||
				(ranges[i+1].MinHigh == ranges[i].MaxHigh && ranges[i+1].MinLow >= ranges[i].MaxLow)

			if !nextGeCurr {
				t.Errorf("128-bit Ranges %d and %d overlap across trunk boundaries!", i, i+1)
			}
		}
	}
}
