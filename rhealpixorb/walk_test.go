package rhealpixorb_test

import (
	"testing"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
	"github.com/sixy6e/go-rhealpix/rhealpixorb"
)

func TestSTACGeometryToTileDBRangesTopDown_Canberra(t *testing.T) {
	el := rhealpix.NewWGS84()

	// Canberra / ACT BBox: approx [148.7, -35.6, 149.4, -35.1]
	canberraBBox := orb.Bound{
		Min: orb.Point{148.7, -35.6},
		Max: orb.Point{149.4, -35.1},
	}.ToPolygon()

	// decompose down to Level 7 (~10km resolution)
	targetRes := uint8(7)

	compactedCells, err := rhealpixorb.STACGeometryToTileDBRangesTopDown(el, canberraBBox, targetRes)
	if err != nil {
		t.Fatalf("Top-down decomposition failed: %v", err)
	}

	if len(compactedCells) == 0 {
		t.Fatal("expected non-empty cell list for Canberra geometry")
	}

	t.Logf("Canberra BBox generated %d compacted cell ranges at Level %d", len(compactedCells), targetRes)

	// verify all returned cells sit inside valid resolution bounds
	for _, cell := range compactedCells {
		if cell.Resolution() > targetRes {
			t.Errorf("cell resolution %d exceeds target resolution %d", cell.Resolution(), targetRes)
		}
	}
}

func BenchmarkSTACGeometryToTileDBRangesTopDown_Canberra(b *testing.B) {
	el := rhealpix.NewWGS84()
	canberraBBox := orb.Bound{
		Min: orb.Point{148.7, -35.6},
		Max: orb.Point{149.4, -35.1},
	}.ToPolygon()
	targetRes := uint8(7)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := rhealpixorb.STACGeometryToTileDBRangesTopDown(el, canberraBBox, targetRes)
		if err != nil {
			b.Fatalf("Decomposition failed: %v", err)
		}
	}
}
