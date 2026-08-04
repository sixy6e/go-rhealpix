package rhealpixorb_test

import (
	"testing"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
	"github.com/sixy6e/go-rhealpix/rhealpixorb"
)

func TestDeriveRegionCodeSingleRoot(t *testing.T) {
	el := rhealpix.NewWGS84()

	// Canberra area bounding box (wholly inside Equatorial Facet Q / Index 3)
	bound := orb.Bound{
		Min: orb.Point{149.0, -35.5},
		Max: orb.Point{149.5, -35.0},
	}

	regionCode, suids, err := rhealpixorb.DeriveRegionCode64FromGeometry(el, bound, 10)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if len(suids) != 1 {
		t.Errorf("expected 1 root cell region, got %d (%v)", len(suids), suids)
	}

	if regionCode == "" {
		t.Errorf("expected non-empty region_code string")
	}
}

func TestDeriveRegionCodeDeterministicSort(t *testing.T) {
	el := rhealpix.NewWGS84()

	// bounding box crossing equatorial boundary (between Facet P and Facet Q across 0 deg Meridian)
	bound := orb.Bound{
		Min: orb.Point{-1.0, -10.0},
		Max: orb.Point{1.0, -5.0},
	}

	regionCode, suids, err := rhealpixorb.DeriveRegionCode64FromGeometry(el, bound, 8)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if len(suids) < 2 {
		t.Logf("Cross-boundary result: %s (components: %v)", regionCode, suids)
	}
}

func TestCell64ToPolygon(t *testing.T) {
	el := rhealpix.NewWGS84()

	cell, err := rhealpix.ParseCellID64("Q012")
	if err != nil {
		t.Fatalf("unexpected parse error: %v", err)
	}

	// standard 4-corner polygon (segmentCount = 1)
	poly, err := rhealpixorb.Cell64ToPolygon(el, cell, 1)
	if err != nil {
		t.Fatalf("unexpected polygon conversion error: %v", err)
	}

	// 4 corners + 1 closure point = 5 points
	if len(poly) == 0 || len(poly[0]) != 5 {
		t.Errorf("expected closed 5-point ring for segmentCount=1, got %d points", len(poly[0]))
	}

	// densified polygon (segmentCount = 4) -> 4 corners * 4 segments + 1 closure point = 17 points
	polyDensified, err := rhealpixorb.Cell64ToPolygon(el, cell, 4)
	if err != nil {
		t.Fatalf("unexpected densified polygon conversion error: %v", err)
	}

	expectedPoints := (4 * 4) + 1
	if len(polyDensified[0]) != expectedPoints {
		t.Errorf("expected %d points for segmentCount=4, got %d", expectedPoints, len(polyDensified[0]))
	}
}

func TestPolarCellToPolygon(t *testing.T) {
	el := rhealpix.NewWGS84()

	// parse a cell on North Polar Facet N (Facet 0)
	cell, err := rhealpix.ParseCellID64("N0")
	if err != nil {
		t.Fatalf("unexpected parse error for polar cell: %v", err)
	}

	poly, err := rhealpixorb.Cell64ToPolygon(el, cell, 1)
	if err != nil {
		t.Fatalf("unexpected polar polygon conversion error: %v", err)
	}

	if len(poly) == 0 || len(poly[0]) < 4 {
		t.Errorf("expected valid closed polar polygon ring, got %v", poly)
	}
}

func TestPointToCellGroundTruth(t *testing.T) {
	// CellID's (SUID's) were calculated with rhealpixdggs_py using cell_from_point output
	el := rhealpix.NewWGS84()

	tests := []struct {
		name       string
		lon, lat   float64
		targetRes  uint8
		expectedID string
	}{
		{
			name:       "Canberra City (Equatorial Facet Q)",
			lon:        149.13,
			lat:        -35.28,
			targetRes:  10,
			expectedID: "R7852345212",
		},
		{
			name:       "Greenwich Meridian (Facet Boundary P/Q)",
			lon:        0.0,
			lat:        51.4778,
			targetRes:  5,
			expectedID: "N22646",
		},
		{
			name:       "Granule Spanning Equatorial Meridian (P/Q boundary) LL Point",
			lon:        -0.2,
			lat:        10.0,
			targetRes:  8,
			expectedID: "P52288776",
		},
		{
			name:       "Granule Spanning Equatorial Meridian (P/Q boundary) LR Point",
			lon:        -0.2,
			lat:        10.5,
			targetRes:  8,
			expectedID: "P52285416",
		},
		{
			name:       "Granule Spanning Equatorial Meridian (P/Q boundary) UR Point",
			lon:        0.2,
			lat:        10.5,
			targetRes:  8,
			expectedID: "Q30063418",
		},
		{
			name:       "Granule Spanning Equatorial Meridian (P/Q boundary) UL Point",
			lon:        0.2,
			lat:        10.0,
			targetRes:  8,
			expectedID: "Q30066778",
		},
		{
			name:       "Polar Satellite Pass (North Cap) Point 0",
			lon:        0.0,
			lat:        85.0,
			targetRes:  6,
			expectedID: "N442224",
		},
		{
			name:       "Polar Satellite Pass (North Cap) Point 0",
			lon:        45.0,
			lat:        88.0,
			targetRes:  6,
			expectedID: "N441771",
		},
		{
			name:       "Polar Satellite Pass (North Cap) Point 0",
			lon:        90.0,
			lat:        85.0,
			targetRes:  6,
			expectedID: "N440004",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			cell, err := rhealpix.ForwardTransform64(el, tt.lon, tt.lat, tt.targetRes)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if cell.String() != tt.expectedID {
				t.Errorf("mismatch: got %s, expected %s", cell.String(), tt.expectedID)
			}
		})
	}
}

func TestSatelliteGranuleIndexingEdgeCases(t *testing.T) {
	el := rhealpix.NewWGS84()

	t.Run("Granule Spanning Equatorial Meridian (P/Q boundary)", func(t *testing.T) {
		// a mock Sentinel-2 granule crossing 0 degrees longitude
		granuleBound := orb.Bound{
			Min: orb.Point{-0.2, 10.0},
			Max: orb.Point{0.2, 10.5},
		}

		expectedRegionCode := "P5228-Q3006"

		regionCode, suids, err := rhealpixorb.DeriveRegionCode64FromGeometry(el, granuleBound, 8)
		if err != nil {
			t.Fatalf("failed to derive region code: %v", err)
		}

		// must detect multi-facet spanning
		if len(suids) < 2 {
			t.Errorf("expected multi-root SUID across 0 deg meridian, got %d (%s)", len(suids), regionCode)
		}

		// exact ground-truth check
		if regionCode != expectedRegionCode {
			t.Errorf("region code mismatch: got %q, expected %q", regionCode, expectedRegionCode)
		}
	})

	t.Run("Polar Satellite Pass (North Cap)", func(t *testing.T) {
		// swath passing through high latitude / polar cap (Facet N)
		polarPass := orb.LineString{
			orb.Point{0.0, 85.0},
			orb.Point{45.0, 88.0},
			orb.Point{90.0, 85.0},
		}

		expectedRegionCode := "N44"

		regionCode, suids, err := rhealpixorb.DeriveRegionCode64FromGeometry(el, polarPass, 6)
		if err != nil {
			t.Fatalf("failed to derive polar region code: %v", err)
		}

		if len(suids) == 0 || regionCode == "" {
			t.Errorf("failed to produce region code for polar pass")
		}

		// exact ground-truth check
		if regionCode != expectedRegionCode {
			t.Errorf("region code mismatch: got %q, expected %q", regionCode, expectedRegionCode)
		}
	})
}
