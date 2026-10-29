package rhealpixorb_test

import (
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
	"github.com/sixy6e/go-rhealpix/rhealpixorb"
)

func TestCell64ToPolygon(t *testing.T) {
	el := rhealpix.NewWGS84()
	d := rhealpix.NewRhpDGGS(
		rhealpix.WithEllipsoid(el),
		rhealpix.WithNSide(3),
		rhealpix.WithNorthSquare(0),
		rhealpix.WithSouthSquare(0),
	)

	cell, err := rhealpix.ParseCellID64("Q012")
	if err != nil {
		t.Fatalf("unexpected parse error: %v", err)
	}

	// standard 4-corner polygon (segmentCount = 1)
	poly, err := rhealpixorb.Cell64ToPolygon(d, cell, 1)
	if err != nil {
		t.Fatalf("unexpected polygon conversion error: %v", err)
	}

	// 4 corners + 1 closure point = 5 points
	if len(poly) == 0 || len(poly[0]) != 5 {
		t.Errorf("expected closed 5-point ring for segmentCount=1, got %d points", len(poly[0]))
	}

	// densified polygon (segmentCount = 4) -> 4 corners * 4 segments + 1 closure point = 17 points
	polyDensified, err := rhealpixorb.Cell64ToPolygon(d, cell, 4)
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
	d := rhealpix.NewRhpDGGS(
		rhealpix.WithEllipsoid(el),
		rhealpix.WithNSide(3),
		rhealpix.WithNorthSquare(0),
		rhealpix.WithSouthSquare(0),
	)

	// parse a cell on North Polar Facet N (Facet 0)
	cell, err := rhealpix.ParseCellID64("N0")
	if err != nil {
		t.Fatalf("unexpected parse error for polar cell: %v", err)
	}

	poly, err := rhealpixorb.Cell64ToPolygon(d, cell, 1)
	if err != nil {
		t.Fatalf("unexpected polar polygon conversion error: %v", err)
	}

	if len(poly) == 0 || len(poly[0]) < 4 {
		t.Errorf("expected valid closed polar polygon ring, got %v", poly)
	}
}
