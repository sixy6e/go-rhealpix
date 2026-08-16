package rhealpixorb

import (
	"math"
	"testing"

	"github.com/paulmach/orb"
)

func TestDensifyBoundEdges_Antimeridian(t *testing.T) {
	// segment spanning 2 degrees across the Pacific seam: 179.0°E to -179.0°W
	corners := []orb.Point{
		{179.0, 0.0},
		{-179.0, 0.0},
	}

	// 2 samples per edge -> t = 0.0 and t = 0.5
	pts := densifyBoundEdges(corners, 2)

	if len(pts) != 4 {
		t.Fatalf("Expected 4 points, got %d", len(pts))
	}

	midpoint := pts[1] // t = 0.5

	// shortest path across seam places midpoint at ±180.0°
	// naive lerp places midpoint at 0.0° (Null Island / Greenwich)
	if math.Abs(math.Abs(midpoint.X())-180.0) > 1e-6 {
		t.Errorf("FAIL: Expected midpoint near ±180.0°, got longitude %f", midpoint.X())
	}
}

func TestDensifyRingEdges_Antimeridian(t *testing.T) {
	// closed ring segment crossing from 179.0°E to -179.0°W across the seam
	ring := orb.Ring{
		{179.0, 10.0},
		{-179.0, 10.0},
		{-179.0, 20.0},
		{179.0, 20.0},
		{179.0, 10.0}, // closed
	}

	// 2 samples per segment
	pts := densifyRingEdges(ring, 2)

	// first segment midpoint (t = 0.5) must land at ±180.0° longitude
	midpoint := pts[1]

	if math.Abs(math.Abs(midpoint.X())-180.0) > 1e-6 {
		t.Errorf("Expected ring segment midpoint near ±180.0°, got longitude %f", midpoint.X())
	}
}
