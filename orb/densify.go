package rhealpixorb

import (
	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// lerp performs linear interpolation between a and b for parameter t in [0.0, 1.0].
func lerp(a, b, t float64) float64 {
	return a + t*(b-a)
}

// buildDensifiedRing interpolates planar points along valid corner edges in local facet space
// and projects them to Lon/Lat (EPSG:4326).
func buildDensifiedRing(el *rhealpix.Ellipsoid, facet uint8, corners [][2]float64, segmentCount int) orb.Ring {
	if segmentCount < 1 {
		segmentCount = 1
	}

	numCorners := len(corners)
	if numCorners == 0 {
		return orb.Ring{}
	}

	totalPoints := numCorners * segmentCount
	ring := make(orb.Ring, 0, totalPoints+1)

	// walk around the polygon edges in order
	for i := 0; i < numCorners; i++ {
		pStart := corners[i]
		pEnd := corners[(i+1)%numCorners]

		for s := 0; s < segmentCount; s++ {
			t := float64(s) / float64(segmentCount)
			xLocal := lerp(pStart[0], pEnd[0], t)
			yLocal := lerp(pStart[1], pEnd[1], t)

			lon, lat := facetLocalToLonLat(el, facet, xLocal, yLocal)
			ring = append(ring, orb.Point{lon, lat})
		}
	}

	// close ring
	if len(ring) > 0 {
		ring = append(ring, ring[0])
	}

	return ring
}

// densifyBoundEdges samples intermediate points along 2D coordinate edges.
func densifyBoundEdges(corners []orb.Point, samplesPerEdge int) []orb.Point {
	if samplesPerEdge < 1 {
		samplesPerEdge = 1
	}

	n := len(corners)
	if n == 0 {
		return nil
	}

	pts := make([]orb.Point, 0, n*samplesPerEdge)

	for i := 0; i < n; i++ {
		p1 := corners[i]
		p2 := corners[(i+1)%n]

		for s := 0; s < samplesPerEdge; s++ {
			t := float64(s) / float64(samplesPerEdge)
			x := lerp(p1.X(), p2.X(), t)
			y := lerp(p1.Y(), p2.Y(), t)
			pts = append(pts, orb.Point{x, y})
		}
	}
	return pts
}
