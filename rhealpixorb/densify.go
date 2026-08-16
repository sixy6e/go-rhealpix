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

// densifyBoundEdges samples intermediate points along 2D coordinate edges,
// handling antimeridian wrapping cleanly.
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

		lon1, lat1 := p1.X(), p1.Y()
		lon2, lat2 := p2.X(), p2.Y()

		// shortest longitude delta across the ±180° antimeridian seam
		dLon := lon2 - lon1
		if dLon > 180.0 {
			dLon -= 360.0
		} else if dLon < -180.0 {
			dLon += 360.0
		}

		for s := 0; s < samplesPerEdge; s++ {
			t := float64(s) / float64(samplesPerEdge)

			x := lon1 + t*dLon
			y := lerp(lat1, lat2, t)

			// wrap back into [-180, 180] domain
			if x > 180.0 {
				x -= 360.0
			} else if x < -180.0 {
				x += 360.0
			}

			pts = append(pts, orb.Point{x, y})
		}
	}
	return pts
}

// densifyRingEdges interpolates N points along each segment of a polygon ring,
// handling antimeridian wrapping cleanly.
func densifyRingEdges(ring orb.Ring, pointsPerSegment int) []orb.Point {
	if len(ring) < 2 {
		return ring
	}
	if pointsPerSegment < 1 {
		pointsPerSegment = 1
	}

	densified := make([]orb.Point, 0, len(ring)*pointsPerSegment)

	for i := 0; i < len(ring)-1; i++ {
		p1 := ring[i]
		p2 := ring[i+1]

		lon1, lat1 := p1.X(), p1.Y()
		lon2, lat2 := p2.X(), p2.Y()

		// shortest longitude delta across the ±180° antimeridian seam
		dLon := lon2 - lon1
		if dLon > 180.0 {
			dLon -= 360.0
		} else if dLon < -180.0 {
			dLon += 360.0
		}

		for j := 0; j < pointsPerSegment; j++ {
			t := float64(j) / float64(pointsPerSegment)

			x := lon1 + t*dLon
			y := lerp(lat1, lat2, t)

			// wrap back into [-180, 180] domain
			if x > 180.0 {
				x -= 360.0
			} else if x < -180.0 {
				x += 360.0
			}

			densified = append(densified, orb.Point{x, y})
		}
	}
	return densified
}
