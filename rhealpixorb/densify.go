package rhealpixorb

import (
	"math"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// lerp performs linear interpolation between a and b for parameter t in [0.0, 1.0].
func lerp(a, b, t float64) float64 {
	return a + t*(b-a)
}

// buildDensifiedRing interpolates planar points along valid corner edges in local facet space
// and projects them to Lon/Lat (EPSG:4326).
func buildDensifiedRing(d *rhealpix.RhpDGGS, facet uint8, corners [][2]float64, segmentCount int) orb.Ring {
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

			lon, lat := facetLocalToLonLat(d, facet, xLocal, yLocal)
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

// DensifyRingByStep samples intermediate points along a geographic ring's segments
// ensuring no segment exceeds maxStepDeg (in degrees). It preserves antimeridian wrapping.
// Similar to densifyRingEdges without the bias of n-points per segment.
// e.g. short segment and long segment both using the same number of points
// one being heavily sampled in a short space, and one getting under sampled
// over a long space.
func DensifyRingByStep(ring orb.Ring, maxStepDeg float64) orb.Ring {
	if len(ring) < 2 {
		return ring
	}
	if maxStepDeg <= 0 {
		maxStepDeg = 0.05 // default ~5km step
	}

	// estimate capacity
	densified := make(orb.Ring, 0, len(ring)*2)

	for i := 0; i < len(ring)-1; i++ {
		p1 := ring[i]
		p2 := ring[i+1]

		lon1, lat1 := p1.X(), p1.Y()
		lon2, lat2 := p2.X(), p2.Y()

		// shortest longitude delta across ±180° seam
		dLon := lon2 - lon1
		if dLon > 180.0 {
			dLon -= 360.0
		} else if dLon < -180.0 {
			dLon += 360.0
		}
		dLat := lat2 - lat1

		// calculate total angular Euclidean step distance in degrees
		dist := math.Hypot(dLon, dLat)

		// calculate required segment count for this specific edge
		steps := 1
		if dist > maxStepDeg {
			steps = int(math.Ceil(dist / maxStepDeg))
		}

		for j := 0; j < steps; j++ {
			t := float64(j) / float64(steps)

			x := lon1 + t*dLon
			y := lerp(lat1, lat2, t)

			// wrap longitude back into [-180, 180] domain
			if x > 180.0 {
				x -= 360.0
			} else if x < -180.0 {
				x += 360.0
			}

			densified = append(densified, orb.Point{x, y})
		}
	}

	// always append the final closing point of the original ring
	if len(ring) > 0 {
		densified = append(densified, ring[len(ring)-1])
	}

	return densified
}

func DensifyAndSplitAtBeta0(d *rhealpix.RhpDGGS, ring orb.Ring, maxStepDeg float64) orb.Ring {
	if len(ring) < 2 {
		return ring
	}
	if maxStepDeg <= 0 {
		maxStepDeg = 0.02 // ~2km step size in degrees
	}

	phi0 := d.GetEllipsoid().Phi0 * rhealpix.RadToDeg // Authalic seam latitude in degrees

	densified := make(orb.Ring, 0, len(ring)*2)

	for i := 0; i < len(ring)-1; i++ {
		p1 := ring[i]
		p2 := ring[i+1]

		lon1, lat1 := p1.X(), p1.Y()
		lon2, lat2 := p2.X(), p2.Y()

		// shortest longitude delta across ±180° antimeridian seam
		dLon := lon2 - lon1
		if dLon > 180.0 {
			dLon -= 360.0
		} else if dLon < -180.0 {
			dLon += 360.0
		}
		dLat := lat2 - lat1

		// check if the segment crosses the North (+Phi0) or South (-Phi0) seam
		crossesNorthSeam := (lat1 < phi0 && lat2 > phi0) || (lat1 > phi0 && lat2 < phi0)
		crossesSouthSeam := (lat1 < -phi0 && lat2 > -phi0) || (lat1 > -phi0 && lat2 < -phi0)

		if crossesNorthSeam || crossesSouthSeam {
			targetPhi := phi0
			if crossesSouthSeam {
				targetPhi = -phi0
			}

			// linear fraction t where segment intersects targetPhi
			tSeam := (targetPhi - lat1) / dLat
			lonSeam := lon1 + tSeam*dLon

			// wrap seam longitude into [-180, 180]
			if lonSeam > 180.0 {
				lonSeam -= 360.0
			} else if lonSeam < -180.0 {
				lonSeam += 360.0
			}

			seamPt := orb.Point{lonSeam, targetPhi}

			// Sub-segment 1: p1 -> seamPt (densifySegment includes p1, excludes seamPt)
			densified = append(densified, densifySegment(p1, seamPt, maxStepDeg)...)
			// Sub-segment 2: seamPt -> p2 (densifySegment includes seamPt, excludes p2)
			densified = append(densified, densifySegment(seamPt, p2, maxStepDeg)...)
			continue
		}

		// standard segment densification (includes p1, excludes p2)
		densified = append(densified, densifySegment(p1, p2, maxStepDeg)...)
	}

	// append the final closing point
	if len(ring) > 0 {
		densified = append(densified, ring[len(ring)-1])
	}

	return densified
}

// densifySegment interpolates points starting at p1 up to (but not including) p2.
// Points p1 and p2 must ALREADY be unwrapped relative to each other.
func densifySegment(p1, p2 orb.Point, maxStepDeg float64) []orb.Point {
	dLon := p2.X() - p1.X()
	dLat := p2.Y() - p1.Y()

	distDeg := math.Hypot(dLon, dLat)
	if distDeg < 1e-12 {
		// ignore zero-length segments to prevent duplicate points
		return nil
	}

	steps := 1
	if distDeg > maxStepDeg {
		steps = int(math.Ceil(distDeg / maxStepDeg))
	}

	pts := make([]orb.Point, 0, steps)
	for s := 0; s < steps; s++ {
		t := float64(s) / float64(steps)
		x := p1.X() + t*dLon
		y := p1.Y() + t*dLat

		// keep x in continuous unwrapped space! Do not wrap back into [-180, 180].
		pts = append(pts, orb.Point{x, y})
	}
	return pts
}

// redo using new projection workflow
// densifyRing inserts intermediate vertices on segments longer than maxStepDeg
func densifyRing(ring orb.Ring, maxStepDeg float64) orb.Ring {
	if len(ring) < 2 {
		return ring
	}

	densified := make(orb.Ring, 0, len(ring)*2)
	for i := 0; i < len(ring)-1; i++ {
		p1, p2 := ring[i], ring[i+1]
		densified = append(densified, densifySegment(p1, p2, maxStepDeg)...)
	}
	if len(ring) > 0 {
		densified = append(densified, ring[len(ring)-1])
	}
	return densified
}

// densifyPolygon samples intermediate points along every ring of an orb.Polygon
// ensuring no segment exceeds maxStepDeg (in degrees).
func densifyPolygon(poly orb.Polygon, maxStepDeg float64) orb.Polygon {
	if len(poly) == 0 {
		return poly
	}

	densified := make(orb.Polygon, len(poly))
	for i, ring := range poly {
		densified[i] = densifyRing(ring, maxStepDeg)
	}
	return densified
}

// DensifyWGS84Ring inserts intermediate sample points along WGS84 ring edges
// so non-linear polar cap projections preserve curved boundary shapes and quadrant seams.
func DensifyWGS84Ring(ring orb.Ring, pointsPerSegment int) orb.Ring {
	if len(ring) < 2 || pointsPerSegment <= 1 {
		return ring
	}

	densified := make(orb.Ring, 0, (len(ring)-1)*pointsPerSegment+1)

	for i := 0; i < len(ring)-1; i++ {
		p1 := ring[i]
		p2 := ring[i+1]

		for s := 0; s < pointsPerSegment; s++ {
			t := float64(s) / float64(pointsPerSegment)
			x := p1[0] + t*(p2[0]-p1[0])
			y := p1[1] + t*(p2[1]-p1[1])
			densified = append(densified, orb.Point{x, y})
		}
	}

	densified = append(densified, ring[len(ring)-1])
	return densified
}
