package rhealpixorb

import (
	"fmt"
	"math"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// FacetSubGeometry holds a WGS84 geometry clipped to a specific rHEALPix base facet.
type FacetSubGeometry struct {
	FacetID  uint8
	Geometry orb.Geometry
}

// UnwrapPolygon converts longitude jumps (>180°) into a continuous coordinate series
// across the anti-meridian without altering ring topology or orientation.
func UnwrapPolygon(poly orb.Polygon) orb.Polygon {
	if len(poly) == 0 {
		return poly
	}
	out := DeepClonePolygon(poly)

	for r := range out {
		ring := out[r]
		if len(ring) < 2 {
			continue
		}
		for i := 1; i < len(ring); i++ {
			diff := ring[i][0] - ring[i-1][0]
			if diff > 180.0 {
				ring[i][0] -= 360.0
			} else if diff < -180.0 {
				ring[i][0] += 360.0
			}
		}
	}
	return out
}

// DeepClonePolygon allocates fresh memory buffers for all rings and points.
func DeepClonePolygon(p orb.Polygon) orb.Polygon {
	if p == nil {
		return nil
	}
	newPoly := make(orb.Polygon, len(p))
	for i, ring := range p {
		newRing := make(orb.Ring, len(ring))
		copy(newRing, ring)
		newPoly[i] = newRing
	}
	return newPoly
}

// normaliseLonDeg wraps unwrapped longitudes back to [-180, 180) for domain checks.
func normaliseLonDeg(lon float64) float64 {
	lon = math.Mod(lon+180.0, 360.0)
	if lon < 0 {
		lon += 360.0
	}
	return lon - 180.0
}

// facetIntersectsBound performs a geodetic bounding box check against
// the 6 rHEALPix base facet domains using the ellipsoid's exact authalic cap latitude.
func facetIntersectsBound(facetID uint8, geomBound orb.Bound, capLat float64) bool {
	minY := geomBound.Min.Y()
	maxY := geomBound.Max.Y()
	minX := normaliseLonDeg(geomBound.Min.X())
	maxX := normaliseLonDeg(geomBound.Max.X())

	// handle bounding box longitude wrapping span
	if geomBound.Max.X()-geomBound.Min.X() >= 360.0 {
		minX = -180.0
		maxX = 180.0
	}

	switch facetID {
	case 0: // North Polar Cap ('N' / 0)
		return maxY >= capLat

	case 1: // Equatorial Facet O / 1 [-180°, -90°), [-pi, -pi/2)
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX < -90.0 && maxX >= -180.0

	case 2: // Equatorial Facet P / 2 [-90°, 0°), [-pi/2, 0)
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX < 0.0 && maxX >= -90.0

	case 3: // Equatorial Facet Q / 3 [0°, 90°), [0, pi/2)
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX < 90.0 && maxX >= 0.0

	case 4: // Equatorial Facet R / 4 [90°, 180°], [pi/2, pi]
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX <= 180.0 && maxX >= 90.0

	case 5: // South Polar Cap ('S' / 5)
		return minY <= -capLat

	default:
		return false
	}
}

// ensurePolygonCCW deep-clones a polygon, forces exterior rings to CCW (1),
// interior rings to CW (-1), and strips degenerate geometries.
func ensurePolygonCCW(poly orb.Polygon) orb.Polygon {
	if len(poly) == 0 {
		return nil
	}
	out := make(orb.Polygon, 0, len(poly))

	for i, ring := range poly {
		if len(ring) < 4 { // minimum closed ring needs 4 points (3 unique + closure)
			continue
		}
		newRing := make(orb.Ring, len(ring))
		copy(newRing, ring)

		// ensure closed ring topology
		if newRing[0] != newRing[len(newRing)-1] {
			newRing = append(newRing, newRing[0])
		}

		// ring 0: Exterior (MUST BE CCW)
		// rings 1..N: Interior Holes (MUST BE CW)
		if i == 0 {
			if newRing.Orientation() != orb.CCW {
				newRing.Reverse()
			}
		} else {
			if newRing.Orientation() != orb.CW {
				newRing.Reverse()
			}
		}

		out = append(out, newRing)
	}

	if len(out) == 0 {
		return nil
	}
	return out
}

// splitPolygonAtLatitude performs exact edge-intersection slicing along an authalic latitude parallel,
// maintaining unwrapped longitudes without generating cartesian bounding box slivers.
func splitPolygonAtLatitude(poly orb.Polygon, splitLat float64) (northPoly, southPoly orb.Polygon) {
	for _, ring := range poly {
		var nRing, sRing orb.Ring

		for i := 0; i < len(ring)-1; i++ {
			p1 := ring[i]
			p2 := ring[i+1]

			lat1, lat2 := p1[1], p2[1]

			// detect authalic boundary crossing edge
			if (lat1 > splitLat && lat2 < splitLat) || (lat1 < splitLat && lat2 > splitLat) {
				t := (splitLat - lat1) / (lat2 - lat1)
				lonInterp := p1[0] + t*(p2[0]-p1[0]) // interpolate along unwrapped longitude
				intersectionPt := orb.Point{lonInterp, splitLat}

				if lat1 >= splitLat {
					nRing = append(nRing, p1, intersectionPt)
					sRing = append(sRing, intersectionPt)
				} else {
					sRing = append(sRing, p1, intersectionPt)
					nRing = append(nRing, intersectionPt)
				}
			} else {
				if lat1 >= splitLat {
					nRing = append(nRing, p1)
				} else {
					sRing = append(sRing, p1)
				}
			}
		}

		// close rings & filter degenerate slices
		if len(nRing) >= 3 {
			if nRing[0] != nRing[len(nRing)-1] {
				nRing = append(nRing, nRing[0])
			}
			if len(nRing) >= 4 {
				northPoly = append(northPoly, nRing)
			}
		}
		if len(sRing) >= 3 {
			if sRing[0] != sRing[len(sRing)-1] {
				sRing = append(sRing, sRing[0])
			}
			if len(sRing) >= 4 {
				southPoly = append(southPoly, sRing)
			}
		}
	}

	return ensurePolygonCCW(northPoly), ensurePolygonCCW(southPoly)
}

// ClipGeometryToFacets partitions a WGS84 geometry across the 6 rHEALPix base facets
// using authalic edge splitting and longitude unwrapping.
func ClipGeometryToFacets(el *rhealpix.Ellipsoid, geom orb.Geometry) ([]FacetSubGeometry, error) {
	if geom == nil {
		return nil, fmt.Errorf("nil input geometry")
	}

	// convert input geometry into orb.Polygon or orb.MultiPolygon
	var polys []orb.Polygon
	switch g := geom.(type) {
	case orb.Polygon:
		polys = append(polys, UnwrapPolygon(g))
	case orb.MultiPolygon:
		for _, p := range g {
			polys = append(polys, UnwrapPolygon(p))
		}
	default:
		return nil, fmt.Errorf("unsupported geometry type %T", geom)
	}

	// compute exact authalic cap latitude for configured ellipsoid
	phi0 := math.Asin(2.0 / 3.0)
	capLatRad := el.AuthLatInverse(phi0)
	capLatDeg := capLatRad * (180.0 / math.Pi)

	var results []FacetSubGeometry

	for _, poly := range polys {
		poly = ensurePolygonCCW(poly)
		if len(poly) == 0 {
			continue
		}

		geomBound := poly.Bound()

		// Case A: geometry straddles the South Cap authalic seam (e.g. 55GDP in Tasmania)
		if geomBound.Min.Y() < -capLatDeg && geomBound.Max.Y() > -capLatDeg {
			northPiece, southPiece := splitPolygonAtLatitude(poly, -capLatDeg)

			if len(northPiece) > 0 {
				facetID := getEquatorialFacetForPoly(northPiece)
				results = append(results, FacetSubGeometry{FacetID: facetID, Geometry: northPiece})
			}
			if len(southPiece) > 0 {
				results = append(results, FacetSubGeometry{FacetID: 5, Geometry: southPiece})
			}
			continue
		}

		// Case B: geometry straddles the North Cap authalic seam
		if geomBound.Min.Y() < capLatDeg && geomBound.Max.Y() > capLatDeg {
			northPiece, southPiece := splitPolygonAtLatitude(poly, capLatDeg)

			if len(northPiece) > 0 {
				results = append(results, FacetSubGeometry{FacetID: 0, Geometry: northPiece})
			}
			if len(southPiece) > 0 {
				facetID := getEquatorialFacetForPoly(southPiece)
				results = append(results, FacetSubGeometry{FacetID: facetID, Geometry: southPiece})
			}
			continue
		}

		// Case C: Single-facet or non-cap-straddling geometries
		if geomBound.Max.Y() <= capLatDeg && geomBound.Min.Y() >= -capLatDeg {
			// Strictly inside Equatorial Belt -> resolve correct facet (1=O, 2=P, 3=Q, 4=R)
			facetID := getEquatorialFacetForPoly(poly)
			results = append(results, FacetSubGeometry{
				FacetID:  facetID,
				Geometry: poly,
			})
			continue
		}

		// Case D: (fallback) for polar cap single-facet geometries
		for facetID := uint8(0); facetID < 6; facetID++ {
			if facetIntersectsBound(facetID, geomBound, capLatDeg) {
				results = append(results, FacetSubGeometry{
					FacetID:  facetID,
					Geometry: poly,
				})
				break
			}
		}
	}

	return results, nil
}

// getEquatorialFacetForPoly resolves which equatorial facet (1=O, 2=P, 3=Q, 4=R)
// a polygon belongs to based on its normalised centroid longitude.
func getEquatorialFacetForPoly(poly orb.Polygon) uint8 {
	if len(poly) == 0 || len(poly[0]) == 0 {
		return 1
	}
	centerLon := normaliseLonDeg(poly[0][0][0])

	switch {
	case centerLon < -90.0:
		return 1 // Facet O
	case centerLon < 0.0:
		return 2 // Facet P
	case centerLon < 90.0:
		return 3 // Facet Q
	default:
		return 4 // Facet R
	}
}
