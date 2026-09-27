package rhealpixorb

import (
	"fmt"
	"math"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/clip"
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

func testNS(d *rhealpix.RhpDGGS, latRad float64) {
	beta0 := d.GetEllipsoid().Beta0
	aulat := d.GetEllipsoid().AuthLat(latRad)
	deg := aulat * rhealpix.RadToDeg
	if math.Abs(aulat) <= beta0 {
		fmt.Printf("Point is in equatorial region: lat: %v, authlat: %v, beta0: %v\n", latRad*rhealpix.RadToDeg, deg, beta0*rhealpix.RadToDeg)
	} else {
		fmt.Printf("Point is in polar region: lat: %v, authlat: %v, beta0: %v\n", latRad*rhealpix.RadToDeg, deg, beta0*rhealpix.RadToDeg)
	}
}

// ClipGeometryToFacets slices an input geographic geometry across the rHEALPix facet boundaries
// (authalic cap parallels and equatorial meridians) and projects each piece into planar space.
func ClipGeometryToFacets(d *rhealpix.RhpDGGS, geom orb.Geometry) ([]FacetSubGeometry, error) {
	if geom == nil {
		return nil, fmt.Errorf("nil input geometry")
	}

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

	capLatDeg := d.GetEllipsoid().Phi0 * rhealpix.RadToDeg
	var results []FacetSubGeometry

	for _, poly := range polys {
		poly = ensurePolygonCCW(poly)
		if len(poly) == 0 {
			continue
		}

		geomBound := poly.Bound()

		// --- Case A: Straddles North Cap Seam (+Phi0) ---
		if geomBound.Min.Y() < capLatDeg && geomBound.Max.Y() > capLatDeg {
			northGeo, southGeo := splitPolygonAtLatitude(poly, capLatDeg)

			// process North Polar Cap piece (Facet 0)
			if len(northGeo) > 0 {
				planarNorth := ForwardProjectPolygonForcedRegion(d, northGeo, rhealpix.RegionNorthPolar)
				if planarNorth != nil {
					results = append(results, FacetSubGeometry{FacetID: 0, Geometry: planarNorth})
				}
			}

			// process Equatorial piece (Facets 1..4)
			if len(southGeo) > 0 {
				processEquatorialPiece(d, southGeo, &results)
			}
			continue
		}

		// --- Case B: Straddles South Cap Seam (-Phi0) ---
		if geomBound.Min.Y() < -capLatDeg && geomBound.Max.Y() > -capLatDeg {
			northGeo, southGeo := splitPolygonAtLatitude(poly, -capLatDeg)

			if len(southGeo) > 0 {
				planarSouth := ForwardProjectPolygonForcedRegion(d, southGeo, rhealpix.RegionSouthPolar)
				if planarSouth != nil {
					results = append(results, FacetSubGeometry{FacetID: 5, Geometry: planarSouth})
				}
			}

			if len(northGeo) > 0 {
				processEquatorialPiece(d, northGeo, &results)
			}
			continue
		}

		// --- Case C: Completely inside North Polar Cap (Facet 0) ---
		if geomBound.Min.Y() >= capLatDeg {
			planarPoly := ForwardProjectPolygonForcedRegion(d, poly, rhealpix.RegionNorthPolar)
			if planarPoly != nil {
				results = append(results, FacetSubGeometry{FacetID: 0, Geometry: planarPoly})
			}
			continue
		}

		// --- Case D: Completely inside South Polar Cap (Facet 5) ---
		if geomBound.Max.Y() <= -capLatDeg {
			planarPoly := ForwardProjectPolygonForcedRegion(d, poly, rhealpix.RegionSouthPolar)
			if planarPoly != nil {
				results = append(results, FacetSubGeometry{FacetID: 5, Geometry: planarPoly})
			}
			continue
		}

		// --- Case E: Entirely inside Equatorial Belt (Facets 1..4) ---
		processEquatorialPiece(d, poly, &results)
	}

	return results, nil
}

// processEquatorialPiece partitions an equatorial polygon across longitude meridians and projects each piece.
func processEquatorialPiece(d *rhealpix.RhpDGGS, poly orb.Polygon, results *[]FacetSubGeometry) {
	eqPieces := partitionEquatorialPolygon(d, poly)
	for _, piece := range eqPieces {
		eqPoly, ok := piece.Geometry.(orb.Polygon)
		if !ok || len(eqPoly) == 0 {
			continue
		}
		planarEq := ForwardProjectPolygonForcedRegion(d, eqPoly, rhealpix.RegionEquatorial)
		if planarEq != nil {
			*results = append(*results, FacetSubGeometry{FacetID: piece.FacetID, Geometry: planarEq})
		}
	}
}

// SplitAtSectorMeridians splits segments crossing polar sector meridians
// (-180°, -90°, 0°, +90°, +180°) when the segment lies within the polar caps.
func SplitAtSectorMeridians(d *rhealpix.RhpDGGS, ring orb.Ring) orb.Ring {
	if len(ring) < 2 {
		return ring
	}

	phi0 := d.GetEllipsoid().Phi0 * rhealpix.RadToDeg
	sectorMeridians := []float64{-180.0, -90.0, 0.0, 90.0, 180.0}

	splitRing := make(orb.Ring, 0, len(ring)*2)

	for i := 0; i < len(ring)-1; i++ {
		p1 := ring[i]
		p2 := ring[i+1]

		lon1, lat1 := p1.X(), p1.Y()
		lon2, lat2 := p2.X(), p2.Y()

		splitRing = append(splitRing, p1)

		// only check sector boundary crossings if in the North or South polar caps
		inPolarCap := (lat1 > phi0 && lat2 > phi0) || (lat1 < -phi0 && lat2 < -phi0)
		if !inPolarCap {
			continue
		}

		dLon := lon2 - lon1
		if dLon > 180.0 {
			dLon -= 360.0
		} else if dLon < -180.0 {
			dLon += 360.0
		}

		if math.Abs(dLon) < 1e-9 {
			continue
		}

		dLat := lat2 - lat1

		// find any sector meridian crossed by this segment
		for _, mer := range sectorMeridians {
			// check if mer falls strictly between lon1 and lon2 (accounting for antimeridian wrap)
			crossed := false
			var t float64

			// standard span
			if lon1 < mer && lon2 > mer {
				t = (mer - lon1) / dLon
				crossed = true
			} else if lon1 > mer && lon2 < mer {
				t = (mer - lon1) / dLon
				crossed = true
			}

			if crossed && t > 1e-6 && t < (1.0-1e-6) {
				latIntersect := lat1 + t*dLat
				intersectPt := orb.Point{mer, latIntersect}
				splitRing = append(splitRing, intersectPt)
			}
		}
	}

	// always preserve the final closing point
	if len(ring) > 0 {
		splitRing = append(splitRing, ring[len(ring)-1])
	}

	return splitRing
}

func ClipPlanarPolygonToBox(poly orb.Polygon, bounds rhealpix.BaseFacetBounds) orb.Polygon {
	if len(poly) == 0 {
		return nil
	}

	// create an orb.Bound for the facet
	b := orb.Bound{
		Min: orb.Point{bounds.MinX(), bounds.MinY()},
		Max: orb.Point{bounds.MaxX(), bounds.MaxY()},
	}

	// clip the polygon using orb/clip
	clipped := clip.Polygon(b, poly)
	if len(clipped) == 0 {
		return nil
	}

	return clipped
}

// clipRingAgainstLine clips an orb.Ring against an axis-aligned boundary line (X or Y).
func clipRingAgainstLine(ring orb.Ring, edgeVal float64, isX, isGreater bool) orb.Ring {
	if len(ring) == 0 {
		return nil
	}

	var output orb.Ring

	isInside := func(pt orb.Point) bool {
		val := pt[0]
		if !isX {
			val = pt[1]
		}
		if isGreater {
			return val >= edgeVal
		}
		return val <= edgeVal
	}

	intersect := func(p1, p2 orb.Point) orb.Point {
		if isX {
			t := (edgeVal - p1[0]) / (p2[0] - p1[0])
			y := p1[1] + t*(p2[1]-p1[1])
			return orb.Point{edgeVal, y}
		}
		t := (edgeVal - p1[1]) / (p2[1] - p1[1])
		x := p1[0] + t*(p2[0]-p1[0])
		return orb.Point{x, edgeVal}
	}

	// remove duplicated closing point during processing loop
	pts := ring
	if len(pts) > 1 && pts[0] == pts[len(pts)-1] {
		pts = pts[:len(pts)-1]
	}

	s := pts[len(pts)-1]

	for _, e := range pts {
		eInside := isInside(e)
		sInside := isInside(s)

		if eInside {
			if !sInside {
				output = append(output, intersect(s, e))
			}
			output = append(output, e)
		} else if sInside {
			output = append(output, intersect(s, e))
		}
		s = e
	}

	return output
}

// getEquatorialFacetForPoly resolves which equatorial facet (1=O, 2=P, 3=Q, 4=R)
// a polygon belongs to based on its normalised centroid longitude.
func getEquatorialFacetForPoly(poly orb.Polygon) uint8 {
	if len(poly) == 0 || len(poly[0]) == 0 {
		return 1
	}
	centreLon := normaliseLonDeg(poly[0][0][0])

	switch {
	case centreLon < -90.0:
		return 1 // Facet O
	case centreLon < 0.0:
		return 2 // Facet P
	case centreLon < 90.0:
		return 3 // Facet Q
	default:
		return 4 // Facet R
	}
}

// partitionEquatorialPolygon splits an equatorial polygon across longitude meridian meridian boundaries (O, P, Q, R)
// (-90°, 0°, +90°) and assigns each resulting piece to its correct equatorial facet (1=O, 2=P, 3=Q, 4=R).
// It also accounts for for d.Ellipsoid().Lon0 rotation and assigns each resulting piece to its correct facet.
func partitionEquatorialPolygon(d *rhealpix.RhpDGGS, poly orb.Polygon) []FacetSubGeometry {
	if len(poly) == 0 {
		return nil
	}

	lon0Deg := d.GetEllipsoid().Lon0 * rhealpix.RadToDeg

	// base meridian seams relative to Lon0:
	// Facet O (-180..-90), Facet P (-90..0), Facet Q (0..90), Facet R (90..180)
	meridians := []float64{
		normaliseLonDeg(lon0Deg - 90.0),
		normaliseLonDeg(lon0Deg),
		normaliseLonDeg(lon0Deg + 90.0),
	}

	currentPolys := []orb.Polygon{poly}

	// slice sequentially across the 3 equatorial meridian seams (-90°, 0°, +90°)
	for _, m := range meridians {
		var nextPolys []orb.Polygon
		for _, p := range currentPolys {
			left, right := splitPolygonAtLongitude(p, m)
			if len(left) > 0 {
				nextPolys = append(nextPolys, left)
			}
			if len(right) > 0 {
				nextPolys = append(nextPolys, right)
			}
		}
		currentPolys = nextPolys
	}

	var results []FacetSubGeometry
	for _, p := range currentPolys {
		if len(p) == 0 {
			continue
		}
		facetID := resolveEquatorialFacetFromCentroid(d, p)
		results = append(results, FacetSubGeometry{FacetID: facetID, Geometry: p})
	}

	return results
}

// resolveEquatorialFacetFromCentroid determines the equatorial facet (1=O, 2=P, 3=Q, 4=R)
// for a sliced polygon piece by evaluating its centroid relative to Lon0.
func resolveEquatorialFacetFromCentroid(d *rhealpix.RhpDGGS, poly orb.Polygon) uint8 {
	lon0Deg := d.GetEllipsoid().Lon0 * rhealpix.RadToDeg
	centroid := poly.Bound().Center()

	// normalise longitude relative to Lon0 into [-180, 180)
	normLon := normaliseLonDeg(centroid.X() - lon0Deg)

	switch {
	case normLon < -90.0:
		return 1 // Facet O
	case normLon < 0.0:
		return 2 // Facet P
	case normLon < 90.0:
		return 3 // Facet Q
	default:
		return 4 // Facet R
	}
}

// splitPolygonAtLongitude cuts an orb.Polygon along a target longitude meridian.
// Returns left (lon <= targetLon) and right (lon >= targetLon) sub-polygons.
func splitPolygonAtLongitude(poly orb.Polygon, targetLon float64) (orb.Polygon, orb.Polygon) {
	if len(poly) == 0 {
		return nil, nil
	}

	var leftRings, rightRings []orb.Ring

	for _, ring := range poly {
		if len(ring) < 3 {
			continue
		}

		leftRing, rightRing := splitRingAtLongitude(ring, targetLon)
		if len(leftRing) >= 4 {
			leftRings = append(leftRings, leftRing)
		}
		if len(rightRing) >= 4 {
			rightRings = append(rightRings, rightRing)
		}
	}

	var leftPoly, rightPoly orb.Polygon
	if len(leftRings) > 0 {
		leftPoly = ensurePolygonCCW(orb.Polygon(leftRings))
	}
	if len(rightRings) > 0 {
		rightPoly = ensurePolygonCCW(orb.Polygon(rightRings))
	}

	return leftPoly, rightPoly
}

// splitRingAtLongitude splits an individual orb.Ring along a vertical meridian line.
func splitRingAtLongitude(ring orb.Ring, targetLon float64) (orb.Ring, orb.Ring) {
	var leftRing, rightRing orb.Ring

	numPts := len(ring)
	for i := 0; i < numPts-1; i++ {
		p1 := ring[i]
		p2 := ring[i+1]

		lon1, lat1 := p1[0], p1[1]
		lon2, lat2 := p2[0], p2[1]

		// add current point to appropriate side(s)
		if lon1 <= targetLon {
			leftRing = append(leftRing, p1)
		}
		if lon1 >= targetLon {
			rightRing = append(rightRing, p1)
		}

		// check if segment crosses the target longitude line
		if (lon1 < targetLon && lon2 > targetLon) || (lon1 > targetLon && lon2 < targetLon) {
			// linear interpolation for intersection latitude
			t := (targetLon - lon1) / (lon2 - lon1)
			intersectLat := lat1 + t*(lat2-lat1)
			intersectPt := orb.Point{targetLon, intersectLat}

			leftRing = append(leftRing, intersectPt)
			rightRing = append(rightRing, intersectPt)
		}
	}

	// close left ring if needed
	if len(leftRing) > 0 && (leftRing[0][0] != leftRing[len(leftRing)-1][0] || leftRing[0][1] != leftRing[len(leftRing)-1][1]) {
		leftRing = append(leftRing, leftRing[0])
	}

	// close right ring if needed
	if len(rightRing) > 0 && (rightRing[0][0] != rightRing[len(rightRing)-1][0] || rightRing[0][1] != rightRing[len(rightRing)-1][1]) {
		rightRing = append(rightRing, rightRing[0])
	}

	return leftRing, rightRing
}

// func ForwardProjectPointForFacet(d *rhealpix.RhpDGGS, facetID uint8, pt orb.Point) orb.Point {
func ForwardProjectPointForFacet(d *rhealpix.RhpDGGS, facetID uint8, pt orb.Point) (xRhp, yRhp float64) {
	lonRad := pt.X() * rhealpix.DegToRad
	latRad := pt.Y() * rhealpix.DegToRad

	capLatRad := d.GetEllipsoid().Phi0

	// clamp ONLY vertices sitting right on the seam boundary (within ~1cm tolerance)
	// DO NOT clamp interior vertices (latRad > capLatRad)!
	if facetID == 0 {
		if math.Abs(latRad-capLatRad) < rhealpix.EPS_SEAM_RAD || latRad < capLatRad {
			latRad = capLatRad
		}
	} else if facetID == 5 {
		if math.Abs(latRad-(-capLatRad)) < rhealpix.EPS_SEAM_RAD || latRad > -capLatRad {
			latRad = -capLatRad
		}
	} else if facetID >= 1 && facetID <= 4 {
		if math.Abs(latRad-capLatRad) < rhealpix.EPS_SEAM_RAD {
			latRad = capLatRad
		} else if math.Abs(latRad-(-capLatRad)) < rhealpix.EPS_SEAM_RAD {
			latRad = -capLatRad
		}
	}

	// ForwardProjectRadians handles NormaliseLongitude and AuthLat internally!
	xRhp, yRhp, _ = d.ForwardProjectRadians(lonRad, latRad)

	// hard-lock Y to exact facet boundary baseline ONLY when latitude is on the seam
	if latRad == capLatRad {
		if facetID == 0 {
			bounds0, _ := d.FacetBounds(0)
			yRhp = bounds0.MinY() // hard-lock Y to 0.7853981633974483
		} else if facetID >= 1 && facetID <= 4 {
			boundsEq, _ := d.FacetBounds(facetID)
			yRhp = boundsEq.MaxY()
		}
	} else if latRad == -capLatRad {
		if facetID == 5 {
			bounds5, _ := d.FacetBounds(5)
			yRhp = bounds5.MaxY()
		} else if facetID >= 1 && facetID <= 4 {
			boundsEq, _ := d.FacetBounds(facetID)
			yRhp = boundsEq.MinY()
		}
	}

	// return orb.Point{xRhp, yRhp}
	return xRhp, yRhp
}
