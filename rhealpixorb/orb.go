package rhealpixorb

import (
	"fmt"
	"math"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// CleanCellSeamVertex aligns vertex longitude signs to match the cell's
// topological quadrant (Eastern vs Western) along the 180° antimeridian.
func CleanCellSeamVertex(pt [2]float64, isEasternCell bool) [2]float64 {
	lon, lat := pt[0], pt[1]
	// lock antimeridian boundary edges cleanly to prevent GIS wrapping artifacts
	if math.Abs(math.Abs(lon)-180.0) < 1e-10 {
		if isEasternCell {
			lon = 180.0
		} else {
			lon = -180.0
		}
	}
	return [2]float64{lon, lat}
}

// IsEasternCell Determines whether a cell ring belongs to an Eastern sector (>0° longitude)
// based on its non-seam interior vertices.
func IsEasternCell(facet uint8, vertices [][2]float64) bool {
	if facet == 4 {
		return true
	}
	if facet == 1 {
		return false
	}

	var sumNonSeamLon float64
	var nonSeamCount int
	for _, pt := range vertices {
		if math.Abs(math.Abs(pt[0])-180.0) >= 1e-10 {
			sumNonSeamLon += pt[0]
			nonSeamCount++
		}
	}

	if nonSeamCount > 0 {
		return (sumNonSeamLon / float64(nonSeamCount)) > 0
	}
	return false
}

// CellBoundary returns the geographic boundary polygon ring for any cell.
// n is the number of points per edge (n >= 2). n=2 returns standard 4 corners.
func CellBoundary[T rhealpix.CellID](d *rhealpix.RhpDGGS, id T, n int) ([][2]float64, error) {
	if n < 2 {
		n = 2
	}

	shape := rhealpix.EllipsoidalCellShape(d, id)

	// for standard Quads and Caps at n=2, use direct 4-corner un-projection
	if n == 2 && (shape == rhealpix.ShapeQuad || shape == rhealpix.ShapeCap) {
		return rhealpix.EllipsoidalVertices(d, id, false)
	}

	// retrieve standard 4 geographic corners to determine quadrant orientation (isEastern)
	corners, err := rhealpix.EllipsoidalVertices(d, id, false)
	if err != nil {
		return nil, fmt.Errorf("failed fetching ellipsoidal vertices for orientation: %w", err)
	}

	var facet uint8
	switch v := any(id).(type) {
	case rhealpix.CellID64:
		facet = v.Facet()
	case rhealpix.CellID128:
		facet = v.Facet()
	}
	isEastern := IsEasternCell(facet, corners)

	// get exact 4 planar corners [UL, UR, LR, LL]
	planarVerts := rhealpix.CellPlanarVertices(d, id)

	// sample planar perimeter by interpolating along the 4 edges
	planarRing := make([][2]float64, 0, 4*(n-1))
	for i := 0; i < 4; i++ {
		p1 := planarVerts[i]
		p2 := planarVerts[(i+1)%4]

		for j := 0; j < n-1; j++ {
			t := float64(j) / float64(n-1)
			px := p1[0] + t*(p2[0]-p1[0])
			py := p1[1] + t*(p2[1]-p1[1])
			planarRing = append(planarRing, [2]float64{px, py})
		}
	}

	// rotate planar ring so it starts at true Geographic NW
	nwIdx, err := rhealpix.NWIndex(d, id)
	if err != nil {
		return nil, err
	}

	rotIdx := (n - 1) * nwIdx
	rotatedPlanar := make([][2]float64, len(planarRing))
	for k := 0; k < len(planarRing); k++ {
		rotatedPlanar[k] = planarRing[(rotIdx+k)%len(planarRing)]
	}

	// un-project rotated planar points and snap antimeridian seam vertices
	geoRing := make([][2]float64, len(rotatedPlanar))
	for k, p := range rotatedPlanar {
		lonRad, latRad, err := d.InverseProjectRadians(p[0], p[1])
		if err != nil {
			return nil, fmt.Errorf("point (%f, %f) for cell %v: %w", p[0], p[1], id, err)
		}
		if d.GetEllipsoid() != nil {
			latRad = d.GetEllipsoid().AuthLatInverse(latRad)
		}

		pt := [2]float64{lonRad * rhealpix.RadToDeg, latRad * rhealpix.RadToDeg}
		geoRing[k] = CleanCellSeamVertex(pt, isEastern)
	}

	// normalise antimeridian longitude jumps across the ring
	return NormaliseRing(geoRing), nil
}

// NormaliseRing locks longitude signs so consecutive points stay within 180° span.
func NormaliseRing(ring [][2]float64) [][2]float64 {
	if len(ring) < 2 {
		return ring
	}

	norm := make([][2]float64, len(ring))
	norm[0] = ring[0]

	for k := 1; k < len(ring); k++ {
		lon := ring[k][0]
		prevLon := norm[k-1][0]

		for lon-prevLon > 180.0 {
			lon -= 360.0
		}
		for lon-prevLon < -180.0 {
			lon += 360.0
		}
		norm[k] = [2]float64{lon, ring[k][1]}
	}

	return norm
}

// CellToPolygon converts any rHEALPix cell (CellID64 or CellID128) into a closed orb.Polygon in geodetic coordinates.
// segmentCount defines how many segments each cell edge is subdivided into (segmentCount <= 1 means 1 segment per edge, i.e., corners only).
func CellToPolygon[T rhealpix.CellID](d *rhealpix.RhpDGGS, id T, segmentCount int) (orb.Polygon, error) {
	var isZero bool
	switch v := any(id).(type) {
	case rhealpix.CellID64:
		isZero = v.IsZero()
	case rhealpix.CellID128:
		isZero = v.IsZero()
	}

	if isZero {
		return nil, fmt.Errorf("zero cell ID")
	}

	// n points per edge = segmentCount + 1 (e.g. segmentCount=1 => n=2 points per edge, corners only)
	pointsPerEdge := segmentCount + 1
	if segmentCount <= 1 {
		pointsPerEdge = 2
	}

	// generate boundary points (interpolated in planar space, then unprojected and ring-normalised)
	boundary, err := CellBoundary(d, id, pointsPerEdge)
	if err != nil {
		return nil, fmt.Errorf("failed to compute cell boundary: %w", err)
	}

	// convert to orb.Ring and close the polygon ring
	ring := make(orb.Ring, len(boundary)+1)
	for i, pt := range boundary {
		ring[i] = orb.Point{pt[0], pt[1]}
	}
	ring[len(boundary)] = ring[0] // close ring

	return orb.Polygon{ring}, nil
}

// Cell64ToPolygon is a concrete helper for 64-bit cell IDs.
func Cell64ToPolygon(d *rhealpix.RhpDGGS, id rhealpix.CellID64, segmentCount int) (orb.Polygon, error) {
	return CellToPolygon[rhealpix.CellID64](d, id, segmentCount)
}

// Cell128ToPolygon is a concrete helper for 128-bit cell IDs.
func Cell128ToPolygon(d *rhealpix.RhpDGGS, id rhealpix.CellID128, segmentCount int) (orb.Polygon, error) {
	return CellToPolygon[rhealpix.CellID128](d, id, segmentCount)
}

// Level0FacetBound returns the exact geographic bounding envelope for a Level 0 base facet (0..5).
func Level0FacetBound(d *rhealpix.RhpDGGS, facet uint8) orb.Bound {
	bounds := d.MustFacetBounds(facet)
	capLatDeg := d.GetEllipsoid().Phi0 * rhealpix.RadToDeg

	switch facet {
	case 0: // Facet N (North Cap covers full lon spectrum down to Phi0)
		return orb.Bound{Min: orb.Point{-180.0, capLatDeg}, Max: orb.Point{180.0, 90.0}}

	case 5: // Facet S (South Cap covers full lon spectrum up to -Phi0)
		return orb.Bound{Min: orb.Point{-180.0, -90.0}, Max: orb.Point{180.0, -capLatDeg}}

	default: // Equatorial Facets 1..4 (O, P, Q, R)
		// unproject planar X bounds to geographic longitude
		lonMinRad, _, _ := d.InverseProjectRadians(bounds.MinX(), 0.0)
		lonMaxRad, _, _ := d.InverseProjectRadians(bounds.MaxX(), 0.0)

		return orb.Bound{
			Min: orb.Point{lonMinRad * rhealpix.RadToDeg, -capLatDeg},
			Max: orb.Point{lonMaxRad * rhealpix.RadToDeg, capLatDeg},
		}
	}
}

// CellToBound calculates the orb.Bound envelope for any CellID (CellID64 or CellID128).
// For Level 0..2, it uses edge-densified polygons to capture global projection edge curvature.
// For Level 3+, it calculates the envelope directly from the cell's geographic vertices.
func CellToBound[T rhealpix.CellID](d *rhealpix.RhpDGGS, id T) (orb.Bound, error) {
	var isZero bool
	var res uint8
	var facet uint8

	switch v := any(id).(type) {
	case rhealpix.CellID64:
		isZero = v.IsZero()
		res = v.Resolution()
		facet = v.Facet()
	case rhealpix.CellID128:
		isZero = v.IsZero()
		res = v.Resolution()
		facet = v.Facet()
	}

	if isZero {
		return orb.Bound{}, fmt.Errorf("zero cell ID")
	}

	// Level 0 base facets: use exact geodetic cap boundaries
	if res == 0 {
		return Level0FacetBound(d, facet), nil
	}

	// coarse levels (Level 1..2): densify edges to handle non-linear boundary curvature
	if res <= 2 {
		poly, err := CellToPolygon(d, id, 8)
		if err != nil {
			return orb.Bound{}, err
		}
		return poly.Bound(), nil
	}

	// deep sub-cells (Level 3+)
	var vertices [][2]float64
	var err error

	switch v := any(id).(type) {
	case rhealpix.CellID64:
		vertices, err = d.GeographicVertices64(v)
	case rhealpix.CellID128:
		vertices, err = d.GeographicVertices128(v)
	}
	if err != nil {
		return orb.Bound{}, fmt.Errorf("failed to get geographic vertices for bound: %w", err)
	}

	minX, maxX := 180.0, -180.0
	minY, maxY := 90.0, -90.0

	for _, pt := range vertices {
		lon, lat := pt[0], pt[1]

		if lon < minX {
			minX = lon
		}
		if lon > maxX {
			maxX = lon
		}
		if lat < minY {
			minY = lat
		}
		if lat > maxY {
			maxY = lat
		}
	}

	return orb.Bound{
		Min: orb.Point{minX, minY},
		Max: orb.Point{maxX, maxY},
	}, nil
}

// Cell64ToBound calculates the orb.Bound envelope for a CellID64 cell.
func Cell64ToBound(d *rhealpix.RhpDGGS, id rhealpix.CellID64) (orb.Bound, error) {
	return CellToBound[rhealpix.CellID64](d, id)
}

// Cell128ToBound calculates the orb.Bound envelope for a CellID128 cell.
func Cell128ToBound(d *rhealpix.RhpDGGS, id rhealpix.CellID128) (orb.Bound, error) {
	return CellToBound[rhealpix.CellID128](d, id)
}

// extractGeometryVerticesDensify adds a little more robustness for geometries
// that cross the anti-meridian by adding more sample points (densifying)
// across the edges.
func extractGeometryVerticesDensify(geom orb.Geometry) ([]orb.Point, error) {
	if geom == nil {
		return nil, fmt.Errorf("nil geometry provided")
	}

	var vertices []orb.Point
	switch g := geom.(type) {
	case orb.Point:
		vertices = []orb.Point{g}
	case orb.MultiPoint:
		vertices = g
	case orb.LineString:
		vertices = densifyBoundEdges(g, 10)
	case orb.Polygon:
		if len(g) > 0 {
			vertices = densifyBoundEdges(g[0], 10)
		}
	case orb.MultiPolygon:
		for _, poly := range g {
			if len(poly) > 0 {
				vertices = append(vertices, densifyBoundEdges(poly[0], 10)...)
			}
		}
	case orb.Bound:
		if g.Min.X() > g.Max.X() {
			// split Antimeridian crossing bounds into East and West boxes
			b1Corners := []orb.Point{g.Min, {180.0, g.Min.Y()}, {180.0, g.Max.Y()}, {g.Min.X(), g.Max.Y()}}
			b2Corners := []orb.Point{{-180.0, g.Min.Y()}, {g.Max.X(), g.Min.Y()}, g.Max, {-180.0, g.Max.Y()}}
			vertices = append(densifyBoundEdges(b1Corners, 10), densifyBoundEdges(b2Corners, 10)...)
		} else {
			corners := []orb.Point{g.Min, {g.Max.X(), g.Min.Y()}, g.Max, {g.Min.X(), g.Max.Y()}}
			vertices = densifyBoundEdges(corners, 10)
		}
	default:
		b := geom.Bound()
		return extractGeometryVerticesDensify(b)
	}

	if len(vertices) == 0 {
		return nil, fmt.Errorf("geometry contains no vertices")
	}
	return vertices, nil
}

func extractGeometryVertices(geom orb.Geometry) ([]orb.Point, error) {
	if geom == nil {
		return nil, fmt.Errorf("nil geometry provided")
	}

	var vertices []orb.Point
	switch g := geom.(type) {
	case orb.Point:
		vertices = []orb.Point{g}
	case orb.MultiPoint:
		vertices = g
	case orb.LineString:
		vertices = g
	case orb.Polygon:
		if len(g) > 0 {
			vertices = g[0]
		}
	case orb.MultiPolygon:
		for _, poly := range g {
			if len(poly) > 0 {
				vertices = append(vertices, poly[0]...)
			}
		}
	case orb.Bound:
		if g.Min.X() > g.Max.X() {
			// crosses antimeridian: split into East and West bounding boxes
			vertices = []orb.Point{
				// East Box: [Min.X, 180]
				g.Min,
				{180.0, g.Min.Y()},
				{180.0, g.Max.Y()},
				{g.Min.X(), g.Max.Y()},
				// West Box: [-180, Max.X]
				{-180.0, g.Min.Y()},
				{g.Max.X(), g.Min.Y()},
				g.Max,
				{-180.0, g.Max.Y()},
			}
		} else {
			vertices = []orb.Point{
				g.Min,
				{g.Max.X(), g.Min.Y()},
				g.Max,
				{g.Min.X(), g.Max.Y()},
			}
		}
	default:
		b := geom.Bound()
		return extractGeometryVertices(b)
	}

	if len(vertices) == 0 {
		return nil, fmt.Errorf("geometry contains no vertices")
	}
	return vertices, nil
}

func facetLocalToLonLat(d *rhealpix.RhpDGGS, facet uint8, xLocal, yLocal float64) (float64, float64) {
	// convert local facet coordinates [0, 1] x [0, 1] to global planar radians
	xRhp, yRhp := d.FacetLocalToRadians(facet, xLocal, yLocal)

	// scale planar radians to projected metres on WGS84
	// inverse project planar radians to relative lonRad and authalic latitude (beta)
	// (InverseProjectRadians automatically applies northSquare, southSquare, and Lon0)
	lonRad, beta, _ := d.InverseProjectRadians(xRhp, yRhp)

	// convert authalic latitude back to geodetic latitude & convert to degrees
	latRad := d.GetEllipsoid().AuthLatInverse(beta)

	// inverse project metres to true geodetic Lon/Lat (degrees)
	return lonRad * rhealpix.RadToDeg, latRad * rhealpix.RadToDeg
}
