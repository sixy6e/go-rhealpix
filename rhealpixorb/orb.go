package rhealpixorb

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// CellToPolygon converts any rHEALPix cell (CellID64 or CellID128) into a closed orb.Polygon in geodetic coordinates.
// If segmentCount <= 1, it returns the fundamental cell boundary (3 points for Dart, 4 points for Quad/Cap/SkewQuad).
// If segmentCount > 1, it densifies each cell edge into segmentCount linear steps in geographic space.
func CellToPolygon[T rhealpix.CellID](d *rhealpix.RhpDGGS, id T, segmentCount int) (orb.Polygon, error) {
	var isZero bool
	var facet uint8

	switch v := any(id).(type) {
	case rhealpix.CellID64:
		isZero = v.IsZero()
		facet = v.Facet()
	case rhealpix.CellID128:
		isZero = v.IsZero()
		facet = v.Facet()
	}

	if isZero {
		return nil, fmt.Errorf("zero cell ID")
	}

	// retrieve topologically aligned, GIS-valid vertices starting at true NW (Dart trimmed to 3 points)
	var vertices [][2]float64
	var err error

	switch v := any(id).(type) {
	case rhealpix.CellID64:
		vertices, err = d.GeographicVertices64(v)
	case rhealpix.CellID128:
		vertices, err = d.GeographicVertices128(v)
	}
	if err != nil {
		return nil, fmt.Errorf("failed to compute geographic vertices: %w", err)
	}

	numVerts := len(vertices)
	if numVerts == 0 {
		return nil, fmt.Errorf("empty vertex set returned")
	}

	// helper to snap seam vertices for Facet 4 (R) from -180.0 -> +180.0
	// (antimeridian)
	cleanPoint := func(pt [2]float64) [2]float64 {
		lon, lat := pt[0], pt[1]
		if facet == 4 && math.Abs(lon - -180.0) < 1e-11 {
			lon = 180.0
		}
		return [2]float64{lon, lat}
	}

	// standard non-densified polygon output (segmentCount <= 1)
	if segmentCount <= 1 {
		ring := make(orb.Ring, 0, numVerts+1)
		for _, pt := range vertices {
			cleaned := cleanPoint(pt)
			ring = append(ring, orb.Point{cleaned[0], cleaned[1]})
		}
		// close ring
		ring = append(ring, ring[0])
		return orb.Polygon{ring}, nil
	}

	// densified edge interpolation along authentic topological boundary edges
	totalPts := numVerts * segmentCount
	ring := make(orb.Ring, 0, totalPts+1)

	for i := 0; i < numVerts; i++ {
		start := cleanPoint(vertices[i])
		end := cleanPoint(vertices[(i+1)%numVerts])

		lon1, lat1 := start[0], start[1]
		lon2, lat2 := end[0], end[1]

		// interpolate intermediate points along current boundary edge
		for j := 0; j < segmentCount; j++ {
			t := float64(j) / float64(segmentCount)
			interpLon := lon1 + t*(lon2-lon1)
			interpLat := lat1 + t*(lat2-lat1)
			ring = append(ring, orb.Point{interpLon, interpLat})
		}
	}

	// close ring for GeoJSON standard compliance
	ring = append(ring, ring[0])

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

// DeriveRegionCode64FromGeometry calculates the canonical region_code string for any orb.Geometry footprint.
func DeriveRegionCode64FromGeometry(d *rhealpix.RhpDGGS, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	vertices, err := extractGeometryVertices(geom)
	if err != nil {
		return "", nil, err
	}

	facetGroups := make(map[uint8][]rhealpix.CellID64)
	for _, pt := range vertices {
		cell, err := d.ForwardTransform64(pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, err
		}
		facetGroups[cell.Facet()] = append(facetGroups[cell.Facet()], cell)
	}

	var suids []string
	for _, groupCells := range facetGroups {
		lcaCell := groupCells[0]
		for _, cell := range groupCells[1:] {
			var lcaErr error
			lcaCell, lcaErr = rhealpix.CommonAncestor64(lcaCell, cell)
			if lcaErr != nil {
				lcaCell = cell
			}
		}
		suids = append(suids, lcaCell.String())
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
}

// DeriveRegionCode128FromGeometry calculates the canonical region_code string using 128-bit precision.
func DeriveRegionCode128FromGeometry(d *rhealpix.RhpDGGS, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	vertices, err := extractGeometryVertices(geom)
	if err != nil {
		return "", nil, err
	}

	facetGroups := make(map[uint8][]rhealpix.CellID128)
	for _, pt := range vertices {
		cell, err := d.ForwardTransform128(pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, err
		}
		facetGroups[cell.Facet()] = append(facetGroups[cell.Facet()], cell)
	}

	var suids []string
	for _, groupCells := range facetGroups {
		lcaCell := groupCells[0]
		for _, cell := range groupCells[1:] {
			var lcaErr error
			lcaCell, lcaErr = rhealpix.CommonAncestor128(lcaCell, cell)
			if lcaErr != nil {
				lcaCell = cell
			}
		}
		suids = append(suids, lcaCell.String())
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
}

// DeriveRegionCode automatically selects 64-bit or 128-bit derivation based on resolution depth.
func DeriveRegionCode(d *rhealpix.RhpDGGS, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	if targetRes <= rhealpix.MaxResolution64 {
		return DeriveRegionCode64FromGeometry(d, geom, targetRes)
	}
	return DeriveRegionCode128FromGeometry(d, geom, targetRes)
}

// DeriveRegionCodeAtRes calculates a fixed-resolution spatial cover string for any orb.Geometry footprint.
// Each decomposed cell in the string is joined by a hyphen (-).
// Output examples for target resolution level 3 (1 facet + 3 levels):
// O003-O006
// S005-S006-S008
// S000-S002-S006-S008
// Output examples for target resolution level 2 (1 facet + 2 levels):
// S03-S04
// S03-S04-S06-S07
// O07
// S03-S04-S06-S07
func DeriveRegionCodeAtRes(
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8, // e.g. 3 or 4 for Landsat
) (string, []string, error) {
	// walk/decompose geometry up to targetRes
	cells, err := STACGeometryToTileDBRangesTopDown128(d, geom, targetRes)
	// cells, err := STACGeometryToTileDBRangesTopDown(el, geom, targetRes)  // testing equivalencce
	if err != nil {
		return "", nil, err
	}

	// extract unique cell string representations at targetRes
	suidSet := make(map[string]struct{})
	for _, cell := range cells {
		suidSet[cell.String()] = struct{}{}
	}

	suids := make([]string, 0, len(suidSet))
	for suid := range suidSet {
		suids = append(suids, suid)
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
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

// DeriveCentroidRegionCode returns the single rHEALPix cell SUID for the geometry's centroid at targetRes.
func DeriveCentroidRegionCode(
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8,
) (string, error) {
	if geom == nil {
		return "", fmt.Errorf("geometry is nil")
	}

	// compute WGS84 bounding box centroid
	center := geom.Bound().Center()

	// projection from (Lon, Lat) to CellID128 (rhealpix space)
	cellID, err := d.ForwardTransform128(center.X(), center.Y(), targetRes)
	if err != nil {
		return "", fmt.Errorf("failed projecting centroid to cell: %w", err)
	}

	return cellID.String(), nil
}

// DeriveBoundingCornerCodes returns the unique rHEALPix cell IDs (max 4, or max 8 if crossing the antimeridian)
// containing the corners of the geometry's WGS84 bounding box at targetRes.
// Each decomposed cell in the string is joined by a hyphen (-).
func DeriveBoundingCornerCodes(
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8,
) (string, []string, error) {
	if geom == nil {
		return "", nil, fmt.Errorf("geometry is nil")
	}

	bound := geom.Bound()
	var corners []orb.Point

	if bound.Min.X() > bound.Max.X() {
		// Antimeridian crossing: split into East [Min.X, 180] and West [-180, Max.X] corner sets
		corners = []orb.Point{
			bound.Min,
			{180.0, bound.Min.Y()},
			{180.0, bound.Max.Y()},
			{bound.Min.X(), bound.Max.Y()},
			{-180.0, bound.Min.Y()},
			{bound.Max.X(), bound.Min.Y()},
			bound.Max,
			{-180.0, bound.Max.Y()},
		}
	} else {
		corners = []orb.Point{
			bound.Min,
			{bound.Max.X(), bound.Min.Y()},
			bound.Max,
			{bound.Min.X(), bound.Max.Y()},
		}
	}

	cellSet := make(map[string]struct{}, len(corners))

	for _, pt := range corners {
		cellID, err := d.ForwardTransform128(pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, fmt.Errorf("failed projecting corner point %v: %w", pt, err)
		}
		cellSet[cellID.String()] = struct{}{}
	}

	suids := make([]string, 0, len(cellSet))
	for suid := range cellSet {
		suids = append(suids, suid)
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
}
