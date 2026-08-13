package rhealpixorb

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// Cell64ToPolygon converts a CellID64 into an orb.Polygon in WGS84 coordinates.
// It handles polar shape classification (Quads, Caps, Darts, and Skew Quads),
// applying vertex re-ordering and dart trimming for clean GIS geometries.
//
// segmentCount defines how many linear sub-segments to create per cell edge.
// Use segmentCount = 1 for standard 4-corner polygons, or segmentCount = 8..16
// for low-resolution / polar cells that curve significantly on a globe.
func Cell64ToPolygon(el *rhealpix.Ellipsoid, id rhealpix.CellID64, segmentCount int) (orb.Polygon, error) {
	if id.IsZero() {
		return nil, fmt.Errorf("zero cell ID")
	}

	facet := id.Facet()
	shape := rhealpix.CellShape64(id)

	// calculate local planar square extent [xMin, xMax], [yMin, yMax] in [0, 1]
	xMin, xMax, yMin, yMax := cellPlanarExtent64(id)
	width := xMax - xMin

	// initial planar corner vertices: UL, UR, LR, LL
	planarCorners := [4][2]float64{
		{xMin, yMax}, // 0: Top-Left (UL)
		{xMax, yMax}, // 1: Top-Right (UR)
		{xMax, yMin}, // 2: Bottom-Right (LR)
		{xMin, yMin}, // 3: Bottom-Left (LL)
	}

	// re-order planar vertices for topological alignment (Skew Quads & Darts)
	alignedCorners := alignVertices(planarCorners, shape, facet, width, el.RA)

	// extract valid, non-collapsed vertices
	validCorners := make([][2]float64, 0, 4)
	for i, pt := range alignedCorners {
		// trim degenerate/collapsed seam vertex on Darts
		if shape == rhealpix.ShapeDart {
			if (facet == 0 && i == 2) || (facet == 5 && i == 1) {
				continue // skip collapsed vertex
			}
		}
		validCorners = append(validCorners, pt)
	}

	// construct ring with edge densification in planar space, then project
	ring := buildDensifiedRing(el, facet, validCorners, segmentCount)
	return orb.Polygon{ring}, nil
}

// Cell128ToPolygon converts a CellID128 into an orb.Polygon with polar shape handling.
func Cell128ToPolygon(el *rhealpix.Ellipsoid, id rhealpix.CellID128, segmentCount int) (orb.Polygon, error) {
	if id.IsZero() {
		return nil, fmt.Errorf("zero cell ID")
	}

	facet := id.Facet()
	shape := rhealpix.CellShape128(id)

	xMin, xMax, yMin, yMax := cellPlanarExtent128(id)
	width := xMax - xMin

	planarCorners := [4][2]float64{
		{xMin, yMax},
		{xMax, yMax},
		{xMax, yMin},
		{xMin, yMin},
	}

	alignedCorners := alignVertices(planarCorners, shape, facet, width, el.RA)

	validCorners := make([][2]float64, 0, 4)
	for i, pt := range alignedCorners {
		if shape == rhealpix.ShapeDart {
			if (facet == 0 && i == 2) || (facet == 5 && i == 1) {
				continue
			}
		}
		validCorners = append(validCorners, pt)
	}

	ring := buildDensifiedRing(el, facet, validCorners, segmentCount)
	return orb.Polygon{ring}, nil
}

// CellToPolygon converts any rHEALPix cell (CellID64 or CellID128) into an orb.Polygon.
func CellToPolygon[T rhealpix.CellID64 | rhealpix.CellID128](el *rhealpix.Ellipsoid, id T, segmentCount int) (orb.Polygon, error) {
	switch cell := any(id).(type) {
	case rhealpix.CellID64:
		return Cell64ToPolygon(el, cell, segmentCount)
	case rhealpix.CellID128:
		return Cell128ToPolygon(el, cell, segmentCount)
	default:
		return nil, fmt.Errorf("unsupported cell ID type")
	}
}

// Cell64ToBound calculates the orb.Bound envelope for a CellID64 cell.
// For low resolutions (Level 0..2), it uses densified edges to handle global projection curvature.
// For deeper resolutions (Level 3+), it directly projects the 4 corner vertices for maximum speed.
func Cell64ToBound(el *rhealpix.Ellipsoid, id rhealpix.CellID64) (orb.Bound, error) {
	if id.IsZero() {
		return orb.Bound{}, fmt.Errorf("zero cell ID")
	}

	res := id.Resolution()

	// for root/coarse levels, use edge densification to capture curved global boundaries
	if res <= 2 {
		poly, err := Cell64ToPolygon(el, id, 8)
		if err != nil {
			return orb.Bound{}, err
		}
		return poly.Bound(), nil
	}

	// for deep sub-cells (Level 3+), directly project the 4 planar corners
	facet := id.Facet()
	xMin, xMax, yMin, yMax := cellPlanarExtent64(id)

	// local planar corners: UL, UR, LR, LL
	corners := [4][2]float64{
		{xMin, yMax},
		{xMax, yMax},
		{xMax, yMin},
		{xMin, yMin},
	}

	// this might seem counter intuitive, but initialising min to the max value
	// avoids potential global wrap arounds
	minX, maxX := 180.0, -180.0
	minY, maxY := 90.0, -90.0

	for _, pt := range corners {
		lon, lat := facetLocalToLonLat(el, facet, pt[0], pt[1])
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

// Cell128ToBound calculates the orb.Bound envelope for a CellID128 cell.
// For low resolutions (Level 0..2), it uses densified edges to handle global projection curvature.
// For deeper resolutions (Level 3+), it directly projects the 4 corner vertices for maximum speed.
func Cell128ToBound(el *rhealpix.Ellipsoid, id rhealpix.CellID128) (orb.Bound, error) {
	if id.IsZero() {
		return orb.Bound{}, fmt.Errorf("zero cell ID")
	}

	res := id.Resolution()

	// for root/coarse levels, use edge densification to capture curved global boundaries
	if res <= 2 {
		poly, err := Cell128ToPolygon(el, id, 8)
		if err != nil {
			return orb.Bound{}, err
		}
		return poly.Bound(), nil
	}

	// for deep sub-cells (Level 3+), directly project the 4 planar corners
	facet := id.Facet()
	xMin, xMax, yMin, yMax := cellPlanarExtent128(id)

	// local planar corners: UL, UR, LR, LL
	corners := [4][2]float64{
		{xMin, yMax},
		{xMax, yMax},
		{xMax, yMin},
		{xMin, yMin},
	}

	// this might seem counter intuitive, but initialising min to the max value
	// avoids potential global wrap arounds
	minX, maxX := 180.0, -180.0
	minY, maxY := 90.0, -90.0

	for _, pt := range corners {
		lon, lat := facetLocalToLonLat(el, facet, pt[0], pt[1])
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

// DeriveRegionCode64FromGeometry calculates the canonical region_code string for any orb.Geometry footprint.
func DeriveRegionCode64FromGeometry(el *rhealpix.Ellipsoid, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	vertices, err := extractGeometryVertices(geom)
	if err != nil {
		return "", nil, err
	}

	facetGroups := make(map[uint8][]rhealpix.CellID64)
	for _, pt := range vertices {
		cell, err := rhealpix.ForwardTransform64(el, pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, err
		}
		facetGroups[cell.Facet()] = append(facetGroups[cell.Facet()], cell)
	}

	var suids []string
	for _, groupCells := range facetGroups {
		var lcaCell rhealpix.CellID64
		lcaCell = groupCells[0]
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
func DeriveRegionCode128FromGeometry(el *rhealpix.Ellipsoid, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	vertices, err := extractGeometryVertices(geom)
	if err != nil {
		return "", nil, err
	}

	facetGroups := make(map[uint8][]rhealpix.CellID128)
	for _, pt := range vertices {
		cell, err := rhealpix.ForwardTransform128(el, pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, err
		}
		facetGroups[cell.Facet()] = append(facetGroups[cell.Facet()], cell)
	}

	var suids []string
	for _, groupCells := range facetGroups {
		var lcaCell rhealpix.CellID128
		lcaCell = groupCells[0]
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
func DeriveRegionCode(el *rhealpix.Ellipsoid, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	if targetRes <= rhealpix.MaxResolution64 {
		return DeriveRegionCode64FromGeometry(el, geom, targetRes)
	}
	return DeriveRegionCode128FromGeometry(el, geom, targetRes)
}

// --- Internal Helpers ---

func cellPlanarExtent64(id rhealpix.CellID64) (xMin, xMax, yMin, yMax float64) {
	xMin, xMax = 0.0, 1.0
	yMin, yMax = 0.0, 1.0
	res := id.Resolution()
	rawID := uint64(id)

	for i := uint8(0); i < res; i++ {
		shift := 52 - (i * 4)
		subCell := uint8((rawID >> shift) & rhealpix.SubCellMask)

		xIdx := float64(subCell % 3)
		yIdx := float64(subCell / 3)

		spanX := (xMax - xMin) / 3.0
		spanY := (yMax - yMin) / 3.0

		xMin = xMin + (xIdx * spanX)
		xMax = xMin + spanX
		yMin = yMin + (yIdx * spanY)
		yMax = yMin + spanY
	}
	return
}

func cellPlanarExtent128(id rhealpix.CellID128) (xMin, xMax, yMin, yMax float64) {
	xMin, xMax = 0.0, 1.0
	yMin, yMax = 0.0, 1.0
	res := id.Resolution()

	for i := uint8(0); i < res; i++ {
		var subCell uint8
		if i < 14 {
			shift := 52 - (i * 4)
			subCell = uint8((id.High >> shift) & rhealpix.SubCellMask)
		} else {
			shift := 60 - ((i - 14) * 4)
			subCell = uint8((id.Low >> shift) & rhealpix.SubCellMask)
		}

		xIdx := float64(subCell % 3)
		yIdx := float64(subCell / 3)

		spanX := (xMax - xMin) / 3.0
		spanY := (yMax - yMin) / 3.0

		xMin = xMin + (xIdx * spanX)
		xMax = xMin + spanX
		yMin = yMin + (yIdx * spanY)
		yMax = yMin + spanY
	}

	// now we may get data/floating point drift, maybe noticeable at >= level 25 ...
	// TODO: need to investigate a potential alternate approach
	return
}

// alignVertices re-orders vertices for Skew Quads and Darts so that geographic NW is preserved.
func alignVertices(corners [4][2]float64, shape rhealpix.CellShape, facet uint8, width, authalicRadius float64) [4][2]float64 {
	if shape == rhealpix.ShapeQuad || shape == rhealpix.ShapeCap {
		return corners
	}

	idx := 0

	if shape == rhealpix.ShapeSkewQuad {
		// Calculate nucleus (cell center in local planar space)
		centerX := corners[3][0] + (width / 2.0)
		centerY := corners[3][1] + (width / 2.0)

		xpt := centerX * (2.0 * math.Pi)
		ypt := centerY * (2.0 * math.Pi)
		eps := 1e-15

		var triangleNum int
		if facet == 0 { // North Pole ('N')
			l1 := xpt - (-3.0 * math.Pi / 4.0)
			l2 := -xpt + (-3.0 * math.Pi / 4.0)
			if ypt < l1-eps && ypt >= l2-eps {
				triangleNum = 1
			} else if ypt >= l1-eps && ypt > l2+eps {
				triangleNum = 2
			} else if ypt > l1+eps && ypt <= l2+eps {
				triangleNum = 3
			} else {
				triangleNum = 0
			}
			idx = (4 - triangleNum) % 4
		} else if facet == 5 { // South Pole ('S')
			l1 := xpt - (-3.0 * math.Pi / 4.0)
			l2 := -xpt + (-3.0 * math.Pi / 4.0)
			if ypt <= l1+eps && ypt > l2+eps {
				triangleNum = 1
			} else if ypt < l1-eps && ypt <= l2+eps {
				triangleNum = 2
			} else if ypt >= l1-eps && ypt < l2-eps {
				triangleNum = 3
			} else {
				triangleNum = 0
			}
			idx = triangleNum % 4
		}
	} else if shape == rhealpix.ShapeDart {
		// find most poleward vertex index along Y
		maxVal := -1.0
		for i, pt := range corners {
			val := math.Abs(pt[1] - 0.5) // distance from center line
			if val > maxVal {
				maxVal = val
				idx = i
			}
		}
		if facet == 5 { // South Pole
			idx = (idx + 1) % 4
		}
	}

	// rotate corner array starting at calculated NW index
	var rotated [4][2]float64
	for i := 0; i < 4; i++ {
		rotated[i] = corners[(idx+i)%4]
	}
	return rotated
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

func facetLocalToLonLat(el *rhealpix.Ellipsoid, facet uint8, xLocal, yLocal float64) (float64, float64) {
	// convert local facet coordinates [0, 1] x [0, 1] to global planar radians
	xRhp, yRhp := rhealpix.FacetLocalToRadians(facet, xLocal, yLocal)

	// scale planar radians to projected meters on WGS84
	xMeters := xRhp * el.RA
	yMeters := yRhp * el.RA

	// inverse project meters to true geodetic Lon/Lat (degrees)
	return el.InverseProject(xMeters, yMeters)
}
