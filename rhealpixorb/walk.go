package rhealpixorb

import (
	"fmt"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// PlanarBox represents an axis-aligned bounding box (AABB) in normalised [0, 1] x [0, 1] local facet space.
type PlanarBox struct {
	MinX, MaxX float64
	MinY, MaxY float64
}

// Subdivide3x3 splits a PlanarBox into a 3x3 grid of 9 child boxes matching the rHEALPix child order (0..8).
// Uses exact fraction multiplication to prevent floating-point boundary drift at high resolutions.
func (pb PlanarBox) Subdivide3x3() [9]PlanarBox {
	var children [9]PlanarBox
	spanX := pb.MaxX - pb.MinX
	spanY := pb.MaxY - pb.MinY

	idx := 0
	for row := 0; row < 3; row++ { // top-to-bottom row ordering
		y0 := pb.MinY + (float64(row) * spanY / 3.0)
		y1 := pb.MinY + (float64(row+1) * spanY / 3.0)

		for col := 0; col < 3; col++ { // left-to-right column ordering
			x0 := pb.MinX + (float64(col) * spanX / 3.0)
			x1 := pb.MinX + (float64(col+1) * spanX / 3.0)

			children[idx] = PlanarBox{MinX: x0, MaxX: x1, MinY: y0, MaxY: y1}
			idx++
		}
	}
	return children
}

// ToBound converts a PlanarBox to an orb.Bound.
func (pb PlanarBox) ToBound() orb.Bound {
	return orb.Bound{
		Min: orb.Point{pb.MinX, pb.MinY},
		Max: orb.Point{pb.MaxX, pb.MaxY},
	}
}

type BoxRelation int

const (
	BoxOutside BoxRelation = iota
	BoxInside
	BoxStraddles
)

// PlanarSegment holds pre-computed axis bounds for a single edge segment.
type PlanarSegment struct {
	MinX, MaxX float64
	MinY, MaxY float64
	P1, P2     orb.Point
}

// PreparedGeometry wraps complex polygons or multipolygons (e.g. Landsat 7 SLC-Off with holes/gaps)
// into flat segment buffers and bounding envelopes for zero-allocation tree walking.
// (initial attempts would calculate bounds, for the same input object,
// for every cell, for every resolution ...).
type PreparedGeometry struct {
	Polygons []orb.Polygon
	Bound    orb.Bound
	Segments []PlanarSegment
}

// NewPreparedGeometry prepares an orb.Geometry for faster spatial tree walking.
func NewPreparedGeometry(geom orb.Geometry) PreparedGeometry {
	if geom == nil {
		return PreparedGeometry{}
	}

	var polys []orb.Polygon
	switch g := geom.(type) {
	case orb.Polygon:
		polys = append(polys, g)
	case orb.MultiPolygon:
		polys = append(polys, g...)
	}

	bound := geom.Bound()
	var segments []PlanarSegment

	for _, poly := range polys {
		for _, ring := range poly {
			for i := 0; i < len(ring)-1; i++ {
				p1, p2 := ring[i], ring[i+1]

				minX, maxX := p1[0], p2[0]
				if minX > maxX {
					minX, maxX = maxX, minX
				}

				minY, maxY := p1[1], p2[1]
				if minY > maxY {
					minY, maxY = maxY, minY
				}

				segments = append(segments, PlanarSegment{
					MinX: minX,
					MaxX: maxX,
					MinY: minY,
					MaxY: maxY,
					P1:   p1,
					P2:   p2,
				})
			}
		}
	}

	return PreparedGeometry{
		Polygons: polys,
		Bound:    bound,
		Segments: segments,
	}
}

// pointInRing performs an inlined crossings ray-casting test.
func pointInRing(pt orb.Point, ring orb.Ring) bool {
	x, y := pt[0], pt[1]
	inside := false
	for i, j := 0, len(ring)-1; i < len(ring); j, i = i, i+1 {
		xi, yi := ring[i][0], ring[i][1]
		xj, yj := ring[j][0], ring[j][1]

		intersect := ((yi > y) != (yj > y)) &&
			(x < (xj-xi)*(y-yi)/(yj-yi)+xi)
		if intersect {
			inside = !inside
		}
	}
	return inside
}

// pointInMultiPolygon checks point containment across multi-polygons and inner holes.
func pointInMultiPolygon(pt orb.Point, polys []orb.Polygon) bool {
	for _, poly := range polys {
		if len(poly) == 0 {
			continue
		}

		// must be inside the outer shell ring (ring 0)
		if !pointInRing(pt, poly[0]) {
			continue
		}

		// must NOT fall inside any hole ring (rings 1..N)
		inHole := false
		for i := 1; i < len(poly); i++ {
			if pointInRing(pt, poly[i]) {
				inHole = true
				break
			}
		}

		if !inHole {
			return true // point falls cleanly inside valid data area
		}
	}

	return false
}

// EvaluateBoxRelation checks box relation against a PreparedGeometry.
func EvaluateBoxRelation(box PlanarBox, prep PreparedGeometry) BoxRelation {
	// instant AABB Envelope Rejection:
	// discard non-overlapping sub-trees at Res 0-3 before doing deeper evaluations!
	if box.MaxX < prep.Bound.Min[0] || box.MinX > prep.Bound.Max[0] ||
		box.MaxY < prep.Bound.Min[1] || box.MinY > prep.Bound.Max[1] {
		return BoxOutside
	}

	// evaluate 4 corners + centre point using stack array
	midX := (box.MinX + box.MaxX) / 2.0
	midY := (box.MinY + box.MaxY) / 2.0

	corners := [5]orb.Point{
		{box.MinX, box.MaxY},
		{box.MaxX, box.MaxY},
		{box.MaxX, box.MinY},
		{box.MinX, box.MinY},
		{midX, midY},
	}

	insideCount := 0
	for _, pt := range corners {
		if pointInMultiPolygon(pt, prep.Polygons) {
			insideCount++
		}
	}

	if insideCount == 5 {
		// verify no outer OR inner hole segment cuts through cell
		if !intersectsAnySegment(prep.Segments, box) {
			return BoxInside
		}
		return BoxStraddles
	} else if insideCount == 0 {
		if intersectsAnySegment(prep.Segments, box) {
			return BoxStraddles
		}
		return BoxOutside
	}

	return BoxStraddles
}

// segmentIntersectsBox performs a Liang-Barsky line clipping test to check
// if directed segment P1->P2 intersects or touches the PlanarBox.
// https://en.wikipedia.org/wiki/Liang%E2%80%93Barsky_algorithm
// https://www.geeksforgeeks.org/computer-graphics/liang-barsky-algorithm/
// https://gamedev.stackexchange.com/questions/112528/liang-barsky-line-clipping-algorithm
func segmentIntersectsBox(p1, p2 orb.Point, box PlanarBox) bool {
	// if either endpoint is strictly inside the box, it intersects
	if p1[0] >= box.MinX && p1[0] <= box.MaxX && p1[1] >= box.MinY && p1[1] <= box.MaxY {
		return true
	}
	if p2[0] >= box.MinX && p2[0] <= box.MaxX && p2[1] >= box.MinY && p2[1] <= box.MaxY {
		return true
	}

	dx := p2[0] - p1[0]
	dy := p2[1] - p1[1]

	tMin := 0.0
	tMax := 1.0

	clip := func(p, q float64) bool {
		if p == 0 {
			return q >= 0
		}
		r := q / p
		if p < 0 {
			if r > tMax {
				return false
			}
			if r > tMin {
				tMin = r
			}
		} else {
			if r < tMin {
				return false
			}
			if r < tMax {
				tMax = r
			}
		}
		return tMin <= tMax
	}

	// test against [Left, Right, Bottom, Top] box boundaries
	if !clip(-dx, p1[0]-box.MinX) ||
		!clip(dx, box.MaxX-p1[0]) ||
		!clip(-dy, p1[1]-box.MinY) ||
		!clip(dy, box.MaxY-p1[1]) {
		return false
	}

	return true
}

// intersectsAnySegment checks if any pre-computed polygon edge actually crosses the PlanarBox.
func intersectsAnySegment(segments []PlanarSegment, box PlanarBox) bool {
	for i := range segments {
		seg := &segments[i]

		if seg.MaxX < box.MinX || seg.MinX > box.MaxX ||
			seg.MaxY < box.MinY || seg.MinY > box.MaxY {
			continue
		}

		if segmentIntersectsBox(seg.P1, seg.P2, box) {
			return true
		}
	}
	return false
}

// projectGeometryToPlanar maps a WGS84 geometry to normalised [0, 1] x [0, 1] facet local planar space.
func projectGeometryToPlanar(el *rhealpix.Ellipsoid, facetID uint8, geom orb.Geometry) orb.Geometry {
	if geom == nil {
		return nil
	}

	switch g := geom.(type) {
	case orb.Polygon:
		proj := projectPolygonToPlanar(el, facetID, g)
		return ensurePolygonCCW(proj)
	case orb.MultiPolygon:
		var mp orb.MultiPolygon
		for _, poly := range g {
			proj := projectPolygonToPlanar(el, facetID, poly)
			ccw := ensurePolygonCCW(proj)
			if len(ccw) > 0 {
				mp = append(mp, ccw)
			}
		}
		if len(mp) == 1 {
			return mp[0]
		}
		return mp
	default:
		return nil
	}
}

// projectPolygonToPlanar maps a WGS84 polygon to normalised [0, 1] x [0, 1] facet local planar space.
func projectPolygonToPlanar(el *rhealpix.Ellipsoid, facetID uint8, poly orb.Polygon) orb.Polygon {
	projected := make(orb.Polygon, len(poly))

	for i, ring := range poly {
		projectedRing := make(orb.Ring, len(ring))
		for j, pt := range ring {
			lon, lat := pt[0], pt[1]

			xLoc, yLoc, err := rhealpix.LonLatToPlanar(el, facetID, lon, lat)
			if err != nil {
				projectedRing[j] = orb.Point{lon, lat}
				continue
			}

			projectedRing[j] = orb.Point{xLoc, yLoc}
		}
		projected[i] = projectedRing
	}

	return projected
}

// WalkPlanarCell recursively walks the rHEALPix nonary tree in local planar space.
// Similar methodology to Uber's H3 PolyFill. But simpler as we have squares, not overlapping
// hexagons and pentagons.
// TODO; test dart shaped cells.
func WalkPlanarCell(
	cell rhealpix.CellID64,
	box PlanarBox,
	prep PreparedGeometry,
	targetRes uint8,
	results []rhealpix.CellID64,
) []rhealpix.CellID64 {
	relation := EvaluateBoxRelation(box, prep)

	switch relation {
	case BoxOutside:
		return results

	case BoxInside:
		return append(results, cell)

	case BoxStraddles:
		if cell.Resolution() >= targetRes {
			return append(results, cell)
		}

		childBoxes := box.Subdivide3x3()

		for i := uint8(0); i < 9; i++ {
			childCell, err := cell.Child(i)
			if err != nil {
				continue
			}

			results = WalkPlanarCell(childCell, childBoxes[i], prep, targetRes, results)
		}
	}

	return results
}

// STACGeometryToTileDBRangesTopDown performs top-down nonary tree spatial decomposition
// of a WGS84 geometry into a compacted set of rHEALPix cell IDs up to targetRes.
//
// It clips the input geometry across rHEALPix ellipsoid facets, projects each facet
// sub-geometry into normalised [0, 1] x [0, 1] local planar space, and constructs a
// PreparedGeometry for zero-allocation point-in-polygon and Liang-Barsky edge checks.
// https://en.wikipedia.org/wiki/Liang%E2%80%93Barsky_algorithm
// https://www.geeksforgeeks.org/computer-graphics/liang-barsky-algorithm/
//
// The resulting cells are compacted (merging nine full child cells into their
// parent cell recursively) to minimise the final range footprint before 1D database
// indexing in TileDB.
// The approach of top down planar decomposition and compaction, is similar to Uber's H3 PolyFill.
func STACGeometryToTileDBRangesTopDown(
	el *rhealpix.Ellipsoid,
	geom orb.Geometry,
	targetRes uint8,
) ([]rhealpix.CellID64, error) {
	facetSubGeoms, err := ClipGeometryToFacets(el, geom)
	if err != nil {
		return nil, fmt.Errorf("failed clipping geometry to facets: %w", err)
	}

	// pre-allocate rawCells slice capacity (128 entries) to avoid dynamic growth allocations
	rawCells := make([]rhealpix.CellID64, 0, 128)
	rootBox := PlanarBox{MinX: 0.0, MaxX: 1.0, MinY: 0.0, MaxY: 1.0}

	for _, sub := range facetSubGeoms {
		rootCell, err := rhealpix.PackCellID64(sub.FacetID, 0, nil)
		if err != nil {
			return nil, err
		}

		planarGeom := projectGeometryToPlanar(el, sub.FacetID, sub.Geometry)
		if planarGeom == nil {
			continue
		}

		prep := NewPreparedGeometry(planarGeom)
		rawCells = WalkPlanarCell(rootCell, rootBox, prep, targetRes, rawCells)
	}

	compacted, err := rhealpix.Compact(rawCells)
	if err != nil {
		return nil, fmt.Errorf("failed compacting decomposed cells: %w", err)
	}

	return compacted, nil
}

// WalkPlanarCell128 recursively walks the rHEALPix nonary tree in local planar space using 128-bit cell IDs.
// Similar methodology to Uber's H3 PolyFill. But simpler as we have squares, not overlapping
// hexagons and pentagons.
// TODO; test dart shaped cells.
func WalkPlanarCell128(
	cell rhealpix.CellID128,
	box PlanarBox,
	prep PreparedGeometry,
	targetRes uint8,
	results []rhealpix.CellID128,
) []rhealpix.CellID128 {
	relation := EvaluateBoxRelation(box, prep)

	switch relation {
	case BoxOutside:
		return results

	case BoxInside:
		return append(results, cell)

	case BoxStraddles:
		if cell.Resolution() >= targetRes {
			return append(results, cell)
		}

		childBoxes := box.Subdivide3x3()

		for i := uint8(0); i < 9; i++ {
			childCell, err := cell.Child(i)
			if err != nil {
				continue
			}

			results = WalkPlanarCell128(childCell, childBoxes[i], prep, targetRes, results)
		}
	}

	return results
}

// STACGeometryToTileDBRangesTopDown128 performs top-down nonary tree spatial decomposition
// of a WGS84 geometry into a compacted set of 128-bit rHEALPix cell IDs up to targetRes.
//
// It clips the input geometry across rHEALPix ellipsoid facets, projects each facet
// sub-geometry into normalised [0, 1] x [0, 1] local planar space, and constructs a
// PreparedGeometry for zero-allocation point-in-polygon and Liang-Barsky edge checks.
// https://en.wikipedia.org/wiki/Liang%E2%80%93Barsky_algorithm
// https://www.geeksforgeeks.org/computer-graphics/liang-barsky-algorithm/
//
// The resulting cells are compacted (merging nine full child cells into their
// parent cell recursively) to minimise the final range footprint before 1D database
// indexing in TileDB.
// The approach of top down planar decomposition and compaction, is similar to Uber's H3 PolyFill.
func STACGeometryToTileDBRangesTopDown128(
	el *rhealpix.Ellipsoid,
	geom orb.Geometry,
	targetRes uint8,
) ([]rhealpix.CellID128, error) {
	facetSubGeoms, err := ClipGeometryToFacets(el, geom)
	if err != nil {
		return nil, fmt.Errorf("failed clipping geometry to facets: %w", err)
	}

	// pre-allocate rawCells slice capacity (128 entries) to avoid dynamic growth allocations
	rawCells := make([]rhealpix.CellID128, 0, 128)
	rootBox := PlanarBox{MinX: 0.0, MaxX: 1.0, MinY: 0.0, MaxY: 1.0}

	for _, sub := range facetSubGeoms {
		rootCell, err := rhealpix.PackCellID128(sub.FacetID, 0, nil)
		if err != nil {
			return nil, err
		}

		planarGeom := projectGeometryToPlanar(el, sub.FacetID, sub.Geometry)
		if planarGeom == nil {
			continue
		}

		prep := NewPreparedGeometry(planarGeom)
		rawCells = WalkPlanarCell128(rootCell, rootBox, prep, targetRes, rawCells)
	}

	compacted, err := rhealpix.Compact(rawCells)
	if err != nil {
		return nil, fmt.Errorf("failed compacting decomposed cells: %w", err)
	}

	return compacted, nil
}
