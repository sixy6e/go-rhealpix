package rhealpixorb

import (
	"github.com/paulmach/orb"
	"github.com/paulmach/orb/planar"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// ============================================================================
// QUERY RANGE DECOMPOSITION (Bounding Boxes -> 1D Subtree Query Ranges)
// ============================================================================

// decomposeCell recursively decomposes an rHEALPix tree branch against a query BBox.
// Fully enclosed cells emit 1D [Min, Max] SubtreeRanges immediately, short-circuiting recursion.
func decomposeCell(
	el *rhealpix.Ellipsoid,
	cell rhealpix.CellID64,
	bound orb.Bound,
	targetRes uint8,
	ranges *[]Uint64Range,
) {
	currentRes := cell.Resolution()

	cellBound, err := Cell64ToBound(el, cell)
	if err != nil {
		return
	}

	// discard if cell does not intersect target bounding box
	if !boundsIntersect(cellBound, bound) {
		return
	}

	// MAX DELTA: only short-circuit on full containment if within 2 levels of targetRes.
	// coarser cells (e.g. Res 5 when target is 12) MUST keep decomposing.
	maxDelta := uint8(2)
	isNearTarget := (targetRes >= maxDelta) && (currentRes >= targetRes-maxDelta)

	// SHORT-CIRCUIT: reached targetRes OR (cell is fully enclosed AND near targetRes) -> Emit Range & STOP
	// don't want any infinite recursion or even getting into nanometre precision,
	// nor millions of children
	if currentRes >= targetRes || (boundsContains(bound, cellBound) && isNearTarget) {
		minCell, maxCell := cell.SubtreeRange(targetRes)
		*ranges = append(*ranges, Uint64Range{
			Min: uint64(minCell),
			Max: uint64(maxCell),
		})
		return // CRITICAL: stop recursing down this branch!
	}

	// straddles boundary -> recurse into 9 sub-cells
	for digit := uint8(0); digit < 9; digit++ {
		childCell, err := cell.Child(digit)
		if err != nil {
			continue
		}
		decomposeCell(el, childCell, bound, targetRes, ranges)
	}
}

// decomposeCell128 recursively decomposes an rHEALPix tree branch against a query BBox.
// Fully enclosed cells emit 1D [Min, Max] SubtreeRanges immediately, short-circuiting recursion.
func decomposeCell128(
	el *rhealpix.Ellipsoid,
	cell rhealpix.CellID128,
	bound orb.Bound,
	targetRes uint8,
	ranges *[]Uint128Range,
) {
	currentRes := cell.Resolution()

	cellBound, err := Cell128ToBound(el, cell)
	if err != nil {
		return
	}

	// discard if cell does not intersect target bounding box
	if !boundsIntersect(cellBound, bound) {
		return
	}

	// MAX DELTA: only short-circuit on full containment if within 2 levels of targetRes.
	maxDelta := uint8(2)
	isNearTarget := (targetRes >= maxDelta) && (currentRes >= targetRes-maxDelta)

	// SHORT-CIRCUIT: reached targetRes OR (cell is fully enclosed AND near targetRes) -> Emit Range & STOP
	// don't want any infinite recursion or even getting into nanometre precision,
	// nor millions of children
	if currentRes >= targetRes || (boundsContains(bound, cellBound) && isNearTarget) {
		minCell, maxCell := cell.SubtreeRange(targetRes)
		*ranges = append(*ranges, NewUint128Range(minCell, maxCell))
		return // CRITICAL: stop recursing down this branch!
	}

	// straddles boundary -> recurse into 9 sub-cells
	for digit := uint8(0); digit < 9; digit++ {
		childCell, err := cell.Child(digit)
		if err != nil {
			continue
		}
		decomposeCell128(el, childCell, bound, targetRes, ranges)
	}
}

// ============================================================================
// INGESTION FOOTPRINT DECOMPOSITION (Polygons -> Multi-Res Array Keys)
// ============================================================================

// decomposePolygon recursively decomposes an rHEALPix tree branch against a STAC polygon footprint.
// Or any polygon for that matter ...
// Fully enclosed coarse cells are retained directly in rawCells without generating leaf children.
// TODO; probably not needed anymore, decomposeCell has reworked the workflow.
// If needed, will need to rework it.
func decomposePolygon(
	el *rhealpix.Ellipsoid,
	cell rhealpix.CellID64,
	poly orb.Polygon,
	polyBound orb.Bound,
	targetRes uint8,
	rawCells *[]rhealpix.CellID64,
) {
	currentRes := cell.Resolution()

	cellBound, err := Cell64ToBound(el, cell)
	if err != nil {
		return
	}

	// quick discard against polygon bounding box
	if !boundsIntersect(cellBound, polyBound) {
		return
	}

	cellCenter := orb.Point{
		(cellBound.Min.X() + cellBound.Max.X()) / 2.0,
		(cellBound.Min.Y() + cellBound.Max.Y()) / 2.0,
	}

	// base case: reached target resolution
	if currentRes >= targetRes {
		if planar.PolygonContains(poly, cellCenter) {
			*rawCells = append(*rawCells, cell)
		}
		return
	}

	// SHORT-CIRCUIT: coarse parent fully inside interior -> EMIT PARENT & STOP
	if boundsContains(polyBound, cellBound) && planar.PolygonContains(poly, cellCenter) {
		*rawCells = append(*rawCells, cell)
		return // stops recursing into millions of children!
	}

	// straddles boundary -> recurse into 9 sub-cells
	// facet, res, path, err := rhealpix.DecodeCellID64(cell)
	// if err != nil {
	// 	return
	// }

	// straddles boundary -> recurse into 9 sub-cells
	for digit := uint8(0); digit < 9; digit++ {
		// childPath := append(append([]uint8(nil), path...), digit)
		// childCell, err := rhealpix.PackCellID64(facet, res+1, childPath)
		childCell, err := cell.Child(digit)
		if err != nil {
			continue
		}
		decomposePolygon(el, childCell, poly, polyBound, targetRes, rawCells)
	}
}

// decomposePolygon128 recursively decomposes an rHEALPix tree branch against a polygon footprint.
func decomposePolygon128(
	el *rhealpix.Ellipsoid,
	cell rhealpix.CellID128,
	poly orb.Polygon,
	polyBound orb.Bound,
	targetRes uint8,
	rawCells *[]rhealpix.CellID128,
) {
	currentRes := cell.Resolution()

	cellBound, err := Cell128ToBound(el, cell)
	if err != nil {
		return
	}

	// quick discard against polygon bounding box
	if !boundsIntersect(cellBound, polyBound) {
		return
	}

	cellCenter := orb.Point{
		(cellBound.Min.X() + cellBound.Max.X()) / 2.0,
		(cellBound.Min.Y() + cellBound.Max.Y()) / 2.0,
	}

	// base case: reached target resolution
	if currentRes >= targetRes {
		if planar.PolygonContains(poly, cellCenter) {
			*rawCells = append(*rawCells, cell)
		}
		return
	}

	// SHORT-CIRCUIT: coarse parent fully inside interior -> EMIT PARENT & STOP
	if boundsContains(polyBound, cellBound) && planar.PolygonContains(poly, cellCenter) {
		*rawCells = append(*rawCells, cell)
		return // stops recursing into millions of children!
	}

	// straddles boundary -> recurse into 9 sub-cells
	for digit := uint8(0); digit < 9; digit++ {
		childCell, err := cell.Child(digit)
		if err != nil {
			continue
		}
		decomposePolygon128(el, childCell, poly, polyBound, targetRes, rawCells)
	}
}

// ============================================================================
// SPATIAL PREDICATE HELPERS
// ============================================================================

// boundsIntersect returns true if axis-aligned bounds b1 and b2 overlap.
func boundsIntersect(b1, b2 orb.Bound) bool {
	return !(b1.Max.X() < b2.Min.X() || b1.Min.X() > b2.Max.X() ||
		b1.Max.Y() < b2.Min.Y() || b1.Min.Y() > b2.Max.Y())
}

// boundsContains returns true if outer bound completely encloses inner bound.
func boundsContains(outer, inner orb.Bound) bool {
	return outer.Min.X() <= inner.Min.X() && outer.Max.X() >= inner.Max.X() &&
		outer.Min.Y() <= inner.Min.Y() && outer.Max.Y() >= inner.Max.Y()
}

// isCellFullyInsidePolygon tests if all 4 corners and the center of cellBound sit strictly inside poly.
func isCellFullyInsidePolygon(poly orb.Polygon, cellBound orb.Bound) bool {
	corners := []orb.Point{
		cellBound.Min,
		{cellBound.Max.X(), cellBound.Min.Y()},
		cellBound.Max,
		{cellBound.Min.X(), cellBound.Max.Y()},
		{(cellBound.Min.X() + cellBound.Max.X()) / 2.0, (cellBound.Min.Y() + cellBound.Max.Y()) / 2.0},
	}

	for _, pt := range corners {
		if !planar.PolygonContains(poly, pt) {
			return false
		}
	}

	// ensure no polygon boundary segment crosses through the cell
	return !intersectsPolygonRing(poly, cellBound)
}

// intersectsPolygonRing checks if any exterior ring segment of poly intersects cellBound.
func intersectsPolygonRing(poly orb.Polygon, cellBound orb.Bound) bool {
	if len(poly) == 0 {
		return false
	}

	outerRing := poly[0]
	cellPolygon := orb.Polygon{
		orb.Ring{
			cellBound.Min,
			{cellBound.Max.X(), cellBound.Min.Y()},
			cellBound.Max,
			{cellBound.Min.X(), cellBound.Max.Y()},
			cellBound.Min,
		},
	}

	// check if any vertex of the ring lies inside the cell
	for _, pt := range outerRing {
		if cellBound.Min.X() <= pt.X() && pt.X() <= cellBound.Max.X() &&
			cellBound.Min.Y() <= pt.Y() && pt.Y() <= cellBound.Max.Y() {
			return true
		}
	}

	// quick centre distance check for boundary overlap
	cellCenter := orb.Point{
		(cellBound.Min.X() + cellBound.Max.X()) / 2.0,
		(cellBound.Min.Y() + cellBound.Max.Y()) / 2.0,
	}
	return planar.PolygonContains(cellPolygon, cellCenter)
}
