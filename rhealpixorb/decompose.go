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

// decomposePolygon recursively decomposes an rHEALPix tree branch against a polygon footprint.
// Or any polygon for that matter ...
// Fully enclosed coarse cells are retained directly in rawCells without generating leaf children.
// Note: STACGeometryToTileDBRangesTopDown in walk.go is preferred for ingestion pipelines
// as it performs exact planar [0, 1] x [0, 1] decomposition.
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

	// base case: reached target resolution
	if currentRes >= targetRes {
		cellCenter := orb.Point{
			(cellBound.Min.X() + cellBound.Max.X()) / 2.0,
			(cellBound.Min.Y() + cellBound.Max.Y()) / 2.0,
		}
		if planar.PolygonContains(poly, cellCenter) {
			*rawCells = append(*rawCells, cell)
		}
		return
	}

	// SHORT-CIRCUIT: coarse parent fully inside interior (strictly enclosed) -> EMIT PARENT & STOP
	if isCellFullyInsidePolygon(poly, cellBound) {
		*rawCells = append(*rawCells, cell)
		return // stops recursing into millions of children!
	}

	// straddles boundary -> recurse into 9 sub-cells
	for digit := uint8(0); digit < 9; digit++ {
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

	// base case: reached target resolution
	if currentRes >= targetRes {
		cellCenter := orb.Point{
			(cellBound.Min.X() + cellBound.Max.X()) / 2.0,
			(cellBound.Min.Y() + cellBound.Max.Y()) / 2.0,
		}
		if planar.PolygonContains(poly, cellCenter) {
			*rawCells = append(*rawCells, cell)
		}
		return
	}

	// SHORT-CIRCUIT: coarse parent fully inside interior (strictly enclosed) -> EMIT PARENT & STOP
	if isCellFullyInsidePolygon(poly, cellBound) {
		*rawCells = append(*rawCells, cell)
		return
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

// boundsIntersect returns true if axis-aligned bounds b1 and b2 overlap,
// handling unwrapped longitude ranges gracefully.
func boundsIntersect(b1, b2 orb.Bound) bool {
	b1MinX, b1MaxX := normaliseLonDeg(b1.Min.X()), normaliseLonDeg(b1.Max.X())
	b2MinX, b2MaxX := normaliseLonDeg(b2.Min.X()), normaliseLonDeg(b2.Max.X())

	// handle full 360-degree span
	if b1.Max.X()-b1.Min.X() >= 360.0 || b2.Max.X()-b2.Min.X() >= 360.0 {
		b1MinX, b1MaxX = -180.0, 180.0
		b2MinX, b2MaxX = -180.0, 180.0
	}

	latOverlap := !(b1.Max.Y() < b2.Min.Y() || b1.Min.Y() > b2.Max.Y())
	lonOverlap := !(b1MaxX < b2MinX || b1MinX > b2MaxX)

	return latOverlap && lonOverlap
}

// boundsContains returns true if outer bound completely encloses inner bound.
func boundsContains(outer, inner orb.Bound) bool {
	oMinX, oMaxX := normaliseLonDeg(outer.Min.X()), normaliseLonDeg(outer.Max.X())
	iMinX, iMaxX := normaliseLonDeg(inner.Min.X()), normaliseLonDeg(inner.Max.X())

	if outer.Max.X()-outer.Min.X() >= 360.0 {
		oMinX, oMaxX = -180.0, 180.0
	}

	latInside := outer.Min.Y() <= inner.Min.Y() && outer.Max.Y() >= inner.Max.Y()
	lonInside := oMinX <= iMinX && oMaxX >= iMaxX

	return latInside && lonInside
}

// isCellFullyInsidePolygon tests if all 4 corners and center of cellBound sit strictly inside poly,
// AND verifies that no polygon boundary segments cut through the cell.
func isCellFullyInsidePolygon(poly orb.Polygon, cellBound orb.Bound) bool {
	if len(poly) == 0 {
		return false
	}

	corners := []orb.Point{
		cellBound.Min,
		{cellBound.Max.X(), cellBound.Min.Y()},
		cellBound.Max,
		{cellBound.Min.X(), cellBound.Max.Y()},
		{(cellBound.Min.X() + cellBound.Max.X()) / 2.0, (cellBound.Min.Y() + cellBound.Max.Y()) / 2.0},
	}

	// ALL 5 sample points must be inside the outer shell (ring 0)
	for _, pt := range corners {
		if !planar.PolygonContains(orb.Polygon{poly[0]}, pt) {
			return false
		}
	}

	// NONE of the sample points may sit inside interior holes (rings 1..N)
	for i := 1; i < len(poly); i++ {
		holePoly := orb.Polygon{poly[i]}
		for _, pt := range corners {
			if planar.PolygonContains(holePoly, pt) {
				return false
			}
		}
	}

	// ensure no polygon ring segment crosses through the cell
	return !intersectsPolygonRing(poly, cellBound)
}

// intersectsPolygonRing checks if any ring segment of poly intersects cellBound.
func intersectsPolygonRing(poly orb.Polygon, cellBound orb.Bound) bool {
	if len(poly) == 0 {
		return false
	}

	minX, maxX := cellBound.Min.X(), cellBound.Max.X()
	minY, maxY := cellBound.Min.Y(), cellBound.Max.Y()

	for _, ring := range poly {
		for i := 0; i < len(ring)-1; i++ {
			p1, p2 := ring[i], ring[i+1]

			// AABB Segment overlap pre-filter
			pMinX, pMaxX := p1[0], p2[0]
			if pMinX > pMaxX {
				pMinX, pMaxX = pMaxX, pMinX
			}
			pMinY, pMaxY := p1[1], p2[1]
			if pMinY > pMaxY {
				pMinY, pMaxY = pMaxY, pMinY
			}

			if pMaxX < minX || pMinX > maxX || pMaxY < minY || pMinY > maxY {
				continue
			}

			// point inside box check
			if minX <= p1[0] && p1[0] <= maxX && minY <= p1[1] && p1[1] <= maxY {
				return true
			}
		}
	}

	return false
}
