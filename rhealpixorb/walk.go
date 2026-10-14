package rhealpixorb

import (
	"fmt"
	"math"

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

// Subdivide splits a PlanarBox into an N x N grid of child boxes matching
// the rHEALPix child ordering (0 .. N^2-1) in Top-Left (UL) space.
func (pb PlanarBox) Subdivide(nSide uint8) []PlanarBox {
	n := int(nSide)
	totalChildren := n * n
	children := make([]PlanarBox, totalChildren)

	stepX := (pb.MaxX - pb.MinX) / float64(n)
	stepY := (pb.MaxY - pb.MinY) / float64(n)

	idx := 0
	for row := 0; row < n; row++ { // Top-to-Bottom (Row ordering)
		y0 := pb.MinY + float64(row)*stepY
		y1 := pb.MinY + float64(row+1)*stepY

		for col := 0; col < n; col++ { // Left-to-Right (Column ordering)
			x0 := pb.MinX + float64(col)*stepX
			x1 := pb.MinX + float64(col+1)*stepX

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
			if len(ring) < 3 {
				continue
			}

			// ensure ring is closed for segment generation
			closedRing := ring
			if ring[0] != ring[len(ring)-1] {
				closedRing = make(orb.Ring, len(ring)+1)
				copy(closedRing, ring)
				closedRing[len(ring)] = ring[0]
			}

			for i := 0; i < len(closedRing)-1; i++ {
				p1, p2 := closedRing[i], closedRing[i+1]

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
func projectGeometryToPlanar(d *rhealpix.RhpDGGS, facetID uint8, geom orb.Geometry) orb.Geometry {
	if geom == nil {
		return nil
	}

	switch g := geom.(type) {
	case orb.Polygon:
		proj := projectPolygonToPlanar(d, facetID, g)
		return ensurePolygonCCW(proj)
	case orb.MultiPolygon:
		var mp orb.MultiPolygon
		for _, poly := range g {
			proj := projectPolygonToPlanar(d, facetID, poly)
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

// projectPolygonToPlanar maps a WGS84 polygon to normalised [0, 1] x [0, 1] facet local planar space.
func projectPolygonToPlanar(d *rhealpix.RhpDGGS, facetID uint8, poly orb.Polygon) orb.Polygon {
	// densify polar cap facets (0=N, 5=S) to capture non-linear parallel curvature and 0° seam apexes
	if facetID == 0 || facetID == 5 {
		densifiedPoly := make(orb.Polygon, len(poly))
		for i, ring := range poly {
			densifiedPoly[i] = DensifyWGS84Ring(ring, 16)
		}
		poly = densifiedPoly
	}

	projected := make(orb.Polygon, len(poly))

	for i, ring := range poly {
		projectedRing := make(orb.Ring, len(ring))
		for j, pt := range ring {
			lon, lat := pt[0], pt[1]

			xLoc, yLoc, err := d.LonLatToPlanar(facetID, lon, lat)
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

// ForwardProjectPolygon projects geographic points directly through RhpDGGS.
// TODO: cleanup; this is replaced with ForwardProjectPolygonForcedRegion
// to enforce topological correctness for geometry structures.
func ForwardProjectPolygon(d *rhealpix.RhpDGGS, poly orb.Polygon) orb.Polygon {
	if len(poly) == 0 {
		return nil
	}

	projected := make(orb.Polygon, len(poly))

	for i, ring := range poly {

		// project directly via standard 3-stage ForwardProjectRadians
		projectedRing := make(orb.Ring, 0, len(ring))
		for _, pt := range ring {
			lonRad := pt.X() * rhealpix.DegToRad
			latRad := pt.Y() * rhealpix.DegToRad

			// pure projection call: HealpixSphere -> Seam Lock -> combineTriangles
			xRhp, yRhp, _ := d.ForwardProjectRadians(lonRad, latRad)
			projectedRing = append(projectedRing, orb.Point{xRhp, yRhp})
		}

		projected[i] = projectedRing
	}

	return projected
}

// ForwardProjectPolygonForcedRegion projects geographic polygon points through RhpDGGS
// using forced regional context, maintaining seam alignment across all regions.
func ForwardProjectPolygonForcedRegion(d *rhealpix.RhpDGGS, poly orb.Polygon, dom rhealpix.GeographicFacetDomain) orb.Polygon {
	if len(poly) == 0 {
		return nil
	}

	projected := make(orb.Polygon, len(poly))

	for i, ring := range poly {
		projectedRing := make(orb.Ring, 0, len(ring))
		for _, pt := range ring {
			lonDeg := pt.X()
			latRad := pt.Y() * rhealpix.DegToRad
			isOnSeam := math.Abs(math.Abs(lonDeg)-180.0) < 1e-11

			var lonRad float64
			if isOnSeam {
				// lock input longitude sign based on domain context
				if dom.FacetID == 4 || dom.MaxLon == 180.0 {
					lonRad = math.Pi // Eastern sector (Sector R) -> +180°
				} else if dom.FacetID == 1 || dom.MinLon == 180.0 || dom.MinLon == -270.0 {
					lonRad = -math.Pi // Western / Extended Wrap sector (Sector O) -> -180°
				} else {
					if lonDeg >= 0 {
						lonRad = math.Pi
					} else {
						lonRad = -math.Pi
					}
				}
			} else {
				lonRad = lonDeg * rhealpix.DegToRad
			}

			// forward project forcing the region a point is to be placed within
			xRhp, yRhp, _ := d.ForwardProjectForcedRegion(lonRad, latRad, dom.Region)

			// HARD LOCK: enforce planar X_RHP alignment ONLY on Equatorial Facets (1 and 4)
			if isOnSeam {
				if dom.FacetID == 4 {
					xRhp = math.Pi
				} else if dom.FacetID == 1 {
					xRhp = -math.Pi
				}
			}

			projectedRing = append(projectedRing, orb.Point{xRhp, yRhp})
		}
		projected[i] = projectedRing
	}

	return projected
}

func PartitionPolygonToCells(d *rhealpix.RhpDGGS, inputGeom orb.Geometry, targetRes uint8) ([]rhealpix.CellID64, error) {
	facetSubGeoms, err := ClipGeometryToFacets(d, inputGeom)
	if err != nil {
		return nil, err
	}

	rawCells := make([]rhealpix.CellID64, 0, 128)
	rootBox := PlanarBox{MinX: 0.0, MaxX: 1.0, MinY: 0.0, MaxY: 1.0}

	for _, sub := range facetSubGeoms {

		rootCell, err := rhealpix.PackCellID64(sub.FacetID, 0, nil)
		if err != nil {
			return nil, err
		}

		// transform global planar radians -> Local [0, 1] Facet UV
		localGeom := PlanarToLocalUVGeometry(d, sub.FacetID, sub.Geometry)

		// prepare local geometry & perform quadtree walk inside [0, 1] x [0, 1]
		prep := NewPreparedGeometry(localGeom)
		rawCells = WalkPlanarCellGeneric(d, rootCell, rootBox, prep, targetRes, rawCells)
	}

	return rawCells, nil
}

func PartitionPolygonToCells128(d *rhealpix.RhpDGGS, inputGeom orb.Geometry, targetRes uint8) ([]rhealpix.CellID128, error) {
	facetSubGeoms, err := ClipGeometryToFacets(d, inputGeom)
	if err != nil {
		return nil, err
	}

	rawCells := make([]rhealpix.CellID128, 0, 128)
	rootBox := PlanarBox{MinX: 0.0, MaxX: 1.0, MinY: 0.0, MaxY: 1.0}

	for _, sub := range facetSubGeoms {

		rootCell, err := rhealpix.PackCellID128(sub.FacetID, 0, nil)
		if err != nil {
			return nil, err
		}

		// transform global planar radians -> Local [0, 1] Facet UV
		localGeom := PlanarToLocalUVGeometry(d, sub.FacetID, sub.Geometry)

		// prepare local geometry & perform quadtree walk inside [0, 1] x [0, 1]
		prep := NewPreparedGeometry(localGeom)
		rawCells = WalkPlanarCellGeneric128(d, rootCell, rootBox, prep, targetRes, rawCells)
	}

	return rawCells, nil
}

// PlanarToLocalUVGeometry transforms an orb.Geometry from global planar radians [x, y]
// into local unit square coordinates [0, 1] x [0, 1] for a given base facet.
func PlanarToLocalUVGeometry(d *rhealpix.RhpDGGS, facetID uint8, geom orb.Geometry) orb.Geometry {
	if geom == nil {
		return nil
	}

	switch g := geom.(type) {
	case orb.Polygon:
		return transformPolygonUV(d, facetID, g)

	case orb.MultiPolygon:
		mp := make(orb.MultiPolygon, len(g))
		for i, poly := range g {
			mp[i] = transformPolygonUV(d, facetID, poly)
		}
		return mp

	default:
		return nil
	}
}

func transformPolygonUV(d *rhealpix.RhpDGGS, facetID uint8, poly orb.Polygon) orb.Polygon {
	out := make(orb.Polygon, len(poly))
	for r, ring := range poly {
		localRing := make(orb.Ring, len(ring))
		for i, pt := range ring {
			// pt[0] = xRhp, pt[1] = yRhp in authalic planar radians
			u, v := d.RadiansToFacetLocal(facetID, pt[0], pt[1])
			localRing[i] = orb.Point{u, v}
		}
		out[r] = localRing
	}
	return out
}

// WalkPlanarCell recursively walks the rHEALPix nonary tree in local planar space.
// Similar methodology to Uber's H3 PolyFill. But simpler as we have squares, not overlapping
// hexagons and pentagons.
// Only works for 3x3 sided case.
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

// WalkPlanarCellGeneric recursively walks the rHEALPix nonary tree in local planar space.
// Similar methodology to Uber's H3 PolyFill. But simpler as we have squares, not overlapping
// hexagons and pentagons.
// Works for generic n-sided case, eg 2x2 or 3x3.
func WalkPlanarCellGeneric(
	d *rhealpix.RhpDGGS,
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

		// dynamically subdivide into N x N children
		childBoxes := box.Subdivide(d.Nside())
		numChildren := uint8(len(childBoxes))

		for i := uint8(0); i < numChildren; i++ {
			childCell, err := cell.Child(i)
			if err != nil {
				continue
			}

			results = WalkPlanarCellGeneric(d, childCell, childBoxes[i], prep, targetRes, results)
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
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8,
) ([]rhealpix.CellID64, error) {
	rawCells, err := PartitionPolygonToCells(d, geom, targetRes)
	if err != nil {
		return nil, fmt.Errorf("failed PartitionPolygonToCells: %w", err)
	}

	// fmt.Println("len(rawCells): ", len(rawCells))
	compacted, err := rhealpix.Compact(rawCells)
	if err != nil {
		return nil, fmt.Errorf("failed compacting decomposed cells: %w", err)
	}

	return compacted, nil
}

// WalkPlanarCell128 recursively walks the rHEALPix nonary tree in local planar space using 128-bit cell IDs.
// Similar methodology to Uber's H3 PolyFill. But simpler as we have squares, not overlapping
// hexagons and pentagons.
// Only works for 3x3 sided case.
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

// WalkPlanarCellGeneric128 recursively walks the rHEALPix nonary tree in local planar space using 128-bit cell IDs.
// Similar methodology to Uber's H3 PolyFill. But simpler as we have squares, not overlapping
// hexagons and pentagons.
// Works for generic n-sided case, eg 2x2 or 3x3.
func WalkPlanarCellGeneric128(
	d *rhealpix.RhpDGGS,
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

		// dynamically subdivide into N x N children
		childBoxes := box.Subdivide(d.Nside())
		numChildren := uint8(len(childBoxes))

		for i := uint8(0); i < numChildren; i++ {
			childCell, err := cell.Child(i)
			if err != nil {
				continue
			}

			results = WalkPlanarCellGeneric128(d, childCell, childBoxes[i], prep, targetRes, results)
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
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8,
) ([]rhealpix.CellID128, error) {
	rawCells, err := PartitionPolygonToCells128(d, geom, targetRes)
	if err != nil {
		return nil, fmt.Errorf("failed PartitionPolygonToCells: %w", err)
	}

	// fmt.Println("len(rawCells): ", len(rawCells))
	compacted, err := rhealpix.Compact(rawCells)
	if err != nil {
		return nil, fmt.Errorf("failed compacting decomposed cells: %w", err)
	}

	return compacted, nil
}
