package rhealpixorb

import (
	"fmt"
	"math"
	"sort"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/geojson"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// Default maxGap parameter to bridge small space-filling curve gaps when querying
const DefaultQueryMaxGap uint64 = 50000000

// Uint64Range represents a 1D min/max boundary for database queries.
type Uint64Range struct {
	Min uint64
	Max uint64
}

// Uint128Range represents a compound 1D min/max boundary for 128-bit database queries
// across two uint64 dimensions ("cell_high" and "cell_low").
type Uint128Range struct {
	MinHigh uint64
	MinLow  uint64
	MaxHigh uint64
	MaxLow  uint64
}

// NewUint128Range constructs a Uint128Range directly from Min and Max CellID128 values.
func NewUint128Range(minCell, maxCell rhealpix.CellID128) Uint128Range {
	return Uint128Range{
		MinHigh: minCell.High,
		MinLow:  minCell.Low,
		MaxHigh: maxCell.High,
		MaxLow:  maxCell.Low,
	}
}

// MinCell returns the lower bound as a native CellID128.
func (r Uint128Range) MinCell() rhealpix.CellID128 {
	return rhealpix.CellID128{High: r.MinHigh, Low: r.MinLow}
}

// MaxCell returns the upper bound as a native CellID128.
func (r Uint128Range) MaxCell() rhealpix.CellID128 {
	return rhealpix.CellID128{High: r.MaxHigh, Low: r.MaxLow}
}

// BoundingBoxToTileDBRanges translates a WGS84 geographic bounding box
// into a set of 1D uint64 ranges using top down spatial tree decomposition.
// Safely handles bounding boxes that cross the 180° Antimeridian.
func BoundingBoxToTileDBRanges(el *rhealpix.Ellipsoid, bound orb.Bound, targetRes uint8) ([]Uint64Range, error) {
	// handle antimeridian crossing (Min.X > Max.X) by splitting into two bounds (East and West boxes)
	if bound.Min.X() > bound.Max.X() {
		b1 := orb.Bound{
			Min: orb.Point{bound.Min.X(), bound.Min.Y()},
			Max: orb.Point{180.0, bound.Max.Y()},
		}
		b2 := orb.Bound{
			Min: orb.Point{-180.0, bound.Min.Y()},
			Max: orb.Point{bound.Max.X(), bound.Max.Y()},
		}

		r1, err := BoundingBoxToTileDBRanges(el, b1, targetRes)
		if err != nil {
			return nil, err
		}
		r2, err := BoundingBoxToTileDBRanges(el, b2, targetRes)
		if err != nil {
			return nil, err
		}
		return MergeRangesWithGap(append(r1, r2...), DefaultQueryMaxGap), nil
	}

	var ranges []Uint64Range

	// decompose across base facets (0..5)
	for facet := uint8(0); facet < 6; facet++ {
		rootCell, err := rhealpix.PackCellID64(facet, 0, nil)
		if err != nil {
			continue
		}

		// check base facet envelope first; skip if non-overlapping
		rootBound, err := Cell64ToBound(el, rootCell)
		if err == nil && !boundsIntersect(rootBound, bound) {
			continue
		}

		decomposeCell(el, rootCell, bound, targetRes, &ranges)
	}

	if len(ranges) == 0 {
		return nil, nil
	}

	// consolidate adjacent or near-adjacent intervals
	return MergeRangesWithGap(ranges, DefaultQueryMaxGap), nil
}

// BoundingBoxToTileDBRanges128 translates a WGS84 geographic bounding box
// into a set of 128-bit compound ranges suitable for high-resolution TileDB queries.
// Safely handles bounding boxes that cross the 180° Antimeridian.
func BoundingBoxToTileDBRanges128(el *rhealpix.Ellipsoid, bound orb.Bound, targetRes uint8) ([]Uint128Range, error) {
	// handle antimeridian crossing (Min.X > Max.X) by splitting into two bounds (East and West boxes)
	if bound.Min.X() > bound.Max.X() {
		b1 := orb.Bound{
			Min: orb.Point{bound.Min.X(), bound.Min.Y()},
			Max: orb.Point{180.0, bound.Max.Y()},
		}
		b2 := orb.Bound{
			Min: orb.Point{-180.0, bound.Min.Y()},
			Max: orb.Point{bound.Max.X(), bound.Max.Y()},
		}

		r1, err := BoundingBoxToTileDBRanges128(el, b1, targetRes)
		if err != nil {
			return nil, err
		}
		r2, err := BoundingBoxToTileDBRanges128(el, b2, targetRes)
		if err != nil {
			return nil, err
		}
		return MergeRanges128(append(r1, r2...)), nil
	}

	// extract 4 corners of the bounding box
	corners := []orb.Point{
		bound.Min,
		{bound.Max.X(), bound.Min.Y()},
		bound.Max,
		{bound.Min.X(), bound.Max.Y()},
	}

	// sample dense points along the edges
	sampledPoints := densifyBoundEdges(corners, 10)

	// forward project all points to CellID128 at targetRes
	cellMap := make(map[rhealpix.CellID128]bool)
	for _, pt := range sampledPoints {
		cell, err := rhealpix.ForwardTransform128(el, pt.X(), pt.Y(), targetRes)
		if err != nil {
			return nil, fmt.Errorf("failed forward project 128 (%f, %f): %w", pt.X(), pt.Y(), err)
		}
		cellMap[cell] = true
	}

	// extract slice for compacting
	rawCells := make([]rhealpix.CellID128, 0, len(cellMap))
	for c := range cellMap {
		rawCells = append(rawCells, c)
	}

	// compact the cell set using generic Compact[CellID128]
	compacted, err := rhealpix.Compact(rawCells)
	if err != nil {
		compacted = rawCells
	}

	// convert each compacted cell into a 128-bit SubtreeRange [Min, Max]
	ranges := make([]Uint128Range, 0, len(compacted))
	for _, c := range compacted {
		minCell, maxCell := c.SubtreeRange(targetRes)
		ranges = append(ranges, NewUint128Range(minCell, maxCell))
	}

	return MergeRanges128(ranges), nil
}

// KRingToTileDBRanges calculates the K-Ring neighbours around an origin cell at resolution R
// and converts them into a compacted set of 1D uint64 ranges for TileDB.
func KRingToTileDBRanges(originCell rhealpix.CellID64, k int) ([]Uint64Range, error) {
	// generate K-Ring cells
	cells, err := rhealpix.KRing(originCell, k)
	if err != nil {
		return nil, fmt.Errorf("failed to compute k-ring: %w", err)
	}

	// compact the cell set to minimise range overhead
	compacted, err := rhealpix.Compact(cells)
	if err != nil {
		return nil, fmt.Errorf("compaction failed on k-ring cells: %w", err)
	}

	// convert each cell to SubtreeRange [Min, Max]
	targetRes := originCell.Resolution()
	ranges := make([]Uint64Range, 0, len(compacted))
	for _, c := range compacted {
		minCell, maxCell := c.SubtreeRange(targetRes)
		ranges = append(ranges, Uint64Range{
			Min: uint64(minCell),
			Max: uint64(maxCell),
		})
	}

	return MergeRangesWithGap(ranges, DefaultQueryMaxGap), nil
}

// KRingToTileDBRanges128 calculates the K-Ring neighbours around a 128-bit origin cell
// and converts them into a compacted set of compound ranges for high-resolution TileDB queries.
func KRingToTileDBRanges128(originCell rhealpix.CellID128, k int) ([]Uint128Range, error) {
	// generate K-Ring cells
	cells, err := rhealpix.KRing(originCell, k)
	if err != nil {
		return nil, fmt.Errorf("failed to compute k-ring 128: %w", err)
	}

	// compact the cell set
	compacted, err := rhealpix.Compact(cells)
	if err != nil {
		return nil, fmt.Errorf("compaction failed on k-ring 128 cells: %w", err)
	}

	// convert each cell to 128-bit SubtreeRange [Min, Max]
	targetRes := originCell.Resolution()
	ranges := make([]Uint128Range, 0, len(compacted))
	for _, c := range compacted {
		minCell, maxCell := c.SubtreeRange(targetRes)
		ranges = append(ranges, NewUint128Range(minCell, maxCell))
	}

	return MergeRanges128(ranges), nil
}

// MergeRanges sorts and merges strictly overlapping or contiguous 64-bit uint64 ranges (maxGap = 1).
func MergeRanges(ranges []Uint64Range) []Uint64Range {
	return MergeRangesWithGap(ranges, 1)
}

// MergeRangesWithGap sorts and merges overlapping or near-adjacent Uint64Ranges.
// If the gap between range[i].Max and range[i+1].Min is <= maxGap, they are combined,
// provided both ranges share the exact same top-level SUID geographic trunk.
func MergeRangesWithGap(ranges []Uint64Range, maxGap uint64) []Uint64Range {
	if len(ranges) <= 1 {
		return ranges
	}

	// sort ranges by Min bound
	sort.Slice(ranges, func(i, j int) bool {
		return ranges[i].Min < ranges[j].Min
	})

	// merge overlapping or adjacent ranges
	merged := make([]Uint64Range, 0, len(ranges))
	current := ranges[0]

	for i := 1; i < len(ranges); i++ {
		next := ranges[i]

		// extract top 12 bits (Facet [3 bits] + Res Header [5 bits] + Level 1 Digit [4 bits])
		// shift by 52 (64 - 12 = 52) to isolate the SUID trunk.
		currentTrunk := current.Min >> 52
		nextTrunk := next.Min >> 52

		// safe check against uint64 overflow when adding maxGap
		canReach := false
		if current.Max <= math.MaxUint64-maxGap {
			canReach = (next.Min <= current.Max+maxGap)
		} else {
			canReach = (next.Min <= current.Max)
		}

		// merge ONLY if ranges belong to the SAME SUID trunk AND fall within maxGap
		if currentTrunk == nextTrunk && canReach {
			if next.Max > current.Max {
				current.Max = next.Max
			}
		} else {
			merged = append(merged, current)
			current = next
		}
	}

	merged = append(merged, current)
	return merged
}

// compare128 returns -1 if a < b, 0 if a == b, and 1 if a > b.
func compare128(aHigh, aLow, bHigh, bLow uint64) int {
	if aHigh < bHigh {
		return -1
	}
	if aHigh > bHigh {
		return 1
	}
	if aLow < bLow {
		return -1
	}
	if aLow > bLow {
		return 1
	}
	return 0
}

// add128 adds maxGap to a 128-bit integer (high, low), properly handling carry across words.
func add128(high, low, gap uint64) (outHigh, outLow uint64) {
	outLow = low + gap
	outHigh = high
	if outLow < low { // Overflow occurred in Low word
		outHigh++
	}
	return outHigh, outLow
}

// MergeRanges128 sorts and merges strictly overlapping or contiguous 128-bit ranges (maxGap = 1).
func MergeRanges128(ranges []Uint128Range) []Uint128Range {
	return MergeRangesWithGap128(ranges, 1)
}

// MergeRangesWithGap128 sorts and merges overlapping or near-adjacent Uint128Ranges.
// Combines ranges if the gap between current.Max and next.Min is <= maxGap,
// provided both ranges share the exact same top-level SUID geographic trunk and resolution header.
func MergeRangesWithGap128(ranges []Uint128Range, maxGap uint64) []Uint128Range {
	if len(ranges) <= 1 {
		return ranges
	}

	// sort ranges by 128-bit Min bound (MinHigh, then MinLow)
	sort.Slice(ranges, func(i, j int) bool {
		if ranges[i].MinHigh != ranges[j].MinHigh {
			return ranges[i].MinHigh < ranges[j].MinHigh
		}
		return ranges[i].MinLow < ranges[j].MinLow
	})

	// merge overlapping or adjacent ranges
	merged := make([]Uint128Range, 0, len(ranges))
	current := ranges[0]

	for i := 1; i < len(ranges); i++ {
		next := ranges[i]

		// extract top 12 bits of High word (Facet [3 bits] + Res Header [5 bits] + Level 1 Digit [4 bits])
		// shift by 52 (64 - 12 = 52) to isolate SUID trunk + resolution header level.
		currentTrunk := current.MaxHigh >> 52
		nextTrunk := next.MinHigh >> 52

		// compute (current.Max + maxGap) with 128-bit overflow safety
		reachHigh := current.MaxHigh
		reachLow := current.MaxLow

		if math.MaxUint64-current.MaxLow < maxGap {
			// Low word overflows into High word
			if current.MaxHigh < math.MaxUint64 {
				reachHigh = current.MaxHigh + 1
				reachLow = maxGap - (math.MaxUint64 - current.MaxLow) - 1
			} else {
				// cap at absolute 128-bit limit
				reachLow = math.MaxUint64
			}
		} else {
			reachLow = current.MaxLow + maxGap
		}

		// check if next.Min <= (current.Max + maxGap) in 128-bit comparison
		canReach := (next.MinHigh < reachHigh) ||
			(next.MinHigh == reachHigh && next.MinLow <= reachLow)

		// merge ONLY if ranges share the exact same SUID trunk & resolution level AND fall within maxGap
		if currentTrunk == nextTrunk && canReach {
			// update current.Max to max(current.Max, next.Max)
			if next.MaxHigh > current.MaxHigh || (next.MaxHigh == current.MaxHigh && next.MaxLow > current.MaxLow) {
				current.MaxHigh = next.MaxHigh
				current.MaxLow = next.MaxLow
			}
		} else {
			merged = append(merged, current)
			current = next
		}
	}

	merged = append(merged, current)
	return merged
}

// STACGeometryToTileDBRanges translates a STAC GeoJSON Geometry into a solid,
// hierarchically compacted set of 1D TileDB ranges using top down planar decomposition.
// The approach of top down planar decomposition and compaction, is similar to Uber's H3 PolyFill.
func STACGeometryToTileDBRanges(
	el *rhealpix.Ellipsoid,
	geom geojson.Geometry,
	targetRes uint8,
) ([]Uint64Range, error) {

	g := geom.Geometry()
	if g == nil {
		return nil, fmt.Errorf("geometry is nil")
	}

	// run top down planar decomposition and compaction
	compactedCells, err := STACGeometryToTileDBRangesTopDown(el, g, targetRes)
	if err != nil {
		return nil, fmt.Errorf("failed top down geometry decomposition: %w", err)
	}

	if len(compactedCells) == 0 {
		return nil, nil
	}

	// convert compacted CellID64 elements into 1D Uint64Ranges
	ranges := make([]Uint64Range, 0, len(compactedCells))
	for _, c := range compactedCells {
		// minCell, maxCell := c.SubtreeRangeMax()
		minCell, maxCell := c.SubtreeRange(targetRes)
		ranges = append(ranges, Uint64Range{
			Min: uint64(minCell),
			Max: uint64(maxCell),
		})
	}

	// consolidate adjacent and near adjacent 1D ranges using DefaultQueryMaxGap
	return MergeRangesWithGap(ranges, DefaultQueryMaxGap), nil
}

// STACGeometryToTileDBRanges128 translates a STAC GeoJSON Geometry into a solid,
// hierarchically compacted set of 1D TileDB ranges using top down planar decomposition.
// The approach of top down planar decomposition and compaction, is similar to Uber's H3 PolyFill.
func STACGeometryToTileDBRanges128(
	el *rhealpix.Ellipsoid,
	geom geojson.Geometry,
	targetRes uint8,
) ([]Uint128Range, error) {

	g := geom.Geometry()
	if g == nil {
		return nil, fmt.Errorf("geometry is nil")
	}

	// run top down planar decomposition and compaction
	compactedCells, err := STACGeometryToTileDBRangesTopDown128(el, g, targetRes)
	if err != nil {
		return nil, fmt.Errorf("failed top down geometry decomposition: %w", err)
	}

	if len(compactedCells) == 0 {
		return nil, nil
	}

	// convert compacted CellID128 elements into 1D Uint128Ranges
	ranges := make([]Uint128Range, 0, len(compactedCells))
	for _, c := range compactedCells {
		// minCell, maxCell := c.SubtreeRangeMax()
		minCell, maxCell := c.SubtreeRange(targetRes)
		ranges = append(ranges, NewUint128Range(minCell, maxCell))
	}

	// consolidate adjacent and near adjacent 1D ranges using DefaultQueryMaxGap
	return MergeRangesWithGap128(ranges, DefaultQueryMaxGap), nil
}

// STACGeometryToCompactedCells decomposes a GeoJSON geometry into a solid,
// hierarchically compacted set of rHEALPix cell IDs down to targetRes using
// 64 bit cell construction.
func STACGeometryToCompactedCells(
	el *rhealpix.Ellipsoid,
	geom geojson.Geometry,
	targetRes uint8,
) ([]rhealpix.CellID64, error) {
	g := geom.Geometry()
	if g == nil {
		return nil, fmt.Errorf("geometry is nil")
	}

	// run top-down planar decomposition and compaction
	compactedCells, err := STACGeometryToTileDBRangesTopDown(el, g, targetRes)
	if err != nil {
		return nil, fmt.Errorf("failed top down geometry decomposition: %w", err)
	}

	return compactedCells, nil
}

// STACGeometryToCompactedCells128 decomposes a GeoJSON geometry into a solid,
// hierarchically compacted set of rHEALPix cell IDs down to targetRes using
// 128 bit cell construction.
func STACGeometryToCompactedCells128(
	el *rhealpix.Ellipsoid,
	geom geojson.Geometry,
	targetRes uint8,
) ([]rhealpix.CellID128, error) {
	g := geom.Geometry()
	if g == nil {
		return nil, fmt.Errorf("geometry is nil")
	}

	// run top-down planar decomposition and compaction
	compactedCells, err := STACGeometryToTileDBRangesTopDown128(el, g, targetRes)
	if err != nil {
		return nil, fmt.Errorf("failed top down geometry decomposition: %w", err)
	}

	return compactedCells, nil
}
