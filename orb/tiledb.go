package rhealpixorb

import (
	"fmt"
	"sort"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

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

// BoundingBoxToTileDBRanges translates a WGS84 geographic bounding box
// into a set of 1D uint64 ranges suitable for TileDB range-filtering.
// Safely handles bounding boxes that cross the 180° Antimeridian.
func BoundingBoxToTileDBRanges(el *rhealpix.Ellipsoid, bound orb.Bound, targetRes uint8) ([]Uint64Range, error) {
	// handle antimeridian crossing (Min.X > Max.X) by splitting into two bounds
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
		return MergeRanges(append(r1, r2...)), nil
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

	// forward project all points to CellID64 at targetRes
	cellMap := make(map[rhealpix.CellID64]bool)
	for _, pt := range sampledPoints {
		cell, err := rhealpix.ForwardTransform64(el, pt.X(), pt.Y(), targetRes)
		if err != nil {
			return nil, fmt.Errorf("failed forward project (%f, %f): %w", pt.X(), pt.Y(), err)
		}
		cellMap[cell] = true
	}

	// extract slice for compacting
	rawCells := make([]rhealpix.CellID64, 0, len(cellMap))
	for c := range cellMap {
		rawCells = append(rawCells, c)
	}

	// compact the cell set
	compacted, err := rhealpix.Compact(rawCells)
	if err != nil {
		return nil, fmt.Errorf("compaction failed: %w", err)
	}

	// convert each compacted cell into a SubtreeRange [Min, Max]
	ranges := make([]Uint64Range, 0, len(compacted))
	for _, c := range compacted {
		minCell, maxCell := c.SubtreeRange()
		ranges = append(ranges, Uint64Range{
			Min: uint64(minCell),
			Max: uint64(maxCell),
		})
	}

	return MergeRanges(ranges), nil
}

// BoundingBoxToTileDBRanges128 translates a WGS84 geographic bounding box
// into a set of 128-bit compound ranges suitable for high-resolution TileDB queries.
// Safely handles bounding boxes that cross the 180° Antimeridian.
func BoundingBoxToTileDBRanges128(el *rhealpix.Ellipsoid, bound orb.Bound, targetRes uint8) ([]Uint128Range, error) {
	// handle antimeridian crossing (Min.X > Max.X) by splitting into two bounds
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
		return nil, fmt.Errorf("compaction failed: %w", err)
	}

	// convert each compacted cell into a 128-bit SubtreeRange [Min, Max]
	ranges := make([]Uint128Range, 0, len(compacted))
	for _, c := range compacted {
		minCell, maxCell := c.SubtreeRange()
		ranges = append(ranges, Uint128Range{
			MinHigh: minCell.High,
			MinLow:  minCell.Low,
			MaxHigh: maxCell.High,
			MaxLow:  maxCell.Low,
		})
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
	ranges := make([]Uint64Range, 0, len(compacted))
	for _, c := range compacted {
		minCell, maxCell := c.SubtreeRange()
		ranges = append(ranges, Uint64Range{
			Min: uint64(minCell),
			Max: uint64(maxCell),
		})
	}

	return MergeRanges(ranges), nil
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
	ranges := make([]Uint128Range, 0, len(compacted))
	for _, c := range compacted {
		minCell, maxCell := c.SubtreeRange()
		ranges = append(ranges, Uint128Range{
			MinHigh: minCell.High,
			MinLow:  minCell.Low,
			MaxHigh: maxCell.High,
			MaxLow:  maxCell.Low,
		})
	}

	return MergeRanges128(ranges), nil
}

// MergeRanges sorts and merges overlapping or contiguous 64-bit uint64 ranges.
func MergeRanges(ranges []Uint64Range) []Uint64Range {
	if len(ranges) <= 1 {
		return ranges
	}

	// sort ranges by Min bound
	sort.Slice(ranges, func(i, j int) bool {
		return ranges[i].Min < ranges[j].Min
	})

	// merge overlapping or adjacent ranges
	merged := make([]Uint64Range, 0, len(ranges))
	merged = append(merged, ranges[0])

	for _, curr := range ranges[1:] {
		last := &merged[len(merged)-1]
		if curr.Min <= last.Max+1 { // overlapping or contiguous
			if curr.Max > last.Max {
				last.Max = curr.Max
			}
		} else {
			merged = append(merged, curr)
		}
	}
	return merged
}

// MergeRanges128 sorts and merges overlapping or contiguous 128-bit Uint128Ranges.
// Compares MinHigh first, then MinLow.
func MergeRanges128(ranges []Uint128Range) []Uint128Range {
	if len(ranges) <= 1 {
		return ranges
	}

	// sort by High bits first, then Low bits
	sort.Slice(ranges, func(i, j int) bool {
		if ranges[i].MinHigh == ranges[j].MinHigh {
			return ranges[i].MinLow < ranges[j].MinLow
		}
		return ranges[i].MinHigh < ranges[j].MinHigh
	})

	// merge overlapping or adjacent 128-bit spans
	merged := make([]Uint128Range, 0, len(ranges))
	merged = append(merged, ranges[0])

	for _, curr := range ranges[1:] {
		last := &merged[len(merged)-1]

		// check if High regions match and Low ranges overlap/touch
		if curr.MinHigh == last.MaxHigh && curr.MinLow <= last.MaxLow+1 {
			if curr.MaxLow > last.MaxLow {
				last.MaxLow = curr.MaxLow
			}
		} else {
			merged = append(merged, curr)
		}
	}
	return merged
}
