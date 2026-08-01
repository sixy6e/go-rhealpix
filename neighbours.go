package rhealpix

import (
	"fmt"
	"math/big"
)

// NeighbourSet holds the 8 surrounding topological cells (4 edge-adjacent, 4 corner-adjacent) for any CellID type.
type NeighbourSet[T CellID] struct {
	N, NE, E, SE, S, SW, W, NW T
}

// Slice returns all non-zero neighbour cell IDs as a slice.
func (n NeighbourSet[T]) Slice() []T {
	cells := make([]T, 0, 8)
	for _, c := range []T{n.N, n.NE, n.E, n.SE, n.S, n.SW, n.W, n.NW} {
		// Use type assertion / zero check for T
		var anyC any = c
		switch v := anyC.(type) {
		case CellID64:
			if !v.IsZero() {
				cells = append(cells, c)
			}
		case CellID128:
			if !v.IsZero() {
				cells = append(cells, c)
			}
		}
	}
	return cells
}

// Neighbours returns the 8 immediate topological neighbour cells for any CellID (64-bit or 128-bit).
func Neighbours[T CellID](id T) (NeighbourSet[T], error) {
	var facet, res uint8
	var path []uint8
	var err error

	// decode
	var anyID any = id
	switch v := anyID.(type) {
	case CellID64:
		if v.IsZero() {
			return NeighbourSet[T]{}, fmt.Errorf("zero cell ID")
		}
		facet, res, path, err = DecodeCellID64(v)
	case CellID128:
		if v.IsZero() {
			return NeighbourSet[T]{}, fmt.Errorf("zero cell ID")
		}
		facet, res, path, err = DecodeCellID128(v)
	}

	if err != nil {
		return NeighbourSet[T]{}, err
	}

	// compute 3^res grid size
	gridSize := big.NewInt(1)
	three := big.NewInt(3)
	for i := uint8(0); i < res; i++ {
		gridSize.Mul(gridSize, three)
	}

	// path digits -> 2D local grid coordinates (row, col)
	row, col := pathToGridCoords(path, res)

	// offsets for N, NE, E, SE, S, SW, W, NW
	offsets := [8][2]int64{
		{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1},
	}

	var result [8]T
	zero := big.NewInt(0)

	// calculate each neighbour position
	for i, off := range offsets {
		nRow := new(big.Int).Add(row, big.NewInt(off[0]))
		nCol := new(big.Int).Add(col, big.NewInt(off[1]))
		nFacet := facet

		// Intra-facet vs Cross-facet check
		inRow := nRow.Cmp(zero) >= 0 && nRow.Cmp(gridSize) < 0
		inCol := nCol.Cmp(zero) >= 0 && nCol.Cmp(gridSize) < 0

		if !inRow || !inCol {
			nFacet, nRow, nCol = wrapFacetCoords(facet, nRow, nCol, gridSize)
		}

		if nFacet != 255 {
			nPath := gridCoordsToPath(nRow, nCol, res)

			// encode
			switch any(id).(type) {
			case CellID64:
				enc, _ := PackCellID64(nFacet, res, nPath)
				result[i] = any(enc).(T)
			case CellID128:
				enc, _ := PackCellID128(nFacet, res, nPath)
				result[i] = any(enc).(T)
			}
		}
	}

	return NeighbourSet[T]{
		N:  result[0],
		NE: result[1],
		E:  result[2],
		SE: result[3],
		S:  result[4],
		SW: result[5],
		W:  result[6],
		NW: result[7],
	}, nil
}

// Convenient type aliases / wrappers for explicit non-generic functions
func Neighbours64(id CellID64) (NeighbourSet[CellID64], error) {
	return Neighbours(id)
}

func Neighbours128(id CellID128) (NeighbourSet[CellID128], error) {
	return Neighbours(id)
}

// --- Shared Grid Coordinate Helpers ---

func pathToGridCoords(path []uint8, res uint8) (row, col *big.Int) {
	row = big.NewInt(0)
	col = big.NewInt(0)
	three := big.NewInt(3)

	for i := uint8(0); i < res; i++ {
		digit := int64(path[i])
		r := big.NewInt(digit / 3)
		c := big.NewInt(digit % 3)

		row.Mul(row, three).Add(row, r)
		col.Mul(col, three).Add(col, c)
	}
	return
}

func gridCoordsToPath(row, col *big.Int, res uint8) []uint8 {
	path := make([]uint8, res)
	rWork := new(big.Int).Set(row)
	cWork := new(big.Int).Set(col)

	three := big.NewInt(3)
	remR := new(big.Int)
	remC := new(big.Int)

	for i := int(res) - 1; i >= 0; i-- {
		rWork.DivMod(rWork, three, remR)
		cWork.DivMod(cWork, three, remC)

		digit := uint8(remR.Int64()*3 + remC.Int64())
		path[i] = digit
	}
	return path
}

func wrapFacetCoords(facet uint8, row, col, gridSize *big.Int) (outFacet uint8, outRow, outCol *big.Int) {
	outRow = new(big.Int).Set(row)
	outCol = new(big.Int).Set(col)
	outFacet = facet

	zero := big.NewInt(0)
	one := big.NewInt(1)

	// Equatorial Facet Wrapping (Facets 1..4)
	if facet >= 1 && facet <= 4 {
		if col.Cmp(zero) < 0 {
			outCol.Add(col, gridSize)
			outFacet = ((facet - 2 + 4) % 4) + 1
			return
		}
		if col.Cmp(gridSize) >= 0 {
			outCol.Sub(col, gridSize)
			outFacet = (facet % 4) + 1
			return
		}
		if row.Cmp(gridSize) >= 0 {
			outFacet = 0
			outRow.Set(col)
			diff := new(big.Int).Sub(row, gridSize)
			diff.Add(diff, one)
			outCol.Sub(gridSize, one).Sub(outCol, diff)
			return
		}
		if row.Cmp(zero) < 0 {
			outFacet = 5
			outRow.Set(col)
			outCol.Neg(row).Sub(outCol, one)
			return
		}
	}

	// Polar Facet Wrapping (Facet 0: North, Facet 5: South)
	gridMinusOne := new(big.Int).Sub(gridSize, one)

	if facet == 0 { // North
		if row.Cmp(gridSize) >= 0 { // Top edge -> Facet 3 (Q)
			outFacet = 3
			outRow.Set(gridMinusOne)
			outCol.Set(col)
		} else if col.Cmp(gridSize) >= 0 { // Right edge -> Facet 4 (R)
			outFacet = 4
			outRow.Set(gridMinusOne)
			outCol.Set(row)
		} else if row.Cmp(zero) < 0 { // Bottom edge -> Facet 1 (O)
			outFacet = 1
			outRow.Set(gridMinusOne)
			outCol.Set(col)
		} else if col.Cmp(zero) < 0 { // Left edge -> Facet 2 (P)
			outFacet = 2
			outRow.Set(gridMinusOne)
			outCol.Set(row)
		}
		return
	}

	if facet == 5 { // South
		if row.Cmp(zero) < 0 { // Bottom edge -> Facet 1 (O)
			outFacet = 1
			outRow.Set(zero)
			outCol.Set(col)
		} else if col.Cmp(gridSize) >= 0 { // Right edge -> Facet 4 (R)
			outFacet = 4
			outRow.Set(zero)
			outCol.Set(row)
		} else if row.Cmp(gridSize) >= 0 { // Top edge -> Facet 3 (Q)
			outFacet = 3
			outRow.Set(zero)
			outCol.Set(col)
		} else if col.Cmp(zero) < 0 { // Left edge -> Facet 2 (P)
			outFacet = 2
			outRow.Set(zero)
			outCol.Set(row)
		}
		return
	}

	return 255, nil, nil
}
