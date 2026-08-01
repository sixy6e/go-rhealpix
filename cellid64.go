package rhealpix

import "fmt"

const MaxResolution64 uint8 = 14

// CellID64 represents a 64-bit bit-packed rHEALPix cell identifier.
// Supports resolution levels 0 through 14 (~100m–1km global accuracy).
//
// Layout (MSB to LSB):
//   - Bits 61..63 (3 bits): Base Facet (0..5)
//   - Bits 56..60 (5 bits): Resolution Depth (0..14)
//   - Bits 0..55  (56 bits): Levels 1..14 Sub-cells (14 x 4-bit nibbles: 0..8)
type CellID64 uint64

// PackCellID64 builds a CellID64 from facet, resolution depth, and sub-cell path slice.
func PackCellID64(facet uint8, res uint8, path []uint8) (CellID64, error) {
	if err := ValidateCellParams(facet, res, MaxResolution64, len(path)); err != nil {
		return 0, err
	}

	id := ((uint64(facet) & FacetMask) << FacetShift) | ((uint64(res) & ResMask) << ResShift)

	for i := 0; i < int(res); i++ {
		subCell := path[i]
		if subCell > 8 {
			return 0, fmt.Errorf("invalid sub-cell index %d at depth %d (must be 0..8)", subCell, i+1)
		}
		shift := 52 - (i * 4)
		id |= (uint64(subCell) & SubCellMask) << shift
	}

	return CellID64(id), nil
}

// Facet returns the base face index (0..5).
func (id CellID64) Facet() uint8 {
	return uint8((uint64(id) >> FacetShift) & FacetMask)
}

// Resolution returns the current resolution depth (0..14).
func (id CellID64) Resolution() uint8 {
	return uint8((uint64(id) >> ResShift) & ResMask)
}

// IsZero returns true if the cell ID is uninitialised.
func (id CellID64) IsZero() bool {
	return id == 0
}

// SubtreeRange calculates the [Min, Max] 64-bit integer range enclosing ALL child cells.
// Essential for direct TileDB, RocksDB, or Postgres B-Tree 1D spatial range queries.
func (id CellID64) SubtreeRange() (CellID64, CellID64) {
	res := id.Resolution()
	minBound := uint64(id)

	if res == MaxResolution64 {
		return id, id
	}

	// unused trailing path bits: (14 - res) * 4
	unusedBits := (14 - res) * 4
	pathMask := (uint64(1) << unusedBits) - 1

	// set max resolution header (14) so maxBound numerically encloses all deeper child keys
	facetBits := (uint64(id.Facet()) & FacetMask) << FacetShift
	maxResBits := (uint64(MaxResolution64) & ResMask) << ResShift
	existingPathBits := uint64(id) & 0x00FFFFFFFFFFFFFF

	maxBound := facetBits | maxResBits | existingPathBits | pathMask

	return CellID64(minBound), CellID64(maxBound)
}

// Parent returns the parent cell at targetLevel by masking out lower level nibbles.
func (id CellID64) Parent(targetLevel uint8) (CellID64, error) {
	curRes := id.Resolution()
	if targetLevel > curRes {
		return 0, fmt.Errorf("target level %d higher than current level %d", targetLevel, curRes)
	}
	if targetLevel == curRes {
		return id, nil
	}

	// shift out the trailing unused 4-bit nibbles
	shift := (14 - targetLevel) * 4
	mask := uint64(0xFFFFFFFFFFFFFFFF) << shift

	facetBits := (uint64(id.Facet()) & FacetMask) << FacetShift
	resBits := (uint64(targetLevel) & ResMask) << ResShift
	pathBits := (uint64(id) & mask) & 0x00FFFFFFFFFFFFFF

	return CellID64(facetBits | resBits | pathBits), nil
}

// CommonAncestor64 calculates the lowest common ancestor (LCA) between two 64-bit cells.
func CommonAncestor64(a, b CellID64) (CellID64, error) {
	if a.Facet() != b.Facet() {
		return 0, fmt.Errorf("cells belong to different base facets (%d vs %d)", a.Facet(), b.Facet())
	}

	minRes := a.Resolution()
	if b.Resolution() < minRes {
		minRes = b.Resolution()
	}

	lcaLevel := uint8(0)
	for lvl := uint8(1); lvl <= minRes; lvl++ {
		pA, _ := a.Parent(lvl)
		pB, _ := b.Parent(lvl)
		if pA != pB {
			break
		}
		lcaLevel = lvl
	}

	return a.Parent(lcaLevel)
}

// Uint64 returns the raw 64-bit unsigned integer register.
func (id CellID64) Uint64() uint64 {
	return uint64(id)
}

// Hex returns a 16-character zero-padded hexadecimal representation (e.g. "0x6080120000000000").
func (id CellID64) Hex() string {
	return fmt.Sprintf("0x%016X", uint64(id))
}

// Binary returns the full 64-bit binary string.
func (id CellID64) Binary() string {
	return fmt.Sprintf("%064b", uint64(id))
}

// DebugString returns a formatted breakdown of the cell's bit layout.
func (id CellID64) DebugString() string {
	facet := id.Facet()
	res := id.Resolution()

	// extract sub-cell path digits
	digits := make([]uint8, res)
	rawID := uint64(id)
	for i := uint8(0); i < res; i++ {
		shift := 52 - (i * 4)
		digits[i] = uint8((rawID >> shift) & SubCellMask)
	}

	return fmt.Sprintf(
		"CellID64[%s] | Facet: %d | Res: %d | Digits: %v | Hex: 0x%016X",
		id.String(), facet, res, digits, uint64(id),
	)
}

// DecodeCellID64 unpacks a 64-bit CellID64 into its base facet (0..5),
// target resolution, and sub-cell digit path (0..8 for each level).
func DecodeCellID64(cellID CellID64) (facet uint8, res uint8, path []uint8, err error) {
	if cellID == 0 {
		return 0, 0, nil, ErrInvalidCellID
	}

	// extract Facet (Top 3 bits: bits 61-63)
	facet = uint8((uint64(cellID) >> FacetShift) & FacetMask)
	if facet > 5 {
		return 0, 0, nil, ErrInvalidCellID
	}

	// extract Resolution (Next 5 bits: bits 56-60)
	res = uint8((uint64(cellID) >> ResShift) & ResMask)
	if res > MaxResolution64 {
		return 0, 0, nil, ErrResolutionExceeded
	}

	path = make([]uint8, res)

	// extract 4-bit digits starting from bit 52 downwards
	bitOffset := 52
	for r := uint8(0); r < res; r++ {
		digit := uint8((uint64(cellID) >> bitOffset) & SubCellMask)
		if digit > 8 {
			return 0, 0, nil, ErrInvalidCellID
		}
		path[r] = digit
		bitOffset -= 4
	}

	return facet, res, path, nil
}

// UnpackCellID64 is an alias for DecodeCellID64.
func UnpackCellID64(cellID CellID64) (facet uint8, res uint8, path []uint8, err error) {
	return DecodeCellID64(cellID)
}
