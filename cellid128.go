package rhealpix

import "fmt"

const MaxResolution128 uint8 = 30

// CellID128 represents a 128-bit bit-packed rHEALPix cell identifier.
// Supports resolution levels 0 through 30 (~sub-millimeter global accuracy).
//
// ============================================================================
// GLOBAL 128-BIT LAYOUT (MSB to LSB across 128 bits):
// ============================================================================
//
//	High Word (Bits 64..127):
//	  - Bits 125..127 (3 bits)  : Base Facet (0..5)
//	  - Bits 120..124 (5 bits)  : Resolution Depth (0..30)
//	  - Bits 64..119  (56 bits) : Levels 1..14 Sub-cells (14 x 4-bit nibbles: 0..8)
//	Low Word (Bits 0..63):
//	  - Bits 0..63    (64 bits) : Levels 15..30 Sub-cells (16 x 4-bit nibbles: 0..8)
//
// ============================================================================
// WORD-RELATIVE BIT OFFSETS:
// ============================================================================
//
//	HIGH WORD (64 bits):
//	  [ 63..61 ] Facet ID (0..5, 3 bits)
//	  [ 60..56 ] Resolution Depth (R = 0..30, 5 bits)
//	  [ 55..52 ] Level 1 Path Digit (0..8, 4 bits)
//	  [ 51..48 ] Level 2 Path Digit (0..8, 4 bits)
//	    ...
//	  [  3..0  ] Level 14 Path Digit (0..8, 4 bits)
//
//	LOW WORD (64 bits):
//	  [ 63..60 ] Level 15 Path Digit (0..8, 4 bits)
//	  [ 59..56 ] Level 16 Path Digit (0..8, 4 bits)
//	    ...
//	  [  3..0  ] Level 30 Path Digit (0..8, 4 bits)
//
// ============================================================================
type CellID128 struct {
	High uint64
	Low  uint64
}

// PackCellID128 builds a CellID128 from facet, resolution depth, and sub-cell path slice.
func PackCellID128(facet uint8, res uint8, path []uint8) (CellID128, error) {
	if err := ValidateCellParams(facet, res, MaxResolution128, len(path)); err != nil {
		return CellID128{}, err
	}

	high := ((uint64(facet) & FacetMask) << FacetShift) | ((uint64(res) & ResMask) << ResShift)
	var low uint64

	for i := 0; i < int(res); i++ {
		subCell := path[i]
		if subCell > 8 {
			return CellID128{}, fmt.Errorf("invalid sub-cell index %d at depth %d (must be 0..8)", subCell, i+1)
		}

		nibble := uint64(subCell) & SubCellMask

		if i < 14 {
			shift := 52 - (i * 4)
			high |= nibble << shift
		} else {
			shift := 60 - ((i - 14) * 4)
			low |= nibble << shift
		}
	}

	return CellID128{High: high, Low: low}, nil
}

// Facet returns the base face index (0..5).
func (id CellID128) Facet() uint8 {
	return uint8((id.High >> FacetShift) & FacetMask)
}

// Resolution returns the current resolution depth (0..30).
func (id CellID128) Resolution() uint8 {
	return uint8((id.High >> ResShift) & ResMask)
}

// IsZero returns true if the cell ID is uninitialised.
func (id CellID128) IsZero() bool {
	return id.Facet() > 5 || id.Resolution() > MaxResolution128
}

// SubtreeRange calculates the [Min, Max] 128-bit integer range enclosing ALL child cells
// down to targetRes. Essential for direct TileDB, RocksDB, or B-Tree 1D spatial range queries.
// If the cell is already at or below targetRes, Min and Max are identical [id, id].
func (id CellID128) SubtreeRange(targetRes uint8) (CellID128, CellID128) {
	res := id.Resolution()

	// If already a leaf at or below targetRes (or absolute MaxResolution128), return exact single key
	if res >= targetRes || res >= MaxResolution128 {
		return id, id
	}

	minBound := id
	maxBound := id

	// Fill trailing 4-bit nibbles from current 'res' up to 'targetRes-1' with max sub-cell digit 8
	for r := res; r < targetRes; r++ {
		if r < 14 {
			// Levels 1..14 affect High word (bits 52 down to 0)
			s := 52 - (r * 4)
			maxBound.High |= (uint64(8) << s)
		} else {
			// Levels 15..30 affect Low word (bits 60 down to 0)
			s := 60 - ((r - 14) * 4)
			maxBound.Low |= (uint64(8) << s)
		}
	}

	return minBound, maxBound
}

// SubtreeRangeMax calculates the [Min, Max] range down to absolute MaxResolution128.
func (id CellID128) SubtreeRangeMax() (CellID128, CellID128) {
	return id.SubtreeRange(MaxResolution128)
}

// Equal checks if two 128-bit cell IDs are identical.
func (id CellID128) Equal(other CellID128) bool {
	return id.High == other.High && id.Low == other.Low
}

// Parent returns the parent cell at targetLevel by masking out lower level nibbles.
func (id CellID128) Parent(targetLevel uint8) (CellID128, error) {
	curRes := id.Resolution()
	if targetLevel > curRes {
		return CellID128{}, fmt.Errorf("target level %d higher than current level %d", targetLevel, curRes)
	}
	if targetLevel == curRes {
		return id, nil
	}

	parent := CellID128{}

	// reconstruct high word Header (facet + new target resolution)
	facetBits := (uint64(id.Facet()) & FacetMask) << FacetShift
	resBits := (uint64(targetLevel) & ResMask) << ResShift
	parent.High = facetBits | resBits

	// mask path nibbles based on target Level
	if targetLevel <= 14 {
		// target parent resides entirely within high word (levels 1..14)
		// low word becomes completely zeroed out
		if targetLevel > 0 {
			shift := (14 - targetLevel) * 4
			mask := ^uint64(0) << shift
			parent.High |= (id.High & 0x00FFFFFFFFFFFFFF) & mask
		}
		// parent.Low remains 0
	} else {
		// target parent includes all 14 high levels, plus a subset of low levels (levels 15..30)
		parent.High |= (id.High & 0x00FFFFFFFFFFFFFF)

		shift := (30 - targetLevel) * 4
		mask := ^uint64(0) << shift
		parent.Low = id.Low & mask
	}

	return parent, nil
}

// CommonAncestor128 calculates the lowest common ancestor (LCA) between two 128-bit cells.
func CommonAncestor128(a, b CellID128) (CellID128, error) {
	if a.Facet() != b.Facet() {
		return CellID128{}, fmt.Errorf("cells belong to different base facets (%d vs %d)", a.Facet(), b.Facet())
	}

	minRes := a.Resolution()
	if b.Resolution() < minRes {
		minRes = b.Resolution()
	}

	lcaLevel := uint8(0)

	// compare levels 1..14 in high word (bits 55..0)
	highDiff := (a.High ^ b.High) & 0x00FFFFFFFFFFFFFF
	if highDiff == 0 {
		// all high levels match! LCA is at least level 14.
		lcaLevel = 14
		if lcaLevel > minRes {
			lcaLevel = minRes
		}

		// if minRes extends into low word (levels 15..30), check low word
		if minRes > 14 {
			lowDiff := a.Low ^ b.Low
			for lvl := uint8(15); lvl <= minRes; lvl++ {
				shift := (30 - lvl) * 4
				if (lowDiff >> shift) != 0 {
					break
				}
				lcaLevel = lvl
			}
		}
	} else {
		// divergence happens in high word (levels 1..14)
		for lvl := uint8(1); lvl <= minRes && lvl <= 14; lvl++ {
			shift := (14 - lvl) * 4
			if (highDiff >> shift) != 0 {
				break
			}
			lcaLevel = lvl
		}
	}

	return a.Parent(lcaLevel)
}

// Uint64s returns the raw high and low 64-bit unsigned integer registers.
func (id CellID128) Uint64s() (high uint64, low uint64) {
	return id.High, id.Low
}

// HighWord returns the upper 64-bit word of the 128-bit register.
func (id CellID128) HighWord() uint64 {
	return id.High
}

// LowWord returns the lower 64-bit word of the 128-bit register.
func (id CellID128) LowWord() uint64 {
	return id.Low
}

// Hex returns a 32-character zero-padded hexadecimal representation.
func (id CellID128) Hex() string {
	return fmt.Sprintf("0x%016X%016X", id.High, id.Low)
}

// Binary returns the full 128-bit binary string across high and low words.
func (id CellID128) Binary() string {
	return fmt.Sprintf("%064b%064b", id.High, id.Low)
}

// DebugString returns a formatted breakdown of the 128-bit cell's bit layout across High and Low words.
func (id CellID128) DebugString() string {
	facet := id.Facet()
	res := id.Resolution()

	digits := make([]uint8, res)
	for i := uint8(0); i < res; i++ {
		if i < 14 {
			shift := 52 - (i * 4)
			digits[i] = uint8((id.High >> shift) & SubCellMask)
		} else {
			shift := 60 - ((i - 14) * 4)
			digits[i] = uint8((id.Low >> shift) & SubCellMask)
		}
	}

	return fmt.Sprintf(
		"CellID128[%s] | Facet: %d | Res: %d | Digits: %v | Hex: 0x%016X%016X",
		id.String(), facet, res, digits, id.High, id.Low,
	)
}

// DecodeCellID128 unpacks a 128-bit CellID128 into its base facet (0..5),
// target resolution, and sub-cell digit path (0..8 for each level).
func DecodeCellID128(cellID CellID128) (facet uint8, res uint8, path []uint8, err error) {
	if cellID.High == 0 && cellID.Low == 0 {
		return 0, 0, nil, ErrInvalidCellID
	}

	// extract Facet (Top 3 bits of High word: bits 61..63)
	facet = uint8((cellID.High >> FacetShift) & FacetMask)
	if facet > 5 {
		return 0, 0, nil, ErrInvalidCellID
	}

	// extract Resolution (5 bits of High uint64: bits 56..60)
	res = uint8((cellID.High >> ResShift) & ResMask)
	if res > MaxResolution128 {
		return 0, 0, nil, ErrResolutionExceeded
	}

	path = make([]uint8, res)

	// extract 4-bit path digits across High and Low uint64 words
	for r := uint8(0); r < res; r++ {
		var digit uint8
		if r < 14 {
			shift := 52 - (r * 4)
			digit = uint8((cellID.High >> shift) & SubCellMask)
		} else {
			shift := 60 - ((r - 14) * 4)
			digit = uint8((cellID.Low >> shift) & SubCellMask)
		}

		if digit > 8 {
			return 0, 0, nil, ErrInvalidCellID
		}
		path[r] = digit
	}

	return facet, res, path, nil
}

// UnpackCellID128 is an alias for DecodeCellID128.
func UnpackCellID128(cellID CellID128) (facet uint8, res uint8, path []uint8, err error) {
	return DecodeCellID128(cellID)
}

// Child returns the nth child cell (0..8) at the next resolution level (res + 1).
// Operates purely bitwise without allocating slices.
func (id CellID128) Child(subCellIdx uint8) (CellID128, error) {
	if subCellIdx > 8 {
		return CellID128{}, fmt.Errorf("invalid sub-cell index %d: must be 0..8", subCellIdx)
	}

	res := id.Resolution()
	if res >= MaxResolution128 {
		return CellID128{}, fmt.Errorf("cannot derive child beyond max resolution %d", MaxResolution128)
	}

	nextRes := res + 1

	child := CellID128{
		High: id.High,
		Low:  id.Low,
	}

	// mask out top 8 header bits (bits 63..56) while preserving all 56 path bits
	const pathMask = uint64(0x00FFFFFFFFFFFFFF)

	// update resolution header in high word
	facetBits := (uint64(id.Facet()) & FacetMask) << FacetShift
	resBits := (uint64(nextRes) & ResMask) << ResShift
	child.High = (child.High & pathMask) | facetBits | resBits

	// insert new 4-bit nibble at current level `res`
	nibble := (uint64(subCellIdx) & SubCellMask)

	if res < 14 {
		// level falls in high word (levels 1..14)
		// level 1 (res=0) -> shift 52 (bits 55..52)
		// level 14 (res=13) -> shift 0 (bits 3..0)
		shift := 52 - (res * 4)
		child.High |= (nibble << shift)
	} else {
		// level falls in Low word (levels 15..30)
		// level 15 (res=14) -> shift 60 (bits 63..60)
		// level 30 (res=29) -> shift 0 (bits 3..0)
		shift := 60 - ((res - 14) * 4)
		child.Low |= (nibble << shift)
	}

	return child, nil
}
