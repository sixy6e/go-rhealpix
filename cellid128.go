package rhealpix

import "fmt"

const MaxResolution128 uint8 = 30

// CellID128 represents a 128-bit bit-packed rHEALPix cell identifier.
// Supports resolution levels 0 through 30 (~sub-millimeter global accuracy).
//
// Layout:
//   - High: [Facet: 3b] [Res: 5b] [Levels 1..14 Sub-cells: 56b (14 x 4b)]
//   - Low:  [Levels 15..30 Sub-cells: 64b (16 x 4b)]
//
// HIGH WORD (64 bits):
//
//	[ 63..61 ] Facet ID (Base Level 0)
//	[ 60..57 ] Padding
//	[ 56..53 ] Resolution Depth (R = 0..30)
//	[ 52..49 ] Level 1 Path Digit  (0..8)
//	[ 48..45 ] Level 2 Path Digit  (0..8)
//	  ...
//	[  3..0  ] Level 14 Path Digit (0..8)  <-- End of High Word
//
//	LOW WORD (64 bits):
//	[ 60..57 ] Level 15 Path Digit (0..8)
//	  ...
//	[  0     ] Level 30 Path Digit
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
			// levels 1..14 go into High word (bits 52 down to 0)
			shift := 52 - (i * 4)
			high |= nibble << shift
		} else {
			// levels 15..30 go into Low word (bits 60 down to 0)
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
	return id.High == 0 && id.Low == 0
}

// Equal checks if two 128-bit cell IDs are identical.
func (id CellID128) Equal(other CellID128) bool {
	return id.High == other.High && id.Low == other.Low
}

// SubtreeRange calculates the [Min, Max] 128-bit key range enclosing ALL child cells.
// Works seamlessly across compound (High, Low) database dimensions or fixed 16-byte slices.
func (id CellID128) SubtreeRange() (CellID128, CellID128) {
	res := id.Resolution()
	minBound := id

	if res == MaxResolution128 {
		return id, id
	}

	var maxHigh, maxLow uint64

	// update resolution in High word to MaxResolution128 (30)
	facetBits := (uint64(id.Facet()) & FacetMask) << FacetShift
	maxResBits := (uint64(MaxResolution128) & ResMask) << ResShift

	if res <= 14 {
		// target cell is coarse (Res <= 14). Unused bits span remaining High nibbles and ALL of Low.
		unusedHighBits := (14 - res) * 4
		highMask := (uint64(1) << unusedHighBits) - 1

		existingHighPath := id.High & 0x00FFFFFFFFFFFFFF
		maxHigh = facetBits | maxResBits | existingHighPath | highMask
		maxLow = ^uint64(0) // fill Low completely (0xFFFFFFFFFFFFFFFF)
	} else {
		// high word path is fixed; set max res header and fill remaining unused bits in Low word
		existingHighPath := id.High & 0x00FFFFFFFFFFFFFF
		maxHigh = facetBits | maxResBits | existingHighPath

		unusedLowBits := (30 - res) * 4
		if unusedLowBits > 0 {
			lowMask := (uint64(1) << unusedLowBits) - 1
			maxLow = id.Low | lowMask
		} else {
			maxLow = id.Low
		}
	}

	return minBound, CellID128{High: maxHigh, Low: maxLow}
}

// Parent returns the parent cell at targetLevel by masking out lower-level nibbles.
func (id CellID128) Parent(targetLevel uint8) (CellID128, error) {
	curRes := id.Resolution()
	if targetLevel > curRes {
		return CellID128{}, fmt.Errorf("target level %d higher than current level %d", targetLevel, curRes)
	}
	if targetLevel == curRes {
		return id, nil
	}

	if targetLevel <= 14 {
		// target level is completely within the High word.
		// low word becomes 0 because no digits remain at depth > 14.
		shift := (14 - targetLevel) * 4
		mask := uint64(0xFFFFFFFFFFFFFFFF) << shift

		facetBits := (uint64(id.Facet()) & FacetMask) << FacetShift
		resBits := (uint64(targetLevel) & ResMask) << ResShift
		pathBits := (id.High & mask) & 0x00FFFFFFFFFFFFFF

		return CellID128{
			High: facetBits | resBits | pathBits,
			Low:  0,
		}, nil
	}

	// target level > 14: High word is fully preserved (except resolution header bits),
	// and Low word has its tail nibbles masked out.
	shift := (30 - targetLevel) * 4
	mask := uint64(0xFFFFFFFFFFFFFFFF) << shift

	// update resolution in High word
	resBits := (uint64(targetLevel) & ResMask) << ResShift
	highBits := (id.High &^ (uint64(ResMask) << ResShift)) | resBits

	return CellID128{
		High: highBits,
		Low:  id.Low & mask,
	}, nil
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
	for lvl := uint8(1); lvl <= minRes; lvl++ {
		pA, _ := a.Parent(lvl)
		pB, _ := b.Parent(lvl)
		if !pA.Equal(pB) {
			break
		}
		lcaLevel = lvl
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

	// extract Facet (Top 3 bits of High uint64: bits 61-63)
	facet = uint8((cellID.High >> FacetShift) & FacetMask)
	if facet > 5 {
		return 0, 0, nil, ErrInvalidCellID
	}

	// extract Resolution (next 6 bits of High uint64: bits 55-60)
	res = uint8((cellID.High >> 55) & 0x3F)
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
