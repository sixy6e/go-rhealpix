package rhealpix

import (
	"errors"
	"fmt"
)

// Bit-masking constants shared across cell ID implementations
const (
	FacetShift = 61
	ResShift   = 56

	// 3 bits (Base facets 0..5)
	FacetMask uint64 = 0x7

	// 5 bits (Resolution levels 0..31)
	ResMask uint64 = 0x1F

	// 4 bits (Sub-cell digits 0..8)
	SubCellMask uint64 = 0xF
)

var (
	ErrInvalidCellID      = errors.New("invalid or unassigned cell id")
	ErrResolutionExceeded = errors.New("resolution exceeds maximum capacity")
)

// CellShape represents the topological geometry type of an rHEALPix cell in geographic space.
type CellShape uint8

const (
	// Standard equatorial or planar interior quad
	ShapeQuad CellShape = 0

	// Polar cap center cell (nested under sub-cell 4)
	ShapeCap CellShape = 1

	// Triangular seam cell along polar cardinal axes
	ShapeDart CellShape = 2

	// Polar cell straddling a triangle boundary requiring vertex rotation
	ShapeSkewQuad CellShape = 3
)

// String returns the human-readable name of the cell shape.
func (s CellShape) String() string {
	switch s {
	case ShapeQuad:
		return "quad"
	case ShapeCap:
		return "cap"
	case ShapeDart:
		return "dart"
	case ShapeSkewQuad:
		return "skew_quad"
	default:
		return "unknown"
	}
}

// CellID is a generic constraint satisfied by both 64-bit and 128-bit cell identifiers.
// This allows high-level algorithms to operate polymorphically at compile time.
type CellID interface {
	CellID64 | CellID128
}

// Cell defines common behavioral methods for cell types, returning T for hierarchy methods.
type Cell[T any] interface {
	Facet() uint8
	Resolution() uint8
	IsZero() bool
	Parent(targetLevel uint8) (T, error)
	String() string
}

// SubtreeRanger is an optional interface implemented by cell IDs that support
// calculating 1D continuous integer bounds for spatial database range scans.
type SubtreeRanger[T CellID] interface {
	SubtreeRange() (minBound T, maxBound T)
}

// --- Topology Classification Functions ---

// CellShape64 determines the geographic shape classification of a CellID64.
func CellShape64(id CellID64) CellShape {
	facet := id.Facet()
	// equatorial facets (O, P, Q, R -> Facets 1..4) are always standard quads
	if facet >= 1 && facet <= 4 {
		return ShapeQuad
	}

	res := id.Resolution()
	if res == 0 {
		return ShapeQuad
	}

	rawID := uint64(id)

	// Cap check: all digits must be 4
	isCap := true
	for i := uint8(0); i < res; i++ {
		shift := 52 - (i * 4)
		if uint8((rawID>>shift)&SubCellMask) != 4 {
			isCap = false
			break
		}
	}
	if isCap {
		return ShapeCap
	}

	// Dart check 1: all digits in main diagonal {0, 4, 8}
	isDart1 := true
	for i := uint8(0); i < res; i++ {
		shift := 52 - (i * 4)
		d := uint8((rawID >> shift) & SubCellMask)
		if d != 0 && d != 4 && d != 8 {
			isDart1 = false
			break
		}
	}
	if isDart1 {
		return ShapeDart
	}

	// Dart check 2: all digits in anti-diagonal {2, 4, 6}
	isDart2 := true
	for i := uint8(0); i < res; i++ {
		shift := 52 - (i * 4)
		d := uint8((rawID >> shift) & SubCellMask)
		if d != 2 && d != 4 && d != 6 {
			isDart2 = false
			break
		}
	}
	if isDart2 {
		return ShapeDart
	}

	// fallback: must be a Skew Quad
	return ShapeSkewQuad
}

// CellShape128 determines the geographic shape classification of a CellID128.
func CellShape128(id CellID128) CellShape {
	facet := id.Facet()
	// equatorial facets (O, P, Q, R -> Facets 1..4) are always standard quads
	if facet >= 1 && facet <= 4 {
		return ShapeQuad
	}

	res := id.Resolution()
	if res == 0 {
		return ShapeQuad
	}

	// helper to extract a 4-bit sub-cell digit from 128-bit ID
	getDigit := func(i uint8) uint8 {
		if i < 14 {
			shift := 52 - (i * 4)
			return uint8((id.High >> shift) & SubCellMask)
		}
		shift := 60 - ((i - 14) * 4)
		return uint8((id.Low >> shift) & SubCellMask)
	}

	// Cap check: all digits must be 4
	isCap := true
	for i := uint8(0); i < res; i++ {
		if getDigit(i) != 4 {
			isCap = false
			break
		}
	}
	if isCap {
		return ShapeCap
	}

	// Dart check 1: all digits in main diagonal {0, 4, 8}
	isDart1 := true
	for i := uint8(0); i < res; i++ {
		d := getDigit(i)
		if d != 0 && d != 4 && d != 8 {
			isDart1 = false
			break
		}
	}
	if isDart1 {
		return ShapeDart
	}

	// Dart check 2: all digits in anti-diagonal {2, 4, 6}
	isDart2 := true
	for i := uint8(0); i < res; i++ {
		d := getDigit(i)
		if d != 2 && d != 4 && d != 6 {
			isDart2 = false
			break
		}
	}
	if isDart2 {
		return ShapeDart
	}

	// fallback: must be a Skew Quad
	return ShapeSkewQuad
}

// --- Generic Utility Functions ---

// IsSameFacet checks if two cells share the same Level-0 base face.
func IsSameFacet[T CellID](a, b T) bool {
	var anyA any = a
	var anyB any = b
	switch vA := anyA.(type) {
	case CellID64:
		return vA.Facet() == anyB.(CellID64).Facet()
	case CellID128:
		return vA.Facet() == anyB.(CellID128).Facet()
	default:
		return false
	}
}

// CompareResolution returns -1 if a is coarser than b, 1 if finer, or 0 if at the same level.
func CompareResolution[T CellID](a, b T) int {
	var resA, resB uint8
	var anyA any = a
	var anyB any = b

	switch vA := anyA.(type) {
	case CellID64:
		resA = vA.Resolution()
		resB = anyB.(CellID64).Resolution()
	case CellID128:
		resA = vA.Resolution()
		resB = anyB.(CellID128).Resolution()
	}

	if resA < resB {
		return -1
	}
	if resA > resB {
		return 1
	}
	return 0
}

// ValidateCellParams validates common rHEALPix cell construction inputs.
func ValidateCellParams(facet uint8, res uint8, maxRes uint8, pathLen int) error {
	if facet > 5 {
		return fmt.Errorf("invalid base facet: %d (must be 0..5)", facet)
	}
	if res > maxRes {
		return fmt.Errorf("resolution depth %d exceeds max %d for this identifier type", res, maxRes)
	}
	if pathLen < int(res) {
		return fmt.Errorf("path length (%d) is shorter than target resolution (%d)", pathLen, res)
	}
	return nil
}

// Compact compresses a slice of any CellID type (CellID64 or CellID128)
// by recursively collapsing complete 9-child sets into their parent cell.
// This is similar in design to Uber's H3 compactCells algorithm, but substantially
// simpler due to how rHEALPix cells are constructed compared to H3.
// H3 doesn't perfectly subdivide into smaller hexagons.
// rHEALPix parent cells are the strict bounds of its 9 child cells.
// H3 children cross parent edges (not strictly contained within).
func Compact[T interface {
	CellID
	Cell[T]
}](cells []T) ([]T, error) {
	if len(cells) == 0 {
		return nil, nil
	}

	// deduplicate input set and find max resolution directly via Cell[T] interface
	currentSet := make(map[T]bool, len(cells))
	maxRes := uint8(0)
	for _, c := range cells {
		if c.IsZero() {
			continue
		}
		currentSet[c] = true
		if res := c.Resolution(); res > maxRes {
			maxRes = res
		}
	}

	// base level-0 facets cannot be compacted further
	if maxRes == 0 {
		result := make([]T, 0, len(currentSet))
		for c := range currentSet {
			result = append(result, c)
		}
		return result, nil
	}

	var compactedResult []T

	// process bottom-up from highest resolution down to Level 1
	for res := maxRes; res > 0; res-- {
		parentMap := make(map[T][]T)
		var nextSet []T

		// group current resolution cells by their parent at res-1
		for c := range currentSet {
			if c.Resolution() == res {
				parent, err := c.Parent(res - 1)
				if err != nil {
					return nil, err
				}
				parentMap[parent] = append(parentMap[parent], c)
			} else {
				// coarser cells stay in the pool for higher loop passes
				nextSet = append(nextSet, c)
			}
		}

		// reset active tracking set for the next pass
		currentSet = make(map[T]bool)
		for _, c := range nextSet {
			currentSet[c] = true
		}

		// check parent counts: 9 children = Collapse! < 9 children = Keep leaves
		for parent, children := range parentMap {
			if len(children) == 9 {
				// full 3x3 block present -> promote parent to active set
				currentSet[parent] = true
			} else {
				// incomplete branch -> add leaf children to final result
				compactedResult = append(compactedResult, children...)
			}
		}

		if len(currentSet) == 0 {
			break
		}
	}

	// add remaining top level or uncompacted parent cells to final output
	for c := range currentSet {
		compactedResult = append(compactedResult, c)
	}

	return compactedResult, nil
}
