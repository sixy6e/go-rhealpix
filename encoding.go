package rhealpix

import (
	"fmt"
	"strings"
)

// Base facet character mapping for facets 0..5 per the rHEALPix specification (Gibb 2016):
// Facet 0: North Pole N
// Facet 1: Equatorial O [-180°, -90°]
// Facet 2: Equatorial P [-90°, 0°]
// Facet 3: Equatorial Q [0°, 90°]
// Facet 4: Equatorial R [90°, 180°]
// Facet 5: South Pole S
var RootFacetChars = [6]byte{'N', 'O', 'P', 'Q', 'R', 'S'}

// --- CellID64 Encoding & Decoding ---
// String returns the canonical string representation of CellID64 (e.g., "Q012", "N").
// Implements fmt.Stringer.
func (id CellID64) String() string {
	if id.IsZero() {
		return ""
	}

	facet := id.Facet()
	res := id.Resolution()

	var sb strings.Builder
	sb.Grow(1 + int(res))

	if facet < uint8(len(RootFacetChars)) {
		sb.WriteByte(RootFacetChars[facet])
	} else {
		// handle corrupt sub-cell digits e.g. Q01?2 instead of panicking
		sb.WriteByte('?')
	}

	rawID := uint64(id)
	for i := 0; i < int(res); i++ {
		shift := 52 - (i * 4)
		subCell := (rawID >> shift) & SubCellMask

		if subCell <= 8 {
			sb.WriteByte(byte('0' + subCell))
		} else {
			// handle corrupt sub-cell digits e.g. Q01?2 instead of panicking
			sb.WriteByte('?')
		}
	}

	return sb.String()
}

// ParseCellID64 parses a canonical cell ID string (e.g., "Q012", "N85") back into CellID64.
func ParseCellID64(s string) (CellID64, error) {
	facet, res, path, err := parseSUIDPayload(s, MaxResolution64)
	if err != nil {
		return 0, err
	}
	return PackCellID64(facet, res, path)
}

// --- CellID128 Encoding & Decoding ---
// String returns the canonical string representation of CellID128.
// Implements fmt.Stringer.
func (id CellID128) String() string {
	if id.IsZero() {
		return ""
	}

	facet := id.Facet()
	res := id.Resolution()

	var sb strings.Builder
	sb.Grow(1 + int(res))

	if facet < uint8(len(RootFacetChars)) {
		sb.WriteByte(RootFacetChars[facet])
	} else {
		// handle corrupt sub-cell digits e.g. Q01?2 instead of panicking
		sb.WriteByte('?')
	}

	for i := uint8(0); i < res; i++ {
		var subCell uint64
		if i < 14 {
			shift := 52 - (i * 4)
			subCell = (id.High >> shift) & SubCellMask
		} else {
			shift := 60 - ((i - 14) * 4)
			subCell = (id.Low >> shift) & SubCellMask
		}

		if subCell <= 8 {
			sb.WriteByte(byte('0' + subCell))
		} else {
			// handle corrupt sub-cell digits e.g. Q01?2 instead of panicking
			sb.WriteByte('?')
		}
	}

	return sb.String()
}

// ParseCellID128 parses a canonical cell ID string back into CellID128.
func ParseCellID128(s string) (CellID128, error) {
	facet, res, path, err := parseSUIDPayload(s, MaxResolution128)
	if err != nil {
		return CellID128{}, err
	}
	return PackCellID128(facet, res, path)
}

// --- Shared Internal Parser ---
func parseSUIDPayload(s string, maxAllowedRes uint8) (uint8, uint8, []uint8, error) {
	if len(s) == 0 {
		return 0, 0, nil, fmt.Errorf("empty cell ID string")
	}

	// parse root region character
	rootChar := s[0]
	facet := uint8(255)
	for i, c := range RootFacetChars {
		if c == rootChar {
			facet = uint8(i)
			break
		}
	}
	if facet == 255 {
		return 0, 0, nil, fmt.Errorf("invalid root facet character '%c' in SUID %q", rootChar, s)
	}

	// parse sub-cell path digits
	res := uint8(len(s) - 1)
	if res > maxAllowedRes {
		return 0, 0, nil, fmt.Errorf("resolution depth %d exceeds max %d", res, maxAllowedRes)
	}

	path := make([]uint8, res)
	for i := 0; i < int(res); i++ {
		ch := s[i+1]
		if ch < '0' || ch > '8' {
			return 0, 0, nil, fmt.Errorf("invalid sub-cell digit '%c' at position %d", ch, i+1)
		}
		path[i] = ch - '0'
	}

	return facet, res, path, nil
}
