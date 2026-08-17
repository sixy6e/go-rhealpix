package rhealpix_test

import (
	"reflect"
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
)

// --- SUID Parsing & Stringer Tests ---

func TestSUIDRoundTrip64(t *testing.T) {
	tests := []string{"Q", "Q0", "Q012", "N85", "S000"}

	for _, original := range tests {
		cell, err := rhealpix.ParseCellID64(original)
		if err != nil {
			t.Fatalf("ParseCellID64(%q) unexpected error: %v", original, err)
		}

		if got := cell.String(); got != original {
			t.Errorf("SUID mismatch: got %q, want %q", got, original)
		}
	}
}

func TestSUIDRoundTrip128(t *testing.T) {
	// tests deep resolution depths beyond level 14
	tests := []string{"Q01234567880123", "R01234567880123456788012"}

	for _, original := range tests {
		cell, err := rhealpix.ParseCellID128(original)
		if err != nil {
			t.Fatalf("ParseCellID128(%q) unexpected error: %v", original, err)
		}

		if got := cell.String(); got != original {
			t.Errorf("128-bit SUID mismatch: got %q, want %q", got, original)
		}
	}
}

// --- Common Ancestor & Range Tests ---

func TestCommonAncestor64(t *testing.T) {
	c1, _ := rhealpix.ParseCellID64("Q0123")
	c2, _ := rhealpix.ParseCellID64("Q0128")

	lca, err := rhealpix.CommonAncestor64(c1, c2)
	if err != nil {
		t.Fatalf("unexpected LCA error: %v", err)
	}

	want := "Q012"
	if lca.String() != want {
		t.Errorf("LCA mismatch: got %q, want %q", lca.String(), want)
	}
}

// OLD; keep around as could use it later to test that we don't try searching higher than ingested
// but also if it gets resolved at a future date, bring the test back.
// func TestSubtreeRange64(t *testing.T) {
// 	cell, err := rhealpix.ParseCellID64("Q012")
// 	if err != nil {
// 		t.Fatalf("ParseCellID64 unexpected error: %v", err)
// 	}
//
// 	minBound, maxBound := cell.SubtreeRange(rhealpix.MaxResolution64)
//
// 	if minBound.Uint64() >= maxBound.Uint64() {
// 		t.Errorf("invalid SubtreeRange bounds: min %s >= max %s", minBound.Hex(), maxBound.Hex())
// 	}
//
// 	// test child cell
// 	childCell, err := rhealpix.ParseCellID64("Q0128")
// 	if err != nil {
// 		t.Fatalf("ParseCellID64 child error: %v", err)
// 	}
//
// 	if childCell.Uint64() < minBound.Uint64() || childCell.Uint64() > maxBound.Uint64() {
// 		t.Errorf("child cell %s (0x%x) fell outside range [%s, %s]",
// 			childCell.String(), childCell.Uint64(), minBound.Hex(), maxBound.Hex())
// 	}
// }

func TestSubtreeRange64(t *testing.T) {
	cell, err := rhealpix.ParseCellID64("Q012")
	if err != nil {
		t.Fatalf("ParseCellID64 unexpected error: %v", err)
	}

	// parse child cell (Level 4)
	childCell, err := rhealpix.ParseCellID64("Q0128")
	if err != nil {
		t.Fatalf("ParseCellID64 child error: %v", err)
	}

	// target resolution must match the child cell's resolution level (Res 4)
	targetRes := childCell.Resolution()
	minBound, maxBound := cell.SubtreeRange(targetRes)

	if minBound.Uint64() >= maxBound.Uint64() {
		t.Errorf("invalid SubtreeRange bounds: min %s >= max %s", minBound.Hex(), maxBound.Hex())
	}

	if childCell.Uint64() < minBound.Uint64() || childCell.Uint64() > maxBound.Uint64() {
		t.Errorf("child cell %s (0x%x) fell outside range [%s, %s]",
			childCell.String(), childCell.Uint64(), minBound.Hex(), maxBound.Hex())
	}
}

// OLD; keep around as could use it later to test that we don't try searching higher than ingested
// but also if it gets resolved at a future date, bring the test back.
//	func TestSubtreeRange128(t *testing.T) {
//		cell, _ := rhealpix.ParseCellID128("Q01234567890123")
//		minBound, maxBound := cell.SubtreeRange(rhealpix.MaxResolution128)
//
//		if minBound.High > maxBound.High || (minBound.High == maxBound.High && minBound.Low >= maxBound.Low) {
//			t.Errorf("invalid SubtreeRange128 bounds: min %s >= max %s", minBound.Hex(), maxBound.Hex())
//		}
//
//		childCell, _ := rhealpix.ParseCellID128("Q012345678901238")
//		if childCell.High < minBound.High || childCell.High > maxBound.High {
//			t.Errorf("child cell High word %x fell outside [%x, %x]", childCell.High, minBound.High, maxBound.High)
//		}
//	}

func TestSubtreeRange128(t *testing.T) {
	// Level 14 cell using valid digits (0..8)
	cell, err := rhealpix.ParseCellID128("Q01234567801234")
	if err != nil {
		t.Fatalf("ParseCellID128 unexpected error: %v", err)
	}

	// Level 15 child cell
	childCell, err := rhealpix.ParseCellID128("Q012345678012348")
	if err != nil {
		t.Fatalf("ParseCellID128 child error: %v", err)
	}

	// target resolution must match childCell's resolution level (Res 15)
	targetRes := childCell.Resolution()
	minBound, maxBound := cell.SubtreeRange(targetRes)

	// validate min < max in 128-bit space
	minIsLess := minBound.High < maxBound.High || (minBound.High == maxBound.High && minBound.Low < maxBound.Low)
	if !minIsLess {
		t.Errorf("invalid SubtreeRange128 bounds: min %s >= max %s", minBound.Hex(), maxBound.Hex())
	}

	// validate childCell falls inside [minBound, maxBound] across both High and Low words
	childGeMin := childCell.High > minBound.High || (childCell.High == minBound.High && childCell.Low >= minBound.Low)
	childLeMax := childCell.High < maxBound.High || (childCell.High == maxBound.High && childCell.Low <= maxBound.Low)

	if !childGeMin || !childLeMax {
		t.Errorf("child cell %s (Hex: %s) fell outside range [%s, %s]",
			childCell.String(), childCell.Hex(), minBound.Hex(), maxBound.Hex())
	}
}

// --- Compaction Tests ---

func TestCompact64_Complete9Children(t *testing.T) {
	// generate all 9 children for parent "Q012"
	children := make([]rhealpix.CellID64, 9)
	suids := []string{"Q0120", "Q0121", "Q0122", "Q0123", "Q0124", "Q0125", "Q0126", "Q0127", "Q0128"}

	for i, s := range suids {
		c, err := rhealpix.ParseCellID64(s)
		if err != nil {
			t.Fatalf("ParseCellID64(%q) error: %v", s, err)
		}
		children[i] = c
	}

	compacted, err := rhealpix.Compact(children)
	if err != nil {
		t.Fatalf("Compact unexpected error: %v", err)
	}

	if len(compacted) != 1 {
		t.Fatalf("expected 1 compacted cell, got %d", len(compacted))
	}

	if got := compacted[0].String(); got != "Q012" {
		t.Errorf("compacted result mismatch: got %q, want \"Q012\"", got)
	}
}

func TestCompact64_IncompleteSet(t *testing.T) {
	// 8 children (missing Q0128) -> should NOT compact
	suids := []string{"Q0120", "Q0121", "Q0122", "Q0123", "Q0124", "Q0125", "Q0126", "Q0127"}

	children := make([]rhealpix.CellID64, len(suids))
	for i, s := range suids {
		c, err := rhealpix.ParseCellID64(s)
		if err != nil {
			t.Fatalf("ParseCellID64(%q) error: %v", s, err)
		}
		children[i] = c
	}

	compacted, err := rhealpix.Compact(children)
	if err != nil {
		t.Fatalf("Compact unexpected error: %v", err)
	}

	if len(compacted) != 8 {
		t.Fatalf("expected 8 uncompacted cells, got %d", len(compacted))
	}
}

func TestCompact64_CascadingMultiLevel(t *testing.T) {
	// Generate all 81 grandchildren for "N0" at resolution 3 (N000 through N088)
	// Should cascade compact: 81 (Res 3) -> 9 (Res 2) -> 1 (Res 1: "N0")
	var children []rhealpix.CellID64

	digits := []string{"0", "1", "2", "3", "4", "5", "6", "7", "8"}
	for _, d1 := range digits {
		for _, d2 := range digits {
			suid := "N0" + d1 + d2
			c, err := rhealpix.ParseCellID64(suid)
			if err != nil {
				t.Fatalf("ParseCellID64(%q) error: %v", suid, err)
			}
			children = append(children, c)
		}
	}

	compacted, err := rhealpix.Compact(children)
	if err != nil {
		t.Fatalf("Compact unexpected error: %v", err)
	}

	if len(compacted) != 1 {
		t.Fatalf("expected 1 cascaded parent cell, got %d", len(compacted))
	}

	if got := compacted[0].String(); got != "N0" {
		t.Errorf("compacted result mismatch: got %q, want \"N0\"", got)
	}
}

func TestCompact128_HighResolution(t *testing.T) {
	base := "P012345678601234"
	children := make([]rhealpix.CellID128, 9)
	digits := []string{"0", "1", "2", "3", "4", "5", "6", "7", "8"}

	for i, d := range digits {
		suid := base + d
		c, err := rhealpix.ParseCellID128(suid)
		if err != nil {
			t.Fatalf("ParseCellID128(%q) error: %v", suid, err)
		}
		children[i] = c
	}

	compacted, err := rhealpix.Compact(children)
	if err != nil {
		t.Fatalf("Compact unexpected error: %v", err)
	}

	if len(compacted) != 1 {
		t.Fatalf("expected 1 compacted cell, got %d", len(compacted))
	}

	if got := compacted[0].String(); got != base {
		t.Errorf("compacted result mismatch: got %q, want %q", got, base)
	}
}

func TestCompact64_MixedInput(t *testing.T) {
	// 9 children of Q012 + 1 independent cell O5
	suids := []string{
		"Q0120", "Q0121", "Q0122", "Q0123", "Q0124", "Q0125", "Q0126", "Q0127", "Q0128",
		"O5",
	}

	var input []rhealpix.CellID64
	for _, s := range suids {
		c, err := rhealpix.ParseCellID64(s)
		if err != nil {
			t.Fatalf("ParseCellID64(%q) error: %v", s, err)
		}
		input = append(input, c)
	}

	compacted, err := rhealpix.Compact(input)
	if err != nil {
		t.Fatalf("Compact unexpected error: %v", err)
	}

	if len(compacted) != 2 {
		t.Fatalf("expected 2 compacted cells, got %d", len(compacted))
	}

	wantSUIDs := map[string]bool{"Q012": true, "O5": true}
	gotSUIDs := map[string]bool{
		compacted[0].String(): true,
		compacted[1].String(): true,
	}

	if !reflect.DeepEqual(gotSUIDs, wantSUIDs) {
		t.Errorf("compacted output mismatch: got %v, want %v", gotSUIDs, wantSUIDs)
	}
}

// --- Test that Cell Facets (Level-0) pack and unpack correctly (roundtrip) ---

func TestPackCellID64_RootFacets(t *testing.T) {
	t.Run("pack Level 0 root cells for all facets with nil path", func(t *testing.T) {
		for facet := uint8(0); facet < 6; facet++ {
			cellID, err := rhealpix.PackCellID64(facet, 0, nil)
			if err != nil {
				t.Fatalf("failed to pack root cell for facet %d: %v", facet, err)
			}

			if cellID.Facet() != facet {
				t.Errorf("facet mismatch: got %d, want %d", cellID.Facet(), facet)
			}

			if cellID.Resolution() != 0 {
				t.Errorf("resolution mismatch: got %d, want 0", cellID.Resolution())
			}

			if cellID.IsZero() {
				t.Errorf("root cell ID for facet %d should not be zero", facet)
			}
		}
	})

	t.Run("pack Level 0 root cells with empty path slice", func(t *testing.T) {
		for facet := uint8(0); facet < 6; facet++ {
			cellID, err := rhealpix.PackCellID64(facet, 0, []uint8{})
			if err != nil {
				t.Fatalf("failed to pack root cell for facet %d with empty slice: %v", facet, err)
			}

			if cellID.Facet() != facet {
				t.Errorf("facet mismatch: got %d, want %d", cellID.Facet(), facet)
			}
		}
	})

	t.Run("reject invalid facet IDs at level 0", func(t *testing.T) {
		invalidFacets := []uint8{6, 7, 255}
		for _, facet := range invalidFacets {
			_, err := rhealpix.PackCellID64(facet, 0, nil)
			if err == nil {
				t.Errorf("expected error when packing invalid facet %d, got nil", facet)
			}
		}
	})

	t.Run("unpack facet across deeper resolution paths", func(t *testing.T) {
		path := []uint8{0, 1, 2, 3, 4, 5, 6, 7, 8} // Level 9 cell
		for facet := uint8(0); facet < 6; facet++ {
			cellID, err := rhealpix.PackCellID64(facet, uint8(len(path)), path)
			if err != nil {
				t.Fatalf("failed to pack deep cell for facet %d: %v", facet, err)
			}

			if cellID.Facet() != facet {
				t.Errorf("facet mismatch at deep res: got %d, want %d", cellID.Facet(), facet)
			}
		}
	})
}

// ============================================================================
// BIT-PACKING TESTS
// ============================================================================

func TestFacet0RootCellRoundTrip(t *testing.T) {
	t.Run("64-bit R0", func(t *testing.T) {
		cell, err := rhealpix.PackCellID64(0, 0, nil)
		if err != nil {
			t.Fatalf("failed to pack Facet 0 Res 0: %v", err)
		}

		if cell.IsZero() {
			t.Fatalf("expected IsZero() to be false for valid R0 cell")
		}

		facet, res, path, err := rhealpix.DecodeCellID64(cell)
		if err != nil {
			t.Fatalf("DecodeCellID64 unexpected error: %v", err)
		}
		if facet != 0 || res != 0 || len(path) != 0 {
			t.Errorf("got facet=%d res=%d path=%v, want facet=0 res=0 path=[]", facet, res, path)
		}
	})

	t.Run("128-bit R0", func(t *testing.T) {
		cell, err := rhealpix.PackCellID128(0, 0, nil)
		if err != nil {
			t.Fatalf("failed to pack 128-bit Facet 0 Res 0: %v", err)
		}

		if cell.IsZero() {
			t.Fatalf("expected IsZero() to be false for valid 128-bit R0 cell")
		}

		facet, res, path, err := rhealpix.DecodeCellID128(cell)
		if err != nil {
			t.Fatalf("DecodeCellID128 unexpected error: %v", err)
		}
		if facet != 0 || res != 0 || len(path) != 0 {
			t.Errorf("got facet=%d res=%d path=%v, want facet=0 res=0 path=[]", facet, res, path)
		}
	})
}

func TestStaleNibbleScrubbing(t *testing.T) {
	digits := []uint8{1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4}

	t.Run("64-bit Stale Nibble Cleaning", func(t *testing.T) {
		deepCell, err := rhealpix.PackCellID64(4, 12, digits)
		if err != nil {
			t.Fatalf("failed to pack deep cell: %v", err)
		}

		parent, err := deepCell.Parent(5)
		if err != nil {
			t.Fatalf("failed to get parent(5): %v", err)
		}

		child, err := parent.Child(8)
		if err != nil {
			t.Fatalf("failed to derive child(8): %v", err)
		}

		facet, res, path, err := rhealpix.DecodeCellID64(child)
		if err != nil {
			t.Fatalf("failed decoding child: %v", err)
		}

		expectedPath := []uint8{1, 2, 3, 4, 5, 8}
		if facet != 4 || res != 6 || !reflect.DeepEqual(path, expectedPath) {
			t.Errorf("got facet=%d res=%d path=%v, want facet=4 res=6 path=%v", facet, res, path, expectedPath)
		}
	})

	t.Run("128-bit Stale Nibble Cleaning", func(t *testing.T) {
		deepCell, err := rhealpix.PackCellID128(4, 12, digits)
		if err != nil {
			t.Fatalf("failed to pack 128-bit deep cell: %v", err)
		}

		parent, err := deepCell.Parent(5)
		if err != nil {
			t.Fatalf("failed to get parent(5): %v", err)
		}

		child, err := parent.Child(8)
		if err != nil {
			t.Fatalf("failed to derive child(8): %v", err)
		}

		facet, res, path, err := rhealpix.DecodeCellID128(child)
		if err != nil {
			t.Fatalf("failed decoding 128-bit child: %v", err)
		}

		expectedPath := []uint8{1, 2, 3, 4, 5, 8}
		if facet != 4 || res != 6 || !reflect.DeepEqual(path, expectedPath) {
			t.Errorf("got facet=%d res=%d path=%v, want facet=4 res=6 path=%v", facet, res, path, expectedPath)
		}
	})
}

func TestWordBoundaryCrossing128(t *testing.T) {
	path14 := []uint8{1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5}
	cell14, err := rhealpix.PackCellID128(2, 14, path14)
	if err != nil {
		t.Fatalf("failed packing level 14 cell: %v", err)
	}

	cell15, err := cell14.Child(7)
	if err != nil {
		t.Fatalf("failed deriving child at level 15: %v", err)
	}

	facet, res, path, err := rhealpix.DecodeCellID128(cell15)
	if err != nil {
		t.Fatalf("failed decoding level 15 cell: %v", err)
	}

	expectedPath := append(path14, 7)
	if facet != 2 || res != 15 || !reflect.DeepEqual(path, expectedPath) {
		t.Errorf("got facet=%d res=%d path=%v, want facet=2 res=15 path=%v", facet, res, path, expectedPath)
	}

	// verify low word receives the 15th digit at shift 60
	expectedLow := uint64(7) << 60
	if cell15.Low != expectedLow {
		t.Errorf("got Low=0x%X, want 0x%X", cell15.Low, expectedLow)
	}
}

// ============================================================================
// SUBTREE RANGE INVARIANT TESTS
// ============================================================================

func TestSubtreeRangeHeaderAlignment(t *testing.T) {
	t.Run("64-bit Header Alignment", func(t *testing.T) {
		parent, _ := rhealpix.PackCellID64(4, 5, []uint8{1, 2, 3, 4, 5})
		targetRes := uint8(12)
		minB, maxB := parent.SubtreeRange(targetRes)

		if minB.Resolution() != targetRes || maxB.Resolution() != targetRes {
			t.Fatalf("got minRes=%d maxRes=%d, want targetRes=%d", minB.Resolution(), maxB.Resolution(), targetRes)
		}
	})

	t.Run("128-bit Header Alignment", func(t *testing.T) {
		parent, _ := rhealpix.PackCellID128(4, 5, []uint8{1, 2, 3, 4, 5})
		targetRes := uint8(20)
		minB, maxB := parent.SubtreeRange(targetRes)

		if minB.Resolution() != targetRes || maxB.Resolution() != targetRes {
			t.Fatalf("got minRes=%d maxRes=%d, want targetRes=%d", minB.Resolution(), maxB.Resolution(), targetRes)
		}
	})
}

func TestSubtreeRangeTrueChildEnclosure(t *testing.T) {
	t.Run("64-bit Enclosure", func(t *testing.T) {
		parent, _ := rhealpix.PackCellID64(4, 5, []uint8{1, 2, 3, 4, 5})
		targetRes := uint8(12)
		minB, maxB := parent.SubtreeRange(targetRes)

		curr := parent
		for r := uint8(5); r < targetRes; r++ {
			curr, _ = curr.Child(r % 9)
		}

		if curr.Uint64() < minB.Uint64() || curr.Uint64() > maxB.Uint64() {
			t.Fatalf("Descendant %s fell outside SubtreeRange [0x%X, 0x%X]", curr.Hex(), minB.Uint64(), maxB.Uint64())
		}
	})

	t.Run("128-bit Enclosure", func(t *testing.T) {
		parent, _ := rhealpix.PackCellID128(4, 5, []uint8{1, 2, 3, 4, 5})
		targetRes := uint8(20)
		minB, maxB := parent.SubtreeRange(targetRes)

		curr := parent
		for r := uint8(5); r < targetRes; r++ {
			curr, _ = curr.Child(r % 9)
		}

		childGeMin := curr.High > minB.High || (curr.High == minB.High && curr.Low >= minB.Low)
		childLeMax := curr.High < maxB.High || (curr.High == maxB.High && curr.Low <= maxB.Low)

		if !childGeMin || !childLeMax {
			t.Fatalf("Descendant %s fell outside 128-bit SubtreeRange [%s, %s]", curr.Hex(), minB.Hex(), maxB.Hex())
		}
	})
}
