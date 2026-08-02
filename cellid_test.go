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

func TestSubtreeRange64(t *testing.T) {
	cell, err := rhealpix.ParseCellID64("Q012")
	if err != nil {
		t.Fatalf("ParseCellID64 unexpected error: %v", err)
	}

	minBound, maxBound := cell.SubtreeRange()

	if minBound.Uint64() >= maxBound.Uint64() {
		t.Errorf("invalid SubtreeRange bounds: min %s >= max %s", minBound.Hex(), maxBound.Hex())
	}

	// test child cell
	childCell, err := rhealpix.ParseCellID64("Q0128")
	if err != nil {
		t.Fatalf("ParseCellID64 child error: %v", err)
	}

	if childCell.Uint64() < minBound.Uint64() || childCell.Uint64() > maxBound.Uint64() {
		t.Errorf("child cell %s fell outside range [%s, %s]", childCell.String(), minBound.Hex(), maxBound.Hex())
	}
}

func TestSubtreeRange128(t *testing.T) {
	cell, _ := rhealpix.ParseCellID128("Q01234567890123")
	minBound, maxBound := cell.SubtreeRange()

	if minBound.High > maxBound.High || (minBound.High == maxBound.High && minBound.Low >= maxBound.Low) {
		t.Errorf("invalid SubtreeRange128 bounds: min %s >= max %s", minBound.Hex(), maxBound.Hex())
	}

	childCell, _ := rhealpix.ParseCellID128("Q012345678901238")
	if childCell.High < minBound.High || childCell.High > maxBound.High {
		t.Errorf("child cell High word %x fell outside [%x, %x]", childCell.High, minBound.High, maxBound.High)
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
	base := "P012345678901234"
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
