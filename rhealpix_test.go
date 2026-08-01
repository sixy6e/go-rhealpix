package rhealpix_test

import (
	"math/rand"
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
)

// --- Benchmarks ---

func BenchmarkForwardBatchParallel64_1M(b *testing.B) {
	el := rhealpix.NewWGS84()
	numPoints := 1_000_000

	lons := make([]float64, numPoints)
	lats := make([]float64, numPoints)
	out := make([]rhealpix.CellID64, numPoints)

	rng := rand.New(rand.NewSource(42))
	for i := 0; i < numPoints; i++ {
		lons[i] = (rng.Float64() * 360.0) - 180.0
		lats[i] = (rng.Float64() * 180.0) - 90.0
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = rhealpix.ForwardBatchParallel64(el, lons, lats, 12, out)
	}
}

func BenchmarkForwardBatchParallel128_1M(b *testing.B) {
	el := rhealpix.NewWGS84()
	numPoints := 1_000_000

	lons := make([]float64, numPoints)
	lats := make([]float64, numPoints)
	out := make([]rhealpix.CellID128, numPoints)

	rng := rand.New(rand.NewSource(42))
	for i := 0; i < numPoints; i++ {
		lons[i] = (rng.Float64() * 360.0) - 180.0
		lats[i] = (rng.Float64() * 180.0) - 90.0
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = rhealpix.ForwardBatchParallel128(el, lons, lats, 24, out)
	}
}

// --- Unit Tests ---

func TestSUIDRoundTrip64(t *testing.T) {
	tests := []string{"Q", "Q0", "Q012", "N85", "S000"}

	for _, original := range tests {
		cell, err := rhealpix.ParseCellID64(original)
		if err != nil {
			t.Fatalf("ParseCellID64(%q) unexpected error: %v", original, err)
		}

		got := cell.String()
		if got != original {
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

		got := cell.String()
		if got != original {
			t.Errorf("128-bit SUID mismatch: got %q, want %q", got, original)
		}
	}
}

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

	t.Logf("Parent Cell: %s | Hex: %s | Bin: %s", cell.String(), cell.Hex(), cell.Binary())
	t.Logf("Min Bound:   Hex: %s | Bin: %s", minBound.Hex(), minBound.Binary())
	t.Logf("Max Bound:   Hex: %s | Bin: %s", maxBound.Hex(), maxBound.Binary())

	if minBound.Uint64() >= maxBound.Uint64() {
		t.Errorf("invalid SubtreeRange bounds: min %s >= max %s", minBound.Hex(), maxBound.Hex())
	}

	// test child cell
	childCell, err := rhealpix.ParseCellID64("Q0128")
	if err != nil {
		t.Fatalf("ParseCellID64 child error: %v", err)
	}

	t.Logf("Child Cell:  %s | Hex: %s | Res: %d", childCell.String(), childCell.Hex(), childCell.Resolution())

	if childCell.Uint64() < minBound.Uint64() || childCell.Uint64() > maxBound.Uint64() {
		t.Errorf(
			"child cell %s (%s) fell outside range [%s, %s]",
			childCell.String(), childCell.Hex(), minBound.Hex(), maxBound.Hex(),
		)
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

func TestCompact64(t *testing.T) {
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
