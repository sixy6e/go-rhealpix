package rhealpix_test

import (
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
)

// Test Generic Neighbours[T] directly
func TestNeighboursGeneric(t *testing.T) {
	t.Run("CellID64", func(t *testing.T) {
		cell, err := rhealpix.ParseCellID64("Q4")
		if err != nil {
			t.Fatalf("parse error: %v", err)
		}

		neighbours, err := rhealpix.Neighbours(cell)
		if err != nil {
			t.Fatalf("neighbours error: %v", err)
		}

		if len(neighbours.Slice()) != 8 {
			t.Errorf("expected 8 neighbours, got %d", len(neighbours.Slice()))
		}
	})

	t.Run("CellID128", func(t *testing.T) {
		cell, err := rhealpix.ParseCellID128("Q4")
		if err != nil {
			t.Fatalf("parse error: %v", err)
		}

		neighbours, err := rhealpix.Neighbours(cell)
		if err != nil {
			t.Fatalf("neighbours error: %v", err)
		}

		if len(neighbours.Slice()) != 8 {
			t.Errorf("expected 8 neighbours, got %d", len(neighbours.Slice()))
		}
	})
}

// Test explicit type wrappers (Neighbours64 & Neighbours128)
func TestNeighboursExplicitWrappers(t *testing.T) {
	t.Run("Neighbours64", func(t *testing.T) {
		cell, err := rhealpix.ParseCellID64("Q4")
		if err != nil {
			t.Fatalf("parse error: %v", err)
		}

		neighbours, err := rhealpix.Neighbours64(cell)
		if err != nil {
			t.Fatalf("Neighbours64 error: %v", err)
		}

		if len(neighbours.Slice()) != 8 {
			t.Errorf("expected 8 neighbours from Neighbours64, got %d", len(neighbours.Slice()))
		}
	})

	t.Run("Neighbours128", func(t *testing.T) {
		cell, err := rhealpix.ParseCellID128("Q4")
		if err != nil {
			t.Fatalf("parse error: %v", err)
		}

		neighbours, err := rhealpix.Neighbours128(cell)
		if err != nil {
			t.Fatalf("Neighbours128 error: %v", err)
		}

		if len(neighbours.Slice()) != 8 {
			t.Errorf("expected 8 neighbours from Neighbours128, got %d", len(neighbours.Slice()))
		}
	})
}

// Test topological edge wrapping across facets
func TestNeighboursFacetBoundaryWrapping(t *testing.T) {
	// West edge of Facet O (1) wraps to East edge of Facet R (4)
	cell, err := rhealpix.ParseCellID64("O0")
	if err != nil {
		t.Fatalf("parse error: %v", err)
	}

	neighbours, err := rhealpix.Neighbours(cell)
	if err != nil {
		t.Fatalf("neighbours error: %v", err)
	}

	if neighbours.W.Facet() != 4 {
		t.Errorf("expected West neighbour of O0 to wrap to Facet 4 (R), got Facet %d", neighbours.W.Facet())
	}
}
