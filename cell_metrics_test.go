package rhealpix_test

import (
	"math"
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
)

func TestCellAreaInvariants(t *testing.T) {
	el := rhealpix.NewWGS84()
	budget := el.AreaErrorBudget(15)

	// total Earth Authalic Surface Area = 4 * pi * R_A^2
	totalEarthArea := 4.0 * math.Pi * el.RA * el.RA

	// at Res 0, there are 6 base facets. Each Res 0 cell should be 1/6th of total Earth area.
	res0Area := el.CellArea(0, false)
	expectedRes0Area := totalEarthArea / 6.0

	// use tolerance from the AreaErrorBudget for Res 0
	tol0 := budget[0].AbsTolerance
	if math.Abs(res0Area-expectedRes0Area) > tol0 {
		t.Errorf("Res 0 area mismatched beyond tolerance: got %f, expected %f (diff: %e, max tol: %e)",
			res0Area, expectedRes0Area, math.Abs(res0Area-expectedRes0Area), tol0)
	}

	// at Res 1, each Res 0 cell splits into 9 sub-cells.
	res1Area := el.CellArea(1, false)
	expectedRes1Area := res0Area / 9.0

	tol1 := budget[1].AbsTolerance
	if math.Abs(res1Area-expectedRes1Area) > tol1 {
		t.Errorf("Res 1 area mismatched beyond tolerance: got %f, expected %f (diff: %e, max tol: %e)",
			res1Area, expectedRes1Area, math.Abs(res1Area-expectedRes1Area), tol1)
	}

	// verify AreaErrorBudget structure
	if len(budget) != 16 {
		t.Errorf("expected 16 resolutions in budget, got %d", len(budget))
	}
}
