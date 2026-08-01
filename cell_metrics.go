package rhealpix

import (
	"math"
)

// AreaErrorEntry holds the theoretical area and floating-point error bounds for a given resolution.
type AreaErrorEntry struct {
	CellAreaM2   float64 // Theoretical equal area per cell in square meters
	AbsTolerance float64 // Absolute tolerance (m²) for float comparisons
	RelTolerance float64 // Relative tolerance (dimensionless)
}

// CellWidth calculates the planar cell width at a given resolution.
// In rHEALPix (Nside = 3), planar width scales down by 3^(-res).
func (el *Ellipsoid) CellWidth(res uint8) float64 {
	return el.RA * (math.Pi / 2.0) * math.Pow(3.0, -float64(res))
}

// CellArea calculates the cell area in meters squared at a given resolution.
// If plane is true, it returns the local planar cell area.
// If plane is false, it returns the true authalic ellipsoidal surface area.
func (el *Ellipsoid) CellArea(res uint8, plane bool) float64 {
	w := el.CellWidth(res)
	planarArea := w * w
	if plane {
		return planarArea
	}
	// scale planar area to true authalic sphere area: 8 / (3 * pi) * planarArea
	return (8.0 / (3.0 * math.Pi)) * planarArea
}

// AreaErrorBudget returns an analytical error budget for cell area equality testing
// from resolution 0 up to maxRes (e.g., MaxResolution64 or MaxResolution128).
//
// Relative tolerance is bounded by 10 * machine epsilon (~2.22e-15).
func (el *Ellipsoid) AreaErrorBudget(maxRes uint8) map[uint8]AreaErrorEntry {
	// machine epsilon for IEEE 754 float64 is 2^(-52) ~ 2.220446049250313e-16
	const float64Epsilon = 1.0 / (1 << 52)
	relTol := 10.0 * float64Epsilon

	budget := make(map[uint8]AreaErrorEntry, maxRes+1)
	for r := uint8(0); r <= maxRes; r++ {
		area := el.CellArea(r, false)
		budget[r] = AreaErrorEntry{
			CellAreaM2:   area,
			AbsTolerance: area * relTol,
			RelTolerance: relTol,
		}
	}
	return budget
}
