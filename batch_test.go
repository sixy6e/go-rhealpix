package rhealpix_test

import (
	"math/rand"
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
)

// --- Benchmarks ---

func BenchmarkForwardBatchParallel64_1M(b *testing.B) {
	targetRes := 12
	el := rhealpix.NewWGS84()
	d := rhealpix.NewRhpDGGS(
		rhealpix.WithEllipsoid(el),
		rhealpix.WithNSide(3),
		rhealpix.WithNorthSquare(0),
		rhealpix.WithSouthSquare(0),
		rhealpix.WithMaxResolution(uint8(targetRes)),
	)
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
		_ = rhealpix.ForwardBatchParallel64(d, lons, lats, uint8(targetRes), out)
	}
}

func BenchmarkForwardBatchParallel128_1M(b *testing.B) {
	targetRes := 24
	el := rhealpix.NewWGS84()
	d := rhealpix.NewRhpDGGS(
		rhealpix.WithEllipsoid(el),
		rhealpix.WithNSide(3),
		rhealpix.WithNorthSquare(0),
		rhealpix.WithSouthSquare(0),
		rhealpix.WithMaxResolution(uint8(targetRes)),
	)
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
		_ = rhealpix.ForwardBatchParallel128(d, lons, lats, uint8(targetRes), out)
	}
}
