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
