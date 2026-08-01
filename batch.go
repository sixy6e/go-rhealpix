package rhealpix

import (
	"fmt"
	"runtime"
	"sync"
)

// ForwardBatchParallel64 converts flat arrays of Lons/Lats into CellID64 keys concurrently.
func ForwardBatchParallel64(el *Ellipsoid, lons, lats []float64, targetRes uint8, out []CellID64) error {
	n := len(lons)
	if n == 0 {
		return nil
	}
	if len(lats) < n || len(out) < n {
		return fmt.Errorf("input/output slice length mismatch: lons(%d), lats(%d), out(%d)", n, len(lats), len(out))
	}

	numWorkers := runtime.GOMAXPROCS(0)
	if numWorkers > n {
		numWorkers = n
	}
	chunkSize := (n + numWorkers - 1) / numWorkers

	var wg sync.WaitGroup
	var errOnce sync.Once
	var firstErr error

	for w := 0; w < numWorkers; w++ {
		start := w * chunkSize
		end := start + chunkSize
		if end > n {
			end = n
		}
		if start >= end {
			continue
		}

		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()

			// compiler BCE (Bounds Check Elimination) hints
			_ = lons[e-1]
			_ = lats[e-1]
			_ = out[e-1]

			for i := s; i < e; i++ {
				id, err := ForwardTransform64(el, lons[i], lats[i], targetRes)
				if err != nil {
					errOnce.Do(func() { firstErr = err })
					return
				}
				out[i] = id
			}
		}(start, end)
	}

	wg.Wait()
	return firstErr
}

// ForwardBatchParallel128 converts flat arrays of Lons/Lats into CellID128 keys concurrently.
func ForwardBatchParallel128(el *Ellipsoid, lons, lats []float64, targetRes uint8, out []CellID128) error {
	n := len(lons)
	if n == 0 {
		return nil
	}
	if len(lats) < n || len(out) < n {
		return fmt.Errorf("input/output slice length mismatch: lons(%d), lats(%d), out(%d)", n, len(lats), len(out))
	}

	numWorkers := runtime.GOMAXPROCS(0)
	if numWorkers > n {
		numWorkers = n
	}
	chunkSize := (n + numWorkers - 1) / numWorkers

	var wg sync.WaitGroup
	var errOnce sync.Once
	var firstErr error

	for w := 0; w < numWorkers; w++ {
		start := w * chunkSize
		end := start + chunkSize
		if end > n {
			end = n
		}
		if start >= end {
			continue
		}

		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()

			// compiler BCE hints
			_ = lons[e-1]
			_ = lats[e-1]
			_ = out[e-1]

			for i := s; i < e; i++ {
				id, err := ForwardTransform128(el, lons[i], lats[i], targetRes)
				if err != nil {
					errOnce.Do(func() { firstErr = err })
					return
				}
				out[i] = id
			}
		}(start, end)
	}

	wg.Wait()
	return firstErr
}
