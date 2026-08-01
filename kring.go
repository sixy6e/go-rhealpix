package rhealpix

import "fmt"

// KRing returns all cells within k steps of the origin cell (including origin).
// k=0 returns [origin]
// k=1 returns origin + 8 immediate neighbors (9 total)
// k=2 returns origin + 1st ring + 2nd ring, etc.
func KRing[T interface {
	CellID
	Cell[T]
}](origin T, k int) ([]T, error) {
	if k < 0 {
		return nil, fmt.Errorf("k must be non-negative")
	}
	if origin.IsZero() {
		return nil, fmt.Errorf("zero origin cell")
	}

	visited := make(map[T]bool)
	visited[origin] = true

	currentRing := []T{origin}

	for step := 0; step < k; step++ {
		var nextRing []T
		for _, cell := range currentRing {
			neighbors, err := Neighbours(cell)
			if err != nil {
				continue
			}

			for _, n := range neighbors.Slice() {
				if !visited[n] {
					visited[n] = true
					nextRing = append(nextRing, n)
				}
			}
		}
		currentRing = nextRing
	}

	// extract all collected unique cells
	result := make([]T, 0, len(visited))
	for c := range visited {
		result = append(result, c)
	}

	return result, nil
}
