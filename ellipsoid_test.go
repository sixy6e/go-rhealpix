package rhealpix_test

import (
	"math"
	"testing"

	rhealpix "github.com/sixy6e/go-rhealpix"
)

func TestWGS84AuthalicRadius(t *testing.T) {
	el := rhealpix.NewWGS84()

	// WGS84 authalic radius R_A is approximately 6,371,007.181 meters
	wantRA := 6371007.181
	if math.Abs(el.RA-wantRA) > 0.01 {
		t.Errorf("authalic radius mismatch: got %f, want ~%f", el.RA, wantRA)
	}
}
