package rhealpix

import (
	"math"
	"testing"
)

func TestComputeCellPath_GroundTruth(t *testing.T) {
	el := NewWGS84()

	// map numerical facet index (0..5) back to canonical character (N, O, P, Q, R, S)
	facetChars := [6]byte{'N', 'O', 'P', 'Q', 'R', 'S'}

	tests := []struct {
		name     string
		lon, lat float64
		res      uint8
		wantSUID string
	}{
		{
			name:     "North Cap (N)",
			lon:      0.0,
			lat:      75.0,
			res:      15,
			wantSUID: "N422244442446644",
		},
		{
			name:     "Facet O",
			lon:      -135.0,
			lat:      0.0,
			res:      15,
			wantSUID: "O444444444444444",
		},
		{
			name:     "Facet P",
			lon:      -45.0,
			lat:      10.0,
			res:      15,
			wantSUID: "P411777777144741",
		},
		{
			name:     "Facet Q",
			lon:      45.0,
			lat:      10.0,
			res:      15,
			wantSUID: "Q411777777144741",
		},
		{
			name:     "Facet R - Pilbara Benchmark",
			lon:      119.6592,
			lat:      -22.0742,
			res:      15,
			wantSUID: "R652226053773451",
		},
		{
			name:     "South Cap (S)",
			lon:      0.0,
			lat:      -75.0,
			res:      15,
			wantSUID: "S488844448440044",
		},
		{
			name:     "Equator Seam",
			lon:      10.0,
			lat:      0.0,
			res:      15,
			wantSUID: "Q343333333333333",
		},
		{
			name:     "Facet Boundary O/P",
			lon:      -90.0,
			lat:      0.0,
			res:      15,
			wantSUID: "P333333333333333",
		},
		{
			name:     "North Pole Apex",
			lon:      0.0,
			lat:      90.0,
			res:      15,
			wantSUID: "N444444444444444",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			facetIdx, path, err := computeCellPath(el, tt.lon, tt.lat, tt.res, 15)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}

			// format SUID string (e.g. "R652226053773451")
			gotSUID := string(facetChars[facetIdx])
			for _, digit := range path {
				gotSUID += string(rune('0' + digit))
			}

			if gotSUID != tt.wantSUID {
				t.Errorf("computeCellPath(%f, %f) =\n  got:  %s\n  want: %s", tt.lon, tt.lat, gotSUID, tt.wantSUID)
			}
		})
	}
}

func TestForwardProject_GroundTruth(t *testing.T) {
	el := NewWGS84()

	// benchmarks sourced directly from rhealpixdggs-py / pyproj (WGS84)
	tests := []struct {
		name     string
		lon, lat float64
		wantX    float64
		wantY    float64
	}{
		{
			name:  "North Cap (N)",
			lon:   0.0,
			lat:   75.0,
			wantX: -13404697.863162834,
			wantY: 11614188.83126379,
		},
		{
			name:  "Facet O",
			lon:   -135.0,
			lat:   0.0,
			wantX: -15011332.016655976,
			wantY: 0.0,
		},
		{
			name:  "Facet P",
			lon:   -45.0,
			lat:   10.0,
			wantX: -5003777.338885325,
			wantY: 1297694.0300918117,
		},
		{
			name:  "Facet Q",
			lon:   45.0,
			lat:   10.0,
			wantX: 5003777.338885325,
			wantY: 1297694.0300918117,
		},
		{
			name:  "Facet R - Pilbara Benchmark",
			lon:   119.6592,
			lat:   -22.0742,
			wantX: 13305510.963314377,
			wantY: -2809845.1733985,
		},
		{
			name:  "South Cap (S)",
			lon:   0.0,
			lat:   -75.0,
			wantX: -13404697.863162834,
			wantY: -11614188.83126379,
		},
		{
			name:  "Equator Seam",
			lon:   10.0,
			lat:   0.0,
			wantX: 1111950.5197522945,
			wantY: 0.0,
		},
		{
			name:  "Facet Boundary O/P",
			lon:   -90.0,
			lat:   0.0,
			wantX: -10007554.67777065,
			wantY: 0.0,
		},
		{
			name:  "North Pole Apex",
			lon:   0.0,
			lat:   90.0,
			wantX: -15011332.016655976,
			wantY: 10007554.67777065,
		},
	}

	const toleranceMeters = 1e-5 // 0.01 mm tolerance

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotX, gotY := el.ForwardProject(tt.lon, tt.lat)

			diffX := math.Abs(gotX - tt.wantX)
			diffY := math.Abs(gotY - tt.wantY)

			if diffX > toleranceMeters || diffY > toleranceMeters {
				t.Errorf("ForwardProject(%f, %f) =\n  got:  (%18.8f, %18.8f)\n  want: (%18.8f, %18.8f)\n  diff: (%.8f m, %.8f m)",
					tt.lon, tt.lat, gotX, gotY, tt.wantX, tt.wantY, diffX, diffY)
			}
		})
	}
}

func TestProjectionRoundTrip(t *testing.T) {
	el := NewWGS84()

	tests := []struct {
		name     string
		lon, lat float64
	}{
		{"North Cap (N)", 0.0, 75.0},
		{"Facet O", -135.0, 0.0},
		{"Facet P", -45.0, 10.0},
		{"Facet Q", 45.0, 10.0},
		{"Facet R - Pilbara", 119.6592, -22.0742},
		{"South Cap (S)", 0.0, -75.0},
		{"Equator Seam", 10.0, 0.0},
		{"Facet Boundary O/P", -90.0, 0.0},
		{"High Precision Point", 151.2093, -33.8688},
	}

	const tolDegrees = 1e-9 // nanodegree accuracy (~0.1 mm on ground)

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			xMeters, yMeters := el.ForwardProject(tt.lon, tt.lat)
			gotLon, gotLat := el.InverseProject(xMeters, yMeters)

			diffLon := math.Abs(gotLon - tt.lon)
			diffLat := math.Abs(gotLat - tt.lat)

			if diffLon > tolDegrees || diffLat > tolDegrees {
				t.Errorf("InverseProject(ForwardProject(%f, %f)) =\n  got:  (%18.9f, %18.9f)\n  want: (%18.9f, %18.9f)\n  diff: (%.9f deg, %.9f deg)",
					tt.lon, tt.lat, gotLon, gotLat, tt.lon, tt.lat, diffLon, diffLat)
			}
		})
	}
}
