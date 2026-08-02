package rhealpix

import "math"

const (
	// WGS84 Constants
	WGS84_A = 6378137.0
	WGS84_F = 1.0 / 298.257223563
)

// Ellipsoid holds pre-computed parameters for high-throughput authalic math.
type Ellipsoid struct {
	A, B, E, E2 float64
	OneMinusE2  float64

	// Total authalic parameter at pole
	QBar float64

	// Authalic radius (R_A)
	RA float64
}

// NewWGS84 creates a pre-computed WGS84 Ellipsoid instance.
func NewWGS84() *Ellipsoid {
	a := WGS84_A
	f := WGS84_F
	b := a * (1.0 - f)
	e2 := (a*a - b*b) / (a * a)
	e := math.Sqrt(e2)
	oneMinusE2 := 1.0 - e2

	// qBar = (1 - e^2) * ( (sin(pi/2)/(1 - e^2 sin^2(pi/2))) - (1/(2e))*ln(...) )
	qBar := oneMinusE2 * (1.0/oneMinusE2 - (1.0/(2.0*e))*math.Log((1.0-e)/(1.0+e)))
	rA := a * math.Sqrt(qBar/2.0)

	return &Ellipsoid{
		A:          a,
		B:          b,
		E:          e,
		E2:         e2,
		OneMinusE2: oneMinusE2,
		QBar:       qBar,
		RA:         rA,
	}
}
