package rhealpix

import "math"

const (
	// WGS84 Default Constants
	WGS84_A = 6378137.0
	WGS84_F = 1.0 / 298.257223563

	// Authalic sphere boundary constants (PROJ parity)
	AuthCapSinBeta0 = 2.0 / 3.0
	AuthCapBeta0    = 0.7297276562269663 // math.Asin(2.0 / 3.0)
	AuthCapBetaDeg  = 41.810314895778596 // asin(2/3) in degrees
)

// Ellipsoid holds pre-computed parameters for high-throughput authalic math.
type Ellipsoid struct {
	A, B, F, E, E2 float64
	OneMinusE2     float64

	// Third flattening n = (a-b)/(a+b)
	N float64

	// Lon0 is the central meridian offset in radians (default: 0.0)
	Lon0 float64

	// Total authalic parameter at pole
	QBar float64

	// Authalic radius (R_A)
	RA float64

	// Phi0 is the geodetic cap boundary latitude in radians (~41.93785° / 0.73196 rad for WGS84)
	// separating the equatorial region from the polar caps.
	// FYI: ~41.81031° on a sphere and ~41.93785° on an ellipsoid.
	// i.e. converting authalic latitude to geodectic latitude
	Phi0 float64

	// Standard model name identifier (e.g., "WGS84", "GRS80", "SPHERE")
	Name string

	// Pre-computed Authalic Latitude (beta_0) at the cap boundary
	// Beta0 = arcsin(2/3) ~ 0.729727658 rad (~41.81031 deg)
	// Useful for direct authalic space tests before converting to geodetic Phi0.
	Beta0 float64

	// Power series coefficients used by AuthLat
	authCoeffsForward [6]float64
	authCoeffsInverse [6]float64
}

// EllipsoidOption defines a functional option for configuring an Ellipsoid.
type EllipsoidOption func(*ellipsoidConfig)

type ellipsoidConfig struct {
	a    float64
	f    float64
	lon0 float64
}

// WithSemiMajorAxis sets the semi-major axis (a) in metres.
func WithSemiMajorAxis(a float64) EllipsoidOption {
	return func(cfg *ellipsoidConfig) {
		cfg.a = a
	}
}

// WithFlattening sets the flattening (f).
func WithFlattening(f float64) EllipsoidOption {
	return func(cfg *ellipsoidConfig) {
		cfg.f = f
	}
}

// WithLon0 sets the central meridian rotation offset in radians.
// For example, 129° East (129.0 * math.Pi / 180.0) re-centres the grid.
// e.g. lon_0 = 129° (or 174° - 45°) rotates the DGGS about the polar axis
// so that New Zealand sits in the middle of an equatorial facet.
func WithLon0(lon0Rad float64) EllipsoidOption {
	return func(cfg *ellipsoidConfig) {
		cfg.lon0 = lon0Rad
	}
}

// NewEllipsoid constructs an Ellipsoid using functional options.
// Defaults to WGS84 with lon0 = 0.0 if options are omitted.
func NewEllipsoid(opts ...EllipsoidOption) *Ellipsoid {
	cfg := &ellipsoidConfig{
		a:    WGS84_A,
		f:    WGS84_F,
		lon0: 0.0,
	}

	for _, opt := range opts {
		if opt != nil {
			opt(cfg)
		}
	}

	a := cfg.a
	f := cfg.f
	b := a * (1.0 - f)
	e2 := (a*a - b*b) / (a * a)
	e := math.Sqrt(e2)
	oneMinusE2 := 1.0 - e2

	// third flattening n
	n := f / (2.0 - f)

	// authalic qBar parameter at pole
	// qBar = (1 - e^2) * ( (sin(pi/2)/(1 - e^2 sin^2(pi/2))) - (1/(2e))*ln(...) )
	// qBar := oneMinusE2 * (1.0/oneMinusE2 - (1.0/(2.0*e))*math.Log((1.0-e)/(1.0+e)))
	// fixed QBar using PROJ/Snyder atanh for RA authalic radius (equivalent to PROJ and rhealpixdggs-py)
	qBar := 1.0 + (oneMinusE2 * math.Atanh(e) / e)
	rA := a * math.Sqrt(qBar/2.0)

	el := &Ellipsoid{
		A:                 a,
		B:                 b,
		F:                 f,
		E:                 e,
		E2:                e2,
		OneMinusE2:        oneMinusE2,
		N:                 n,
		Lon0:              cfg.lon0,
		QBar:              qBar,
		RA:                rA,
		authCoeffsForward: authLatCoefficients(n, false),
		authCoeffsInverse: authLatCoefficients(n, true),
	}

	// pre-compute the exact geodetic cap transition latitude (phi_0) for this ellipsoid:
	// phi_0 = AuthLatInverse( asin(2/3) )
	authalicCapBeta := math.Asin(2.0 / 3.0)
	el.Beta0 = authalicCapBeta
	el.Phi0 = el.AuthLatInverse(authalicCapBeta)

	return el
}

// NewWGS84 creates a default WGS84 Ellipsoid instance.
func NewWGS84() *Ellipsoid {
	return NewEllipsoid()
}

// NormaliseLongitude adjusts longitude (in radians) relative to Lon0 into [-pi, pi).
func (e *Ellipsoid) NormaliseLongitude(lonRad float64) float64 {
	return NormaliseLongitudeRad(lonRad - e.Lon0)
}

// AuthLat computes authalic latitude (beta in radians) from geodetic latitude (phi in radians).
// Uses a direct formula which isn't as accurate as the power series coefficients method.
func (el *Ellipsoid) AuthLatDirect(phi float64) float64 {
	sinPhi := math.Sin(phi)

	// quick check for polar extremities
	if sinPhi >= 0.99999999 {
		return HalfPi
	}
	if sinPhi <= -0.99999999 {
		return -HalfPi
	}

	// handle spherical case (E == 0) safely
	if el.E == 0 {
		return phi
	}

	sin2 := sinPhi * sinPhi
	q := el.OneMinusE2 * (sinPhi/(1.0-el.E2*sin2) - (1.0/(2.0*el.E))*math.Log((1.0-el.E*sinPhi)/(1.0+el.E*sinPhi)))
	sinBeta := q / el.QBar

	// clamp ratio to [-1.0, 1.0] to protect math.Asin from rounding precision overflow
	if sinBeta > 1.0 {
		sinBeta = 1.0
	} else if sinBeta < -1.0 {
		sinBeta = -1.0
	}

	return math.Asin(sinBeta)
}

func (e *Ellipsoid) GeodeticToAuthalic(phi float64) float64 {
	return e.AuthLat(phi)
}

// authLatCoefficients computes the 6 Fourier coefficients matching rhealpixdggs-py verbatim.
// source: https://arxiv.org/pdf/2212.05818v2
// The power-series coefficients of sin(2 phi), ..., sin(12 phi) used by (el *Ellipsoid) AuthLat
// for third flattening n. Equations A19 (forward) and A20 (inverse).
// Ideal for when ellipsoidal flattenings are small, or according to the paper, (f <= 1/150).
// For reference, the flattening for WGS84 is 1.0 / 298.257223563
func authLatCoefficients(n float64, inverse bool) [6]float64 {
	n2 := n * n
	n3 := n2 * n
	n4 := n3 * n
	n5 := n4 * n
	n6 := n5 * n

	if !inverse {
		// forward series: phi -> beta (Eq. A19)
		c2 := n * (-4.0/3.0 + n*(-4.0/45.0+n*(88.0/315.0+n*(538.0/4725.0+n*(20824.0/467775.0+n*(-44732.0/2837835.0))))))
		c4 := n2 * (34.0/45.0 + n*(8.0/105.0+n*(-2482.0/14175.0+n*(-37192.0/467775.0+n*(-12467764.0/212837625.0)))))
		c6 := n3 * (-1532.0/2835.0 + n*(-898.0/14175.0+n*(54968.0/467775.0+n*(100320856.0/1915538625.0))))
		c8 := n4 * (6007.0/14175.0 + n*(24496.0/467775.0+n*(-5884124.0/70945875.0)))
		c10 := n5 * (-23356.0/66825.0 + n*(-839792.0/19348875.0))
		c12 := n6 * (570284222.0 / 1915538625.0)

		return [6]float64{c2, c4, c6, c8, c10, c12}
	}

	// inverse series: beta -> phi (Eq. A20)
	c2 := n * (4.0/3.0 + n*(4.0/45.0+n*(-16.0/35.0+n*(-2582.0/14175.0+n*(60136.0/467775.0+n*(28112932.0/212837625.0))))))
	c4 := n2 * (46.0/45.0 + n*(152.0/945.0+n*(-11966.0/14175.0+n*(-21016.0/51975.0+n*(251310128.0/638512875.0)))))
	c6 := n3 * (3044.0/2835.0 + n*(3802.0/14175.0+n*(-94388.0/66825.0+n*(-8797648.0/10945935.0))))
	c8 := n4 * (6059.0/4725.0 + n*(41072.0/93555.0+n*(-1472637812.0/638512875.0)))
	c10 := n5 * (768272.0/467775.0 + n*(455935736.0/638512875.0))
	c12 := n6 * (4210684958.0 / 1915538625.0)

	return [6]float64{c2, c4, c6, c8, c10, c12}
}

// AuthLat converts geodetic latitude (phi) -> authalic latitude (beta) in radians.
func (el *Ellipsoid) AuthLat(phi float64) float64 {
	if el.E == 0 {
		return phi
	}
	c := el.authCoeffsForward
	return phi + (c[0]*math.Sin(2*phi) +
		c[1]*math.Sin(4*phi) +
		c[2]*math.Sin(6*phi) +
		c[3]*math.Sin(8*phi) +
		c[4]*math.Sin(10*phi) +
		c[5]*math.Sin(12*phi))
}

// AuthLatInverse converts authalic latitude (beta) -> geodetic latitude (phi) in radians.
func (el *Ellipsoid) AuthLatInverse(beta float64) float64 {
	if el.E == 0 {
		return beta
	}
	c := el.authCoeffsInverse
	return beta + (c[0]*math.Sin(2*beta) +
		c[1]*math.Sin(4*beta) +
		c[2]*math.Sin(6*beta) +
		c[3]*math.Sin(8*beta) +
		c[4]*math.Sin(10*beta) +
		c[5]*math.Sin(12*beta))
}
