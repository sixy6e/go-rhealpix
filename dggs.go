package rhealpix

import (
	"fmt"
	"math"
)

// Standard planar domain constants for an authalic sphere of radius R_A = 1
// i.e. working in a dimensionless unit coordinate system where the radius
// is normalised to 1.0.
const (
	Pi         = math.Pi
	TwoPi      = 2.0 * math.Pi
	HalfPi     = math.Pi / 2.0
	QuarterPi  = math.Pi / 4.0
	ThreeQtrPi = 3.0 * math.Pi / 4.0

	// Default polar cap alignment positions (0 = directly above/below Facet O)
	DefaultNorthSquare = 0
	DefaultSouthSquare = 0

	// N-sides relates to either a 2x2 cell configuration or a 3x3 cell configuration
	DefaultNsides = 3

	// MaxSupportedResolution represents the maximum subdivision depth supported
	// by the 128-bit cell encoding structure (Level 30 = nanometre precision).
	MaxSupportedResolution uint8 = 30

	// Angle conversion constants
	RadToDeg = 180.0 / math.Pi
	DegToRad = math.Pi / 180.0

	// EPS_DOMAIN: Used for floating-point domain containment (pnpoly checks & interval bounds)
	// Matches PROJ's #define EPS 1e-15
	EPS_DOMAIN = 1e-15

	// EPS_SEAM_RAD: Radian tolerance for snapping geodetic latitude to the authalic seam
	// Prevents double-precision roundtrip noise between AuthLat and AuthLatInverse (~1cm on Earth)
	EPS_SEAM_RAD = 1e-11

	// EPS_SEAM_DEG: Degree tolerance for snapping geographic split points (~1cm on Earth)
	EPS_SEAM_DEG = 1e-7

	// EPS_CLIP_BOUND: Planar bounding box inflation tolerance
	// Prevents orb/clip from stripping interior vertices sitting directly on facet boundaries
	EPS_CLIP_BOUND = 1e-12

	// Micro-fuzz margin (1e-10) to ensure vertices sitting directly on
	// projection boundaries pass domain validation despite floating-point noise.
	EPS_IMAGE = 1e-10

	// Epsilon for planar radians and absorbs float64 arithmetic noise (~0.06mm on Earth)
	EPS_PLANAR = 1e-11
)

// BaseFacetBounds defines the exact global planar bounding box for each base facet [0..5]
// in authalic radius R_A = 1 unit space.
type BaseFacetBounds struct {
	ULX, ULY float64 // Upper-Left corner (x, y)
	Width    float64 // Planar width
	Height   float64 // Planar height
}

// RhpDGGS encapsulates the full configuration, ellipsoid parametres,
// and layout bounds for an rHEALPix DGGS structure.
type RhpDGGS struct {
	ellipsoid   *Ellipsoid
	nSide       uint8 // Subdivision factor (e.g., 3 for nonary)
	northSquare uint8 // Alignment offset [0..3]
	southSquare uint8 // Alignment offset [0..3]

	maxArealResolution float64 // Targeted max area in m²
	maxResolution      uint8   // Derived max grid depth level

	baseFacetBounds [6]BaseFacetBounds

	// Area in m² of a Level 0 base facet (4 * Pi * R_A^2 / 6 = 2/3 * Pi * R_A^2)
	baseFacetArea float64

	// Planar extent width in metres at Level 0: W = R_A * (Pi / 2)
	// Eliminates redundant planar-to-metre scaling calculations in hot loops.
	baseFacetWidthMetres float64

	// Total number of cells at maxResolution: 6 * (nSide^(2 * maxResolution))
	// Critical for spatial range partitioning, key space validation, and TileDB dimension sizing.
	totalCellsAtMaxRes uint64 // or big.Int for level > 15
}

// Option defines a functional option for configuring RhpDGGS.
type RhpOption func(*RhpDGGS)

func WithEllipsoid(e *Ellipsoid) RhpOption {
	return func(r *RhpDGGS) {
		r.ellipsoid = e
	}
}

// WithNSide defines each planar cell to have nSide x nSide child cells.
// An integer of size at least 2.
func WithNSide(nSide uint8) RhpOption {
	return func(r *RhpDGGS) {
		if nSide < 2 {
			panic("nSide must be at least 2")
		}
		r.nSide = nSide
	}
}

// WithNorthSquare sets the north polar cap alignment. 0 is directly above Facet O.
func WithNorthSquare(square uint8) RhpOption {
	return func(r *RhpDGGS) {
		r.northSquare = square % 4 // 0, 1, 2 or 3
	}
}

// WithSouthSquare sets the south polar cap alignment. 0 is directly below Facet O.
func WithSouthSquare(square uint8) RhpOption {
	return func(r *RhpDGGS) {
		r.southSquare = square % 4 // 0, 1, 2 or 3
	}
}

// WithMaxArealResolution sets max ground area per cell (in m²) and computes maxResolution depth.
func WithMaxArealResolution(sqMetres float64) RhpOption {
	return func(r *RhpDGGS) {
		if sqMetres <= 0 {
			return
		}
		r.maxArealResolution = sqMetres
		r.maxResolution = CalculateMaxResolution(r.ellipsoid.RA, r.nSide, sqMetres)
	}
}

// WithMaxResolution set the target resolution for cell generation.
// e.g. 12 is ideal for 30m Landsat data, and 13 is ideal for 10m Sentinel-2 data.
func WithMaxResolution(level uint8) RhpOption {
	return func(r *RhpDGGS) {
		r.maxResolution = level
	}
}

// NewRhpDGGS constructs a fully initialised, immutable RhpDGGS context.
func NewRhpDGGS(opts ...RhpOption) *RhpDGGS {
	dggs := &RhpDGGS{
		ellipsoid:          NewWGS84(),
		nSide:              3, // default nonary (3x3)
		northSquare:        0,
		southSquare:        0,
		maxArealResolution: 1.0, // default 1 m²
	}

	for _, opt := range opts {
		opt(dggs)
	}

	// compute max resolution if not explicitly set
	if dggs.maxResolution == 0 && dggs.maxArealResolution > 0 {
		dggs.maxResolution = CalculateMaxResolution(dggs.ellipsoid.RA, dggs.nSide, dggs.maxArealResolution)
	}

	// pre-calculate fixed base facet bounds for instantaneous lookups
	for f := uint8(0); f < 6; f++ {
		dggs.baseFacetBounds[f] = GetBaseFacetBounds(f, int(dggs.northSquare), int(dggs.southSquare))
	}

	return dggs
}

// CalculateMaxResolution computes the maximum resolution level required to achieve
// cells with an area no larger than maxArealResolution (in m²).
func CalculateMaxResolution(ra float64, nSide uint8, maxArealResolution float64) uint8 {
	if maxArealResolution <= 0 {
		return 0
	}

	baseArea := ra * ra * (2.0 * math.Pi / 3.0)
	ratio := baseArea / maxArealResolution
	if ratio <= 1.0 {
		return 0
	}

	logN := math.Log(float64(nSide))
	val := math.Log(ratio) / (2.0 * logN)

	// epsilon prevents float64 floating-point precision noise from over-shooting level boundaries
	const eps = 1e-11
	level := math.Ceil(val - eps)

	if level <= 0 {
		return 0
	}

	res := uint8(level)
	if res > MaxSupportedResolution {
		return MaxSupportedResolution
	}

	return res
}

// GetBaseFacetBounds returns the exact planar extent for a base facet given N_square and S_square alignment.
func GetBaseFacetBounds(facet uint8, northSquare, southSquare int) BaseFacetBounds {
	switch facet {
	case 0: // Facet N (North Cap)
		ulX := -Pi + float64(northSquare)*HalfPi
		return BaseFacetBounds{ULX: ulX, ULY: ThreeQtrPi, Width: HalfPi, Height: HalfPi}
		// return BaseFacetBounds{ULX: ulX, ULY: ThreeQtrPi, Width: Pi, Height: HalfPi} // very wrong, just attempting to correct for triangular shaped overshoots and undershoots.

	case 1: // Facet O (Equatorial)
		return BaseFacetBounds{ULX: -Pi, ULY: QuarterPi, Width: HalfPi, Height: HalfPi}

	case 2: // Facet P (Equatorial)
		return BaseFacetBounds{ULX: -HalfPi, ULY: QuarterPi, Width: HalfPi, Height: HalfPi}

	case 3: // Facet Q (Equatorial)
		return BaseFacetBounds{ULX: 0.0, ULY: QuarterPi, Width: HalfPi, Height: HalfPi}

	case 4: // Facet R (Equatorial)
		return BaseFacetBounds{ULX: HalfPi, ULY: QuarterPi, Width: HalfPi, Height: HalfPi}

	case 5: // Facet S (South Cap)
		ulX := -Pi + float64(southSquare)*HalfPi
		return BaseFacetBounds{ULX: ulX, ULY: -QuarterPi, Width: HalfPi, Height: HalfPi}

	default:
		panic("invalid facet ID")
	}
}

func (d *RhpDGGS) GetEllipsoid() *Ellipsoid {
	return d.ellipsoid
}

func (d *RhpDGGS) Nside() uint8 {
	return d.nSide
}

// FacetBounds returns the pre-calculated planar extent bounds for a base facet [0..5].
func (d *RhpDGGS) FacetBounds(facet uint8) (BaseFacetBounds, error) {
	if facet > 5 {
		return BaseFacetBounds{}, fmt.Errorf("invalid facet ID: %d (must be 0..5)", facet)
	}
	return d.baseFacetBounds[facet], nil
}

// MustFacetBounds returns the pre-calculated planar extent bounds for a base facet [0..5].
// It panics if the facet ID is out of bounds [0..5].
func (d *RhpDGGS) MustFacetBounds(facet uint8) BaseFacetBounds {
	if facet > 5 {
		panic(fmt.Sprintf("invalid facet ID: %d (must be 0..5)", facet))
	}
	return d.baseFacetBounds[facet]
}

// --- Convenience Helpers on BaseFacetBounds ---

// MinX returns the western planar bound (ULX).
func (b BaseFacetBounds) MinX() float64 { return b.ULX }

// MaxX returns the eastern planar bound (ULX + Width).
func (b BaseFacetBounds) MaxX() float64 { return b.ULX + b.Width }

// MinY returns the southern planar bound (ULY - Height).
func (b BaseFacetBounds) MinY() float64 { return b.ULY - b.Height }

// MaxY returns the northern planar bound (ULY).
func (b BaseFacetBounds) MaxY() float64 { return b.ULY }

// Centroid returns the centre point [x, y] of the base facet in authalic planar space.
func (b BaseFacetBounds) Centroid() [2]float64 {
	return [2]float64{
		b.ULX + (b.Width / 2.0),
		b.ULY - (b.Height / 2.0),
	}
}

// RhpSpec represents the serializable configuration of an RhpDGGS instance.
type RhpSpec struct {
	// Grid Topology & Alignment
	NSide         uint8 `json:"nside"`
	NorthSquare   uint8 `json:"north_square"`
	SouthSquare   uint8 `json:"south_square"`
	MaxResolution uint8 `json:"max_resolution"`

	// Pre-computed Geometric Metrics
	BaseFacetArea        float64 `json:"base_facet_area_m2"`
	BaseFacetWidthMetres float64 `json:"base_facet_width_m"`

	// Ellipsoid Parameters
	EllipsoidName  string  `json:"ellipsoid_name,omitempty"`
	SemiMajorAxis  float64 `json:"semi_major_axis"`
	Flattening     float64 `json:"flattening"`
	Lon0           float64 `json:"lon0_rad"`
	AuthalicRadius float64 `json:"authalic_radius_m"`
	Phi0           float64 `json:"phi0_rad"`
}

// Spec returns the serializable grid parameters for metadata storage.
func (d *RhpDGGS) Spec() RhpSpec {
	el := d.GetEllipsoid()
	return RhpSpec{
		NSide:                d.nSide,
		NorthSquare:          d.northSquare,
		SouthSquare:          d.southSquare,
		MaxResolution:        d.maxResolution,
		BaseFacetArea:        d.baseFacetArea,
		BaseFacetWidthMetres: d.baseFacetWidthMetres,
		EllipsoidName:        el.Name,
		SemiMajorAxis:        el.A,
		Flattening:           el.F,
		Lon0:                 el.Lon0,
		AuthalicRadius:       el.RA,
		Phi0:                 el.Phi0,
	}
}
