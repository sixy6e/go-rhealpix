package rhealpix

import (
	"fmt"
	"math"
)

// NormaliseLongitudeRad normalises any longitude in radians into [-pi, pi).
func NormaliseLongitudeRad(lonRad float64) float64 {
	return lonRad - math.Floor((lonRad+Pi)/TwoPi)*TwoPi
}

// NormaliseLonDeg normalises any longitude in degrees into [-180, 180).
func NormaliseLonDeg(lonDeg float64) float64 {
	return lonDeg - math.Floor((lonDeg+180.0)/360.0)*360.0
}

// ForwardProject converts geodetic longitude and latitude (in degrees)
// to planar rHEALPix projected coordinates (x, y in metres) on the ellipsoid.
func (d *RhpDGGS) ForwardProject(lonDeg, latDeg float64) (xMetres, yMetres float64, facet uint8) {
	lonRad := lonDeg * DegToRad
	latRad := latDeg * DegToRad

	// projected dimensionless planar coordinates on unit sphere
	xRhp, yRhp, facet := d.ForwardProjectRadians(lonRad, latRad)

	// scale by Authalic Radius (R_A) to get physical ground metres
	return xRhp * d.ellipsoid.RA, yRhp * d.ellipsoid.RA, facet
}

// ForwardProjectRadians takes lon and lat as radians and projects them to the
// rHEALPix projection. Returned coordinates are in radians relative to the facet they lie within.
func (d *RhpDGGS) ForwardProjectRadians(lonRad, latRad float64) (xRhp float64, yRhp float64, facet uint8) {
	lon := d.ellipsoid.NormaliseLongitude(lonRad)
	beta := d.ellipsoid.AuthLat(latRad)

	xRhp, yRhp = d.HealpixSphere(lon, beta)

	xRhp, yRhp, facet = d.CombineTrianglesWithFacet(xRhp, yRhp, false, lon, RegionAuto)
	return xRhp, yRhp, facet
}

// ForwardProjectForcedRegion projects a coordinate while enforcing a specific regional math branch
// (RegionNorthPolar, RegionEquatorial, RegionSouthPolar). This bypasses auto-routing to maintain
// topological integrity and seam alignment during polygon transformation.
func (d *RhpDGGS) ForwardProjectForcedRegion(lonRad, latRad float64, region CellRegion) (float64, float64, uint8) {
	lon := d.ellipsoid.NormaliseLongitude(lonRad)
	beta := d.ellipsoid.AuthLat(latRad)

	// seam check: if point sits on/near either seam, force exact Beta0 with matching sign
	if math.Abs(math.Abs(beta)-d.ellipsoid.Beta0) < EPS_SEAM_RAD {
		sign := 1.0
		if beta < 0 {
			sign = -1.0
		}

		switch region {
		case RegionNorthPolar:
			x, y := d.healpixSpherePolar(lon, d.ellipsoid.Beta0)
			return d.CombineTrianglesWithFacet(x, y, false, lon, region)

		case RegionSouthPolar:
			x, y := d.healpixSpherePolar(lon, -d.ellipsoid.Beta0)
			return d.CombineTrianglesWithFacet(x, y, false, lon, region)

		case RegionEquatorial:
			x, y := d.healpixSphereEquatorial(lon, sign*d.ellipsoid.Beta0)
			return d.CombineTrianglesWithFacet(x, y, false, lon, region)
		}
	}

	// unforced region (fallback for interior vertices)
	x, y := d.HealpixSphere(lon, beta)
	return d.CombineTrianglesWithFacet(x, y, false, lon, region)
}

// ForwardTransform64 converts Lon/Lat (degrees) into a 64-bit CellID64 (up to Res 14).
func (d *RhpDGGS) ForwardTransform64(lonDeg, latDeg float64, targetRes uint8) (CellID64, error) {
	facet, path, err := d.computeCellPath(lonDeg, latDeg, targetRes, MaxResolution64)
	if err != nil {
		return 0, err
	}
	return PackCellID64(facet, targetRes, path)
}

// ForwardTransform128 converts Lon/Lat (degrees) into a 128-bit CellID128 (up to Res 30).
func (d *RhpDGGS) ForwardTransform128(lonDeg, latDeg float64, targetRes uint8) (CellID128, error) {
	facet, path, err := d.computeCellPath(lonDeg, latDeg, targetRes, MaxResolution128)
	if err != nil {
		return CellID128{}, err
	}
	return PackCellID128(facet, targetRes, path)
}

// computeCellPath converts Lon/Lat degrees to base facet and sub-cell (e.g. 2x2, 3x3) digit path.
func (d *RhpDGGS) computeCellPath(lonDeg, latDeg float64, targetRes uint8, maxAllowedRes uint8) (uint8, []uint8, error) {
	if targetRes > maxAllowedRes {
		return 0, nil, ValidateCellParams(0, targetRes, maxAllowedRes, int(targetRes))
	}

	lonRad := lonDeg * DegToRad
	latRad := latDeg * DegToRad

	// identify base facet (0..5) and local planar coordinates [0, 1]
	facet, xLocal, yLocal := d.IdentifyBaseFacet(lonRad, latRad)

	// traverse hierarchy using configured nSide (subdivision factor)
	nSideFloat := float64(d.nSide)
	path := make([]uint8, targetRes)

	for r := uint8(0); r < targetRes; r++ {
		xIdx := int(xLocal * nSideFloat)
		yIdx := int(yLocal * nSideFloat)

		// clamp boundary precision artifacts
		if xIdx >= int(d.nSide) {
			xIdx = int(d.nSide) - 1
		} else if xIdx < 0 {
			xIdx = 0
		}
		if yIdx >= int(d.nSide) {
			yIdx = int(d.nSide) - 1
		} else if yIdx < 0 {
			yIdx = 0
		}

		subCell := uint8(yIdx*int(d.nSide) + xIdx)
		path[r] = subCell

		// scale local coordinates down into selected sub-cell
		xLocal = (xLocal * nSideFloat) - float64(xIdx)
		yLocal = (yLocal * nSideFloat) - float64(yIdx)
	}

	return facet, path, nil
}

// IdentifyBaseFacet determines which base facet (0..5: N, O, P, Q, R, S)
// a given longitude and authalic latitude (in radians) belongs to,
// returning the facet index (0..5) and normalised local coordinates (xLocal, yLocal in [0, 1])
// measured from the Upper-Left (North-West) corner of the facet.
//
// Direct Go port of rhealpixdggs-py (healpix_sphere + combine_triangles)
// with default north_square = 0 and south_square = 0.
//
// Topology:
//
//	+---+
//	| N |           <- North Polar Cap (N) attached above Facet O
//	+---+---+---+---+
//	| O | P | Q | R | <- Equatorial Belt (Facets O, P, Q, R)
//	+---+---+---+---+
//	| S |           <- South Polar Cap (S) attached below Facet O
//	+---+
func (d *RhpDGGS) IdentifyBaseFacet(lonRad, latRad float64) (uint8, float64, float64) {
	normLonRad := d.ellipsoid.NormaliseLongitude(lonRad)
	xRhp, yRhp, facet := d.ForwardProjectRadians(normLonRad, latRad)

	xLocal, yLocal := d.RadiansToFacetLocal(facet, xRhp, yRhp)
	return facet, xLocal, yLocal
}

func clampUnit(v float64) float64 {
	if v < 0.0 {
		return 0.0
	}
	if v > 1.0 {
		return 1.0
	}
	return v
}

// InverseProject converts planar rHEALPix projected coordinates (x, y in metres)
// back to geodetic longitude and latitude (in degrees) on the ellipsoid.
func (d *RhpDGGS) InverseProject(xMetres, yMetres float64) (lonDeg, latDeg float64, err error) {
	// scale down metres to authalic unit sphere radians (dimensionless unit sphere radians)
	xRhp := xMetres / d.ellipsoid.RA
	yRhp := yMetres / d.ellipsoid.RA

	// inverse project radians to geodetic lonRad and latRad
	lonRad, latRad, err := d.InverseProjectRadians(xRhp, yRhp)
	if err != nil {
		return 0, 0, err
	}

	return lonRad * RadToDeg, latRad * RadToDeg, nil
}

// FacetLocalToRadians transforms local unit square coordinates [0, 1] x [0, 1]
// back to global rHEALPix planar radians (xRhp, yRhp) for a given base facet.
func (d *RhpDGGS) FacetLocalToRadians(facet uint8, xLocal, yLocal float64) (xRhp, yRhp float64) {
	b := d.baseFacetBounds[facet]

	xRhp = b.ULX + (xLocal * b.Width)
	yRhp = b.ULY - (yLocal * b.Height)

	return xRhp, yRhp
}

// RadiansToFacetLocal transforms global planar radians (xRhp, yRhp) into
// local unit square coordinates [0, 1] x [0, 1] for a given base facet.
func (d *RhpDGGS) RadiansToFacetLocal(facet uint8, xRhp, yRhp float64) (xLocal, yLocal float64) {
	b := d.baseFacetBounds[facet]

	// xLocal: distance right from ULX divided by facet width
	xLocal = (xRhp - b.ULX) / b.Width

	// yLocal: distance down from ULY divided by facet height
	yLocal = (b.ULY - yRhp) / b.Height

	return clampUnit(xLocal), clampUnit(yLocal)
}

// LonLatToPlanar converts geodetic longitude and latitude (in degrees)
// specifically FOR the target base facet (0..5), returning normalised local
// planar coordinates (xLocal, yLocal in [0, 1]) measured from the Upper-Left (NW) corner.
//
// Longitudes outside [-180, 180] (unwrapped) are automatically normalised relative
// to the target facet's domain.
func (d *RhpDGGS) LonLatToPlanar(facetID uint8, lonDeg, latDeg float64) (xLocal, yLocal float64, err error) {
	if facetID > 5 {
		return 0, 0, fmt.Errorf("invalid facet ID %d: must be 0..5", facetID)
	}

	lonRad := lonDeg * DegToRad
	latRad := latDeg * DegToRad

	capLatRad := d.ellipsoid.Phi0

	// clamp geodetic latRad near cap boundary to prevent branch flapping
	switch facetID {
	case 0: // North Cap
		if latRad < capLatRad || math.Abs(latRad-capLatRad) < EPS_SEAM_RAD {
			latRad = capLatRad
		}
	case 5: // South Cap
		if latRad > -capLatRad || math.Abs(latRad-(-capLatRad)) < EPS_SEAM_RAD {
			latRad = -capLatRad
		}
	default: // Equatorial Facets (1..4)
		if math.Abs(latRad-capLatRad) < EPS_SEAM_RAD {
			latRad = capLatRad
		} else if math.Abs(latRad-(-capLatRad)) < EPS_SEAM_RAD {
			latRad = -capLatRad
		}
	}

	// ForwardProjectRadians handles normalisation and AuthLat conversion internally
	xRhp, yRhp, _ := d.ForwardProjectRadians(lonRad, latRad)
	xLocal, yLocal = d.RadiansToFacetLocal(facetID, xRhp, yRhp)

	return clampUnit(xLocal), clampUnit(yLocal), nil
}

// Triangle returns the HEALPix polar triangle number (0..3, or -1 for equatorial)
// and the CellRegion ('RegionNorthPolar', 'RegionSouthPolar', 'RegionEquatorial') for a point (x, y)
// given in global authalic planar radians.
//
// PARAMETERS:
//   - x, y: Coordinates in global authalic planar radians [-pi, pi] for X, [-pi/2, pi/2] for Y.
//   - inverse: If false, computes forward mapping assuming (x, y) is in HEALPix space.
//     If true, computes inverse mapping assuming (x, y) is in rHEALPix space (rearranged polar caps).
//
// top-level HEALPix vs rHEALPix POLAR CAP LAYOUTS:
//
//  1. Forward HEALPix Projection Space (inverse = false):
//     The 4 northern and 4 southern polar cap regions exist as triangular fins (0, 1, 2, 3):
//
// //          *       *       *       *
// //        * 0 *   * 1 *   * 2 *   * 3 *     <- North Polar Triangles
// //	   *-------*-------*-------*-------*
// //	   |   1   |   2   |   3   |   4   |   <- Equatorial Facets (1..4)
// //	   *-------*-------*-------*-------*
// //        * 0 *   * 1 *   * 2 *   * 3 *     <- South Polar Triangles
// //          *       *       *       *
//
//  2. Inverse rHEALPix Projection Space (inverse = true):
//     The triangular fins are rearranged into square polar cap facets (Facet 0 North, Facet 5 South).
//     The triangles rotate around the configured fixed square positions (northSquare / southSquare).
//
//     For example, if northSquare = 1 and southSquare = 3:
//
//     North Polar Square (Facet 0):
//     *-------*-------*
//     |       | * 3 * |       L2: -x + (-3pi/4 + (northSquare+1)pi/2)
//     |   0   | 0 * 2 |       L1:  x - (-3pi/4 + (northSquare-1)pi/2)
//     |       | * 1 * |
//     *-------*-------*
//
//     South Polar Square (Facet 5):
//     *-------*-------*
//     |       | * 3 * |       L2: -x + (-3pi/4 + (southSquare-1)pi/2)
//     |   5   | 2 * 0 |       L1:  x - (-3pi/4 + (southSquare+1)pi/2)
//     |       | * 1 * |
//     *-------*-------*
//
// OUTPUT:
//
//	Returns (triangle_number, region), where triangle_number is 0, 1, 2, or 3 for polar regions,
//	or -1 if the region is RegionEquatorial.
func (d *RhpDGGS) Triangle(x, y float64, inverse bool) (int8, CellRegion) {
	// determine regional belt by planar latitude boundaries (+/- Pi/4)
	var region CellRegion
	if y > QuarterPi {
		region = RegionNorthPolar
	} else if y < -QuarterPi {
		region = RegionSouthPolar
	} else {
		region = RegionEquatorial
	}

	if region == RegionEquatorial {
		return -1, RegionEquatorial
	}

	ns := float64(d.northSquare)
	ss := float64(d.southSquare)

	if !inverse {
		return TrianglePolarSector(x), region
	}

	// inverse rHEALPix domain mapping (triangle identification in polar cap squares)
	// EPS_DOMAIN: fuzz to avoid some rounding errors (direct port from Python)

	if region == RegionNorthPolar {
		L1 := x - (-ThreeQtrPi + (ns-1.0)*HalfPi)
		L2 := -x + (-ThreeQtrPi + (ns+1.0)*HalfPi)

		if y < L1-EPS_DOMAIN && y >= L2-EPS_DOMAIN {
			return int8((d.northSquare + 1) % 4), RegionNorthPolar
		} else if y >= L1-EPS_DOMAIN && y > L2+EPS_DOMAIN {
			return int8((d.northSquare + 2) % 4), RegionNorthPolar
		} else if y > L1+EPS_DOMAIN && y <= L2+EPS_DOMAIN {
			return int8((d.northSquare + 3) % 4), RegionNorthPolar
		}
		return int8(d.northSquare), RegionNorthPolar
	}

	// RegionSouthPolar
	L1 := x - (-ThreeQtrPi + (ss+1.0)*HalfPi)
	L2 := -x + (-ThreeQtrPi + (ss-1.0)*HalfPi)

	if y <= L1+EPS_DOMAIN && y > L2+EPS_DOMAIN {
		return int8((d.southSquare + 1) % 4), RegionSouthPolar
	} else if y < L1-EPS_DOMAIN && y <= L2+EPS_DOMAIN {
		return int8((d.southSquare + 2) % 4), RegionSouthPolar
	} else if y >= L1-EPS_DOMAIN && y < L2-EPS_DOMAIN {
		return int8((d.southSquare + 3) % 4), RegionSouthPolar
	}
	return int8(d.southSquare), RegionSouthPolar
}

// TrianglePolarSector identifies the cap triangle sector (0, 1, 2, or 3)
// based on planar X longitude in standard HEALPix space.
func TrianglePolarSector(x float64) int8 {
	switch {
	case x < -HalfPi:
		return 0
	case x < 0.0:
		return 1
	case x < HalfPi:
		return 2
	default:
		return 3
	}
}

// healpixSpherePolar handles pure Polar Cap math for both North and South caps.
func (d *RhpDGGS) healpixSpherePolar(lonRad, beta float64) (x, y float64) {
	// radius contraction factor (must use Abs so southern negative sines don't invert)
	sigma := math.Sqrt(3.0 * (1.0 - math.Abs(math.Sin(beta))))

	capNum := int(math.Floor(2.0*lonRad/Pi + 2.0))
	if capNum < 0 {
		capNum = 0
	} else if capNum > 3 {
		capNum = 3
	}

	lamc := -ThreeQtrPi + (HalfPi * float64(capNum))
	x = lamc + (lonRad-lamc)*sigma

	sign := 1.0
	if beta < 0 {
		sign = -1.0
	}
	y = sign * QuarterPi * (2.0 - sigma)

	return x, y
}

// healpixSphereEquatorial handles pure Equatorial Belt math.
func (d *RhpDGGS) healpixSphereEquatorial(lonRad, beta float64) (x, y float64) {
	return lonRad, ThreeQtrPi / 2.0 * math.Sin(beta) // 3*Pi/8 * sin(beta)
}

// HealpixSphere auto-routes coordinates to the correct domain based on beta vs Beta0 thresholds.
func (d *RhpDGGS) HealpixSphere(lonRad, beta float64) (x, y float64) {
	capLat := d.ellipsoid.Beta0

	if math.Abs(beta) > capLat {
		x, y = d.healpixSpherePolar(lonRad, beta)
		if beta > 0 {
			return x, y
		}
		return x, y
	}

	x, y = d.healpixSphereEquatorial(lonRad, beta)
	return x, y
}

// applyCapRotation performs the 2D affine rotation and translation
// for the polar cap triangles, matching rhealpixdggs-py reference behavior.
func (d *RhpDGGS) applyCapRotation(x, y float64, c int8, region CellRegion, inverse bool) (float64, float64) {
	northSquare := int8(d.northSquare % 4)
	southSquare := int8(d.southSquare % 4)

	// determine sign for Y-axis scaling depending on North vs South polar region
	signY := 1.0
	if region == RegionSouthPolar {
		signY = -1.0
	}

	// centre of the polar triangle apex in standard HEALPix space: tc
	// tc = (-3*pi/4 + c*pi/2, signY * pi/2)
	tcX := -ThreeQtrPi + float64(c)*HalfPi
	tcY := signY * HalfPi

	var pole int8
	var rotK int8

	if !inverse {
		// --- FORWARD TRANSFORMATION ---
		if region == RegionNorthPolar {
			pole = northSquare
			rotK = (c - pole) % 4
		} else {
			pole = southSquare
			rotK = (-(c - pole)) % 4
		}
		if rotK < 0 {
			rotK += 4
		}

		// u = (-3*pi/4 + pole*pi/2, signY * pi/2)
		uX := -ThreeQtrPi + float64(pole)*HalfPi
		uY := tcY

		// translate point relative to tc: (v - tc)
		dx := x - tcX
		dy := y - tcY

		// apply 2D matrix rotation R[rotK]
		rx, ry := rotate2D(dx, dy, rotK)

		// translate to destination apex u
		return rx + uX, ry + uY
	}

	// --- INVERSE TRANSFORMATION ---
	if region == RegionNorthPolar {
		pole = northSquare
		rotK = (-(c - pole)) % 4
	} else {
		pole = southSquare
		rotK = (c - pole) % 4
	}
	if rotK < 0 {
		rotK += 4
	}

	// target apex in rHEALPix space: u
	uX := -ThreeQtrPi + float64(pole)*HalfPi
	uY := tcY

	// translate point relative to u: (v - u)
	dx := x - uX
	dy := y - uY

	// apply inverse 2D matrix rotation R[rotK]
	rx, ry := rotate2D(dx, dy, rotK)

	// translate back to standard HEALPix apex tc
	return rx + tcX, ry + tcY
}

// rotate2D applies discrete 90-degree rotations (k = 0, 1, 2, 3)
func rotate2D(x, y float64, k int8) (float64, float64) {
	switch k % 4 {
	case 0: // 0 degrees
		return x, y
	case 1: // 90 degrees CCW: (-y, x)
		return -y, x
	case 2: // 180 degrees: (-x, -y)
		return -x, -y
	case 3: // 270 degrees CCW: (y, -x)
		return y, -x
	default:
		return x, y
	}
}

// CombineTriangles is the "Triangle Rearrangement" (rHEALPix Facet Assembly)
func (d *RhpDGGS) CombineTriangles(x, y float64, inverse bool) (float64, float64) {
	c, region := d.Triangle(x, y, inverse)
	if region == RegionEquatorial {
		return x, y
	}

	// apply 2D matrix rotation + vector translation (matching Python reference)
	return d.applyCapRotation(x, y, c, region, inverse)
}

// CombineTrianglesWithFacet performs triangle rearrangement and returns the base facet ID (0..5).
// If forcedRegion is RegionAuto (0), it auto-detects the region via Triangle().
func (d *RhpDGGS) CombineTrianglesWithFacet(x, y float64, inverse bool, normLonRad float64, forcedRegion CellRegion) (float64, float64, uint8) {
	region := forcedRegion
	var c int8

	if region == RegionAuto {
		c, region = d.Triangle(x, y, inverse)
	} else if region != RegionEquatorial {
		c = TrianglePolarSector(x)
	}

	switch region {
	case RegionNorthPolar:
		rx, ry := d.applyCapRotation(x, y, c, region, inverse)
		return rx, ry, 0 // Facet 0 (N)

	case RegionSouthPolar:
		rx, ry := d.applyCapRotation(x, y, c, region, inverse)
		return rx, ry, 5 // Facet 5 (S)

	default: // RegionEquatorial
		var facet uint8
		switch {
		case normLonRad < -HalfPi:
			facet = 1 // Facet O
		case normLonRad < 0.0:
			facet = 2 // Facet P
		case normLonRad < HalfPi:
			facet = 3 // Facet Q
		default:
			facet = 4 // Facet R
		}
		return x, y, facet
	}
}

// InverseProjectRadians converts global rHEALPix planar coordinates (in radians on unit sphere)
// back to longitude and latitude (in geodetic radians).
func (d *RhpDGGS) InverseProjectRadians(x, y float64) (float64, float64, error) {
	// domain guard
	if !d.InRhealpixImage(x, y) {
		return 0, 0, fmt.Errorf("point (%f, %f) is outside the valid rHEALPix domain", x, y)
	}

	// uncombine polar triangles (rHEALPix -> HEALPix)
	hX, hY := d.CombineTriangles(x, y, true)

	// inverse HEALPix projection (HEALPix -> Authalic lonRad, beta)
	lonRad, beta, err := d.inverseHealpixSphere(hX, hY)
	if err != nil {
		return 0, 0, err
	}

	// convert authalic latitude (beta) -> geodetic latitude (phi)
	latRad := d.ellipsoid.AuthLatInverse(beta)

	return d.ellipsoid.NormaliseLongitude(lonRad), latRad, nil
}

// inverseHealpixSphere converts standard HEALPix planar coordinates (x, y)
// into authalic spherical coordinates (lonRad, beta).
func (d *RhpDGGS) inverseHealpixSphere(x, y float64) (float64, float64, error) {
	absY := math.Abs(y)

	// --- Equatorial Region: |y| <= Pi/4 ---
	if absY <= QuarterPi {
		lon := x
		beta := math.Asin((8.0 * y) / (3.0 * Pi))
		return lon, beta, nil
	}

	// --- Polar Region (Non-Pole): Pi/4 < |y| < Pi/2 ---
	if absY < HalfPi {
		// cap number assignment clamped strictly to [0, 3]
		capNum := math.Floor((2.0*x)/Pi + 2.0)
		if capNum < 0 {
			capNum = 0
		} else if capNum > 3 {
			capNum = 3
		}

		xc := -ThreeQtrPi + (HalfPi * capNum)
		tau := 2.0 - (4.0 * absY / Pi)

		lon := xc + (x-xc)/tau

		// handle floating-point rounding errors in longitude along polar boundaries
		if lon < -Pi {
			lon = -Pi
		} else if lon > Pi {
			lon = Pi
		}

		sign := 1.0
		if y < 0 {
			sign = -1.0
		}
		beta := sign * math.Asin(1.0-(tau*tau)/3.0)

		return lon, beta, nil
	}

	// --- Poles: |y| >= Pi/2 ---
	if math.Abs(absY-HalfPi) < 1e-12 {
		sign := 1.0
		if y < 0 {
			sign = -1.0
		}
		return -Pi, sign * HalfPi, nil
	}

	return 0, 0, fmt.Errorf("y coordinate %f outside HEALPix domain [-pi/2, pi/2]", y)
}

// InHealpixImage returns true if (x, y) lies inside (or on the boundary of)
// the unarranged HEALPix projection image.
func InHealpixImage(x, y float64) bool {
	absY := math.Abs(y)

	// hard bounding box check
	if math.Abs(x) > Pi+EPS_IMAGE || absY > HalfPi+EPS_IMAGE {
		return false
	}

	// Equatorial belt check (|y| <= pi/4 + eps)
	if absY <= QuarterPi+EPS_IMAGE {
		return true
	}

	// Polar cap triangles check
	// |x - xc| + |y| <= pi/2 + eps
	capNum := math.Floor((2.0*x)/Pi + 2.0)
	if capNum < 0 {
		capNum = 0
	} else if capNum > 3 {
		capNum = 3
	}

	xc := -ThreeQtrPi + (HalfPi * capNum)
	return (math.Abs(x-xc) + absY) <= HalfPi+EPS_IMAGE
}

// InRhealpixImage returns true if (x, y) lies inside (or on the boundary of)
// the assembled rHEALPix projection image.
func (d *RhpDGGS) InRhealpixImage(x, y float64) bool {
	absY := math.Abs(y)

	// hard bounding box check: X in [-pi, pi], Y in [-3pi/4, 3pi/4]
	if math.Abs(x) > Pi+EPS_IMAGE || absY > ThreeQtrPi+EPS_IMAGE {
		return false
	}

	// main rectangular body: X in [-pi, pi] and Y in [-pi/4, pi/4]
	if absY <= QuarterPi+EPS_IMAGE {
		return true
	}

	ns := float64(d.northSquare % 4)
	ss := float64(d.southSquare % 4)

	// North Polar Cap Square (Facet 0):
	// Y in [pi/4, 3pi/4] and X in [-pi + ns*pi/2, -pi + (ns+1)*pi/2]
	if y > QuarterPi {
		minX := -Pi + (ns * HalfPi) - EPS_IMAGE
		maxX := -Pi + ((ns + 1.0) * HalfPi) + EPS_IMAGE
		return x >= minX && x <= maxX
	}

	// South Polar Cap Square (Facet 5):
	// Y in [-3pi/4, -pi/4] and X in [-pi + ss*pi/2, -pi + (ss+1)*pi/2]
	if y < -QuarterPi {
		minX := -Pi + (ss * HalfPi) - EPS_IMAGE
		maxX := -Pi + ((ss + 1.0) * HalfPi) + EPS_IMAGE
		return x >= minX && x <= maxX
	}

	return false
}
