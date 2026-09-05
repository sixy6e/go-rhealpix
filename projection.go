package rhealpix

import (
	"fmt"
	"math"
)

// NormaliseLongitudeRad normalises any longitude in radians into the canonical [-pi, pi) range.
func NormaliseLongitudeRad(lonRad float64) float64 {
	pi := math.Pi
	lon := math.Mod(lonRad+pi, 2.0*pi)
	if lon < 0 {
		lon += 2.0 * pi
	}
	return lon - pi
}

// AuthLat computes authalic latitude (beta in radians) from geodetic latitude (phi in radians).
func (el *Ellipsoid) AuthLat(phi float64) float64 {
	sinPhi := math.Sin(phi)

	// quick check for polar extremities
	if sinPhi >= 0.99999999 {
		return math.Pi / 2.0
	}
	if sinPhi <= -0.99999999 {
		return -math.Pi / 2.0
	}

	// handle spherical case (E == 0) safely
	if el.E == 0 {
		return phi
	}

	sin2 := sinPhi * sinPhi
	q := el.OneMinusE2 * (sinPhi/(1.0-el.E2*sin2) - (1.0/(2.0*el.E))*math.Log((1.0-el.E*sinPhi)/(1.0+el.E*sinPhi)))
	sinBeta := q / el.QBar

	if sinBeta > 1.0 {
		sinBeta = 1.0
	} else if sinBeta < -1.0 {
		sinBeta = -1.0
	}

	return math.Asin(sinBeta)
}

// ForwardProject converts geodetic longitude and latitude (in degrees)
// to planar rHEALPix projected coordinates (x, y in meters) on the ellipsoid.
func (el *Ellipsoid) ForwardProject(lonDeg, latDeg float64) (xMeters, yMeters float64) {
	lonRad := lonDeg * (math.Pi / 180.0)
	latRad := latDeg * (math.Pi / 180.0)

	// compute authalic latitude (beta)
	beta := el.AuthLat(latRad)

	// projected coordinates in radians on unit sphere
	xRhp, yRhp := ForwardProjectRadians(lonRad, beta)

	// scale by Authalic Radius (R_A) to get meters
	return xRhp * el.RA, yRhp * el.RA
}

// ForwardProjectRadians converts longitude and authalic latitude (in radians)
// into global rHEALPix planar coordinates (in radians on the unit authalic sphere).
func ForwardProjectRadians(lonRad, beta float64) (xRhp, yRhp float64) {
	pi := math.Pi
	halfPi := pi / 2.0
	quarterPi := pi / 4.0
	threeQuarterPi := 3.0 * pi / 4.0
	phi0 := math.Asin(2.0 / 3.0)

	lon := NormaliseLongitudeRad(lonRad)

	// exact Poles
	if beta >= halfPi-1e-15 {
		return -threeQuarterPi, halfPi
	}
	if beta <= -halfPi+1e-15 {
		return -threeQuarterPi, -halfPi
	}

	// Equatorial Region
	if math.Abs(beta) <= phi0 {
		xHp := lon
		yHp := (3.0 * pi / 8.0) * math.Sin(beta)
		return xHp, yHp
	}

	// Polar Caps
	sigma := math.Sqrt(3.0 * (1.0 - math.Abs(math.Sin(beta))))

	capNumber := int(math.Floor(2.0*lon/pi + 2.0))
	if capNumber >= 4 {
		capNumber = 3
	} else if capNumber < 0 {
		capNumber = 0
	}

	lamc := -threeQuarterPi + (halfPi * float64(capNumber))
	xHp := lamc + (lon-lamc)*sigma

	if beta > 0 {
		yHp := quarterPi * (2.0 - sigma)

		c := capNumber
		tcX := -threeQuarterPi + float64(c)*halfPi
		tcY := halfPi

		dx := xHp - tcX
		dy := yHp - tcY

		var rx, ry float64
		switch c {
		case 0:
			rx, ry = dx, dy
		case 1:
			rx, ry = -dy, dx
		case 2:
			rx, ry = -dx, -dy
		case 3:
			rx, ry = dy, -dx
		}

		return rx + (-threeQuarterPi), ry + halfPi
	}

	// South Polar Cap
	yHp := -quarterPi * (2.0 - sigma)

	c := capNumber
	tcX := -threeQuarterPi + float64(c)*halfPi
	tcY := -halfPi

	dx := xHp - tcX
	dy := yHp - tcY

	rotIdx := (4 - c) % 4
	var rx, ry float64
	switch rotIdx {
	case 0:
		rx, ry = dx, dy
	case 1:
		rx, ry = -dy, dx
	case 2:
		rx, ry = -dx, -dy
	case 3:
		rx, ry = dy, -dx
	}

	return rx + (-threeQuarterPi), ry + (-halfPi)
}

// ForwardTransform64 converts Lon/Lat (degrees) into a 64-bit CellID64 (up to Res 14).
func ForwardTransform64(el *Ellipsoid, lonDeg, latDeg float64, targetRes uint8) (CellID64, error) {
	facet, path, err := computeCellPath(el, lonDeg, latDeg, targetRes, MaxResolution64)
	if err != nil {
		return 0, err
	}
	return PackCellID64(facet, targetRes, path)
}

// ForwardTransform128 converts Lon/Lat (degrees) into a 128-bit CellID128 (up to Res 30).
func ForwardTransform128(el *Ellipsoid, lonDeg, latDeg float64, targetRes uint8) (CellID128, error) {
	facet, path, err := computeCellPath(el, lonDeg, latDeg, targetRes, MaxResolution128)
	if err != nil {
		return CellID128{}, err
	}
	return PackCellID128(facet, targetRes, path)
}

// computeCellPath converts Lon/Lat degrees to base facet and 3x3 sub-cell digit path.
func computeCellPath(el *Ellipsoid, lonDeg, latDeg float64, targetRes uint8, maxAllowedRes uint8) (uint8, []uint8, error) {
	if targetRes > maxAllowedRes {
		return 0, nil, ValidateCellParams(0, targetRes, maxAllowedRes, int(targetRes))
	}

	lonRad := lonDeg * (math.Pi / 180.0)
	latRad := latDeg * (math.Pi / 180.0)

	// convert to authalic latitude
	beta := el.AuthLat(latRad)

	// identify base facet (0..5) and top-down local planar coordinates (X, Y in [0, 1])
	facet, xLocal, yLocal := IdentifyBaseFacet(lonRad, beta)

	// traverse 3x3 hierarchy to collect sub-cell indices
	path := make([]uint8, targetRes)

	for r := uint8(0); r < targetRes; r++ {
		// divide local square into 3x3 grid
		xIdx := int(xLocal * 3.0)
		yIdx := int(yLocal * 3.0)

		// clamp edge boundary precision artifacts
		if xIdx > 2 {
			xIdx = 2
		} else if xIdx < 0 {
			xIdx = 0
		}
		if yIdx > 2 {
			yIdx = 2
		} else if yIdx < 0 {
			yIdx = 0
		}

		subCell := uint8(yIdx*3 + xIdx)
		path[r] = subCell

		// scale local coordinates down into the selected sub-cell
		xLocal = (xLocal * 3.0) - float64(xIdx)
		yLocal = (yLocal * 3.0) - float64(yIdx)
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
func IdentifyBaseFacet(lonRad, beta float64) (uint8, float64, float64) {
	pi := math.Pi
	halfPi := pi / 2.0
	quarterPi := pi / 4.0

	normLonRad := NormaliseLongitudeRad(lonRad)
	xRhp, yRhp := ForwardProjectRadians(normLonRad, beta)

	// Facet N
	if yRhp > quarterPi {
		xLocal := (xRhp - (-pi)) / halfPi
		yLocal := ((3.0 * pi / 4.0) - yRhp) / halfPi
		return 0, clampUnit(xLocal), clampUnit(yLocal)
	}

	// Facet S
	if yRhp < -quarterPi {
		xLocal := (xRhp - (-pi)) / halfPi
		yLocal := (-quarterPi - yRhp) / halfPi
		return 5, clampUnit(xLocal), clampUnit(yLocal)
	}

	// Equatorial Facets (O, P, Q, R)
	var facet uint8
	var lonMin float64

	switch {
	case normLonRad < -halfPi: // [-pi, -pi/2) -> Facet O
		facet = 1
		lonMin = -pi
	case normLonRad < 0.0: // [-pi/2, 0) -> Facet P
		facet = 2
		lonMin = -halfPi
	case normLonRad < halfPi: // [0, pi/2) -> Facet Q
		facet = 3
		lonMin = 0.0
	default: // [pi/2, pi] -> Facet R
		facet = 4
		lonMin = halfPi
	}

	xLocal := (normLonRad - lonMin) / halfPi
	yLocal := 0.5 * (1.0 - (yRhp / quarterPi))

	return facet, clampUnit(xLocal), clampUnit(yLocal)
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

// AuthLatInverse converts authalic latitude (beta in radians) back to geodetic latitude (phi in radians).
func (el *Ellipsoid) AuthLatInverse(beta float64) float64 {
	sinBeta := math.Sin(beta)

	if sinBeta >= 0.99999999 {
		return math.Pi / 2.0
	}
	if sinBeta <= -0.99999999 {
		return -math.Pi / 2.0
	}

	if el.E == 0 {
		return beta
	}

	// target authalic q value
	qTarget := sinBeta * el.QBar

	// Newton-Raphson solver for phi
	phi := beta // initial guess
	for i := 0; i < 5; i++ {
		sinPhi := math.Sin(phi)
		cosPhi := math.Cos(phi)
		sin2 := sinPhi * sinPhi

		if math.Abs(cosPhi) < 1e-12 {
			break
		}

		// q(phi)
		q := el.OneMinusE2 * (sinPhi/(1.0-el.E2*sin2) - (1.0/(2.0*el.E))*math.Log((1.0-el.E*sinPhi)/(1.0+el.E*sinPhi)))

		// dq/dphi = 2 * (1 - e^2) * cos(phi) / (1 - e^2 * sin^2(phi))^2
		denom := 1.0 - el.E2*sin2
		dqdphi := 2.0 * el.OneMinusE2 * cosPhi / (denom * denom)

		diff := q - qTarget
		if math.Abs(diff) < 1e-14 {
			break
		}

		phi -= diff / dqdphi
	}

	return phi
}

// InverseProject converts planar rHEALPix coordinates (x, y in meters)
// back to geodetic longitude and latitude (in degrees).
func (el *Ellipsoid) InverseProject(xMeters, yMeters float64) (lonDeg, latDeg float64) {
	// scale down meters to authalic unit sphere radians
	xRhp := xMeters / el.RA
	yRhp := yMeters / el.RA

	lonRad, beta := InverseProjectRadians(xRhp, yRhp)
	latRad := el.AuthLatInverse(beta)

	return lonRad * (180.0 / math.Pi), latRad * (180.0 / math.Pi)
}

// InverseProjectRadians converts global rHEALPix planar coordinates (in radians on unit sphere)
// back to longitude and authalic latitude (in radians).
func InverseProjectRadians(xRhp, yRhp float64) (lonRad, beta float64) {
	pi := math.Pi
	halfPi := pi / 2.0
	quarterPi := pi / 4.0
	threeQuarterPi := 3.0 * pi / 4.0
	eps := 1e-15

	var xHp, yHp float64

	if yRhp > quarterPi {
		// North Polar Region (north_square = 0, at least until configuration is supported)
		// calculate diagonal bounding lines L1 and L2
		L1 := xRhp - (-threeQuarterPi - halfPi)
		L2 := -xRhp + (-threeQuarterPi + halfPi)

		var c int
		if yRhp < L1-eps && yRhp >= L2-eps {
			c = 1
		} else if yRhp >= L1-eps && yRhp > L2+eps {
			c = 2
		} else if yRhp > L1+eps && yRhp <= L2+eps {
			c = 3
		} else {
			c = 0
		}

		uX := -threeQuarterPi
		uY := halfPi

		rx := xRhp - uX
		ry := yRhp - uY

		// apply inverse rotation ROTATE[-(c - north_square)]
		// for north_square = 0, this is ROTATE[-c]
		rotIdx := (4 - (c % 4)) % 4
		var dx, dy float64
		switch rotIdx {
		case 0: // 0 deg
			dx, dy = rx, ry
		case 1: // +90 deg CCW
			dx, dy = -ry, rx
		case 2: // 180 deg
			dx, dy = -rx, -ry
		case 3: // +270 deg CCW
			dx, dy = ry, -rx
		}

		tcX := -threeQuarterPi + float64(c)*halfPi
		tcY := halfPi

		xHp = dx + tcX
		yHp = dy + tcY

	} else if yRhp < -quarterPi {
		// South Polar Region (south_square = 0, at least until configuration is supported)
		L1 := xRhp - (-threeQuarterPi + halfPi)
		L2 := -xRhp + (-threeQuarterPi - halfPi)

		var c int
		if yRhp <= L1+eps && yRhp > L2+eps {
			c = 1
		} else if yRhp < L1-eps && yRhp <= L2+eps {
			c = 2
		} else if yRhp >= L1-eps && yRhp < L2-eps {
			c = 3
		} else {
			c = 0
		}

		uX := -threeQuarterPi
		uY := -halfPi

		rx := xRhp - uX
		ry := yRhp - uY

		// apply inverse rotation ROTATE[c - south_square]
		// for south_square = 0, this is ROTATE[c]
		rotIdx := (c%4 + 4) % 4
		var dx, dy float64
		switch rotIdx {
		case 0:
			dx, dy = rx, ry
		case 1:
			dx, dy = -ry, rx
		case 2:
			dx, dy = -rx, -ry
		case 3:
			dx, dy = ry, -rx
		}

		tcX := -threeQuarterPi + float64(c)*halfPi
		tcY := -halfPi

		xHp = dx + tcX
		yHp = dy + tcY

	} else {
		// Equatorial Region
		xHp = xRhp
		yHp = yRhp
	}

	// invert HEALPix spherical projection (healpix_sphere_inverse)
	if math.Abs(yHp) <= quarterPi {
		// Equatorial Belt
		lonRad = xHp
		sinBeta := (8.0 / (3.0 * pi)) * yHp
		if sinBeta > 1.0 {
			sinBeta = 1.0
		} else if sinBeta < -1.0 {
			sinBeta = -1.0
		}
		beta = math.Asin(sinBeta)
	} else {
		// Polar Cap
		sigma := 2.0 - (4.0 * math.Abs(yHp) / pi)

		sinBeta := 1.0 - (sigma * sigma / 3.0)
		if yHp < 0 {
			sinBeta = -sinBeta
		}
		if sinBeta > 1.0 {
			sinBeta = 1.0
		} else if sinBeta < -1.0 {
			sinBeta = -1.0
		}
		beta = math.Asin(sinBeta)

		capNumber := int(math.Floor(2.0*xHp/pi + 2.0))
		if capNumber >= 4 {
			capNumber = 3
		} else if capNumber < 0 {
			capNumber = 0
		}

		lamc := -threeQuarterPi + (halfPi * float64(capNumber))

		if sigma > 1e-15 {
			lonRad = lamc + (xHp-lamc)/sigma
		} else {
			lonRad = lamc
		}
	}

	return NormaliseLongitudeRad(lonRad), beta
}

// FacetLocalToRadians converts a base facet index (0..5) and local coordinates [0, 1] x [0, 1]
// back to global rHEALPix planar radians.
func FacetLocalToRadians(facet uint8, xLocal, yLocal float64) (xRhp, yRhp float64) {
	pi := math.Pi
	halfPi := pi / 2.0
	quarterPi := pi / 4.0
	threeQuarterPi := 3.0 * pi / 4.0

	switch facet {
	case 0: // Facet N
		xRhp = -pi + (xLocal * halfPi)
		yRhp = threeQuarterPi - (yLocal * halfPi)
	case 1: // Facet O
		xRhp = -pi + (xLocal * halfPi)
		yRhp = quarterPi - (yLocal * halfPi)
	case 2: // Facet P
		xRhp = -halfPi + (xLocal * halfPi)
		yRhp = quarterPi - (yLocal * halfPi)
	case 3: // Facet Q
		xRhp = 0.0 + (xLocal * halfPi)
		yRhp = quarterPi - (yLocal * halfPi)
	case 4: // Facet R
		xRhp = halfPi + (xLocal * halfPi)
		yRhp = quarterPi - (yLocal * halfPi)
	case 5: // Facet S
		xRhp = -pi + (xLocal * halfPi)
		yRhp = -quarterPi - (yLocal * halfPi)
	}

	return xRhp, yRhp
}

// RadiansToFacetLocal converts global planar radians (xRhp, yRhp)
// directly into local [0, 1] x [0, 1] coordinates for a specific base facet.
func RadiansToFacetLocal(facet uint8, xRhp, yRhp float64) (xLocal, yLocal float64) {
	pi := math.Pi
	halfPi := pi / 2.0
	quarterPi := pi / 4.0
	threeQuarterPi := 3.0 * pi / 4.0

	switch facet {
	case 0: // Facet N
		xLocal = (xRhp - (-pi)) / halfPi
		yLocal = (threeQuarterPi - yRhp) / halfPi
	case 1: // Facet O
		xLocal = (xRhp - (-pi)) / halfPi
		yLocal = (quarterPi - yRhp) / halfPi
	case 2: // Facet P
		xLocal = (xRhp - (-halfPi)) / halfPi
		yLocal = (quarterPi - yRhp) / halfPi
	case 3: // Facet Q
		xLocal = (xRhp - 0.0) / halfPi
		yLocal = (quarterPi - yRhp) / halfPi
	case 4: // Facet R
		xLocal = (xRhp - halfPi) / halfPi
		yLocal = (quarterPi - yRhp) / halfPi
	case 5: // Facet S
		xLocal = (xRhp - (-pi)) / halfPi
		yLocal = (-quarterPi - yRhp) / halfPi
	}

	return clampUnit(xLocal), clampUnit(yLocal)
}

// LonLatToPlanar converts geodetic longitude and latitude (in degrees)
// specifically FOR the target base facet (0..5), returning normalised local
// planar coordinates (xLocal, yLocal in [0, 1]) measured from the Upper-Left (NW) corner.
//
// Longitudes outside [-180, 180] (unwrapped) are automatically normalised relative
// to the target facet's domain.
func LonLatToPlanar(el *Ellipsoid, facetID uint8, lonDeg, latDeg float64) (xLocal, yLocal float64, err error) {
	if facetID > 5 {
		return 0, 0, fmt.Errorf("invalid facet ID %d: must be 0..5", facetID)
	}

	lonRad := lonDeg * (math.Pi / 180.0)
	latRad := latDeg * (math.Pi / 180.0)

	normLonRad := NormaliseLongitudeRad(lonRad)
	beta := el.AuthLat(latRad)

	pi := math.Pi
	halfPi := pi / 2.0

	// Equatorial Facets (O=1, P=2, Q=3, R=4)
	if facetID >= 1 && facetID <= 4 {
		var lonMin float64
		switch facetID {
		case 1: // Facet O: [-pi, -pi/2)
			lonMin = -pi
		case 2: // Facet P: [-pi/2, 0)
			lonMin = -halfPi
		case 3: // Facet Q: [0, pi/2)
			lonMin = 0.0
		case 4: // Facet R: [pi/2, pi)
			lonMin = halfPi
		}

		_, yRhp := ForwardProjectRadians(normLonRad, beta)

		xLocal = (normLonRad - lonMin) / halfPi
		yLocal = 0.5 * (1.0 - (yRhp / (pi / 4.0)))

		return clampUnit(xLocal), clampUnit(yLocal), nil
	}

	// Polar Cap Facets (N=0, S=5)
	detectedFacet, xLoc, yLoc := IdentifyBaseFacet(normLonRad, beta)
	if detectedFacet == facetID {
		return xLoc, yLoc, nil
	}

	xRhp, yRhp := ForwardProjectRadians(normLonRad, beta)
	xLoc, yLoc = RadiansToFacetLocal(facetID, xRhp, yRhp)

	return clampUnit(xLoc), clampUnit(yLoc), nil
}
