package rhealpix

import (
	"math"
)

// AreaErrorEntry holds the theoretical area and floating-point error bounds for a given resolution.
type AreaErrorEntry struct {
	CellAreaM2   float64 // Theoretical equal area per cell in square metres
	AbsTolerance float64 // Absolute tolerance (m²) for float comparisons
	RelTolerance float64 // Relative tolerance (dimensionless)
}

// CellWidth calculates the planar cell width in metres at a given resolution level.
// Planar width scales down by nSide^(-res) relative to the root facet width.
func (d *RhpDGGS) CellWidth(res uint8) float64 {
	// root facet planar width in metres = R_A * (pi / 2)
	baseWidth := d.ellipsoid.RA * HalfPi
	return baseWidth * math.Pow(float64(d.nSide), -float64(res))
}

// CellArea calculates the cell area in metres squared at a given resolution level.
// If plane is true, it returns the local planar cell area.
// If plane is false, it returns the true authalic ellipsoidal surface area.
func (d *RhpDGGS) CellArea(res uint8, plane bool) float64 {
	w := d.CellWidth(res)
	planarArea := w * w
	if plane {
		return planarArea
	}
	// scale planar area to true authalic sphere area: (8 / (3 * pi)) * planarArea
	return (8.0 / (3.0 * Pi)) * planarArea
}

// AreaErrorBudget returns an analytical error budget for cell area equality testing
// from resolution 0 up to maxRes (defaults to d.maxResolution if maxRes is 0).
//
// Relative tolerance is bounded by 10 * machine epsilon (~2.22e-15).
func (d *RhpDGGS) AreaErrorBudget(maxRes uint8) map[uint8]AreaErrorEntry {
	if maxRes == 0 {
		maxRes = d.maxResolution
	}

	// machine epsilon for IEEE 754 float64 is 2^(-52) ~ 2.220446049250313e-16
	const float64Epsilon = 1.0 / (1 << 52)
	relTol := 10.0 * float64Epsilon

	budget := make(map[uint8]AreaErrorEntry, maxRes+1)
	for r := uint8(0); r <= maxRes; r++ {
		area := d.CellArea(r, false)
		budget[r] = AreaErrorEntry{
			CellAreaM2:   area,
			AbsTolerance: area * relTol,
			RelTolerance: relTol,
		}
	}
	return budget
}

// CellPlanarWidthMetres returns the uniform planar cell edge length in metres at level res.
func (d *RhpDGGS) CellPlanarWidthMetres(res uint8) float64 {
	baseWidthMetres := d.ellipsoid.RA * HalfPi
	return baseWidthMetres * math.Pow(float64(d.nSide), -float64(res))
}

// CellPlanarWidthRadians returns the uniform planar cell edge length in radians at level res.
func (d *RhpDGGS) CellPlanarWidthRadians(res uint8) float64 {
	return HalfPi * math.Pow(float64(d.nSide), -float64(res))
}

// CellPlanarFractionalExtent returns normalised unit space bounds [0.0, 1.0] within the base facet.
// It accumulates integer grid coordinates to prevent float drift at deep resolutions.
func CellPlanarFractionalExtent[T CellID](d *RhpDGGS, id T) (uMinX, uMaxX, uMinY, uMaxY float64) {
	var anyID any = id
	var res uint8
	switch v := anyID.(type) {
	case CellID64:
		res = v.Resolution()
	case CellID128:
		res = v.Resolution()
	}

	if res == 0 {
		// Level 0 base facet bounds
		return 0.0, 1.0, 0.0, 1.0
	}

	// rawID := uint64(id)
	var colIdx, rowIdx uint64 = 0, 0
	nside := uint64(d.nSide)

	// accumulate discrete integer grid coordinates across resolution levels
	for i := uint8(0); i < res; i++ {
		// shift := 52 - (i * 4)
		// subCell := uint8((rawID >> shift) & SubCellMask)
		subCell := extractSubCellDigit(id, i)

		col := uint64(subCell % d.nSide)
		row := uint64(subCell / d.nSide)

		colIdx = colIdx*nside + col
		rowIdx = rowIdx*nside + row
	}

	scale := math.Pow(float64(d.nSide), float64(res))

	// normalised fractional offsets within the base facet [0.0, 1.0]
	uMinX = float64(colIdx) / scale
	uMaxX = float64(colIdx+1) / scale
	uMinY = float64(rowIdx) / scale
	uMaxY = float64(rowIdx+1) / scale

	return uMinX, uMaxX, uMinY, uMaxY
}

// CellPlanarExtent64 calculates the global planar bounding box (xMin, xMax, yMin, yMax)
// in authalic radians for a given 64-bit or 128-bit CellID.
func CellPlanarExtent[T CellID](d *RhpDGGS, id T) (xMin, xMax, yMin, yMax float64) {
	var facet uint8
	switch v := any(id).(type) {
	case CellID64:
		facet = v.Facet()
	case CellID128:
		facet = v.Facet()
	}

	uMinX, uMaxX, uMinY, uMaxY := CellPlanarFractionalExtent(d, id)
	b := d.baseFacetBounds[facet]

	// map fractional unit space directly to global planar radians:

	// X expands Right from ULX
	xMin = b.ULX + (uMinX * b.Width)
	xMax = b.ULX + (uMaxX * b.Width)

	// Y expands Down from ULY
	yMax = b.ULY - (uMinY * b.Height) // Northern / Top edge
	yMin = b.ULY - (uMaxY * b.Height) // Southern / Bottom edge

	return xMin, xMax, yMin, yMax
}

func CellPlanarVertices[T CellID](d *RhpDGGS, id T) [4][2]float64 {
	xMin, xMax, yMin, yMax := CellPlanarExtent(d, id)

	// local planar corners: UL, UR, LR, LL
	vertices := [4][2]float64{
		{xMin, yMax}, // UL
		{xMax, yMax}, // UR
		{xMax, yMin}, // LR
		{xMin, yMin}, // LL
	}

	return vertices
}

// Nucleus returns the central interior point of a cell in global planar radians [x, y].
// Matches Python's c.nucleus(plane=True).
func Nucleus[T CellID](d *RhpDGGS, id T) [2]float64 {
	xMin, xMax, yMin, yMax := CellPlanarExtent(d, id)
	return [2]float64{(xMin + xMax) / 2.0, (yMin + yMax) / 2.0}
}

// EllipsoidalVertices returns the geographic coordinates ([lon, lat] in degrees)
// of a cell in clockwise order starting from the NW vertex on the configured ellipsoid.
// If trimDart is true, the collapsed non-vertex point on Dart cells is omitted (returning 3 vertices).
func EllipsoidalVertices[T CellID](d *RhpDGGS, id T, trimDart bool) ([][2]float64, error) {
	var facet uint8
	switch v := any(id).(type) {
	case CellID64:
		facet = v.Facet()
	case CellID128:
		facet = v.Facet()
	}

	planarVerts := CellPlanarVertices(d, id) // [UL, UR, LR, LL]
	shape := EllipsoidalCellShape(d, id)

	// get index of true NW vertex in planar space directly
	nwIdx, err := NWIndex(d, id)
	if err != nil {
		return nil, err
	}

	// unproject planar authalic radians back to geographic degrees
	// starting at true NW
	geoVerts := make([][2]float64, 4)
	for i := 0; i < 4; i++ {
		pv := planarVerts[(nwIdx+i)%4]
		lonRad, latRad, err := d.InverseProjectRadians(pv[0], pv[1])
		if err != nil {
			return nil, err
		}

		// apply ellipsoid authalic latitude inversion
		if d.ellipsoid != nil {
			latRad = d.ellipsoid.AuthLatInverse(latRad)
		}

		geoVerts[i] = [2]float64{lonRad * RadToDeg, latRad * RadToDeg}
	}

	// trim collapsed point for Dart cells (returning 3 vertices)
	if trimDart && shape == ShapeDart {
		trimmed := make([][2]float64, 3)
		if FacetRegion(facet) == RegionNorthPolar {
			trimmed[0], trimmed[1], trimmed[2] = geoVerts[0], geoVerts[1], geoVerts[3]
		} else {
			trimmed[0], trimmed[1], trimmed[2] = geoVerts[0], geoVerts[2], geoVerts[3]
		}
		return trimmed, nil
	}

	return geoVerts, nil
}

// GeographicVertices is an alias for EllipsoidalVertices, providing standard GIS-friendly terminology.
func GeographicVertices[T CellID](d *RhpDGGS, id T) ([][2]float64, error) {
	return EllipsoidalVertices(d, id, true)
}

// GeographicVertices64 returns GIS-valid geographic coordinates for a 64-bit cell ID.
func (d *RhpDGGS) GeographicVertices64(id CellID64) ([][2]float64, error) {
	return GeographicVertices(d, id)
}

// GeographicVertices128 returns GIS-valid geographic coordinates for a 128-bit cell ID.
func (d *RhpDGGS) GeographicVertices128(id CellID128) ([][2]float64, error) {
	return GeographicVertices(d, id)
}

// EllipsoidalVertices64 returns geographic coordinates with explicit control over dart trimming.
func (d *RhpDGGS) EllipsoidalVertices64(id CellID64, trimDart bool) ([][2]float64, error) {
	return EllipsoidalVertices(d, id, trimDart)
}

// EllipsoidalVertices128 returns geographic coordinates with explicit control over dart trimming.
func (d *RhpDGGS) EllipsoidalVertices128(id CellID128, trimDart bool) ([][2]float64, error) {
	return EllipsoidalVertices(d, id, trimDart)
}

// EllipsoidalCellShape64 determines the shape of a 64-bit cell ID.
func (d *RhpDGGS) EllipsoidalCellShape64(id CellID64) CellShape {
	return EllipsoidalCellShape(d, id)
}

// EllipsoidalCellShape128 determines the shape of a 128-bit cell ID.
func (d *RhpDGGS) EllipsoidalCellShape128(id CellID128) CellShape {
	return EllipsoidalCellShape(d, id)
}

// CellRegion returns whether the cell is Equatorial, North Polar, or South Polar.
func (d *RhpDGGS) CellRegion(facet uint8) CellRegion {
	return FacetRegion(facet)
}

func (d *RhpDGGS) CellRegion64(id CellID64) CellRegion {
	return FacetRegion(id.Facet())
}

func (d *RhpDGGS) CellRegion128(id CellID128) CellRegion {
	return FacetRegion(id.Facet())
}

// NWVertex returns the true geographical Northwest vertex of a cell.
// If plane is true, it returns the point in global authalic planar radians [x, y].
// If plane is false, it returns the point in geographic degrees [lon, lat] on the configured ellipsoid.
func NWVertex[T CellID](d *RhpDGGS, id T, plane bool) ([2]float64, error) {
	planarVerts := CellPlanarVertices(d, id) // [UL, UR, LR, LL]
	shape := EllipsoidalCellShape(d, id)

	var planarNW [2]float64

	switch shape {
	case ShapeQuad, ShapeCap:
		// Northwest vertex is the upper-left vertex (UL)
		planarNW = planarVerts[0]

	case ShapeSkewQuad:
		// determine rotation offset based on the Collignon triangle of the cell's nucleus
		// xNuc, yNuc := Nucleus(d, id)
		nuc := Nucleus(d, id)
		triangle, region := d.Triangle(nuc[0], nuc[1], true)

		if region == RegionNorthPolar {
			ns := int(d.northSquare)
			i := (int(triangle) - ns) % 4
			if i < 0 {
				i += 4
			}
			// Python v[-i] -> Go v[(4 - i) % 4]
			idx := (4 - i) % 4
			planarNW = planarVerts[idx]
		} else {
			ss := int(d.southSquare)
			i := (int(triangle) - ss) % 4
			if i < 0 {
				i += 4
			}
			// Python v[i] -> Go v[i]
			planarNW = planarVerts[i]
		}

	case ShapeDart:
		// find the poleward vertex in planar space.
		// the pole is the center of the base polar facet (Facet 0 or 5).
		var facet uint8
		switch v := any(id).(type) {
		case CellID64:
			facet = v.Facet()
		case CellID128:
			facet = v.Facet()
		}

		b := d.baseFacetBounds[facet]
		centre := b.Centroid()
		// cx := b.ULX + (b.Width / 2.0)
		// cy := b.ULY - (b.Height / 2.0)
		cx, cy := centre[0], centre[1]

		// find vertex with smallest Chebyshev distance max(|x - cx|, |y - cy|) to pole
		minIdx := 0
		minDist := math.MaxFloat64
		for j := 0; j < 4; j++ {
			dx := math.Abs(planarVerts[j][0] - cx)
			dy := math.Abs(planarVerts[j][1] - cy)
			dist := math.Max(dx, dy)
			if dist < minDist {
				minDist = dist
				minIdx = j
			}
		}

		region := FacetRegion(facet)
		if region == RegionNorthPolar {
			// Northwest vertex is the poleward vertex
			planarNW = planarVerts[minIdx]
		} else {
			// Northwest vertex is one step clockwise from poleward vertex
			planarNW = planarVerts[(minIdx+1)%4]
		}
	}

	if plane {
		return planarNW, nil
	}

	// unproject planar radians to geographic degrees (lon, lat)
	lonRad, beta, err := d.InverseProjectRadians(planarNW[0], planarNW[1])
	if err != nil {
		return [2]float64{}, err
	}
	latRad := d.ellipsoid.AuthLatInverse(beta)

	return [2]float64{lonRad * RadToDeg, latRad * RadToDeg}, nil
}

// NWIndex returns the slice index (0, 1, 2, or 3) within planarVerts [UL, UR, LR, LL]
// that corresponds to the true NW vertex of the cell.
func NWIndex[T CellID](d *RhpDGGS, id T) (int, error) {
	shape := EllipsoidalCellShape(d, id)

	switch shape {
	case ShapeQuad, ShapeCap:
		return 0, nil // Upper-left (UL)

	case ShapeSkewQuad:
		nuc := Nucleus(d, id)
		triangle, region := d.Triangle(nuc[0], nuc[1], true)

		if region == RegionNorthPolar {
			ns := int(d.northSquare)
			i := (int(triangle) - ns) % 4
			if i < 0 {
				i += 4
			}
			return (4 - i) % 4, nil // Python v[-i] -> Go v[(4 - i) % 4]
		} else {
			ss := int(d.southSquare)
			i := (int(triangle) - ss) % 4
			if i < 0 {
				i += 4
			}
			return i, nil // Python v[i] -> Go v[i]
		}

	case ShapeDart:
		planarVerts := CellPlanarVertices(d, id)
		var facet uint8
		switch v := any(id).(type) {
		case CellID64:
			facet = v.Facet()
		case CellID128:
			facet = v.Facet()
		}

		b := d.baseFacetBounds[facet]
		centre := b.Centroid()
		cx, cy := centre[0], centre[1]

		minIdx := 0
		minDist := math.MaxFloat64
		for j := 0; j < 4; j++ {
			dx := math.Abs(planarVerts[j][0] - cx)
			dy := math.Abs(planarVerts[j][1] - cy)
			dist := math.Max(dx, dy)
			if dist < minDist {
				minDist = dist
				minIdx = j
			}
		}

		if FacetRegion(facet) == RegionNorthPolar {
			return minIdx, nil
		}
		return (minIdx + 1) % 4, nil
	}

	return 0, nil
}
