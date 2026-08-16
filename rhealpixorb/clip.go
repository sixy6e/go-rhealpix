package rhealpixorb

import (
	"fmt"
	"math"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/clip"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// FacetSubGeometry holds a WGS84 geometry clipped to a specific rHEALPix base facet.
type FacetSubGeometry struct {
	FacetID  uint8
	Geometry orb.Geometry
}

// facetIntersectsBound performs a geodetic bounding box check against
// the 6 rHEALPix base facet domains.
func facetIntersectsBound(facetID uint8, geomBound orb.Bound) bool {
	// WGS84 authalic transition latitude: math.Asin(2/3) * (180 / Pi) ≈ 41.8103149°
	capLat := math.Asin(2.0/3.0) * (180.0 / math.Pi)

	minY := geomBound.Min.Y()
	maxY := geomBound.Max.Y()
	minX := geomBound.Min.X()
	maxX := geomBound.Max.X()

	switch facetID {
	case 0: // North Polar Cap ('N' / 0)
		return maxY >= capLat

	case 1: // Equatorial Facet O / 1 [-180°, -90°]
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX <= -90.0 && maxX >= -180.0

	case 2: // Equatorial Facet P / 2 [-90°, 0°]
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX <= 0.0 && maxX >= -90.0

	case 3: // Equatorial Facet Q / 3 [0°, 90°]
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX <= 90.0 && maxX >= 0.0

	case 4: // Equatorial Facet R / 4 [90°, 180°]
		if maxY < -capLat || minY > capLat {
			return false
		}
		return minX <= 180.0 && maxX >= 90.0

	case 5: // South Polar Cap ('S' / 5)
		return minY <= -capLat

	default:
		return false
	}
}

// getFacetClipBound returns an approximate WGS84 bounding box for coarse filtering.
// Facet IDs match rHEALPix specification (0=N, 1..4=O,P,Q,R, 5=S).
func getFacetClipBound(facetID uint8) orb.Bound {
	// WGS84 authalic transition latitude: math.Asin(2/3) * (180 / Pi) ≈ 41.8103149°
	capLat := math.Asin(2.0/3.0) * (180.0 / math.Pi)

	switch facetID {
	case 0: // North Cap ('N')
		return orb.Bound{Min: orb.Point{-180.0, capLat}, Max: orb.Point{180.0, 90.0}}
	case 1: // Facet O
		return orb.Bound{Min: orb.Point{-180.0, -capLat}, Max: orb.Point{-90.0, capLat}}
	case 2: // Facet P
		return orb.Bound{Min: orb.Point{-90.0, -capLat}, Max: orb.Point{0.0, capLat}}
	case 3: // Facet Q
		return orb.Bound{Min: orb.Point{0.0, -capLat}, Max: orb.Point{90.0, capLat}}
	case 4: // Facet R
		return orb.Bound{Min: orb.Point{90.0, -capLat}, Max: orb.Point{180.0, capLat}}
	case 5: // South Cap ('S')
		return orb.Bound{Min: orb.Point{-180.0, -90.0}, Max: orb.Point{180.0, -capLat}}
	default:
		return orb.Bound{}
	}
}

// ClipGeometryToFacets splits an arbitrary WGS84 geometry across the 6 rHEALPix base facets.
func ClipGeometryToFacets(el *rhealpix.Ellipsoid, geom orb.Geometry) ([]FacetSubGeometry, error) {
	if geom == nil {
		return nil, fmt.Errorf("nil input geometry")
	}

	geomBound := geom.Bound()
	var results []FacetSubGeometry

	for facetID := uint8(0); facetID < 6; facetID++ {
		// facet domain intersection check
		if !facetIntersectsBound(facetID, geomBound) {
			continue
		}

		// rectangular bounding box for polygon clipping
		clipBound := getFacetClipBound(facetID)

		// clip input geometry against exact facet box
		clipped := clip.Geometry(clipBound, geom)
		if clipped != nil {
			results = append(results, FacetSubGeometry{
				FacetID:  facetID,
				Geometry: clipped,
			})
		}
	}

	return results, nil
}
