package rhealpixorb

import (
	"fmt"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/clip"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// FacetSubGeometry holds a WGS84 geometry clipped to a specific rHEALPix base facet.
type FacetSubGeometry struct {
	FacetID  uint8
	Geometry orb.Geometry
}

// ClipGeometryToFacets splits an arbitrary WGS84 geometry across the 6 rHEALPix base facets.
// Uses Cell64ToPolygon with segmentCount = 16 to capture geodetic boundary curvature accurately.
func ClipGeometryToFacets(el *rhealpix.Ellipsoid, geom orb.Geometry) ([]FacetSubGeometry, error) {
	if geom == nil {
		return nil, fmt.Errorf("nil input geometry")
	}

	var results []FacetSubGeometry

	for facetID := uint8(0); facetID < 6; facetID++ {
		// construct root cell ID for base facet Level 0
		rootCell, err := rhealpix.PackCellID64(facetID, 0, nil)
		if err != nil {
			return nil, fmt.Errorf("failed to build root cell for facet %d: %w", facetID, err)
		}

		// densify facet envelope with 16 sub-segments per edge
		facetPoly, err := Cell64ToPolygon(el, rootCell, 16)
		if err != nil {
			return nil, fmt.Errorf("failed to build polygon for facet %d: %w", facetID, err)
		}

		// get the WGS84 bounding envelope for the facet
		facetBound := facetPoly.Bound()

		// quick check: if input geometry bound doesn't even touch the facet bound, skip early
		if !geom.Bound().Intersects(facetBound) {
			continue
		}

		// clip input geometry using orb/clip against the facet bound
		clipped := clip.Geometry(facetBound, geom)
		if clipped != nil {
			results = append(results, FacetSubGeometry{
				FacetID:  facetID,
				Geometry: clipped,
			})
		}

	}

	return results, nil
}
