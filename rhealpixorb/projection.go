package rhealpixorb

import (
	"math"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// projectGeometryToPlanar maps a WGS84 geometry to normalised [0, 1] x [0, 1] facet local planar space.
func projectGeometryToPlanar(d *rhealpix.RhpDGGS, facetID uint8, geom orb.Geometry) orb.Geometry {
	if geom == nil {
		return nil
	}

	switch g := geom.(type) {
	case orb.Polygon:
		proj := projectPolygonToPlanar(d, facetID, g)
		return ensurePolygonCCW(proj)
	case orb.MultiPolygon:
		var mp orb.MultiPolygon
		for _, poly := range g {
			proj := projectPolygonToPlanar(d, facetID, poly)
			ccw := ensurePolygonCCW(proj)
			if len(ccw) > 0 {
				mp = append(mp, ccw)
			}
		}
		if len(mp) == 1 {
			return mp[0]
		}
		return mp
	default:
		return nil
	}
}

// projectPolygonToPlanar maps a WGS84 polygon to normalised [0, 1] x [0, 1] facet local planar space.
func projectPolygonToPlanar(d *rhealpix.RhpDGGS, facetID uint8, poly orb.Polygon) orb.Polygon {
	// densify polar cap facets (0=N, 5=S) to capture non-linear parallel curvature and 0° seam apexes
	if facetID == 0 || facetID == 5 {
		densifiedPoly := make(orb.Polygon, len(poly))
		for i, ring := range poly {
			densifiedPoly[i] = DensifyWGS84Ring(ring, 16)
		}
		poly = densifiedPoly
	}

	projected := make(orb.Polygon, len(poly))

	for i, ring := range poly {
		projectedRing := make(orb.Ring, len(ring))
		for j, pt := range ring {
			lon, lat := pt[0], pt[1]

			xLoc, yLoc, err := d.LonLatToPlanar(facetID, lon, lat)
			if err != nil {
				projectedRing[j] = orb.Point{lon, lat}
				continue
			}

			projectedRing[j] = orb.Point{xLoc, yLoc}
		}
		projected[i] = projectedRing
	}

	return projected
}

// ForwardProjectPolygon projects geographic points directly through RhpDGGS.
// TODO: cleanup; this is replaced with ForwardProjectPolygonForcedRegion
// to enforce topological correctness for geometry structures.
func ForwardProjectPolygon(d *rhealpix.RhpDGGS, poly orb.Polygon) orb.Polygon {
	if len(poly) == 0 {
		return nil
	}

	projected := make(orb.Polygon, len(poly))

	for i, ring := range poly {

		// project directly via standard 3-stage ForwardProjectRadians
		projectedRing := make(orb.Ring, 0, len(ring))
		for _, pt := range ring {
			lonRad := pt.X() * rhealpix.DegToRad
			latRad := pt.Y() * rhealpix.DegToRad

			// pure projection call: HealpixSphere -> Seam Lock -> combineTriangles
			xRhp, yRhp, _ := d.ForwardProjectRadians(lonRad, latRad)
			projectedRing = append(projectedRing, orb.Point{xRhp, yRhp})
		}

		projected[i] = projectedRing
	}

	return projected
}

// ForwardProjectPolygonForcedRegion projects geographic polygon points through RhpDGGS
// using forced regional context, maintaining seam alignment across all regions.
func ForwardProjectPolygonForcedRegion(d *rhealpix.RhpDGGS, poly orb.Polygon, dom rhealpix.GeographicFacetDomain) orb.Polygon {
	if len(poly) == 0 {
		return nil
	}

	projected := make(orb.Polygon, len(poly))

	for i, ring := range poly {
		projectedRing := make(orb.Ring, 0, len(ring))
		for _, pt := range ring {
			lonDeg := pt.X()
			latRad := pt.Y() * rhealpix.DegToRad
			isOnSeam := math.Abs(math.Abs(lonDeg)-180.0) < 1e-11

			var lonRad float64
			if isOnSeam {
				// lock input longitude sign based on domain context
				if dom.FacetID == 4 || dom.MaxLon == 180.0 {
					lonRad = math.Pi // Eastern sector (Sector R) -> +180°
				} else if dom.FacetID == 1 || dom.MinLon == 180.0 || dom.MinLon == -270.0 {
					lonRad = -math.Pi // Western / Extended Wrap sector (Sector O) -> -180°
				} else {
					if lonDeg >= 0 {
						lonRad = math.Pi
					} else {
						lonRad = -math.Pi
					}
				}
			} else {
				lonRad = lonDeg * rhealpix.DegToRad
			}

			// forward project forcing the region a point is to be placed within
			xRhp, yRhp, _ := d.ForwardProjectForcedRegion(lonRad, latRad, dom.Region)

			// HARD LOCK: enforce planar X_RHP alignment ONLY on Equatorial Facets (1 and 4)
			if isOnSeam {
				if dom.FacetID == 4 {
					xRhp = math.Pi
				} else if dom.FacetID == 1 {
					xRhp = -math.Pi
				}
			}

			projectedRing = append(projectedRing, orb.Point{xRhp, yRhp})
		}
		projected[i] = projectedRing
	}

	return projected
}

// PlanarToLocalUVGeometry transforms an orb.Geometry from global planar radians [x, y]
// into local unit square coordinates [0, 1] x [0, 1] for a given base facet.
func PlanarToLocalUVGeometry(d *rhealpix.RhpDGGS, facetID uint8, geom orb.Geometry) orb.Geometry {
	if geom == nil {
		return nil
	}

	switch g := geom.(type) {
	case orb.Polygon:
		return transformPolygonUV(d, facetID, g)

	case orb.MultiPolygon:
		mp := make(orb.MultiPolygon, len(g))
		for i, poly := range g {
			mp[i] = transformPolygonUV(d, facetID, poly)
		}
		return mp

	default:
		return nil
	}
}

func transformPolygonUV(d *rhealpix.RhpDGGS, facetID uint8, poly orb.Polygon) orb.Polygon {
	out := make(orb.Polygon, len(poly))
	for r, ring := range poly {
		localRing := make(orb.Ring, len(ring))
		for i, pt := range ring {
			// pt[0] = xRhp, pt[1] = yRhp in authalic planar radians
			u, v := d.RadiansToFacetLocal(facetID, pt[0], pt[1])
			localRing[i] = orb.Point{u, v}
		}
		out[r] = localRing
	}
	return out
}
