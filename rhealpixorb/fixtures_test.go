package rhealpixorb_test

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/paulmach/orb"
	"github.com/paulmach/orb/encoding/wkt"
	"github.com/paulmach/orb/geojson"
	rhealpix "github.com/sixy6e/go-rhealpix"
	"github.com/sixy6e/go-rhealpix/rhealpixorb"
)

type FixtureTestCase struct {
	Name      string
	Center    orb.Point
	WithHole  bool
	TargetRes uint8
}

// MakeSyntheticS2Footprint creates a synthetic 100km x 100km WGS84 polygon
// centreed at (centreLon, centreLat), with optional interior hole.
func MakeSyntheticS2Footprint(centreLon, centreLat float64, withHole bool) orb.Polygon {
	halfLatSpan := 0.45 // ~50 km in latitude

	// adjust longitude width for meridian convergence away from equator
	cosLat := math.Cos(centreLat * rhealpix.DegToRad)
	if math.Abs(cosLat) < 1e-4 {
		cosLat = 1e-4
	}
	halfLonSpan := 0.45 / cosLat

	minX := centreLon - halfLonSpan
	maxX := centreLon + halfLonSpan
	minY := centreLat - halfLatSpan
	maxY := centreLat + halfLatSpan

	// exterior ring (counter-clockwise)
	outer := orb.Ring{
		{minX, minY},
		{maxX, minY},
		{maxX, maxY},
		{minX, maxY},
		{minX, minY},
	}

	if !withHole {
		return orb.Polygon{outer}
	}

	// optional interior hole (clockwise) ~20km x 20km in centre
	hLon := halfLonSpan * 0.2
	hLat := halfLatSpan * 0.2
	inner := orb.Ring{
		{centreLon - hLon, centreLat - hLat},
		{centreLon - hLon, centreLat + hLat},
		{centreLon + hLon, centreLat + hLat},
		{centreLon + hLon, centreLat - hLat},
		{centreLon - hLon, centreLat - hLat},
	}

	return orb.Polygon{outer, inner}
}

func GetSyntheticTestMatrix() []FixtureTestCase {
	capLat := math.Asin(2.0/3.0) * rhealpix.RadToDeg // ~41.8103149°

	return []FixtureTestCase{
		// --- Pure In-Facet Granules (6 Cases) ---
		{Name: "01_Facet_N_North_Cap", Center: orb.Point{0.0, 65.0}},
		{Name: "02_Facet_O_Equatorial_West", Center: orb.Point{-135.0, 0.0}},
		{Name: "03_Facet_P_Equatorial_MidWest", Center: orb.Point{-45.0, 0.0}},
		{Name: "04_Facet_Q_Equatorial_MidEast", Center: orb.Point{45.0, 0.0}},
		{Name: "05_Facet_R_Equatorial_East", Center: orb.Point{135.0, 0.0}},
		{Name: "06_Facet_S_South_Cap", Center: orb.Point{0.0, -65.0}},

		// --- Equatorial-Polar Cap Straddles (8 Cases) ---
		{Name: "07_Straddle_O_N", Center: orb.Point{-135.0, capLat}},
		{Name: "08_Straddle_P_N", Center: orb.Point{-45.0, capLat}},
		{Name: "09_Straddle_Q_N", Center: orb.Point{45.0, capLat}},
		{Name: "10_Straddle_R_N", Center: orb.Point{135.0, capLat}},
		{Name: "11_Straddle_O_S", Center: orb.Point{-135.0, -capLat}},
		{Name: "12_Straddle_P_S", Center: orb.Point{-45.0, -capLat}},
		{Name: "13_Straddle_Q_S", Center: orb.Point{45.0, -capLat}},
		{Name: "14_Straddle_R_S_55GDP_Style", Center: orb.Point{146.5, -capLat}},

		// --- Equatorial Inter-Facet Longitude Seams (4 Cases) ---
		{Name: "15_Seam_O_P_Minus90", Center: orb.Point{-90.0, 0.0}},
		{Name: "16_Seam_P_Q_PrimeMeridian", Center: orb.Point{0.0, 0.0}},
		{Name: "17_Seam_Q_R_Plus90", Center: orb.Point{90.0, 0.0}},
		{Name: "18_Seam_R_O_AntiMeridian", Center: orb.Point{180.0, 0.0}},

		// --- 3-Facet Tri-Junction Corners (8 Cases) ---
		{Name: "19_Corner_OP_N", Center: orb.Point{-90.0, capLat}},
		{Name: "20_Corner_PQ_N", Center: orb.Point{0.0, capLat}},
		{Name: "21_Corner_QR_N", Center: orb.Point{90.0, capLat}},
		{Name: "22_Corner_RO_N", Center: orb.Point{180.0, capLat}},
		{Name: "23_Corner_OP_S", Center: orb.Point{-90.0, -capLat}},
		{Name: "24_Corner_PQ_S", Center: orb.Point{0.0, -capLat}},
		{Name: "25_Corner_QR_S", Center: orb.Point{90.0, -capLat}},
		{Name: "26_Corner_RO_S", Center: orb.Point{180.0, -capLat}},

		// --- Quad Points & Extremes (3 Cases) ---
		{Name: "27_QuadPoint_North_West", Center: orb.Point{-180.0, capLat}},
		{Name: "28_QuadPoint_South_East", Center: orb.Point{180.0, -capLat}},
		{Name: "29_Unwrapped_AntiMeridian_Cross", Center: orb.Point{180.05, 10.0}},

		// --- Singularities & Hole Topology (2 Cases) ---
		{Name: "30_North_Pole_Enclosure", Center: orb.Point{0.0, 89.6}},
		{Name: "31_Seam_R_S_With_Hole", Center: orb.Point{146.5, -capLat}, WithHole: true},
	}
}

func InspectCellWKT(d *rhealpix.RhpDGGS, cellIDStr string) {
	// parse cell ID string (e.g. N4444444444482)
	cellID, err := rhealpix.ParseCellID64(cellIDStr)
	if err != nil {
		fmt.Printf("Error parsing CellID %s: %v\n", cellIDStr, err)
		return
	}

	// raw geographic vertices from un-projection
	rawVerts, err := d.GeographicVertices64(cellID)
	if err != nil {
		fmt.Printf("Error computing raw vertices: %v\n", err)
		return
	}

	planarVerts := rhealpix.CellPlanarVertices(d, cellID)
	ra := d.GetEllipsoid().RA

	fmt.Printf("\n=== Cell ID: %s ===\n", cellIDStr)
	fmt.Printf("Facet: %d, Level: %d\n", cellID.Facet(), cellID.Resolution())
	fmt.Println("Raw Geographic Vertices:")
	for i, pt := range rawVerts {
		fmt.Printf("  Pt %d: Lon = %.14f, Lat = %.14f\n", i, pt[0], pt[1])
	}
	fmt.Println("Planar Vertices:")
	for i, pt := range planarVerts {
		fmt.Printf("  Pt %d: x = %.14f, y = %.14f\n", i, pt[0], pt[1])
		fmt.Printf("  Pt %d: x = %.14f, y = %.14f\n", i, pt[0]*ra, pt[1]*ra)
	}

	fmt.Println("Cell Shape")
	fmt.Println(rhealpix.EllipsoidalCellShape(d, cellID))

	// CellToPolygon formatted output
	poly, err := rhealpixorb.CellToPolygon(d, cellID, 4)
	if err != nil {
		fmt.Printf("Error generating CellToPolygon: %v\n", err)
		return
	}

	fmt.Println("CellToPolygon WKT:")
	str := wkt.MarshalString(poly)
	fmt.Println(str)
}

func ExportComparisonGeoJSON() error {
	fc := geojson.NewFeatureCollection()

	// Python reference polygon (from rhealpixdggs-py)
	pyPoly := orb.Polygon{
		orb.Ring{
			{-77.14285713436743, 89.99979551794536},
			{-64.28571428119999, 89.99979551794536},
			{-59.999999997078895, 89.99973709433294},
			{-69.99999999438026, 89.99973709433294},
			{-77.14285713436743, 89.99979551794536},
		},
	}
	fPy := geojson.NewFeature(pyPoly)
	fPy.Properties["source"] = "python_reference"
	fPy.Properties["cell_id"] = "N4444444444482"
	fc.Append(fPy)

	// Go generated polygon (from go-rhealpix)
	goPoly := orb.Polygon{
		orb.Ring{
			{-77.14285714503936, 89.99979643134944},
			{-64.28571428702361, 89.99979643134944},
			{-60.000000000512635, 89.99973826871036},
			{-70.0000000008544, 89.99973826871036},
			{-77.14285714503936, 89.99979643134944},
		},
	}
	fGo := geojson.NewFeature(goPoly)
	fGo.Properties["source"] = "go_implementation"
	fGo.Properties["cell_id"] = "N4444444444482"
	fc.Append(fGo)

	data, err := json.MarshalIndent(fc, "", "  ")
	if err != nil {
		return err
	}

	return os.WriteFile("N4444444444482_comparison.geojson", data, 0644)
}

func InspectDartCell(d *rhealpix.RhpDGGS, cellIDStr string) {
	cellID, err := rhealpix.ParseCellID64(cellIDStr)
	if err != nil {
		fmt.Printf("Error parsing CellID %s: %v\n", cellIDStr, err)
		return
	}

	shape := rhealpix.EllipsoidalCellShape(d, cellID)
	nwIdx, _ := rhealpix.NWIndex(d, cellID)
	planarVerts := rhealpix.CellPlanarVertices(d, cellID)

	fmt.Printf("\n=== Cell ID: %s ===\n", cellIDStr)
	fmt.Printf("Facet: %d, Resolution: %d, Shape: %s, NWIndex: %d\n",
		cellID.Facet(), cellID.Resolution(), shape, nwIdx)

	fmt.Println("Planar Vertices [UL, UR, LR, LL]:")
	for i, pv := range planarVerts {
		fmt.Printf("  [%d]: X = %.14f, Y = %.14f\n", i, pv[0], pv[1])
	}

	// full 4-point geographic before trimming
	geoVerts4, _ := rhealpix.EllipsoidalVertices(d, cellID, false)
	fmt.Println("4-Point Geographic Vertices (Ordered from NWIndex):")
	for i, g := range geoVerts4 {
		fmt.Printf("  Pt %d: Lon = %.14f, Lat = %.14f\n", i, g[0], g[1])
	}

	// trimmed 3-point geographic
	geoVerts3, _ := rhealpix.EllipsoidalVertices(d, cellID, true)
	fmt.Println("Trimmed Geographic Vertices:")
	for i, g := range geoVerts3 {
		fmt.Printf("  Pt %d: Lon = %.14f, Lat = %.14f\n", i, g[0], g[1])
	}
}

func TestGenerateAndExportSyntheticFixtures(t *testing.T) {
	el := rhealpix.NewWGS84()
	tests := GetSyntheticTestMatrix()

	outputDir := filepath.Join(".", "test_output_geojson")
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		t.Fatalf("failed creating output directory: %v", err)
	}

	targetRes := uint8(13)
	d := rhealpix.NewRhpDGGS(
		rhealpix.WithEllipsoid(el),
		rhealpix.WithNSide(3),
		rhealpix.WithNorthSquare(0),
		rhealpix.WithSouthSquare(0),
		rhealpix.WithMaxResolution(targetRes),
	)

	// c, _ := rhealpix.ParseCellID64("R55555")
	// cvert := rhealpix.CellPlanarVertices(d, c)
	// wgs84verts, _ := rhealpix.GeographicVertices(d, c)
	// fmt.Printf("R5555 planar vertices: %v\n", cvert)
	// fmt.Printf("R5555 geographic vertices: %v\n", wgs84verts)

	// InspectCellWKT(d, "N4444444444482")
	// ExportComparisonGeoJSON()
	// InspectCellWKT(d, "N44440220")
	// InspectCellWKT(d, "N44442")
	// InspectDartCell(d, "N44442")

	for _, tc := range tests {
		t.Run(tc.Name, func(t *testing.T) {
			inputPoly := MakeSyntheticS2Footprint(tc.Center.X(), tc.Center.Y(), tc.WithHole)
			// ================ Unwrap Logging =====================
			// unwrapped := rhealpixorb.UnwrapPolygon(inputPoly, 0.0)
			// for i, pt := range unwrapped[0] {
			// 	fmt.Printf("Vertex %d: Lon = %f, Lat = %f\n", i, pt.X(), pt.Y())
			// }

			// export synthetic input geometry
			inputFC := geojson.NewFeatureCollection()
			inputFeat := geojson.NewFeature(inputPoly)
			inputFeat.Properties["test_name"] = tc.Name
			inputFeat.Properties["has_hole"] = tc.WithHole
			inputFC.Append(inputFeat)

			inputPath := filepath.Join(outputDir, fmt.Sprintf("%s_input.geojson", tc.Name))
			if err := writeGeoJSONFile(inputPath, inputFC); err != nil {
				t.Fatalf("failed exporting input GeoJSON for %s: %v", tc.Name, err)
			}

			// clip/partition geometry into facet sub-geometries
			// results, err := rhealpixorb.ClipGeometryToFacets(d, inputFeat.Geometry)
			// if err != nil {
			// 	t.Fatalf("ClipGeometryToFacets failed: %v", err)
			// }

			// ==================== BOUNDS LOGGING ====================
			// for i, res := range results {
			// 	bounds := res.Geometry.Bound()
			// 	t.Logf("Piece %d -> Facet %d | Bounds: Lon [%.4f, %.4f], Lat [%.4f, %.4f]",
			// 		i, res.FacetID, bounds.Min.X(), bounds.Max.X(), bounds.Min.Y(), bounds.Max.Y())
			// }

			// run top-down planar decomposition & compaction
			gjsn := geojson.NewGeometry(inputPoly)
			compactedCells, err := rhealpixorb.STACGeometryToCompactedCells(d, *gjsn, targetRes)
			if err != nil {
				t.Fatalf("Decomposition failed for %s: %v", tc.Name, err)
			}

			if len(compactedCells) == 0 {
				t.Fatalf("FAIL: 0 compacted cells produced for valid polygon %s", tc.Name)
			}

			counts := make(map[uint8]uint64)
			for _, c := range compactedCells {
				counts[c.Facet()] += 1
			}
			fmt.Println("Facet counts: ", counts)

			// export compacted cell geometries
			cellsPath := filepath.Join(outputDir, fmt.Sprintf("%s_cells.geojson", tc.Name))
			if err := rhealpixorb.ExportCellsToGeoJSON(d, compactedCells, cellsPath, 16); err != nil {
				t.Fatalf("failed exporting cell GeoJSON for %s: %v", tc.Name, err)
			}

			// assert reasonableness bounds
			// 500000 -> just a guesstimate
			if len(compactedCells) > 500000 {
				t.Errorf("Cell explosion detected in %s: got %d cells", tc.Name, len(compactedCells))
			}

			t.Logf("EXPORTED: %-35s -> Input: %s | Cells: %s (%d cells)",
				tc.Name, filepath.Base(inputPath), filepath.Base(cellsPath), len(compactedCells))
		})
	}
}

func writeGeoJSONFile(path string, fc *geojson.FeatureCollection) error {
	rawJSON, err := fc.MarshalJSON()
	if err != nil {
		return err
	}
	var prettyJSON bytes.Buffer
	if err := json.Indent(&prettyJSON, rawJSON, "", "  "); err != nil {
		return err
	}
	return os.WriteFile(path, prettyJSON.Bytes(), 0644)
}
