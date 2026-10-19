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

// Set saveGeoJSONDisk to true when debugging to write GeoJSON files to disk.
const saveGeoJSONDisk = false

type FixtureTestCase struct {
	Name              string
	Centre            orb.Point
	WithHole          bool
	TargetRes         uint8
	ExpectedCellCount int
}

// MakeSyntheticS2Footprint creates a synthetic 100km x 100km WGS84 polygon
// centred at (centreLon, centreLat), with optional interior hole.
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
		{Name: "01_Facet_N_North_Cap", Centre: orb.Point{0.0, 65.0}, ExpectedCellCount: 156837},
		{Name: "02_Facet_O_Equatorial_West", Centre: orb.Point{-135.0, 0.0}, ExpectedCellCount: 44627},
		{Name: "03_Facet_P_Equatorial_MidWest", Centre: orb.Point{-45.0, 0.0}, ExpectedCellCount: 44627},
		{Name: "04_Facet_Q_Equatorial_MidEast", Centre: orb.Point{45.0, 0.0}, ExpectedCellCount: 44627},
		{Name: "05_Facet_R_Equatorial_East", Centre: orb.Point{135.0, 0.0}, ExpectedCellCount: 44627},
		{Name: "06_Facet_S_South_Cap", Centre: orb.Point{0.0, -65.0}, ExpectedCellCount: 156832},

		// --- Equatorial-Polar Cap Straddles (8 Cases) ---
		{Name: "07_Straddle_O_N", Centre: orb.Point{-135.0, capLat}, ExpectedCellCount: 144987},
		{Name: "08_Straddle_P_N", Centre: orb.Point{-45.0, capLat}, ExpectedCellCount: 144987},
		{Name: "09_Straddle_Q_N", Centre: orb.Point{45.0, capLat}, ExpectedCellCount: 144987},
		{Name: "10_Straddle_R_N", Centre: orb.Point{135.0, capLat}, ExpectedCellCount: 144987},
		{Name: "11_Straddle_O_S", Centre: orb.Point{-135.0, -capLat}, ExpectedCellCount: 144987},
		{Name: "12_Straddle_P_S", Centre: orb.Point{-45.0, -capLat}, ExpectedCellCount: 144987},
		{Name: "13_Straddle_Q_S", Centre: orb.Point{45.0, -capLat}, ExpectedCellCount: 144987},
		{Name: "14_Straddle_R_S_55GDP_Style", Centre: orb.Point{146.5, -capLat}, ExpectedCellCount: 116629},

		// --- Equatorial Inter-Facet Longitude Seams (4 Cases) ---
		{Name: "15_Seam_O_P_Minus90", Centre: orb.Point{-90.0, 0.0}, ExpectedCellCount: 85896},
		{Name: "16_Seam_P_Q_PrimeMeridian", Centre: orb.Point{0.0, 0.0}, ExpectedCellCount: 85896},
		{Name: "17_Seam_Q_R_Plus90", Centre: orb.Point{90.0, 0.0}, ExpectedCellCount: 85896},
		{Name: "18_Seam_R_O_AntiMeridian", Centre: orb.Point{180.0, 0.0}, ExpectedCellCount: 85896},

		// --- 3-Facet Tri-Junction Corners (8 Cases) ---
		{Name: "19_Corner_OP_N", Centre: orb.Point{-90.0, capLat}, ExpectedCellCount: 132947},
		{Name: "20_Corner_PQ_N", Centre: orb.Point{0.0, capLat}, ExpectedCellCount: 132947},
		{Name: "21_Corner_QR_N", Centre: orb.Point{90.0, capLat}, ExpectedCellCount: 132947},
		{Name: "22_Corner_RO_N", Centre: orb.Point{180.0, capLat}, ExpectedCellCount: 132947},
		{Name: "23_Corner_OP_S", Centre: orb.Point{-90.0, -capLat}, ExpectedCellCount: 132947},
		{Name: "24_Corner_PQ_S", Centre: orb.Point{0.0, -capLat}, ExpectedCellCount: 132947},
		{Name: "25_Corner_QR_S", Centre: orb.Point{90.0, -capLat}, ExpectedCellCount: 132947},
		{Name: "26_Corner_RO_S", Centre: orb.Point{180.0, -capLat}, ExpectedCellCount: 132947},

		// --- Quad Points & Extremes (3 Cases) ---
		{Name: "27_QuadPoint_North_West", Centre: orb.Point{-180.0, capLat}, ExpectedCellCount: 132947},
		{Name: "28_QuadPoint_South_East", Centre: orb.Point{180.0, -capLat}, ExpectedCellCount: 132947},
		{Name: "29_Unwrapped_AntiMeridian_Cross", Centre: orb.Point{180.05, 10.0}, ExpectedCellCount: 120530},

		// --- Singularities & Hole Topology (2 Cases) ---
		{Name: "30_North_Pole_Enclosure", Centre: orb.Point{0.0, 89.6}, ExpectedCellCount: 109409},
		{Name: "31_Seam_R_S_With_Hole", Centre: orb.Point{146.5, -capLat}, WithHole: true, ExpectedCellCount: 137006},
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

	fmt.Println("Cell Shape:")
	fmt.Println(rhealpix.EllipsoidalCellShape(d, cellID))

	// CellToPolygon formatted output
	poly, err := rhealpixorb.CellToPolygon(d, cellID, 4)
	if err != nil {
		fmt.Printf("Error generating CellToPolygon: %v\n", err)
		return
	}

	fmt.Println("CellToPolygon WKT:")
	fmt.Println(wkt.MarshalString(poly))
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
	if saveGeoJSONDisk {
		if err := os.MkdirAll(outputDir, 0755); err != nil {
			t.Fatalf("failed creating output directory: %v", err)
		}
	}

	targetRes := uint8(13)
	d := rhealpix.NewRhpDGGS(
		rhealpix.WithEllipsoid(el),
		rhealpix.WithNSide(3),
		rhealpix.WithNorthSquare(0),
		rhealpix.WithSouthSquare(0),
		rhealpix.WithMaxResolution(targetRes),
	)

	for _, tc := range tests {
		t.Run(tc.Name, func(t *testing.T) {
			inputPoly := MakeSyntheticS2Footprint(tc.Centre.X(), tc.Centre.Y(), tc.WithHole)

			inputPath := filepath.Join(outputDir, fmt.Sprintf("%s_input.geojson", tc.Name))
			if saveGeoJSONDisk {
				// export synthetic input geometry
				inputFC := geojson.NewFeatureCollection()
				inputFeat := geojson.NewFeature(inputPoly)
				inputFeat.Properties["test_name"] = tc.Name
				inputFeat.Properties["has_hole"] = tc.WithHole
				inputFC.Append(inputFeat)

				if err := writeGeoJSONFile(inputPath, inputFC); err != nil {
					t.Fatalf("failed exporting input GeoJSON for %s: %v", tc.Name, err)
				}
			}

			// run top-down planar decomposition & compaction
			gjsn := geojson.NewGeometry(inputPoly)
			compactedCells, err := rhealpixorb.STACGeometryToCompactedCells(d, *gjsn, targetRes)
			if err != nil {
				t.Fatalf("Decomposition failed for %s: %v", tc.Name, err)
			}

			if len(compactedCells) == 0 {
				t.Fatalf("FAIL: 0 compacted cells produced for valid polygon %s", tc.Name)
			}

			cellsPath := filepath.Join(outputDir, fmt.Sprintf("%s_cells.geojson", tc.Name))
			if saveGeoJSONDisk {
				// export compacted cell geometries
				if err := rhealpixorb.ExportCellsToGeoJSON(d, compactedCells, cellsPath, 0); err != nil {
					t.Fatalf("failed exporting cell GeoJSON for %s: %v", tc.Name, err)
				}
			}

			// cell count assertion
			if tc.ExpectedCellCount > 0 && len(compactedCells) != tc.ExpectedCellCount {
				t.Errorf("%s cell count mismatch: got %d cells, want %d",
					tc.Name, len(compactedCells), tc.ExpectedCellCount)
			}

			if saveGeoJSONDisk {
				fmt.Printf("EXPORTED: %-35s -> Input: %s | Cells: %s (%d cells)\n",
					tc.Name, filepath.Base(inputPath), filepath.Base(cellsPath), len(compactedCells))
			} else {
				fmt.Printf("PASSED: %-35s -> %d compacted cells\n", tc.Name, len(compactedCells))
			}
			os.Stdout.Sync()
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
