package rhealpixorb

import (
	// "encoding/json"
	"bytes"
	"encoding/json"
	"fmt"
	"os"

	"github.com/paulmach/orb/geojson"
	rhealpix "github.com/sixy6e/go-rhealpix"
)

// ExportRangesToGeoJSON exports the Min and Max endpoint cells of each Uint64Range
// to a GeoJSON FeatureCollection file for QGIS / GeoJSON.io inspection.
func ExportRangesToGeoJSON(el *rhealpix.Ellipsoid, ranges []Uint64Range, outputPath string) error {
	fc := geojson.NewFeatureCollection()

	for i, r := range ranges {
		// min boundary cell
		minID := rhealpix.CellID64(r.Min)
		if !minID.IsZero() {
			minPoly, err := Cell64ToPolygon(el, minID, 1)
			if err == nil {
				feat := geojson.NewFeature(minPoly)
				feat.Properties["range_index"] = i
				feat.Properties["bound"] = "MIN"
				feat.Properties["cell_id_dec"] = r.Min
				feat.Properties["cell_id_hex"] = fmt.Sprintf("0x%X", r.Min)
				feat.Properties["span_cells"] = r.Max - r.Min
				fc.Append(feat)
			}
		}

		// max boundary cell
		maxID := rhealpix.CellID64(r.Max)
		if !maxID.IsZero() {
			maxPoly, err := Cell64ToPolygon(el, maxID, 1)
			if err == nil {
				feat := geojson.NewFeature(maxPoly)
				feat.Properties["range_index"] = i
				feat.Properties["bound"] = "MAX"
				feat.Properties["cell_id_dec"] = r.Max
				feat.Properties["cell_id_hex"] = fmt.Sprintf("0x%X", r.Max)
				feat.Properties["span_cells"] = r.Max - r.Min
				fc.Append(feat)
			}
		}
	}

	rawJSON, err := fc.MarshalJSON()
	if err != nil {
		return fmt.Errorf("failed marshaling range GeoJSON: %w", err)
	}

	return os.WriteFile(outputPath, rawJSON, 0644)
}

// ExportRangesToGeoJSON128 exports the Min and Max endpoint cells of each Uint128Range
// to a GeoJSON FeatureCollection file for QGIS / GeoJSON.io inspection.
func ExportRangesToGeoJSON128(el *rhealpix.Ellipsoid, ranges []Uint128Range, outputPath string) error {
	fc := geojson.NewFeatureCollection()

	for i, r := range ranges {
		minID := r.MinCell()
		if !minID.IsZero() {
			minPoly, err := Cell128ToPolygon(el, minID, 1)
			if err == nil {
				feat := geojson.NewFeature(minPoly)
				feat.Properties["range_index"] = i
				feat.Properties["bound"] = "MIN"
				feat.Properties["cell_id_dec"] = minID.String()
				feat.Properties["cell_id_hex"] = fmt.Sprintf("0x%016X%016X", r.MinHigh, r.MinLow)
				fc.Append(feat)
			}
		}

		maxID := r.MaxCell()
		if !maxID.IsZero() {
			maxPoly, err := Cell128ToPolygon(el, maxID, 1)
			if err == nil {
				feat := geojson.NewFeature(maxPoly)
				feat.Properties["range_index"] = i
				feat.Properties["bound"] = "MAX"
				feat.Properties["cell_id_dec"] = maxID.String()
				feat.Properties["cell_id_hex"] = fmt.Sprintf("0x%016X%016X", r.MaxHigh, r.MaxLow)
				fc.Append(feat)
			}
		}
	}

	rawJSON, err := fc.MarshalJSON()
	if err != nil {
		return fmt.Errorf("failed marshaling 128-bit range GeoJSON: %w", err)
	}

	return os.WriteFile(outputPath, rawJSON, 0644)
}

// ExportCellsToGeoJSON converts a slice of CellID64 into a GeoJSON FeatureCollection file
// with cell metadata properties for inspection in QGIS / ogrinfo.
func ExportCellsToGeoJSON(el *rhealpix.Ellipsoid, cells []rhealpix.CellID64, outputPath string) error {
	fc := geojson.NewFeatureCollection()

	for _, cell := range cells {
		poly, err := Cell64ToPolygon(el, cell, 2)
		if err != nil {
			return fmt.Errorf("failed converting cell %d to polygon: %w", cell, err)
		}

		feature := geojson.NewFeature(poly)

		minRange, maxRange := cell.SubtreeRangeMax()

		// convert uint64 integers to strings to prevent JSON float precision/parser overflow
		feature.Properties["cell_id_str"] = cell.String()
		feature.Properties["cell_id_uint64"] = fmt.Sprintf("%d", uint64(cell))
		feature.Properties["cell_id_hex"] = fmt.Sprintf("0x%X", uint64(cell))
		feature.Properties["resolution"] = cell.Resolution()
		feature.Properties["facet"] = cell.Facet()
		feature.Properties["range_min"] = fmt.Sprintf("%d", uint64(minRange))
		feature.Properties["range_max"] = fmt.Sprintf("%d", uint64(maxRange))

		fc.Append(feature)
	}

	rawJSON, err := fc.MarshalJSON()
	if err != nil {
		return fmt.Errorf("failed encoding GeoJSON: %w", err)
	}

	var prettyJSON bytes.Buffer
	if err := json.Indent(&prettyJSON, rawJSON, "", "  "); err != nil {
		return fmt.Errorf("failed formatting GeoJSON indent: %w", err)
	}

	if err := os.WriteFile(outputPath, prettyJSON.Bytes(), 0644); err != nil {
		return fmt.Errorf("failed writing file %s: %w", outputPath, err)
	}

	// fmt.Printf("Exported %d compacted cells to %s\n", len(cells), outputPath)
	return nil
}

// ExportCellsToGeoJSON128 converts a slice of CellID128 into a GeoJSON FeatureCollection file
// with cell metadata properties for inspection in QGIS / ogrinfo.
func ExportCellsToGeoJSON128(el *rhealpix.Ellipsoid, cells []rhealpix.CellID128, outputPath string) error {
	fc := geojson.NewFeatureCollection()

	for _, cell := range cells {
		poly, err := Cell128ToPolygon(el, cell, 2)
		if err != nil {
			return fmt.Errorf("failed converting cell %s to polygon: %w", cell.String(), err)
		}

		feature := geojson.NewFeature(poly)

		minRange, maxRange := cell.SubtreeRangeMax()

		feature.Properties["cell_id_str"] = cell.String()
		feature.Properties["cell_id_high"] = fmt.Sprintf("%d", cell.High)
		feature.Properties["cell_id_low"] = fmt.Sprintf("%d", cell.Low)
		feature.Properties["cell_id_hex"] = fmt.Sprintf("0x%016X%016X", cell.High, cell.Low)
		feature.Properties["resolution"] = cell.Resolution()
		feature.Properties["facet"] = cell.Facet()
		feature.Properties["range_min_str"] = minRange.String()
		feature.Properties["range_max_str"] = maxRange.String()

		fc.Append(feature)
	}

	rawJSON, err := fc.MarshalJSON()
	if err != nil {
		return fmt.Errorf("failed encoding GeoJSON: %w", err)
	}

	var prettyJSON bytes.Buffer
	if err := json.Indent(&prettyJSON, rawJSON, "", "  "); err != nil {
		return fmt.Errorf("failed formatting GeoJSON indent: %w", err)
	}

	if err := os.WriteFile(outputPath, prettyJSON.Bytes(), 0644); err != nil {
		return fmt.Errorf("failed writing file %s: %w", outputPath, err)
	}

	// fmt.Printf("Exported %d compacted 128-bit cells to %s\n", len(cells), outputPath)
	return nil
}
