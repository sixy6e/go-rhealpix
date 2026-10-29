package rhealpixorb

import (
	"fmt"
	"sort"
	"strings"

	"github.com/paulmach/orb"
	rhealpix "github.com/sixy6e/go-rhealpix" // Assumes RootFacetChars is exported here
)

// RegionTags holds the formatted spatial identifiers for filenames and metadata.
type RegionTags struct {
	FilenameTag string   `json:"filename_tag"` // e.g., "P081", "P08-Q12", or "N-O-P-Q-R-S"
	MetadataTag []string `json:"metadata_tag"` // e.g., ["P081", "P082"], ["P08", "Q12"], or ["N", "O", "P", "Q", "R", "S"]
}

// DeriveRegionTagsFromCompactCells generates canonical, glob-friendly filename and metadata
// spatial tags from a compacted slice of rHEALPix cell identifiers (CellID64 or CellID128).
//
// The function groups cells by base facet (N, O, P, Q, R, S) and branches by facet coverage:
//   - Single Facet (Local/Regional):
//     FilenameTag: Top 1 coarsest cell (e.g. "P081")
//     MetadataTag: Top 2 coarsest cells (e.g. ["P081", "P082"])
//   - Multi-Facet (2–5 Facets / Seam Straddles):
//     FilenameTag: Top 1 coarsest cell per facet joined by hyphen (e.g. "P08-Q12")
//     MetadataTag: Top 1 coarsest cell per facet slice (e.g. ["P08", "Q12"])
//   - Global Coverage (All 6 Base Facets Present):
//     FilenameTag: Full canonical base facet sequence "N-O-P-Q-R-S"
//     MetadataTag: Slice of all base facets ["N", "O", "P", "Q", "R", "S"]
//
// Output examples:
//
//	Input examples given here are for visual purposes.
//	It would normally be either []CellID64 or or []CellID128.
//	Single-facet scene:
//	  Input:  ["P081", "P082", "P083"]  (string version of packed uint64/uint128 Cell)
//	  Return: RegionTags{FilenameTag: "P081", MetadataTag: ["P081", "P082"]}
//	Facet-straddling swath:
//	  Input:  ["P081", "Q120", "Q121"]  (string version of packed uint64/uint128 Cell)
//	  Return: RegionTags{FilenameTag: "P081-Q120", MetadataTag: ["P081", "Q120"]}
//	Global model dataset:
//	  Input:  ["N", "O", "P", "Q", "R", "S"]  (string version of packed uint64/uint128 Cell)
//	  Return: RegionTags{FilenameTag: "N-O-P-Q-R-S", MetadataTag: ["N", "O", "P", "Q", "R", "S"]}
func DeriveRegionTagsFromCompactCells[T interface {
	rhealpix.CellID
	rhealpix.Cell[T]
}](compactCells []T) (RegionTags, error) {
	if len(compactCells) == 0 {
		return RegionTags{}, fmt.Errorf("compacted cell list cannot be empty")
	}

	// group cells directly by their uint8 facet index (0..5)
	facetMap := make(map[uint8][]T)
	for _, cell := range compactCells {
		if cell.IsZero() {
			continue
		}
		facet := cell.Facet()
		facetMap[facet] = append(facetMap[facet], cell)
	}

	if len(facetMap) == 0 {
		return RegionTags{}, fmt.Errorf("compacted cell list contains only zero values")
	}

	// sort cells within each facet by coarsest resolution first (lowest Resolution() number)
	for facet := range facetMap {
		sort.Slice(facetMap[facet], func(i, j int) bool {
			resI := facetMap[facet][i].Resolution()
			resJ := facetMap[facet][j].Resolution()
			if resI == resJ {
				// lexicographical tie-breaker
				return facetMap[facet][i].String() < facetMap[facet][j].String()
			}
			return resI < resJ
		})
	}

	// extract sorted list of present facet indices (0..5)
	var presentFacets []uint8
	for facet := range facetMap {
		presentFacets = append(presentFacets, facet)
	}
	sort.Slice(presentFacets, func(i, j int) bool {
		return presentFacets[i] < presentFacets[j]
	})

	// scenario A: global coverage (all 6 base facets present: N, O, P, Q, R, S)
	if len(presentFacets) == 6 {
		allFacets := make([]string, len(rhealpix.RootFacetChars))
		for i, b := range rhealpix.RootFacetChars {
			allFacets[i] = string(b)
		}

		return RegionTags{
			FilenameTag: strings.Join(allFacets, "-"), // "N-O-P-Q-R-S" (Gibb 2016 canonical order)
			MetadataTag: allFacets,                    // ["N", "O", "P", "Q", "R", "S"]
		}, nil
	}

	// scenario B: single-facet dataset (in-facet local/regional scene)
	if len(presentFacets) == 1 {
		facetCells := facetMap[presentFacets[0]]

		filenameTag := facetCells[0].String() // top 1 coarsest cell SUID string

		var metadataTag []string
		if len(facetCells) >= 2 {
			metadataTag = []string{facetCells[0].String(), facetCells[1].String()} // top 2 coarsest cells
		} else {
			metadataTag = []string{facetCells[0].String()}
		}

		return RegionTags{
			FilenameTag: filenameTag,
			MetadataTag: metadataTag,
		}, nil
	}

	// scenario C: multi-facet dataset (2 to 5 facets - seam straddle or wide swath)
	var filenameParts []string
	var metadataTag []string

	for _, facet := range presentFacets {
		coarsestCellStr := facetMap[facet][0].String() // top 1 coarsest per facet
		filenameParts = append(filenameParts, coarsestCellStr)
		metadataTag = append(metadataTag, coarsestCellStr)
	}

	return RegionTags{
		FilenameTag: strings.Join(filenameParts, "-"), // e.g. "P08-Q12"
		MetadataTag: metadataTag,                      // ["P08", "Q12"]
	}, nil
}

// DeriveCentroidRegionCode returns the single rHEALPix cell SUID for the geometry's centroid at targetRes.
func DeriveCentroidRegionCode(
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8,
) (string, error) {
	if geom == nil {
		return "", fmt.Errorf("geometry is nil")
	}

	// compute WGS84 bounding box centroid
	center := geom.Bound().Center()

	// projection from (Lon, Lat) to CellID128 (rhealpix space)
	cellID, err := d.ForwardTransform128(center.X(), center.Y(), targetRes)
	if err != nil {
		return "", fmt.Errorf("failed projecting centroid to cell: %w", err)
	}

	return cellID.String(), nil
}

// DeriveBoundingCornerCodes returns the unique rHEALPix cell IDs (max 4, or max 8 if crossing the antimeridian)
// containing the corners of the geometry's WGS84 bounding box at targetRes.
// Each decomposed cell in the string is joined by a hyphen (-).
func DeriveBoundingCornerCodes(
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8,
) (string, []string, error) {
	if geom == nil {
		return "", nil, fmt.Errorf("geometry is nil")
	}

	bound := geom.Bound()
	var corners []orb.Point

	if bound.Min.X() > bound.Max.X() {
		// Antimeridian crossing: split into East [Min.X, 180] and West [-180, Max.X] corner sets
		corners = []orb.Point{
			bound.Min,
			{180.0, bound.Min.Y()},
			{180.0, bound.Max.Y()},
			{bound.Min.X(), bound.Max.Y()},
			{-180.0, bound.Min.Y()},
			{bound.Max.X(), bound.Min.Y()},
			bound.Max,
			{-180.0, bound.Max.Y()},
		}
	} else {
		corners = []orb.Point{
			bound.Min,
			{bound.Max.X(), bound.Min.Y()},
			bound.Max,
			{bound.Min.X(), bound.Max.Y()},
		}
	}

	cellSet := make(map[string]struct{}, len(corners))

	for _, pt := range corners {
		cellID, err := d.ForwardTransform128(pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, fmt.Errorf("failed projecting corner point %v: %w", pt, err)
		}
		cellSet[cellID.String()] = struct{}{}
	}

	suids := make([]string, 0, len(cellSet))
	for suid := range cellSet {
		suids = append(suids, suid)
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
}

// DeriveRegionCode64FromGeometry calculates the canonical region_code string for any orb.Geometry footprint.
func DeriveRegionCode64FromGeometry(d *rhealpix.RhpDGGS, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	vertices, err := extractGeometryVertices(geom)
	if err != nil {
		return "", nil, err
	}

	facetGroups := make(map[uint8][]rhealpix.CellID64)
	for _, pt := range vertices {
		cell, err := d.ForwardTransform64(pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, err
		}
		facetGroups[cell.Facet()] = append(facetGroups[cell.Facet()], cell)
	}

	var suids []string
	for _, groupCells := range facetGroups {
		lcaCell := groupCells[0]
		for _, cell := range groupCells[1:] {
			var lcaErr error
			lcaCell, lcaErr = rhealpix.CommonAncestor64(lcaCell, cell)
			if lcaErr != nil {
				lcaCell = cell
			}
		}
		suids = append(suids, lcaCell.String())
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
}

// DeriveRegionCode128FromGeometry calculates the canonical region_code string using 128-bit precision.
func DeriveRegionCode128FromGeometry(d *rhealpix.RhpDGGS, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	vertices, err := extractGeometryVertices(geom)
	if err != nil {
		return "", nil, err
	}

	facetGroups := make(map[uint8][]rhealpix.CellID128)
	for _, pt := range vertices {
		cell, err := d.ForwardTransform128(pt.X(), pt.Y(), targetRes)
		if err != nil {
			return "", nil, err
		}
		facetGroups[cell.Facet()] = append(facetGroups[cell.Facet()], cell)
	}

	var suids []string
	for _, groupCells := range facetGroups {
		lcaCell := groupCells[0]
		for _, cell := range groupCells[1:] {
			var lcaErr error
			lcaCell, lcaErr = rhealpix.CommonAncestor128(lcaCell, cell)
			if lcaErr != nil {
				lcaCell = cell
			}
		}
		suids = append(suids, lcaCell.String())
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
}

// DeriveRegionCode automatically selects 64-bit or 128-bit derivation based on resolution depth.
func DeriveRegionCode(d *rhealpix.RhpDGGS, geom orb.Geometry, targetRes uint8) (string, []string, error) {
	if targetRes <= rhealpix.MaxResolution64 {
		return DeriveRegionCode64FromGeometry(d, geom, targetRes)
	}
	return DeriveRegionCode128FromGeometry(d, geom, targetRes)
}

// DeriveRegionCodeAtRes calculates a fixed-resolution spatial cover string for any orb.Geometry footprint.
// Each decomposed cell in the string is joined by a hyphen (-).
// Output examples for target resolution level 3 (1 facet + 3 levels):
// O003-O006
// S005-S006-S008
// S000-S002-S006-S008
// Output examples for target resolution level 2 (1 facet + 2 levels):
// S03-S04
// S03-S04-S06-S07
// O07
// S03-S04-S06-S07
func DeriveRegionCodeAtRes(
	d *rhealpix.RhpDGGS,
	geom orb.Geometry,
	targetRes uint8, // e.g. 3 or 4 for Landsat
) (string, []string, error) {
	// walk/decompose geometry up to targetRes
	cells, err := STACGeometryToTileDBRangesTopDown128(d, geom, targetRes)
	// cells, err := STACGeometryToTileDBRangesTopDown(el, geom, targetRes)  // testing equivalencce
	if err != nil {
		return "", nil, err
	}

	// extract unique cell string representations at targetRes
	suidSet := make(map[string]struct{})
	for _, cell := range cells {
		suidSet[cell.String()] = struct{}{}
	}

	suids := make([]string, 0, len(suidSet))
	for suid := range suidSet {
		suids = append(suids, suid)
	}

	sort.Strings(suids)
	return strings.Join(suids, "-"), suids, nil
}
