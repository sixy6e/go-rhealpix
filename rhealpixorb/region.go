package rhealpixorb

import (
	"fmt"
	"sort"
	"strings"

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
