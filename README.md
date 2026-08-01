# go-rhealpix

A high-performance Go implementation of the **rHEALPix Discrete Global Grid System (DGGS)**.

`go-rhealpix` provides memory-efficient spatial indexing, hierarchical cell bit-packing, bounding box/polygon range generation, and lossless spatial compaction for global geospatial datasets (e.g., satellite imagery telemetry, point clouds, and TileDB spatial dimensions).

---

## Overview & Lineage

## Overview, Lineage & References

This library is an idiomatic, high-performance Go port and extension inspired by the reference Python implementation and formal specification by Dr. Dmitry Raichev and Manaaki Whenua – Landcare Research:

* 🐍 **Python Reference Implementation:** [`manaakiwhenua/rhealpixdggs-py`](https://github.com/manaakiwhenua/rhealpixdggs-py)
* 📄 **Original Preprint Paper:** [The rHEALPix Discrete Global Grid System (Raichev, 2014)](https://raichev.net/files/rhealpix_dggs_preprint.pdf)
* 🔬 **IOP Publishing Paper:** [The rHEALPix Discrete Global Grid System (IOP Conf. Series: Earth and Environmental Science)](https://iopscience.iop.org/article/10.1088/1755-1315/34/1/012012/pdf)
* 📊 **Official Landcare Research Paper:** [rHEALPix Discrete Global Grid System Data Store](https://datastore.landcareresearch.co.nz/dataset/rhealpix-discrete-global-grid-system)

### Scope & Initial Implementation Constraints

To optimize for maximum computational performance, simple bit manipulation, and database index friendliness, **this implementation currently focuses on the standard baseline rHEALPix DGGS configuration**:

* **Aperture:** $3 \times 3$ (Aperture 9 nonary subdivision).
* **Polar Cap Alignment:** `north_square = 0` and `south_square = 0`.
  * The North ($N$) and South ($S$) polar caps attach directly to **Facet O** (Equatorial Facet 0).
* **Ellipsoid:** WGS84 (EPSG:4326 authalic coordinates).

---

## Key Features

* **Compact Integer Encoding:**
  * **`CellID64`:** Ultra-fast 64-bit unsigned integer encoding supporting **Resolutions 0 through 14**.
  * **`CellID128`:** 128-bit unsigned integer (`High`/`Low` pair) supporting **Resolutions 0 through 30** for millimeter-level sub-meter precision.
* **Lossless Spatial Compaction:**
  * Generic `Compact[T]` function that recursively collapses complete $3 \times 3$ (9-child) cell blocks into their parent cell bottom-up across resolutions.
* **Database & Cloud Querying (TileDB / S3):**
  * `BoundingBoxToTileDBRanges` & `KRingToTileDBRanges` for $O(\log N)$ 1D range queries over $1\text{D}$ space-filling curves.
* **Vector Integration:**
  * Native conversion to/from [`orb`](https://github.com/paulmach/orb) geometries with boundary densification along facet edges.

---

## Resolution, Scale, and Metrics

Below are the exact grid metrics derived from the reference WGS84 ellipsoid ($10,007,554.68\text{ m}$ base facet width):

| Resolution ($R$) | Cell Width | Cell Area | Storage Type | Typical Use Case |
| :---: | :---: | :---: | :---: | :--- |
| **0** | $\sim 10,007.55\text{ km}$ | $\sim 100,151,148\text{ km}^2$ | `CellID64` | Base Facets ($N, O, P, Q, R, S$) |
| **1** | $\sim 3,335.85\text{ km}$ | $\sim 11,127,905\text{ km}^2$ | `CellID64` | Continental scale |
| **5** | $\sim 41.18\text{ km}$ | $\sim 1,696.07\text{ km}^2$ | `CellID64` | Major metropolitan regions |
| **8** | $\sim 1.525\text{ km}$ | $\sim 2.326\text{ km}^2$ | `CellID64` | Sentinel-2 / Landsat scene blocks |
| **12** | $\sim 18.83\text{ m}$ | $\sim 354.6\text{ m}^2$ | `CellID64` | Suburban building footprints |
| **13** | $\sim 6.28\text{ m}$ | $\sim 39.4\text{ m}^2$ | `CellID64` | Fine urban detail |
| **14** | $\mathbf{\sim 2.09\text{ m}}$ | $\mathbf{\sim 4.38\text{ m}^2}$ | **`CellID64` (Max)** | **Parking space / High-res aerial** |
| **15** | $\sim 69.7\text{ cm}$ | $\sim 0.486\text{ m}^2$ | `CellID128` | Sub-meter drone scans |
| **20** | $\sim 2.87\text{ cm}$ | $\sim 8.24\text{ cm}^2$ | `CellID128` | Precision agriculture / Surveying |
| **30** | $\sim 0.028\text{ mm}$ | $\sim 0.0008\text{ mm}^2$ | `CellID128` | High-density LiDAR point clouds |

> **Note on 64-bit Limit:** `CellID64` packs the 3-bit Facet ID, 4-bit Resolution depth, and 14 sub-cell path digits into 63 bits, capping out at **Resolution 14 ($\approx 2.09\text{ meters}$ width)**. Higher resolutions automatically use `CellID128`.

---

## Topology

```text
   +---+
   | N |           <- North Polar Cap (N) attached above Facet O
   +---+---+---+---+
   | O | P | Q | R |   <- Equatorial Belt (Facets O, P, Q, R)
   +---+---+---+---+
   | S |           <- South Polar Cap (S) attached below Facet O
   +---+
```
