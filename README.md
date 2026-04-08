# CrIS_VIIRS_collocation-master

Collocation library for matching CrIS (Cross-track Infrared Sounder) L1B sounder
footprints with VIIRS (Visible Infrared Imaging Radiometer Suite) 750-m imager
pixels. Supports both the **SNPP** (Suomi NPP) and **JPSS-1** (NOAA-20)
satellites.

This repository implements the sounder/imager spatial matching algorithm
(Wang et al. 2016) used to produce the following MEaSUREs WVCC products for
NASA GES DISC:

- **SNPP CrIS-VIIRS 750-m Matchup Indexes V1** (DOI: [10.5067/MEASURES/WVCC/DATA211](https://doi.org/10.5067/MEASURES/WVCC/DATA211))
- **JPSS-1 CrIS-VIIRS 750-m Matchup Indexes V1** (DOI: [10.5067/MEASURES/WVCC/DATA212](https://doi.org/10.5067/MEASURES/WVCC/DATA212))

## Contents

```
CrIS_VIIRS_collocation-master/
├── code_test_QY.py       # Entry point: call_match_cris_viirs() - end-to-end collocation + NetCDF writer
├── code_test.py          # Legacy test harness (SNPP VIIRS HDF5 format)
├── geo_QY.py             # Core geometry + NetCDF readers for NASA CrIS L1B / VIIRS VNP03MOD/VJ103MOD
├── post_process/         # Post-production utilities (metadata insertion, quality checks)
├── total_nc_size.py      # Utility: aggregate output volume
└── README.md
```

This repo is a library — it does not run production directly. The production
runner lives in [`matchup_pge`](https://github.com/JPL-WVCC/matchup_pge). See
`parallel_run_matchup.py` for how these modules are invoked.

## Algorithm overview

The collocation uses the method described in Wang et al. (2016):

1. **Read CrIS geolocation** from the L1B NetCDF4 file
   (`geo_QY.read_nasa_cris_geo`): lat, lon, satellite azimuth/range/zenith,
   observation time (TAI93), and longwave radiance.
2. **Read VIIRS geolocation** from VNP03MOD (SNPP) or VJ103MOD (J1) NetCDF4
   files (`geo_QY.read_nasa_viirs_geo`) spanning 3 consecutive days: lat, lon,
   satellite azimuth/range/zenith, height, time.
3. **Convert to ECEF** (Earth-Centered Earth-Fixed) coordinates via
   `geo_QY.LLA2ECEF`, then compute the sounder FOV position in ENU via
   `geo_QY.RAE2ENU` and back to ECEF via `geo_QY.ENU2ECEF`.
4. **KD-tree spatial matching** (`pykdtree.kdtree.KDTree` inside
   `geo_QY.match_cris_viirs_QY`) finds VIIRS pixels within each CrIS footprint,
   filtered by a 900-second time tolerance.
5. **Write output** as NetCDF4 in `code_test_QY.call_match_cris_viirs`, with
   spacecraft-aware metadata (SNPP vs J1 DOI, SHORT_NAME, TITLE) and per-FOV
   indexes into the VIIRS granules.

## Output format

Output NetCDF4 files contain indexes pointing from each CrIS FOV to the
matching VIIRS pixels. Dimensions:

- `GranuleCount_ImagerPixel` — total matched VIIRS pixels
- `sounder_atrack`, `sounder_xtrack` — CrIS scan line / footprint
- `sounder_fov` — 9 fields of view per CrIS footprint (3x3 array)

Spacecraft-specific metadata:

| Attribute              | SNPP                              | JPSS-1                            |
|------------------------|-----------------------------------|-----------------------------------|
| `SHORT_NAME`           | `SNPP_CrIS_VIIRS750m_IND`         | `J1_CrIS_VIIRS750m_IND`           |
| `TITLE`                | SNPP CrIS-VIIRS 750-m Matchup...  | JPSS-1 CrIS-VIIRS 750-m Matchup...|
| `IDENTIFIER_PRODUCT_DOI` | `10.5067/MEASURES/WVCC/DATA211` | `10.5067/MEASURES/WVCC/DATA212`   |

The underlying collocation algorithm is identical for SNPP and J1 — only the
output metadata and input file patterns differ.

Each output directory also contains a `.dataset.json` and `.met.json` sidecar
for GES DISC ingest.

## Input data

### CrIS L1B

| Satellite | Path on AIRS/SIPS                                       | File pattern                             |
|-----------|---------------------------------------------------------|------------------------------------------|
| SNPP      | `/peate_archive/NPPOps/snpp/gdisc/2/{YYYY}/{MM}/{DD}/crisl1b/` | `SNDR.SNPP.CRIS.*.L1B.std.v*.nc`         |
| JPSS-1    | `/peate_archive/NPPOps/jpss1/gdisc/3/{YYYY}/{MM}/{DD}/crisl1b/` | `SNDR.J1.CRIS.*.L1B.std.v03_*.nc` (v3)   |

**Note:** JPSS-1 CrIS L1B v3 uses a different `time_coverage_start` format
(`YYYY-MM-DDTHH:MM:SS.00Z`) than v2 (`YYYY-MM-DDTHH:MM:SSZ`). The production
runner in `matchup_pge` handles both; this library consumes whatever the
runner passes in.

### VIIRS geolocation

| Satellite | Product    | LAADS Collection | Typical location                   |
|-----------|------------|------------------|------------------------------------|
| SNPP      | VNP03MOD   | 5110             | `/raid15/leipan/VIIRS/VNP03MOD/{YYYY}/{DOY}/` |
| JPSS-1    | VJ103MOD   | 5201             | `{user_home}/measures/VIIRS/VJ103MOD/{YYYY}/{DOY}/` |

VIIRS data is not pre-staged on AIRS/SIPS. Use
[`matchup_pge/scripts/download_viirs.py`](https://github.com/JPL-WVCC/matchup_pge/blob/master/scripts/download_viirs.py)
to fetch from LAADS DAAC.

## Dependencies

- Python 3.8+
- `numpy`, `netCDF4`, `pykdtree`

```bash
pip install -r requirements.txt
```

Or use the pre-built environment on AIRS/SIPS:

```bash
export PATH=/home/leipan/anaconda3/bin:$PATH
```

## Usage

This library is invoked from the production runner in `matchup_pge`. Direct
usage example:

```python
import sys
sys.path.append('/path/to/CrIS_VIIRS_collocation-master')

from code_test_QY import call_match_cris_viirs

cris_files = ['SNDR.J1.CRIS.20240115T0000.m06.g001.L1B.std.v03_08.G.240115070053.nc']
viirs_files = ['VJ103MOD.A2024015.0000.021.2024015155121.nc', ...]
product_root = '/path/to/output/J1_CrIS_VIIRS/'
spacecraft = 'J1'   # or 'SNPP'

call_match_cris_viirs(cris_files, viirs_files, product_root, spacecraft)
```

## Performance note

For each CrIS granule, the production runner opens ~800 VIIRS files (3 days x
~288 files/day) to read `time_coverage_start/end` for time filtering before
invoking this library. At production scale (~240 CrIS granules/day) that is
~200,000 VIIRS file opens per day of production. Running against remote data
(OPeNDAP, S3) is therefore impractical — local disk access is required.

## References

- Wang, L., B. Chen, Q. Yue, and E. Fetzer (2016), *Combining AIRS, MODIS, and
  CALIOP observations to identify contrail and contrail-induced cirrus in the
  upper troposphere*, J. Geophys. Res. Atmos., 121, 11,547–11,568.
- SNPP product landing page:
  https://disc.gsfc.nasa.gov/datasets/SNPP_CrIS_VIIRS_IND_1/summary
- J1 product landing page:
  https://disc.gsfc.nasa.gov/datasets/J1_CrIS_VIIRS_IND_1/summary

## See also

- [`matchup_pge`](https://github.com/JPL-WVCC/matchup_pge) — production runner and VIIRS download tooling
- [`AIRS_MODIS_collocation-master`](https://github.com/JPL-WVCC/AIRS_MODIS_collocation-master) — sibling library for Aqua AIRS-MODIS collocation
