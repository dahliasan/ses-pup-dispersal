---
title: "Data Access"
layout: default
---

# Data Access

## Complete Datasets
The following datasets are available through Dryad: [DOI: pending]

### Tracking Data
1. **Juvenile Seal Tracks** (`tracks_processed_12h.rds`)
   - 49 individual seals
   - 14,162 locations
   - 12-hour intervals
   - ARGOS locations processed using state-space models
   - Variables: id, date, lon, lat, flag

2. **Adult Female Tracks** (`adult_female_locs.csv`)
   - 67 individual females
   - 101,905 locations
   - 6-hour intervals
   - ARGOS and GLS locations
   - Variables: id, date, lon, lat, type, month, season

### Particle Tracking Data
3. **Deep Particle Traces** (`particle-trace-186.13m.csv`)
   - 48 particle trajectories
   - 5,603 locations
   - Depth: 186.13m
   - Variables: id, date, lon, lat

4. **Surface Particle Traces** (`currently_particleTrace.csv`)
   - 48 particle trajectories
   - 11,164 locations
   - Depth: surface
   - Variables: id, date, lon, lat

### Metadata
5. **Dataset Summary** (`dataset_metadata_summary.csv`)
   - Summary statistics for each dataset
   - Temporal coverage
   - Number of records
   - Dataset descriptions

## Data Processing
- All seal tracks processed using `foieGras` R package
- Particle traces generated using `currently` R package
- Quality control flags indicate suspicious tracks
- All coordinates in decimal degrees (WGS84)

## Citation
If you use this data, please cite:
```bibtex
@article{foo2024seal,
  title={},
  author={Foo, Dahlia},
  journal={},
  year={2024}
}
```

## Contact
For questions about the data, please contact [your contact information or institutional email] 