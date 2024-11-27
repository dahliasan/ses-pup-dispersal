---
title: "Southern Elephant Seal Pup Dispersal Analysis"
format:
  html:
    toc: true
---

# Code Documentation

This repository contains the code used for analyzing southern elephant seal pup dispersal patterns from Macquarie Island.

## Analysis Scripts

### [Dispersal Analysis](code/dispersal-analysis.qmd)
Main analysis including:
- Track preprocessing and interpolation
- Movement pattern analysis
- Survival modeling
- Statistical analysis

### [Data Export](code/export-data.qmd)
Data preparation and export including:
- Seal tracking data (12-hour intervals)
- Adult female locations (6-hour intervals)
- Particle tracking at surface and 186.13m depth
- Quality control and filtering

### [Utility Functions](code/functions.qmd)
Helper functions for:
- Track processing
- Departure date calculations
- Survival analysis
- Data preprocessing

## [Data Access](data.qmd)
- Complete datasets available on Dryad
- Data documentation and metadata
- File descriptions and variables
- Processing methods