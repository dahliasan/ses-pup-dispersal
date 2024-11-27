---
title: Southern Elephant Seal Pup Dispersal Analysis
layout: default
---

# Code Documentation

This repository contains the code used for analyzing southern elephant seal pup dispersal patterns from Macquarie Island.

## Analysis Scripts

### [Dispersal Analysis](code/dispersal-analysis.md)
Main analysis including:
- Track preprocessing and interpolation
- Movement pattern analysis
- Survival modeling
- Statistical analysis

### [Data Export](code/export-data.md)
Data preparation and export including:
- Seal tracking data (12-hour intervals)
- Adult female locations (6-hour intervals)
- Particle tracking at surface and 186.13m depth
- Quality control and filtering

### [Utility Functions](code/functions.md)
Helper functions for:
- Track processing
- Departure date calculations
- Survival analysis
- Data preprocessing

## Data Structure
The analysis uses four main datasets:
1. Juvenile seal tracks (ARGOS, 12-hour intervals)
2. Adult female tracks (ARGOS/GLS, 6-hour intervals)
3. Particle traces at 186.13m depth
4. Surface particle traces

All tracks were processed using state-space models (`foieGras` R package) and interpolated to regular time intervals.