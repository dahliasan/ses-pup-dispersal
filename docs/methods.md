Methods

Data Collection and Preprocessing

Sixty-nine weaned elephant seal pups at Macquarie Island (54°30' S, 158°57' E) were tagged with ARGOS satellite-relayed data loggers (SRDL, Sea Mammal Research Unit) during their post-weaning fast in December 1995, 1998, 1999 and 2000. These SRDLs provided ARGOS location estimates which are categorised into different classes based on accuracy (~0.5-10 km) [37]. The pups were also flipper tagged (Jumbo Rototag, Dalton Supplies, Henley-on-Thames, UK) and permanently marked by hot iron branding for future identification as part of a long-term demographic study [2,38,39]. Researchers made resights of marked animals every year at Macquarie Island. The resight dataset spanned from 23 November 1951 - 1 October 2014. Details on how seals were captured, weighed, and handled can be found in Hindell et al. [17].

It is important to note that the weaned pups in our dataset were not selected randomly. Instead, the aim was to get a range of weaner sizes [40]. While this may lead to some bias in our results, given our relatively large sample size and the biological context of this research, we believe our results are robust.

We also obtained 67 adult female post-moult winter migration tracks from Macquarie Island which were collected from February-October in 2000-2005 and 2010 [41]. Out of the 67 tracks, 30 were obtained by SRDLs and the remaining 37 by geolocation light loggers (GLS; from time-depth recorders, Wildlife Computers, Redmond, USA). GLS tags are accurate to ~70-120 km for southern elephant seals [42,43]. Some of these tracks did not start at the colony and were removed from further analyses.

The seal tracking data was processed to provide locations at 12-hour intervals. Concurrently, we utilized particle tracking data to represent ocean currents in the study area. Both datasets were preprocessed and matched temporally and spatially using custom R functions.

Processing tracks

Only tracked data from weaned seals that left the colony and provided > 10 days of ARGOS data were kept for further analysis. Additionally, some seals had missing weaning mass in the dataset (n = 3) – these were also removed from further analyses. The at-sea locations were filtered using a correlated random walk state-space model with a 4 m s⁻¹ max velocity threshold via the R package `aniMotum`.

Particle Trace Generation

We compared two datasets for generating particle traces: the Copernicus Global Ocean Physics Reanalysis (186m depth currents) and the Global Ocean Gridded SSALTO/DUACS Sea Surface Height L4 (surface currents). Preliminary analysis showed that seal tracks aligned more closely with surface particle traces (SI Fig. 1), which also provided more complete data (SI Fig. 2). Consequently, we used the surface current dataset for our analysis, allowing for a more comprehensive representation of ocean currents encountered by seals during dispersal.

Dispersal Analysis

We conducted a comprehensive dispersal analysis to investigate the relationship between seal movements and ocean currents. For each seal and corresponding particle track, we calculated bearings between each sequential location using the `geosphere` package in R. This provided us with a series of bearings representing the direction of movement for both seals and particles over time. We then computed cumulative circular correlations between these seal and particle bearings using the `circular` package to assess the similarity in movement patterns.

To identify critical periods in the seals' dispersal, we determined the maximum correlation and the time to reach this maximum for each individual. We used the median of these times to maximum correlation across all seals as the critical period. For each seal, we calculated the mean bearing during this critical period. We tested whether the mean maximum correlation differed significantly from zero using a one-sample t-test. Additionally, we examined the relationship between the time to maximum correlation and the maximum correlation value using a correlation test.

To account for individual variability and temporal patterns, we fitted a Generalized Additive Mixed Model (GAMM) using the `mgcv` package. The model included a smooth term for days since the start of the journey and a random effect for individual seals.

We compared each seal's mean bearing during the critical period to the overall mean particle bearing to determine if seals were following the prevailing currents. The overall mean particle bearing was calculated using all particle traces in the study area, providing a representation of the general current direction across the entire region. A seal was considered to be "following" if its mean bearing was within 45 degrees of this overall particle mean bearing. This classification method allowed us to assess whether seals were generally aligning their movement with the dominant current patterns in the region, and was later used as a factor in our survival analysis.

Survival Analysis

Survival data for the seals was obtained from the long-term resight dataset. We conducted survival analyses for two periods: the first trip and the first year post-weaning. We used binomial Generalized Linear Models (GLMs) to investigate the relationship between survival and various predictors, including whether the seal was following currents, weaning mass, birth year, and environmental variables such as sea surface temperature, sea surface height anomaly, eddy kinetic energy, terrain ruggedness index, slope, sea surface temperature gradient, ice coverage, distance to ice edge, chlorophyll concentration, and chlorophyll gradient.

Model selection was performed using the dredge function from the `MuMIn` package, which fits all possible combinations of predictor variables. We considered models with a delta AICc (Akaike Information Criterion corrected for small sample sizes) of 2 or less to be the best-fitting models. We then used model averaging to account for model selection uncertainty and to obtain robust parameter estimates and relative variable importance.

Statistical Analysis and Visualization

All statistical analyses were conducted in R (version 4.4.0). We used the `tidyverse` suite of packages for data manipulation and visualization, `sf` for spatial data handling, `patchwork` for combining multiple plots, and `rnaturalearth` for mapping. We generated various plots to visualize the dispersal patterns, including bearings over time, spatial trajectories, and cumulative correlations.

For the survival analysis, we used the `broom` package to extract and summarize model results. We calculated variable importance using the `sw` function from the `MuMIn` package.

All results, including summary statistics, test results, model outputs, and plots, were saved to a date-stamped output directory for reproducibility and further analysis.