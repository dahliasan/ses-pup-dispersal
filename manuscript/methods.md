Methods

Data collection

Sixty-nine weaned elephant seal pups at Macquarie Island (54°30′ S, 158°57′ E) were tagged with ARGOS satellite-relayed data loggers (SRDLs, Sea Mammal Research Unit) during their post-weaning fast in December 1995, 1998, 1999, and 2000. These SRDLs provided ARGOS location estimates categorized into different classes based on accuracy (~0.5–10 km) [37]. The pups were also flipper-tagged (Jumbo Rototag, Dalton Supplies, Henley-on-Thames, UK) and permanently marked by hot iron branding for future identification as part of a long-term demographic study [2,38,39]. Researchers made resights of marked animals every year at Macquarie Island, with the resight dataset spanning from 23 November 1951 to 1 October 2014. Details on how seals were captured, weighed, and handled can be found in Hindell et al. [17].

It is important to note that the weaned pups in our dataset were not selected randomly. Instead, the aim was to obtain a range of weaner sizes [40]. While this may lead to some bias in our results, given our relatively large sample size and the biological context of this research, we believe that our results are robust.

We also obtained 67 adult female post-moult winter migration tracks from Macquarie Island, collected from February to October in 2000–2005 and 2010 [41]. Of these, 30 tracks were obtained by SRDLs and the remaining 37 by geolocation light loggers (GLS; from time-depth recorders, Wildlife Computers, Redmond, USA). GLS tags are accurate to approximately 70–120 km for southern elephant seals [42,43]. Some tracks that did not start at the colony were removed from further analyses.

Processing tracks

Only tracked data from weaned seals that left the colony and provided more than 10 days of ARGOS data were retained for further analysis. The at-sea locations provided by ARGOS were filtered using a correlated random walk state-space model with a 4 m s⁻¹ maximum velocity threshold via the R package aniMotum [45]. We used a 1-day time step to provide daily location estimates, simplifying the dataset while still capturing the essential movement patterns of each seal. Adult female tracks, which were initially received in 6-hour timesteps, were reinterpolated to 1-day intervals using the same methods to ensure consistency in data handling across all datasets.

Particle trace generation

To investigate the influence of ocean currents on seal movements, we generated particle traces corresponding to each seal's first foraging trip. We utilized the currently R package [47] to simulate the paths that Lagrangian particles would take if released from the seals' starting locations and dates. Particle traces were created using two datasets: the Copernicus Global Ocean Physics Reanalysis (providing subsurface currents at approximately 186 meters depth) and the Global Ocean Gridded SSALTO/DUACS Sea Surface Height Level 4 (providing surface currents). Preliminary analyses indicated that seal tracks aligned more closely with particle traces generated from surface currents (see Supplementary Figures 1 and 2). Additionally, the surface current dataset offered more complete spatial and temporal coverage. Consequently, we used the surface current dataset for our primary analysis, as it provided a more comprehensive representation of the ocean currents encountered by the seals during dispersal.

Dispersal analysis

We conducted a comprehensive dispersal analysis to examine the relationship between seal movements and ocean currents. For each seal and its corresponding particle trace, we calculated bearings between sequential locations using the geosphere R package. This yielded a series of bearings representing the direction of movement over time for both seals and particles.

To compare seal and particle movements, we calculated the angular difference between the seal's bearing and the particle's bearing at each time point. We also computed cumulative mean bearings for both seals and particles over time. These methods allowed us to assess both the immediate alignment of seal movements with currents and the overall trend of their travel directions in relation to prevailing currents throughout their journeys.

To classify seals as "following" or "not following" the currents, we compared the overall mean bearing of each seal to that of the corresponding particle trace. Seals with a mean bearing within 45 degrees of the particle trace's mean bearing were classified as "following" currents. This 45-degree threshold balances strict alignment with currents and natural deviations due to factors such as foraging behavior and individual preferences.

First foraging trip and first-year survival

A seal was considered to have survived its first foraging trip if there was either a complete track of this trip or if the seal was seen again at the colony during the resight period, indicating that the tag failed before it returned to the colony during its first foraging trip. Seals were considered to have survived their first year of life if they were resighted at least once after one year of being tagged until the end of the study in 2014. The probability of a seal not returning to the colony due to tag failure alone is considered low, as tag failures are rare and the seals have strong site fidelity to their natal colonies [46].

Environmental variables

We extracted environmental variables along the seals' tracks to include in our analyses. These variables included sea surface temperature (SST), sea surface height anomaly (SSHA), eddy kinetic energy (EKE), terrain ruggedness index (TRI), slope, SST gradient, ice coverage, distance to ice margin, chlorophyll concentration (chl), and chlorophyll gradient (chlgrad). These variables were averaged over the duration of each seal's first foraging trip.

Statistical analysis

Angular data were analyzed using the circular R package [48]. To test whether mean bearings differed between groups of seals, we used Watson's two-sample test of homogeneity [49]. We also used the Rayleigh test of uniformity to determine if there was a preferred direction of travel within certain groups.

Logistic regression models were fitted to test two response variables: survival of the first foraging trip and first-year survival. These response variables were analyzed against predictor variables including whether the seal was following currents, weaning mass, birth year, and the environmental variables mentioned above. Seals with missing weaning mass were excluded from the survival analysis.

Model selection was performed using the dredge function from the MuMIn R package [50], fitting all possible combinations of predictor variables. Models with a delta AICc (Akaike Information Criterion corrected for small sample sizes) of less than or equal to 2 were considered the best-fitting models. We then used model averaging to account for model selection uncertainty and to obtain robust parameter estimates and relative variable importance.

Model assumptions were checked by simulating residuals of the global model using the DHARMa R package [51], assessing issues such as overdispersion, zero-inflation, and outliers.

All statistical analyses were conducted in R (version 4.4.0). We used the tidyverse suite of packages for data manipulation and visualization, sf for spatial data handling, and rnaturalearth for mapping. Results, including summary statistics, test results, model outputs, and plots, were saved to a date-stamped output directory for reproducibility and further analysis.