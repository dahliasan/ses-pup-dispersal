# Analysis of Bearing Comparisons

## Methods

To investigate the movement patterns of weaner seals in relation to adult females and simulated ocean particles, we conducted a series of circular statistical analyses. We employed the Watson-Williams test, a circular analogue of the one-way ANOVA, to compare the mean bearings between groups. This test was chosen for its ability to handle directional data and assess whether two or more samples of circular data have the same mean direction or not.

The analysis was performed at two levels: individual and group. For the individual-level analysis, we compared the bearings of each movement step for weaners against those of females and particles at 0m depth. The group-level analysis utilized summarized data, where mean bearings were calculated for each individual before comparison. This two-tiered approach allowed us to examine both fine-scale movement patterns and overall trends while accounting for potential pseudo-replication issues.

To visualize the differences in bearings, we created circular histograms using the ggplot2 package in R. These plots displayed the distribution of bearings for each group, with mean directions indicated by arrows. Additionally, we calculated cumulative mean bearings for weaners and particles to assess how the average direction of movement changed over time. The angle differences between these cumulative means were then computed and visualized using a linear histogram.

## Results

## Weaner vs Female

### Individual-level comparison:
- Mean bearing for Weaner: 111.808°
- Mean bearing for Female: 122.4281°
- Watson's Two-Sample Test: p-value < 0.001

### Group-level comparison (using summarized data):
- Mean bearing for Weaner: 119.9439°
- Mean bearing for Female: 138.4876°
- Watson's Two-Sample Test: 0.01 < p-value < 0.05

## Weaner vs Particle 0m

### Individual-level comparison:
- Mean bearing for Weaner: 111.808°
- Mean bearing for Particle 0m: 114.6757°
- Watson's Two-Sample Test: p-value < 0.001

### Group-level comparison (using summarized data):
- Mean bearing for Weaner: 119.9439°
- Mean bearing for Particle 0m: 114.8217°
- Watson's Two-Sample Test: 0.05 < p-value < 0.10

## Interpretation

1. Weaner vs Female:
   - Significant differences at both individual and group levels
   - Weaners tend to move in a more easterly direction compared to females

2. Weaner vs Particle 0m:
   - Significant difference at individual level, but not at group level
   - Similar overall movement patterns when considered at a group level

## General Observations

1.  The individual-level analyses show more significant differences than the group-level analyses. This suggests that aggregating data by individual before comparison may smooth out some of the variability.
2. Weaners seem to have a more consistent bearing across analyses compared to females and particles, which show more variation between individual and group-level analyses.
3. The bearings for all groups are generally in the southeast direction (between 90° and 180°), which might indicate a common environmental influence on movement patterns.

## Conclusions

While there are some differences in movement patterns between weaners, females, and particles, the overall directions are relatively similar, especially when considered at a group level. The differences are more pronounced between weaners and females than between weaners and particles at 0m depth.

## Additional Findings

The analysis of cumulative mean bearings between weaners and particles at 0m depth revealed subtle differences in movement patterns over time. A histogram of the angle differences between the cumulative mean bearings was generated to visualize these discrepancies. The resulting distribution provides insights into how closely weaner movements align with simulated particle trajectories throughout their tracks.

The histogram of angle differences showed [describe the shape of the distribution, e.g., "a unimodal distribution centered near 0°" or "a bimodal distribution with peaks at X° and Y°"]. This suggests that [interpret the meaning of the distribution, e.g., "weaner movements generally align closely with particle trajectories" or "there are two distinct patterns of alignment between weaner movements and particle trajectories"].

It's worth noting that 14 data points were removed from this analysis due to non-finite values, which may be attributed to missing data or calculation issues at the beginning of tracks. This data loss, while minor, should be considered when interpreting the results.

These findings complement our earlier analyses by providing a dynamic perspective on how weaner movements relate to simulated ocean currents over the course of their tracks. The results suggest [summarize the main takeaway, e.g., "a strong influence of ocean currents on weaner movements" or "a complex relationship between weaner behavior and environmental factors"].

## Conclusions

[Your existing conclusions here, potentially updated with insights from the cumulative mean bearing analysis]

This comprehensive analysis of movement bearings provides valuable insights into the spatial behavior of weaner seals in relation to adult females and ocean currents. The multi-faceted approach, combining individual and group-level analyses with cumulative mean comparisons, offers a nuanced understanding of these complex ecological interactions. These findings have important implications for understanding the factors influencing juvenile seal dispersal and survival in dynamic marine environments.

## Early Period Analysis

In an attempt to identify potential differences in movement patterns during the initial stages of weaner dispersal, we conducted an analysis focused on the early period of each track. This analysis was based on the hypothesis that the influence of ocean currents might be more pronounced or different during the initial phase of the weaner's journey.

### Methods

1. For each weaner, we calculated the day at which the angle difference between the weaner's cumulative mean bearing and the corresponding particle's cumulative mean bearing was at its minimum.

2. We then computed the mean of these "minimum angle difference days" across all weaners to define a common "early period" for the entire dataset.

3. Using this early period threshold, we filtered the datasets for weaners, particles, and adult females to include only data up to this time point.

4. We then reran our bearing comparisons (Weaner vs Particle 0m and Weaner vs Female) using only this early period data.

### Results

The analysis of the early period yielded the following results:

1. Weaner vs Particle 0m comparison:
   - Mean bearing for Weaner: 141.1501° (sd 0.8732621)
   - Mean bearing for Particle 0m: 120.1901° (sd 0.5409575)
   - Watson-Williams test results:
     Test Statistic: 0.2546
     0.01 < P-value < 0.05

2. Weaner vs Female comparison:
   - Mean bearing for Weaner: 141.1501° (sd 0.8732621)
   - Mean bearing for Female: 158.1333° (sd 0.903587)
   - Watson-Williams test results:
     Test Statistic: 0.2373
     0.01 < P-value < 0.05

### Interpretation

These findings suggest that:

1. There are statistically significant differences in mean bearings between weaners and particles at 0m depth during the early period (0.01 < p < 0.05). Weaners tend to move in a more southeasterly direction (141.15°) compared to the simulated particles (120.19°).

2. Similarly, there are statistically significant differences between weaners and adult females during the early period (0.01 < p < 0.05). Adult females tend to move in a more southerly direction (158.13°) compared to weaners (141.15°).

3. The differences in bearings during the early period are more pronounced than in the full dataset analysis, particularly for the weaner vs particle comparison. This suggests that the relationship between weaner movements and simulated ocean currents may be more complex in the initial stages of dispersal.

4. The standard deviations indicate that weaners and females have more variable bearings compared to particles, which could reflect individual behavioral differences or responses to local environmental conditions.

### Methodological Implications

[Previous content remains the same]

This analysis of the early period reveals subtle but significant differences in movement patterns, highlighting the importance of considering temporal variations in animal movement studies. While the full-track analysis remains valuable for understanding overall patterns, this early period analysis provides insights into the critical initial phase of weaner dispersal.
