# hydrology-trends
## Quantifying trends in discharge in the United States and Canada ##

Trend analysis has been conducted for stations with long records of daily discharge (1960-present) in the United States (United States Geological Survey, USGS) and Canada (Water Survey of Canada, WSC).


# Introduction to trend analysis approach #

**Challenges to estimating changing flood and drought frequency in the United States and Canada**

The seasonality and magnitude of river discharge varies widely in the United States and Canada due to regional and local variations in precipitation timing and intensity, snowpack water storage, topography, and landcover. 

Coherent estimates of changing extreme discharge on a continental scale are hindered by this variability, as rivers in different regions may be changing at different rates or during different seasons. We use unsupervised K-Means clustering to group rivers into with similar annual flood timing, location, and elevation. These groups are called "hydro-regions". Rivers in each hydro-region are analyzed together, allowing a more coherent picture of changes in discharge to emerge.

<p align="center">
  <img src="/hydro-trends-figures/hclusters-fig1.png" width="85%" >
</p>

In addition, significant changes in the timing and magnitude of precipitation have been observed in the United States and Canada, but these changes may not align with conditions necessary to generate a flood or drought. As a result, analysis of annual maximum and/or minimum flow may not show significant changes, despite fundamental changes occurring. To address this challenge, we not only analyze annual changes in high-flow and low-flow frequency, but also quantify changes occurring in each season.

**Significant changes in seasonal extreme discharge have occurred since the mid-20th century**

Seasonal analysis of changes extreme high- and low-flow event frequency show increased likelihood of extreme events in many regions in the United States.

<p align="center">
  <img src="/hydro-trends-figures/hclusters-fig2.png" width="85%" >
</p>









