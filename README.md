# Atqasuk and Utqiaġvik flower count data and code repository

## Content
This repository contains the code and data necessary to replicate data analysis, figures and tables in the manuscript:

Limits on phenological response to high temperature in the Arctic

## Contacts
Sarah C. Elmendorf; sarah.elmendorf [at] colorado.edu


Robert D. Hollister; Robert Hollister hollistr [at] gvsu.edu

## Data usage guidelines and license 
### Data
The name of Barrow, AK, changed to Utqiaġvik, AK in 2016. These data are derived from a database that predates the name transition, and thus abbreviations and raw data files usually reference the outdated name (Barrow). The Inupiaq name (Utqiaġvik) should be used in any derivative work or publications.

Data used in the analyses, which represent a snapshot of ongoing data collection at the time of analyses and in some instances have had further processing or qa/qc applied, can be found within:

```
/data/BARROW_ATQASUK
```

Please note that raw data for ongoing collections for these datasets (archived in a slightly different format) can be found at:

Robert Hollister and Katlyn Betway-May. 2022. Air and soil temperatures and soil moisture in the International Tundra Experiment (ITEX) plots at Utqiaġvik and Atqasuk, Alaska. Arctic Data Center. doi:10.18739/A2V40K12R, version: urn:uuid:02022a31-97b5-4178-b692-6d2a77c120eb.

Robert Hollister and Katlyn Betway-May. 2022. Plant phenology and performance in the International Tundra Experiment (ITEX) plots at Utqiaġvik and Atqasuk, Alaska. Arctic Data Center. urn:uuid:59f94e1f-93b2-4f8f-b68b-cbf260cd61fc.

### Code
All code provided for data preparation and analysis is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

The main script to reproduce the analyses and figures can be found in:
```
/scripts/analyze_inflflowercts.R
```
which sources script 
```
/scripts/calculate_hourly_climate.R as an input to the downstream analyses
```
