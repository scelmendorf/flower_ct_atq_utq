


# -- SETUP ----
# clear workspace
rm(list = ls())
lapply(paste("package:", names(sessionInfo()$otherPkgs), sep = ""), detach, character.only = TRUE, unload = TRUE)
# so that all the joins dont give warnings
options(stringsAsFactors = FALSE)
options(warn = 1)

# libraries
library(tidyverse)
library(NMOF) # for grid search
library(ggrepel) # for plots
library(scales) # for plots
library(cowplot) # for plots

# gamm4 also required by coded by name in funcs
# library (dplyr)
# library (stringr)
# library (data.table)
# NMOF
# Metrics
# scales
# viridis
# slider
# grid
# gridExtra
# egg


# options -
# set to false if have already run once and just want to tweak figs
# the overlap section takes a long time to run
run_overap <- TRUE

# -- READ INFLOCT DATA ----
# inflocts as dumped from sql db lacks headers so is read in 2 lines

flct <- read.delim("data/BARROW_ATQASUK/tblInfloCountsInfl.txt", sep = "\t", header = TRUE)
names(flct) <- read.delim("data/BARROW_ATQASUK/tblInfloCounts_headers.txt", sep = "\t", header = TRUE) %>% names()

# name decoder
taxon <- read.csv("data/BARROW_ATQASUK/taxon.csv", header = TRUE)

# fix a few counts of 0.5, set to 0 as poisson model (and reality) doesn't
# understand a flower count of 0.5. This happens once in
# PBISB, 2014, AD, CTL
# and presumably is a typo
flct <- flct %>%
  mutate(
    AvgOfAvgOfnumResponse =
      ifelse(AvgOfAvgOfnumResponse == 0.5, 0, AvgOfAvgOfnumResponse)
  ) %>%
  # use modern site name Utqiavuk
  mutate(strSitCom = case_when(
    strSitCom == "BD" ~ "UD",
    strSitCom == "BW" ~ "UW",
    TRUE ~ strSitCom
  ))


# -- READ PHENOLOGY FIRST EVENT DATA FOR LEADING 0s ----
# We will assume given fairly frequent revisit frequency at the beginning
# of the season that 2 d before the first flower there were 0 flowers
# at the site. This allows us to add an extra "no flowers yet" to some
# subsite/years where the first flower count recorded is a positive number
# and is consistent with the field methods used for the dataPhenology tables

phen <- read.csv("data/BARROW_ATQASUK/dataPhenology.csv", header = TRUE)

# rename like flct
phen %<>% mutate(
  type = stringr::str_sub(strIdPl, 6, 7),
  strSexx = str_sub(strIdPl, 5, 5),
  strGSpp = str_sub(strIdPl, 1, 4),
  strGSppSex = str_sub(strIdPl, 1, 5),
  strSitCom = str_sub(strIdLo, 1, 2),
  strTrea = str_sub(strIdLo, 3, 3)
) %>%
  rename(numYear = strYear) %>%
  # use only total plot counts, single individuals might not be flowering
  # but what's relevant here is if there's NO flowers anywhere
  filter(type == "TP") %>%
  select(numYear, strSexx, strGSpp, strGSppSex, strTrea, numJulFl, strIdLo, strSitCom) %>%
  # remove a few typo outliers
  filter(!is.na(numJulFl) & numJulFl > 100 & numJulFl < 300) %>%
  # assume 2d before first count, 0 flowers
  mutate(numJulian = numJulFl - 2, AvgOfAvgOfnumResponse = 0) %>%
  select(-numJulFl) %>%
  # use modern site name Utqiavuk
  mutate(strSitCom = case_when(
    strSitCom == "BD" ~ "UD",
    strSitCom == "BW" ~ "UW",
    TRUE ~ strSitCom
  ))

# -- READ SNOWFREE DATA ----
all_sf <- list()
for (s in c("ATQASUK.AD.RH", "ATQASUK.AW.RH", "BARROW.BD.RH", "BARROW.BW.RH")) {
  all_sf[[s]] <- read.csv(paste0("data/BARROW_ATQASUK/", s, ".sf.csv"))
}
all_sf %<>% data.table::rbindlist(., fill = TRUE) %>%
  mutate(subsite = gsub("^ATQASUK\\.|^BARROW\\.|\\.RH$", "", subsite)) %>%
  rename(strIdLo = plot, numYear = year, doy_sf = doy) %>%
  # use modern site name Utqiavuk
  mutate(subsite = case_when(
    subsite == "BD" ~ "UD",
    subsite == "BW" ~ "UW",
    TRUE ~ subsite
  ))

# -- READ AND INFILL HOURLY CLIMATE DATA DATA ----
# this script will generate an object "infilled_hourly" which is a dataframe
# of infilled hourly data. Very few hours needed to be infilled, just used
# linear interpolation of hours before/after or 24hrs before/after

source("scripts/calculate_hourly_climate.R")
infilled_hourly <- infilled_hourly %>%
  # use modern site name Utqiavuk
  mutate(strSitCom = case_when(
    strSitCom == "BD" ~ "UD",
    strSitCom == "BW" ~ "UW",
    TRUE ~ strSitCom
  ))

# -- ADD LEADING 0s from first even to flct ----
# determine which ones have the first ct as something >0

need_leading_0s <- flct %>%
  ungroup() %>%
  group_by(strSitCom, strTrea, numYear, strGSppSex, numJulian) %>%
  summarize(ct = sum(AvgOfAvgOfnumResponse)) %>%
  ungroup() %>%
  filter(ct == 0 & !is.na(ct)) %>%
  group_by(numYear, strSitCom, strTrea, strGSppSex) %>%
  summarize(start_0 = min(numJulian)) %>%
  full_join(
    .,
    flct %>%
      ungroup() %>%
      group_by(strSitCom, strTrea, numYear, strGSppSex, numJulian) %>%
      summarize(ct = sum(AvgOfAvgOfnumResponse)) %>%
      ungroup() %>%
      filter(ct != 0 & !is.na(ct)) %>%
      group_by(numYear, strSitCom, strTrea, strGSppSex) %>%
      summarize(start_any = min(numJulian))
  ) %>%
  rowwise() %>%
  mutate(probs = ifelse(start_any < start_0, TRUE, FALSE)) %>%
  filter(probs == TRUE)

# infill a count of 0 for all subplots 2 d before first event
# using the individual plot surveys for flowering pres/abs
phen <- phen %>%
  ungroup() %>%
  inner_join(
    .,
    need_leading_0s %>%
      ungroup() %>%
      select(strSitCom, strTrea, numYear, strGSppSex) %>%
      distinct()
  ) %>%
  group_by(strSitCom, strTrea, numYear, strGSppSex) %>%
  summarize(numJulian = min(numJulian))

phen <- phen %>%
  # add to all the plots in that year
  left_join(
    .,
    flct %>%
      ungroup() %>%
      select(
        strSitCom, strTrea, numYear, strIdLo, strGSpp, strSexx,
        strGSppSex, strTreaYear
      ) %>%
      distinct()
  ) %>%
  mutate(AvgOfAvgOfnumResponse = 0)


# append these rows to the noninfilled flct data, some renaming required
flct <- bind_rows(
  flct, phen
) %>%
  rename(subsite = strSitCom) %>%
  mutate(strTrea = ifelse(strTrea == "C", "CTL", "OTC")) %>%
  distinct()

# View some summary stats
sum_stats <- flct %>%
  group_by(
    subsite,
    strTrea, numYear, strGSpp, strSexx, strGSppSex
  ) %>%
  summarize(
    tot = sum(AvgOfAvgOfnumResponse),
    nplot = length(unique(strIdLo))
  ) %>%
  left_join(., flct %>%
    group_by(
      subsite,
      strTrea, numYear, strGSpp, strSexx, strGSppSex
    ) %>%
    summarize(tot = sum(AvgOfAvgOfnumResponse)) %>%
    ungroup() %>%
    group_by(
      subsite,
      strTrea, strGSpp, strSexx, strGSppSex
    ) %>%
    summarize(mean_flct = mean(tot))) %>%
  mutate(perc_of_avg_yr = tot / mean_flct)


# -- FILTER to min required totals ----
# require at least 10 total flowers in a year*sub to make a sensible curve
flct <- flct %>%
  anti_join(., flct %>%
    group_by(
      subsite,
      strTrea, numYear, strGSpp, strSexx, strGSppSex
    ) %>%
    summarize(tot = sum(AvgOfAvgOfnumResponse)) %>%
    filter(tot < 10))


# -- SUPPLMENTAL cleaning ----
# year*spp specific comments from R. Hollister on partial years or where
# notes/data indicate incomplete sampling
flct <- flct %>%
  filter(numYear != 2001) %>%
  filter(!(subsite == "AW" & strGSpp == "CRAR" & numYear == 2011)) %>%
  filter(!(subsite == "BD" & strGSpp == "CTET" & numYear %in% c(2001, 2015))) %>%
  filter(!(subsite == "BD" & strGSpp == "SLAE" & numYear %in% c(2014, 2018))) %>%
  filter(!(subsite == "BW" & strGSpp == "DFIS" & numYear %in% c(2000, 2001, 2018)))


# remove 2nd curves where RH says flowers
# were missed in earlier surveys based on inspection/field notes
flct <- flct %>%
  filter(!(subsite == "AD" & numYear == 2016 & strGSpp == "DLAP")) %>%
  filter(!(subsite == "AW" & numYear == 1999 & strGSpp == "CAQU" & strJulian > 200)) %>%
  filter(!(subsite == "AW" & numYear == 2011 & strGSpp == "CAQU" & strJulian > 200)) %>%
  filter(!(subsite == "AW" & numYear == 2016 & strGSpp == "CAQU" & strJulian > 205)) %>%
  filter(!(subsite == "AW" & numYear == 2017 & strGSpp == "CAQU")) %>%
  filter(!(subsite == "AW" & numYear == 1999 & strGSpp == "EANG" & strJulian > 200)) %>%
  filter(!(subsite == "AW" & numYear == 2011 & strGSpp == "EANG" & strJulian > 200)) %>%
  filter(!(subsite == "AW" & numYear == 2012 & strGSpp == "EANG" & strJulian > 200))

# -- CALCULATE trt*subsites*year mean SF date ----
# because curves are fit at the total ct per subsite level
# plots w/o counts are a combination of spp not present and 0 flowering
# individuals but ok to aggregate to a subsite scale and assume a constant counting
# intensity across space and time as that's the protocol even if not all the 0
# plots were entered

sf <- all_sf %>%
  group_by(subsite, treatment, numYear) %>%
  summarize(
    doy_sf = round(mean(doy_sf))
  ) %>%
  rename(strTrea = treatment)


# -- MERGE phen with hourly temp and sf data----
# 2019 is missing some late season met data so omitted
# note also based in inspection of clim data a few yrs have
# obvious sensor probs for airT data:
# - but these yrs not in flwcts anyhow
# 1994
# 1996 AD/AW
# 2004

# calculate some basic summary stats of the temperatures
# helpful to interpreting differences among sites, and used later in plots
temp_summry_stats <- infilled_hourly %>%
  filter(!numYear %in% c(1994, 2004, 2019)) %>%
  filter(!(numYear == 1996 & strSitCom %in% c("AD", "AW"))) %>%
  filter(numMonth %in% c(6, 7, 8))

# Calculate 95% of hottest hours for reference line in figues
AD_thresh <- quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "AD" & temp_summry_stats$strTrea == "C"], 0.95)
AW_thresh <- quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "AW" & temp_summry_stats$strTrea == "C"], 0.95)
UD_thresh <- quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "UD" & temp_summry_stats$strTrea == "C"], 0.95)
Uw_thresh <- quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "UW" & temp_summry_stats$strTrea == "C"], 0.95)

flct <- flct %>%
  group_by(subsite, strTrea, numYear, strGSpp, strSexx, strGSppSex, numJulian) %>%
  summarize(
    ct = sum(AvgOfAvgOfnumResponse),
    # per RH protocol - always 24 plots counted
    plot_ct = 24
  ) %>%
  ungroup() %>%
  # for join to hourly temp data, we will assume the counts occurred
  # at noon so sum GDDs up to noon on the day of counting
  mutate(numHr = 12) %>%
  full_join(
    .,
    # join to temp data that has been subset to years with flowering data
    infilled_hourly %>%
      ungroup() %>%
      rename(subsite = strSitCom, meanT = AvgofnumTemp) %>%
      mutate(strTrea = ifelse(strTrea == "C", "CTL", "OTC")) %>%
      select(subsite, numYear, numJulian, numHr, strTrea, meanT) %>%
      # add in species, sex
      right_join(., flct %>% select(subsite, strTrea, numYear, strGSpp, strSexx, strGSppSex) %>% distinct())
  ) %>%
  # get rid of climate only years
  filter(!is.na(strGSpp)) %>%
  # add in mean doy sf over sampled plots
  left_join(., sf) %>%
  # make sf binary, 1 if has snow on it, else 0
  # this is really snowfull, not snowfree...I should have renamed
  mutate(sf = ifelse(doy_sf >= numJulian, 1, 0)) %>%
  # 121 is april 30, gets rid of some NAs in clim
  # that won't matter to the gdd sums but cause probs in summing
  # 2019 no complete clim data yet
  filter(numJulian > 121 & numYear != 2019) %>% # min sf doy over all treats
  arrange(subsite, strTrea, numYear, numJulian, numHr, strGSpp, strSexx)

# -- DEFINE functions to calculated DD sums in different ways

# using this method
# DDs only accumulate when there is no snow
# when no snow, the accumulate as 0 when temp is negative,
# temp when 0<temp<threshold
# threshold when temp exceeds threshold

DD_func <- function(temp, thresh, sf) {
  if (temp < 0 | sf == 1) {
    val <- 0
  }
  else if (temp >= 0 & temp < thresh & sf == 0) {
    val <- temp
  } else if (temp >= thresh & sf == 0) {
    val <- thresh
  }
  return(val)
}

# -- DEFINE data subsets for analysis

full_season_sampling_ctl <- flct %>%
  ungroup() %>%
  filter(strTrea == "CTL") %>%
  filter(numYear != 1998) %>% # dropped not enough counts in this year
  group_by(subsite, strTrea, strGSpp, strSexx) %>%
  summarize(yearct = length(unique(numYear))) %>%
  filter(yearct > 9) %>%
  select(-yearct) %>%
  left_join(., flct %>% filter(numYear != 1998)) %>%
  droplevels()

full_season_sampling_otc <- flct %>%
  ungroup() %>%
  filter(strTrea == "OTC") %>%
  filter(numYear != 1998) %>%
  # keep only spp used in the ctl fitting
  full_join(., full_season_sampling_ctl %>%
    ungroup() %>%
    select(subsite, strGSpp, strSexx) %>%
    distinct())

# backdating cts to midpoints
full_season_sampling_midpoint_ctl <-
  full_season_sampling_ctl %>%
  ungroup() %>%
  filter(!is.na(ct)) %>%
  group_by(subsite, strTrea, strGSpp, strSexx, strGSppSex, numYear) %>%
  mutate(interval = c(7, diff(numJulian))) %>%
  rowwise() %>%
  mutate(numJulian_backdate = numJulian - round(interval / 2), numHr = 12) %>%
  select(subsite, strTrea, strGSpp, strSexx, strGSppSex, numYear, numJulian_backdate, numHr, ct) %>%
  rename(numJulian = numJulian_backdate) %>%
  ungroup() %>%
  full_join(
    ., # add to clim data
    full_season_sampling_ctl %>%
      ungroup() %>%
      select(
        subsite, strTrea, strGSpp, strSexx, numYear, strGSppSex,
        numJulian, numHr, meanT, sf
      )
  )

full_season_sampling_midpoint_otc <-
  full_season_sampling_otc %>%
  ungroup() %>%
  filter(!is.na(ct)) %>%
  group_by(subsite, strTrea, strGSpp, strSexx, strGSppSex, numYear) %>%
  mutate(interval = c(7, diff(numJulian))) %>%
  rowwise() %>%
  mutate(numJulian_backdate = numJulian - round(interval / 2), numHr = 12) %>%
  select(subsite, strTrea, strGSpp, strSexx, strGSppSex, numYear, numJulian_backdate, numHr, ct) %>%
  rename(numJulian = numJulian_backdate) %>%
  ungroup() %>%
  full_join(
    ., # add to clim data
    full_season_sampling_otc %>%
      ungroup() %>%
      select(
        subsite, strTrea, strGSpp, strSexx, numYear, strGSppSex,
        numJulian, numHr, meanT, sf
      )
  )

# summary stats for reporting in paper, nyrs, n subsites, nyears
length(unique(full_season_sampling_midpoint_ctl$numYear))


# -- DEFINE function to estimate overlap####
# to avoid recoding, this can do a lot of different things
# not set up with a lot of warnings so you can definitely set incompatible
# parameters and get errors without helpful warnings
# based on the overlap package function overlap
# but here using splines instead of just curve fitting, to account for
# sometimes incomplete sampling at the end of the year

overlap_stats_GDD <- function(x, thresh = NULL, method = "DD_func",
                              opt = "AIC",
                              s4only = FALSE,
                              predictionmode = FALSE,
                              input_model = NULL,
                              scale_effort = TRUE,
                              scale_input = TRUE,
                              ...) {
  # @param x A named list, each list element will be compared to all others
  # each element should be a dataframe wtih columns
  # numJulian (if method=doy); numHr, ct (count of flowers on that date),
  # meanT (hourly temps), sf (0 if has snow, 1 if snowfree) - only
  # need meanT and sf if/when they are needed for the method selected
  # need continuous data - no gaps in hourly temp/sf if the temp models

  # @param thresh real scalar. Max temperature at which GDD hours will accumulate
  # see DD func above to understand how this works
  # @param method string scalar. Method for accumulating degree days. Must be one
  # of 'DD_func'  # or 'doy)
  # @param opt AIC of OV, AIC used currently to make gamm models, OV will just plot overlap
  # @param s4only logical scalar. Return the poisson gam only (retain S4 class for subsequent function)
  # @param predictionmode logical scalar whether to just use it to make predicted cts for an input
  # dataset based on an input_model
  # @param input_model - should be an s4 class output from a poisson GAM with the relevant inputs
  # @param scale_effort assume cts need to be modeled accounting for different intercensus interval
  # @param scale_input, normalize the predictor variables before running model?

  if (is.null(names(x))) {
    names(x) <- paste("Y", 1:length(x), sep = "")
  }
  dd_raw <- list()
  DD <- list()

  OV <- OV_fitted <- DD_fitted <- xpoints <- xpoints_fitted <- COMPTITLE <- FUNC <- FUNC_FITTED <- poisfit <- NULL

  # calc boundaries for fitting:
  if (opt == "OV") {
    if (method == "doy") {
      boundaries <- lapply(x, function(y) {
        mutate(y %>% rowwise(), DDs = NA) %>%
          ungroup() %>%
          arrange(numJulian, numHr) %>%
          mutate(DD_sum = numJulian) %>%
          filter(!is.na(ct)) %>%
          summarize(min = min(DD_sum), max = max(DD_sum))
      }) %>%
        data.table::rbindlist(.) %>%
        summarise(max = max(max), min = min(min))
    }
    if (method == "DD_func") {
      boundaries <- lapply(x, function(y) {
        mutate(y %>% rowwise(), DDs = DD_func(temp = meanT, thresh = thresh, sf = sf)) %>%
          ungroup() %>%
          arrange(numJulian, numHr) %>%
          mutate(DD_sum = cumsum(DDs)) %>%
          filter(!is.na(ct)) %>%
          summarize(min = min(DD_sum), max = max(DD_sum))
      }) %>%
        data.table::rbindlist(.) %>%
        summarise(max = max(max), min = min(min))
    }
    from <- boundaries$min
    to <- boundaries$max
  }

  # loop over set to fit densities:
  for (j in 1:length(x)) {
    if (method == "doy") {
      iter <- x[[j]] %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          DDs = NA,
          group = names(x)[j]
        ) %>%
        ungroup() %>%
        arrange(numJulian, numHr) %>%
        mutate(DD_sum = numJulian) %>%
        filter(!is.na(ct)) %>%
        rename(x = DD_sum, y = ct)
      if (scale_effort) {
        iter <- iter %>%
          mutate(effort = c(7, diff(numJulian)))
      } else {
        iter <- iter %>%
          mutate(effort = 1)
      }
    }
    if (method == "DD_func") {
      iter <- x[[j]] %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          DDs = DD_func(temp = meanT, thresh = thresh, sf = sf),
          group = names(x)[j]
        ) %>%
        ungroup() %>%
        arrange(numJulian, numHr) %>%
        mutate(DD_sum = cumsum(DDs)) %>%
        filter(!is.na(ct)) %>%
        rename(x = DD_sum, y = ct)
      if (scale_effort) {
        iter <- iter %>%
          mutate(effort = c(7, diff(numJulian)))
      } else {
        iter <- iter %>%
          mutate(effort = 1)
      }
    }

    # just make a normalized one of the raw data, proportion at each
    # sampling point
    iter$ystd <- (iter$y / iter$effort) / sum(iter$y / iter$effort)
    if (opt == "OV") {
      f_raw <- approxfun(iter$x, iter$y / iter$effort, ties = "ordered")
      # normalize
      normcurve <-
        f_raw(seq(from, to)) / sum(f_raw(seq(from, to)), na.rm = TRUE)
      normcurve[is.na(normcurve)] <- 0
      FUNC <- c(FUNC, list(normcurve))
    }
    dd_raw[[j]] <- iter
  }
  if (opt == "OV") {
    # calculate overlap from this:
    for (i1 in 1:(length(x) - 1)) {
      for (i2 in (i1 + 1):(length(x))) {
        comptitle <- paste0(names(x)[i1], "-", names(x)[i2])
        # calc overlap based on linear approx
        dd2 <- data.frame(
          x = seq(from, to), y1 = FUNC[[i1]],
          y2 = FUNC[[i2]]
        )
        # set to 0 outside range
        dd2[is.na(dd2)] <- 0
        # find min and max over eachpoint in the curve
        dd2$ovy <- apply(
          dd2[, c("y1", "y2")],
          1, min
        )
        dd2$ally <- apply(dd2[, c("y1", "y2")],
          1, max,
          na.rm = TRUE
        )
        dd2$dominance <- NA
        dd2$dominance <- ifelse(dd2$y1 > dd2$y2, 1, 2) # which is bigger
        dd2$dominance[dd2$y1 == dd2$y2] <- NA
        dd2$k <- comptitle
        OV <- c(OV, sum(dd2$ovy, na.rm = TRUE) / sum(dd2$ally,
          na.rm = TRUE
        ))
        dd2 <- dd2[order(dd2$x), ]
        CHANGE <- dd2 %>% filter(!is.na(dominance))
        if (nrow(CHANGE) > 1) {
          CHANGE <- CHANGE$x[which(CHANGE$dominance[2:nrow(CHANGE)] !=
            CHANGE$dominance[1:(nrow(CHANGE) - 1)])]
          xpoints <- c(xpoints, list(CHANGE))
        } else {
          xpoints <- c(xpoints, list(NA))
        }
        DD <- rbind(DD, dd2)
        COMPTITLE <- c(COMPTITLE, comptitle)
      }
    }

    names(xpoints) <- names(OV) <- COMPTITLE
  }
  dd_raw <- dd_raw %>% data.table::rbindlist()
  dd_raw$fGroup <- factor(dd_raw$group)
  if (opt == "AIC" & !predictionmode) {
    poisfit <- gamm4::gamm4(y ~ s(x) + offset(log(effort)),
      random = ~ (1 | fGroup),
      family = poisson(), data = dd_raw, weights = NULL,
      subset = NULL, knots = NULL, drop.unused.levels = TRUE,
      REML = FALSE, control = NULL, start = NULL, verbose = 0L
    )
    AIC <- AIC(poisfit$mer)
  }

  if (s4only & opt == "AIC") {
    return(poisfit)
  } else if (predictionmode) {
    newd1 <- dd_raw
    newd1$fGroup <- NULL ## remove predictor - otherwise you get probs with checks
    # if (pkg == 'gamm4'){
    preds <- predict(input_model$gam, newd1,
      type = "response"
    )
    # }
    return(dd_raw %>% mutate(fitted = preds))
  } else {
    return(list(
      DD = data.frame(DD), dd_raw = data.frame(dd_raw), OV = OV,
      DD_fitted = DD_fitted, OV = OV,
      OV_fitted = OV_fitted, xpoints = xpoints, xpoints_fitted = xpoints_fitted,
      AIC = ifelse(!is.null(AIC), AIC, NA)
    )) # ,
  }
} # end function definition


# -- RUN the overlaps function to estimate overlap####
# you can pretty much skip all this if you don't want to remake the overlaps
# it does take a while to go
if (run_overap) {
  out_AIC_spline_gamm4 <- list()

  options(warn = 2) # to fail on warnings
  # sometimes won't converge
  # otherwise will just skip some
  # seems to now work ok w the error handling
  # iterate over everythign

  # running with a max of 20deg as realistically there are hardly any hours above that
  # and including up to 40 does not change the results at all

  # first just fit with all years to get the optimal threshold

  # use gamm4 rather than mcgv sce per the advice that's ok to compare aics more
  # strictly advice

  # note there would be a lot of ways to make this more efficient, e.g.
  # don't define the function inside the loop, parallelize, calc DD's outside
  # but using the "set computer to run and ignore" method of time-saving, left
  # as is. Have patience...

  for (s in unique(full_season_sampling_midpoint_ctl$subsite)) {
    df <- full_season_sampling_midpoint_ctl %>% filter(subsite == s)
    for (spp in unique(df$strGSpp)) {
      df1 <- df %>%
        filter(strGSpp == spp) %>%
        ungroup()
      for (sex in unique(df1$strSexx)) {
        df2 <- df1 %>%
          filter(strSexx == sex) %>%
          mutate(numYear = as.factor(numYear)) %>%
          droplevels()
        df2 <- base::split(df2, df2$numYear)
        # every once in a while the mcgv gam won't solve
        # set AIC to Inf for these, those are not optimal
        min.overlap_aic_spline_gamm4 <- function(par) {
          res <- tryCatch(overlap_stats_GDD(
            x = df2, thresh = par, method = "DD_func",
            opt = "AIC", scale_effort = TRUE,
            scale_input = FALSE
          ),
          error = function(err) {
            NULL
          }
          )
          if (!is.null(res)) {
            return(res$AIC)
          } else {
            return(Inf)
          }
        }

        # Iterate AIC over all values 0.5C to 20C by 0.2 for thresholds
        # note this takes an extremely long time
        out_AIC_spline_gamm4[[paste(s, spp, sex, sep = "_")]] <- tryCatch(NMOF::gridSearch(min.overlap_aic_spline_gamm4, list(seq(0.5, 20, by = 0.2))),
          error = function(err) {
            NA
          }
        )
      }
    }
  }

  # save to avoid rerunning slow code above
  saveRDS(out_AIC_spline_gamm4, "calculations/inflo_out_AIC_spline_gamm4.Rdata")

  # define spp that have local minima for plotting###############
  # list to DF the AIC models, then plot


  out_AIC_spline_gamm4 <- lapply(out_AIC_spline_gamm4, function(x) {
    if (all(is.na(x))) {
      NULL
    }
    else {
      data.frame(
        thresh = unlist(x$levels), y = x$values,
        min_AIC = x$minfun,
        best_thresh = x$minlevels
      )
    }
  }) %>%
    data.table::rbindlist(., idcol = "sub_spp_sex") %>%
    separate(sub_spp_sex,
      sep = "_", remove = FALSE,
      into = c("sub", "spp", "sex")
    ) %>%
    unite("spp_sex", spp, sex, remove = FALSE) %>%
    left_join(., taxon) %>%
    mutate(
      subname =
        case_when(
          sub == "UW" ~ "Utqiagvik - Wet Meadow",
          sub == "UD" ~ "Utqiagvik - Dry Heath",
          sub == "AD" ~ "Atqasuk - Dry Heath",
          sub == "AW" ~ "Atqasuk - Wet Meadow",
          TRUE ~ "NA"
        )
    )

  # these determind byviz inspection after making plots, but added here to alter
  # linetype for final figs

  has_local_minima_spline_for_plotting <- tibble(
    sub_spp_sex = c(
      "AD_LCON_B",
      "AD_PBIS_B",
      "AW_PVIV_B",
      "AW_EANG_B",
      "UD_CTET_B",
      "UD_LARC_B",
      "UD_LCON_B",
      "UD_PHUL_B",
      "UD_PHYP_B",
      "UD_SLAE_B",
      "UD_SPUN_B",
      "UD_SROT_F",
      "UD_SROT_M",
      "UW_CSTA_B",
      "UW_DFIS_B",
      "UW_ETRI_B",
      "UW_HPAU_B",
      "UW_LARC_B",
      "UW_SHIE_B"
    ),
    linetype = "solid"
  )

  # makes plots S3
  for (subsite in unique(out_AIC_spline_gamm4$sub)) {
    # nice title for plots
    subname <- case_when(
      subsite == "UD" ~ expression(paste("Utqia", dot(g), "vik - Dry")),
      subsite == "UW" ~ expression(paste("Utqia", dot(g), "vik - Wet")),
      subsite == "AD" ~ expression("Atqasuk - Dry"),
      subsite == "AW" ~ expression("Atqasuk - Wet")
    )

    panel_lab <- case_when(
      subsite == "UD" ~ "c",
      subsite == "UW" ~ "d",
      subsite == "AD" ~ "a",
      subsite == "AW" ~ "b"
    )

    minfunc_plot <- ggplot(
      out_AIC_spline_gamm4 %>% filter(sub == subsite) %>%
        left_join(., has_local_minima_spline_for_plotting) %>%
        mutate(lineformat = ifelse(is.na(linetype), "dotdash", linetype)) %>%
        mutate(
          y = ifelse(!is.finite(y), NA, y),
          sci_name = ifelse(sci_name == "Eriophorum angustifolium ssp. triste",
            "E. angustifolium ssp. triste", sci_name
          )
        ),
      aes(x = thresh, y = y)
    ) +
      geom_line() +
      facet_wrap(~sci_name, scales = "free") +
      geom_vline(aes(xintercept = best_thresh, linetype = lineformat), color = "red") +
      scale_linetype_identity() +
      geom_hline(aes(yintercept = min_AIC),
        color = "blue",
        alpha = 0.4
      ) +
      xlab("Maximum temperature threshold used in model (°C)") +
      ylab("AIC") +
      theme(strip.text.x = element_text(size = 8)) +
      ggtitle(subname)

    minfunc_plot <- cowplot::plot_grid(minfunc_plot,
      ncol = 1,
      labels = panel_lab, label_x = 0, label_y = 1
    )

    if (subsite != "AW") {
      ggsave(minfunc_plot,
        file = paste0("plots/", subsite, "_inflo_gamm_spline_minfunc_plot.jpg"),
        width = 6, height = 6, dpi = 1200
      )
    } else { # adjust aspect ratio for 6 rather than 7-9 panel plot
      ggsave(minfunc_plot,
        file = paste0("plots/", subsite, "_inflo_gamm_spline_minfunc_plot.jpg"),
        width = 6, height = 4.2, dpi = 1200
      )
    }
  }


  # -- LOO CV for CIs ----
  # drop out one year for sequentially
  # takes a long time to run

  out_AIC_jack_spline_gamm4 <- list()
  for (s in unique(full_season_sampling_midpoint_ctl$subsite)) {
    df <- full_season_sampling_midpoint_ctl %>% filter(subsite == s)
    for (spp in unique(df$strGSpp)) {
      df1 <- df %>% filter(strGSpp == spp)
      for (sex in unique(df1$strSexx)) {
        df2 <- df1 %>%
          filter(strSexx == sex) %>%
          mutate(numYear = factor(numYear)) %>%
          droplevels()
        df2 <- base::split(df2, df2$numYear)
        for (l in 1:length(df2)) {
          min.overlap_aic <- function(par) {
            res <- tryCatch(overlap_stats_GDD(
              x = df2[-l], thresh = par, method = "DD_func",
              opt = "AIC",
              scale_effort = TRUE, scale_input = FALSE
            ),
            error = function(err) {
              NULL
            }
            )
            if (!is.null(res)) {
              return(res$AIC)
            } else {
              return(Inf)
            }
          }
          out_AIC_jack_spline_gamm4[[paste(s, spp, sex, names(df2)[[l]], sep = "_")]] <- tryCatch(NMOF::gridSearch(min.overlap_aic, list(seq(0.5, 20, by = 0.2))),
            error = function(err) {
              NA
            }
          )
        }
      }
    }
  }

  # save the results so if wanting to tweak figures not nec to rerun
  saveRDS(out_AIC_jack_spline_gamm4, "calculations/inflo_out_AIC_jack_spline_gamm4.Rdata")
} # end run overlap section


# -- READ IN OVERLAP OUTPUTS ----
# or you could just read them in to tweak the plots
out_AIC_spline_gamm4 <- readRDS("calculations/inflo_out_AIC_spline_gamm4.Rdata")
out_AIC_jack_spline_gamm4 <- readRDS("calculations/inflo_out_AIC_jack_spline_gamm4.Rdata")


# -- USE CTL PARAMS TO MODEL OTC AND CTL CTS ----
# take best params from AIC on ctls,
# calculate overlap stats on paired OTC vs ctl years
# fit all with ctl fitting set
# make predictions of proportion that show
# up on each date
# ctls and otcs both


best_mods_spline_gamm4 <- lapply(out_AIC_spline_gamm4, function(x) {
  if (all(is.na(x))) {
    NULL
  }
  else {
    data.frame(
      thresh = unlist(x$levels), y = x$values,
      min_AIC = x$minfun,
      best_thresh = x$minlevels
    )
  }
}) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex") %>%
  select(sub_spp_sex, best_thresh) %>%
  distinct()

# poisson model predictions without ranefs for years
# then take proportions of total for year counted on each date to normalize
otc_ctl_overlap_best_gamm4 <- list()
otc_ctl_overlap_default_gamm4 <- list()
otc_ctl_overlap_doy_gamm4 <- list()
ctl_ctl_overlap_best_gamm4 <- list()
ctl_ctl_overlap_default_gamm4 <- list()
ctl_ctl_overlap_doy_gamm4 <- list()

# loop over all combinations
for (s in unique(full_season_sampling_midpoint_ctl$subsite)) {
  df <- full_season_sampling_midpoint_ctl %>% filter(subsite == s)
  for (spp in unique(df$strGSpp)) {
    df1 <- df %>%
      filter(strGSpp == spp) %>%
      ungroup()
    for (sex in unique(df1$strSexx)) {
      df2_ctl <- df1 %>%
        filter(strSexx == sex & strTrea == "CTL") %>%
        mutate(numYear = factor(numYear)) %>%
        droplevels()
      df2_ctl <- base::split(df2_ctl, df2_ctl$numYear)
      df2_otc <- full_season_sampling_midpoint_otc %>%
        ungroup() %>%
        filter(subsite == s & strSexx == sex & strGSpp == spp)
      df2_otc <- base::split(df2_otc, df2_otc$numYear)
      par_spline_gamm4 <- best_mods_spline_gamm4 %>%
        filter(sub_spp_sex == paste0(c(s, spp, sex), collapse = "_")) %>%
        pull(best_thresh)


      predmod_best_gamm4 <- overlap_stats_GDD(
        x = df2_ctl, thresh = par_spline_gamm4, method = "DD_func",
        opt = "AIC", s4only = TRUE, # pkg='gamm4',
        scale_effort = TRUE,
        scale_input = FALSE
      ) # , gaussian=FALSE)

      df2_ctl_best_gamm4 <- overlap_stats_GDD(
        x = df2_ctl, thresh = par_spline_gamm4, method = "DD_func",
        opt = "AIC", s4only = FALSE, predictionmode = TRUE,
        input_model = predmod_best_gamm4,
        scale_effort = TRUE,
        scale_input = FALSE
      )

      df2_otc_best_gamm4 <- overlap_stats_GDD(
        x = df2_otc, thresh = par_spline_gamm4, method = "DD_func",
        opt = "AIC", s4only = FALSE, predictionmode = TRUE,
        input_model = predmod_best_gamm4,
        scale_effort = TRUE,
        scale_input = FALSE
      )

      # default, set thresh to 40C which is max so no temps above here
      predmod_default_gamm4 <- overlap_stats_GDD(
        x = df2_ctl, thresh = 40, method = "DD_func",
        opt = "AIC", s4only = TRUE,
        scale_effort = TRUE,
        scale_input = FALSE
      )

      df2_ctl_default_gamm4 <- overlap_stats_GDD(
        x = df2_ctl, thresh = 40, method = "DD_func",
        opt = "AIC", s4only = FALSE, predictionmode = TRUE,
        input_model = predmod_default_gamm4,
        scale_effort = TRUE,
        scale_input = FALSE
      )
      df2_otc_default_gamm4 <- overlap_stats_GDD(
        x = df2_otc, thresh = 40, method = "DD_func",
        opt = "AIC", s4only = FALSE, predictionmode = TRUE,
        input_model = predmod_default_gamm4,
        scale_effort = TRUE,
        scale_input = FALSE
      )

      # doy
      predmod_doy_gamm4 <- overlap_stats_GDD(
        x = df2_ctl, thresh = NULL, method = "doy",
        opt = "AIC", s4only = TRUE,
        scale_effort = TRUE,
        scale_input = FALSE
      )

      df2_ctl_doy_gamm4 <- overlap_stats_GDD(
        x = df2_ctl, thresh = NULL, method = "doy",
        opt = "AIC", s4only = FALSE, predictionmode = TRUE,
        input_model = predmod_doy_gamm4,
        scale_effort = TRUE,
        scale_input = FALSE
      )

      df2_otc_doy_gamm4 <- overlap_stats_GDD(
        x = df2_otc, thresh = NULL, method = "doy",
        opt = "AIC", s4only = FALSE, predictionmode = TRUE,
        input_model = predmod_doy_gamm4,
        scale_effort = TRUE,
        scale_input = FALSE
      )

      # rmse on proportion of total that show up at each sampling point - wmax
      tst_ctl_best_gamm4 <- df2_ctl_best_gamm4 %>%
        group_by(numYear) %>%
        mutate(fitted_prop = fitted / sum(fitted))

      ctl_ctl_overlap_best_gamm4[[paste0(c(s, spp, sex), collapse = "_")]] <-
        Metrics::rmse(tst_ctl_best_gamm4$ystd, tst_ctl_best_gamm4$fitted_prop)

      tst_otc_best_gamm4 <- df2_otc_best_gamm4 %>%
        group_by(numYear) %>%
        mutate(fitted_prop = fitted / sum(fitted))

      otc_ctl_overlap_best_gamm4[[paste0(c(s, spp, sex), collapse = "_")]] <-
        Metrics::rmse(tst_otc_best_gamm4$ystd, tst_otc_best_gamm4$fitted_prop)

      # rmse on proportion of total that show up at each sampling point - default
      tst_ctl_default_gamm4 <- df2_ctl_default_gamm4 %>%
        group_by(numYear) %>%
        mutate(fitted_prop = fitted / sum(fitted))

      ctl_ctl_overlap_default_gamm4[[paste0(c(s, spp, sex), collapse = "_")]] <-
        Metrics::rmse(tst_ctl_default_gamm4$ystd, tst_ctl_default_gamm4$fitted_prop)

      tst_otc_default_gamm4 <- df2_otc_default_gamm4 %>%
        group_by(numYear) %>%
        mutate(fitted_prop = fitted / sum(fitted))

      otc_ctl_overlap_default_gamm4[[paste0(c(s, spp, sex), collapse = "_")]] <-
        Metrics::rmse(tst_otc_default_gamm4$ystd, tst_otc_default_gamm4$fitted_prop)


      # rmse on proportion of total that show up at each sampling point - doy
      tst_ctl_doy_gamm4 <- df2_ctl_doy_gamm4 %>%
        group_by(numYear) %>%
        mutate(fitted_prop = fitted / sum(fitted))

      ctl_ctl_overlap_doy_gamm4[[paste0(c(s, spp, sex), collapse = "_")]] <-
        Metrics::rmse(tst_ctl_doy_gamm4$ystd, tst_ctl_doy_gamm4$fitted_prop)

      tst_otc_doy_gamm4 <- df2_otc_doy_gamm4 %>%
        group_by(numYear) %>%
        mutate(fitted_prop = fitted / sum(fitted))
      otc_ctl_overlap_doy_gamm4[[paste0(c(s, spp, sex), collapse = "_")]] <-
        Metrics::rmse(tst_otc_doy_gamm4$ystd, tst_otc_doy_gamm4$fitted_prop)
    }
  }
}


### make ctl predicted by ctl model
ctl_ctl_overlap_best_gamm4 <- lapply(ctl_ctl_overlap_best_gamm4, function(x) data.frame(x)) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex") %>%
  rename(best = x)
ctl_ctl_overlap_default_gamm4 <- lapply(ctl_ctl_overlap_default_gamm4, function(x) data.frame(x)) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex") %>%
  rename(default = x)
ctl_ctl_overlap_doy_gamm4 <- lapply(ctl_ctl_overlap_doy_gamm4, function(x) data.frame(x)) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex") %>%
  rename(doy = x)

ctl_ctl_overlap_pois_gamm4 <- full_join(ctl_ctl_overlap_best_gamm4, ctl_ctl_overlap_default_gamm4) %>%
  full_join(., ctl_ctl_overlap_doy_gamm4) %>%
  separate(., sub_spp_sex,
    into = c("sub", "spp", "sex"), sep = "_",
    remove = FALSE
  ) %>%
  unite("sub_spp_sex", sub, spp, sex, remove = FALSE) %>%
  filter(sub_spp_sex %in% has_local_minima_spline_for_plotting$sub_spp_sex)

# -- PLOT EXAMPLE UD CTET FITTED 3 WAYS ----
# Makes Fig S2
# plot example spp to show method
UD_CTET <- full_season_sampling_midpoint_ctl %>%
  filter(strGSpp == "CTET" & strSexx == "B" & subsite == "UD") %>%
  mutate(numYear = factor(numYear)) %>%
  droplevels()
UD_CTET <- base::split(UD_CTET, UD_CTET$numYear)
par_best <- best_mods_spline_gamm4$best_thresh[best_mods_spline_gamm4$sub_spp_sex == "UD_CTET_B"]

# calc GDDs 2 ways + doy for comparison viz
# gdd_max
UD_CTET_gdd_max <- overlap_stats_GDD(UD_CTET,
  method = "DD_func", thresh = par_best,
  predictionmode = FALSE,
  opt = "AIC",
  s4only = FALSE,
  input_model = NULL,
  scale_effort = TRUE,
  scale_input = TRUE
)

UD_CTET_gdd <- overlap_stats_GDD(UD_CTET,
  method = "DD_func", thresh = 40,
  predictionmode = FALSE,
  opt = "AIC",
  s4only = FALSE,
  input_model = NULL,
  scale_effort = TRUE,
  scale_input = TRUE
)


UD_CTET_doy <- overlap_stats_GDD(UD_CTET, # thresh=par,
  method = "doy",
  predictionmode = FALSE,
  boundaries = NULL, opt = "AIC",
  s4only = FALSE,
  input_model = NULL,
  scale_effort = TRUE,
  scale_input = TRUE
)

max_example <- ggplot(UD_CTET_gdd_max$dd_raw %>%
  mutate(year = group), aes(
  x = x / 24, y = (y / effort) + 1, group = year,
  color = year
)) + # ,
  geom_line() +
  # ylab ('new flowers/day + 1')+
  ylab(NULL) +
  xlab(expression(GDD[max])) +
  # xlab('GDD_max') +
  scale_y_continuous(
    trans = scales::log_trans(),
    breaks = c(1, 10, 100, 1000)
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "D") +
  # viridis::scale_color_viridis(discrete = TRUE)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gdd_example <- ggplot(UD_CTET_gdd$dd_raw %>% mutate(year = group), aes(
  x = x / 24, y = (y / effort) + 1, group = year,
  color = year
)) + # ,
  # linetype=group))+
  geom_line() +
  # ylab ('new flowers/day + 1')+
  ylab(NULL) +
  xlab("GDD") +
  scale_y_continuous(
    trans = scales::log_trans(),
    breaks = c(1, 10, 100, 1000)
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "D") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

doy_example <- ggplot(UD_CTET_doy$dd_raw %>% mutate(year = group), aes(
  x = x, y = (y / effort) + 1, group = year,
  color = year
)) +
  geom_line() +
  ylab(NULL) +
  xlab("DOY") +
  scale_y_continuous(
    trans = scales::log_trans(),
    breaks = c(1, 10, 100, 1000)
  ) +
  viridis::scale_color_viridis(discrete = TRUE, option = "D") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# https://cran.r-project.org/web/packages/lemon/vignettes/legends.html#shared-legend-across-multiple-plot
UD_CTET <- cowplot::plot_grid(doy_example + theme(legend.position = "none"),
  gdd_example + theme(legend.position = "none"),
  max_example + theme(legend.position = "none"),
  ncol = 1,
  labels = c("a", "b", "c"), label_x = 0, label_y = 1,
  hjust = -4.5, vjust = 2
)

UD_CTET <- cowplot::plot_grid(UD_CTET, cowplot::get_legend(doy_example),
  ncol = 2, rel_widths = c(1, .1)
)

UD_CTET <- gridExtra::grid.arrange(gridExtra::arrangeGrob(UD_CTET,
  left = grid::textGrob("(flowers censused +1) / intercensus period", rot = 90, vjust = 1)
))


jpeg("plots/UD_CTET_example.jpeg",
  width = 8, height = 6, units = "in",
  res = 600
)
grid::grid.draw(UD_CTET)
dev.off()

# -- PLOT PREDICTED PROPORTIONS ----
ctl_vs_ctl_plot_pois_basic_gamm4 <-
  ggplot(
    data = ctl_ctl_overlap_pois_gamm4 %>% pivot_longer(.,
      cols = c("best", "default", "doy"),
      names_to = "method",
      values_to = "rmse"
    ) %>%
      unite("rep", sub, spp, sex) %>%
      mutate(
        method =
          case_when(
            method == "best" ~ "Max Temperature GDD",
            method == "default" ~ "GDD (no Max Temperature)",
            method == "doy" ~ "DOY",
            TRUE ~ "NA"
          )
      ),
    aes(x = method, y = rmse, fill = method)
  ) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  labs(y = "RMSE (per species-location)", x = NULL) +
  # Removing legends
  guides(fill = FALSE, color = FALSE) +
  scale_fill_grey() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 7)) +
  scale_x_discrete(
    labels = c(expression(DOY), expression(GDD), expression(GDD[max]))
  ) +
  ggtitle("Predicted vs actual proportion of flowers in each doy \n in ambient plots using 3 different methods")

otc_ctl_overlap_best_gamm4 <- lapply(otc_ctl_overlap_best_gamm4, function(x) data.frame(x)) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex") %>%
  rename(best = x)
otc_ctl_overlap_default_gamm4 <- lapply(otc_ctl_overlap_default_gamm4, function(x) data.frame(x)) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex") %>%
  rename(default = x)
otc_ctl_overlap_doy_gamm4 <- lapply(otc_ctl_overlap_doy_gamm4, function(x) data.frame(x)) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex") %>%
  rename(doy = x)

otc_ctl_overlap_pois_gamm4 <- full_join(otc_ctl_overlap_best_gamm4, otc_ctl_overlap_default_gamm4) %>%
  full_join(., otc_ctl_overlap_doy_gamm4) %>%
  separate(., sub_spp_sex,
    into = c("sub", "spp", "sex"), sep = "_",
    remove = FALSE
  ) %>%
  unite("sub_spp_sex", sub, spp, sex, remove = FALSE) %>%
  filter(sub_spp_sex %in% has_local_minima_spline_for_plotting$sub_spp_sex) %>%
  unite("spp_sex", spp, sex) %>%
  left_join(., taxon)

otc_vs_ctl_plot_pois_basic_gamm4 <-
  ggplot(
    data = otc_ctl_overlap_pois_gamm4 %>% pivot_longer(.,
      cols = c("best", "default", "doy"),
      names_to = "method",
      values_to = "rmse"
    ) %>%
      mutate(
        method =
          case_when(
            method == "best" ~ "Max Temperature GDD",
            method == "default" ~ "GDD (no Max Temperature)",
            method == "doy" ~ "DoY",
            TRUE ~ "NA"
          )
      ),
    aes(x = method, y = rmse, fill = method)
  ) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  labs(y = "RMSE (per species/location combination)", x = NULL) +
  # Removing legends
  guides(fill = FALSE, color = FALSE) +
  scale_fill_grey() +
  scale_x_discrete(
    labels = c(expression(DOY), expression(GDD), expression(GDD[max]))
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  theme(axis.text.x = element_text(size = 10, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Predicted vs actual proportion of flowers in each doy \n in warmed plots using 3 different methods")

# makes Fig 3
# combine above into 2 panel plot
joint_panel_pois_infloct_validation_spline_basic_gamm4 <- egg::ggarrange(
  ctl_vs_ctl_plot_pois_basic_gamm4 +
    ggtitle("ambient") +
    ylim(0.00, 0.27) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(size = 10, face = "bold", color = "black"),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 10)
    ),
  otc_vs_ctl_plot_pois_basic_gamm4 +
    ggtitle("warmed") +
    ylim(0.00, 0.27) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 10, face = "bold", color = "black"),
      axis.text.y = element_text(size = 9)
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ),
  nrow = 1, labels = c("a", "b"), padding = 1,
  label.args = list(gp = grid::gpar(
    cex =
      1.2, fontsize = 7, fontface = "bold", x = unit(1, "line")
  ), hjust = -0.5, vjust = 1.5)
)

ggsave(joint_panel_pois_infloct_validation_spline_basic_gamm4,
  file = "plots/joint_panel_pois_infloct_validation_spline_basic_gamm4.jpg",
  width = 90, height = 90, units = "mm"
)

# -- PLOT SPP*SITE THRESHOLDS AND CIs ####
jack_CI_spline_gamm4 <- lapply(out_AIC_jack_spline_gamm4, function(x) {
  if (all(is.na(x))) {
    NULL
  }
  else {
    data.frame(
      thresh = unlist(x$levels), y = x$values,
      min_AIC = x$minfun,
      best_thresh = x$minlevels
    )
  }
}) %>%
  data.table::rbindlist(., idcol = "sub_spp_sex_yr") %>%
  separate(., sub_spp_sex_yr,
    into = c("sub", "spp", "sex", "yr"),
    remove = FALSE
  ) %>%
  select(sub, spp, sex, best_thresh, yr) %>%
  distinct() %>%
  unite(., "spp_sex", spp, sex) %>%
  group_by(sub, spp_sex) %>%
  summarize(
    maxT_median = median(best_thresh),
    maxT_mean = mean(best_thresh),
    q95 = quantile(best_thresh, 0.95),
    q05 = quantile(best_thresh, 0.05),
    min = min(best_thresh),
    max = max(best_thresh),
    ct = dplyr::n()
  ) %>%
  unite("sub_spp_sex", sub, spp_sex, remove = FALSE) %>%
  left_join(., data.frame(sub_spp_sex = has_local_minima_spline_for_plotting$sub_spp_sex, has_local_minima_spline = 1)) %>%
  left_join(., taxon) %>%
  mutate(
    site = ifelse(sub %in% c("UW", "UD"), "Utqiagvik", "Atqasuk"),
    commtype = ifelse(sub %in% c("AD", "UD"), "Dry Heath", "Wet Meadow")
  ) %>%
  mutate(relationship_to_max_temperature = ifelse(is.na(has_local_minima_spline), "non-saturating", "saturating"))

# change factor levels for plotting
df <- jack_CI_spline_gamm4 %>% mutate(
  commtype = ifelse(
    commtype == "Wet Meadow", "Wet", "Dry"
  ),
  site = factor(site)
)
levels(df$site) <- c("Atqasuk", expression(paste("Utqia", dot(g), "vik")))

# shorten for vertical fit
df <- df %>%
  mutate(sci_name = ifelse(sci_name == "Eriophorum angustifolium ssp. triste",
    "E. angustifolium ssp. triste",
    sci_name
  ))

# add ref line to illustrate
summer_temps <- tibble(
  commtype = c("Dry", "Wet", "Dry", "Wet"),
  site = factor(c("Atqasuk", "Atqasuk", 'paste("Utqia", dot(g), "vik")', 'paste("Utqia", dot(g), "vik")'),
    levels = (c("Atqasuk", expression(paste("Utqia", dot(g), "vik"))))
  ),
  sum_temp = c(
    quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "AD" & temp_summry_stats$strTrea == "C"], 0.95),
    quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "AW" & temp_summry_stats$strTrea == "C"], 0.95),
    quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "UD" & temp_summry_stats$strTrea == "C"], 0.95),
    quantile(temp_summry_stats$AvgofnumTemp[temp_summry_stats$strSitCom == "UW" & temp_summry_stats$strTrea == "C"], 0.95)
  )
)

# Makes Fig 2
# plot thresholds by spp and subsite
sumryplot_spline_gamm4 <- ggplot(df, aes(
  x = sci_name, y = maxT_median, ymin = min, ymax = max, # group=sub,
  fill = relationship_to_max_temperature,
  color = relationship_to_max_temperature
)) +
  geom_hline(
    data = summer_temps,
    aes(yintercept = sum_temp), alpha = 0.3,
    linetype = "solid", size = 0.5
  ) +
  geom_col(position = position_dodge()) +
  geom_errorbar(position = position_dodge(), color = "black", width = 0.5) +
  facet_grid(site ~ commtype,
    scales = "free_x",
    labeller = label_parsed
  ) +
  scale_fill_manual(values = c("#D3D3D3", "#4D4D4D")) +
  scale_color_manual(values = c("#D3D3D3", "#4D4D4D")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(strip.text.y = element_text(size = 10)) +
  # https://stackoverflow.com/questions/37554118/ggplot-inserting-space-before-degree-symbol-on-axis-label
  labs(y = expression("Temperature threshold (" * degree * C * ")"), x = NULL) +
  guides(
    fill = guide_legend("relationship to \n max temp "),
    color = guide_legend("relationship to \n max temp ")
  ) +
  theme(axis.text.x = element_text(face = "italic", size = 7)) +
  theme(legend.text = element_text(size = 8))

ggsave(sumryplot_spline_gamm4,
  file = paste0("plots/sumryplot_inflo_maxT_midpoint_spline_gamm4.jpg"),
  width = 180, height = 90, units = "mm"
)

# -- SUMMARY TABLES AND STATS ####
# make table of plotyrs
sumry_table_plotyrs <- full_season_sampling_midpoint_ctl %>%
  filter(!is.na(ct) & ct > 0 & numYear != 1998 & numYear != 2019) %>%
  select(subsite, strGSppSex, numYear) %>%
  distinct() %>%
  group_by(subsite, strGSppSex) %>%
  summarize(nyr = dplyr::n()) %>%
  left_join(taxon %>%
    mutate(strGSppSex = gsub("_", "", spp_sex))) %>%
  filter(!is.na(sci_name)) %>% # remove spp things that are only in OTC
  select(-strGSppSex, -spp_sex) %>%
  pivot_wider(names_from = subsite, values_from = nyr) %>%
  arrange(sci_name)

write.csv(sumry_table_plotyrs, "tables/sumry_table_plotyrs.csv", row.names = FALSE, na = "")

# summary table of maxct by year
sumry_table_maxct <- full_season_sampling_midpoint_ctl %>%
  group_by(numYear, strGSppSex, subsite) %>%
  summarize(max_ct = max(ct, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(., names_from = numYear, values_from = max_ct) %>%
  left_join(taxon %>%
    mutate(strGSppSex = gsub("_", "", spp_sex))) %>%
  arrange(sci_name, subsite) %>%
  select(-spp_sex) %>%
  select(strGSppSex, sci_name, subsite, everything())

write.csv(sumry_table_maxct, "tables/max_ct_by_spp_by_year.csv",
  row.names = FALSE
)

# -- PLOT GDD accumulation typical ----
# Makes Fig 1
# temp plots
# figure out avg doy with non0 flcts
sumry <- full_season_sampling_midpoint_ctl %>%
  filter(!is.na(ct) & ct > 0 & numYear != 1998 & numYear != 2019) %>%
  # just grab paired OTC data
  bind_rows(full_season_sampling_midpoint_otc %>%
    inner_join(
      .,
      full_season_sampling_midpoint_ctl %>%
        select(subsite, strGSppSex) %>%
        distinct()
    )) %>%
  filter(!is.na(ct) & ct > 0 & numYear != 1998 & numYear != 2019) %>%
  group_by(subsite, strTrea, strGSppSex) %>%
  summarize(numJulian = mean(weighted.mean(numJulian, weight = ct))) %>%
  left_join(taxon %>%
    mutate(strGSppSex = gsub("_", "", spp_sex))) %>%
  filter(!is.na(sci_name)) # remove spp things that are only in OTC

# filter down to temps actually used
temps <- full_season_sampling_midpoint_ctl %>%
  bind_rows(full_season_sampling_midpoint_otc) %>%
  filter(numYear != 1998 & numYear != 2019) %>%
  select(-strGSpp, -strSexx, -strGSppSex, -ct) %>%
  distinct()
# forcings illustration - calc cc for all values 0-20

DD_by_day <- function(df, n) {
  varname <- paste("DD", n, sep = ".")
  mutate(df %>% rowwise(), !!varname := DD_func(temp = meanT, thresh = n, sf = sf))
}

for (b in (c(0, 2, 4, 6, 10, 20))) {
  temps <- full_join(
    temps,
    DD_by_day(temps, b)
  )
}
# cumsum
DD_cum_sum <- temps %>%
  group_by(subsite, strTrea, numYear) %>%
  arrange(numJulian, numHr) %>%
  mutate_at(vars(contains("DD")), .funs = list(cat = ~ cumsum(.))) %>%
  arrange(subsite, strTrea, numYear, numJulian, numHr)

# average over yrs
DD_cum_sum <- DD_cum_sum %>%
  ungroup() %>%
  group_by(subsite, strTrea, numJulian) %>%
  summarise_at(vars(contains("_cat")), .funs = list(
    ave = ~ mean(.),
    mn = ~ min(.),
    mx = ~ max(.)
  )) %>%
  pivot_longer(., -c("subsite", "strTrea", "numJulian"),
    names_to = "threshold", values_to = "DD_sum"
  ) %>%
  mutate(
    how =
      case_when(
        grepl("ave", threshold) ~ "ave",
        grepl("mn", threshold) ~ "mn",
        grepl("mx", threshold) ~ "mx",
        TRUE ~ "NA"
      )
  ) %>%
  mutate(threshold = gsub("[^0-9]+", "", threshold)) %>%
  mutate(threshold = as.numeric(threshold)) %>%
  mutate(
    strTrea = ifelse(strTrea == "OTC", "warmed", "ambient"),
    subsite = factor(subsite)
  )

levels(DD_cum_sum$subsite) <- c(
  "Atqasuk - Dry",
  "Atqasuk - Wet",
  expression(paste("Utqia", dot(g), "vik - Dry")),
  expression(paste("Utqia", dot(g), "vik - Wet"))
)

sci_labs <- sumry %>%
  mutate(sci_name = ifelse(sci_name == "Eriophorum angustifolium ssp. triste",
    "E. angustifolium ssp. triste", sci_name
  )) %>%
  mutate(
    strTrea = ifelse(strTrea == "OTC", "warmed", "ambient"),
    subsite = factor(subsite)
  )
levels(sci_labs$subsite) <- c(
  "Atqasuk - Dry",
  "Atqasuk - Wet",
  expression(paste("Utqia", dot(g), "vik - Dry")),
  expression(paste("Utqia", dot(g), "vik - Wet"))
)

DD_sum_plot_ave <- ggplot(DD_cum_sum %>% filter(how == "ave" & numJulian <= 235), aes(x = numJulian, y = DD_sum / 24)) +
  geom_line(aes(group = threshold, color = threshold)) +
  facet_grid(strTrea ~ subsite, labeller = label_parsed) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_point(data = sci_labs, aes(x = numJulian, y = 0)) +
  geom_text_repel(
    data = sci_labs, aes(x = numJulian, y = 0, label = sci_name),
    nudge_y = -25,
    direction = "x",
    angle = -90,
    vjust = 1,
    segment.size = 0.1,
    segment.color = "grey0",
    size = 2
  ) +
  ylim(-400, 800) +
  theme_bw() +
  xlim(122, 255) +
  geom_text(
    data = DD_cum_sum %>% filter(how == "ave" & numJulian == 235),
    aes(
      label = threshold, colour = threshold,
      x = numJulian + 10, y = DD_sum / 24
    ), size = 3
  ) +
  theme(
    plot.margin = unit(c(1, 3, 1, 1), "lines"),
    text = element_text(size = 9),
    strip.text = element_text(size = 9)
  ) +
  guides(color = FALSE) + # ggtitle("Growing degree accumulation in an average year")+
  labs(
    x = "day of year (DOY)",
    y = expression("GDD"[max] * " \n (sum of hourly temperatures above 0° C \n and capped at the maximum temperature threshold)")
  )


DD_sum_plot_ave <- DD_sum_plot_ave +
  ylab("") + # empty label
  # Tweak the margins (push the label down by forcing a wider top margin)
  theme(axis.title.y = element_text(
    size = 10, # also adjust text size if needed
    margin = margin(
      t = 0, r = 3, b = 0, l = 0,
      unit = "mm"
    )
  ))

# multiline y axis
line_1 <- expression("GDD"[max])
line_2 <- "(sum of hourly temperatures above 0 °C"
line_3 <- "and capped at the maximum temperature threshold)/24"

# Call cowplot::draw_label two times to plot two lines of text
DD_sum_plot_ave <- ggdraw(DD_sum_plot_ave) +
  draw_label(line_1, x = 0.015, y = 0.5, angle = 90, size = 9) + # use relative coordinates for positioning
  draw_label(line_2, x = 0.030, y = 0.5, angle = 90, size = 9) +
  draw_label(line_3, x = 0.045, y = 0.5, angle = 90, size = 9)

ggsave(DD_sum_plot_ave,
  file = paste0("plots/DD_sum_plot_ave.jpg"),
  width = 180, height = 180, units = "mm", dpi = 600
)

# -- PLOT EXAMPLE TEMP DIFS OTC VS CTL MIN/MAX/MEAN ----

difs <- infilled_hourly %>%
  filter(numYear == 2008) %>% # pick a year to show
  group_by(strSitCom, strTrea, date) %>%
  summarize(
    meanT = mean(AvgofnumTemp), maxT = max(AvgofnumTemp),
    minT = min(AvgofnumTemp), .groups = "drop"
  ) %>%
  pivot_wider(.,
    names_from = "strTrea",
    values_from = c("meanT", "maxT", "minT")
  ) %>%
  mutate(
    dif_max = maxT_E - maxT_C,
    dif_min = minT_E - minT_C,
    dif_mean = meanT_E - meanT_C
  ) %>%
  select(strSitCom, date, dif_max, dif_min, dif_mean) %>%
  pivot_longer(cols = c(dif_max, dif_min, dif_mean), names_to = "param") %>%
  group_by(strSitCom, param) %>%
  arrange(date) %>%
  mutate(param_roll = slider::slide_dbl(value, mean, .before = 3, .after = 3)) %>%
  mutate(
    site = ifelse(strSitCom %in% c("UW", "UD"), "Utqiagvik", "Atqasuk"),
    commtype = ifelse(strSitCom %in% c("AD", "UD"), "Dry", "Wet")
  ) %>%
  mutate(param = case_when(
    param == "dif_max" ~ "Max Temperature",
    param == "dif_min" ~ "Min Temperature",
    param == "dif_mean" ~ "Mean Temperature",
    TRUE ~ "NA"
  )) %>%
  filter(date >= "2008-07-01" & date < "2008-08-01")

# for nice plots
difs$site <- factor(difs$site)
levels(difs$site) <- c("Atqasuk", expression(paste("Utqia", dot(g), "vik")))

use_theme <- function() {
  theme_bw() +
    theme(
      axis.text.y = element_text(size = 4),
      axis.text.x = element_text(angle = -90),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
}

# Makes Fig S1
# plot
temp_plot <- ggplot(
  difs,
  aes(x = date, y = value, color = param, group = param)
) +
  geom_point(alpha = 0.2) +
  geom_line(aes(x = date, y = param_roll)) +
  facet_grid(site ~ commtype, labeller = label_parsed) +
  labs(y = "Temperature difference (°C) \n(warmed-ambient)", x = NULL) +
  use_theme() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "#4D4D4D", "blue")) +
  theme(strip.text.x = element_text(size = 10)) +
  theme(strip.text.y = element_text(size = 10)) +
  scale_x_date(
    labels = date_format("%m-%d"),
    limits = c(as.Date("2008-07-01"), as.Date("2008-07-31")),
    expand = c(0, 0)
  ) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10))

ggsave(temp_plot, file = "plots/temp_plot.jpg", width = 180, height = 90, units = "mm")
