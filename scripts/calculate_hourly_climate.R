# script to infill hourly temperature
# based on Keith Jennings script here:
# https://github.com/NWTlter/NWT_climate_infilling/blob/master/hrly_met/r_code/nwt_hrly_qc_infill_TEMP_V2.R

# Temporal, single station infilling (gap <= 72 h) based on Liston and Elder (2006) — aka, MicroMet
# e.g. <1 hour interpolate, 2-24 - mean of day before and after
# otherwise skip

# -- read and process temperature data from repo ----
dataTpHrcanopyTr <- read.csv("data/BARROW_ATQASUK/dataTpHrcanopyTr.csv")

# remove missing value code, rename y col to AvgofnumTemp so rest of code works
dataTpHrcanopyTr <- dataTpHrcanopyTr %>%
  mutate(date = as.Date(numJulian, origin = paste0(numYear, "-01-01")) - 1) %>%
  mutate(AvgofnumTemp = ifelse(numTp == -6999, NA, numTp))
dataTpHrcanopyTr <- dataTpHrcanopyTr %>% arrange(date, numHr, strSitCom, strTrea)

# raw data - all have 24 hours per day
range(dataTpHrcanopyTr %>%
  group_by(date, strSitCom, strTrea, numYear, numJulian) %>%
  summarize(ct = dplyr::n()) %>%
  pull(ct))

# find gaps
gaps <- dataTpHrcanopyTr %>%
  group_by(strSitCom, strTrea, date) %>%
  filter(!is.na(AvgofnumTemp)) %>%
  summarize(ct = dplyr::n()) %>%
  filter(ct < 24) %>%
  left_join(., dataTpHrcanopyTr) %>%
  arrange(strSitCom, strTrea, date, numHr)

dataTpHrcanopyTr$gapID <- NA

# number gaps by ID
dataTpHrcanopyTr <- dataTpHrcanopyTr %>%
  arrange(strSitCom, strTrea, date, numHr)
counter <- 1
for (i in 1:nrow(dataTpHrcanopyTr)) {
  if (is.na(dataTpHrcanopyTr$AvgofnumTemp[i])) {
    dataTpHrcanopyTr$gapID[i] <- counter
    if (!is.na(dataTpHrcanopyTr$AvgofnumTemp[i + 1]) | i == nrow(dataTpHrcanopyTr)) {
      counter <- counter + 1
    }
  }
}

# make gaps unique to subsite and year just in case any are at the boundaries
dataTpHrcanopyTr <- dataTpHrcanopyTr %>%
  mutate(gapID = paste0(gapID, strSitCom, strTrea, numYear))

range(dataTpHrcanopyTr$AvgofnumTemp, na.rm = TRUE)

fillType <- dataTpHrcanopyTr %>%
  filter(!is.na(gapID)) %>%
  group_by(gapID) %>%
  summarize(ct = dplyr::n()) %>%
  ungroup() %>%
  mutate(
    filltype =
      case_when(
        ct == 1 ~ "INTERP",
        ct >= 2 & ct <= 24 ~ "AVG",
        ct > 24 ~ "ARIMA"
      )
  )

dataTpHrcanopyTr <- dataTpHrcanopyTr %>%
  left_join(., fillType)

# cts are all <24 hours or very multi day, in which case arima not useful
# so just use the INTERP and AVG methods:
for (i in 1:nrow(dataTpHrcanopyTr)) {
  if (is.na(dataTpHrcanopyTr$AvgofnumTemp[i])) {
    if (dataTpHrcanopyTr$filltype[i] == "INTERP") {
      dataTpHrcanopyTr$AvgofnumTemp[i] <-
        mean(dataTpHrcanopyTr$AvgofnumTemp[i - 1], dataTpHrcanopyTr$AvgofnumTemp[i + 1])
    }
    if (dataTpHrcanopyTr$filltype[i] == "AVG" &
      dataTpHrcanopyTr$strSitCom[i - 24] ==
        dataTpHrcanopyTr$strSitCom[i + 24]) {
      dataTpHrcanopyTr$AvgofnumTemp[i] <-
        mean(dataTpHrcanopyTr$AvgofnumTemp[i - 24], dataTpHrcanopyTr$AvgofnumTemp[i + 24],
          na.rm = TRUE
        )
    }
  }
}

# rename to infilled for the subsequent code chunk
infilled_hourly <- dataTpHrcanopyTr

# cleanup
rm(dataTpHrcanopyTr)
