## code to prepare `BQdata` dataset goes here

# Data taken from https://academic.oup.com/restud/advance-article/doi/10.1093/restud/rdz028/5490841#supplementary-data

# Parameters ####
MAKE_PLOTS = FALSE

# Loading packages and raw data ####
pkgs <- c("tidyverse", "xts")
void <- lapply(pkgs, library, character.only = TRUE)


tbq = read_delim("inst/extdata/bqdata.csv",
                 delim = ";",
                 escape_double = FALSE,
                 col_types = cols(DATE = col_date(format = "%Y-%m-%d")),
                 trim_ws = TRUE,)

tbq %>% names()
tbq = tbq %>%
  rename(date = DATE,
         unemp = LHMUR)

if(MAKE_PLOTS){
  # Plot as xts
  tbq_xts = xts(tbq %>% select(-date) %>% as.matrix(), order.by = tbq$date)
  plot(tbq_xts)

  # Plot with ggplot2
  tbq_gg = tbq %>%
    gather(key = "type", value = "value", -date)

  p <- tbq_gg %>%
    ggplot(aes(date, value, color = type)) + geom_line()

  plotly::ggplotly(p)

}

# Transform to real GDP, then obtain growth rate ####
#
tbq = tbq %>%
  mutate(realGDP = GNP/GD87) %>%
  select(-GNP, -GD87) %>%
  mutate(lag_realGDP = dplyr::lag(realGDP)) %>%
  mutate(rGDPgrowth = (realGDP - lag_realGDP)/lag_realGDP) %>%
  select(-lag_realGDP)

if(MAKE_PLOTS){
  # Plot again
  tbq_gg = tbq %>%
    gather(key = "type", value = "value", -date)

  p <- tbq_gg %>%
    ggplot(aes(date, value, color = type)) + geom_line()

  plotly::ggplotly(p)
}

# Separate into two periods, demean real GDP growth for different periods ####

tbq = tbq %>%
  mutate(regime = c(rep(1, 104), rep(2, nrow(tbq)-104))) %>%
  group_by(regime) %>%
  select(-realGDP) %>%
  mutate(rGDPgrowth = 100*rGDPgrowth) %>%
  mutate(mean_rGDPgrowth = mean(rGDPgrowth, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(trend_lin = 1:nrow(tbq))


# Detrend unemployment ####
lm_unemp = lm(unemp ~ trend_lin, data = tbq)
#lm_unemp %>% str(1)

tbq %>%
  mutate(rGDPgrowth_demeaned = rGDPgrowth - mean_rGDPgrowth) %>%
  mutate(trend_unemp = lm_unemp$fitted.values) %>%
  mutate(unemp_detrended = unemp - lm_unemp$fitted.values) -> tbq

if(MAKE_PLOTS){
  p = tbq %>%
    select(date, rGDPgrowth, mean_rGDPgrowth, unemp, trend_unemp) %>%
    filter(!is.na(rGDPgrowth)) %>%
    gather(key = "type", value = "value", -date) %>%
    ggplot(aes(date, value, color = type)) + geom_line()
  plotly::ggplotly(p)

  p = tbq %>%
    select(date, rGDPgrowth_demeaned, unemp_detrended) %>%
    filter(!is.na(rGDPgrowth_demeaned)) %>%
    gather(key = "type", value = "value", -date) %>%
    ggplot(aes(date, value, color = type)) + geom_line()
  plotly::ggplotly(p)
}

# Save data  ####
BQdata = tbq %>%
  select(date, rGDPgrowth_demeaned, unemp_detrended) %>%
  filter(!is.na(rGDPgrowth_demeaned))

BQdata_ts = ts(BQdata %>% select(-date) %>% as.matrix(),
               frequency = 4, start = c(1948, 2))

BQdata_xts = xts::xts(BQdata %>% select(-date) %>% as.matrix(),
                      order.by = BQdata$date)

usethis::use_data(BQdata, BQdata_ts, BQdata_xts, overwrite = TRUE)
