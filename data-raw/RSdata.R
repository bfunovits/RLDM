## code to prepare `RSdata` dataset goes here

pkgs <- c("tidyverse", "lubridate")
void <- lapply(pkgs, library, character.only = TRUE)

# 1) Get Ramey data ####

# Warnings because NAs are coded as dots (.).
# They are imported as NAs, so that is fine
data <- read_csv("inst/extdata/govdat3908.csv",
                 col_types = cols(nwbus = col_double(),
                                  pbus = col_double(), rcndsv = col_double(),
                                  romerexog = col_double(), rsl = col_double(),
                                  tothoursces = col_double(), wara = col_double())) %>%
  filter(quarter>=1947) # like in Ramey, important that coefficients for intercept and trends coincide

# Adding linear and quadratic trend
data = data %>%
  mutate(trend_lin = 1:nrow(data),
         trend_quadr = (1:nrow(.))^2)

# 2) All Data for VAR ####

data_var1 = data %>%
  mutate(rwbus = nwbus/pbus) %>%
  select(quarter, rgov, rgdp, rcndsv, rinvfx, tothoursces, rwbus) %>%
  mutate_all(list(log = log)) %>%
  select(quarter, ends_with("_log")) %>%
  select(-quarter_log)

data_var = data_var1 %>%
  left_join(data %>% select(quarter, amtbr, starts_with("trend_")), by = "quarter")

# date column of correct type
data_var = data_var %>%
  mutate(date = (quarter %% 1)*4 + 1) %>%
  mutate(quarter = as.integer(quarter)) %>%
  mutate(date = paste0(quarter, " Q", date)) %>%
  mutate(date = yq(date)) %>%
  select(-quarter)

# 3) Separate in exogonous and endogenous ####

# Transform tibble (data.frame) to matrix and replace NAs
data_mat = data_var %>% select(-date, -starts_with("trend_")) %>% as.matrix()
data_mat[is.na(data_mat)] = 0
col_names = dimnames(data_mat)[[2]]

# Matrix of exogenous variables: Linear and quadratic trend
data_exog = data_var %>% select(starts_with("trend_")) %>% as.matrix()
data_exog[is.na(data_exog)] = 0

# 4) Orthogonalize the data w.r.t. the exogenous variables ####
B_ex = solve(t(data_exog) %*% data_exog) %*% t(data_exog) %*% data_mat
RSdata = data_mat - data_exog %*% B_ex

# 5) Transform to tibble, ts, xts
RSdata_ts = ts(RSdata,
               frequency = 4, start = c(1947, 1))

RSdata = bind_cols(data_var %>% select(date),
                   RSdata %>% as_tibble())

RSdata_xts = xts::xts(RSdata %>% select(-date) %>% as.matrix(),
                  order.by = RSdata$date)

# 5) Save data ####
usethis::use_data(RSdata, RSdata_ts, RSdata_xts, overwrite = TRUE)
