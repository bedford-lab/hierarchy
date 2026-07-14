# Libraries ----
library(tidyverse) 

# Paths ----
setwd("~/Dropbox (Personal)/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data/")

# ==================================================
# Load ----
# ==================================================

md <- read.csv("SDH_Cohorts - Mice.csv", stringsAsFactors = TRUE)
bw <- read.csv("predicted_weights.csv",  stringsAsFactors = TRUE)
tc <- read.csv("testes_clean.csv",       stringsAsFactors = TRUE)

# ==================================================
# Standardize dates ----
# ==================================================

to_date_mdy2 <- function(x) as.Date(as.character(x), format = "%m/%d/%y")

md <- md %>% mutate(DOB = to_date_mdy2(DOB), Sac_date = to_date_mdy2(Sac_date))
bw <- bw %>% mutate(DOB = to_date_mdy2(DOB))
tc <- tc %>% mutate(DOB = to_date_mdy2(DOB), Sac_date = to_date_mdy2(Sac_date))

# ==================================================
# Join metadata + weights + testes ----
# ==================================================

common_md_bw <- Reduce(intersect, list(names(md), names(bw)))

md_bw <- md %>%
  full_join(bw %>% select(-c(2:10)), by = common_md_bw)  # drop X28–X140

common_md_bw_tc <- Reduce(intersect, list(names(md_bw), names(tc)))

md_bw_tc <- md_bw %>% full_join(tc, by = common_md_bw_tc) %>% filter(!is.na(Experimenter))

# ==================================================
# Export ----
# ==================================================
write.csv(md_bw_tc, "meta_data.csv", row.names = FALSE)

# ---- Weaning statistics ----
df <- md_bw_tc %>%
  select(Cohort, MouseID, Treatment, Wean_age, Wean_weight, Pred_Wean_Weight, Wean_date) %>%
  mutate(Wean_date = as.Date(as.character(Wean_date), format = "%m/%d/%y"))

mean(df$Wean_age, na.rm = TRUE)
median(df$Wean_age, na.rm = TRUE)

mean_pairwise_diff <- function(x) {
  x <- na.omit(x)
  if (length(x) < 2) return(NA_real_)
  mean(abs(dist(x)))
}

cohort_summary <- df %>%
  group_by(Cohort, Treatment) %>%
  summarise(
    n_mice = dplyr::n(),
    same_wean_date = n_distinct(Wean_date, na.rm = TRUE) == 1,
    # Mean pairwise age difference
    mean_age_diff = mean_pairwise_diff(Wean_age),
    # Mean pairwise WEAN weight difference
    mean_wean_weight_diff = mean_pairwise_diff(Wean_weight),
    # Mean pairwise PREDICTED weight difference
    mean_pred_wean_weight_diff = mean_pairwise_diff(Pred_Wean_Weight),
    .groups = "drop"
  )

plot(cohort_summary$mean_wean_weight_diff, cohort_summary$mean_pred_wean_weight_diff)

col_means <- cohort_summary %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
