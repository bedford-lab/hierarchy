# ---- Libraries ----
suppressPackageStartupMessages({
library(Rmisc); library(data.table); library(tidyverse); library(lme4); library(lmerTest);
library(emmeans); library(EloRating); library(igraph); library(ggraph); 
library(stringr); library(purrr); library(broom); library(patchwork); library(car)
})

# ---- Paths ----
setwd("~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data")
plotdir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/warm_spot"

theme_nolabels <- theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        strip.text = element_blank(), legend.position = "none")

# ==================================================
# LUT + Event logs
# ==================================================
lut <- read.csv("SDH_WSA - LUT.csv", stringsAsFactors = TRUE)

lut %>%
  summarise(n_files = n_distinct(file_name),
            n_cohorts = n_distinct(Cohort),
            n_mice = n_distinct(MouseID))

lut <- lut %>%
  mutate(
    Age        = case_when(str_detect(Treatment, "10W$") ~ "10W",
                           str_detect(Treatment, "20W$") ~ "20W",
                           TRUE ~ NA_character_),
    Geno       = if_else(str_detect(Treatment, "MUT"), "MUT", "WT"),
    Background = if_else(str_detect(Treatment, "BalbC"), "BalbC", "C57"),
    Group      = paste(Background, Geno, sep = "_"),
    Age        = factor(Age, levels = c("10W", "20W")),
    Geno       = factor(Geno, levels = c("WT", "MUT")),
    Background = factor(Background, levels = c("BalbC", "C57")),
    Group      = factor(Group, levels = c("BalbC_WT", "C57_WT", "C57_MUT"))
  )

lut %>% group_by(Group) %>% summarise(n_cohorts = n_distinct(Cohort), .groups = "drop")

singleTrial_cohorts <- lut %>%
  distinct(Cohort, Trial) %>%
  dplyr::count(Cohort, name = "n_trials") %>%
  filter(n_trials == 1) %>%
  pull(Cohort) %>%
  droplevels()

singleTrial_cohorts

# ==================================================
# Read + basic checks
# ==================================================
list_of_files <- list.files("WSA_eventLogs", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

df <- list_of_files %>%
  setNames(list_of_files) %>%
  map_df(~ readr::read_csv(.x, col_types = cols(), col_names = TRUE), .id = "file_name") %>%
  drop_na(Behavior, Subject) %>%
  mutate(file_name = str_remove(file_name, "^WSA_eventLogs/")) %>%
  droplevels()

trial_dur <- df %>%
  group_by(file_name) %>%
  summarise(TrialDur = max(Time_Relative_sf, na.rm = TRUE), .groups = "drop")

hist(trial_dur$TrialDur, breaks = 100)

# ==================================================
# Merge split trials (P1/P2)
# ==================================================
two_parts_base <- lut %>%
  filter(Part %in% c("Part1", "Part2")) %>%
  mutate(base_name = str_remove(file_name, "_P[12]( - Event Logs)?\\.csv$")) %>%
  distinct(base_name)

df2 <- df %>%
  mutate(base_name = str_remove(file_name, "_P[12]( - Event Logs)?\\.csv$"),
         part = case_when(str_detect(file_name, "P1") ~ "P1",
                          str_detect(file_name, "P2") ~ "P2",
                          TRUE ~ NA_character_))

merged_parts <- df2 %>%
  semi_join(two_parts_base, by = "base_name") %>%
  drop_na(part) %>%
  group_by(base_name) %>%
  mutate(offset = max(Time_Relative_sf[part == "P1"], na.rm = TRUE) + 5,
         Time_Relative_sf = if_else(part == "P2", Time_Relative_sf + offset, Time_Relative_sf),
         file_name = str_c(base_name, " - Event Logs.csv")) %>%
  ungroup() %>%
  select(-base_name, -part, -offset)

df <- df2 %>% filter(is.na(part)) %>% select(-base_name, -part) %>% bind_rows(merged_parts)

n_distinct(df$file_name)  # 144

# ==================================================
# LUT: normalize filenames + trial-level metadata
# ==================================================
lut_merged <- lut %>%
  mutate(file_name = case_when(
    Part %in% c("Part1", "Part2") & str_detect(file_name, "_P[12] - Event Logs\\.csv$") ~
      sub("_P[12] - Event Logs\\.csv$", " - Event Logs.csv", file_name),
    TRUE ~ file_name
  ))

lut_trials <- lut_merged %>%
  group_by(MouseID, TailMark, Trial, Date, Observer, file_name) %>%
  summarise(Cohort = first(Cohort),
            Age = first(Age), Group = first(Group), Sex = first(Sex),
            Notes = paste(unique(na.omit(Notes)), collapse = "; "),
            .groups = "drop")

n_distinct(lut_trials$file_name)  # 144

double_trials <- lut_trials %>%
  filter(Notes == "DoubleObserved") %>%
  group_by(Trial, Date, Cohort, Age, Group, Sex) %>%
  summarise(n_obs = n_distinct(Observer), .groups = "drop") %>%
  filter(n_obs > 1)

# ==================================================
# Timestamp template + state filling
# ==================================================
df <- df %>% mutate(Time_Relative_sf = round(Time_Relative_sf, 1))
timestamp_tbl <- tibble(Time_Relative_sf = seq(0, 20 * 60, by = 0.1))

TrialAnalysis <- function(trial_name) {
  trial_df <- filter(df, file_name == trial_name)
  subs <- unique(trial_df$Subject)
  
  StateFiller <- function(subj) {
    f <- trial_df %>% filter(Subject == subj) %>% select(Time_Relative_sf, Subject, Behavior, Event_Type)
    
    timestamp_tbl %>%
      left_join(f, by = "Time_Relative_sf") %>%
      mutate(State = if_else(is.na(Behavior) & is.na(Event_Type), NA_character_,
                             paste(Behavior, Event_Type, sep = "_"))) %>%
      fill(State, .direction = "down") %>%
      fill(Subject, .direction = "downup") %>%
      mutate(
        Plate = if_else(State %in% c("2 paws_State start", "4 paws_State start"), "ON", "OFF"),
        Qualifier = case_when(State == "2 paws_State start" ~ "Partial",
                              State == "4 paws_State start" ~ "Full",
                              TRUE ~ NA_character_),
        rleid = data.table::rleid(Plate),
        runLength = ave(seq_along(Plate), rleid, FUN = length)
      ) %>%
      filter(runLength != 1) %>%  # drop 0.1 s blips
      mutate(file_name = trial_name) %>%
      select(Time_Relative_sf, Subject, Plate, Qualifier, file_name)
  }
  
  map_dfr(subs, StateFiller)
}

trials <- unique(df$file_name)
result <- map_dfr(trials, TrialAnalysis)

dups <- result %>%
  summarise(n = n(), .by = c(file_name, Time_Relative_sf, Subject)) %>%
  filter(n > 1L)

dups2 <- result %>%
  semi_join(dups, by = c("file_name", "Time_Relative_sf", "Subject")) %>%
  arrange(file_name, Subject, Time_Relative_sf)

# ==================================================
# Wide format + group state
# ==================================================
# Rule: if OFF and ON exist at same timestamp, treat as ON (new state wins)
result_dt <- as.data.table(result); result_dt[, is_on := Plate == "ON"]
setorder(result_dt, file_name, Time_Relative_sf, Subject, is_on)

result_clean <- result_dt[, .SD[.N], by = .(file_name, Time_Relative_sf, Subject)][, is_on := NULL] %>%
  as_tibble()

sf_wide <- result_clean %>%
  pivot_wider(id_cols = c(file_name, Time_Relative_sf),
              names_from = Subject,
              values_from = c(Plate, Qualifier)) %>%
  mutate(groupState = case_when(
    Plate_Base == "OFF" & Plate_Middle == "OFF" & Plate_Tip == "OFF" ~ "All_OFF",
    Plate_Base == "ON"  & Plate_Middle == "OFF" & Plate_Tip == "OFF" ~ "Base_Only",
    Plate_Base == "OFF" & Plate_Middle == "ON"  & Plate_Tip == "OFF" ~ "Middle_Only",
    Plate_Base == "OFF" & Plate_Middle == "OFF" & Plate_Tip == "ON"  ~ "Tip_Only",
    Plate_Base == "ON"  & Plate_Middle == "ON"  & Plate_Tip == "OFF" ~ "BaseMiddle",
    Plate_Base == "ON"  & Plate_Middle == "OFF" & Plate_Tip == "ON"  ~ "BaseTip",
    Plate_Base == "OFF" & Plate_Middle == "ON"  & Plate_Tip == "ON"  ~ "MiddleTip",
    Plate_Base == "ON"  & Plate_Middle == "ON"  & Plate_Tip == "ON"  ~ "All_ON",
    TRUE ~ NA_character_
  ))

# ==================================================
# Cumulative warm spot duration (ranked trajectories ± SE)
# ==================================================
trial_len <- 1200
bin_w <- 30
is_on <- function(x) x %in% c("ON", "On", "on", 1, TRUE)

sf_cum <- sf_wide %>%
  filter(between(Time_Relative_sf, 0, trial_len)) %>%
  arrange(file_name, Time_Relative_sf) %>%
  group_by(file_name) %>%
  mutate(dt = pmax(0, Time_Relative_sf - lag(Time_Relative_sf, default = first(Time_Relative_sf)))) %>%
  ungroup() %>%
  pivot_longer(starts_with("Plate_"), names_to = "Subject", values_to = "plate_state") %>%
  mutate(Subject = str_remove(Subject, "^Plate_"), on_plate = is_on(plate_state)) %>%
  group_by(file_name, Subject) %>%
  mutate(cum_onplate_sec = cumsum(dt * as.numeric(on_plate))) %>%
  ungroup()

sf_cum_30 <- sf_cum %>%
  mutate(bin = round(Time_Relative_sf / bin_w) * bin_w) %>%
  group_by(file_name, Subject, bin) %>%
  slice_min(abs(Time_Relative_sf - bin), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(file_name, Subject, time_sec = bin, cum_onplate_sec,
            TailMark = dplyr::recode(Subject, Base = "B", Middle = "M", Tip = "T")) %>%
  left_join(select(lut_trials, file_name, Date, Trial, TailMark, Cohort, MouseID, Age, Group, Sex, Notes),
            by = c("file_name", "TailMark"))

id_vars <- c("Cohort", "Date", "Trial", "Subject", "time_sec", "Group", "Age", "Sex", "MouseID")

# ---- Take average of DoubleObserved trials ----
sf_cum_30_resolved <- bind_rows(
  sf_cum_30 %>% filter(Notes != "DoubleObserved") %>%
    select(all_of(id_vars), cum_onplate_sec) %>%
    mutate(n_obs = 1),
  sf_cum_30 %>% filter(Notes == "DoubleObserved") %>%
    group_by(across(all_of(id_vars))) %>%
    summarise(cum_onplate_sec = mean(cum_onplate_sec, na.rm = TRUE),
              n_obs = n(), .groups = "drop")
)

final_ranked <- sf_cum_30_resolved %>%
  group_by(Cohort, Date, Trial, MouseID) %>%
  slice_max(time_sec, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(Cohort, Date, Trial) %>%
  mutate(Rank = factor(dense_rank(desc(cum_onplate_sec)), levels = 1:3, ordered = TRUE)) %>%
  ungroup() %>%
  select(Cohort, Date, Trial, MouseID, Rank)

sf_rank_ts <- sf_cum_30_resolved %>% left_join(final_ranked, by = c("Cohort", "Date", "Trial", "MouseID")) %>%
  mutate(cum_onplate_min = cum_onplate_sec/60)

rank_summary_time <- sf_rank_ts %>%
  group_by(Group, Sex, Rank, time_sec, Cohort) %>%
  summarise(cum_onplate_min = mean(cum_onplate_min, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Sex, Rank, time_sec) %>%
  summarise(mean_cum = mean(cum_onplate_min, na.rm = TRUE),
            sd_cum = sd(cum_onplate_min, na.rm = TRUE),
            n_cohorts = n(),
            se_cum = sd_cum / sqrt(n_cohorts),
            lower = mean_cum - se_cum, upper = mean_cum + se_cum,
            .groups = "drop")

ggplot(rank_summary_time, aes(time_sec/60, mean_cum, color = Rank, fill = Rank, linetype = Sex)) +
  facet_wrap(~ Group) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 0.5) +
  theme_nolabels
ggsave(file.path(plotdir, "WS_CumDuration.pdf"), width = 7, height = 3.5, units = "cm")

sf_rank_ts20 <- sf_rank_ts %>% filter(time_sec == 1200) %>% droplevels()

m1 <- lmer(cum_onplate_sec ~ Rank * Group * Sex + Age + (1 | MouseID), data = sf_rank_ts20); summary(m1)
Anova(m1, type = 3)
pairs(emmeans(m1, ~ Sex))
pairs(emmeans(m1, ~ Rank))
pairs(emmeans(m1, ~ Group))
pairs(emmeans(m1, ~ Group | Sex))

# ==================================================
# Per-mouse per-file summary (one value per Observed Trial)
# ==================================================
sf_sum <- result_clean %>%
  group_by(file_name, Subject) %>%
  mutate(rleid = data.table::rleid(Plate),
         runLength = ave(seq_along(Plate), rleid, FUN = length)) %>%
  select(-Time_Relative_sf, -Qualifier) %>%
  filter(Plate == "ON") %>%
  distinct() %>%
  mutate(TailMark = dplyr::recode(Subject, Base = "B", Middle = "M", Tip = "T"))

sf_sum_lut <- left_join(sf_sum, lut_trials) %>% ungroup() 

sf_sum_lut <- sf_sum_lut %>%
  group_by(Cohort, Date, Trial, Subject, MouseID, Observer) %>%
  summarise(numBouts = n(),
            TotDur = sum(runLength),
            meanBout_dur = mean(runLength),
            medBout_dur = median(runLength),
            maxBout_dur = max(runLength),
            .groups = "drop")

mouse_info <- lut %>% select(Cohort, Treatment, Group, Sex, Age, MouseID, TailMark) %>% distinct()
cohort_info <- lut %>% select(Cohort, Group, Sex, Age) %>% distinct()

sf_sum_lut <- sf_sum_lut %>% full_join(mouse_info)

write.csv(sf_sum_lut %>% dplyr::select(Cohort, Trial, Observer, Group, Sex, Age, MouseID, TotDur), 
          "WS_allTrials.csv", row.names = FALSE)

# ---- Trial 1 vs. Trial 2 ----
sf_sum_wide <- sf_sum_lut %>% group_by(Cohort, MouseID, Group, Sex, Age, Trial) %>%
  summarise(TotDur = mean(TotDur), .groups = "drop") %>% 
  pivot_wider(names_from = Trial, values_from = TotDur, names_prefix = "Trial_")

ggplot(sf_sum_wide, aes(x=Trial_1/600, y=Trial_2/600)) + 
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  theme_nolabels
ggsave(file.path(plotdir, "WS_T1T2.pdf"), width = 3, height = 3, units = "cm")

sf_sum_wide %>% group_by(Cohort) %>%
  mutate(rank1 = min_rank(desc(Trial_1)),
         rank2 = min_rank(desc(Trial_2))) %>% ungroup() %>%
  mutate(rankDiff = abs(rank2 - rank1)) %>%
  filter(!is.na(rankDiff)) %>%
  count(rankDiff) %>%
  mutate(prop = n / sum(n), rankDiff = factor(rankDiff, levels = c(0, 1, 2)))
 
cor.test(sf_sum_wide$Trial_1, sf_sum_wide$Trial_2) # 0.4312807

c1 <- lm(Trial_2 ~ Trial_1 * Group * Sex + Age, data = sf_sum_wide); summary(c1)
Anova(c1, type = 3)

c1a <- lm(Trial_2 ~ Trial_1 * Sex + Age, data = sf_sum_wide %>% filter(Group == "BalbC_WT")); summary(c1a)
c1b <- lm(Trial_2 ~ Trial_1 * Sex + Age, data = sf_sum_wide %>% filter(Group == "C57_WT")); summary(c1b)
c1c <- lm(Trial_2 ~ Trial_1 * Sex, data = sf_sum_wide %>% filter(Group == "C57_MUT")); summary(c1c)

cor_result <- sf_sum_wide %>%
  group_by(Group) %>%
  dplyr::summarise(
    r = cor(Trial_1, Trial_2, use = "complete.obs"),
    p = cor.test(Trial_1, Trial_2)$p.value,
    n = n())

# ---- Repeatability or Intraclass Correlation Coefficient (ICC) ----
icc_df <- sf_sum_lut %>% group_by(Group, Sex) %>% tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data, ~ lmer(TotDur ~ Trial + (1 | MouseID), data = .x)),
    icc = purrr::map(model, ~ performance::icc(.x, ci = 0.95, iterations = 5000)))

icc_plot <- icc_df %>% dplyr::mutate(
    ICC     = map_dbl(icc, ~ .x$ICC_adjusted[1]),
    CI_low  = map_dbl(icc, ~ .x$ICC_adjusted[2]),
    CI_high = map_dbl(icc, ~ .x$ICC_adjusted[3])) %>% dplyr::select(Group, Sex, ICC, CI_low, CI_high)

ggplot(icc_plot, aes(x = Group, y = ICC, color = Sex)) +
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high),
    position = position_dodge(width = .9), width = 0) +
  geom_point(position = position_dodge(width = .9), size = 3) +
  ylim(0, .8) + theme_nolabels
ggsave(file.path(plotdir, "WS_ICC.pdf"), width = 3, height = 3, units = "cm")

# ---- Function to calculate ICC
get_icc <- function(dat, covariate) {
  form <- if (covariate == "Group") TotDur ~ Trial + Group + (1|MouseID) else TotDur ~ Trial + Sex + (1|MouseID)
  m <- lmer(form, data = dat)
  vc <- as.data.frame(VarCorr(m))
  vc$vcov[vc$grp == "MouseID"] / sum(vc$vcov)
}

set.seed(123); B <- 5000

# ---- Sex: F vs M, controlling for Group
icc_F <- get_icc(filter(sf_sum_lut, Sex == "F"), "Group")
icc_M <- get_icc(filter(sf_sum_lut, Sex == "M"), "Group")
obs_sex <- icc_F - icc_M

mice <- sf_sum_lut %>% distinct(MouseID, Sex, Group)

perm_sex <- replicate(B, {
  perm_labels <- mice %>% group_by(Group) %>% mutate(Sex_perm = sample(Sex)) %>% ungroup()
  d <- sf_sum_lut %>% select(-Sex) %>% left_join(select(perm_labels, MouseID, Sex_perm), by = "MouseID")
  F <- d %>% filter(Sex_perm == "F") %>% rename(Sex = Sex_perm)
  M <- d %>% filter(Sex_perm == "M") %>% rename(Sex = Sex_perm)
  get_icc(F, "Group") - get_icc(M, "Group")
})

p_sex <- (sum(abs(perm_sex) >= abs(obs_sex)) + 1) / (B + 1)
c(ICC_F = icc_F, ICC_M = icc_M, difference = obs_sex, p = p_sex)

# ---- BalbC_WT vs C57_WT, controlling for Sex
d_BC <- sf_sum_lut %>% filter(Group %in% c("BalbC_WT", "C57_WT"))
icc_Balb <- get_icc(filter(d_BC, Group == "BalbC_WT"), "Sex")
icc_WT <- get_icc(filter(d_BC, Group == "C57_WT"), "Sex")
obs_BC <- icc_Balb - icc_WT

mouse_BC <- d_BC %>% distinct(MouseID, Group, Sex)

perm_BC <- replicate(B, {
  labs <- mouse_BC %>% group_by(Sex) %>% mutate(Group_perm = sample(Group)) %>% ungroup()
  d <- d_BC %>% select(-Group) %>% left_join(select(labs, MouseID, Group_perm), by = "MouseID")
  B <- d %>% filter(Group_perm == "BalbC_WT") %>% rename(Group = Group_perm)
  W <- d %>% filter(Group_perm == "C57_WT") %>% rename(Group = Group_perm)
  get_icc(B, "Sex") - get_icc(W, "Sex")
})

p_BC <- (sum(abs(perm_BC) >= abs(obs_BC)) + 1) / (B + 1)
c(ICC_BalbC = icc_Balb, ICC_C57WT = icc_WT, difference = obs_BC, p = p_BC)

# ---- C57_WT vs C57_MUT, controlling for Sex
d_WM <- sf_sum_lut %>% filter(Group %in% c("C57_WT", "C57_MUT"))
icc_WT2 <- get_icc(filter(d_WM, Group == "C57_WT"), "Sex")
icc_MUT <- get_icc(filter(d_WM, Group == "C57_MUT"), "Sex")
obs_WM <- icc_WT2 - icc_MUT

mouse_WM <- d_WM %>% distinct(MouseID, Group, Sex)

perm_WM <- replicate(B, {
  labs <- mouse_WM %>% group_by(Sex) %>% mutate(Group_perm = sample(Group)) %>% ungroup()
  d <- d_WM %>% select(-Group) %>% left_join(select(labs, MouseID, Group_perm), by = "MouseID")
  W <- d %>% filter(Group_perm == "C57_WT") %>% rename(Group = Group_perm)
  M <- d %>% filter(Group_perm == "C57_MUT") %>% rename(Group = Group_perm)
  get_icc(W, "Sex") - get_icc(M, "Sex")
})

p_WM <- (sum(abs(perm_WM) >= abs(obs_WM)) + 1) / (B + 1)
c(ICC_C57WT = icc_WT2, ICC_C57MUT = icc_MUT, difference = obs_WM, p = p_WM)

# ---- BalbC_WT vs C57_MUT, controlling for Sex
d_BM <- sf_sum_lut %>% filter(Group %in% c("BalbC_WT", "C57_MUT"))
icc_Balb2 <- get_icc(filter(d_BM, Group == "BalbC_WT"), "Sex")
icc_MUT2 <- get_icc(filter(d_BM, Group == "C57_MUT"), "Sex")
obs_BM <- icc_Balb2 - icc_MUT2

mouse_BM <- d_BM %>% distinct(MouseID, Group, Sex)

perm_BM <- replicate(B, {
  labs <- mouse_BM %>% group_by(Sex) %>% mutate(Group_perm = sample(Group)) %>% ungroup()
  d <- d_BM %>% select(-Group) %>% left_join(select(labs, MouseID, Group_perm), by = "MouseID")
  Bc <- d %>% filter(Group_perm == "BalbC_WT") %>% rename(Group = Group_perm)
  M <- d %>% filter(Group_perm == "C57_MUT") %>% rename(Group = Group_perm)
  get_icc(Bc, "Sex") - get_icc(M, "Sex")
})

p_BM <- (sum(abs(perm_BM) >= abs(obs_BM)) + 1) / (B + 1)
c(ICC_BalbC = icc_Balb2, ICC_C57MUT = icc_MUT2, difference = obs_BM, p = p_BM)

p.adjust(c(BalbC_vs_C57WT = p_BC, C57WT_vs_C57MUT = p_WM, BalbC_vs_C57MUT = p_BM), method = "holm")

# ---- Take average of DoubleObserved trials ----
sf_sum_lut_resolved <- sf_sum_lut %>% 
  group_by(Cohort, Date, Trial, Subject, MouseID, Group, Sex, Age) %>%
  summarise(across(c(numBouts, TotDur, meanBout_dur, medBout_dur, maxBout_dur),
                   ~ mean(.x, na.rm = TRUE)), n = n(), 
            .groups = "drop")

sf_sum_lut_resolved_avg <- sf_sum_lut_resolved %>% 
  group_by(Cohort, MouseID, Group, Sex, Age) %>%
  summarise(across(c(numBouts, TotDur, meanBout_dur, medBout_dur, maxBout_dur),
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>% mutate(across(7:10, ~ .x / 10)) # convert to seconds

# ---- Plot TotDur ----
TotDur_se <- summarySE(sf_sum_lut_resolved_avg, measurevar = "TotDur", groupvars = c("Group", "Sex"), na.rm = TRUE)

ggplot(sf_sum_lut_resolved_avg, aes(Sex, TotDur, color = Sex)) +
  facet_wrap(~ Group) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = TotDur_se, aes(ymin = TotDur - se, ymax = TotDur + se),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  theme_nolabels
ggsave(file.path(plotdir, "WS_TotDur.pdf"), width = 3.5, height = 3.5, units = "cm")

# ---- TotDur (all mice; per Trial) ----
TotDur_All <- sf_sum_lut %>% 
  group_by(Cohort, Trial, Observer, Group, Sex, Age) %>% 
  summarise(TotDur_All = (sum(TotDur))/600, .groups = "drop") %>%
  mutate(TotDur_All = case_when(TotDur_All > 20 ~ 20, .default = TotDur_All)) %>%
  group_by(Cohort, Trial, Group, Sex, Age) %>%
  summarise(TotDur_All = mean(TotDur_All), .groups = "drop")

summarySE(TotDur_All, measurevar = "TotDur_All", na.rm = TRUE)

# ---- Plot numBouts ----
numBouts_se <- summarySE(sf_sum_lut_resolved_avg, measurevar = "numBouts", groupvars = c("Group", "Sex"), na.rm = TRUE)

ggplot(sf_sum_lut_resolved_avg, aes(Sex, numBouts, color = Sex)) +
  facet_wrap(~ Group) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = numBouts_se, aes(ymin = numBouts - se, ymax = numBouts + se),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  theme_nolabels
ggsave(file.path(plotdir, "WS_numBouts.pdf"), width = 3.5, height = 3.5, units = "cm")

m2 <- lmer(numBouts ~ Group * Sex + Age + (1 | MouseID), data = sf_sum_lut); summary(m2)
pairs(emmeans(m2, ~ Group))
pairs(emmeans(m2, ~ Group | Sex)) 

# ==================================================
# Inter-observer consistency (MW vs RP vs TR)
# ==================================================
consistency <- sf_sum_lut %>%
  group_by(Date, Cohort, MouseID, Trial, Subject) %>%
  filter(n_distinct(Observer) >= 2) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Date, Cohort, MouseID, Trial, Subject),
              names_from = Observer,
              values_from = c(numBouts, TotDur, meanBout_dur, medBout_dur, maxBout_dur),
              names_sep = "_")

observers <- consistency %>%
  select(matches("_(MW|RP|TR)$")) %>%
  names() %>%
  str_extract("(?<=_)\\w+$") %>%
  unique()

metrics <- c("numBouts", "TotDur", "meanBout_dur", "medBout_dur", "maxBout_dur")
observer_pairs <- as_tibble(t(combn(observers, 2)), .name_repair = "minimal") %>%
  setNames(c("obs1", "obs2"))

compute_absolute_consistency <- function(df, metric, obs1, obs2) {
  col1 <- paste0(metric, "_", obs1); col2 <- paste0(metric, "_", obs2)
  if (!(col1 %in% names(df) && col2 %in% names(df))) return(tibble())

  df2 <- df %>%
    filter(!is.na(.data[[col1]]), !is.na(.data[[col2]])) %>%
    mutate(abs_diff = abs(.data[[col1]] - .data[[col2]]),
           max_val = pmax(.data[[col1]], .data[[col2]]),
           abs_consistency = 1 - abs_diff / max_val)

  tibble(metric = metric, obs1 = obs1, obs2 = obs2, n = nrow(df2),
         mean_abs_consistency = mean(df2$abs_consistency, na.rm = TRUE),
         median_abs_consistency = median(df2$abs_consistency, na.rm = TRUE),
         sd_abs_consistency = sd(df2$abs_consistency, na.rm = TRUE))
}

consistency_absolute <- map_dfr(metrics, \(m) {
  map_dfr(seq_len(nrow(observer_pairs)), \(i) {
    compute_absolute_consistency(consistency, m, observer_pairs$obs1[i], observer_pairs$obs2[i])
  })
}) %>% filter(!is.na(mean_abs_consistency))

consistency_absolute %>% group_by(metric) %>%
  summarise(mean_abs_consistency = mean(mean_abs_consistency))

# ==================================================
# Displacement transitions (groupState) + Elo-format events
# ==================================================
file_info <- lut_trials %>% select(file_name, Cohort, Trial, Observer) %>% distinct()

sf_wide_filled <- sf_wide %>%
  left_join(file_info, by = "file_name") %>%
  mutate(Time_Relative_sf = as.numeric(Time_Relative_sf),
         groupState = if_else(Time_Relative_sf == 0, "All_OFF", groupState),
         Time_Relative_sf = if_else(Trial == 2, Time_Relative_sf + 1200, Time_Relative_sf)) %>%
  group_by(file_name) %>%
  arrange(Time_Relative_sf, .by_group = TRUE) %>%
  fill(groupState, .direction = "updown") %>%
  ungroup()

sf_collapsed <- sf_wide_filled %>%
  arrange(Cohort, Trial, Observer, Time_Relative_sf) %>%
  group_by(Cohort, Trial, Observer) %>%
  mutate(rleid = data.table::rleid(groupState)) %>%
  group_by(Cohort, Trial, Observer, rleid) %>%
  mutate(runLength = n()) %>%
  ungroup() %>%
  mutate(state2 = if_else(groupState == "All_OFF" & runLength < 30, NA_character_, groupState)) %>%
  group_by(Cohort, Observer) %>%
  fill(state2, .direction = "downup") %>%
  ungroup()

transitions <- sf_collapsed %>%
  arrange(Cohort, Trial, Observer, Time_Relative_sf) %>%
  group_by(Cohort, Trial, Observer) %>%
  mutate(PrevState = lag(state2), PrevTime = lag(Time_Relative_sf)) %>%
  ungroup() %>%
  filter(!is.na(PrevState), !is.na(state2), state2 != PrevState) %>%
  transmute(Cohort, Trial, Observer, Time_trans = Time_Relative_sf, Pre = PrevState, Post = state2)

sf_trim_wide <- transitions %>%
  mutate(Transition = case_when(
    Pre %in% c("BaseMiddle", "Middle_Only") & Post %in% c("Base_Only")   ~ "B displace M",
    Pre %in% c("BaseTip",    "Tip_Only")    & Post %in% c("Base_Only")   ~ "B displace T",
    Pre %in% c("BaseMiddle", "Base_Only")   & Post %in% c("Middle_Only") ~ "M displace B",
    Pre %in% c("MiddleTip",  "Tip_Only")    & Post %in% c("Middle_Only") ~ "M displace T",
    Pre %in% c("BaseTip",    "Base_Only")   & Post %in% c("Tip_Only")    ~ "T displace B",
    Pre %in% c("MiddleTip",  "Middle_Only") & Post %in% c("Tip_Only")    ~ "T displace M",
    TRUE ~ NA_character_
  ))

ws <- sf_trim_wide %>%
  mutate(parts = str_split(as.character(Transition), " displace ", n = 2),
         winner = vapply(parts, `[`, character(1), 1),
         loser  = vapply(parts, `[`, character(1), 2)) %>%
  select(-Pre, -Post, -parts) %>%
  drop_na(Transition)

ws_elo <- ws %>%
  mutate(Cohort = as.character(Cohort), Observer = as.character(Observer)) %>%
  group_by(Cohort, Observer) %>%
  arrange(Time_trans, .by_group = TRUE) %>%
  mutate(FakeTrial_Date = as.Date("2000-01-01") + row_number() - 1L) %>%
  ungroup()

# ==================================================
# Hierarchy metrics (Landau h, David's scores)
# ==================================================
get_mat <- function(df) EloRating::creatematrix(winners = df$winner, losers = df$loser)

metrics <- ws_elo %>%
  group_by(Cohort, Observer) %>%
  summarise(mat = list(get_mat(pick(everything()))), .groups = "drop") %>%
  mutate(landau = map(mat, ~ EloRating::h.index(.x, loops = 1000)),
         ds     = map(mat, ~ EloRating::DS(.x)))

cohort_unpacked <- metrics %>%
  transmute(Cohort, Observer,
            landau_h      = map_dbl(landau, ~ .x$value[.x$variable == "h index"]),
            landau_hprime = map_dbl(landau, ~ .x$value[.x$variable == "h' index"]),
            landau_p      = map_dbl(landau, ~ .x$value[.x$variable == "p right"]))

metrics_mouse <- metrics %>%
  select(Cohort, Observer, ds) %>%
  unnest(ds) %>%
  rename(TailMark = ID, DS_raw = DS, DS_norm = normDS)

metrics_all_avg <- metrics_mouse %>%
  left_join(cohort_unpacked, by = c("Cohort", "Observer")) %>%
  group_by(Cohort, TailMark) %>%
  summarise(DS_raw = mean(DS_raw), DS_norm = mean(DS_norm),
            landau_h = mean(landau_h), landau_hprime = mean(landau_hprime), landau_p = mean(landau_p),
            n_obs = n(), .groups = "drop")

# ==================================================
# Per-mouse summary (one value per two consecutive trials)
# ==================================================
result_clean_con <- result_clean %>%
  left_join(file_info, by = "file_name") %>%
  mutate(Time_Relative_sf = if_else(Trial == 2, Time_Relative_sf + 1200, Time_Relative_sf))

sf_sum_con <- result_clean_con %>%
  group_by(Cohort, Observer, Subject) %>%
  arrange(Time_Relative_sf, .by_group = TRUE) %>%
  mutate(rleid = data.table::rleid(Plate),
         runLength = ave(seq_along(Plate), rleid, FUN = length)) %>%
  select(-Time_Relative_sf, -Qualifier) %>%
  filter(Plate == "ON") %>%
  distinct() %>%
  summarise(numBouts = n(),
            TotDur = sum(runLength),
            meanBout_dur = mean(runLength),
            medBout_dur = median(runLength),
            maxBout_dur = max(runLength),
            .groups = "drop") %>%
  mutate(TailMark = dplyr::recode(Subject, Base = "B", Middle = "M", Tip = "T"),
         across(c(TotDur, meanBout_dur, medBout_dur, maxBout_dur), ~ .x / 10))  # 0.1 s -> sec

# ---- Take average of DoubleObserved trials ----
sf_sum_con_avg <- sf_sum_con %>%
  group_by(Cohort, Subject, TailMark) %>%
  summarise(numBouts = mean(numBouts),
            TotDur = mean(TotDur),
            meanBout_dur = mean(meanBout_dur),
            medBout_dur = median(medBout_dur),
            maxBout_dur = max(maxBout_dur),
            n_obs = n(), .groups = "drop")

sf_sum_con_avg <- sf_sum_con_avg %>% full_join(mouse_info)

# ---- Plot medBout_dur ----
medBout_se <- summarySE(sf_sum_con_avg, measurevar = "medBout_dur", groupvars = c("Group", "Sex"), na.rm = TRUE)

ggplot(sf_sum_con_avg, aes(Sex, medBout_dur, color = Sex)) +
  facet_wrap(~ Group) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = medBout_se, aes(ymin = medBout_dur - se, ymax = medBout_dur + se),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  scale_y_continuous(trans = "log10") +
  theme_nolabels
ggsave(file.path(plotdir, "WS_BoutDuration.pdf"), width = 3.5, height = 3.5, units = "cm")

m3 <- lm(medBout_dur ~ Group * Sex + Age, data = sf_sum_con_avg); summary(m3)
pairs(emmeans(m3, ~ Group))
pairs(emmeans(m3, ~ Group | Sex))

# ---- Summarize bout data ----
# Rule: Take mean of Trial 1 and Trial 2 for numBouts and TotDur (sf_sum_lut_resolved_avg)
# Rule: Take grand mean/median/maximum across Trials 1 and 2 for meanBout_dur, medBout_dur, and maxBout_dur (sf_sum_con_avg)

mouse_summary <- sf_sum_lut_resolved_avg %>% select(-meanBout_dur, -medBout_dur, -maxBout_dur) %>% 
  full_join(sf_sum_con_avg %>% select(-numBouts, -TotDur))

# ==================================================
# Network analysis (displacement edges)
# ==================================================
dsp_edges <- ws_elo %>%
  count(Cohort, Trial, Observer, winner, loser, Transition, name = "n_displacements")

Tot_dsp_edges <- dsp_edges %>%
  group_by(Cohort, Trial, Observer) %>%
  summarise(Tot_dsp = sum(n_displacements, na.rm = TRUE), .groups = "drop") %>%
  left_join(cohort_info)

# ---- Plot Total Displacements per Cohort ----
dsp_se <- summarySE(Tot_dsp_edges, measurevar = "Tot_dsp", groupvars = c("Group", "Sex"), na.rm = TRUE)

ggplot(dsp_se, aes(Sex, Tot_dsp, fill = Sex)) +
  facet_grid(~ Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  geom_errorbar(aes(ymin = Tot_dsp - se, ymax = Tot_dsp + se),
                position = position_dodge(.9), width = 0) +
  theme_nolabels
ggsave(file.path(plotdir, "WS_TotDsp.pdf"), width = 3.5, height = 3.5, units = "cm")

m4 <- lmer(Tot_dsp ~ Group * Sex + Age + (1 | Cohort), data = Tot_dsp_edges); summary(m4)
pairs(emmeans(m4, ~ Group))
pairs(emmeans(m4, ~ Group | Sex))

# ---- Make Network Graphs ----
nodes <- sf_sum_lut %>%
  mutate(across(c(Cohort, Trial, Observer, TailMark, MouseID), ~ trimws(as.character(.x)))) %>%
  select(Cohort, Trial, Observer, TailMark, MouseID, TotDur)

edges <- dsp_edges %>%
  mutate(across(c(Cohort, Trial, Observer, winner, loser), ~ trimws(as.character(.x)))) %>%
  left_join(nodes %>% transmute(Cohort, Trial, Observer, TailMark, from = MouseID),
    by = c("Cohort", "Trial", "Observer", "winner" = "TailMark")) %>%
  left_join(nodes %>% transmute(Cohort, Trial, Observer, TailMark, to = MouseID),
    by = c("Cohort", "Trial", "Observer", "loser" = "TailMark")) %>%
  transmute(Cohort, Trial, Observer, from, to, weight = n_displacements)

out_plot_dir <- file.path(plotdir, "network_plots")
key_cols <- c("Cohort","Trial","Observer")

nodes <- nodes %>% mutate(across(all_of(key_cols), as.character), MouseID = as.character(MouseID))
edges <- edges %>% mutate(across(all_of(key_cols), as.character),
                          from = as.character(from), to = as.character(to), weight = as.numeric(weight))

scale_to <- function(x, to = c(12, 32)){x <- as.numeric(x); r <- range(x, na.rm = TRUE)
if(!is.finite(r[1])||diff(r) == 0)return(rep(mean(to), length(x)))
(x-r[1])/diff(r)*diff(to)+to[1]}

centralization_from <- function(x){
  x <- as.numeric(x); x <- x[is.finite(x)]
  n <- length(x); xmax <- max(x); denom <- (n-1)*xmax
  if(denom == 0)return(0); sum(xmax-x)/denom}

nodes_inst <- nodes %>% group_by(Cohort, Trial, Observer, MouseID) %>%
  summarise(TotDur = sum(TotDur, na.rm = TRUE), .groups = "drop")

instances <- inner_join(edges %>% distinct(Cohort, Trial, Observer),
                        nodes_inst %>% distinct(Cohort, Trial, Observer),
                        by = key_cols) %>%
  arrange(Cohort, Trial, Observer) %>%
  mutate(instance_id = paste(Cohort, Trial, Observer, sep = "__"))

for(i in seq_len(nrow(instances))){
  k <- instances[i, ]
  vdf <- nodes_inst %>% filter(Cohort == k$Cohort, Trial == k$Trial, Observer == k$Observer) %>%
    transmute(name = MouseID, TotDur)
  edf <- edges %>% filter(Cohort == k$Cohort, Trial == k$Trial, Observer == k$Observer) %>%
    select(from, to, weight)
  g <- graph_from_data_frame(edf, directed = TRUE, vertices = vdf)
  
  V(g)$size <- scale_to(sqrt(V(g)$TotDur), to = c(10, 22))
  E(g)$width <- if(ecount(g)>0)scale_to(E(g)$weight, to = c(1.5, 8))else 1.5
  
  w <- if(ecount(g)>0)E(g)$weight else numeric(0)
  out_str <- strength(g, mode = "out", weights = w); in_str  <- strength(g, mode = "in",  weights = w); net_str <- out_str - in_str
  vord <- names(net_str)[order(-net_str, names(net_str))]
  vids <- match(vord, V(g)$name)
  coords <- layout_in_circle(g, order = vids); coords <- cbind(-coords[, 2], coords[, 1])
  
  png(file.path(out_plot_dir, paste0("NetworkPlot_", k$instance_id, ".png")),
      width = 800, height = 800, res = 150)
  plot(g, layout = coords,
       main = paste("Cohort", k$Cohort, "| Trial", k$Trial, "| Observer", k$Observer),
       vertex.label = V(g)$name, vertex.label.cex = 0.8,
       edge.arrow.size = 0.6, edge.curved = 0.25)
  dev.off()
}

instance_stats <- pmap_dfr(instances, function(Cohort, Trial, Observer, instance_id){
  vdf <- nodes_inst %>% filter(Cohort ==!!Cohort, Trial ==!!Trial, Observer ==!!Observer) %>%
    transmute(name = MouseID, TotDur)
  edf <- edges %>% filter(Cohort ==!!Cohort, Trial ==!!Trial, Observer ==!!Observer) %>%
    select(from, to, weight)
  g <- graph_from_data_frame(edf, directed = TRUE, vertices = vdf)
  
  w <- if(ecount(g)>0)E(g)$weight else numeric(0)
  out_str <- strength(g, mode = "out", weights = w); in_str  <- strength(g, mode = "in",  weights = w); net_str <- out_str - in_str
  vord <- names(net_str)[order(-net_str, names(net_str))]
  Aord <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))[vord, vord]
  
  downward <- Aord[1, 2]+Aord[1, 3]+Aord[2, 3]
  upward <- Aord[2, 1]+Aord[3, 1]+Aord[3, 2]
  linearity_index <- downward/(downward+upward+1e-9)
  
  tibble(Cohort = Cohort, Trial = Trial, Observer = Observer, instance_id = instance_id,
         n_nodes = vcount(g), n_edges = ecount(g), total_weight = sum(w),
         density = edge_density(g, loops = FALSE),
         reciprocity = reciprocity(g, ignore.loops = TRUE),
         out_centralization_w = centralization_from(out_str),
         in_centralization_w = centralization_from(in_str),
         despotism = max(out_str)/(sum(out_str)+1e-9),
         linearity_index = linearity_index)
})
instance_stats <- full_join(instance_stats, cohort_info)

node_stats <- pmap_dfr(instances, function(Cohort, Trial, Observer, instance_id){
  vdf <- nodes_inst %>% filter(Cohort ==!!Cohort, Trial ==!!Trial, Observer ==!!Observer) %>%
    transmute(name = MouseID, TotDur)
  edf <- edges %>% filter(Cohort ==!!Cohort, Trial ==!!Trial, Observer ==!!Observer) %>%
    select(from, to, weight)
  g <- graph_from_data_frame(edf, directed = TRUE, vertices = vdf)
  
  out_strength <- strength(g, mode = "out", weights = E(g)$weight)
  in_strength  <- strength(g, mode = "in",  weights = E(g)$weight)
  
  tibble(Cohort = Cohort, Trial = Trial, Observer = Observer, instance_id = instance_id,
         MouseID = V(g)$name, TotDur = V(g)$TotDur,
         out_strength = out_strength, in_strength = in_strength,
         net_strength = out_strength-in_strength,
         pagerank = page_rank(g, directed = TRUE, weights = E(g)$weight)$vector)
})
node_stats <- full_join(node_stats, mouse_info)

net_stats <- full_join(node_stats, instance_stats); write.csv(net_stats, "net_stats.csv", row.names = FALSE)

net_stats_avg <- net_stats %>%
  group_by(Cohort, MouseID, Group, Sex, Age) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop") %>% select(-TotDur, -n_nodes, -density, -reciprocity)

# ==================================================
# Network Analysis Statistics
# ==================================================
cohort_metrics <- c("out_centralization_w", "in_centralization_w", "despotism", "linearity_index")

for (v in cohort_metrics) {
  p <- ggplot(instance_stats, aes(Sex, .data[[v]], fill = Sex)) +
    geom_boxplot() +
    geom_point() +
    facet_wrap(~ Group) +
    labs(x = NULL, y = v, title = v) +
    theme_classic()
  print(p)
}

plot(instance_stats$despotism, instance_stats$out_centralization_w)

m5 <- lmer(out_centralization_w ~ Group * Sex + Age + (1 | Cohort), instance_stats); summary(m5)
pairs(emmeans(m5, ~ Sex)) # NS
pairs(emmeans(m5, ~ Group)) # NS
pairs(emmeans(m5, ~ Group | Sex)) # NS

m6 <- lmer(in_centralization_w ~ Group * Sex + Age + (1 | Cohort), instance_stats); summary(m6)
pairs(emmeans(m6, ~ Sex)) # NS
pairs(emmeans(m6, ~ Group)) # NS
pairs(emmeans(m6, ~ Group | Sex)) # NS

m7 <- lmer(despotism ~ Group * Sex + Age + (1 | Cohort), instance_stats); summary(m7)
pairs(emmeans(m7, ~ Sex)) # NS
pairs(emmeans(m7, ~ Group)) # NS
pairs(emmeans(m7, ~ Group | Sex)) # NS

m8 <- lmer(linearity_index ~ Group * Sex + Age + (1 | Cohort), instance_stats); summary(m8)
pairs(emmeans(m8, ~ Sex)) # NS
pairs(emmeans(m8, ~ Group)) # NS
pairs(emmeans(m8, ~ Group | Sex)) # NS

# ---- Net Displacement and TotDur ----
ggplot(net_stats, aes(x = net_strength, y = TotDur, group = Sex, color = Sex)) +
  facet_wrap(~ Group) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(x = "out_strength - in_strength", y = "Total Duration (s)")

m9 <- lmer(TotDur ~ net_strength * Group * Sex + Age + (1 | MouseID), net_stats); summary(m9)
emm_slope <- emtrends(m9, specs = c("Group", "Sex"), var = "net_strength"); summary(emm_slope)

# ==================================================
# Make custom Network Graphs
# ==================================================
instance_ids <- c("C26_2_RP", "C72_2_RP", "C2_MW.R1_2_MW", "C50_2_RP")

nodes2 <- nodes %>% mutate(instance_id = paste(Cohort, Trial, Observer, sep = "_"))
edges2 <- edges %>% mutate(instance_id = paste(Cohort, Trial, Observer, sep = "_"))

nodes_inst <- nodes2 %>% group_by(instance_id, MouseID) %>%
  summarise(TotDur = sum(TotDur, na.rm = TRUE), .groups = "drop")

nodes_sel <- nodes_inst %>% filter(instance_id %in% instance_ids)
edges_sel <- edges2 %>% filter(instance_id %in% instance_ids)

global_node_range <- range(nodes_sel$TotDur, na.rm = TRUE)
global_edge_range <- range(edges_sel$weight, na.rm = TRUE)

graphs <- map(instance_ids, function(id) {
  vdf <- nodes_inst %>% filter(instance_id == id) %>% transmute(name = MouseID, TotDur)
  edf <- edges2 %>% filter(instance_id == id) %>% transmute(from, to, weight)
  graph_from_data_frame(edf, directed = TRUE, vertices = vdf)
})

plot_list <- imap(graphs, function(g, i)
  ggraph(g, layout = "circle") +
    geom_edge_fan(aes(edge_width = weight), alpha = 0.8) +
    geom_node_point(aes(size = TotDur), shape = 21, stroke = 0.4) +
    scale_edge_width(range = c(0.3, 2.5), limits = global_edge_range,
                     name = "Edge weight", breaks = pretty(global_edge_range, n = 3)) +
    scale_size(range = c(3, 15), limits = global_node_range,
               name = "Total Duration", breaks = pretty(global_node_range, n = 3)) +
    ggtitle(instance_ids[i]) +
    theme_void() +
    theme(plot.title = element_text(size = 9, hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7))
)

wrap_plots(plot_list, ncol = 2)
ggsave(file.path(plotdir, "WS_NetworkGraphs.pdf"), width = 16, height = 12, units = "cm")

# ==================================================
# Hierarchy classification from TotDur gaps
# ==================================================
mouse_summary_full <- full_join(mouse_summary, metrics_all_avg)

TD_plot <- mouse_summary_full %>%
  filter(!is.na(TotDur)) %>%
  group_by(Cohort) %>%
  arrange(desc(TotDur), .by_group = TRUE) %>%
  mutate(rank_TD = row_number()) %>%
  ungroup()

ggplot(TD_plot, aes(rank_TD, TotDur, group = Cohort)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  facet_grid(Sex ~ Group) +
  scale_x_continuous(breaks = 1:3, labels = c("Top", "Mid", "Low")) +
  theme_classic() +
  labs(x = "Within-cohort TotDur rank", y = "TotDur")

global_dists <- mouse_summary_full %>%
  select(Cohort, MouseID, TotDur) %>%
  group_by(Cohort) %>%
  filter(n() == 3) %>%
  summarise(dist = list(as.numeric(dist(TotDur))), .groups = "drop") %>%
  unnest(dist)

cohort_dists <- mouse_summary_full %>%
  select(Cohort, MouseID, TotDur) %>%
  group_by(Cohort) %>%
  filter(n() == 3) %>%
  summarise(TD_sorted = list(sort(TotDur, decreasing = TRUE)), .groups = "drop") %>%
  mutate(d1 = map_dbl(TD_sorted, 1), d2 = map_dbl(TD_sorted, 2), d3 = map_dbl(TD_sorted, 3),
         delta_12 = d1 - d2, delta_23 = d2 - d3, delta_13 = d1 - d3)

cuts <- quantile(cohort_dists$delta_13, c(0.25, 0.75), na.rm = TRUE)
small_cut <- unname(cuts[1])

primary_classified <- cohort_dists %>%
  mutate(hierarchy_type = case_when(
    delta_13 <  small_cut ~ "weak",
    delta_13 >= small_cut ~ "hierarchical",
    TRUE ~ NA_character_
  ))

cohort_classified <- primary_classified %>%
  mutate(hierarchy_type = case_when(
    hierarchy_type != "hierarchical" ~ hierarchy_type,
    d2 <= (d3 + delta_13 / 3)        ~ "despotic",
    d2 >= (d1 - delta_13 / 3)        ~ "co-dominant",
    TRUE                             ~ "transitive"
  ))

mouse_summary_full <- mouse_summary_full %>%
  left_join(select(cohort_classified, Cohort, hierarchy_type), by = "Cohort")

cohort_level <- mouse_summary_full %>%
  distinct(Cohort, Group, Sex, hierarchy_type) %>%
  mutate(hierarchy_type = factor(hierarchy_type,
                                 levels = c("despotic", "transitive", "co-dominant", "weak")))

ggplot(cohort_level, aes(Sex, fill = hierarchy_type)) +
  geom_bar(position = "fill") +
  facet_wrap(~ Group) +
  theme_nolabels
ggsave(file.path(plotdir, "WS_HierarchyType.pdf"), width = 3, height = 2.5, units = "cm") 

fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "BalbC_WT"))) # NS
fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "C57_WT")))   # NS
fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "C57_MUT")))  # NS

fisher.test(xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "F"))) # NS
fisher.test(xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "M"))) # NS

# ---- Combine C57_WT and C57_MUT ----
cohort_level <- cohort_level %>% 
  mutate(Background = if_else(str_detect(Group, "BalbC"), "BalbC", "C57"))

fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Background == "C57"))) # 0.008901

# ---- Proportions ----
hierarchy_props <- cohort_level %>% 
  dplyr::count(Group, Sex, hierarchy_type) %>%
  group_by(Group, Sex) %>% mutate(percent = 100 * n / sum(n)) %>% ungroup()

cohort_level %>% 
  dplyr::count(Background, Sex, hierarchy_type) %>%
  group_by(Background, Sex) %>% mutate(percent = 100 * n / sum(n)) %>% ungroup()

# ==================================================
# Net displacements by Rank (total duration)
# ==================================================
mouse_summary_full <- mouse_summary_full %>% full_join(net_stats_avg) %>%
  group_by(Cohort) %>%
  mutate(TotDur_rank = factor(dense_rank(desc(TotDur)), levels = 1:3, ordered = TRUE),
         DS_rank = factor(dense_rank(desc(DS_norm)), levels = 1:3, ordered = TRUE)) %>%
  ungroup()

ggplot(mouse_summary_full, aes(Sex, net_strength, fill = TotDur_rank)) +
  facet_wrap(~ Group) +
  geom_boxplot(position = position_dodge(.7), width =.7, outlier.shape = NA) +
  theme_nolabels
ggsave(file.path(plotdir, "WS_netStrength.pdf"), width = 7, height = 3.5, units = "cm")

m10 <- lm(net_strength ~ TotDur_rank * Group * Sex + Age, data = mouse_summary_full); summary(m10)
Anova(m10, type = 3)
pairs(emmeans(m10, ~ TotDur_rank))
pairs(emmeans(m10, ~ TotDur_rank | Sex | Group))

m10a <- lm(net_strength ~ TotDur_rank + Age, data = mouse_summary_full %>% filter(Group == "BalbC_WT" & Sex == "M")); summary(m10a); Anova(m10a, type = 3) # NS
m10b <- lm(net_strength ~ TotDur_rank + Age, data = mouse_summary_full %>% filter(Group == "C57_WT" & Sex == "M")); summary(m10b); Anova(m10b, type = 3)   # 0.009962 **
m10c <- lm(net_strength ~ TotDur_rank, data = mouse_summary_full %>% filter(Group == "C57_MUT" & Sex == "M")); summary(m10c); Anova(m10c, type = 3)        # NS

null_zero <- test(emmeans(m10, ~ Group * Sex * TotDur_rank), null = 0)

# ---- Combine Ranks 2 & 3 ----
rankBin <- mouse_summary_full %>% mutate(Rank_binary = ifelse(TotDur_rank == 1, "Dominant", "Subordinate"))

m11 <- lm(net_strength ~ Rank_binary * Group * Sex + Age, data = rankBin); summary(m11)
Anova(m11, type = 3)
pairs(emmeans(m11, ~ Rank_binary | Sex))

# ==================================================
# Export: final mouse summary
# ==================================================
mouse_summary_full_clean <- mouse_summary_full %>%
  select(Cohort, MouseID, Treatment, Group, Age, Sex, 
         numBouts, TotDur, meanBout_dur, medBout_dur, maxBout_dur, DS_norm, out_strength, in_strength, net_strength, pagerank, hierarchy_type)

mouse_summary_full_clean %>%
  rename_with(~ paste0(.x, "_WS"), .cols = 7:17) %>%
  write_csv("warm_spot_clean.csv")
