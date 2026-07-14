# ---- Libraries ----
suppressPackageStartupMessages({
  library(Rmisc); library(data.table); library(tidyverse); library(lme4); library(lmerTest)
  library(emmeans); library(lubridate); library(broom.mixed); library(purrr); library(brms)
  library(tidybayes); library(kSamples); library(EloRating); library(ggbeeswarm)
})

# ---- Paths + theme ----
setwd("~/Dropbox (Personal)/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data")
plotdir <- "~/Dropbox (Personal)/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/tube_test"

theme_nolabels <- theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        strip.text = element_blank(), legend.position = "none")

# ---- Load ----
tt <- read.csv("SDH_TT.csv", stringsAsFactors = TRUE) # 3828

# ==================================================
# Clean + recode
# ==================================================
tt <- tt %>%
  mutate(
    across(c(MouseID, PartnerID),
           ~ as.integer(gsub("[C\\.]", "", gsub("b", "0", .x))),
           .names = "{.col}_Clean"),
    Age        = factor(case_when(str_detect(Treatment, "10W$") ~ "10W",
                                  str_detect(Treatment, "20W$") ~ "20W",
                                  TRUE ~ NA_character_),
                        levels = c("10W", "20W")),
    Geno       = factor(if_else(str_detect(Treatment, "MUT"), "MUT", "WT"), levels = c("WT", "MUT")),
    Background = factor(if_else(str_detect(Treatment, "BalbC"), "BalbC", "C57"), levels = c("C57", "BalbC")),
    Group      = factor(paste(Background, Geno, sep = "_"),
                        levels = c("BalbC_WT", "C57_WT", "C57_MUT")),
    Trial_perMouse = Trial,
    TRIAL          = paste(Date, TrialStart),
    DateTime       = as.POSIXct(paste(Date, TrialStart), format = "%Y-%m-%d %I:%M:%S %p")
  ) %>%
  arrange(Cohort, DateTime) %>%
  group_by(Cohort) %>% mutate(Trial_perCohort = dense_rank(TRIAL)) %>% ungroup()

tt %>% group_by(Group, Sex, Age) %>% summarise(n_cohorts = n_distinct(Cohort), .groups = "drop")

# ==================================================
# Pair-level dataset (1 row per pair per trial)
# ==================================================
pairs <- tt %>%
  mutate(Pair = if_else(MouseID_Clean > PartnerID_Clean,
                        paste(MouseID_Clean, PartnerID_Clean, sep = "_"),
                        paste(PartnerID_Clean, MouseID_Clean, sep = "_"))) %>%
  select(-c(MouseID, MouseID_Clean, PartnerID, PartnerID_Clean, Trial, Trial_perMouse, Side, WL,
            SingleHouse_start, TrialStart, TrialEnd)) %>%
  distinct() %>% droplevels() # 1914

# ==================================================
# Duration (pair level)
# ==================================================
dur_se <- pairs %>%
  group_by(Group, Sex) %>%
  summarise(n = sum(!is.na(Duration)),
            median = median(Duration, na.rm = TRUE),
            Q1 = quantile(Duration, 0.25, na.rm = TRUE),
            Q3 = quantile(Duration, 0.75, na.rm = TRUE),
            IQR = IQR(Duration, na.rm = TRUE),
            .groups = "drop")

ggplot(pairs, aes(Sex, Duration, color = Sex)) +
  facet_wrap(~ Group) +
  geom_beeswarm(method = "hex", size = 0.5, alpha = 0.5) +
  geom_pointrange(data = dur_se, aes(y = median, ymin = Q3, ymax = Q1),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  scale_y_continuous(trans = "log10") +
  theme_nolabels
ggsave(file.path(plotdir, "TT_Duration.pdf"), width = 7, height = 3.5, units = "cm")

m1 <- lmer(Duration ~ Group * Sex + Age + (1 | Pair), data = pairs); summary(m1)
pairs(emmeans(m1, ~ Group))
pairs(emmeans(m1, ~ Sex | Group)) # NS

m1a <- lmer(Duration ~ Sex * Age + (1 | Pair), data = pairs %>% filter(Group == "BalbC_WT")); summary(m1a) # NS
m1b <- lmer(Duration ~ Sex * Age + (1 | Pair), data = pairs %>% filter(Group == "C57_WT")); summary(m1b)
m1c <- lmer(Duration ~ Sex + (1 | Pair), data = pairs %>% filter(Group == "C57_MUT")); summary(m1c)        # NS

pairs(emmeans(m1b, ~ Sex)) # NS
pairs(emmeans(m1b, ~ Sex | Age))
pairs(emmeans(m1b, ~ Age | Sex))

summarySE(pairs, measurevar = "Duration", groupvars = c("Group", "Sex", "Age"), na.rm = TRUE)

pairs %>% group_by(Group) %>%
  summarise(n = sum(!is.na(Duration)),
            median = median(Duration, na.rm = TRUE),
            Q1 = quantile(Duration, 0.25, na.rm = TRUE),
            Q3 = quantile(Duration, 0.75, na.rm = TRUE),
            IQR = IQR(Duration, na.rm = TRUE),
            .groups = "drop")

# ---- Add Trial number per Pair ----
pairs <- pairs %>%
  arrange(Group, Cohort, Pair, DateTime) %>%
  group_by(Pair) %>%
  mutate(Outcome_num = if_else(Outcome == "LR", 1, 0),
         Trial_perPair = dense_rank(TRIAL)) %>%
  ungroup()

# ==================================================
# Outcome (LR vs WA) (pair level)
# ==================================================
outcome <- pairs %>%
  count(Cohort, Group, Age, Sex, Pair, Outcome, name = "n") %>%
  pivot_wider(names_from = Outcome, values_from = n) %>%
  mutate(across(c(WA, LR, `NA`), ~ replace_na(.x, 0L)),
         n_trials = LR + WA + `NA`,
         por_WA = WA / (n_trials))

outcome_filt <- outcome %>% filter(n_trials > 2) %>% droplevels()

por_se <- summarySE(outcome_filt, measurevar = "por_WA", groupvars = c("Group", "Sex"), na.rm = TRUE)

ggplot(outcome_filt, aes(Sex, por_WA, color = Sex)) +
  facet_wrap(~ Group) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = por_se, aes(ymin = por_WA - se, ymax = por_WA + se),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  theme_nolabels
ggsave(file.path(plotdir, "TT_Outcome.pdf"), width = 3.5, height = 3.5, units = "cm")

m2 <- lm(por_WA ~ Group * Sex + Age, data = outcome_filt); summary(m2)
pairs(emmeans(m2, ~ Group))
pairs(emmeans(m2, ~ Sex | Group)) # NS

m2a <- lm(por_WA ~ Sex * Age, data = outcome_filt %>% filter(Group == "BalbC_WT")); summary(m2a) # NS
m2b <- lm(por_WA ~ Sex * Age, data = outcome_filt %>% filter(Group == "C57_WT")); summary(m2b)   # NS
m2c <- lm(por_WA ~ Sex, data = outcome_filt %>% filter(Group == "C57_MUT")); summary(m2c)        # NS

summarySE(outcome_filt, measurevar = "por_WA", groupvars = c("Group", "Sex", "Age"), na.rm = TRUE)

# # ==================================================
# # Winner & loser effects (Bayesian; brms)
# # ==================================================
# tt_streak <- tt %>%
#   arrange(MouseID, Trial_perMouse) %>%
#   group_by(MouseID) %>%
#   mutate(
#     WL_prev = lag(WL),
#     win_streak  = accumulate(WL_prev, \(p, x) if (is.na(x)) p else if (x == 1) p + 1 else 0, .init = 0)[-1],
#     lose_streak = accumulate(WL_prev, \(p, x) if (is.na(x)) p else if (x == 0) p + 1 else 0, .init = 0)[-1]
#   ) %>% ungroup()
# 
# priors <- c(prior(normal(0, 1.5), class = "Intercept"),
#             prior(normal(0, 1.0), class = "b"),
#             prior(exponential(1), class = "sd"))
# 
# tt_win  <- tt_streak %>% filter(!is.na(WL), lose_streak == 0, !(win_streak == 0 & lose_streak == 0))
# tt_lose <- tt_streak %>% filter(!is.na(WL), win_streak  == 0, !(win_streak == 0 & lose_streak == 0))
# 
# win_model <- brm(WL ~ win_streak * Group * Sex + (1 | MouseID) + (1 | PartnerID),
#                  data = tt_win, family = bernoulli(), prior = priors,
#                  chains = 4, iter = 4000, warmup = 1000,
#                  control = list(adapt_delta = 0.99, max_treedepth = 12))
# 
# lose_model <- brm(WL ~ lose_streak * Group * Sex + (1 | MouseID) + (1 | PartnerID),
#                   data = tt_lose, family = bernoulli(), prior = priors,
#                   chains = 4, iter = 4000, warmup = 1000,
#                   control = list(adapt_delta = 0.99, max_treedepth = 12))
# 
# summary(win_model); summary(lose_model)
# 
# ce_win  <- conditional_effects(win_model,  "win_streak:Sex",
#                                conditions = list(Group = c("BalbC_WT", "C57_WT", "C57_MUT")),
#                                resolution = 100, plot = FALSE)[[1]]
# ce_lose <- conditional_effects(lose_model, "lose_streak:Sex",
#                                conditions = list(Group = c("BalbC_WT", "C57_WT", "C57_MUT")),
#                                resolution = 100, plot = FALSE)[[1]]
# 
# plot_ce <- function(ce, xvar, ylab_txt) {
#   ggplot(ce, aes(.data[[xvar]], estimate__, color = Sex, fill = Sex)) +
#     geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.25, color = NA) +
#     geom_line(linewidth = 1) +
#     facet_wrap(~ Group, nrow = 1) +
#     coord_cartesian(ylim = c(0, 1)) +
#     theme_classic(base_size = 14) +
#     labs(x = str_replace_all(xvar, "_", " "), y = ylab_txt)
# }
# 
# plot_ce(ce_win, "win_streak", "Predicted P(win)")
# plot_ce(ce_lose, "lose_streak", "Predicted P(lose)")
# 
# ggplot(ce_win, aes(win_streak, estimate__, color = Sex, fill = Sex)) +
#   geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.25, color = NA) +
#   geom_line(linewidth = 0.5) +
#   facet_wrap(~ Group) +
#   coord_cartesian(ylim = c(0, 1)) +
#   theme_nolabels
# ggsave(file.path(plotdir, "TT_WinStreak.pdf"), width = 7, height = 3, units = "cm")
# 
# # ---- Compare winner effect magnitude ----
# emm_win  <- emtrends(win_model,  specs = c("Group", "Sex"), var = "win_streak"); summary(emm_win)
# emm_lose <- emtrends(lose_model, specs = c("Group", "Sex"), var = "lose_streak"); summary(emm_lose)
# 
# win_df  <- as.data.frame(emm_win)  %>% mutate(effect = "Winner effect (win_streak)", estimate = win_streak.trend)
# lose_df <- as.data.frame(emm_lose) %>% mutate(effect = "Loser effect (lose_streak)", estimate = lose_streak.trend)
# 
# forest_df <- bind_rows(win_df, lose_df) %>%
#   mutate(GroupSex = factor(paste(Group, Sex),
#                            levels = c("BalbC_WT F","BalbC_WT M","C57_WT F","C57_WT M","C57_MUT F","C57_MUT M")))
# 
# ggplot(forest_df, aes(estimate, GroupSex, color = Sex)) +
#   geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
#   geom_point(size = 2) +
#   geom_errorbarh(aes(xmin = lower.HPD, xmax = upper.HPD), height = 0.15, linewidth = 0.7) +
#   facet_wrap(~ effect, scales = "free_x") +
#   theme_bw(base_size = 13) +
#   theme(panel.grid.major.y = element_blank(), legend.position = "bottom") +
#   labs(x = "Slope (log-odds per unit streak)", y = NULL,
#        title = "Winner and loser streak effects by Group and Sex",
#        subtitle = "Points = posterior medians; bars = 95% HPD intervals")
# 
# # ==================================================
# # Max streaks (20-trial mice) + AD tests
# # ==================================================
# max_streak_df <- tt_streak %>%
#   group_by(Group, Age, Sex, Cohort, MouseID) %>%
#   summarise(n_trials = max(Trial_perMouse, na.rm = TRUE),
#             max_win_streak  = max(win_streak,  na.rm = TRUE),
#             max_lose_streak = max(lose_streak, na.rm = TRUE),
#             trial_of_max_win  = Trial_perMouse[which(win_streak  == max(win_streak,  na.rm = TRUE))][1],
#             trial_of_max_lose = Trial_perMouse[which(lose_streak == max(lose_streak, na.rm = TRUE))][1],
#             .groups = "drop")
# 
# df20 <- filter(max_streak_df, n_trials == 20)
# 
# plot_den <- function(dat, x, title_txt) {
#   ggplot(dat, aes(.data[[x]], color = Group, fill = Group)) +
#     geom_density(alpha = 0.2, adjust = 1.2) +
#     theme_bw(base_size = 14) +
#     labs(x = x, y = "Density", title = title_txt)
# }
# 
# plot_den(filter(df20, Sex == "F"), "max_win_streak",  "Max win streak (Females; 20 trials)")
# plot_den(filter(df20, Sex == "M"), "max_win_streak",  "Max win streak (Males; 20 trials)")
# plot_den(filter(df20, Sex == "F"), "max_lose_streak", "Max lose streak (Females; 20 trials)")
# plot_den(filter(df20, Sex == "M"), "max_lose_streak", "Max lose streak (Males; 20 trials)")
# 
# df20_F <- filter(df20, Sex == "F"); df20_M <- filter(df20, Sex == "M")
# 
# ad.test(df20_F$max_win_streak[df20_F$Group == "BalbC_WT"], df20_F$max_win_streak[df20_F$Group == "C57_WT"])    # NS
# ad.test(df20_F$max_win_streak[df20_F$Group == "C57_WT"],  df20_F$max_win_streak[df20_F$Group == "C57_MUT"])
# ad.test(df20_F$max_win_streak[df20_F$Group == "BalbC_WT"], df20_F$max_win_streak[df20_F$Group == "C57_MUT"])   # NS
# 
# ad.test(df20_M$max_win_streak[df20_M$Group == "BalbC_WT"], df20_M$max_win_streak[df20_M$Group == "C57_WT"])    # NS
# ad.test(df20_M$max_win_streak[df20_M$Group == "C57_WT"],  df20_M$max_win_streak[df20_M$Group == "C57_MUT"])    # NS
# ad.test(df20_M$max_win_streak[df20_M$Group == "BalbC_WT"], df20_M$max_win_streak[df20_M$Group == "C57_MUT"])   # NS
# 
# ad.test(df20_F$max_lose_streak[df20_F$Group == "BalbC_WT"], df20_M$max_lose_streak[df20_M$Group == "BalbC_WT"]) # NS
# ad.test(df20_F$max_lose_streak[df20_F$Group == "C57_WT"], df20_M$max_lose_streak[df20_M$Group == "C57_WT"])     # NS
# ad.test(df20_F$max_lose_streak[df20_F$Group == "C57_MUT"], df20_M$max_lose_streak[df20_M$Group == "C57_MUT"])   # NS
# 
# ggplot(df20, aes(max_win_streak, color = Sex, fill = Sex)) +
#   geom_density(alpha = 0.2) +
#   facet_wrap(~ Group) +
#   theme_nolabels
# ggsave(file.path(plotdir, "TT_WinStreak_Density.pdf"), width = 7, height = 2, units = "cm")
# 
# summarySE(df20, measurevar = "max_win_streak", groupvars = c("Group", "Sex"), na.rm = TRUE)

# ==================================================
# Win summary (per mouse; across days)
# ==================================================
win_rate <- tt %>%
  group_by(Cohort, Treatment, Group, Age, Sex, MouseID) %>%
  summarise(wins = sum(WL == 1, na.rm = TRUE),
            losses = sum(WL == 0, na.rm = TRUE),
            missing = sum(is.na(WL)),
            nTrials_perMouse = n(),
            .groups = "drop") %>%
  mutate(valid_trials = nTrials_perMouse - missing,
         win_prop = if_else(valid_trials > 0, wins / valid_trials, NA_real_))

tt_sum <- tt %>%
  group_by(Cohort, Date, Treatment, Group, Age, Sex, MouseID) %>%
  summarise(wins = sum(WL, na.rm = TRUE), .groups = "drop") %>%
  group_by(Cohort, Date) %>% mutate(dailyRank = rank(-wins, ties.method = "average")) %>% ungroup()

rank_summary <- tt_sum %>% group_by(Cohort, MouseID) %>% summarise(totWins = sum(wins), .groups = "drop")
mouse_levels <- rank_summary %>% arrange(Cohort, desc(totWins)) %>% pull(MouseID) %>% as.character() %>% unique()

tt_sum <- tt_sum %>%
  mutate(MouseID = factor(MouseID, levels = mouse_levels)) %>%
  group_by(Cohort, MouseID) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(Day = row_number(), mean_dailyRank = mean(dailyRank, na.rm = TRUE), nDays = n()) %>%
  ungroup()

ggplot(tt_sum, aes(Day, dailyRank, group = MouseID, color = Group)) +
  facet_wrap(~ Cohort) +
  geom_line(linewidth = 1) +
  scale_y_reverse(breaks = 1:3, labels = as.character(1:3)) +
  labs(x = "Day", y = "Daily rank") +
  theme_classic()

mouse_summary <- tt_sum %>%
  distinct(Cohort, MouseID, Treatment, Group, Age, Sex, mean_dailyRank, nDays) %>%
  left_join(rank_summary, by = c("Cohort", "MouseID")) %>%
  left_join(win_rate, by = c("Cohort", "Treatment", "Group", "Age", "Sex", "MouseID")) %>%
  group_by(Cohort) %>% mutate(win_prop_Rank = rank(-win_prop, ties.method = "average")) %>% ungroup()

# ==================================================
# Elo rating: data prep
# ==================================================
tt_elo <- tt %>%
  mutate(Date = mdy(as.character(Date)),
         DateTime = parse_date_time(paste(Date, as.character(TrialStart)),
                                    orders = c("Y-m-d I:M:S p", "Y-m-d H:M:S"),
                                    tz = "America/Denver")) %>%
  filter(Side == "L") %>%
  mutate(
    winner = case_when(WL == 1 ~ as.character(MouseID),
                       WL == 0 ~ as.character(PartnerID),
                       TRUE ~ as.character(MouseID)),
    loser  = case_when(WL == 1 ~ as.character(PartnerID),
                       WL == 0 ~ as.character(MouseID),
                       TRUE ~ as.character(PartnerID)),
    draw = is.na(WL)
  ) %>%
  arrange(Cohort, Date, SingleHouse_start, TrialStart) %>%
  select(Cohort, Sex, Treatment, MouseID, PartnerID, Date, TrialStart, DateTime, Outcome, winner, loser, draw) %>%
  group_by(Cohort) %>%
  mutate(FakeTrial_Date = seq.Date(as.Date("2000-01-01"), by = "day", length.out = n())) %>%
  arrange(FakeTrial_Date) %>% ungroup()

seqcheck(winner = tt_elo$winner, loser = tt_elo$loser, Date = tt_elo$FakeTrial_Date)

# ==================================================
# Elo + ISI + stability per cohort
# ==================================================
threshold <- 0.85
Cohorts <- unique(as.character(tt_elo$Cohort))

stability_over_time <- function(elo_obj, threshold = threshold) {
  dates <- elo_obj$stability$date
  traj <- tibble(
    FakeTrial_Date = dates,
    TrialWithinCohort = seq_along(dates),
    S = map_dbl(dates, ~ stab_elo(elo_obj, to = as.character(.x)))
  )
  idx <- which(traj$S >= threshold)
  stable_from_trial <- if (length(idx)) min(idx) else NA_integer_
  n_trials <- nrow(traj)
  
  list(
    traj = traj,
    summary = tibble(
      threshold = threshold,
      stable_from_trial = stable_from_trial,
      stable_from_date  = if (!is.na(stable_from_trial)) traj$FakeTrial_Date[stable_from_trial] else as.Date(NA),
      S_at_stable = if (!is.na(stable_from_trial)) traj$S[stable_from_trial] else NA_real_,
      S_final = last(traj$S),
      n_trials = n_trials,
      S_at_trial10 = if (n_trials > 6) traj$S[9] else NA_real_,
      S_at_trial20 = if (n_trials > 6) traj$S[19] else NA_real_,
      S_at_trial30 = if (n_trials > 6) traj$S[29] else NA_real_
    )
  )
}

elo_perCohort <- function(cohort_id, threshold = threshold) {
  fo <- filter(tt_elo, Cohort == cohort_id)
  
  res <- elo.seq(winner = fo$winner, loser = fo$loser, Date = fo$FakeTrial_Date,
                 draw = fo$draw, runcheck = TRUE, startvalue = 1000)
  
  pdf(file.path(plotdir, "elo_plots", paste0("EloPlot_", cohort_id, ".pdf")))
  eloplot(res, color = FALSE); dev.off()
  
  rtbl <- as.data.frame(ISIranks(ISI(creatematrix(res)))); names(rtbl)[1] <- "MouseID"
  elo_df <- as.data.frame(extract_elo(res)) %>% rownames_to_column("MouseID") %>% rename(Elo = 2)
  
  mouse_metrics <- elo_df %>%
    left_join(rtbl, by = "MouseID") %>%
    mutate(Cohort = cohort_id, Stability_lastTrial = stab_elo(res)) %>%
    select(-any_of(c("rnkg2", "rnkg3", "avg")))
  
  stab <- stability_over_time(res, threshold)
  list(mouse = mouse_metrics,
       traj = stab$traj %>% mutate(Cohort = cohort_id, .before = 1),
       summary = stab$summary %>% mutate(Cohort = cohort_id, .before = 1))
}

elo_list <- lapply(Cohorts, elo_perCohort, threshold = threshold)
elo_mouse_all    <- bind_rows(lapply(elo_list, `[[`, "mouse"))
stab_traj_all    <- bind_rows(lapply(elo_list, `[[`, "traj"))
stab_summary_all <- bind_rows(lapply(elo_list, `[[`, "summary"))

# ---- Strict stability (S >= threshold for all remaining trials) ----
strict_stability_forever <- function(df, threshold = threshold) {
  df <- arrange(df, TrialWithinCohort)
  ok_from_here <- rev(cummin(rev(df$S >= threshold))) > 0
  idx <- which(ok_from_here)[1]
  
  if (is.na(idx)) tibble(Cohort = df$Cohort[1], threshold = threshold,
                         stable_from_trial_strict = NA_integer_,
                         stable_from_date_strict = as.Date(NA),
                         S_at_strict_stable = NA_real_)
  else tibble(Cohort = df$Cohort[1], threshold = threshold,
              stable_from_trial_strict = idx,
              stable_from_date_strict = df$FakeTrial_Date[idx],
              S_at_strict_stable = df$S[idx])
}

stab_strict_all <- stab_traj_all %>%
  group_by(Cohort) %>% group_modify(~ strict_stability_forever(.x, threshold)) %>% ungroup()

stab_summary_all <- full_join(stab_summary_all, stab_strict_all, by = c("Cohort", "threshold"))

# ==================================================
# Add nTrials + Elo + stability to mouse_summary
# ==================================================
nTrials_perCohort <- pairs %>%
  group_by(Cohort) %>%
  summarise(nTrials_perCohort = max(Trial_perCohort, na.rm = TRUE), .groups = "drop")

stab_summary_clean <- stab_summary_all %>%
  select(-threshold, -stable_from_trial, -stable_from_date, -S_at_stable, -stable_from_date_strict)

elo_mouse_all <- elo_mouse_all %>% mutate(MouseID = as.character(MouseID))
mouse_summary <- mouse_summary %>% mutate(MouseID = as.character(MouseID))

mouse_summary_full <- mouse_summary %>%
  left_join(nTrials_perCohort, by = "Cohort") %>%
  left_join(elo_mouse_all, by = c("MouseID", "Cohort")) %>%
  left_join(stab_summary_clean, by = "Cohort") %>%
  group_by(Cohort) %>%
  mutate(EloRank_final = rank(-Elo, ties.method = "first"),
         EloRank_final = if_else(duplicated(Elo) | duplicated(Elo, fromLast = TRUE), NA_real_, EloRank_final)) %>%
  ungroup()

# ---- Plot S(t) per cohort ----
ggplot(stab_traj_all, aes(TrialWithinCohort, S)) +
  geom_line() +
  geom_hline(yintercept = threshold, linetype = "dashed") +
  facet_wrap(~ Cohort, scales = "free_x") +
  labs(x = "Trial within cohort", y = "Stability index S",
       title = sprintf("Hierarchy stability over time (threshold = %.2f)", threshold)) +
  theme_classic()

# ==================================================
# Hierarchy metrics (Landau h, David's scores)
# ==================================================
get_mat <- function(df) EloRating::creatematrix(winners = df$winner, losers = df$loser)

metrics <- tt_elo %>%
  group_by(Cohort) %>%
  summarise(mat = list(get_mat(pick(everything()))), .groups = "drop") %>%
  mutate(landau = map(mat, ~ EloRating::h.index(.x, loops = 1000)),
         ds     = map(mat, ~ EloRating::DS(.x)))

cohort_unpacked <- metrics %>%
  transmute(Cohort,
            landau_h      = map_dbl(landau, ~ .x$value[.x$variable == "h index"]),
            landau_hprime = map_dbl(landau, ~ .x$value[.x$variable == "h' index"]),
            landau_p      = map_dbl(landau, ~ .x$value[.x$variable == "p right"]))

metrics_mouse <- metrics %>%
  select(Cohort, ds) %>%
  unnest(ds) %>%
  rename(MouseID = ID, DS_raw = DS, DS_norm = normDS)

mouse_summary_full <- mouse_summary_full %>%
  full_join(left_join(metrics_mouse, cohort_unpacked, by = "Cohort"), by = c("Cohort", "MouseID")) %>%
  select(-landau_h, -landau_hprime, -landau_p)

# ==================================================
# Hierarchy models (cohort level)
# ==================================================
mouse_summary_full_cohort <- mouse_summary_full %>%
  select(-MouseID, -mean_dailyRank, -totWins, -wins, -losses, -missing, -valid_trials,
         -win_prop, -win_prop_Rank, -Elo, -rnkg1, -EloRank_final, -DS_raw, -DS_norm) %>%
  distinct()

mouse_summary_full_cohort_clean <- filter(mouse_summary_full_cohort, nTrials_perCohort > 6)

cohort_info <- mouse_summary_full_cohort %>%
  select(Cohort, Treatment, Group, Age, Sex, nDays, nTrials_perCohort) %>%
  distinct()

stab_traj_all <- full_join(stab_traj_all, cohort_info)

stab_traj_all_sum <- stab_traj_all %>%
  filter(nTrials_perCohort > 6) %>%
  mutate(TrialWithinCohort = TrialWithinCohort + 1) %>%
  group_by(Group, Sex, TrialWithinCohort) %>%
  summarise(mean_stab = mean(S, na.rm = TRUE),
            sd_stab   = sd(S, na.rm = TRUE),
            n_cohorts = n_distinct(Cohort),
            se_stab   = sd_stab / sqrt(n_cohorts),
            lower = mean_stab - se_stab, upper = mean_stab + se_stab,
            .groups = "drop")

ggplot(filter(stab_traj_all_sum, n_cohorts > 2),
       aes(TrialWithinCohort, mean_stab, color = Sex, fill = Sex)) +
  facet_wrap(~ Group) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
  geom_line(linewidth = 0.5) +
  coord_cartesian(ylim = c(0.6, 1.0)) +
  theme_nolabels
ggsave(file.path(plotdir, "TT_Stability.pdf"), width = 9, height = 3, units = "cm")

m3 <- lmer(S ~ TrialWithinCohort * Group * Sex + Age + (1 | Cohort), stab_traj_all); summary(m3)

emm_slope <- emtrends(m3, specs = c("Group", "Sex"), var = "TrialWithinCohort"); summary(emm_slope)
pairs(emm_slope, by = "Sex") 
pairs(emm_slope, by = "Group") # NS

m3a <- lmer(S ~ Sex * Age + (1 | Cohort), stab_traj_all %>% filter(nTrials_perCohort > 6 & Group == "BalbC_WT")); summary(m3a)
m3b <- lmer(S ~ Sex * Age + (1 | Cohort), stab_traj_all %>% filter(nTrials_perCohort > 6 & Group == "C57_WT")); summary(m3b) # NS
m3c <- lmer(S ~ Sex + (1 | Cohort), stab_traj_all %>% filter(nTrials_perCohort > 6 & Group == "C57_MUT")); summary(m3c)      # NS

pairs(emmeans(m3a, ~ Sex))
pairs(emmeans(m3a, ~ Sex | Age))
pairs(emmeans(m3a, ~ Age | Sex))

summarySE(stab_traj_all %>% filter(nTrials_perCohort == 30), 
                                   measurevar = "S", groupvars = c("Group", "Sex", "Age"), na.rm = TRUE)

# ---- LMMs on trial-wise stability S ----
combos <- stab_traj_all %>% distinct(Group, Sex)

fit_one_lmm <- function(g, s) {
  dat <- filter(stab_traj_all, nDays > 2, Group == g, Sex == s)
  tidy(lmer(S ~ TrialWithinCohort + (1 | Cohort), data = dat), effects = "fixed") %>%
    mutate(Group = g, Sex = s, .before = 1)
}

lmm_results <- map2_dfr(combos$Group, combos$Sex, fit_one_lmm)

stab_summary_all <- stab_summary_all %>%
  left_join(cohort_info, by = "Cohort") %>%
  distinct()

# ---- Make plots ----
mouse_summary_full_cohort_clean <- mouse_summary_full_cohort_clean %>%
  mutate(stable_from_trial_strict =
      if_else(is.na(stable_from_trial_strict), 30L, stable_from_trial_strict + 1L))

plot_df <- mouse_summary_full_cohort_clean %>%
  select(Group, Sex, S_at_trial10, S_at_trial20, S_at_trial30) %>%
  pivot_longer(starts_with("S_at_trial"), names_to = "Trial", values_to = "Stability")

plot_summary <- plot_df %>% group_by(Group, Sex, Trial) %>%
  summarise(mean = mean(Stability, na.rm=TRUE), 
            se = sd(Stability, na.rm=TRUE)/sqrt(n()), 
            .groups="drop")

ggplot(plot_summary, aes(Sex, mean, fill = Trial)) +
  facet_wrap(~ Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(.9), width = 0) +
  coord_cartesian(ylim = c(0.7, 1.0)) +
  theme_nolabels
ggsave(file.path(plotdir, "TT_StabilityTrial.pdf"), width = 7, height = 3.5, units = "cm")

# ---- Stability at Trial 10 ----
m4 <- lm(S_at_trial10 ~ Group * Sex + Age, data = mouse_summary_full_cohort_clean); summary(m4)
pairs(emmeans(m4, ~ Group))       # NS
pairs(emmeans(m4, ~ Sex | Group))  
pairs(emmeans(m4, ~ Group | Sex)) # NS

# ---- Stability at Trial 20 ----
m5 <- lm(S_at_trial20 ~ Group * Sex + Age, data = mouse_summary_full_cohort_clean); summary(m5)
pairs(emmeans(m5, ~ Group)) 
pairs(emmeans(m5, ~ Sex | Group)) 
pairs(emmeans(m5, ~ Group | Sex))

# ---- Stability at Trial 30 ----
m6 <- lm(S_at_trial30 ~ Group * Sex + Age, data = mouse_summary_full_cohort_clean); summary(m6)
pairs(emmeans(m6, ~ Group))       # NS
pairs(emmeans(m6, ~ Sex | Group)) # NS  
pairs(emmeans(m6, ~ Group | Sex))

ggplot(mouse_summary_full_cohort_clean, aes(S_at_trial30, color = "green", fill = "green")) +
    geom_density(alpha = 0.2) +
    theme_nolabels
ggsave(file.path(plotdir, "TT_Stab30_Density.pdf"), width = 3, height = 3, units = "cm")

summarySE(mouse_summary_full_cohort_clean, measurevar = "S_at_trial30", na.rm = TRUE)

mouse_summary_full_cohort_clean %>% group_by(Group, Sex, Age) %>%
  summarise(n_cohorts = n_distinct(Cohort), .groups = "drop")

# ==================================================
# Hierarchy classification (DS-based; cohort-level)
# ==================================================
ds_plot <- mouse_summary_full %>%
  filter(!is.na(DS_norm)) %>%
  group_by(Cohort) %>%
  arrange(desc(DS_norm), .by_group = TRUE) %>%
  mutate(rank_DS = row_number()) %>%
  ungroup()

ggplot(ds_plot %>% filter(nTrials_perCohort > 6),
       aes(rank_DS, DS_norm, group = Cohort)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  facet_grid(Sex ~ Group) +
  scale_x_continuous(breaks = 1:3, labels = c("Top", "Mid", "Low")) +
  theme_classic() +
  labs(x = "Within-cohort DS rank", y = "DS_norm")

global_dists <- mouse_summary_full %>%
  select(Cohort, MouseID, DS_norm) %>%
  group_by(Cohort) %>%
  filter(n() == 3) %>%
  summarise(dist = list(as.numeric(dist(DS_norm))), .groups = "drop") %>%
  unnest(dist)

hist(global_dists$dist)

cohort_dists <- mouse_summary_full %>%
  select(Cohort, MouseID, DS_norm) %>%
  group_by(Cohort) %>%
  filter(n() == 3) %>%
  summarise(ds_sorted = list(sort(DS_norm, decreasing = TRUE)), .groups = "drop") %>%
  mutate(d1 = map_dbl(ds_sorted, 1), d2 = map_dbl(ds_sorted, 2), d3 = map_dbl(ds_sorted, 3),
         delta_12 = d1 - d2, delta_23 = d2 - d3, delta_13 = d1 - d3)

cohort_dists <- full_join(cohort_dists, cohort_info)

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
  left_join(select(cohort_classified, Cohort, hierarchy_type), by = "Cohort") %>%
  group_by(Cohort) %>%
  mutate(DS_rank = dense_rank(desc(DS_norm))) %>% ungroup()

mouse_summary_full %>% count(Group, Sex, hierarchy_type) %>%
  group_by(Group, Sex) %>% mutate(total = sum(n)) %>%
  slice_max(n, n = 1, with_ties = TRUE) %>% ungroup()

cohort_level <- mouse_summary_full %>%
  distinct(Cohort, Group, Sex, hierarchy_type, nTrials_perCohort) %>%
  mutate(hierarchy_type = factor(hierarchy_type,
                                 levels = c("despotic", "transitive", "co-dominant", "weak")))

ggplot(cohort_level, aes(Sex, fill = hierarchy_type)) +
  geom_bar(position = "fill") +
  facet_wrap(~ Group) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(y = "Proportion of cohorts", fill = "Hierarchy type")

ggplot(cohort_level, aes(Sex, fill = hierarchy_type)) +
  geom_bar(position = "fill") +
  facet_wrap(~ Group) +
  theme_nolabels
ggsave(file.path(plotdir, "TT_HierarchyType.pdf"), width = 3, height = 2.5, units = "cm")

fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "BalbC_WT"))) # NS
fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "C57_WT")))   # NS
fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "C57_MUT")))  # NS

fisher.test(xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "F"))) # 0.02353
fisher.test(xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "M"))) # NS

groups <- levels(droplevels(cohort_level$Group))
pairs_grp <- combn(groups, 2, simplify = FALSE)

posthoc_fisher_F <- map_dfr(pairs_grp, \(g) {
  tab <- xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "F", Group %in% g))
  ft <- fisher.test(tab)
  tibble(contrast = paste(g, collapse = " vs "), p.value = ft$p.value)
})

# ---- Proportions ----
hierarchy_props <- cohort_level %>% 
  dplyr::count(Group, Sex, hierarchy_type) %>%
  group_by(Group, Sex) %>% mutate(percent = 100 * n / sum(n)) %>% ungroup()

cohort_level %>% mutate(Geno = if_else(grepl("WT", Group), "WT", "MUT")) %>%
  dplyr::count(Geno, Sex, hierarchy_type) %>%
  group_by(Geno, Sex) %>% mutate(percent = 100 * n / sum(n)) %>% ungroup()

# ---- Representative Elo plots for each hierarchy type ----
target_cohorts <- c("C16", "C24", "C33", "C74")

elo_cmat_subset <- lapply(target_cohorts, \(cohort_id)
                          tt_elo %>% filter(Cohort == cohort_id) %>%
                            { elo.seq(.$winner, .$loser, Date = .$FakeTrial_Date, draw = .$draw,
                                      runcheck = TRUE, startvalue = 1000)$cmat } %>%
                            as.data.frame() %>% mutate(Cohort = cohort_id, trial = row_number(), .before = 1) %>%
                            pivot_longer(-c(Cohort, trial), names_to = "MouseID", values_to = "Elo")
)

elo_cmat_all <- bind_rows(elo_cmat_subset)

elo_cmat_all <- elo_cmat_all %>%
  left_join(ds_plot %>% filter(Cohort %in% target_cohorts) %>%
              select(Cohort, MouseID, mean_dailyRank, win_prop, rnkg1, EloRank_final, DS_norm, rank_DS),
            by = c("Cohort", "MouseID")) %>%
  mutate(rank_DS = factor(rank_DS, levels = c(1,2,3)),
         Cohort  = factor(Cohort, levels = c("C24", "C33","C16", "C74")))

ggplot(elo_cmat_all, aes(trial, Elo, color = rank_DS)) +
  facet_wrap(~ Cohort, nrow = 1) +
  geom_line(linewidth = 0.5) +
  theme_nolabels
ggsave(file.path(plotdir, "TT_Elo.pdf"), width = 9, height = 3, units = "cm")

sample_size <- cohort_info %>% group_by(Group, Sex) %>%
  summarise(n_cohorts = n_distinct(Cohort), .groups = "drop")

# ==================================================
# Export: final mouse summary
# ==================================================
mouse_summary_full_clean <- mouse_summary_full %>%
  select(Cohort, MouseID, Treatment, Group, Age, Sex, wins, losses, win_prop, mean_dailyRank, Elo, DS_norm, hierarchy_type)

mouse_summary_full_clean %>%
  rename_with(~ paste0(.x, "_TT"), .cols = 7:13) %>%
  write_csv("tube_test_clean.csv")
