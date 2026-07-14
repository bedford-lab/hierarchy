# Libraries ----
library(tidyverse)
library(minpack.lm)  # nlsLM

# Paths ----
setwd("~/Dropbox (Personal)/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data/")
plotdir <- "~/Dropbox (Personal)/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/weight_curves/"

# ==================================================
# Load + basic checks ----
# ==================================================

wc <- read.csv("SDH_WeightCurves - Weight.csv", stringsAsFactors = TRUE)

n_distinct(wc$Cohort)   # 67
n_distinct(wc$MouseID)  # 201

# P0 reference weights (JAX): C57 ~1.3g; Balb ~1.5g
bw <- wc %>%
  select(MouseID, Cohort, Sex, Treatment, DOB) %>%
  distinct() %>%
  mutate(
    Date = DOB, Age_days = 0,
    Weight = case_when(
      str_detect(Treatment, "C57|Shank") ~ 1.3,
      str_detect(Treatment, "Balb")      ~ 1.5,
      TRUE                               ~ NA_real_
    )
  )

# ==================================================
# Exploratory plots ----
# ==================================================

ggplot(filter(wc, str_detect(Treatment, "10W")), aes(Age_days, Weight, group = MouseID)) +
  geom_line() +
  geom_vline(xintercept = 42, linetype = "dashed", color = "red") +   # early weight
  geom_vline(xintercept = 63, linetype = "dashed", color = "blue") +  # late weight
  facet_wrap(vars(Treatment, Sex, Cohort))

ggplot(filter(wc, str_detect(Treatment, "20W")), aes(Age_days, Weight, group = MouseID)) +
  geom_line() +
  geom_vline(xintercept = 42,  linetype = "dashed", color = "red") +  # early weight
  geom_vline(xintercept = 105, linetype = "dashed", color = "blue") + # late weight
  facet_wrap(vars(Treatment, Sex, Cohort))

# ==================================================
# Logistic fits per mouse (nlsLM) + parametric CI simulation ----
# ==================================================

df <- wc %>% mutate(MouseID = as.character(MouseID))

target_ages <- tibble(Age_days = c(28, 42, 56, 70, 84, 98, 112, 126, 140))
n_sim <- 1000

fit_summary      <- tibble()
all_predictions  <- tibble()
fit_success      <- character()
fit_fail         <- character()

logistic_fun <- function(K, r, t0, age) K / (1 + exp(-r * (age - t0)))

for (mouse in unique(df$MouseID)) {
  df_sub <- filter(df, MouseID == mouse)
  
  if (nrow(df_sub) < 4 || all(df_sub$Weight <= 0, na.rm = TRUE)) {
    fit_fail <- c(fit_fail, mouse)
    next
  }
  
  ok <- try({
    fit <- nlsLM(
      Weight ~ K / (1 + exp(-r * (Age_days - t0))),
      data   = df_sub,
      start  = list(K = max(df_sub$Weight, na.rm = TRUE), r = 0.05, t0 = 50),
      control = nls.lm.control(maxiter = 200)
    )
    
    coefs <- coef(fit)
    fit_success <- c(fit_success, mouse)
    
    fit_summary <- bind_rows(
      fit_summary,
      tibble(MouseID = mouse, K = coefs[["K"]], r = coefs[["r"]], t0 = coefs[["t0"]])
    )
    
    vcov_mat <- tryCatch(vcov(fit), error = function(e) NULL)
    if (is.null(vcov_mat)) stop("vcov failed")
    
    param_sim <- MASS::mvrnorm(n_sim, mu = coefs, Sigma = vcov_mat)
    
    # Predictions at target ages (+ 95% CI from parametric simulation)
    pred_vals <- predict(fit, newdata = target_ages)
    sim_preds <- sapply(
      target_ages$Age_days,
      \(age) logistic_fun(param_sim[, "K"], param_sim[, "r"], param_sim[, "t0"], age)
    )
    
    pred_df <- tibble(
      MouseID    = mouse,
      Age_days   = target_ages$Age_days,
      Predicted  = pred_vals,
      CI_lower   = apply(sim_preds, 2, quantile, probs = 0.025, na.rm = TRUE),
      CI_upper   = apply(sim_preds, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
    
    all_predictions <- bind_rows(all_predictions, pred_df)
    
    # Full curve + ribbon
    plot_ages <- tibble(Age_days = seq(min(df_sub$Age_days), max(df_sub$Age_days), length.out = 200))
    sim_curve <- sapply(
      plot_ages$Age_days,
      \(age) logistic_fun(param_sim[, "K"], param_sim[, "r"], param_sim[, "t0"], age)
    )
    
    plot_df <- tibble(
      Age_days = plot_ages$Age_days,
      Fit      = predict(fit, newdata = plot_ages),
      CI_lower = apply(sim_curve, 2, quantile, probs = 0.025, na.rm = TRUE),
      CI_upper = apply(sim_curve, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
    
    mouse_safe <- gsub("[^A-Za-z0-9]", "_", mouse)
    
    p <- ggplot(df_sub, aes(Age_days, Weight)) +
      geom_point() +
      geom_ribbon(data = plot_df, aes(y = Fit, ymin = CI_lower, ymax = CI_upper),
                  fill = "lightblue", alpha = 0.4) +
      geom_line(data = plot_df, aes(y = Fit), color = "blue", linewidth = 1) +
      ggtitle(paste("Mouse", mouse)) +
      theme_minimal()
    
    ggsave(file.path(plotdir, paste0("logistic_fit_", mouse_safe, ".pdf")),
           plot = p, width = 6, height = 4)
    
    TRUE
  }, silent = TRUE)
  
  if (!isTRUE(ok)) fit_fail <- c(fit_fail, mouse)
}

n_distinct(all_predictions$MouseID)  # 201
table(all_predictions$MouseID)

# ==================================================
# Export model parameters + predictions ----
# ==================================================

write.csv(fit_summary,     "logistic_model_parameters.csv",     row.names = FALSE)
write.csv(all_predictions, "logistic_predictions_with_ci.csv",  row.names = FALSE)

# ==================================================
# Build predicted-weight summary table ----
# ==================================================

all_predictions_clean <- all_predictions %>%
  filter(Predicted > 5) %>%
  select(MouseID:Predicted) %>%
  group_by(MouseID) %>%
  pivot_wider(names_from = Age_days, values_from = Predicted) %>%
  ungroup()

all_predictions_clean <- full_join(
  all_predictions_clean,
  wc %>% select(MouseID:DOB) %>% distinct(),
  by = "MouseID"
)

all_predictions_clean <- all_predictions_clean %>%
  mutate(
    Pred_Sac_Weight = case_when(
      str_detect(Treatment, "10W") ~ `70`,
      str_detect(Treatment, "20W") ~ `140`,
      TRUE ~ NA_real_
    ),
    Pred_Wean_Weight = `28`,
    Pred_Mat_Weight  = `56`
  )

wc_summary <- all_predictions_clean %>%
  select(MouseID, Cohort, Sex, Treatment, Pred_Wean_Weight, Pred_Mat_Weight, Pred_Sac_Weight) %>%
  group_by(Cohort, Sex, Treatment) %>%
  mutate(
    Trio_max_WeanWeight = max(Pred_Wean_Weight, na.rm = TRUE),
    Rel_WeanWeight      = Pred_Wean_Weight / Trio_max_WeanWeight,
    Trio_max_SacWeight  = max(Pred_Sac_Weight,  na.rm = TRUE),
    Rel_SacWeight       = Pred_Sac_Weight  / Trio_max_SacWeight,
    Trio_max_MatWeight  = max(Pred_Mat_Weight, na.rm = TRUE),
    Rel_MatWeight       = Pred_Mat_Weight / Trio_max_MatWeight
  ) %>%
  ungroup()

fo <- full_join(all_predictions_clean, wc_summary, by = c("MouseID","Cohort","Sex","Treatment",
                                                          "Pred_Wean_Weight","Pred_Mat_Weight","Pred_Sac_Weight"))

str(fo)
write.csv(fo, "predicted_weights.csv", row.names = FALSE)

# ---- Make Plots ----
theme_nolabels <- theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        strip.text = element_blank(), legend.position = "none")

wc <- wc %>% mutate(
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

wc_ranked <- wc %>% left_join(all_predictions_clean %>% select(Cohort, MouseID, Pred_Sac_Weight)) %>%
  group_by(Cohort) %>%
  mutate(Rank = factor(dense_rank(desc(Pred_Sac_Weight)), levels = 1:3, ordered = TRUE)) %>%
  ungroup()

wc_ranked_summary <- wc_ranked %>%
  group_by(Group, Sex, Rank, Age_days) %>%
  summarise(mean_wt = mean(Weight, na.rm = TRUE),
            sd_wt = sd(Weight, na.rm = TRUE),
            n_cohorts = n(),
            se_wt = sd_wt / sqrt(n_cohorts),
            lower = mean_wt - se_wt, upper = mean_wt + se_wt,
            .groups = "drop")

plot_df <- wc_ranked_summary %>%
  filter(Age_days >= 21,
         (Group == "C57_MUT" & Age_days <= 84) | (Group != "C57_MUT" & Age_days <= 154))

ggplot(plot_df, aes(Age_days/7, mean_wt, color = Rank, fill = Rank, linetype = Sex)) +
  facet_grid(~ Group) +
  scale_x_continuous(breaks = seq(3, 21, by = 3)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = TRUE, linewidth = 0.5) +
  theme_nolabels
ggsave(file.path(plotdir, "Summary_WeightCurves.pdf"), width = 9, height = 3, units = "cm")

wc_ranked_clean <- wc_ranked %>% 
  select(MouseID, Cohort, Rank, Group, Sex, Age, Pred_Sac_Weight) %>%
  distinct()

m1 <- lmer(Pred_Sac_Weight ~ Rank * Group * Sex + Age + (1 | Cohort), data = wc_ranked_clean); summary(m1)
Anova(m1, type = 3)
pairs(emmeans(m1, ~ Rank | Group | Sex))