# ============================================================
# Generalist vs. specialist simulation
# ============================================================
library(Rmisc); library(tidyverse); library(EloRating); library(ggtern); library(qvalue)
library(car); library(emmeans); library(MCMCpack); library(ggridges)

setwd("~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data")
plotdir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/comp"

theme_nolabels <- theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        strip.text = element_blank(), legend.position = "none")

# ---- Load data ----
tern_df      <- read.csv("tern_df.csv")
WS_allTrials <- read.csv("WS_allTrials.csv")
VS_allTrials <- read.csv("VS_allTrials.csv")
tt_elo       <- read.csv("tt_elo.csv") %>% filter(Cohort != "C32")

# ============================================================
# 1. Estimate common latent dominance (g)
# ============================================================

base <- tern_df %>%
  mutate(MouseID = as.character(MouseID),
         g = rowMeans(cbind(z_DS_norm_TT, z_TotDur_WS, z_totVol_overall_VS), na.rm = TRUE),
         g = as.numeric(scale(g)))

# ============================================================
# 2. Estimate assay-specific sampling variation
# ============================================================

# ---- Warm Spot (WS): average observers within trial ----
ws_trial <- WS_allTrials %>%
  mutate(MouseID = as.character(MouseID)) %>%
  group_by(Cohort, MouseID, Trial) %>%
  dplyr::summarise(WS = mean(TotDur, na.rm = TRUE) / 10, .groups = "drop") %>%
  left_join(base %>% dplyr::select(MouseID, g), by = "MouseID")

ws_mouse <- ws_trial %>% group_by(MouseID) %>%
  dplyr::summarise(WS_mean = mean(WS, na.rm = TRUE), g = first(g), .groups = "drop")
ws_fit <- lm(WS_mean ~ g, data = ws_mouse)

sigma_WS <- ws_trial %>% group_by(MouseID) %>%
  mutate(resid = WS - mean(WS, na.rm = TRUE)) %>% ungroup() %>%
  dplyr::summarise(sigma = sd(resid, na.rm = TRUE)) %>% pull(sigma)

# ---- Void Spot (VS): model positive volumes on log scale ----
vs_trial <- VS_allTrials %>%
  mutate(MouseID = as.character(MouseID), nonzero = totVol_overall > 0) %>%
  left_join(base %>% dplyr::select(MouseID, g), by = "MouseID")

vs_pos <- vs_trial %>% filter(nonzero) %>% mutate(logVS = log(totVol_overall))
vs_mouse <- vs_pos %>% group_by(MouseID) %>%
  dplyr::summarise(logVS_mean = mean(logVS, na.rm = TRUE), g = first(g), .groups = "drop")
vs_fit <- lm(logVS_mean ~ g, data = vs_mouse)

sigma_VS <- vs_pos %>% group_by(MouseID) %>%
  mutate(resid = logVS - mean(logVS, na.rm = TRUE)) %>% ungroup() %>%
  dplyr::summarise(sigma = sd(resid, na.rm = TRUE)) %>% pull(sigma)
p_zero_VS <- mean(!vs_trial$nonzero)

# ---- Tube Test (TT): estimate effect of latent dominance on contest outcome ----
tt_sim_base <- tt_elo %>%
  mutate(winner = as.character(winner), loser = as.character(loser)) %>%
  left_join(base %>% dplyr::select(MouseID, g), by = c("winner" = "MouseID")) %>% dplyr::rename(g_winner = g) %>%
  left_join(base %>% dplyr::select(MouseID, g), by = c("loser" = "MouseID")) %>% dplyr::rename(g_loser = g) %>%
  drop_na(g_winner, g_loser)

tt_bt <- bind_rows(
  tt_sim_base %>% transmute(y = 1, dg = g_winner - g_loser),
  tt_sim_base %>% transmute(y = 0, dg = g_loser - g_winner)
)
tt_fit <- glm(y ~ 0 + dg, family = binomial, data = tt_bt)
beta_TT <- coef(tt_fit)[["dg"]]

# ============================================================
# 3. Simulation functions
# ============================================================

simulate_TT <- function(tt, base, beta_TT) {
  d <- tt %>%
    mutate(winner = as.character(winner), loser = as.character(loser)) %>%
    left_join(base %>% dplyr::select(MouseID, g), by = c("winner" = "MouseID")) %>% dplyr::rename(g1 = g) %>%
    left_join(base %>% dplyr::select(MouseID, g), by = c("loser" = "MouseID")) %>% dplyr::rename(g2 = g) %>%
    drop_na(g1, g2)
  animal1_wins <- rbinom(nrow(d), 1, plogis(beta_TT * (d$g1 - d$g2)))
  d %>%
    mutate(winner_sim = if_else(animal1_wins == 1, winner, loser),
           loser_sim = if_else(animal1_wins == 1, loser, winner)) %>%
    group_by(Cohort) %>%
    dplyr::summarise(mat = list(EloRating::creatematrix(winners = winner_sim, losers = loser_sim)), .groups = "drop") %>%
    mutate(ds = map(mat, ~ EloRating::DS(.x))) %>% dplyr::select(Cohort, ds) %>% unnest(ds) %>%
    transmute(Cohort, MouseID = as.character(ID), DS_norm_TT = normDS)
}

simulate_WS <- function(ws_trial, ws_fit, sigma_WS) {
  d <- ws_trial
  d$mu <- predict(ws_fit, newdata = d)
  d$WS_sim <- pmax(0, rnorm(nrow(d), d$mu, sigma_WS))
  d %>% group_by(Cohort, MouseID) %>%
    dplyr::summarise(TotDur_WS = mean(WS_sim), .groups = "drop")
}

simulate_VS <- function(vs_trial, vs_fit, sigma_VS, p_zero_VS) {
  d <- vs_trial
  d$mu_log <- predict(vs_fit, newdata = d)
  d$zero_sim <- rbinom(nrow(d), 1, p_zero_VS)
  d$VS_sim <- if_else(d$zero_sim == 1, 0, exp(rnorm(nrow(d), d$mu_log, sigma_VS)))
  d %>% group_by(Cohort, MouseID) %>%
    dplyr::summarise(totVol_overall_VS = if_else(any(VS_sim > 0),
                                                 mean(VS_sim[VS_sim > 0]), 0), .groups = "drop")
}

min_max_norm <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

make_ternary <- function(d) {
  d %>%
    mutate(TT_mm = min_max_norm(DS_norm_TT), WS_mm = min_max_norm(TotDur_WS),
           VS_mm = min_max_norm(totVol_overall_VS), tot_dom = TT_mm + WS_mm + VS_mm,
           TT_dom = TT_mm/tot_dom, WS_dom = WS_mm/tot_dom, VS_dom = VS_mm/tot_dom,
           R = 1 - 3*pmin(TT_dom, WS_dom, VS_dom)) %>%
    filter(is.finite(TT_dom), is.finite(WS_dom), is.finite(VS_dom))
}

simulate_generalists <- function() {
  simulate_TT(tt_elo, base, beta_TT) %>%
    full_join(simulate_WS(ws_trial, ws_fit, sigma_WS), by = c("Cohort", "MouseID")) %>%
    full_join(simulate_VS(vs_trial, vs_fit, sigma_VS, p_zero_VS), by = c("Cohort", "MouseID")) %>%
    drop_na(DS_norm_TT, TotDur_WS, totVol_overall_VS) %>% make_ternary()
}

# ============================================================
# 4. Simulate generalist null distribution
# ============================================================

set.seed(123); nsim <- 10000
null_results <- map_dfr(seq_len(nsim), \(i) simulate_generalists() %>% transmute(sim = i, MouseID, R))

quantile(null_results$R, c(.50, .90, .95, .975, .99), na.rm = TRUE)
R_cutoff <- quantile(null_results$R, .90, na.rm = TRUE)

# ============================================================
# 5. Empirical specialization magnitude and direction
# ============================================================

# R measures specialization magnitude; archetypes are used only to label direction.
archetypes <- tribble(
  ~strategy, ~TT, ~WS, ~VS,
  "TT", 1, 0, 0, "WS", 0, 1, 0, "VS", 0, 0, 1,
  "TT+WS", .5, .5, 0, "TT+VS", .5, 0, .5, "WS+VS", 0, .5, .5
)

get_strategy <- function(dat) {
  crossing(dat, archetypes) %>%
    mutate(dx = TT - 1/3, dy = WS - 1/3, dz = VS - 1/3,
           lambda = pmax(0, pmin(1, ((TT_dom - 1/3)*dx + (WS_dom - 1/3)*dy + (VS_dom - 1/3)*dz) /
                                   (dx^2 + dy^2 + dz^2))),
           error = sqrt((TT_dom - (1/3 + lambda*dx))^2 +
                        (WS_dom - (1/3 + lambda*dy))^2 +
                        (VS_dom - (1/3 + lambda*dz))^2)) %>%
    group_by(MouseID) %>% slice_min(error, n = 1, with_ties = FALSE) %>% ungroup() %>%
    dplyr::select(MouseID, strategy)
}

tern_df <- tern_df %>%
  mutate(MouseID = as.character(MouseID), R = 1 - 3*pmin(TT_dom, WS_dom, VS_dom)) %>%
  left_join(get_strategy(.), by = "MouseID") %>%
  mutate(evidence_global = if_else(R > R_cutoff, "Specialist", "Not distinguishable from generalist"))

summary(tern_df$R)
quantile(tern_df$R, c(.50, .75, .90, .95, 1), na.rm = TRUE)
mean(tern_df$R > R_cutoff, na.rm = TRUE)

# ============================================================
# 6. Mouse-specific null inference
# ============================================================

empirical_p <- tern_df %>% dplyr::select(MouseID, R) %>% dplyr::rename(R_emp = R) %>%
  left_join(null_results, by = "MouseID") %>% group_by(MouseID, R_emp) %>%
  dplyr::summarise(p_null = (sum(R >= R_emp, na.rm = TRUE) + 1) / (sum(!is.na(R)) + 1), .groups = "drop")

mouse_cutoffs <- null_results %>% group_by(MouseID) %>%
  dplyr::summarise(R50 = quantile(R, .50, na.rm = TRUE), 
                   R90 = quantile(R, .90, na.rm = TRUE),
                   R95 = quantile(R, .95, na.rm = TRUE), .groups = "drop")

tern_df <- tern_df %>%
  left_join(empirical_p, by = "MouseID") %>% left_join(mouse_cutoffs, by = "MouseID") %>%
  mutate(specialist_sig = R > R90,
         specialization = if_else(specialist_sig, "Specialist", "Not distinguishable from generalist"))

summary(tern_df$p_null)
tern_df %>% count(specialist_sig) %>% mutate(proportion = n/sum(n))

# ============================================================
# 7. Population proportion of generalists vs. specialists
# ============================================================

qobj <- qvalue(tern_df$p_null)
pi_generalist <- qobj$pi0
pi_specialist <- 1 - pi_generalist
pi_generalist; pi_specialist

# ============================================================
# 8. Bootstrap uncertainty in specialist proportion
# ============================================================

set.seed(123); B <- 5000
boot_pi0 <- replicate(B, {
  p_boot <- sample(tern_df$p_null, nrow(tern_df), replace = TRUE)
  tryCatch(qvalue(p_boot)$pi0, error = function(e) NA_real_)
})
boot_specialist <- 1 - boot_pi0
quantile(boot_specialist, c(.025, .50, .975), na.rm = TRUE)
mean(boot_specialist, na.rm = TRUE); sd(boot_specialist, na.rm = TRUE)

# ============================================================
# 9. Test whether specialists exceed chance expectation
# ============================================================

# Each simulated mouse is classified against its own 90th-percentile null threshold.
null_class <- null_results %>%
  left_join(mouse_cutoffs %>% dplyr::select(MouseID, R90), by = "MouseID") %>%
  mutate(specialist = R > R90)

null_prop <- null_class %>% group_by(sim) %>%
  dplyr::summarise(prop_specialist = mean(specialist, na.rm = TRUE), .groups = "drop")

emp_prop <- mean(tern_df$specialist_sig, na.rm = TRUE)
p_specialist <- (sum(null_prop$prop_specialist >= emp_prop) + 1) / (nrow(null_prop) + 1)

emp_prop
summary(null_prop$prop_specialist)
quantile(null_prop$prop_specialist, c(.025, .50, .975))
p_specialist

# ============================================================
# 10. Group/Sex differences in specialization magnitude
# ============================================================

tern_df <- tern_df %>% mutate(Group = factor(Group), Sex = factor(Sex))
m1 <- lm(R ~ Group * Sex, data = tern_df); summary(m1); Anova(m1, type = 3)
pairs(emmeans(m1, ~ Sex)); pairs(emmeans(m1, ~ Group))
pairs(emmeans(m1, ~ Group | Sex)); pairs(emmeans(m1, ~ Sex | Group))

ggplot(tern_df, aes(x = R, y = Group_Sex, fill = Group_Sex)) +
  geom_density_ridges(scale = 1.4) +
  theme_nolabels
ggsave(file.path(plotdir, "R_Distribution.pdf"), width = 7, height = 3, units = "cm")

summarySE(tern_df, measurevar = "R", na.rm = TRUE)

R <- summarySE(tern_df, groupvars = c("Group", "Sex"), measurevar = "R", na.rm = TRUE) %>%
  mutate(lower = R - ci, upper = R + ci); R

R$Group_Sex <- factor(paste(R$Group, R$Sex),
                            levels = c("C57_MUT M", "C57_MUT F", "C57_WT M", "C57_WT F", "BalbC_WT M", "BalbC_WT F"))

tern_df$Group_Sex <- factor(paste(tern_df$Group, tern_df$Sex),
  levels = c("C57_MUT M", "C57_MUT F", "C57_WT M", "C57_WT F", "BalbC_WT M", "BalbC_WT F"))

ggplot(tern_df, aes(x = R, y = Group_Sex, color = Sex)) +
  geom_point(position = position_jitter(width = 0, height = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = R_se, aes(xmin = lower, xmax = upper),
    orientation = "y", color = "black", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_nolabels
ggsave(file.path(plotdir, "R_Pointrange.pdf"), width = 7, height = 3, units = "cm")

# ============================================================
# 11. Diagnostic plots
# ============================================================

# Null ternary distribution
ggtern(test_null <- simulate_generalists(), aes(x = TT_dom, z = WS_dom, y = VS_dom)) +
  geom_point(alpha = 0.5) + theme_bw()

# Empirical vs. simulated-generalist R
ggplot() +
  geom_density(data = null_results, aes(R), linewidth = 1) +
  geom_density(data = tern_df, aes(R), linewidth = 1, linetype = 2) +
  geom_vline(xintercept = R_cutoff, linetype = 3) + theme_classic() +
  labs(x = "Radial specialization index (R)", y = "Density")

# Mouse-specific null p-values
ggplot(tern_df, aes(p_null)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  scale_x_continuous(breaks = seq(0, 1, .1), limits = c(0, 1)) + theme_classic() +
  labs(x = "Mouse-specific generalist-null p-value", y = "Number of animals")

# ============================================================
# 12. Radial plots: strategy direction x median R
# ============================================================

angles <- tribble(
  ~strategy, ~angle,
  "VS", 0, "WS+VS", 60, "WS", 120, "TT+WS", 180, "TT", 240, "TT+VS", 300
)

radial_R_df <- tern_df %>%
  group_by(Group, Sex, strategy) %>%
  dplyr::summarise(median_R = median(R, na.rm = TRUE), .groups = "drop") %>%
  right_join(expand_grid(distinct(tern_df, Group, Sex), angles), by = c("Group", "Sex", "strategy")) %>%
  mutate(Group = factor(Group, levels = c("BalbC_WT", "C57_WT", "C57_MUT")))

ggplot(radial_R_df %>% filter(Sex == "M"), aes(angle, median_R)) +
  geom_segment(aes(xend = angle, y = 0, yend = median_R), linewidth = 1, lineend = "round",
               arrow = arrow(length = unit(.1, "in"))) +
  facet_grid(~Group) + coord_polar() +
  scale_x_continuous(breaks = angles$angle, labels = angles$strategy, limits = c(0, 360)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(.25, .5, .75, 1)) + theme_minimal()
ggsave(file.path(plotdir, "RadialPlot_M.pdf"), width = 21, height = 9, units = "cm")

ggplot(radial_R_df %>% filter(Sex == "F"), aes(angle, median_R)) +
  geom_segment(aes(xend = angle, y = 0, yend = median_R), linewidth = 1, lineend = "round",
               arrow = arrow(length = unit(.1, "in"))) +
  facet_grid(~Group) + coord_polar() +
  scale_x_continuous(breaks = angles$angle, labels = angles$strategy, limits = c(0, 360)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(.25, .5, .75, 1)) + theme_minimal()
ggsave(file.path(plotdir, "RadialPlot_F.pdf"), width = 21, height = 9, units = "cm")

# ============================================================
# 13. Simulated ternary plots
# ============================================================

set.seed(123)
n <- 200

mu_gen  <- c(1/3, 1/3, 1/3)
mu_TT   <- c(.6, .2, .2)
mu_TTWS <- c(.4, .4, .2)

k_gen  <- 100
k_TT   <- 45
k_TTWS <- 45

sim_df <- bind_rows(
  as_tibble(rdirichlet(n, mu_gen*k_gen), .name_repair = ~c("TT_dom","WS_dom","VS_dom")) %>%
    mutate(strategy = "Generalist"),
  as_tibble(rdirichlet(n, mu_TT*k_TT), .name_repair = ~c("TT_dom","WS_dom","VS_dom")) %>%
    mutate(strategy = "TT"),
  as_tibble(rdirichlet(n, mu_TTWS*k_TTWS), .name_repair = ~c("TT_dom","WS_dom","VS_dom")) %>%
    mutate(strategy = "TT+WS")
) %>%
  mutate(MouseID = paste0("Sim", row_number())) %>%
  dplyr::select(MouseID, TT_dom, WS_dom, VS_dom, strategy)

sim_df <- sim_df %>%
  mutate(tot_dom = TT_dom + WS_dom + VS_dom)

for (i in unique(sim_df$strategy)) {
  
  p <- ggtern(sim_df %>% filter(strategy == i),
              aes(x = TT_dom, z = WS_dom, y = VS_dom)) +
    stat_density_tern(
      geom = "polygon", n = 400, h = 1.3, bins = 12,
      aes(fill = after_stat(level))
    ) +
    scale_fill_viridis_c(option = "plasma", direction = 1, trans = "sqrt") +
    theme_nolabels
  
  print(p)
  
  ggsave(
    filename = file.path(plotdir, paste0("TernaryPlot_", i, ".pdf")),
    plot = p, width = 3, height = 3, units = "cm"
  )
}

