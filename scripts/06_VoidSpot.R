# Libraries ----
library(Rmisc); library(data.table); library(tidyverse); library(lme4); library(lmerTest);
library(emmeans); library(broom.mixed); library(MASS); library(caret); library(kSamples); 
library(pheatmap); library(ggbeeswarm); library(car); library(diptest); library(mclust)

# Paths ----
setwd("~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data")
plotdir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/void_spot"
vsa_dir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data/VSA_outputs"

theme_nolabels <- theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        strip.text = element_blank(), legend.position = "none")

# ==================================================
# LUT ----
# ==================================================
lut <- read.csv("SDH_VSA - LUT.csv", stringsAsFactors = TRUE) %>%
  dplyr::mutate(
    Age        = factor(case_when(str_detect(Treatment, "10W$") ~ "10W",
                                  str_detect(Treatment, "20W$") ~ "20W", TRUE ~ NA_character_),
                        levels = c("10W", "20W")),
    Geno       = factor(if_else(str_detect(Treatment, "MUT"), "MUT", "WT"), levels = c("WT", "MUT")),
    Background = factor(if_else(str_detect(Treatment, "BalbC"), "BalbC", "C57"), levels = c("BalbC", "C57")),
    Group      = factor(paste(Background, Geno, sep = "_"),
                        levels = c("BalbC_WT", "C57_WT", "C57_MUT"))
  )

lut %>% summarise(n_files = n_distinct(file_name),
                  n_cohorts = n_distinct(Cohort),
                  n_mice = n_distinct(MouseID))

lut %>% dplyr::count(Group, name = "n_files") %>%
  left_join(lut %>% group_by(Group) %>% summarise(n_cohorts = n_distinct(Cohort), .groups = "drop"),
            by = "Group")

ctrl_excl_cohorts <- c(
  "C3_MW.R2","C5_MW.R1","C5_MW.R2","C6_MW.R1","C6_MW.R2",
  "C7","C12","C15","C15b","C19","C23","C31","C32","C37","C40")

# ==================================================
# FIJI import + geometry-based edge/middle classifier ----
# ==================================================
read_fiji <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(Spot = 1) %>%
    dplyr::mutate(file_name = basename(path))
}

edge_mid <- function(dat, x_shift = 0.75, y_shift = 1.50) {
  dat <- dat %>% dplyr::mutate(X0 = X - min(X, na.rm = TRUE), Y0 = Y - min(Y, na.rm = TRUE))
  Xmax <- max(dat$X0, na.rm = TRUE); Ymax <- max(dat$Y0, na.rm = TRUE)
  x <- sqrt(((Xmax * Ymax) / 2 * Xmax) / Ymax)
  y <- (Ymax * x) / Xmax
  x_min <- (Xmax - x) / 2; x_max <- Xmax - x_min - x_shift
  y_min <- (Ymax - y) / 2; y_max <- Ymax - y_min - y_shift
  dat %>%
    dplyr::mutate(
      X_pos = if_else(between(X0, x_min, x_max), "middle", "edge"),
      Y_pos = if_else(between(Y0, y_min, y_max), "middle", "edge"),
      pos   = factor(if_else(X_pos == "middle" & Y_pos == "middle", "middle", "edge"),
                     levels = c("edge", "middle"))
    )
}

files <- list.files(vsa_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

df <- map_dfr(files, read_fiji) %>%
  dplyr::mutate(
    across(c(Spot, Area, X, Y), as.numeric),
    TrialType = case_when(str_detect(file_name, "BL|nostim|PBS") ~ "BL", TRUE ~ "VSA"),
    Vol = 12.7722 * Area
  ) %>%
  filter(is.na(Label), Area >= 0.05, !is.na(TrialType)) %>%
  dplyr::select(-Label) %>%
  group_by(TrialType) %>%
  group_modify(~ edge_mid(.x)) %>%
  ungroup()

n_distinct(df$file_name) # 689

# ==================================================
# QC plots: Stim vs No Stim trials ----
# =================================================
plot_zone <- function(type) {
  ggplot(filter(df, TrialType == type), aes(X, Y, color = pos, size = Vol)) +
    geom_point() +
    coord_fixed() +
    scale_color_brewer(palette = "Dark2") +
    scale_size(range = c(0, 10)) +
    theme_classic() +
    labs(title = paste(type, "Trials"))
}

plot_zone("BL")
plot_zone("VSA")

# ==================================================
# Per-paper summaries (voids ≥ 1 uL) ----
# =================================================
df1 <- df %>% filter(Vol >= 1)

sum_overall <- df1 %>%
  group_by(file_name) %>%
  summarise(
    numVoids_overall = n(), minVol_overall = min(Vol), maxVol_overall = max(Vol),
    meanVol_overall = mean(Vol), medVol_overall = median(Vol), totVol_overall = sum(Vol),
    .groups = "drop"
  )

sum_by_pos <- df1 %>%
  group_by(file_name, pos) %>%
  summarise(
    numVoids = n(), minVol = min(Vol), maxVol = max(Vol),
    meanVol = mean(Vol), medVol = median(Vol), totVol = sum(Vol),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = pos,
              values_from = c(numVoids, minVol, maxVol, meanVol, medVol, totVol)) %>%
  dplyr::mutate(across(-file_name, ~ replace_na(.x, 0)))

vs_sum <- sum_overall %>%
  full_join(sum_by_pos, by = "file_name") %>%
  dplyr::mutate(file_name = str_remove(file_name, "\\.csv$"))

# ==================================================
# Merge with LUT ----
# ==================================================
vsa <- vs_sum %>%
  full_join(lut, by = "file_name") %>%
  dplyr::mutate(
    across(numVoids_overall:totVol_middle, ~ replace_na(.x, 0)),
    Stim_bin = factor(if_else(Stim %in% c("60µL C57 F urine", "60µL C57 M urine"), "S", "N"),
                      levels = c("N", "S")),
    TrueZero = factor(case_when(
      numVoids_overall == 1 & numVoids_middle == 1 & !is.na(Stim) & totVol_overall < 110 ~ "Y",
      TRUE ~ TrueZero), levels = c("N","Y")
      )
  )

# ---- Trial 1 vs. Trial 2 ----
vsa_wide <- vsa %>% filter(!Cohort %in% ctrl_excl_cohorts) %>%
  group_by(Cohort, MouseID, Group, Sex, Age, Stim_bin, Trial) %>%
  summarise(TotVol = mean(totVol_overall), .groups = "drop") %>%
  group_by(Cohort, MouseID, Group, Sex, Age, Stim_bin) %>% arrange(Trial, .by_group = TRUE) %>%
  dplyr::mutate(trial_within_bin = row_number()) %>% filter(trial_within_bin <= 2) %>% ungroup() %>%
  dplyr::select(-Trial) %>%
  pivot_wider(names_from = trial_within_bin, values_from = TotVol, names_prefix = "Trial_")

ggplot(vsa_wide, aes(x=Trial_1, y=Trial_2, group = Stim_bin, color = Stim_bin)) + 
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_T1T2.pdf"), width = 3, height = 3, units = "cm")

vsa_wide_B <- vsa_wide %>% filter(Stim_bin == "N"); sum(is.na(vsa_wide_B))
vsa_wide_U <- vsa_wide %>% filter(Stim_bin == "S"); sum(is.na(vsa_wide_U))

cor.test(vsa_wide_B$Trial_1, vsa_wide_B$Trial_2) # 0.5429664
cor.test(vsa_wide_U$Trial_1, vsa_wide_U$Trial_2) # 0.7211486

c1 <- lm(Trial_2 ~ Trial_1 * Group * Sex * Stim_bin + Age, data = vsa_wide); summary(c1)
Anova(c1, type = 3)

c1a <- lm(Trial_2 ~ Trial_1 * Sex * Stim_bin + Age, data = vsa_wide %>% filter(Group == "BalbC_WT")); summary(c1a)
c1b <- lm(Trial_2 ~ Trial_1 * Sex * Stim_bin + Age, data = vsa_wide %>% filter(Group == "C57_WT")); summary(c1b)
c1c <- lm(Trial_2 ~ Trial_1 * Sex * Stim_bin, data = vsa_wide %>% filter(Group == "C57_MUT")); summary(c1c)

cor_result <- vsa_wide %>%
  group_by(Group, Sex) %>%
  summarise(
    r = cor(Trial_1, Trial_2, use = "complete.obs"),
    p = cor.test(Trial_1, Trial_2)$p.value,
    n = n())

# ==================================================
# Exploratory mixed models across VSA metrics (cols 2:19) ----
# ==================================================
vsa_clean <- vsa %>% filter(!Cohort %in% ctrl_excl_cohorts)
resp_vars <- names(vsa_clean)[2:19]

fit_one <- function(y) {
  f <- reformulate(c("Group*Sex", "Stim_bin", "Age", "(1|Cohort/MouseID)"), response = y)
  tidy(lmer(f, data = vsa_clean, REML = FALSE), effects = "fixed") %>%
    dplyr::mutate(response = y, .before = 1)
}

res_fixed <- map_dfr(resp_vars, possibly(fit_one, otherwise = tibble()))

write.csv(vsa_clean %>% filter(Stim_bin == "S") %>%
            dplyr::select(Cohort, Trial, Group, Sex, Age, MouseID, totVol_overall), "VS_allTrials.csv", row.names = FALSE)

# ---- Repeatability or Intraclass Correlation Coefficient (ICC) ----
vsa_icc <- vsa_clean %>% filter(Stim_bin == "S") %>%
  dplyr::mutate(Trial = factor(Trial))

icc_df <- vsa_icc %>% group_by(Group, Sex) %>% tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data, ~ lmer(totVol_overall ~ Trial + (1 | MouseID), data = .x)),
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
ggsave(file.path(plotdir, "VS_ICC.pdf"), width = 3, height = 3, units = "cm")

# ---- Function to calculate ICC
get_icc <- function(dat, covariate) {
  form <- if (covariate == "Group") totVol_overall ~ Trial + Group + (1|MouseID) else totVol_overall ~ Trial + Sex + (1|MouseID)
  m <- lmer(form, data = dat)
  vc <- as.data.frame(VarCorr(m))
  vc$vcov[vc$grp == "MouseID"] / sum(vc$vcov)
}

set.seed(123); B <- 5000

# ---- Sex: F vs M, controlling for Group
icc_F <- get_icc(filter(vsa_icc, Sex == "F"), "Group")
icc_M <- get_icc(filter(vsa_icc, Sex == "M"), "Group")
obs_sex <- icc_F - icc_M

mice <- vsa_icc %>% distinct(MouseID, Sex, Group)

perm_sex <- replicate(B, {
  perm_labels <- mice %>% group_by(Group) %>% dplyr::mutate(Sex_perm = sample(Sex)) %>% ungroup()
  d <- vsa_icc %>% dplyr::select(-Sex) %>% left_join(dplyr::select(perm_labels, MouseID, Sex_perm), by = "MouseID")
  F <- d %>% filter(Sex_perm == "F") %>% rename(Sex = Sex_perm)
  M <- d %>% filter(Sex_perm == "M") %>% rename(Sex = Sex_perm)
  get_icc(F, "Group") - get_icc(M, "Group")
})

p_sex <- (sum(abs(perm_sex) >= abs(obs_sex)) + 1) / (B + 1)
c(ICC_F = icc_F, ICC_M = icc_M, difference = obs_sex, p = p_sex)

# ---- BalbC_WT vs C57_WT, controlling for Sex
d_BC <- vsa_icc %>% filter(Group %in% c("BalbC_WT", "C57_WT"))
icc_Balb <- get_icc(filter(d_BC, Group == "BalbC_WT"), "Sex")
icc_WT <- get_icc(filter(d_BC, Group == "C57_WT"), "Sex")
obs_BC <- icc_Balb - icc_WT

mouse_BC <- d_BC %>% distinct(MouseID, Group, Sex)

perm_BC <- replicate(B, {
  labs <- mouse_BC %>% group_by(Sex) %>% dplyr::mutate(Group_perm = sample(Group)) %>% ungroup()
  d <- d_BC %>% dplyr::select(-Group) %>% left_join(dplyr::select(labs, MouseID, Group_perm), by = "MouseID")
  B <- d %>% filter(Group_perm == "BalbC_WT") %>% rename(Group = Group_perm)
  W <- d %>% filter(Group_perm == "C57_WT") %>% rename(Group = Group_perm)
  get_icc(B, "Sex") - get_icc(W, "Sex")
})

p_BC <- (sum(abs(perm_BC) >= abs(obs_BC)) + 1) / (B + 1)
c(ICC_BalbC = icc_Balb, ICC_C57WT = icc_WT, difference = obs_BC, p = p_BC)

# ---- C57_WT vs C57_MUT, controlling for Sex
d_WM <- vsa_icc %>% filter(Group %in% c("C57_WT", "C57_MUT"))
icc_WT2 <- get_icc(filter(d_WM, Group == "C57_WT"), "Sex")
icc_MUT <- get_icc(filter(d_WM, Group == "C57_MUT"), "Sex")
obs_WM <- icc_WT2 - icc_MUT

mouse_WM <- d_WM %>% distinct(MouseID, Group, Sex)

perm_WM <- replicate(B, {
  labs <- mouse_WM %>% group_by(Sex) %>% dplyr::mutate(Group_perm = sample(Group)) %>% ungroup()
  d <- d_WM %>% dplyr::select(-Group) %>% left_join(dplyr::select(labs, MouseID, Group_perm), by = "MouseID")
  W <- d %>% filter(Group_perm == "C57_WT") %>% rename(Group = Group_perm)
  M <- d %>% filter(Group_perm == "C57_MUT") %>% rename(Group = Group_perm)
  get_icc(W, "Sex") - get_icc(M, "Sex")
})

p_WM <- (sum(abs(perm_WM) >= abs(obs_WM)) + 1) / (B + 1)
c(ICC_C57WT = icc_WT2, ICC_C57MUT = icc_MUT, difference = obs_WM, p = p_WM)

# ---- BalbC_WT vs C57_MUT, controlling for Sex
d_BM <- vsa_icc %>% filter(Group %in% c("BalbC_WT", "C57_MUT"))
icc_Balb2 <- get_icc(filter(d_BM, Group == "BalbC_WT"), "Sex")
icc_MUT2 <- get_icc(filter(d_BM, Group == "C57_MUT"), "Sex")
obs_BM <- icc_Balb2 - icc_MUT2

mouse_BM <- d_BM %>% distinct(MouseID, Group, Sex)

perm_BM <- replicate(B, {
  labs <- mouse_BM %>% group_by(Sex) %>% dplyr::mutate(Group_perm = sample(Group)) %>% ungroup()
  d <- d_BM %>% dplyr::select(-Group) %>% left_join(dplyr::select(labs, MouseID, Group_perm), by = "MouseID")
  Bc <- d %>% filter(Group_perm == "BalbC_WT") %>% rename(Group = Group_perm)
  M <- d %>% filter(Group_perm == "C57_MUT") %>% rename(Group = Group_perm)
  get_icc(Bc, "Sex") - get_icc(M, "Sex")
})

p_BM <- (sum(abs(perm_BM) >= abs(obs_BM)) + 1) / (B + 1)
c(ICC_BalbC = icc_Balb2, ICC_C57MUT = icc_MUT2, difference = obs_BM, p = p_BM)

p.adjust(c(BalbC_vs_C57WT = p_BC, C57WT_vs_C57MUT = p_WM, BalbC_vs_C57MUT = p_BM), method = "holm")

# ==================================================
# Summary + plots (violin + jitter + mean ± SE) ----
# ==================================================
plot_violin <- function(y) {
  se_df <- summarySE(vsa_clean, measurevar = y, groupvars = c("Group", "Sex"), na.rm = TRUE)
  
  ggplot(vsa_clean, aes(Sex, .data[[y]])) +
    facet_wrap(vars(Group)) +
    geom_violin() +
    geom_point(aes(fill = Group), colour = "black", pch = 21,
               position = position_jitter(0.1), size = 1) +
    geom_pointrange(data = se_df,
                    aes(ymin = .data[[y]] - se, ymax = .data[[y]] + se),
                    size = 0.25) +
    theme_classic() +
    labs(x = NULL, y = y)
}

plot_violin("totVol_overall")
plot_violin("numVoids_overall")
plot_violin("meanVol_overall")
plot_violin("medVol_overall")

# ==================================================
# Hurdle models + paired spaghetti plots ----
# ==================================================
hist(vsa_clean$numVoids_overall)
hist(vsa_clean$meanVol_overall)
hist(vsa_clean$totVol_overall)

vsa_clean <- vsa_clean %>% 
  dplyr::mutate(across(2:19, ~ if_else(TrueZero == "Y", 0, .x)),
         TrueZero_num = ifelse(TrueZero == "Y", 0, 1))

resp_se <- summarySE(vsa_clean, measurevar = "TrueZero_num", 
                     groupvars = c("Group", "Sex", "Stim_bin"), na.rm = TRUE)

ggplot(resp_se, aes(Sex, TrueZero_num, fill = Stim_bin)) +
  facet_wrap(~ Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  geom_errorbar(aes(ymin = TrueZero_num - se, ymax = TrueZero_num + se),
                position = position_dodge(.9), width = 0) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_ResponseRate.pdf"), width = 5, height = 3.5, units = "cm")

m1 <- glmer(TrueZero ~ Group + Sex + Stim_bin + Age + (1 | MouseID), family = binomial, data = vsa_clean); summary(m1)
pairs(emmeans(m1, ~ Group))
pairs(emmeans(m1, ~ Stim_bin))

# ---- Effect size ----
resp_se_wide <- resp_se %>% 
  pivot_wider(id_cols = c("Group", "Sex"),
              names_from = "Stim_bin", values_from = "TrueZero_num") %>%
  dplyr::mutate(Diff = S/N)

# ---- Spaghetti Plots (includes True Zero trials) ----
vsa_spag <- vsa_clean %>%
  group_by(Cohort, MouseID, Group, Sex, Age, Stim_bin) %>%
  summarise(totVol = mean(totVol_overall),
            numVoids = mean(numVoids_overall),
            meanVol = mean(meanVol_overall),
            medVol = mean(medVol_overall),
            .groups = "drop") %>%
  add_count(MouseID, name = "paired_obs")

plot_spag <- function(y, ylab_txt) {
  ggplot(filter(vsa_spag, paired_obs == 2),
         aes(Stim_bin, .data[[y]], color = Sex, group = MouseID)) +
    facet_wrap(vars(Sex, Group), scales = "free_y") +
    geom_line(alpha = 0.8) +
    geom_point(size = 2) +
    theme_classic() +
    labs(x = "Urine stimulus", y = ylab_txt)
}

plot_spag("numVoids", "No. void spots")
plot_spag("meanVol", "Mean void volume (uL)")
plot_spag("totVol", "Total urine volume (uL)")

vsa_clean_nz <- vsa_clean %>% filter(TrueZero == "N")
hist(log10(vsa_clean_nz$totVol_overall))

m2 <- lmer(log10(totVol_overall) ~ (Group + Sex + Stim_bin)^2 + Age + (1 | MouseID), data = vsa_clean_nz); summary(m2)
Anova(m2, type = 3)                     
pairs(emmeans(m2, ~ Group))             # BalbC_WT diff
pairs(emmeans(m2, ~ Group | Sex))       # BalbC_WT diff
pairs(emmeans(m2, ~ Sex | Group))       # Not C57_MUT

ggplot(vsa_spag, aes(Stim_bin, totVol, color = Sex, group = MouseID)) +
  geom_line(linewidth = 0.25) +
  geom_point(size = 0.5) +
  facet_wrap(~ Group + Sex, nrow = 1) +
  scale_x_discrete(expand = c(.05, .05)) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_TotalVolume.pdf"), width = 7, height = 3.5, units = "cm")

# ==================================================
# Group-stratified models (Sex * Stim within Group) ----
# ==================================================
fit_by_group <- function(y, g, add_age = TRUE) {
  rhs <- paste0("Sex*Stim_bin", if (add_age) " + Age" else "", " + (1|Cohort/MouseID)")
  lmer(as.formula(paste(y, rhs, sep = " ~ ")), data = filter(vsa_clean, Group == g))
}

summary(fit_by_group("totVol_overall", "BalbC_WT"))
summary(fit_by_group("totVol_overall", "C57_WT")) # Sex:Stim
summary(fit_by_group("totVol_overall", "C57_MUT", add_age = FALSE))
summary(fit_by_group("numVoids_overall", "BalbC_WT"))
summary(fit_by_group("numVoids_overall", "C57_WT")) # Sex:Stim
summary(fit_by_group("numVoids_overall", "C57_MUT", add_age = FALSE))
summary(fit_by_group("meanVol_overall", "BalbC_WT")) # Sex:Stim
summary(fit_by_group("meanVol_overall", "C57_WT"))
summary(fit_by_group("meanVol_overall", "C57_MUT", add_age = FALSE))

# ==================================================
# SMI (Scent Mark Index): LDA classifier + LD1-based index ----
# ==================================================
# LDA features + labeled training set ----
vsa <- vsa %>%
  dplyr::mutate(
    numVoids_ratio = numVoids_middle / (numVoids_edge + numVoids_middle),
    totVol_ratio   = totVol_middle   / (totVol_edge   + totVol_middle)
  )

ld <- vsa %>%
  transmute(
    file_name, TrueZero,
    numVoids_overall, meanVol_overall, maxVol_overall, totVol_overall,
    numVoids_ratio, totVol_ratio
  )

pred_vars <- names(ld)[3:8]
pred_vars

cls <- read_csv("VoidPattern_Classifier - Classifier_NEW.csv", show_col_types = FALSE)
table(cls$Class)

cls <- cls %>%
  right_join(ld, by = "file_name") %>%
  filter(TrueZero == "N") %>%
  transmute(
    file_name, Class = factor(Class),
    across(all_of(pred_vars), ~ as.numeric(scale(.x)))
  )

train <- filter(cls, !is.na(Class))
test  <- filter(cls,  is.na(Class))

fit_lda <- function(method = c("moment", "t"), dat = train) {
  method <- match.arg(method)
  lda(x = dat %>% dplyr::select(all_of(pred_vars)), grouping = dat$Class, method = method)
}

add_preds <- function(mod, dat) {
  if (nrow(dat) == 0) return(tibble())
  p <- predict(mod, newdata = dat %>% dplyr::select(all_of(pred_vars)))
  bind_cols(as_tibble(p$x), dat,
            tibble(Predicted_Class = p$class),
            as_tibble(p$posterior))
}

cm <- function(pred, obs) confusionMatrix(table(predicted = pred, observed = obs))

m_moment <- fit_lda("moment")
m_robust <- fit_lda("t")

train_robust <- add_preds(m_robust, train)
cm(train_robust$Predicted_Class, train$Class)

cv_loo <- lda(Class ~ ., data = train %>% dplyr::select(-file_name), CV = TRUE)
cm(cv_loo$class, train$Class)

LDA_output <- bind_rows(
  train_robust %>% dplyr::mutate(Set = "Train"),
  add_preds(m_robust, test) %>% dplyr::mutate(Set = "Test")
) %>%
  dplyr::mutate(Set = factor(Set, levels = c("Train", "Test")),
         Class = factor(Class))

# LDA diagnostic plots (saved) ----
p_ellipse <- ggplot(LDA_output, aes(LD1, LD2, color = Predicted_Class)) +
  geom_point(alpha = 0.5) +
  stat_ellipse(aes(fill = Predicted_Class), alpha = 0.2, geom = "polygon") +
  theme_classic()

p_hist_ds <- ggplot(filter(LDA_output, Predicted_Class %in% c("D", "S")),
                    aes(LD1, fill = Predicted_Class, colour = Predicted_Class)) +
  geom_histogram(alpha = 0.5, position = "dodge") +
  theme_classic()

p_den_ld1 <- ggplot(LDA_output, aes(LD1, fill = Predicted_Class, colour = Predicted_Class)) +
  geom_density(alpha = 0.5) +
  theme_classic()

p_den_ld2 <- ggplot(LDA_output, aes(LD2, fill = Predicted_Class, colour = Predicted_Class)) +
  geom_density(alpha = 0.5) +
  theme_classic()

p_ellipse
p_hist_ds
p_den_ld1
p_den_ld2

ggplot(LDA_output, aes(LD1, LD2, color = Predicted_Class)) +
  geom_point(size = 0.5, alpha = 0.5) +
  stat_ellipse(aes(fill = Predicted_Class), alpha = 0.2, geom = "polygon") +
  theme_nolabels
ggsave(file.path(plotdir, "VS_Ellipse.pdf"), width = 3, height = 3, units = "cm")

# Merge LD scores back to VSA + define SMI (LD1-based) ----
vsa_NEW <- vsa_clean %>%
  left_join(LDA_output %>% dplyr::select(file_name, LD1, LD2, Predicted_Class, Set, Class), by = "file_name") %>%
  dplyr::mutate(Discord = case_when(
      Set == "Train" & Class == Predicted_Class ~ "N",
      Set == "Train" & Class != Predicted_Class ~ "Y",
      TRUE ~ NA_character_),
    SMI = {smi0 <- -LD1; smi0 <- smi0 + abs(min(smi0, na.rm = TRUE)) + 1
    }
  )

ggplot(filter(vsa_NEW, SMI >= 1), aes(SMI, fill = Predicted_Class, colour = Predicted_Class)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 14, by = 4)) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_Density.pdf"), width = 3, height = 3, units = "cm")

sum(is.na(vsa_NEW$LD1)) # 203
sum(is.na(vsa_NEW$SMI)) # 203

ggplot(vsa_NEW, aes(SMI, fill = Sex, colour = Sex)) +
  facet_wrap(~ Group) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_SMI_Sex.pdf"), width = 5, height = 3, units = "cm")

m3 <- lmer(SMI ~ Group * Sex + Stim_bin + Age + (1 | MouseID), data = vsa_NEW); summary(m3)
Anova(m3, type = 3)                     
pairs(emmeans(m3, ~ Group))             # BalbC_WT diff
pairs(emmeans(m3, ~ Sex | Group))       # BalbC_WT only

vsa_NEW %>% 
  filter(file_name %in% c("20231121_BL2_C28_2530_female_nostim", "20231121_BL2_C29_2531_male_nostim", "20240527_BL2_C71_3394_male_nostim")) %>%
  dplyr::select("file_name", "Group", "Sex", "Age", "MouseID", "totVol_overall", "LD1", "SMI", "Predicted_Class")

ggplot(vsa_NEW, aes(x=SMI, y=totVol_overall, group = Sex, color = Sex)) + 
  facet_grid(Sex ~ Group) +
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm", se = TRUE)

plot_summary <- vsa_NEW %>% group_by(Group, Sex, Predicted_Class) %>%
  summarise(mean = mean(totVol_overall, na.rm=TRUE), 
            n = n(), 
            se = sd(totVol_overall, na.rm=TRUE)/sqrt(n()),
            .groups = "drop") %>%
  filter(!is.na(Predicted_Class))

ggplot(plot_summary, aes(Sex, mean, fill = Predicted_Class)) +
  facet_wrap(~ Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(.9), width = 0)

# ----Test for bimodality ----
dip_results <- vsa_NEW %>% filter(!is.na(SMI)) %>%
  group_by(Group, Sex) %>%
  summarise(n = n(), dip_obj = list(dip.test(SMI)), 
            dip = dip_obj[[1]]$statistic, p_value = dip_obj[[1]]$p.value, .groups = "drop")

mclust_results <- vsa_NEW %>%
  filter(!is.na(SMI)) %>% group_by(Group, Sex) %>%
  summarise(n = n(), fit = list(Mclust(SMI, G = 1:3)), G_selected = fit[[1]]$G,
    model = fit[[1]]$modelName, BIC = max(fit[[1]]$bic), logLik = fit[[1]]$loglik,
    uncertainty = mean(fit[[1]]$uncertainty), .groups = "drop")

vsa_NEW <- vsa_NEW %>%
  dplyr::mutate(Predicted_Class = factor(Predicted_Class, levels = c("D", "N", "S")))

vsa_NEW_noNA <- vsa_NEW %>% filter(!is.na(Predicted_Class)) 

ggplot(vsa_NEW_noNA, aes(Sex, fill = Predicted_Class)) +
  geom_bar(position = "fill") +
  facet_wrap(~ Group) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_PredictedClass.pdf"), width = 3, height = 2.5, units = "cm") 

fisher.test(xtabs(~ Predicted_Class + Group, data = filter(vsa_NEW_noNA, Sex == "F")), workspace = 2e8)
fisher.test(xtabs(~ Predicted_Class + Group, data = filter(vsa_NEW_noNA, Sex == "M")), workspace = 2e8)

groups <- levels(droplevels(vsa_NEW_noNA$Group))
pairs_grp <- combn(groups, 2, simplify = FALSE)

posthoc_fisher_F <- map_dfr(pairs_grp, \(g) {
  tab <- xtabs(~ Predicted_Class + Group, data = filter(vsa_NEW_noNA, Sex == "F", Group %in% g))
  ft <- fisher.test(tab)
  tibble(contrast = paste(g, collapse = " vs "), p.value = ft$p.value)
})

posthoc_fisher_M <- map_dfr(pairs_grp, \(g) {
  tab <- xtabs(~ Predicted_Class + Group, data = filter(vsa_NEW_noNA, Sex == "M", Group %in% g))
  ft <- fisher.test(tab)
  tibble(contrast = paste(g, collapse = " vs "), p.value = ft$p.value)
})

fisher.test(xtabs(~ Predicted_Class + Sex, data = filter(vsa_NEW_noNA, Group == "BalbC_WT")), workspace = 2e8)
fisher.test(xtabs(~ Predicted_Class + Sex, data = filter(vsa_NEW_noNA, Group == "C57_WT")), workspace = 2e8)
fisher.test(xtabs(~ Predicted_Class + Sex, data = filter(vsa_NEW_noNA, Group == "C57_MUT")), workspace = 2e8)

# ---- Proportions ----
Class_props <- vsa_NEW_noNA %>% 
  dplyr::count(Group, Sex, Predicted_Class) %>%
  group_by(Group, Sex) %>% dplyr::mutate(percent = round(100 * n / sum(n), digits = 1)) %>% ungroup()

# ==================================================
# Hierarchy classification from totVol (mouse → cohort) ----
# ==================================================
meta_vars <- c("Cohort", "MouseID", "Treatment", "Group", "Age", "Sex")

mouse_summary <- vsa_NEW %>%
  left_join(ld %>% dplyr::select(file_name, numVoids_ratio, totVol_ratio), by = "file_name") %>%
  filter(Stim_bin == "S") %>%
  group_by(across(all_of(meta_vars))) %>%
  summarise(
    avgLD1           = mean(LD1, na.rm = TRUE),
    avgLD2           = mean(LD2, na.rm = TRUE),
    avgSMI           = mean(SMI, na.rm = TRUE),
    numVoids_overall = mean(numVoids_overall, na.rm = TRUE),
    meanVol_overall  = mean(meanVol_overall, na.rm = TRUE),
    maxVol_overall   = mean(maxVol_overall, na.rm = TRUE),
    totVol_overall   = mean(totVol_overall, na.rm = TRUE),
    numVoids_ratio   = mean(numVoids_ratio, na.rm = TRUE),
    totVol_ratio     = mean(totVol_ratio, na.rm = TRUE),
    nTrial_stim      = n(),
    .groups = "drop"
  ) %>%
  full_join(lut %>% dplyr::count(Cohort, MouseID, Treatment, Group, Age, Sex, name = "nTrials"), by = meta_vars) %>%
  filter(!Cohort %in% ctrl_excl_cohorts) %>%
  dplyr::mutate(across(avgSMI:nTrial_stim, ~ replace_na(.x, 0))) %>%
  group_by(Cohort) %>%
  dplyr::mutate(SMI_Rank = factor(dense_rank(desc(avgSMI)))) %>%
  ungroup()

global_dists <- mouse_summary %>%
  dplyr::select(Cohort, MouseID, totVol_overall) %>%
  group_by(Cohort) %>%
  filter(n() == 3) %>%
  summarise(dist = list(as.numeric(dist(totVol_overall))), .groups = "drop") %>%
  unnest(dist)

cohort_dists <- mouse_summary %>%
  dplyr::select(Cohort, MouseID, totVol_overall) %>%
  group_by(Cohort) %>%
  filter(n() == 3) %>%
  summarise(ds_sorted = list(sort(totVol_overall, decreasing = TRUE)), .groups = "drop") %>%
  dplyr::mutate(
    d1 = map_dbl(ds_sorted, 1),
    d2 = map_dbl(ds_sorted, 2),
    d3 = map_dbl(ds_sorted, 3),
    delta_12 = d1 - d2,
    delta_23 = d2 - d3,
    delta_13 = d1 - d3
  )

cuts <- quantile(cohort_dists$delta_13, c(0.25, 0.75), na.rm = TRUE)
small_cut <- unname(cuts[1])
large_cut <- unname(cuts[2])

primary_classified <- cohort_dists %>%
  dplyr::mutate(hierarchy_type = case_when(
    delta_13 <  small_cut ~ "weak",
    delta_13 >= small_cut ~ "hierarchical",
    TRUE ~ NA_character_
  ))

cohort_classified <- primary_classified %>%
  dplyr::mutate(hierarchy_type = case_when(
    hierarchy_type != "hierarchical" ~ hierarchy_type,
    d2 <= (d3 + delta_13 / 3)        ~ "despotic",
    d2 >= (d1 - delta_13 / 3)        ~ "co-dominant",
    TRUE                             ~ "transitive"
  ))

mouse_summary_full <- mouse_summary %>%
  left_join(dplyr::select(cohort_classified, Cohort, hierarchy_type), by = "Cohort")

cohort_level <- mouse_summary_full %>%
  distinct(Cohort, Group, Sex, hierarchy_type) %>%
  dplyr::mutate(hierarchy_type = factor(hierarchy_type,
                                 levels = c("despotic", "transitive", "co-dominant", "weak")))

ggplot(cohort_level, aes(Sex, fill = hierarchy_type)) +
  geom_bar(position = "fill") +
  facet_wrap(~ Group) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_HierarchyType.pdf"), width = 3, height = 2.5, units = "cm") 

fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "BalbC_WT"))) # NS
fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "C57_WT")))   # 0.01401
fisher.test(xtabs(~ hierarchy_type + Sex, data = filter(cohort_level, Group == "C57_MUT")))  # NS

fisher.test(xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "F"))) # 0.01743
fisher.test(xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "M"))) # NS

groups <- levels(droplevels(cohort_level$Group))
pairs_grp <- combn(groups, 2, simplify = FALSE)

posthoc_fisher_F2 <- map_dfr(pairs_grp, \(g) {
  tab <- xtabs(~ hierarchy_type + Group, data = filter(cohort_level, Sex == "F", Group %in% g))
  ft <- fisher.test(tab)
  tibble(contrast = paste(g, collapse = " vs "), p.value = ft$p.value)
})

# ---- Proportions ----
hierarchy_props <- cohort_level %>% 
  dplyr::count(Group, Sex, hierarchy_type) %>%
  group_by(Group, Sex) %>% dplyr::mutate(percent = 100 * n / sum(n)) %>% ungroup()

# ---- Plot totVol by Rank ----
mouse_summary_full <- mouse_summary_full %>% 
  group_by(Cohort) %>%
  dplyr::mutate(totVol_Rank = factor(dense_rank(desc(totVol_overall)))) %>%
  ungroup()

ggplot(mouse_summary_full, aes(Sex, totVol_overall, fill = totVol_Rank)) +
  facet_wrap(~ Group) +
  geom_boxplot(position = position_dodge(.7), width =.7, outlier.shape = NA) +
  #coord_cartesian(ylim = c(0, 650)) +
  theme_nolabels
ggsave(file.path(plotdir, "VS_totVol_Rank.pdf"), width = 7, height = 3.5, units = "cm")

m4 <- lm(totVol_overall ~ totVol_Rank * Group * Sex + Age, data = mouse_summary_full); summary(m4)
Anova(m4, type = 3)
pairs(emmeans(m4, ~ totVol_Rank))
pairs(emmeans(m4, ~ totVol_Rank | Sex | Group))

# ==================================================
# Export: final mouse summary
# ==================================================
mouse_summary_full_clean <- mouse_summary_full %>%
  dplyr::select(Cohort, MouseID, Treatment, Group, Age, Sex, 
         avgLD1, avgLD2, avgSMI, numVoids_overall, meanVol_overall, maxVol_overall, totVol_overall, numVoids_ratio, totVol_ratio, hierarchy_type)

mouse_summary_full_clean %>%
  rename_with(~ paste0(.x, "_VS"), .cols = 7:16) %>%
  write_csv("void_spot_clean.csv")
