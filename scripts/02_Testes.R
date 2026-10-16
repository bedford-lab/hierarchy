# Libraries ----
library(Rmisc); library(tidyverse); library(plotly); library(stringr); library(emmeans)

# Paths ----
setwd("~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data/")
plotdir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/testes"

theme_nolabels <- theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        strip.text = element_blank(), legend.position = "none")

# ==================================================
# Load ----
# ==================================================

testes <- read.csv("SDH_Testes - Testes.csv", stringsAsFactors = TRUE)
bw     <- read.csv("predicted_weights.csv", stringsAsFactors = TRUE)

outliers <- c("1068", "2409", "2127", "2160", "2220", "2221")

# ==================================================
# 1) Left vs right testis weight ----
# ==================================================

p_testes <- ggplot(testes, aes(Weight_LT, Weight_RT, color = Treatment)) + geom_point()
ggplotly(p_testes)

summary(lm(Weight_RT ~ Weight_LT, data = testes))

testes_no_out <- testes %>% filter(!MouseID %in% outliers)
summary(lm(Weight_RT ~ Weight_LT, data = testes_no_out))

# ==================================================
# 2) Testis weight vs FIJI size measurements ----
# ==================================================

testes_long <- testes %>%
  pivot_longer(cols = Weight_LT:Minor_RT, names_to = "measure", values_to = "value") %>%
  mutate(
    Testis   = if_else(str_detect(measure, "LT"), "Left", "Right"),
    measure  = str_remove(measure, "_(LT|RT)$")
  ) %>%
  pivot_wider(names_from = measure, values_from = value)

testes_long_no_out <- testes_long %>% filter(!MouseID %in% outliers)

resp_vars <- c("Area", "Perim", "Major", "Minor")
lm_fiji_results <- list()

for (resp in resp_vars) {
  fit <- lm(reformulate("Weight", response = resp), data = testes_long_no_out)
  lm_fiji_results[[resp]] <- fit
  print(summary(fit))
}

# ==================================================
# 3) Merge testes with predicted weights ----
# ==================================================

df <- testes %>%
  select(Round, MouseID, Cohort, Sex, Treatment, DOB, Sac_weight, Sac_date, Testes_weight) %>%
  full_join(bw %>% select(MouseID, Cohort, Sex, Treatment, DOB, Pred_Sac_Weight),
            by = c("MouseID","Cohort","Sex","Treatment","DOB")) %>%
  filter(Sex == "M")

df <- df %>% mutate(
  Geno       = if_else(str_detect(Treatment, "MUT"), "MUT", "WT"),
  Background = if_else(str_detect(Treatment, "BalbC"), "BalbC", "C57"),
  Group      = factor(paste(Background, Geno, sep = "_"),
                      levels = c("BalbC_WT", "C57_WT", "C57_MUT"))) %>%
  filter(!stringr::str_detect(Treatment, "Ctrl"))

nrow(testes)  # 118
nrow(df)      # 102

df_no_out <- df %>% filter(!MouseID %in% outliers)

# ==================================================
# 4) Testis weight vs actual vs predicted sac weight ----
# ==================================================

pred_vars <- c("Sac_weight", "Pred_Sac_Weight")
lm_pred_results <- list()

for (pred in pred_vars) {
  fit <- lm(reformulate(pred, response = "Testes_weight"), data = df_no_out)
  lm_pred_results[[pred]] <- fit
  print(summary(fit))
}

plot(df$Pred_Sac_Weight, df$Sac_weight)
cor.test(df$Pred_Sac_Weight, df$Sac_weight) # 0.94

df <- df %>% mutate(
    Testes_weight_adj = case_when(
      !is.na(Sac_weight) ~ (Testes_weight / Sac_weight) * 100,
      !is.na(Pred_Sac_Weight) ~ (Testes_weight / Pred_Sac_Weight) * 100,
      TRUE ~ NA_real_
      )
  )

# ---- Plot Testes_weight ----
weight_ci <- summarySE(df, measurevar = "Testes_weight", groupvars = "Group", na.rm = TRUE)

ggplot(df, aes(Group, Testes_weight)) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = weight_ci, aes(ymin = Testes_weight - ci, ymax = Testes_weight + ci),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  theme_nolabels
ggsave(file.path(plotdir, "Testes_weight.pdf"), width = 2.5, height = 3, units = "cm")

m1 <- lm(Testes_weight ~ Group, df); summary(m1)
pairs(emmeans(m1, ~ Group)) 

# ---- Plot Testes_weight_adj ----
weight_adj_ci <- summarySE(df, measurevar = "Testes_weight_adj", groupvars = "Group", na.rm = TRUE)

ggplot(df, aes(Group, Testes_weight_adj)) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = weight_adj_ci, aes(ymin = Testes_weight_adj - ci, ymax = Testes_weight_adj + ci),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  theme_nolabels
ggsave(file.path(plotdir, "Testes_weight_adj.pdf"), width = 2.5, height = 3, units = "cm")

m2 <- lm(Testes_weight_adj ~ Group, df); summary(m2)
pairs(emmeans(m2, ~ Group))  

# ---- Plot Sac_weight ----
sac_weight_ci <- summarySE(df, measurevar = "Sac_weight", groupvars = "Group", na.rm = TRUE)

ggplot(df, aes(Group, Sac_weight)) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  geom_pointrange(data = sac_weight_ci, aes(ymin = Sac_weight - ci, ymax = Sac_weight + ci),
                  size = 1, fatten = 0.5, position = position_nudge(x = 0.4),
                  color = "black") +
  theme_nolabels
ggsave(file.path(plotdir, "Sac_weight_Males.pdf"), width = 2.5, height = 3, units = "cm")

m3 <- lm(Sac_weight ~ Group, df); summary(m3)
pairs(emmeans(m3, ~ Group)) 

# ==================================================
# Export ----
# ==================================================

write.csv(df, "testes_clean.csv", row.names = FALSE)
