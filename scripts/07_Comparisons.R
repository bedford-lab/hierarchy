# Libraries / paths ----
library(Rmisc); library(tidyverse); library(lme4); library(lmerTest); library(emmeans); library(caret)
library(ordinal); library(purrr); library(irr); library(broom); library(pheatmap); library(scales); 
library(umap); library(car); library(vcd); library(Rtsne); library(networkD3); library(webshot2);
library(ggtern); library(viridis)

setwd("~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/data")
plotdir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/comp"

# Setup ----
theme_nolabels <- theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        strip.text=element_blank(), legend.position="none")

ctrl_excl <- c("C3_MW.R2","C5_MW.R1","C5_MW.R2","C6_MW.R1","C6_MW.R2",
               "C7","C12","C15","C15b","C19","C23","C31","C37","C40")

weight_vars <- c("Pred_Wean_Weight","Pred_Mat_Weight","Pred_Sac_Weight","Sac_weight",
                 "Testes_weight","Testes_weight_adj")

cont_vars <- c("wins_TT","losses_TT","win_prop_TT","mean_dailyRank_TT","Elo_TT","DS_norm_TT",
               "numBouts_WS","TotDur_WS","meanBout_dur_WS","medBout_dur_WS","maxBout_dur_WS",
               "DS_norm_WS","out_strength_WS","in_strength_WS","net_strength_WS","pagerank_WS",
               "avgSMI_VS",
               "numVoids_overall_VS","meanVol_overall_VS","maxVol_overall_VS","totVol_overall_VS",
               "numVoids_ratio_VS","totVol_ratio_VS")

categories <- c("_VS","_TT","_WS")
cat_pattern <- paste(categories, collapse="|")

# Read / merge ----
md <- read.csv("meta_data.csv", stringsAsFactors=TRUE)
tt <- read.csv("tube_test_clean.csv", stringsAsFactors=TRUE)
ws <- read.csv("warm_spot_clean.csv", stringsAsFactors=TRUE)
vs <- read.csv("void_spot_clean.csv", stringsAsFactors=TRUE)

md_trim <- md %>%
  mutate(Sac_weight = if_else(is.na(Sac_weight), Pred_Sac_Weight, Sac_weight)) %>%
  filter(!Cohort %in% ctrl_excl) %>%
  dplyr::select(Cohort, MouseID, Treatment, Sex, all_of(weight_vars)) %>%
  distinct()

join_cols <- c("Cohort","MouseID","Treatment","Group","Age","Sex")

all <- tt %>%
  full_join(ws, by=join_cols) %>%
  full_join(vs, by=join_cols) %>%
  mutate(Group = factor(Group, levels = c("BalbC_WT", "C57_WT", "C57_MUT")))

all_cont <- all %>%
  dplyr::select(Cohort, MouseID, Treatment, Group, Age, Sex, all_of(cont_vars)) %>%
  left_join(md_trim, by=c("Cohort","MouseID","Treatment","Sex"))

resp_cols <- c(cont_vars, weight_vars)

# Helpers ----
safe_cor <- function(x,y,method="pearson"){
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if(n<3 || sd(x)==0 || sd(y)==0)
    return(tibble(n_pair=n,r=NA_real_,p=NA_real_))
  ct <- suppressWarnings(cor.test(x,y,method=method))
  tibble(n_pair=n,r=unname(ct$estimate),p=ct$p.value)
}

all_by_all <- function(df,vars){
  expand_grid(var1=vars,var2=vars) %>%
    filter(var1 < var2) %>%
    mutate(out = map2(var1,var2,~safe_cor(df[[.x]],df[[.y]]))) %>%
    unnest(out)
}

corr_mat <- function(df,vars,min_finite=3){
  v <- vars[map_lgl(vars,~sum(is.finite(df[[.x]]))>=min_finite)]
  if(length(v)<2) return(NULL)
  C <- cor(df[,v],use="pairwise.complete.obs")
  ok <- rowSums(is.finite(C))>0
  C <- C[ok,ok,drop=FALSE]
  if(nrow(C)<2) NULL else C
}

# Correlations ----
cor_results <- all_cont %>%
  group_by(Group,Sex) %>%
  group_modify(~{
    v <- resp_cols
    if(.y$Sex[[1]]=="F") v <- setdiff(v,c("Testes_weight","Testes_weight_adj"))
    all_by_all(.x,v)
  }) %>%
  ungroup() %>%
  group_by(Group,Sex) %>%
  mutate(p_fdr = p.adjust(p,"fdr")) %>%
  ungroup()

# Cross-category significant ----
all_cor <- cor_results %>%
  mutate(cat1=str_extract(var1,cat_pattern),
         cat2=str_extract(var2,cat_pattern)) %>%
  filter(p<0.05, !is.na(cat1), !is.na(cat2), cat1!=cat2,
         !str_detect(var1,regex("eight",TRUE)),
         !str_detect(var2,regex("eight",TRUE)))

# Variable importance ----
var_scores <- cor_results %>%
  mutate(cat1=str_extract(var1,cat_pattern),
         cat2=str_extract(var2,cat_pattern),
         weight=abs(r)*(p<0.05)) %>%
  filter(!is.na(cat1), !is.na(cat2), cat1!=cat2) %>%
  pivot_longer(c(var1,var2), names_to="side", values_to="variable") %>%
  mutate(category=str_extract(variable,cat_pattern)) %>%
  group_by(category,variable) %>%
  summarise(sum_weight=sum(weight,na.rm=TRUE),
            n_sig=sum(p<0.05,na.rm=TRUE), .groups="drop")

best_vars <- var_scores %>%
  group_by(category) %>%
  slice_max(sum_weight,n=3,with_ties=FALSE) %>%
  ungroup()

# Best pairs ----
best_pairs <- cor_results %>%
  mutate(cat1=str_extract(var1,cat_pattern),
         cat2=str_extract(var2,cat_pattern),
         pair=map2_chr(cat1,cat2,~paste(sort(c(.x,.y)),collapse="_")),
         sig_fdr=p_fdr<0.1,
         score=if_else(sig_fdr,abs(r),0),
         Sig=if_else(p<0.05,"*","")) %>%
  filter(!is.na(cat1), !is.na(cat2), cat1!=cat2) %>%
  group_by(Group,Sex,pair) %>%
  slice_max(score,n=1,with_ties=FALSE) %>%
  ungroup() %>%
  dplyr::select(Group,Sex,pair,var1,var2,r,p,p_fdr,Sig,score)

# Heatmaps ----
all_cont %>%
  group_by(Group,Sex) %>%
  group_walk(~{
    v <- resp_cols
    if(.y$Sex[[1]]=="F") v <- setdiff(v,c("Testes_weight","Testes_weight_adj"))
    C <- corr_mat(.x,v)
    if(!is.null(C))
      pheatmap(C,
               clustering_distance_rows=as.dist(1-abs(C)),
               clustering_distance_cols=as.dist(1-abs(C)),
               main=paste(.y$Group[[1]],.y$Sex[[1]],sep=" | "))
  })

# ==================================================
# Hierarchy Type Heatmaps ----
# ==================================================
hier_levels <- c("despotic", "transitive", "co-dominant", "weak")

type <- all %>% 
  dplyr::select(Cohort, Group, Sex, starts_with("hierarchy_type")) %>%
  distinct() %>% drop_na() %>%
  mutate(across(starts_with("hierarchy_type"),
                ~ factor(.x, levels = hier_levels))) %>%
  mutate(consistency = apply(
      dplyr::select(., hierarchy_type_TT, hierarchy_type_WS, hierarchy_type_VS), 1,
      \(x) c(`1` = 3, `2` = 2, `3` = 0)[as.character(n_distinct(x))]))

table(type$consistency)

make_confusion <- function(df, x, y) {
  df %>%
    dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>%
    dplyr::count(Group, Sex, .data[[x]], .data[[y]], name = "n") %>%
    group_by(Group, Sex, .data[[x]]) %>%
    mutate(prop = n / sum(n)) %>% ungroup() %>%
    complete(Group, Sex, !!sym(x), !!sym(y), fill = list(n = 0, prop = 0))
}

save_confusion_plot <- function(df, x, y, filename) {
  p <- ggplot(df, aes(forcats::fct_rev(.data[[x]]), forcats::fct_rev(.data[[y]]), fill = prop)) +
    geom_tile(color = NA) +
    facet_grid(Sex ~ Group) +
    scale_fill_gradientn(colours = c("black", viridisLite::viridis(100, option = "magma")),
      limits = c(0, 1), labels = percent, name = "Proportion") + coord_equal() + theme_nolabels
  print(p)
  ggsave(file.path(plotdir, filename),
         plot = p,
         width = 7, height = 3, units = "cm")
}

tt_ws <- make_confusion(type, "hierarchy_type_TT", "hierarchy_type_WS")
tt_vs <- make_confusion(type, "hierarchy_type_TT", "hierarchy_type_VS")
ws_vs <- make_confusion(type, "hierarchy_type_WS", "hierarchy_type_VS")

save_confusion_plot(tt_ws, "hierarchy_type_TT", "hierarchy_type_WS", "TT_WS.pdf")
save_confusion_plot(tt_vs, "hierarchy_type_TT", "hierarchy_type_VS", "TT_VS.pdf")
save_confusion_plot(ws_vs, "hierarchy_type_WS", "hierarchy_type_VS", "WS_VS.pdf")

tt_ws_diag <- tt_ws %>% filter(hierarchy_type_TT == hierarchy_type_WS)

# Statistics ----
cohorts <- all %>% distinct(Cohort, Group, Sex,
           hierarchy_type_TT, hierarchy_type_WS, hierarchy_type_VS) %>% filter(Cohort != "C32") %>%
  mutate(Agreement = ((hierarchy_type_TT==hierarchy_type_WS) + (hierarchy_type_TT==hierarchy_type_VS) + 
                        (hierarchy_type_WS==hierarchy_type_VS))/3)

set.seed(123); B <- 5000
agree_sum <- cohorts %>% group_by(Group,Sex) %>% nest() %>% mutate(
  Agreement = map_dbl(data, ~mean(.x$Agreement)),
  boot = map(data, ~replicate(B, mean(sample(.x$Agreement, replace = TRUE)))),
  low = map_dbl(boot, ~quantile(.x, .025)), high = map_dbl(boot, ~quantile(.x, .975))) %>%
  dplyr::select(-data, -boot)

ggplot(agree_sum, aes(x = Group, y = Agreement, color = Sex)) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = .9), width = 0) +
  geom_point(position = position_dodge(width = .9), size = 3) +
  ylim(0, .6) + theme_nolabels
ggsave(file.path(plotdir, "Agreement.pdf"), width = 3, height = 3, units = "cm")

# Hierarchy type consistency ----
tests <- c("TT", "WS", "VS")
hier_types <- c("weak", "co-dominant", "transitive", "despotic")

hier_summary <- cohorts %>%
  pivot_longer(c(hierarchy_type_TT, hierarchy_type_WS, hierarchy_type_VS),
               names_to = "test", values_to = "hierarchy_type") %>%
  mutate(
    test = dplyr::recode(test, hierarchy_type_TT = "TT", hierarchy_type_WS = "WS", hierarchy_type_VS = "VS"),
    test = factor(test, levels = tests),
    hierarchy_type = factor(hierarchy_type, levels = hier_types)
  ) %>%
  count(Sex, Group, test, hierarchy_type, name = "n") %>%
  right_join(expand_grid(
    Sex = unique(cohorts$Sex), Group = unique(cohorts$Group),
    test = factor(tests, levels = tests),
    hierarchy_type = factor(hier_types, levels = hier_types)
  ), by = c("Sex", "Group", "test", "hierarchy_type")) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(Sex, Group, test) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  ungroup()

ggplot(filter(hier_summary, Sex == "F"),
       aes(hierarchy_type, proportion, fill = test)) +
  facet_wrap(~Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  coord_cartesian(ylim = c(0, 0.7)) +
  theme_nolabels
ggsave(file.path(plotdir, "HierarchyType_Prop_F.pdf"), width = 9, height = 3.5, units = "cm")

ggplot(filter(hier_summary, Sex == "M"),
       aes(hierarchy_type, proportion, fill = test)) +
  facet_wrap(~Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  coord_cartesian(ylim = c(0, 0.7)) +
  theme_nolabels
ggsave(file.path(plotdir, "HierarchyType_Prop_M.pdf"), width = 9, height = 3.5, units = "cm")

hier_summary %>%
  filter(hierarchy_type == "transitive") %>% group_by(test) %>%
  slice_max(proportion, n = 1, with_ties = TRUE) %>% ungroup() %>%
  select(test, Group, Sex, n, total, proportion)

# ==================================================
# Targeted correlations ----
# ==================================================
targ_vars <- c("DS_norm_TT","TotDur_WS","totVol_overall_VS")

all_cont <- all_cont %>%
  mutate(across(all_of(targ_vars), ~ as.numeric(scale(.)), .names = "z_{.col}"))

hist(all_cont$z_DS_norm_TT)
hist(all_cont$z_TotDur_WS)
hist(all_cont$z_totVol_overall_VS)

# Warm spot vs Tube test ----
m1 <- lmer(z_TotDur_WS ~ z_DS_norm_TT * Group * Sex + Age + (1 | Cohort), data = all_cont); summary(m1); Anova(m1, type = 3)
emm_slope <- emtrends(m1, specs = c("Group", "Sex"), var = "z_DS_norm_TT"); summary(emm_slope) # C57_MUT F
summary(emm_slope, infer = TRUE) # 0.0263

# Void spot vs Tube test ----
m2 <- lmer(z_totVol_overall_VS ~ z_DS_norm_TT * Group * Sex + Age + (1 | Cohort), data = all_cont); summary(m2); Anova(m2, type = 3) 
emm_slope <- emtrends(m2, specs = c("Group", "Sex"), var = "z_DS_norm_TT"); summary(emm_slope)
summary(emm_slope, infer = TRUE) # NS

# Void spot vs Warm spot ----
m3 <- lmer(z_totVol_overall_VS ~ z_TotDur_WS * Group * Sex + Age + (1 | Cohort), data = all_cont); summary(m3); Anova(m3, type = 3) 
emm_slope <- emtrends(m3, specs = c("Group", "Sex"), var = "z_TotDur_WS"); summary(emm_slope) # BalbC_WT M
summary(emm_slope, infer = TRUE) # 0.0006

# Z-score Plots ----
ggplot(all_cont, aes(x = z_DS_norm_TT, y = z_TotDur_WS, group = Sex, color = Sex)) +
  facet_wrap(vars(Group), scales = "free_x") +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  theme_nolabels
ggsave(file.path(plotdir, "TT_WS_lm.pdf"), width = 7, height = 3, units = "cm")

ggplot(all_cont, aes(x = z_DS_norm_TT, y = z_totVol_overall_VS, group = Sex, color = Sex)) +
  facet_wrap(vars(Group), scales = "free_x") +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  theme_nolabels
ggsave(file.path(plotdir, "TT_VS_lm.pdf"), width = 7, height = 3, units = "cm")

ggplot(all_cont, aes(x = z_TotDur_WS, y = z_totVol_overall_VS, group = Sex, color = Sex)) +
  facet_wrap(vars(Group), scales = "free_x") +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  theme_nolabels
ggsave(file.path(plotdir, "WS_VS_lm.pdf"), width = 7, height = 3, units = "cm")

# Weight correlations ----
weight_cor <- cor_results %>%
  filter((str_detect(var1,regex("eight",TRUE)) |
            str_detect(var2,regex("eight",TRUE))) & p<0.05)

# Wean Weight (4 weeks) ----
m4 <- lmer(z_DS_norm_TT ~ Pred_Wean_Weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m4)
emm_slope <- emtrends(m4, specs = c("Group", "Sex"), var = "Pred_Wean_Weight"); summary(emm_slope) # NS

m5 <- lmer(z_TotDur_WS ~ Pred_Wean_Weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m5)
emm_slope <- emtrends(m5, specs = c("Group", "Sex"), var = "Pred_Wean_Weight"); summary(emm_slope) # NS

m6 <- lmer(z_totVol_overall_VS ~ Pred_Wean_Weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m6)
emm_slope <- emtrends(m6, specs = c("Group", "Sex"), var = "Pred_Wean_Weight"); summary(emm_slope) # NS

# Maturity Weight (8 weeks) ----
m7 <- lmer(z_DS_norm_TT ~ Pred_Mat_Weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m7)
emm_slope <- emtrends(m7, specs = c("Group", "Sex"), var = "Pred_Mat_Weight"); summary(emm_slope) # NS

m8 <- lmer(z_TotDur_WS ~ Pred_Mat_Weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m8)
emm_slope <- emtrends(m8, specs = c("Group", "Sex"), var = "Pred_Mat_Weight"); summary(emm_slope) # NS

m9 <- lmer(z_totVol_overall_VS ~ Pred_Mat_Weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m9)
emm_slope <- emtrends(m9, specs = c("Group", "Sex"), var = "Pred_Mat_Weight"); summary(emm_slope) # NS

# Sac Weight ----
all_cont <- all_cont %>% mutate(z_Sac_weight = scale(Sac_weight),
                                z_Testes_weight_adj = scale(Testes_weight_adj))

m10 <- lmer(z_DS_norm_TT ~ z_Sac_weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m10)
emm_slope <- emtrends(m10, specs = c("Group", "Sex"), var = "z_Sac_weight"); summary(emm_slope) # NS

m11 <- lmer(z_TotDur_WS ~ z_Sac_weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m11)
emm_slope <- emtrends(m11, specs = c("Group", "Sex"), var = "z_Sac_weight"); summary(emm_slope) # NS

m12 <- lmer(z_totVol_overall_VS ~ Sac_weight * Group * Sex + Age + (1|Cohort), data = all_cont); summary(m12)
emm_slope <- emtrends(m12, specs = c("Group", "Sex"), var = "Sac_weight"); summary(emm_slope) # BalbC_WT M
summary(emm_slope, infer = TRUE)
pairs(emm_slope, by = "Sex")   
pairs(emm_slope, by = "Group")

ggplot(all_cont, aes(x = z_Sac_weight, y = z_totVol_overall_VS, group = Sex, color = Sex)) +
  facet_wrap(vars(Group), scales = "free_x") +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  theme_nolabels
ggsave(file.path(plotdir, "totVol_Sac_weight.pdf"), width = 7, height = 3, units = "cm")

ggplot(all_cont %>% filter(Group == "BalbC_WT"), 
       aes(x = Sac_weight, y = z_totVol_overall_VS, group = Sex, color = Sex)) +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  xlim(15.8, 31.2) +
  theme_nolabels
ggsave(file.path(plotdir, "totVol_Sac_weight_BalbC_WT.pdf"), width = 3, height = 3, units = "cm")

# Testes Weight ----
m13 <- lmer(z_DS_norm_TT ~ z_Testes_weight_adj * Group + (1|Cohort), data = all_cont); summary(m13)
emm_slope <- emtrends(m13, specs = "Group", var = "z_Testes_weight_adj"); summary(emm_slope) # NS

m14 <- lmer(z_TotDur_WS ~ z_Testes_weight_adj * Group + (1|Cohort), data = all_cont); summary(m14)
emm_slope <- emtrends(m13, specs = "Group", var = "z_Testes_weight_adj"); summary(emm_slope) # NS

m15 <- lmer(z_totVol_overall_VS ~ z_Testes_weight_adj * Group + (1|Cohort), data = all_cont); summary(m15)
emm_slope <- emtrends(m15, specs = "Group", var = "z_Testes_weight_adj"); summary(emm_slope) # BalbC_WT
summary(emm_slope, infer = TRUE)

ggplot(all_cont, aes(x = z_Testes_weight_adj, y = z_totVol_overall_VS, group = Sex, color = Sex)) +
  facet_wrap(vars(Group), scales = "free_x") +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  theme_nolabels
ggsave(file.path(plotdir, "totVol_Testes.pdf"), width = 7, height = 3, units = "cm")

# ==================================================
# Consistency ----
# ==================================================
z_vars <- c("z_DS_norm_TT", "z_TotDur_WS", "z_totVol_overall_VS")

z_long <- all_cont %>%
  filter(Cohort != "C32") %>%
  dplyr::select(MouseID, Cohort, Group, Sex, Age, all_of(z_vars)) %>%
  pivot_longer(cols = all_of(z_vars), names_to = "Assay", values_to = "z_score") %>%
  mutate(
    Assay = dplyr::recode(Assay,
      z_DS_norm_TT = "TT",
      z_TotDur_WS = "WS",
      z_totVol_overall_VS = "VS"))

consistency <- z_long %>% group_by(MouseID, Cohort, Group, Sex, Age) %>%
  summarise(SD = sd(z_score, na.rm = TRUE),
            range = max(z_score, na.rm = TRUE) - min(z_score, na.rm = TRUE),
            .groups = "drop")

m16 <- lmer(SD ~ Group * Sex + Age + (1 | Cohort), data = consistency); summary(m16)
pairs(emmeans(m16, ~ Sex)) # NS
pairs(emmeans(m16, ~ Group)) # NS
pairs(emmeans(m16, ~ Group | Sex))

# ==================================================
# Rank consistency ----
# ==================================================
all_cont_clean_Rank <- all_cont %>% filter(Cohort != "C32") %>%
  dplyr::select(Cohort:Sex, all_of(targ_vars)) %>%
  group_by(Cohort) %>%
  mutate(across(all_of(targ_vars), ~ {
    r <- rank(-., ties.method = "min", na.last = "keep")
    tied <- duplicated(.) | duplicated(., fromLast = TRUE)
    if_else(tied & r > 1, sum(!is.na(.)), r)
  }, .names = "{.col}_rank")) %>% 
  ungroup()

all_cont_clean_Rank <- all_cont_clean_Rank %>%
  mutate(TT_WS = abs(DS_norm_TT_rank - TotDur_WS_rank),
         TT_VS = abs(DS_norm_TT_rank - totVol_overall_VS_rank),
         WS_VS = abs(TotDur_WS_rank - totVol_overall_VS_rank)) %>%
  mutate(n_zero = rowSums(across(c(TT_WS, TT_VS, WS_VS), ~ .x == 0)))

all_cont_clean_Rank %>% group_by(Group, Sex) %>%
  summarise(
    N = n(),
    n_consistent = sum(n_zero == 3),
    percent_consistent = 100 * n_consistent / N,
    .groups = "drop")

table(all_cont_clean_Rank$n_zero) # 19/201 = 9.5%

# ---- Plot rank shifts ----
plotdir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/comp/sankey_plots"

sankey_input <- all_cont_clean_Rank %>%
  transmute(Group, Sex,
    TT = factor(paste0("TT_", DS_norm_TT_rank), levels = c("TT_1","TT_2","TT_3")),
    WS = factor(paste0("WS_", TotDur_WS_rank), levels = c("WS_1","WS_2","WS_3")),
    VS = factor(paste0("VS_", totVol_overall_VS_rank), levels = c("VS_1","VS_2","VS_3")),
    TT_end = factor(paste0("TT_end_", DS_norm_TT_rank), levels = c("TT_end_1","TT_end_2","TT_end_3")))

make_sankey <- function(dat, group_sel, sex_sel) {
  d <- dat %>% filter(Group == group_sel, Sex == sex_sel)
  links <- bind_rows(d %>% count(source = TT, target = WS, name = "value"),
                     d %>% count(source = WS, target = VS, name = "value"),
                     d %>% count(source = VS, target = TT_end, name = "value"))
  node_order <- c("TT_1","TT_2","TT_3", "WS_1","WS_2","WS_3", "VS_1","VS_2","VS_3", "TT_end_1","TT_end_2","TT_end_3")
  nodes <- tibble(name = node_order) %>% filter(name %in% c(as.character(links$source), as.character(links$target))) %>%
    mutate(label = sub("TT_end_", "TT_", name))
  links <- links %>% mutate(source = match(as.character(source), nodes$name) - 1,
                            target = match(as.character(target), nodes$name) - 1)
  sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target",
                Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30,
                sinksRight = FALSE, iterations = 0)
}

make_sankey(sankey_input, "BalbC_WT", "F")
make_sankey(sankey_input, "BalbC_WT", "M")

clean_name <- function(x) stringr::str_replace_all(as.character(x), "[^A-Za-z0-9]+", "_")
panels <- sankey_input %>% distinct(Group, Sex) %>% arrange(Group, Sex)

purrr::walk2(panels$Group, panels$Sex, \(g, s) {
  p <- make_sankey(sankey_input, g, s)
  base <- paste0("sankey_", clean_name(g), "_", clean_name(s))
  html_file <- file.path(plotdir, paste0(base, ".html"))
  pdf_file  <- file.path(plotdir, paste0(base, ".pdf"))
  htmlwidgets::saveWidget(widget = p, file = html_file, selfcontained = TRUE)
  webshot2::webshot(url = html_file, file = pdf_file,
    vwidth = 1200, vheight = 600, delay = 1)
})

# ---- Plot min/max normalized data ----
plotdir <- "~/Dropbox/UWyo/Projects/SocialHierarchy/Fall25_SDH/SDH/plots/comp"

min_max_norm <- function(x) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

df <- all_cont_clean_Rank %>% 
  group_by(Group, Sex) %>%
  mutate(across(c(DS_norm_TT, TotDur_WS, totVol_overall_VS), min_max_norm, .names = "{.col}_mm")) %>%
  rename(
    TT_mm = DS_norm_TT_mm,
    WS_mm = TotDur_WS_mm,
    VS_mm = totVol_overall_VS_mm) %>%
  mutate(tot_dom_mm = rowSums(across(c(TT_mm, WS_mm, VS_mm)))) %>% ungroup() %>%
  group_by(Cohort) %>%
  mutate(tot_dom_mm_rank = rank(-tot_dom_mm, ties.method = "min", na.last = "keep")) %>% ungroup()

df <- df %>% left_join(all_cont %>% dplyr::select(MouseID, z_DS_norm_TT, z_TotDur_WS, z_totVol_overall_VS)) %>%
  mutate(tot_dom_z = rowSums(across(c(z_DS_norm_TT, z_TotDur_WS, z_totVol_overall_VS)))) %>%
  group_by(Cohort) %>%
  mutate(tot_dom_z_rank = rank(-tot_dom_z, ties.method = "min", na.last = "keep"),
         tot_dom_mm_rank_diff = tot_dom_mm_rank - tot_dom_z_rank) %>% ungroup()

df_long <- df %>% group_by(Cohort) %>% 
  dplyr::select(-c(DS_norm_TT, TotDur_WS, totVol_overall_VS, DS_norm_TT_rank, TotDur_WS_rank, totVol_overall_VS_rank)) %>%
  pivot_longer(cols = c(TT_mm, WS_mm, VS_mm), names_to = "test") %>%
  mutate(
    test = factor(test, levels = c("TT_mm", "WS_mm", "VS_mm")),
    tot_dom__mm_rank = factor(tot_dom_mm_rank, levels = c(1, 2, 3)))

df_summary <- df_long %>%
  group_by(Group, Sex, tot_dom_mm_rank, test) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(sum(!is.na(value))),
            .groups = "drop")

ggplot(filter(df_summary, Sex == "F"),
       aes(tot_dom_mm_rank, mean, fill = test)) +
  facet_wrap(~Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(.9), width = 0) +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_nolabels
ggsave(file.path(plotdir, "Dominance_byRank_F.pdf"), width = 9, height = 3.5, units = "cm")

ggplot(filter(df_summary, Sex == "M"),
       aes(tot_dom_mm_rank, mean, fill = test)) +
  facet_wrap(~Group) +
  geom_bar(stat = "identity", position = position_dodge(.9), width =.9) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(.9), width = 0) +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_nolabels
ggsave(file.path(plotdir, "Dominance_byRank_M.pdf"), width = 9, height = 3.5, units = "cm")

m18_F <- lmer(value ~ Group * tot_dom_mm_rank * test + (1|MouseID), data = df_long %>% filter(Sex == "F")); summary(m18_F); Anova(m18_F, type = 3)
pairs(emmeans(m18_F, ~ test | Group | tot_dom_mm_rank)) # Pairwise contrasts
joint_tests(m18_F, by = c("Group", "tot_dom_mm_rank")) # Overall assay effect

m18_M <- lmer(value ~ Group * tot_dom_mm_rank * test + (1|MouseID), data = df_long %>% filter(Sex == "M")); summary(m18_M); Anova(m18_M, type = 3)
pairs(emmeans(m18_M, ~ test | Group | tot_dom_mm_rank)) # Pairwise contrasts
joint_tests(m18_M, by = c("Group", "tot_dom_mm_rank")) # Overall assay effect

jt_F <- as.data.frame(joint_tests(m18_F, by = c("Group", "tot_dom_mm_rank"))) %>% mutate(Sex = "F")
jt_M <- as.data.frame(joint_tests(m18_M, by = c("Group", "tot_dom_mm_rank"))) %>% mutate(Sex = "M")

jt <- bind_rows(jt_F, jt_M) %>% group_by(Sex) %>%
  mutate(p_holm = p.adjust(p.value, method = "holm")) %>%
  ungroup()

jt_plot <- jt %>%
  mutate(sig = case_when(p_holm < 0.001 ~ "***", 
                         p_holm < 0.01 ~ "**", 
                         p_holm < 0.05 ~ "*", TRUE ~ ""))

sig_df <- df_summary %>%
  group_by(Sex, Group, tot_dom_mm_rank) %>%
  dplyr::summarise(
    y = max(mean + se, na.rm = TRUE) + 0.08,
    .groups = "drop") %>%
  left_join(
    jt_plot %>% dplyr::select(Sex, Group, tot_dom_mm_rank, sig),
    by = c("Sex", "Group", "tot_dom_mm_rank"))

# ---- Make Ternary plots ----
tern_df <- df %>% 
dplyr::select(-tot_dom_z, -tot_dom_z_rank, -tot_dom_mm_rank_diff) %>%
  mutate(TT_dom = TT_mm/tot_dom_mm,
         WS_dom = WS_mm/tot_dom_mm,
         VS_dom = VS_mm/tot_dom_mm,
         tot = TT_dom + WS_dom + VS_dom)

write.csv(tern_df, "tern_df.csv", row.names = FALSE)

ggtern(tern_df %>% filter(Sex == "F"),
       aes(x = TT_dom, z = WS_dom, y = VS_dom)) +
  stat_density_tern(geom = "polygon", n = 200, h = 0.5, base = "identity",
                    aes(fill = after_stat(level))) +
  geom_point(color = "black", size = 0.5) +
  facet_grid(~Group) +
  scale_fill_viridis_c(option = "plasma", direction = 1)
ggsave(file.path(plotdir, "TernaryPlot_F.pdf"), width = 16, height = 4, units = "cm")

ggtern(tern_df %>% filter(Sex == "M"),
       aes(x = TT_dom, z = WS_dom, y = VS_dom)) +
  stat_density_tern(geom = "polygon", n = 200, h = 0.5, base = "identity",
                    aes(fill = after_stat(level))) +
  geom_point(color = "black", size = 0.5) +
  facet_grid(~Group) +
  scale_fill_viridis_c(option = "plasma", direction = 1)
ggsave(file.path(plotdir, "TernaryPlot_M.pdf"), width = 16, height = 4, units = "cm")

h_common <- mean(c(MASS::bandwidth.nrd(tern_df$TT_dom),
                   MASS::bandwidth.nrd(tern_df$WS_dom),
                   MASS::bandwidth.nrd(tern_df$VS_dom)))

plot_combos <- expand_grid(Sex = levels(tern_df$Sex), Group = levels(tern_df$Group))

purrr::pwalk(plot_combos, function(Sex, Group) {
  i <- Sex
  j <- Group
  
  p <- ggtern(tern_df %>% filter(Sex == i, Group == j),
              aes(x = TT_dom, z = WS_dom, y = VS_dom)) +
    stat_density_tern(geom = "polygon", n = 400, h = h_common, base = "identity",
                      aes(fill = after_stat(level))) +
    geom_point(color = "black", size = 0.5) +
    scale_fill_viridis_c(option = "plasma", direction = 1, trans = "sqrt",
                         name = "Density") +
    theme_nolabels
  
  print(p)
  ggsave(file.path(plotdir, paste0("TernaryPlot_", j, "_", i, ".pdf")),
         plot = p, width = 3, height = 3, units = "cm")
})