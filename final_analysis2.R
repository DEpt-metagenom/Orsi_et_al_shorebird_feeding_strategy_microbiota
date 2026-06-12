###############################################################################
# Microbiota shows major difference in case of two shorebird species with different feeding strategy
# Ákos Őrsi, Levente Laczkó, Renáta Bőkényné Tóth, Csongor Freytag, Pál Tóth, Gábor Simay, Nándor Szabó, Gábor Kardos, Ádám Lovas-Kiss
################################################################################

############################
# 0) Working directories
############################

# Windows
#setwd("C:/Users/orsia/Dropbox/Microbiome ms/stat")

# Mac
setwd("/Users/lovas-kissadam/Library/CloudStorage/Dropbox/Ákos dropbox/Microbiome_TRIGLA_GALGAL/stat")

getwd()


############################
# 1) Packages
############################

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
  library(stringr)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(forcats)
  library(readr)
})


############################
# 2) Load data
############################

mvb <- read.csv2("merged_wider_javított2.csv",
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

colnames(mvb) <- gsub("\\[|\\]", "", colnames(mvb))


# Remove samples with zero total abundance
mvb$.__tot <- rowSums(mvb[, -c(1:2)], na.rm = TRUE)
mvb <- subset(mvb, .__tot > 0, select = -.__tot)

# Clean bird names
mvb$BIRD <- str_trim(mvb$BIRD)

# Community matrix
comm_raw <- mvb[, -c(1:2)]

# Make sure all taxa columns are numeric
comm_raw <- as.data.frame(
  lapply(comm_raw, function(x) suppressWarnings(as.numeric(as.character(x))))
)

# Replace NA values with zero
comm_raw[is.na(comm_raw)] <- 0

# Fourth-root transformed community matrix
# This is the same transformation you used previously before Bray-Curtis
comm <- comm_raw^0.25

# Bird factor with full species names
mvb$Bird_species <- factor(
  mvb$BIRD,
  levels = c("GALGAL", "TRIGLA"),
  labels = c("Gallinago gallinago", "Tringa glareola")
)


###############################################################################
# SECTION A — PERMANOVA, dispersion test, SIMPER and PCoA
###############################################################################

############################
# A1) Bray-Curtis dissimilarity
############################

bray_dist <- vegdist(comm, method = "bray")


############################
# A2) PERMANOVA
############################

pmv <- adonis2(
  bray_dist ~ Bird_species,
  data = mvb,
  permutations = 999
)

print(pmv)

write.csv(
  as.data.frame(pmv),
  "PERMANOVA_BrayCurtis_results.csv"
)


############################
# A3) Test for homogeneity of multivariate dispersion
############################



bd <- betadisper(
  bray_dist,
  group = mvb$Bird_species
)

print(bd)
print(anova(bd))
print(TukeyHSD(bd))

capture.output(
  {
    print(bd)
    print(anova(bd))
    print(TukeyHSD(bd))
  },
  file = "Betadisper_BrayCurtis_results.txt"
)


############################
# A4) SIMPER analysis
############################

sim <- simper(
  comm,
  group = mvb$Bird_species,
  permutations = 999
)

sim_sum <- summary(sim)

full_df <- do.call(
  rbind,
  lapply(names(sim_sum), function(comp) {
    
    df <- as.data.frame(sim_sum[[comp]])
    df$Taxon <- rownames(sim_sum[[comp]])
    df$Comparison <- comp
    rownames(df) <- NULL
    
    df[, c("Comparison", "Taxon",
           setdiff(names(df), c("Comparison", "Taxon")))]
  })
)

write.csv(
  full_df,
  "SIMPER_full_results.csv",
  row.names = FALSE
)

sig_list <- lapply(sim, function(x) {
  df <- as.data.frame(x)
  df[df$p <= 0.05, ]
})

sig_names <- names(Filter(function(dd) nrow(dd) > 0, sig_list))

if (length(sig_names) > 0) {
  
  sig_df <- do.call(
    rbind,
    lapply(sig_names, function(name) {
      cbind(Comparison = name, sig_list[[name]])
    })
  )
  
  write.csv(
    sig_df,
    "SIMPER_significant_results.csv",
    row.names = FALSE
  )
}


############################
# A5) PCoA ordination — Bray-Curtis
############################

pcoa_res <- cmdscale(
  bray_dist,
  eig = TRUE,
  k = 2
)

pcoa_df <- data.frame(
  sample = mvb$sample,
  Bird_species = mvb$Bird_species,
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2]
)

# Percentage variance explained
eig_vals <- pcoa_res$eig
var_exp <- round(
  100 * eig_vals / sum(eig_vals[eig_vals > 0]),
  1
)

# PERMANOVA label
R2 <- round(pmv$R2[1], 3)
P_value <- signif(pmv$`Pr(>F)`[1], 3)

pcoa_label <- paste0(
  "PERMANOVA\n",
  "R² = ", R2,
  "\np = ", P_value
)

p_pcoa <- ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = Bird_species
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(linewidth = 0.8, linetype = 2) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = pcoa_label,
    hjust = 1.1,
    vjust = 1.4,
    size = 4
  ) +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    colour = "Bird species"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right"
  )

p_pcoa

ggsave(
  "Figure_PCoA_BrayCurtis.png",
  p_pcoa,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  "Figure_PCoA_BrayCurtis.pdf",
  p_pcoa,
  width = 7,
  height = 5
)


###############################################################################
# SECTION B — Diversity indices and differences
###############################################################################

df_div <- read.csv2(
  "merged_wider_javított2.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(df_div) <- gsub("\\[|\\]", "", colnames(df_div))


if (!"sample" %in% names(df_div)) {
  df_div$sample <- paste0("sample", seq_len(nrow(df_div)))
}

df_div <- df_div %>%
  mutate(
    BIRD = str_trim(BIRD),
    sample = str_trim(sample)
  ) %>%
  relocate(BIRD, sample)

taxa_cols <- setdiff(names(df_div), c("BIRD", "sample"))

X <- as.data.frame(df_div[, taxa_cols, drop = FALSE])

X_num <- as.data.frame(
  lapply(X, function(col) {
    suppressWarnings(as.numeric(as.character(col)))
  })
)

X_num[is.na(X_num)] <- 0

# Remove zero-abundance samples
df_div$.__tot <- rowSums(X_num, na.rm = TRUE)
keep_rows <- df_div$.__tot > 0

df_div <- df_div[keep_rows, ]
X_num <- X_num[keep_rows, ]

df_div$.__tot <- NULL

row_fun <- function(v) {
  
  s <- sum(v)
  
  if (s <= 0) {
    return(
      c(
        shan = NA_real_,
        lambda = NA_real_,
        simpson = NA_real_,
        inv = NA_real_,
        rich = 0
      )
    )
  }
  
  p <- v / s
  lambda <- sum(p^2)
  
  c(
    shan = vegan::diversity(v, index = "shannon"),
    lambda = lambda,
    simpson = 1 - lambda,
    inv = 1 / lambda,
    rich = vegan::specnumber(v)
  )
}

M <- t(apply(X_num, 1, row_fun)) %>%
  as.data.frame()

diversity_df <- df_div %>%
  bind_cols(M) %>%
  rename(
    shannon = shan,
    simpson_1mD = simpson,
    invsimpson = inv,
    lambda = lambda,
    richness = rich
  )

div2 <- diversity_df %>%
  filter(BIRD %in% c("GALGAL", "TRIGLA")) %>%
  mutate(
    Bird_species = recode(
      BIRD,
      "GALGAL" = "Gallinago gallinago",
      "TRIGLA" = "Tringa glareola"
    ),
    Bird_species = factor(
      Bird_species,
      levels = c("Gallinago gallinago", "Tringa glareola")
    )
  )


############################
# B1) Wilcoxon tests and effect sizes
############################

effects <- bind_rows(
  rstatix::wilcox_effsize(div2, shannon ~ Bird_species) %>%
    mutate(metric = "Shannon"),
  
  rstatix::wilcox_effsize(div2, simpson_1mD ~ Bird_species) %>%
    mutate(metric = "Simpson (1 − D)"),
  
  rstatix::wilcox_effsize(div2, invsimpson ~ Bird_species) %>%
    mutate(metric = "Inverse Simpson (1/D)"),
  
  rstatix::wilcox_effsize(div2, richness ~ Bird_species) %>%
    mutate(metric = "Species richness")
) %>%
  select(metric, effsize, magnitude)

tests_tbl <- bind_rows(
  rstatix::wilcox_test(div2, shannon ~ Bird_species) %>%
    mutate(metric = "Shannon"),
  
  rstatix::wilcox_test(div2, simpson_1mD ~ Bird_species) %>%
    mutate(metric = "Simpson (1 − D)"),
  
  rstatix::wilcox_test(div2, invsimpson ~ Bird_species) %>%
    mutate(metric = "Inverse Simpson (1/D)"),
  
  rstatix::wilcox_test(div2, richness ~ Bird_species) %>%
    mutate(metric = "Species richness")
) %>%
  rstatix::add_significance() %>%
  mutate(alternative = "two.sided") %>%
  select(
    metric,
    group1,
    group2,
    n1,
    n2,
    statistic,
    p,
    p.signif,
    alternative
  ) %>%
  left_join(
    effects %>% rename(r = effsize),
    by = "metric"
  ) %>%
  mutate(
    W = round(statistic, 1),
    p = signif(p, 3),
    r = round(r, 3),
    stars = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(
    metric,
    group1,
    group2,
    n1,
    n2,
    W,
    p,
    stars,
    r,
    magnitude
  )

print(tests_tbl)

write.csv(
  tests_tbl,
  "Alpha_diversity_Wilcoxon_effect_sizes.csv",
  row.names = FALSE
)


############################
# B2) Diversity figure with n, p and effect size
############################

label_tbl <- tests_tbl %>%
  mutate(
    label = paste0(
      "\np = ", format(p, digits = 3, scientific = TRUE),
      "\nr = ", r
    )
  ) %>%
  select(metric, label)

mk_div_plot <- function(var, ylab, metric_name) {
  
  lab <- label_tbl %>%
    filter(metric == metric_name) %>%
    pull(label)
  
  ggplot(
    div2,
    aes(
      x = Bird_species,
      y = .data[[var]],
      fill = Bird_species
    )
  ) +
    geom_boxplot(
      width = 0.6,
      alpha = 0.85,
      outlier.shape = 21
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.45,
      size = 1.8
    ) +
    annotate(
      "text",
      x = 1.5,
      y = Inf,
      vjust = 1.25,
      label = lab,
      size = 3.5
    ) +
    labs(
      x = NULL,
      y = ylab
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 20,
        hjust = 1
      )
    )
}

p1 <- mk_div_plot(
  "shannon",
  "Shannon diversity",
  "Shannon"
)

p2 <- mk_div_plot(
  "simpson_1mD",
  "Simpson (1 − D)",
  "Simpson (1 − D)"
)

p3 <- mk_div_plot(
  "invsimpson",
  "Inverse Simpson (1/D)",
  "Inverse Simpson (1/D)"
)

p4 <- mk_div_plot(
  "richness",
  "Species richness",
  "Species richness"
)

diversity_panel <- (p1 | p2) / (p3 | p4)

diversity_panel

ggsave(
  "Figure_alpha_diversity_n_p_effectsize.png",
  diversity_panel,
  width = 9,
  height = 7,
  dpi = 600
)

ggsave(
  "Figure_alpha_diversity_n_p_effectsize.pdf",
  diversity_panel,
  width = 9,
  height = 7
)


###############################################################################
# SECTION C — Within-species Bray-Curtis dissimilarity + permutation test
###############################################################################

dist_mat <- as.matrix(bray_dist)
grp <- mvb$Bird_species

within_stats <- do.call(
  rbind,
  lapply(unique(grp), function(sp) {
    
    idx <- which(grp == sp)
    dvals <- dist_mat[idx, idx]
    dvals <- dvals[upper.tri(dvals)]
    
    data.frame(
      Species = sp,
      Mean_within_dissimilarity = mean(dvals),
      SD_within_dissimilarity = sd(dvals),
      N_pairs = length(dvals)
    )
  })
)

print(within_stats)

write.csv(
  within_stats,
  "Within_species_BrayCurtis_dissimilarity.csv",
  row.names = FALSE
)

set.seed(123)

nperm <- 999

perm_results <- do.call(
  rbind,
  lapply(unique(grp), function(sp) {
    
    idx <- which(grp == sp)
    
    obs <- mean(as.dist(dist_mat[idx, idx]))
    
    perm_means <- replicate(nperm, {
      rand_idx <- sample(seq_len(nrow(comm)), length(idx))
      mean(as.dist(dist_mat[rand_idx, rand_idx]))
    })
    
    p_value <- mean(perm_means <= obs)
    
    data.frame(
      Species = sp,
      Obs_mean = obs,
      P_value = p_value
    )
  })
)

print(perm_results)

write.csv(
  perm_results,
  "Within_species_permutation_test.csv",
  row.names = FALSE
)


###############################################################################
# SECTION D — Relative abundance stacked bar plots
###############################################################################

df <- read.csv2(
  "merged_wider_javított2.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    BIRD = str_trim(BIRD),
    sample = str_trim(sample),
    Bird_species = recode(
      BIRD,
      "GALGAL" = "Gallinago gallinago",
      "TRIGLA" = "Tringa glareola",
      .default = BIRD
    )
  ) %>%
  mutate(
    across(
      -c(BIRD, sample, Bird_species),
      ~ suppressWarnings(as.numeric(as.character(.)))
    )
  )

colnames(df) <- gsub("\\[|\\]", "", colnames(df))


df[is.na(df)] <- 0

long_df <- df %>%
  pivot_longer(
    cols = -c(BIRD, sample, Bird_species),
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  group_by(sample) %>%
  mutate(
    total_abundance = sum(Abundance, na.rm = TRUE),
    RelAbundance = ifelse(
      total_abundance > 0,
      Abundance / total_abundance,
      0
    )
  ) %>%
  ungroup() %>%
  filter(total_abundance > 0) %>%
  select(-total_abundance)


############################
# D1) Simplified stacked bar plot: Top 10 taxa + Other
############################

top_taxa <- long_df %>%
  group_by(Taxon) %>%
  summarise(
    mean_abundance = mean(RelAbundance, na.rm = TRUE),
    prevalence = mean(RelAbundance > 0),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Taxon)

long_top10 <- long_df %>%
  mutate(
    Taxon_group = if_else(
      Taxon %in% top_taxa,
      Taxon,
      "Other"
    )
  ) %>%
  group_by(Bird_species, sample, Taxon_group) %>%
  summarise(
    RelAbundance = sum(RelAbundance),
    .groups = "drop"
  )

sample_levels <- long_top10 %>%
  distinct(Bird_species, sample) %>%
  arrange(Bird_species, sample) %>%
  pull(sample)

long_top10 <- long_top10 %>%
  mutate(
    sample = factor(sample, levels = sample_levels),
    Taxon_group = fct_relevel(Taxon_group, "Other", after = Inf)
  )

p_bar_top10 <- ggplot(
  long_top10,
  aes(
    x = sample,
    y = RelAbundance,
    fill = Taxon_group
  )
) +
  geom_col(width = 0.9) +
  facet_grid(
    ~ Bird_species,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_x_discrete(
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Taxon"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_bar_top10

ggsave(
  "Figure_stacked_bar_top10_taxa.png",
  p_bar_top10,
  width = 10,
  height = 5.5,
  dpi = 600
)

ggsave(
  "Figure_stacked_bar_top10_taxa.pdf",
  p_bar_top10,
  width = 10,
  height = 5.5
)


############################
# D2) Mean relative abundance plot by bird species
# This is an even simpler summary figure
############################

mean_taxa_top10 <- long_df %>%
  mutate(
    Taxon_group = if_else(
      Taxon %in% top_taxa,
      Taxon,
      "Other"
    )
  ) %>%
  group_by(Bird_species, sample, Taxon_group) %>%
  summarise(
    RelAbundance = sum(RelAbundance),
    .groups = "drop"
  ) %>%
  group_by(Bird_species, Taxon_group) %>%
  summarise(
    Mean_rel_abundance = mean(RelAbundance),
    SE = sd(RelAbundance) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Taxon_group = fct_reorder(Taxon_group, Mean_rel_abundance)
  )

p_bar_mean <- ggplot(
  mean_taxa_top10,
  aes(
    x = Bird_species,
    y = Mean_rel_abundance,
    fill = Taxon_group
  )
) +
  geom_col(width = 0.75) +
  scale_y_continuous(
    labels = percent_format(),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Mean relative abundance",
    fill = "Taxon"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    legend.position = "right"
  )

p_bar_mean

ggsave(
  "Figure_mean_relative_abundance_top10_taxa.png",
  p_bar_mean,
  width = 8,
  height = 5,
  dpi = 600
)

ggsave(
  "Figure_mean_relative_abundance_top10_taxa.pdf",
  p_bar_mean,
  width = 8,
  height = 5
)



###############################################################################
# SECTION E — Core microbiome per bird
###############################################################################

df_core <- read.csv2(
  "merged_wider_javított2.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    BIRD = str_trim(BIRD),
    sample = str_trim(sample),
    Bird_species = recode(
      BIRD,
      "GALGAL" = "Gallinago gallinago",
      "TRIGLA" = "Tringa glareola",
      .default = BIRD
    )
  ) %>%
  mutate(
    across(
      -c(BIRD, sample, Bird_species),
      ~ suppressWarnings(as.numeric(as.character(.)))
    )
  )

colnames(df_core) <- gsub("\\[|\\]", "", colnames(df_core))


df_core[is.na(df_core)] <- 0

long_full <- df_core %>%
  pivot_longer(
    cols = -c(BIRD, sample, Bird_species),
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  mutate(
    Taxon = str_trim(Taxon)
  ) %>%
  group_by(sample) %>%
  mutate(
    total_abundance = sum(Abundance),
    RelAbundance = ifelse(
      total_abundance > 0,
      Abundance / total_abundance,
      0
    )
  ) %>%
  ungroup() %>%
  select(-total_abundance)

prevalence_threshold <- 0.80
abundance_threshold <- 0.001

prev_summary <- long_full %>%
  mutate(
    present = RelAbundance >= abundance_threshold
  ) %>%
  group_by(Bird_species, Taxon) %>%
  summarise(
    prevalence = mean(present),
    mean_abundance = mean(RelAbundance),
    .groups = "drop"
  )

core_microbiome <- prev_summary %>%
  filter(prevalence >= prevalence_threshold) %>%
  arrange(
    Bird_species,
    desc(prevalence),
    desc(mean_abundance)
  )

write.csv(
  core_microbiome,
  "Core_microbiome_per_bird.csv",
  row.names = FALSE
)

plot_df <- core_microbiome %>%
  group_by(Bird_species) %>%
  arrange(desc(mean_abundance), .by_group = TRUE) %>%
  mutate(
    Taxon = factor(Taxon, levels = rev(unique(Taxon)))
  ) %>%
  ungroup()

gg_core <- ggplot(
  plot_df,
  aes(
    x = Taxon,
    y = mean_abundance,
    fill = prevalence
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = percent(prevalence, accuracy = 1)),
    hjust = -0.15,
    size = 3.3,
    color = "black"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_gradient(
    limits = c(prevalence_threshold, 1),
    low = "#9ecae1",
    high = "#3b5bdb",
    oob = scales::squish,
    name = "Prevalence"
  ) +
  facet_wrap(
    ~ Bird_species,
    ncol = 1,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "Mean relative abundance"
  ) +
  theme_cowplot(font_size = 12) +
  theme(
    plot.margin = margin(5.5, 30, 5.5, 5.5),
    strip.text = element_text(face = "bold", size = 14)
  )

gg_core

ggsave(
  "core_microbiome_per_bird.pdf",
  gg_core,
  width = 180,
  height = 160,
  units = "mm"
)

ggsave(
  "core_microbiome_per_bird.png",
  gg_core,
  width = 1600,
  height = 1400,
  units = "px",
  dpi = 300
)


###############################################################################
# END OF SCRIPT
###############################################################################