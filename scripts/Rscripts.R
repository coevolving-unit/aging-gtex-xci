setwd("~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/COEVOLviNG/XCI_escape_aging")

library(biomaRt)
library(data.table)
library(dplyr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggplot2)
library(ggrepel)
library(egg)
library(rstatix)
library(readxl)
library(stringr)
library(paletteer)


#### 1. Load & merge data ####

d <- setDT(read.csv("GTEX_ASE_table_for_DeCasien.csv")) # Gylemo et al. medRxiv (WES only)
setnames(d, c(3, 4), c("tissue", "external_gene_name"))

g <- setDT(read.csv("gylemo.csv")) # Gylemo et al. eLife (WGS+WES)
setnames(g, c(3, 6), c("external_gene_name", "tissue"))

g[, participant := fcase(
  individual == "UPIC",    "GTEX-UPIC",
  individual == "nmXCI-2", "GTEX-ZZPU",
  individual == "nmXCI-1", "GTEX-13PLJ",
  default = "nmXCI"
)]
g$age = ifelse(g$participant == "GTEX-UPIC", 21, NA)
g$age = ifelse(g$participant == "GTEX-13PLJ", 60, g$age)
g$age = ifelse(g$participant == "GTEX-ZZPU", 50, g$age)

sample_lookup <- unique(
  d[!is.na(sampleID), .(participant, tissue, sampleID)])
g[sample_lookup, sampleID := i.sampleID, on = .(participant, tissue)]
table(is.na(g$sampleID))

cols_to_add <- c("external_gene_name", "participant", "tissue","sampleID",
                 "refCount", "altCount", "totalCount", "allelic_expression")
g_keys   <- unique(g[, .(external_gene_name, sampleID)])
d_keys   <- unique(d[, .(external_gene_name, sampleID)])
missing  <- g[g_keys[!d_keys, on = names(g_keys)],
              on = names(g_keys), ..cols_to_add]
age_map  <- unique(d[!is.na(age), .(participant, age)])
missing  <- age_map[missing, on = "participant"]
d        <- rbindlist(list(d, missing), fill = TRUE)

table(is.na(d$sampleID))
View(table(d$tissue, d$external_gene_name))

#### end ####

#### 2. Gene → chromosome annotation #####

# hsap <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 115)
# 
# gene2chrom <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
#   filters    = "external_gene_name",
#   values     = unique(d$external_gene_name),
#   mart       = hsap
# )
# saveRDS(gene2chrom, file = 'gene2chrom.rds')

gene2chrom = readRDS('gene2chrom.rds')
gene2chrom <- subset(gene2chrom, chromosome_name %in% c("X", "Y"))

multi_ensembl <- gene2chrom |>
  group_by(external_gene_name) |>
  summarise(
    n_ensembl  = n_distinct(ensembl_gene_id, na.rm = TRUE),
    n_chrom    = n_distinct(chromosome_name, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  filter(n_ensembl > 1 | n_chrom > 1) |>
  filter(external_gene_name != "CA5BP1")   # two Ensembl IDs, both on X

gene2chrom <- setDT(subset(gene2chrom, chromosome_name == "X"))
gene2chrom <- gene2chrom[, .SD[1], by = external_gene_name]
d <- merge(d, gene2chrom, by = "external_gene_name", all.x = TRUE)

# Recover ENSG-named rows
ensg_mask <- is.na(d$ensembl_gene_id) & startsWith(d$external_gene_name, "ENSG")
table(ensg_mask)
d$ensembl_gene_id[ensg_mask] <- d$external_gene_name[ensg_mask]

# ensg_extra <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
#   filters    = "ensembl_gene_id",
#   values     = unique(d$external_gene_name[ensg_mask]),
#   mart       = hsap
# )
#saveRDS(ensg_extra, file = 'ensg_extra.rds')

ensg_extra = readRDS('ensg_extra.rds')
id2chrom <- setNames(ensg_extra$chromosome_name, ensg_extra$ensembl_gene_id)
mask <- d$ensembl_gene_id %in% names(id2chrom)
d$chromosome_name[mask] <- id2chrom[d$ensembl_gene_id[mask]]

d <- subset(d, chromosome_name == "X")
d[, one_minus_AE := 1-allelic_expression]

#### end ####

#### 3. PAR & published escape status ####

balaton     <- setDT(read.csv("balaton.csv"))
setnames(balaton, 1, "external_gene_name")
balaton = subset(balaton, external_gene_name != 'Centromere')
table(balaton$Balaton.consensus.calls)
# Discordant          E   Mostly E   Mostly S  Mostly VE    No call        PAR          S         VE 
# 44         29         26        131         10        514         22        331         37

tukiainen   <- setDT(read_excel("Tukiainen_escape_genes.xlsx"))
setnames(tukiainen, "Gene name", "external_gene_name")
table(tukiainen$`Reported XCI status`)
# Escape Inactive  Unknown Variable 
# 82      392      120       89 
table(tukiainen$`Possible new call for escape (support from several analyses)`)
# x 
# 7 
tukiainen$combined_status = ifelse(is.na(tukiainen$`Possible new call for escape (support from several analyses)`) == TRUE, tukiainen$`Reported XCI status`, "Escape")
table(tukiainen$combined_status)
# Escape Inactive  Unknown Variable 
# 89      386      119       89 

merged <- merge(balaton, tukiainen, by = "external_gene_name", all = TRUE)

merged <- merged %>%
  mutate(
    escape_both     = Balaton.consensus.calls == "E" & `Reported XCI status` == "Escape",
    escape_either   = Balaton.consensus.calls %in% c("E", "Mostly E") | `Reported XCI status` == "Escape",
    
    silenced_both     = Balaton.consensus.calls == "S" & `Reported XCI status` == "Inactive",
    silenced_either   = Balaton.consensus.calls %in% c("S", "Mostly S") | `Reported XCI status` == "Inactive",
    
    variable_both     = Balaton.consensus.calls == "VE" & `Reported XCI status` == "Variable",
    variable_either   = Balaton.consensus.calls %in% c("VE", "Mostly VE", "Discordant") | `Reported XCI status` == "Variable",
    
    PAR = Balaton.consensus.calls == "PAR" | Region == 'PAR'
  )

colSums(merged[, c(
  "escape_both",   "escape_either",
  "silenced_both", "silenced_either",
  "variable_both", "variable_either", 
  "PAR"
)], na.rm = TRUE)

par_genes        <- merged$external_gene_name[which(merged$PAR)] # N = 27

escape_both      <- merged$external_gene_name[which(merged$escape_both)] # N = 22
escape_either    <- merged$external_gene_name[which(merged$escape_either)] # N = 94

silenced_both    <- merged$external_gene_name[which(merged$silenced_both)] # N = 263
silenced_either  <- merged$external_gene_name[which(merged$silenced_either)] # N = 504

variable_both    <- merged$external_gene_name[which(merged$variable_both)] # N = 26
variable_either  <- merged$external_gene_name[which(merged$variable_either)] # N = 130

d[, PAR := ifelse(external_gene_name %in% par_genes, "PAR", "NPX")]

d[, pl := fcase(
  PAR == "PAR",                              "PAR",
  external_gene_name %in% silenced_both,    "silenced",
  external_gene_name %in% escape_both,      "escape",
  external_gene_name %in% variable_both,    "variable",
  default = "other"
)]

table(silenced_either %in% escape_either)
table(silenced_either %in% variable_either)
table(escape_either %in% variable_either)

d[, pl2 := fcase(
  PAR == "PAR",                              "PAR",
  external_gene_name %in% variable_either,    "variable", 
  external_gene_name %in% escape_either,      "escape",
  external_gene_name %in% silenced_either,    "silenced",
  default = "other"
)]

#### end ####

#### 4. Individual- and tissue-level skew stats ####
ind_stats <- d[external_gene_name %in% silenced_both, .(
  AE_center = median(one_minus_AE, na.rm = TRUE),
  AE_IQR    = IQR(one_minus_AE,    na.rm = TRUE)
), by = .(participant, tissue)]

global_m  <- median(ind_stats$AE_center)
global_sd <- sd(ind_stats$AE_center)
ind_stats[, m1  := AE_center <= global_m - 1.0 * global_sd]
ind_stats[, m15 := AE_center <= global_m - 1.5 * global_sd]
ind_stats[, m2 := AE_center <= global_m - 2.0 * global_sd]
ind_stats[, m25  := AE_center <= global_m - 2.5 * global_sd]
ind_stats[, m3  := AE_center <= global_m - 3.0 * global_sd]
table(ind_stats$m1) 
table(ind_stats$m15) 
table(ind_stats$m2) 
table(ind_stats$m25) 
table(ind_stats$m3) 

all_tissues <- levels(factor(ind_stats$tissue))
n_tissues   <- length(all_tissues)

tissue_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#7570B3",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#1B9E77", "#D95F02", "#FF7F00", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "blue",
  "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#FCCDE5", "#BC80BD"
)[1:n_tissues]

tissue_shapes <- c(0:25, 0, 1, 2, 3, 4, 5)[1:n_tissues]

ind_plot <- ggplot(ind_stats, aes(x = reorder(participant, AE_center, FUN = median), y = AE_center)) +
  geom_boxplot() +
  geom_point(
    data = ind_stats[m3 == TRUE],
    aes(color = tissue, shape = tissue),
    size = 2
  ) +
  geom_hline(yintercept = global_m, linetype = 'dashed') +
  geom_hline(yintercept = global_m - 1.0 * global_sd, linetype = 'dashed', color = 'blue') +
  geom_hline(yintercept = global_m - 1.5 * global_sd, linetype = 'dashed', color = 'deeppink') +
  geom_hline(yintercept = global_m - 2.0 * global_sd, linetype = 'dashed', color = 'green') +
  geom_hline(yintercept = global_m - 2.5 * global_sd, linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = global_m - 3.0 * global_sd, linetype = 'dashed', color = 'red') +
  scale_shape_manual(values = c(0:25, 0, 1, 2)) +
  scale_color_manual(values = tissue_colors) +
  theme_article() +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(
    color = guide_legend(ncol = 2),
    shape = guide_legend(ncol = 2)
  )


dh <- ind_stats[d, on = c("participant", "tissue")]
dh$group = ifelse(dh$age < 50, "young", "old")

length(unique(subset(dh, m3 == TRUE)$participant))
length(unique(subset(dh, m3 == TRUE)$tissue))
table(unique(subset(dh, m3 == TRUE)[,c('sampleID','tissue')])$tissue)
length(unique(subset(dh, m3 == TRUE)$sampleID))
ma = median(subset(dh, m3 == TRUE)$age)
ma
View(unique(subset(dh, m3 == TRUE)[,c('participant','sampleID','age')]))

# Age plot
age_plot <- ggplot(
  unique(subset(dh, m3 == TRUE)[, c('participant', 'age', 'group')]),
  aes(x = age, fill = group)
) +
  geom_histogram(
    boundary  = 0,
    color     = "white",
    linewidth = 0.3,
    alpha     = 0.85,
    position  = "identity"
  ) +
  geom_vline(
    xintercept = 50,
    linetype   = "dashed",
    linewidth  = 0.6,
    color      = "grey30"
  ) +
  scale_fill_manual(
    values = c("young" = "#a00000", "old" = "#298c8c"),
    labels = c("young" = "Young", "old" = "Old")
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Age (years)", y = "Number of participants", fill = NULL) +
  theme_article() +
  theme(
    legend.position    = "left",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )
ggsave("age_m3.pdf", age_plot, height = 4, width = 7)

# m3 individuals - by gene group
pl2_colors <- c(
  "PAR"      = "#73d2de",
  "variable" = "#8f2d56",
  "escape"   = "#ffbc42",
  "silenced" = "#d81159"
)

m3_ind_plot <- ggplot(subset(dh, m3 == TRUE & pl2 != "other"), aes(
  x    = reorder(pl2, one_minus_AE, FUN = median),
  y    = one_minus_AE,
  fill = pl2,
  color = pl2
)) +
  geom_violin(alpha = 0.25, linewidth = 0.4, trim = TRUE) +
  geom_boxplot(
    width         = 0.38,
    alpha         = 0.7,
    linewidth     = 0.4,
    outlier.shape = 1       
  ) +
  scale_fill_manual(values  = pl2_colors) +
  scale_color_manual(values = pl2_colors) +
  scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, 0.1),
  ) +
  labs(
    x    = NULL,
    y    = "1 - Allelic Expression"
  ) +
  theme_article() +
  theme(
    legend.position  = "none",          
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.title       = element_text(size = 11, face = "bold", margin = margin(b = 8))
  )
ggsave("m3_AE_byGroup.pdf", m3_ind_plot, height = 3, width = 6)

#### end ####

#### 5. prior consensus analysis ####

length(silenced_both) # 263
length(variable_either) # 130

# ----- 5.1 pie charts --------------
m3_sil = subset(dh, m3 == TRUE & external_gene_name %in% silenced_both)
m3_var = subset(dh, m3 == TRUE & external_gene_name %in% variable_either)

hist(m3_sil$age)
hist(m3_var$age)
unique(m3_sil$age)[order(unique(m3_sil$age))]
unique(m3_var$age)[order(unique(m3_var$age))]
unique(m3_sil[,c('participant','age')])
unique(m3_var[,c('participant','age')])
table(m3_sil$participant)
table(m3_var$participant)

m3_sil$currentescape = ifelse(m3_sil$one_minus_AE > 0.6, 'escaping', 'silenced') #modified to 1-AE
m3_var$currentescape = ifelse(m3_var$one_minus_AE > 0.6, 'escaping', 'silenced') #modified to 1-AE
m3_sil$decade = str_sub(m3_sil$age, 1, 1)
m3_var$decade = str_sub(m3_var$age, 1, 1)
table(unique(m3_sil[,c('participant','decade')])$decade)
table(unique(m3_var[,c('participant','decade')])$decade)


#run this twice: once with m3-silenced and once with m3-variable
gene_escape_groups_sil <- m3_sil %>%
  mutate(
    age_group2 = case_when(
      decade %in% c("2", "4") ~ "young",
      decade %in% c("5", "6") ~ "old",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(external_gene_name),
    !is.na(age_group2),
    !is.na(currentescape)
  ) %>%
  mutate(
    is_escaping = currentescape == "escaping"
  ) %>%
  group_by(external_gene_name, age_group2) %>%
  summarise(
    escaping_in_group = any(is_escaping, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(external_gene_name) %>%
  summarise(
    groups_with_data = paste(sort(unique(age_group2)), collapse = ","),
    escaping_groups = paste(
      sort(unique(age_group2[escaping_in_group])),
      collapse = ","
    ),
    n_groups_with_data = n_distinct(age_group2),
    n_escaping_groups = sum(escaping_in_group),
    escape_group = case_when(
      n_escaping_groups == 0 ~ "silenced in all samples",
      TRUE ~ paste0("escaping in ", escaping_groups)
    ),
    .groups = "drop"
  )

gene_escape_groups_var <- m3_var %>%
  mutate(
    age_group2 = case_when(
      decade %in% c("2", "4") ~ "young",
      decade %in% c("5", "6") ~ "old",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(external_gene_name),
    !is.na(age_group2),
    !is.na(currentescape)
  ) %>%
  mutate(
    is_escaping = currentescape == "escaping"
  ) %>%
  group_by(external_gene_name, age_group2) %>%
  summarise(
    escaping_in_group = any(is_escaping, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(external_gene_name) %>%
  summarise(
    groups_with_data = paste(sort(unique(age_group2)), collapse = ","),
    escaping_groups = paste(
      sort(unique(age_group2[escaping_in_group])),
      collapse = ","
    ),
    n_groups_with_data = n_distinct(age_group2),
    n_escaping_groups = sum(escaping_in_group),
    escape_group = case_when(
      n_escaping_groups == 0 ~ "silenced in all samples",
      TRUE ~ paste0("escaping in ", escaping_groups)
    ),
    .groups = "drop"
  )

gene_escape_groups_sil = subset(gene_escape_groups_sil, n_groups_with_data > 1)
gene_escape_groups_var = subset(gene_escape_groups_var, n_groups_with_data > 1)
length(unique(gene_escape_groups_sil$external_gene_name))
length(unique(gene_escape_groups_var$external_gene_name)) # 104 (silenced) 40 (variable)


pie_silenced <- gene_escape_groups_sil %>%
  count(escape_group) %>%
  mutate(prop = n / sum(n), dataset = "Silenced")

pie_variable <- gene_escape_groups_var %>%
  count(escape_group) %>%
  mutate(prop = n / sum(n), dataset = "Variable")

pie_df <- bind_rows(pie_silenced, pie_variable) %>%
  group_by(dataset) %>%
  mutate(
    ymax    = cumsum(n),
    ymin    = lag(ymax, default = 0),
    label_y = (ymin + ymax) / 2,
    label   = paste0(n, "\n(", scales::percent(prop, accuracy = 0.1), ")")
  ) %>%
  ungroup()

pie_p <- ggplot(pie_df, aes(x = "", y = n, fill = escape_group)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(y = label_y, label = label),
            size = 3, color = "white", fontface = "bold", lineheight = 0.9) +
  coord_polar(theta = "y") +
  facet_wrap(~dataset, scales = "free") +
  scale_fill_paletteer_d("fishualize::Antennarius_commerson") +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12)) +
  labs(fill = "Escape pattern", title = "Escape status across age groups")
ggsave("pie_charts.pdf", pie_p, height = 4, width = 8)

# ----- 5.2 combined plots -----
# -- Silenced genes --
heat_df <- m3_sil %>%
  semi_join(gene_escape_groups_sil, by = "external_gene_name") %>%
  left_join(
    gene_escape_groups_sil %>%
      dplyr::select(external_gene_name, escape_group),
    by = "external_gene_name"
  )
length(unique(heat_df$external_gene_name))
heat_df$currentescape = ifelse(heat_df$one_minus_AE > 0.6, 'escaping', 'silenced')

# order samples by age
sample_order <- heat_df %>%
  distinct(sampleID, age) %>%
  arrange(age, sampleID) %>%
  pull(sampleID)

# order genes grouped by escape_group
gene_order <- heat_df %>%
  distinct(external_gene_name, escape_group) %>%
  arrange(escape_group, external_gene_name) %>%
  pull(external_gene_name)

heat_df <- heat_df %>%
  mutate(
    sampleID = factor(sampleID, levels = sample_order),
    external_gene_name = factor(external_gene_name, levels = rev(gene_order))
  )

sample_ann <- m3_sil %>%
  distinct(sampleID, participant, tissue, age) %>%
  filter(sampleID %in% sample_order) %>%
  mutate(sampleID = factor(sampleID, levels = sample_order))
sample_ann$group = ifelse(sample_ann$age >= 50, 'old', 'young')

heat_df <- heat_df %>%
  mutate(currentescape = factor(currentescape, levels = c("silenced", "escaping")))

escaping_genes <- heat_df %>%
  filter(currentescape == "escaping") %>%
  pull(external_gene_name) %>%
  unique()

heat_df <- heat_df %>%
  filter(external_gene_name %in% escaping_genes) %>%
  mutate(external_gene_name = droplevels(external_gene_name))

sample_ann <- sample_ann %>%
  mutate(
    participant = factor(participant),
    tissue      = factor(tissue),
    group       = factor(group)
  )

tissue_levels <- levels(sample_ann$tissue)
tissue_cols <- setNames(
  as.character(Polychrome::createPalette(
    length(tissue_levels),
    seedcolors = c("#FF0000", "#00FF00", "#0000FF")
  )),
  tissue_levels
)

participant_cols <- setNames(
  qualitative_hcl(n = nlevels(sample_ann$participant), palette = "Dynamic"),
  levels(sample_ann$participant)
)

group_cols <- setNames(
  qualitative_hcl(n = nlevels(sample_ann$group), palette = "Set 2"),
  levels(sample_ann$group)
)

currentescape_cols <- c(
  "silenced"  = "darkblue",
  "escaping"  = "lightgreen"
)

p_heat <- ggplot(heat_df, aes(x = sampleID, y = external_gene_name, color = currentescape)) +
  geom_point(size = 2, alpha = 1) +
  facet_grid(escape_group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = currentescape_cols, name = "Escape status", na.value = "grey80") +
  theme_article() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.text.y = element_text(angle = 0)) +
  labs(x = NULL, y = "Gene")

p_participant <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = participant)) +
  geom_tile() +
  scale_fill_manual(values = participant_cols, name = "Participant") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), legend.key.size = unit(0.25, "cm")) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))

p_age <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = age)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Age") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

p_age_group <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = group_cols, name = "Group") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

p_tissue <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = tissue)) +
  geom_tile() +
  scale_fill_manual(values = tissue_cols, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), legend.key.size = unit(0.25, "cm")) +
  guides(fill = guide_legend(ncol = 6, byrow = TRUE))

gene_escape_counts <- heat_df %>%
  filter(currentescape == "escaping") %>%
  group_by(external_gene_name) %>%
  summarise(
    n_escaping_participants = n_distinct(participant),
    n_escaping_participant_tissue = n_distinct(paste(participant, tissue, sep = "__")),
    .groups = "drop"
  )

gene_counts_df <- heat_df %>%
  distinct(external_gene_name, escape_group) %>%
  left_join(gene_escape_counts, by = "external_gene_name") %>%
  mutate(
    n_escaping_participants       = tidyr::replace_na(n_escaping_participants, 0L),
    n_escaping_participant_tissue = tidyr::replace_na(n_escaping_participant_tissue, 0L),
    external_gene_name = factor(external_gene_name, levels = levels(heat_df$external_gene_name))
  )

p_part_count <- ggplot(gene_counts_df, aes(y = external_gene_name, x = n_escaping_participants)) +
  geom_col(fill = "steelblue", width = 0.8) +
  facet_grid(escape_group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "none") +
  labs(x = "# participants")

p_pt_count <- ggplot(gene_counts_df, aes(y = external_gene_name, x = n_escaping_participant_tissue)) +
  geom_col(fill = "firebrick", width = 0.8) +
  facet_grid(escape_group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "none") +
  labs(x = "# participant × tissue")

top_panel <- p_heat + p_part_count + p_pt_count + plot_layout(widths = c(20, 3, 3))

silenced_combined_plot <- (top_panel / p_participant / p_age / p_age_group / p_tissue) +
  plot_layout(heights = c(20, 1, 1, 1, 1))
ggsave("silenced_combined_plot.pdf", silenced_combined_plot, height = 14, width = 9)


# -- Variable genes --
heat_df <- m3_var %>%
  semi_join(gene_escape_groups_var, by = "external_gene_name") %>%
  left_join(
    gene_escape_groups_var %>%
      dplyr::select(external_gene_name, escape_group),
    by = "external_gene_name"
  )
length(unique(heat_df$external_gene_name))
heat_df$currentescape = ifelse(heat_df$one_minus_AE > 0.6, 'escaping', 'silenced')

sample_order <- heat_df %>%
  distinct(sampleID, age) %>%
  arrange(age, sampleID) %>%
  pull(sampleID)

gene_order <- heat_df %>%
  distinct(external_gene_name, escape_group) %>%
  arrange(escape_group, external_gene_name) %>%
  pull(external_gene_name)

heat_df <- heat_df %>%
  mutate(
    sampleID = factor(sampleID, levels = sample_order),
    external_gene_name = factor(external_gene_name, levels = rev(gene_order))
  )

sample_ann <- m3_var %>%
  distinct(sampleID, participant, tissue, age) %>%
  filter(sampleID %in% sample_order) %>%
  mutate(sampleID = factor(sampleID, levels = sample_order))
sample_ann$group = ifelse(sample_ann$age >= 50, 'old', 'young')

heat_df <- heat_df %>%
  mutate(currentescape = factor(currentescape, levels = c("silenced", "escaping")))

escaping_genes <- heat_df %>%
  filter(currentescape == "escaping") %>%
  pull(external_gene_name) %>%
  unique()

heat_df <- heat_df %>%
  filter(external_gene_name %in% escaping_genes) %>%
  mutate(external_gene_name = droplevels(external_gene_name))

sample_ann <- sample_ann %>%
  mutate(
    participant = factor(participant),
    tissue      = factor(tissue),
    group       = factor(group)
  )

tissue_levels <- levels(sample_ann$tissue)
tissue_cols <- setNames(
  as.character(Polychrome::createPalette(
    length(tissue_levels),
    seedcolors = c("#FF0000", "#00FF00", "#0000FF")
  )),
  tissue_levels
)

participant_cols <- setNames(
  qualitative_hcl(n = nlevels(sample_ann$participant), palette = "Dynamic"),
  levels(sample_ann$participant)
)

group_cols <- setNames(
  qualitative_hcl(n = nlevels(sample_ann$group), palette = "Set 2"),
  levels(sample_ann$group)
)

p_heat <- ggplot(heat_df, aes(x = sampleID, y = external_gene_name, color = currentescape)) +
  geom_point(size = 2, alpha = 1) +
  facet_grid(escape_group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = currentescape_cols, name = "Escape status", na.value = "grey80") +
  theme_article() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.text.y = element_text(angle = 0)) +
  labs(x = NULL, y = "Gene")

p_participant <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = participant)) +
  geom_tile() +
  scale_fill_manual(values = participant_cols, name = "Participant") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), legend.key.size = unit(0.25, "cm")) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))

p_age <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = age)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Age") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

p_age_group <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = group_cols, name = "Group") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

p_tissue <- ggplot(sample_ann, aes(x = sampleID, y = 1, fill = tissue)) +
  geom_tile() +
  scale_fill_manual(values = tissue_cols, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), legend.key.size = unit(0.25, "cm")) +
  guides(fill = guide_legend(ncol = 6, byrow = TRUE))

gene_escape_counts <- heat_df %>%
  filter(currentescape == "escaping") %>%
  group_by(external_gene_name) %>%
  summarise(
    n_escaping_participants = n_distinct(participant),
    n_escaping_participant_tissue = n_distinct(paste(participant, tissue, sep = "__")),
    .groups = "drop"
  )

gene_counts_df <- heat_df %>%
  distinct(external_gene_name, escape_group) %>%
  left_join(gene_escape_counts, by = "external_gene_name") %>%
  mutate(
    n_escaping_participants       = tidyr::replace_na(n_escaping_participants, 0L),
    n_escaping_participant_tissue = tidyr::replace_na(n_escaping_participant_tissue, 0L),
    external_gene_name = factor(external_gene_name, levels = levels(heat_df$external_gene_name))
  )

p_part_count <- ggplot(gene_counts_df, aes(y = external_gene_name, x = n_escaping_participants)) +
  geom_col(fill = "steelblue", width = 0.8) +
  facet_grid(escape_group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "none") +
  labs(x = "# participants")

p_pt_count <- ggplot(gene_counts_df, aes(y = external_gene_name, x = n_escaping_participant_tissue)) +
  geom_col(fill = "firebrick", width = 0.8) +
  facet_grid(escape_group ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "none") +
  labs(x = "# participant × tissue")

top_panel <- p_heat + p_part_count + p_pt_count + plot_layout(widths = c(20, 3, 3))

variable_combined_plot <- (top_panel / p_participant / p_age / p_age_group / p_tissue) +
  plot_layout(heights = c(20, 1, 1, 1, 1))
ggsave("variable_combined_plot.pdf", variable_combined_plot, height = 8.9, width = 9)

#### end ####

#### 6. matched tissue ####

m3 = subset(dh, m3 == TRUE)
m3$currentescape = ifelse(m3$one_minus_AE >= 0.6, 'escaping', 'silenced')

gene_tissue_counts <- m3 %>%
  group_by(external_gene_name, tissue, group) %>%
  summarise(
    n_participants = n_distinct(participant),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = n_participants,
    values_fill = 0
  ) %>%
  rename(
    n_young = young,
    n_old = old
  )

gene_tissue_counts_use = subset(gene_tissue_counts, n_old > 0 & n_young > 0)
gene_tissue_counts_use$key = paste(gene_tissue_counts_use$external_gene_name, gene_tissue_counts_use$tissue)

m3$key = paste(m3$external_gene_name, m3$tissue)
m3n = subset(m3, key %in% gene_tissue_counts_use$key)

# status per gene x tissue x PAR key across young/old
pl <- m3n %>%
  mutate(group = factor(group, levels = c("young", "old"))) %>%
  group_by(external_gene_name, tissue, key, group, PAR) %>%
  summarise(
    max = max(one_minus_AE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    currentescape = ifelse(max >= 0.6, 'escaping', 'silenced')
  )

# classify switch pattern
switch_df <- pl %>%
  dplyr::select(external_gene_name, key, group, currentescape) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = currentescape
  ) %>%
  mutate(
    switch_pattern = case_when(
      young == old ~ "no switch",
      young == "escaping" & old == "silenced" ~ "escaping in young only",
      young == "silenced" & old == "escaping" ~ "escaping in old only",
      TRUE ~ "incomplete"
    )
  ) %>%
  dplyr::select(key, external_gene_name, switch_pattern)

# Genes that escape in young only
young_only_genes <- switch_df %>%
  filter(switch_pattern == "escaping in young only") %>%
  pull(external_gene_name) %>%
  unique()

# Genes that escape in old only
old_only_genes <- switch_df %>%
  filter(switch_pattern == "escaping in old only") %>%
  pull(external_gene_name) %>%
  unique()

young_only_genes
old_only_genes

pl <- pl %>%
  left_join(switch_df %>% dplyr::select(-external_gene_name), by = "key") %>%
  mutate(
    status_change = switch_pattern != "no switch",
    switch_pattern = factor(
      switch_pattern,
      levels = c(
        "no switch",
        "escaping in young only",
        "escaping in old only",
        "incomplete"
      )
    )
  )

lines_true <- pl %>% filter(status_change == TRUE)
lines_false <- pl %>% filter(status_change == FALSE)

matched_tissue <- ggplot(pl, aes(x = group, y = max)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(PAR ~ tissue) +
  geom_hline(yintercept = 0.6, linetype = "dashed") +
  geom_line(
    data = lines_false,
    aes(group = key),
    color = "gray70",
    linewidth = 0.7,
    alpha = 0.6
  ) +
  geom_line(
    data = lines_true,
    aes(group = key, color = switch_pattern),
    linewidth = 1,
    alpha = 0.8
  ) +
  geom_point(
    aes(color = switch_pattern, shape = currentescape),
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "no switch" = "gray70",
      "escaping in young only" = "darkblue",
      "escaping in old only" = "red",
      "incomplete" = "black"
    )
  ) +
  scale_shape_manual(
    values = c(
      "escaping" = 16,
      "silenced" = 17
    )
  ) +
  theme_article(base_size = 9) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Minimum allelic expression",
    color = "Switch pattern",
    shape = "Status"
  )
ggsave("matched_tissue.pdf", matched_tissue, height = 4, width = 11)

#### end ####

#### 7. variance partitioning ####

get_var_parts <- function(df) {
  
  df <- df %>%
    filter(
      !is.na(allelic_expression),
      !is.na(age),
      !is.na(AE_center),
      !is.na(AE_IQR),
      !is.na(totalCount),
      !is.na(participant),
      !is.na(tissue)
    )
  
  if (
    nrow(df) < 5 ||
    n_distinct(df$participant) < 3 ||
    n_distinct(df$tissue) < 2
  ) {
    return(tibble())
  }
  
  fit <- tryCatch(
    lmer(
      allelic_expression ~
        age +
        AE_center +
        AE_IQR +
        totalCount +
        (1 | participant) +
        (1 | tissue),
      data = df,
      REML = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(tibble())
  
  vc <- as.data.frame(VarCorr(fit))
  
  get_vcov <- function(term) {
    out <- vc %>%
      filter(grp == term) %>%
      pull(vcov)
    
    if (length(out) == 0) NA_real_ else out[1]
  }
  
  var_participant <- get_vcov("participant")
  var_tissue      <- get_vcov("tissue")
  var_residual    <- get_vcov("Residual")
  
  X <- model.matrix(fit)
  beta <- fixef(fit)
  
  fixed_terms <- c("age", "AE_center", "AE_IQR", "totalCount")
  
  fixed_vars <- sapply(fixed_terms, function(term) {
    if (term %in% colnames(X) && term %in% names(beta)) {
      var(X[, term] * beta[term], na.rm = TRUE)
    } else {
      NA_real_
    }
  })
  
  total_var <- sum(
    fixed_vars,
    var_participant,
    var_tissue,
    var_residual,
    na.rm = TRUE
  )
  
  tibble(
    predictor = c(
      fixed_terms,
      "participant",
      "tissue",
      "residual"
    ),
    variance = c(
      fixed_vars,
      var_participant,
      var_tissue,
      var_residual
    ),
    variance_fraction = variance / total_var,
    gene_group = first(na.omit(df$pl2))
  )
}

# Run variance part all samples  
varpart_df_all <- dh %>%
  group_by(external_gene_name) %>%
  group_modify(~ get_var_parts(.x)) %>%
  ungroup()

ordered_predictors <- c("AE_center", "participant", "tissue", "age", "totalCount", "AE_IQR", "residual")

varpart_df_all$predictor <- factor(varpart_df_all$predictor, levels = ordered_predictors)
base_pal <- as.character(paletteer::paletteer_d("MetBrewer::Archambault"))
pal <- colorRampPalette(base_pal)(length(ordered_predictors))
names(pal) <- ordered_predictors

outliers <- varpart_df_all %>%
  group_by(predictor) %>%
  mutate(
    q1 = quantile(variance_fraction, 0.25, na.rm = TRUE),
    q3 = quantile(variance_fraction, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr
  ) %>%
  ungroup() %>%
  filter(variance_fraction < lower | variance_fraction > upper)

mean_df <- varpart_df_all %>%
  group_by(predictor) %>%
  summarise(
    mean = mean(variance_fraction, na.rm = TRUE),
    .groups = "drop"
  )

varpart_all <- ggplot(varpart_df_all, aes(x = predictor, y = variance_fraction, fill = predictor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = pal) +
  geom_point(
    data = outliers,
    aes(shape = gene_group),
    size = 2,
    alpha = 0.9
  ) +
  geom_text(
    data = mean_df,
    aes(x = predictor, y = -0.02, label = sprintf("%.3f", mean)),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_shape_manual(
    values = c(
      "PAR"      = 15,
      "escape"   = 16,
      "other"    = 4,
      "silenced" = 18,
      "variable" = 17
    ),
    name = "Gene group"
  ) +
  guides(fill = "none") +
  theme_bw() +
  ylab("Proportion of variance explained") +
  xlab(NULL) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    legend.position = "right"
  )
ggsave("varpart_all.pdf", varpart_all, height = 7, width = 12)

# Run variance part for m3 samples only 
varpart_df_m3 <- dh %>%
  filter(m3 == TRUE) %>%
  group_by(external_gene_name) %>%
  group_modify(~ get_var_parts(.x)) %>%
  ungroup()

ordered_predictors <- c("tissue", "participant", "AE_IQR", "AE_center", "age", "totalCount", "residual")

varpart_df_m3$predictor <- factor(varpart_df_m3$predictor, levels = ordered_predictors)
base_pal <- as.character(paletteer::paletteer_d("MetBrewer::Archambault"))
pal <- colorRampPalette(base_pal)(length(ordered_predictors))
names(pal) <- ordered_predictors

outliers <- varpart_df_m3 %>%
  group_by(predictor) %>%
  mutate(
    q1 = quantile(variance_fraction, 0.25, na.rm = TRUE),
    q3 = quantile(variance_fraction, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr
  ) %>%
  ungroup() %>%
  filter(variance_fraction < lower | variance_fraction > upper)

mean_df <- varpart_df_m3 %>%
  group_by(predictor) %>%
  summarise(
    mean = mean(variance_fraction, na.rm = TRUE),
    .groups = "drop"
  )

varpart_m3 <- ggplot(varpart_df_m3, aes(x = predictor, y = variance_fraction, fill = predictor)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = pal) +
  geom_point(
    data = outliers,
    aes(shape = gene_group),
    size = 2,
    alpha = 0.9
  ) +
  geom_text(
    data = mean_df,
    aes(x = predictor, y = -0.02, label = sprintf("%.3f", mean)),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_shape_manual(
    values = c(
      "PAR"      = 15,
      "escape"   = 16,
      "other"    = 4,
      "silenced" = 18,
      "variable" = 17
    ),
    name = "Gene group"
  ) +
  guides(fill = "none") +
  theme_bw() +
  ylab("Proportion of variance explained") +
  xlab(NULL) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    legend.position = "right"
  )
ggsave("varpart_m3.pdf", varpart_m3, height = 7, width = 12)

#### end ####

#### 8. linear modelling ####

fit_gene_models <- function(data, cutoff,
                            use_AE_center, use_IQR, use_totalCount,
                            use_participant_re, use_tissue_re) {
  
  # fixed-effects 
  rhs_fixed <- "age"
  if (use_AE_center)  rhs_fixed <- c(rhs_fixed, "AE_center")
  if (use_IQR)        rhs_fixed <- c(rhs_fixed, "AE_IQR")
  if (use_totalCount) rhs_fixed <- c(rhs_fixed, "totalCount")
  
  # random-effects
  re_terms <- c(
    if (use_participant_re) "(1 | participant)",
    if (use_tissue_re)      "(1 | tissue)"
  )
  use_lmer <- length(re_terms) > 0
  
  # apply cutoff
  data <- data[get(cutoff) == TRUE]
  
  data[!is.na(one_minus_AE)] |>
    split(by = "external_gene_name") |>
    lapply(function(g) {
      
      gene       <- g$external_gene_name[1]
      n_obs      <- nrow(g)
      n_ind      <- length(unique(g$participant))
      n_ind_tiss <- nrow(unique(g[, .(participant, tissue)]))
      
      if (n_ind <= 3) return(data.table(
        external_gene_name = gene, n_obs = n_obs, n_ind = n_ind,
        n_ind_tiss = n_ind_tiss, model = "skipped_n_ind<=3", term = "age"
      ))
      
      g[, one_minus_AE := scale(one_minus_AE)]
      g[, age          := scale(age)]
      if (use_AE_center)  g[, AE_center  := scale(AE_center)]
      if (use_IQR)        g[, AE_IQR     := scale(AE_IQR)]
      if (use_totalCount) g[, totalCount := scale(totalCount)]
      
      fml <- paste(
        "one_minus_AE ~",
        paste(rhs_fixed, collapse = " + "),
        if (use_lmer) paste("+", paste(re_terms, collapse = " + "))
      )
      
      fit <- tryCatch(
        if (use_lmer) lmerTest::lmer(as.formula(fml), data = g, REML = FALSE)
        else          lm(as.formula(fml), data = g),
        error = function(e) NULL
      )
      
      if (is.null(fit)) return(data.table(
        external_gene_name = gene, n_obs = n_obs, n_ind = n_ind,
        n_ind_tiss = n_ind_tiss, model = "failed", term = "age"
      ))
      
      model_note <- if (use_lmer) {
        if (isSingular(fit)) "lmer_singular" else "lmer"
      } else "lm"
      
      res <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
      as.data.table(res)[, `:=`(
        external_gene_name = gene,
        n_obs              = n_obs,
        n_ind              = n_ind,
        n_ind_tiss         = n_ind_tiss,
        model              = model_note,
        BIC                = BIC(fit),
        AIC                = AIC(fit)
      )]
    }) |>
    rbindlist(fill = TRUE)
}

# ── Factorial grid ─────────────────────────────────────────────
pl_lookup <- unique(d[, .(external_gene_name, pl2)])

grid <- expand.grid(
  cutoff             = c("m25", "m3"),
  use_AE_center      = c(TRUE, FALSE),
  use_IQR            = c(TRUE, FALSE),
  use_participant_re = c(TRUE, FALSE),
  use_tissue_re      = c(TRUE, FALSE),
  use_totalCount     = c(TRUE, FALSE),
  stringsAsFactors   = FALSE
)

all_results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  cfg <- grid[i, ]
  cat(sprintf(
    "[%d/%d] cutoff=%-3s  AE_center=%s  IQR=%s  participant_re=%s  tissue_re=%s  totalCount=%s\n",
    i, nrow(grid),
    cfg$cutoff, cfg$use_AE_center, cfg$use_IQR,
    cfg$use_participant_re, cfg$use_tissue_re, cfg$use_totalCount
  ))
  
  raw <- fit_gene_models(
    data               = copy(dh),
    cutoff             = cfg$cutoff,
    use_AE_center      = cfg$use_AE_center,
    use_IQR            = cfg$use_IQR,
    use_participant_re = cfg$use_participant_re,
    use_tissue_re      = cfg$use_tissue_re,
    use_totalCount     = cfg$use_totalCount
  )
  
  res <- raw[grepl("^age$", term) | model %in% c("failed", "skipped_n_ind<=3")]
  
  res[, p.adj := NA_real_]
  valid <- !is.nan(res$p.value) & !is.na(res$p.value)
  res[valid, p.adj := p.adjust(p.value, method = "BH")]
  setorder(res, p.adj)
  
  res <- merge(res, pl_lookup, by = "external_gene_name", all.x = TRUE)
  
  res[, `:=`(
    cutoff             = cfg$cutoff,
    use_AE_center      = cfg$use_AE_center,
    use_IQR            = cfg$use_IQR,
    use_participant_re = cfg$use_participant_re,
    use_tissue_re      = cfg$use_tissue_re,
    use_totalCount     = cfg$use_totalCount,
    direction          = fcase(
      p.value < 0.05 & estimate > 0, "more escaped with age",
      p.value < 0.05 & estimate < 0, "more silenced with age",
      default = "NS"
    )
  )]
  
  all_results[[i]] <- res
}

results <- rbindlist(all_results, fill = TRUE)


# AIC rank per gene x cutoff
setDT(results)
results[!is.na(AIC) & is.finite(AIC),
        aic_rank := rank(AIC, ties.method = "min"),
        by = .(external_gene_name, cutoff)]

results[, delta_aic := AIC - min(AIC[is.finite(AIC)], na.rm = TRUE),
        by = .(external_gene_name, cutoff)]

results[, best_aic := fcase(
  delta_aic <= 2, TRUE,
  delta_aic > 2, FALSE,
  default = NA)]

write_csv(results, "age_escape_results.csv")

#### 8.1 shape by gene category ####

m3_res <- results %>%
  filter(cutoff == "m3")

df_plot <- m3_res %>%
  filter(!is.na(estimate)) %>%
  group_by(external_gene_name) %>%
  filter(any(p.value < 0.05)) %>%
  ungroup() %>%
  mutate(
    sig = factor(
      ifelse(p.value < 0.05, "p < 0.05", "p >= 0.05"),
      levels = c("p < 0.05", "p >= 0.05")
    ),
    pl2 = factor(pl2)
  )

gene_order <- df_plot %>%
  group_by(external_gene_name) %>%
  summarise(max_est = max(estimate, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_est)) %>%
  pull(external_gene_name)

df_plot <- df_plot %>%
  mutate(external_gene_name = factor(external_gene_name, levels = rev(gene_order)))

df_plot$conf.low  <- df_plot$estimate - df_plot$std.error
df_plot$conf.high <- df_plot$estimate + df_plot$std.error

df_plot <- df_plot %>%
  filter(!is.na(p.value))

m3_p <- ggplot(df_plot, aes(
  x     = estimate,
  y     = external_gene_name,
  xmin  = conf.low,
  xmax  = conf.high,
  color = best_aic == "TRUE" & !is.na(best_aic),
  alpha = sig,
  shape = pl2
)) +
  geom_pointrange(data = ~ subset(., is.na(best_aic) | best_aic != "TRUE")) +
  geom_pointrange(data = ~ subset(., !is.na(best_aic) & best_aic == "TRUE")) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    name   = "Best model"
  ) +
  scale_alpha_manual(
    values = c("p >= 0.05" = 0.15, "p < 0.05" = 1),
    name   = "Significant"
  ) +
  scale_shape_manual(
    values = c("PAR" = 15, "escape" = 16, "other" = 4, "silenced" = 18, "variable" = 17),
    name = "Gene group"
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Coefficient", y = "Gene", title = "m3") +
  theme_bw()
ggsave("32models_m3_by_shape.pdf", m3_p, height = 10, width = 6.5)

#### end ####

#### 9. consistency ####

# ── 1. Summarise each analysis per gene ──────────────────────────────────────

# Analysis 1: among best-AIC models, proportion with sig positive age effect
a1 <- as.data.table(results)[best_aic == TRUE,
                             .(prop_sig_pos = mean(estimate > 0 & p.value < 0.05, na.rm = TRUE),
                               prop_sig_neg = mean(estimate < 0 & p.value < 0.05, na.rm = TRUE),
                               n_best       = .N),
                             by = external_gene_name]
a1[, a1_direction := fcase(
  prop_sig_pos > 0 & prop_sig_neg == 0, "more escape with age",
  prop_sig_neg > 0 & prop_sig_pos == 0, "less escape with age",
  prop_sig_pos > 0 & prop_sig_neg > 0,  "mixed", # some best AIC go one way, some go the other
  default = "no effect" #no best AIC has a significant effect in either direction
)]

# a1 based on all models
a1_all <- as.data.table(results)[p.value < 0.05,
                                 .(prop_sig_pos_all = mean(estimate > 0, na.rm = TRUE),
                                   prop_sig_neg_all = mean(estimate < 0, na.rm = TRUE),
                                   n_all            = .N),
                                 by = external_gene_name]
a1_all[, a1_direction_all := fcase(
  prop_sig_pos_all > 0 & prop_sig_neg_all == 0, "more escape with age",
  prop_sig_neg_all > 0 & prop_sig_pos_all == 0, "less escape with age",
  prop_sig_pos_all > 0 & prop_sig_neg_all > 0,  "mixed",
  default = "no effect"
)]

# Analysis 2: across tissues, any "escaping in old only" or "escaping in young only"
pl_dt <- as.data.table(pl)
a2 <- pl_dt[, .(
  n_old_only   = sum(switch_pattern == "escaping in old only",   na.rm = TRUE),
  n_young_only = sum(switch_pattern == "escaping in young only", na.rm = TRUE),
  n_no_switch  = sum(switch_pattern == "no switch",              na.rm = TRUE)
), by = external_gene_name]

a2[, a2_direction := fcase(
  n_old_only > 0 & n_young_only == 0, "more escape with age",
  n_young_only > 0 & n_old_only == 0, "less escape with age",
  n_old_only > 0 & n_young_only > 0,  "mixed",
  default = "no switch"
)]

# Analysis 3: gene_escape_groups_sil label
a3 <- as.data.table(gene_escape_groups_sil)[, .(
  external_gene_name = external_gene_name,
  a3_direction = fcase(
    escape_group == "escaping in old",          "more escape with age",
    escape_group == "escaping in young",        "less escape with age",
    escape_group == "escaping in old,young","mixed",
    default = "silenced"
  )
)]

# Analysis 4: gene_escape_groups_var label
a4 <- as.data.table(gene_escape_groups_var)[, .(
  external_gene_name = external_gene_name,
  a4_direction = fcase(
    escape_group == "escaping in old",          "more escape with age",
    escape_group == "escaping in young",        "less escape with age",
    escape_group == "escaping in old,young","mixed",
    default = "silenced"
  )
)]

# ── 2. Merge all four ────────────────────────────────────────────────────────
consistency_df <- Reduce(
  function(x, y) merge(x, y, by = "external_gene_name", all = TRUE),
  list(
    as.data.frame(a1)[, c("external_gene_name", "a1_direction", "prop_sig_pos")],
    as.data.frame(a1_all)[, c("external_gene_name", "a1_direction_all")],
    as.data.frame(a2)[, c("external_gene_name", "a2_direction")],
    a3,
    a4
  )
)
consistency_df %>%
  dplyr::select(external_gene_name, a1_direction, a1_direction_all, a2_direction, a3_direction, a4_direction) %>%
  write_csv("analysis_consistency.csv")


# ── 3. Count genes showing more/less escape with age ─────────────────────────
more_escape <- consistency_df %>%
  filter(if_any(c(a1_direction, a2_direction, a3_direction, a4_direction),
                ~ . == "more escape with age")) %>%
  pull(external_gene_name) %>%
  unique()

less_escape <- consistency_df %>%
  filter(if_any(c(a1_direction, a2_direction, a3_direction, a4_direction),
                ~ . == "less escape with age")) %>%
  pull(external_gene_name) %>%
  unique()

# any mixed sign across analyses
any_mixed <- consistency_df %>%
  filter(if_any(c(a1_direction, a2_direction, a3_direction, a4_direction),
                ~ . == "mixed")) %>%
  pull(external_gene_name) %>%
  unique()

# inconsistent in direction across analyses
inconsistent <- intersect(more_escape, less_escape)

# any mixed OR inconsistent direction
both_escape <- unique(c(any_mixed, inconsistent))

# remove from more and less escape
more_escape <- setdiff(more_escape, both_escape)
less_escape <- setdiff(less_escape, both_escape)

length(more_escape)   
length(less_escape)   
length(both_escape) 

# ── 4. Upset plot: overlap of genes showing "more escape with age" ────────────
upset_dt <- consistency_df %>%
  filter(external_gene_name %in% more_escape) %>%
  mutate(
    A1 = a1_direction == "more escape with age" & !is.na(a1_direction),
    A2 = a2_direction == "more escape with age" & !is.na(a2_direction),
    A3 = a3_direction == "more escape with age" & !is.na(a3_direction),
    A4 = a4_direction == "more escape with age" & !is.na(a4_direction)
  ) %>%
  dplyr::select(external_gene_name, A1, A2, A3, A4)

upset_dt <- as.data.table(upset_dt)
upset_dt[, analyses := lapply(seq_len(.N), function(i) {
  v <- c()
  if (A1[i]) v <- c(v, "A1: lmer")
  if (A2[i]) v <- c(v, "A2: tissue switch")
  if (A3[i]) v <- c(v, "A3: group escape")
  if (A4[i]) v <- c(v, "A4: group escape variable")
  v
})]

upset_dt_sub <- upset_dt %>%
  filter(!(A1 == "FALSE" & A2 == "FALSE" & A3 == "FALSE" & A4 == "FALSE"))

upset_p <- ggplot(upset_dt_sub, aes(x = analyses)) +
  geom_bar(fill = "#d62728", alpha = 0.8) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 3.5) +
  scale_x_upset(n_intersections = 20) +
  scale_y_continuous(limits = c(0, 20)) +
  labs(x = NULL, y = "Number of genes") +
  theme_bw(base_size = 10)
ggsave("upset_plot_4groups.pdf", upset_p, height = 4, width = 11)

#### end ####

#### 10. positional enrichment ####

# ── 0. Background: all X chromosome genes ────────────────────
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 115)

x_genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position"),
  filters    = "chromosome_name",
  values     = "X",
  mart       = ensembl
) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  dplyr::select(external_gene_name = hgnc_symbol, start_position, end_position) %>%
  mutate(chromosome_name = "chrX")

# ── 1. Bin setup ─────────────────────────────────────────────
bin_size      <- 10e6
chrX_length   <- 156040895
PAR1_boundary <- 2.8  # Mb

bins <- data.frame(
  start = seq(1, chrX_length, by = bin_size),
  end   = pmin(seq(1, chrX_length, by = bin_size) + bin_size - 1, chrX_length)
)
bins$bin_id <- paste0("bin_", seq_len(nrow(bins)))

assign_bin <- function(pos) {
  i <- which(bins$start <= pos & bins$end >= pos)
  if (length(i) == 1) return(bins$bin_id[i])
  return(NA_character_)
}

# ── 2. Background bin counts ──────────────────────────────────
bg <- x_genes %>%
  filter(!is.na(start_position), !is.na(end_position)) %>%
  mutate(
    midpoint = (start_position + end_position) / 2,
    bin      = sapply(midpoint, assign_bin)
  )

bg_counts   <- bg %>% group_by(bin) %>% summarise(genes_in_bin = n(), .groups = "drop")
genes_total <- nrow(bg)

# ── 3. Fisher's enrichment helper ────────────────────────────
run_enrichment <- function(gene_vec, label) {
  hits <- bg %>% filter(external_gene_name %in% gene_vec)
  hits_total <- nrow(hits)
  if (hits_total == 0) { message("No genes found for: ", label); return(NULL) }
  
  hit_counts <- hits %>% group_by(bin) %>% summarise(hits_in_bin = n(), .groups = "drop")
  
  left_join(bg_counts, hit_counts, by = "bin") %>%
    mutate(hits_in_bin = replace_na(hits_in_bin, 0)) %>%
    rowwise() %>%
    mutate(
      mat = list(matrix(c(
        hits_in_bin,
        genes_in_bin - hits_in_bin,
        hits_total   - hits_in_bin,
        genes_total  - genes_in_bin - (hits_total - hits_in_bin)
      ), nrow = 2, byrow = TRUE)),
      fisher_p   = fisher.test(mat, alternative = "greater")$p.value,
      odds_ratio = fisher.test(mat, alternative = "greater")$estimate,
      binom_p    = binom.test(
        x           = hits_in_bin,
        n           = hits_total,
        p           = genes_in_bin / genes_total,
        alternative = "greater"
      )$p.value
    ) %>%
    ungroup() %>%
    dplyr::select(-mat) %>%
    mutate(
      neg_log10_p       = -log10(fisher_p),
      neg_log10_binom_p = -log10(binom_p),
      sig_fisher        = fisher_p < 0.05,
      sig_binom         = binom_p  < 0.05,
      sig               = fisher_p < 0.05,
      label             = label,
      n_hits            = hits_total
    ) %>%
    left_join(bins, by = c("bin" = "bin_id"))
}

# ── 4. Gene sets ─────────────────────────────────────────────
gene_sets <- list(
  union_more = more_escape,
  union_less = less_escape,
  union_all  = unique(c(more_escape, less_escape))
)

set_labels <- c(
  union_more = "Union: More escape with age",
  union_less = "Union: Less escape with age",
  union_all  = "Union: All XCI-changing"
)

# ── 5. Run enrichment tests ───────────────────────────────────
all_results <- mapply(
  run_enrichment,
  gene_vec = gene_sets,
  label    = set_labels[names(gene_sets)],
  SIMPLIFY = FALSE
) %>%
  Filter(Negate(is.null), .) %>%
  bind_rows()

# ── 6. Plot function ──────────────────────────────────────────
plot_enrichment <- function(df, title_str) {
  ggplot(df) +
    geom_rect(aes(
      xmin = start / 1e6,
      xmax = end   / 1e6,
      ymin = 0,
      ymax = neg_log10_p,
      fill = sig
    )) +
    scale_fill_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
      guide  = "none"
    ) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_vline(xintercept = PAR1_boundary,
               linetype = "dotted", color = "black", linewidth = 0.5) +
    annotate("text", x = PAR1_boundary, y = Inf,
             label = "PAR1", hjust = -0.15, vjust = 1.5, size = 2.5, color = "gray30") +
    labs(
      title    = title_str,
      subtitle = paste0("n = ", unique(df$n_hits), " genes"),
      x        = "Position on chrX (Mb)",
      y        = expression(-log[10](italic(p)))
    ) +
    scale_x_continuous(breaks = seq(0, 155, by = 20), limits = c(0, 157)) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(size = 9, face = "bold"),
      plot.subtitle    = element_text(size = 8, color = "gray40"),
      axis.title       = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
}

# ── 7. Union panel ────────────────────────────────────────────
panel_union <- (
  plot_enrichment(all_results %>% filter(label == set_labels["union_more"]), set_labels["union_more"]) /
    plot_enrichment(all_results %>% filter(label == set_labels["union_less"]), set_labels["union_less"]) /
    plot_enrichment(all_results %>% filter(label == set_labels["union_all"]),  set_labels["union_all"])
) +
  plot_annotation(
    title    = "Positional Enrichment of XCI Changes with Age — Union Across Analyses",
    subtitle = "Filled = p < 0.05 | Dashed = significance threshold | Dotted = PAR1 boundary",
    theme    = theme(plot.title = element_text(size = 12, face = "bold"))
  )
ggsave("pos_enrichment.pdf", panel_union, height = 6, width = 4.5)

#### end ####

#### 12.1 Functional enrichment ####
library(gprofiler2)
library(ggplot2)
library(dplyr)
library(stringr)

# ── 1. Gene sets ─────────────────────────────────────────────
genes_aged  <- more_escape
genes_young <- less_escape

cat("Genes escaping in aged: ",  length(genes_aged),  "\n")
cat("Genes escaping in young:", length(genes_young), "\n")
# ── 2. Run gProfiler ─────────────────────────────────────────
sources_of_interest <- c("CORUM", "GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "MIRNA", "TF", "HP")

run_gost_paper <- function(gene_vec, label) {
  if (length(gene_vec) == 0) { message("Empty set: ", label); return(NULL) }
  message("Running gProfiler: ", label, " (n = ", length(gene_vec), ")")
  res <- gost(
    query             = gene_vec,
    organism          = "hsapiens",
    ordered_query     = FALSE,
    correction_method = "g_SCS",
    user_threshold    = 0.1,
    sources           = sources_of_interest,
    evcodes           = TRUE
  )
  if (is.null(res) || nrow(res$result) == 0) {
    message("No significant results for: ", label); return(NULL)
  }
  res$result %>% mutate(group = label)
}

gp_aged  <- run_gost_paper(genes_aged,  "Aged")
gp_young <- run_gost_paper(genes_young, "Young")

gp_all <- bind_rows(gp_aged, gp_young) %>%
  mutate(
    group          = factor(group, levels = c("Young", "Aged")),
    neg_log10_p    = -log10(p_value),
    term_name_wrap = str_wrap(term_name, width = 15)
  )

# ── 3. Filter sources and order terms ────────────────────────
sources_show <- c("CORUM", "GO:BP", "GO:CC", "GO:MF", "HP", "MIRNA")
gp_plot <- gp_all %>%
  filter(source %in% sources_show) %>%
  mutate(source = factor(source, levels = sources_show))

term_order <- gp_plot %>%
  group_by(source, term_name_wrap) %>%
  summarise(min_p = min(p_value), .groups = "drop") %>%
  arrange(source, min_p) %>%
  pull(term_name_wrap)

gp_plot <- gp_plot %>%
  mutate(term_name_wrap = factor(term_name_wrap, levels = unique(term_order)))

# ── 4. Plot ───────────────────────────────────────────────────
fig4_b <- ggplot(gp_plot, aes(
  x     = term_name_wrap,
  y     = group,
  size  = intersection_size,
  color = p_value
)) +
  geom_point(alpha = 0.9) +
  facet_grid(
    . ~ source,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_color_gradient(
    low      = "#08519c",
    high     = "#9ecae1",
    name     = "g:SCS",
    limits   = c(0, 0.1),
    breaks   = c(0.1, 0.05, 0.01),
    labels   = c("0.1", "0.05", "0.01")
  ) +
  scale_size_continuous(
    name   = "Number of Genes",
    range  = c(3, 14),
    breaks = c(5, 10, 20, 40, 60)
  ) +
  scale_y_discrete(limits = c("Aged", "Young")) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_text(size = 8, angle = 0, hjust = 0.5, vjust = 1,
                                      lineheight = 0.85),
    axis.text.y        = element_text(size = 11),
    strip.background   = element_rect(fill = "gray92", color = NA),
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.box         = "horizontal",
    panel.spacing.x    = unit(0.3, "lines")
  ) +
  guides(
    size  = guide_legend(order = 1, override.aes = list(color = "gray30")),
    color = guide_colorbar(order = 2, barwidth = 5, barheight = 0.8,
                           title.position = "top", direction = "horizontal",
                           reverse = TRUE)
  )
ggsave("functional_enrichment_gSCS_0.1_new.pdf", fig4_b, width = 14, height = 4)
#### End ####

#### 13 position on X chromosome ####

genes_export <- bg %>%
  filter(external_gene_name %in% c(more_escape, less_escape, both_escape)) %>%
  mutate(
    direction = case_when(
      external_gene_name %in% both_escape   ~ "escape in both",
      external_gene_name %in% more_escape   ~ "more escape with age",
      external_gene_name %in% less_escape   ~ "less escape with age"
    )
  ) %>%
  arrange(start_position) %>%
  mutate(
    midpoint = (start_position + end_position) / 2,
    x_label  = seq(1, 156040895, length.out = n())
  )
write_csv(genes_export, "genes_export.csv")

all_x <- bg %>%
  mutate(midpoint = (start_position + end_position) / 2)

consensus <- bg %>%
  filter(external_gene_name %in% c(silenced_either, escape_either)) %>%
  mutate(
    status = case_when(
      external_gene_name %in% silenced_either & external_gene_name %in% escape_either ~ "both",
      external_gene_name %in% silenced_either ~ "silenced",
      external_gene_name %in% escape_either   ~ "escape"
    )
  )
write.csv(consensus, "consensus_genes.csv", row.names = FALSE)

# ── 2. Colors ─────────────────────────────────────────────────
dir_colors <- c(
  "more escape with age" = "#d62728",
  "less escape with age" = "#1f77b4",
  "escape in both" = "orange"
)
status_colors <- c(
  "silenced" = "#2ca02c",
  "escape" = "#ff69b4"
)

marker_colors <- dir_colors[genes_123$direction]

#  ── 3. GRanges objects ────────────────────────────────────────
genes_gr <- regioneR::toGRanges(data.frame(
  chr   = "chrX",
  start = genes_123$start_position,
  end   = genes_123$end_position
))
genes_gr$gene_name <- genes_123$external_gene_name
genes_gr$direction <- genes_123$direction

all_x_gr <- regioneR::toGRanges(data.frame(
  chr   = "chrX",
  start = all_x$start_position,
  end   = all_x$end_position
))

all_x <- all_x %>%
  mutate(midpoint = (start_position + end_position) / 2)

consensus_silenced <- subset(consensus, status == "silenced" & !is.na(start_position)) %>%
  mutate(midpoint = (start_position + end_position) / 2)

consensus_escape <- subset(consensus, status == "escape" & !is.na(start_position)) %>%
  mutate(midpoint = (start_position + end_position) / 2)

# ── 4. Plot ───────────────────────────────────────────────────
pdf("gene_positions_X.pdf", width = 20, height = 6)

kp <- plotKaryotype(
  genome      = "hg38",
  chromosomes = "chrX",
  plot.type   = 3,
  zoom        = "chrX:1-156040895"
)

# PAR1 boundary
kpAbline(kp, v = 2800000, col = "black", lty = 2, lwd = 1)

# Vertical stem from true position to elbow
kpSegments(kp,
           chr  = "chrX",
           x0   = genes_123$midpoint,
           x1   = genes_123$midpoint,
           y0   = 0.05,
           y1   = 0.3,
           col  = marker_colors,
           lwd  = 0.8,
           r0   = 0, r1 = 1
)

# Diagonal from elbow to evenly spaced label position
kpSegments(kp,
           chr  = "chrX",
           x0   = genes_123$midpoint,
           x1   = genes_123$x_label,
           y0   = 0.3,
           y1   = 0.55,
           col  = marker_colors,
           lwd  = 0.8,
           r0   = 0, r1 = 1
)

# Gene labels at evenly spaced positions
kpText(kp,
       chr    = "chrX",
       x      = genes_123$x_label,
       y      = 0.6,
       col    = marker_colors,
       labels = genes_123$external_gene_name,
       srt    = 90,
       cex    = 0.5,
       r0     = 0, r1 = 1
)

# ── data.panel = 2: bottom tracks ────────────────────────────

# All coding & noncoding X genes
kpSegments(kp,
           chr        = "chrX",
           x0         = all_x$midpoint,
           x1         = all_x$midpoint,
           y0         = 0, y1 = 1,
           col        = "grey60",
           lwd        = 0.5,
           r0         = 0.55, r1 = 1,
           data.panel = 2
)

# Consensus silenced genes
kpSegments(kp,
           chr        = "chrX",
           x0         = consensus_silenced$midpoint,
           x1         = consensus_silenced$midpoint,
           y0         = 0, y1 = 1,
           col        = status_colors["silenced"],
           lwd        = 0.5,
           r0         = 0, r1 = 0.45,
           data.panel = 2
)

# Consensus escape genes
kpSegments(kp,
           chr        = "chrX",
           x0         = consensus_escape$midpoint,
           x1         = consensus_escape$midpoint,
           y0         = 0, y1 = 1,
           col        = status_colors["escape"],
           lwd        = 0.5,
           r0         = 0, r1 = 0.45,
           data.panel = 2
)

# Row labels
kpAddLabels(kp,
            labels     = "all coding &\nnoncoding genes",
            r0         = 0.55, r1 = 1,
            data.panel = 2,
            cex        = 0.6
)
kpAddLabels(kp,
            labels     = "consensus escape",
            r0         = 0, r1 = 0.45,
            data.panel = 2,
            cex        = 0.6
)

kpAddBaseNumbers(kp)

legend("bottomright",
       legend = c(names(dir_colors), names(status_colors)),
       fill   = c(dir_colors, status_colors),
       border = NA,
       bty    = "n",
       cex    = 0.8
)

dev.off()

#### end ####



