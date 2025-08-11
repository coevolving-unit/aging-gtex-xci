library(dplyr)
library(ggplot2)
library(egg)
library(karyoploteR)
library(regioneR)
library(scales)
library(dplyr)
library(GenomicRanges)
library(biomaRt)

`%!in%` = Negate(`%in%`)

#### load ####

# UPIC = young
# nmXCI-1 = 13PLJ = aged
# nmXCI-2 = ZZPU = aged

# hsap = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
# map = getBM(mart = hsap, attributes = c('ensembl_gene_id', 'external_gene_name','chromosome_name','start_position','end_position'))
# saveRDS(map, file = 'map.rds')

map = readRDS('map.rds')

ase = read.csv('gylemo.csv')
colnames(ase)[3] = 'external_gene_name'
colnames(ase)[6] = 'tissue'
ase = merge(ase, map, by = 'external_gene_name', all.x = T) 
ase$tissue = as.factor(ase$tissue)
levels(ase$tissue) = c('AdiposeSubcutaneous','AdiposeVisceralOmentum','AdrenalGland','ArteryAorta',
                       'ArteryCoronary','ArteryTibial','BrainAnteriorcingulatecortexBA24','BrainCaudatebasalganglia',
                       'BrainCerebellarHemisphere','BrainSpinalcordcervicalc1','Breast','ColonTransverse','EsophagusMucosa','EsophagusMuscularis',
                       'HeartAtrialAppendage','HeartLeftVentricle','KidneyCortex','Liver','Lung','MuscleSkeletal',
                       'Ovary','Pancreas','SkinSunExposedLowerleg','SkinNotSunExposedSuprapubic','Spleen','Stomach',
                       'Thyroid','Uteris','Vagina','WholeBlood')
ase$current_escape = ase$new_category

## limit to shared tissues

table(ase$individual, ase$tissue)

upic = data.frame(table(ase$individual, ase$tissue))
upic = subset(upic, Var1 == 'UPIC' & Freq > 0)

tissue_indiv_counts <- ase %>%
  distinct(individual, tissue) %>%
  group_by(tissue) %>%
  summarise(n_individuals = n_distinct(individual)) %>%
  ungroup()
View(tissue_indiv_counts)

keep = subset(tissue_indiv_counts, tissue %in% upic$Var2 & n_individuals > 1)$tissue
keep

asenow = subset(ase, tissue %in% keep)
asenow$group = ifelse(asenow$individual == 'UPIC', 'young', 'old')
asenow$key = paste(asenow$external_gene_name, asenow$tissue)
asenow$status = ifelse(asenow$allelic_expression<0.4, 'escape', 'silenced')

#### end ####

#### compare age groups ####

asenow_escape = unique(asenow[,c('key','tissue','allelic_expression','external_gene_name','status','group','PAR')])

genes_in_both <- asenow_escape %>%
  group_by(tissue, external_gene_name) %>%
  summarise(n_individuals = n_distinct(group), .groups = "drop") %>%  # assuming 'group' indicates individual group (e.g., young/old)
  filter(n_individuals == 2) %>%  # keep only genes present in BOTH individuals for that tissue
  dplyr::select(tissue, external_gene_name)

asenow_filtered <- asenow_escape %>%
  inner_join(genes_in_both, by = c("tissue", "external_gene_name"))

unique_gene_counts <- asenow_filtered %>%
  distinct(tissue, external_gene_name) %>%
  count(tissue, name = "n_genes") %>%
  arrange(desc(n_genes))

print(unique_gene_counts)

status_switch_genes <- asenow_filtered %>%
  group_by(tissue, external_gene_name) %>%
  summarise(n_status = n_distinct(status), .groups = "drop") %>%
  filter(n_status > 1)

asenow_filtered <- asenow_filtered %>%
  mutate(status_change = paste(tissue, external_gene_name) %in%
           paste(status_switch_genes$tissue, status_switch_genes$external_gene_name))
asenow_filtered$group = factor(asenow_filtered$group, levels = c('young','old'))
levels(asenow_filtered$group)

lines_true <- asenow_filtered %>% filter(status_change == TRUE)
lines_false <- asenow_filtered %>% filter(status_change == FALSE)
ggplot(asenow_filtered, aes(x = group, y = 1 - allelic_expression)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(PAR~tissue) +
  geom_hline(yintercept = 1 - 0.4, linetype = "dashed") +
  geom_line(data = lines_true, aes(group = key), color = "red", linewidth = 0.5, alpha = 0.7) +
  geom_line(data = lines_false, aes(group = key), color = "gray70", linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(color = status, shape = status_change), size = 2) +
  scale_color_manual(values = c("escape" = "darkblue", "silenced" = "darkgrey")) +
  labs(y = "1 - Allelic Expression", x = "Group") +
  theme_article()

library(tidyr)
status_compare <- asenow_filtered %>%
  dplyr::select(external_gene_name, tissue, group, status) %>%
  distinct() %>%
  pivot_wider(names_from = group, values_from = status) %>%
  filter(!is.na(young) & !is.na(old)) %>%
  mutate(
    comparison = case_when(
      young == "silenced" & old == "silenced" ~ "silence_both",
      young == "escape" & old == "escape" ~ "escape_both",
      young == "silenced" & old == "escape" ~ "silence_young_escape_old",
      young == "escape" & old == "silenced" ~ "escape_young_silence_old",
      TRUE ~ "other"
    )
  )

combination_counts <- status_compare %>%
  count(comparison)

gene_counts <- status_compare %>%
  group_by(comparison) %>%
  summarise(unique_genes = n_distinct(external_gene_name))

dim(asenow_filtered)[1]/2
length(unique(asenow_filtered$external_gene_name))
length(unique(asenow_filtered$tissue))
print(combination_counts)
print(gene_counts)
length(unique(subset(status_compare, comparison == 'silence_young_escape_old' | comparison == 'escape_young_silence_old')$external_gene_name))

escape_young_silence_old_df <- status_compare %>%
  filter(comparison == "escape_young_silence_old") %>%
  dplyr::select(external_gene_name, tissue, young, old)
escape_young_silence_old_df
length(unique(escape_young_silence_old_df$external_gene_name))
table(escape_young_silence_old_df$tissue)
table(escape_young_silence_old_df$external_gene_name)

silence_young_escape_old_df <- status_compare %>%
  filter(comparison == "silence_young_escape_old") %>%
  dplyr::select(external_gene_name, tissue, young, old)
silence_young_escape_old_df
table(silence_young_escape_old_df$tissue)
table(silence_young_escape_old_df$external_gene_name)
length(unique(silence_young_escape_old_df$external_gene_name))

combo = rbind(escape_young_silence_old_df, silence_young_escape_old_df)
table(combo$tissue)
table(combo$external_gene_name)
length(unique(combo$external_gene_name))

mod = aov(formula = 1-allelic_expression ~ group + tissue + PAR + external_gene_name, data = asenow_filtered)
summary(mod)
TukeyHSD(mod)

silence_both <- status_compare %>%
  filter(comparison == "silence_both") %>%
  dplyr::select(external_gene_name, tissue, young, old)
silence_both$key = paste(silence_both$external_gene_name, silence_both$tissue)

mod = aov(formula = 1-allelic_expression ~ group + tissue + external_gene_name, 
          data = subset(asenow_filtered, key %in% silence_both$key & PAR != 'PAR'))
summary(mod)
TukeyHSD(mod)

#### end ####

#### plot along X chr ####

ase_coords <- ase %>%
  dplyr::select(external_gene_name, chromosome_name, start_position, end_position) %>%
  filter(chromosome_name == "X") %>%
  distinct() %>%
  group_by(external_gene_name) %>%
  dplyr::slice(1) %>%
  ungroup() 

asenow = subset(ase, tissue %in% keep)
asenow$group = ifelse(asenow$individual == 'UPIC', 'young', 'old')
asenow$key = paste(asenow$external_gene_name, asenow$tissue)
asenow$status = ifelse(asenow$allelic_expression<0.4, 'escape', 'silenced')
asenow_escape = unique(asenow[,c('key','tissue','allelic_expression','external_gene_name','status','group','PAR')])
asenow_filtered <- asenow_escape
status_switch_genes <- asenow_filtered %>%
  group_by(tissue, external_gene_name) %>%
  summarise(n_status = n_distinct(status), .groups = "drop") %>%
  filter(n_status > 1)
asenow_filtered <- asenow_filtered %>%
  mutate(status_change = paste(tissue, external_gene_name) %in%
           paste(status_switch_genes$tissue, status_switch_genes$external_gene_name))
asenow_filtered$group = factor(asenow_filtered$group, levels = c('young','old'))

# Get unique status per group

status_compare_full <- status_collapsed %>%
  pivot_wider(names_from = group, values_from = status)

# Add comparison categories, even when only young or old is present
status_compare <- status_compare_full %>%
  mutate(
    comparison = case_when(
      !is.na(young) & !is.na(old) & young == "silenced" & old == "silenced" ~ "silence_both",
      !is.na(young) & !is.na(old) & young == "escape" & old == "escape" ~ "escape_both",
      !is.na(young) & !is.na(old) & young == "silenced" & old == "escape" ~ "silence_young_escape_old",
      !is.na(young) & !is.na(old) & young == "escape" & old == "silenced" ~ "escape_young_silence_old",
      !is.na(young) & is.na(old) & young == "silenced" ~ "silenced_young_only",
      !is.na(young) & is.na(old) & young == "escape" ~ "escape_young_only",
      is.na(young) & !is.na(old) & old == "silenced" ~ "silenced_old_only",
      is.na(young) & !is.na(old) & old == "escape" ~ "escape_old_only",
      TRUE ~ "other"
    )
  )

table(status_compare$comparison)

# 2. Merge with status_compare to annotate with comparison and tissue
plot_df <- left_join(status_compare, ase_coords, by = "external_gene_name")

plot_df <- plot_df %>%
  filter(tissue %in% unique(ase$tissue))

# 3. Format for GRanges
plot_df <- plot_df %>%
  mutate(chromosome_name = paste0("chr", chromosome_name),
         label = paste0(external_gene_name, "_", tissue))

# Rename comparisons (update your data first)
plot_df <- plot_df %>%
  mutate(comparison = case_when(
    comparison == "silence_both" ~ "Silence (Both)",
    comparison == "escape_both" ~ "Escape (Both)",
    comparison == "escape_young_silence_old" ~ "Escape (Young), Silence (Old)",
    comparison == "silence_young_escape_old" ~ "Silence (Young), Escape (Old)",
    TRUE ~ comparison
  ))
table(plot_df$comparison)
plot_df$comparison = factor(plot_df$comparison, levels = c('Silence (Both)','silenced_old_only','silenced_young_only',
                                                           'Escape (Both)','escape_old_only','escape_young_only',
                                                           'Escape (Young), Silence (Old)','Silence (Young), Escape (Old)'))
levels(plot_df$comparison)

fixed_width <- 100000
plot_df <- plot_df %>%
  mutate(mid = round((start_position + end_position) / 2),
         start_position = mid - fixed_width / 2,
         end_position = mid + fixed_width / 2)

# Define tissue color palette
unique(plot_df$tissue)
tissue_colors <- c(
  "AdiposeSubcutaneous" = "#FFC857",          # warm yellow
  "AdiposeVisceralOmentum" = "#FF8C42",       # orange
  "AdrenalGland" = "#C8553D",                 # brick
  "ArteryAorta" = "#8B1812",                  # (existing)
  "ArteryCoronary" = "#C97A41",               # (existing)
  "ArteryTibial" = "#D1495B",                 # dusty red
  "Breast" = "#ED6A5A",                       # coral
  "EsophagusMuscularis" = "#DCB465",          # (existing)
  "HeartAtrialAppendage" = "#7E909A",         # cool gray
  "HeartLeftVentricle" = "#525E75",           # blue-gray
  "Liver" = "#4C6A92",                        # (existing)
  "Lung" = "#A77A76",                         # (existing)
  "MuscleSkeletal" = "#7FB685",               # muted green
  "Ovary" = "#BC6C25",                        # burnt sienna
  "Pancreas" = "#244C51",                     # (existing)
  "SkinNotSunExposedSuprapubic" = "#9C89B8",  # lavender gray
  "SkinSunExposedLowerleg" = "#6A0572",       # deep purple
  "Spleen" = "#386641",                       # deep forest green
  "Stomach" = "#555463",                      # (existing)
  "Thyroid" = "#8A5C74",                      # (existing)
  "KidneyCortex" = "#588157",                # moss green
  "BrainAnteriorcingulatecortexBA24" = "#023047",  # dark teal
  "BrainCaudatebasalganglia" = "#126782",          # steel blue
  "BrainCerebellarHemisphere" = "#8199A3",         # slate
  "BrainSpinalcordcervicalc1" = "#4B3B47",         # deep plum
  "ColonTransverse" = "#AB4967",              # dusty rose
  "EsophagusMucosa" = "#F28482",              # salmon
  "Uteris" = "#A26769",                       # warm mauve
  "Vagina" = "#B977A6",                       # (existing)
  "WholeBlood" = "#A33B3B"                    # (existing)
)



# Organize plotting data
comparisons <- levels(plot_df$comparison)
plotting_list <- list()
track_map <- list()

track_index <- 1
for (comp in comparisons) {
  tissues <- unique(plot_df$tissue[plot_df$comparison == comp])
  for (tissue in tissues) {
    subset_df <- plot_df %>% filter(comparison == comp & tissue == !!tissue)
    plotting_list[[paste(comp, tissue, sep = "_")]] <- subset_df
    track_map[[paste(comp, tissue, sep = "_")]] <- list(
      comparison = comp,
      tissue = tissue,
      track = track_index
    )
    track_index <- track_index + 1
  }
}

# Gene coordinates for labels (optional)
label_coords <- plot_df %>%
  distinct(chromosome_name, start_position, end_position, external_gene_name, comparison) %>%
  mutate(chromosome_name = as.character(chromosome_name)) %>%
  as.data.frame()
genes_coord_all <- regioneR::toGRanges(label_coords)
View(as.data.frame(genes_coord_all))

# tads
tads <- data.table::fread('/data/A-172_GSE147123_tad.bed', header = FALSE)
tads = subset(tads, V1 == 'chrX')
tads <- tads %>% arrange(V2)
tads_regions <- tads %>%
  mutate(start = V3,          # use end of current boundary as start of TAD
         end = lead(V2)) %>%  # use start of next boundary as end of TAD
  filter(!is.na(end)) %>%     # drop the last row (no following boundary)
  filter(end > start) %>%     # filter out invalid ranges
  transmute(chr = V1, start, end)
tads_gr <- regioneR::toGRanges(tads_regions)

# Plot
kp <- plotKaryotype(
  genome = "hg38", 
  chromosomes = "chrX", 
  plot.type = 3, 
  zoom = "chrX:1-156040895"
)
kpRect(kp,
       data = tads_gr,
       y0 = 0, y1 = 1,
       col = "#D3D3D3",
       border = NA,
       r0 = 0.05, r1 = 0.15)
kpAddBaseNumbers(kp)

kpPlotMarkers(kp, 
              data = subset(genes_coord_all, comparison %in% 
                              c('Escape (Young), Silence (Old)','Silence (Young), Escape (Old)')), 
              label.dist = 0.000001,
              adjust.label.position	= T,
              marker.parts = c(0.3,0.03,0.1),
              #ymax=5,
              r1 = 0.5,
              cex = 0.5,
              labels = subset(genes_coord_all, comparison %in% 
                                c('Escape (Young), Silence (Old)','Silence (Young), Escape (Old)'))$external_gene_name,
              text.orientation = "vertical")
View(as.data.frame(subset(genes_coord_all, comparison %in% 
                            c('Escape (Young), Silence (Old)','Silence (Young), Escape (Old)'))))

kpAbline(kp, v = 2800000, col = "black", lty = 2, lwd = 1)

# Plot each tissue-comparison track

for (name in names(plotting_list)) {
  df <- plotting_list[[name]]
  track <- track_map[[name]]$track
  tissue <- track_map[[name]]$tissue
  kpRect(kp, chr = df$chromosome_name, x0 = df$start_position, x1 = df$end_position, y0 = 0, y1 = 1,
         r0 = autotrack(track, track_index), data.panel = 2,
         col = alpha(tissue_colors[[tissue]], 0.8), border = NA)
}

# Draw comparison boxes and labels
comparison_bounds <- aggregate(unlist(lapply(track_map, function(x) x$track)), 
                               by = list(comparison = sapply(track_map, function(x) x$comparison)), 
                               FUN = range)

for (i in 1:nrow(comparison_bounds)) {
  comp <- comparison_bounds$comparison[i]
  bounds <- comparison_bounds$x[i, ]
  r0r1 <- autotrack(bounds[1], track_index, margin = 0.002)
  r0r1_end <- autotrack(bounds[2], track_index, margin = 0.002)
  
  kpRect(kp, chr = "chrX", x0 = 0, x1 = 156040895, y0 = 0, y1 = 1,
         r0 = r0r1$r0, r1 = r0r1_end$r1, data.panel = 2,
         border = "black", col = NA, lwd = 1.2)
  
  kpAddLabels(kp, labels = comp, r0 = r0r1$r0, r1 = r0r1_end$r1,
              data.panel = 2, pos = 2, srt = 0, cex = 0.75)
}

#### end ####

#### positional enrichment ####

agedep = as.data.frame(subset(genes_coord_all, comparison %in% c('Escape (Young), Silence (Old)','Silence (Young), Escape (Old)')))
agedep <- unique(agedep[, c("external_gene_name", "start", "end")])
agedep$midpoint <- (agedep$start + agedep$end) / 2
agedep$external_gene_name
all = unique(ase$external_gene_name)
all = all[which(all %!in% agedep$external_gene_name)]
all
silenced = as.data.frame(subset(genes_coord_all, comparison %in% c('Silence (Both)')))
unique(silenced$external_gene_name)
dim(agedep)

bin_size <- 10e6
chrX_length <- 156040895  
bins <- data.frame(
  start = seq(1, chrX_length, by = bin_size),
  end = pmin(seq(1, chrX_length, by = bin_size) + bin_size - 1, chrX_length)
)
bins$bin_id <- paste0("bin_", seq_len(nrow(bins)))
assign_bin <- function(pos) {
  bin_index <- which(bins$start <= pos & bins$end >= pos)
  if (length(bin_index) == 1) return(bins$bin_id[bin_index])
  return(NA)
}
agedep$bin <- sapply(agedep$midpoint, assign_bin)
bin_counts <- agedep %>%
  group_by(bin) %>%
  summarise(k = n(), .groups = "drop")
n_total <- nrow(agedep)
n_bins <- nrow(bins)
p_null <- 1 / n_bins
bin_counts$expected <- n_total * p_null
bin_counts$binomial_p <- mapply(function(k) {
  binom.test(k, n_total, p_null, alternative = "greater")$p.value
}, bin_counts$k)
bin_counts <- left_join(bin_counts, bins, by = c("bin" = "bin_id"))

ggplot(bin_counts, aes(x = start / 1e6, y = -log10(binomial_p))) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "Position on chrX (Mb)", y = "-log10(p-value)", title = "Enrichment of Age-Dependent Genes") +
  theme_minimal()

#### end ####


