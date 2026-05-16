########################################
##
## Di Nezio et al.
## Microbial diversity of groundwater related environments in Lanzarote
## Ecological analyses
## Script written by Alejandro Martinez and Francesco Di Nezio
##
########################################

library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)

###### PART 1 - ARRANGING DATA FILE --------------------------------------------
project_dir <- "C:/Users/ASUS/Desktop/lastversionscripts"

data_dir <- file.path(project_dir, "rawdata2")
fun_dir  <- file.path(data_dir, "functions")

### cleaning community matrix

# comm <- read.csv2("sequence_tab.csv")
# taxonomy <- read.csv2("taxonomy.csv")
# names <- readxl::read_xlsx("names.xlsx")
# names$Sample <- gsub("sum_","",names$Sample, fixed = T)
# 
# comm.names <- as.data.frame(t(comm))
# comm.names$taxonomy <- rownames(comm.names)
# 
# colnames(comm.names) <- comm.names[1,]
# comm.names <- comm.names[-1,]
# 
# comm.names <- merge(comm.names, taxonomy[c("X","ID")], by = "X")
# row.names(comm.names) <- comm.names$ID
# comm.names <- comm.names[,c(-1,-32)]
# 
# comm.names[] <- lapply(comm.names, function(x) {
#   if (is.character(x)) gsub(" ", "", x) else x
# })
# 
# comm.names[] <- lapply(comm.names, function(x) as.numeric(as.character(x)))
# 
# names$Sample %in% colnames(comm.names)
# 
# name_vector <- setNames(names$Name, names$Sample)
# colnames(comm.names) <- ifelse(colnames(comm.names) %in% names(name_vector),
#                          name_vector[colnames(comm.names)],
#                          colnames(comm.names))
# 
# write.csv2(comm.names, "sequence_tab_asvs.csv")
# 
# 
# 
# stations <- read.csv2("Metadat.csv", sep = ";", dec = ".")

  
###### 1.1 - Arrange taxonomy-free community matrix ----------------------------

comm <- read.csv2(file.path(data_dir, "sequence_tab_asvs.csv"), header = 1)
taxonomy <- read.csv2(file.path(data_dir, "taxonomy.csv"))
stations <- read.csv2(file.path(data_dir, "Metadat.csv"))
rownames(stations) <- stations$Sample
traits <- read.csv2(file.path(data_dir, "pathogenicity_review_v2.csv"))

#### 1.1.1 - Clean community matrix --------------------------------------------

colnames(comm) <- comm[1,]
comm <- comm[-1,]

row.names(comm) <- comm[,1]
comm <- comm[,-1]

comm[] <- lapply(comm, function(x) {
  if (is.character(x)) gsub(" ", "", x) else x
})

comm[] <- lapply(comm, function(x) as.numeric(as.character(x)))
str(comm)

#### 1.1.2 - Filtering out ASVs corresponding to wrong groups ------------------

taxonomy[is.na(taxonomy)] <- "Unclassified"

# --- BEFORE filtering ---

# total reads per ASV
asv_reads <- colSums(comm, na.rm = TRUE)

# subset per category
archaea_ids       <- taxonomy[taxonomy$Kingdom == "Archaea", "ID"]
mito_ids          <- taxonomy[taxonomy$Family == "Mitochondria", "ID"]
chloro_ids        <- taxonomy[taxonomy$Order == "Chloroplast", "ID"]
unclassified_ids  <- taxonomy[taxonomy$Phylum == "Unclassified", "ID"]

# calculate reads
reads_archaea      <- sum(asv_reads[names(asv_reads) %in% archaea_ids], na.rm = TRUE)
reads_mito         <- sum(asv_reads[names(asv_reads) %in% mito_ids], na.rm = TRUE)
reads_chloro       <- sum(asv_reads[names(asv_reads) %in% chloro_ids], na.rm = TRUE)
reads_unclassified <- sum(asv_reads[names(asv_reads) %in% unclassified_ids], na.rm = TRUE)

total_reads <- sum(asv_reads, na.rm = TRUE)

# summary table
reads_removed <- data.frame(
  Category = c("Archaea", "Mitochondria", "Chloroplast", "Unclassified"),
  Reads = c(reads_archaea, reads_mito, reads_chloro, reads_unclassified)
)

reads_removed$Percentage <- 100 * reads_removed$Reads / total_reads

reads_removed

archaea <- taxonomy[ which (taxonomy$Kingdom == "Archaea"),"ID"]
mitochondria <- taxonomy[ which (taxonomy$Family == "Mitochondria"),"ID"]
chloroplast <- taxonomy[ which (taxonomy$Order == "Chloroplast"),"ID"]
wrong.seq <- taxonomy[ which (taxonomy$Phylum == "Unclassified"),"ID"]

drop.seqs <- c(archaea,mitochondria,chloroplast,wrong.seq)

"%ni%" <- Negate("%in%")
comm <- comm[ , which (colnames(comm) %ni% drop.seqs)]

rm(archaea,mitochondria,chloroplast,wrong.seq, drop.seqs)

#### 1.1.3 - Deal with rare community matrix -----------------------------------

### In the taxonomy-free matrix, we filter out sequences present in less than 
### two samples and totaling less than three reads.

comm <- comm[,
  colSums(comm > 0) >= 2 & colSums(comm) >= 3
]

# comm.rare <- comm[,
#              colSums(comm > 0) < 2 & colSums(comm) < 3
# ]

###### 1.2 - Complete the stations matrix --------------------------------------

#### 1.2.1 - Calculate richness per station ------------------------------------

comm.rich <- comm
comm.rich[comm.rich > 0] <- 1

richness <- as.data.frame(rowSums(comm.rich))
colnames(richness) <- "richness"
richness$Sample <- rownames(richness)

stations <- merge(stations, richness, by = "Sample")
rownames(stations) <- stations$Sample

rm(comm.rich,richness)

#### 1.2.2 -  Calculate shannon diversity per station ----------

d.shannonH <- as.data.frame(vegan::diversity(comm, index = "shannon"))
colnames(d.shannonH) <- "shannon"
d.shannonH$Sample <- rownames(d.shannonH)

stations <- merge(stations, d.shannonH, by = "Sample")
rownames(stations) <- stations$Sample

rm(d.shannonH)

#### 1.2.3 - Add number of reads ----------

reads <- read.csv(file.path(data_dir,"number_of_reads.csv"))
names <- as.data.frame(readxl::read_xlsx(file.path(data_dir,"names.xlsx")))

names$Sample <- gsub("sum_","", names$Sample, fixed = T)
reads <- merge(reads, names, by = "Sample")

stations <- merge(stations, reads[c("Name","Reads_PostChimera")], by.x = "Sample", by.y = "Name")
colnames(stations)[colnames(stations) == "Reads_PostChimera"] <- "reads"

rm(reads, names)

#### 1.2.4 -  Count buildings and measure m of road ----------

source(file.path(fun_dir, "GIS_covariates.R"))

library(sf)
library(osmdata)
library(units)

stations$buildings <- count_buildings_osm(latitudes = as.numeric(stations$Latitude),
                                         longitudes = as.numeric(stations$Longitude),
                                         radius = 1000)

# stations$buildings <- c(69,69,15,475,0,0,68,4,
#                        141,136,136,3,0,29,0,2,
#                        0,2,2,2,0,2,0,0,344,0,2,
#                        2,2,2)

stations$roads <- calculate_road_lengths(latitudes = as.numeric(stations$Latitude),
                                          longitudes = as.numeric(stations$Longitude),
                                          radius = 1000)

# stations$roads <-  c(8579.325,8579.325,5925.543,75528.494,5717.383,
#                     5182.547,7890.656,6479.438,11556.677,10361.297,
#                     10361.297,10195.618,10256.365,16542.026,6820.342,
#                     6413.168,0,6857.950,0,9609.712,9609.712,9609.712,
#                    4451.792,4400.386,24303.961,0,8980.558,7864.168,
#                     6461.287,7449.642)

##### 1.3 - Complete taxonomy matrix ----------------------------------------

#### 1.3.1 -  Add taxon name -------------

taxonomy[is.na(taxonomy)] <- "Unclassified"

taxonomy <- taxonomy %>%
  mutate(
    Taxon = case_when(
      Genus != "Unclassified" & Species != "Unclassified" ~ paste(Genus, Species),
      Genus != "Unclassified" & Species == "Unclassified" ~ Genus,
      Genus == "Unclassified" & Species == "Unclassified" ~ coalesce(
        ifelse(Family != "Unclassified", Family, NA_character_),
        ifelse(Order != "Unclassified", Order, NA_character_),
        ifelse(Class != "Unclassified", Class, NA_character_),
        ifelse(Phylum != "Unclassified", Phylum, NA_character_)
      )
    )
  )

#### 1.3.2 -  Save clean taxonomy -------------

taxonomy.checked <- taxonomy[ which (taxonomy$ID %in% colnames(comm)),]

comm.names <- as.data.frame(t(comm))
comm.names$ID <- row.names(comm.names)
taxonomy.checked <- merge(taxonomy.checked,comm.names, by ="ID", all = F)

colnames(taxonomy.checked)[2] <- "Sequence"

write.table(taxonomy.checked, "taxonomy.checked.csv", dec = ".", sep = ";", row.names = FALSE)

#### 1.3.3 -  Save ASVs sequences as fasta -------------

source(file.path(fun_dir, "dataframe2fas.R"))

dataframe2fas(taxonomy.checked[c("ID","Sequence")], "taxonomy.check.fasta")

rm(comm.names)

###### 1.3.4 - Delete samples with little reads --------------------------------

drop <- c("CSP","TDA_MA")

stations <- stations[ which (stations$Sample %ni% drop),]
comm <- comm[ which (rownames(comm) %ni% drop),]

stations <- stations[match(rownames(comm), stations$Sample), ]
rownames(stations) <- stations$Sample
stopifnot(all(rownames(comm) == stations$Sample))

rm(drop)

library(indicspecies)
library(vegan)
library(tibble)

########## 1.3.5 - Generate a taxonomically-ranked community matrix ----------

source(file.path(fun_dir, "summarytaxonomy.R"))

comm.tax <- as.data.frame(t(comm))
comm.tax$asvs <- row.names(comm.tax)

comm.tax <- merge(taxonomy, comm.tax, by.x = "ID", by.y = "asvs", all.y = TRUE)

comm.tax <- comm.tax %>%
  mutate(
    Taxon = case_when(
      Genus != "Unclassified" & Species != "Unclassified" ~ paste(Genus, Species),
      Genus != "Unclassified" & Species == "Unclassified" ~ Genus,
      Genus == "Unclassified" & Species == "Unclassified" ~ coalesce(
        ifelse(Family != "Unclassified", Family, NA_character_),
        ifelse(Order != "Unclassified", Order, NA_character_),
        ifelse(Class != "Unclassified", Class, NA_character_),
        ifelse(Phylum != "Unclassified", Phylum, NA_character_)
      )
    )
  )

comm.tax.list <- summarize_by_taxonomic_rank(comm.tax)

###### 1.3.6 - Indicator taxa analysis -------------

# Genus matrix: samples in rows, genera in columns
genus_mat <- comm.tax.list$Genus

# keep only samples present in stations
genus_mat <- genus_mat[rownames(genus_mat) %in% stations$Sample, , drop = FALSE]
stations_ind <- stations[match(rownames(genus_mat), stations$Sample), ]

# optional: exclude habitats with too few replicates (e.g. pond if only one sample)
tab_groups <- table(stations_ind$Type)
keep_groups <- names(tab_groups)[tab_groups >= 2]

keep <- stations_ind$Type %in% keep_groups
genus_mat <- genus_mat[keep, , drop = FALSE]
stations_ind <- stations_ind[keep, , drop = FALSE]

group <- droplevels(factor(stations_ind$Type))

# remove very rare genera
genus_mat <- genus_mat[, colSums(genus_mat > 0) >= 2, drop = FALSE]

# run IndVal
set.seed(123)
indval_genus <- multipatt(x = genus_mat, cluster = group, func = "IndVal.g",
  duleg = TRUE, control = how(nperm = 9999))

# extract results
indval_sign <- as.data.frame(indval_genus$sign) %>%
  tibble::rownames_to_column("Genus")

# identify the habitat associated with each genus
indval_sign <- indval_sign %>%
  mutate(
    Habitat = case_when(
      s.cave == 1 ~ "cave",
      s.pool == 1 ~ "pool",
      s.salt == 1 ~ "salt",
      s.sea  == 1 ~ "sea",
      s.well == 1 ~ "well",
      TRUE ~ NA_character_
    )
  )

# build Genus -> Family lookup from taxonomy
genus_family_lookup <- taxonomy %>%
  select(Genus, Family) %>%
  filter(!is.na(Genus), Genus != "Unclassified") %>%
  distinct() %>%
  group_by(Genus) %>%
  summarise(
    Family = paste(sort(unique(Family[!is.na(Family) & Family != "Unclassified"])), collapse = "; "),
    .groups = "drop"
  )

# add Family and prevalence / mean abundance summaries
genus_pa <- as.data.frame(genus_mat > 0) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Genus", values_to = "Present")

genus_abund <- as.data.frame(genus_mat) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Genus", values_to = "Abundance")

sample_habitat <- stations_ind %>%
  select(Sample, Type) %>%
  rename(Habitat = Type)

prev_tbl <- genus_pa %>%
  left_join(sample_habitat, by = "Sample") %>%
  group_by(Habitat, Genus) %>%
  summarise(prevalence = mean(Present, na.rm = TRUE),
    .groups = "drop")

mean_abund_tbl <- genus_abund %>%
  left_join(sample_habitat, by = "Sample") %>%
  group_by(Habitat, Genus) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop")

# final table
indval_genus_tbl <- indval_sign %>%
  rename(IndVal = stat, p_value = p.value) %>%
  left_join(genus_family_lookup, by = "Genus") %>%
  left_join(prev_tbl, by = c("Habitat", "Genus")) %>%
  left_join(mean_abund_tbl, by = c("Habitat", "Genus")) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"),
    Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family)) %>%
  arrange(p_adj, desc(IndVal))

# significant indicators
indval_genus_sig <- indval_genus_tbl %>%
  filter(p_adj < 0.1) %>%
  arrange(Habitat, p_adj, desc(IndVal))

writexl::write_xlsx(indval_genus_tbl, "indval_genus_all.xlsx")
writexl::write_xlsx(indval_genus_sig, "indval_genus_significant.xlsx")

##### 1.4 - Definition of functional groups ------------------------------------

pot_pathogens <- traits[ which (traits$Pathogenicity == "2" |
                              traits$Pathogenicity == "3"),]
nrow(pot_pathogens)

pot_pathogens.asv <- taxonomy[ which(taxonomy$Taxon %in% pot_pathogens$Taxon ), ]
comm.pathogens <- comm[ , which ( colnames(comm) %in%  pot_pathogens.asv$ID)]
rich.pathogens <- as.data.frame(rowSums(comm.pathogens > 0))
rich.pathogens$Sample <- rownames(comm.pathogens)
colnames(rich.pathogens)[1] <- "richness.pathogens"

shannonH.patho <- as.data.frame(vegan::diversity(comm.pathogens, index = "shannon"))
shannonH.patho$Sample <- rownames(shannonH.patho)
colnames(shannonH.patho)[1] <- "shannon.pathogens"

abund.pathogens <- as.data.frame(rowSums(comm.pathogens))
abund.pathogens$Sample <- rownames(abund.pathogens)
colnames(abund.pathogens)[1] <- "abundance.pathogens"

#
waste <- traits[ which (traits$Origin == "W"),]

waste.asv <- taxonomy[ which(taxonomy$Taxon %in% waste$Taxon ), ]
comm.waste <- comm[ , which ( colnames(comm) %in%  waste.asv$ID)]
rich.waste <- as.data.frame(rowSums(comm.waste > 0))
rich.waste$Sample <- rownames(comm.waste)
colnames(rich.waste)[1] <- "richness.waste"

shannonH.waste <- as.data.frame(vegan::diversity(comm.waste, index = "shannon"))
shannonH.waste$Sample <- rownames(shannonH.waste)
colnames(shannonH.waste)[1] <- "shannon.waste"

abund.waste <- as.data.frame(rowSums(comm.waste))
abund.waste$Sample <- rownames(abund.waste)
colnames(abund.waste)[1] <- "abund.waste"

#
environmental <- traits[ which (traits$Origin == "E"),]
nrow(environmental)

environmental.asv <- taxonomy[ which(taxonomy$Taxon %in% environmental$Taxon ), ]
comm.environmental <- comm[ , which ( colnames(comm) %in%  environmental.asv$ID)]
rich.environment <- as.data.frame(rowSums(comm.environmental > 0))
rich.environment$Sample <- rownames(comm.environmental)
colnames(rich.environment)[1] <- "richness.environment"

shannonH.environment <- as.data.frame(vegan::diversity(comm.environmental, index = "shannon"))
shannonH.environment$Sample <- rownames(shannonH.environment)
colnames(shannonH.environment)[1] <- "shannon.environment"

abund.environment <- as.data.frame(rowSums(comm.environmental))
abund.environment$Sample <- rownames(abund.environment)
colnames(abund.environment)[1] <- "abund.environment"

stations <- merge(stations, abund.pathogens, by = "Sample")
stations <- merge(stations, rich.pathogens, by = "Sample")
stations <- merge(stations, shannonH.patho, by = "Sample")
stations <- merge(stations, abund.waste, by = "Sample")
stations <- merge(stations, rich.waste, by = "Sample")
stations <- merge(stations, shannonH.waste, by = "Sample")
stations <- merge(stations, abund.environment, by = "Sample")
stations <- merge(stations, rich.environment, by = "Sample")
stations <- merge(stations, shannonH.environment, by = "Sample")

rm(pot_pathogens, pot_pathogens.asv, rich.pathogens, waste, waste.asv, rich.waste, 
   rich.environment, environmental, environmental.asv, shannonH.environment, 
   shannonH.patho, shannonH.waste, abund.environment, abund.pathogens, abund.waste)

###### GOAL 1: DESCRIPTIVE PART ------------------------------------------------

########## G 1.1. Functional groups overall ------------------------------------

# --- palette for overlap-aware categories 

cols_combo <- c(
  "Environmental only"                      = "#2c7bb6",
  "Environmental + Potentially Pathogenic"  = "#5e3c99",
  "Human-related only"                      = "#fdae61",
  "Human-related + Potentially Pathogenic"  = "#d95f02",
  "Potentially Pathogenic (origin unknown)" = "#d7191c",
  "Unassigned"                              = "grey80")

# --- map traits to ASVs; keep Pathogenic orthogonal to origin

traits_by_taxon <- traits %>%
  group_by(Taxon) %>%
  summarise(is_patho = any(Pathogenicity %in% c("2","3"), na.rm = TRUE),
    env_flag = any(Origin == "E", na.rm = TRUE),
    ant_flag = any(Origin == "W", na.rm = TRUE),
    .groups  = "drop") %>%
  mutate(origin = dplyr::case_when(
      env_flag & !ant_flag ~ "Environmental",
      ant_flag & !env_flag ~ "Anthropogenic",
      TRUE                 ~ NA_character_)) %>%
  select(Taxon, is_patho, origin, env_flag, ant_flag)

ambig <- traits_by_taxon %>% filter(env_flag & ant_flag)
if (nrow(ambig) > 0) message("Ambiguous origin for ", nrow(ambig), " taxa; origin set to NA.")


cat_map <- taxonomy.checked %>%
  select(ID, Taxon) %>%
  left_join(traits_by_taxon, by = "Taxon", relationship = "many-to-one")

dup_ids <- cat_map %>% count(ID) %>% filter(n > 1)
if (nrow(dup_ids) > 0) {
  stop("Join duplicated ", nrow(dup_ids), " ASV IDs. Investigate taxonomy or trait duplicates.")
}

cat_map <- cat_map %>%
  mutate(
    combo = case_when(
      origin == "Environmental"  &  is_patho ~ "Environmental + Potentially Pathogenic",
      origin == "Environmental"  & !is_patho ~ "Environmental only",
      origin == "Anthropogenic"  &  is_patho ~ "Human-related + Potentially Pathogenic",
      origin == "Anthropogenic"  & !is_patho ~ "Human-related only",
      is_patho & is.na(origin)                 ~ "Potentially Pathogenic (origin unknown)",
      TRUE                                     ~ "Unassigned"))

comm_long <- as.data.frame(comm) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ID", values_to = "reads") %>%
  left_join(cat_map, by = "ID") %>%
  mutate(
    origin2 = forcats::fct_na_value_to_level(origin),
    combo   = forcats::fct_relevel(combo, names(cols_combo)))

########### G 1.2. Functional groups per sample --------------------------------

hab_order <- c("cave","pool","salt","sea","well","pond")

sample_order <- c("CUL1","CUL2","JDA_pool1","JDA_pool2","JDA_pool3",
                  "TDA_entrada","TDA_sima","TDA_LE",
                  "CHL","CHP","CHS","CHZ","MBB","MBS",
                  "CBC1","CBC2","CBL","CHG","JDA_beach1","JDA_beach2",
                  "CSS1","CSS2",
                  "FWH","FWM","GAF","SHW","TAB",
                  "FUC")

## --- Abundance relative to total community reads

reads_sample_combo <- comm_long %>%
  group_by(Sample) %>%
  mutate(total_reads = sum(reads, na.rm = TRUE)) %>%
  group_by(Sample, combo, total_reads) %>%
  summarise(n_reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  mutate(p = n_reads / total_reads) %>%
  left_join(stations[, c("Sample","Type")], by = "Sample") %>%
  mutate(
    Type  = factor(Type, levels = hab_order),
    combo = forcats::fct_relevel(combo, names(cols_combo))) %>%
  group_by(Type) %>%
  arrange(Sample, .by_group = TRUE) %>%
  mutate(
    Type = forcats::fct_recode(Type,
                               "cave" = "cave",
                               "anc.p" = "pool",
                               "sltw" = "salt",
                               "sea" = "sea",
                               "well" = "well",
                               "spr" = "pond"),
    Sample_ord = factor(Sample, levels = sample_order)) %>%
  ungroup()

p_hab_combo <- ggplot(reads_sample_combo, aes(x = Sample_ord, y = p, fill = combo)) +
  geom_col(color = "white", width = 0.95) +
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  scale_y_continuous( limits= c(0,1), labels = scales::label_percent(accuracy = 1)) +  
  scale_fill_manual(values = cols_combo, name = "") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 11),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(6, "mm")) +  
  labs(x = NULL, y = NULL, title = "a") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 11, hjust = 1),
        axis.text.y = element_text(angle = 90, vjust = 0.5, size = 11, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        strip.text.x = element_text(size = 10, angle = 45))

# write.csv2(reads_sample_combo, "functional_abundance_persample.csv")

# ---- Richness 
reads_sample_combo_presence <- comm_long %>%
  filter(reads > 0) %>%
  distinct(Sample, ID, combo) %>%   # garantisce ASV uniche
  group_by(Sample, combo) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  left_join(stations[, c("Sample","Type")], by = "Sample") %>%
  mutate(
    Type = factor(Type, levels = hab_order),
    Type = forcats::fct_recode(Type,
                               "cave" = "cave",
                               "anc.p" = "pool",
                               "sltw" = "salt",
                               "sea" = "sea",
                               "well" = "well",
                               "spr" = "pond"),
    combo = forcats::fct_relevel(combo, names(cols_combo)),
    Sample_ord = factor(Sample, levels = sample_order))


p_hab_rich <- ggplot(reads_sample_combo_presence, aes(x = Sample_ord, y = p, fill = combo)) +
  geom_col(color = "white", width = 0.95) +
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  scale_y_continuous( limits= c(0,1), labels = scales::label_percent(accuracy = 1)) +  
  scale_fill_manual(values = cols_combo, name = "") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 11),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(6, "mm")) +  
    labs(x = NULL, y = NULL, title = "b") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11, hjust = 1),
        axis.text.y = element_text(angle = 90, vjust = 0.5, size = 11, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        strip.text.x = element_text(size = 10, angle = 45,
                                    margin = margin(t = 2, r = 0, b = 1, l = 0)))

# combine in multipanel figure
library(patchwork)
Figure_3 <- (p_hab_combo + p_hab_rich) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", plot.margin = margin(5.5, 15, 5.5, 15))

ggsave("Figure_3.pdf", plot = Figure_3, device = cairo_pdf, width = 250,
  height = 150, units = "mm")

ggsave("Figure_3.tiff", plot = Figure_3, width = 11, height = 8, units = "in",
  dpi = 600, compression = "lzw")

ggsave("Figure_3.png", plot = Figure_3, width = 11, height = 8, units = "in",
  dpi = 600)

# write.csv2(reads_sample_combo_presence, "functional_richness_persample.csv")

###### G 1.3 Taxonomic Community composition graphs ---------------------------

source(file.path(fun_dir, "taxonomy-plot.R"))

library(pals)  
library(grid)

comm_long.tax <- prep_comm_long(comm, taxonomy.checked)

comm_long.patho <- prep_comm_long(comm.pathogens, taxonomy.checked)
comm_long.anthro <- prep_comm_long(comm.waste, taxonomy.checked)
comm_long.envi <- prep_comm_long(comm.environmental, taxonomy.checked)

groups <- c("cave","cave","cave","cave","cave",
           "cave","cave","cave",
           "anc.p","anc.p","anc.p","anc.p","anc.p","anc.p",
           "sea","sea","sea","sea","sea","sea",
           "sltw","sltw",
           "well","well","well","well","well",
           "spr")

#### G 1.3.1 All bacterial community graphs --------------------------

# Recalculate relative abundance against total community
total_reads_sample <- as.data.frame(comm) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(total_reads = rowSums(dplyr::across(-Sample), na.rm = TRUE)) %>%
  dplyr::select(Sample, total_reads)

comm_long.tax <- comm_long.tax %>%
  dplyr::left_join(total_reads_sample, by = "Sample") %>%
  dplyr::mutate(rel_abund = abund / total_reads) %>%
  dplyr::select(-total_reads)

# Sample -> habitat mapping
groups_named <- groups
names(groups_named) <- sample_order

# Keep these labels exactly as used in your plots
hab_order <- c("cave", "anc.p", "sltw", "sea", "well", "spr")

# Load heatmap function
source(file.path(fun_dir, "taxrank_heatmap.R"))

# Build Family and Genus panels
hm_family <- taxrank_heatmap(
  df = comm_long.tax,
  tax_rank = "Family",
  panel_title = "a",
  n_show = 10,
  show_unclassified = FALSE,
  fill_breaks = c(0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.075, 0.1, 0.25, 0.5),
  text_white_threshold = 0.01
)

hm_genus <- taxrank_heatmap(
  df = comm_long.tax,
  tax_rank = "Genus",
  panel_title = "b",
  n_show = 10,
  show_unclassified = FALSE,
  fill_breaks = c(0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.075, 0.1, 0.25, 0.5),
  text_white_threshold = 0.01
)

# Combine panels
library(patchwork)

Figure_2 <- (hm_family | hm_genus) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "right")

Figure_2


# Save

ggsave("Figure_2.pdf", plot = Figure_2, width = 11, height = 8,
  units = "in", device = cairo_pdf)

ggsave("Figure_2.tiff", plot = Figure_2, width = 11, height = 8,
  units = "in", dpi = 600, compression = "lzw")

ggsave("Figure_2.png", plot = Figure_2, width = 12, height = 8,
  units = "in", dpi = 600)

#### G 1.3.2 Genus -level graphs --------------------------

write.csv2(comm_long,"taxonomy_all_persample.csv")
write.csv2(comm_long.patho,"taxonomy_pathogenes_persample.csv")
write.csv2(comm_long.anthro,"taxonomy_anthropogenic_persample.csv")
write.csv2(comm_long.envi,"taxonomy_environmental_persample.csv")

# Recalculate relative abundance against total community

total_reads_sample <- as.data.frame(comm) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(total_reads = rowSums(across(-Sample), na.rm = TRUE)) %>%
  select(Sample, total_reads)

comm_long.patho <- comm_long.patho %>%
  left_join(total_reads_sample, by = "Sample") %>%
  mutate(rel_abund = abund / total_reads) %>%
  select(-total_reads)

comm_long.anthro <- comm_long.anthro %>%
  left_join(total_reads_sample, by = "Sample") %>%
  mutate(rel_abund = abund / total_reads) %>%
  select(-total_reads)

comm_long.envi <- comm_long.envi %>%
  left_join(total_reads_sample, by = "Sample") %>%
  mutate(rel_abund = abund / total_reads) %>%
  select(-total_reads)

# Sample -> habitat mapping

groups_named <- groups
names(groups_named) <- sample_order

hab_order <- c("cave","anc.p","sltw","sea","well","spr")

# Shared color scale limit for one common legend

global_max <- quantile(
  c(comm_long.patho$rel_abund,
    comm_long.anthro$rel_abund,
    comm_long.envi$rel_abund),
  0.95, na.rm = TRUE)

# Heatmap function
source(file.path(fun_dir, "habitat_heatmap.R"))

# Build panels

hm_patho <- habitat_heatmap(
  comm_long.patho,
  panel_title = "a",
  n_show = 10,
  show_unclassified = FALSE)

hm_anthro <-habitat_heatmap(
  comm_long.anthro,
  panel_title = "b",
  n_show = 10,
  show_unclassified = FALSE)

hm_envi <- habitat_heatmap(
  comm_long.envi,
  panel_title = "c",
  n_show = 10,
  show_unclassified = FALSE)

Figure_S2 <- (hm_patho | hm_anthro | hm_envi) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Save
ggsave("Figure_S2.pdf", plot = Figure_S2, width = 14, height = 9,
  units = "in", device = cairo_pdf)

ggsave("Figure_S2.tiff", plot = Figure_S2, width = 14, height = 9,
  units = "in", dpi = 600, compression = "lzw")

ggsave("Figure_S2.png", plot = Figure_S2, width = 14,
  height = 9, units = "in", dpi = 600)

########## GOAL 2 - INFERENTIAL PART  ------------------------------------------

###### G 2.1 - Total bacterial communities in Lanzarote ------------------------

###### G 2.1.1a - ASV richness global dataset -----------------------------------

psych::pairs.panels(stations[c(4:12,14:17)]) 

mod.rich <- MASS::glm.nb(richness ~ Type + reads + Distance.sea + log(buildings + 1),
                       data = stations)

performance::check_model(mod.rich) 
performance::check_overdispersion(mod.rich) 
car::Anova(mod.rich)
summary(mod.rich)

emm_type <- emmeans::emmeans(mod.rich, ~ Type)
tukey_type <- pairs(emm_type, adjust = "tukey")
summary(tukey_type)

###### G 2.1.1b - Shannon diversity global dataset -----------------------------

mod.shannon <- lm(shannon ~ Type + reads + Distance.sea + log(buildings + 1),
                  data = stations)

performance::check_model(mod.shannon)
car::Anova(mod.shannon)
summary(mod.shannon)

emm_type_shannon <- emmeans::emmeans(mod.shannon, ~ Type)
tukey_type_shannon <- pairs(emm_type_shannon, adjust = "tukey")
summary(tukey_type_shannon)

###### G 2.1.2a - ASV richness visualization ----------------------------------------

# --- Custom palette ---
my_palette <- c(
  cave = "#4B4B4B",   
  pool = "#1B9E77",   
  salt = "#E6AB02",   
  sea  = "#1F78B4",   
  well = "#A6761D",   
  pond = "#66C2A5")

stations <- stations |>
  mutate(Type = fct_relevel(Type, c("cave","pool","salt","sea","well","pond")))

# --- Plot ---

# Count observations per habitat
stations2 <- stations %>%
  group_by(Type) %>%
  mutate(n_type = sum(!is.na(richness))) %>%
  ungroup()

# Separate datasets
stations_box <- stations2 %>% filter(n_type > 4)
stations_dot <- stations2 %>% filter(n_type <= 4)

p_rich.all <- ggplot() +
  # Boxplots only for habitats with > 4 observations
  geom_boxplot(
    data = stations_box,
    aes(x = Type, y = richness, fill = Type),
    width = 0.6,
    color = "black",
    outlier.shape = NA
  ) +
  # Jittered observations on top of boxplots
  geom_jitter(
    data = stations_box,
    aes(x = Type, y = richness, fill = Type),
    width = 0.12,
    size = 2.2,
    shape = 21,
    color = "black",
    alpha = 0.8
  ) +
  # Dot plot for habitats with <= 4 observations
  geom_jitter(
    data = stations_dot,
    aes(x = Type, y = richness, fill = Type),
    width = 0.12,
    size = 2.5,
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_manual(values = my_palette) +
  scale_x_discrete(labels = c(
    cave = "cave",
    pool = "anc.p.",
    salt = "sltw",
    sea = "sea",
    well = "well",
    pond = "spr"
  )) +
  labs(x = "", y = "Richness", title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )

p_rich.all


###### G 2.1.2b - Shannon diversity visualization ------------------------------

(p_shannon.all <- ggplot(stations, aes(x = Type, y = shannon, color = Type)) +
  geom_point(size = 3, position = position_jitter(width = 0, height = 0)) + 
  stat_summary(fun = mean, geom = "crossbar", linetype = "solid",
               width = 0.6, fatten = 2, color = "black") +
  scale_color_manual(values = my_palette) +
  scale_x_discrete(labels = c(
    cave = "cave",
    pool = "anc.p",
    salt = "sltw",
    sea  = "sea",
    well = "well",
    pond = "spr"
  )) +
  labs(x = "", y = "Shannon diversity", title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ))

###### G 2.1.3 - ASV community composition analyses global dataset ----------------------------------------

# check for multivariate dispersion  
dist.bc <- vegdist(comm, method = "bray")

bd <- betadisper(dist.bc, stations$Type)
permutest(bd, permutations = 9999)

beta.tax <- BAT::beta(comm, abund = T)

mean(beta.tax$Btotal)
mean(beta.tax$Brepl / beta.tax$Btotal)
mean(beta.tax$Brich / beta.tax$Btotal)

source(file.path(fun_dir, "beta_plot2.R"))

(pb.all <- beta_density_plot(comm, stations$Type, beta.tax,
                              fill_values = c("Within habitat" = "#1b9e77",
                                              "Across habitats" = "#d95f02"),
                              include_total = TRUE,
                              include_nestedness = FALSE,
                              facet_cols = 2,
                              title = NULL,
                              alpha = 0.6) +
                              theme_classic() +
                              theme(axis.text.x = element_text(size = 12, 
                                                  angle = 45, hjust = 1),
                                      axis.text.y  = element_text(size = 12),
                                      axis.title.y = element_text(size = 14),
                                      axis.title.x = element_text(size = 12),
                                      strip.text.x = element_text(size = 14),
                                      strip.background = element_blank()))

## Total beta diversity
model.Btotal <- beta.tax$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxT <- vegan::adonis2(model.Btotal, stations, permutations = 9999, by = "terms"))

## Beta nestness
model.Bnest <- beta.tax$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi <- vegan::adonis2(model.Bnest, stations, permutations = 9999, by = "terms"))

## Beta turnover
model.Bturn <- beta.tax$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe <- vegan::adonis2(model.Bturn, stations, permutations = 9999, by = "terms"))

rm(model.Btotal,model.Bnest,model.Bturn,perm.taxRi,perm.taxRe,perm.taxT)

###### G 2.1.4 - ASV nMDS visualization ----------------------------------------

source(file.path(fun_dir, "plot_nmds.R"))

library(vegan)

lvl_order <- c("cave","pool","salt","sea","well","pond")
disp_lab <- c(cave="cave",pool="anc.p",salt='sltw',sea="sea",well="well",pond='spr')

stations_ord <- stations |> mutate(Type = fct_relevel(Type, lvl_order))
rownames(stations_ord) <- stations_ord$Sample

nmds1 <- plot_nmds(
  comm, stations_ord, my_palette,
  plot_title = NULL,
  point_size = 3,
  label_size = 0,
  label_vjust = -1
) +
  scale_color_manual(
    values = my_palette,
    labels = disp_lab[names(my_palette)],
    name = "Habitat"
  ) +
  theme(
    axis.text.x = element_text(size = 12, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 12),
    strip.text.x = element_text(size = 14)
  )

######### G 2.1.5 - OVERALL FIGURE ---------------------------------------------
# Layout: row1 = A | B ; row2 = C spans both columns
design <- "
ab
cc
"

(fig_small <-
    ((p_rich.all + guides(fill = "none", color = "none")) +
       nmds1 +
       pb.all) +
    plot_layout(design = design, guides = "collect", heights = c(1, 1.3)) +
    plot_annotation(tag_levels = "a") &
    theme_classic() &
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = ggplot2::margin(b = 4)),
      strip.text = element_text(size = 14),
      strip.background = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 10),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0, 1)))

# PDF
ggsave(filename = "Figure_4.pdf", plot = fig_small, device = cairo_pdf, width = 11,
  height = 8.5, units = "in")

# TIFF
ggsave(filename = "Figure_4.tiff", plot = fig_small, width = 11, height = 8.5,
  units = "in", dpi = 600, compression = "lzw")

# PNG
ggsave(filename = "Figure_4.png", plot = fig_small, width = 11, height = 8.5,
  units = "in", dpi = 600)

###### G 2.2 - Functional bacterial groups in Lanzarote ----------------------------------------

library(glmmTMB)

###### G 2.2.1a - Pathogens abundance and richness  -------------------------------------------------------

stations$Type <- relevel(factor(stations$Type), ref = "sea")

#
mod.patho.abund.prob <- glmmTMB(cbind(abundance.pathogens, reads - abundance.pathogens) ~ 
                              Type + Distance.sea + log(buildings + 1),
                      family=betabinomial(link = "logit"), data=stations)

performance::check_model(mod.patho.abund.prob) 
performance::check_overdispersion(mod.patho.abund.prob) 
car::Anova(mod.patho.abund.prob)
summary(mod.patho.abund.prob)

emm_type.patho.abund.prob <- emmeans::emmeans(mod.patho.abund.prob, ~ Type)
tukey_type.patho.abund.prob <- pairs(emm_type.patho.abund.prob, adjust = "tukey")
summary(tukey_type.patho.abund.prob)

rm(emm_type.patho.abund.prob, tukey_type.patho.abund.prob)

#
mod.patho.prob <- glm(cbind(richness.pathogens, richness - richness.pathogens) ~ 
                        Type + Distance.sea + reads + log(buildings+1),
                      family=binomial(logit), data=stations)

performance::check_model(mod.patho.prob) 
performance::check_overdispersion(mod.patho.prob) 
car::Anova(mod.patho.prob)
summary(mod.patho.prob)

emm_type.patho.prob <- emmeans::emmeans(mod.patho.prob, ~ Type)
tukey_type.patho.prob <- pairs(emm_type.patho.prob, adjust = "tukey")
summary(tukey_type.patho.prob)

rm(emm_type.patho.prob, tukey_type.patho.prob)

###### G 2.2.1b - Anthropogenic abundance and richness  ---------------------------------------------------

mod.waste.abund.prob <- glmmTMB(cbind(abund.waste, reads - abund.waste) ~ 
                                  Type + Distance.sea + log(buildings + 1),
                                family=betabinomial(link = "logit"), data=stations)

performance::check_model(mod.waste.abund.prob) 
performance::check_overdispersion(mod.waste.abund.prob) 
car::Anova(mod.waste.abund.prob)
summary(mod.waste.abund.prob)

emm_type.waste.abund.prob <- emmeans::emmeans(mod.waste.abund.prob, ~ Type)
tukey_type.waste.abund.prob <- pairs(emm_type.waste.abund.prob, adjust = "tukey")
summary(tukey_type.waste.abund.prob)

rm(emm_type.waste.abund.prob, tukey_type.waste.abund.prob)

#
mod.waste.prob <- glm(cbind(richness.waste, richness - richness.waste) ~ 
                        Type + reads + Distance.sea + log(buildings + 1),
                      family=binomial(logit), data=stations)

car::Anova(mod.waste.prob)
summary(mod.waste.prob)

emm_type.waste.prob <- emmeans::emmeans(mod.waste.prob, ~ Type)
tukey_type.waste.prob <- pairs(emm_type.waste.prob, adjust = "tukey")
summary(tukey_type.waste.prob)

rm(emm_type.waste.prob, tukey_type.waste.prob)

###### G 2.2.1c - Environment abundance and richness ------------------------------------

mod.envi.abund.prob <- glmmTMB(cbind(abund.environment, reads - abund.environment) ~ 
                                  Type + Distance.sea + log(buildings + 1),
                                family=betabinomial(link = "logit"), data=stations)

performance::check_model(mod.envi.abund.prob) 
performance::check_overdispersion(mod.envi.abund.prob) 
car::Anova(mod.envi.abund.prob)
summary(mod.envi.abund.prob)

emm_type.envi.abund.prob <- emmeans::emmeans(mod.envi.abund.prob, ~ Type)
tukey_type.envi.abund.prob <- pairs(emm_type.envi.abund.prob, adjust = "tukey")
summary(tukey_type.envi.abund.prob)

#
mod.env.prob <- glm(cbind(richness.environment, richness - richness.environment) ~ 
                      Type + reads +  Distance.sea + log(buildings+1),
                    family=binomial(logit), data=stations)

car::Anova(mod.env.prob)
summary(mod.env.prob)

emm_type.env.prob <- emmeans::emmeans(mod.env.prob, ~ Type)
tukey_type.env.prob <- pairs(emm_type.env.prob, adjust = "tukey")
summary(tukey_type.env.prob)

###### G 2.2.1d - Heatmap difference richness and abundace per habitat ----------------------------------------------

source(file.path(fun_dir, "heatmap.R"))

h.abund_patho <- mk_heat(mod.patho.abund.prob,
                         title = "",
                         show_x = FALSE,
                         show_y = TRUE)

h.abund_waste <- mk_heat(mod.waste.abund.prob,
                         title = "",
                         show_x = FALSE,
                         show_y = FALSE)

h.abund_env <- mk_heat(mod.envi.abund.prob,
                       title = "",
                       show_x = FALSE,
                       show_y = FALSE)

h.rich_patho <- mk_heat(mod.patho.prob,
                        show_x = TRUE,
                        show_y = TRUE)

h.rich_waste <- mk_heat(mod.waste.prob,
                        show_x = TRUE,
                        show_y = FALSE)

h.rich_env <- mk_heat(mod.env.prob,
                      show_x = TRUE,
                      show_y = FALSE)

###### G 2.2.2a - Beta diversity pathogens  ---------------------------------

beta.patho <- BAT::beta(comm.pathogens, abund=T) 
print(beta.patho)

mean(beta.patho$Btotal)
mean(beta.patho$Brepl / beta.patho$Btotal)
mean(beta.patho$Brich / beta.patho$Btotal)


(pb.patho <- beta_density_plot(comm.pathogens, stations$Type, beta.patho,
                             fill_values = c("Within habitat" = "#1b9e77",
                                             "Across habitats" = "#d95f02"),
                             include_total = TRUE,
                             include_nestedness = FALSE,
                             title = "",
                             facet_cols = 1,
                             alpha = 0.6))

## Total beta diversity
model.Btotal.patho <- beta.patho$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxT.patho <- vegan::adonis2(model.Btotal.patho, stations, permutations = 9999, by = "terms"))

## Beta nestness
model.Bnest.patho <- beta.patho$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi.patho <- vegan::adonis2(model.Bnest.patho, stations, permutations = 9999, by = "terms"))

## Beta turnover
model.Bturn.patho <- beta.patho$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe.patho <- vegan::adonis2(model.Bturn.patho, stations, permutations = 9999, by = "terms"))

rm(model.Btotal.patho,model.Bnest.patho,model.Bturn.patho,perm.taxRi.patho,perm.taxRe.patho,perm.taxT.patho)

###### G 2.2.2b - Beta diversity Anthropogenic  ---------

beta.waste <- BAT::beta(comm.waste, abund=T) 

mean(beta.waste$Btotal)
mean(beta.waste$Brepl / beta.waste$Btotal)
mean(beta.waste$Brich / beta.waste$Btotal)


(pb.waste <- beta_density_plot(comm.waste, stations$Type, beta.waste,
                             fill_values = c("Within habitat" = "#1b9e77",
                                             "Across habitats" = "#d95f02"),
                             include_total = TRUE,
                             include_nestedness = FALSE,
                             title = "",
                             facet_cols = 1,
                             alpha = 0.6))

## Total beta diversity
model.Btotal.waste <- beta.waste$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxT.waste <- vegan::adonis2(model.Btotal.waste, stations, permutations = 9999, by = "terms"))

## Beta nestness
model.Bnest.waste <- beta.waste$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi.waste <- vegan::adonis2(model.Bnest.waste, stations, permutations = 9999, by = "terms"))

## Beta turnover
model.Bturn.waste <- beta.waste$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe.waste <- vegan::adonis2(model.Bturn.waste, stations, permutations = 9999, by = "terms"))

rm(model.Btotal.waste,model.Bnest.waste,model.Bturn.waste,perm.taxRi.waste,perm.taxRe.waste,perm.taxT.waste)

###### G 2.2.2c - Beta diversity Environmental  ---------

beta.envi <- BAT::beta(comm.environmental, abund=T)

mean(beta.envi$Btotal)
mean(beta.envi$Brepl / beta.envi$Btotal)
mean(beta.envi$Brich / beta.envi$Btotal)

(pb.envi <- beta_density_plot(comm.environmental, stations$Type, beta.envi,
                               fill_values = c("Within habitat" = "#1b9e77",
                                               "Across habitats" = "#d95f02"),
                               include_total = TRUE,
                               include_nestedness = FALSE,
                              title = "",
                               facet_cols = 1,
                               alpha = 0.6))

## Environmental
model.Btotal.envi <- beta.envi$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxT.envi <- vegan::adonis2(model.Btotal.envi, stations, permutations = 9999, by = "terms"))

## Beta nestness
model.Bnest.envi <- beta.envi$Brich ~  Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRi.envi <- vegan::adonis2(model.Bnest.envi, stations, permutations = 9999, by = "terms"))

## Beta turnover
model.Bturn.envi <- beta.envi$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads
(perm.taxRe.envi <- vegan::adonis2(model.Bturn.envi, stations, permutations = 9999, by = "terms"))

rm(model.Btotal.envi,model.Bnest.envi,model.Bturn.envi,perm.taxRi.envi,perm.taxRe.envi,perm.taxT.envi)

###### G 2.2.2d - ASV nMDS visualization  -----------------------------------------------------

library(vegan)

lvl_order <- c("cave","pool","salt","sea","well","pond")
disp_lab <- c(cave="cave",pool="anc.p",salt='sltw',sea="sea",well="well",pond='spr')

stations_ord <- stations |> mutate(Type = fct_relevel(Type, lvl_order))
rownames(stations_ord) <- stations_ord$Sample

my_palette <- c(
  cave = "#4B4B4B",   
  pool = "#1B9E77",   
  salt = "#E6AB02",   
  sea  = "#1F78B4",   
  well = "#A6761D",   
  pond = "#66C2A5")

# Plot all
nmds.p <- plot_nmds(
  comm.pathogens,
  stations_ord,
  my_palette,
  "",
  point_size = 2.5,
  label_size = 0,
  label_vjust = -1
) +
  scale_colour_manual(
    values = my_palette[levels(stations_ord$Type)],
    breaks = levels(stations_ord$Type),
    labels = disp_lab[levels(stations_ord$Type)],
    name = "Habitat"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    strip.text.x = element_text(size = 14)
  )

(nmds.w <- plot_nmds(comm.waste,
                     stations_ord,
                    my_palette,
                    "",point_size = 2.5,
                    label_size = 0,
                    label_vjust = -1)+
    scale_colour_manual(
      values = my_palette[levels(stations_ord$Type)],
      breaks = levels(stations_ord$Type),
      labels = disp_lab[levels(stations_ord$Type)],
      name = "Habitat"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      strip.text.x = element_text(size = 14)
    ))

(nmds.e <- plot_nmds(comm.environmental,
                     stations_ord,
                    my_palette,
                    "",
                    point_size = 2.5,
                    label_size = 0,
                    label_vjust = -1)+
    scale_colour_manual(
      values = my_palette[levels(stations_ord$Type)],
      breaks = levels(stations_ord$Type),
      labels = disp_lab[levels(stations_ord$Type)],
      name = "Habitat"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 12, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      strip.text.x = element_text(size = 14)
    ))
                    
########## G 2.2.3 - FIGURE FUNCTIONAL  -------------------------------------

# -------- Richness --------
p_rich <-
  (h.rich_patho | h.rich_waste | h.rich_env)

# -------- Abundance --------
p_abund <-
  (h.abund_patho | h.abund_waste | h.abund_env)

# -------- Beta-diversity --------
p_beta <-
  (pb.patho | pb.waste | pb.envi)

# -------- NMDS --------
p_nmds <-
  (nmds.p | nmds.w | nmds.e)

# -------- COMBINE ALL --------
fig_all <-
  (p_beta /
    p_rich /
      p_abund /
      p_nmds
  ) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0, 1)
  )

fig_all

# PNG
ggsave("Figure_5.png", plot = fig_all, width = 12, height = 18, units = "in", dpi = 600)

# TIFF
ggsave("Figure_5.tiff", plot = fig_all, width = 18, height = 11, units = "cm", dpi = 600,
  compression = "lzw")

########## G 2.3 - Classification  -----------------------------------------------------

########## G 2.3.1 - RF Entire community-------------------------------------------

library(randomForest)
library(caret)
library(compositions)

source(file.path(fun_dir, "randomforest_onebyone.R"))

rownames(stations) <- stations$Sample

rf.all <- run_rf_one_vs_rest_by_rank(
  comm.tax.list = comm.tax.list,
  stations      = stations,
  target_col    = "Type",
  ranks         = c("Family","Genus","Species"),
  ntree         = 10000,
  balance       = TRUE,
  plot          = TRUE
)

########## G 2.3.2 - RF Pathogens community-------------------------------------------

comm.pathogens.tax <- as.data.frame(t(comm.pathogens))
comm.pathogens.tax$ID <- row.names(comm.pathogens.tax)
comm.pathogens.tax <- merge(taxonomy, comm.pathogens.tax, by = "ID")

comm.patho.list <- summarize_by_taxonomic_rank(comm.pathogens.tax,
                                               c("Family", "Genus", "Species"),
                                               min_reads = 20)

rf.pathogens <- run_rf_one_vs_rest_by_rank(
  comm.tax.list = comm.patho.list,
  stations      = stations,
  target_col    = "Type",
  ranks         = c("Family","Genus","Species"),
  ntree         = 10000,
  balance       = TRUE,
  plot          = TRUE
)

########## G 2.3.3 - RF Anthropogenic community -------------------------------------------

comm.waste.tax <- as.data.frame(t(comm.waste))
comm.waste.tax$ID <- row.names(comm.waste.tax)
comm.waste.tax <- merge(taxonomy, comm.waste.tax, by = "ID")

comm.waste.list <- summarize_by_taxonomic_rank(comm.waste.tax,
                                               c("Family", "Genus", "Species"),
                                               min_reads = 20)

rf.waste <- run_rf_one_vs_rest_by_rank(
                comm.tax.list = comm.waste.list,
                stations      = stations,
                target_col    = "Type",
                ranks         = c("Family","Genus","Species"),
                ntree         = 10000,
                balance       = TRUE,
                plot          = TRUE
              )
               
########## G 2.3.4 -  RF Environmental -------------------------------------------

comm.environmental.tax <- as.data.frame(t(comm.environmental))
comm.environmental.tax$ID <- row.names(comm.environmental.tax)
comm.environmental.tax <- merge(taxonomy, comm.environmental.tax, by = "ID")

comm.environmental.list <- summarize_by_taxonomic_rank(comm.environmental.tax,
                                                       c("Family", "Genus", "Species"),
                                                       min_reads = 20)

rf.environmental <- run_rf_one_vs_rest_by_rank(
                        comm.tax.list = comm.environmental.list,
                        stations      = stations,
                        target_col    = "Type",
                        ranks         = c("Family","Genus","Species"),
                        ntree         = 10000,
                        balance       = TRUE,
                        plot          = TRUE
                      )
                        
########### G 2.3.5 -  RF Visualization -------------------------------------------

library(tibble)
library(purrr)

col_low  <- "#d7191c" 
col_mid  <- "white"
col_high <- "#2c7bb6" 

hab_order   <- c("cave","pool","salt","sea","well")   
rank_order  <- c("Family","Genus","Species")
group_order <- c("All","Pathogens","Anthropogenic","Environmental")

# Build tidy metrics table
m_all <- bind_rows(
  bind_metrics(rf.all,           "All"),
  bind_metrics(rf.pathogens,     "Pathogens"),
  bind_metrics(rf.waste,         "Anthropogenic"),
  bind_metrics(rf.environmental, "Environmental")
) %>%
  # keep known habitats and drop pond in one go
  filter(habitat %in% c(hab_order, "pond")) %>%
  filter(habitat != "pond") %>%
  mutate(
    group    = factor(group,   levels = group_order),
    habitat  = factor(habitat, levels = hab_order),
    rank     = factor(rank,    levels = rank_order),
    Accuracy = 1 - OOB
  ) %>%
  # drop rows where rank not in rank_order (e.g., Order)
  filter(!is.na(rank)) %>%
  arrange(group, habitat, rank)


# AUC heatmap
(p_auc <- ggplot(m_all, aes(habitat, rank, fill = AUC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(AUC), "NA", sprintf("%.2f", AUC))), size = 3) +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0.85, limits = c(0.6, 1),
                       breaks = c(0, 0.85, 1), labels = c("0", "0.85", "1"),
                       oob = scales::squish, name = "AUC") +
  labs(title = "A. ROC AUC", x = NULL, y = NULL) +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank()))


# Accuracy heatmap
(p_acc <- ggplot(m_all, aes(habitat, rank, fill = Accuracy)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Accuracy)), size = 3) +
  scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high,
                       midpoint = 0.80, limits = c(0.6, 1),
                       breaks = c(0, 0.80, 1), labels = c("0", "0.80", "1"),
                       oob = scales::squish, name = "Accuracy") +
  labs(title = "B. Accuracy (1 − OOB)", x = NULL, y = NULL) +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank()))

final.rf <-
  (p_auc / p_acc) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 9),
    axis.text       = element_text(size = 9),
    strip.text      = element_text(size = 9)  # facet labels
  )

final.rf

##### SUPPLEMENTARY ANALYSES --------------------------------------- 

###### S 1 - Dendrograms  -------------------------------------------

# Set up a 2x2 plotting layout
par(mfrow = c(2, 2))

# Plot the four dendrograms
plot(hclust(beta.tax$Btotal, method = "average"), main = "All Taxa", xlab = "", sub = "")
plot(hclust(beta.patho$Btotal, method = "average"), main = "Pathogenic Taxa", xlab = "", sub = "")
plot(hclust(beta.waste$Btotal, method = "average"), main = "Wastewater Taxa", xlab = "", sub = "")
plot(hclust(beta.envi$Btotal, method = "average"), main = "Environmental Taxa", xlab = "", sub = "")

rm(comm.environmental,comm.environmental.tax,comm.environmental.list,
   comm.pathogens,comm.pathogens.tax,comm.patho.list,
   comm.waste,comm.waste.tax,comm.waste.list)

########## S 2 - Richness within the La Corona  -----------------------------------------------------
 
 ## Test for a gradient of antropogenic pollution in the cave
 
 stations.cave <- stations[ which (stations$Sample %in% c("TDA_sima",
                                                          "TDA_entrada",
                                                          "TDA_LE",
                                                          "JDA_pool1",
                                                          "JDA_pool2",
                                                          "JDA_pool3",
                                                          "CUL1",
                                                          "CUL2",
                                                          "CHL",
                                                          "CHZ") ),]
 
 
 # stations.cave$Distance.center <- c(5250,1800,250,50,0,0,0,0,50,250)
 
 mod.rich.cave <- MASS::glm.nb(richness ~ Light + reads +  Distance.sea,
                               data = stations.cave)
 
 performance::check_model(mod.rich.cave) 
 performance::check_overdispersion(mod.rich.cave) 
 car::Anova(mod.rich.cave)
 summary(mod.rich.cave)
 
 
 emm_light <- emmeans::emmeans(mod.rich.cave, ~ Light)
 tukey_light <- pairs(emm_light, adjust = "tukey")
 summary(tukey_light)
 
 boxplot(richness ~ Light, data = stations.cave)
 
 plot(richness ~ Distance.sea, data = stations.cave)
 
 
 rm(mod.rich.cave, emm_light, tukey_light)
 
########## S 2.1 - Composition within the La Corona  -----------------------------------------------------
 
 comm.cave <- comm[ which (rownames(comm) %in% stations.cave$Sample),]
 comm.cave <- comm.cave[ , colSums(comm.cave) > 0 ]
 
 beta.tax.cave <- BAT::beta(comm.cave, abund=T) 
 
 mean(beta.tax.cave$Btotal)
 mean(beta.tax.cave$Brepl / beta.tax.cave$Btotal)
 mean(beta.tax.cave$Brich / beta.tax.cave$Btotal)
 
 (pb.cave <- beta_density_plot(comm.cave, stations.cave$Type, beta.tax.cave,
                               fill_values = c("Within habitat" = "#1b9e77",
                                               "Across habitats" = "#d95f02"),
                               alpha = 0.6))
 
 ## Total beta diversity
 model.Btotal.cave <- beta.tax.cave$Btotal ~ Light + reads +  Distance.sea
 (perm.taxT.cave <- vegan::adonis2(model.Btotal.cave, stations.cave, permutations = 9999, by = "terms"))
 
 ## Beta nestness
 model.Bnest.cave <- beta.tax.cave$Brich ~  Light + reads +  Distance.sea
 (perm.taxRi.cave <- vegan::adonis2(model.Bnest.cave, stations.cave, permutations = 9999, by = "terms"))
 
 ## Beta turnover
 model.Bturn.cave <- beta.tax.cave$Brepl ~ Light + reads +  Distance.sea
 (perm.taxRe.cave <- vegan::adonis2(model.Bturn.cave, stations.cave, permutations = 9999, by = "terms"))
 
 rm(model.Btotal.cave,model.Bnest.cave,model.Bturn.cave,perm.taxRi.cave,perm.taxRe.cave,perm.taxT.cave)
 
 rm(comm.cave,stations.cave)
 
########## S 3 - Phylogenetic diversity  -----------------------------------------------------
 

########## S 3.1 - Calculate phylogenetic tree  -----------------------------------------------------

  
 # ## Read alignment
 # alignment.all <- ape::read.FASTA("phylogenetic_diversity/taxonomy.check_Mafft.fasta")
 # 
 # ## Calculate distance matrix
 # dist <- ape::dist.dna(alignment.all, pairwise.deletion = TRUE) 
 # 
 # ## Compute tree
 # tree <- ape::njs(dist)
 # tree <- phangorn::midpoint(tree)
 # ape::write.tree(tree,"phylogenetic_diversity/Lanzarote_16s.tre")
 # 
 # ### Annotate tree
 # taxonomy$label.long <- paste0(taxonomy$Phylum,"|",
 #                               taxonomy$Class,"|",
 #                               taxonomy$Order, "|",
 #                               taxonomy$Family,"|",
 #                               taxonomy$ID)
 # 
 # tree.names <- suppressWarnings(phylotools::sub.taxa.label(tree, taxonomy[c("ID", "label.long")]))
 # ape::write.tree(tree.names,"phylogenetic_diversity/Lanzarote_16s_names.tre")
 
 
########## S 3.2 - Calculate phylogenetic diversity  -----------------------------------------------------
 

 tree <- ape::read.tree("Lanzarote_16s.tre")
 
 comm_pa <- (comm > 0) * 1

 pd_sample <- BAT::alpha(comm_pa, tree = tree) 
 pd_sample <- as.data.frame(pd_sample)
 colnames(pd_sample) <- "phylo.diver"
 pd_sample$Sample <- rownames(pd_sample)
 
 stations <- merge(stations, pd_sample, by = "Sample")
 rownames(stations) <- stations$Sample
 
 rm(pd_sample)
 
 
 ########## S 3.3 - Phylogenetic richness  -----------------------------------------------------
 
 plot(stations$phylo.diver ~ stations$richness)
 
 mod.phylo <- glm(phylo.diver ~ Type + reads + Distance.sea + log(buildings + 1),
                          data = stations,
                          family = Gamma(link="log"))
 
 performance::check_model(mod.phylo) 
 
 
 #performance::check_overdispersion(mod.phylo) 
 car::Anova(mod.phylo)
 summary(mod.phylo)
 
 emm_type <- emmeans::emmeans(mod.phylo, ~ Type)
 tukey_type <- pairs(emm_type, adjust = "tukey")
 summary(tukey_type)
 
 ########## S 3.4 - Phylogenetic composition  -----------------------------------------------------
 
 beta.phylo <- BAT::beta(comm_pa, tree = tree)


 mean(beta.phylo$Btotal)
 mean(beta.phylo$Brepl / beta.phylo$Btotal)
 mean(beta.phylo$Brich / beta.phylo$Btotal)
 
 (pb.phylo <- beta_density_plot(comm_pa, stations$Type, beta.phylo,
                                fill_values = c("Within habitat" = "#1b9e77",
                                                "Across habitats" = "#d95f02"),
                                alpha = 0.6))
 
 
 ## Total beta diversity
 model.Bphtotal <- beta.phylo$Btotal ~ Type + Distance.sea + log(buildings + 1) + reads + richness
 (perm.phyT <- vegan::adonis2(model.Bphtotal, stations, permutations = 9999, by = "terms"))
 
 ## Beta nestness
 model.Bphnest <- beta.phylo$Brich ~  Type + Distance.sea + log(buildings + 1) + reads + richness
 (perm.phyRi <- vegan::adonis2(model.Bphnest, stations, permutations = 9999, by = "terms"))
 
  ## Beta turnover
 model.Bphturn <- beta.phylo$Brepl ~ Type + Distance.sea + log(buildings + 1) + reads + richness
 (perm.phyRe <- vegan::adonis2(model.Bphturn, stations, permutations = 9999, by = "terms"))
 
 rm(model.Bphtotal, model.Bphnest, model.Bphturn, perm.phyRi, perm.phyRe, perm.phyT) 
 
  
 
 
 
 
 
