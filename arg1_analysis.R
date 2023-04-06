library(limma)
library(tidyverse)
library(illuminaio)
library(magrittr)

# read in data ------------------------------------------------------------
files <- list.files(path = "~/My Drive/analysis_for_people/erin/arg1_patient/raw_data",
                    pattern = "idat",
                    full.names = T)

df <- read.idat(idatfiles = files, 
                bgxfile = "~/My Drive/analysis_for_people/erin/arg1_patient/raw_data/HumanHT-12_V4_0_R2_15002873_B")


# create metadata ---------------------------------------------------------
meta <- data.frame(files = files,
                   sample = 1:length(files),
                   status = rep(c("healthy", "arg1_def"), each = 6),
                   stimulation = rep(c("none", "cd3", "cd3_46"), times = 4))


# normalisation and filtering ---------------------------------------------
boxplot(log2(df$E), range = 0, ylab = "log2 intensity", las = 2,
        names = paste0("sample", 1:12))

pvals <- detectionPValues(df)
df$other$Detection <- pvals

norm_data <- neqc(df)
expressed <- rowSums(norm_data$other$Detection < 0.05) >= 2 # any probes expressed in >=2 samples are kept
norm_data <- norm_data[expressed, ]


# export the normalised data ----------------------------------------------
sample_names <- paste(meta$status, meta$stimulation, sep = "_")
norm_output <- norm_data$E %>%
  set_colnames(sample_names)

# get gene annotations
annots <- norm_data$genes
annots$Array_Address_Id <- as.character(annots$Array_Address_Id)

# combine gene annotations and expression file
norm_output <- norm_output %>%
  data.frame() %>%
  rownames_to_column("Array_Address_Id") %>%
  full_join(annots, by = "Array_Address_Id") %>%
  select(-c(Probe_Id, Array_Address_Id)) %>%
  select(Symbol, tidyselect::peek_vars()) %>%
  rename("gene" = "Symbol") %>%
  janitor::clean_names()

# export expression file
write.csv(norm_output, "~/My Drive/analysis_for_people/erin/arg1_patient/normalised_data.csv",
          row.names = F)

# check the normalised data -----------------------------------------------
boxplot(norm_data$E, range = 0, ylab = "log2 intensity", las = 2,
        names = paste0("sample", 1:12))

# export data for gsea ----------------------------------------------------
norm_output %>%
  filter(gene != "") %>%
  select(gene, healthy_cd3_46, healthy_cd3_46_1, arg1_def_cd3_46, arg1_def_cd3_46_1) %>%
  mutate(DESCRIPTION = 1) %>%
  rename(NAME = gene) %>%
  select(NAME, DESCRIPTION, tidyselect::peek_vars()) %>%
  write.table("~/My Drive/analysis_for_people/erin/arg1_patient/gsea/stim_gsea.txt",
              quote = F, sep = "\t", row.names = F)


# visualise pathways ------------------------------------------------------
stim_gsea <- rbind(
  read.delim("~/My Drive/analysis_for_people/erin/arg1_patient/gsea/stim_analysis.Gsea.1680798280344/gsea_report_for_healthy_1680798280344.tsv") %>%
    select(NAME, NES, FDR.q.val) %>%
    mutate(up = "healthy"),
  
  read.delim("~/My Drive/analysis_for_people/erin/arg1_patient/gsea/stim_analysis.Gsea.1680798280344/gsea_report_for_arg1_1680798280344.tsv") %>%
    select(NAME, NES, FDR.q.val) %>%
    mutate(up = "arg1")
)

stim_gsea %>% arrange(FDR.q.val) %>% head(40)

highlight_paths <- stim_gsea %>%
  filter(NAME %in% c("REACTOME_INTERLEUKIN_10_SIGNALING", "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES"))
highlight_sig <- stim_gsea %>% filter(FDR.q.val < 0.05)

ggplot(stim_gsea, aes(NAME, abs(NES))) +
  geom_hline(yintercept = 2, col = "gray50", linetype = "dashed") +
  geom_point(shape = 21, fill = "#cfe3d5", size = 2, alpha = 0.4, stroke = 0.1) +
  geom_point(data = highlight_sig, shape = 21, fill = "#fae3e1", size = 3, stroke = 0.3) +
  geom_point(data = highlight_paths, shape = 21, fill = "#34c6eb", size = 4, stroke = 0.3) +
  ggrepel::geom_label_repel(data = highlight_paths, 
                            label = gsub("_", " ", highlight_paths$NAME),
                            size = 2.2,
                            nudge_y = ifelse(highlight_paths$NAME == "REACTOME_INTERLEUKIN_10_SIGNALING", -0.4, 0.6), 
                            nudge_x = ifelse(highlight_paths$NAME == "REACTOME_INTERLEUKIN_10_SIGNALING", 400, 200),
                            point.padding = 0.55,
                            arrow = arrow(angle = 20, length = unit(2, "mm"))) +
  labs(x = "Pathway", y = "Absolute NES") +
  scale_x_discrete(expand = c(0.03, 0.03)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = "~/My Drive/analysis_for_people/erin/arg1_patient/gsea_output.png", 
       dpi = 600, width = 3.5, height = 3)









