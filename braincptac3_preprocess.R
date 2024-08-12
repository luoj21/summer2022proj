setwd("/Users/jasonluo/Documents/NIH/Summer2022_Work/NIH_Summer_2022/braincptac3_repo")
library(tidyverse)

cptac3_gdc_sample <- read_tsv("gdc_sample_sheet.2022-06-04.tsv")
cptac3_gdc_sample1 <- cptac3_gdc_sample %>%
  select(`File ID`, `File Name`, `Project ID`, `Sample ID`, `Sample Type`, `Case ID`) %>%
  mutate(`Case ID` = gsub('\\,.*', '', `Case ID`)) %>%
  mutate(`Case ID`, case_id = if_else(grepl('Primary Tumor', `Sample Type`)
                                      , gsub('$', '-Tu', `Case ID`), gsub('$', '-No', `Case ID`))) %>%
  select(-c('Case ID', 'Sample Type'))

write_csv(cptac3_gdc_sample1, file = "braincptac3SAMPLE.csv")


cptac3_gdc_manifest <- read_tsv("cptac3_manifest.tsv")
cptac3_gdc_manifest1 <- cptac3_gdc_manifest %>%
  drop_na(gene_name, gene_type) %>%
  rename(file_id_name = `04a95ed6-a66b-4a95-9b12-da3b48ab7ff6/0ab80f6a-8139-42a9-886c-e8eb5cc9de88.rna_seq.augmented_star_gene_counts.tsv:gene_id`) %>%
  mutate(EID = gsub('.*\\:', '', file_id_name )) %>%
  mutate(file_id_name = gsub('\\:.*', '', file_id_name)) %>%
  mutate(file_name  = gsub('.*\\/', '', file_id_name)) %>%
  rename(file_id = file_id_name) %>%
  mutate(file_id = gsub('\\/.*', '', file_id)) %>%
  select(file_id, file_name, EID, gene_name, gene_type, fpkm_unstranded)

write_csv(cptac3_gdc_manifest1, file = "braincptac3RNAmanifest.csv")


cptac3_proteome <- read_tsv("CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv")
cptac3_proteome1_LR <- cptac3_proteome %>%
  subset(Gene != "Mean" & Gene != "Median" & Gene != "StdDev") %>%
  select(-ends_with("Unshared Log Ratio")) %>%
  select(contains("Gene") | contains("Log Ratio")) %>%
  select(-ends_with("NCBIGeneID")) %>%
  pivot_longer(!Gene, names_to = "patientID", values_to = "value") %>%
  mutate(patientID = gsub('Log Ratio', '', patientID)) %>%
  mutate(patientID = trimws(as.character(patientID)))

proteome_biospecimen <- read_tsv("PDC_biospecimen_manifest_06082022_102638.tsv")
proteome_biospecimen1 <- proteome_biospecimen %>%
  select(`Aliquot Submitter ID`, `Case Submitter ID` ,`Sample Type`) %>%
  rename(case_id = `Case Submitter ID`, aliquot = `Aliquot Submitter ID`, sample_type = `Sample Type`) %>%
  mutate(case_id = if_else(grepl('Primary Tumor', sample_type), gsub('$', '-Tu', case_id), gsub('$', '-No', case_id))) %>%
  select(-c(sample_type)) %>%
  mutate(aliquot = trimws(as.character(aliquot)))

cptac3_proteome_joined_shared <- cptac3_proteome1_LR %>%
  left_join(proteome_biospecimen1, by = c("patientID" = "aliquot")) %>%
  select(-c(patientID)) %>%
  select(Gene, case_id, value) %>%
  na.omit

which(duplicated(cptac3_proteome_joined_shared))
write_csv(cptac3_proteome_joined_shared, file = "braincptac3prodata_shared.csv")

manifest_RNA <- read_csv("braincptac3RNAmanifest.csv")
sample_sheet <-  read_csv("braincptac3SAMPLE.csv")

braincptac3RNA <- manifest_RNA %>%
  left_join(sample_sheet, by = c("file_id" = "File ID", "file_name" = "File Name")) %>%
  filter(gene_type == "protein_coding") %>%
  select(gene_name, case_id, fpkm_unstranded) %>%
  mutate(Tumor = rep(c("GBM")))


which(duplicated(braincptac3RNA))
write_csv(braincptac3RNA, file = "braincptac3RNA.csv")



