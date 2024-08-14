setwd("/Users/jasonluo/Documents/NIH/Summer2022_Work/NIH_Summer_2022/braincptac3_repo") # change to appropriate folder
library(tidyverse)
library(ggpubr)
library(skimr)
library(EnvStats)
library(fuzzyjoin)
library(stringr)
library(survival)
library(survminer)

manifest_RNA <- read_csv("braincptac3RNAmanifest.csv")
RNA <- read_csv("braincptac3RNA.csv")
protein_shared <- read_csv("braincptac3prodata_shared.csv")
sample_sheet <-  read_csv("braincptac3SAMPLE.csv")
brain_clinical <- read_tsv("brain_clinical.tsv")

gene_subset <- protein_shared %>%
  filter(Gene %in% c("RCC1", "RANGRF", "RANGAP1", "RANBP2", "RANBP1")) %>%
  mutate(type = if_else(grepl('-Tu', case_id), "Primary Tumor", "Solid Tissue Normal"))

## Boxplots
ggplot(gene_subset, aes(Gene, value, fill = type)) + 
  geom_boxplot(color = "black", alpha = 0.5, width = 1.1) + 
  facet_grid(~ Gene, scale = "free_x") + 
  geom_jitter(aes(color = type, shape = type), size = 1, alpha=0.9) +
  ylab("Protein Abundance (Shared log2 Ratio)") +
  stat_compare_means(method = "wilcox.test", aes(label = paste0("          p = ", ..p.format..)), size = 4) +
  theme_bw() +
  #ggtitle("Protein Abundance Distribution in Brain Cancer Patients") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        axis.title = element_text(size = 15), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 17),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)
  ) +
  stat_n_text(aes(type), size = 4)
 
 
## Distribution of protein abundances
ggplot(gene_subset, aes(value, color = type)) + 
   geom_density() + 
   facet_wrap(~ Gene, scale = "free_x") + 
   theme(strip.text.x = element_text(size = 10, face = "bold"),
         legend.text = element_text(size = 12), 
         legend.position = "top",
         axis.title.x = element_text(face = "bold")) +
   xlab("Protein Abundance (Shared log2 Ratio)")
 
 
## Linear Regression
filtered_manifest_RNA <- manifest_RNA %>%
   left_join(sample_sheet, by = c("file_id" = "File ID", "file_name" = "File Name")) %>%
   select(fpkm_unstranded, case_id, gene_name) %>%
   inner_join(protein_shared, by = c("case_id", "gene_name" = "Gene")) %>%
   mutate(type = if_else(grepl('-Tu', case_id), "Primary Tumor", "Solid Tissue Normal"))
 
filtered_manifest_RNA <- filtered_manifest_RNA %>%
   fuzzy_left_join(brain_clinical, match_fun = str_detect, by = c("case_id" = "case_submitter_id")) %>%
   select(fpkm_unstranded, gene_name, case_submitter_id, value, gender) %>%
   drop_na()
 
ggplot(subset(filtered_manifest_RNA, gene_name %in% c("RANBP1", "RANGRF")), aes(fpkm_unstranded, value, color = gene_name)) +
   geom_point(aes(shape = gene_name)) + 
   geom_smooth(method = "lm", se = TRUE) + 
   stat_regline_equation() + 
   stat_cor(aes(label = ..r.label..), label.x = 18) + 
   facet_wrap(~ gender, scale = "free_y") + 
   labs(x = "mRNA Abundance (FPKM)", y = "Protein Abundance (Shared log2 Ratio)") + 
   theme_bw() + 
   theme(strip.text.x = element_text(size = 15, face = "bold"),
         axis.text = element_text(size = 10))
 
 
## Survival Analysis
## Base survival curve based on gender
brain_survfit1 <- brain_clinical %>%
   select(case_submitter_id, gender, vital_status, days_to_last_follow_up, days_to_death) %>%
   filter(!grepl("Not Reported", vital_status)) %>%
   mutate(time_recent_visit = if_else(days_to_death > days_to_last_follow_up, days_to_death, days_to_last_follow_up)) %>%
   mutate(time_recent_visit = as.numeric(time_recent_visit)) %>%
   mutate(time_recent_visit = time_recent_visit / 30) %>%
   mutate(vital_status = if_else(vital_status == "Alive", 0, 1))
 
kmfit1 <- survfit(Surv(time_recent_visit, vital_status) ~ gender, data = brain_survfit1, type = "kaplan-meier")
 
ggsurvplot(kmfit1, data = brain_survfit1, conf.int = TRUE, risk.table = TRUE, pval = TRUE) + 
   xlab("Time in Months") +
   ggtitle("Brain Cancer Survival Curve")

## Survival curve based on gender and protein expression
brain_survfit2 <- protein_shared %>%
  mutate(tissue_type = if_else(grepl('-Tu', case_id), "Primary Tumor", "Solid Tissue Normal")) %>%
  mutate(case_id = gsub('.{3}$', '', case_id)) %>%
  left_join(brain_clinical, by = c("case_id" = "case_submitter_id")) %>%
  select(case_id, gender, Gene, value, vital_status, days_to_last_follow_up, days_to_death, tissue_type) %>%
  filter(!grepl("Not Reported", vital_status)) %>%
  mutate(time_recent_visit = if_else(days_to_death > days_to_last_follow_up, days_to_death, days_to_last_follow_up)) %>%
  mutate(time_recent_visit = as.numeric(time_recent_visit)) %>%
  mutate(time_recent_visit = time_recent_visit / 30) %>%
  mutate(vital_status = if_else(vital_status == "Alive", 0, 1)) %>%
  drop_na(gender) %>%
  mutate(value_status = case_when(value > median(subset(., gender == "male")$value) & gender == "male" ~ 'high-expression-male',
                                  value < median(subset(., gender == "male")$value) & gender == "male" ~ 'low-expression-male',
                                  value > median(subset(., gender == "female")$value) & gender == "female" ~ 'high-expression-female',
                                  value < median(subset(., gender == "female")$value) & gender == "female" ~ 'low-expression-female')) %>%
  mutate(case_id = if_else(tissue_type == "Primary Tumor", gsub('$', '-Tu', case_id),gsub('$', '-No', case_id)))


gene_survival_curve <- function(gene_name, states = c("low-expression-male", "high-expression-male")) { 
 
  brain_survfit2_subset <- brain_survfit2 %>%
    filter(Gene == gene_name)

  kmfit2 <- survfit(Surv(time_recent_visit, vital_status) ~ value_status, data = brain_survfit2_subset, type = "kaplan-meier")
  
  
  brain_survfit3 <- brain_survfit2_subset %>%
    filter(value_status %in% states) 
  
  kmfit3 <- survfit(Surv(time_recent_visit, vital_status) ~ value_status, data = brain_survfit3, type = "kaplan-meier")
  plot <- ggsurvplot(kmfit3, data = brain_survfit3, conf.int = TRUE, risk.table = TRUE, pval = TRUE) +
    xlab("Time in Months") +
    ggtitle(paste0(gene_name, " Gene Survival Curve For ", states))
  
  log_rank_test <- survdiff(Surv(time_recent_visit, vital_status) ~ value_status, data = brain_survfit3)
  p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
  
  print(paste0(p_value))
  return(plot)
}


gene_survival_curve("RANGRF", c("low-expression-female", "high-expression-female")) 


## Combined survival curves
kmfit_full <- survfit(Surv(time_recent_visit, vital_status) ~ value_status, 
                  data = brain_survfit2 %>% filter(Gene == "RANGRF"), # change this to whatever gene to analyze
                  type = "kaplan-meier")
p1 <- ggsurvplot(kmfit_full, 
                 linetype =  c("dashed", "solid", "dashed", "solid"), 
                 data = brain_survfit2, 
                 conf.int = FALSE, 
                 risk.table = TRUE,
                 risk.table.col = "strata",
                 risk.table.y.text = FALSE,
                 pval = FALSE, 
                 legend.title= "") + 
  xlab("Time in Months") +
  ggtitle("RANGRF Gene Brain Cancer Survival Curve")

p1$plot <- p1$plot + 
  annotate("text", x = 5, y = 0.25, label = "Male: p = 0.0057 \n Female: p = 0.87", size = 5) + # manually inputted pvalues
  theme(legend.text = element_text(size = 11))

p1

 