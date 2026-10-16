source("Scripts/00_packages_functions.R")

cohort6 <- readRDS("Scratch_data/cohort6.Rds")
length(unique(cohort6$ALF_PE)) #

# n= and % individuals with only one creatinine value in relevant time period

cohort6 %>% 
  summarise(n = n(), 
            n1creat = sum(n_creat == 1),
            pc1creat = n1creat/n()*100, 
            morecreat = sum(n_creat >1),
            pcmorecreat =morecreat/n()*100 ) %>% 
  write_csv("Outputs/n_creat.csv")

# baseline data table 1 ----

base <- c("age", "sex_name", "egfr", "ckd_cat", "mm_count", "n_cancers", "m_stage", "tx_smry")
cat_base <- c("sex_name", "ckd_cat", "m_stage", "tx_smry")

# sex-stratified baseline table  ----
table1 <- CreateTableOne(vars = base, 
                             factorVars = cat_base, 
                             strata = c("cancer_type", "sex_name"),
                             addOverall = TRUE,
                             data = cohort6)
print(table1, 
      nonnormal = c("age", "egfr", "mm_count", "n_cancers"),
      contDigits = 1) %>% 
  write.csv("Outputs/table1.csv")
  