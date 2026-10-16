source("Scripts/00_packages_functions.R")

mmcount_wide <- readRDS("Scratch_data/mmcount_wide.Rds")
cohort4 <- readRDS("Scratch_data/cohort4.Rds")

## add comorbidities to cohort information ----

cohort5 <- cohort4 %>% 
  left_join(mmcount_wide %>% select(ALF_PE, MIX_DT, mm_count))

saveRDS(cohort5, "Scratch_data/cohort5.Rds")
