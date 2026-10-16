source("Scripts/00_packages_functions.R")

cohort5 <- readRDS("Scratch_data/cohort5.Rds")
kf_on_krt <- readRDS("Processed_data/kf_on_krt.Rds") %>% 
  mutate(renal_tx = if_else(grepl("N165|T861|Z940|M012|M013|M014|M015|M018|M019|M084|M172|M174|M178|M179", krt_code), 1, 0))

## 
cohort6 <- cohort5 %>% 
  left_join(kf_on_krt %>% select(ALF_PE, krt_date, krt_code, renal_tx)) %>% 
  mutate(DEATH_DT = as.Date(DEATH_DT, format = "%Y-%m-%d")) %>% 
  mutate(sex_name = if_else(sex == 1, "Male", "Female")) %>% 
  mutate(age10 = age/10) %>% 
  mutate(cancer_type = case_when(
    grepl("C50", SITE_CD) ~ "Breast",
    grepl("C61", SITE_CD) ~  "Prostate",
    grepl("C18|C19|C20", SITE_CD) ~ "Colorectal", 
    grepl("C33|C34", SITE_CD) ~ "Lung and bronchus", 
    grepl("C43", SITE_CD) ~ "Melanoma", 
    grepl("C64|C65|C66|C67|C68", SITE_CD) ~ "Renal/urinary tract"
  )) %>% 
  mutate(
    t_stage = case_when(
      DIAGSTAGE_HIST_T == "0" ~ "0",
      grepl("1", DIAGSTAGE_HIST_T) ~ "1",
      grepl("2", DIAGSTAGE_HIST_T) ~ "2",
      grepl("3", DIAGSTAGE_HIST_T) ~ "3",
      grepl("4", DIAGSTAGE_HIST_T) ~ "4",
      TRUE ~ "Unknown"),
    n_stage = case_when(
      DIAGSTAGE_HIST_N == "0" ~ "0",
      grepl("1", DIAGSTAGE_HIST_N) ~ "1",
      grepl("2", DIAGSTAGE_HIST_N) ~ "2",
      grepl("3", DIAGSTAGE_HIST_N) ~ "3",
      TRUE ~ "Unknown"), 
    m_stage = case_when(
      DIAGSTAGE_HIST_M == "0" ~ "0",
      grepl("1", DIAGSTAGE_HIST_M) ~ "1",
      TRUE ~ "Unknown")) %>% 
  mutate_at(.vars = vars(c(DT_TREAT_START_DT, RTSTART_DT, SURGDATE_OPER, DEATH_DT)),
            .funs = list(time = ~ as.numeric(. - MIX_DT))) %>% 
  mutate(surv1y = if_else(DEATH_DT_time >= 365, 1L, 0L), 
         surv5y = if_else(DEATH_DT_time >= 1826, 1L, 0L)) %>% 
  mutate(sact = if_else(is.na(sact), 0, sact), 
         rtx = if_else(is.na(rtx), 0, rtx), 
         surg = if_else(is.na(surg), 0, surg)) %>% 
  mutate(tx_smry = case_when(
    sact == 0 & rtx == 0 & surg == 0 ~ "None", 
    sact == 1 & rtx == 0 & surg == 0 ~ "SACT only",
    sact == 0 & rtx == 1 & surg == 0 ~ "RTx only",
    sact == 0 & rtx == 0 & surg == 1 ~ "Surg only",
    sact == 1 & rtx == 1 & surg == 0 ~ "Treatment combination",
    sact == 1 & rtx == 0 & surg == 1 ~ "Treatment combination", 
    sact == 0 & rtx == 1 & surg == 1 ~ "Treatment combination", 
    sact == 1 & rtx == 1 & surg == 1 ~ "Treatment combination")) %>% 
  mutate(any_tx = case_when(
    sact == 0 & rtx == 0 & surg == 0 ~ "None", 
    sact == 1 & rtx == 0 & surg == 0 ~ "Treatment",
    sact == 0 & rtx == 1 & surg == 0 ~ "Treatment",
    sact == 0 & rtx == 0 & surg == 1 ~ "Treatment",
    sact == 1 & rtx == 1 & surg == 0 ~ "Treatment",
    sact == 1 & rtx == 0 & surg == 1 ~ "Treatment", 
    sact == 0 & rtx == 1 & surg == 1 ~ "Treatment", 
    sact == 1 & rtx == 1 & surg == 1 ~ "Treatment")) %>% 
  mutate(
    krt = if_else(!is.na(krt_code), 1, 0),
    krt_pre = if_else(krt == 1 & krt_date < MIX_DT, 1, 0)) %>% 
  filter(krt_pre != 1) %>% # remove individuals on maintenance KRT 
  mutate(ckd_cat = factor(case_when(
    #krt_pre == 1 ~ "KRT", # remove individuals on maintenance KRT 
    egfr <= 15 ~ "<=15",
    egfr > 15 & egfr <= 30 ~ ">15-30",
    egfr > 30 & egfr <= 45 ~ ">30-45",
    egfr > 45 & egfr <= 60 ~ ">45-60",
    egfr > 60 & egfr <= 90 ~ ">60-90",
    egfr > 90 & egfr <= 105 ~ ">90-105",
    egfr > 105 ~ ">105"), 
    levels = c(
      #"KRT", # remove individuals on maintenance KRT 
      "<=15", ">15-30", ">30-45", ">45-60", ">60-90", ">90-105", ">105"))) %>% 
  mutate(#egfr = if_else(krt_pre == 1 & renal_tx == 0, 1, egfr), # remove individuals on maintenance KRT 
         egfr10 = -egfr/10) %>%
  arrange(MIX_DT) %>% 
  group_by(ALF_PE) %>% 
  mutate(n_cancers = row_number() - 1) %>% 
  ungroup() %>% 
  select(ALF_PE, MIX_DT, SITE_CD, 
         age, age10, sex, sex_name, n_creat, scr, egfr, egfr10, ckd_cat, mm_count, n_cancers, 
         sact, rtx, surg, surv1y, surv5y, any_tx, tx_smry, cancer_type, t_stage, n_stage, m_stage,
         DT_TREAT_START_DT_time, RTSTART_DT_time, SURGDATE_OPER_time) 

## check for implausible cancers or small numbers
cohort6 %>% 
  group_by(sex_name, cancer_type) %>% 
  summarise(n = n())

## remove females with "prostate cancer" (can't explain) and males with breast cancer (due to small numbers) 
cohort6 <- cohort6 %>% 
  mutate(tofilter = if_else(sex == 0 & cancer_type == "Prostate" | sex == 1 & cancer_type == "Breast", 1, 0)) %>% 
  filter(tofilter == 0)

cohort6 %>% 
  group_by(sex_name, cancer_type) %>%  
  summarise(n = sum(any_tx == "Treatment"), 
            pc = n/n()*100)

testid <- cohort6 %>% 
  filter(any_tx == "None") %>% 
  select(ALF_PE, MIX_DT, SITE_CD)
testid %>%  head()

saveRDS(cohort6, "Scratch_data/cohort6.Rds")
