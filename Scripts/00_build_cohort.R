source("Scripts/00_packages_functions.R")

login <- getLogin('')
channel <- odbcConnect('PR_SAIL',login[1],login[2])

## Available SAIL tables
xmn <- sqlTables(channel)
xmn_us <- xmn[stringr::str_detect(xmn$TABLE_SCHEM, "1214"), ]
rm(xmn)
write.csv(xmn_us, "\\\\dfssail/UKSeRPSAIL/shared/1214 - The Relationship Between Cancer and Kidney Disease/SQL_tables_2025.csv")
xmn_us

#### Inspect tables ----

xmn_us <- read.csv("\\\\dfssail/UKSeRPSAIL/shared/1214 - The Relationship Between Cancer and Kidney Disease/SQL_tables_2025.csv")
cnis <- xmn_us %>% 
  filter(grepl("CNIS", TABLE_NAME))
cnis %>%  distinct()

#X TABLE_CAT TABLE_SCHEM                         TABLE_NAME TABLE_TYPE REMARKS
#1  231827        NA   SAIL1214V        CNIS_COMPLICATIONS_20240908       VIEW      NA
#2  231828        NA   SAIL1214V            CNIS_DIAGNOSIS_20240908       VIEW      NA
#3  231829        NA   SAIL1214V    CNIS_DISEASERECURRENCE_20240908       VIEW      NA
#4  231830        NA   SAIL1214V          CNIS_DRUGTHERAPY_20240908       VIEW      NA
#5  231831        NA   SAIL1214V       CNIS_INVESTIGATIONS_20240908       VIEW      NA
#6  231832        NA   SAIL1214V            CNIS_KEYWORKER_20240908       VIEW      NA
#7  231833        NA   SAIL1214V                  CNIS_MDT_20240908       VIEW      NA
#8  231834        NA   SAIL1214V         CNIS_OBSERVATIONS_20240908       VIEW      NA
#9  231835        NA   SAIL1214V            CNIS_PATHOLOGY_20240908       VIEW      NA
#10 231836        NA   SAIL1214V              CNIS_PATIENT_20240908       VIEW      NA
#11 231837        NA   SAIL1214V         CNIS_RADIOTHERAPY_20240908       VIEW      NA
#12 231838        NA   SAIL1214V CNIS_RADIOTHERAPY_ENDDATE_20240908       VIEW      NA
#13 231839        NA   SAIL1214V             CNIS_REFERRAL_20240908       VIEW      NA
#14 231840        NA   SAIL1214V              CNIS_SURGERY_20240908       VIEW      NA

## cohort should be limited to those:
# 1. with cancer diagnosis from 2006-01-01 (when creatinine values were standardised)
# 2. with cancer diagnosis before 2023-09-01 (to allow 1 year of survival follow up)
# 3. with creatinine value at time of cancer diagnosis (within 1 year before or 3 months after)

## cancer tables ----
# Replace * with variable names, separated by comma as required

sqlQuery(channel, 
        "SELECT * 
         FROM SAIL1214V.CNIS_DIAGNOSIS_20240908
         LIMIT 10")

cancers <- sqlQuery(channel, 
         "SELECT * 
        FROM SAIL1214V.CNIS_DIAGNOSIS_20240908;")

# now keep cancer data for diagnoses in correct time frame
# limit to breast, prostate, colorectal, lung, melanoma, renal tract cancers
# filter to first diagnosis at each site for each individual

cancers_2006 <- cancers %>% 
  mutate(SITE_CD = str_extract(SITE_CD, "^[A-Za-z]+\\d{2}")) %>% 
  arrange(MIX_DT) %>% 
  group_by(ALF_PE, SITE_CD) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  filter(MIX_DT >= as.Date("2006-01-01") & MIX_DT <= as.Date("2023-09-01")) %>% 
  filter(grepl("C50|C61|C18|C19|C20|C33|C34|C43|C64|C65|C66|C67|C68", SITE_CD)) 

saveRDS(cancers_2006_test, "Scratch_data/cancers_2006.Rds")

# list of individual patient IDs for those with cancer record in CNIS
ca_pt <- cancers_2006 %>% 
  select(ALF_PE) %>% 
  distinct()

length(unique(ca_pt$ALF_PE)) # 197,196 unique individuals

## death ----

death <- sqlQuery(channel, 
                  "SELECT ALF_PE, DEATH_DT, DEATHCAUSE_DIAG_UNDERLYING_CD
         FROM SAIL1214V.ADDE_DEATHS_20240901;")
death <- death %>% 
  arrange(DEATH_DT) %>% 
  group_by(ALF_PE) %>% 
  slice_head(n=1) %>% 
  ungroup()

cancers_death <- cancers_2006 %>% 
  left_join(death, by = c("ALF_PE")) %>% 
  mutate(DEATH_DT = if_else(is.na(DEATH_DT), as.Date("2024-09-03"), DEATH_DT)) %>% 
  filter(DEATH_DT >= MIX_DT)  # remove individuals who had a death date before cancer diagnosis date 

saveRDS(cancers_death, "Scratch_data/cancers_death.Rds")

## creatinine values ---- 

# extract all from GP records
gp_scr = sqlQuery(channel, 
                    "SELECT P.ALF_PE, P.WOB, P.GNDR_CD, EVENT_VAL AS CREAT, EVENT_DT, EVENT_CD
        FROM SAIL1214V.WLGP_GP_EVENT_CLEANSED_20240401 E
        LEFT JOIN SAIL1214V.WLGP_PATIENT_ALF_CLEANSED_20240401 P
        ON E.PRAC_CD_PE = P.PRAC_CD_PE
        AND E.LOCAL_NUM_PE = P.LOCAL_NUM_PE
        WHERE EVENT_CD IN ('44J33', '44J30', '44J3z', '44J3.', 'EMISREQ|44J3', '44JF.', '44J32', '44JC.', '44JD.', '44J31', '4Q40.', 
        '9iB..', '9iC..', '9iD2.', '9iD..', '9iD3.', '9iD0.', '9iD1.', '9iD4.', '9S2..', '9S3..', '9S47.', '9S42.', '9S48.', '9S4..', '9S46.', '9S44.', '9S45.', '9S41.', '9S43','9S5..', '9SA3.', '9SA5.', '9SG..');")

# then limit only to creatinine values in patients with a diagnosis of cancer in CNIS
## limit to individuals with cancer and creatinine values >= 2006 
# if multiple creatinine measures were taken on the same day, calculate mean result
# creatinine values must be within 1 year before or up to 90 days after cancer diagnosis

creat <- gp_scr %>%
  filter(ALF_PE %in% ca_pt$ALF_PE) %>% 
  filter(EVENT_DT >= as.Date("2005-01-01")) %>% # filter creatinine on or after 01/01/2005 (to allow values for those with cancer in 2006)
  filter(!is.na(CREAT) & CREAT >= 20 & CREAT <= 1000) %>% # filter implausible creatinine values
  group_by(ALF_PE, EVENT_DT) %>% 
  mutate(CREAT = mean(CREAT)) %>% # calculate mean creatinine if 2 or more values taken on same day
  slice_head(n = 1) %>% # only select 1 record from each date 
  ungroup()

saveRDS(creat, "Scratch_data/creat.Rds")

## join creatinine values for these patients ----
# filter only to individuals >= 18 years
# filter to those with creatinine values within 1 year before or 3 months after cancer diagnosis

cohort1 <- cancers_death %>% 
  left_join(creat %>% select(ALF_PE, EVENT_DT, CREAT, WOB, GNDR_CD), "ALF_PE", relationship = "many-to-many") %>% 
  rename(scr_date = EVENT_DT, 
         scr = CREAT,
         dob = WOB) %>% 
  mutate(sex = if_else(GNDR_CD == 1, 1, 0), # 1 = Male, 0 = Female
         age = as.numeric(MIX_DT - dob)/365.25) %>%
  filter(age >= 18) %>% 
  filter(as.numeric(scr_date - MIX_DT) >= -365 & as.numeric(scr_date - MIX_DT) <= 90) %>% 
  mutate(diff = abs(as.numeric(scr_date - MIX_DT))) %>% 
  group_by(ALF_PE, MIX_DT, SITE_CD) %>% 
  mutate(n_creat = n()) %>% 
  ungroup() %>% 
  arrange(diff) %>% 
  group_by(ALF_PE, MIX_DT, SITE_CD) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  mutate(
    scr_mgdl = scr * 0.0113,
    age_cr = as.numeric(scr_date - dob)/365.25,
    ethnicity = 0) %>% 
  mutate(egfr = CKDEpi.creat(scr_mgdl, sex, age_cr, ethnicity))
length(unique(cohort1$ALF_PE))  # 124,042 selected cancers in 129,607 individuals

cohort1 %>% count() # 129,607 records 
cohort1 %>% filter(n_creat == 1) %>% count() # 60,582 records 
# 60,582 of 136,182 eligible individuals had only 1 creatinine value in the required time scale

saveRDS(cohort1, "Scratch_data/cohort1.Rds")

## drug treatments ---- 

drug <- sqlQuery(channel, 
         "SELECT * 
         FROM SAIL1214V.CNIS_DRUGTHERAPY_20240908
         WHERE DT_TREAT_START_DT > '2005-01-01'; ")

drug1 <- drug %>% 
  filter(ALF_PE %in% cohort1$ALF_PE) %>% 
  select(ALF_PE, MIX_DT, SITE_CD, DT_TREAT_TYPE, DT_TREAT_START_DT, DT_TREAT_END_DT, DT_TREAT_INTENT) %>% 
  mutate(SITE_CD = str_extract(SITE_CD, "^[A-Za-z]+\\d{2}")) %>% 
  group_by(ALF_PE, MIX_DT, SITE_CD) %>% 
  mutate(sact = 1) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  select(ALF_PE, MIX_DT, SITE_CD, sact, DT_TREAT_START_DT, DT_TREAT_TYPE)
  
# cohort1 <- readRDS("Scratch_data/cohort1.Rds")
cohort2 <- cohort1 %>% 
  mutate(SITE_CD = str_extract(SITE_CD, "^[A-Za-z]+\\d{2}")) %>% 
  left_join(drug1, by = c("ALF_PE", "MIX_DT", "SITE_CD"))

length(unique(cohort2$ALF_PE)) # 124,042  

## radiotherapy ----

rt <- sqlQuery(channel, 
                "SELECT * 
         FROM SAIL1214V.CNIS_RADIOTHERAPY_20240908
               WHERE RTSTART_DT > '2005-01-01';")

rt1 <- rt %>%
  filter(ALF_PE %in% cohort1$ALF_PE) %>% 
  select(ALF_PE, MIX_DT, SITE_CD, RTNUMFRACTIONSPLANNED, RTSTART_DT, RT_TREAT_INTENT) %>% 
  mutate(SITE_CD = str_extract(SITE_CD, "^[A-Za-z]+\\d{2}")) %>% 
  group_by(ALF_PE, MIX_DT, SITE_CD) %>% 
  mutate(rtx = 1) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  select(ALF_PE, MIX_DT, SITE_CD, rtx, RTSTART_DT, RT_TREAT_INTENT)

cohort3 <- cohort2 %>% 
  left_join(rt1, by = c("ALF_PE", "MIX_DT", "SITE_CD"))

length(unique(cohort3$ALF_PE))  # 124,042 individuals

## surgery ----

surg <- sqlQuery(channel, 
               "SELECT * 
         FROM SAIL1214V.CNIS_SURGERY_20240908; ") # SURGDATE_OPER is a character - need to convert to date first

surg1 <- surg %>% 
  filter(ALF_PE %in% cohort1$ALF_PE) %>% 
  select(ALF_PE, MIX_DT, SITE_CD, SURGDATE_OPER, SURG_ADMIS_DT, SURG_DISCH_DT, SURGREASONINTENTION, SURGMODEOP_URGENCYSURG, SURGPOSTOPDEATH, SURGLAPAROSCOPY) %>% 
  mutate(SITE_CD = str_extract(SITE_CD, "^[A-Za-z]+\\d{2}")) %>% 
  group_by(ALF_PE, MIX_DT, SITE_CD) %>% 
  mutate(surg = 1) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  select(ALF_PE, MIX_DT, SITE_CD, surg, SURGDATE_OPER, SURGREASONINTENTION) %>% 
  mutate(SURGDATE_OPER = as.Date(SURGDATE_OPER, format = "%d/%m/%Y")) %>% 
  filter(SURGDATE_OPER > as.Date("2005-01-01"))

cohort4 <- cohort3 %>% 
  left_join(surg1, by = c("ALF_PE", "MIX_DT", "SITE_CD"))

length(unique(cohort4$ALF_PE))  # 124,042

saveRDS(cohort4, "Scratch_data/cohort4.Rds")

rm(cohort1, cohort2, cohort3, drug, drug1, rt, rt1, surg, surg1)


