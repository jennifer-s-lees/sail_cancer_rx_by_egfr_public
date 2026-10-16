source("Scripts/00_packages_functions.R")

login <- getLogin('')
channel <- odbcConnect('PR_SAIL',login[1],login[2])

## identify kidney failure ----

krt_icd = sqlQuery(channel, 
                   "SELECT ALF_PE, ADMIS_DT, D.DIAG_CD_1234 AS EVENT_CODE FROM SAIL1214V.PEDW_SPELL_20240909 S
INNER JOIN SAIL1214V.PEDW_EPISODE_20240909 E
ON S.PROV_UNIT_CD = E.PROV_UNIT_CD
AND S.SPELL_NUM_PE = E.SPELL_NUM_PE
left JOIN SAIL1214V.PEDW_DIAG_20240909 D
ON S.PROV_UNIT_CD = D.PROV_UNIT_CD
AND S.SPELL_NUM_PE = D.SPELL_NUM_PE
and e.epi_num = d.epi_num
where d.diag_cd_1234 IN ('E853', 'N165', 'T824', 'T861', 'Y602', 'Y612', 'Y622', 'Y841', 'Z490', 'Z491', 'Z492', 'Z940', 'Z992')"
)

krt_maintenance_icd = sqlQuery(channel, 
                               "SELECT ALF_PE, ADMIS_DT, D.DIAG_CD_1234 AS EVENT_CODE FROM SAIL1214V.PEDW_SPELL_20240909 S
INNER JOIN SAIL1214V.PEDW_EPISODE_20240909 E
ON S.PROV_UNIT_CD = E.PROV_UNIT_CD
AND S.SPELL_NUM_PE = E.SPELL_NUM_PE
left JOIN SAIL1214V.PEDW_DIAG_20240909 D
ON S.PROV_UNIT_CD = D.PROV_UNIT_CD
AND S.SPELL_NUM_PE = D.SPELL_NUM_PE
and e.epi_num = d.epi_num
where d.diag_cd_1234 IN ('N165', 'T861', 'Z492', 'Z940')"
)

ckd5_icd = sqlQuery(channel, 
                    "SELECT ALF_PE, ADMIS_DT, D.DIAG_CD_1234 AS EVENT_CODE FROM SAIL1214V.PEDW_SPELL_20240909 S
INNER JOIN SAIL1214V.PEDW_EPISODE_20240909 E
ON S.PROV_UNIT_CD = E.PROV_UNIT_CD
AND S.SPELL_NUM_PE = E.SPELL_NUM_PE
left JOIN SAIL1214V.PEDW_DIAG_20240909 D
ON S.PROV_UNIT_CD = D.PROV_UNIT_CD
AND S.SPELL_NUM_PE = D.SPELL_NUM_PE
and e.epi_num = d.epi_num
where d.diag_cd_1234 IN ('N165', 'T861', 'Z940', 'E853', 'N180', 'N185', 'Q601')"
)

krt_opcs = sqlQuery(channel, 
                    "SELECT ALF_PE, ADMIS_DT, O.OPER_CD AS EVENT_CODE FROM SAIL1214V.PEDW_SPELL_20240909 S
INNER JOIN SAIL1214V.PEDW_EPISODE_20240909 E
ON S.PROV_UNIT_CD = E.PROV_UNIT_CD
AND S.SPELL_NUM_PE = E.SPELL_NUM_PE
left JOIN SAIL1214V.PEDW_OPER_20240909 O
ON S.PROV_UNIT_CD = O.PROV_UNIT_CD
AND S.SPELL_NUM_PE = O.SPELL_NUM_PE
and e.epi_num = O.epi_num
where O.OPER_CD IN ('M012', 'M013', 'M014', 'M015', 'M018', 'M019', 'M084', 'M174', 'M178', 'M179', 'X401', 'X402', 'X403', 'X404', 'X405', 'X406', 'X407', 'X408', 'X409', 'X411', 'X412', 'X418', 'X419', 'X421', 'X428', 'X429', 'X431')"
)

krt_maintenance_opcs = sqlQuery(channel, 
                                "SELECT ALF_PE, ADMIS_DT, O.OPER_CD AS EVENT_CODE FROM SAIL1214V.PEDW_SPELL_20240909 S
INNER JOIN SAIL1214V.PEDW_EPISODE_20240909 E
ON S.PROV_UNIT_CD = E.PROV_UNIT_CD
AND S.SPELL_NUM_PE = E.SPELL_NUM_PE
left JOIN SAIL1214V.PEDW_OPER_20240909 O
ON S.PROV_UNIT_CD = O.PROV_UNIT_CD
AND S.SPELL_NUM_PE = O.SPELL_NUM_PE
and e.epi_num = O.epi_num
where O.OPER_CD IN ('M012', 'M013', 'M014', 'M015', 'M018', 'M019', 'M084', 'M174', 'M178', 'M179', 'X402', 'X405', 'X406', 'X411', 'X412')"
)

ckd5_opcs = sqlQuery(channel, 
                     "SELECT ALF_PE, ADMIS_DT, O.OPER_CD AS EVENT_CODE FROM SAIL1214V.PEDW_SPELL_20240909 S
INNER JOIN SAIL1214V.PEDW_EPISODE_20240909 E
ON S.PROV_UNIT_CD = E.PROV_UNIT_CD
AND S.SPELL_NUM_PE = E.SPELL_NUM_PE
left JOIN SAIL1214V.PEDW_OPER_20240909 O
ON S.PROV_UNIT_CD = O.PROV_UNIT_CD
AND S.SPELL_NUM_PE = O.SPELL_NUM_PE
and e.epi_num = O.epi_num
where O.OPER_CD IN ('M012', 'M013', 'M014', 'M015', 'M018', 'M019', 'M084', 'M174', 'M178', 'M179', 'X402', 'X405', 'X406', 'X411', 'X412', 'M172', 'M023', 'L741', 'L742', 'L743', 'L744', 'L745', 'L746', 'L748', 'L749')"
)

krt = bind_rows(krt_icd, krt_opcs) %>%
  distinct(ALF_PE, ADMIS_DT, .keep_all = T) %>%
  rename(krt_date = ADMIS_DT) %>%
  rename(krt_code = EVENT_CODE) %>% 
  filter(ALF_PE %in% cohort5_alfpe$ALF_PE)
krt_maintenance = bind_rows(krt_maintenance_icd, krt_maintenance_opcs) %>%
  distinct(ALF_PE, ADMIS_DT, .keep_all = T) %>%
  rename(krt_date = ADMIS_DT) %>%
  rename(krt_maintenance_code = EVENT_CODE) %>% 
  filter(ALF_PE %in% cohort5_alfpe$ALF_PE)
ckd5 = bind_rows(ckd5_icd, ckd5_opcs) %>%
  distinct(ALF_PE, ADMIS_DT, .keep_all = T) %>%
  rename(ckd5_date = ADMIS_DT) %>%
  rename(ckd5_code = EVENT_CODE) %>% 
  filter(ALF_PE %in% cohort5_alfpe$ALF_PE)

saveRDS(krt, "Processed_data/krt.Rds")
saveRDS(krt_maintenance, "Processed_data/krt_maintenance.Rds")
saveRDS(ckd5, "Processed_data/ckd5.Rds")

## read in datasets ----

krt <- readRDS("Processed_data/krt.Rds")
krt_maintenance <- readRDS("Processed_data/krt_maintenance.Rds")
ckd5 <- readRDS("Processed_data/ckd5.Rds")

kf_on_krt = krt %>%
  inner_join(ckd5, by="ALF_PE", relationship = "many-to-many") %>% 
  mutate(gap = as.numeric(ckd5_date-krt_date)) %>%
  filter(ckd5_date<krt_date | between(gap, -366, 366)) %>%
  group_by(ALF_PE) %>%
  filter(krt_date == min(krt_date)) %>%
  ungroup() %>%
  distinct(ALF_PE, .keep_all = T) %>%
  bind_rows(krt_maintenance) %>%
  group_by(ALF_PE) %>%
  filter(krt_date == min(krt_date)) %>%
  ungroup() %>%
  distinct(ALF_PE, .keep_all = T)

saveRDS(kf_on_krt, "Processed_data/kf_on_krt.Rds")
