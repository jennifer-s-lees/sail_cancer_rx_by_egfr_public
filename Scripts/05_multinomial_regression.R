source("Scripts/00_packages_functions.R")

cohort6 <- readRDS("Scratch_data/cohort6.Rds")
cohort6$cancer_type <- factor(cohort6$cancer_type, levels = c("Breast", "Prostate", "Colorectal", "Lung and bronchus", "Melanoma", "Renal/urinary tract"))

cohort6_m <- cohort6 %>% filter(sex_name == "Male")
cohort6_m$cancer_type <- factor(cohort6_m$cancer_type, levels = c("Prostate", "Colorectal", "Lung and bronchus", "Melanoma", "Renal/urinary tract"))
cohort6_f <- cohort6 %>% filter(sex_name == "Female")
cohort6_f$cancer_type <- factor(cohort6_f$cancer_type, levels = c("Breast", "Colorectal", "Lung and bronchus", "Melanoma", "Renal/urinary tract"))

# Assess treatment type by eGFR, adjusted for other factors ----

cohort6 <- cohort6 %>% 
  mutate(tx_smry = FctTreatmentSmry(tx_smry),
         any_tx = FctAnyTreatment(any_tx))

## Multinomial models ----

#### Decide whether to use eGFR as log transformed, quadratic term or spline ----

mod_notrans = multinom(tx_smry ~ egfr * age + sex + cancer_type + mm_count + n_cancers + m_stage, data = cohort6)
mod_log = multinom(tx_smry ~ log(egfr) * age + sex + cancer_type + mm_count + n_cancers + m_stage, data = cohort6)
mod_quad = multinom(tx_smry ~ I(egfr^2) * age + sex + cancer_type + mm_count + n_cancers + m_stage, data = cohort6)
mod_spline = multinom(tx_smry ~ ns(egfr, df = 5) * age + sex + cancer_type + mm_count + n_cancers + m_stage, data = cohort6)

# the lowest AIC is for mod_spline --> spline term best captures relationship with egfr
# this also allows creation of nice plots to visualise the relationship between egfr and treatment type
AIC(mod_notrans, mod_log, mod_quad, mod_spline)

# likelihood ratio test - including spline terms significantly improves model
anova(mod_notrans, mod_spline, test = "Chisq") # p<0.001 

mod_multi <- multinom(tx_smry ~ ns(egfr, df = 5) + sex_name + age + mm_count + n_cancers
                      + cancer_type 
                      + m_stage
                      , # models do not converge when adjusted also for t and n stages
                      data = cohort6) 

mod_multi_int <- multinom(tx_smry ~ ns(egfr, df = 5) * sex_name + age + mm_count + n_cancers
                      + cancer_type 
                      + m_stage
                      , # models do not converge when adjusted also for t and n stages
                      data = cohort6) 

AIC(mod_multi, mod_multi_int) # interaction model has lower AIC
multill1 <- logLik(mod_multi)
multill2 <- logLik(mod_multi_int)
multiLRstat <- 2 * (as.numeric(multill2) - as.numeric(multill1))
multidf <- attr(multill2, "df") - attr(multill1, "df")
multipval <- pchisq(multiLRstat, multidf, lower.tail = FALSE) 
multipval # p<0.00001


# check and interpret sex interaction for treatment type
mod_multi_int %>% 
  tidy() %>% 
  mutate(or = exp(estimate), 
         lower = exp(estimate - 1.96*std.error), 
         upper = exp(estimate + 1.96*std.error)) %>% 
  filter(grepl("egfr|:sex", term)) %>% 
  mutate(sigp = if_else(p.value <0.01, 1, 0)) %>% 
  view()

# check and interpret sex interaction for any treatment
mod_multi_int2 <- multinom(any_tx ~ ns(egfr, df = 5) * sex_name + age + mm_count + n_cancers
                          + cancer_type 
                          + m_stage
                          , # models do not converge when adjusted also for t and n stages
                          data = cohort6) 

mod_multi_int2 %>% 
  tidy() %>% 
  mutate(or = exp(estimate), 
         lower = exp(estimate - 1.96*std.error), 
         upper = exp(estimate + 1.96*std.error)) %>% 
  filter(grepl("egfr|:sex", term)) %>% 
  mutate(sigp = if_else(p.value <0.05, 1, 0)) %>% 
  view()


# there are significant sex*egfr interactions for range of treatment choices
# remainder of analyses therefore stratified by sex

#### male ----

# test value of including age interaction ----

mod_male <- multinom(tx_smry ~ egfr10 + age10 + mm_count + n_cancers
                      + cancer_type 
                      + m_stage, # models do not converge when adjusted also for t and n stages
                      data = cohort6 %>% filter(sex_name == "Male"))

mod_male_int <- multinom(tx_smry ~ egfr10 * age10 + mm_count + n_cancers
                     + cancer_type 
                     + m_stage, # models do not converge when adjusted also for t and n stages
                     data = cohort6 %>% filter(sex_name == "Male"))

AIC(mod_male, mod_male_int)
malell1 <- logLik(mod_male)
malell2 <- logLik(mod_male_int)
maleLRstat <- 2 * (as.numeric(malell2) - as.numeric(malell1))
maledf <- attr(malell2, "df") - attr(malell1, "df")
malepval <- pchisq(maleLRstat, maledf, lower.tail = FALSE) 
malepval # p<0.00001 - age interaction included in presented models

#### female ----

# test value of including age interaction ----

mod_female <- multinom(tx_smry ~ egfr10 + age10 + mm_count + n_cancers
                     + cancer_type 
                     + m_stage, # models do not converge when adjusted also for t and n stages
                     data = cohort6)
mod_female %>% tidy()

mod_female_int <- multinom(tx_smry ~ egfr10 * age10 + mm_count + n_cancers
                         + cancer_type 
                         + m_stage, # models do not converge when adjusted also for t and n stages
                         data = cohort6) 
mod_female_int %>% tidy()

AIC(mod_female, mod_female_int)
femalell1 <- logLik(mod_female)
femalell2 <- logLik(mod_female_int)
femaleLRstat <- 2 * (as.numeric(femalell2) - as.numeric(femalell1))
femaledf <- attr(femalell2, "df") - attr(femalell1, "df")
femalepval <- pchisq(femaleLRstat, femaledf, lower.tail = FALSE) 
femalepval # p<0.00001

## create multinomial models using spline terms for egfr and including age interaction term ----

## combined cohort to compare sexes ----

adjmod_spline = multinom(
  tx_smry ~ ns(egfr, df = 5) * sex_name * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6)

# this version of code uses marginaleffects "prediction" function to allow generation of splines with prediction intervals
pred_long <- predictions(
  adjmod_spline, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    sex_name = c("Male", "Female"), # when plotting splines, will therefore compare females against males
    age = mean(cohort6$age, na.rm = TRUE),
    cancer_type = levels(cohort6$cancer_type), 
    mm_count = 0, 
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_long, file = "Outputs/pred_long.csv")

## male ----

adjmod_spline_m = multinom(
  tx_smry ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m)

# this version of code uses marginaleffects "prediction" function to allow generation of splines with prediction intervals
pred_long_m <- predictions(
  adjmod_spline_m, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_m$age, na.rm = TRUE),
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0, 
    n_cancers = 0,
    m_stage = "0"
    ),
  type = "probs"
)

write_csv(pred_long_m, file = "Outputs/pred_long_m.csv")

## female ----

adjmod_spline_f = multinom(
  tx_smry ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_f)

pred_long_f <- predictions(
  adjmod_spline_f, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_f$age, na.rm = TRUE),
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_long_f, file = "Outputs/pred_long_f.csv")

## for both male and female models, check for age*egfr spline term interaction with treatment type ----
# do this with age tested for each additional 10 years of aging

int_adjmod_spline_m = multinom(
  tx_smry ~ ns(egfr, df = 5) * age10 + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m)%>% 
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%   
  filter(grepl(":age", term)) 

int_adjmod_spline_f = multinom(
  tx_smry ~ ns(egfr, df = 5) * age10 + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_f) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%   
  filter(grepl(":age", term)) 

txsmry_age_egfr_interaction <- 
  bind_rows(int_adjmod_spline_m %>% mutate(sex = "Male"), int_adjmod_spline_f %>% mutate(sex = "Female")) %>% 
  mutate(OR = round(estimate, 2)) %>% 
  mutate(`95% Confidence Interval` = paste0(round(conf.low, 2), " - ", round(conf.high, 2)),
         p.value = round(p.value, 3)) %>% 
  rename(`Treatment type` = y.level) %>% 
  select(sex, `Treatment type`, term, OR, `95% Confidence Interval`, p.value)
txsmry_age_egfr_interaction %>%  view()

write.csv(txsmry_age_egfr_interaction, "Outputs/txsmry_age_egfr_interaction.csv") 

## create logistic regression models using spline terms for egfr ----

## create combined cohort to compare female vs. male 

logmod_spline = multinom(
  any_tx ~ ns(egfr, df = 5) * sex_name * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6)

pred_logmod <- predictions(
  logmod_spline, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6$age, na.rm = TRUE),
    sex_name = c("Male", "Female"), # compare females to males
    cancer_type = levels(cohort6$cancer_type), 
    mm_count = 0,
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod, file = "Outputs/pred_logmod.csv")

## run in males

logmod_spline_m = multinom(
  any_tx ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m)

pred_logmod_m <- predictions(
  logmod_spline_m, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_m$age, na.rm = TRUE),
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0,
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_m, file = "Outputs/pred_logmod_m.csv")

logmod_spline_f = multinom(
  any_tx ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_f)

pred_logmod_f <- predictions(
  logmod_spline_f, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_f$age, na.rm = TRUE),
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_f, file = "Outputs/pred_logmod_f.csv")

## for both male and female models, check for age *egfr spline term interaction with no treatment ----
# do this with age tested for each additional 10 years of aging

int_logmod_spline_m = multinom(
  any_tx ~ ns(egfr, df = 5) * age10 + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%   
  filter(grepl(":age", term)) 

int_logmod_spline_f = multinom(
  any_tx ~ ns(egfr, df = 5) * age10 + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m) %>% 
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%   
  filter(grepl(":age", term)) 

anytx_age_egfr_interaction <- 
  bind_rows(int_logmod_spline_m %>% mutate(sex = "Male"), int_logmod_spline_f %>% mutate(sex = "Female")) %>% 
  mutate(OR = round(estimate, 2)) %>% 
  mutate(`95% Confidence Interval` = paste0(round(conf.low, 2), " - ", round(conf.high, 2)),
         p.value = round(p.value, 3)) %>% 
  select(sex, term, OR, `95% Confidence Interval`, p.value)
anytx_age_egfr_interaction %>%  view()

write.csv(anytx_age_egfr_interaction, "Outputs/anytx_age_egfr_interaction.csv") 

# Run models using spline terms for 1 year survival ----

## create combined cohort to compare female vs. male ---- 

logmod_spline_surv1y = multinom(
  surv1y ~ ns(egfr, df = 5) * sex_name * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6)

pred_logmod_surv1y <- predictions(
  logmod_spline_surv1y, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6$age, na.rm = TRUE),
    sex_name = c("Male", "Female"), # compare females to males
    cancer_type = levels(cohort6$cancer_type), 
    mm_count = 0,
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_surv1y, file = "Outputs/pred_logmod_surv1y.csv")

## run in males ----

logmod_spline_m_surv1y = multinom(
  surv1y ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m)

pred_logmod_m_surv1y <- predictions(
  logmod_spline_m_surv1y, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_m$age, na.rm = TRUE),
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0,
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_m_surv1y, file = "Outputs/pred_logmod_m_surv1y.csv")

# run in females ----

logmod_spline_f_surv1y = multinom(
  surv1y ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_f)

pred_logmod_f_surv1y <- predictions(
  logmod_spline_f_surv1y, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_f$age, na.rm = TRUE),
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_f_surv1y, file = "Outputs/pred_logmod_f_surv1y.csv")

# Run models using spline terms for 5 year survival ----

## create combined cohort to compare female vs. male ---- 

logmod_spline_surv5y = multinom(
  surv5y ~ ns(egfr, df = 5) * sex_name * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6)

pred_logmod_surv5y <- predictions(
  logmod_spline_surv5y, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6$age, na.rm = TRUE),
    sex_name = c("Male", "Female"), # compare females to males
    cancer_type = levels(cohort6$cancer_type), 
    mm_count = 0, 
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_surv5y, file = "Outputs/pred_logmod_surv5y.csv")

## run in males ----

logmod_spline_m_surv5y = multinom(
  surv5y ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m)

pred_logmod_m_surv5y <- predictions(
  logmod_spline_m_surv5y, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_m$age, na.rm = TRUE),
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0,
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_m_surv5y, file = "Outputs/pred_logmod_m_surv5y.csv")

# run in females ----

logmod_spline_f_surv5y = multinom(
  surv5y ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_f)

pred_logmod_f_surv5y <- predictions(
  logmod_spline_f_surv5y, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = mean(cohort6_f$age, na.rm = TRUE),
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    n_cancers = 0,
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_f_surv5y, file = "Outputs/pred_logmod_f_surv5y.csv")
