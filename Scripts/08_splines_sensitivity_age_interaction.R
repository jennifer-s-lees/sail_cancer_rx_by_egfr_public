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

## Multinomial models holding age constant at younger and older ages ----
young <- 45
old  <- 85

# Treatment summary ----
## male ----

adjmod_spline_m = multinom(
  tx_smry ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m)
cohort6_m$age %>%  hist()

# this version of code uses marginaleffects "prediction" function to allow generation of splines with prediction intervals
pred_long_m_young <- predictions(
  adjmod_spline_m, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = young,
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
    ),
  type = "probs"
)

pred_long_m_old <- predictions(
  adjmod_spline_m, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = old,
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_long_m_young, file = "Outputs/pred_long_m_young.csv")
write_csv(pred_long_m_old, file = "Outputs/pred_long_m_old.csv")

## female ----

adjmod_spline_f = multinom(
  tx_smry ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_f)
cohort6_f$age %>%  hist()

# this version of code uses marginaleffects "prediction" function to allow generation of splines with prediction intervals
pred_long_f_young <- predictions(
  adjmod_spline_f, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = young,
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
  ),
  type = "probs"
)

pred_long_f_old <- predictions(
  adjmod_spline_f, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = old,
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_long_f_young, file = "Outputs/pred_long_f_young.csv")
write_csv(pred_long_f_old, file = "Outputs/pred_long_f_old.csv")

# now draw the splines ----
## young ----

pred_long_m_young <- read_csv("Outputs/pred_long_m_young.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))
pred_long_m_old <- read_csv("Outputs/pred_long_m_old.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))
pred_long_f_young <- read_csv("Outputs/pred_long_f_young.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))
pred_long_f_old <- read_csv("Outputs/pred_long_f_old.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))

# male
male_plot_young <- ggplot(pred_long_m_young, aes(x = egfr, y = estimate, color = tx_smry)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tx_smry), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_brewer(palette = "Spectral") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male - age held at 45 years", color = "", fill = "") 

# female
female_plot_young <- ggplot(pred_long_f_young, aes(x = egfr, y = estimate, color = tx_smry)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tx_smry), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_brewer(palette = "Spectral") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female - age held at 45 years", color = "", fill = "") 

comb_plot_young <- ggarrange(male_plot_young, female_plot_young,
                       nrow = 2, 
                       common.legend = TRUE)
comb_plot_young 

ggsave(comb_plot_young, 
       file = "Outputs/spline_plot_tx_smry_by_egfr_age_45.png",
       height = 8, width = 16)

## old ----

# male
male_plot_old <- ggplot(pred_long_m_old, aes(x = egfr, y = estimate, color = tx_smry)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tx_smry), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_brewer(palette = "Spectral") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male - age held at 85 years", color = "", fill = "") 

# female
female_plot_old <- ggplot(pred_long_f_old, aes(x = egfr, y = estimate, color = tx_smry)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tx_smry), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_brewer(palette = "Spectral") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female - age held at 85 years", color = "", fill = "") 

comb_plot_old <- ggarrange(male_plot_old, female_plot_old,
                             nrow = 2, 
                             common.legend = TRUE)
comb_plot_old 

ggsave(comb_plot_old, 
       file = "Outputs/spline_plot_tx_smry_by_egfr_age_85.png",
       height = 8, width = 16)




# Any treatment  ----
## male ----

logmod_spline_m = multinom(
  any_tx ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_m)
cohort6_m$age %>%  hist()

# this version of code uses marginaleffects "prediction" function to allow generation of splines with prediction intervals
pred_logmod_spline_m_young <- predictions(
  logmod_spline_m, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = young,
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
  ),
  type = "probs"
)

pred_logmod_spline_m_old <- predictions(
  logmod_spline_m, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = old,
    cancer_type = levels(cohort6_m$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_spline_m_young, file = "Outputs/pred_logmod_spline_m_young.csv")
write_csv(pred_logmod_spline_m_old, file = "Outputs/pred_logmod_spline_m_old.csv")

## female ----

logmod_spline_f = multinom(
  any_tx ~ ns(egfr, df = 5) * age + cancer_type + mm_count + n_cancers + m_stage, 
  data = cohort6_f)
cohort6_f$age %>%  hist()

# this version of code uses marginaleffects "prediction" function to allow generation of splines with prediction intervals
pred_logmod_spline_f_young <- predictions(
  logmod_spline_f, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = young,
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
  ),
  type = "probs"
)

pred_logmod_spline_f_old <- predictions(
  logmod_spline_f, 
  newdata = datagrid(
    egfr = seq(0, 150, length = 200),
    age = old,
    cancer_type = levels(cohort6_f$cancer_type), 
    mm_count = 0, 
    m_stage = "0"
  ),
  type = "probs"
)

write_csv(pred_logmod_spline_f_young, file = "Outputs/pred_logmod_spline_f_young.csv")
write_csv(pred_logmod_spline_f_old, file = "Outputs/pred_logmod_spline_f_old.csv")

# now draw the splines ----
## young ----

pred_logmod_spline_m_young <- read_csv("Outputs/pred_logmod_spline_m_young.csv") %>% 
  mutate(any_tx = FctAnyTreatment(group),
         cancer_type = FctCancerType(cancer_type))
pred_logmod_spline_m_old <- read_csv("Outputs/pred_logmod_spline_m_old.csv") %>% 
  mutate(any_tx = FctAnyTreatment(group),
         cancer_type = FctCancerType(cancer_type))
pred_logmod_spline_f_young <- read_csv("Outputs/pred_logmod_spline_f_young.csv") %>% 
  mutate(any_tx = FctAnyTreatment(group),
         cancer_type = FctCancerType(cancer_type))
pred_logmod_spline_f_old <- read_csv("Outputs/pred_logmod_spline_f_old.csv") %>% 
  mutate(any_tx = FctAnyTreatment(group),
         cancer_type = FctCancerType(cancer_type))

# male
male_plot_logmod_young <- ggplot(pred_logmod_spline_m_young, aes(x = egfr, y = estimate, color = any_tx)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = any_tx), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = "royalblue") +
  scale_fill_manual(values = "royalblue") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  guides(colour = "none", fill = "none") +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male - age held at 45 years", color = "", fill = "") 

# female
female_plot_logmod_young <- ggplot(pred_logmod_spline_f_young, aes(x = egfr, y = estimate, color = any_tx)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = any_tx), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = "hotpink") +
  scale_fill_manual(values = "hotpink") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  guides(colour = "none", fill = "none") +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female - age held at 45 years", color = "", fill = "") 

comb_plot_logmod_young <- ggarrange(male_plot_logmod_young, female_plot_logmod_young,
                             nrow = 2, 
                             common.legend = TRUE)
comb_plot_logmod_young 

ggsave(comb_plot_logmod_young, 
       file = "Outputs/spline_plot_any_tx_by_egfr_age_45.png",
       height = 8, width = 16)

## old ----

# male
male_plot_logmod_old <- ggplot(pred_logmod_spline_m_old, aes(x = egfr, y = estimate, color = any_tx)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = any_tx), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = "royalblue") +
  scale_fill_manual(values = "royalblue") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  guides(colour = "none", fill = "none") +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male - age held at 85 years", color = "", fill = "") 

# female
female_plot_logmod_old <- ggplot(pred_logmod_spline_f_old, aes(x = egfr, y = estimate, color = any_tx)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = any_tx), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30, 60), linetype = "dashed") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = "hotpink") +
  scale_fill_manual(values = "hotpink") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  guides(colour = "none", fill = "none") +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female - age held at 85 years", color = "", fill = "") 

comb_plot_logmod_old <- ggarrange(male_plot_logmod_old, female_plot_logmod_old,
                                    nrow = 2, 
                                    common.legend = TRUE)
comb_plot_logmod_old 

ggsave(comb_plot_logmod_old, 
       file = "Outputs/spline_plot_logmod_any_tx_by_egfr_age_85.png",
       height = 8, width = 16)

# combine male groups and compare age

comb_plot_logmod_m <- ggarrange(male_plot_logmod_young, male_plot_logmod_old,
                                  nrow = 2, 
                                  common.legend = TRUE)
comb_plot_logmod_m

ggsave(comb_plot_logmod_m, 
       file = "Outputs/spline_plot_logmod_any_tx_by_egfr_m.png",
       height = 8, width = 16)

# combine female groups and compare age

# combine female groups and compare age

comb_plot_logmod_f <- ggarrange(female_plot_logmod_young, female_plot_logmod_old,
                                nrow = 2, 
                                common.legend = TRUE)
comb_plot_logmod_f

ggsave(comb_plot_logmod_f, 
       file = "Outputs/spline_plot_logmod_any_tx_by_egfr_f.png",
       height = 8, width = 16)


