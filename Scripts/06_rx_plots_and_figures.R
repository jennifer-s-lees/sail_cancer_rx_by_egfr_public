source("Scripts/00_packages_functions.R")

cohort6 <- readRDS("Scratch_data/cohort6.Rds")

## BAR PLOTS of treatment types ----

# review counts of people who received any treatment
any_tx_counts_m <- cohort6 %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>% 
  select(sex_name, ckd_cat, any_tx, cancer_type) %>% 
  mutate(any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type)) %>% 
  FilterSmallGroups(ckd_cat, any_tx, cancer_type) %>% 
  group_by(sex_name, ckd_cat, any_tx, cancer_type) %>% 
  summarise(count = n()) 
any_tx_counts_f <- cohort6 %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>% 
  select(sex_name, ckd_cat, any_tx, cancer_type) %>% 
  mutate(any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type)) %>% 
  FilterSmallGroups(ckd_cat, any_tx, cancer_type) %>% 
  group_by(sex_name, ckd_cat, any_tx, cancer_type) %>% 
  summarise(count = n()) 

any_tx_counts <- bind_rows(
  any_tx_counts_m,
  any_tx_counts_f
  )
any_tx_counts %>% filter(count < 10) %>% glimpse() # no small groups

# save counts 
write_csv(any_tx_counts, "Outputs/any_tx_counts.csv")

# by Any treatment
any_tx_plot_m <- any_tx_counts %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>% 
  ggplot(aes(x = ckd_cat, y = count, fill = any_tx)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("lightblue", "royalblue")) +
  facet_wrap(~ cancer_type, nrow = 1) +
  labs(x = "eGFR category (mL/min/1.73m2)", y = "Percentage", fill = "", title = "Male") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45))
any_tx_plot_m

any_tx_plot_f <- any_tx_counts %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>% 
  ggplot(aes(x = ckd_cat, y = count, fill = any_tx)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("pink", "hotpink")) +
  facet_wrap(~ cancer_type, nrow = 1) +
  labs(x = "eGFR category (mL/min/1.73m2)", y = "Percentage", fill = "", title = "Female") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45))
any_tx_plot_f

any_tx_plot <- ggarrange(any_tx_plot_m, any_tx_plot_f, 
                         nrow = 2, 
                         common.legend = TRUE)
any_tx_plot

ggsave("Outputs/any_tx_plot.png",
       any_tx_plot,
       height = 8, width = 16)

# by Treatment summary
tx_smry_counts_m <- cohort6 %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>% 
  FilterSmallGroups(ckd_cat, tx_smry, cancer_type) %>% 
  group_by(sex_name, ckd_cat, tx_smry, cancer_type) %>% 
  summarise(count = n(), .groups = "drop") 
tx_smry_counts_f <- cohort6 %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>% 
  FilterSmallGroups(ckd_cat, tx_smry, cancer_type) %>% 
  group_by(sex_name, ckd_cat, tx_smry, cancer_type) %>% 
  summarise(count = n(), .groups = "drop") 
tx_smry_counts <- bind_rows(tx_smry_counts_m, tx_smry_counts_f)
tx_smry_counts %>% filter(count < 10) %>% glimpse() #no groups with counts <10

# save file that include treatment counts for tx_smry_plot below
tx_smry_counts %>% write_csv("Outputs/tx_smry_counts.csv")

# males
tx_smry_plot_m <- tx_smry_counts %>%   
  mutate(cancer_type = FctCancerType(cancer_type)) %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>% 
  ggplot(aes(x = ckd_cat, y = count, fill = tx_smry)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Spectral") +
  facet_wrap(~ cancer_type, nrow = 1) +
  labs(x = "eGFR category (mL/min/1.73m2)", y = "Percentage", fill = "", title = "Male") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45))
tx_smry_plot_m

tx_smry_plot_f <- tx_smry_counts %>% 
  mutate(cancer_type = FctCancerType(cancer_type)) %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>% 
  ggplot(aes(x = ckd_cat, y = count, fill = tx_smry)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Spectral") +
  facet_wrap(~ cancer_type, nrow = 1) +
  labs(x = "eGFR category (mL/min/1.73m2)", y = "Percentage", fill = "", title = "Female") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45))
tx_smry_plot_f

tx_smry_plot <- ggarrange(tx_smry_plot_m, tx_smry_plot_f, 
                          nrow = 2, 
                          common.legend = TRUE)
tx_smry_plot

ggsave("Outputs/tx_smry_plot.png",
       tx_smry_plot,
       height = 8, width = 16)

## BAR PLOTS of treatment types by METASTATIC DISEASE ----

# create table of counts for plot
tx_meta_counts_m <- cohort6 %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>% 
  FilterSmallGroups(any_tx, cancer_type, m_stage) %>% 
  select(sex_name, ckd_cat, any_tx, cancer_type, m_stage) %>% 
  mutate(any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type),
         m_stage = FctMetastatic(m_stage)) %>% 
  group_by(sex_name, any_tx, cancer_type, m_stage) %>% 
  summarise(count = n()) #
tx_meta_counts_f <- cohort6 %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>% 
  FilterSmallGroups(any_tx, cancer_type, m_stage) %>% 
  select(sex_name, ckd_cat, any_tx, cancer_type, m_stage) %>% 
  mutate(any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type),
         m_stage = FctMetastatic(m_stage)) %>% 
  group_by(sex_name, any_tx, cancer_type, m_stage) %>% 
  summarise(count = n()) #
tx_meta_counts <- bind_rows(
  tx_meta_counts_m,
  tx_meta_counts_f
)
tx_meta_counts %>% filter(count < 10) %>%  glimpse() # no small groups

# save file with counts for tx_meta_plot
tx_meta_counts %>%  write_csv("Outputs/tx_meta_counts.csv")

# male plot
tx_meta_plot_m <- tx_meta_counts %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>% 
  ggplot(aes(x = m_stage, y = count, fill = any_tx)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("lightblue", "royalblue")) +
  facet_wrap(~ cancer_type, nrow = 1) +
  labs(x = "Metastatic disease status", y = "Percentage", fill = "", title = "Male") +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme_bw()
tx_meta_plot_m

# female plot
tx_meta_plot_f <- tx_meta_counts %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>% 
  ggplot(aes(x = m_stage, y = count, fill = any_tx)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("pink", "hotpink")) +
  facet_wrap(~ cancer_type, nrow = 1) +
  labs(x = "Metastatic disease status", y = "Percentage", fill = "", title = "Female") +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme_bw()
tx_meta_plot_f

tx_meta_plot <- ggarrange(tx_meta_plot_m, tx_meta_plot_f, 
                          nrow = 2, 
                          common.legend = FALSE)
tx_meta_plot

ggsave("Outputs/tx_meta_plot.png",
       tx_meta_plot,
       height = 10, width = 18)

### VIOLIN plots by age at treatment ----

## plot of age/treatment

age_rx_m <- cohort6 %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>%   
  mutate(age_cat = cut(age, breaks = 3), 
         any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type)) %>% 
  group_by(any_tx, cancer_type) %>% 
  mutate(lower = quantile(age, 0.025),
         upper = quantile (age, 0.975)) %>% 
  filter(age >= lower & age <= upper) %>% 
  ggplot(aes(x = any_tx, y = age, fill = any_tx)) + 
  geom_violin(alpha = 0.8, trim = TRUE) + 
  geom_boxplot(width = 0.1, coef = 0, outlier.shape = NA) +
  scale_fill_manual(values = c("lightblue", "royalblue")) +
  facet_wrap(~cancer_type, nrow = 1) + 
  labs(x = "", y = "Age (years)", fill = "", title = "Male") + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme_bw()

age_rx_f <- cohort6 %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>%   
  mutate(age_cat = cut(age, breaks = 3), 
         any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type)) %>% 
  group_by(any_tx, cancer_type) %>% 
  mutate(lower = quantile(age, 0.025),
         upper = quantile (age, 0.975)) %>% 
  filter(age >= lower & age <= upper) %>% 
  ggplot(aes(x = any_tx, y = age, fill = any_tx)) + 
  geom_violin(alpha = 0.8, trim = TRUE) + 
  geom_boxplot(width = 0.1, coef = 0, outlier.shape = NA) +
  scale_fill_manual(values = c("pink", "hotpink")) +
  facet_wrap(~cancer_type, nrow = 1) + 
  labs(x = "", y = "Age (years)", fill = "", title = "Female") + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme_bw()

age_rx <- ggarrange(age_rx_m, age_rx_f, 
                    nrow = 2, 
                    common.legend = TRUE)
age_rx

ggsave("Outputs/age_rx.png",
       age_rx,
       height = 10, width = 16)

### VIOLIN plots by comorbidity count at treatment ----

mm_rx_m <- cohort6 %>% 
  filter(sex_name == "Male" & cancer_type != "Breast") %>%   
  mutate(any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type)) %>% 
  group_by(any_tx, cancer_type) %>% 
  mutate(upper = quantile(mm_count, 0.975, na.rm = TRUE),
         lower = quantile(mm_count, 0.025, na.rm = TRUE)) %>% 
  filter(mm_count >= lower & mm_count <= upper) %>% 
  ggplot(aes(x = any_tx, y = mm_count, fill = any_tx)) + 
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.1, coef = 0, outlier.shape = NA) +
  scale_fill_manual(values = c("lightblue", "royalblue")) +
  facet_wrap(~cancer_type, nrow = 1) + 
  labs(x = "", y = "Number of comorbidities", fill = "", title = "Male") + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme_bw()

mm_rx_f <- cohort6 %>% 
  filter(sex_name == "Female" & cancer_type != "Prostate") %>%   
  mutate(any_tx = FctAnyTreatment(any_tx),
         cancer_type = FctCancerType(cancer_type)) %>% 
  group_by(any_tx, cancer_type) %>% 
  mutate(upper = quantile(mm_count, 0.975, na.rm = TRUE),
         lower = quantile(mm_count, 0.025, na.rm = TRUE)) %>% 
  filter(mm_count >= lower & mm_count <= upper) %>% 
  ggplot(aes(x = any_tx, y = mm_count, fill = any_tx)) + 
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.1, coef = 0, outlier.shape = NA) +
  scale_fill_manual(values = c("pink", "hotpink")) +
  facet_wrap(~cancer_type, nrow = 1) + 
  labs(x = "", y = "Number of comorbidities", fill = "", title = "Female") + 
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45)) + 
  theme_bw()

mm_rx <- ggarrange(mm_rx_m, mm_rx_f, 
                   nrow = 2, 
                   common.legend = TRUE)
mm_rx

ggsave("Outputs/mm_rx.png",
       mm_rx,
       height = 10, width = 16)

## SPLINES of any treatment by egfr ----

pred_logmod_m <- read_csv("Outputs/pred_logmod_m.csv") %>% 
  mutate(cancer_type = FctCancerType(cancer_type))
pred_logmod_f <- read_csv("Outputs/pred_logmod_f.csv") %>% 
  mutate(cancer_type = FctCancerType(cancer_type))

# male
male_logplot <- ggplot(pred_logmod_m, aes(x = egfr, y = estimate, color = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = "royalblue") +
  scale_fill_manual(values = "royalblue") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male", color = "", fill = "") +
  theme_bw()
male_logplot

# female
female_logplot <- ggplot(pred_logmod_f, aes(x = egfr, y = estimate, color = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = "hotpink") +
  scale_fill_manual(values = "hotpink") +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female", colour = "", fill = "") +
  theme_bw()
female_logplot

comb_logplot <- ggarrange(male_logplot, female_logplot,
                          nrow = 2, 
                          common.legend = TRUE)
comb_logplot 

ggsave(comb_logplot, 
       file = "Outputs/spline_plot_any_tx_by_egfr.png",
       height = 8, width = 16)

# female vs. male plot

pred_logmod <- read_csv("Outputs/pred_logmod.csv") %>% 
  filter(cancer_type != "Breast" & cancer_type != "Prostate") |>
  mutate(cancer_type = FctCancerType(cancer_type))

female_vs_male_logplot <- ggplot(pred_logmod, aes(x = egfr, y = estimate, color = sex_name)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = sex_name), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = c("hotpink", "royalblue")) +
  scale_fill_manual(values = c("hotpink", "royalblue")) +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female vs. Male", fill = "", color = "") +
  theme_bw()
female_vs_male_logplot

ggsave(female_vs_male_logplot, 
       file = "Outputs/female_vs_male_anytx_plot.png",
       height = 4, width = 16)

## SPLINES of treatment type by egfr ----

pred_long_m <- read_csv("Outputs/pred_long_m.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))
pred_long_f <- read_csv("Outputs/pred_long_f.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))

# male
male_plot <- ggplot(pred_long_m, aes(x = egfr, y = estimate, color = tx_smry)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tx_smry), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_brewer(palette = "Spectral") +
  #scale_colour_manual(values = BlueSpectrum) + 
  #scale_fill_manual(values = BlueSpectrum) +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male", color = "", fill = "") +
  theme_bw()

# female
female_plot <- ggplot(pred_long_f, aes(x = egfr, y = estimate, color = tx_smry)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tx_smry), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_brewer(palette = "Spectral") +
  #scale_colour_manual(values = PinkSpectrum) + 
  #scale_fill_manual(values = PinkSpectrum) +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female", color = "", fill = "") + 
  theme_bw()

comb_plot <- ggarrange(male_plot, female_plot,
                       nrow = 2, 
                       common.legend = TRUE)
comb_plot 

ggsave(comb_plot, 
       file = "Outputs/spline_plot_tx_smry_by_egfr.png",
       height = 8, width = 16)

## SPLINES of 1 year survival by egfr ----

pred_logmod_m_surv1y <- read_csv("Outputs/pred_logmod_m_surv1y.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type), 
         sex_name = "Male")
pred_logmod_f_surv1y <- read_csv("Outputs/pred_logmod_f_surv1y.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type), 
         sex_name = "Female")

# male
male_plot_surv1y <- ggplot(pred_logmod_m_surv1y, aes(x = egfr, y = estimate)) +
  geom_line(size = 1.2, colour = "royalblue") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "royalblue") +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male") +
  theme_bw()

# female
female_plot_surv1y <- ggplot(pred_logmod_f_surv1y, aes(x = egfr, y = estimate)) +
  geom_line(size = 1.2, colour = "hotpink") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "hotpink") +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female") +
  theme_bw()

comb_plot_surv1y <- ggarrange(male_plot_surv1y, female_plot_surv1y,
                       nrow = 2, 
                       common.legend = TRUE)
comb_plot_surv1y 

ggsave(comb_plot_surv1y, 
       file = "Outputs/spline_plot_surv1y_by_egfr.png",
       height = 8, width = 16)


## SPLINES of 5 year survival by egfr ----

pred_logmod_m_surv5y <- read_csv("Outputs/pred_logmod_m_surv5y.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type), 
         sex_name = "Male")
pred_logmod_f_surv5y <- read_csv("Outputs/pred_logmod_f_surv5y.csv") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type), 
         sex_name = "Female")

# male
male_plot_surv5y <- ggplot(pred_logmod_m_surv5y, aes(x = egfr, y = estimate)) +
  geom_line(size = 1.2, colour = "royalblue") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "royalblue") +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Male") +
  theme_bw()

# female
female_plot_surv5y <- ggplot(pred_logmod_f_surv5y, aes(x = egfr, y = estimate)) +
  geom_line(size = 1.2, colour = "hotpink") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "hotpink") +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  ylim(c(0, 1)) +
  facet_wrap(~cancer_type, nrow = 1) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female") +
  theme_bw()

comb_plot_surv5y <- ggarrange(male_plot_surv5y, female_plot_surv5y,
                              nrow = 2, 
                              common.legend = TRUE)
comb_plot_surv5y 

ggsave(comb_plot_surv5y, 
       file = "Outputs/spline_plot_surv5y_by_egfr.png",
       height = 8, width = 16)


## SPLINES of male vs. female comparison for treatment type by eGFR ----

pred_long <- read_csv("Outputs/pred_long.csv") %>% 
  filter(cancer_type != "Breast" & cancer_type != "Prostate") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))

# female vs. male plot
female_vs_male_tx_smry_plot <- ggplot(pred_long, aes(x = egfr, y = estimate, color = sex_name)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = sex_name), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = c("hotpink", "royalblue")) +
  scale_fill_manual(values = c("hotpink", "royalblue")) +
  ylim(c(0, 1)) +
  facet_grid(tx_smry ~ cancer_type) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female vs. Male", color = "", fill = "") 

female_vs_male_tx_smry_plot

ggsave(female_vs_male_tx_smry_plot, 
       file = "Outputs/female_vs_male_tx_smry_plot.png",
       height = 12, width = 16)


## SPLINES of male vs. female comparison for 1 year survival by eGFR ----

pred_logmod_surv1y <- read_csv("Outputs/pred_logmod_surv1y.csv") %>% 
  filter(cancer_type != "Breast" & cancer_type != "Prostate") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))

# female vs. male plot
female_vs_male_surv1y_plot <- ggplot(pred_logmod_surv1y, aes(x = egfr, y = estimate, color = sex_name)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = sex_name), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = c("hotpink", "royalblue")) +
  scale_fill_manual(values = c("hotpink", "royalblue")) +
  ylim(c(0, 1)) +
  facet_wrap(~ cancer_type) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female vs. Male", color = "", fill = "") 

female_vs_male_surv1y_plot

ggsave(female_vs_male_surv1y_plot, 
       file = "Outputs/female_vs_male_surv1y_plot.png",
       height = 12, width = 16)

## SPLINES of male vs. female comparison for 5 year survival by eGFR ----

pred_logmod_surv5y <- read_csv("Outputs/pred_logmod_surv5y.csv") %>% 
  filter(cancer_type != "Breast" & cancer_type != "Prostate") %>% 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))

# female vs. male plot
female_vs_male_surv5y_plot <- ggplot(pred_logmod_surv5y, aes(x = egfr, y = estimate, color = sex_name)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = sex_name), alpha = 0.2, colour = NA) +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_vline(xintercept = c(60), linetype = "dashed") +
  geom_vline(xintercept = c(90), linetype = "solid") +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150)) +
  scale_colour_manual(values = c("hotpink", "royalblue")) +
  scale_fill_manual(values = c("hotpink", "royalblue")) +
  ylim(c(0, 1)) +
  facet_wrap(~ cancer_type) +
  labs(x = "eGFR (mL/min/1.73m2)", y = "Adjusted probability", title = "Female vs. Male", color = "", fill = "") 

female_vs_male_surv5y_plot

ggsave(female_vs_male_surv5y_plot, 
       file = "Outputs/female_vs_male_surv5y_plot.png",
       height = 12, width = 16)

