source("Scripts/00_packages_functions.R")

library(ggpubr)

## SPLINES of treatment type by egfr ----

pred_long_m <- read_csv("Outputs/pred_long_m.csv") |> 
  mutate(tx_smry = FctTreatmentSmry(group),
         cancer_type = FctCancerType(cancer_type))
pred_long_f <- read_csv("Outputs/pred_long_f.csv") |> 
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
       file = "Outputs/spline_plot_tx_smry_by_egfr_for_wcn.png",
       device = ragg::agg_png,
       width = 6.4, height = 4.8, units = "in", dpi = 100)
