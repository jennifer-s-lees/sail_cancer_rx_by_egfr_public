library(tidyverse)
library(lme4)
library(nephro)
library(RODBC)
library(tcltk)
library(data.table)
library(broom)
library(RColorBrewer)
library(tableone)
library(ggpubr)
library(nnet)
library(splines)
library(reshape2)
library(marginaleffects)
library(boot)


#### Login to SAIL and identify tables ----
getLogin <- function(userName=''){
  require(tcltk);
  wnd <- tktoplevel();
  tclVar(userName) -> user;
  tclVar("") -> passVar;
  
  tkgrid(tklabel(wnd, text="Username:"));
  tkgrid(tkentry(wnd, textvariable=user) -> passBox);
  
  tkgrid(tklabel(wnd, text="Password:"));
  tkgrid(tkentry(wnd, textvariable=passVar, show="*") -> passBox);
  
  tkbind(passBox, "<Return>",function()
    tkdestroy(wnd));
  
  tkwait.window(wnd);
  password <- tclvalue(passVar);
  userName <- tclvalue(user);
  return(c(userName, password));
}


# Cancer types etc

FctCancerType <- function(x) {
  factor(x, 
         levels = c("Breast", "Prostate", "Colorectal", "Lung and bronchus", "Melanoma", "Renal/urinary tract"))
}

FctCKDCat <- function(x) {
  factor(x, 
         levels = c("KRT", "<=15", ">15-30", ">30-45", ">45-60", ">60-90", ">90-105", ">105"))
}

## change factor levels for tx_smry ----

FctAnyTreatment <- function(x) {
  factor(x, 
         levels = c("None", "Treatment"),
         labels = c("None", "Treatment"))
}

FctTreatmentSmry <- function(x) {
  factor(x, 
         levels = c("None", "RTx.only", "SACT.only", "Surg.only", "Treatment.combination",
                    "RTx only", "SACT only", "Surg only", "Treatment combination"),
         labels = c("None", "Radiotherapy only", "SACT only", "Surgery only", "Treatment combination",
                    "Radiotherapy only", "SACT only", "Surgery only", "Treatment combination"))
}

FctMetastatic <- function(x) {
  factor(x,
         levels = c("1", "0", "Unknown"),
         labels = c("Metastatic", "Not metastatic", "Unknown"))
}

FilterSmallGroups <- function(data, ..., threshold = 10) {
  # capture grouping variables
  group_vars <- rlang::enquos(...)
  # count totals for each combination of the grouping variables
  combo_totals <- data %>% 
    dplyr::count(!!!group_vars) %>% 
    dplyr::filter(n >= threshold)
  # keep only combinations that meet the threshold
  data %>% 
    dplyr::semi_join(combo_totals, by = purrr::map_chr(group_vars, rlang::as_name))
}

## Colour spectra ---

PinkSpectrum <- c("black", "lightpink", "hotpink", "maroon4", "plum4")
BlueSpectrum <- c("black", "paleturquoise3", "lightseagreen", "royalblue", "midnightblue")


