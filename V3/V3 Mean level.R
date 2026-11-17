### START OF CODE ###

#R version 4.5.1
#Load / Install Packages
library(readxl)   #version 1.4.5; Wickham & Bryan (2025)
library(metafor)  #version 4.8.0; Viechtbauer (2010)
library(dplyr)    #version 1.1.4; Wickham et al. (2023)
library(tidyr)    #version 1.3.1; Wickham et al. (2024)
library(lmerTest) #version 3.1.3; Kuznetsova et al. (2017)
library(ggplot2)  #version 3.5.2; Wickham (2016)
library(stringr)  #version 1.6.0; Wickham (2025)

##Becker's d and Var Calculation & Multilevel meta ----
#Load and clean data
dat <- 
  readxl::read_xlsx("STAnalysis7.xlsx", sheet = "Rumination") %>% 
  dplyr::mutate(
    
    #prep nesting variables by making them into factors; Konstantopoulos (2011); Viechtbauer (2010)
    `Paper ID`  = as.factor(`Paper ID`),
    `Sample ID` = as.factor(`Sample ID`),
    `Meta ID`   = as.factor(`Meta ID`),
    
    #fix numeric columns; Wickham et al. (2023)
    dplyr::across(
      c(
        MeanT1, MeanT2, MeanT3, 
        SDT1, SDT2, SDT3,
        SampleSize,
        rT1T2, rT1T3, rT2T3
      ),
      as.numeric
    ),
    
    #calc Beckers d; Becker (1988)
    beckersdT1T2 = (MeanT2 - MeanT1) / SDT1,
    beckersdT1T3 = (MeanT3 - MeanT1) / SDT1,
    beckersdT2T3 = (MeanT3 - MeanT2) / SDT2,
    
    #calc variance of Beckers d; Morris & DeShon (2002)
    varbeckersdT1T2 = (2 * (1 - rT1T2) / SampleSize) + 
      (beckersdT1T2^2 / (2 * SampleSize)),
    
    varbeckersdT1T3 = (2 * (1 - rT1T3) / SampleSize) + 
      (beckersdT1T3^2 / (2 * SampleSize)),
    
    varbeckersdT2T3 = (2 * (1 - rT2T3) / SampleSize) + 
      (beckersdT2T3^2 / (2 * SampleSize))
  )

#Reshape data to long format and set up moderators for analysis later; Wickham et al. (2024); Konstantopoulos (2011); Viechtbauer (2010)
dat_long <- dat %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("beckersd"),
    names_to = "comparison",
    values_to = "d"
  ) %>%
  dplyr::mutate(
    
    #Match each Becker's d to its variance; Becker (1988); Morris & DeShon (2002)
    var = dplyr::case_when(
      comparison == "beckersdT1T2" ~ varbeckersdT1T2,
      comparison == "beckersdT1T3" ~ varbeckersdT1T3,
      comparison == "beckersdT2T3" ~ varbeckersdT2T3,
      TRUE ~ NA_real_
    ),
    
    #Clean comparison label and extract timepoints; Wickham et al. (2023)
    comp_clean = sub("^beckersd", "", comparison),
    t1 = as.numeric(sub("T([0-9])T([0-9])", "\\1", comp_clean)),
    t2 = as.numeric(sub("T([0-9])T([0-9])", "\\2", comp_clean)),
    lag = t2 - t1,
    
    #Keep original clustering structure for multilevel models; Konstantopoulos (2011); Viechtbauer (2010)
    `Paper ID`  = as.factor(`Paper ID`),
    `Sample ID` = as.factor(`Sample ID`),
    `Meta ID`   = as.factor(`Meta ID`),
    
    #Moderators; Thompson & Higgins (2002); Viechtbauer (2010)
    QA_mod = dplyr::case_when(
      Qatest %in% c("HIGH") ~ 1,
      Qatest %in% c("LOW") ~ 0,
      TRUE ~ NA_real_
    ),
    Country_mod = trimws(as.character(Country)),
    Country_mod = dplyr::na_if(Country_mod, "NA"),
    Country_mod = as.factor(Country_mod),
    
    #Numeric moderators; Thompson & Higgins (2002); Viechtbauer (2010)
    AgeMeanT1 = suppressWarnings(as.numeric(AgeMeanT1)),
    AgeSDT1   = suppressWarnings(as.numeric(AgeSDT1))
  ) %>%
  dplyr::filter(!is.na(d), !is.na(var))

#descriptives; Wickham et al. (2023)
dat_long %>%
  dplyr::distinct(`Paper ID`, `Sample ID`, SampleSize, AgeMeanT1, AgeSDT1) %>%
  dplyr::summarise(
    k = nrow(dat_long),  #total effect sizes
    
    N = sum(SampleSize, na.rm = TRUE),  #total participants (unique samples)
    
    mean_age = mean(AgeMeanT1, na.rm = TRUE),
    sd_age   = sd(AgeMeanT1, na.rm = TRUE),
    
    min_age = min(AgeMeanT1, na.rm = TRUE),
    max_age = max(AgeMeanT1, na.rm = TRUE)
  )

#d_year ----
#Convert Becker's d into standardised yearly change; Becker (1988); Morris & DeShon (2002)
dat_long <- dat_long %>%
  dplyr::mutate(
    d_year     = d / lag,
    var_d_year = var / (lag^2)
  )

#Pooled d_year ----
#Multilevel random-effects meta-analysis; Konstantopoulos (2011); Viechtbauer (2010)
pooled_d_year <- metafor::rma.mv(
  yi     = d_year,
  V      = var_d_year,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`, #between-study, within-study/between-sample, and within-sample clustering
  data   = dat_long,
  method = "REML"
) 
summary(pooled_d_year)

##Moderators ----
#Mean Age; Thompson & Higgins (2002); Viechtbauer (2010)
dat_long <- dat_long %>%
  dplyr::mutate(
    AgeMeanT1_c = scale(AgeMeanT1, center = TRUE, scale = FALSE)[, 1]
  )

m_age <- metafor::rma.mv(
  yi     = d_year,
  V      = var_d_year,
  mods   = ~ AgeMeanT1_c,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long,
  method = "REML"
)
summary(m_age)

#Check multicollinearity between age variables; Thompson & Higgins (2002)
cor(dat_long$AgeMeanT1, dat_long$AgeSDT1, use = "complete.obs") #Since cor=0.98, we will not test AgeSDT1 due to multicollinearity (unable to uniquely attribute variance to each predictor)

#QA moderator; Thompson & Higgins (2002); Viechtbauer (2010)
m_QA <- metafor::rma.mv(
  yi     = d_year,
  V      = var_d_year,
  mods   = ~ QA_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(QA_mod)),
  method = "REML"
)
summary(m_QA)

#Exploratory country moderator; Thompson & Higgins (2002); Viechtbauer (2010)
m_country <- metafor::rma.mv(
  yi     = d_year,
  V      = var_d_year,
  mods   = ~ Country_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(Country_mod)),
  method = "REML"
)

summary(m_country)

##Sensitivity analysis / Pub Bias Tests ----
#Egger's test for small-study effects; Egger et al. (1997); Sterne & Egger (2001); Kuznetsova et al. (2017)
dat_long <- dat_long %>%
  dplyr::mutate(
    sei_year       = sqrt(var_d_year),
    precision_year = 1 / sei_year
  )

egger_mixed <- lmerTest::lmer(
  I(d_year / sei_year) ~ 1 + precision_year + (1 | `Paper ID`/`Sample ID`/`Meta ID`),
  data = dat_long
)

summary(egger_mixed, correlation = FALSE)

#Leave-one-out sensitivity analysis; Viechtbauer & Cheung (2010)
study_ids <- unique(dat_long$`Paper ID`)

loo_results <- data.frame(
  PaperID = study_ids,
  k       = NA,
  est     = NA,
  se      = NA,
  zval    = NA,
  pval    = NA
)

for (i in seq_along(study_ids)) {
  dat_loo <- subset(dat_long, `Paper ID` != study_ids[i])
  
  fit_loo <- try(
    metafor::rma.mv(
      yi     = d_year,
      V      = var_d_year,
      random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
      data   = dat_loo,
      method = "REML"
    ),
    silent = TRUE
  )
  
  if (!inherits(fit_loo, "try-error")) {
    loo_results$k[i]    <- nrow(dat_loo)
    loo_results$est[i]  <- fit_loo$b[1, 1]
    loo_results$se[i]   <- fit_loo$se[1]
    loo_results$zval[i] <- fit_loo$zval[1]
    loo_results$pval[i] <- fit_loo$pval[1]
  }
}

full_est <- pooled_d_year$b[1, 1]

loo_results$diff <- loo_results$est - full_est
loo_results_sorted <- loo_results[order(abs(loo_results$diff), decreasing = TRUE), ]

print(loo_results_sorted)

## END OF CODE ###
