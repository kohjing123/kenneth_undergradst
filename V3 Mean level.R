### START OF CODE ###

#R version 4.5.1
#Load / Install Packages
library(readxl) #version 1.4.5
library(metafor) #version 4.8.0
library(dplyr) #version 4.8.0
library(tidyr) #version 1.3.1
library(lmerTest) #version 3.1.3
library(ggplot2) #version 3.5.2
library(stringr) #version 1.6.0

##Becker's d and Var Calculation & Multilevel meta ----
#Load and clean data
dat <- 
  readxl::read_xlsx("STAnalysisNEW.xlsx", sheet = "Rumination") %>% 
  dplyr::mutate(
    
    #prep nesting variables
    `Paper ID`  = as.factor(`Paper ID`),
    `Sample ID` = as.factor(`Sample ID`),
    `Meta ID`   = as.factor(`Meta ID`),
    
    #fix numeric columns 
    dplyr::across(
      c(
        MeanT1, MeanT2, MeanT3, 
        SDT1, SDT2, SDT3,
        SampleSize,
        rT1T2, rT1T3, rT2T3
      ),
      as.numeric
    ),
    
    #calc beckers d
    beckersdT1T2 = (MeanT2 - MeanT1) / SDT1,
    beckersdT1T3 = (MeanT3 - MeanT1) / SDT1,
    beckersdT2T3 = (MeanT3 - MeanT2) / SDT2,
    
    #calc variance of beckers d
    varbeckersdT1T2 = (2 * (1 - rT1T2) / SampleSize) + 
      (beckersdT1T2^2 / (2 * SampleSize)),
    
    varbeckersdT1T3 = (2 * (1 - rT1T3) / SampleSize) + 
      (beckersdT1T3^2 / (2 * SampleSize)),
    
    varbeckersdT2T3 = (2 * (1 - rT2T3) / SampleSize) + 
      (beckersdT2T3^2 / (2 * SampleSize))
  )

#Reshape data to long format and set up moderators for analysis later
dat_long <- dat %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("beckersd"),
    names_to = "comparison",
    values_to = "d"
  ) %>%
  dplyr::mutate(
    
    #Match each Becker's d to its variance
    var = dplyr::case_when(
      comparison == "beckersdT1T2" ~ varbeckersdT1T2,
      comparison == "beckersdT1T3" ~ varbeckersdT1T3,
      comparison == "beckersdT2T3" ~ varbeckersdT2T3,
      TRUE ~ NA_real_
    ),
    
    #Clean comparison label and extract timepoints
    comp_clean = sub("^beckersd", "", comparison),
    t1 = as.numeric(sub("T([0-9])T([0-9])", "\\1", comp_clean)),
    t2 = as.numeric(sub("T([0-9])T([0-9])", "\\2", comp_clean)),
    lag = t2 - t1,
    
    #Keep original clustering structure for multilevel models
    `Paper ID`  = as.factor(`Paper ID`),
    `Sample ID` = as.factor(`Sample ID`),
    `Meta ID`   = as.factor(`Meta ID`),
    
    #Moderators
    QA_mod = dplyr::case_when(
      Qatest %in% c("HIGH", "High", "high") ~ 1,
      Qatest %in% c("LOW", "Low", "low") ~ 0,
      TRUE ~ NA_real_
    ),
    WestEast_mod = dplyr::case_when(
      Region %in% c("America", "Europe", "Oceania") ~ 1,
      Region == "Asia" ~ 0,
      TRUE ~ NA_real_
    ),
    
    #Numeric moderators
    AgeMeanT1 = suppressWarnings(as.numeric(AgeMeanT1)),
    AgeSDT1   = suppressWarnings(as.numeric(AgeSDT1))
  ) %>%
  dplyr::filter(!is.na(d), !is.na(var))

#descriptives
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
dat_long <- dat_long %>%
  dplyr::mutate(
    d_year     = d / lag,
    var_d_year = var / (lag^2)
  )

#Pooled d_year ----
pooled_d_year <- metafor::rma.mv(
  yi     = d_year,
  V      = var_d_year,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`, #between-study, within-study/between-sample, and within-sample clustering
  data   = dat_long,
  method = "REML"
) 
summary(pooled_d_year)

##Scatterplot of d over mean age with cumulative d line graph For now Optional ----
dat_long <- dat_long %>%
  dplyr::mutate(sampleMeanAge = AgeMeanT1) %>%
  dplyr::filter(!is.na(sampleMeanAge)) %>%
  dplyr::arrange(sampleMeanAge) %>%
  dplyr::mutate(yi_cumulative = as.numeric(cummean(d_year)))

#Scatterplot
ggplot(dat_long, aes(x = sampleMeanAge)) +
  geom_line(aes(y = yi_cumulative)) +
  geom_point(aes(y = d_year), size = 2.5) +
  scale_x_continuous("Sample Mean Age", breaks = seq(0, 55, 5)) +
  labs(
    y = "Becker's d per year",
    title = "Rumination"
  ) +
  theme(panel.grid.minor.x = element_blank())

##Moderators ----
#Mean Age
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
cor(dat_long$AgeMeanT1, dat_long$AgeSDT1, use = "complete.obs") #Since cor=0.98, we will not test AgeSDT1 due to multicollinearity (unable to uniquely attribute variance to each predictor)

#QA
m_QA <- metafor::rma.mv(
  yi     = d_year,
  V      = var_d_year,
  mods   = ~ QA_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(QA_mod)),
  method = "REML"
)
summary(m_QA)

#Exploratories
m_westeast <- metafor::rma.mv(
  yi     = d_year,
  V      = var_d_year,
  mods   = ~ WestEast_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(WestEast_mod)),
  method = "REML"
)
summary(m_westeast)

##Sensitivity analysis / Pub Bias Tests ----
#Egger's test
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

#Leave-one-out
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