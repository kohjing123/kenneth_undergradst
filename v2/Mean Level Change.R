### START OF CODE ###

#R version 4.5.1
#Load / Install Packages
library(readxl) #version 1.4.5
library(metafor) #version 4.8.0
library(dplyr) #version 4.8.0
library(lmerTest) #version 3.1.3

#Load and clean data
dat <- 
  readxl::read_xlsx("STAnalysisNew.xlsx", sheet = "Acceptance") %>% 
  dplyr::mutate(
    # fix column formats
    dplyr::across(
      c(
        MeanT1, MeanT2, MeanT3, 
        SDT1, SDT2, SDT3,
        SampleSize
      ),
      as.numeric
    ),
    # calc adjustment factor
    adjFactor = 1 - (3 / (4 * (SampleSize - 1) - 1)),
    # calc multiplication factor
    mulFactor = (SampleSize - 1) / (SampleSize - 3),
    # calc beckers d
    beckersdT1T2 = (MeanT2 - MeanT1) / SDT1,
    beckersdT1T3 = (MeanT3 - MeanT1) / SDT1,
    beckersdT2T3 = (MeanT3 - MeanT2) / SDT2,
    # calc variance of beckers d
    varbeckersdT1T2 = ((1 + SampleSize * beckersdT1T2^2) / SampleSize) * (mulFactor) -
      (beckersdT1T2^2 / (adjFactor * (SampleSize - 1))^2),
    varbeckersdT1T3 = ((1 + SampleSize * beckersdT1T3^2) / SampleSize) * (mulFactor) -
      (beckersdT1T3^2 / (adjFactor * (SampleSize - 1))^2),
    varbeckersdT2T3 = ((1 + SampleSize * beckersdT2T3^2) / SampleSize) * (mulFactor) -
      (beckersdT2T3^2 / (adjFactor * (SampleSize - 1))^2)
  )

#Reshape data to long format and set up moderators for analysis later
dat_long <- dat %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("beckersd"),
    names_to = "comparison",
    values_to = "d"
  ) %>%
  dplyr::mutate(
    var = dplyr::case_when(
      comparison == "beckersdT1T2" ~ varbeckersdT1T2,
      comparison == "beckersdT1T3" ~ varbeckersdT1T3,
      comparison == "beckersdT2T3" ~ varbeckersdT2T3,
      TRUE ~ NA_real_
    ),
    comp_clean = stringr::str_remove(comparison, "^beckersd"),
    t1 = as.numeric(stringr::str_match(comp_clean, "T(\\d)T(\\d)")[, 2]),
    t2 = as.numeric(stringr::str_match(comp_clean, "T(\\d)T(\\d)")[, 3]),
    lag = t2 - t1,
    StudyID = paste(Author, `Year Collected`, sep = "_"),
    QA_mod = dplyr::case_when(Qatest == "High" ~ 1, TRUE ~ 0),
    Inv_mod = dplyr::case_when(`Reports Invariance?` == "YES" ~ 1, TRUE ~ 0),
    FullInv_mod = dplyr::case_when(AchievedFullScalarInvariance == "YES" ~ 1, TRUE ~ 0),
    WestEast_mod = dplyr::case_when(
      Region %in% c("America", "Europe", "Oceania") ~ 1,
      Region == "Asia" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(d), !is.na(var))

#Multilevel meta across timepoints
m_lag <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
base::summary(m_lag)

#Moderators
#Control age mean 
m_lag_age <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * AgeMeanT1, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
base::summary(m_lag_age)

#Control age + sd
m_lag_age_sd <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * AgeMeanT1 + lag * AgeSDT1, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
base::summary(m_lag_age_sd)

#QA moderator 
m_QA <- metafor::rma.mv(
  yi = d, V = var, 
  mods = ~ lag * QA_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
base::summary(m_QA)

#Reported Invariance moderator
m_Inv <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * Inv_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
base::summary(m_Inv)

#Full Invariance moderator
m_FullInv <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * FullInv_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
base::summary(m_FullInv)

#West v East moderator
m_westeast <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * WestEast_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
base::summary(m_westeast)

#Sensitivity analysis
#Funnel-plot
metafor::funnel(
  m_lag,
  main = "Funnel Plot (Becker’s d)",
  xlab = "Standardized Mean Change (Becker’s d)",
  ylab = "Standard Error of Becker’s d",
  pch = 21, bg = "skyblue", cex = 1.1
)

#Egger's test
  #Compute standard error and precision
  dat_long <- dat_long %>%
  mutate(
    sei = sqrt(var),
    precision = 1 / sei
  )
  #Mixed-effects Egger regression
  egger_mixed <- lmerTest::lmer(
  I(d / sei) ~ 1 + precision + (1 | StudyID),
  data = dat_long)
  #Intercept
  summary(egger_mixed, correlation = FALSE)
  
#Leave-one-out
  #Get all unique study IDs
  study_ids <- unique(dat_long$StudyID)
  #Create empty results table
  loo_results <- data.frame(
    StudyID = study_ids,
    k = NA, est = NA, se = NA, zval = NA, pval = NA
  )
  #Loop through studies, drop one at a time
  for (i in seq_along(study_ids)) {
    dat_loo <- subset(dat_long, StudyID != study_ids[i])
  #Try fitting the same model without that study
  fit_loo <- try(
      metafor::rma.mv(
        yi = d,
        V = var,
        mods = ~ lag,
        random = ~ 1 | StudyID,
        data = dat_loo,
        method = "REML"
      ),
      silent = TRUE
    )
   #If model runs successfully, store results
    if (!inherits(fit_loo, "try-error")) {
      loo_results$k[i]    <- nrow(dat_loo)
      loo_results$est[i]  <- fit_loo$b[1, 1]
      loo_results$se[i]   <- fit_loo$se[1]
      loo_results$zval[i] <- fit_loo$zval[1]
      loo_results$pval[i] <- fit_loo$pval[1]
    }
  }
  #Compare to full model estimate
  full_est <- m_lag$b[1, 1]
  loo_results$diff <- loo_results$est - full_est
  #Sort by largest effect change
  loo_results_sorted <- loo_results[order(abs(loo_results$diff), decreasing = TRUE), ]
  print(loo_results_sorted)

### END OF CODE ###



