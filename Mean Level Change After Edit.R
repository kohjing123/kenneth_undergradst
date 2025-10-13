### START OF CODE ###

#Packages
library(readxl)
library(metafor)
library(dplyr)
library(tidyr)
library(stringr)
library(clubSandwich)

#Load and clean data, and calculate beckers d and intermediate step
dat <- 
  read_xlsx("STAnalysisNew.xlsx", sheet = "Acceptance") %>% 
  mutate(
    MeanT1 = as.numeric(MeanT1), 
    MeanT2 = as.numeric(MeanT2),
    MeanT3 = as.numeric(MeanT3),
    SDT1   = as.numeric(SDT1),
    SDT2   = as.numeric(SDT2),
    SDT3   = as.numeric(SDT3),
    SampleSize = as.numeric(SampleSize),
    beckersdT1T2 = (MeanT2-MeanT1)/SDT1,
    beckersdT1T3 = (MeanT3-MeanT1)/SDT1,
    beckersdT2T3 = (MeanT3-MeanT2)/SDT2,
    intermediate = 1 - (3/(4*(SampleSize-1)-1))
  )

#Calculate beckers d variance
dat <- dat %>% mutate(
  varbeckersdT1T2 = ((1 + SampleSize*beckersdT1T2^2)/SampleSize)*((SampleSize-1)/(SampleSize-3)) -
    (beckersdT1T2^2/(intermediate*(SampleSize-1))^2),
  varbeckersdT1T3 = ((1 + SampleSize*beckersdT1T3^2)/SampleSize)*((SampleSize-1)/(SampleSize-3)) -
    (beckersdT1T3^2/(intermediate*(SampleSize-1))^2),
  varbeckersdT2T3 = ((1 + SampleSize*beckersdT2T3^2)/SampleSize)*((SampleSize-1)/(SampleSize-3)) -
    (beckersdT2T3^2/(intermediate*(SampleSize-1))^2)
) %>% filter(
  is.na(beckersdT1T2) | (beckersdT1T2 > -2 & beckersdT1T2 < 2),
  is.na(beckersdT1T3) | (beckersdT1T3 > -2 & beckersdT1T3 < 2),
  is.na(beckersdT2T3) | (beckersdT2T3 > -2 & beckersdT2T3 < 2)
) 

#Reshape data to long format and set up moderators for analysis later
dat_long <- dat %>%
  pivot_longer(
    cols = starts_with("beckersd"),
    names_to = "comparison",
    values_to = "d"
  ) %>%
  mutate(
    var = case_when(
      comparison == "beckersdT1T2" ~ varbeckersdT1T2,
      comparison == "beckersdT1T3" ~ varbeckersdT1T3,
      comparison == "beckersdT2T3" ~ varbeckersdT2T3,
      TRUE ~ NA_real_
    ),
    comp_clean = stringr::str_remove(comparison, "^beckersd"),
    t1 = as.numeric(stringr::str_match(comp_clean, "T(\\d)T(\\d)")[,2]),
    t2 = as.numeric(stringr::str_match(comp_clean, "T(\\d)T(\\d)")[,3]),
    lag = t2 - t1,
    StudyID = paste(Author, `Year Collected`, sep = "_"),
    QA_mod  = ifelse(Qatest == "High", 1, 0),
    Inv_mod = ifelse(`Reports Invariance?` == "YES", 1, 0),
    FullInv_mod = ifelse(AchievedFullScalarInvariance == "YES", 1, 0),
    WestEast_mod = ifelse(Region %in% c("America", "Europe", "Oceania"), 1,
                          ifelse(Region == "Asia", 0, NA_real_))
  ) %>%
  filter(!is.na(d), !is.na(var))

#Multilevel meta across timepoints
m_lag <- rma.mv(yi = d, 
                V = var, 
                mods = ~ lag, 
                random = ~ 1 | StudyID, 
                data = dat_long)
summary(m_lag)

#Moderators
#Control age mean 
m_lag_age <- rma.mv(yi = d, 
                    V = var, 
                    mods = ~ lag + AgeMeanT1, 
                    random = ~ 1 | StudyID, 
                    data = dat_long)
summary(m_lag_age)

#Control age + sd
m_lag_age_sd <- rma.mv(yi = d, 
                       V = var, 
                       mods = ~ lag + AgeMeanT1 + AgeSDT1, 
                       random = ~ 1 | StudyID, 
                       data = dat_long)
summary(m_lag_age_sd)

#QA moderator 
m_QA  <- rma.mv(yi = d, V = var, 
                mods = ~ lag * QA_mod, 
                random = ~ 1 | StudyID, 
                data = dat_long)
summary(m_QA)

#Reported Invariance moderator
m_Inv <- rma.mv(yi = d, 
                V = var, 
                mods = ~ lag * Inv_mod, 
                random = ~ 1 | StudyID, 
                data = dat_long)
summary(m_INV)

#Full Invariance moderator
m_FullInv <- rma.mv(yi = d, 
                V = var, 
                mods = ~ lag * FullInv_mod, 
                random = ~ 1 | StudyID, 
                data = dat_long)
summary(m_FullInv)

#USA vs others moderator
m_westeast <- rma.mv(yi = d, 
                V = var, mods = ~ lag * WestEast_mod, 
                random = ~ 1 | StudyID, 
                data = dat_long)
summary(m_USA)

#Sensitivity analysis
#Cluster-robust SEs
cr2 <- function(mod) clubSandwich::coef_test(mod, vcov = vcovCR(mod, cluster = dat_long$StudyID, type = "CR2"))
print(cr2(m_lag))
print(cr2(m_QA))
print(cr2(m_Inv))
print(cr2(m_FullInv))
print(cr2(m_westeast))
print(cr2(m_lag_age))
print(cr2(m_lag_age_sd))

extract_cr2 <- function(model, model_name, coef_name = "lag") {
  res <- tryCatch({
    cr <- clubSandwich::coef_test(model, vcov = clubSandwich::vcovCR(model, cluster = dat_long$StudyID, type = "CR2"))
    if (coef_name %in% rownames(cr)) {
      data.frame(
        Model = model_name,
        Coefficient = coef_name,
        Estimate = cr[coef_name, "beta"],
        SE = cr[coef_name, "SE"],
        CI_lb = cr[coef_name, "CI_L"],
        CI_ub = cr[coef_name, "CI_U"],
        p = cr[coef_name, "p_Satt"],
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }, error = function(e) {
    message(paste("Error extracting", model_name, ":", e$message))
    NULL
  })
  return(res)
}

#Bootstrapping for lag
set.seed(123)
B <- 1000
studies <- unique(dat_long$StudyID)
boot_beta <- numeric(B)
for (b in 1:B) {
  s_ids <- sample(studies, length(studies), replace = TRUE)
  dat_b <- do.call(rbind, lapply(s_ids, function(id) dat_long[dat_long$StudyID == id, , drop=FALSE]))
  fit_b <- try(rma.mv(yi = d, V = var, mods = ~ lag, random = ~ 1 | StudyID, data = dat_b), silent=TRUE)
  boot_beta[b] <- if (inherits(fit_b,"try-error")) NA_real_ else as.numeric(coef(fit_b)["lag"])
}
boot_beta <- boot_beta[is.finite(boot_beta)]
boot_ci <- quantile(boot_beta, c(.025,.5,.975), na.rm=TRUE)
print(boot_ci)

#Publication bias using PET-PEESE
dat_long$sei <- sqrt(dat_long$var)
  #PET
  pet <- rma.mv(yi = d, V = var, mods = ~ sei, random = ~ 1 | StudyID, data = dat_long)
  print(clubSandwich::coef_test(pet, vcov = clubSandwich::vcovCR(pet, cluster = dat_long$StudyID, type = "CR2")))
  #If PET slope is significant (> 0), run PEESE
  peese <- rma.mv(yi = d, V = var, mods = ~ I(sei^2), random = ~ 1 | StudyID, data = dat_long)
  print(clubSandwich::coef_test(peese, vcov = clubSandwich::vcovCR(peese, cluster = dat_long$StudyID, type = "CR2")))

#Combine all results into a summary table
ml_summary_table <- do.call(rbind, list(
  extract_cr2(m_lag, "Lag"),
  extract_cr2(m_lag_age, "Lag + AgeMeanT1"),
  extract_cr2(m_lag_age_sd, "Lag + AgeMeanT1 + AgeSDT1"),
  extract_cr2(m_QA, "Lag * QA"),
  extract_cr2(m_Inv, "Lag * Inv"),
  extract_cr2(m_FullInv, "Lag * FullInv"),
  extract_cr2(m_westeast, "Lag * WestEast"),
  extract_cr2(pet, "PET (Lag ~ SE)"),
  extract_cr2(peese, "PEESE (Lag ~ SE^2)")
  ))
  
print(ml_summary_table)
write.csv(ml_summary_table, "meanlevelchange.csv", row.names = FALSE)
  

### END OF CODE ###
