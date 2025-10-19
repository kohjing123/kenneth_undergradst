### START OF CODE ###

#R version 4.5.1
#Load / Install Packages
library(readxl) #version 1.4.5
library(metafor) #version 4.8.0
library(dplyr) #version 4.8.0
library(tidyr) #version 1.3.1
library(stringr) #version 1.5.1
library(clubSandwich) #version 0.6.1

#Load and clean data, and calculate beckers d and intermediate step
dat <- 
  readxl::read_xlsx("STAnalysisNew.xlsx", sheet = "Acceptance") %>% 
  dplyr::mutate(
    MeanT1 = as.numeric(MeanT1), 
    MeanT2 = as.numeric(MeanT2),
    MeanT3 = as.numeric(MeanT3),
    SDT1   = as.numeric(SDT1),
    SDT2   = as.numeric(SDT2),
    SDT3   = as.numeric(SDT3),
    SampleSize = as.numeric(SampleSize),
    beckersdT1T2 = (MeanT2 - MeanT1) / SDT1,
    beckersdT1T3 = (MeanT3 - MeanT1) / SDT1,
    beckersdT2T3 = (MeanT3 - MeanT2) / SDT2,
    intermediate = 1 - (3 / (4 * (SampleSize - 1) - 1))
  )

#Calculate beckers d variance
dat <- dat %>% 
  dplyr::mutate(
    varbeckersdT1T2 = ((1 + SampleSize * beckersdT1T2^2) / SampleSize) * ((SampleSize - 1) / (SampleSize - 3)) -
      (beckersdT1T2^2 / (intermediate * (SampleSize - 1))^2),
    varbeckersdT1T3 = ((1 + SampleSize * beckersdT1T3^2) / SampleSize) * ((SampleSize - 1) / (SampleSize - 3)) -
      (beckersdT1T3^2 / (intermediate * (SampleSize - 1))^2),
    varbeckersdT2T3 = ((1 + SampleSize * beckersdT2T3^2) / SampleSize) * ((SampleSize - 1) / (SampleSize - 3)) -
      (beckersdT2T3^2 / (intermediate * (SampleSize - 1))^2)
  ) %>% 
  dplyr::filter(
    is.na(beckersdT1T2) | (beckersdT1T2 > -2 & beckersdT1T2 < 2),
    is.na(beckersdT1T3) | (beckersdT1T3 > -2 & beckersdT1T3 < 2),
    is.na(beckersdT2T3) | (beckersdT2T3 > -2 & beckersdT2T3 < 2)
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
summary(m_lag)

#Moderators
#Control age mean 
m_lag_age <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag + AgeMeanT1, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
summary(m_lag_age)

#Control age + sd
m_lag_age_sd <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag + AgeMeanT1 + AgeSDT1, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
summary(m_lag_age_sd)

#QA moderator 
m_QA <- metafor::rma.mv(
  yi = d, V = var, 
  mods = ~ lag * QA_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
summary(m_QA)

#Reported Invariance moderator
m_Inv <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * Inv_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
summary(m_Inv)

#Full Invariance moderator
m_FullInv <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * FullInv_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
summary(m_FullInv)

#USA vs others moderator
m_westeast <- metafor::rma.mv(
  yi = d, 
  V = var, 
  mods = ~ lag * WestEast_mod, 
  random = ~ 1 | StudyID, 
  data = dat_long
)
summary(m_westeast)

#Sensitivity analysis
#Cluster-robust SEs
cr2 <- function(mod) clubSandwich::coef_test(mod, vcov = clubSandwich::vcovCR(mod, cluster = dat_long$StudyID, type = "CR2"))
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
  dat_b <- do.call(rbind, lapply(s_ids, function(id) dat_long[dat_long$StudyID == id, , drop = FALSE]))
  fit_b <- try(
    metafor::rma.mv(yi = d, V = var, mods = ~ lag, random = ~ 1 | StudyID, data = dat_b),
    silent = TRUE
  )
  boot_beta[b] <- if (inherits(fit_b, "try-error")) NA_real_ else as.numeric(coef(fit_b)["lag"])
}
boot_beta <- boot_beta[is.finite(boot_beta)]
boot_ci <- stats::quantile(boot_beta, c(.025, .5, .975), na.rm = TRUE)
print(boot_ci)

#Publication bias using PET-PEESE
dat_long$sei <- sqrt(dat_long$var)

#PET
pet <- metafor::rma.mv(yi = d, V = var, mods = ~ sei, random = ~ 1 | StudyID, data = dat_long)
print(clubSandwich::coef_test(pet, vcov = clubSandwich::vcovCR(pet, cluster = dat_long$StudyID, type = "CR2")))

#PEESE
peese <- metafor::rma.mv(yi = d, V = var, mods = ~ I(sei^2), random = ~ 1 | StudyID, data = dat_long)
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
                                 
utils::write.csv(ml_summary_table, "meanlevelchange.csv", row.names = FALSE)

### END OF CODE ###
