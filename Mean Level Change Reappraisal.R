#Do What? 
1. Beckers d
2. Variance of beckers d
3. metafor rma #skip escalc cus alr doing in 1 and 2
4. Moderator analysis w metafor here #eventually doing with both metasem and mean level change

### START OF CODE ###

#Packages
library(readxl)
library(metafor)
library(dplyr)
library(tidyr)
library(stringr)

#Load and clean data, and calculate beckers d and intermediate step
dat <- 
  read_xlsx("STAnalysisNew.xlsx", sheet = "Reappraisal") %>% 
  mutate(
    MeanT1 = as.numeric(MeanT1), 
    MeanT2 = as.numeric(MeanT2),
    MeanT3 = as.numeric(MeanT3),
    MeanT4 = as.numeric(MeanT4),
    MeanT5 = as.numeric(MeanT5),
    MeanT6 = as.numeric(MeanT6),
    SDT1 = as.numeric(SDT1),
    SDT2 = as.numeric(SDT2),
    SDT3 = as.numeric(SDT3),
    SDT4 = as.numeric(SDT4),
    SDT5 = as.numeric(SDT5),
    SDT6 = as.numeric(SDT6),
    SampleSize = as.numeric(SampleSize),
    beckersdT1T2 = (MeanT2-MeanT1)/SDT1, #beckersd=+-2 is not realistic. Check DE again
    beckersdT1T3 = (MeanT3-MeanT1)/SDT1,
    beckersdT1T4 = (MeanT4-MeanT1)/SDT1,
    beckersdT1T5 = (MeanT5-MeanT1)/SDT1,
    beckersdT1T6 = (MeanT6-MeanT1)/SDT1,
    beckersdT2T3 = (MeanT3-MeanT2)/SDT2,
    beckersdT2T4 = (MeanT4-MeanT2)/SDT2,
    beckersdT2T5 = (MeanT5-MeanT2)/SDT2,
    beckersdT2T6 = (MeanT6-MeanT2)/SDT2,
    beckersdT3T4 = (MeanT4-MeanT3)/SDT3,
    beckersdT3T5 = (MeanT5-MeanT3)/SDT3,
    beckersdT3T6 = (MeanT6-MeanT3)/SDT3,
    beckersdT4T5 = (MeanT5-MeanT4)/SDT4,
    beckersdT4T6 = (MeanT6-MeanT4)/SDT4,
    beckersdT5T6 = (MeanT6-MeanT5)/SDT5,
    intermediate = 1 - (3/(4*(SampleSize-1)-1)))

#calculate beckers d variance
dat <- dat %>% mutate(
  varbeckersdT1T2 = (((1+SampleSize*beckersdT1T2^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     -(beckersdT1T2^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT1T3 = (((1 + SampleSize*beckersdT1T3^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT1T3^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT1T4 = (((1 + SampleSize*beckersdT1T4^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT1T4^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT1T5 = (((1 + SampleSize*beckersdT1T5^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT1T5^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT1T6 = (((1 + SampleSize*beckersdT1T6^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT1T6^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT2T3 = (((1 + SampleSize*beckersdT2T3^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT2T3^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT2T4 = (((1 + SampleSize*beckersdT2T4^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT2T4^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT2T5 = (((1 + SampleSize*beckersdT2T5^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT2T5^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT2T6 = (((1 + SampleSize*beckersdT2T6^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT2T6^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT3T4 = (((1 + SampleSize*beckersdT3T4^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT3T4^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT3T5 = (((1 + SampleSize*beckersdT3T5^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT3T5^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT3T6 = (((1 + SampleSize*beckersdT3T6^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT3T6^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT4T5 = (((1 + SampleSize*beckersdT4T5^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT4T5^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT4T6 = (((1 + SampleSize*beckersdT4T6^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT4T6^2/(intermediate*(SampleSize-1))^2)),
  varbeckersdT5T6 = (((1 + SampleSize*beckersdT5T6^2)/SampleSize)*((SampleSize-1)/(SampleSize-3))
                     - (beckersdT5T6^2/(intermediate*(SampleSize-1))^2))
) %>% filter(is.na(beckersdT1T2) | (beckersdT1T2 > -2 & beckersdT1T2 <  2),
             is.na(beckersdT1T3) | (beckersdT1T3 > -2 & beckersdT1T3 <  2),
             is.na(beckersdT1T4) | (beckersdT1T4 > -2 & beckersdT1T4 <  2),
             is.na(beckersdT1T5) | (beckersdT1T5 > -2 & beckersdT1T5 <  2),
             is.na(beckersdT1T6) | (beckersdT1T6 > -2 & beckersdT1T6 <  2),
             is.na(beckersdT2T3) | (beckersdT2T3 > -2 & beckersdT2T3 <  2),
             is.na(beckersdT2T4) | (beckersdT2T4 > -2 & beckersdT2T4 <  2),
             is.na(beckersdT2T5) | (beckersdT2T5 > -2 & beckersdT2T5 <  2),
             is.na(beckersdT2T6) | (beckersdT2T6 > -2 & beckersdT2T6 <  2),
             is.na(beckersdT3T4) | (beckersdT3T4 > -2 & beckersdT3T4 <  2),
             is.na(beckersdT3T5) | (beckersdT3T5 > -2 & beckersdT3T5 <  2),
             is.na(beckersdT3T6) | (beckersdT3T6 > -2 & beckersdT3T6 <  2),
             is.na(beckersdT4T5) | (beckersdT4T5 > -2 & beckersdT4T5 <  2),
             is.na(beckersdT4T6) | (beckersdT4T6 > -2 & beckersdT4T6 <  2),
             is.na(beckersdT5T6) | (beckersdT5T6 > -2 & beckersdT5T6 <  2))

#T1–T2
dat2_T1T2 <- metafor::rma(yi = beckersdT1T2, vi = varbeckersdT1T2, data = dat)
summary(dat2_T1T2)

#T1–T3
dat2_T1T3 <- metafor::rma(yi = beckersdT1T3, vi = varbeckersdT1T3, data = dat)
summary(dat2_T1T3)

#T1–T4
dat2_T1T4 <- metafor::rma(yi = beckersdT1T4, vi = varbeckersdT1T4, data = dat)
summary(dat2_T1T4)

#T1–T5
dat2_T1T5 <- metafor::rma(yi = beckersdT1T5, vi = varbeckersdT1T5, data = dat)
summary(dat2_T1T5)

#T1–T6
dat2_T1T6 <- metafor::rma(yi = beckersdT1T6, vi = varbeckersdT1T6, data = dat)
summary(dat2_T1T6)

#T2–T3
dat2_T2T3 <- metafor::rma(yi = beckersdT2T3, vi = varbeckersdT2T3, data = dat)
summary(dat2_T2T3)

#T2–T4
dat2_T2T4 <- metafor::rma(yi = beckersdT2T4, vi = varbeckersdT2T4, data = dat)
summary(dat2_T2T4)

#T2–T5
dat2_T2T5 <- metafor::rma(yi = beckersdT2T5, vi = varbeckersdT2T5, data = dat)
summary(dat2_T2T5)

#T2–T6
dat2_T2T6 <- metafor::rma(yi = beckersdT2T6, vi = varbeckersdT2T6, data = dat)
summary(dat2_T2T6)

#T3–T4
dat2_T3T4 <- metafor::rma(yi = beckersdT3T4, vi = varbeckersdT3T4, data = dat)
summary(dat2_T3T4)

#T3–T5
dat2_T3T5 <- metafor::rma(yi = beckersdT3T5, vi = varbeckersdT3T5, data = dat)
summary(dat2_T3T5)

#T3–T6
dat2_T3T6 <- metafor::rma(yi = beckersdT3T6, vi = varbeckersdT3T6, data = dat)
summary(dat2_T3T6)

#T4–T5
dat2_T4T5 <- metafor::rma(yi = beckersdT4T5, vi = varbeckersdT4T5, data = dat)
summary(dat2_T4T5)

#T4–T6
dat2_T4T6 <- metafor::rma(yi = beckersdT4T6, vi = varbeckersdT4T6, data = dat)
summary(dat2_T4T6)

#T5–T6
dat2_T5T6 <- metafor::rma(yi = beckersdT5T6, vi = varbeckersdT5T6, data = dat)
summary(dat2_T5T6)
                          
#Sequential Change for Meta-regression

#Gather all Becker's d and variance columns into long format
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
      comparison == "beckersdT1T4" ~ varbeckersdT1T4,
      comparison == "beckersdT1T5" ~ varbeckersdT1T5,
      comparison == "beckersdT1T6" ~ varbeckersdT1T6,
      comparison == "beckersdT2T3" ~ varbeckersdT2T3,
      comparison == "beckersdT2T4" ~ varbeckersdT2T4,
      comparison == "beckersdT2T5" ~ varbeckersdT2T5,
      comparison == "beckersdT2T6" ~ varbeckersdT2T6,
      comparison == "beckersdT3T4" ~ varbeckersdT3T4,
      comparison == "beckersdT3T5" ~ varbeckersdT3T5,
      comparison == "beckersdT3T6" ~ varbeckersdT3T6,
      comparison == "beckersdT4T5" ~ varbeckersdT4T5,
      comparison == "beckersdT4T6" ~ varbeckersdT4T6,
      comparison == "beckersdT5T6" ~ varbeckersdT5T6,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(d), !is.na(var)) %>%
  mutate(
    #Extract numeric timepoints from comparison names
    t1 = as.numeric(str_extract(comparison, "(?<=T)\\d(?=T)")),
    t2 = as.numeric(str_extract(comparison, "(?<=T\\d)T\\d")),
    lag = t2 - t1
  )

#Check structure
head(dat_long)

#Clean comparison labels and extract timepoints
dat_long <- dat_long %>%
  mutate(
    comp_clean = str_remove(comparison, "^beckersd"),  # remove prefix
    t1 = as.numeric(str_sub(comp_clean, 2, 2)),
    t2 = as.numeric(str_sub(comp_clean, 4, 4)),
    lag = t2 - t1
  )

#Quick check
dat_long %>% select(comparison, comp_clean, d, var, t1, t2, lag) %>% head(10)

#Meta-regression: Becker's d ~ time lag
res_lag <- rma(yi = d, vi = var, mods = ~ lag, method = "REML", data = dat_long)
summary(res_lag)

#Predicted values across lag lengths 1–5
pred_lag <- predict(res_lag, newmods = 1:5)
print(pred_lag)

#Plot predictions with 95% CI
plot(1:5, pred_lag$pred, type = "b", pch = 19,
     ylim = range(pred_lag$ci.lb, pred_lag$ci.ub),
     xlab = "Lag Length (Timepoints Apart)",
     ylab = "Pooled Becker's d",
     main = "Predicted Mean-Level Change by Lag")
arrows(1:5, pred_lag$ci.lb, 1:5, pred_lag$ci.ub,
       angle = 90, code = 3, length = 0.05)

#Moderator Analysis for Sequential Change

#Make sure moderators are in dat_long (join them from original dat — which has QA, Invariance, USA_mod)
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
      comparison == "beckersdT1T4" ~ varbeckersdT1T4,
      comparison == "beckersdT1T5" ~ varbeckersdT1T5,
      comparison == "beckersdT1T6" ~ varbeckersdT1T6,
      comparison == "beckersdT2T3" ~ varbeckersdT2T3,
      comparison == "beckersdT2T4" ~ varbeckersdT2T4,
      comparison == "beckersdT2T5" ~ varbeckersdT2T5,
      comparison == "beckersdT2T6" ~ varbeckersdT2T6,
      comparison == "beckersdT3T4" ~ varbeckersdT3T4,
      comparison == "beckersdT3T5" ~ varbeckersdT3T5,
      comparison == "beckersdT3T6" ~ varbeckersdT3T6,
      comparison == "beckersdT4T5" ~ varbeckersdT4T5,
      comparison == "beckersdT4T6" ~ varbeckersdT4T6,
      comparison == "beckersdT5T6" ~ varbeckersdT5T6,
      TRUE ~ NA_real_
    ),
    #Extract numbers safely
    t1 = as.numeric(str_match(comparison, "T(\\d)T(\\d)")[,2]),
    t2 = as.numeric(str_match(comparison, "T(\\d)T(\\d)")[,3]),
    lag = t2 - t1,
    #Encode moderators directly
    QA_mod  = ifelse(Qatest == "High", 1, 0),
    Inv_mod = ifelse(`Reports Invariance?` == "YES", 1, 0),
    USA_mod = ifelse(Country == "USA", 1, 0)
  ) %>%
  filter(!is.na(d), !is.na(var), !is.na(lag))

#QA
res_lag_QA <- rma(yi = d, vi = var, mods = ~ lag * QA_mod, method = "REML", data = dat_long)
summary(res_lag_QA)

#Inv
res_lag_INV <- rma(yi = d, vi = var, mods = ~ lag * Inv_mod, method = "REML", data = dat_long)
summary(res_lag_INV)

#USA v Others
res_lag_USA <- rma(yi = d, vi = var, mods = ~ lag * USA_mod, method = "REML", data = dat_long)
summary(res_lag_USA)

#Full Sensitivity and Bias Tests
  
#All pairwise comparisons
pairs <- c("T1T2","T1T3","T1T4","T1T5","T1T6",
           "T2T3","T2T4","T2T5","T2T6",
           "T3T4","T3T5","T3T6",
           "T4T5","T4T6",
           "T5T6")

#Storage
results <- list()
summary_table <- data.frame(
  Pair = character(),
  k = integer(),
  FEM_est = numeric(), FEM_ciLB = numeric(), FEM_ciUB = numeric(),
  REM_est = numeric(), REM_ciLB = numeric(), REM_ciUB = numeric(),
  I2 = numeric(),
  Egger_p = numeric(),
  TrimFill_k = integer(), TrimFill_est = numeric(),
  stringsAsFactors = FALSE
)

#Loop
for (p in pairs) {
  d_col <- paste0("beckersd", p)
  v_col <- paste0("varbeckersd", p)
  
  dat_p <- dat %>% filter(!is.na(.data[[d_col]]), !is.na(.data[[v_col]]))
  k <- nrow(dat_p)
  
  if (k >= 2) {
    cat("\n=============================\n")
    cat("Comparison:", p, " (k =", k, ")\n")
    
    #FEM v REM
    res_fem <- tryCatch(rma(yi = dat_p[[d_col]], vi = dat_p[[v_col]], method="FE"), error=function(e) NULL)
    res_rem <- tryCatch(rma(yi = dat_p[[d_col]], vi = dat_p[[v_col]], method="REML"), error=function(e) NULL)
    
    #Egger’s test
    egger_p <- tryCatch({
      reg <- regtest(res_rem, model="lm")
      reg$pval
    }, error=function(e) NA)
    
    #Trim and fill
    trimfill_k <- NA; trimfill_est <- NA
    taf <- tryCatch({
      tf <- trimfill(res_rem)
      trimfill_k <- tf$k
      trimfill_est <- tf$b
      tf
    }, error=function(e) NULL)
    
    #Sensitivity Analysis, Leave-one-out LOO and Influence inf
    loo <- tryCatch(leave1out(res_rem), error=function(e) NULL)
    inf <- tryCatch(influence(res_rem), error=function(e) NULL)
    
    #FEMvREM, Egger's test, Trim and fill, LOO and inf AGAIN but only if k >= 5
    if (k >= 5 && !is.null(res_rem)) {
      funnel(res_rem, main = paste("Funnel:", p))
      if (!is.null(taf)) funnel(taf, main = paste("Trim-and-Fill:", p))
      forest(res_rem,
             slab = paste(dat_p$Author, dat_p$`Year Collected`),
             xlab = "Becker's d", main = paste("Forest Plot:", p))
      plot(inf); title(main = paste("Influence Diagnostics:", p))
    }
    
    #Save detailed results
    results[[p]] <- list(
      FEM=res_fem, REM=res_rem, Egger_p=egger_p,
      TrimFill=taf, LOO=loo, Influence=inf
    )
    
    #Save summary row
    summary_table <- rbind(summary_table, data.frame(
      Pair = p,
      k = k,
      FEM_est = if (!is.null(res_fem)) res_fem$b else NA,
      FEM_ciLB = if (!is.null(res_fem)) res_fem$ci.lb else NA,
      FEM_ciUB = if (!is.null(res_fem)) res_fem$ci.ub else NA,
      REM_est = if (!is.null(res_rem)) res_rem$b else NA,
      REM_ciLB = if (!is.null(res_rem)) res_rem$ci.lb else NA,
      REM_ciUB = if (!is.null(res_rem)) res_rem$ci.ub else NA,
      I2 = if (!is.null(res_rem)) res_rem$I2 else NA,
      Egger_p = egger_p,
      TrimFill_k = trimfill_k,
      TrimFill_est = trimfill_est
    ))
    
  } else {
    cat("Skipping", p, "(not enough studies)\n")
    results[[p]] <- "Not enough studies"
    summary_table <- rbind(summary_table, data.frame(
      Pair=p, k=k, FEM_est=NA, FEM_ciLB=NA, FEM_ciUB=NA,
      REM_est=NA, REM_ciLB=NA, REM_ciUB=NA, I2=NA,
      Egger_p=NA, TrimFill_k=NA, TrimFill_est=NA
    ))
  }
}

#Clean output
print(summary_table)
write.csv(summary_table, "meta_sensitivity_summary_reappraisal.csv", row.names=FALSE)

### END OF CODE ###
