### START OF CODE ###

#R Version 4.51
#Load libraries/ Install Packages
library(metaSEM) #version 1.5.0
library(readxl) #version 1.4.5
library(dplyr) #version 1.1.4
library(metafor) #version 4.8.0
library(clubSandwich) #version 0.6.1

### START OF CODE ###

#Load libraries / Install Packages
library(metaSEM)
library(readxl)
library(dplyr)
library(metafor)
library(clubSandwich)

#Read in your data
df <- readxl::read_xlsx("STAnalysisNEW.xlsx", sheet = "Reappraisal", na = c("NA", "N/A", "—", "-", ".", ""))

#Filter out rows missing SampleSize or if missing T1T2
df2 <- df %>%
  dplyr::filter(!is.na(SampleSize)) %>%
  dplyr::filter(rowSums(!is.na(dplyr::select(., dplyr::starts_with("rT")))) > 0)

#Restrict to correlation columns for T1–T6
cor_cols <- c("rT1T2", "rT1T3", "rT2T3")

#Timepoints we care about
all_tps <- paste0("T", 1:3)

#Function to build a 3x3 correlation matrix for one row
row_to_cor_matrix <- function(row) {
  mat <- diag(1, length(all_tps))
  dimnames(mat) <- list(all_tps, all_tps)
  
  for (p in cor_cols) {
    val <- row[[p]]
    if (!is.na(val)) {
      matches <- regmatches(p, gregexpr("T[0-9]+", p))[[1]]
      t1 <- matches[1]
      t2 <- matches[2]
      mat[t1, t2] <- val
      mat[t2, t1] <- val
    }
  }
  mat
}

#Build the list of correlation matrices
cor_matrices <- lapply(
  seq_len(nrow(df2)),
  function(i) row_to_cor_matrix(df2[i, , drop = TRUE])
)

#Extract the list of sample sizes
n_list <- df2$SampleSize

#Stage 1: pool study-specific correlation matrices into one (REM)
re_stage1 <- metaSEM::tssem1(
  cor.analysis = TRUE,
  Cov          = cor_matrices,
  n            = n_list,
  method       = "REM"
)

summary(re_stage1)

#Stage 2: Hypothesised path model
hypothesised_path_model <- metaSEM::lavaan2RAM(
  'T3 ~ T2
  T2 ~ T1
  T1 ~~ 1*T1',
  obs.variables = all_tps
)

#Fit the model
re_stage2 <- metaSEM::tssem2(
  re_stage1,
  RAM = hypothesised_path_model
)
summary(re_stage2)

#Moderator Analysis using Invariance
#Clean up Reports Invariance? column
df2 <- df2 %>%
  dplyr::mutate(
    Invariance = dplyr::case_when(
      `Reports Invariance?` == "YES" ~ 1,
      TRUE ~ 0
    )
  )

#Quick check
table(df2$Invariance, useNA = "ifany")

#Invariance as a matrix
inv <- as.matrix(df2[,"Invariance"])

#Stage 1 with Invariance as moderator 
re_stage1_inv <- tryCatch({
  metaSEM::tssem1(
    cor.analysis = TRUE,
    Cov          = cor_matrices,
    n            = n_list,
    method       = "REM",
    x            = inv 
  )
}, error = function(e) {
  message("Invariance model failed: ", e$message)
  NULL
})

if (!is.null(re_stage1_inv)) summary(re_stage1_inv)

#Moderator Analysis using AchievedFullScalarInvariance 
df2 <- df2 %>%
  dplyr::mutate(
    Fullinv = dplyr::case_when(
      AchievedFullScalarInvariance == "YES" ~ 1,
      TRUE ~ 0
    )
  )

#Quick check
table(df2$Fullinv, useNA = "ifany")

#AchievedFullScalarInvariance
fullinv <- as.matrix(df2[, "Fullinv"])

#Stage 1 with AchievedFullScalarInvariance as moderator
re_stage1_scalar_inv <- tryCatch({
  metaSEM::tssem1(
    cor.analysis = TRUE,
    Cov          = cor_matrices,
    n            = n_list,
    method       = "REM",
    x            = fullinv 
  )
}, error = function(e) {
  message("Full Scalar Invariance model failed: ", e$message)
  NULL
})

if (!is.null(re_stage1_scalar_inv)) summary(re_stage1_scalar_inv)

#Moderator Analysis with QA 
df2 <- df2 %>%
  dplyr::mutate(
    QA_mod = dplyr::case_when(
      Qatest %in% c("HIGH") ~ 1,
      TRUE ~ 0
    )
  )

#Quick check
table(df2$QA_mod, useNA = "ifany")

#QA as a matrix
qa <- as.matrix(df2[, "QA_mod"])

#Stage 1 with QA as moderator
re_stage1_qa <- tryCatch({
  metaSEM::tssem1(
    cor.analysis = TRUE,
    Cov          = cor_matrices,
    n            = n_list,
    method       = "REM",
    x            = qa
  )
}, error = function(e) {
  message("QA model failed: ", e$message)
  NULL
})

if (!is.null(re_stage1_qa)) summary(re_stage1_qa)

#Moderator Analysis using West v East
df2 <- df2 %>%
  dplyr::mutate(
    WestEast = dplyr::case_when(
      Region %in% c("America", "Europe", "Oceania") ~ 1,
      TRUE ~ 0
    )
  )

#Quick Check
table(df2$WestEast, useNA = "ifany")

#West v east as a matrix
west_east <- as.matrix(df2[, "WestEast"])

#Stage 1 with West v East as moderator  
re_stage1_west_east <- tryCatch({
  metaSEM::tssem1(
    cor.analysis = TRUE,
    Cov          = cor_matrices,
    n            = n_list,
    method       = "REM",
    x            = west_east
  )
}, error = function(e) {
  message("West vs East model failed: ", e$message)
  NULL
})

if (!is.null(re_stage1_west_east)) summary(re_stage1_west_east)

#Sensitivity Analysis
cor_cols <- c("rT1T2", "rT1T3", "rT2T3")

#Ensure StudyID is defined
df2 <- df2 %>%
  dplyr::mutate(StudyID = paste(Author, `Year Collected`, sep = "_"))

#Store results
results_all <- list()

for (cc in cor_cols) {
  df_es <- df2 %>%
    dplyr::filter(!is.na(.data[[cc]])) %>%
    dplyr::mutate(
      FisherZ = atanh(.data[[cc]]),
      VarZ = 1 / (SampleSize - 3),
      Comparison = cc
    )
  
  if (nrow(df_es) > 2) {
    #Multilevel model
    res_mv <- metafor::rma.mv(
      yi = FisherZ,
      V = VarZ,
      random = ~ 1 | StudyID/Comparison,
      data = df_es,
      method = "REML"
    )
    print(summary(res_mv))
    
    #Funnel Plot Inspection for Pub Bias Test
    metafor::funnel(res_mv, main = "Funnel Plot (Multilevel Model)")
    
    #Cluster-robust SEs robustness check
    cr2_res <- clubSandwich::coef_test(res_mv, vcov = "CR2", cluster = df_es$StudyID)
    print(cr2_res)
    
    #Manual leave-one-out robustness check
    loo_df <- do.call(rbind, lapply(unique(df_es$StudyID), function(id) {
      dat_loo <- df_es[df_es$StudyID != id, ]
      fit_loo <- try(
        metafor::rma.mv(yi = FisherZ, V = VarZ, random = ~ 1 | StudyID/Comparison, data = dat_loo, method = "REML"),
        silent = TRUE
      )
      if (!inherits(fit_loo, "try-error")) {
        return(data.frame(LeftOut = id, Beta = coef(fit_loo)[1]))
      } else {
        return(data.frame(LeftOut = id, Beta = NA))
      }
    }))
    print(loo_df)
    
    #Save results
    results_all[[cc]] <- list(model = res_mv, cr2 = cr2_res, loo = loo_df)
    
  } else {
    cat("Not enough studies for:", cc, "\n")
    results_all[[cc]] <- "Not enough studies"
  }
}

#Build summary table from all results
summary_table <- data.frame(
  Correlation = character(),
  k = numeric(),
  Estimate = numeric(),
  SE = numeric(),
  CI_lb = numeric(),
  CI_ub = numeric(),
  stringsAsFactors = FALSE
)

for (cc in cor_cols) {
  if (!is.character(results_all[[cc]])) {
    mod <- results_all[[cc]]$model
    summary_table <- rbind(summary_table, data.frame(
      Correlation = cc,
      k = mod$k,
      Estimate = mod$b[1],
      SE = mod$se,
      CI_lb = mod$ci.lb,
      CI_ub = mod$ci.ub
    ))
  }
}

print(summary_table)

utils::write.csv(summary_table, "rankedorder.csv", row.names = FALSE)

### END OF CODE ###
