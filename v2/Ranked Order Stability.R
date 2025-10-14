### START OF CODE ###

#Load libraries/ Install Packages
library(metaSEM)
library(readxl)
library(dplyr)
library(metafor)
library(clubSandwich)

#Read in your data
df <- read_xlsx("STAnalysisNEW.xlsx", sheet = "Reappraisal", na = c("NA", "N/A", "—", "-", ".", ""))

#Filter out rows missing SampleSize or if missing T1T2
df2 <- df %>%
  filter(!is.na(SampleSize)) %>%
  filter(rowSums(!is.na(select(., starts_with("rT")))) > 0)

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
re_stage1 <- tssem1(
  cor.analysis = TRUE,
  Cov          = cor_matrices,
  n            = n_list,
  method       = "REM"
)

summary(re_stage1)

#Stage 2: Hypothesised path model
hypothesised_path_model <- lavaan2RAM(
  'T3 ~ T2
  T2 ~ T1
  T1 ~~ 1*T1',
  obs.variables = all_tps
)

#Fit the model
re_stage2 <- tssem2(
  re_stage1,
  RAM = hypothesised_path_model
)
summary(re_stage2)

#Moderator Analysis using Invariance
#Clean up Reports Invariance? column
df2$Invariance <- ifelse(df2$`Reports Invariance?` == "YES", 1, 0)

#Quick check
table(df2$Invariance, useNA = "ifany")

#Invariance as a matrix
inv <- as.matrix(df2[,"Invariance"])

#Stage 1 with Invariance as moderator 
re_stage1_inv <- tryCatch({
  tssem1(
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
#Clean up  AchievedFullScalarInvariance column
df2$Fullinv <- ifelse(df2$AchievedFullScalarInvariance == "YES", 1, 0)

#Quick check
table(df2$Fullinv, useNA = "ifany")

#AchievedFullScalarInvariance
fullinv <- as.matrix(df2[, "Fullinv"])

#Stage 1 with AchievedFullScalarInvariance as moderator
re_stage1_scalar_inv <- tryCatch({
  tssem1(
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
#Clean up QA moderator
df2$QA_mod <- ifelse(df2$Qatest %in% c("HIGH"), 1, 0)

#Quick check
table(df2$QA_mod, useNA = "ifany")

#QA as a matrix
qa <- as.matrix(df2[, "QA_mod"])

#Stage 1 with QA as moderator
re_stage1_qa <- tryCatch({
  tssem1(
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
#Create West v East moderator
df2$WestEast <- ifelse(df2$Region %in% c("America", "Europe", "Oceania"), 1, 0)

#Quick Check
table(df2$WestEast, useNA = "ifany")

#West v east as a matrix
west_east <- as.matrix(df2[, "WestEast"])

#Stage 1 with West v East as moderator  
re_stage1_west_east <- tryCatch({
  tssem1(
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
# Define the correlation columns
cor_cols <- c("rT1T2", "rT1T3", "rT2T3")

# Ensure StudyID is defined
df2$StudyID <- paste(df2$Author, df2$`Year Collected`, sep = "_")

# Store results
results_all <- list()

for (cc in cor_cols) {
  df_es <- df2 %>%
    filter(!is.na(.data[[cc]])) %>%
    mutate(
      FisherZ = atanh(.data[[cc]]),
      VarZ = 1 / (SampleSize - 3),
      Comparison = cc
    )
  
  if (nrow(df_es) > 2) {
    cat("\n=============================\n")
    cat("Testing correlation:", cc, "\n")
    
    #Multilevel model
    res_mv <- rma.mv(
      yi = FisherZ,
      V = VarZ,
      random = ~ 1 | StudyID/Comparison,
      data = df_es,
      method = "REML"
    )
    print(summary(res_mv))
    
    #Funnel Plot Inspection for Pub Bias Test
    funnel(res_mv, main = "Funnel Plot (Multilevel Model)")
    
    #Cluster-robust SEs robustness check
    cr2_res <- coef_test(res_mv, vcov = "CR2", cluster = df_es$StudyID)
    print(cr2_res)
    
    #Manual leave-one-out robustness check
    loo_df <- do.call(rbind, lapply(unique(df_es$StudyID), function(id) {
      dat_loo <- df_es[df_es$StudyID != id, ]
      fit_loo <- try(
        rma.mv(yi = FisherZ, V = VarZ, random = ~ 1 | StudyID/Comparison, data = dat_loo, method = "REML"),
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

# Build summary table from all results
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
write.csv(summary_table, "rankedorder.csv", row.names = FALSE)


### END OF CODE ###

