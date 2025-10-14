### START OF CODE ###

#Load libraries/ Install Packages
library("metaSEM")
library("readxl")
library("dplyr")
library("metafor")

#Read in your data
df <- read_excel("STAnalysisNEW.xlsx", sheet = "Reappraisal", na = c("NA", "N/A", "—", "-", ".", ""))

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

#Stage 2: Hypothesised path model (sequential chain T1→T2→T3→T4→T5→T6)
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

#Moderator Analaysis using Invariance

#Clean up Reports Invariance? column
df2$Invariance <- ifelse(df2$`Reports Invariance?` == "YES", 1, 0)

#Quick check
table(df2$Invariance, useNA = "ifany")

#Invariance as a matrix
inv <- as.matrix(df2[,"Invariance"])

#Stage 1 with Invariance as moderator 
re_stage1_inv <- tssem1(
  cor.analysis = TRUE,
  Cov          = cor_matrices,
  n            = n_list,
  method       = "REM",
  x            = inv 
)

summary(re_stage1_inv)

#Moderator Analysis with QA 

#Clean up QA moderator
df2$QA_mod <- ifelse(df2$Qatest %in% c("HIGH"), 1, 0)

#Quick check
table(df2$QA_mod, useNA = "ifany")

#QA as a matrix
qa <- as.matrix(df2[, "QA_mod"])

#Stage 1 with QA as moderator 
re_stage1_qa <- tssem1(
  cor.analysis = TRUE,
  Cov          = cor_matrices,
  n            = n_list,
  method       = "REM",
  x            = qa
)

summary(re_stage1_qa)

#Moderator Analysis using country

#Create USA vs Others moderator
df2$USA_mod <- ifelse(df2$Country == "USA", 1, 0)

#Quick Check
table(df2$USA_mod, useNA = "ifany")

# Convert to matrix
usa <- as.matrix(df2[, "USA_mod"])

#Stage 1 with USA as moderator
re_stage1_usa <- tssem1(
  cor.analysis = TRUE,
  Cov          = cor_matrices,
  n            = n_list,
  method       = "REM",
  x            = usa
)

summary(re_stage1_usa)

#Sensitivity Analysis

#Create columns
cor_cols <- c("rT1T2", "rT1T3", "rT2T3")

results_all <- list()

for (cc in cor_cols) {
  df_es <- df2 %>%
    filter(!is.na(.data[[cc]])) %>%
    mutate(
      FisherZ = atanh(.data[[cc]]),
      VarZ    = 1 / (SampleSize - 3)
    )
  
  if (nrow(df_es) > 2) {
    cat("\n=============================\n")
    cat("Correlation:", cc, "\n")
    
    #Random-effects model
    res <- rma(yi = FisherZ, vi = VarZ, method = "REML", data = df_es)
    print(summary(res))

    #Egger’s regression test
    egger <- regtest(res, model = "lm")
    print(egger)
    
    #Funnel plot
    funnel(res, main = paste("Funnel Plot:", cc))
    
    #Leave-one-out
    loo <- leave1out(res)
    print(loo)
    
    #Influence diagnostics
    inf <- influence(res)
    title(main = paste("Influence Diagnostics:", cc))
    
    #Trim-and-fill
    taf <- trimfill(res)
    print(summary(taf))
    funnel(taf, main = paste("Trim-and-Fill:", cc))
    
    #Save results
    results_all[[cc]] <- list(model=res, egger=egger, loo=loo, inf=inf, taf=taf)
    
  } else {
    cat("Not enough studies for:", cc, "\n")
    results_all[[cc]] <- "Not enough studies"
  }
}

summary_table <- data.frame(
  Correlation = character(),
  k = numeric(),
  Estimate = numeric(),
  SE = numeric(),
  CI_lb = numeric(),
  CI_ub = numeric(),
  Egger_p = numeric(),
  TrimFill_Est = numeric(),
  stringsAsFactors = FALSE
)

for (cc in cor_cols) {
  df_es <- df2 %>%
    filter(!is.na(.data[[cc]])) %>%
    mutate(FisherZ = atanh(.data[[cc]]),
           VarZ = 1 / (SampleSize - 3))
  
  if (nrow(df_es) > 2) {
    res <- rma(yi = FisherZ, vi = VarZ, method = "REML", data = df_es)
    egger <- regtest(res, model = "lm")
    taf <- trimfill(res)
    
    #Store summary stats
    summary_table <- rbind(summary_table, data.frame(
      Correlation = cc,
      k = res$k,
      Estimate = res$b[1],
      SE = res$se,
      CI_lb = res$ci.lb,
      CI_ub = res$ci.ub,
      Egger_p = egger$pval,
      TrimFill_Est = taf$b[1]
    ))
  }
}

print(summary_table)
write.csv(summary_table, "rankedordersummary_reppraisal.csv", row.names=FALSE)

### END OF CODE ###

