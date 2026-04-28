### START OF CODE ###

#R Version 4.51
#Load / Install Packages
library(readxl) #version 1.4.5
library(dplyr) #version 1.1.4
library(metafor) #version 4.8.0
library(lmerTest) #version 3.1.3
library(tidyr) #version 1.3.1

##Rank Order Stability: Breit-style correlation meta-analysis ----
#Read data
df <- readxl::read_xlsx(
  "STAnalysisNEW.xlsx",
  sheet = "Rumination",
  na = c("NA", "N/A", "—", "-", ".", "")
)

#Correlation columns
cor_cols <- c("rT1T2", "rT1T3", "rT2T3")

#Prepare data
df2 <- df %>%
  dplyr::mutate(
    SampleSize = suppressWarnings(as.numeric(SampleSize)),
    AgeMeanT1  = suppressWarnings(as.numeric(AgeMeanT1)),
    AgeSDT1    = suppressWarnings(as.numeric(AgeSDT1)),
    
    `Paper ID`  = as.factor(`Paper ID`),
    `Sample ID` = as.factor(`Sample ID`),
    `Meta ID`   = as.factor(`Meta ID`)
  ) %>%
  dplyr::filter(!is.na(SampleSize)) %>%
  dplyr::filter(rowSums(!is.na(dplyr::select(., dplyr::all_of(cor_cols)))) > 0)

#Reshape correlations to long format
dat_long <- df2 %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(cor_cols),
    names_to = "comparison",
    values_to = "r"
  ) %>%
  dplyr::filter(!is.na(r)) %>%
  dplyr::mutate(
    r = suppressWarnings(as.numeric(r)),
    
    comp_clean = sub("^r", "", comparison),
    t1 = as.numeric(sub("T([0-9]+)T([0-9]+)", "\\1", comp_clean)),
    t2 = as.numeric(sub("T([0-9]+)T([0-9]+)", "\\2", comp_clean)),
    
    lag = t2 - t1,
    lag_c = lag - 1,
    
    FisherZ = atanh(r),
    VarZ = 1 / (SampleSize - 3),
    SEZ = sqrt(VarZ)
  ) %>%
  dplyr::filter(
    !is.na(r),
    r > -1,
    r < 1,
    is.finite(FisherZ),
    is.finite(VarZ),
    VarZ > 0
  )

#Pooled rank-order stability
pooled_rank <- metafor::rma.mv(
  yi     = FisherZ,
  V      = VarZ,
  mods   = ~ lag_c,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long,
  method = "REML"
)

summary(pooled_rank)

#Back-transform pooled estimate to Pearson's r
tanh(coef(pooled_rank)[1])

##Moderation ----
#Mean age
dat_long <- dat_long %>%
  dplyr::mutate(
    AgeMeanT1_c = scale(AgeMeanT1, center = TRUE, scale = FALSE)[, 1]
  )

m_age <- metafor::rma.mv(
  yi     = FisherZ,
  V      = VarZ,
  mods   = ~ lag_c + AgeMeanT1_c,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long,
  method = "REML"
)

summary(m_age)

#Study quality
dat_long <- dat_long %>%
  dplyr::mutate(
    QA_mod = dplyr::case_when(
      Qatest %in% c("HIGH", "High", "high") ~ 1,
      Qatest %in% c("LOW", "Low", "low") ~ 0,
      TRUE ~ NA_real_
    )
  )

m_QA <- metafor::rma.mv(
  yi     = FisherZ,
  V      = VarZ,
  mods = ~ lag_c + QA_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(QA_mod)),
  method = "REML"
)

summary(m_QA)

#Exploratory: West vs East
dat_long <- dat_long %>%
  dplyr::mutate(
    WestEast_mod = dplyr::case_when(
      Region %in% c("America", "Europe", "Oceania") ~ 1,
      Region == "Asia" ~ 0,
      TRUE ~ NA_real_
    )
  )

m_westeast <- metafor::rma.mv(
  yi     = FisherZ,
  V      = VarZ,
  mods = ~ lag_c + WestEast_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(WestEast_mod)),
  method = "REML"
)

summary(m_westeast)

##Sensitivity Analysis / Pub Bias Tests ----
#Egger's test
dat_long <- dat_long %>%
  dplyr::mutate(
    precision = 1 / SEZ
  )

egger_rank <- lmerTest::lmer(
  I(FisherZ / SEZ) ~ 1 + precision + (1 | `Paper ID`/`Sample ID`/`Meta ID`),
  data = dat_long
)

summary(egger_rank, correlation = FALSE)

#Leave-one-out by Paper ID
paper_ids <- unique(dat_long$`Paper ID`)

loo_results <- data.frame(
  PaperID = paper_ids,
  k       = NA,
  est_z   = NA,
  est_r   = NA,
  se      = NA,
  zval    = NA,
  pval    = NA
)

for (i in seq_along(paper_ids)) {
  dat_loo <- subset(dat_long, `Paper ID` != paper_ids[i])
  
  fit_loo <- try(
    metafor::rma.mv(
      yi     = FisherZ,
      V      = VarZ,
      mods   = ~ lag_c,
      random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
      data   = dat_loo,
      method = "REML"
    ),
    silent = TRUE
  )
  
  if (!inherits(fit_loo, "try-error")) {
    loo_results$k[i]     <- nrow(dat_loo)
    loo_results$est_z[i] <- fit_loo$b[1, 1]
    loo_results$est_r[i] <- tanh(fit_loo$b[1, 1])
    loo_results$se[i]    <- fit_loo$se[1]
    loo_results$zval[i]  <- fit_loo$zval[1]
    loo_results$pval[i]  <- fit_loo$pval[1]
  }
}

full_est_z <- pooled_rank$b[1, 1]

loo_results$diff_z <- loo_results$est_z - full_est_z
loo_results$diff_r <- loo_results$est_r - tanh(full_est_z)

loo_results_sorted <- loo_results[order(abs(loo_results$diff_z), decreasing = TRUE), ]

print(loo_results_sorted)

## -----
### END OF CODE ###



