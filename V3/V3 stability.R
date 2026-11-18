### START OF CODE ###

#R Version 4.5.1
#Load / Install Packages
library(readxl)   #version 1.4.5; Wickham & Bryan (2025)
library(dplyr)    #version 1.1.4; Wickham et al. (2023)
library(metafor)  #version 4.8.0; Viechtbauer (2010)
library(lmerTest) #version 3.1.3; Kuznetsova et al. (2017)
library(tidyr)    #version 1.3.1; Wickham et al. (2024)

##Rank Order Stability: Breit-style correlation meta-analysis ----
#Read data; Wickham & Bryan (2025)
df <- readxl::read_xlsx(
  "STAnalysis7.xlsx",
  sheet = "Suppression",
  na = c("NA", "N/A", "—", "-", ".", "")
)

#Correlation columns
cor_cols <- c("rT1T2", "rT1T3", "rT2T3")

#Prepare data; Wickham et al. (2023)
df2 <- df %>%
  dplyr::mutate(
    SampleSize = suppressWarnings(as.numeric(SampleSize)),
    AgeMeanT1  = suppressWarnings(as.numeric(AgeMeanT1)),
    AgeSDT1    = suppressWarnings(as.numeric(AgeSDT1)),
    
    #prep nesting variables by making them into factors; Konstantopoulos (2011); Viechtbauer (2010)
    `Paper ID`  = as.factor(`Paper ID`),
    `Sample ID` = as.factor(`Sample ID`),
    `Meta ID`   = as.factor(`Meta ID`)
  ) %>%
  dplyr::filter(!is.na(SampleSize)) %>%
  dplyr::filter(rowSums(!is.na(dplyr::select(., dplyr::all_of(cor_cols)))) > 0)

#Reshape correlations to long format; Wickham et al. (2024); Konstantopoulos (2011); Viechtbauer (2010)
dat_long <- df2 %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(cor_cols),
    names_to = "comparison",
    values_to = "r"
  ) %>%
  dplyr::filter(!is.na(r)) %>%
  dplyr::mutate(
    r = suppressWarnings(as.numeric(r)),
    
    #Clean comparison label and extract timepoints; Wickham et al. (2023)
    comp_clean = sub("^r", "", comparison),
    t1 = as.numeric(sub("T([0-9]+)T([0-9]+)", "\\1", comp_clean)),
    t2 = as.numeric(sub("T([0-9]+)T([0-9]+)", "\\2", comp_clean)),
    
    #Time lag centred at one year for interpretation of intercept; Viechtbauer (2010); Thompson & Higgins (2002)
    lag = t2 - t1,
    lag_c = lag - 1,
    
    #Fisher's z transformation and sampling variance for correlations; Borenstein et al. (2009); Cox (2008)
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

#Pooled rank-order stability; Konstantopoulos (2011); Viechtbauer (2010)
pooled_rank <- metafor::rma.mv(
  yi     = FisherZ,
  V      = VarZ,
  mods   = ~ lag_c,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long,
  method = "REML"
)

summary(pooled_rank)

#Back-transform pooled Fisher's z estimate to Pearson's r; Cox (2008)
tanh(coef(pooled_rank)[1])

##Moderation ----
#Mean age moderator; Thompson & Higgins (2002); Viechtbauer (2010)
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

#Study quality moderator; Thompson & Higgins (2002); Viechtbauer (2010)
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
  mods   = ~ lag_c + QA_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(QA_mod)),
  method = "REML"
)

summary(m_QA)

#Exploratory country moderator; Thompson & Higgins (2002); Viechtbauer (2010)
dat_long <- dat_long %>%
  dplyr::mutate(
    Country_mod = trimws(as.character(Country)),
    Country_mod = dplyr::na_if(Country_mod, "NA"),
    Country_mod = as.factor(Country_mod)
  )

m_country <- metafor::rma.mv(
  yi     = FisherZ,
  V      = VarZ,
  mods   = ~ lag_c + Country_mod,
  random = ~ 1 | `Paper ID`/`Sample ID`/`Meta ID`,
  data   = dat_long %>% dplyr::filter(!is.na(Country_mod)),
  method = "REML"
)

summary(m_country)

##Sensitivity Analysis / Pub Bias Tests ----
#Egger's test for small-study effects; Egger et al. (1997); Sterne & Egger (2001); Kuznetsova et al. (2017)
dat_long <- dat_long %>%
  dplyr::mutate(
    precision = 1 / SEZ
  )

egger_rank <- lmerTest::lmer(
  I(FisherZ / SEZ) ~ 1 + precision + (1 | `Paper ID`/`Sample ID`/`Meta ID`),
  data = dat_long
)

summary(egger_rank, correlation = FALSE)

#Leave-one-out sensitivity analysis by Paper ID; Viechtbauer & Cheung (2010)
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

for (i in seq_along(paper_ids))

full_est_z <- pooled_rank$b[1, 1]

loo_results$diff_z <- loo_results$est_z - full_est_z

#Back-transform full-model Fisher's z estimate to Pearson's r before computing r-scale difference; Cox (2008)
loo_results$diff_r <- loo_results$est_r - tanh(full_est_z)

loo_results_sorted <- loo_results[order(abs(loo_results$diff_z), decreasing = TRUE), ]

print(loo_results_sorted)

### END OF CODE ###

packages <- c(
  "readxl",
  "dplyr",
  "metafor",
  "lmerTest",
  "tidyr"
)

data.frame(
  Package = packages,
  Version = sapply(packages, function(x) {
    as.character(packageVersion(x))
  }),
  row.names = NULL
)