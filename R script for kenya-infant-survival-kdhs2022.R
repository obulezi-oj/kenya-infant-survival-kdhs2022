# ==============================================================================
# TITLE: Born Unequal: Structural Stratification, Institutional Exclusion, and
#        Infant Survival Across Kenyan Counties
# DATA:  2022 Kenya Demographic and Health Survey (KDHS) - Births Recode (KEBR8CFL)
# ==============================================================================

# 0. SET WORKING DIRECTORY & OUTPUT FOLDER -------------------------------------
work_dir <- "C:/Users/Dr. O. J. Obulezi/Documents/R projects/Dr Adimchi"
setwd(work_dir)

out_dir <- file.path(work_dir, "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir)


# 1. INSTALL & LOAD REQUIRED PACKAGES ------------------------------------------
required_pkgs <- c("haven", "tidyverse", "survey", "srvyr", "DHS.rates", 
                   "broom", "scales", "forcats", "glmnet", "ranger", "xgboost")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(haven)      # Read DHS .DTA / .SAV files
library(tidyverse)  # Data manipulation & ggplot2
library(survey)     # Complex survey design & GLM
library(srvyr)      # Tidyverse survey wrapper
library(DHS.rates)  # Official DHS synthetic-cohort mortality rates
library(broom)      # Tidy regression outputs
library(scales)     # Axis & label formatting
library(forcats)    # Factor ordering
library(glmnet)     # Elastic Net regularized regression
library(ranger)     # Random Forest ensemble
library(xgboost)    # Gradient Boosted Trees ensemble

options(survey.lonely.psu = "adjust")
set.seed(2026)


# 2. UNZIP (IF NEEDED) & LOAD DATASET AUTOMATICALLY ----------------------------
zip_files <- list.files(work_dir, pattern = "\\.zip$", full.names = TRUE)
if (length(zip_files) > 0) {
  for (zf in zip_files) {
    unzip(zf, exdir = file.path(work_dir, "extracted_data"), overwrite = FALSE)
  }
}

dta_file <- list.files(work_dir, pattern = "(?i)KEBR8CFL\\.dta$", recursive = TRUE, full.names = TRUE)
sav_file <- list.files(work_dir, pattern = "(?i)KEBR8CFL\\.sav$", recursive = TRUE, full.names = TRUE)
dat_file <- list.files(work_dir, pattern = "(?i)KEBR8CFL\\.dat$", recursive = TRUE, full.names = TRUE)
dct_file <- list.files(work_dir, pattern = "(?i)KEBR8CFL\\.dct$", recursive = TRUE, full.names = TRUE)

if (length(dta_file) > 0) {
  br_raw <- read_dta(dta_file[1])
} else if (length(sav_file) > 0) {
  br_raw <- read_sav(sav_file[1])
} else if (length(dat_file) > 0 && length(dct_file) > 0) {
  dct_lines <- readLines(dct_file[1], warn = FALSE)
  m <- str_match(dct_lines, "_column\\((\\d+)\\)\\s+(\\S+)\\s+(\\S+)\\s+%(\\d+)(?:\\.\\d+)?([sf])")
  m <- m[!is.na(m[, 1]), , drop = FALSE]
  starts <- as.integer(m[, 2]); varnames <- tolower(m[, 4]); widths <- as.integer(m[, 5])
  ends <- starts + widths - 1
  col_types <- paste(ifelse(m[, 6] == "s", "c", "d"), collapse = "")
  br_raw <- read_fwf(dat_file[1], col_positions = fwf_positions(starts, ends, col_names = varnames),
                     col_types = col_types, show_col_types = FALSE)
} else {
  stop("No KEBR8CFL dataset found in: ", work_dir)
}

names(br_raw) <- tolower(names(br_raw))


# 3. FEATURE ENGINEERING FOR SOCIOLOGICAL & CRIMINOLOGICAL ANALYSIS ------------
to_factor_safe <- function(x, fallback_levels = NULL, fallback_labels = NULL) {
  if (is.labelled(x)) as_factor(x)
  else if (!is.null(fallback_levels)) factor(x, levels = fallback_levels, labels = fallback_labels)
  else factor(x)
}

br <- br_raw %>%
  mutate(
    wt         = v005 / 1e6,
    psu        = v001,
    strata     = v022,
    county     = to_factor_safe(v024),
    wealth     = to_factor_safe(v190, 1:5, c("Poorest", "Poorer", "Middle", "Richer", "Richest")),
    education  = to_factor_safe(v106, 0:3, c("No education", "Primary", "Secondary", "Higher")),
    residence  = to_factor_safe(v025, 1:2, c("Urban", "Rural")),
    child_sex  = to_factor_safe(b4, 1:2, c("Male", "Female")),
    
    # Household & Family Structure Indicators
    female_head   = if_else(v151 == 2, 1, 0),
    partner_union = factor(if_else(v501 %in% c(1, 2), "In union", "Single/Widowed/Divorced")),
    birth_order   = bord,
    children_born = v201,
    months_ago    = v008 - b3,
    
    # Life-Course & Spacing Indicators
    mother_age_birth = (b3 - v011) / 12,
    teen_mother      = factor(if_else(mother_age_birth < 20, "Under 20", "20 and older")),
    birth_interval   = if_else(is.na(b11), 48, as.numeric(b11)),
    short_interval   = factor(case_when(
      is.na(b11) ~ "First birth",
      b11 < 24   ~ "Short (<24 mos)",
      TRUE       ~ "Optimal (24+ mos)"
    )),
    
    # Institutional Healthcare Access (recent births module)
    home_birth = case_when(
      m15 %in% c(11, 12, 96) ~ 1,
      m15 >= 20 & m15 < 90   ~ 0,
      TRUE                   ~ NA_real_
    ),
    delivery_cat = factor(case_when(
      m15 >= 20 & m15 < 90   ~ "Health Facility",
      m15 %in% c(11, 12, 96) ~ "Home / Non-Facility",
      TRUE                   ~ "Outside 3-Yr Module"
    ), levels = c("Health Facility", "Home / Non-Facility", "Outside 3-Yr Module")),
    
    inadequate_anc = case_when(
      m14 < 4             ~ 1,
      m14 >= 4 & m14 < 90 ~ 0,
      TRUE                ~ NA_real_
    ),
    anc_cat = factor(case_when(
      m14 >= 4 & m14 < 90 ~ "Adequate (4+ visits)",
      m14 < 4             ~ "Inadequate (<4 visits)",
      TRUE                ~ "Outside 3-Yr Module"
    ), levels = c("Adequate (4+ visits)", "Inadequate (<4 visits)", "Outside 3-Yr Module")),
    
    # Infant Mortality Outcome (<12 months, restricted to births 12-59 mos prior)
    infant_death_5yr = case_when(
      months_ago >= 12 & months_ago < 60 & b5 == 0 & b7 < 12 ~ 1,
      months_ago >= 12 & months_ago < 60                     ~ 0,
      TRUE                                                   ~ NA_real_
    )
  )

kdhs_design <- br %>%
  as_survey_design(ids = psu, strata = strata, weights = wt, nest = TRUE)


# ==============================================================================
# PART A: DESCRIPTIVE, DEMOGRAPHIC & UNSUPERVISED ML TYPOLOGY (TABLES 1-4)
# ==============================================================================

# TABLE 1: Weighted & Unweighted Sample Characteristics
make_desc_rows <- function(var_col, label_prefix) {
  br %>%
    filter(!is.na(.data[[var_col]])) %>%
    group_by(Category = as.character(.data[[var_col]])) %>%
    summarise(
      Variable      = label_prefix,
      Unweighted_N  = n(),
      Weighted_N    = round(sum(wt), 1),
      .groups       = "drop"
    ) %>%
    mutate(Weighted_Pct = round(Weighted_N / sum(Weighted_N) * 100, 2)) %>%
    select(Variable, Category, Unweighted_N, Weighted_N, Weighted_Pct)
}

table01_sample <- bind_rows(
  make_desc_rows("wealth", "Household Wealth Quintile"),
  make_desc_rows("education", "Maternal Education"),
  make_desc_rows("residence", "Type of Residence"),
  make_desc_rows("teen_mother", "Maternal Age at Birth"),
  make_desc_rows("short_interval", "Preceding Birth Interval"),
  make_desc_rows("delivery_cat", "Place of Delivery (Recent Births)"),
  make_desc_rows("anc_cat", "Antenatal Care Utilization")
)
write_csv(table01_sample, file.path(out_dir, "Table01_Sample_Characteristics_Weighted.csv"))

# TABLE 2: Official DHS Synthetic-Cohort Mortality Rates (National, Wealth, Residence)
format_chmort <- function(res, group_name = "Subgroup") {
  if (is.list(res) && !is.data.frame(res)) {
    lbls <- names(res)
    if (is.null(lbls)) lbls <- paste0("Group_", seq_along(res))
    bind_rows(lapply(seq_along(res), function(i) {
      df <- rownames_to_column(as.data.frame(res[[i]]), var = "Indicator")
      df[[group_name]] <- lbls[i]
      relocate(df, all_of(group_name), Indicator)
    }))
  } else {
    rownames_to_column(as.data.frame(res), var = "Indicator") %>%
      mutate(!!group_name := "National") %>%
      relocate(all_of(group_name), Indicator)
  }
}

table02_mortality <- bind_rows(
  format_chmort(chmort(br, Period = 60), "Stratum"),
  format_chmort(chmort(br, Class = "wealth", Period = 60), "Stratum"),
  format_chmort(chmort(br, Class = "residence", Period = 60), "Stratum")
)
write_csv(table02_mortality, file.path(out_dir, "Table02_National_and_Subgroup_Mortality.csv"))

# TABLE 3: 47-County Structural Inequality League Table
county_disparities <- kdhs_design %>%
  group_by(county) %>%
  summarise(
    unweighted_births = n(),
    pct_poorest       = survey_mean(v190 == 1, na.rm = TRUE) * 100,
    pct_no_education  = survey_mean(v106 == 0, na.rm = TRUE) * 100,
    pct_rural         = survey_mean(v025 == 2, na.rm = TRUE) * 100,
    pct_teen_birth    = survey_mean(mother_age_birth < 20, na.rm = TRUE) * 100,
    pct_short_spacing = survey_mean(b11 < 24, na.rm = TRUE) * 100,
    pct_home_birth    = survey_mean(home_birth, na.rm = TRUE) * 100,
    pct_low_anc       = survey_mean(inadequate_anc, na.rm = TRUE) * 100,
    infant_deaths_1k  = survey_mean(infant_death_5yr, na.rm = TRUE) * 1000
  ) %>%
  arrange(desc(pct_home_birth))

write_csv(county_disparities, file.path(out_dir, "Table03_County_Structural_Inequality_47.csv"))

# TABLE 4: Unsupervised Machine Learning — County Structural Vulnerability Regimes (PCA + K-Means)
county_features <- county_disparities %>%
  select(county, pct_poorest, pct_no_education, pct_rural, pct_teen_birth,
         pct_short_spacing, pct_home_birth, pct_low_anc, infant_deaths_1k)

scaled_mat <- scale(county_features %>% select(-county))
pca_res    <- prcomp(scaled_mat, center = FALSE, scale. = FALSE)
km_res     <- kmeans(scaled_mat, centers = 4, nstart = 50)

county_clusters <- county_features %>%
  mutate(
    Cluster_ID = km_res$cluster,
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )

cluster_ranks <- county_clusters %>%
  group_by(Cluster_ID) %>%
  summarise(mean_exclusion = mean(pct_home_birth + pct_poorest)) %>%
  arrange(desc(mean_exclusion)) %>%
  mutate(Structural_Regime = c(
    "Regime 1: Acute Structural Neglect",
    "Regime 2: Rural Institutional Deficit",
    "Regime 3: Transitional Mixed Vulnerability",
    "Regime 4: Institutionalized Advantage"
  ))

table04_clusters <- county_clusters %>%
  left_join(cluster_ranks %>% select(Cluster_ID, Structural_Regime), by = "Cluster_ID") %>%
  arrange(Structural_Regime, desc(pct_home_birth))

write_csv(table04_clusters, file.path(out_dir, "Table04_County_Structural_Clusters_ML.csv"))


# ==============================================================================
# PART B: ECONOMETRIC SURVEY LOGISTIC REGRESSION (TABLES 5-6)
# ==============================================================================

clean_or_terms <- function(term_vec) {
  case_when(
    term_vec == "(Intercept)"                    ~ "(Intercept)",
    term_vec == "wealthPoorer"                   ~ "Wealth: Poorer (vs. Poorest)",
    term_vec == "wealthMiddle"                   ~ "Wealth: Middle (vs. Poorest)",
    term_vec == "wealthRicher"                   ~ "Wealth: Richer (vs. Poorest)",
    term_vec == "wealthRichest"                  ~ "Wealth: Richest (vs. Poorest)",
    term_vec == "educationPrimary"               ~ "Education: Primary (vs. None)",
    term_vec == "educationSecondary"             ~ "Education: Secondary (vs. None)",
    term_vec == "educationHigher"                ~ "Education: Higher (vs. None)",
    term_vec == "residenceRural"                 ~ "Residence: Rural (vs. Urban)",
    term_vec == "home_birth"                     ~ "Delivery: Home / Non-Facility",
    term_vec == "teen_motherUnder 20"            ~ "Maternal Age: Under 20 (vs. 20+)",
    term_vec == "short_intervalOptimal (24+ mos)"~ "Spacing: Optimal 24+ mos (vs. First)",
    term_vec == "short_intervalShort (<24 mos)"  ~ "Spacing: Short <24 mos (vs. First)",
    term_vec == "female_head"                    ~ "Household: Female-Headed",
    term_vec == "child_sexFemale"                ~ "Child Sex: Female (vs. Male)",
    TRUE                                         ~ term_vec
  )
}

model_home_birth <- svyglm(
  home_birth ~ wealth + education + residence + teen_mother + short_interval + female_head,
  design = kdhs_design,
  family = quasibinomial()
)

model_infant_death <- svyglm(
  infant_death_5yr ~ wealth + education + residence + home_birth + teen_mother + short_interval + child_sex,
  design = kdhs_design,
  family = quasibinomial()
)

table05_or_home <- tidy(model_home_birth, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(term = clean_or_terms(term)) %>%
  select(term, aOR = estimate, conf.low, conf.high, p.value)

table06_or_death <- tidy(model_infant_death, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(term = clean_or_terms(term)) %>%
  select(term, aOR = estimate, conf.low, conf.high, p.value)

write_csv(table05_or_home,  file.path(out_dir, "Table05_Survey_Logistic_HomeBirth_aOR.csv"))
write_csv(table06_or_death, file.path(out_dir, "Table06_Survey_Logistic_InfantDeath_aOR.csv"))


# ==============================================================================
# PART C: SUPERVISED ML, DEEP LEARNING & STACKED ENSEMBLE (TABLES 7-10)
# ==============================================================================
ml_df <- br %>%
  filter(!is.na(infant_death_5yr)) %>%
  left_join(table04_clusters %>% select(county, Structural_Regime), by = "county") %>%
  mutate(
    Structural_Regime = factor(Structural_Regime),
    y = as.integer(infant_death_5yr)
  ) %>%
  select(y, wt, psu, wealth, education, residence, Structural_Regime,
         mother_age_birth, birth_interval, birth_order, children_born,
         female_head, partner_union, child_sex, delivery_cat, anc_cat) %>%
  drop_na()

# Cluster-Grouped 75/25 Train-Test Split
unique_psus <- unique(ml_df$psu)
train_psus  <- sample(unique_psus, size = floor(0.75 * length(unique_psus)))

train_df <- ml_df %>% filter(psu %in% train_psus)
test_df  <- ml_df %>% filter(!psu %in% train_psus)

psu_fold_map <- tibble(
  psu  = unique(train_df$psu),
  fold = sample(rep(1:5, length.out = length(unique(train_df$psu))))
)
train_df <- train_df %>% left_join(psu_fold_map, by = "psu")

feature_cols <- c("wealth", "education", "residence", "Structural_Regime",
                  "mother_age_birth", "birth_interval", "birth_order", "children_born",
                  "female_head", "partner_union", "child_sex", "delivery_cat", "anc_cat")

mm_formula <- as.formula(paste("~", paste(feature_cols, collapse = " + "), "- 1"))
X_all      <- model.matrix(mm_formula, data = bind_rows(train_df, test_df))
X_train    <- X_all[1:nrow(train_df), , drop = FALSE]
X_test     <- X_all[(nrow(train_df) + 1):nrow(X_all), , drop = FALSE]

build_aligned_X <- function(df) {
  mat <- model.matrix(mm_formula, data = df)
  aligned <- matrix(0, nrow = nrow(mat), ncol = ncol(X_train))
  colnames(aligned) <- colnames(X_train)
  common_cols <- intersect(colnames(mat), colnames(X_train))
  aligned[, common_cols] <- mat[, common_cols, drop = FALSE]
  aligned
}

col_means <- colMeans(X_train)
col_sds   <- apply(X_train, 2, sd)
col_sds[col_sds == 0] <- 1
X_train_sc <- scale(X_train, center = col_means, scale = col_sds)
X_test_sc  <- scale(X_test,  center = col_means, scale = col_sds)

y_train <- train_df$y
w_train <- train_df$wt
y_test  <- test_df$y
w_test  <- test_df$wt

# --- Custom Vectorized 3-Layer Deep Neural Network (MLP: 64 -> 32 -> 16 -> 1) ---
relu_mat <- function(z) {
  out <- pmax(z, 0)
  dim(out) <- dim(z)
  out
}

sigmoid_vec <- function(z) {
  1 / (1 + exp(-pmin(pmax(z, -20), 20)))
}

train_deep_mlp <- function(X, y, w, epochs = 35, batch_size = 256, lr = 0.003, l2 = 1e-4) {
  X <- as.matrix(X)
  d_in <- ncol(X)
  h1 <- 64; h2 <- 32; h3 <- 16
  
  W1 <- matrix(rnorm(d_in * h1, 0, sqrt(2 / d_in)), d_in, h1); b1 <- rep(0, h1)
  W2 <- matrix(rnorm(h1 * h2, 0, sqrt(2 / h1)), h1, h2);     b2 <- rep(0, h2)
  W3 <- matrix(rnorm(h2 * h3, 0, sqrt(2 / h2)), h2, h3);     b3 <- rep(0, h3)
  W4 <- matrix(rnorm(h3 * 1,  0, sqrt(2 / h3)), h3, 1);      b4 <- 0
  
  w_norm <- w / mean(w)
  n <- nrow(X)
  
  for (ep in seq_len(epochs)) {
    idx <- sample(n)
    for (start_i in seq(1, n, by = batch_size)) {
      b_idx <- idx[start_i:min(start_i + batch_size - 1, n)]
      Xb <- X[b_idx, , drop = FALSE]
      yb <- y[b_idx]
      wb <- w_norm[b_idx]
      m  <- length(b_idx)
      
      # Forward pass preserving 2D matrix dimensions
      Z1 <- sweep(Xb %*% W1, 2, b1, "+"); A1 <- relu_mat(Z1)
      Z2 <- sweep(A1 %*% W2, 2, b2, "+"); A2 <- relu_mat(Z2)
      Z3 <- sweep(A2 %*% W3, 2, b3, "+"); A3 <- relu_mat(Z3)
      Z4 <- sweep(A3 %*% W4, 2, b4, "+"); P  <- as.vector(sigmoid_vec(Z4))
      
      # Backpropagation
      dZ4 <- matrix((P - yb) * wb / m, m, 1)
      dW4 <- t(A3) %*% dZ4 + l2 * W4; db4 <- sum(dZ4)
      
      dA3 <- dZ4 %*% t(W4); dZ3 <- dA3 * (Z3 > 0)
      dW3 <- t(A2) %*% dZ3 + l2 * W3; db3 <- colSums(dZ3)
      
      dA2 <- dZ3 %*% t(W3); dZ2 <- dA2 * (Z2 > 0)
      dW2 <- t(A1) %*% dZ2 + l2 * W2; db2 <- colSums(dZ2)
      
      dA1 <- dZ2 %*% t(W2); dZ1 <- dA1 * (Z1 > 0)
      dW1 <- t(Xb) %*% dZ1 + l2 * W1; db1 <- colSums(dZ1)
      
      # Parameter updates
      W4 <- W4 - lr * dW4; b4 <- b4 - lr * db4
      W3 <- W3 - lr * dW3; b3 <- b3 - lr * db3
      W2 <- W2 - lr * dW2; b2 <- b2 - lr * db2
      W1 <- W1 - lr * dW1; b1 <- b1 - lr * db1
    }
  }
  list(W1 = W1, b1 = b1, W2 = W2, b2 = b2, W3 = W3, b3 = b3, W4 = W4, b4 = b4)
}

predict_deep_mlp <- function(model, X) {
  X <- as.matrix(X)
  A1 <- relu_mat(sweep(X  %*% model$W1, 2, model$b1, "+"))
  A2 <- relu_mat(sweep(A1 %*% model$W2, 2, model$b2, "+"))
  A3 <- relu_mat(sweep(A2 %*% model$W3, 2, model$b3, "+"))
  as.vector(sigmoid_vec(sweep(A3 %*% model$W4, 2, model$b4, "+")))
}

# --- Generate Out-of-Fold Predictions for Stacked Super-Learner Ensemble ---
oof_preds <- matrix(0, nrow = nrow(train_df), ncol = 5)
colnames(oof_preds) <- c("Logistic", "ElasticNet", "RandomForest", "XGBoost", "DeepMLP")

glm_formula <- as.formula(paste("y ~", paste(feature_cols, collapse = " + ")))

for (k in 1:5) {
  tr_idx <- which(train_df$fold != k)
  va_idx <- which(train_df$fold == k)
  
  m_glm <- glm(glm_formula, data = train_df[tr_idx, ], weights = wt, family = quasibinomial())
  oof_preds[va_idx, "Logistic"] <- predict(m_glm, newdata = train_df[va_idx, ], type = "response")
  
  m_enet <- cv.glmnet(X_train_sc[tr_idx, ], y_train[tr_idx], weights = w_train[tr_idx],
                      family = "binomial", alpha = 0.5, nfolds = 5)
  oof_preds[va_idx, "ElasticNet"] <- as.vector(predict(m_enet, X_train_sc[va_idx, ], s = "lambda.min", type = "response"))
  
  rf_tr <- train_df[tr_idx, c("y", feature_cols)] %>% mutate(y = factor(y))
  m_rf  <- ranger(y ~ ., data = rf_tr, case.weights = w_train[tr_idx],
                  probability = TRUE, num.trees = 250, min.node.size = 30, seed = 2026 + k)
  oof_preds[va_idx, "RandomForest"] <- predict(m_rf, data = train_df[va_idx, feature_cols])$predictions[, "1"]
  
  dtr <- xgb.DMatrix(data = X_train[tr_idx, ], label = y_train[tr_idx], weight = w_train[tr_idx])
  m_xgb <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "logloss",
                                   max_depth = 3, eta = 0.04, subsample = 0.8, colsample_bytree = 0.8),
                     data = dtr, nrounds = 140, verbose = 0)
  oof_preds[va_idx, "XGBoost"] <- predict(m_xgb, X_train[va_idx, ])
  
  m_mlp <- train_deep_mlp(X_train_sc[tr_idx, ], y_train[tr_idx], w_train[tr_idx], epochs = 30)
  oof_preds[va_idx, "DeepMLP"] <- predict_deep_mlp(m_mlp, X_train_sc[va_idx, ])
}

clip_p <- function(p) pmin(pmax(p, 1e-5), 1 - 1e-5)
qlogis_mat <- qlogis(clip_p(oof_preds))
meta_fit <- glm(y_train ~ qlogis_mat, weights = w_train, family = quasibinomial())

# --- Fit Full Training Models & Evaluate on Held-Out Cluster Test Set ---
full_glm  <- glm(glm_formula, data = train_df, weights = wt, family = quasibinomial())
full_enet <- cv.glmnet(X_train_sc, y_train, weights = w_train, family = "binomial", alpha = 0.5, foldid = train_df$fold)

rf_full_df <- train_df[, c("y", feature_cols)] %>% mutate(y = factor(y))
full_rf    <- ranger(y ~ ., data = rf_full_df, case.weights = w_train,
                     probability = TRUE, importance = "permutation", num.trees = 400, min.node.size = 30, seed = 2026)

dtrain_full <- xgb.DMatrix(data = X_train, label = y_train, weight = w_train)
full_xgb    <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "logloss",
                                       max_depth = 3, eta = 0.04, subsample = 0.8, colsample_bytree = 0.8),
                         data = dtrain_full, nrounds = 150, verbose = 0)

full_mlp <- train_deep_mlp(X_train_sc, y_train, w_train, epochs = 40)

test_preds <- tibble(
  Logistic     = predict(full_glm, newdata = test_df, type = "response"),
  ElasticNet   = as.vector(predict(full_enet, X_test_sc, s = "lambda.min", type = "response")),
  RandomForest = predict(full_rf, data = test_df[, feature_cols])$predictions[, "1"],
  XGBoost      = predict(full_xgb, X_test),
  DeepMLP      = predict_deep_mlp(full_mlp, X_test_sc)
)

test_qlogis <- qlogis(clip_p(as.matrix(test_preds)))
meta_coefs  <- coef(meta_fit)
meta_coefs[is.na(meta_coefs)] <- 0
test_preds$StackedEnsemble <- plogis(as.vector(cbind(1, test_qlogis) %*% meta_coefs))

calc_weighted_metrics <- function(y_true, p_pred, weights) {
  p_pred <- clip_p(p_pred)
  ord    <- order(p_pred, decreasing = TRUE)
  y_s    <- y_true[ord]; w_s <- weights[ord]
  
  cum_tp <- cumsum(w_s * y_s)
  cum_fp <- cumsum(w_s * (1 - y_s))
  tpr    <- c(0, cum_tp / sum(w_s * y_s))
  fpr    <- c(0, cum_fp / sum(w_s * (1 - y_s)))
  auc    <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  
  brier   <- sum(weights * (y_true - p_pred)^2) / sum(weights)
  logloss <- -sum(weights * (y_true * log(p_pred) + (1 - y_true) * log(1 - p_pred))) / sum(weights)
  
  top10_cut <- quantile(p_pred, 0.90)
  top10_capture <- sum(weights[y_true == 1 & p_pred >= top10_cut]) / sum(weights[y_true == 1]) * 100
  
  list(AUC = round(auc, 4), Brier = round(brier, 5), LogLoss = round(logloss, 4),
       TopDecile_Death_Capture_Pct = round(top10_capture, 2),
       ROC_Curve = tibble(FPR = fpr, TPR = tpr))
}

eval_list <- lapply(names(test_preds), function(m_name) {
  res <- calc_weighted_metrics(y_test, test_preds[[m_name]], w_test)
  tibble(Model = m_name, Weighted_ROC_AUC = res$AUC, Weighted_LogLoss = res$LogLoss,
         Weighted_Brier = res$Brier, Top10_Risk_Death_Capture_Pct = res$TopDecile_Death_Capture_Pct)
})

# TABLE 7: Out-of-Sample ML, Deep Learning & Ensemble Benchmark
table07_benchmark <- bind_rows(eval_list) %>% arrange(desc(Weighted_ROC_AUC))
write_csv(table07_benchmark, file.path(out_dir, "Table07_ML_DL_Ensemble_Model_Comparison.csv"))

# TABLE 8: Cross-Model Explainable AI (XAI) Feature Importance
rf_imp <- tibble(Feature = names(full_rf$variable.importance),
                 RF_Permutation_Imp = as.numeric(full_rf$variable.importance))

base_mlp_loss <- calc_weighted_metrics(y_test, test_preds$DeepMLP, w_test)$LogLoss
mlp_imp_vals <- sapply(feature_cols, function(fc) {
  perm_df <- test_df
  perm_df[[fc]] <- sample(perm_df[[fc]])
  X_perm <- build_aligned_X(perm_df)
  X_perm_sc <- scale(X_perm, center = col_means, scale = col_sds)
  perm_p <- predict_deep_mlp(full_mlp, X_perm_sc)
  max(0, calc_weighted_metrics(y_test, perm_p, w_test)$LogLoss - base_mlp_loss)
})

base_xgb_loss <- calc_weighted_metrics(y_test, test_preds$XGBoost, w_test)$LogLoss
xgb_imp_vals <- sapply(feature_cols, function(fc) {
  perm_df <- test_df
  perm_df[[fc]] <- sample(perm_df[[fc]])
  X_perm <- build_aligned_X(perm_df)
  perm_p <- predict(full_xgb, X_perm)
  max(0, calc_weighted_metrics(y_test, perm_p, w_test)$LogLoss - base_xgb_loss)
})

safe_norm <- function(v) {
  v_pos <- pmax(0, v)
  mx <- max(v_pos, na.rm = TRUE)
  if (is.na(mx) || mx <= 0) rep(0, length(v_pos)) else v_pos / mx * 100
}

clean_feature_names <- function(f_vec) {
  recode(f_vec,
         "wealth"            = "Household Wealth Quintile",
         "education"         = "Maternal Education",
         "residence"         = "Urban / Rural Residence",
         "Structural_Regime" = "County Structural Regime",
         "mother_age_birth"  = "Maternal Age at Birth",
         "birth_interval"    = "Preceding Birth Interval",
         "birth_order"       = "Child Birth Order",
         "children_born"     = "Total Children Ever Born",
         "female_head"       = "Female-Headed Household",
         "partner_union"     = "Marital / Union Status",
         "child_sex"         = "Sex of Child",
         "delivery_cat"      = "Place of Delivery",
         "anc_cat"           = "Antenatal Care Adequacy"
  )
}

table08_xai <- rf_imp %>%
  mutate(
    XGB_Permutation_Imp = xgb_imp_vals[Feature],
    MLP_Permutation_Imp = mlp_imp_vals[Feature],
    RF_Norm  = safe_norm(RF_Permutation_Imp),
    XGB_Norm = safe_norm(XGB_Permutation_Imp),
    MLP_Norm = safe_norm(MLP_Permutation_Imp),
    Ensemble_Consensus_Score = round((RF_Norm + XGB_Norm + MLP_Norm) / 3, 2),
    Feature_Label = clean_feature_names(Feature)
  ) %>%
  arrange(desc(Ensemble_Consensus_Score))

write_csv(table08_xai %>% select(Feature = Feature_Label, RF_Permutation_Imp, XGB_Permutation_Imp,
                                 MLP_Permutation_Imp, RF_Norm, XGB_Norm, MLP_Norm, Ensemble_Consensus_Score),
          file.path(out_dir, "Table08_Ensemble_Feature_Importance_XAI.csv"))

predict_stacked_ensemble <- function(newdata) {
  X_new    <- build_aligned_X(newdata)
  X_new_sc <- scale(X_new, center = col_means, scale = col_sds)
  p_mat <- cbind(
    predict(full_glm, newdata = newdata, type = "response"),
    as.vector(predict(full_enet, X_new_sc, s = "lambda.min", type = "response")),
    predict(full_rf, data = newdata[, feature_cols])$predictions[, "1"],
    predict(full_xgb, X_new),
    predict_deep_mlp(full_mlp, X_new_sc)
  )
  ql <- qlogis(clip_p(p_mat))
  plogis(as.vector(cbind(1, ql) %*% meta_coefs))
}

# TABLE 9: Intersectional Risk Matrix (Wealth x Education x Residence via Stacked Ensemble)
intersectional_grid <- expand_grid(
  wealth    = levels(ml_df$wealth),
  education = levels(ml_df$education),
  residence = levels(ml_df$residence)
) %>%
  mutate(
    wealth            = factor(wealth, levels = levels(ml_df$wealth)),
    education         = factor(education, levels = levels(ml_df$education)),
    residence         = factor(residence, levels = levels(ml_df$residence)),
    Structural_Regime = factor("Regime 2: Rural Institutional Deficit", levels = levels(ml_df$Structural_Regime)),
    mother_age_birth  = median(ml_df$mother_age_birth),
    birth_interval    = median(ml_df$birth_interval),
    birth_order       = median(ml_df$birth_order),
    children_born     = median(ml_df$children_born),
    female_head       = 0,
    partner_union     = factor("In union", levels = levels(ml_df$partner_union)),
    child_sex         = factor("Male", levels = levels(ml_df$child_sex)),
    delivery_cat      = factor("Health Facility", levels = levels(ml_df$delivery_cat)),
    anc_cat           = factor("Adequate (4+ visits)", levels = levels(ml_df$anc_cat))
  )

intersectional_grid$Predicted_IMR_per_1000 <- round(predict_stacked_ensemble(intersectional_grid) * 1000, 2)
table09_intersectional <- intersectional_grid %>%
  select(Residence = residence, Wealth = wealth, Education = education, Predicted_IMR_per_1000) %>%
  arrange(desc(Predicted_IMR_per_1000))

write_csv(table09_intersectional, file.path(out_dir, "Table09_Intersectional_Risk_Matrix.csv"))

# TABLE 10: Criminological / Public Health Risk Concentration by Predicted Decile
test_scored <- test_df %>%
  mutate(
    Ensemble_Risk = test_preds$StackedEnsemble,
    Risk_Decile   = ntile(Ensemble_Risk, 10)
  )

table10_concentration <- test_scored %>%
  group_by(Risk_Decile) %>%
  summarise(
    Unweighted_Births   = n(),
    Weighted_Births     = round(sum(wt), 1),
    Observed_IMR_1k     = round(sum(wt * y) / sum(wt) * 1000, 2),
    Mean_Predicted_1k   = round(sum(wt * Ensemble_Risk) / sum(wt) * 1000, 2),
    Weighted_Deaths     = sum(wt * y)
  ) %>%
  arrange(desc(Risk_Decile)) %>%
  mutate(
    Share_of_All_Deaths_Pct = round(Weighted_Deaths / sum(Weighted_Deaths) * 100, 2),
    Cumulative_Deaths_Pct   = round(cumsum(Share_of_All_Deaths_Pct), 2)
  ) %>%
  select(-Weighted_Deaths)

write_csv(table10_concentration, file.path(out_dir, "Table10_Risk_Decile_Calibration_Concentration.csv"))


# ==============================================================================
# PART D: 10 PUBLICATION-READY 300-DPI FIGURES (NO TITLES, SUBTITLES, OR CAPTIONS)
# ==============================================================================
base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_blank(),
    plot.subtitle    = element_blank(),
    plot.caption     = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# FIGURE 1: Class Gradient in Institutional Exclusion (Bar Chart + Aligned Frequency Polygon)
wealth_grad <- kdhs_design %>%
  group_by(wealth) %>%
  summarise(
    Home_Birth_Pct = survey_mean(home_birth, na.rm = TRUE) * 100,
    IMR_per_1k     = survey_mean(infant_death_5yr, na.rm = TRUE) * 1000
  )

fig01 <- ggplot(wealth_grad, aes(x = wealth, y = Home_Birth_Pct, group = 1)) +
  geom_col(fill = "#2B5C8F", width = 0.58) +
  geom_line(color = "#C0392B", linewidth = 1.2) +
  geom_point(color = "#C0392B", size = 3.5) +
  geom_text(aes(label = sprintf("%.1f%%", Home_Birth_Pct)), vjust = -0.8, fontface = "bold", size = 3.6) +
  scale_y_continuous(limits = c(0, max(wealth_grad$Home_Birth_Pct) * 1.15),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title    = "",
    subtitle = "",
    x        = "Household Wealth Quintile",
    y        = "Non-Facility / Home Births (%)"
  ) + base_theme
ggsave(file.path(out_dir, "Figure01_Class_Gradient_Exclusion_Mortality.png"), fig01, width = 9, height = 6, dpi = 300)

# FIGURE 2: Geography as Destiny (Top 10 vs. Bottom 10 Counties)
county_extremes <- bind_rows(
  county_disparities %>% slice_head(n = 10) %>% mutate(Tier = "Highest Exclusion (Top 10 Counties)"),
  county_disparities %>% slice_tail(n = 10) %>% mutate(Tier = "Highest Access (Bottom 10 Counties)")
)

fig02 <- ggplot(county_extremes, aes(x = pct_home_birth, y = reorder(county, pct_home_birth), color = Tier)) +
  geom_segment(aes(x = 0, xend = pct_home_birth, y = county, yend = county), linewidth = 0.85) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Highest Exclusion (Top 10 Counties)" = "#C0392B",
                                "Highest Access (Bottom 10 Counties)" = "#1F618D")) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "",
    subtitle = "",
    x        = "Non-Facility / Home Births (%)",
    y        = NULL,
    color    = NULL
  ) + base_theme + theme(legend.position = "top")
ggsave(file.path(out_dir, "Figure02_County_Divide_Top_vs_Bottom.png"), fig02, width = 10, height = 7, dpi = 300)

# FIGURE 3: Unsupervised ML — PCA Map of County Structural Regimes
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
fig03 <- ggplot(table04_clusters, aes(x = PC1, y = PC2, color = Structural_Regime)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_text(aes(label = county), vjust = -0.6, size = 2.9, check_overlap = TRUE) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "",
    subtitle = "",
    x        = sprintf("Principal Component 1: Structural Deprivation Axis (%.1f%% variance)", var_exp[1]),
    y        = sprintf("Principal Component 2: Demographic & Spacing Axis (%.1f%% variance)", var_exp[2]),
    color    = NULL
  ) + base_theme + theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "Figure03_County_Structural_Typology_PCA.png"), fig03, width = 10.5, height = 7.5, dpi = 300)

# FIGURE 4: Structural Regime Profiles Heatmap
regime_profile <- table04_clusters %>%
  group_by(Structural_Regime) %>%
  summarise(
    `Poorest Quintile (%)`   = mean(pct_poorest),
    `No Education (%)`       = mean(pct_no_education),
    `Rural Residence (%)`    = mean(pct_rural),
    `Teen Births (%)`        = mean(pct_teen_birth),
    `Short Spacing (%)`      = mean(pct_short_spacing),
    `Home Births (%)`        = mean(pct_home_birth),
    `Inadequate ANC (%)`     = mean(pct_low_anc),
    `Infant Deaths (/1,000)` = mean(infant_deaths_1k)
  ) %>%
  pivot_longer(-Structural_Regime, names_to = "Indicator", values_to = "Mean_Value")

fig04 <- ggplot(regime_profile, aes(x = Indicator, y = fct_rev(Structural_Regime), fill = Mean_Value)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", Mean_Value)), color = "black", fontface = "bold", size = 3.5) +
  scale_fill_gradient(low = "#EBF5FB", high = "#E74C3C") +
  labs(
    title    = "",
    subtitle = "",
    x        = NULL,
    y        = NULL,
    fill     = "Mean Rate"
  ) + base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "Figure04_Cluster_Structural_Profiles.png"), fig04, width = 10.5, height = 5.5, dpi = 300)

# FIGURE 5: Forest Plot of Adjusted Odds Ratios (Survey Logistic Regression)
or_plot_df <- bind_rows(
  table05_or_home %>% mutate(Outcome = "Outcome A: Non-Facility (Home) Birth"),
  table06_or_death %>% mutate(Outcome = "Outcome B: Infant Death (<12 Months)")
) %>%
  filter(term != "(Intercept)")

fig05 <- ggplot(or_plot_df, aes(x = aOR, y = term, xmin = conf.low, xmax = conf.high, color = Outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_errorbarh(height = 0.25, linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_log10() +
  facet_wrap(~Outcome, scales = "free_x") +
  scale_color_manual(values = c("#1F618D", "#C0392B")) +
  labs(
    title    = "",
    subtitle = "",
    x        = "Adjusted Odds Ratio (log scale)",
    y        = NULL
  ) + base_theme + theme(legend.position = "none")
ggsave(file.path(out_dir, "Figure05_Forest_Plot_Adjusted_Odds_Ratios.png"), fig05, width = 11.5, height = 6.5, dpi = 300)

# FIGURE 6: Out-of-Sample ROC Curves (Econometric, ML, Deep Learning & Stacked Ensemble)
roc_df <- bind_rows(lapply(names(test_preds), function(m_name) {
  res <- calc_weighted_metrics(y_test, test_preds[[m_name]], w_test)
  res$ROC_Curve %>% mutate(Model = sprintf("%s (AUC = %.3f)", m_name, res$AUC))
}))

fig06 <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title    = "",
    subtitle = "",
    x        = "False Positive Rate (1 - Specificity)",
    y        = "True Positive Rate (Sensitivity)",
    color    = NULL
  ) + base_theme + theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "Figure06_ML_DL_Ensemble_ROC_Curves.png"), fig06, width = 8.5, height = 6.8, dpi = 300)

# FIGURE 7: Explainable AI (XAI) Consensus Feature Importance
xai_long <- table08_xai %>%
  select(Feature_Label, `Random Forest` = RF_Norm, `XGBoost` = XGB_Norm, `Deep Neural Net` = MLP_Norm) %>%
  pivot_longer(-Feature_Label, names_to = "Algorithm", values_to = "Normalized_Importance")

fig07 <- ggplot(xai_long, aes(x = Normalized_Importance, y = reorder(Feature_Label, Normalized_Importance), fill = Algorithm)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_manual(values = c("Random Forest" = "#2E86C1", "XGBoost" = "#E67E22", "Deep Neural Net" = "#8E44AD")) +
  labs(
    title    = "",
    subtitle = "",
    x        = "Normalized Permutation Importance (0-100 Scale)",
    y        = NULL,
    fill     = NULL
  ) + base_theme + theme(legend.position = "top")
ggsave(file.path(out_dir, "Figure07_XAI_Feature_Importance_Comparison.png"), fig07, width = 9.5, height = 6.5, dpi = 300)

# FIGURE 8: Intersectional Structural Risk Heat-Matrix (Wealth x Education x Residence)
fig08 <- ggplot(table09_intersectional, aes(x = Wealth, y = fct_rev(Education), fill = Predicted_IMR_per_1000)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", Predicted_IMR_per_1000)), fontface = "bold", size = 3.6) +
  facet_wrap(~Residence) +
  scale_fill_gradient(low = "#EAFAF1", high = "#C0392B") +
  labs(
    title    = "",
    subtitle = "",
    x        = "Household Wealth Quintile",
    y        = "Maternal Education",
    fill     = "IMR (/1,000)"
  ) + base_theme + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "Figure08_Intersectional_Heat_Matrix.png"), fig08, width = 10.5, height = 5.5, dpi = 300)

# FIGURE 9: Non-Linear Life-Course Risk Curve (Maternal Age at Birth)
age_grid <- tibble(mother_age_birth = seq(15, 44, by = 1)) %>%
  mutate(
    wealth            = factor("Middle", levels = levels(ml_df$wealth)),
    education         = factor("Primary", levels = levels(ml_df$education)),
    residence         = factor("Rural", levels = levels(ml_df$residence)),
    Structural_Regime = factor("Regime 2: Rural Institutional Deficit", levels = levels(ml_df$Structural_Regime)),
    birth_interval    = 36,
    birth_order       = 2,
    children_born     = 2,
    female_head       = 0,
    partner_union     = factor("In union", levels = levels(ml_df$partner_union)),
    child_sex         = factor("Male", levels = levels(ml_df$child_sex)),
    delivery_cat      = factor("Health Facility", levels = levels(ml_df$delivery_cat)),
    anc_cat           = factor("Adequate (4+ visits)", levels = levels(ml_df$anc_cat))
  )
age_grid$Predicted_IMR_1k <- predict_stacked_ensemble(age_grid) * 1000

y_top <- max(age_grid$Predicted_IMR_1k) * 1.12
fig09 <- ggplot(age_grid, aes(x = mother_age_birth, y = Predicted_IMR_1k)) +
  geom_area(fill = "#2B5C8F", alpha = 0.15) +
  geom_line(color = "#2B5C8F", linewidth = 1.3) +
  geom_vline(xintercept = c(20, 35), linetype = "dashed", color = "#C0392B") +
  annotate("text", x = 17.5, y = y_top * 0.96, label = "Adolescent\nRisk Zone", color = "#C0392B", size = 3.3) +
  annotate("text", x = 39.5, y = y_top * 0.96, label = "Advanced Maternal\nAge Risk Zone", color = "#C0392B", size = 3.3) +
  scale_y_continuous(limits = c(0, y_top)) +
  labs(
    title    = "",
    subtitle = "",
    x        = "Mother's Age at Childbirth (Years)",
    y        = "Predicted Infant Mortality Rate (per 1,000 Live Births)"
  ) + base_theme
ggsave(file.path(out_dir, "Figure09_LifeCourse_NonLinear_Risk_Curve.png"), fig09, width = 9, height = 5.5, dpi = 300)

# FIGURE 10: Criminological / Public Health Risk Concentration Curve (Cumulative Gains)
gains_df <- table10_concentration %>%
  arrange(desc(Risk_Decile)) %>%
  mutate(Pop_Share_Pct = row_number() * 10) %>%
  bind_rows(tibble(Pop_Share_Pct = 0, Cumulative_Deaths_Pct = 0)) %>%
  arrange(Pop_Share_Pct)

fig10 <- ggplot(gains_df, aes(x = Pop_Share_Pct, y = Cumulative_Deaths_Pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray55") +
  geom_line(color = "#8E44AD", linewidth = 1.3) +
  geom_point(color = "#8E44AD", size = 3) +
  geom_text(data = gains_df %>% filter(Pop_Share_Pct %in% c(10, 20, 30, 50)),
            aes(label = sprintf("%.1f%%", Cumulative_Deaths_Pct)), vjust = -0.8, fontface = "bold", size = 3.5) +
  scale_x_continuous(breaks = seq(0, 100, 10), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(breaks = seq(0, 100, 10), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "",
    subtitle = "",
    x        = "Cumulative Share of Births Ranked from Highest to Lowest Structural Risk (%)",
    y        = "Cumulative Share of Preventable Infant Deaths (%)"
  ) + base_theme
ggsave(file.path(out_dir, "Figure10_Criminological_Risk_Concentration_Lift.png"), fig10, width = 8.5, height = 6, dpi = 300)

cat("\nCompleted! All 10 CSV tables and 10 title-free 300-DPI PNG figures saved in:\n", out_dir, "\n")