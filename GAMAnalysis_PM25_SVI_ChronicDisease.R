# =============================================================================
# TWO-MODEL GAM ANALYSIS WITH INTERACTION ANALYSIS
# =============================================================================
#
# Title: Examining the Relationship Between PM2.5 Exposure, Social
#        Vulnerability, and Chronic Disease Prevalence Using Generalized
#        Additive Models (GAMs)
#
# Author: Emily Myers
# Institution: Villanova University
# Advisor: Dr. Peleg Kremer
# Committee Members: Dr. Kabin Shakya, Dr. Yimin Zhang, Dr. Bonnie Henderson
# Last Updated: March 25, 2026
#
# =============================================================================
# SCRIPT OVERVIEW
# =============================================================================
#
# R Version Required: >= 4.0.0
# Primary Package: mgcv (>= 1.8-40)
#
# This script fits two GAMs to answer the research questions:
#
#   MODEL 1 (RQ1): PM2.5 as the OUTCOME -- does SVI predict PM2.5 exposure?
#     1A: pm25 ~ s(Theme1) + s(Theme2) + s(Theme3) + s(Theme4) + controls
#     1B: pm25 ~ s(svi_city_mean) + s(svi_within) + controls (with spatial smooth)
#           Overall SVI decomposed into between-city and within-city components.
#           Model 1C (no-spatial baseline) removed -- attenuation % between
#           models with different covariate structures are misleading (Yimin).
#     Family: inverse.gaussian(link = "1/mu^2")
#
#     Geographic control (Model 1): s(longitude, latitude, k=50)
#       -- k=50 ensures the spatial smooth captures the broad atmospheric 
#          gradients of PM2.5, preventing spatial autocorrelation from
#          inflating the significance of SVI predictors.
#
#   MODEL 2 (RQ2 & RQ3): Disease prevalence as the OUTCOME
#     disease ~ s(PM2.5) + s(SVI) + ti(PM2.5, SVI) + controls
#     Family: quasibinomial(link = "logit")
#     
#     Geographic control (Model 2): s(longitude, latitude, k=50)
#       -- Increased from k=30 to k=50 to better account for localized 
#          clustering in health data and ensure a flat variogram (F > 0.95).
#
#
# Research Questions:
#   RQ1: Do communities with a higher SVI disproportionately experience
#        higher levels of PM2.5 air pollution?
#   RQ2: Do higher levels of PM2.5 increase the prevalence of chronic disease?
#   RQ3: Do communities with a high SVI have increased prevalence of chronic
#        disease when exposed to PM2.5 air pollution?
#
# Interaction Analysis (replaces Bhatia-Lin amplification ratios):
#   - Predicted prevalence tables (3x3 and 2x2 grids)
#   - Difference-in-differences interaction contrasts
#   - Marginal effect curves (PM2.5 effect at different SVI levels)
#   - Interaction surface heatmaps
#
# =============================================================================
# STATISTICAL METHODS
# =============================================================================
#
# MODEL 1 -- Family: inverse.gaussian(link = "1/mu^2")
#   PM2.5 is continuous, positive, right-skewed concentration data.
#   Variance proportional to mu^3, appropriate for concentration measurements.
#
# MODEL 2 -- Family: quasibinomial(link = "logit")
#   Disease prevalence is proportion data bounded [0, 1].
#   Data shows underdispersion (Var << p(1-p)).
#   Quasibinomial estimates dispersion parameter from data.
#   Logit link ensures predictions remain within [0, 1].
#
# Geographic Control Strategy (informed by Section 16 sensitivity analysis):
#   Model 1: s(lon, lat, k=20) -- k=10 exhausted df (statistician feedback);
#             k=20 balances flexibility with risk of signal absorption
#   Model 2: s(lon, lat, k=30) -- needed for disease models (F=0.91-1.04)
#   Climate region fixed effects dropped from both -- adding them on top of
#   the spatial smooth does not improve variogram flatness but may absorb
#   real predictor signal (see sensitivity results).
#
# Smoothing Method: REML (more conservative than ML or GCV)
#
# Predictor Transformations:
#   - PM2.5: log-transformed then standardized
#   - SVI: standardized
#   - Population density: log-transformed then standardized
#   - Time: standardized (centered at 2016)
#
# SVI Decomposition (Model 1):
#   svi_city_mean = city-average SVI (between-city component)
#   svi_within    = tract SVI - city mean (within-city component)
#   Allows testing whether the SVI-PM2.5 relationship operates at
#   the city scale vs. the neighborhood scale.
#
# =============================================================================
# OUTPUT STRUCTURE
# =============================================================================
#
# All outputs saved to: [base_path]/Thesis_Results_TwoModel/
#
#   /Models/      - Fitted model objects (.rds)
#   /Tables/      - Summary statistics, effect sizes, prevalence tables (.csv)
#   /Figures/     - Publication-quality plots (.pdf)
#   /Diagnostics/ - Diagnostic plots (.pdf)
#
# =============================================================================

# =============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# =============================================================================
library(magrittr)
suppressPackageStartupMessages({
  library(mgcv)
  library(tidyverse)
  library(data.table)
  library(viridis)
  library(patchwork)
  library(scales)
  library(gridExtra)
  library(car)
  library(gstat)
  library(corrplot)
  library(lmtest)
  library(ggrepel)   # city labels on bubble plot
  library(stringr)   # string cleaning for city names
  library(vegan)     # PERMANOVA
})
# Helper 1: Spatial Residual Analysis
compute_variogram <- function(model, data, n_sample = 5000) {
  dev_resid <- residuals(model, type = "deviance")
  df <- data.frame(resid = dev_resid, lon = data$longitude, lat = data$latitude)
  if (nrow(df) > n_sample) {
    set.seed(SEED)
    df <- df[sample(nrow(df), n_sample), ]
  }
  gstat::variogram(resid ~ 1, locations = ~lon + lat, data = df)
}

# Helper 2: Statistical Validation Metric
flatness_ratio <- function(vario) {
  short_sv <- mean(vario$gamma[1:2])
  long_sv <- mean(tail(vario$gamma, 3))
  return(short_sv / long_sv)
}

# Helper 3: Prediction Grid Generator (defined in Section 13 after data_model is available)

# ---- PATH CONFIGURATION ----
# NOTE FOR USERS:
#   Save this script in the same folder as your data file
#   (MyersThesis_CensusDataFull.csv). The script will automatically detect
#   its own location and use that as the base path.
#   All output folders (Models, Tables, Figures, Diagnostics) will be created
#   inside a "Thesis_Results_Full_Final" subfolder next to this script.

base_path    <- dirname(rstudioapi::getSourceEditorContext()$path)
results_path <- file.path(base_path, "Thesis_Results_Full_Final")

DISEASES <- c("asthma", "diabetes", "copd", "chd", "cancer")
DISEASE_LABELS <- c("Asthma", "Diabetes", "COPD", "CHD", "Cancer")
names(DISEASE_LABELS) <- DISEASES

# SVI Variants: Overall + 4 CDC Themes
SVI_VARIANTS <- list(
  overall = list(raw = "svi_overall", scaled = "svi_scaled",    label = "Overall SVI"),
  theme1  = list(raw = "rpl_theme1",  scaled = "svi_t1_scaled", label = "Theme 1: Socioeconomic"),
  theme2  = list(raw = "rpl_theme2",  scaled = "svi_t2_scaled", label = "Theme 2: Household Char."),
  theme3  = list(raw = "rpl_theme3",  scaled = "svi_t3_scaled", label = "Theme 3: Minority/Language"),
  theme4  = list(raw = "rpl_theme4",  scaled = "svi_t4_scaled", label = "Theme 4: Housing/Transp.")
)

SEED <- 12345
set.seed(SEED)

print_section <- function(title) cat("\n", rep("=", 80), "\n", title, "\n", rep("=", 80), "\n\n", sep = "")

dirs <- c("Models", "Tables", "Figures", "Diagnostics/Model1", "Diagnostics/Model2", "Sensitivity")

# =============================================================================
# SECTION 2: DATA LOADING AND PREPARATION (FIXED)
# =============================================================================
print_section("SECTION 2: DATA LOADING AND PREPARATION")

data_file <- file.path(base_path, "MyersThesis_CensusDataFull.csv")
data_raw <- fread(data_file)

# 1. FORCE LOWERCASE
names(data_raw) <- tolower(names(data_raw))

# 2. REMOVE DUPLICATE COLUMNS CREATED BY TOLOWER
# This fixes the "placename at locations 3 and 311" error by keeping only the first occurrence
data_raw <- data_raw[, !duplicated(names(data_raw)), with = FALSE]


# -----------------------------------------------------------------------------
# INTEGRATED STEP: AGGREGATE CENSUS TRACT POPULATIONS TO CITY LEVEL
# -----------------------------------------------------------------------------
cat("Aggregating tract populations up to the city level...\n")

city_populations <- data_raw %>%
  # Filter out rows where either the city name or population is missing
  filter(!is.na(placename), !is.na(population)) %>%
  
  # Group by the standardized city column
  group_by(placename) %>%
  
  # Compute total populations and keep track of total census tracts per city
  summarise(
    city_total_population = sum(population, na.rm = TRUE),
    total_census_tracts   = n(),
    .groups = "drop"
  ) %>%
  
  # Arrange alphabetically for presentation consistency
  arrange(placename)

# Save the compiled CSV down to your predefined results folder structure
write.csv(city_populations, 
          file.path(results_path, "Tables", "City_Level_Populations.csv"), 
          row.names = FALSE)

cat("Success! City population benchmarks saved to the Tables folder.\n\n")
# -----------------------------------------------------------------------------


data_raw$tract2010 <- as.character(data_raw$tract2010)

# 3. DEFINE STATIC GEOGRAPHIC COLUMNS
geo_cols <- c("placename", "state", "st_abbr", "county", "climate_region", "longitude", "latitude")

# 4. IDENTIFY YEARLY COLUMNS (Prefix-Only)
all_year_cols <- names(data_raw)[grepl("^x(16|18|20|22)_", names(data_raw))]

# ---- 2A: Deduplicate split tracts ----
cat("Handling duplicate tracts (tracts split across multiple cities)...\n")

data_clean <- data_raw %>%
  group_by(tract2010) %>%
  summarise(
    # Pull in all yearly data (x16_..., x18_..., etc.)
    across(all_of(all_year_cols), first),
    
    # Handle population (Note: ensure lowercase 'population' matches names(data_raw))
    population = sum(population, na.rm = TRUE),
    pop_density = ifelse(sum(population, na.rm = TRUE) > 0,
                         sum(pop_density * population, na.rm = TRUE) / sum(population, na.rm = TRUE),
                         mean(pop_density, na.rm = TRUE)),
    
    # Handle geographic identifiers
    state = names(sort(table(state), decreasing = TRUE))[1],
    st_abbr = names(sort(table(st_abbr), decreasing = TRUE))[1],
    county = names(sort(table(county), decreasing = TRUE))[1],
    climate_region = names(sort(table(climate_region), decreasing = TRUE))[1],
    longitude = mean(longitude, na.rm = TRUE),
    latitude  = mean(latitude, na.rm = TRUE),
    placename = names(sort(table(placename), decreasing = TRUE))[1]
    
  )  

# ---- 2B: Reshape to long format ----
SVI_THEMES_VEC <- c("rpl_theme1", "rpl_theme2", "rpl_theme3", "rpl_theme4")
vars_to_reshape <- c(DISEASES, "mean_pm25", "svi_overall", SVI_THEMES_VEC)

# Identify specifically which x16/x18/x20/x22 columns to pivot
year_cols_subset <- names(data_clean)[grepl(paste0("_(", paste(vars_to_reshape, collapse="|"), ")$"), names(data_clean))]

cat("Reshaping data to long format...\n")
data_long <- data_clean %>%
  dplyr::select(all_of(c("tract2010", "population", "pop_density", geo_cols, year_cols_subset))) %>%
  pivot_longer(cols = all_of(year_cols_subset), names_to = "var_year", values_to = "value") %>%
  mutate(
    year = as.integer(gsub("x(\\d{2})_.*", "20\\1", var_year)),
    variable = gsub("x\\d{2}_", "", var_year)
  ) %>%
  dplyr::select(-var_year) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  rename(GEOID = tract2010, pm25 = mean_pm25) %>%
  mutate(time = year - 2016)

cat(sprintf("  Data Ready: %d rows (Tract-Years)\n", nrow(data_long)))

# =============================================================================
# SECTION 3: DATA EXPLORATION AND CLEANING
# =============================================================================

print_section("SECTION 3: DATA EXPLORATION AND CLEANING")

while (!is.null(dev.list())) dev.off()

# ----------------------------------------------------------------------------
# STATISTICIAN FIX: Do NOT drop disease NAs at this stage.
# Model 1 does not use disease variables at all -- dropping disease NAs here
# would unnecessarily shrink the Model 1 dataset.
# Model 2 filters disease NAs disease-by-disease inside the modeling loop,
# so each model only loses rows where *that* disease is missing.
# ----------------------------------------------------------------------------

# Base complete cases for NON-disease predictors only (used for Model 1)
data_complete_m1 <- data_long %>%
  drop_na(all_of(c("pm25", "svi_overall", SVI_THEMES_VEC, "population",
                   "pop_density", "longitude", "latitude", "climate_region", "time")))

# For visualizations and exploration, also keep a version that requires at least
# one non-zero disease value (used downstream in Sections 5, 13, 14, 15)
# This retains rows with some disease missing (NA) but drops all-zero rows
data_complete <- data_long %>%
  drop_na(all_of(c("pm25", "svi_overall", SVI_THEMES_VEC, "population",
                   "pop_density", "longitude", "latitude", "climate_region", "time"))) %>%
  # Only remove rows where ALL disease values are simultaneously 0
  # (indicates a data entry error, not a legitimate missing value for one disease)
  filter(!(asthma == 0 & diabetes == 0 & copd == 0 & chd == 0 & cancer == 0) |
           is.na(asthma) | is.na(diabetes) | is.na(copd) | is.na(chd) | is.na(cancer))

# Distribution diagnostics
tryCatch({
  pdf(file.path(results_path, "Diagnostics", "01_initial_distributions.pdf"),
      width = 12, height = 10)
  par(mfrow = c(3, 3))
  
  hist(data_complete$pm25, main = "PM2.5 Distribution (Raw)", xlab = "ug/m3", col = "steelblue")
  hist(data_complete$svi_overall, main = "Overall SVI", xlab = "Score", col = "coral")
  hist(data_complete$rpl_theme1, main = "Theme 1 (Socioeconomic)", col = "coral1")
  hist(data_complete$rpl_theme2, main = "Theme 2 (Household)", col = "coral2")
  hist(data_complete$rpl_theme3, main = "Theme 3 (Minority)", col = "coral3")
  hist(data_complete$rpl_theme4, main = "Theme 4 (Housing)", col = "coral4")
  
}, error = function(e) cat("ERROR IN PLOTTING: ", e$message, "\n"),
finally = { if (!is.null(dev.list())) dev.off() })

# =============================================================================
# SECTION 4: DATA FILTERING AND TRANSFORMATION
# =============================================================================

print_section("SECTION 4: DATA FILTERING AND TRANSFORMATION")

cat("Filters applied:\n")
cat("  - Remove PM2.5 <= 0 (physically impossible)\n")
cat("  - Remove SVI < -100 (extreme outliers, likely data errors)\n")
cat("  - Remove population = 0\n")
cat("  - Remove rows where ALL five diseases are simultaneously zero\n")
cat("  NOTE: Rows with a single disease missing/zero are retained.\n")
cat("  Model 1 uses data_model_m1 (no disease filter at all).\n")
cat("  Model 2 filters per-disease inside the modeling loop.\n\n")

# Model 1 base data: no disease filtering whatsoever
data_model_m1 <- data_complete_m1 %>%
  filter(pm25 > 0 & svi_overall >= -100 & population > 0)

# Shared model 2 base: disease columns converted to proportions,
# all-zero rows dropped. Per-disease NAs handled inside the loop.
data_model <- data_complete %>%
  filter(pm25 > 0 & svi_overall >= -100 & population > 0) %>%
  mutate(across(all_of(DISEASES), ~./100))

cat(sprintf("  Model 1 base dataset: %d rows\n", nrow(data_model_m1)))
cat(sprintf("  Model 2 base dataset: %d rows\n\n", nrow(data_model)))

# ---- Variable transformations (applied to both datasets) ----
cat("Transforming and Scaling Variables (including Themes 1-4)...\n")

transform_vars <- function(df) {
  df %>%
    mutate(
      pm25_scaled        = as.numeric(scale(log(pm25))),
      svi_scaled         = as.numeric(scale(svi_overall)),
      svi_t1_scaled      = as.numeric(scale(rpl_theme1)),
      svi_t2_scaled      = as.numeric(scale(rpl_theme2)),
      svi_t3_scaled      = as.numeric(scale(rpl_theme3)),
      svi_t4_scaled      = as.numeric(scale(rpl_theme4)),
      pop_density_scaled = as.numeric(scale(log(pop_density + 1))),
      time_scaled        = as.numeric(scale(time))
    ) %>%
    mutate(climate_region = factor(climate_region))
}

data_model_m1 <- transform_vars(data_model_m1)
data_model     <- transform_vars(data_model)

# Scaling parameters for back-transformation in interaction analysis
# Derived from data_model (Model 2 base) for consistency with Model 2 predictions
pm25_log_mean <- mean(log(data_model$pm25))
pm25_log_sd   <- sd(log(data_model$pm25))
svi_mean      <- mean(data_model$svi_overall)
svi_sd        <- sd(data_model$svi_overall)

pm25_to_scaled <- function(x) (log(x) - pm25_log_mean) / pm25_log_sd
svi_to_scaled  <- function(x) (x - svi_mean) / svi_sd

# Scaling parameters for all SVI themes
svi_params <- list(
  overall = list(m = mean(data_model$svi_overall), s = sd(data_model$svi_overall)),
  t1      = list(m = mean(data_model$rpl_theme1),  s = sd(data_model$rpl_theme1)),
  t2      = list(m = mean(data_model$rpl_theme2),  s = sd(data_model$rpl_theme2)),
  t3      = list(m = mean(data_model$rpl_theme3),  s = sd(data_model$rpl_theme3)),
  t4      = list(m = mean(data_model$rpl_theme4),  s = sd(data_model$rpl_theme4))
)

cat("  Scaling parameters stored for all 4 Themes and Overall SVI.\n")
cat("\n--- DESCRIPTIVE STATISTICS (NATIONAL) ---\n")
data_model %>%
  select(pm25, svi_overall, population, all_of(tolower(DISEASES))) %>%
  summary() %>%
  print()

# =============================================================================
# INTEGRATED STEP: GENERATE APPENDIX TABLE 1 (498 CITIES SUMMARY)
# =============================================================================
cat("Generating Appendix Table 1 (City-level summary and EJ Quadrants)...\n")

# 1. Establish national medians from the long-format dataset as objective benchmarks
pm25_national_median <- median(data_model$pm25, na.rm = TRUE)
svi_national_median  <- median(data_model$svi_overall, na.rm = TRUE)

cat(sprintf("  National Medians for EJ Quadrants -> PM2.5: %.2f, SVI: %.3f\n", 
            pm25_national_median, svi_national_median))

# 2. Extract city-level averages and merge with our corrected city populations
appendix_table_1 <- data_model %>%
  group_by(placename, state) %>%
  summarise(
    # Calculate true temporal/spatial averages across tract-years
    average_pm25 = mean(pm25, na.rm = TRUE),
    average_svi  = mean(svi_overall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Bring in our exact corrected total populations calculated in Section 2
  left_join(
    city_populations %>% select(placename, city_total_population), 
    by = "placename"
  ) %>%
  
  # 3. Apply the EJ Quadrant assignment logic based on national thresholds
  mutate(
    ej_quadrant = case_when(
      average_pm25 >= pm25_national_median & average_svi >= svi_national_median ~ "High Pollution / High Vulnerability",
      average_pm25 >= pm25_national_median & average_svi <  svi_national_median ~ "High Pollution / Low Vulnerability",
      average_pm25 <  pm25_national_median & average_svi >= svi_national_median ~ "Low Pollution / High Vulnerability",
      average_pm25 <  pm25_national_median & average_svi <  svi_national_median ~ "Low Pollution / Low Vulnerability",
      TRUE ~ "Unclassified"
    )
  ) %>%
  
  # 4. Clean column names and structure to match your committee requirements
  select(
    City = placename,
    State = state,
    Population = city_total_population,
    `Average PM2.5` = average_pm25,
    `Average SVI` = average_svi,
    `EJ Quadrant` = ej_quadrant
  ) %>%
  
  # Sort alphabetically by State, then by City
  arrange(State, City)

# 5. Format character strings to Title Case for polished table presentation
# (Fixes any lowercase text styling artifacts)
appendix_table_1 <- appendix_table_1 %>%
  mutate(
    City  = stringr::str_to_title(City),
    State = stringr::str_to_title(State)
  )

# 6. Save the table down as a clean CSV file to your Tables directory
output_appendix_path <- file.path(results_path, "Tables", "Appendix_Table1_City_EJ_Quadrants.csv")
write.csv(appendix_table_1, output_appendix_path, row.names = FALSE)

cat(sprintf("Success! Appendix Table 1 generated with %d unique cities.\n", nrow(appendix_table_1)))
cat("Saved cleanly to:", output_appendix_path, "\n\n")


# =============================================================================
# SECTION 5: DISEASE DISTRIBUTION ANALYSIS
# =============================================================================

print_section("SECTION 5: DISEASE DISTRIBUTION ANALYSIS")

tryCatch({
  pdf(file.path(results_path, "Diagnostics", "02_disease_distributions.pdf"),
      width = 12, height = 10)
  par(mfrow = c(3, 2))
  
  disease_summaries <- list()
  for (disease in DISEASES) {
    disease_data <- data_model[[disease]]
    cat(sprintf("%s: Mean = %.4f (%.2f%%), Var/Mean = %.4f\n",
                toupper(disease), mean(disease_data, na.rm = TRUE),
                mean(disease_data, na.rm = TRUE) * 100,
                var(disease_data, na.rm = TRUE) / mean(disease_data, na.rm = TRUE)))
    
    hist(disease_data, main = paste(toupper(disease), "Distribution"),
         xlab = "Prevalence (proportion)", col = viridis(1), breaks = 30)
    
    disease_summaries[[disease]] <- data.frame(
      Disease = toupper(disease),
      Mean = mean(disease_data, na.rm = TRUE),
      Var = var(disease_data, na.rm = TRUE),
      Var_Mean_Ratio = var(disease_data, na.rm = TRUE) / mean(disease_data, na.rm = TRUE)
    )
  }
  
  disease_summary_table <- do.call(rbind, disease_summaries)
  write.csv(disease_summary_table,
            file.path(results_path, "Tables", "disease_summaries.csv"),
            row.names = FALSE)
}, error = function(e) cat("ERROR: ", e$message, "\n"),
finally = { if (!is.null(dev.list())) dev.off() })

# =============================================================================
# SECTION 6: MULTICOLLINEARITY CHECK
# =============================================================================

print_section("SECTION 6: MULTICOLLINEARITY CHECK")

cor_vars <- c("pm25_scaled", "svi_scaled", "pop_density_scaled",
              "time_scaled", "longitude", "latitude")
cor_matrix <- cor(data_model[, cor_vars, drop = FALSE], use = "complete.obs")

pdf(file.path(results_path, "Diagnostics", "04_correlation_matrix.pdf"),
    width = 10, height = 10)
corrplot::corrplot(cor_matrix, method = "color", type = "upper",
                   addCoef.col = "black", tl.col = "black",
                   title = "Predictor Correlation Matrix", mar = c(0, 0, 2, 0))
dev.off()

high_cor <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("WARNING: High correlations detected (|r| > 0.7):\n")
  for (i in 1:nrow(high_cor)) {
    cat(sprintf("  %s <-> %s: r = %.3f\n",
                rownames(cor_matrix)[high_cor[i, 1]],
                colnames(cor_matrix)[high_cor[i, 2]],
                cor_matrix[high_cor[i, 1], high_cor[i, 2]]))
  }
} else {
  cat("No problematic multicollinearity detected (all |r| < 0.7)\n")
}

# =============================================================================
# SECTION 7: SAVE PREPARED DATA
# =============================================================================

print_section("SECTION 7: SAVE PREPARED DATA")

write.csv(data_model,
          file.path(results_path, "Tables", "data_prepared_for_modeling.csv"),
          row.names = FALSE)
cat(sprintf("  Saved: %d rows x %d columns\n", nrow(data_model), ncol(data_model)))

# =============================================================================
# SECTION 8: MODEL 1 (RQ1) -- SVI THEMES & SCALE DECOMPOSITION
# =============================================================================

print_section("SECTION 8: MODEL 1 (RQ1) -- SVI THEMES & SCALE DECOMPOSITION")

# ---- 8A: THEMATIC ANALYSIS (MODEL 1A) ----
# Identifies which specific social factors (Themes 1-4) drive exposure.
cat("--- Part A: Fitting Model 1A (Thematic Analysis) ---\n")

m1a_formula <- as.formula(paste0(
  "pm25 ~ s(rpl_theme1, k=10) + s(rpl_theme2, k=10) + ",
  "s(rpl_theme3, k=10) + s(rpl_theme4, k=10) + ",
  "s(pop_density_scaled, k=15) + ", 
  "s(longitude, latitude, k=50) + ", 
  "s(time, k=4)"
))

# Filter for complete cases
m1a_data <- data_model_m1 %>%
  filter(complete.cases(pm25, rpl_theme1, rpl_theme2, rpl_theme3, rpl_theme4, 
                        pop_density_scaled, longitude, latitude, time))

m1a_fit <- gam(m1a_formula, data = m1a_data, 
               method = "REML", family = inverse.gaussian(link = "1/mu^2"))

saveRDS(m1a_fit, file.path(results_path, "Models", "model1a_thematic.rds"))


# ---- 8B: GEOGRAPHIC SCALE ANALYSIS (MODEL 1B) ----
# Determines if disparities operate at the City scale or Neighborhood scale.
cat("\n--- Part B: Fitting Model 1B (Geographic Scale Analysis) ---\n")

city_stats <- data_model_m1 %>%
  group_by(placename) %>%
  summarise(svi_city_mean = mean(svi_overall, na.rm = TRUE), .groups = "drop")

data_model_m1 <- data_model_m1 %>%
  left_join(city_stats, by = "placename") %>%
  mutate(svi_within = svi_overall - svi_city_mean)

m1b_fit <- gam(pm25 ~ s(svi_city_mean, k=10) + s(svi_within, k=10) +
                 s(pop_density_scaled, k=8) + s(longitude, latitude, k=50) + s(time, k=4),
               data = data_model_m1, 
               method = "REML", 
               family = inverse.gaussian(link = "1/mu^2"))

saveRDS(m1b_fit, file.path(results_path, "Models", "model1b_scale_decomposition.rds"))


# ---- 8C: EFFECT EXTRACTION (For Table 4 & Forest Plot) ----

# Helper to calculate 10th-90th percentile predicted change
calc_diff <- function(model, data, var) {
  p10 <- quantile(data[[var]], 0.1, na.rm=TRUE)
  p90 <- quantile(data[[var]], 0.9, na.rm=TRUE)
  nd1 <- data[1,]; nd2 <- data[1,]
  # Set all vars to median/reference
  for(v in names(data)) if(is.numeric(data[[v]])) { nd1[[v]] <- median(data[[v]], na.rm=T); nd2[[v]] <- nd1[[v]] }
  nd1[[var]] <- p10; nd2[[var]] <- p90
  diff <- predict(model, newdata=nd2, type="response") - predict(model, newdata=nd1, type="response")
  return(as.numeric(diff))
}

# Extract Theme Effects (from 1A)
theme_effects <- data.frame(
  Theme = c("Socioeconomic", "Household/Disability", "Minority/Language", "Housing/Transp"),
  Effect_Size = c(calc_diff(m1a_fit, m1a_data, "rpl_theme1"),
                  calc_diff(m1a_fit, m1a_data, "rpl_theme2"),
                  calc_diff(m1a_fit, m1a_data, "rpl_theme3"),
                  calc_diff(m1a_fit, m1a_data, "rpl_theme4"))
)

# Extract Scale Effects (from 1B) for Table 4
table4_results <- data.frame(
  Scale = c("Between-City (Scale)", "Within-City (Scale)"),
  Predicted_Change = c(calc_diff(m1b_fit, data_model_m1, "svi_city_mean"),
                       calc_diff(m1b_fit, data_model_m1, "svi_within"))
)
write.csv(table4_results, file.path(results_path, "Tables", "Table4_Scale_Effects.csv"), row.names=FALSE)


# ---- 8D: VISUALIZATION (Figure 4 & Forest Plot) ----

# 1. Forest Plot for Themes
p_forest <- ggplot(theme_effects, aes(x = reorder(Theme, Effect_Size), y = Effect_Size)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 4, color = "#e67e22") +
  coord_flip() +
  labs(title = "Influence of SVI Themes on PM2.5", y = "Predicted Increase in PM2.5 (ug/m3)", x = "") +
  theme_minimal()

ggsave(file.path(results_path, "Figures", "Model1_Theme_ForestPlot.tiff"), p_forest, width=8, height=5)

# 2. Figure 4: Scale Decomposition
tiff(file.path(results_path, "Figures", "Figure4_Scale_Decomposition.tiff"), width=12, height=6, units="in", res=300)
par(mfrow = c(1, 2))
plot(m1b_fit, select = 1, main = "Between-City SVI Effect", xlab = "City Mean SVI", 
     shade = TRUE, shade.col="lightblue", shift = coef(m1b_fit)[1], seWithMean = TRUE)
plot(m1b_fit, select = 2, main = "Within-City SVI Effect", xlab = "Neighborhood SVI Deviation", 
     shade = TRUE, shade.col="lightblue", shift = coef(m1b_fit)[1], seWithMean = TRUE)
dev.off()
# =============================================================================
# SECTION 9: MODEL 1 DIAGNOSTICS
# =============================================================================
print_section("SECTION 9: MODEL 1 DIAGNOSTICS")

# Start PDF for all diagnostic plots
# This creates the file you will show your committee if they ask about model fit.
pdf(file.path(results_path, "Diagnostics/Model1", "model1_diagnostics_report.pdf"), 
    width = 12, height = 10)

# ---- 9A: STANDARD GAM CHECK (REFINED SCALE MODEL 1B) ----
# This proves that our k-values (basis dimensions) are high enough 
# for the new Overall SVI decomposition.
cat("Running standard GAM checks for the new Model 1B...\n")
par(mfrow = c(2, 2))
gam.check(m1b_fit) 


# ---- 9B: CONCURVITY CHECK (THEMATIC MODEL 1A) ----
# Proves the 4 SVI themes are independent enough to trust their individual effects.
cat("Checking Concurvity for Model 1A (Themes)...\n")
mc <- concurvity(m1a_fit, full = TRUE)

plot.new()
text(0.5, 0.9, "Concurvity Check (Model 1A: Themes)", cex=1.5, font=2)
# We take the "worst" case estimate (the maximum) to be conservative for the defense.
text(0.5, 0.7, paste("Worst Concurvity (Estimate):", round(max(mc[3,]), 3)), cex=1.2)
text(0.5, 0.5, "Note: Values < 0.8 indicate themes are sufficiently independent.", cex=1.1)


# ---- 9C: SPATIAL AUTOCORRELATION (VARIOGRAM) ----
# CRITICAL DEFENSE PROOF: This shows that the s(lon, lat) smooth 
# successfully removed geographic bias from the scale model.
cat("Generating residual variogram for the new Model 1B...\n")
resids_1b <- residuals(m1b_fit, type = "deviance")

set.seed(SEED)
# Use data_model_m1 -- the dataset m1b_fit was fitted on.
# (data_model is the Model 2 dataset and has a different row count.)
samp_idx <- sample(1:nrow(data_model_m1), min(5000, nrow(data_model_m1)))
vario_data <- data.frame(
  res = resids_1b[samp_idx],
  lon = data_model_m1$longitude[samp_idx],
  lat = data_model_m1$latitude[samp_idx]
)

v_model <- gstat::variogram(res ~ 1, locations = ~lon + lat, data = vario_data)
par(mfrow = c(1, 1))
plot(v_model$dist, v_model$gamma, type="b", pch=19, col="#2c3e50",
     main="Residual Variogram (Scale Model 1B)", 
     xlab="Distance (Geographic)", ylab="Semivariance")
# The red dashed line represents "white noise" (no spatial bias)
abline(h = var(vario_data$res), col="#e74c3c", lty=2, lwd=2)


# ---- 9D: OBSERVED VS. PREDICTED COMPARISON ----
# Visual proof that the models are capturing the variance in PM2.5 correctly.
# Both models fitted on data_model_m1; use it here to align row counts correctly.
cat("Generating Observed vs. Predicted plots...\n")
par(mfrow = c(1, 2))

# Model 1A (Thematic) -- fitted on m1a_data (complete-cases subset of data_model_m1)
plot(m1a_data$pm25, fitted(m1a_fit), 
     main="1A: Thematic Model Fit", xlab="Observed PM2.5", ylab="Predicted PM2.5",
     pch=16, cex=0.4, col=rgb(0.1, 0.1, 0.8, 0.2))
abline(0, 1, col="red", lwd=2)

# Model 1B (Scale - Overall SVI decomposition) -- fitted on data_model_m1
plot(data_model_m1$pm25, fitted(m1b_fit), 
     main="1B: Scale Model Fit", xlab="Observed PM2.5", ylab="Predicted PM2.5",
     pch=16, cex=0.4, col=rgb(0.1, 0.8, 0.1, 0.2))
abline(0, 1, col="red", lwd=2)


# ---- 9E: RESIDUALS VS PREDICTORS ----
# Ensures there is no "shape" left in the errors relative to SVI.
# Use data_model_m1 to match the rows of m1b_fit residuals.
par(mfrow = c(1, 2))
plot(data_model_m1$svi_overall, resids_1b, 
     xlab="Overall SVI", ylab="Deviance Residuals", main="Residuals vs SVI (1B)")
abline(h=0, col="red", lty=2)

plot(data_model_m1$pop_density_scaled, resids_1b, 
     xlab="Pop Density (Scaled)", ylab="Deviance Residuals", main="Residuals vs Density")
abline(h=0, col="red", lty=2)

dev.off()
cat("Section 9 Complete. Updated Diagnostics saved to PDF.\n")

print_section("SECTION 10: EXPORT FINAL MODEL STATISTICS")
# =============================================================================
# SECTION 10: MASTER EXPORT FOR THESIS TABLES
# =============================================================================
library(broom)

# --- 10A: MODEL 1A - THEMATIC DRIVERS (For Table 3) ---
if(exists("m1a_fit")) {
  # 1. Detailed Statistics (F-stats, EDF, P-values)
  m1a_tidy <- broom::tidy(m1a_fit)
  
  # 2. Add the "Effect Size" (Predicted Change) we calculated in Section 8C
  # We match them by theme name to keep everything in one spreadsheet
  m1a_final_table <- m1a_tidy %>%
    mutate(Effect_Label = case_when(
      term == "s(rpl_theme1)" ~ "Socioeconomic",
      term == "s(rpl_theme2)" ~ "Household/Disability",
      term == "s(rpl_theme3)" ~ "Minority/Language",
      term == "s(rpl_theme4)" ~ "Housing/Transp",
      TRUE ~ term
    )) %>%
    left_join(theme_effects, by = c("Effect_Label" = "Theme"))
  
  write.csv(m1a_final_table, 
            file.path(results_path, "Tables", "Model1A_Final_Thematic_Results.csv"), 
            row.names = FALSE)
  
  # 3. Model Fit (R-squared / AIC)
  write.csv(broom::glance(m1a_fit), 
            file.path(results_path, "Tables", "Model1A_Fit_Metrics.csv"), 
            row.names = FALSE)
}

# --- 10B: MODEL 1B - SCALE DECOMPOSITION (For Table 4) ---
if(exists("m1b_fit")) {
  # 1. Detailed Statistics (F-stats, EDF, P-values)
  m1b_tidy <- broom::tidy(m1b_fit)
  
  # 2. Add the Between/Within predicted changes
  m1b_final_table <- m1b_tidy %>%
    mutate(Scale_Label = case_when(
      term == "s(svi_city_mean)" ~ "Between-City (Scale)",
      term == "s(svi_within)"    ~ "Within-City (Scale)",
      TRUE ~ term
    )) %>%
    left_join(table4_results, by = c("Scale_Label" = "Scale"))
  
  write.csv(m1b_final_table, 
            file.path(results_path, "Tables", "Model1B_Final_Scale_Results.csv"), 
            row.names = FALSE)
  
  # 3. Model Fit (R-squared / AIC)
  write.csv(broom::glance(m1b_fit), 
            file.path(results_path, "Tables", "Model1B_Fit_Metrics.csv"), 
            row.names = FALSE)
}

cat("\nDone! Check your 'Tables' folder for the Final_Results CSVs.\n")

# This joins your calculated effect sizes with the actual model significance
m1a_final_clean <- tidy(m1a_fit) %>%
  mutate(Theme = case_when(
    term == "s(rpl_theme1)" ~ "Socioeconomic",
    term == "s(rpl_theme2)" ~ "Household/Disability",
    term == "s(rpl_theme3)" ~ "Minority/Language",
    term == "s(rpl_theme4)" ~ "Housing/Transp",
    TRUE ~ term
  )) %>%
  left_join(theme_effects, by = "Theme")

# Now view the version that should have all the numbers
print(m1a_final_clean)

# =============================================================================
# SECTION 10: PERMANOVA -- REGIONAL DIFFERENCES
# =============================================================================
print_section("SECTION 10: PERMANOVA -- REGIONAL DIFFERENCES")

tryCatch({
  perm_vars <- c(DISEASES, "svi_overall", "pm25", "climate_region")
  df_perm <- na.omit(data_model[, perm_vars])
  
  set.seed(SEED)
  df_sample <- df_perm[sample(nrow(df_perm), min(nrow(df_perm), 2000)), ]
  
  dist_mat <- dist(scale(as.matrix(df_sample[, c(DISEASES, "svi_overall", "pm25")])),
                   method = "euclidean")
  perm_res <- vegan::adonis2(dist_mat ~ climate_region, data = df_sample, permutations = 999)
  
  cat("PERMANOVA (n=2000 sample, 999 permutations):\n")
  print(perm_res)
  write.csv(as.data.frame(perm_res), 
            file.path(results_path, "Tables", "permanova_results.csv"))
}, error = function(e) cat("PERMANOVA skipped: ", e$message, "\n"))

# =============================================================================
# SECTION 11 & 12: MODEL 2 -- DISEASE PREVALENCE & INTEGRATED DIAGNOSTICS
# =============================================================================
print_section("SECTION 11 & 12: MODEL 2 & DIAGNOSTICS")

SVI_MODEL2_VARIANTS <- list(
  overall = list(col = "svi_scaled",    label = "Overall SVI"),
  theme1  = list(col = "svi_t1_scaled", label = "Theme 1: Socioeconomic"),
  theme2  = list(col = "svi_t2_scaled", label = "Theme 2: Household Char."),
  theme3  = list(col = "svi_t3_scaled", label = "Theme 3: Minority/Language"),
  theme4  = list(col = "svi_t4_scaled", label = "Theme 4: Housing/Transp.")
)
for (v_name in names(SVI_MODEL2_VARIANTS)) {
  v_info <- SVI_MODEL2_VARIANTS[[v_name]]
  s_col <- v_info$col
  
  cat(sprintf("\n--- STARTING ANALYSIS FOR: %s ---\n", v_info$label))
  
  for (disease in DISEASES) {
    cat(sprintf("Fitting %s x %s...", toupper(disease), v_info$label))
    
    # 1. Define Formula
    m2_formula <- as.formula(sprintf(
      "%s ~ s(pm25_scaled, k=10) + s(%s, k=10) + 
           ti(pm25_scaled, %s, k=c(5,5)) + 
           s(pop_density_scaled, k=15) + 
           s(longitude, latitude, k=50) + 
           s(time_scaled, k=4)", 
      disease, s_col, s_col
    ))
    
    # 2. Prep Data
    m2_data <- data_model %>%
      filter(complete.cases(.data[[disease]], pm25_scaled, .data[[s_col]], 
                            pop_density_scaled, longitude, latitude, time_scaled))
    
    # 3. Fit Model
    m2_fit <- tryCatch({
      gam(m2_formula, data = m2_data, method = "REML", family = quasibinomial())
    }, error = function(e) {
      cat(" Error: ", e$message, "\n")
      return(NULL)
    })
    
    # 4. Save and Plot
    if (!is.null(m2_fit)) {
      saveRDS(m2_fit, file.path(results_path, "Models", 
                                sprintf("model2_%s_%s.rds", disease, v_name)))
      
      diag_dir <- file.path(results_path, "Diagnostics/Model2", v_name)
      if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
      
      # --- FIX: USE A SAFE STRING FOR COLOR ---
      pdf(file.path(diag_dir, sprintf("surface_%s_%s.pdf", disease, v_name)))
      
      # Using "topo" (Topographic) colors which are robust and clear
      vis.gam(m2_fit, view = c("pm25_scaled", s_col), 
              theta = 35, phi = 35, 
              color = "topo", # This is the specific fix
              type = "response",
              main = paste("Risk Surface:", toupper(disease)))
      
      dev.off()
    }
  }
}
# =============================================================================
# SECTION 13: INTERACTION ANALYSIS & VISUALIZATION (DiD & MARGINAL EFFECTS)
# =============================================================================
print_section("SECTION 13: INTERACTION ANALYSIS")

THEMES <- c("overall", "theme1", "theme2", "theme3", "theme4")
all_theme_results <- list()

for (t in THEMES) {
  cat(sprintf("\n--- Processing SVI %s ---\n", toupper(t)))
  
  target_col <- if(t == "overall") "svi_overall" else paste0("rpl_", t)
  if(!target_col %in% names(data_model)) next
  
  # Calculate scaling
  svi_mean <- mean(data_model[[target_col]], na.rm = TRUE)
  svi_sd   <- sd(data_model[[target_col]], na.rm = TRUE)
  svi_to_sc <- function(x) (x - svi_mean) / svi_sd
  svi_sc_2x2 <- svi_to_sc(quantile(data_model[[target_col]], c(0.10, 0.90), na.rm = TRUE))
  
  theme_contrast_list <- list()
  
  for (disease in DISEASES) {
    model_name <- sprintf("model2_%s_%s.rds", tolower(disease), t)
    model_path <- file.path(results_path, "Models", model_name)
    if (!file.exists(model_path)) next
    model <- readRDS(model_path)
    
    # Grid for DiD Plot
    model_svi_var <- all.vars(formula(model))[grepl("svi", all.vars(formula(model))) & grepl("scaled", all.vars(formula(model)))]
    # Build 2x2 prediction grid inline (make_newdata_fixed was undefined)
    pm_sc_2x2 <- quantile(data_model$pm25_scaled, c(0.10, 0.90), na.rm = TRUE)
    nd2 <- expand.grid(pm25_scaled = pm_sc_2x2, svi_val = svi_sc_2x2,
                       pop_density_scaled = 0,
                       longitude = median(data_model$longitude, na.rm = TRUE),
                       latitude  = median(data_model$latitude,  na.rm = TRUE),
                       time_scaled = 0)
    colnames(nd2)[colnames(nd2) == "svi_val"] <- model_svi_var
    
    p2  <- predict(model, newdata = nd2, type = "response", se.fit = TRUE)
    nd2$fit <- p2$fit * 100
    nd2$se  <- p2$se.fit * 100
    nd2$svi_label <- ifelse(nd2[[model_svi_var]] == min(nd2[[model_svi_var]]), "Low SVI (p10)", "High SVI (p90)")
    nd2$pm25_label <- ifelse(nd2$pm25_scaled == min(nd2$pm25_scaled), "Low PM2.5", "High PM2.5")
    
    # --- NEW: DiD INTERACTION PLOT ---
    p_did <- ggplot(nd2, aes(x = pm25_label, y = fit, group = svi_label, color = svi_label)) +
      geom_line(size = 1.2) + geom_point(size = 3) +
      geom_errorbar(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se), width = 0.1) +
      scale_color_manual(values = c("Low SVI (p10)" = "#3182bd", "High SVI (p90)" = "#d73027")) +
      labs(title = paste("Interaction Penalty:", DISEASE_LABELS[disease]),
           subtitle = paste("Theme:", t), y = "Predicted Prevalence (%)", x = "") +
      theme_minimal()
    
    ggsave(file.path(results_path, "Figures", sprintf("DiD_%s_%s.pdf", disease, t)), p_did, width = 6, height = 5)
    
    # expand.grid row order (pm varies fastest):
    #   [1] pm10+svi10  [2] pm90+svi10  [3] pm10+svi90  [4] pm90+svi90
    # DiD = (high_pm|high_svi - low_pm|high_svi) - (high_pm|low_svi - low_pm|low_svi)
    #     = (row4 - row3) - (row2 - row1)  [equivalent to code below]
    theme_contrast_list[[disease]] <- data.frame(
      Disease = toupper(disease), Theme = t,
      Interaction_DiD = (nd2$fit[4] - nd2$fit[2]) - (nd2$fit[3] - nd2$fit[1])
    )
    rm(model); gc()
  }
  all_theme_results[[t]] <- bind_rows(theme_contrast_list)
}
# Bind all theme results into a single summary (full_contrast_summary was undefined)
full_contrast_summary <- bind_rows(all_theme_results)
write.csv(full_contrast_summary, file.path(results_path, "Tables", "Interaction_DiD_Summary.csv"), row.names = FALSE)
# =============================================================================
# SECTION 13.5: MARGINAL RESPONSE CURVES (SVI effect at PM2.5 levels)
# =============================================================================
print_section("SECTION 13.5: MARGINAL SVI CURVES")

# 1. Define PM2.5 thresholds for comparison (10th, 50th, and 90th percentiles)
pm_levels <- quantile(data_model$pm25_scaled, c(0.1, 0.5, 0.9), na.rm = TRUE)
pm_labels <- c("Low PM2.5 (p10)", "Median PM2.5 (p50)", "High PM2.5 (p90)")

# 2. Setup SVI sequence for the X-axis (0 to 1)
# Recompute svi_mean/svi_sd explicitly from svi_overall here.
# The loop in Section 13 overwrites these with theme-specific values;
# the last iteration leaves them set to Theme 4, not overall SVI.
svi_mean_overall <- mean(data_model$svi_overall, na.rm = TRUE)
svi_sd_overall   <- sd(data_model$svi_overall,   na.rm = TRUE)
svi_seq_raw    <- seq(0, 1, length.out = 100)
svi_seq_scaled <- (svi_seq_raw - svi_mean_overall) / svi_sd_overall

for (disease in DISEASES) {
  # Load the Overall SVI model for this disease
  model_path <- file.path(results_path, "Models", sprintf("model2_%s_overall.rds", disease))
  if (!file.exists(model_path)) next
  model <- readRDS(model_path)
  
  # 3. Build a prediction grid
  # We vary SVI and PM2.5 while holding other factors at their median (0 for scaled vars)
  pred_grid_marginal <- expand.grid(
    svi_scaled = svi_seq_scaled,
    pm25_scaled = pm_levels,
    pop_density_scaled = 0,
    longitude = median(data_model$longitude, na.rm = TRUE),
    latitude = median(data_model$latitude, na.rm = TRUE),
    time_scaled = 0
  )
  
  # Map back raw SVI for plotting and add labels
  pred_grid_marginal$svi_raw <- pred_grid_marginal$svi_scaled * svi_sd_overall + svi_mean_overall
  pred_grid_marginal$pm_level <- factor(pred_grid_marginal$pm25_scaled, 
                                        levels = pm_levels, labels = pm_labels)
  
  # 4. Generate Predictions (Response scale: 0-1)
  preds <- predict(model, newdata = pred_grid_marginal, type = "response", se.fit = TRUE)
  pred_grid_marginal$fit <- preds$fit * 100
  pred_grid_marginal$se  <- preds$se.fit * 100
  
  # 5. Create the Visualization
  p_marginal <- ggplot(pred_grid_marginal, aes(x = svi_raw, y = fit, color = pm_level, fill = pm_level)) +
    # 95% Confidence Ribbons
    geom_ribbon(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se), alpha = 0.1, color = NA) +
    # Marginal Mean Lines
    geom_line(size = 1.2) +
    # Styling
    scale_color_viridis_d(option = "magma", end = 0.8) +
    scale_fill_viridis_d(option = "magma", end = 0.8) +
    labs(title = paste("Vulnerability Amplification:", DISEASE_LABELS[disease]),
         subtitle = "Relationship between SVI and Prevalence at Different Pollution Levels",
         x = "Social Vulnerability Index (SVI)", 
         y = "Predicted Prevalence (%)",
         color = "PM2.5 Level", fill = "PM2.5 Level") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  
  # Save to Figures folder
  ggsave(file.path(results_path, "Figures", sprintf("Marginal_SVI_Curves_%s.pdf", disease)), 
         p_marginal, width = 8, height = 6)
  
  rm(model); gc()
}

# =============================================================================
# SECTION 14: RESULTS VISUALIZATION (FINAL FIXED VERSION)
# =============================================================================
print_section("SECTION 14: RESULTS VISUALIZATION")

# ---- 14.0: DATA SCHEMA ALIGNMENT ----

# Using your specific names: 'placename' and 'svi_overall'
actual_city_col <- "placename"
actual_svi_col  <- "svi_overall"

cat(sprintf("Using '%s' for City and '%s' for SVI grouping...\n", actual_city_col, actual_svi_col))

data_model <- data_model %>%
  group_by(!!sym(actual_city_col)) %>%
  mutate(
    svi_city_mean = mean(!!sym(actual_svi_col), na.rm = TRUE),
    svi_within = !!sym(actual_svi_col) - svi_city_mean
  ) %>%
  ungroup()

# Load Model 1B
m1b_path <- file.path(results_path, "Models", "model1b_scale_decomposition.rds")
if(file.exists(m1b_path)){
  m1b_fit <- readRDS(m1b_path)
} else {
  stop("Model file 'model1b_scale_decomposition.rds' not found in the Models folder.")
}

# Ensure Section 13 data is available
if (exists("full_contrast_summary")) { contrast_summary <- full_contrast_summary }


# ---- 14.1: RQ1 NEIGHBORHOOD EXPOSURE CURVE (FIXED) ----

# 1. Create Prediction Grid
mean_city_svi  <- median(data_model$svi_city_mean, na.rm = TRUE)
svi_within_seq <- seq(min(data_model$svi_within, na.rm = TRUE), 
                      max(data_model$svi_within, na.rm = TRUE), 
                      length.out = 200)

pred_grid_svi <- data.frame(
  svi_within         = svi_within_seq,
  svi_city_mean      = mean_city_svi,
  pop_density_scaled = 0,
  longitude          = median(data_model$longitude, na.rm = TRUE),
  latitude           = median(data_model$latitude, na.rm = TRUE),
  time               = median(data_model$time, na.rm = TRUE)
  # Note: m1b_fit uses raw 'time', not 'time_scaled'; time_scaled removed.
)

# 2. Predict (Family: inverse.gaussian | Link: 1/mu^2)
# The 'm1b_fit' model should now find the 'time' column and run correctly
preds <- predict(m1b_fit, newdata = pred_grid_svi, type = "link", se.fit = TRUE)

# 3. Back-Transform
safe_fit <- pmax(preds$fit, 0.0001)
pred_grid_svi$fit <- 1 / sqrt(safe_fit)
pred_grid_svi$upr <- 1 / sqrt(pmax(preds$fit - (1.96 * preds$se.fit), 0.0001))
pred_grid_svi$lwr <- 1 / sqrt(pmax(preds$fit + (1.96 * preds$se.fit), 0.0001))
pred_grid_svi$absolute_svi <- pred_grid_svi$svi_within + mean_city_svi

# 4. Plot RQ1
p_rq1 <- ggplot(pred_grid_svi, aes(x = absolute_svi, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "#3182bd", alpha = 0.2) +
  geom_line(color = "#3182bd", size = 1.2) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "RQ1: Neighborhood-Scale Exposure Curve",
       subtitle = "Prediction based on local SVI deviations (Within-City Effect)",
       x = "Social Vulnerability Index (SVI)", 
       y = expression(Predicted ~ PM[2.5] ~ (mu*g/m^3))) +
  theme_minimal()

print(p_rq1)
# =============================================================================
# SECTION 14.3: SINGLE OVERALL BUBBLE PLOT (NON-FACETED)
# =============================================================================

print_section("SECTION 14.3: NATIONAL EJ QUADRANT")

library(patchwork)
library(cowplot)
library(ggrepel)
library(scales)

# 1. Prepare City Summary Data
city_summary <- data_model %>%
  mutate(placename = str_to_title(str_trim(placename)),
         state = str_to_title(str_trim(state))) %>%
  group_by(placename, state) %>%
  summarise(avg_pm25 = mean(pm25, na.rm = TRUE),
            avg_svi = mean(svi_overall, na.rm = TRUE),
            total_pop = sum(population, na.rm = TRUE), .groups = "drop") %>%
  mutate(quadrant = {
    med_pm25 <- median(avg_pm25, na.rm = TRUE)
    med_svi <- median(avg_svi, na.rm = TRUE)
    case_when(
      avg_pm25 >= med_pm25 & avg_svi >= med_svi ~ "High Pollution / High Vulnerability",
      avg_pm25 >= med_pm25 & avg_svi < med_svi  ~ "High Pollution / Low Vulnerability",
      avg_pm25 < med_pm25  & avg_svi >= med_svi ~ "Low Pollution / High Vulnerability",
      avg_pm25 < med_pm25  & avg_svi < med_svi  ~ "Low Pollution / Low Vulnerability"
    )
  })

# 2. Define Case Study Cities for Labeling
target_cities <- c("Camden", "San Bernardino", "Fresno", "Trenton", "Everett",
                   "New York", "Naperville", "Boulder", "Fishers", "Clovis", "Mount Pleasant")
case_studies_plot <- city_summary %>% filter(placename %in% target_cities)

med_pm25_threshold <- median(city_summary$avg_pm25, na.rm = TRUE)
med_svi_threshold  <- median(city_summary$avg_svi, na.rm = TRUE)

# 3. Build the Main Plot Object
p_national <- ggplot(city_summary, aes(x = avg_pm25, y = avg_svi)) +
  # Background Quadrant Shading
  annotate("rect", xmin = med_pm25_threshold, xmax = Inf, ymin = med_svi_threshold, ymax = Inf,
           fill = "#d73027", alpha = 0.04) +
  annotate("rect", xmin = 0, xmax = med_pm25_threshold, ymin = 0, ymax = med_svi_threshold,
           fill = "#084594", alpha = 0.04) +
  # Data Points
  geom_point(aes(size = total_pop, color = quadrant), alpha = 0.45) +
  # Highlight Case Studies
  geom_point(data = case_studies_plot, aes(size = total_pop),
             shape = 21, color = "black", fill = NA, stroke = 1.2) +
  geom_label_repel(data = case_studies_plot,
                   aes(label = paste0(placename, ", ", state)),
                   size = 3.8, fontface = "bold", box.padding = 0.6,
                   segment.color = "black", max.overlaps = 35) +
  # Threshold Lines
  geom_vline(xintercept = med_pm25_threshold, linetype = "dotted") +
  geom_hline(yintercept = med_svi_threshold, linetype = "dotted") +
  # Scales & Styling
  scale_color_manual(values = c(
    "High Pollution / High Vulnerability" = "#d73027",
    "High Pollution / Low Vulnerability"  = "#f46d43",
    "Low Pollution / High Vulnerability"  = "#b3c5f4",
    "Low Pollution / Low Vulnerability"   = "#084594"),
    name = "EJ Quadrant") +
  scale_size_continuous(range = c(2, 16), labels = comma, name = "Total Population") +
  labs(x = expression(paste("Avg. PM2.5 (", mu, "g/m"^3, ")")),
       y = "Avg. SVI (0-1)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(face = "bold"))

# 4. Apply the Layout Fix (The "No-Cutoff" Strategy)
p_national_no_legend <- p_national + theme(legend.position = "none")

# Extract legend specifically
# get_legend() was removed in cowplot >= 1.1.3; extract legend via ggplot2
# by building a legend-only plot and grabbing its grob safely.
legend_shared <- tryCatch(
  cowplot::get_legend(
    p_national +
      theme(legend.position = "bottom",
            legend.box = "vertical",
            legend.spacing.y = unit(0.3, "cm"),
            legend.text = element_text(size = 10))
  ),
  error = function(e) {
    # Fallback: use ggplotGrob to extract the legend grob
    g <- ggplotGrob(p_national +
                      theme(legend.position = "bottom",
                            legend.box = "vertical"))
    leg_idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    if (length(leg_idx)) g$grobs[[leg_idx[1]]] else grid::nullGrob()
  }
)

# Assemble with Patchwork/Cowplot
# Use cowplot::plot_grid throughout -- it handles both ggplot objects and
# grobs (like legend_shared). patchwork wrap_plots() does not accept grobs.
title_grob <- grid::textGrob(
  "National Environmental Justice Landscape
Aggregated City-Level Pollution Exposure vs. Social Vulnerability",
  gp = grid::gpar(fontsize = 16, fontface = "bold"), just = "center"
)
final_plot <- cowplot::plot_grid(
  title_grob,
  p_national_no_legend,
  legend_shared,
  ncol = 1,
  rel_heights = c(0.08, 1, 0.2)
)

# 5. Save Output
ggsave(file.path(results_path, "Figures", "National_EJ_Quadrant_Final.pdf"),
       final_plot, width = 14, height = 12, dpi = 300)

cat("SUCCESS: National EJ Quadrant saved without legend cutoffs.\n")

# --- SECTION 15: CLEANED CASE STUDIES ---
# NOTE: city_combined, city_palette, city_levels, and health_matrix must be
# constructed from your case study city data before this section runs.
# Wrapped in tryCatch so a missing object does not crash the overnight run.
tryCatch({
  if (!exists("city_combined") || !exists("city_palette") ||
      !exists("city_levels")   || !exists("health_matrix")) {
    stop("city_combined / city_palette / city_levels / health_matrix are not defined. ",
         "Build these objects from your case study subset before running Section 15.")
  }
  
  tiff(file.path(results_path, "Figures", "Case_Studies_Comparison_Final.tiff"),
       width = 12, height = 13, units = "in", res = 300)
  
  layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE), heights=c(4.5, 4.5, 1.2))
  
  # Panel A: PM2.5
  boxplot(pm25 ~ City, data = city_combined, main = "PM2.5 Exposure",
          ylab = expression(PM[2.5] ~ (mu*g/m^3)), col = city_palette[city_levels])
  mtext("A", side = 3, adj = -0.1, line = 1, cex = 2, font = 2)
  
  # Panel B: SVI
  boxplot(svi_overall ~ City, data = city_combined, main = "Social Vulnerability Index",
          ylab = "SVI (0-1)", col = city_palette[city_levels])
  mtext("B", side = 3, adj = -0.1, line = 1, cex = 2, font = 2)
  
  # Panel C: EJ Profile
  plot(city_combined$pm25, city_combined$svi_overall,
       col = city_palette[as.character(city_combined$City)],
       main = "Local Environmental Justice Profile",
       xlab = expression(PM[2.5] ~ (mu*g/m^3)), ylab = "Social Vulnerability Index")
  mtext("C", side = 3, adj = -0.1, line = 1, cex = 2, font = 2)
  
  # Panel D: Health
  barplot(health_matrix, beside = TRUE, col = city_palette[city_levels],
          main = "Chronic Disease Prevalence", ylab = "Mean Prevalence (%)")
  mtext("D", side = 3, adj = -0.1, line = 1, cex = 2, font = 2)
  
  dev.off()
  cat("Section 15 case study figure saved.\n")
}, error = function(e) {
  if (!is.null(dev.list())) dev.off()
  cat("WARNING: Section 15 skipped --", e$message, "\n")
})
# =============================================================================
# SECTION 16: SPATIAL SENSITIVITY (THE ROBUSTNESS CHECK)
# =============================================================================
print_section("SECTION 16: SPATIAL SENSITIVITY")

# We prove that adding a spatial smooth s(lon, lat) successfully removed 
# spatial autocorrelation from the residuals across all SVI themes.

SVI_TARGETS <- c("svi_overall", "rpl_theme1", "rpl_theme2", "rpl_theme3", "rpl_theme4")
pop_col <- if("pop_density_scaled" %in% names(data_model)) "pop_density_scaled" else "pop_density"

sens_results_list <- list()

for (target in SVI_TARGETS) {
  cat(sprintf("Testing Spatial Robustness for: %s...\n", target))
  
  # Model A: Baseline (No spatial control)
  f_A <- as.formula(paste("pm25 ~", target, "+", pop_col, "+ time_scaled"))
  
  # Model C: Spatial (With neighborhood-scale smooth)
  f_C <- as.formula(paste("pm25 ~", target, "+", pop_col, "+ time_scaled + s(longitude, latitude, k=30)"))
  
  # Fit models using bam() for speed
  m_A <- bam(f_A, data = data_model, family = gaussian(link = "log"))
  m_C <- bam(f_C, data = data_model, family = gaussian(link = "log"))
  
  # Calculate Flatness Ratios (Proving the variogram is now "flat" or "white noise")
  v_A <- flatness_ratio(compute_variogram(m_A, data_model))
  v_C <- flatness_ratio(compute_variogram(m_C, data_model))
  
  # Clean Metric Name for Table
  clean_name <- str_to_title(str_replace_all(target, "rpl_|_", " "))
  
  sens_results_list[[target]] <- data.frame(
    SVI_Metric = clean_name,
    NonSpatial_Flatness = round(v_A, 3),
    Spatial_Flatness = round(v_C, 3),
    Autocorr_Removed = ifelse(v_C > 0.95, "COMPLETE", "PARTIAL")
  )
  
  rm(m_A, m_C); gc() # Clear memory
}

# --- LABELED VARIOGRAM PROOF VISUAL ---
# We use the Overall SVI model to create the visual proof for the defense
target <- "svi_overall"
f_A <- as.formula(paste("pm25 ~", target, "+", pop_col, "+ time_scaled"))
f_C <- as.formula(paste("pm25 ~", target, "+", pop_col, "+ time_scaled + s(longitude, latitude, k=30)"))

m_A <- bam(f_A, data = data_model, family = gaussian(link = "log"))
m_C <- bam(f_C, data = data_model, family = gaussian(link = "log"))

pdf(file.path(results_path, "Sensitivity", "Spatial_Correction_Proof.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 3))

# Panel A: Baseline Model
# gstat's plot.variogram() opens its own device and doesn't play nicely with
# mtext(). Extract the values and use base R plot() instead so mtext() works.
v_proof_A <- compute_variogram(m_A, data_model)
plot(v_proof_A$dist, v_proof_A$gamma,
     type = "b", pch = 19, col = "#2c3e50",
     main = "Non-Spatial Model (Clustered Errors)",
     xlab = "Distance", ylab = "Semivariance")
mtext("A", side = 3, adj = -0.1, line = 1.5, cex = 2, font = 2)

# Panel B: Spatial Model
v_proof_C <- compute_variogram(m_C, data_model)
plot(v_proof_C$dist, v_proof_C$gamma,
     type = "b", pch = 19, col = "#27ae60",
     main = "Spatial Smooth Model (White Noise/Flat)",
     xlab = "Distance", ylab = "Semivariance")
mtext("B", side = 3, adj = -0.1, line = 1.5, cex = 2, font = 2)

dev.off()

# Save the final proof table with cleaned column names
spatial_robustness_table <- bind_rows(sens_results_list)
colnames(spatial_robustness_table) <- c("SVI Metric", "Baseline Flatness", "Spatial Flatness", "Status")

write.csv(spatial_robustness_table, 
          file.path(results_path, "Tables", "Spatial_Sensitivity_Summary.csv"), 
          row.names = FALSE)

print(spatial_robustness_table)

cat("\n--- END OF ANALYSIS SCRIPT ---\n")

# =============================================================================
# SECTION 17: CONSOLIDATED THESIS TABLES (RESULTS EXPORT)
# =============================================================================
cat("\n--- GENERATING MASTER RESULTS TABLES ---\n")

# Helper function to clean names (e.g., "svi_overall" -> "Social Vulnerability Index")
clean_thesis_names <- function(x) {
  x <- str_replace_all(x, "svi_overall", "Social Vulnerability Index")
  x <- str_replace_all(x, "rpl_theme", "SVI Theme ")
  x <- str_replace_all(x, "_prev", "")
  x <- str_replace_all(x, "_", " ")
  return(str_to_title(x))
}

# --- MASTER TABLE 1: POLLUTION EXPOSURE (RQ1 & RQ2) ---
# Pulling results from both the Thematic and Scale Decomposition models
m1a <- readRDS(file.path(results_path, "Models", "model1a_thematic.rds"))
m1b <- readRDS(file.path(results_path, "Models", "model1b_scale_decomposition.rds"))

table1_results <- bind_rows(
  as.data.frame(summary(m1a)$s.table) %>% mutate(Model = "Thematic (National)", Term = rownames(.)),
  as.data.frame(summary(m1b)$s.table) %>% mutate(Model = "Scale (Within/Between)", Term = rownames(.))
) %>%
  mutate(Term = clean_thesis_names(Term)) %>%
  select(Model, Term, EDF = edf, `P-Value` = `p-value`)

write.csv(table1_results, file.path(results_path, "Tables", "Master_Results_Model1_Exposure.csv"), row.names = FALSE)


# --- MASTER TABLE 2: HEALTH INTERACTIONS (RQ3) ---
# Pulling the interaction effects for every disease and vulnerability theme
m2_list <- list()

for (v_name in names(SVI_MODEL2_VARIANTS)) {
  for (disease in DISEASES) {
    path <- file.path(results_path, "Models", sprintf("model2_%s_%s.rds", disease, v_name))
    if (file.exists(path)) {
      m <- readRDS(path)
      # Extract only the Interaction Term: ti(pm25, SVI)
      s_table <- as.data.frame(summary(m)$s.table)
      int_term <- s_table[grepl("ti\\(", rownames(s_table)), ]
      
      if (nrow(int_term) > 0) {
        m2_list[[paste(disease, v_name)]] <- data.frame(
          Outcome = clean_thesis_names(disease),
          Vulnerability_Metric = clean_thesis_names(v_name),
          Interaction_EDF = round(int_term$edf, 2),
          P_Value = format.pval(int_term$`p-value`, digits = 3, eps = 0.001),
          Model_R2 = round(summary(m)$r.sq, 3)
        )
      }
    }
  }
}

table2_results <- bind_rows(m2_list)
write.csv(table2_results, file.path(results_path, "Tables", "Master_Results_Model2_Interactions.csv"), row.names = FALSE)

cat("Done! Master tables saved in the /Tables folder.\n")
# =============================================================================
# SECTION 18: PUBLICATION-READY TABLES (GT FORMATTING)
# =============================================================================
library(gt)
library(gtsummary)

cat("\n--- FORMATTING TABLES FOR PUBLICATION ---\n")

# --- 1. Master Table: Environmental Justice Exposure (Model 1) ---
# This table combines the national thematic effects and neighborhood scale effects
tab1_pub <- table1_results %>%
  gt(groupname_col = "Model") %>%
  tab_header(
    title = "Table 1. Associations Between Social Vulnerability and PM2.5 Exposure",
    subtitle = "Results from Generalized Additive Models (GAMs) with Spatial Controls"
  ) %>%
  cols_label(
    Term = "Vulnerability Predictor",
    EDF = "Est. Degrees of Freedom (EDF)",
    `P-Value` = "p-value"
  ) %>%
  fmt_number(columns = EDF, decimals = 2) %>%
  # Apply significance stars logic
  # `P-Value` is character after select(); compare numerically via a helper col
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `P-Value`,
      rows    = as.numeric(`P-Value`) < 0.05
    )
  ) %>%
  tab_source_note(source_note = "Note: Bolded rows indicate statistical significance at p < 0.05.")

# --- 2. Master Table: Health Interaction Effects (Model 2) ---
tab2_pub <- table2_results %>%
  gt(groupname_col = "Outcome") %>%
  tab_header(
    title = "Table 2. Interaction Effects of PM2.5 and SVI on Chronic Disease",
    subtitle = "Thematic Analysis of Vulnerability Amplification"
  ) %>%
  cols_label(
    Vulnerability_Metric = "SVI Dimension",
    Interaction_EDF = "Interaction Complexity (EDF)",
    P_Value = "Significance (p)",
    Model_R2 = "Model R-squared"
  ) %>%
  tab_options(row_group.font.weight = "bold") %>%
  tab_style(
    style = cell_fill(color = "gray95"),
    locations = cells_body(rows = seq(1, nrow(table2_results), by = 2))
  )

# --- SAVE AS WORD DOCUMENTS ---
gtsave(tab1_pub, file.path(results_path, "Tables", "Pub_Table1_Exposure.docx"))
gtsave(tab2_pub, file.path(results_path, "Tables", "Pub_Table2_Interactions.docx"))

cat("Publication-ready Word documents saved in /Tables.\n")