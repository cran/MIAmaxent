## ----setup, include = FALSE---------------------------------------------------
# Check if required packages are available
has_suggested <- requireNamespace("sf", quietly = TRUE) &&
  requireNamespace("disdat", quietly = TRUE) &&
  requireNamespace("purrr", quietly = TRUE) &&
  requireNamespace("tidyr", quietly = TRUE) &&
  requireNamespace("tibble", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  out.width = "100%",
  eval = has_suggested
)

## ----check-packages, echo = FALSE, results = "asis", eval = TRUE--------------
if (!has_suggested) {
  cat("**Note:** This vignette requires the following packages to run: `sf`, `disdat`, `purrr`, `tidyr`, `tibble`, and `ggplot2`. Please install them to view the full vignette output.\n\n")
  cat("```r\n")
  cat('install.packages(c("sf", "disdat", "purrr", "tidyr", "tibble", "ggplot2"))\n')
  cat("```\n\n")
}

## ----libraries, message=FALSE-------------------------------------------------
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)
library(terra)
library(sf)
library(disdat)
library(MIAmaxent)

## ----load-data----------------------------------------------------------------
# Choose a species (change this to model a different species!)
selected_species <- "swi01" # `disdat` has anonymized species codes: swi01 to swi30

# Define predictor variables
predictors <- c(
  "bcc", "calc", "ccc", "ddeg", "nutri", "pday",
  "precyy", "sfroyy", "slope", "sradyy", "swb", "tavecc", "topo"
)

# Extract presence-only (PO) data for the selected species
swi_po <- disPo("SWI") # function from `disdat` package
po_records <- swi_po[swi_po$spid == selected_species, ]
cat("Number of presence records:", nrow(po_records), "\n")

# Background data (10,000 random locations)
background <- disBg("SWI") # function from `disdat` package

# Combine presence and background into training data frame
# RV = 1 for presences, NA for background
train_data <- rbind(
  data.frame(RV = 1, po_records[, predictors]),
  data.frame(RV = NA, background[, predictors])
)

cat("Training presences:", sum(train_data$RV == 1, na.rm = TRUE), "\n")

## ----plot-po, fig.height=4----------------------------------------------------
studyarea <- disBorder("SWI")
po_sf <- st_as_sf(po_records, coords = c("x", "y"), crs = 21781)
ggplot() +
  geom_sf(data = studyarea, fill = "lightgrey") +
  geom_sf(data = po_sf, color = "blue", size = 0.5) +
  labs(
    title = paste("Presence records for species", selected_species),
    x = "Easting (m)", y = "Northing (m)"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
# Rename columns to more informative names
col_name_map <- c(
  RV = "RV",
  bcc = "broadleaf.cover",
  calc = "bedrock.calcareous",
  ccc = "conifer.cover",
  ddeg = "growing.degree.days",
  nutri = "soil.nutrients",
  pday = "precip.days",
  precyy = "annual.precip",
  sfroyy = "summer.frost.freq",
  slope = "slope",
  sradyy = "solar.radiation",
  swb = "water.balance",
  tavecc = "temp.coldest.month",
  topo = "topographic.position"
)
names(train_data) <- col_name_map[names(train_data)]

# Convert the categorical variable to factor
train_data$bedrock.calcareous <- as.factor(train_data$bedrock.calcareous)

## ----single-fop, fig.height=4-------------------------------------------------
plotFOP(train_data, "temp.coldest.month")

## ----fop-data-extraction, echo=FALSE, fig.show='hide'-------------------------
# Identify continuous variables (exclude RV and categorical variables)
continuous_vars <- setdiff(names(train_data), c("RV", "bedrock.calcareous"))

# Extract FOP data for all continuous variables
fop_results <- map(continuous_vars, \(x) plotFOP(train_data, x))
names(fop_results) <- continuous_vars

## ----fop-combine-data, echo=FALSE---------------------------------------------
# Combine FOP data into a single data frame
fop_combined <- map2_dfr(
  fop_results, continuous_vars,
  \(x, y) mutate(x$FOPdata, variable = y)
)

# Prepare raw training data for density plots
train_data_long <- train_data %>%
  select(all_of(continuous_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

## ----fop-faceted-plot, fig.height=10, echo=FALSE------------------------------
# Create a faceted plot to compare FOP across all continuous variables
ggplot() +
  # Grey density showing data distribution (scaled to fit on same axis)
  geom_density(
    data = train_data_long, aes(x = value, y = after_stat(scaled)),
    fill = "grey90", color = NA
  ) +
  # FOP points and loess line
  geom_point(data = fop_combined, aes(x = intEV, y = intRV), size = 1.5) +
  geom_line(data = fop_combined, aes(x = intEV, y = loess), color = "red") +
  facet_wrap(~variable, scales = "free_x", ncol = 3) +
  labs(
    x = "Environmental variable value",
    y = "Frequency of Observed Presence (FOP)"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 8))

## ----derive-vars--------------------------------------------------------------
# Create derived variables using multiple transformation types
# L = Linear, M = Monotonic, D = Deviation (unimodal),
# HF = Forward hinge, HR = Reverse hinge, T = Threshold
train_DVs <- deriveVars(train_data,
  transformtype = c("L", "M", "D", "HF", "HR", "T"),
  allsplines = TRUE,
  quiet = TRUE
)

## ----examine-dvs--------------------------------------------------------------
# Example: some of the DVs created for temperature
head(names(train_DVs$dvdata$temp.coldest.month))

## ----select-dvs---------------------------------------------------------------
# Select the best DVs for each individual EV
train_DV_select <- selectDVforEV(train_DVs$dvdata,
  alpha = 0.001,
  quiet = TRUE
)

## ----examine-selected-dvs-----------------------------------------------------
# Example: DVs selected for temperature
names(train_DV_select$dvdata$temp.coldest.month)

## ----select-evs---------------------------------------------------------------
# Select which EVs (represented by their best DVs) to include
train_model <- selectEV(train_DV_select$dvdata,
  alpha = 0.001,
  interaction = FALSE,
  quiet = TRUE
)

# View the selection trail, only the top model from each round
selection_trail <- train_model$selection[!duplicated(train_model$selection$round), ]
knitr::kable(selection_trail, digits = c(0, 0, 0, 3, 1, 0, 3)) # Print the data.frame as a table

## ----final-model--------------------------------------------------------------
# View the final model
train_model$selectedmodel$formula

## ----response-curves, fig.height=4--------------------------------------------
# Get the most important variable (first selected)
top_EV <- train_model$selection$variables[1]

# Plot response curve for the most important variable
plotResp(
  train_model$selectedmodel,
  train_DVs$transformations,
  top_EV
)

## ----compare-fop-response, fig.show='hold', fig.height=4----------------------
plotFOP(train_data, top_EV)

## ----all-responses, include=FALSE---------------------------------------------
evs <- (train_model$selectedmodel$coefficients)[-1] |>
  names() |>
  sub("_.*$", "", x = _) |>
  unique()
for (ev in evs[-1]) {
  plotResp(
    train_model$selectedmodel,
    train_DVs$transformations,
    ev
  )
}
train_model$selectedmodel$formula

## ----load-data-pa-------------------------------------------------------------
# Extract presence-absence (PA) data
swi_pa <- disPa("SWI") # function from `disdat` package

# Merge with environmental data by location
swi_env <- disEnv("SWI") # function from `disdat` package
pa_merged <- merge(swi_pa, swi_env)

# Create validation data frame with RV and environmental variables (parallel to train_data)
validation_data <- data.frame(
  RV = pa_merged[, selected_species],
  pa_merged[, predictors]
)

# Rename columns to match training data
names(validation_data) <- col_name_map[names(validation_data)]

# Convert the categorical variable to factor
validation_data$bedrock.calcareous <- as.factor(validation_data$bedrock.calcareous)

cat("Validation presences:", sum(validation_data$RV == 1), "\n")
cat("Validation absences:", sum(validation_data$RV == 0), "\n")

## ----plot-pa, fig.height=4----------------------------------------------------
pa_sf <- st_as_sf(swi_pa, coords = c("x", "y"), crs = 21781)
# Convert to factor for discrete colors (presence/absence)
pa_sf[[selected_species]] <- as.factor(pa_sf[[selected_species]])
ggplot() +
  geom_sf(data = studyarea, fill = "lightgrey") +
  geom_sf(data = pa_sf, aes(color = .data[[selected_species]]), size = 0.5) +
  scale_color_manual(
    values = c("0" = "red", "1" = "blue"),
    labels = c("0" = "Absence", "1" = "Presence"),
    name = "Record type"
  ) +
  labs(
    title = paste("Presence-absence records for species", selected_species),
    x = "Easting (m)", y = "Northing (m)"
  ) +
  theme_minimal()

## ----validation-auc, fig.height=4---------------------------------------------
# Calculate validation AUC
auc_result <- testAUC(
  model = train_model$selectedmodel,
  transformations = train_DVs$transformations,
  data = validation_data
)

## ----reconstruct-spatial, include=FALSE---------------------------------------
# To create spatial predictions directly, we need raster data
# Here we reconstruct the (incomplete) raster from coordinates, and metadata instead

# Create projection data  (similar to validation data but with coordinates)
proj_data <- data.frame(
  x = pa_merged$x,
  y = pa_merged$y,
  pa_merged[, predictors]
)

# Rename columns to match training data
col_name_map <- c(
  x = "x",
  y = "y",
  col_name_map
)
names(proj_data) <- col_name_map[names(proj_data)]

# Convert the categorical variable to factor
proj_data$bedrock.calcareous <- as.factor(proj_data$bedrock.calcareous)

# SWI
# Prepared by Niklaus E. Zimmermann and Antoine Guisan
# 13 variables:  see SWI\01_metadata_SWI_environment.csv
# Coordinate reference system: Transverse, spheroid Bessel (note all data has had a constant shift applied to it)
# EPSG:21781
# Units: meters
# Raster cell size: 100m

# Create a raster to fill in with predictions
# Here we assume env_data has columns 'x' and 'y' for coordinates
coordinates <- proj_data[, c("x", "y")]
env_raster <- rast(
  xmin = min(coordinates$x),
  xmax = max(coordinates$x),
  ymin = min(coordinates$y),
  ymax = max(coordinates$y),
  resolution = 100,
  crs = "EPSG:21781"
)

# Create a simple prediction map using the environmental data
predictions <- projectModel(
  model = train_model$selectedmodel,
  transformations = train_DVs$transformations,
  data = proj_data
)

# Rasterize the predictions
pred_raster <- rasterize(as.matrix(coordinates), env_raster,
  values = predictions$output$PRO
)
plot(pred_raster, main = "Predicted Probability of Occurrence")

## ----auc-po, include=FALSE----------------------------------------------------
testAUC(
  model = train_model$selectedmodel,
  transformations = train_DVs$transformations,
  data = train_data
)

## ----bv-setup-----------------------------------------------------------------
# Create a tibble to store results for three alpha levels
model_comparison <- tibble(
  alpha = c(1e-1, 1e-3, 1e-9)
)

## ----bv-split-----------------------------------------------------------------
set.seed(123) # For reproducibility

# Create indices for 70/30 split
n_total <- nrow(train_data)
train_idx <- sample(1:n_total, size = floor(0.7 * n_total))

# Split the data into training and test sets
train_70 <- train_data[train_idx, ]
test_30 <- train_data[-train_idx, ]

# Derive variables for the training set
train_DVs_70 <- deriveVars(train_70,
  transformtype = c("L", "M", "D", "HF", "HR", "T"),
  allsplines = TRUE,
  quiet = TRUE
)

## ----bv-fit-models------------------------------------------------------------
# Add model fitting results using map
model_comparison <- model_comparison %>%
  mutate(
    # Step 1: Select DVs for each EV
    dv_selection = map(alpha, \(a) selectDVforEV(train_DVs_70$dvdata,
      alpha = a,
      quiet = TRUE
    )),
    # Step 2: Select EVs for final model
    ev_selection = map2(
      dv_selection, alpha,
      \(dv, a) selectEV(dv$dvdata,
        alpha = a,
        interaction = FALSE,
        quiet = TRUE
      )
    ),
    # Extract number of DVs in final model
    n_dvs = map_int(ev_selection, \(ev) length(ev$selectedmodel$coefficients) - 1)
  )

# Display intermediate results
model_comparison %>%
  select(alpha, n_dvs)

## ----bv-evaluate, fig.show='hide', warning=FALSE------------------------------
# Add AUC calculations
model_comparison <- model_comparison %>%
  mutate(
    # AUC on 30% test set (internal evaluation)
    auc_test = map_dbl(
      ev_selection,
      \(ev) testAUC(
        model = ev$selectedmodel,
        transformations = train_DVs_70$transformations,
        data = test_30
      )
    ),
    # AUC on independent validation set (external evaluation)
    auc_validation = map_dbl(
      ev_selection,
      \(ev) testAUC(
        model = ev$selectedmodel,
        transformations = train_DVs_70$transformations,
        data = validation_data
      )
    )
  )

# Display final results table
results_table <- model_comparison %>%
  select(alpha, n_dvs, auc_test, auc_validation)

knitr::kable(results_table, digits = c(10, 0, 3, 3))

