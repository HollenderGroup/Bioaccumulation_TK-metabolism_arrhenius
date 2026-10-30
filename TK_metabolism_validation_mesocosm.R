###############################################################################
###############################################################################
##
##    Temperature corrected TK-metabolims (TK-metabolism-temp)
##
##   VALIDATION OF THE CALIBRATED TK-METABOLISM MODEL AGAINST THE MESOCOSM
##   STUDY (MesoTK1) WITH ARRHENIUS TEMPERATURE CORRECTION
##
##   parent compound (diflufenican)  --->  BTP (AE B107137)
##
##   Companion script to TK_metabolism_model.R (calibration) and
##   TK_arrhenius.R (Arrhenius analysis of the calibrated rate constants).
##   Johannes Raths, 2026. R script: October 2026.
##
###############################################################################
###############################################################################
##
##  WHAT THIS SCRIPT DOES
##  ---------------------
##  The rate constants calibrated in the laboratory (AE0 + AE1 pooled, full TK
##  design, 20.25 degC) are transferred to the mesocosm, where the water
##  temperature changes over the day (logger every 5 min). Every rate constant
##  follows the measured water temperature T(t) with the Arrhenius equation
##
##        k(T) = A * exp(-T_A / T)                T in K, T_A = Ea/R in K
##
##  The model is then run BLIND (nothing is fitted) with the measured water
##  concentrations as exposure, and compared with the measured tissue
##  concentrations of diflufenican and AE B107137 in A. aquaticus.
##
##  Which T_A is used for which rate constant (SECTION 1.6):
##     ku               : T_A of ku
##     ke, km and kem   : T_A of k1 = ke + km (total parent loss), because the
##                        T_A of ke, km and kem could not be estimated reliably
##
##  TWO MODEL VERSIONS (switch in SECTION 1.7; run one or both)
##  ------------------------------------------------------------
##  "reference" : the AE0+AE1 rate constants are kept exactly and only scaled
##                with temperature:
##                    k(T) = k_ref * exp(-T_A * (1/T - 1/T_ref))
##                (= the equation above with A = k_ref * exp(T_A / T_ref)).
##                At 20.25 degC the model is identical to the calibration.
##  "intercept" : the Arrhenius regression lines of ALL calibrations
##                (AE0+AE1, AE2, AE4) are used, i.e. A AND T_A come from the
##                Arrhenius fit; no reference temperature:
##                    k(T) = A * exp(-T_A / T)
##                ku     : the weighted Arrhenius line of ku (T_A of ku)
##                ke, km : the weighted Arrhenius line of k1 (T_A of k1),
##                         split into ke and km with the AE0+AE1 ratios
##                         ke/k1 and km/k1
##                kem    : a weighted Arrhenius line through the kem
##                         estimates with the slope fixed to T_A of k1
##                (A is calculated from the rate-constant file with the same
##                 weights as in TK_arrhenius.R, w = 1/SE^2 of ln k; for ku and
##                 k1 this reproduces lnA in Arrhenius_results.xlsx exactly.)
##  "none"      : optional reference run WITHOUT temperature correction
##                (AE0+AE1 rate constants at all temperatures); useful as a
##                baseline for the AIC comparison.
##
##  THE MODEL (identical to TK_metabolism_model.R / BYOM derivatives.m)
##  --------------------------------------------------------------------
##      dCp/dt = ku(T) * Cw_p(t)  -  ke(T) * Cp  -  km(T) * Cp
##      dCm/dt = km(T) * Cp       -  kem(T) * Cm
##   Cp, Cm  internal concentration of parent and BTP     [umol/kg ww]
##   Cw_p    parent concentration in the water             [umol/L]
##           mean of the replicates, linear interpolation between the
##           sampling times, constant after the last one (as in the calibration)
##   T(t)    water temperature of the logger, linear interpolation   [K]
##   Uptake of the BTP from the water is not modelled (kum was not calibrated;
##   same as in the calibration).
##   All calculations are in molar units (mass balance parent -> BTP); plots
##   and exports are in ug/L and ug/kg ww (as measured), the export also in umol.
##
##   Solver: because the rate constants change with T(t), the analytical
##   solution of the calibration script cannot be used. The ODEs are solved
##   with a fixed-step Runge-Kutta (4th order) on a grid that contains every
##   logger time, every water and tissue sampling time and has steps of at
##   most STEP_MAX_D (default 0.01 d = 14 min; with k*dt < 0.15 the error is
##   < 1e-6). All Monte Carlo draws are solved at once (vectorised), which
##   makes the CI bands fast. Optionally the best estimate is checked with the
##   deSolve ODE solver (lsoda) (SECTION 1.9).
##
##  UNCERTAINTY (95% CI bands; three switches in SECTION 1.8)
##  ---------------------------------------------------------
##   Monte Carlo: every draw combines the sources that are switched on, the
##   model is solved for every draw, and the CI band is the 2.5-97.5% range.
##   Each version gets the uncertainty of ITS OWN parameter source:
##   "reference" (k = k_ref * exp(-T_A (1/T - 1/T_ref))):
##     CI_TA_SE          T_A(ku) and T_A(k1) ~ normal(mean, SE), independent;
##                       the lines rotate around T_ref
##     CI_REF_PARAMETERS k_ref = one parameter set (ku, ke, km, kem) drawn at
##                       random from the AE0+AE1 parameter sets of the 95%
##                       region (parameter-space explorer of the calibration).
##                       Whole sets are drawn, so the correlations (e.g. km and
##                       kem) are kept.
##   "intercept" (k = A * exp(-T_A / T), Arrhenius regression lines):
##     CI_TA_SE          slope: T_A ~ normal(mean, SE); the lines rotate around
##                       their weighted centroid (where the slope and the
##                       height of a weighted regression line are independent)
##     CI_LINE_HEIGHT    height of each line at its centroid ~ normal(0, SE),
##                       SE = sigma_hat * sqrt(1 / sum(w)) as in lm() of
##                       TK_arrhenius.R (the grey band of the Arrhenius figure).
##                       The AE0+AE1 CIs are NOT added here: they are already
##                       part of the regression (no double counting).
##   "none": only CI_REF_PARAMETERS (no temperature correction).
##   Not included: residual (animal-to-animal) variability, i.e. the bands are
##   confidence bands of the mean model, not prediction intervals.
##
##  QUALITY OF FIT / VALIDITY (per version, compound and treatment + pooled)
##  ------------------------------------------------------------------------
##   Only tissue samples of SPECIES (e.g. "A. aquaticus") are used; other
##   entries in the Species column (e.g. "A. aquaticus (dead)") are ignored.
##   R2     1 - SSQ/SST on the transformed scale (square root, as BYOM and the
##          calibration script). Can be negative in a validation (model worse
##          than the mean of the data).
##   MLL    (n/2) * log(SSQ) + n/2 (residual variance profiled out, as BYOM)
##   AIC    2 * MLL + 2 * k with k = 0 (nothing is fitted in a validation).
##          Only meaningful for COMPARING the model versions on the same data
##          (lower = better; see sheet "version_comparison").
##   TK analogues of the three validation criteria of EFSA PPR Panel (2018),
##   EFSA Journal 16(8):5377, section 4.1.4.5 (pp. 53-54), defined there for
##   GUTS survival models:
##   PPC    posterior predictive check: % of the observations inside the 95%
##          CI band of the model. EFSA allows confidence instead of credible
##          intervals ("PPC-like" check for frequentist fits); < 50% indicates
##          poor model performance (EFSA). Here the band contains parameter
##          uncertainty only (no residual error).
##   NRMSE  RMSE / mean of the observations x 100% (untransformed; same
##          definition as EFSA); should not exceed 50% (EFSA).
##   PE     prediction error at the end of the exposure (TK analogue of the
##          survival-probability prediction error SPPE of EFSA):
##          (predicted - observed) / observed x 100% at the LAST tissue
##          sampling time of each treatment (mean of the replicates;
##          positive = overprediction = conservative). Pooled = mean over the
##          treatments. EFSA gives NO formal cut-off for the SPPE; the +/-20%
##          used here (PE_MAX_ABS) follows the worked GUTS example of EFSA
##          (SPPE "within +/-20%" of the observed survival, section 4.1.4.5).
##
##  OUTPUT (folder OUTPUT_DIR in the working directory, one sub-folder per run)
##   <RUN_NAME>_<version>_figure*.png           figures per version and layout
##        _figure.png                   2 x 2 panels, treatments overlaid
##        _figure_treatment_rows.png    3 x 4 panels (rows = treatments)
##        _figure_treatment_columns.png 4 x 3 panels (columns = treatments)
##   <RUN_NAME>_validation_results.xlsx        quality of fit, rate constants,
##                                              modelled time courses, ...
##   <RUN_NAME>_log.txt                        settings and console output
##
##  HOW TO USE
##   1. Adjust SECTION 1 (files, switches, figure settings).
##   2. Run the whole script (Ctrl+Shift+Enter / source).
##   You do NOT need to change anything below SECTION 1.
##
##  Packages: ggplot2, patchwork, readxl, writexl (required); deSolve,
##            rstudioapi (optional)
##
###############################################################################


###############################################################################
###############################################################################
##                                                                           ##
##   SECTION 1 :  USER SETTINGS   (this is the only part you need to edit)   ##
##                                                                           ##
###############################################################################
###############################################################################

## ---------------------------------------------------------------------------
## 1.1  Files and folders
## ---------------------------------------------------------------------------
# Working directory: all input files are in this folder (use "/" not "\")
setwd("C:/Onedrive/OneDrive - Eawag/Offline Working/Johannes_Projects/Paper_007_DFF/01_Rscript_others/04_Mesocosm")

# File names in the working directory (or a full path).
# "" = choose the file with a file browser instead.
DATA_FILE          <- "raths_dff_SI_B01.xlsx"         # measured concentrations (long format, SI B01)
TEMPERATURE_FILE   <- "raths_dff_SI_B04.xlsx"         # temperature loggers of the mesocosm (SI B04)
RATE_CONSTANT_FILE <- "rate_constants_AE0-AE4.xlsx"   # calibrated rate constants
ARRHENIUS_FILE     <- "Arrhenius_results.xlsx"        # results of TK_arrhenius.R

DATA_SHEET    <- 1         # sheet (number or name) of DATA_FILE
LOGGER_SHEET  <- 1         # sheet of TEMPERATURE_FILE

OUTPUT_DIR    <- "TK_validation_output"   # relative name = sub-folder of the working directory
RUN_NAME      <- "DFF_MesoTK1"            # prefix of all output files


## ---------------------------------------------------------------------------
## 1.2  Measured concentrations: columns and data selection
## ---------------------------------------------------------------------------
COL_EXPERIMENT <- "Experiment code"
COL_MATRIX     <- "Matrix"
COL_SPECIES    <- "Species"
COL_TREATMENT  <- "Temperature [C] / Treatment"   # tolerant: also matches [°C]
COL_TIME       <- "Time [d]"                      # time since start, in DAYS
COL_REPLICATE  <- "Replicate"
COL_PARENT     <- "Diflufenican Concentration"
COL_BTP        <- "AE B107137 Concentration"

EXPERIMENT      <- "MesoTK1"           # experiment code of the mesocosm
EXPOSURE_MATRIX <- "Medium"            # water samples (ug/L)
ORGANISM_MATRIX <- "Tissue"            # tissue samples (ug/kg ww)
SPECIES         <- "A. aquaticus"   # tissue samples used (other species entries,
# e.g. "A. aquaticus (dead)", are ignored)

# Treatments: code in the data (and in LOGGER_COLUMNS) = label in the figure
TREATMENTS <- c(Amb = "Ambient (Amb)",
                E   = "Elevated (E)",
                HW  = "Heatwave (HW)")

# Single tissue values excluded from the ENTIRE analysis (figures, quality of
# fit, exports): treatment, time [d], replicate, compound ("parent", "btp" or
# "both"). Default: the two outliers of HW at 3.07 d (replicate 1: diflufenican
# 264 ug/kg; replicate 2: AE B107137 171 ug/kg). Empty table = keep everything.
EXCLUDE_TISSUE <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
  treatment  time       replicate  compound
  HW         3.066667   1          parent
  HW         3.066667   2          btp
")

# Values below the LOQ: text entries ("<LOQ", "n.d.", ...) are set to this value;
# empty cells and MISSING_CODES are missing (not used)
BELOW_LOQ_VALUE <- 0
MISSING_CODES   <- c("NA", "n.a.", "na", "-", "missing")


## ---------------------------------------------------------------------------
## 1.3  Temperature loggers
## ---------------------------------------------------------------------------
# Column NUMBERS (not names) of the exposure time [d] and the temperature
# [degC] of each treatment. In raths_dff_SI_B04.xlsx the header row is
# "#, Time, Exposure time [d], Ambient, Elevated, HW, Column1, Column2", but
# the content is:  3-5 = exposure time of Amb, E, HW (each relative to the
# start of its own mesocosm: 09:03, 09:13, 09:23) and 6-8 = temperature of
# Amb, E, HW.
LOGGER_COLUMNS <- list(Amb = c(time = 3, temp = 6),
                       E   = c(time = 4, temp = 7),
                       HW  = c(time = 5, temp = 8))
MAX_LOGGER_GAP_H <- 1      # warn if the logger has a gap longer than this (hours)


## ---------------------------------------------------------------------------
## 1.4  Compounds
## ---------------------------------------------------------------------------
PARENT_NAME       <- "Diflufenican"   # full names (figure titles, tables)
PARENT_MOLAR_MASS <- 394.29           # g/mol, C19H11F5N2O2
BTP_NAME          <- "AE B107137"
BTP_MOLAR_MASS    <- 283.20           # g/mol, C13H8F3NO3


## ---------------------------------------------------------------------------
## 1.5  Calibrated TK parameters (reference)
## ---------------------------------------------------------------------------
# Row of RATE_CONSTANT_FILE (column "model_run") with the reference parameters
# ku, ke, km, kem and its temperature T_C (= T_ref of version "reference").
REFERENCE_RUN <- "AE0_AE1_pooled"     # AE0 | 20.1 + AE1 | 20.4 -> 20.25 degC

# Initial internal concentrations [ug/kg ww] (MesoTK1 has no tissue data at t = 0)
X0_PARENT <- 0
X0_BTP    <- 0


## ---------------------------------------------------------------------------
## 1.6  Arrhenius temperatures T_A
## ---------------------------------------------------------------------------
# Which T_A is used for which rate constant:
TA_FOR <- c(ku = "ku", ke = "k1", km = "k1", kem = "k1")

# "file"   : T_A and SE from ARRHENIUS_FILE (sheet arrhenius_results:
#            ku 5686.9 +/- 946.8 K, k1 6382.5 +/- 1754.8 K)  (DEFAULT)
# "manual" : the values below (rounded values for cross-validation)
TA_SOURCE <- "file"
TA_MANUAL <- data.frame(rate = c("ku", "k1"),
                        T_A  = c(5690, 6380),     # K
                        SE   = c( 950, 1750))     # K


## ---------------------------------------------------------------------------
## 1.7  Model versions
## ---------------------------------------------------------------------------
# "reference" : AE0+AE1 rate constants scaled from T_ref with T_A
# "intercept" : Arrhenius regression lines (A and T_A), no reference temperature
# "none"      : no temperature correction (optional baseline for the AIC)
# One, two or all three, e.g. c("reference", "intercept")
MODEL_VERSIONS <- c("reference", "intercept")

# Version "intercept": which calibrated estimates define the line (A) of each
# rate constant (the slope is always T_A of TA_FOR):
#   the rate itself (e.g. kem = "kem") : weighted line through its own estimates
#   another rate (ke = "k1")           : line of that rate x AE0+AE1 ratio
#                                        (e.g. ke = k1_line * ke_ref / k1_ref)
#   "reference"                        : anchored at AE0+AE1 (as version "reference")
INTERCEPT_FROM <- c(ku = "ku", ke = "k1", km = "k1", kem = "kem")
# Calibrations used for these lines ("all" = every row of RATE_CONSTANT_FILE,
# or a vector of model_run names, e.g. c("AE0_AE1_pooled", "AE2 | 15.6", ...))
INTERCEPT_RUNS <- "all"


## ---------------------------------------------------------------------------
## 1.8  Uncertainty (Monte Carlo, CI bands)
## ---------------------------------------------------------------------------
# Sources of uncertainty in the CI bands (TRUE = included):
CI_TA_SE          <- TRUE   # SE of T_A (ku and k1)        versions "reference", "intercept"
CI_REF_PARAMETERS <- TRUE   # AE0+AE1 parameter sets       versions "reference", "none"
CI_LINE_HEIGHT    <- TRUE   # SE of the Arrhenius lines    version  "intercept"

# AE0+AE1 parameter sets (only for CI_REF_PARAMETERS): output of
# TK_metabolism_model.R for the pooled AE0+AE1 run (version with the
# "parameter_sample" sheet).
#   "auto" : the NEWEST file matching PARAMETER_SAMPLE_PATTERN anywhere in the
#            working directory or its sub-folders (e.g. TK_output_for_CIs/...)
#   or a file name / path, e.g. "TK_output_for_CIs/AE0_AE1_2026-10-18_230512/
#            DFF_AE0_AE1_pooled_parameters_quality_of_fit.xlsx"
PARAMETER_SAMPLE_FILE    <- "auto"
PARAMETER_SAMPLE_PATTERN <- "AE0_AE1_pooled_parameters_quality_of_fit\\.xlsx$"
PARAMETER_SAMPLE_SHEET <- "parameter_sample"

N_MC        <- 1000     # number of Monte Carlo draws (1000 is plenty; 200 for a quick look)
CI_LEVEL    <- 0.95     # width of the CI band
RANDOM_SEED <- 42       # reproducible draws
MC_CHUNK    <- 100      # draws solved together (progress bar step; memory)


## ---------------------------------------------------------------------------
## 1.9  Numerics
## ---------------------------------------------------------------------------
STEP_MAX_D         <- 0.01   # largest Runge-Kutta step [d] (logger times are
# always grid points, i.e. 5-min steps where logged)
OUTPUT_DT          <- 0.01   # time resolution of the exported model curves [d]
CHECK_WITH_DESOLVE <- TRUE   # compare the best estimate with deSolve::lsoda
# (only if the package is installed)


## ---------------------------------------------------------------------------
## 1.10  Quality-of-fit criteria
## ---------------------------------------------------------------------------
DATA_TRANSFORM <- 0.5   # for R2 / MLL / AIC: 0.5 = square root (as calibration),
# 0 = log, 1 = none
PPC_MIN        <- 50    # acceptable if >= PPC_MIN % of the observations in the CI band
NRMSE_MAX      <- 50    # acceptable if NRMSE <= NRMSE_MAX %
PE_MAX_ABS     <- 20    # acceptable if -PE_MAX_ABS <= PE <= +PE_MAX_ABS %
# (EFSA 2018 gives no cut-off for the SPPE; 20% follows
# the worked GUTS example, section 4.1.4.5)


## ---------------------------------------------------------------------------
## 1.11  Figure
## ---------------------------------------------------------------------------
# Treatment colours (as in the temperature-logger figure)
TREATMENT_COLOURS <- c(Amb = "#3B75AF", E = "#EF8636", HW = "#B4182D")
BAND_ALPHA  <- 0.22       # transparency of the CI bands
GRID_COL    <- "#E6E6E6"  # horizontal grid lines (as Arrhenius figure)
BASE_SIZE   <- 8          # font size [pt] (Arrhenius figure: pointsize 8)
# Layouts (one PNG each; one, two or all three):
#   "overview"          : 2 x 2 panels a-d, treatments overlaid
#   "treatment_rows"    : 3 x 4 panels, rows = treatments, columns = a-d
#   "treatment_columns" : 4 x 3 panels, rows = a-d, columns = treatments
FIGURE_LAYOUTS <- c("overview", "treatment_rows", "treatment_columns")
FIG_SIZE_MM <- list(overview          = c(width = 180, height = 140),   # 180 mm =
                    treatment_rows    = c(width = 180, height = 150),   # double column
                    treatment_columns = c(width = 180, height = 210))
# Axis ranges (NULL = automatic). Values outside a range stay in all
# calculations; points, lines and bands are only cut off at the panel edge.
AXIS_LIMITS <- list(water  = c(0, 1.1),    # water concentration [ug/L]
                    temp   = c(15, 35),    # water temperature [degC]
                    parent = c(0, 200),    # diflufenican in tissue [ug/kg ww]
                    btp    = c(0, 110))    # AE B107137 in tissue [ug/kg ww]
# Labelled (major) ticks, NULL = automatic; minor ticks lie halfway in between
AXIS_BREAKS <- list(water = seq(0, 1, by = 0.25),
                    btp   = seq(0, 100, by = 25))
# Version-specific ranges / ticks: replace the entry above for that version only
AXIS_LIMITS_VERSION <- list(reference = list(btp = c(0, 160)))
AXIS_BREAKS_VERSION <- list(reference = list(btp = seq(0, 150, by = 50)))
X_LIMITS <- c(Amb = 14, E = 10, HW = 3)    # time axis 0 - X_LIMITS [d] per treatment
# (overview layout: 0 - largest value)
MINOR_TICKS <- TRUE       # unlabelled ticks between the labelled ones (ggplot2 >= 3.5.0)
# AE B107137 in the water (dotted line, triangles; measured, not used by the model)
#   "auto": plotted unless it is below BTP_WATER_MIN_FRACTION of the parent
#           at every sampling time; TRUE / FALSE: always / never
PLOT_BTP_WATER         <- "auto"
POINT_ALPHA    <- 0.8     # transparency of the measured points (overlaps look darker)
COLOURED_STRIPS <- TRUE   # treatment names in their colour (needs package ggtext)
BTP_WATER_MIN_FRACTION <- 0.05
SHOW_QC_IN_PANELS <- TRUE # treatment layouts: NRMSE, PPC, PE of each treatment
# printed in the tissue panels
PNG_DPI     <- 600
SHOW_TREF_LINE <- TRUE    # dashed T_ref line in the temperature panel
# (version "reference" only)
CAPTION_WRAP <- 140       # characters per line of the figure caption


###############################################################################
###############################################################################
##                                                                           ##
##   NOTHING TO CHANGE BELOW THIS LINE                                       ##
##                                                                           ##
###############################################################################
###############################################################################


###############################################################################
## SECTION 2 : Packages and helper functions
###############################################################################

## 2.1 Packages -------------------------------------------------------------------
for (p in c("ggplot2", "patchwork", "readxl", "writexl")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Package '", p, "' is needed. Install it with: install.packages(\"", p, "\")")
}
library(ggplot2)
library(patchwork)
HAVE_DESOLVE <- requireNamespace("deSolve", quietly = TRUE)
if (MINOR_TICKS && utils::packageVersion("ggplot2") < "3.5.0") {
  message("Minor ticks need ggplot2 >= 3.5.0 (installed: ", utils::packageVersion("ggplot2"),
          "); figures without minor ticks. Update with install.packages(\"ggplot2\").")
  MINOR_TICKS <- FALSE
}

## 2.2 Console + log file -----------------------------------------------------------
LOG <- character(0)
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  LOG <<- c(LOG, msg)
  invisible(msg)
}
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(signif(x, d), format = "fg", digits = d, flag = "#"))

## 2.3 Choose a file (path or file browser) -------------------------------------------
choose_file <- function(path, caption) {
  if (nzchar(path)) {
    if (!file.exists(path))
      stop("File not found: ", normalizePath(path, mustWork = FALSE),
           "\n  Working directory is: ", getwd())
    return(normalizePath(path))
  }
  if (!interactive()) stop("No path given for: ", caption, " (and R is not interactive).")
  f <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    f <- tryCatch(rstudioapi::selectFile(caption = caption,
                                         filter = "Excel files (*.xlsx *.xls)"),
                  error = function(e) NULL)
  } else {
    cat(caption, "(the dialog may open behind other windows) ...\n")
    f <- file.choose()
  }
  if (is.null(f) || !nzchar(f)) stop("No file selected for: ", caption)
  normalizePath(f)
}

## 2.4 Tolerant column names, text to numbers ------------------------------------------
norm_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
find_col <- function(df, wanted) {
  nms <- names(df)
  i <- which(nms == wanted)
  if (!length(i)) i <- which(norm_name(nms) == norm_name(wanted))
  if (!length(i)) stop("Column '", wanted, "' not found. Available columns: ",
                       paste0("'", nms, "'", collapse = ", "), "\n  Adjust SECTION 1.2.")
  nms[i[1]]
}
to_num <- function(x) {
  s   <- trimws(as.character(x))
  out <- suppressWarnings(as.numeric(gsub(",", ".", s)))
  miss <- is.na(x) | s == "" | s %in% MISSING_CODES
  out[is.na(out) & !miss] <- BELOW_LOQ_VALUE        # "<LOQ" etc.
  out[miss] <- NA
  out
}

## 2.5 Small numeric helpers ------------------------------------------------------------
tr_data <- function(x) {                       # transformation for R2 / MLL
  if (DATA_TRANSFORM == 0) log(pmax(x, 1e-10)) else pmax(x, 0)^DATA_TRANSFORM
}
tw_mean <- function(t, y) {                    # time-weighted mean (trapezoid)
  if (length(t) < 2) return(y[1])
  sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2) / diff(range(t))
}
RATE_NAMES <- c("ku", "ke", "km", "kem")
RATE_UNITS <- c(ku = "L/kg ww/d", ke = "1/d", km = "1/d", kem = "1/d")
VERSION_LABELS <- c(reference = "Arrhenius anchored at the AE0+AE1 calibration (T_ref)",
                    intercept = "Arrhenius regression lines (A and T_A, no T_ref)",
                    none      = "No temperature correction (AE0+AE1 rate constants)")
STATE_NAMES <- c(parent = PARENT_NAME, btp = BTP_NAME)
MM          <- c(parent = PARENT_MOLAR_MASS, btp = BTP_MOLAR_MASS)

## 2.6 Check the settings ---------------------------------------------------------------
if (!length(MODEL_VERSIONS) || !all(MODEL_VERSIONS %in% names(VERSION_LABELS)))
  stop("MODEL_VERSIONS must be one or more of: ", paste(names(VERSION_LABELS), collapse = ", "))
MODEL_VERSIONS <- unique(MODEL_VERSIONS)
if (!all(RATE_NAMES %in% names(TA_FOR)))       stop("TA_FOR needs ku, ke, km and kem.")
if (!all(RATE_NAMES %in% names(INTERCEPT_FROM))) stop("INTERCEPT_FROM needs ku, ke, km and kem.")
if (!all(names(TREATMENTS) %in% names(LOGGER_COLUMNS)))
  stop("Every treatment in TREATMENTS needs an entry in LOGGER_COLUMNS.")
if (!all(names(TREATMENTS) %in% names(TREATMENT_COLOURS)))
  stop("Every treatment in TREATMENTS needs a colour in TREATMENT_COLOURS.")
if (N_MC < 10) stop("N_MC should be at least 10 (better 1000).")
if (!all(TR %in% names(X_LIMITS))) stop("X_LIMITS needs a value for every treatment.")
if (!all(c(names(AXIS_LIMITS_VERSION), names(AXIS_BREAKS_VERSION)) %in% names(VERSION_LABELS)))
  stop("AXIS_LIMITS_VERSION / AXIS_BREAKS_VERSION: names must be model versions (",
       paste(names(VERSION_LABELS), collapse = ", "), ").")
if (!all(c("water", "temp", "parent", "btp") %in% names(AXIS_LIMITS)))
  stop("AXIS_LIMITS needs the entries water, temp, parent and btp (NULL = automatic).")
if (!length(FIGURE_LAYOUTS) || !all(FIGURE_LAYOUTS %in% names(FIG_SIZE_MM)))
  stop("FIGURE_LAYOUTS must be one or more of: ", paste(names(FIG_SIZE_MM), collapse = ", "))
TR <- names(TREATMENTS)


###############################################################################
## SECTION 3 : Read the input files
###############################################################################

DATA_FILE          <- choose_file(DATA_FILE,          "Select the measured concentrations (SI B01)")
TEMPERATURE_FILE   <- choose_file(TEMPERATURE_FILE,   "Select the temperature loggers (SI B04)")
RATE_CONSTANT_FILE <- choose_file(RATE_CONSTANT_FILE, "Select the calibrated rate constants")
ARRHENIUS_FILE     <- if (TA_SOURCE == "file")
  choose_file(ARRHENIUS_FILE, "Select the Arrhenius results") else ARRHENIUS_FILE
USE_PSAMPLE <- CI_REF_PARAMETERS && any(c("reference", "none") %in% MODEL_VERSIONS)
if (USE_PSAMPLE && identical(PARAMETER_SAMPLE_FILE, "auto")) {
  hits <- list.files(getwd(), pattern = PARAMETER_SAMPLE_PATTERN, recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("^~\\$", basename(hits))]                 # skip Excel lock files
  if (!length(hits))
    stop("No file matching '", PARAMETER_SAMPLE_PATTERN, "' found in ", getwd(), " or its sub-folders.\n",
         "  Run TK_metabolism_model.R for the pooled AE0+AE1 model first, give the path in\n",
         "  PARAMETER_SAMPLE_FILE, or set CI_REF_PARAMETERS <- FALSE.")
  PARAMETER_SAMPLE_FILE <- hits[which.max(file.mtime(hits))]     # newest
  if (length(hits) > 1) message(length(hits), " parameter-sample files found; using the newest.")
}
PARAMETER_SAMPLE_FILE <- if (USE_PSAMPLE)
  choose_file(PARAMETER_SAMPLE_FILE, "Select the AE0+AE1 parameter sets (..._parameters_quality_of_fit.xlsx)") else
    PARAMETER_SAMPLE_FILE

# Output folder: <OUTPUT_DIR>/<RUN_NAME>_<date_time>
if (!grepl("^(/|~|[A-Za-z]:|\\\\\\\\)", OUTPUT_DIR))
  OUTPUT_DIR <- file.path(getwd(), OUTPUT_DIR)
OUT <- file.path(OUTPUT_DIR, paste0(RUN_NAME, "_", format(Sys.time(), "%Y-%m-%d_%H%M%S")))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
out_file <- function(...) file.path(OUT, paste0(RUN_NAME, "_", ...))

say("=====================================================================")
say(" TK model validation: ", EXPERIMENT, "  (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")")
say("=====================================================================")
say("Data          : ", DATA_FILE)
say("Temperatures  : ", TEMPERATURE_FILE)
say("Rate constants: ", RATE_CONSTANT_FILE)
if (TA_SOURCE == "file") say("Arrhenius     : ", ARRHENIUS_FILE)
if (USE_PSAMPLE)        say("Parameter sets: ", PARAMETER_SAMPLE_FILE)
say("Output        : ", OUT)

## 3.1 Measured concentrations -------------------------------------------------------
raw <- as.data.frame(readxl::read_excel(DATA_FILE, sheet = DATA_SHEET, guess_max = 1e6))
C <- sapply(c(exp = COL_EXPERIMENT, mat = COL_MATRIX, spe = COL_SPECIES, trt = COL_TREATMENT,
              tim = COL_TIME, rep = COL_REPLICATE, par = COL_PARENT, btp = COL_BTP),
            function(w) find_col(raw, w))
DAT <- data.frame(experiment = trimws(as.character(raw[[C["exp"]]])),
                  matrix     = trimws(as.character(raw[[C["mat"]]])),
                  species    = trimws(as.character(raw[[C["spe"]]])),
                  treatment  = trimws(as.character(raw[[C["trt"]]])),
                  time       = to_num(raw[[C["tim"]]]),
                  replicate  = as.character(raw[[C["rep"]]]),
                  parent_ug  = to_num(raw[[C["par"]]]),
                  btp_ug     = to_num(raw[[C["btp"]]]),
                  stringsAsFactors = FALSE)
DAT <- DAT[DAT$experiment == EXPERIMENT & DAT$treatment %in% TR & !is.na(DAT$time), ]
if (!nrow(DAT)) stop("No data for experiment '", EXPERIMENT, "' and treatments ",
                     paste(TR, collapse = ", "), ".")
DAT$parent_umol <- DAT$parent_ug / PARENT_MOLAR_MASS    # ug / (g/mol) = umol
DAT$btp_umol    <- DAT$btp_ug    / BTP_MOLAR_MASS

# water (exposure): mean of the replicates per time point
WATER <- DAT[DAT$matrix == EXPOSURE_MATRIX & !is.na(DAT$parent_ug), ]
EXPO  <- aggregate(cbind(parent_ug, parent_umol) ~ treatment + time, data = WATER, FUN = mean)
EXPO  <- EXPO[order(EXPO$treatment, EXPO$time), ]
# AE B107137 in the water (mean per time point; only for the figure)
WATER_BTP <- WATER[!is.na(WATER$btp_ug), ]
EXPO_BTP  <- if (nrow(WATER_BTP))
  aggregate(btp_ug ~ treatment + time, data = WATER_BTP, FUN = mean) else
    data.frame(treatment = character(0), time = numeric(0), btp_ug = numeric(0))
EXPO_BTP  <- EXPO_BTP[order(EXPO_BTP$treatment, EXPO_BTP$time), ]

# tissue: only SPECIES
TISSUE <- DAT[DAT$matrix == ORGANISM_MATRIX & DAT$species == SPECIES, ]
TISSUE <- TISSUE[order(TISSUE$treatment, TISSUE$time), ]
n_other <- sum(DAT$matrix == ORGANISM_MATRIX & DAT$species != SPECIES)

say("")
say("Data per treatment (water samples / tissue samples of ", SPECIES, "):")
T_END <- c()
for (tr in TR) {
  nw <- sum(WATER$treatment == tr); nt <- sum(TISSUE$treatment == tr)
  if (!nw) stop("No water data for treatment ", tr)
  if (!nt) stop("No tissue data of '", SPECIES, "' for treatment ", tr)
  T_END[tr] <- max(c(WATER$time[WATER$treatment == tr], TISSUE$time[TISSUE$treatment == tr]))
  say(sprintf("  %-4s %-16s %3d / %3d   (0 - %s d)", tr, TREATMENTS[[tr]], nw, nt, fmt(T_END[tr], 4)))
}
if (n_other) say("  (", n_other, " tissue samples with another Species entry are ignored)")

# single values excluded from the entire analysis (EXCLUDE_TISSUE)
if (nrow(EXCLUDE_TISSUE)) for (i in seq_len(nrow(EXCLUDE_TISSUE))) {
  x   <- EXCLUDE_TISSUE[i, ]
  hit <- TISSUE$treatment == x$treatment & abs(TISSUE$time - x$time) < 1e-3 &
    TISSUE$replicate == as.character(x$replicate)
  if (!any(hit)) { say("  ! EXCLUDE_TISSUE row ", i, " matches no tissue sample"); next }
  sts <- if (x$compound == "both") c("parent", "btp") else x$compound
  if (!all(sts %in% c("parent", "btp"))) stop("EXCLUDE_TISSUE: compound must be parent, btp or both.")
  for (st in sts) {
    say(sprintf("  Excluded: %s, t = %s d, replicate %s, %s = %s ug/kg ww", x$treatment, fmt(x$time, 4),
                x$replicate, STATE_NAMES[[st]], fmt(TISSUE[hit, paste0(st, "_ug")][1], 4)))
    TISSUE[hit, paste0(st, "_ug")] <- NA; TISSUE[hit, paste0(st, "_umol")] <- NA
  }
}

# plot AE B107137 in the water?
btp_frac <- if (nrow(EXPO_BTP)) {
  mm <- merge(EXPO_BTP, EXPO, by = c("treatment", "time"))
  mm$btp_ug / pmax(mm$parent_ug, 1e-12)
} else 0
SHOW_BTP_WATER <- if (identical(PLOT_BTP_WATER, "auto")) any(btp_frac >= BTP_WATER_MIN_FRACTION) else
  isTRUE(PLOT_BTP_WATER) && nrow(EXPO_BTP) > 0
say("  ", BTP_NAME, " in the water: ", fmt(100 * min(btp_frac), 2), " - ", fmt(100 * max(btp_frac), 3),
    "% of the parent -> ", if (SHOW_BTP_WATER) "plotted (dotted line)" else "not plotted")

## 3.2 Temperature loggers -----------------------------------------------------------------
lg <- as.data.frame(readxl::read_excel(TEMPERATURE_FILE, sheet = LOGGER_SHEET, guess_max = 1e6))
LOGGER <- list(); TEMP_DF <- list()
say("")
say("Water temperature during the exposure (logger):")
for (tr in TR) {
  cc <- LOGGER_COLUMNS[[tr]]
  if (max(cc) > ncol(lg)) stop("LOGGER_COLUMNS of ", tr, " point to a column that does not exist.")
  d <- data.frame(time = to_num(lg[[cc[["time"]]]]), T_C = to_num(lg[[cc[["temp"]]]]))
  d <- d[!is.na(d$time) & !is.na(d$T_C), ]
  d <- d[order(d$time), ]
  d <- d[!duplicated(d$time), ]
  if (nrow(d) < 2) stop("No temperature data for treatment ", tr)
  LOGGER[[tr]] <- d
  w <- d[d$time >= 0 & d$time <= T_END[tr], ]
  TEMP_DF[[tr]] <- data.frame(treatment = tr, w)
  if (min(d$time) > 0 || max(d$time) < T_END[tr])
    say(sprintf("  ! %s: logger covers %s - %s d, exposure 0 - %s d: temperature held constant outside",
                tr, fmt(min(d$time), 4), fmt(max(d$time), 4), fmt(T_END[tr], 4)))
  gaps <- diff(c(max(0, min(d$time)), w$time, min(T_END[tr], max(d$time)))) * 24
  if (any(gaps > MAX_LOGGER_GAP_H))
    say(sprintf("  ! %s: logger gap of up to %s h during the exposure (linear interpolation)",
                tr, fmt(max(gaps), 2)))
  say(sprintf("  %-4s mean %s degC (min %s, max %s), %d logger values",
              tr, fmt(mean(w$T_C), 3), fmt(min(w$T_C), 3), fmt(max(w$T_C), 3), nrow(w)))
}
TEMP_DF <- do.call(rbind, TEMP_DF)

## 3.3 Calibrated rate constants (reference) ------------------------------------------------
rc <- as.data.frame(readxl::read_excel(RATE_CONSTANT_FILE, sheet = 1))
rc <- rc[!is.na(rc$model_run), ]
need <- c("model_run", "T_C", as.vector(outer(c("ku", "k1", "ke", "km", "kem"), c("", "_lo", "_hi"), paste0)))
if (!all(need %in% names(rc)))
  stop("RATE_CONSTANT_FILE needs the columns: ", paste(setdiff(need, names(rc)), collapse = ", "))
ref <- rc[rc$model_run == REFERENCE_RUN, ]
if (nrow(ref) != 1) stop("REFERENCE_RUN '", REFERENCE_RUN, "' not found (once) in column model_run.")
REF    <- c(ku = ref$ku, ke = ref$ke, km = ref$km, kem = ref$kem)
REF    <- c(REF, k1 = REF[["ke"]] + REF[["km"]])
T_REF_C <- ref$T_C
T_REF_K <- T_REF_C + 273.15
say("")
say("Reference parameters (", REFERENCE_RUN, ", T_ref = ", fmt(T_REF_C, 4), " degC = ",
    fmt(T_REF_K, 5), " K):")
for (r in RATE_NAMES) say(sprintf("  %-4s = %s %s", r, fmt(REF[[r]], 5), RATE_UNITS[[r]]))
if (ref$ke <= ref$ke_lo * 1.01)
  say("  (ke is at the lower bound of its search range: elimination of the parent",
      " is dominated by biotransformation)")

## 3.4 Arrhenius temperatures ---------------------------------------------------------------
TA_RATES <- unique(unname(TA_FOR))
if (TA_SOURCE == "file") {
  sh <- readxl::excel_sheets(ARRHENIUS_FILE)
  ar <- as.data.frame(readxl::read_excel(ARRHENIUS_FILE,
                                         sheet = if ("arrhenius_results" %in% sh) "arrhenius_results" else 1))
  TA <- data.frame(rate = ar$rate, T_A = ar$T_A_K, SE = ar$SE_T_A)
} else if (TA_SOURCE == "manual") {
  TA <- TA_MANUAL
} else stop("TA_SOURCE must be \"file\" or \"manual\".")
TA <- TA[match(TA_RATES, TA$rate), ]
if (any(is.na(TA$T_A))) stop("T_A missing for: ", paste(TA_RATES[is.na(TA$T_A)], collapse = ", "))
rownames(TA) <- TA$rate
say("")
say("Arrhenius temperatures (", TA_SOURCE, "):")
for (r in RATE_NAMES)
  say(sprintf("  %-4s uses T_A of %-3s = %s +/- %s K (Ea = %s kJ/mol)", r, TA_FOR[[r]],
              fmt(TA[TA_FOR[[r]], "T_A"], 5), fmt(TA[TA_FOR[[r]], "SE"], 4),
              fmt(TA[TA_FOR[[r]], "T_A"] * 8.314 / 1000, 3)))


## 3.5 AE0+AE1 parameter sets of the 95% region (CI_REF_PARAMETERS) ------------------
PSAMPLE <- NULL
if (USE_PSAMPLE) {
  sh <- readxl::excel_sheets(PARAMETER_SAMPLE_FILE)
  if (!PARAMETER_SAMPLE_SHEET %in% sh)
    stop("Sheet '", PARAMETER_SAMPLE_SHEET, "' not found in ", basename(PARAMETER_SAMPLE_FILE), ".\n",
         "  Run TK_metabolism_model.R (version with the parameter_sample export) for the pooled\n",
         "  AE0+AE1 model with FIT_METHOD = \"parspace\", or set CI_REF_PARAMETERS <- FALSE.")
  PSAMPLE <- as.data.frame(readxl::read_excel(PARAMETER_SAMPLE_FILE, sheet = PARAMETER_SAMPLE_SHEET))
  if (!all(RATE_NAMES %in% names(PSAMPLE)) || nrow(PSAMPLE) < 2)
    stop("The parameter sample in ", basename(PARAMETER_SAMPLE_FILE), " has no sets or lacks ",
         paste(RATE_NAMES, collapse = ", "), ".")
  PSAMPLE <- PSAMPLE[stats::complete.cases(PSAMPLE[, RATE_NAMES]), ]
  say("")
  say("AE0+AE1 parameter sets of the 95% region: ", nrow(PSAMPLE), " sets (", basename(PARAMETER_SAMPLE_FILE), ")")
  for (r in RATE_NAMES)
    say(sprintf("  %-4s %s - %s %s  (set 1 = best fit: %s; %s: %s)", r, fmt(min(PSAMPLE[[r]]), 4),
                fmt(max(PSAMPLE[[r]]), 4), RATE_UNITS[[r]], fmt(PSAMPLE[[r]][1], 5),
                REFERENCE_RUN, fmt(REF[[r]], 5)))
  dev_best <- max(abs(unlist(PSAMPLE[1, RATE_NAMES]) / REF[RATE_NAMES] - 1))
  if (dev_best > 0.02)
    say("  ! The best set of the sample differs by up to ", fmt(100 * dev_best, 2), "% from ", REFERENCE_RUN,
        " in RATE_CONSTANT_FILE: are both from the same calibration run?")
}


###############################################################################
## SECTION 4 : Temperature dependence of the rate constants
###############################################################################
##  Every rate constant r is written as
##      k_r(T) = k_anchor * exp(-T_A * (1/T - x_anchor))
##  with the anchor point (x_anchor = 1/T_anchor, k_anchor):
##    reference : (1/T_ref, k_ref)
##    intercept : (weighted mean 1/T, weighted geometric mean of the
##                 calibrated estimates) = point of the weighted Arrhenius line
##                 that does not depend on the slope; x ratio for ke, km
##    none      : (1/T_ref, k_ref) and T_A = 0
##  Monte Carlo: T_A draws rotate the line around the anchor point; the
##  AE0+AE1 parameter sets ("reference", "none") or the SE of the line height
##  ("intercept") move k_anchor.

## 4.1 Anchor of the weighted Arrhenius line of one rate (version "intercept") ---------
runs <- if (identical(INTERCEPT_RUNS, "all")) rc$model_run else INTERCEPT_RUNS
if (!all(runs %in% rc$model_run))
  stop("INTERCEPT_RUNS not found: ", paste(setdiff(runs, rc$model_run), collapse = ", "))
rcI <- rc[rc$model_run %in% runs, ]
TREF_TXT  <- format(round(T_REF_C + 1e-9, 1), nsmall = 1)   # label, e.g. 20.25 -> "20.3"
#   L: rate whose calibrated estimates define the line; TA_L: its slope (T_A)
#   Returns the anchor (centroid) and the SE of ln k at the centroid:
#   SE = sigma_hat * sqrt(1/sum(w)), sigma_hat^2 = sum(w * resid^2) / df, with
#   df = n - 2 if the slope was fitted to these estimates (ku, k1), else n - 1
#   (slope taken from another rate, e.g. kem with T_A of k1). Same as lm().
line_anchor <- function(L, TA_L) {
  y  <- log(rcI[[L]])
  se <- (log(rcI[[paste0(L, "_hi")]]) - log(rcI[[paste0(L, "_lo")]])) / (2 * 1.96)
  w  <- 1 / se^2                                      # as in TK_arrhenius.R
  x  <- 1 / (rcI$T_C + 273.15)
  xc <- sum(w * x) / sum(w); yc <- sum(w * y) / sum(w)
  res <- y - (yc - TA_L * (x - xc))
  df  <- length(y) - if (L %in% TA_RATES) 2 else 1
  c(k = exp(yc), x = xc, se_ln = sqrt(sum(w * res^2) / df / sum(w)), n = length(y))
}

## 4.2 Anchor table per model version ---------------------------------------------------
make_anchors <- function(version) {
  out <- data.frame(rate = RATE_NAMES, TA_rate = unname(TA_FOR[RATE_NAMES]),
                    k_anchor = NA_real_, x_anchor = NA_real_, use_TA = 1, how = "",
                    line = "reference", SE_ln_anchor = 0, stringsAsFactors = FALSE)
  rownames(out) <- RATE_NAMES
  for (r in RATE_NAMES) {
    src <- INTERCEPT_FROM[[r]]
    if (version %in% c("reference", "none") || src == "reference") {
      out[r, c("k_anchor", "x_anchor")] <- c(REF[[r]], 1 / T_REF_K)
      out[r, "how"] <- paste0(REFERENCE_RUN, " value at T_ref = ", fmt(T_REF_C, 4), " degC")
    } else {
      if (!src %in% c("ku", "k1", "ke", "km", "kem")) stop("INTERCEPT_FROM: unknown value '", src, "'")
      a <- line_anchor(src, TA[TA_FOR[[r]], "T_A"])
      share <- REF[[r]] / REF[[src]]                  # 1 if the rate has its own line
      out[r, c("k_anchor", "x_anchor", "SE_ln_anchor")] <- c(a[["k"]] * share, a[["x"]], a[["se_ln"]])
      out[r, "line"] <- src
      out[r, "how"] <- if (src == r)
        paste0("weighted Arrhenius line through the ", r, " estimates (slope = T_A of ", TA_FOR[[r]], ")") else
          paste0("weighted Arrhenius line of ", src, " x ", REFERENCE_RUN, " ratio ", r, "/", src,
                 " = ", fmt(share, 4))
    }
    if (version == "none") { out[r, "use_TA"] <- 0; out[r, "how"] <- paste0(out[r, "how"], "; no T correction") }
  }
  out$T_anchor_C <- 1 / out$x_anchor - 273.15
  out$T_A        <- TA[out$TA_rate, "T_A"] * out$use_TA
  out$SE_T_A     <- TA[out$TA_rate, "SE"]  * out$use_TA
  out$lnA        <- log(out$k_anchor) + out$T_A * out$x_anchor     # k = A exp(-T_A/T)
  out
}
ANCHORS <- lapply(MODEL_VERSIONS, make_anchors)
names(ANCHORS) <- MODEL_VERSIONS

say("")
for (v in MODEL_VERSIONS) {
  say("Version \"", v, "\": ", VERSION_LABELS[[v]])
  a <- ANCHORS[[v]]
  for (r in RATE_NAMES)
    say(sprintf("  %-4s anchor %s %s at %s degC, lnA = %s%s  (%s)", r, fmt(a[r, "k_anchor"], 4),
                RATE_UNITS[[r]], fmt(a[r, "T_anchor_C"], 4), fmt(a[r, "lnA"], 5),
                if (a[r, "line"] != "reference") paste0(", SE(ln k) = ", fmt(a[r, "SE_ln_anchor"], 3)) else "",
                a[r, "how"]))
}

## 4.3 Rate constants for a vector of draws at one temperature ---------------------------
#   anc : anchor table of the version
#   D   : draws, list(TA = list(ku = .., k1 = ..), KA = list(ku = .., ke = .., km = .., kem = ..))
#         TA = T_A per T_A source, KA = k_anchor per rate (vectors of equal length)
#   invT: 1/T (scalar)
rates_at <- function(anc, D, invT) {
  k <- vector("list", 4); names(k) <- RATE_NAMES
  for (r in RATE_NAMES) {
    i <- which(anc$rate == r)
    k[[r]] <- D$KA[[r]] * exp(-D$TA[[anc$TA_rate[i]]] * anc$use_TA[i] * (invT - anc$x_anchor[i]))
  }
  k
}
subset_draws <- function(D, idx) list(TA = lapply(D$TA, `[`, idx), KA = lapply(D$KA, `[`, idx))

## 4.4 Best estimate and Monte Carlo draws per version -----------------------------------
best_draw <- function(v) {
  anc <- ANCHORS[[v]]
  list(TA = setNames(lapply(TA_RATES, function(r) TA[r, "T_A"]), TA_RATES),
       KA = setNames(lapply(RATE_NAMES, function(r) anc[r, "k_anchor"]), RATE_NAMES))
}
make_draws <- function(v) {
  anc <- ANCHORS[[v]]; D <- best_draw(v); src <- character(0); set_id <- NA
  D$TA <- lapply(D$TA, rep, N_MC); D$KA <- lapply(D$KA, rep, N_MC)
  if (CI_TA_SE && v != "none") {                         # slope: T_A ~ normal(mean, SE)
    D$TA <- setNames(lapply(TA_RATES, function(r) pmax(rnorm(N_MC, TA[r, "T_A"], TA[r, "SE"]), 0)), TA_RATES)
    src <- c(src, "T_A (SE)")
  }
  if (CI_REF_PARAMETERS && v %in% c("reference", "none")) {   # AE0+AE1 parameter sets
    set_id <- sample.int(nrow(PSAMPLE), N_MC, replace = TRUE)
    for (r in RATE_NAMES) if (anc[r, "line"] == "reference") D$KA[[r]] <- PSAMPLE[[r]][set_id]
    src <- c(src, paste0(REFERENCE_RUN, " parameter sets"))
  }
  if (CI_LINE_HEIGHT && v == "intercept") {                   # height of the Arrhenius lines
    for (L in setdiff(unique(anc$line), "reference")) {
      rr <- anc$rate[anc$line == L]
      dl <- rnorm(N_MC, 0, anc$SE_ln_anchor[anc$line == L][1])   # one shift per line
      for (r in rr) D$KA[[r]] <- D$KA[[r]] * exp(dl)
    }
    src <- c(src, "Arrhenius line height (SE)")
  }
  list(D = D, sources = src, set_id = set_id)
}
set.seed(RANDOM_SEED)
DRAWS <- lapply(MODEL_VERSIONS, make_draws); names(DRAWS) <- MODEL_VERSIONS
HAS_BAND <- sapply(DRAWS, function(d) length(d$sources) > 0)
say("")
say("Uncertainty in the CI bands (", N_MC, " Monte Carlo draws):")
for (v in MODEL_VERSIONS)
  say(sprintf("  %-10s %s", v, if (HAS_BAND[[v]]) paste(DRAWS[[v]]$sources, collapse = " + ") else "none (no CI band)"))


###############################################################################
## SECTION 5 : Model and solver
###############################################################################

## 5.1 Forcing functions per treatment: water concentration and temperature ----------
CW_FUN <- list(); TK_FUN <- list(); GRID <- list()
for (tr in TR) {
  e <- EXPO[EXPO$treatment == tr, ]
  CW_FUN[[tr]] <- approxfun(e$time, e$parent_umol, rule = 2)     # constant after last
  TK_FUN[[tr]] <- approxfun(LOGGER[[tr]]$time, LOGGER[[tr]]$T_C + 273.15, rule = 2)
  # integration grid: logger times, sampling times, output times, max step
  t_out  <- sort(unique(round(c(seq(0, T_END[tr], by = OUTPUT_DT), T_END[tr],
                                TISSUE$time[TISSUE$treatment == tr], e$time), 7)))
  lt     <- LOGGER[[tr]]$time
  t_grid <- sort(unique(round(c(t_out, lt[lt > 0 & lt < T_END[tr]]), 7)))
  big    <- which(diff(t_grid) > STEP_MAX_D)          # refine long steps
  if (length(big))
    t_grid <- sort(unique(round(c(t_grid, unlist(lapply(big, function(i)
      seq(t_grid[i], t_grid[i + 1], length.out = ceiling((t_grid[i + 1] - t_grid[i]) / STEP_MAX_D) + 1)))), 7)))
  GRID[[tr]] <- list(t = t_grid, out = t_out, rec = match(t_grid, t_out))
}

## 5.2 Runge-Kutta 4 for all draws at once ---------------------------------------------
#   Returns Cp and Cm [umol/kg ww] at the output times (rows) for every draw (columns).
simulate <- function(tr, anc, D) {
  g   <- GRID[[tr]]; tg <- g$t; h <- diff(tg); N <- length(D$KA[[1]])
  tm  <- head(tg, -1) + h / 2
  iT0 <- 1 / TK_FUN[[tr]](tg); iTm <- 1 / TK_FUN[[tr]](tm)
  cw0 <- CW_FUN[[tr]](tg);     cwm <- CW_FUN[[tr]](tm)
  deriv <- function(Cp, Cm, k, cw)
    list(p = k$ku * cw - (k$ke + k$km) * Cp,
         m = k$km * Cp - k$kem * Cm)
  Cp <- rep(X0_PARENT / PARENT_MOLAR_MASS, N)
  Cm <- rep(X0_BTP / BTP_MOLAR_MASS, N)
  P  <- matrix(NA_real_, length(g$out), N); M <- P
  if (!is.na(g$rec[1])) { P[g$rec[1], ] <- Cp; M[g$rec[1], ] <- Cm }
  k0 <- rates_at(anc, D, iT0[1])
  for (i in seq_along(h)) {
    kh <- rates_at(anc, D, iTm[i])
    k1 <- rates_at(anc, D, iT0[i + 1])
    s  <- h[i]
    d1 <- deriv(Cp, Cm, k0, cw0[i])
    d2 <- deriv(Cp + s / 2 * d1$p, Cm + s / 2 * d1$m, kh, cwm[i])
    d3 <- deriv(Cp + s / 2 * d2$p, Cm + s / 2 * d2$m, kh, cwm[i])
    d4 <- deriv(Cp + s * d3$p,     Cm + s * d3$m,     k1, cw0[i + 1])
    Cp <- Cp + s / 6 * (d1$p + 2 * d2$p + 2 * d3$p + d4$p)
    Cm <- Cm + s / 6 * (d1$m + 2 * d2$m + 2 * d3$m + d4$m)
    j  <- g$rec[i + 1]
    if (!is.na(j)) { P[j, ] <- Cp; M[j, ] <- Cm }
    k0 <- k1
  }
  list(t = g$out, Cp = P, Cm = M)
}

## 5.3 Check with deSolve (best estimate) ----------------------------------------------
check_desolve <- function(tr, anc, best, Db) {
  fn <- function(t, y, parms) {
    k <- rates_at(anc, Db, 1 / TK_FUN[[tr]](t))
    list(c(k$ku * CW_FUN[[tr]](t) - (k$ke + k$km) * y[1],
           k$km * y[1] - k$kem * y[2]))
  }
  y0  <- c(X0_PARENT / PARENT_MOLAR_MASS, X0_BTP / BTP_MOLAR_MASS)
  sol <- deSolve::lsoda(y0, GRID[[tr]]$out, fn, NULL, rtol = 1e-9, atol = 1e-12,
                        hmax = STEP_MAX_D)
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  c(parent = rel(best$Cp[, 1], sol[, 2]), btp = rel(best$Cm[, 1], sol[, 3]))
}


###############################################################################
## SECTION 6 : Run all versions and treatments (with progress bar)
###############################################################################
say("")
say("Simulations: ", length(MODEL_VERSIONS), " version(s) x ", length(TR), " treatments, ",
    N_MC, " Monte Carlo draws each ...")
alpha  <- (1 - CI_LEVEL) / 2
chunks <- split(seq_len(N_MC), ceiling(seq_len(N_MC) / MC_CHUNK))
mc_versions <- MODEL_VERSIONS[HAS_BAND[MODEL_VERSIONS]]
n_ticks <- max(1, length(mc_versions) * length(TR) * length(chunks))
pb <- txtProgressBar(min = 0, max = n_ticks, style = 3); tick <- 0
t_start <- Sys.time()

MODEL <- list(); CHECK <- list()
for (v in MODEL_VERSIONS) {
  anc <- ANCHORS[[v]]; Db <- best_draw(v)
  for (tr in TR) {
    best <- simulate(tr, anc, Db)
    if (v %in% mc_versions) {
      Pmc <- NULL; Mmc <- NULL
      for (ch in chunks) {
        s   <- simulate(tr, anc, subset_draws(DRAWS[[v]]$D, ch))
        Pmc <- cbind(Pmc, s$Cp); Mmc <- cbind(Mmc, s$Cm)
        tick <- tick + 1; setTxtProgressBar(pb, tick)
      }
      qP <- apply(Pmc, 1, quantile, probs = c(alpha, 1 - alpha), names = FALSE)
      qM <- apply(Mmc, 1, quantile, probs = c(alpha, 1 - alpha), names = FALSE)
    } else {                                           # no uncertainty source: no band
      qP <- rbind(best$Cp[, 1], best$Cp[, 1]); qM <- rbind(best$Cm[, 1], best$Cm[, 1])
    }
    t  <- best$t
    kb <- rates_at(anc, Db, 1 / TK_FUN[[tr]](t))
    MODEL[[paste(v, tr)]] <- data.frame(
      version = v, treatment = tr, time_d = t,
      T_water_C = TK_FUN[[tr]](t) - 273.15,
      Cw_parent_ug_L = CW_FUN[[tr]](t) * PARENT_MOLAR_MASS,
      parent_ug_kg = best$Cp[, 1] * PARENT_MOLAR_MASS,
      parent_ug_kg_lo = qP[1, ] * PARENT_MOLAR_MASS, parent_ug_kg_hi = qP[2, ] * PARENT_MOLAR_MASS,
      btp_ug_kg = best$Cm[, 1] * BTP_MOLAR_MASS,
      btp_ug_kg_lo = qM[1, ] * BTP_MOLAR_MASS, btp_ug_kg_hi = qM[2, ] * BTP_MOLAR_MASS,
      Cw_parent_umol_L = CW_FUN[[tr]](t),
      parent_umol_kg = best$Cp[, 1], parent_umol_kg_lo = qP[1, ], parent_umol_kg_hi = qP[2, ],
      btp_umol_kg = best$Cm[, 1], btp_umol_kg_lo = qM[1, ], btp_umol_kg_hi = qM[2, ],
      ku = kb$ku, ke = kb$ke, km = kb$km, kem = kb$kem,
      stringsAsFactors = FALSE)
    if (CHECK_WITH_DESOLVE && HAVE_DESOLVE) CHECK[[paste(v, tr)]] <- check_desolve(tr, anc, best, Db)
  }
}
setTxtProgressBar(pb, n_ticks); close(pb)
MODEL <- do.call(rbind, MODEL); rownames(MODEL) <- NULL
say("Done in ", format(round(difftime(Sys.time(), t_start, units = "secs"), 1)), ".")
if (length(CHECK)) {
  dev <- max(unlist(CHECK))
  say("Check with deSolve::lsoda: max. relative deviation of the best estimate = ",
      formatC(dev, format = "e", digits = 1), if (dev < 1e-4) " (ok)" else "  ! check STEP_MAX_D")
} else if (CHECK_WITH_DESOLVE) say("(deSolve not installed: solver check skipped)")


###############################################################################
## SECTION 7 : Quality of fit and validity criteria
###############################################################################

## 7.1 Observed vs predicted (all tissue samples) ---------------------------------------
OVP <- list()
for (v in MODEL_VERSIONS) for (tr in TR) {
  m <- MODEL[MODEL$version == v & MODEL$treatment == tr, ]
  o <- TISSUE[TISSUE$treatment == tr, ]
  i <- match(round(o$time, 7), round(m$time_d, 7))
  for (st in c("parent", "btp")) {
    OVP[[paste(v, tr, st)]] <- data.frame(
      version = v, treatment = tr, compound = STATE_NAMES[[st]], state = st,
      time_d = o$time, replicate = o$replicate,
      observed_ug_kg  = o[[paste0(st, "_ug")]],
      predicted_ug_kg = m[[paste0(st, "_ug_kg")]][i],
      CI_lo_ug_kg     = m[[paste0(st, "_ug_kg_lo")]][i],
      CI_hi_ug_kg     = m[[paste0(st, "_ug_kg_hi")]][i],
      observed_umol_kg  = o[[paste0(st, "_umol")]],
      predicted_umol_kg = m[[paste0(st, "_umol_kg")]][i],
      stringsAsFactors = FALSE)
  }
}
OVP <- do.call(rbind, OVP); rownames(OVP) <- NULL
OVP$inside_CI <- OVP$observed_ug_kg >= OVP$CI_lo_ug_kg - 1e-12 &
  OVP$observed_ug_kg <= OVP$CI_hi_ug_kg + 1e-12

## 7.2 Criteria for one subset ---------------------------------------------------------
metrics <- function(d, v) {
  o <- d$observed_umol_kg; p <- pmax(d$predicted_umol_kg, 0); n <- length(o)
  ssq <- sum((tr_data(o) - tr_data(p))^2)
  sst <- sum((tr_data(o) - mean(tr_data(o)))^2)
  mll <- if (ssq > 0) n / 2 * log(ssq) + n / 2 else NA
  # PE at the last sampling time of each treatment in the subset
  sp <- sapply(split(d, d$treatment), function(s) {
    tl <- max(s$time_d); sel <- abs(s$time_d - tl) < 1e-9
    mo <- mean(s$observed_umol_kg[sel])
    c(sppe = if (mo > 0) 100 * (mean(s$predicted_umol_kg[sel]) - mo) / mo else NA, t = tl)
  })
  ppc <- if (!HAS_BAND[[v]]) NA else 100 * mean(d$inside_CI)
  nrmse <- if (mean(o) > 0) 100 * sqrt(mean((p - o)^2)) / mean(o) else NA
  data.frame(n = n, R2 = 1 - ssq / sst, MLL = mll, AIC = 2 * mll,
             NRMSE_pct = nrmse, PPC_pct = ppc,
             PE_pct = mean(sp["sppe", ]),
             PE_time_d = if (ncol(sp) == 1) sp["t", 1] else NA,
             NRMSE_ok = !is.na(nrmse) & nrmse <= NRMSE_MAX,
             PPC_ok = if (is.na(ppc)) NA else ppc >= PPC_MIN,
             PE_ok = !is.na(mean(sp["sppe", ])) & abs(mean(sp["sppe", ])) <= PE_MAX_ABS)
}

## 7.3 Table per version x compound x treatment (+ pooled) --------------------------------
QOF <- list()
for (v in MODEL_VERSIONS) {
  for (st in c("parent", "btp")) {
    d <- OVP[OVP$version == v & OVP$state == st & !is.na(OVP$observed_umol_kg), ]
    for (sc in c(TR, "all")) {
      s <- if (sc == "all") d else d[d$treatment == sc, ]
      if (!nrow(s)) next
      QOF[[length(QOF) + 1]] <- data.frame(
        version = v, compound = STATE_NAMES[[st]],
        treatment = if (sc == "all") "all (pooled)" else sc, metrics(s, v),
        stringsAsFactors = FALSE)
    }
  }
}
QOF <- do.call(rbind, QOF); rownames(QOF) <- NULL
QOF$all_criteria_met <- QOF$NRMSE_ok & QOF$PE_ok & (is.na(QOF$PPC_ok) | QOF$PPC_ok)

## 7.4 Version comparison (both compounds, pooled) -------------------------------------
VCOMP <- do.call(rbind, lapply(MODEL_VERSIONS, function(v) {
  q <- QOF[QOF$version == v & QOF$treatment == "all (pooled)", ]
  data.frame(version = v, description = VERSION_LABELS[[v]],
             n_total = sum(q$n), MLL_total = sum(q$MLL), AIC_total = sum(q$AIC),
             AIC_parent = q$AIC[q$compound == PARENT_NAME], AIC_btp = q$AIC[q$compound == BTP_NAME],
             R2_parent = q$R2[q$compound == PARENT_NAME], R2_btp = q$R2[q$compound == BTP_NAME],
             NRMSE_parent = q$NRMSE_pct[q$compound == PARENT_NAME],
             NRMSE_btp = q$NRMSE_pct[q$compound == BTP_NAME],
             stringsAsFactors = FALSE)
}))
VCOMP$delta_AIC <- VCOMP$AIC_total - min(VCOMP$AIC_total)

## 7.5 Console summary ---------------------------------------------------------------------
say("")
say(sprintf("Quality of fit (PPC ok >= %g%%, NRMSE ok <= %g%%, PE ok within +/-%g%%; PE > 0 = overprediction)",
            PPC_MIN, NRMSE_MAX, PE_MAX_ABS))
say(sprintf("  %-10s %-13s %-13s %4s %7s %8s %7s %7s %8s", "version", "compound", "treatment",
            "n", "R2", "NRMSE%", "PPC%", "PE%", "AIC"))
for (i in seq_len(nrow(QOF))) with(QOF[i, ], say(sprintf(
  "  %-10s %-13s %-13s %4d %7s %7s%s %6s%s %6s%s %8s", version, compound, treatment, n,
  fmt(R2, 3), fmt(NRMSE_pct, 3), if (isTRUE(NRMSE_ok)) " " else "!",
  fmt(PPC_pct, 3), if (isFALSE(PPC_ok)) "!" else " ", fmt(PE_pct, 3),
  if (isTRUE(PE_ok)) " " else "!", fmt(AIC, 4))))
if (length(MODEL_VERSIONS) > 1) {
  say("")
  say("Version comparison (AIC of both compounds, pooled; lower = better):")
  for (i in seq_len(nrow(VCOMP))) say(sprintf("  %-10s AIC = %s   dAIC = %s", VCOMP$version[i],
                                              fmt(VCOMP$AIC_total[i], 5), fmt(VCOMP$delta_AIC[i], 3)))
}


###############################################################################
## SECTION 8 : Figure (one 4-panel figure per version)
###############################################################################

## 8.1 Theme (layout as the Arrhenius figure) --------------------------------------------
theme_tk <- function() {
  theme_classic(base_size = BASE_SIZE, base_family = "sans") +
    theme(panel.border       = element_rect(fill = NA, colour = "black", linewidth = 0.45),
          axis.line          = element_blank(),
          panel.grid.major.y = element_line(colour = GRID_COL, linewidth = 0.3),
          axis.ticks         = element_line(colour = "black", linewidth = 0.35),
          axis.ticks.length  = unit(1, "mm"),
          axis.text          = element_text(colour = "black", size = BASE_SIZE + 1),
          axis.text.x        = element_text(margin = margin(t = 2.5)),
          axis.text.y        = element_text(margin = margin(r = 2.5)),
          axis.title         = element_text(colour = "black", size = BASE_SIZE + 1.5),
          plot.title         = element_text(face = "bold", size = BASE_SIZE + 2, hjust = 0,
                                            margin = margin(b = 2)),
          plot.subtitle      = element_text(colour = "grey20", size = BASE_SIZE, hjust = 0,
                                            margin = margin(b = 5)),
          plot.title.position = "panel",
          plot.tag           = element_text(face = "bold", size = BASE_SIZE + 3),
          plot.margin        = margin(4, 8, 4, 4),
          legend.title       = element_text(colour = "grey35", size = BASE_SIZE + 0.5),
          legend.margin      = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "mm"),
          legend.text        = element_text(size = BASE_SIZE + 1.5),
          legend.key         = element_blank(),
          legend.key.height  = unit(5, "mm"),
          legend.key.width   = unit(7, "mm"),
          legend.spacing.y   = unit(1, "mm"),
          strip.background   = element_blank(),
          strip.text         = element_text(face = "bold", size = BASE_SIZE + 1),
          strip.placement    = "outside",
          panel.spacing      = unit(3, "mm"))
}
LAB  <- TREATMENTS[TR]
COLS <- TREATMENT_COLOURS[TR]
sc_col  <- scale_colour_manual("Treatment", values = COLS, labels = LAB, breaks = TR, limits = TR)
sc_fill <- scale_fill_manual("Treatment",   values = COLS, labels = LAB, breaks = TR, limits = TR)

## Axes: fixed ranges from AXIS_LIMITS / X_LIMITS; values outside are kept (oob_keep)
## and only cut off at the panel edge. Time ticks: every 2 d (minor 1 d) for
## ranges > 6 d, else every 1 d (minor 0.5 d).
brk_x   <- function(l) seq(0, l[2], by = if (diff(l) > 6) 2 else 1)
minor_x <- function(l) seq(0, l[2], by = if (diff(l) > 6) 1 else 0.5)
sc_x <- scale_x_continuous("Time [d]", limits = c(0, max(X_LIMITS[TR])), breaks = brk_x,
                           minor_breaks = minor_x, oob = scales::oob_keep,
                           expand = expansion(mult = c(0.02, 0.02)))
sc_x_free <- scale_x_continuous("Time [d]", breaks = brk_x, minor_breaks = minor_x,
                                expand = expansion(mult = c(0.02, 0.03)))
x_blank <- function()             # sets 0 - X_LIMITS per treatment (free time axes)
  geom_blank(data = data.frame(treatment = factor(rep(TR, 2), levels = TR),
                               time = c(rep(0, length(TR)), unname(X_LIMITS[TR]))),
             aes(x = time), inherit.aes = FALSE)
or_else <- function(a, b) if (is.null(a)) b else a
sc_y <- function(q, lab, top = 0.06, v = NULL) {
  lim <- or_else(if (!is.null(v)) AXIS_LIMITS_VERSION[[v]][[q]], AXIS_LIMITS[[q]])
  brk <- or_else(if (!is.null(v)) AXIS_BREAKS_VERSION[[v]][[q]], or_else(AXIS_BREAKS[[q]], waiver()))
  if (is.null(lim)) {
    scale_y_continuous(lab, limits = if (q == "temp") NULL else c(0, NA), breaks = brk,
                       expand = expansion(mult = c(if (q == "temp") 0.06 else 0, top)))
  } else {
    scale_y_continuous(lab, limits = lim, breaks = brk, oob = scales::oob_keep,
                       expand = expansion(mult = c(0.02, 0.02)))
  }
}
minor_ticks <- function()
  if (MINOR_TICKS) list(guides(x = guide_axis(minor.ticks = TRUE), y = guide_axis(minor.ticks = TRUE)),
                        theme(axis.minor.ticks.length = rel(0.5))) else NULL
LAB_W  <- expression(italic(C)[w] ~ "[\u00b5g" ~ L^-1 * "]")
LAB_T  <- "Water temperature [\u00b0C]"
LAB_CT <- expression(italic(C)[tissue] ~ "[\u00b5g" ~ kg^-1 ~ "ww]")

qof_txt <- function(v, st) {
  q <- QOF[QOF$version == v & QOF$compound == STATE_NAMES[[st]] & QOF$treatment == "all (pooled)", ]
  if (!nrow(q)) return(NULL)
  paste0("R\u00b2 = ", sprintf("%.2f", q$R2), ", NRMSE = ", sprintf("%.0f", q$NRMSE_pct), "%",
         if (!is.na(q$PPC_pct)) paste0(", PPC = ", sprintf("%.0f", q$PPC_pct), "%") else "",
         "\nPE = ", sprintf("%+.0f", q$PE_pct), "%  (pooled)")
}

## Water data in long format: parent (model input, solid, circles) and, if
## SHOW_BTP_WATER, AE B107137 (measured only, dotted, triangles)
CMP_LEVELS <- c(PARENT_NAME, BTP_NAME)
water_lines <- function(v) {
  m <- MODEL[MODEL$version == v, ]
  d <- data.frame(treatment = m$treatment, time = m$time_d, y = m$Cw_parent_ug_L, compound = PARENT_NAME)
  if (SHOW_BTP_WATER)
    d <- rbind(d, data.frame(treatment = EXPO_BTP$treatment, time = EXPO_BTP$time, y = EXPO_BTP$btp_ug,
                             compound = BTP_NAME))
  d$compound <- factor(d$compound, levels = CMP_LEVELS); d
}
water_points <- function() {
  d <- data.frame(treatment = WATER$treatment, time = WATER$time, y = WATER$parent_ug, compound = PARENT_NAME)
  if (SHOW_BTP_WATER)
    d <- rbind(d, data.frame(treatment = WATER_BTP$treatment, time = WATER_BTP$time, y = WATER_BTP$btp_ug,
                             compound = BTP_NAME))
  d$compound <- factor(d$compound, levels = CMP_LEVELS); d
}
sc_lt <- scale_linetype_manual("Water", values = setNames(c("solid", "11"), CMP_LEVELS), drop = TRUE)
sc_sh <- scale_shape_manual("Water", values = setNames(c(21, 24), CMP_LEVELS), drop = TRUE)
water_layers <- function(v, pt_size = 1.9) list(
  geom_line(data = water_lines(v), aes(time, y, colour = treatment, linetype = compound), linewidth = 0.7),
  geom_point(data = water_points(), aes(time, y, fill = treatment, shape = compound),
             colour = "white", size = pt_size, stroke = 0.35, alpha = POINT_ALPHA),
  sc_lt, sc_sh)
WATER_NOTE <- paste0(PARENT_NAME, ": solid line (model input) and circles",
                     if (SHOW_BTP_WATER) paste0("; ", BTP_NAME, ": dotted line and triangles (measured, not used by the model)") else "")

## 8.2 Panels of the overview layout -------------------------------------------------------
panel_water <- function(v) {
  ggplot() + water_layers(v) +
    sc_col + sc_fill + sc_x + sc_y("water", LAB_W, v = v) +
    labs(title = "Water concentration", subtitle = "Points: measured; lines: linear interpolation\n") +
    theme_tk() + minor_ticks() +
    guides(colour = "none", fill = "none",
           linetype = guide_legend(order = 2), shape = guide_legend(order = 2,
                                                                    override.aes = list(fill = "grey45", colour = "grey30", size = 2)))
}
panel_temp <- function(v) {
  p <- ggplot(TEMP_DF, aes(time, T_C, colour = treatment)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    sc_col + sc_x + sc_y("temp", LAB_T, v = v) +
    labs(title = "Water temperature",
         subtitle = if (v == "reference" && SHOW_TREF_LINE)
           paste0("Logger data (5 min)\ndashed: T_ref = ", TREF_TXT, " \u00b0C") else
             "Logger data (5 min)\n") +
    theme_tk() + minor_ticks() + theme(legend.position = "none")
  if (v == "reference" && SHOW_TREF_LINE)
    p <- p + geom_hline(yintercept = T_REF_C, linetype = "22", colour = "grey30", linewidth = 0.4)
  p
}
panel_tissue <- function(v, st, legend = FALSE) {
  m <- MODEL[MODEL$version == v, ]
  m$y <- m[[paste0(st, "_ug_kg")]]; m$lo <- m[[paste0(st, "_ug_kg_lo")]]; m$hi <- m[[paste0(st, "_ug_kg_hi")]]
  o <- TISSUE; o$y <- o[[paste0(st, "_ug")]]; o <- o[!is.na(o$y), ]
  p <- ggplot() +
    geom_ribbon(data = m, aes(time_d, ymin = lo, ymax = hi, fill = treatment), alpha = BAND_ALPHA) +
    geom_line(data = m, aes(time_d, y, colour = treatment), linewidth = 0.7) +
    geom_point(data = o, aes(time, y, fill = treatment), shape = 21,
               colour = "white", size = 1.9, stroke = 0.35, alpha = POINT_ALPHA) +
    sc_col + sc_fill + sc_x + sc_y(st, LAB_CT, v = v) +
    labs(title = paste(STATE_NAMES[[st]], "in tissue"), subtitle = qof_txt(v, st)) +
    theme_tk() + minor_ticks() +
    guides(colour = guide_legend(order = 1), fill = guide_legend(order = 1))
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

## 8.3 Panels of the treatment layouts: one facet per treatment -------------------------
#   dir = "rows"    : facets stacked (rows = treatments), strip on the right
#   dir = "columns" : facets side by side (columns = treatments), strip on top
#   Every treatment has its own time axis (0 - X_LIMITS); y axes as AXIS_LIMITS.
fac <- function(d) { d$treatment <- factor(d$treatment, levels = TR); d }
USE_GGTEXT <- COLOURED_STRIPS && requireNamespace("ggtext", quietly = TRUE)
if (COLOURED_STRIPS && !USE_GGTEXT)
  say("(Package 'ggtext' not installed: treatment names in black. install.packages(\"ggtext\"))")
STRIP_LAB <- if (USE_GGTEXT)
  setNames(sprintf("<span style='color:%s'>%s</span>", TREATMENT_COLOURS[TR], TREATMENTS[TR]), TR) else
    TREATMENTS[TR]
facet_tr <- function(dir)
  facet_wrap(~ treatment, scales = "free_x", labeller = as_labeller(STRIP_LAB),
             ncol = if (dir == "rows") 1 else length(TR),
             strip.position = if (dir == "rows") "right" else "top")
theme_alt <- function()                               # compact rows, coloured strips
  theme(legend.position = "none",
        plot.title  = element_text(face = "bold", size = BASE_SIZE + 2, hjust = 0, margin = margin(b = 2)),
        plot.margin = margin(4, 6, 1, 3),
        strip.text  = if (USE_GGTEXT) ggtext::element_markdown(face = "bold", size = BASE_SIZE + 1) else
          element_text(face = "bold", size = BASE_SIZE + 1))
first_tr <- function(d) { d$treatment <- factor(rep(TR[1], nrow(d)), levels = TR); d }  # first panel only
TXT_COL  <- "grey40"
TXT_SIZE <- (BASE_SIZE - 0.5) / .pt
ok_sym   <- function(ok) ifelse(is.na(ok), "", ifelse(ok, " \u2713", " \u2717"))

# per-treatment criteria with tick (met) or cross (not met), all in grey
qc_labels <- function(v, st) {
  q <- QOF[QOF$version == v & QOF$compound == STATE_NAMES[[st]] & QOF$treatment %in% TR, ]
  data.frame(treatment = factor(q$treatment, levels = TR),
             label = paste0("NRMSE ", sprintf("%.0f", q$NRMSE_pct), "%", ok_sym(q$NRMSE_ok),
                            ifelse(is.na(q$PPC_pct), "",
                                   paste0("\nPPC ", sprintf("%.0f", q$PPC_pct), "%", ok_sym(q$PPC_ok))),
                            "\nPE ", sprintf("%+.0f", q$PE_pct), "%", ok_sym(q$PE_ok)))
}

alt_temp <- function(v, dir, letter) {
  p <- ggplot(fac(TEMP_DF), aes(time, T_C, colour = treatment)) +
    geom_line(linewidth = 0.6, alpha = 0.9) + x_blank() +
    facet_tr(dir) + sc_col + sc_x_free + sc_y("temp", LAB_T, v = v) +
    labs(title = paste0(letter, "   Water temperature")) +
    theme_tk() + minor_ticks() + theme_alt()
  if (v == "reference" && SHOW_TREF_LINE)               # reference temperature (dashed)
    p <- p + geom_hline(yintercept = T_REF_C, linetype = "22", colour = TXT_COL, linewidth = 0.4)
  p
}
alt_water <- function(v, dir, letter) {
  L <- water_layers(v, pt_size = 1.7)
  L[[1]]$data <- fac(L[[1]]$data); L[[2]]$data <- fac(L[[2]]$data)
  p <- ggplot() + L + x_blank() +
    facet_tr(dir) + sc_col + sc_fill + sc_x_free + sc_y("water", LAB_W, v = v) +
    labs(title = paste0(letter, "   Water concentration")) +
    theme_tk() + minor_ticks() + theme_alt()
  # key for the water lines, top right of the first panel
  xr  <- unname(X_LIMITS[TR[1]])
  top <- or_else(AXIS_LIMITS$water[2], max(WATER$parent_ug, na.rm = TRUE))
  key <- data.frame(compound = factor(CMP_LEVELS, levels = CMP_LEVELS), y = top * c(0.95, 0.84))
  if (!SHOW_BTP_WATER) key <- key[1, , drop = FALSE]
  key <- first_tr(key)
  p + geom_segment(data = key, aes(x = 0.50 * xr, xend = 0.62 * xr, y = y, yend = y, linetype = compound),
                   colour = "grey30", linewidth = 0.6, inherit.aes = FALSE) +
    geom_point(data = key, aes(x = 0.56 * xr, y = y, shape = compound), fill = "grey55",
               colour = "white", size = 1.7, stroke = 0.35, inherit.aes = FALSE) +
    geom_text(data = key, aes(x = 0.65 * xr, y = y, label = compound), hjust = 0, vjust = 0.5,
              size = TXT_SIZE, colour = "grey25", inherit.aes = FALSE)
}
alt_tissue <- function(v, st, dir, letter) {
  m <- fac(MODEL[MODEL$version == v, ])
  m$y <- m[[paste0(st, "_ug_kg")]]; m$lo <- m[[paste0(st, "_ug_kg_lo")]]; m$hi <- m[[paste0(st, "_ug_kg_hi")]]
  o <- fac(TISSUE); o$y <- o[[paste0(st, "_ug")]]; o <- o[!is.na(o$y), ]
  p <- ggplot() +
    geom_ribbon(data = m, aes(time_d, ymin = lo, ymax = hi, fill = treatment), alpha = BAND_ALPHA) +
    geom_line(data = m, aes(time_d, y, colour = treatment), linewidth = 0.7) +
    geom_point(data = o, aes(time, y, fill = treatment), shape = 21,
               colour = "white", size = 1.7, stroke = 0.3, alpha = POINT_ALPHA) + x_blank() +
    facet_tr(dir) + sc_col + sc_fill + sc_x_free +
    sc_y(st, LAB_CT, top = if (SHOW_QC_IN_PANELS) 0.40 else 0.06, v = v) +
    labs(title = paste0(letter, "   ", STATE_NAMES[[st]], " in tissue")) +
    theme_tk() + minor_ticks() + theme_alt()
  if (SHOW_QC_IN_PANELS)
    p <- p + geom_text(data = qc_labels(v, st), aes(x = Inf, y = Inf, label = label),
                       hjust = 1.06, vjust = 1.12, size = TXT_SIZE, colour = TXT_COL, lineheight = 0.9)
  p
}

## 8.4 Assemble and save --------------------------------------------------------------------
# PNG with ragg if installed (renders the tick/cross symbols reliably), else cairo
save_fig <- function(fig, file, layout) {
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else
    function(filename, width, height, ...)
      grDevices::png(filename, width = width, height = height, units = "in", res = PNG_DPI, ..., type = "cairo")
  ggsave(paste0(file, ".png"), fig, device = dev, width = FIG_SIZE_MM[[layout]][["width"]],
         height = FIG_SIZE_MM[[layout]][["height"]], units = "mm", dpi = PNG_DPI, bg = "white")
}
annotate_fig <- function(fig, v, extra = "", tags = TRUE) {
  band_txt <- if (!HAS_BAND[[v]]) "no CI band" else
    sprintf("bands: %g%% CI from %s (Monte Carlo, n = %d)", 100 * CI_LEVEL,
            paste(DRAWS[[v]]$sources, collapse = " + "), N_MC)
  fig + plot_annotation(
    tag_levels = if (tags) "a" else NULL,
    title   = paste0("Model validation ", EXPERIMENT, ": ", VERSION_LABELS[[v]]),
    caption = paste(strwrap(paste0(
      "Lines: model prediction (best estimate); ", band_txt, ". ",
      "Points: measured (", SPECIES, ").", extra), CAPTION_WRAP), collapse = "\n"),
    theme = theme(plot.title   = element_text(face = "bold", size = BASE_SIZE + 2, margin = margin(b = 2)),
                  plot.caption = element_text(colour = "grey35", size = BASE_SIZE - 0.5, hjust = 0)))
}
build_figure <- function(v, layout) {
  if (layout == "overview")
    return(annotate_fig((panel_water(v) + panel_temp(v)) /
                          (panel_tissue(v, "parent", legend = TRUE) + panel_tissue(v, "btp")) +
                          plot_layout(guides = "collect"), v))
  dir <- if (layout == "treatment_rows") "rows" else "columns"
  # order: a temperature, b water concentration, c parent, d BTP (letter in the title)
  p <- list(alt_temp(v, dir, "a"), alt_water(v, dir, "b"),
            alt_tissue(v, "parent", dir, "c"), alt_tissue(v, "btp", dir, "d"))
  if (dir == "rows") {                       # treatment names only right of the last column
    for (i in 1:3) p[[i]] <- p[[i]] + theme(strip.text = element_blank())
    fig <- wrap_plots(p, nrow = 1)
  } else {                                   # treatment names only above the first row
    for (i in 2:4) p[[i]] <- p[[i]] + theme(strip.text = element_blank())
    for (i in 1:3) p[[i]] <- p[[i]] + theme(axis.title.x = element_blank())
    fig <- wrap_plots(p, ncol = 1)
  }
  annotate_fig(fig, v, tags = FALSE, extra = paste0(
    if (v == "reference" && SHOW_TREF_LINE)
      paste0(" Dashed line in a: reference temperature of the calibration (", TREF_TXT, " \u00b0C).") else "",
    if (SHOW_QC_IN_PANELS) sprintf(paste0(
      " Numbers: NRMSE, PPC and PE (prediction error at the last sampling time) of each treatment;",
      " \u2713 / \u2717 = criterion met / not met (NRMSE \u2264 %g%%, PPC \u2265 %g%%, |PE| \u2264 %g%%)."),
      NRMSE_MAX, PPC_MIN, PE_MAX_ABS) else ""))
}
FIG_SUFFIX <- c(overview = "", treatment_rows = "_treatment_rows", treatment_columns = "_treatment_columns")
FIGS <- list(); fig_files <- character(0)
for (v in MODEL_VERSIONS) for (lay in FIGURE_LAYOUTS) {
  fig <- build_figure(v, lay)
  FIGS[[paste(v, lay)]] <- fig
  save_fig(fig, out_file(v, "_figure", FIG_SUFFIX[[lay]]), lay)
  fig_files <- c(fig_files, paste0(RUN_NAME, "_", v, "_figure", FIG_SUFFIX[[lay]], ".png"))
  if (lay == FIGURE_LAYOUTS[1]) print(fig)
}
say("")
say("Figures saved: ", paste(fig_files, collapse = ", "))


###############################################################################
## SECTION 9 : Excel output
###############################################################################

## 9.1 Rate constants per version and treatment ----------------------------------------
RATES_TAB <- do.call(rbind, lapply(split(MODEL, list(MODEL$version, MODEL$treatment), drop = TRUE), function(m) {
  t <- m$time_d
  row <- data.frame(version = m$version[1], treatment = m$treatment[1],
                    T_mean_C = tw_mean(t, m$T_water_C), T_min_C = min(m$T_water_C), T_max_C = max(m$T_water_C))
  for (r in RATE_NAMES) {
    row[[paste0(r, "_mean")]] <- tw_mean(t, m[[r]])
    row[[paste0(r, "_min")]]  <- min(m[[r]])
    row[[paste0(r, "_max")]]  <- max(m[[r]])
  }
  k1 <- m$ke + m$km
  row$BCF_parent_mean_L_kg  <- tw_mean(t, m$ku / k1)
  row$t_half_parent_mean_d  <- tw_mean(t, log(2) / k1)
  row$t_half_BTP_mean_d     <- tw_mean(t, log(2) / m$kem)
  row
}))
RATES_TAB <- RATES_TAB[order(match(RATES_TAB$version, MODEL_VERSIONS), match(RATES_TAB$treatment, TR)), ]
REF_ROW <- data.frame(version = "reference parameters", treatment = REFERENCE_RUN,
                      T_mean_C = T_REF_C, T_min_C = T_REF_C, T_max_C = T_REF_C)
for (r in RATE_NAMES) { REF_ROW[[paste0(r, "_mean")]] <- REF[[r]]; REF_ROW[[paste0(r, "_min")]] <- REF[[r]]
REF_ROW[[paste0(r, "_max")]] <- REF[[r]] }
REF_ROW$BCF_parent_mean_L_kg <- REF[["ku"]] / REF[["k1"]]
REF_ROW$t_half_parent_mean_d <- log(2) / REF[["k1"]]
REF_ROW$t_half_BTP_mean_d    <- log(2) / REF[["kem"]]
RATES_TAB <- rbind(REF_ROW, RATES_TAB); rownames(RATES_TAB) <- NULL

## 9.2 Settings and notes --------------------------------------------------------------------
SETTINGS <- data.frame(setting = c(
  "data file", "temperature file", "rate-constant file", "Arrhenius file", "experiment",
  "treatments", "species (tissue)", "model versions",
  "reference run", "T_ref [degC]", "T_A source",
  paste0("T_A of ", TA_RATES, " [K]"), paste0("SE of T_A of ", TA_RATES, " [K]"),
  "T_A used for ku / ke / km / kem", "intercept lines from (ku / ke / km / kem)", "intercept runs",
  "CI: SE of T_A", "CI: AE0+AE1 parameter sets", "CI: SE of the Arrhenius line height",
  "parameter-set file", "Monte Carlo draws", "CI level", "random seed", "max. RK4 step [d]", "output resolution [d]",
  "data transformation (R2, MLL)", "PPC acceptable [%]", "NRMSE acceptable [%]",
  "X0 parent / BTP [ug/kg ww]", "molar mass parent / BTP [g/mol]", "R version", "date"),
  value = c(DATA_FILE, TEMPERATURE_FILE, RATE_CONSTANT_FILE, ARRHENIUS_FILE, EXPERIMENT,
            paste(TR, collapse = ", "), SPECIES, paste(MODEL_VERSIONS, collapse = ", "),
            REFERENCE_RUN, T_REF_C, TA_SOURCE, TA[TA_RATES, "T_A"], TA[TA_RATES, "SE"],
            paste(TA_FOR[RATE_NAMES], collapse = " / "), paste(INTERCEPT_FROM[RATE_NAMES], collapse = " / "),
            paste(runs, collapse = ", "), CI_TA_SE, CI_REF_PARAMETERS, CI_LINE_HEIGHT,
            if (USE_PSAMPLE) PARAMETER_SAMPLE_FILE else "not used", N_MC, CI_LEVEL, RANDOM_SEED, STEP_MAX_D, OUTPUT_DT,
            DATA_TRANSFORM, PPC_MIN, NRMSE_MAX, paste(X0_PARENT, X0_BTP, sep = " / "),
            paste(PARENT_MOLAR_MASS, BTP_MOLAR_MASS, sep = " / "), R.version.string,
            format(Sys.time(), "%Y-%m-%d %H:%M")), stringsAsFactors = FALSE)
ANCH_TAB <- do.call(rbind, lapply(MODEL_VERSIONS, function(v)
  data.frame(version = v, ANCHORS[[v]][, c("rate", "TA_rate", "T_A", "SE_T_A", "k_anchor",
                                           "T_anchor_C", "lnA", "line", "SE_ln_anchor", "how")])))
rownames(ANCH_TAB) <- NULL
NOTES <- data.frame(note = c(
  "Model: dCp/dt = ku(T) Cw_p(t) - (ke(T) + km(T)) Cp;  dCm/dt = km(T) Cp - kem(T) Cm  (molar units)",
  "Arrhenius: k(T) = A exp(-T_A/T) = k_anchor exp(-T_A (1/T - 1/T_anchor)); T = water temperature of the logger (linear interpolation)",
  "reference: anchor = AE0+AE1 rate constants at T_ref; intercept: anchor = weighted Arrhenius line (w = 1/SE^2 of ln k), see sheet anchors",
  paste0("CI band: Monte Carlo (sheet MC_draws); band = 2.5-97.5% of the draws. Sources: ",
         paste(paste0(MODEL_VERSIONS, ": ", sapply(DRAWS, function(d) if (length(d$sources)) paste(d$sources, collapse = " + ") else "none")), collapse = "; ")),
  "T_A (SE): T_A ~ normal(mean, SE), ku and k1 independent; the lines rotate around their anchor (T_ref or weighted centroid)",
  "AE0+AE1 parameter sets: whole sets (ku, ke, km, kem) drawn uniformly from the 95% region of the calibration (keeps correlations)",
  "Arrhenius line height (SE): ln k at the weighted centroid ~ normal(0, sigma_hat * sqrt(1/sum w)) per line (as lm() in TK_arrhenius.R); ke and km share the k1 line",
  "Exposure: mean water concentration of the parent per sampling time, linear interpolation, constant after the last sample",
  paste0("Tissue data: ", SPECIES, " only (other Species entries ignored)"),
  "R2 and MLL on the transformed scale (DATA_TRANSFORM); AIC = 2 MLL (no fitted parameters): compare versions only",
  "NRMSE = RMSE / mean(obs) x 100 (acceptable <= NRMSE_MAX); PPC = % obs inside the CI band (acceptable >= PPC_MIN; CI band only, no residual error)",
  "PE = (pred - obs)/obs x 100 at the last tissue sampling time (mean of replicates; TK analogue of the EFSA SPPE); pooled = mean of the treatments; acceptable within +/- PE_MAX_ABS (no formal EFSA cut-off; +/-20% as in the EFSA worked example)",
  "Criteria: EFSA PPR Panel (2018) EFSA Journal 16(8):5377, section 4.1.4.5: PPC < 50% indicates poor performance; NRMSE should not exceed 50%",
  "model_timecourse: best estimate (mean T_A) with CI band, in ug and umol; ku..kem = temperature-corrected rate constants (best estimate)"))

MC_TAB <- do.call(rbind, lapply(MODEL_VERSIONS, function(v) {
  D <- DRAWS[[v]]$D
  data.frame(version = v, draw = seq_len(N_MC),
             setNames(as.data.frame(D$TA), paste0("T_A_", names(D$TA))),
             setNames(as.data.frame(D$KA), paste0(names(D$KA), "_anchor")),
             parameter_set = DRAWS[[v]]$set_id)
}))

## 9.3 Write ------------------------------------------------------------------------------------
xl_file <- out_file("validation_results.xlsx")
writexl::write_xlsx(list(
  quality_of_fit     = QOF,
  version_comparison = VCOMP,
  rate_constants     = RATES_TAB,
  anchors            = ANCH_TAB,
  model_timecourse   = MODEL,
  obs_vs_pred        = OVP,
  exposure_input     = EXPO,
  temperature_input  = TEMP_DF,
  MC_draws           = MC_TAB,
  settings           = SETTINGS,
  notes              = NOTES), xl_file)
say("Excel saved  : ", basename(xl_file))

writeLines(LOG, out_file("log.txt"))
say("Log saved    : ", RUN_NAME, "_log.txt")
say("All results in: ", OUT)