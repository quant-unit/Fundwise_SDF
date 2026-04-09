# ==============================================================================
# plot_degenerate_alpha.R
# ==============================================================================
# Objective: 
#   Calculate and visualize the objective function landscape (`err.sqr.calc`) 
#   over a 2D parameter grid of Alpha (x-axis) and Market Beta (y-axis).
#
# Idea:
#   By evaluating simulated Cash Flows (from a Data Generating Process) against 
#   historical market factors over a predefined mesh, we can visually detect 
#   "degenerate" regions of Alpha—specifically contrasting the fully scaled 
#   discount logic (max.month = 180) against the unscaled raw lifetime 
#   aggregate discounted cash flows (max.month = 0).
#
# Method:
#   1. Merges simulated private CFs (e.g., base_case_vyp) with public q_factors.
#   2. Rapidly pads CF data to monthly granularity for a single sample.
#   3. Evaluates a L2-Lasso objective function across an (N x N) Alpha/Beta grid.
#   4. Uses base R persp() to emit an academic-quality PDF containing side-by-side
#      topographical 3D charts with continuous color mapping.
# ==============================================================================

# 1. Setup Data Paths and Parameters ------------------------------------------

setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
print(getwd())

# Set working directory to r_project
if (basename(getwd()) == "simulation") {
  setwd("..")
}

# Source the NPV calculation helper
source("helper/getNPVs.R")

# Parameters
cache.folder.tag <- "base_case_vyp"
simulation.filename <- "base_case_vyp/base_case_vyp_simulated_cashflows_EW_VYP.csv"
public.filename <- "q_factors"

# Grid parameters (can be easily increased)
grid_n <- 100
alpha_seq <- seq(-0.01, 0.01, length.out = grid_n)
beta_seq <- seq(0, 2, length.out = grid_n)

# Plot viewing angles (rotate the chart)
# theta: Azimuthal rotation (left/right)
# phi: Colatitude rotation (up/down)
view_theta <- 25
view_phi <- 10  

max.months <- c(0, 180)
lambda <- 0

# True parameters for 'base_case_vyp' DGP
true_alpha <- 0
true_beta <- 1

# 2. Data Loading -------------------------------------------------------------
cat("Loading Public and Private Data...\n")

# Load public data
df.public <- read.csv(paste0("empirical/data_prepared_2026/", public.filename, ".csv"))
colnames(df.public) <- gsub("_World", "", colnames(df.public))
df.public$Date <- as.Date(df.public$Date)

# Normalize dates to last calendar day of month
df.public$Date <- as.Date(format(df.public$Date, "%Y-%m-01")) 
df.public$Date <- seq.Date(df.public$Date[1], by = "month", length.out = nrow(df.public))
df.public$Date <- as.Date(format(df.public$Date + 32, "%Y-%m-01")) - 1

df.public$Alpha <- 1

# Load private data (simulated cash flows)
df.private.cfs <- read.csv(paste0("simulation/data_prepared_sim/", simulation.filename))
df.private.cfs$Date <- as.Date(df.private.cfs$Date)
df.private.cfs$type <- as.factor(as.character(df.private.cfs$type))
df.private.cfs$Fund.ID <- as.factor(paste(df.private.cfs$Fund.ID, df.private.cfs$type, sep = "_"))

# Convert to monthly (padding zeros)
to.monthly <- function(df.ss) {
  max.date <- as.Date("2020-01-01")
  df.m <- data.frame(Date = seq(min(df.ss$Date) + 1, max.date, by = "month") - 1)
  df.ss <- merge(df.ss, df.m, by = "Date", all = TRUE)
  df.ss$CF[is.na(df.ss$CF)] <- 0
  for (col in c("type", "Vintage", "Fund.ID")) {
    if (col %in% colnames(df.ss)) {
      df.ss[, col] <- df.ss[1, col]
    }
  }
  return(df.ss)
}

# Subset to a single sample type *before* padding to save time and memory!
type_to_use <- levels(df.private.cfs$type)[1]
df.private.cfs <- df.private.cfs[df.private.cfs$type == type_to_use, ]

# Crucial: Drop unused levels so split() doesn't create empty data frames
df.private.cfs <- droplevels(df.private.cfs)

cat("Padding to monthly cash flows...\n")
df.private.cfs <- as.data.frame(data.table::rbindlist(lapply(split(df.private.cfs, df.private.cfs$Fund.ID), to.monthly)))

df.private.cfs$type <- as.factor(as.character(df.private.cfs$type))
df.private.cfs$Fund.ID <- as.factor(as.character(df.private.cfs$Fund.ID))

# Merge private and public
df0 <- base::merge(df.private.cfs, df.public, by = "Date", all.x = TRUE)
df0$Fund.ID <- as.factor(df0$Fund.ID)
rm(df.private.cfs)

# Final subset 
df.in <- df0
df.in$Fund.ID <- as.character(df.in$Fund.ID)

cat("Data prepared.\n")

# 3. Objective Function (L2_Lasso, Linear SDF) --------------------------------

# Evaluates the Net Present Value (NPV) logic for the cash flows.
f1 <- function(df.ss, max.month, par0) {
  if (max.month == 0) {
    # NOTE: When max.month = 0, we calculate the TRUE UNSCALED aggregate sum
    # of the discounted cash flows over the fund's entire lifetime. 
    # This intentionally removes the `/ max.month` divisor seen in the main 
    # estim_model_optimized.R script to avoid dividing by 0 (which results in Inf or 1e10).
    return(
      sum(
        df.ss$CF / exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0))))
      ) 
    )
  }
  # For max.month > 0, calculate the appropriately averaged NPV using the C++ helper.
  return(getNPVs(
    df.ss$CF,
    exp(cumsum(log(1 + (as.matrix(df.ss[, names(par0)]) %*% par0)))),
    max.month
  ) / max.month)
}

err.sqr.calc <- function(par, max.month, lambda, df) {
  dfx <- split(df, df$Fund.ID)
  npvs <- sapply(dfx, f1, max.month = max.month, par0 = c("RF" = 1, par))
  if (any(!is.finite(npvs))) {
    return(1e10)
  }
  p <- ifelse(length(par) == 1, par, par[-1])
  return(
    sqrt(sum(npvs^2)) / length(npvs) + lambda / length(p) * sum(abs(p))
  )
}

# 4. Grid Calculation ---------------------------------------------------------
cat("Calculating objective function over Alpha/Beta grid...\n")

calculate_grid_surface <- function(df, alpha_seq, beta_seq, max.month, lambda) {
  z_matrix <- matrix(NA, nrow = length(alpha_seq), ncol = length(beta_seq))
  
  for (i in seq_along(alpha_seq)) {
    for (j in seq_along(beta_seq)) {
      par <- c("MKT" = beta_seq[j], "Alpha" = alpha_seq[i])
      z_matrix[i, j] <- err.sqr.calc(par, max.month, lambda, df)
    }
  }
  return(z_matrix)
}

z_0   <- calculate_grid_surface(df.in, alpha_seq, beta_seq, max.months[1], lambda)
z_180 <- calculate_grid_surface(df.in, alpha_seq, beta_seq, max.months[2], lambda)

cat("Grid calculation completed.\n")

# 5. Plotting and Export ------------------------------------------------------
cat("Generating 3D plot to PDF...\n")

output_pdf <- "simulation/degenerate_alpha_grid.pdf"

pdf(output_pdf, width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))

# Helper to plot one side
plot_surface <- function(z_mat, max_m) {
  # True parameter projection
  true_z <- err.sqr.calc(c("MKT" = true_beta, "Alpha" = true_alpha), max_m, lambda, df.in)
  
  # Ensure valid z limits even if surface is completely flat (e.g. 1e10)
  z_range <- range(z_mat, na.rm = TRUE)
  if (diff(z_range) == 0) {
    z_range <- z_range + c(-1, 1)
  }
  
  # Compute faceted colors based on altitude (z)
  nrz <- nrow(z_mat)
  ncz <- ncol(z_mat)
  zfacet <- (z_mat[-1, -1] + z_mat[-1, -ncz] + z_mat[-nrz, -1] + z_mat[-nrz, -ncz]) / 4
  
  if (diff(range(zfacet, na.rm=TRUE)) == 0) {
    facet_cols <- rep("lightblue", length(zfacet))
  } else {
    nbcol <- 100
    facetcol <- cut(zfacet, nbcol)
    col_palette <- hcl.colors(nbcol, "Spectral", rev=FALSE)
    facet_cols <- col_palette[facetcol]
  }
  
  # base R persp plot
  pmat <- persp(x = alpha_seq, y = beta_seq, z = z_mat, 
          zlim = z_range,
          xlab = "\nAlpha", ylab = "\nBeta (MKT)", zlab = "\nSquared Error",
          main = paste0("Objective Function (Horizon = ", max_m, ")"), 
          theta = view_theta, phi = view_phi, 
          ticktype = "detailed",
          col = facet_cols,
          border = NA, shade = 0.5)
  
  # Map 3D points to 2D using trans3d
  true_pt <- trans3d(true_alpha, true_beta, true_z, pmat)
  min_z <- min(z_mat, na.rm=TRUE)
  floor_pt <- trans3d(true_alpha, true_beta, min_z, pmat)
  
  # Add the true parameter point
  points(true_pt, pch = 19, cex = 1.5, col = "green")
  
  # Add a vertical drop line to the floor to visually anchor the point
  segments(floor_pt$x, floor_pt$y, true_pt$x, true_pt$y, col = "green", lwd = 2, lty = 2)
          
  # Add a text label
  text(true_pt$x, true_pt$y, labels = "True DGP", col = "green", pos = 3, cex = 1.2)
}

plot_surface(z_0, max.months[1])
plot_surface(z_180, max.months[2])

dev.off()

cat("Saved plot to:", output_pdf, "\n")
