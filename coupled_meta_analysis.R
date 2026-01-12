# ============================================================================
# COUPLED META-ANALYSIS FOREST PLOT FOR SENSITIVITY AND SPECIFICITY
# ============================================================================
# This script creates a coupled forest plot with:
# - Data columns: TN, FP, TP, FN
# - Sensitivity and Specificity estimates with 95% CI
# - Pooled estimates with diamonds
# - Heterogeneity statistics (I², τ², Q)
# - Test statistics for pooled effects
# ============================================================================

# 1. Load Required Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, meta, grid, gridExtra, mada)

# 2. Load Data
# ============================================================================
# OPTION A: Load from Excel file (uncomment and modify path as needed)
# df <- read_excel("your_data.xlsx")

# OPTION B: Sample data for demonstration
df <- data.frame(
  ID = paste("Study", 1:6),
  TP = c(45, 23, 89, 56, 12, 67),
  FN = c(5, 7, 11, 4, 3, 8),
  TN = c(50, 20, 80, 55, 15, 60),
  FP = c(10, 5, 20, 10, 2, 15)
)

# If your Excel has pre-calculated estimates, you can use them
# Otherwise, we calculate from TP, FN, TN, FP below

# 3. Data Preparation
# ============================================================================
df$TP <- as.integer(df$TP)
df$FN <- as.integer(df$FN)
df$TN <- as.integer(df$TN)
df$FP <- as.integer(df$FP)

# Calculate totals
df$n_sens <- df$TP + df$FN  # Total diseased
df$n_spec <- df$TN + df$FP  # Total non-diseased

# Calculate estimates and CIs using exact binomial method
calc_ci <- function(events, total) {
  if (is.na(events) || is.na(total) || total <= 0) {
    return(c(est = NA, low = NA, upp = NA))
  }
  result <- binom.test(events, total)
  c(est = events/total,
    low = result$conf.int[1],
    upp = result$conf.int[2])
}

# Calculate sensitivity estimates
sens_ci <- t(mapply(calc_ci, df$TP, df$n_sens))
df$sens_est <- sens_ci[, "est"]
df$sens_low <- sens_ci[, "low"]
df$sens_upp <- sens_ci[, "upp"]

# Calculate specificity estimates
spec_ci <- t(mapply(calc_ci, df$TN, df$n_spec))
df$spec_est <- spec_ci[, "est"]
df$spec_low <- spec_ci[, "low"]
df$spec_upp <- spec_ci[, "upp"]

# Format for display in table
fmt_est_ci <- function(est, low, upp) {
  if (any(is.na(c(est, low, upp)))) return(" ")
  sprintf("%.2f [%.2f, %.2f]", est, low, upp)
}

df$Sens_Txt <- mapply(fmt_est_ci, df$sens_est, df$sens_low, df$sens_upp)
df$Spec_Txt <- mapply(fmt_est_ci, df$spec_est, df$spec_low, df$spec_upp)

# 4. Perform Meta-Analysis (Random Effects)
# ============================================================================
# Using Freeman-Tukey double arcsine transformation for better pooling
# This handles proportions near 0 or 1 better

m_sens <- metaprop(
  event = TP,
  n = n_sens,
  data = df,
  studlab = ID,
  sm = "PFT",           # Freeman-Tukey transformation (better for proportions)
  method.ci = "CP",     # Clopper-Pearson exact CI
  method.tau = "DL",    # DerSimonian-Laird for tau²
  random = TRUE,
  common = FALSE
)

m_spec <- metaprop(
  event = TN,
  n = n_spec,
  data = df,
  studlab = ID,
  sm = "PFT",
  method.ci = "CP",
  method.tau = "DL",
  random = TRUE,
  common = FALSE
)

# 5. Print Summary Statistics
# ============================================================================
cat("\n========== SENSITIVITY META-ANALYSIS ==========\n")
print(summary(m_sens))

cat("\n========== SPECIFICITY META-ANALYSIS ==========\n")
print(summary(m_spec))

# 6. Create Coupled Forest Plot
# ============================================================================

# Add the text columns to the meta objects for display
m_sens$data$TN <- df$TN
m_sens$data$FP <- df$FP
m_sens$data$TP <- df$TP
m_sens$data$FN <- df$FN
m_sens$data$Sens_Txt <- df$Sens_Txt
m_sens$data$Spec_Txt <- df$Spec_Txt

m_spec$data$TN <- df$TN
m_spec$data$FP <- df$FP
m_spec$data$TP <- df$TP
m_spec$data$FN <- df$FN
m_spec$data$Sens_Txt <- df$Sens_Txt
m_spec$data$Spec_Txt <- df$Spec_Txt

# Function to create forest plot panels
create_forest_panel <- function(model, side = "left", main_title = "") {
  gc <- "gray30"  # Color for squares and diamonds

  if (side == "left") {
    # Left panel: Table columns + Sensitivity forest plot
    forest(model,
           # Table columns on the left
           leftcols = c("studlab", "TP", "FN", "TN", "FP", "Sens_Txt", "Spec_Txt"),
           leftlabs = c("Study", "TP", "FN", "TN", "FP",
                        "Sensitivity [95% CI]", "Specificity [95% CI]"),

           # No columns on right (graph only)
           rightcols = FALSE,

           # X-axis settings
           xlab = "Sensitivity",
           xlim = c(0, 1),

           # Appearance
           col.square = gc,
           col.square.lines = gc,
           col.diamond = "darkblue",
           col.diamond.lines = "darkblue",

           # CRITICAL: Force display of pooled estimate and statistics
           overall = TRUE,
           overall.hetstat = TRUE,
           hetstat = TRUE,
           test.overall.random = TRUE,

           # Heterogeneity statistics to display
           print.I2 = TRUE,
           print.tau2 = TRUE,
           print.Q = TRUE,
           print.pval.Q = TRUE,

           # Labels
           text.random = "Pooled Sensitivity",
           smlab = "",

           # Layout settings
           plotwidth = unit(5, "cm"),
           colgap.studlab = unit(3, "mm"),
           colgap.left = unit(2, "mm"),

           # Font sizes
           fs.study = 9,
           fs.heading = 10,
           fs.random = 10,
           fs.hetstat = 9,

           # Spacing
           spacing = 1.0,

           # Backtransform from Freeman-Tukey
           backtransf = TRUE
    )
  } else {
    # Right panel: Specificity forest plot only (no table)
    forest(model,
           # No text columns
           leftcols = FALSE,
           rightcols = FALSE,
           studlab = FALSE,

           # X-axis settings
           xlab = "Specificity",
           xlim = c(0, 1),

           # Appearance
           col.square = gc,
           col.square.lines = gc,
           col.diamond = "darkred",
           col.diamond.lines = "darkred",

           # CRITICAL: Force display of pooled estimate and statistics
           overall = TRUE,
           overall.hetstat = TRUE,
           hetstat = TRUE,
           test.overall.random = TRUE,

           # Heterogeneity statistics
           print.I2 = TRUE,
           print.tau2 = TRUE,
           print.Q = TRUE,
           print.pval.Q = TRUE,

           # Labels
           text.random = "Pooled Specificity",
           smlab = "",

           # Layout settings
           plotwidth = unit(5, "cm"),

           # Font sizes
           fs.study = 9,
           fs.heading = 10,
           fs.random = 10,
           fs.hetstat = 9,

           # Spacing
           spacing = 1.0,

           # Backtransform
           backtransf = TRUE
    )
  }
}

# 7. Generate Output
# ============================================================================

# Option A: Save to PDF (recommended for publication)
pdf("coupled_forest_plot.pdf", width = 14, height = 8)
grid.newpage()

# Capture both panels
g_sens <- grid.grabExpr(create_forest_panel(m_sens, "left"), wrap = TRUE)
g_spec <- grid.grabExpr(create_forest_panel(m_spec, "right"), wrap = TRUE)

# Arrange side by side
grid.arrange(g_sens, g_spec, ncol = 2, widths = c(2.5, 1))
dev.off()

# Option B: Display in R viewer
grid.newpage()
g_sens <- grid.grabExpr(create_forest_panel(m_sens, "left"), wrap = TRUE)
g_spec <- grid.grabExpr(create_forest_panel(m_spec, "right"), wrap = TRUE)
grid.arrange(g_sens, g_spec, ncol = 2, widths = c(2.5, 1))

# 8. Print Pooled Results Summary
# ============================================================================
cat("\n\n============================================================\n")
cat("                    POOLED RESULTS SUMMARY                    \n")
cat("============================================================\n\n")

cat("SENSITIVITY:\n")
cat(sprintf("  Pooled estimate: %.3f (95%% CI: %.3f - %.3f)\n",
            m_sens$TE.random, m_sens$lower.random, m_sens$upper.random))
cat(sprintf("  Heterogeneity: I² = %.1f%%, τ² = %.4f\n", m_sens$I2*100, m_sens$tau2))
cat(sprintf("  Cochran's Q = %.2f, df = %d, p = %.4f\n", m_sens$Q, m_sens$df.Q, m_sens$pval.Q))
cat(sprintf("  Test for overall effect: z = %.2f, p = %.4f\n\n", m_sens$zval.random, m_sens$pval.random))

cat("SPECIFICITY:\n")
cat(sprintf("  Pooled estimate: %.3f (95%% CI: %.3f - %.3f)\n",
            m_spec$TE.random, m_spec$lower.random, m_spec$upper.random))
cat(sprintf("  Heterogeneity: I² = %.1f%%, τ² = %.4f\n", m_spec$I2*100, m_spec$tau2))
cat(sprintf("  Cochran's Q = %.2f, df = %d, p = %.4f\n", m_spec$Q, m_spec$df.Q, m_spec$pval.Q))
cat(sprintf("  Test for overall effect: z = %.2f, p = %.4f\n", m_spec$zval.random, m_spec$pval.random))

cat("\n============================================================\n")
cat("Forest plot saved to: coupled_forest_plot.pdf\n")
cat("============================================================\n")
