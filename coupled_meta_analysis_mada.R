# ============================================================================
# COUPLED META-ANALYSIS FOREST PLOT - USING MADA PACKAGE
# ============================================================================
# The 'mada' package is specifically designed for diagnostic test accuracy
# meta-analysis and creates proper coupled forest plots
# ============================================================================

# 1. Load Required Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, mada, meta, grid, gridExtra)

# 2. Load Data
# ============================================================================
# OPTION A: Load from Excel file
# df <- read_excel("your_data.xlsx")

# OPTION B: Sample data for demonstration
df <- data.frame(
  ID = paste("Study", 1:6),
  TP = c(45, 23, 89, 56, 12, 67),
  FN = c(5, 7, 11, 4, 3, 8),
  TN = c(50, 20, 80, 55, 15, 60),
  FP = c(10, 5, 20, 10, 2, 15)
)

# 3. Data Preparation
# ============================================================================
df$TP <- as.integer(df$TP)
df$FN <- as.integer(df$FN)
df$TN <- as.integer(df$TN)
df$FP <- as.integer(df$FP)

# Create a madad object (diagnostic accuracy data)
# mada expects: TP, FN, FP, TN
dta_data <- data.frame(
  TP = df$TP,
  FN = df$FN,
  FP = df$FP,
  TN = df$TN
)
rownames(dta_data) <- df$ID

# 4. Calculate Estimates with CIs
# ============================================================================
calc_binom_ci <- function(events, total) {
  if (is.na(events) || is.na(total) || total <= 0) {
    return(c(est = NA, low = NA, upp = NA))
  }
  result <- binom.test(events, total)
  c(est = events/total,
    low = result$conf.int[1],
    upp = result$conf.int[2])
}

# Sensitivity
sens_results <- t(mapply(calc_binom_ci, df$TP, df$TP + df$FN))
df$sens_est <- sens_results[, "est"]
df$sens_low <- sens_results[, "low"]
df$sens_upp <- sens_results[, "upp"]

# Specificity
spec_results <- t(mapply(calc_binom_ci, df$TN, df$TN + df$FP))
df$spec_est <- spec_results[, "est"]
df$spec_low <- spec_results[, "low"]
df$spec_upp <- spec_results[, "upp"]

# Format text
df$Sens_Txt <- sprintf("%.2f [%.2f, %.2f]", df$sens_est, df$sens_low, df$sens_upp)
df$Spec_Txt <- sprintf("%.2f [%.2f, %.2f]", df$spec_est, df$spec_low, df$spec_upp)

# 5. Perform Bivariate Meta-Analysis (Reitsma Model)
# ============================================================================
# This is the recommended approach for DTA meta-analysis
# It accounts for correlation between sensitivity and specificity

cat("\n========== BIVARIATE (REITSMA) MODEL ==========\n")
fit_reitsma <- reitsma(dta_data)
print(summary(fit_reitsma))

# 6. Alternative: Univariate Meta-Analysis for each
# ============================================================================
# Using meta package for univariate pooling

m_sens <- metaprop(
  event = TP,
  n = TP + FN,
  data = df,
  studlab = ID,
  sm = "PFT",
  method.ci = "CP",
  method.tau = "DL",
  random = TRUE,
  common = FALSE
)

m_spec <- metaprop(
  event = TN,
  n = TN + FP,
  data = df,
  studlab = ID,
  sm = "PFT",
  method.ci = "CP",
  method.tau = "DL",
  random = TRUE,
  common = FALSE
)

cat("\n========== UNIVARIATE SENSITIVITY ==========\n")
print(summary(m_sens))

cat("\n========== UNIVARIATE SPECIFICITY ==========\n")
print(summary(m_spec))

# 7. Create Coupled Forest Plot using MADA
# ============================================================================
# This creates the standard coupled forest plot for DTA

pdf("mada_coupled_forest.pdf", width = 12, height = 7)
forest(madad(dta_data),
       type = "sens",
       main = "Coupled Forest Plot: Sensitivity and Specificity",
       snames = df$ID)
dev.off()

# 8. Create Custom Coupled Forest Plot with Full Statistics
# ============================================================================

# Add data to meta objects
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

# Create coupled plot function
create_coupled_dta_forest <- function(m_sens, m_spec, df) {

  # Open PDF device
  pdf("coupled_forest_full.pdf", width = 16, height = 9)

  # Set up layout: 3 columns (table, sens plot, spec plot)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(2.5, 1, 1))
  par(mar = c(5, 1, 4, 0))

  n_studies <- nrow(df)

  # --- Panel 1: Data Table ---
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, n_studies + 4))

  # Header row
  y_header <- n_studies + 3
  text(0.3, y_header, "Study", font = 2, adj = 0, cex = 0.9)
  text(1.8, y_header, "TP", font = 2, adj = 0.5, cex = 0.9)
  text(2.6, y_header, "FN", font = 2, adj = 0.5, cex = 0.9)
  text(3.4, y_header, "TN", font = 2, adj = 0.5, cex = 0.9)
  text(4.2, y_header, "FP", font = 2, adj = 0.5, cex = 0.9)
  text(6.0, y_header, "Sens [95% CI]", font = 2, adj = 0.5, cex = 0.9)
  text(8.5, y_header, "Spec [95% CI]", font = 2, adj = 0.5, cex = 0.9)

  # Data rows
  for (i in 1:n_studies) {
    y_pos <- n_studies - i + 2
    text(0.3, y_pos, df$ID[i], adj = 0, cex = 0.85)
    text(1.8, y_pos, df$TP[i], adj = 0.5, cex = 0.85)
    text(2.6, y_pos, df$FN[i], adj = 0.5, cex = 0.85)
    text(3.4, y_pos, df$TN[i], adj = 0.5, cex = 0.85)
    text(4.2, y_pos, df$FP[i], adj = 0.5, cex = 0.85)
    text(6.0, y_pos, df$Sens_Txt[i], adj = 0.5, cex = 0.85)
    text(8.5, y_pos, df$Spec_Txt[i], adj = 0.5, cex = 0.85)
  }

  # Pooled row
  y_pooled <- 0.5
  abline(h = 1.3, lty = 1)
  text(0.3, y_pooled, "Pooled", font = 2, adj = 0, cex = 0.9)
  pooled_sens_txt <- sprintf("%.2f [%.2f, %.2f]",
                              m_sens$TE.random, m_sens$lower.random, m_sens$upper.random)
  pooled_spec_txt <- sprintf("%.2f [%.2f, %.2f]",
                              m_spec$TE.random, m_spec$lower.random, m_spec$upper.random)
  text(6.0, y_pooled, pooled_sens_txt, adj = 0.5, cex = 0.85, font = 2)
  text(8.5, y_pooled, pooled_spec_txt, adj = 0.5, cex = 0.85, font = 2)

  # --- Panel 2: Sensitivity Forest Plot ---
  par(mar = c(5, 0, 4, 0.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, n_studies + 4))

  # Title
  text(0.5, n_studies + 3, "Sensitivity", font = 2, cex = 1)

  # X-axis
  axis(1, at = seq(0, 1, 0.2), cex.axis = 0.8)
  abline(v = 0.5, lty = 2, col = "gray60")

  # Study-level estimates
  for (i in 1:n_studies) {
    y_pos <- n_studies - i + 2

    # Point estimate (square)
    points(df$sens_est[i], y_pos, pch = 15, cex = 1.5, col = "gray30")

    # Confidence interval (line)
    lines(c(df$sens_low[i], df$sens_upp[i]), c(y_pos, y_pos), lwd = 2, col = "gray30")
  }

  # Pooled estimate (diamond)
  y_pooled <- 0.5
  abline(h = 1.3, lty = 1)

  diamond_x <- c(m_sens$lower.random, m_sens$TE.random, m_sens$upper.random, m_sens$TE.random)
  diamond_y <- c(y_pooled, y_pooled + 0.3, y_pooled, y_pooled - 0.3)
  polygon(diamond_x, diamond_y, col = "darkblue", border = "darkblue")

  # Heterogeneity stats at bottom
  het_text <- sprintf("I² = %.1f%%, τ² = %.3f, Q = %.2f (p = %.3f)",
                      m_sens$I2 * 100, m_sens$tau2, m_sens$Q, m_sens$pval.Q)
  test_text <- sprintf("Test: z = %.2f, p = %.4f", m_sens$zval.random, m_sens$pval.random)
  mtext(het_text, side = 1, line = 2.5, cex = 0.7)
  mtext(test_text, side = 1, line = 3.5, cex = 0.7)

  # --- Panel 3: Specificity Forest Plot ---
  par(mar = c(5, 0.5, 4, 1))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, n_studies + 4))

  # Title
  text(0.5, n_studies + 3, "Specificity", font = 2, cex = 1)

  # X-axis
  axis(1, at = seq(0, 1, 0.2), cex.axis = 0.8)
  abline(v = 0.5, lty = 2, col = "gray60")

  # Study-level estimates
  for (i in 1:n_studies) {
    y_pos <- n_studies - i + 2

    # Point estimate (square)
    points(df$spec_est[i], y_pos, pch = 15, cex = 1.5, col = "gray30")

    # Confidence interval (line)
    lines(c(df$spec_low[i], df$spec_upp[i]), c(y_pos, y_pos), lwd = 2, col = "gray30")
  }

  # Pooled estimate (diamond)
  y_pooled <- 0.5
  abline(h = 1.3, lty = 1)

  diamond_x <- c(m_spec$lower.random, m_spec$TE.random, m_spec$upper.random, m_spec$TE.random)
  diamond_y <- c(y_pooled, y_pooled + 0.3, y_pooled, y_pooled - 0.3)
  polygon(diamond_x, diamond_y, col = "darkred", border = "darkred")

  # Heterogeneity stats at bottom
  het_text <- sprintf("I² = %.1f%%, τ² = %.3f, Q = %.2f (p = %.3f)",
                      m_spec$I2 * 100, m_spec$tau2, m_spec$Q, m_spec$pval.Q)
  test_text <- sprintf("Test: z = %.2f, p = %.4f", m_spec$zval.random, m_spec$pval.random)
  mtext(het_text, side = 1, line = 2.5, cex = 0.7)
  mtext(test_text, side = 1, line = 3.5, cex = 0.7)

  dev.off()

  cat("\nCustom coupled forest plot saved to: coupled_forest_full.pdf\n")
}

# Generate the custom plot
create_coupled_dta_forest(m_sens, m_spec, df)

# 9. Summary Output
# ============================================================================
cat("\n\n============================================================\n")
cat("                    POOLED RESULTS SUMMARY                    \n")
cat("============================================================\n\n")

cat("UNIVARIATE RESULTS:\n")
cat("-" , rep("-", 50), "\n", sep = "")
cat("\nSENSITIVITY:\n")
cat(sprintf("  Pooled estimate: %.3f (95%% CI: %.3f - %.3f)\n",
            m_sens$TE.random, m_sens$lower.random, m_sens$upper.random))
cat(sprintf("  Heterogeneity: I² = %.1f%%, τ² = %.4f\n", m_sens$I2*100, m_sens$tau2))
cat(sprintf("  Cochran's Q = %.2f, df = %d, p = %.4f\n", m_sens$Q, m_sens$df.Q, m_sens$pval.Q))
cat(sprintf("  Test for effect: z = %.2f, p = %.4f\n", m_sens$zval.random, m_sens$pval.random))

cat("\nSPECIFICITY:\n")
cat(sprintf("  Pooled estimate: %.3f (95%% CI: %.3f - %.3f)\n",
            m_spec$TE.random, m_spec$lower.random, m_spec$upper.random))
cat(sprintf("  Heterogeneity: I² = %.1f%%, τ² = %.4f\n", m_spec$I2*100, m_spec$tau2))
cat(sprintf("  Cochran's Q = %.2f, df = %d, p = %.4f\n", m_spec$Q, m_spec$df.Q, m_spec$pval.Q))
cat(sprintf("  Test for effect: z = %.2f, p = %.4f\n", m_spec$zval.random, m_spec$pval.random))

cat("\n\nBIVARIATE (REITSMA) RESULTS:\n")
cat("-" , rep("-", 50), "\n", sep = "")
cat("(Recommended for DTA meta-analysis - accounts for correlation)\n")
cat(sprintf("  Pooled Sensitivity: %.3f (95%% CI: %.3f - %.3f)\n",
            fit_reitsma$coefficients["tsens", "Estimate"],
            fit_reitsma$coefficients["tsens", "Estimate"] - 1.96 * fit_reitsma$coefficients["tsens", "Std. Error"],
            fit_reitsma$coefficients["tsens", "Estimate"] + 1.96 * fit_reitsma$coefficients["tsens", "Std. Error"]))
cat(sprintf("  Pooled Specificity: %.3f (95%% CI: %.3f - %.3f)\n",
            fit_reitsma$coefficients["tfpr", "Estimate"],
            fit_reitsma$coefficients["tfpr", "Estimate"] - 1.96 * fit_reitsma$coefficients["tfpr", "Std. Error"],
            fit_reitsma$coefficients["tfpr", "Estimate"] + 1.96 * fit_reitsma$coefficients["tfpr", "Std. Error"]))

cat("\n============================================================\n")
cat("Output files generated:\n")
cat("  1. mada_coupled_forest.pdf    - Standard MADA forest plot\n")
cat("  2. coupled_forest_full.pdf    - Custom coupled plot with all stats\n")
cat("============================================================\n")
