#' @title .ComputeFeatureThreshold
#'
#' @description
#' Computes an automatic threshold on normalized expression a feature to filter
#' cells based on this value. This function has been adapted from the BaderLab
#' MALAT1_threshold R module (see references) and is intented to apply the MALAT1
#' threshold approach developed by BaderLab (see reference).
#' The function estimates the kernel density of a given feature, fits a
#' smoothing spline, locates local maxima and minima, and then uses a quadratic
#' approximation around the main cell peak to derive a robust cutoff.
#'
#' @details
#' MALAT1 is a nuclear-retained lncRNA whose expression correlates strongly with
#' nuclear (intronic) RNA content and splice ratio, making it a convenient proxy
#' for cell quality in droplet-based scRNA-seq.
#' Low MALAT1 values typically correspond to empty droplets, cytosolic debris,
#' or damaged cells, whereas higher values correspond to intact nuclei.
#' This function assumes you provide normalized MALAT1 counts from a single,
#' unintegrated sample, and it attempts to:
#' \itemize{
#'   \item Identify the main MALAT1 peak representing real cells.
#'   \item Find a local minimum separating low-MALAT1 droplets from that peak.
#'   \item Fit a quadratic to the density around the peak and solve for the
#'         left intercept as the threshold.
#' }
#' Robustness safeguards include minimum and rough peak positions (`abs.min`,
#' `rough.max`) and explicit error handling for atypical distributions (e.g.
#' entire sample of low-quality droplets).
#'
#' NOTE: This parameters are set to default values that work well for MALAT1,
#' but can be adjusted for other features.
#'
#' @param counts.vector Numeric vector of normalized feature expression values
#'   for one sample; typically log-normalized counts per cell from an
#'   unintegrated dataset containing multiple cell types. Note that the rest of
#'   the parameters are tuned for MALAT1 and may need adjustment for other features.
#' @param bw.bandwidth Numeric scalar passed to \code{stats::density()} as the
#'   kernel bandwidth; default is \code{0.01}. Smaller values (e.g. \code{0.001})
#'   make the density line less stiff and more closely track the histogram.
#' @param chosen.min Numeric scalar giving the minimum MALAT1 value above which
#'   a density peak is considered a candidate for the real cell peak; default
#'   \code{2} - fine-tuned for MALAT1. This helps ignore spurious peaks at or
#'   near zero due to empty droplets or ambient RNA.
#' @param smooth.spar Numeric smoothing parameter (\code{spar}) passed to
#'   \code{stats::smooth.spline()}, controlling the trade-off between smoothness
#'   and fidelity to the density curve; default \code{1} - fine-tuned for MALAT1.
#' @param abs.min Numeric scalar specifying the absolute minimum allowed
#'   MALAT1 threshold. Protects against thresholds collapsing to zero in
#'   pathological cases; default \code{1} - fine-tned for MALAT1.
#' @param rough.max Numeric scalar giving a rough expected position of the
#'   MALAT1 peak when no local maximum can be found. The closest x-value in the
#'   density to this value is used as a surrogate peak; default \code{6}
#'   - fine-tuned for MALAT1.
#' @param conservative.threshold Numeric scalar specifying a conservative default
#'   threshold to return in case of errors (e.g. no high MALAT1 peaks). Default is \code{2}
#'   - fine-tuned for MALAT1.
#'
#' @import stats
#' @return threshold Numeric scalar: the MALAT1 threshold (left quadratic intercept,
#' constrained to be at least \code{abs.min}). On error, the function returns the
#' numeric value \code{2} as a conservative default threshold.
#'
#' @noRd
#' @keywords internal
#'
#' @references
#' This function has been adapted from the BaderLab MALAT1_threshold R module.
#'
#' Clarke, Bader et al. "MALAT1 expression indicates cell quality in
#' single-cell RNA sequencing data." bioRxiv (2024).
#' \url{https://www.biorxiv.org/content/10.1101/2024.07.14.603469v2}
#'
#' BaderLab MALAT1_threshold GitHub repository:
#' \url{https://github.com/BaderLab/MALAT1_threshold}
#'
#' @examples
#' \dontrun{
#'   # Assume malat1_norm is a numeric vector of normalized MALAT1 counts
#'   threshold.res <- .ComputeFeatureThreshold(
#'     counts.vector    = malat1_norm,
#'     bw.bandwidth     = 0.01,
#'     chosen.min       = 2,
#'     smooth.spar      = 2,
#'     abs.min          = 1,
#'     rough.max        = 3,
#'     conservative.threshold = 2
#'   )
#' }
#'
#' @seealso
#' \itemize{
#'   \item DropletQC: nuclear fraction-based QC for empty droplets and damaged cells.
#'   \item EmptyDrops: ambient RNA-based empty droplet detection.
#' }

.ComputeFeatureThreshold <- function(
  counts.vector,
  bw.bandwidth = 0.01,
  chosen.min = 2,
  smooth.spar = 1,
  abs.min = 1,
  rough.max = 6,
  conservative.threshold = 2
) {
  # Wrap the entire computation in tryCatch to handle atypical input gracefully.
  tryCatch(
    expr = {
      # Compute kernel density estimate of normalized MALAT1 counts over their observed range.
      density.data <- stats::density(
        x = counts.vector, # Use the input MALAT1 counts as data.
        bw = bw.bandwidth, # Apply user-specified bandwidth for smoothness.
        from = min(counts.vector), # Start the density at the minimum observed value.
        to = max(counts.vector) # End the density at the maximum observed value.
      )

      # Fit a smoothing spline to the density curve to obtain a smooth representation.
      fit.spline <- stats::smooth.spline(
        x = density.data$x, # Use the density x-coordinates as spline input.
        y = density.data$y, # Use the density y-values as spline targets.
        spar = smooth.spar # Smoothing parameter controlling fit flexibility.
      )

      # Compute the first derivative of the spline at each x to detect local maxima and minima.
      first.derivative <- stats::predict(
        object = fit.spline, # Use the fitted spline as the function to differentiate.
        x = density.data$x, # Evaluate the derivative at all density x positions.
        deriv = 1 # Request the first derivative.
      )

      # Identify local maxima where the derivative changes sign from positive to negative.
      local.maxima <- density.data$x[
        which(
          diff(sign(first.derivative$y)) == -2 # Sign flip +1 -> -1 indicates a peak.
        )
      ]

      # If no local maxima are found, fall back to the x-value closest to rough.max as a surrogate peak.
      if (length(local.maxima) == 0L) {
        # Compute absolute differences between each x and rough.max.
        abs.diff.rough <- abs(density.data$x - rough.max)
        # Find index of the minimal difference (closest point).
        index.closest.rough <- which(
          abs.diff.rough == min(abs.diff.rough)
        )
        # Use the corresponding x-value as the only local maximum.
        local.maxima <- density.data$x[index.closest.rough]
      }

      # Build data frame for the raw density to use in ggplot2 visualizations.
      density.df <- data.frame(
        x = density.data$x, # Density x-coordinates.
        y = density.data$y # Density y-values.
      )

      # Build data frame for the spline fit evaluated across the density x-grid.
      fit.df <- data.frame(
        x = density.data$x, # Same x grid as density.
        y = stats::predict(fit.spline, density.data$x)$y # Corresponding smoothed y-values.
      )

      # Build data frame for local maxima points.
      local.maxima.df <- data.frame(
        x = local.maxima, # x-positions of maxima.
        y = stats::predict(fit.spline, local.maxima)$y # Spline-predicted y at maxima.
      )

      # Identify local minima where the derivative changes sign from negative to positive.
      local.minima <- density.data$x[
        which(
          diff(sign(first.derivative$y)) == 2 # Sign flip -1 -> +1 indicates a valley.
        )
      ]

      # If no local minima are found, use abs.min as a fallback minimum position.
      if (length(local.minima) == 0L) {
        local.minima <- abs.min # Assign scalar fallback as minima.
      }

      # Select the largest local maximum above chosen.min to represent the real cell MALAT1 peak.
      biggest.y <- max(
        density.data$y[
          density.data$x %in% local.maxima[local.maxima > chosen.min]
        ]
      )

      # Find the x-position (max.index) of this peak in the density.
      max.index <- density.data$x[density.data$y == biggest.y]

      # Restrict local minima to those strictly to the left of the chosen peak.
      local.minima.left <- local.minima[local.minima < max.index]

      # If minima vector to the left is shorter than abs.min (original code’s robustness check), enforce abs.min.
      if (length(local.minima.left) < abs.min) {
        local.minima.left <- abs.min # Replace with scalar minimum threshold.
      }

      # Choose the local minimum closest to the peak from the left as the separating valley (min.index).
      min.index <- local.minima.left[
        (max.index - local.minima.left) == min(max.index - local.minima.left)
      ]

      # Compute delta as the horizontal distance between peak (max.index) and chosen minimum (min.index).
      delta <- max.index - min.index

      # Assemble a full data frame of density x and y values for subsetting and modelling.
      density.full.df <- data.frame(
        x = density.data$x, # All x-coordinates.
        y = density.data$y # All density values.
      )

      # Subset the density to a symmetric window around the peak of width 2 * delta.
      subset.df <- density.full.df[
        density.full.df[, "x"] >= (max.index - delta) &
          density.full.df[, "x"] <= (max.index + delta),
        ,
        drop = FALSE
      ]

      # Fit a quadratic model to the subset: y ~ x + x^2 via poly(x, 2, raw = TRUE).
      quad.model <- stats::lm(
        formula = y ~ stats::poly(x, 2, raw = TRUE), # Quadratic polynomial in x.
        data = subset.df # Use subset around peak for fitting.
      )

      # Extract coefficient estimates (a, b, c) from the quadratic model summary.
      coef.table <- summary(quad.model)$coefficients
      c.intercept <- coef.table[1, 1] # Constant term c.
      b.linear <- coef.table[2, 1] # Linear term coefficient b.
      a.quadratic <- coef.table[3, 1] # Quadratic term coefficient a.

      # Compute the discriminant of the quadratic: b^2 - 4ac.
      discriminant <- b.linear^2 - 4 * a.quadratic * c.intercept

      # Solve for the first x-intercept of the quadratic using the quadratic formula.
      x.intercept1 <- (-b.linear + sqrt(discriminant)) / (2 * a.quadratic)

      # If the computed intercept is negative, enforce abs.min as a lower bound for the threshold.
      if (x.intercept1 < 0) {
        x.intercept1 <- abs.min # Replace negative intercept with minimum allowed.
      }

      return(x.intercept1) # Return only the numeric threshold.
    },

    error = function(e) {
      # On error, emit a diagnostic message explaining common causes (e.g. input not normalized).
      warning(
        "An error occurred: Please make sure you have used a vector of normalized counts as input. ",
        "This may also indicate that you have no high feature peaks, meaning this particular sample may be poor quality (if feature =  'MALAT1'). ",
        "A conservative default threshold of ",
        conservative.threshold,
        " will be returned. Error details: ",
        e$message
      )

      # Return a conservative default threshold value (2 by default) to avoid failing downstream code.
      return(conservative.threshold)
    }
  )
}
