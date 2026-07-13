# Modules for feature threshold computation (feature-based QC thresholds for scRNA-seq in Seurat)

This repository implements robust, feature-based quality-control (QC) thresholds for single-cell RNA‑seq data, centered on the MALAT1-based approach originally developed by BaderLab (see refrences) and extended here with a modern, idiomatic R interface and additional modules.

The main goal is to provide:

- A clean, well‑documented R implementation of MALAT1‐based cell quality thresholds.
- Generalized threshold routines for other nuclear or QC marker genes (“features”) using the same density/spline framework directly on SeuratOjects.

> **Note:** This repository contains derivative work based on code and ideas from the BaderLab `MALAT1_threshold` project and the article “MALAT1 expression indicates cell quality in single-cell RNA sequencing data” by Clarke & Bader (2024).

---

## Background: MALAT1 and cell quality

Single‑cell RNA‑seq experiments routinely capture:

- **Empty droplets**, largely containing ambient cytosolic RNA.
- **Damaged or low‑quality cells**, with depleted nuclear RNA and altered splice ratios.
- **Real cells**, with intact nuclei and appropriate intronic (unspliced) content.

Many QC frameworks rely on **spliced vs unspliced read ratios** or nuclear fraction metrics to detect these problematic droplets, but these approaches require re‑examining all reads and can be computationally expensive.

Clarke & Bader showed that the expression of **MALAT1**, a long non‑coding RNA that is retained in the nucleus and expressed across most cell types, is strongly correlated with splice ratio and nuclear content. In practice:

- **Low MALAT1** expression tends to indicate empty droplets, ambient RNA–dominated droplets, or damaged cells.
- **High MALAT1** expression, with appropriate intronic fraction, tends to indicate good‑quality cells.

Because MALAT1 is a single gene, its expression distribution is easy to visualize and threshold, making it a fast proxy for more expensive splice‑ratio metrics.

---

## Overview of this implementation

This repository contains:

- A direct, documented implementation of the MALAT1‑based thresholding algorithm (density + smoothing spline + quadratic fit around the main peak) adapted from BaderLab and refactored to an idiomatic R style (CamelCase function names, dot.case arguments and variables, explicit `package::function` calls).
- Companion functions to compute thresholds for  MALAT1 or other QC features (e.g., nuclear marker genes, mitochondrial fraction, or custom features), reusing the same density/spline infrastructure directly on the SeuratObject. The functions are designed to be called from a single wrapper function (`CalculateQC()`) from the [SCutils](https://github.com/jmgs7/SCutils) package that computes multiple QC metrics in one go.

---

## Core MALAT1 threshold function

The core function in this repository implements the algorithm described by BaderLab, with the following main steps:

1. **Kernel density estimation**  
   - Compute a smooth kernel density estimate of normalized MALAT1 expression across all cells using `stats::density()` with a user‑defined bandwidth.  
   - The bandwidth controls how “stiff” or “wiggly” the density curve is relative to the histogram.

2. **Spline smoothing and derivative**  
   - Fit a smoothing spline (`stats::smooth.spline()`) to the density curve with parameter `spar` controlling smoothness.  
   - Compute the first derivative of the spline (`stats::predict(..., deriv = 1)`) to locate local maxima and minima in the smoothed density.

3. **Local maxima and minima detection**  
   - Identify **local maxima** (density peaks) where the derivative changes sign from positive to negative.  
   - Identify **local minima** (valleys) where the derivative flips from negative to positive.  
   - If no maxima are found, fall back to using the x‑value closest to a rough expected peak position (e.g. `rough.max = 6`).  
   - If no minima are found, fall back to a robust lower bound (e.g. `abs.min = 1`).

4. **Peak selection and valley to the left**  
   - Among local maxima, select the largest peak above a minimum x‑value (`chosen.min`, typically 2) as the main MALAT1 peak representing real cells.  
   - Restrict local minima to those strictly to the left of this peak, then choose the closest such minimum as the separating valley.

5. **Quadratic fit around the peak**  
   - Take a symmetric window around the selected peak, with width determined by the peak–valley distance (`delta = max.index - min.index`).  
   - Fit a quadratic polynomial (`stats::lm(y ~ poly(x, 2, raw = TRUE))`) to the density within this window, approximating the shape of the peak.

6. **Threshold from quadratic intercept**  
   - Solve the quadratic for its x‑intercepts, and take the **left intercept** as the threshold separating low‑MALAT1 droplets from the main cell distribution.  
   - If the computed intercept is negative, enforce `abs.min` (e.g. 1) as a lower bound to avoid pathological thresholds near zero.

8. **Error handling**  
   - If the computation fails (e.g. no meaningful peak because the sample is dominated by low‑quality droplets, or input not properly normalized), a `tryCatch` block:  
     - Prints an informative message suggesting to inspect the histogram and ensure normalized counts.  
     - Plots a raw histogram for visual diagnosis.  
     - Returns a conservative default threshold (e.g. 2) to avoid breaking downstream code.

This logic is a direct adaptation of BaderLab’s implementation, with only stylistic and documentation changes in this fork.

---

## Generalized feature threshold modules

Beyond MALAT1, the same framework can be used to derive thresholds for other features:

- Nuclear marker genes with high intronic content.
- Mitochondrial fraction (after preprocessing).
- Protocol‑specific QC markers that separate “good” vs “bad” droplets.

## Direct interface with SeuratObjects.

The repository provides a convenience function `.ComputeFeatureThresholdSeurat()` that takes a SeuratObject and computes the threshold for a specified feature (e.g., MALAT1) directly from the normalized counts in the specified assay and layer. This allows users to integrate the thresholding step seamlessly into their single-cell analysis workflow.

The computed threshold value and test results (boolean value if the cell passes the threshold) are stored in the SeuratObject’s metadata for downstream filtering and analysis. Each cell (row) contains the threshold it has been tested against, and a boolean column indicating whether it passes the threshold.

In the case of spliited objects, if the normalized counts are stored in the standar `data.layer` layers of the assay or the user specifies the concrete data layers, it computes the threshold per individual layer. 

---

## Usage

Below is a sketch of how to use the MALAT1 threshold function in a typical single‑cell workflow. Adapt names to match your actual function signatures in this repository.

```r
# Example 1: calculating MALAT1 threshold from a vector of log counts

# malat1_norm: numeric vector of normalized MALAT1 expression per cell
# e.g. from a Seurat object:
# malat1_norm <- FetchData(seurat_obj, vars = "MALAT1")[, "MALAT1"]
threshold.res <- .ComputeFeatureThreshold(
  counts.vector = malat1_norm,
  bw.bandwidth = 0.01,
  chosen.min = 2,
  smooth.spar = 1,
  abs.min  = 1,
  rough.max = 6,
  conservative.threshold = 2
)

# Example 2: applying the MALAT1 threshold to a Seurat object
# seurat_obj: a Seurat object with normalized counts in the "RNA" assay
seurat_obj <- .ComputeFeatureThresholdSeurat(
  SeuObj = seurat_obj,
  assay = "RNA",
  layers = NULL,  # or specify custom layers if not in the default "data" layer
  feature = "MALAT1", # Can be any feature name present in the assay
  # Parameters for the thresholding function (send as aditional arguments `...`)
  bw.bandwidth = 0.01,
  chosen.min = 2,
  smooth.spar = 1,
  abs.min = 1,
  rough.max = 6,
  conservative.threshold = 2
)

# Example 3: Compute multiple QC metrics and the MALAT1 threshold in one go with
# the CalculateQC() wrapper function (from SCutils)
seurat_obj <- SCutils::CalculateQC(
  seurat_obj,
  assay = "RNA",
  layers = NULL,
  perform.cell.cycle.scoring = TRUE,
  perform.MALAT1.test = TRUE
)
# Visualize the MALAT1 distribution and other multiple QC metrics in one single command:
plot.list <- FeatureDensityPlot(
   SeuratObject = seurat_obj,
   features = c("percent.mt", "MALAT1", "nCount_RNA"),
   layers = c(NA, "data", NA),
   vline = c(10, threshold.res$threshold, "both"),
   group.by = "batch",
   plot.title = "QC metrics distribution per batch",
 )
 invisible(lapply(plot.list, print))  # print all plots in the list
```

Typical recommendations:

- Run the threshold per **sample/batch**, on unintegrated normalized counts.
- Always inspect the MALAT1 distribution and applied threshold plots to ensure the automatic threshold aligns with the visual distribution. You can find useful plotting function for Seurat objects in the [SCutils](https://github.com/jmgs7/SCutils) package.
- Combine MALAT1 thresholds with other QC metrics (UMI counts, mitochondrial fraction, complexity...) rather than relying on a single metric.

---

## Installation

This repository is intended to be vendored into the [SCutils](https://github.com/jmgs7/SCutils) package.
For standalone use, clone the repository and source the relevant `.R` files in your analysis scripts:

```r
source("R/ComputeFeatureThresholdSeurat.R")   # adjust path to your installation
.ComputeFeatureThresholdSeurat(...)  # call the function with your Seurat object
```

---

## Licensing, attribution and references.

### Upstream work (BaderLab)

The original MALAT1 threshold code and concept were developed by BaderLab and released under the MIT License in the [`BaderLab/MALAT1_threshold`](https://github.com/BaderLab/MALAT1_threshold) repository.

This repository contains **derivative work** based on that MIT‑licensed upstream code. The original MIT license and copyright notice for BaderLab’s code must be preserved in the distribution of the derived portions.

### This repository (derivative work)

The adaptations in this repository include:

- Refactoring of the original MALAT1 threshold function to a more idiomatic R style (CamelCase function names, dot.case arguments and variables, explicit `package::function` calls).
- Extended, line‑by‑line comments and roxygen documentation.
- Direct interface with SeuratObjects to compute thresholds and store results in metadata.
- Generalization to threshold other features using the same density/spline/quadratic strategy.

> **License for this repository**  
> This derivative work is distributed under the **GNU General Public License, version 3 (GPL‑3.0)**.  
> See the file `LICENSE` in this repository for the full GPL‑3.0 text.
> Because the upstream MIT license is GPL‑compatible, it is legally permissible to incorporate MIT‑licensed code into a GPLv3‑licensed project, provided the MIT notice is retained for the parts derived from that upstream code.
> The original MIT license text for BaderLab’s MALAT1_threshold is included verbatim in a separate file (`LICENSE-BaderLab.txt`) as part of attribution.

### Summary of licensing

- Upstream BaderLab MALAT1_threshold repository: **MIT License**.
- This repository (derivative work, extended modules): **GPL‑3.0**, with MIT attribution retained for upstream‑derived portions.

GPLv3 © José Manuel Gómez Silva

### References

Clarke ZA, Bader GD. *MALAT1 expression indicates cell quality in single‑cell RNA sequencing data.* [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.14.603469v2) (2024). 

Please cite the original paper if you use this method in scientific publications.

---

## Contact and contributions

If you find bugs, have suggestions for improvements, please open an issue or pull request on GitHub.

When contributing code:

- Follow the existing style (based partially on Google's R style guide: CamelCase for exported functions, dot.case for arguments/variables, explicit `package::function` calls, etc).
- Include clear roxygen documentation and tests where appropriate.
- Respect licensing and attribution, particularly for any further upstream code you may incorporate.


