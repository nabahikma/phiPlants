# === Required Libraries ===
suppressMessages({
  library(xcms)
  library(CAMERA)
  library(MSnbase)
  library(tibble)
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
})

# === USER ARGUMENTS ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript save_final_plots_and_report.R <mzml_folder> <feature_ids_file> <output_excel_name>")
}

mzml_folder       <- args[1]
feature_ids_file  <- args[2]   # we pass finalReport.csv here (but TXT also supported)
output_excel_name <- args[3]

# ---- parameters you can tweak ----
GAUSS_BW_SEC     <- 5     # Gaussian bandwidth for EIC smoothing (seconds)
MZ_WINDOW_MS1    <- 0.5   # ± Da around the feature m/z for MS1 spectrum view
MZ_TOL_MS2       <- 0.01  # Da tolerance to consider MS2 matching the precursor
LEGEND_POS       <- "topright"
EIC_PADDING_SEC  <- 30    # final padding already applied in rtmin-30/rtmax+30
# ----------------------------------

# === Load MS data ===
mzml_files <- list.files(mzml_folder, pattern = "[.]mzML$", full.names = TRUE)
if (length(mzml_files) == 0) stop("No .mzML files found in folder: ", mzml_folder)
raw_data <- readMSData(mzml_files, mode = "onDisk")

# Optional xdata.rds if needed later
xdata_rds <- file.path(mzml_folder, "xdata.rds")
if (file.exists(xdata_rds)) {
  filled_data <- readRDS(xdata_rds)
}

# Annotated table (for mz & annotations)
ann_csvs <- list.files(mzml_folder, pattern = "_Annotated[.]csv$", full.names = TRUE)
if (length(ann_csvs) == 0) stop("Annotated CSV not found in folder: ", mzml_folder)
annotated_peaks <- read.csv(ann_csvs[1], check.names = TRUE, stringsAsFactors = FALSE)

# Features matrix (to fetch rtmin/rtmax in seconds)
features_matrix <- NULL
fm_path <- file.path(mzml_folder, "Features_Matrix.csv")
if (file.exists(fm_path)) {
  features_matrix <- read.csv(fm_path, check.names = TRUE, stringsAsFactors = FALSE)
} else {
  cand <- list.files(mzml_folder, pattern = "^Features_Matrix[.]csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(cand) > 0) features_matrix <- read.csv(cand[1], check.names = TRUE, stringsAsFactors = FALSE)
}

# Selected features: accept CSV with feature_id OR a plain TXT list
selected_feature_ids <- NULL
selected_df <- tryCatch(
  read.csv(feature_ids_file, check.names = TRUE, stringsAsFactors = FALSE),
  error = function(e) NULL
)
if (!is.null(selected_df) && "feature_id" %in% names(selected_df)) {
  selected_feature_ids <- unique(as.character(selected_df$feature_id))
} else {
  if (file.exists(feature_ids_file)) {
    ids_txt <- tryCatch(readLines(feature_ids_file, warn = FALSE), error = function(e) character(0))
    ids_txt <- trimws(ids_txt)
    ids_txt <- ids_txt[nchar(ids_txt) > 0]
    if (length(ids_txt) > 0) selected_feature_ids <- unique(as.character(ids_txt))
  }
}
if (is.null(selected_feature_ids) || length(selected_feature_ids) == 0) {
  stop("feature_ids_file does not provide any feature_id values: ", feature_ids_file)
}

# === Helpers ===
first_nonempty <- function(x) {
  x <- x[!is.na(x) & nchar(trimws(x)) > 0]
  if (length(x) == 0) return(NA_character_)
  x[1]
}
pick_annotation_name <- function(row) {
  cols <- names(row)
  val <- NA
  if ("SelectedAnnotation" %in% cols) val <- row[["SelectedAnnotation"]]
  if (is.na(val) && "annotation_candidates" %in% cols) val <- row[["annotation_candidates"]]
  if (is.na(val) && "Candidates" %in% cols) val <- row[["Candidates"]]
  if (is.na(val)) return(NA_character_)
  parts <- unlist(strsplit(as.character(val), "/"))
  first_nonempty(parts)
}
derive_rt_center_sec <- function(row) {
  cols <- names(row)
  if ("rt" %in% cols) {
    v <- suppressWarnings(as.numeric(row[["rt"]]))
    if (!is.na(v)) return(v)                # assume seconds if named 'rt'
  }
  if ("RT" %in% cols) {
    v <- suppressWarnings(as.numeric(row[["RT"]]))
    if (!is.na(v)) return(v * 60)           # minutes -> seconds
  }
  if ("RT.in.min" %in% cols) {
    v <- suppressWarnings(as.numeric(row[["RT.in.min"]]))
    if (!is.na(v)) return(v * 60)           # minutes -> seconds
  }
  alt <- grep("^RT[._]in[._]min$", cols, value = TRUE)
  if (length(alt) > 0) {
    v <- suppressWarnings(as.numeric(row[[alt[1]]]))
    if (!is.na(v)) return(v * 60)
  }
  NA_real_
}
gauss_smooth <- function(x, y, bw) {
  # Gaussian smoothing using stats::ksmooth, robust to NA and non-monotonic x
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 5 || bw <= 0) return(list(x = x, y = y))
  o <- order(x); x <- x[o]; y <- y[o]
  tryCatch({
    s <- ksmooth(x, y, kernel = "normal", bandwidth = bw)
    list(x = s$x, y = s$y)
  }, error = function(e) list(x = x, y = y))
}

# labels & colors per sample (consistent across plots)
sample_labels <- sub("[.]mzML$", "", basename(mzml_files))
sample_cols   <- setNames(rainbow(length(sample_labels)), sample_labels)

final_plots_folder <- file.path(mzml_folder, "Final_Plots")
if (!dir.exists(final_plots_folder)) dir.create(final_plots_folder, recursive = TRUE)

final_table <- data.frame()

# === Main loop per feature ===
for (feature_id_to_plot in selected_feature_ids) {

  feat_row <- annotated_peaks[annotated_peaks$feature_id == feature_id_to_plot, , drop = FALSE]
  if (nrow(feat_row) == 0) next

  mz <- suppressWarnings(as.numeric(feat_row$mz[1]))
  if (is.na(mz)) next

  # RT window: [rtmin-30, rtmax+30] in seconds, from Features_Matrix.csv
  rtrange <- NULL
  center_rt <- NA_real_
  if (!is.null(features_matrix) && all(c("feature_id","rtmin","rtmax") %in% names(features_matrix))) {
    fm_row <- features_matrix[features_matrix$feature_id == feature_id_to_plot, , drop = TRUE]
    if (is.data.frame(fm_row) && nrow(fm_row) >= 1) {
      rtmin <- suppressWarnings(as.numeric(fm_row$rtmin[1]))
      rtmax <- suppressWarnings(as.numeric(fm_row$rtmax[1]))
      if (!is.na(rtmin) && !is.na(rtmax)) {
        rtrange  <- c(max(0, rtmin - EIC_PADDING_SEC), rtmax + EIC_PADDING_SEC)
        center_rt <- (rtmin + rtmax) / 2
      }
    }
  }
  # Fallback if no rtmin/rtmax
  if (is.null(rtrange)) {
    center_rt <- derive_rt_center_sec(feat_row[1, ])
    if (is.na(center_rt)) next
    rtrange <- c(max(0, center_rt - 60), center_rt + 60)
  }
  if (is.na(center_rt)) center_rt <- mean(rtrange)

  # Annotation name for title
  ann_name <- pick_annotation_name(feat_row[1, ])
  ann_label <- if (is.na(ann_name) || ann_name == "" || ann_name == "No Annotation") "" else paste0("  |  ", ann_name)

  # Safe filename
  fname <- if (ann_label == "") paste0(feature_id_to_plot, ".png") else {
    paste0(gsub("[^A-Za-z0-9_\\-]", "_", ann_name), ".png")
  }
  full_plot_path <- file.path(final_plots_folder, fname)

  # === Open device and lay out 3 clean panels ===
  png(full_plot_path, width = 1200, height = 1800)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(3, 1), mar = c(4, 5, 4, 2) + 0.1)

  # ===== Panel 1: EIC (Gaussian smoothed), legend =====
  mzrange <- c(mz - 0.02, mz + 0.02)
  eic_obj <- chromatogram(raw_data, mz = mzrange, rt = rtrange)

  # Precompute smoothed lines and y-limits
  smoothed <- vector("list", length(sample_labels))
  ymax <- 1
  if (length(eic_obj) > 0) {
    for (i in seq_along(sample_labels)) {
      ch <- eic_obj[[i]]
      rt_vals  <- rtime(ch)
      int_vals <- intensity(ch)
      sm <- gauss_smooth(rt_vals, int_vals, GAUSS_BW_SEC)
      smoothed[[i]] <- sm
      if (length(sm$y) > 0 && any(is.finite(sm$y))) {
        ymax <- max(ymax, max(sm$y, na.rm = TRUE))
      }
    }
  }
  plot(NULL, xlim = rtrange, ylim = c(0, ymax),
       xlab = "Retention Time (s)", ylab = "Intensity",
       main = paste0("EIC (Gaussian bw=", GAUSS_BW_SEC, "s): feature ", feature_id_to_plot,
                     "   m/z=", round(mz, 4), ann_label))
  if (length(eic_obj) > 0) {
    for (i in seq_along(sample_labels)) {
      sm <- smoothed[[i]]
      if (is.null(sm) || length(sm$x) == 0) next
      lines(sm$x, sm$y, col = sample_cols[[ sample_labels[i] ]], lwd = 2)
    }
    # legend
    legend(LEGEND_POS, legend = sample_labels,
           col = sample_cols[sample_labels], lwd = 2, bty = "n", cex = 0.9)
  }

  # ===== Panel 2: MS1 near center RT =====
  rt_all <- rtime(raw_data)
  idx <- suppressWarnings(which.min(abs(rt_all - center_rt)))
  if (length(idx) == 1 && is.finite(idx) && idx > 0) {
    ms1_spec <- spectra(raw_data)[[idx]]
    mz_vals  <- mz(ms1_spec)
    int_vals <- intensity(ms1_spec)
    filt <- mz_vals >= (mz - MZ_WINDOW_MS1) & mz_vals <= (mz + MZ_WINDOW_MS1)
    ms1_mz  <- mz_vals[filt]
    ms1_int <- int_vals[filt]
    cols <- ifelse(abs(ms1_mz - mz) < 0.01, "red", "black")
    plot(ms1_mz, ms1_int, type = "h", col = cols, lwd = 2,
         main = paste0("MS1 Spectrum (±", MZ_WINDOW_MS1, " Da) at ~", round(center_rt), " s"),
         xlab = "m/z", ylab = "Intensity")
    abline(v = mz, col = "red", lty = 2)
  } else {
    plot.new(); text(0.5, 0.5, "No MS1 at this RT", cex = 1.3)
  }

  # ===== Panel 3: MS2 closest to center RT (if precursor matches) =====
  ms2_data <- filterMsLevel(raw_data, 2L)
  if (length(ms2_data) > 0) {
    ms2_rt <- rtime(ms2_data)
    j <- suppressWarnings(which.min(abs(ms2_rt - center_rt)))
    if (length(j) == 1 && is.finite(j) && j > 0) {
      spec2 <- spectra(ms2_data)[[j]]
      if (!is.na(precursorMz(spec2)) && abs(precursorMz(spec2) - mz) <= MZ_TOL_MS2) {
        plot(mz(spec2), intensity(spec2), type = "h", col = "blue", lwd = 2,
             main = paste0("MS2 Spectrum (precursor ~", round(precursorMz(spec2), 4), ")"),
             xlab = "m/z", ylab = "Intensity")
        abline(v = mz, col = "red", lty = 2)
      } else {
        plot.new(); text(0.5, 0.5, "No matching MS2 spectrum", cex = 1.3)
      }
    } else {
      plot.new(); text(0.5, 0.5, "No MS2 near RT", cex = 1.3)
    }
  } else {
    plot.new(); text(0.5, 0.5, "No MS2 data available", cex = 1.3)
  }

  # accumulate rows for Excel
  final_table <- bind_rows(final_table, feat_row)
}

# === Save Final Excel Report ===
openxlsx::write.xlsx(final_table, file = file.path(mzml_folder, output_excel_name))
cat("\n✅ Done: Final plots (Gaussian-smoothed EIC) and Excel report saved.\n")
