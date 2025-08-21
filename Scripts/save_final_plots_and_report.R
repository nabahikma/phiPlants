# === Required Libraries ===
suppressMessages({
  # Portable-friendly library paths (from env)
  lib_envs <- unique(Filter(nzchar, c(Sys.getenv("R_LIBS"),
                                      Sys.getenv("R_LIBS_USER"),
                                      Sys.getenv("R_USER_LIBS"))))
  if (length(lib_envs)) .libPaths(c(lib_envs, .libPaths()))

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
feature_ids_file  <- args[2]   # preferred: finalReport.csv (has feature_id); TXT list also accepted
output_excel_name <- args[3]

# ---- KNOBS ----
GAUSS_BW_SEC  <- 5
MZ_WINDOW_MS1 <- 0.5
MZ_TOL_MS2    <- 0.01
EIC_PAD_S     <- 30
LEGEND_POS    <- "topright"

PNG_WIDTH_IN  <- 8
PNG_HEIGHT_IN <- 12
PNG_DPI       <- 600
USE_RAGG      <- TRUE

Y_PAD_FACTOR  <- 1.30
# --------------

# === Load MS data ===
mzml_files <- list.files(mzml_folder, pattern = "[.]mzML$", full.names = TRUE)
if (length(mzml_files) == 0) stop("No .mzML files found in folder: ", mzml_folder)
raw_data <- readMSData(mzml_files, mode = "onDisk")

# Optional xdata.rds
xdata_rds <- file.path(mzml_folder, "xdata.rds")
if (file.exists(xdata_rds)) {
  filled_data <- readRDS(xdata_rds)
}

# Annotated table (for mz & annotations)
ann_csvs <- list.files(mzml_folder, pattern = "_Annotated[.]csv$", full.names = TRUE)
if (length(ann_csvs) == 0) stop("Annotated CSV not found in folder: ", mzml_folder)
annotated_peaks <- read.csv(ann_csvs[1], check.names = TRUE, stringsAsFactors = FALSE)

# Features_Matrix for rtmin/rtmax (seconds)
features_matrix <- NULL
fm_path <- file.path(mzml_folder, "Features_Matrix.csv")
if (file.exists(fm_path)) {
  features_matrix <- read.csv(fm_path, check.names = TRUE, stringsAsFactors = FALSE)
} else {
  cand <- list.files(mzml_folder, pattern = "^Features_Matrix[.]csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(cand) > 0) features_matrix <- read.csv(cand[1], check.names = TRUE, stringsAsFactors = FALSE)
}

# Selected features (CSV with feature_id OR TXT list)
selected_feature_ids <- NULL
selected_df <- tryCatch(
  read.csv(feature_ids_file, check.names = TRUE, stringsAsFactors = FALSE),
  error = function(e) NULL
)
if (!is.null(selected_df) && "feature_id" %in% names(selected_df)) {
  selected_feature_ids <- unique(as.character(selected_df$feature_id))
} else if (file.exists(feature_ids_file)) {
  ids_txt <- tryCatch(readLines(feature_ids_file, warn = FALSE), error = function(e) character(0))
  ids_txt <- trimws(ids_txt); ids_txt <- ids_txt[nchar(ids_txt) > 0]
  if (length(ids_txt) > 0) selected_feature_ids <- unique(as.character(ids_txt))
}
if (is.null(selected_feature_ids) || length(selected_feature_ids) == 0) {
  stop("feature_ids_file does not provide any feature_id values: ", feature_ids_file)
}

# === Helpers ===
safe_num <- function(x) suppressWarnings(as.numeric(x))

first_nonempty <- function(x) {
  x <- x[!is.na(x) & nchar(trimws(x)) > 0]
  if (length(x) == 0) return(NA_character_); x[1]
}
pick_annotation_name <- function(row) {
  cols <- names(row); val <- NA
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
    v <- safe_num(row[["rt"]]); if (!is.na(v)) return(v)
  }
  if ("RT" %in% cols) {
    v <- safe_num(row[["RT"]]); if (!is.na(v)) return(v * 60)
  }
  if ("RT.in.min" %in% cols) {
    v <- safe_num(row[["RT.in.min"]]); if (!is.na(v)) return(v * 60)
  }
  alt <- grep("^RT[._]in[._]min$", cols, value = TRUE)
  if (length(alt) > 0) {
    v <- safe_num(row[[alt[1]]]); if (!is.na(v)) return(v * 60)
  }
  NA_real_
}
gauss_smooth <- function(x, y, bw) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 5 || bw <= 0) return(list(x = x, y = y))
  o <- order(x); x <- x[o]; y <- y[o]
  out <- tryCatch(ksmooth(x, y, kernel = "normal", bandwidth = bw),
                  error = function(e) list(x = x, y = y))
  if (is.list(out) && !is.null(out$x) && !is.null(out$y)) list(x = out$x, y = out$y) else list(x = x, y = y)
}
pad_top <- function(ymax) {
  ymax <- ifelse(is.finite(ymax) && ymax > 0, ymax, 1)
  top  <- ymax * Y_PAD_FACTOR
  max(pretty(c(0, top)))
}
finite_sorted_window <- function(a, b) {
  a <- safe_num(a); b <- safe_num(b)
  if (!is.finite(a) || !is.finite(b)) return(NULL)
  if (a > b) { tmp <- a; a <- b; b <- tmp }
  c(a, b)
}

# Pull rt window for a feature_id (robust)
get_rt_window <- function(fid_chr) {
  # 1) From Features_Matrix.csv
  if (!is.null(features_matrix) && all(c("feature_id","rtmin","rtmax") %in% names(features_matrix))) {
    fm_idx <- which(as.character(features_matrix$feature_id) == fid_chr)
    if (length(fm_idx) >= 1) {
      rtmin <- safe_num(features_matrix$rtmin[fm_idx[1]])
      rtmax <- safe_num(features_matrix$rtmax[fm_idx[1]])
      if (is.finite(rtmin) && is.finite(rtmax)) {
        win <- finite_sorted_window(max(0, rtmin - EIC_PAD_S), rtmax + EIC_PAD_S)
        if (!is.null(win) && diff(win) > 0) return(list(rtrange = win, center = mean(c(rtmin, rtmax))))
      }
    }
  }
  # 2) From finalReport.csv (selected_df), if it carries rtmin/rtmax
  if (!is.null(selected_df) && all(c("rtmin","rtmax") %in% names(selected_df))) {
    sd_idx <- which(as.character(selected_df$feature_id) == fid_chr)
    if (length(sd_idx) >= 1) {
      rtmin <- safe_num(selected_df$rtmin[sd_idx[1]])
      rtmax <- safe_num(selected_df$rtmax[sd_idx[1]])
      if (is.finite(rtmin) && is.finite(rtmax)) {
        win <- finite_sorted_window(max(0, rtmin - EIC_PAD_S), rtmax + EIC_PAD_S)
        if (!is.null(win) && diff(win) > 0) return(list(rtrange = win, center = mean(c(rtmin, rtmax))))
      }
    }
  }
  # 3) From annotated_peaks row as last resort
  ap_idx <- which(as.character(annotated_peaks$feature_id) == fid_chr)
  if (length(ap_idx) >= 1) {
    row <- annotated_peaks[ap_idx[1], , drop = FALSE]
    rt_c <- derive_rt_center_sec(row)
    if (is.finite(rt_c)) {
      win <- finite_sorted_window(max(0, rt_c - 60), rt_c + 60)
      if (!is.null(win) && diff(win) > 0) return(list(rtrange = win, center = rt_c))
    }
  }
  NULL
}

# Colors per sample
sample_labels <- sub("[.]mzML$", "", basename(mzml_files))
sample_cols   <- setNames(rainbow(length(sample_labels)), sample_labels)

# Output
final_plots_folder <- file.path(mzml_folder, "Final_Plots")
if (!dir.exists(final_plots_folder)) dir.create(final_plots_folder, recursive = TRUE)

final_table <- data.frame()

# === Main loop ===
for (fid in selected_feature_ids) {
  fid_chr <- as.character(fid)

  ap_idx <- which(as.character(annotated_peaks$feature_id) == fid_chr)
  if (length(ap_idx) == 0) { message("Skip ", fid_chr, ": not found in annotated table"); next }
  feat_row <- annotated_peaks[ap_idx[1], , drop = FALSE]

  mz <- safe_num(feat_row$mz[1])
  if (!is.finite(mz)) { message("Skip ", fid_chr, ": mz not finite"); next }

  win <- get_rt_window(fid_chr)
  if (is.null(win)) { message("Skip ", fid_chr, ": could not determine finite RT window"); next }
  rtrange  <- win$rtrange
  center_rt <- win$center
  if (!is.finite(center_rt)) center_rt <- mean(rtrange)

  # Annotation for title + filename
  ann_name  <- pick_annotation_name(feat_row[1, ])
  ann_label <- if (is.na(ann_name) || ann_name == "" || ann_name == "No Annotation") "" else paste0("  |  ", ann_name)
  fname <- if (ann_label == "") paste0(fid_chr, ".png") else paste0(gsub("[^A-Za-z0-9_\\-]", "_", ann_name), ".png")
  full_plot_path <- file.path(final_plots_folder, fname)

  # --- open high-res device ---
  if (requireNamespace("ragg", quietly = TRUE) && USE_RAGG) {
    ragg::agg_png(full_plot_path, width = PNG_WIDTH_IN, height = PNG_HEIGHT_IN,
                  units = "in", res = PNG_DPI, scaling = 1)
  } else {
    png(full_plot_path, width = PNG_WIDTH_IN, height = PNG_HEIGHT_IN, units = "in",
        res = PNG_DPI, type = ifelse(.Platform$OS.type == "windows", "windows", "cairo"))
  }
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(3, 1), mar = c(4, 5, 4, 2) + 0.1)

  # ===== Panel 1: EIC (Gaussian-smoothed) =====
  mzrange <- c(mz - 0.02, mz + 0.02)
  eic_obj <- tryCatch(chromatogram(raw_data, mz = mzrange, rt = rtrange),
                      error = function(e) { message("EIC error for ", fid_chr, ": ", e$message); NULL })
  smoothed <- list(); ymax <- 1

  if (!is.null(eic_obj) && length(eic_obj) > 0) {
    n_ch <- length(eic_obj)
    for (i in seq_len(n_ch)) {
      ch <- eic_obj[[i]]
      rt_vals  <- rtime(ch); int_vals <- intensity(ch)
      sm <- gauss_smooth(rt_vals, int_vals, GAUSS_BW_SEC)
      smoothed[[i]] <- sm
      if (length(sm$y) > 0 && any(is.finite(sm$y))) ymax <- max(ymax, max(sm$y, na.rm = TRUE))
    }
  }
  ylim_eic <- c(0, pad_top(ymax))
  plot(NULL, xlim = rtrange, ylim = ylim_eic,
       xlab = "Retention Time (s)", ylab = "Intensity",
       main = paste0("EIC (Gaussian bw=", GAUSS_BW_SEC, "s): feature ", fid_chr,
                     "   m/z=", round(mz, 4), ann_label))
  if (!is.null(eic_obj) && length(eic_obj) > 0) {
    n_ch <- length(eic_obj)
    lbls <- sample_labels[seq_len(min(length(sample_labels), n_ch))]
    for (i in seq_len(n_ch)) {
      sm <- smoothed[[i]]
      if (is.null(sm) || length(sm$x) == 0) next
      col_i <- if (i <= length(lbls)) sample_cols[[ lbls[i] ]] else "gray40"
      lines(sm$x, sm$y, col = col_i, lwd = 2)
    }
    legend(LEGEND_POS, legend = lbls,
           col = sample_cols[lbls], lwd = 2, bty = "n", cex = 0.9)
  }

  # ===== Panel 2: MS1 =====
  rt_all <- rtime(raw_data)
  idx <- suppressWarnings(which.min(abs(rt_all - center_rt)))
  if (length(idx) == 1 && is.finite(idx) && idx > 0) {
    ms1_spec <- spectra(raw_data)[[idx]]
    mz_vals  <- mz(ms1_spec); int_vals <- intensity(ms1_spec)
    filt <- mz_vals >= (mz - MZ_WINDOW_MS1) & mz_vals <= (mz + MZ_WINDOW_MS1)
    ms1_mz  <- mz_vals[filt]; ms1_int <- int_vals[filt]
    cols <- ifelse(abs(ms1_mz - mz) < 0.01, "red", "black")
    ylim_ms1 <- c(0, pad_top(if (length(ms1_int)) max(ms1_int, na.rm = TRUE) else 1))
    plot(ms1_mz, ms1_int, type = "h", col = cols, lwd = 2,
         main = paste0("MS1 Spectrum (±", MZ_WINDOW_MS1, " Da) at ~", round(center_rt), " s"),
         xlab = "m/z", ylab = "Intensity", ylim = ylim_ms1)
    abline(v = mz, col = "red", lty = 2)
  } else {
    plot.new(); text(0.5, 0.5, "No MS1 at this RT", cex = 1.3)
  }

  # ===== Panel 3: MS2 =====
  ms2_data <- filterMsLevel(raw_data, 2L)
  if (length(ms2_data) > 0) {
    ms2_rt <- rtime(ms2_data)
    j <- suppressWarnings(which.min(abs(ms2_rt - center_rt)))
    if (length(j) == 1 && is.finite(j) && j > 0) {
      spec2 <- spectra(ms2_data)[[j]]
      if (!is.na(precursorMz(spec2)) && abs(precursorMz(spec2) - mz) <= MZ_TOL_MS2) {
        ms2_mz <- mz(spec2); ms2_int <- intensity(spec2)
        ylim_ms2 <- c(0, pad_top(if (length(ms2_int)) max(ms2_int, na.rm = TRUE) else 1))
        plot(ms2_mz, ms2_int, type = "h", col = "blue", lwd = 2,
             main = paste0("MS2 Spectrum (precursor ~", round(precursorMz(spec2), 4), ")"),
             xlab = "m/z", ylab = "Intensity", ylim = ylim_ms2)
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

  # accumulate for Excel
  final_table <- bind_rows(final_table, feat_row)
}

# === Save Final Excel Report ===
openxlsx::write.xlsx(final_table, file = file.path(mzml_folder, output_excel_name))
cat("\n✅ Done: Final plots (Gaussian-smoothed EIC) and Excel report saved.\n")
