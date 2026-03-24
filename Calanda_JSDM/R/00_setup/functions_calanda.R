# ==============================================================================
# functions_calanda.R
# Author: Billur Bektas
# Claude (Anthropic) was used to assist with code refactoring, validation, and documentation.
#
# Helper functions for the Calanda JSDM project
#
# Sections:
#   1. DATA PROCESSING: Snow, temperature, mowing, imputation
#   2. PLOTTING: Venn diagram, labeling helpers
#   3. ANALYSIS PLOTTING FUNCTIONS: VP scatter, partial effects, paired plots
# ==============================================================================


# ==============================================================================
# 1. DATA PROCESSING FUNCTIONS
# ==============================================================================

#' Process Copernicus GFSC snow cover data from ZIP files
#'
#' Extracts gap-filled snow cover (GF) raster values at vegetation plot
#' locations from Copernicus GFSC ZIP archives.
#'
#' @param copernicus_dir Directory containing GFSC ZIP files
#' @param output_dir Directory for output files
#' @param veg_coords_path Path to CSV with plot coordinates (plot_id, releve_id, x, y)
#' @param calanda_mask_path Path to Calanda study area shapefile
#' @param force_reprocess If TRUE, reprocess even if output exists
#' @return Path to output CSV file
process_gfsc_data = function(copernicus_dir,
                             output_dir,
                             veg_coords_path,
                             calanda_mask_path,
                             force_reprocess = FALSE) {

  suppressPackageStartupMessages({
    library(tidyverse)
    library(terra)
    library(lubridate)
    library(sf)
  })

  output_file = file.path(output_dir, "gf_vegetation_data.csv")

  if (file.exists(output_file) && !force_reprocess) {
    message("Output file already exists. Set force_reprocess=TRUE to reprocess.")
    return(output_file)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  zip_files = list.files(copernicus_dir, pattern = "GFSC_\\d{8}-007_S1-S2_T32T\\w{2}_V101_\\d+\\.zip$", full.names = TRUE)

  if (length(zip_files) == 0) {
    stop("No matching ZIP files found in ", copernicus_dir)
  }

  if (!file.exists(veg_coords_path)) {
    stop("Vegetation coordinates file not found: ", veg_coords_path)
  }

  if (!file.exists(calanda_mask_path)) {
    if (!grepl("\\.shp$", calanda_mask_path) && !file.exists(paste0(calanda_mask_path, ".shp"))) {
      stop("Calanda mask shapefile not found: ", calanda_mask_path)
    }
  }

  calanda_mask_sf = st_read(calanda_mask_path, quiet = TRUE)

  extract_date = function(filename) {
    date_str = str_extract(basename(filename), "\\d{8}")
    return(ymd(date_str))
  }

  veg_coords = read_csv(veg_coords_path)

  veg_points_sf = st_as_sf(veg_coords,
                            coords = c("x", "y"),
                            crs = 4326)

  safe_create_dir = function(dir_path) {
    if (!dir.exists(dir_path)) {
      disk_info = system2("df", args = c("-k", dirname(dir_path)), stdout = TRUE)
      avail_space = as.numeric(strsplit(disk_info[2], "\\s+")[[1]][4])

      if (avail_space < 100 * 1024) {
        message("Low disk space detected. Cleaning up extraction directory...")
        temp_dirs = list.dirs(output_dir, recursive = FALSE)
        for (dir in temp_dirs) {
          if (grepl("extracted", dir)) {
            unlink(dir, recursive = TRUE)
          }
        }
      }

      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(dir_path)) {
        stop("Failed to create directory: ", dir_path, ". Check disk space and permissions.")
      }
    }
    return(dir_path)
  }

  temp_extract_dir = file.path(tempdir(), "gfsc_extract")
  safe_create_dir(temp_extract_dir)

  results = map_df(zip_files, function(zip_file) {
    tryCatch({
      date = extract_date(zip_file)

      if (dir.exists(temp_extract_dir)) {
        unlink(list.files(temp_extract_dir, full.names = TRUE), recursive = TRUE)
      }

      message("Processing ", basename(zip_file), " (", date, ")")

      message("  Extracting...")
      extracted_files = unzip(zip_file, exdir = temp_extract_dir, list = TRUE)
      unzip(zip_file, exdir = temp_extract_dir)

      gf_file = list.files(temp_extract_dir, pattern = "_GF\\.tif$", recursive = TRUE, full.names = TRUE)

      if (length(gf_file) == 0) {
        warning(paste("No GF.tif file found in", zip_file))
        return(tibble())
      }

      message("  Found GF file: ", basename(gf_file))

      gf_raster = rast(gf_file)

      calanda_crs = st_crs(calanda_mask_sf)

      message("  Reprojecting raster...")
      gf_raster_reproj = project(gf_raster, calanda_crs$wkt)

      if (st_crs(veg_points_sf) != calanda_crs) {
        veg_points_sf = st_transform(veg_points_sf, calanda_crs)
      }

      message("  Extracting values at vegetation points...")
      extracted_values = terra::extract(gf_raster_reproj, veg_points_sf)

      extracted_df = as.data.frame(extracted_values)

      result = veg_coords %>%
        bind_cols(raster_value = extracted_df[[2]]) %>%
        mutate(date = date)

      message("  Cleaning up extracted files...")
      unlink(list.files(temp_extract_dir, full.names = TRUE), recursive = TRUE)

      return(result)
    },
    error = function(e) {
      warning(paste("Error processing", zip_file, ":", e$message))
      unlink(list.files(temp_extract_dir, full.names = TRUE), recursive = TRUE)
      return(tibble())
    })
  })

  message("Final cleanup...")
  unlink(temp_extract_dir, recursive = TRUE)

  if (nrow(results) == 0) {
    warning("No data was extracted from any of the files")
    final_results = veg_coords[0, ] %>%
      mutate(date = as.Date(character()),
             raster_value = numeric())
  } else {
    final_results = results %>%
      select(date, plot_id, releve_id, x, y, raster_value) %>%
      arrange(date, plot_id)
  }

  write_csv(final_results, output_file)

  message("Processed ", n_distinct(final_results$date), " dates with ", nrow(final_results), " total observations")
  message("Data saved to ", output_file)

  if (nrow(final_results) > 0) {
    p = final_results %>%
      ggplot(aes(x = date, y = raster_value, color = factor(plot_id))) +
      geom_line() +
      geom_point() +
      labs(title = "Raster Values Over Time by Plot ID",
           x = "Date",
           y = "Raster Value",
           color = "Plot ID") +
      theme_bw()

    plot_path = file.path(output_dir, "temporal_plot.png")
    ggsave(plot_path, p, width = 10, height = 6)
    message("Plot saved to ", plot_path)
  } else {
    message("No data to plot")
  }

  return(output_file)
}


#' Preprocess snow cover data with BISE correction and interpolation
#'
#' Applies BISE correction for sudden drops/spikes, spline interpolation for
#' missing values, Savitzky-Golay smoothing, and daily interpolation.
#'
#' @param raw_snow_df Data frame with columns: date, raster_value, plot_id, releve_id
#' @param bise_threshold Threshold for BISE correction (default 0.2)
#' @param sliding_window Window size for sliding max (default 3)
#' @param sgfilter_p Savitzky-Golay polynomial order (default 3)
#' @param sgfilter_n Savitzky-Golay window size (default 7)
#' @param expand_time_series Whether to expand to a global date range (default TRUE)
#' @param date_range Optional date range for expansion
#' @return Data frame with daily snow cover values per plot-releve combination
process_snow_data = function(raw_snow_df,
                             bise_threshold = 0.2,
                             sliding_window = 3,
                             sgfilter_p = 3,
                             sgfilter_n = 7,
                             expand_time_series = TRUE,
                             date_range = NULL) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(zoo)
    library(signal)
    library(lubridate)
  })

  required_cols = c("date", "raster_value", "plot_id", "releve_id")
  missing_cols = setdiff(required_cols, colnames(raw_snow_df))
  if (length(missing_cols) > 0) {
    stop("Input dataframe missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  raw_snow_df = raw_snow_df %>%
    mutate(plot_id_releve = paste(plot_id, releve_id, sep = "_")) %>%
    arrange(plot_id_releve, date)

  if (expand_time_series) {
    if (!is.null(date_range)) {
      global_date_range = date_range
      message("Using user-provided date range: ",
              format(global_date_range[1], "%Y-%m-%d"), " to ",
              format(global_date_range[2], "%Y-%m-%d"))
    } else {
      global_date_range = range(raw_snow_df$date, na.rm = TRUE)
      message("Using global date range from data: ",
              format(global_date_range[1], "%Y-%m-%d"), " to ",
              format(global_date_range[2], "%Y-%m-%d"))
    }
  } else {
    global_date_range = NULL
    message("Time series will not be expanded (using location-specific date ranges)")
  }

  plot_releve_combos = unique(raw_snow_df$plot_id_releve)
  message("Processing ", length(plot_releve_combos), " unique plot-releve combinations...")

  all_results = list()

  for (curr_combo in plot_releve_combos) {
    curr_data = raw_snow_df %>%
      filter(plot_id_releve == curr_combo)

    if (nrow(curr_data) < 3 || sum(!is.na(curr_data$raster_value)) < 3) {
      message("Skipping ", curr_combo, " - insufficient data points")
      next
    }

    curr_data = curr_data %>% arrange(date)

    message("Processing ", curr_combo, " with ", nrow(curr_data), " observations")

    snow_raw = curr_data$raster_value
    snow_processed = snow_raw

    # 1. BISE correction to remove sudden drops/spikes
    if (length(snow_raw) > 1) {
      diff_vals = diff(snow_raw)
      decrease_idx = which(diff_vals < 0) + 1
      data_threshold = c(NA, snow_raw[-1] - bise_threshold * diff_vals)

      if (length(snow_raw) >= sliding_window) {
        val_sliding_period = zoo::rollapply(snow_raw,
                                            width = sliding_window,
                                            FUN = max,
                                            fill = NA,
                                            align = "left")

        if (length(decrease_idx) > 0) {
          valid_idx = decrease_idx[decrease_idx <= length(val_sliding_period) &
                                     decrease_idx <= length(data_threshold)]

          if (length(valid_idx) > 0) {
            law_check = val_sliding_period[valid_idx] - data_threshold[valid_idx]
            reject_decrease = valid_idx[which(law_check > 0)]
            snow_processed[reject_decrease] = NA
          }
        }

        max_val = max(snow_raw, na.rm = TRUE)
        if (!is.infinite(max_val)) {
          increase_threshold = bise_threshold * max_val
          increase_idx = which(diff_vals > increase_threshold)

          if (length(increase_idx) > 0) {
            reject_increase = increase_idx[!increase_idx %in% decrease_idx]
            valid_increases = reject_increase[reject_increase < length(snow_processed)]

            if (length(valid_increases) > 0) {
              snow_processed[valid_increases + 1] = NA
            }
          }
        }
      }
    }

    if (sum(!is.na(snow_processed)) < 3) {
      message("  Skipping ", curr_combo, " - insufficient data after BISE correction")
      next
    }

    # 2. Fill start and end gaps
    first_valid = which(!is.na(snow_processed))[1]
    if (!is.na(first_valid) && first_valid > 1) {
      snow_processed[1:(first_valid - 1)] = snow_processed[first_valid]
    }

    last_valid = max(which(!is.na(snow_processed)))
    if (!is.na(last_valid) && last_valid < length(snow_processed)) {
      snow_processed[(last_valid + 1):length(snow_processed)] = snow_processed[last_valid]
    }

    # 3. Interpolate missing values using spline
    if (any(is.na(snow_processed))) {
      if (sum(!is.na(snow_processed)) >= 3) {
        snow_interpolated = tryCatch({
          zoo::na.spline(snow_processed)
        }, error = function(e) {
          message("  Error in spline interpolation for ", curr_combo, ": ", e$message)
          approx(x = which(!is.na(snow_processed)),
                 y = snow_processed[!is.na(snow_processed)],
                 xout = 1:length(snow_processed),
                 rule = 2)$y
        })
      } else {
        message("  Skipping ", curr_combo, " - insufficient data for interpolation")
        next
      }
    } else {
      snow_interpolated = snow_processed
    }

    # 4. Apply Savitzky-Golay filter for smoothing
    if (length(snow_interpolated) >= sgfilter_n) {
      snow_smoothed = tryCatch({
        signal::sgolayfilt(snow_interpolated, p = sgfilter_p, n = sgfilter_n, m = 0)
      }, error = function(e) {
        message("  Error in Savitzky-Golay filter for ", curr_combo, ": ", e$message)
        snow_interpolated
      })
    } else {
      snow_smoothed = snow_interpolated
    }

    # 5. Create daily interpolation with expanded time series
    loc_date_range = range(curr_data$date)

    if (!is.null(global_date_range)) {
      expanded_range = global_date_range
    } else {
      expanded_range = loc_date_range
    }

    dates_seq = seq(expanded_range[1], expanded_range[2], by = "days")

    needs_extrapolation = min(dates_seq) < min(curr_data$date) || max(dates_seq) > max(curr_data$date)

    if (needs_extrapolation) {
      message("  Extrapolating outside observed range for ", curr_combo)

      loc_date_range = range(curr_data$date)

      if (length(curr_data$date) < 4) {
        daily_values = approx(x = as.numeric(curr_data$date),
                              y = snow_smoothed,
                              xout = as.numeric(dates_seq),
                              rule = 2)$y
      } else {
        tryCatch({
          daily_values = spline(x = as.numeric(curr_data$date),
                                y = snow_smoothed,
                                xout = as.numeric(dates_seq),
                                method = "natural")$y
        }, error = function(e) {
          message("  Error in spline extrapolation for ", curr_combo, ": ", e$message)
          daily_values = approx(x = as.numeric(curr_data$date),
                                y = snow_smoothed,
                                xout = as.numeric(dates_seq),
                                rule = 2)$y
        })
      }
    } else {
      if (length(curr_data$date) < 4) {
        daily_values = approx(x = as.numeric(curr_data$date),
                              y = snow_smoothed,
                              xout = as.numeric(dates_seq),
                              rule = 2)$y
      } else {
        daily_values = tryCatch({
          spline(x = as.numeric(curr_data$date),
                 y = snow_smoothed,
                 xout = as.numeric(dates_seq))$y
        }, error = function(e) {
          message("  Error in daily spline interpolation for ", curr_combo, ": ", e$message)
          approx(x = as.numeric(curr_data$date),
                 y = snow_smoothed,
                 xout = as.numeric(dates_seq),
                 rule = 2)$y
        })
      }
    }

    # 6. Ensure values are within valid range (0-100 for percentages)
    daily_values = pmin(pmax(daily_values, 0), 100)

    observed_dates = curr_data$date

    result_df = data.frame(
      date = dates_seq,
      snow_raw = approx(x = as.numeric(curr_data$date),
                        y = curr_data$raster_value,
                        xout = as.numeric(dates_seq),
                        rule = 2)$y,
      snow_processed = approx(x = as.numeric(curr_data$date),
                              y = snow_processed,
                              xout = as.numeric(dates_seq),
                              rule = 2)$y,
      snow_interpolated = approx(x = as.numeric(curr_data$date),
                                 y = snow_interpolated,
                                 xout = as.numeric(dates_seq),
                                 rule = 2)$y,
      snow_smoothed = daily_values,
      plot_id = curr_data$plot_id[1],
      releve_id = curr_data$releve_id[1],
      plot_id_releve = curr_combo,
      is_observed = dates_seq %in% observed_dates
    )

    result_df = result_df %>%
      mutate(
        year = year(date),
        month = month(date),
        day = day(date),
        doy = yday(date)
      )

    all_results[[curr_combo]] = result_df
  }

  if (length(all_results) == 0) {
    stop("No valid data could be processed for any plot-releve combination")
  }

  final_results = bind_rows(all_results)

  message("Completed processing with ",
          n_distinct(final_results$plot_id_releve), " plot-releve combinations and ",
          nrow(final_results), " total daily records")

  return(final_results)
}


#' Calculate snow disappearance date from daily snow cover data
#'
#' Finds the day of year when snow cover drops below 10% for at least 30
#' consecutive days after the winter peak.
#'
#' @param snow Numeric vector of daily snow cover values
#' @param doy Numeric vector of day-of-year values
#' @param month Numeric vector of month values
#' @return Day of year of snow disappearance, or NA
calculate_snow_disappearance = function(snow, doy, month) {
  winter_idx = which(month %in% c(1, 2, 3))
  if (length(winter_idx) == 0) return(NA)

  peak_idx = winter_idx[which.max(snow[winter_idx])]
  if (length(peak_idx) == 0 || is.na(peak_idx)) return(NA)

  summer_cutoff = which(doy >= 200)[1]
  if (is.na(summer_cutoff)) summer_cutoff = length(snow)

  remaining_idx = (peak_idx + 1):min(summer_cutoff, length(snow))
  if (length(remaining_idx) == 0) return(NA)

  low_snow_days = which(snow[remaining_idx] <= 10)

  if (length(low_snow_days) == 0) return(NA)

  consecutive_days = 30
  streaks = rle(diff(low_snow_days) == 1)
  streak_idx = which(streaks$values & streaks$lengths >= (consecutive_days - 1))

  if (length(streak_idx) == 0) {
    snow_melt_idx = remaining_idx[min(low_snow_days)]
  } else {
    start_pos = ifelse(streak_idx[1] > 1,
                       sum(streaks$lengths[1:(streak_idx[1] - 1)]) + 1,
                       1)
    snow_melt_idx = remaining_idx[low_snow_days[start_pos]]
  }

  return(doy[snow_melt_idx])
}


#' Calculate snow metrics per plot-releve combination
#'
#' Computes snow median, sum, days with snow, and snow disappearance day.
#'
#' @param data Data frame from process_snow_data()
#' @param stats_by_year If TRUE, group by year as well (default TRUE)
#' @return Data frame of snow metrics
calculate_snow_metrics = function(data, stats_by_year = TRUE) {
  if (stats_by_year) {
    group_vars = c("plot_id", "releve_id", "plot_id_releve", "year")
  } else {
    group_vars = c("plot_id", "releve_id", "plot_id_releve")
  }

  snow_metrics = data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      snow_median = median(snow_smoothed, na.rm = TRUE),
      snow_sum = sum(snow_smoothed, na.rm = TRUE),
      days_with_snow = sum(snow_smoothed > 5, na.rm = TRUE),
      snow_disappearance_doy = calculate_snow_disappearance(snow_smoothed, doy, month),
      start_date = min(date, na.rm = TRUE),
      end_date = max(date, na.rm = TRUE),
      days_covered = as.numeric(difftime(max(date, na.rm = TRUE),
                                         min(date, na.rm = TRUE),
                                         units = "days")) + 1,
      n_observed = sum(is_observed),
      n_total = n(),
      .groups = "drop"
    )

  return(snow_metrics)
}


#' Process temperature data with interpolation and smoothing
#'
#' Processes raw ECOSTRESS LST data: daily aggregation, spline smoothing,
#' gap filling, and optional time series expansion.
#'
#' @param raw_temp_df Data frame with columns: Date, plot_id_releve, lst
#' @param spline_smoothing Smoothing parameter for smooth.spline (default 0.5)
#' @param expand_time_series Whether to expand to global date range (default TRUE)
#' @param date_range Optional date range for expansion
#' @param temp_col Name of temperature column (default "lst")
#' @return Data frame with daily interpolated temperatures
process_temperature_data = function(raw_temp_df,
                                    spline_smoothing = 0.5,
                                    expand_time_series = TRUE,
                                    date_range = NULL,
                                    temp_col = "lst") {
  suppressPackageStartupMessages({
    library(zoo)
    library(lubridate)
  })

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }

  required_cols = c("Date", "plot_id_releve", temp_col)
  missing_cols = setdiff(required_cols, colnames(raw_temp_df))
  if (length(missing_cols) > 0) {
    stop("Input dataframe missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (is.character(raw_temp_df$Date) || is.factor(raw_temp_df$Date)) {
    raw_temp_df$Date = as.POSIXct(raw_temp_df$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    if (all(is.na(raw_temp_df$Date))) {
      raw_temp_df$Date = as.POSIXct(raw_temp_df$Date, format = "%Y-%m-%d", tz = "UTC")
    }
  }

  raw_temp_df$date_only = as.Date(raw_temp_df$Date)

  raw_temp_df = raw_temp_df[order(raw_temp_df$plot_id_releve, raw_temp_df$Date), ]

  if (expand_time_series) {
    if (!is.null(date_range)) {
      global_date_range = date_range
      message("Using user-provided date range: ",
              format(global_date_range[1], "%Y-%m-%d"), " to ",
              format(global_date_range[2], "%Y-%m-%d"))
    } else {
      global_date_range = range(raw_temp_df$date_only, na.rm = TRUE)
      message("Using global date range from data: ",
              format(global_date_range[1], "%Y-%m-%d"), " to ",
              format(global_date_range[2], "%Y-%m-%d"))
    }
  } else {
    global_date_range = NULL
    message("Time series will not be expanded (using location-specific date ranges)")
  }

  plot_releve_combos = unique(raw_temp_df$plot_id_releve)
  message("Processing ", length(plot_releve_combos), " unique plot-releve combinations...")

  all_results = list()

  for (curr_combo in plot_releve_combos) {
    curr_data = raw_temp_df[raw_temp_df$plot_id_releve == curr_combo, ]

    if (nrow(curr_data) < 3 || sum(!is.na(curr_data[[temp_col]])) < 3) {
      message("Skipping ", curr_combo, " - insufficient data points")
      next
    }

    curr_data = curr_data[order(curr_data$Date), ]

    message("Processing ", curr_combo, " with ", nrow(curr_data), " observations")

    daily_data = NULL
    try({
      daily_data = tapply(curr_data[[temp_col]], curr_data$date_only, mean, na.rm = TRUE)
      daily_data = data.frame(
        date = as.Date(names(daily_data)),
        temp = as.numeric(daily_data)
      )
      daily_data = daily_data[!is.nan(daily_data$temp), ]
    }, silent = TRUE)

    if (is.null(daily_data) || nrow(daily_data) < 3) {
      message("  Using original time points (daily aggregation failed)")
      daily_data = data.frame(
        date = curr_data$date_only,
        temp = curr_data[[temp_col]]
      )
      daily_data = daily_data[!is.na(daily_data$temp), ]
    }

    if (nrow(daily_data) < 3) {
      message("  Skipping ", curr_combo, " - insufficient data after aggregation")
      next
    }

    loc_date_range = range(daily_data$date, na.rm = TRUE)

    if (!is.null(global_date_range)) {
      expanded_range = global_date_range
      needs_extrapolation = TRUE
    } else {
      expanded_range = loc_date_range
      needs_extrapolation = FALSE
    }

    dates_seq = seq(expanded_range[1], expanded_range[2], by = "days")

    interpolated_temps = NULL

    observed_dates = daily_data$date

    if (needs_extrapolation || length(dates_seq) > length(daily_data$date)) {
      message("  Filling gaps in time series for ", curr_combo)

      if (nrow(daily_data) >= 5) {
        try({
          interpolated_temps = smooth.spline(
            x = as.numeric(daily_data$date),
            y = daily_data$temp,
            spar = spline_smoothing
          )

          interpolated_temps = predict(
            interpolated_temps,
            x = as.numeric(dates_seq)
          )$y
        }, silent = TRUE)
      }

      if (is.null(interpolated_temps)) {
        message("  Using linear interpolation/extrapolation")
        interpolated_temps = approx(
          x = as.numeric(daily_data$date),
          y = daily_data$temp,
          xout = as.numeric(dates_seq),
          rule = 2
        )$y
      }
    } else {
      interpolated_temps = daily_data$temp
    }

    result_df = data.frame(
      date = dates_seq,
      temp_raw = approx(
        x = as.numeric(daily_data$date),
        y = daily_data$temp,
        xout = as.numeric(dates_seq),
        rule = 2
      )$y,
      temp_interpolated = interpolated_temps,
      plot_id_releve = curr_combo,
      is_observed = dates_seq %in% observed_dates
    )

    if ("Latitude" %in% colnames(curr_data) && "Longitude" %in% colnames(curr_data)) {
      result_df$Latitude = curr_data$Latitude[1]
      result_df$Longitude = curr_data$Longitude[1]
    }

    result_df$year = lubridate::year(result_df$date)
    result_df$month = lubridate::month(result_df$date)
    result_df$doy = lubridate::yday(result_df$date)

    all_results[[curr_combo]] = result_df
  }

  if (length(all_results) == 0) {
    stop("No valid data could be processed for any plot-releve combination")
  }

  final_results = do.call(rbind, all_results)

  message("Completed processing with ",
          length(all_results), " plot-releve combinations and ",
          nrow(final_results), " total daily records")

  return(final_results)
}


#' Find the first spring date when temperature reaches 1 degrees C
#'
#' Searches March-July for the first day temperature exceeds 1 degrees C.
#'
#' @param dates Date vector
#' @param temps Temperature vector
#' @param year Year to filter
#' @return Day of year, or NA
find_spring_warming_date = function(dates, temps, year) {
  idx = lubridate::month(dates) %in% 3:7 & lubridate::year(dates) == year
  if (sum(idx) == 0) return(NA)

  spring_dates = dates[idx]
  spring_temps = temps[idx]

  ord = order(spring_dates)
  spring_dates = spring_dates[ord]
  spring_temps = spring_temps[ord]

  idx_warm = which(spring_temps >= 1)

  if (length(idx_warm) == 0) {
    return(NA)
  } else {
    warming_date = spring_dates[min(idx_warm)]
    return(lubridate::yday(warming_date))
  }
}


#' Calculate temperature metrics: summer temperature and spring warming date
#'
#' @param data Data frame from process_temperature_data()
#' @return Data frame with summer_temp and spring_warming_doy per plot-releve-year
calculate_temp_metrics = function(data) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for the metrics calculation")
  }

  summer_temp = data %>%
    dplyr::filter(month %in% 6:8) %>%
    dplyr::group_by(plot_id_releve, year) %>%
    dplyr::summarise(
      summer_temp = mean(temp_interpolated, na.rm = TRUE),
      .groups = "drop"
    )

  warming_dates = data %>%
    dplyr::group_by(plot_id_releve, year) %>%
    dplyr::summarise(
      spring_warming_doy = find_spring_warming_date(date, temp_interpolated, first(year)),
      .groups = "drop"
    )

  metrics = dplyr::full_join(summer_temp, warming_dates, by = c("plot_id_releve", "year"))

  if (all(c("Latitude", "Longitude") %in% colnames(data))) {
    lat_lon = data %>%
      dplyr::group_by(plot_id_releve) %>%
      dplyr::summarise(
        Latitude = dplyr::first(Latitude),
        Longitude = dplyr::first(Longitude),
        .groups = "drop"
      )

    metrics = dplyr::left_join(metrics, lat_lon, by = "plot_id_releve")
  }

  return(metrics)
}


#' Impute missing environmental data using miceRanger
#'
#' Performs multiple imputation with random forests, including spatial features
#' (distance to centroid, elevation-distance ratio) and cross-validation.
#'
#' @param data Data frame with environmental variables and Longitude, Latitude, altitude
#' @param variables_to_impute Character vector of column names to impute
#' @param m Number of imputed datasets (default 5)
#' @param maxiter Number of iterations (default 10)
#' @param num_trees Trees per random forest (default 100)
#' @param seed Random seed (default 123)
#' @param validation_fraction Fraction of data to mask for validation (default 0.2)
#' @return List with imputed_data, performance, and mice_model
impute_environmental_data = function(
    data,
    variables_to_impute = c("summer_temp", "fdd", "et_annual", "snow_sum"),
    m = 5,
    maxiter = 10,
    num_trees = 100,
    seed = 123,
    validation_fraction = 0.2
) {

  library(tidyverse)
  library(miceRanger)
  library(sf)

  input_data = data

  data_sf = st_as_sf(input_data, coords = c("Longitude", "Latitude"), crs = 4326)

  data_sf_projected = st_transform(data_sf, 32632)

  centroid = st_centroid(st_combine(data_sf_projected))

  input_data$dist_to_centroid = st_distance(data_sf_projected, centroid) %>% as.numeric()

  input_data$elev_dist_ratio = input_data$altitude / (input_data$dist_to_centroid + 1)

  set.seed(seed)

  original_data = input_data

  validation_sets = list()
  masked_indices = list()

  for (var in variables_to_impute) {
    complete_indices = which(!is.na(input_data[[var]]))

    if (length(complete_indices) > 5) {
      mask_size = floor(length(complete_indices) * validation_fraction)
      mask_idx = sample(complete_indices, mask_size)

      validation_sets[[var]] = input_data[[var]][mask_idx]

      input_data[[var]][mask_idx] = NA

      masked_indices[[var]] = mask_idx
    } else {
      validation_sets[[var]] = NULL
      masked_indices[[var]] = NULL
    }
  }

  numeric_cols = names(input_data)[sapply(input_data, is.numeric)]

  mice_data = input_data %>%
    select(all_of(numeric_cols))

  set.seed(seed)

  mice_model = miceRanger(
    mice_data,
    m = m,
    maxiter = maxiter,
    returnModels = TRUE,
    verbose = FALSE,
    num.trees = num_trees,
    num.threads = parallel::detectCores() - 1
  )

  completed_data = completeData(mice_model)

  result_data = original_data

  for (var in variables_to_impute) {
    imputed_values = sapply(completed_data, function(dataset) {
      dataset[[var]]
    })

    imputed_mean = rowMeans(imputed_values)

    result_data[[paste0(var, "_imputed")]] = imputed_mean

    result_data[[paste0(var, "_final")]] = ifelse(
      is.na(original_data[[var]]),
      imputed_mean,
      original_data[[var]]
    )

    result_data[[paste0(var, "_flag")]] = ifelse(
      is.na(original_data[[var]]),
      "Imputed",
      "Original"
    )
  }

  performance = data.frame(variable = character(),
                           r_squared = numeric(),
                           rmse = numeric(),
                           mae = numeric(),
                           stringsAsFactors = FALSE)

  for (var in variables_to_impute) {
    if (!is.null(validation_sets[[var]])) {
      mask_idx = masked_indices[[var]]

      actual = validation_sets[[var]]

      predicted = result_data[[paste0(var, "_imputed")]][mask_idx]

      r2 = cor(actual, predicted)^2
      rmse = sqrt(mean((actual - predicted)^2))
      mae = mean(abs(actual - predicted))

      performance = rbind(performance,
                          data.frame(variable = var,
                                     r_squared = r2,
                                     rmse = rmse,
                                     mae = mae))
    }
  }

  return(list(
    imputed_data = result_data,
    performance = performance,
    mice_model = mice_model
  ))
}


#' Impute missing functional trait data using miceRanger
#'
#' Performs multiple imputation with random forests for species-level trait data.
#'
#' @param data Data frame with trait columns
#' @param variables_to_impute Character vector of trait column names
#' @param m Number of imputed datasets (default 5)
#' @param maxiter Number of iterations (default 10)
#' @param num_trees Trees per random forest (default 100)
#' @param seed Random seed (default 123)
#' @param validation_fraction Fraction of data to mask for validation (default 0.2)
#' @return List with imputed_data and performance
impute_functional_traits = function(
    data,
    variables_to_impute = c("LA", "LDMC", "N_percent", "SLA", "vegetative_height", "seed_mass"),
    m = 5,
    maxiter = 10,
    num_trees = 100,
    seed = 123,
    validation_fraction = 0.2
) {

  library(tidyverse)
  library(miceRanger)

  input_data = data

  original_data = input_data

  validation_sets = list()
  masked_indices = list()

  set.seed(seed)

  for (var in variables_to_impute) {
    complete_indices = which(!is.na(input_data[[var]]))

    if (length(complete_indices) > 5) {
      mask_size = floor(length(complete_indices) * validation_fraction)
      mask_idx = sample(complete_indices, mask_size)

      validation_sets[[var]] = input_data[[var]][mask_idx]

      input_data[[var]][mask_idx] = NA

      masked_indices[[var]] = mask_idx
    } else {
      validation_sets[[var]] = NULL
      masked_indices[[var]] = NULL
    }
  }

  numeric_cols = names(input_data)[sapply(input_data, is.numeric)]

  mice_data = input_data %>%
    select(all_of(numeric_cols))

  set.seed(seed)

  mice_model = miceRanger(
    mice_data,
    m = m,
    maxiter = maxiter,
    returnModels = TRUE,
    verbose = FALSE,
    num.trees = num_trees,
    num.threads = parallel::detectCores() - 1
  )

  completed_data = completeData(mice_model)

  result_data = original_data

  for (var in variables_to_impute) {
    imputed_values = sapply(completed_data, function(dataset) {
      dataset[[var]]
    })

    imputed_mean = rowMeans(imputed_values)

    result_data[[paste0(var, "_imputed")]] = imputed_mean

    result_data[[paste0(var, "_final")]] = ifelse(
      is.na(original_data[[var]]),
      imputed_mean,
      original_data[[var]]
    )

    result_data[[paste0(var, "_flag")]] = ifelse(
      is.na(original_data[[var]]),
      "Imputed",
      "Original"
    )
  }

  performance = data.frame(
    variable = character(),
    r_squared = numeric(),
    rmse = numeric(),
    stringsAsFactors = FALSE
  )

  for (var in variables_to_impute) {
    if (!is.null(validation_sets[[var]])) {
      predicted_values = result_data[[paste0(var, "_imputed")]][masked_indices[[var]]]
      actual_values = validation_sets[[var]]

      ss_total = sum((actual_values - mean(actual_values))^2)
      ss_residual = sum((actual_values - predicted_values)^2)
      r_squared = 1 - (ss_residual / ss_total)

      rmse = sqrt(mean((actual_values - predicted_values)^2))

      performance = rbind(performance, data.frame(
        variable = var,
        r_squared = r_squared,
        rmse = rmse,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(list(
    imputed_data = result_data,
    performance = performance
  ))
}


#' Calculate community-weighted mean traits
#'
#' Joins a community matrix (sites x species) with species trait data and
#' calculates community-weighted means per site. Optionally log- or
#' sqrt-transforms specified trait columns before computing CWMs.
#'
#' @param community_data Matrix or data frame of sites x species (presence/absence or abundance)
#' @param traits_data Data frame of species traits
#' @param species_col Column name for species in traits_data
#' @param abundance_col If provided, use abundance weighting
#' @param trait_cols Character vector of trait column names (auto-detected if NULL)
#' @param log_traits Character vector of trait column names to log-transform before CWM calculation (default NULL)
#' @param sqrt_traits Character vector of trait column names to sqrt-transform before CWM calculation (default NULL)
#' @return Data frame with CWM traits per site
calculate_community_traits = function(community_data, traits_data, species_col = "species",
                                       abundance_col = NULL, trait_cols = NULL,
                                       log_traits = NULL, sqrt_traits = NULL) {
  site_ids = rownames(community_data)
  community_data = as_tibble(community_data) %>%
    mutate(plot_id_releve = site_ids, .before = 1)

  traits_data = as_tibble(traits_data)

  # Remove index column if present
  if ("...1" %in% colnames(traits_data)) {
    traits_data = traits_data %>% select(-`...1`)
  }

  if (is.null(trait_cols)) {
    trait_cols = setdiff(colnames(traits_data), species_col)
  }

  # Apply log-transformation to specified trait columns
  if (!is.null(log_traits)) {
    traits_data = traits_data %>%
      mutate(across(all_of(log_traits), ~ log(.x)))
  }

  # Apply sqrt-transformation to specified trait columns
  if (!is.null(sqrt_traits)) {
    traits_data = traits_data %>%
      mutate(across(all_of(sqrt_traits), ~ sqrt(.x)))
  }

  species_cols = setdiff(colnames(community_data), "plot_id_releve")

  community_long = community_data %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "species",
      values_to = "value"
    ) %>%
    filter(value > 0)

  if (is.null(abundance_col)) {
    community_long = community_long %>%
      mutate(abundance = 1)
  } else {
    community_long = community_long %>%
      rename(abundance = value)
  }

  community_traits = community_long %>%
    left_join(traits_data, by = c("species" = species_col))

  # Trait coverage: proportion of total abundance with non-NA trait data
  # Computed before transformations would matter (NA pattern is the same)
  trait_coverage = community_traits %>%
    group_by(plot_id_releve) %>%
    summarise(
      across(
        all_of(trait_cols),
        list(
          coverage = ~ sum(abundance[!is.na(.x)]) / sum(abundance)
        )
      )
    )

  community_weighted_means = community_traits %>%
    group_by(plot_id_releve) %>%
    summarise(
      species_richness = n(),
      across(
        all_of(trait_cols),
        list(
          cwm = ~ weighted.mean(.x, w = abundance, na.rm = TRUE)
        )
      )
    ) %>%
    left_join(trait_coverage, by = "plot_id_releve")

  return(community_weighted_means)
}


#' Calculate unweighted community trait means and variances
#'
#' Joins a community matrix with species trait data and calculates
#' unweighted (equal weight per species) mean and variance per site.
#' Optionally log-transforms specified trait columns before computing.
#'
#' @param community_data Matrix or data frame of sites x species (presence/absence or abundance)
#' @param traits_data Data frame of species traits
#' @param species_col Column name for species in traits_data
#' @param trait_cols Character vector of trait column names
#' @param log_traits Character vector of trait column names to log-transform (default NULL)
#' @return Data frame with unweighted mean and variance traits per site
calculate_community_traits_unweighted = function(community_data, traits_data,
                                                  species_col = "species",
                                                  trait_cols, log_traits = NULL) {
  site_ids = rownames(community_data)
  comm_long = as_tibble(community_data) %>%
    mutate(plot_id_releve = site_ids, .before = 1) %>%
    pivot_longer(-plot_id_releve, names_to = "species", values_to = "value") %>%
    filter(value > 0)

  traits_tmp = as_tibble(traits_data)
  if ("...1" %in% colnames(traits_tmp)) traits_tmp = traits_tmp %>% select(-`...1`)

  if (!is.null(log_traits)) {
    traits_tmp = traits_tmp %>%
      mutate(across(all_of(log_traits), ~ log(.x)))
  }

  comm_traits = comm_long %>%
    left_join(traits_tmp, by = c("species" = species_col))

  comm_mean = comm_traits %>%
    group_by(plot_id_releve) %>%
    summarise(
      species_richness = n(),
      across(all_of(trait_cols),
             list(mean = ~ mean(.x, na.rm = TRUE)),
             .names = "{.col}_mean")
    )

  comm_var = comm_traits %>%
    group_by(plot_id_releve) %>%
    summarise(
      across(all_of(trait_cols),
             list(var = ~ var(.x, na.rm = TRUE)),
             .names = "{.col}_var")
    )

  comm_mean %>% left_join(comm_var, by = "plot_id_releve")
}


#' Calculate per-trait species functional distinctiveness
#'
#' For each trait separately: for each site, for each present species,
#' compute the mean absolute difference to all other co-occurring species.
#' Then for each species: take the median across all sites where it occurs.
#'
#' @param Y Binary site x species matrix (rows = sites, cols = species)
#' @param traits_data Data frame with species traits
#' @param species_col Column name for species in traits_data
#' @param trait_cols Character vector of trait column names to use
#' @param log_traits Character vector of trait columns to log-transform (default NULL)
#' @return Tibble with species and one distinctiveness column per trait
calculate_functional_distinctiveness = function(Y, traits_data, species_col = "species",
                                                trait_cols, log_traits = NULL) {
  traits_tmp = as_tibble(traits_data)
  if ("...1" %in% colnames(traits_tmp)) traits_tmp = traits_tmp %>% select(-`...1`)

  # Log-transform specified traits
  if (!is.null(log_traits)) {
    traits_tmp = traits_tmp %>%
      mutate(across(all_of(log_traits), ~ log(.x)))
  }

  # Build trait matrix aligned with Y columns
  sp_names = colnames(Y)
  trait_mat = traits_tmp %>%
    filter(.data[[species_col]] %in% sp_names) %>%
    arrange(match(.data[[species_col]], sp_names))

  trait_vals = as.matrix(trait_mat[, trait_cols])
  rownames(trait_vals) = trait_mat[[species_col]]
  complete = complete.cases(trait_vals)
  trait_vals = trait_vals[complete, , drop = FALSE]
  valid_sp = rownames(trait_vals)

  cat(sprintf("  %d / %d species with complete trait data\n",
              length(valid_sp), length(sp_names)))

  # For each trait, compute per-species median distinctiveness
  all_results = tibble(species = valid_sp)

  for (tr in trait_cols) {
    tr_vec = trait_vals[, tr]

    site_distinct = list()
    for (i in seq_len(nrow(Y))) {
      present = sp_names[Y[i, ] == 1]
      present = intersect(present, valid_sp)
      if (length(present) < 2) next

      vals = tr_vec[present]
      for (j in seq_along(present)) {
        mean_diff = mean(abs(vals[-j] - vals[j]))
        site_distinct[[length(site_distinct) + 1]] = tibble(
          species = present[j],
          diff = mean_diff
        )
      }
    }

    df_tr = bind_rows(site_distinct) %>%
      group_by(species) %>%
      summarise(med_distinct = median(diff, na.rm = TRUE), .groups = "drop")

    col_name = paste0("distinct_", gsub("^Median_", "", tr))
    all_results = all_results %>%
      left_join(df_tr %>% rename(!!col_name := med_distinct), by = "species")

    cat(sprintf("  %s: done\n", tr))
  }

  return(all_results)
}


# ==============================================================================
# 2. PLOTTING FUNCTIONS
# ==============================================================================

#' Abbreviate a species name to 4+4 format (e.g., "Agrostis capillaris" -> "AgroCapi")
#' Handles single-word names and lowercase input via str_to_title.
#'
#' @param x Character vector of species names
#' @return Character vector of abbreviated names
#' @used_by 03_merge_and_assess_traits.R (PCA biplots)
abbrev_species = function(x) {
  parts = str_split(x, "\\s+")
  sapply(parts, function(p) {
    if (length(p) >= 2) {
      paste0(str_to_title(str_sub(p[1], 1, 4)), str_to_title(str_sub(p[2], 1, 4)))
    } else {
      str_sub(p[1], 1, 8)
    }
  })
}


#' Custom Venn diagram for sjSDM ANOVA results
#'
#' Draws overlapping circles showing variance explained by environment,
#' species associations, and space (and their intersections).
#'
#' @param x sjSDM anova result object
#' @param y Unused (kept for compatibility)
#' @param type R-squared type: "McFadden", "Deviance", or "Nagelkerke"
#' @param fractions How to handle shared fractions: "discard", "proportional", "equal"
#' @param cols Colors for environment, associations, spatial circles
#' @param alpha Circle transparency
#' @param env_deviance Unused (kept for compatibility)
#' @param ... Additional arguments (unused)
#' @return Invisible list with Venn values
plot_anova_custom = function(x,
                             y,
                             type = c("McFadden", "Deviance", "Nagelkerke"),
                             fractions = c("discard", "proportional", "equal"),
                             cols = c("#81caf3", "#00bd89", "#d00000"),
                             alpha = 0.15,
                             env_deviance = NULL,
                             ci_data = NULL,
                             ...) {
  fractions = match.arg(fractions)
  lineSeq = 0.3
  nseg = 100
  dr = 1.0
  type = match.arg(type)
  out = list()
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))
  values = x$results
  select_rows =
    if (x$spatial) {
      sapply(c("F_A", "F_B", "F_AB", "F_S", "F_AS", "F_BS", "F_ABS"), function(i) which(values$models == i, arr.ind = TRUE))
    } else {
      sapply(c("F_A", "F_B", "F_AB"), function(i) which(values$models == i, arr.ind = TRUE))
    }
  values = values[select_rows, ]
  col_index =
    switch(type,
           Deviance = 4,
           Nagelkerke = 5,
           McFadden = 6
    )

  # Helper: format value by model name
  fmt_val = function(model_name) {
    row_idx = which(values$models == model_name)
    if (length(row_idx) == 0) return("")
    as.character(round(values[row_idx, col_index], 3))
  }

  # Circle centers — original layout
  cx_env = lineSeq          # 0.3
  cy_env = 1 - lineSeq      # 0.7
  cx_bio = 1 - lineSeq      # 0.7
  cy_bio = 1 - lineSeq      # 0.7
  cx_spa = 0.5
  cy_spa = lineSeq           # 0.3
  r_circ = 1.1 * lineSeq     # 0.33

  graphics::plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), pty = "s", axes = FALSE, xlab = "", ylab = "")
  xx = r_circ * cos(seq(0, 2 * pi, length.out = nseg))
  yy = r_circ * sin(seq(0, 2 * pi, length.out = nseg))

  # Draw circles
  graphics::polygon(xx + cx_env, yy + cy_env, col = add_alpha(cols[1], alpha = alpha), border = "black", lty = 1, lwd = 1)
  graphics::polygon(xx + cx_bio, yy + cy_bio, col = add_alpha(cols[2], alpha = alpha), border = "black", lty = 1, lwd = 1)
  graphics::text(cx_env, cy_env + r_circ - 0.05, labels = "Environment", cex = 1.2)
  graphics::text(cx_bio, cy_bio + r_circ - 0.05, labels = "Species\nAssociations", cex = 1.2)

  if (x$spatial) {
    graphics::polygon(xx + cx_spa, yy + cy_spa, col = add_alpha(cols[3], alpha = alpha), border = "black", lty = 1, lwd = 1)
    graphics::text(0.5, 0.07, labels = "Space", pos = 1, cex = 1.2)

    # --- Text positions for each Venn region ---
    # F_A: Env only — far left
    graphics::text(cx_env - 0.12, cy_env, labels = fmt_val("F_A"), cex = 1.3, font = 2)
    # F_B: Assoc only — far right
    graphics::text(cx_bio + 0.12, cy_bio, labels = fmt_val("F_B"), cex = 1.3, font = 2)
    # F_AB: Env ∩ Assoc, NOT Space — top center
    graphics::text(0.5, cy_env + 0.02, labels = fmt_val("F_AB"), cex = 1.2, font = 2)
    # F_S: Space only — bottom center
    graphics::text(0.5, cy_spa - 0.02, pos = 1, labels = fmt_val("F_S"), cex = 1.3, font = 2)
    # F_AS: Env ∩ Space — bottom-left
    graphics::text(0.3, 0.48, pos = 1, labels = fmt_val("F_AS"), cex = 1.1, font = 2)
    # F_BS: Assoc ∩ Space — bottom-right
    graphics::text(0.7, 0.48, pos = 1, labels = fmt_val("F_BS"), cex = 1.1, font = 2)
    # F_ABS: all three — center
    graphics::text(0.5, 0.55, labels = fmt_val("F_ABS"), cex = 1.2, font = 2)
  } else {
    graphics::text(cx_env - 0.12, cy_env, labels = fmt_val("F_A"), cex = 1.3, font = 2)
    graphics::text(cx_bio + 0.12, cy_bio, labels = fmt_val("F_B"), cex = 1.3, font = 2)
    graphics::text(0.5, cy_env, labels = fmt_val("F_AB"), cex = 1.2, font = 2)
  }
  out$VENN = values
  return(invisible(out))
}


#' Add alpha transparency to color(s)
#'
#' @param col Color string(s)
#' @param alpha Alpha value (0-1)
#' @return Color string(s) with alpha
add_alpha = function(col, alpha = 0.25) {
  apply(sapply(col, grDevices::col2rgb) / 255, 2, function(x) grDevices::rgb(x[1], x[2], x[3], alpha = alpha))
}


# ==============================================================================
# 3. ANALYSIS PLOTTING FUNCTIONS
# ==============================================================================

#' VP scatter plot: env vs codist, colored by spa, sized by r2. Used by 08_variance_partitioning.R.
build_vp_scatter_combined = function(df, label_species = FALSE) {
  # Pearson correlations
  cor_eb = cor.test(df$env, df$codist, method = "pearson")
  cor_es = cor.test(df$env, df$spa, method = "pearson")

  cor_label = paste0(
    "r(Env, Sp.Assoc.) = ", sprintf("%.2f", cor_eb$estimate),
    "\nr(Env, Space) = ", sprintf("%.2f", cor_es$estimate)
  )

  p = ggplot(df, aes(x = env, y = codist, color = spa, size = r2)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "grey85", high = color_spa, name = "Space") +
    scale_size_continuous(range = c(1, 6), name = "Variance explained") +
    annotate("text",
             x = 0.02,
             y = max(df$codist, na.rm = TRUE),
             label = cor_label, hjust = 0, vjust = 1, size = 4.5,
             fontface = "italic") +
    labs(x = "Environment", y = "Species Associations") +
    theme_bw(base_size = 18) +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.position = "right"
    )

  if (label_species && "species_name" %in% names(df)) {
    labeled = bind_rows(
      df %>% slice_max(env, n = 3),
      df %>% slice_max(spa, n = 3),
      df %>% slice_max(codist, n = 3),
      df %>% slice_min(r2, n = 2),
      df %>% slice_max(r2, n = 2)
    ) %>% distinct(species_name, .keep_all = TRUE)
    p = p + geom_text_repel(data = labeled, aes(label = species_name),
                            size = 3.5, fontface = "italic",
                            max.overlaps = 15, segment.color = "grey50",
                            show.legend = FALSE)
  }
  p
}


#' Partial effect plot for one env variable across VP components. Used by 09_environmental_gradient_analysis.R.
build_partial_plot = function(evar, models, df, predictor_cols,
                               vp_components, vp_labels, comp_colors,
                               var_full_names, x_label_prefix = "") {

  pred_parts = list()
  annot_parts = list()

  for (comp in vp_components) {
    best_mod = models[[comp]]
    s = summary(best_mod)
    retained = names(coef(best_mod))[-1]

    relevant = retained[grepl(paste0("^", evar, "$|^I\\(", evar), retained)]
    if (length(relevant) == 0) {
      annot_parts[[comp]] = sprintf("%s: not retained", vp_labels[comp])
      next
    }

    pred_df = as.data.frame(lapply(predictor_cols, function(v) rep(mean(df[[v]]), 200)))
    names(pred_df) = predictor_cols
    pred_df[[evar]] = seq(min(df[[evar]]), max(df[[evar]]), length.out = 200)
    pred_df$y = predict(best_mod, pred_df)
    pred_df$component = comp

    coef_info = as.data.frame(s$coefficients) %>%
      rownames_to_column("term") %>%
      filter(grepl(paste0("^", evar, "$|^I\\(", evar), term)) %>%
      mutate(
        short = ifelse(grepl("\\^2", term), "x\u00B2", "x"),
        sig = ifelse(`Pr(>|t|)` < 0.001, "***",
              ifelse(`Pr(>|t|)` < 0.01, "**",
              ifelse(`Pr(>|t|)` < 0.05, "*", ""))),
        label = sprintf("%s = %.4f%s", short, Estimate, sig)
      )

    # Check if any term for this variable is significant (p < 0.05)
    any_sig = any(coef_info$`Pr(>|t|)` < 0.05)
    pred_df$line_alpha = ifelse(any_sig, 1.0, 0.2)
    pred_parts[[comp]] = pred_df

    annot_parts[[comp]] = sprintf("%s: %s",
                                   vp_labels[comp],
                                   paste(coef_info$label, collapse = ", "))
  }

  df_preds = bind_rows(pred_parts)

  # Points in long format
  df_pts = df %>%
    select(all_of(paste0(vp_components, "_mean_val")), all_of(evar)) %>%
    pivot_longer(cols = all_of(paste0(vp_components, "_mean_val")),
                 names_to = "component", values_to = "value") %>%
    mutate(component = gsub("_mean_val", "", component))

  annot_text = paste(annot_parts[vp_components], collapse = "\n")
  full_name = ifelse(evar %in% names(var_full_names), var_full_names[evar], evar)
  x_lab = if (x_label_prefix != "") paste0(x_label_prefix, full_name) else full_name

  p = ggplot() +
    geom_point(data = df_pts,
               aes(x = .data[[evar]], y = value, color = component),
               size = 0.3, alpha = 0.3)

  if (nrow(df_preds) > 0) {
    p = p + geom_line(data = df_preds,
                      aes(x = .data[[evar]], y = y, color = component,
                          alpha = line_alpha),
                      linewidth = 1.2) +
      scale_alpha_identity()
  }

  p = p +
    scale_color_manual(values = comp_colors, labels = vp_labels, name = NULL) +
    annotate("text",
             x = min(df[[evar]]),
             y = max(df_pts$value, na.rm = TRUE),
             label = annot_text,
             hjust = 0, vjust = 1, size = 4, fontface = "italic") +
    labs(x = x_lab, y = "Variance explained") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 9),
          legend.position = "bottom",
          legend.text = element_text(size = 10))

  return(p)
}


#' Build paired species/community PDF for a set of env variables. Used by 09_environmental_gradient_analysis.R.
build_combined_pdf = function(vars, filename, pdf_width, pdf_height) {

  sp_plots = list()
  comm_plots = list()

  for (v in vars) {
    sp_plots[[v]] = build_partial_plot(
      evar = v, models = beta_models, df = df_merged,
      predictor_cols = beta_cols,
      vp_components = vp_components, vp_labels = vp_labels,
      comp_colors = comp_colors, var_full_names = var_full_names,
      x_label_prefix = "Response to "
    )

    comm_plots[[v]] = build_partial_plot(
      evar = v, models = community_env_models, df = df_wide,
      predictor_cols = env_cols,
      vp_components = vp_components, vp_labels = vp_labels,
      comp_colors = comp_colors, var_full_names = var_full_names,
      x_label_prefix = ""
    )
  }

  n = length(vars)

  # Column headers via ggtitle on first row only
  sp_plots[[1]] = sp_plots[[1]] + ggtitle("A) Species") +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))
  comm_plots[[1]] = comm_plots[[1]] + ggtitle("B) Community") +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0))

  # Interleave: sp1, comm1, sp2, comm2, ...
  # Hide legend on all individual panels; add a shared legend via collect
  paired = list()
  for (i in seq_along(vars)) {
    paired[[2 * i - 1]] = sp_plots[[i]] + theme(legend.position = "none")
    paired[[2 * i]]     = comm_plots[[i]] + theme(legend.position = "none")
  }

  # Create a standalone legend from a dummy plot with all 3 components
  legend_plot = ggplot(tibble(x = 1:3, y = 1:3,
                               component = factor(names(comp_colors), levels = names(comp_colors))),
                        aes(x = x, y = y, color = component)) +
    geom_point(size = 3) +
    scale_color_manual(values = comp_colors, labels = vp_labels, name = NULL) +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 12))
  legend_grob = cowplot::get_legend(legend_plot)

  p_combined = wrap_plots(paired, ncol = 2) +
    plot_annotation(theme = theme(plot.margin = margin(0, 0, 30, 0)))

  # Combine with legend at bottom
  p_final = cowplot::plot_grid(p_combined, legend_grob, ncol = 1,
                                rel_heights = c(1, 0.04))

  pdf(here("Calanda_JSDM", "plot", filename), width = pdf_width, height = pdf_height)
  print(p_final)
  dev.off()
  cat(sprintf("  Saved %s\n", filename))
}


# Global linetype mapping for paired plots
# Species median & Community mean -> solid
# Species plasticity & Community variance -> dashed
# Species distinctiveness -> dotted
global_ltype_values = c(
  "Median / Mean"            = "solid",
  "Plasticity / Variance"    = "dashed",
  "Distinctiveness"          = "dotted"
)

#' Paired partial effect plot (two variables overlaid with different linetypes). Used by 10_functional_post_analysis.R.
build_paired_plot = function(evar1, evar2, label1, label2, trait_name,
                              models, df, predictor_cols,
                              vp_components, vp_labels, comp_colors,
                              ltype1 = "Median / Mean",
                              ltype2 = "Plasticity / Variance",
                              evar3 = NULL, label3 = NULL,
                              ltype3 = "Distinctiveness",
                              annot_labels = NULL) {

  # Short labels for annotations (default to vp_labels if not provided)
  if (is.null(annot_labels)) annot_labels = vp_labels

  # Generate partial predictions for one variable (holding others at mean)
  # Also checks significance: lines for non-significant terms get alpha = 0.2
  get_preds = function(evar) {
    pred_parts = list()
    for (comp in vp_components) {
      best_mod = models[[comp]]
      s = summary(best_mod)
      retained = names(coef(best_mod))[-1]
      relevant = retained[grepl(paste0("^", evar, "$|^I\\(", evar), retained)]
      if (length(relevant) == 0) next

      pred_df = as.data.frame(lapply(predictor_cols, function(v) rep(mean(df[[v]]), 200)))
      names(pred_df) = predictor_cols
      pred_df[[evar]] = seq(min(df[[evar]]), max(df[[evar]]), length.out = 200)
      pred_df$y = predict(best_mod, pred_df)
      pred_df$component = comp
      pred_df$x_val = pred_df[[evar]]

      # Check if any term for this variable is significant
      p_vals = as.data.frame(s$coefficients) %>%
        rownames_to_column("term") %>%
        filter(grepl(paste0("^", evar, "$|^I\\(", evar), term)) %>%
        pull(`Pr(>|t|)`)
      pred_df$line_alpha = ifelse(any(p_vals < 0.05), 1.0, 0.2)

      pred_parts[[comp]] = pred_df
    }
    bind_rows(pred_parts)
  }

  # Annotation text for one variable
  get_annot = function(evar, label) {
    annot_parts = list()
    for (comp in vp_components) {
      best_mod = models[[comp]]
      s = summary(best_mod)
      retained = names(coef(best_mod))[-1]
      relevant = retained[grepl(paste0("^", evar, "$|^I\\(", evar), retained)]
      if (length(relevant) == 0) {
        annot_parts[[comp]] = sprintf("%s: not retained", annot_labels[comp])
        next
      }
      coef_info = as.data.frame(s$coefficients) %>%
        rownames_to_column("term") %>%
        filter(grepl(paste0("^", evar, "$|^I\\(", evar), term)) %>%
        mutate(
          short = ifelse(grepl("\\^2", term), "x\u00B2", "x"),
          sig = ifelse(`Pr(>|t|)` < 0.001, "***",
                ifelse(`Pr(>|t|)` < 0.01, "**",
                ifelse(`Pr(>|t|)` < 0.05, "*", ""))),
          lbl = sprintf("%s = %.4f%s", short, Estimate, sig)
        )
      annot_parts[[comp]] = sprintf("%s: %s", annot_labels[comp],
                                     paste(coef_info$lbl, collapse = ", "))
    }
    paste0(label, ": ", paste(annot_parts[vp_components], collapse = "; "))
  }

  preds1 = get_preds(evar1)
  preds2 = get_preds(evar2)

  annot1 = get_annot(evar1, label1)
  annot2 = get_annot(evar2, label2)
  annot_lines = c(annot1, annot2)

  if (nrow(preds1) > 0) preds1$ltype = ltype1
  if (nrow(preds2) > 0) preds2$ltype = ltype2

  # Optional third variable
  preds3_df = NULL
  if (!is.null(evar3) && !is.null(label3)) {
    preds3 = get_preds(evar3)
    annot3 = get_annot(evar3, label3)
    annot_lines = c(annot_lines, annot3)
    if (nrow(preds3) > 0) {
      preds3$ltype = ltype3
      preds3_df = preds3 %>% select(x_val, y, component, ltype, line_alpha)
    }
  }

  annot_text = paste(annot_lines, collapse = "\n")

  df_all_preds = bind_rows(
    if (nrow(preds1) > 0) preds1 %>% select(x_val, y, component, ltype, line_alpha) else NULL,
    if (nrow(preds2) > 0) preds2 %>% select(x_val, y, component, ltype, line_alpha) else NULL,
    preds3_df
  )

  # Points: VP values vs evar1 (the primary measure: median or mean)
  df_pts = df %>%
    select(all_of(paste0(vp_components, "_mean_val")), all_of(evar1)) %>%
    pivot_longer(cols = all_of(paste0(vp_components, "_mean_val")),
                 names_to = "component", values_to = "value") %>%
    mutate(component = gsub("_mean_val", "", component))

  # y range for annotation
  y_max = max(c(if (nrow(df_all_preds) > 0) df_all_preds$y else NULL,
                df_pts$value),
              na.rm = TRUE)
  x_min = min(df[[evar1]], df[[evar2]])

  # Invisible dummy rows to ensure all linetype levels appear in legend
  dummy_ltypes = tibble(
    x_val = NA_real_, y = NA_real_, component = vp_components[1],
    ltype = names(global_ltype_values)
  )

  p = ggplot() +
    geom_point(data = df_pts,
               aes(x = .data[[evar1]], y = value, color = component),
               size = 0.3, alpha = 0.3)

  if (nrow(df_all_preds) > 0) {
    p = p + geom_line(data = df_all_preds,
                      aes(x = x_val, y = y, color = component, linetype = ltype,
                          alpha = line_alpha),
                      linewidth = 1.2) +
      scale_alpha_identity()
  }

  p = p +
    geom_line(data = dummy_ltypes,
              aes(x = x_val, y = y, linetype = ltype),
              na.rm = TRUE, show.legend = TRUE) +
    scale_color_manual(values = comp_colors, labels = vp_labels, name = NULL) +
    scale_linetype_manual(values = global_ltype_values, name = NULL,
                          drop = FALSE,
                          guide = guide_legend(keywidth = unit(1.8, "cm"))) +
    annotate("text", x = x_min, y = y_max,
             label = annot_text,
             hjust = 0, vjust = 1, size = 3.3, fontface = "italic") +
    labs(x = trait_name, y = "Variance explained") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 9),
          legend.text = element_text(size = 10),
          legend.key.height = unit(0.7, "cm"))

  return(p)
}


# ==============================================================================
# 4. MODEL SUMMARY FUNCTIONS
# ==============================================================================

#' Extract stepwise AIC model results as tidy data frames
#'
#' Given a full lm model and its stepwise-selected best model, extracts:
#'   1. Coefficient table (estimates, SE, t, p, significance stars)
#'   2. BIC step history (which terms were dropped and BIC at each step)
#'   3. Model-level summary (adj R2, F-stat, p, df, AIC, BIC, n)
#'
#' @param full_mod The full lm model before stepwise selection
#' @param best_mod The final lm model after backward stepwise BIC/AIC
#' @param analysis_label Character label for the analysis (e.g., "env_gradient_community")
#' @param component_label Character label for the VP component (e.g., "env", "codist", "spa")
#' @param k_penalty Penalty per parameter for stepwise (default log(n) for BIC; use 2 for AIC)
#' @return List with $coefficients, $bic_steps, $model_summary (all tibbles)
#' @used_by 09_environmental_gradient_analysis.R, 10_functional_post_analysis.R
extract_model_results = function(full_mod, best_mod, analysis_label, component_label,
                                  k_penalty = NULL) {

  if (is.null(k_penalty)) k_penalty = 2  # AIC by default

  s = summary(best_mod)

  # 1. Coefficient table
  coef_tbl = as.data.frame(s$coefficients) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    rename(estimate = Estimate, std_error = `Std. Error`,
           t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    mutate(
      sig = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        p_value < 0.1   ~ ".",
        TRUE            ~ ""
      ),
      analysis = analysis_label,
      component = component_label
    ) %>%
    select(analysis, component, term, estimate, std_error, t_value, p_value, sig)

  # 2. BIC/AIC step history — re-run step with trace to capture steps
  step_log = capture.output(step(full_mod, direction = "backward", trace = 1, k = k_penalty))
  # Parse step lines: "- <term>   <df>  <sum_sq>  <BIC>"
  step_lines = step_log[grepl("^- ", step_log)]
  if (length(step_lines) > 0) {
    bic_steps = tibble(raw = step_lines) %>%
      mutate(
        term_dropped = str_trim(str_extract(raw, "(?<=^- )\\S.*?(?=\\s+\\d)")),
        bic = as.numeric(str_extract(raw, "[\\-\\.0-9]+$"))
      ) %>%
      mutate(
        step = seq_len(n()),
        analysis = analysis_label,
        component = component_label
      ) %>%
      select(analysis, component, step, term_dropped, bic)
  } else {
    bic_steps = tibble(
      analysis = character(), component = character(),
      step = integer(), term_dropped = character(), bic = numeric()
    )
  }

  # 3. Model-level summary
  f_stat = s$fstatistic
  model_p = if (!is.null(f_stat)) pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE) else NA
  model_summary = tibble(
    analysis = analysis_label,
    component = component_label,
    n_obs = nobs(best_mod),
    n_terms_full = length(coef(full_mod)) - 1,
    n_terms_final = length(coef(best_mod)) - 1,
    r_squared = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    f_statistic = if (!is.null(f_stat)) f_stat[1] else NA,
    df_model = if (!is.null(f_stat)) f_stat[2] else NA,
    df_residual = if (!is.null(f_stat)) f_stat[3] else NA,
    model_p_value = model_p,
    aic_full = AIC(full_mod),
    aic_final = AIC(best_mod),
    bic_full = BIC(full_mod),
    bic_final = BIC(best_mod),
    bic_drop = BIC(full_mod) - BIC(best_mod)
  )

  list(coefficients = coef_tbl, bic_steps = bic_steps, model_summary = model_summary)
}
