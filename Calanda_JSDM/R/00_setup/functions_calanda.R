# ==============================================================================
# functions_calanda.R
# Helper functions for the Calanda JSDM project
#
# Sections:
#   1. DATA PROCESSING: Snow, temperature, mowing, imputation
#   2. PLOTTING: Ternary plots, Venn diagram, labeling helpers
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


#' Replace NA values with data from the geographically closest site
#'
#' For sites missing et_annual, finds the nearest site with valid data
#' and copies its value.
#'
#' @param data Data frame with Longitude, Latitude, and et_annual columns
#' @return Data frame with NAs filled from nearest neighbour
replace_with_closest_site_data = function(data) {
  filled_data = data

  na_rows = which(is.na(data$et_annual))

  if (length(na_rows) > 0) {
    all_sites = data %>%
      filter(!is.na(et_annual)) %>%
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

    for (i in na_rows) {
      target_site = data[i, ] %>%
        st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

      distances = st_distance(target_site, all_sites)

      closest_idx = which.min(distances)

      filled_data$et_annual[i] = all_sites$et_annual[closest_idx]
    }
  }

  return(filled_data)
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


#' Extract mowing events from grassland-use intensity rasters
#'
#' Loads multi-year grassland-use rasters, extracts the mowingEvents layer,
#' and optionally extracts values at vegetation plot coordinates.
#'
#' @param raster_pattern Regex pattern to match raster filenames
#' @param coords_path Path to CSV with coordinates (x, y columns in WGS84)
#' @param output_path Path for output CSV
#' @param raster_dir Directory containing raster files
#' @param processed_raster_path Path to save/load processed raster stack
#' @param force_reprocess If TRUE, reprocess even if output exists
#' @return Data frame with extracted values, or SpatRaster if no coords provided
extract_mowing_events = function(raster_pattern = "grassland-use_intensity_.*\\.tif$",
                                 coords_path = NULL,
                                 output_path = "mowing_events_at_veg_points.csv",
                                 raster_dir = ".",
                                 processed_raster_path = "processed_mowing_events.tif",
                                 force_reprocess = FALSE) {

  library(tidyverse)
  library(terra)

  if (file.exists(processed_raster_path) && !force_reprocess) {
    message("Using existing processed raster from: ", processed_raster_path)
    mowing_events_stack = rast(processed_raster_path)
  } else {
    raster_files = list.files(path = raster_dir,
                              pattern = raster_pattern,
                              full.names = TRUE)

    if (length(raster_files) == 0) {
      stop("No raster files found matching the pattern")
    }

    mowing_events_list = list()

    for (i in seq_along(raster_files)) {
      r = rast(raster_files[i])

      if ("mowingEvents" %in% names(r)) {
        mowing_layer = r[["mowingEvents"]]

        year = str_extract(raster_files[i], "\\d{4}")
        names(mowing_layer) = paste0("mowingEvents_", year)

        mowing_events_list[[i]] = mowing_layer
      } else {
        warning(paste("mowingEvents layer not found in", raster_files[i]))
      }
    }

    if (length(mowing_events_list) == 0) {
      stop("No mowingEvents layers found in any of the raster files")
    }

    mowing_events_stack = do.call(c, mowing_events_list)

    writeRaster(mowing_events_stack, filename = processed_raster_path, overwrite = TRUE)
    message("Processed raster saved to: ", processed_raster_path)
  }

  if (!is.null(coords_path)) {
    veg_coords = read_csv(coords_path)

    veg_vector = vect(veg_coords, geom = c("x", "y"), crs = "EPSG:4326")

    if (crs(veg_vector) != crs(mowing_events_stack)) {
      veg_vector = project(veg_vector, crs(mowing_events_stack))
    }

    extracted_values = extract(mowing_events_stack, veg_vector)

    result = bind_cols(
      veg_coords,
      extracted_values %>%
        select(-ID)
    )

    if (!is.null(output_path)) {
      write_csv(result, output_path)
    }

    return(result)
  } else {
    return(mowing_events_stack)
  }
}


#' Calculate community-weighted mean traits
#'
#' Joins a community matrix (sites x species) with species trait data and
#' calculates community-weighted means per site.
#'
#' @param community_data Matrix or data frame of sites x species (presence/absence or abundance)
#' @param traits_data Data frame of species traits
#' @param species_col Column name for species in traits_data
#' @param abundance_col If provided, use abundance weighting
#' @param trait_cols Character vector of trait column names (auto-detected if NULL)
#' @return Data frame with CWM traits per site
calculate_community_traits = function(community_data, traits_data, species_col = "species", abundance_col = NULL, trait_cols = NULL) {
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
    )

  return(community_weighted_means)
}


# ==============================================================================
# 2. PLOTTING FUNCTIONS
# ==============================================================================

#' Label environmental variables for plotting
#'
#' @param var Character vector of variable names
#' @return Factor with human-readable labels
label_env_var = function(var) {
  factor(case_when(
    var == "summer_temp" ~ "Summer temperature",
    var == "fdd" ~ "Freezing degree days",
    var == "et_annual" ~ "Annual evapotranspiration",
    var == "land_use" ~ "Land use intensity"
  ), levels = c("Summer temperature", "Annual evapotranspiration",
                "Freezing degree days", "Land use intensity"))
}

#' Abbreviate a species name (first 4 letters of genus + first 4 of epithet)
#'
#' @param sp_name Full species name (e.g. "Ranunculus acris")
#' @return Abbreviated name (e.g. "RanuAcri")
abbrev_species = function(sp_name) {
  paste0(substring(word(sp_name, 1), 1, 4), toupper(substring(word(sp_name, 2), 1, 1)), substring(word(sp_name, 2), 2, 4))
}


#' Convert degrees to radians
deg2rad = function(deg) {
  return(deg * pi / 180)
}


#' Convert ternary coordinates to Cartesian
get_coords = function(X) {
  x_w = X[1] + sin(deg2rad(30)) * X[2]
  y_w = cos(deg2rad(30)) * X[2]
  return(c(x_w, y_w))
}


#' Species ternary plot (ggplot2) colored by mean altitude
#'
#' Creates a ternary plot of species-level variance components (environment,
#' species associations, space) with points colored by mean altitude and
#' labels for corner and center species.
#'
#' @param res JSDM results object (from sjSDM anova)
#' @param veg Vegetation data frame with plot_id_releve, species, species_cover
#' @param veg_clim Environmental data with plot_id_releve and altitude
#' @param color_env Color for environment axis
#' @param color_codist Color for species associations axis
#' @param color_spa Color for spatial axis
#' @param color_palette Color palette function (default gray gradient)
#' @param alpha Point transparency
#' @param cex Point size multiplier
#' @return ggplot object
plot_tern_species = function(res, veg, veg_clim,
                             color_env = "#81caf3",
                             color_codist = "#00bd89",
                             color_spa = "#d00000",
                             color_palette = NULL,
                             alpha = 0.7,
                             cex = 1.5) {

  species_altitude_ranges = veg_clim %>%
    select(plot_id_releve, altitude) %>%
    distinct() %>%
    left_join(
      veg %>% select(plot_id_releve, species, species_cover),
      by = "plot_id_releve"
    ) %>%
    filter(species_cover > 0) %>%
    group_by(species) %>%
    summarize(
      min_altitude = min(altitude, na.rm = TRUE),
      max_altitude = max(altitude, na.rm = TRUE),
      altitude_range = max_altitude - min_altitude,
      mean_altitude = mean(altitude, na.rm = TRUE),
      .groups = "drop"
    )

  species_data = res$internals$Species %>%
    rownames_to_column("species") %>%
    left_join(species_altitude_ranges, by = "species") %>%
    filter(!is.na(altitude_range))

  species_data = species_data %>%
    mutate(
      total = env + codist + spa,
      env_prop = env / total,
      codist_prop = codist / total,
      spa_prop = spa / total
    )

  species_data$x = NA
  species_data$y = NA
  for (i in 1:nrow(species_data)) {
    total = species_data$env[i] + species_data$codist[i] + species_data$spa[i]
    env_norm = species_data$env[i] / total
    codist_norm = species_data$codist[i] / total
    spa_norm = species_data$spa[i] / total
    coords = get_coords(c(spa_norm, codist_norm, env_norm))
    species_data$x[i] = coords[1]
    species_data$y[i] = coords[2]
  }

  if (is.null(color_palette)) {
    color_palette = colorRampPalette(c("gray90", "gray20"))
  }

  species_data$dist_env = sqrt(species_data$x^2 + species_data$y^2)
  species_data$dist_codist = sqrt((species_data$x - 0.5)^2 + (species_data$y - 0.866)^2)
  species_data$dist_spa = sqrt((species_data$x - 1)^2 + species_data$y^2)

  corner_species = bind_rows(
    species_data %>% filter(env_prop > 0.8) %>% arrange(dist_env) %>% head(3),
    species_data %>% filter(codist_prop > 0.8) %>% arrange(dist_codist) %>% head(3),
    species_data %>% filter(spa_prop > 0.8) %>% arrange(dist_spa) %>% head(3)
  ) %>%
    distinct(species, .keep_all = TRUE)

  center_species = species_data %>%
    filter(
      env_prop >= 0.25 & env_prop <= 0.40,
      codist_prop >= 0.25 & codist_prop <= 0.40,
      spa_prop >= 0.25 & spa_prop <= 0.40
    ) %>%
    mutate(dist_center = sqrt((x - 0.5)^2 + (y - 0.289)^2)) %>%
    arrange(dist_center) %>%
    head(3)

  labeled_species = bind_rows(corner_species, center_species) %>%
    distinct(species, .keep_all = TRUE) %>%
    mutate(label = abbrev_species(species))

  triangle_edges = data.frame(
    x_start = c(0, 1, 0.5),
    y_start = c(0, 0, 0.866),
    x_end = c(1, 0.5, 0),
    y_end = c(0, 0.866, 0),
    edge = c("spa", "codist", "env")
  )

  grid_lines = data.frame()
  grid_labels = data.frame()

  for (val in seq(0.2, 0.8, 0.2)) {
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val / 2, 1 - val / 2),
      y = c(val * 0.866, val * 0.866),
      group = paste0("h", val),
      axis = "codist"
    ))
    label_pos = get_coords(c(1 - val, 0.0, val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1],
      y = label_pos[2] - 0.03,
      label = as.character(1 - val),
      angle = 60
    ))

    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val, 0.5 + val / 2),
      y = c(0, (1 - val) * 0.866),
      group = paste0("l", val),
      axis = "spa"
    ))
    label_pos = get_coords(c(1 - val, val, 0.0))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] + 0.05,
      y = label_pos[2],
      label = as.character(val),
      angle = 0
    ))

    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(1 - val, 0.5 - val / 2),
      y = c(0, (1 - val) * 0.866),
      group = paste0("r", val),
      axis = "env"
    ))
    label_pos = get_coords(c(0, val, 1 - val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] - 0.03,
      y = label_pos[2] + 0.03,
      label = as.character(1 - val),
      angle = -50
    ))
  }

  p = ggplot() +
    geom_segment(data = triangle_edges,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = edge),
                 linewidth = 1.2) +
    scale_color_manual(values = c("env" = color_env, "codist" = color_codist, "spa" = color_spa),
                       guide = "none") +
    geom_line(data = grid_lines, aes(x = x, y = y, group = group, color = axis),
              linewidth = 0.5, alpha = 0.3) +
    geom_text(data = grid_labels, aes(x = x, y = y, label = label, angle = angle),
              size = 5, color = "gray30") +
    geom_point(data = species_data, aes(x = x, y = y, fill = mean_altitude),
               shape = 21, size = cex * 2, alpha = alpha, color = "black", stroke = 0.2) +
    geom_text_repel(data = labeled_species, aes(x = x, y = y, label = label),
                    size = 4, fontface = "italic",
                    force = 10,
                    min.segment.length = 0.01,
                    segment.color = "gray50") +
    scale_fill_gradientn(colors = color_palette(100),
                         name = "Mean\naltitude (m)",
                         guide = guide_colorbar(barwidth = 0.8, barheight = 8)) +
    annotate("text", x = 0, y = 0, label = "Environment",
             color = color_env, fontface = "bold", hjust = 1.2, size = 4.5) +
    annotate("text", x = 0.5, y = 0.95, label = "Species associations",
             color = color_codist, fontface = "bold", hjust = 1, size = 4.5) +
    annotate("text", x = 1, y = 0.01, label = "Space",
             color = color_spa, fontface = "bold", hjust = -0.6, size = 4.5) +
    ggtitle("Species") +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = -0.1, face = "bold", size = 12),
      legend.position = "right",
      plot.margin = margin(20, 10, 30, 40)
    )

  return(p)
}


#' Sites ternary plot (ggplot2) colored by altitude
#'
#' Creates a ternary plot of site-level variance components (environment,
#' species associations, space) with points colored by altitude.
#'
#' @param res JSDM results object (from sjSDM anova)
#' @param veg_clim Environmental data with plot_id_releve and altitude
#' @param color_env Color for environment axis
#' @param color_codist Color for species associations axis
#' @param color_spa Color for spatial axis
#' @param color_palette Color palette function (default gray gradient)
#' @param alpha Point transparency
#' @param cex Point size multiplier
#' @return ggplot object
plot_tern_sites = function(res, veg_clim,
                           color_env = "#81caf3",
                           color_codist = "#00bd89",
                           color_spa = "#d00000",
                           color_palette = NULL,
                           alpha = 0.7,
                           cex = 1.5) {

  library(ggplot2)

  sites_data = res$internals$Sites %>%
    rownames_to_column("plot_id_releve") %>%
    left_join(veg_clim %>% select(plot_id_releve, altitude) %>% distinct(),
              by = "plot_id_releve") %>%
    filter(!is.na(altitude))

  sites_data$x = NA
  sites_data$y = NA
  for (i in 1:nrow(sites_data)) {
    total = sites_data$env[i] + sites_data$codist[i] + sites_data$spa[i]
    env_norm = sites_data$env[i] / total
    codist_norm = sites_data$codist[i] / total
    spa_norm = sites_data$spa[i] / total
    coords = get_coords(c(spa_norm, codist_norm, env_norm))
    sites_data$x[i] = coords[1]
    sites_data$y[i] = coords[2]
  }

  if (is.null(color_palette)) {
    color_palette = colorRampPalette(c("gray90", "gray20"))
  }

  triangle_edges = data.frame(
    x_start = c(0, 1, 0.5),
    y_start = c(0, 0, 0.866),
    x_end = c(1, 0.5, 0),
    y_end = c(0, 0.866, 0),
    edge = c("spa", "codist", "env")
  )

  grid_lines = data.frame()
  grid_labels = data.frame()

  for (val in seq(0.2, 0.8, 0.2)) {
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val / 2, 1 - val / 2),
      y = c(val * 0.866, val * 0.866),
      group = paste0("h", val),
      axis = "codist"
    ))
    label_pos = get_coords(c(1 - val, 0.0, val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1],
      y = label_pos[2] - 0.03,
      label = as.character(1 - val),
      angle = 60
    ))

    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val, 0.5 + val / 2),
      y = c(0, (1 - val) * 0.866),
      group = paste0("l", val),
      axis = "spa"
    ))
    label_pos = get_coords(c(1 - val, val, 0.0))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] + 0.05,
      y = label_pos[2],
      label = as.character(val),
      angle = 0
    ))

    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(1 - val, 0.5 - val / 2),
      y = c(0, (1 - val) * 0.866),
      group = paste0("r", val),
      axis = "env"
    ))
    label_pos = get_coords(c(0, val, 1 - val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] - 0.03,
      y = label_pos[2] + 0.03,
      label = as.character(1 - val),
      angle = -50
    ))
  }

  p = ggplot() +
    geom_segment(data = triangle_edges,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = edge),
                 linewidth = 1.2) +
    scale_color_manual(values = c("env" = color_env, "codist" = color_codist, "spa" = color_spa),
                       guide = "none") +
    geom_line(data = grid_lines, aes(x = x, y = y, group = group, color = axis),
              linewidth = 0.5, alpha = 0.3) +
    geom_text(data = grid_labels, aes(x = x, y = y, label = label, angle = angle),
              size = 5, color = "gray30") +
    geom_point(data = sites_data, aes(x = x, y = y, fill = altitude),
               shape = 21, size = cex * 2, alpha = alpha, color = "black", stroke = 0.2) +
    scale_fill_gradientn(colors = color_palette(100),
                         name = "Altitude\n(m)",
                         guide = guide_colorbar(barwidth = 0.8, barheight = 8)) +
    annotate("text", x = 0, y = 0, label = "Environment",
             color = color_env, fontface = "bold", hjust = 1.2, size = 4.5) +
    annotate("text", x = 0.5, y = 0.95, label = "Species associations",
             color = color_codist, fontface = "bold", hjust = 1, size = 4.5) +
    annotate("text", x = 1, y = 0.01, label = "Space",
             color = color_spa, fontface = "bold", hjust = -0.6, size = 4.5) +
    ggtitle("Sites") +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = -0.1, face = "bold", size = 12),
      legend.position = "right",
      plot.margin = margin(20, 10, 30, 40)
    )

  return(p)
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
  graphics::plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), pty = "s", axes = FALSE, xlab = "", ylab = "")
  xx = 1.1 * lineSeq * cos(seq(0, 2 * pi, length.out = nseg))
  yy = 1.1 * lineSeq * sin(seq(0, 2 * pi, length.out = nseg))
  graphics::polygon(xx + lineSeq, yy + (1 - lineSeq), col = add_alpha(cols[1], alpha = alpha), border = "black", lty = 1, lwd = 1)
  graphics::text(lineSeq - 0.1, (1 - lineSeq), labels = round(values[1, col_index], 3))
  graphics::text(mean(xx + lineSeq), 0.9, labels = "Environmental", pos = 3)
  graphics::polygon(xx + 1 - lineSeq, yy + 1 - lineSeq, col = add_alpha(cols[2], alpha = alpha), border = "black", lty = 1, lwd = 1)
  graphics::text(1 - lineSeq + 0.1, (1 - lineSeq), labels = round(values[2, col_index], 3))
  graphics::text(1 - mean(xx + lineSeq), 0.9, labels = "Associations", pos = 3)
  graphics::text(0.5, (1 - lineSeq), labels = round(values[3, col_index], 3))
  if (x$spatial) {
    graphics::polygon(xx + 0.5, yy + lineSeq, col = add_alpha(cols[3], alpha = alpha), border = "black", lty = 1, lwd = 1)
    graphics::text(0.5, lineSeq + 0.0, pos = 1, labels = round(values[4, col_index], 3))
    graphics::text(0.5, 0.1, labels = "Spatial", pos = 1)
    graphics::text(0.3, 0.5, pos = 1, labels = round(values[5, col_index], 3))
    graphics::text(1 - 0.3, 0.5, pos = 1, labels = round(values[6, col_index], 3))
    graphics::text(0.5, 0.5 + 0.05, labels = round(values[7, col_index], 3))
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
