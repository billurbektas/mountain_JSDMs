bioclim_data = data.frame(
  Variable = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", 
               "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "AI", "etpY","gdd0Y", "gdd3Y", "gdd5Y"),
  Description = c("Annual Mean Temperature",
                  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
                  "Isothermality",
                  "Temperature Seasonality (standard deviation ×100)",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter",
                  "Aridity Index", 
                  "Evapotranspiration",
                  "Sum of growing degree days (above 0 degrees C)",
                  "Sum of growing degree days (above 3 degrees C)",
                  "Sum of growing degree days (above 5 degrees C)"))
match_clim_var = function(var){
bioclim_data = data.frame(
  Variable = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", 
               "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "AI", "etpY","gdd0Y", "gdd3Y", "gdd5Y"),
  Description = c("Annual Mean Temperature",
                  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
                  "Isothermality",
                  "Temperature Seasonality (standard deviation ×100)",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter",
                  "Aridity Index", 
                  "Evapotranspiration",
                  "Sum of growing degree days (above 0 degrees C)",
                  "Sum of growing degree days (above 3 degrees C)",
                  "Sum of growing degree days (above 5 degrees C)"))
  return(bioclim_data$Description[bioclim_data$Variable == var])
}

aucplot =  function(auc) {
  mean.e=signif(mean(auc$e.aucs,na.rm = T),digits = 2)  
  mean.p=signif(mean(auc$p.aucs,na.rm = T),digits = 2) 
  ggplot(auc,aes(x =e.aucs, y =p.aucs))+
    scale_color_viridis(option="D",discrete=TRUE)+
    xlab("explanatory AUC") + 
    ylab("predictive AUC") + 
    xlim(0.3,1)+
    ylim(0.3,1)+
    geom_abline(intercept=0,slope=1,colour='azure3')+
    geom_hline(yintercept = 0.5,colour='azure3')+
    geom_vline(xintercept = 0.5,colour='azure3')+
    geom_hline(yintercept = mean.p,linetype="dotted")+
    geom_vline(xintercept = mean.e,linetype="dotted")+
    labs(fill="Species Index") + 
    geom_point()+ 
    ggplot2::annotate("text",label=mean.p, x = 0.45, y = mean.p+0.02)+
    ggplot2::annotate("text",label=mean.e, x = 0.75, y = 0.4)+
    theme(panel.grid.major =element_blank(), panel.grid.minor  =element_blank(),panel.background = element_blank(),panel.border = element_blank())+
    geom_text_repel(aes(e.aucs, p.aucs, label=sp,size=20,fontface = "italic"),show.legend = FALSE)
}


process_gfsc_data <- function(copernicus_dir, 
                              output_dir,
                              veg_coords_path,
                              calanda_mask_path,
                              force_reprocess = FALSE) {
  
  # Load required packages silently
  suppressPackageStartupMessages({
    require(tidyverse)
    require(terra)
    require(lubridate)
    require(sf)
  })
  
  # Set output file path
  output_file <- file.path(output_dir, "gf_vegetation_data.csv")
  
  # Check if output file already exists
  if (file.exists(output_file) && !force_reprocess) {
    message("Output file already exists. Set force_reprocess=TRUE to reprocess.")
    return(output_file)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get list of all zip files
  zip_files <- list.files(copernicus_dir, pattern = "GFSC_\\d{8}-007_S1-S2_T32T\\w{2}_V101_\\d+\\.zip$", full.names = TRUE)
  
  if (length(zip_files) == 0) {
    stop("No matching ZIP files found in ", copernicus_dir)
  }
  
  # Check if vegetation coordinates file exists
  if (!file.exists(veg_coords_path)) {
    stop("Vegetation coordinates file not found: ", veg_coords_path)
  }
  
  # Check if calanda mask shapefile exists
  # For shapefiles, we need to check if the .shp file exists
  if (!file.exists(calanda_mask_path)) {
    # If direct path doesn't exist, try adding .shp extension if not already present
    if (!grepl("\\.shp$", calanda_mask_path) && !file.exists(paste0(calanda_mask_path, ".shp"))) {
      stop("Calanda mask shapefile not found: ", calanda_mask_path)
    }
  }
  
  # Load the calanda mask shapefile using sf
  calanda_mask_sf <- st_read(calanda_mask_path, quiet = TRUE)
  
  # Function to extract date from filename
  extract_date <- function(filename) {
    # Extract the date portion (YYYYMMDD) from the filename
    date_str <- str_extract(basename(filename), "\\d{8}")
    # Convert to Date object
    return(ymd(date_str))
  }
  
  # Read vegetation coordinates from CSV
  veg_coords <- read_csv(veg_coords_path)
  
  # Create sf object from veg_coords
  # Based on the screenshot, we have plot_id, releve_id, x, y columns
  veg_points_sf <- st_as_sf(veg_coords, 
                            coords = c("x", "y"),
                            crs = 4326)  # Assuming WGS84 coordinates, adjust if needed
  
  # Function to safely create directory and check space
  safe_create_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
      # Check disk space before creating directory
      disk_info <- system2("df", args = c("-k", dirname(dir_path)), stdout = TRUE)
      avail_space <- as.numeric(strsplit(disk_info[2], "\\s+")[[1]][4])
      
      # If less than 100MB available, clean up temporary directories
      if (avail_space < 100 * 1024) {
        message("Low disk space detected. Cleaning up extraction directory...")
        temp_dirs <- list.dirs(output_dir, recursive = FALSE)
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
  
  # Create temp extraction directory
  temp_extract_dir <- file.path(tempdir(), "gfsc_extract")
  safe_create_dir(temp_extract_dir)
  
  # Process each zip file
  results <- map_df(zip_files, function(zip_file) {
    tryCatch({
      # Extract date from filename
      date <- extract_date(zip_file)
      
      # Clear temp directory before processing new file
      if (dir.exists(temp_extract_dir)) {
        unlink(list.files(temp_extract_dir, full.names = TRUE), recursive = TRUE)
      }
      
      message("Processing ", basename(zip_file), " (", date, ")")
      
      # Unzip the file to a temporary directory
      message("  Extracting...")
      extracted_files <- unzip(zip_file, exdir = temp_extract_dir, list = TRUE)
      unzip(zip_file, exdir = temp_extract_dir)
      
      # Find the GF.tif file
      gf_file <- list.files(temp_extract_dir, pattern = "_GF\\.tif$", recursive = TRUE, full.names = TRUE)
      
      # If no GF.tif file found, return empty tibble
      if (length(gf_file) == 0) {
        warning(paste("No GF.tif file found in", zip_file))
        return(tibble())
      }
      
      message("  Found GF file: ", basename(gf_file))
      
      # Read the raster
      gf_raster <- rast(gf_file)
      
      # Get the CRS from calanda mask
      calanda_crs <- st_crs(calanda_mask_sf)
      
      # Reproject the raster to match calanda_mask projection
      message("  Reprojecting raster...")
      gf_raster_reproj <- project(gf_raster, calanda_crs$wkt)
      
      # Transform vegetation points to match the raster CRS
      if (st_crs(veg_points_sf) != calanda_crs) {
        veg_points_sf <- st_transform(veg_points_sf, calanda_crs)
      }
      
      # Extract values from raster
      message("  Extracting values at vegetation points...")
      extracted_values <- terra::extract(gf_raster_reproj, veg_points_sf)
      
      # Convert extracted values to a data frame
      extracted_df <- as.data.frame(extracted_values)
      
      # Combine the original veg_coords with extracted values and add date
      result <- veg_coords %>%
        bind_cols(raster_value = extracted_df[[2]]) %>%
        mutate(date = date)
      
      # Clean up extracted files immediately after processing
      message("  Cleaning up extracted files...")
      unlink(list.files(temp_extract_dir, full.names = TRUE), recursive = TRUE)
      
      return(result)
    },
    error = function(e) {
      warning(paste("Error processing", zip_file, ":", e$message))
      # Clean up any partial extraction
      unlink(list.files(temp_extract_dir, full.names = TRUE), recursive = TRUE)
      return(tibble())
    })
  })
  
  # Clean up temporary extraction directory
  message("Final cleanup...")
  unlink(temp_extract_dir, recursive = TRUE)
  
  # If no results were obtained, return empty dataframe with correct structure
  if (nrow(results) == 0) {
    warning("No data was extracted from any of the files")
    # Create empty dataframe with structure matching veg_coords plus date and raster_value
    final_results <- veg_coords[0, ] %>%
      mutate(date = as.Date(character()),
             raster_value = numeric())
  } else {
    # Clean up and organize final results
    final_results <- results %>%
      # Make sure plot_id, releve_id columns are preserved
      select(date, plot_id, releve_id, x, y, raster_value) %>%
      arrange(date, plot_id)
  }
  
  # Save the results
  write_csv(final_results, output_file)
  
  # Print summary
  message("Processed ", n_distinct(final_results$date), " dates with ", nrow(final_results), " total observations")
  message("Data saved to ", output_file)
  
  # Create a plot of values over time for each point (only if we have data)
  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(final_results) > 0) {
    p <- final_results %>%
      ggplot(aes(x = date, y = raster_value, color = factor(plot_id))) +
      geom_line() +
      geom_point() +
      labs(title = "Raster Values Over Time by Plot ID",
           x = "Date",
           y = "Raster Value",
           color = "Plot ID") +
      theme_minimal()
    
    plot_path <- file.path(output_dir, "temporal_plot.png")
    ggsave(plot_path, p, width = 10, height = 6)
    message("Plot saved to ", plot_path)
  } else if(nrow(final_results) == 0) {
    message("No data to plot")
  } else {
    message("ggplot2 package not available, skipping plot creation")
  }
  
  return(output_file)
}

# Function to preprocess snow cover data with BISE correction and interpolation
process_snow_data <- function(raw_snow_df, 
                              bise_threshold = 0.2,
                              sliding_window = 3,
                              sgfilter_p = 3, 
                              sgfilter_n = 7,
                              expand_time_series = TRUE,
                              date_range = NULL) {
  # Ensure library dependencies are loaded
  suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(zoo)
    require(signal)
    require(lubridate)
  })
  
  # Check inputs
  required_cols <- c("date", "raster_value", "plot_id", "releve_id")
  missing_cols <- setdiff(required_cols, colnames(raw_snow_df))
  if (length(missing_cols) > 0) {
    stop("Input dataframe missing required columns: ", paste(missing_cols, collapse=", "))
  }
  
  # Create unique identifier for plot_id and releve_id combination
  raw_snow_df <- raw_snow_df %>%
    mutate(plot_id_releve = paste(plot_id, releve_id, sep="_")) %>%
    arrange(plot_id_releve, date)
  
  # Define global date range for all time series
  if (expand_time_series) {
    if (!is.null(date_range)) {
      # Use user-provided date range
      global_date_range <- date_range
      message("Using user-provided date range: ", 
              format(global_date_range[1], "%Y-%m-%d"), " to ", 
              format(global_date_range[2], "%Y-%m-%d"))
    } else {
      # Use min/max dates from the entire dataset
      global_date_range <- range(raw_snow_df$date, na.rm = TRUE)
      message("Using global date range from data: ", 
              format(global_date_range[1], "%Y-%m-%d"), " to ", 
              format(global_date_range[2], "%Y-%m-%d"))
    }
  } else {
    global_date_range <- NULL
    message("Time series will not be expanded (using location-specific date ranges)")
  }
  
  # Get all unique plot_id_releve combinations
  plot_releve_combos <- unique(raw_snow_df$plot_id_releve)
  message("Processing ", length(plot_releve_combos), " unique plot-releve combinations...")
  
  # Initialize list to store results
  all_results <- list()
  
  # Process each plot_id/releve_id combination separately
  for (curr_combo in plot_releve_combos) {
    # Extract data for current combination
    curr_data <- raw_snow_df %>% 
      filter(plot_id_releve == curr_combo)
    
    # Skip if less than 3 valid observations (minimum needed for interpolation)
    if (nrow(curr_data) < 3 || sum(!is.na(curr_data$raster_value)) < 3) {
      message("Skipping ", curr_combo, " - insufficient data points")
      next
    }
    
    # Ensure data is sorted by date
    curr_data <- curr_data %>% arrange(date)
    
    message("Processing ", curr_combo, " with ", nrow(curr_data), " observations")
    
    # Extract raw values
    snow_raw <- curr_data$raster_value
    snow_processed <- snow_raw
    
    # 1. BISE correction to remove sudden drops/spikes
    if (length(snow_raw) > 1) {
      # Calculate differences between consecutive observations
      diff_vals <- diff(snow_raw)
      
      # Identify points showing decrease
      decrease_idx <- which(diff_vals < 0) + 1  
      
      # Using the specified slope threshold
      data_threshold <- c(NA, snow_raw[-1] - bise_threshold * diff_vals)
      
      # Forward sliding period
      if (length(snow_raw) >= sliding_window) {
        val_sliding_period <- zoo::rollapply(snow_raw, 
                                             width = sliding_window,
                                             FUN = max, 
                                             fill = NA, 
                                             align = "left")
        
        # Identify points to reject based on decrease
        if (length(decrease_idx) > 0) {
          # Check if indices are valid
          valid_idx <- decrease_idx[decrease_idx <= length(val_sliding_period) & 
                                      decrease_idx <= length(data_threshold)]
          
          if (length(valid_idx) > 0) {
            law_check <- val_sliding_period[valid_idx] - data_threshold[valid_idx]
            reject_decrease <- valid_idx[which(law_check > 0)]
            
            # Apply corrections for sudden decreases
            snow_processed[reject_decrease] <- NA
          }
        }
        
        # Check for sudden increases
        max_val <- max(snow_raw, na.rm = TRUE)
        if (!is.infinite(max_val)) {
          increase_threshold <- bise_threshold * max_val
          increase_idx <- which(diff_vals > increase_threshold)
          
          if (length(increase_idx) > 0) {
            reject_increase <- increase_idx[!increase_idx %in% decrease_idx]
            
            # Check for valid indices
            valid_increases <- reject_increase[reject_increase < length(snow_processed)]
            
            if (length(valid_increases) > 0) {
              # Apply corrections for sudden increases
              snow_processed[valid_increases + 1] <- NA
            }
          }
        }
      }
    }
    
    # Skip if too many NA values after correction
    if (sum(!is.na(snow_processed)) < 3) {
      message("  Skipping ", curr_combo, " - insufficient data after BISE correction")
      next
    }
    
    # 2. Fill start and end gaps if they exist
    first_valid <- which(!is.na(snow_processed))[1]
    if (!is.na(first_valid) && first_valid > 1) {
      snow_processed[1:(first_valid-1)] <- snow_processed[first_valid]
    }
    
    last_valid <- max(which(!is.na(snow_processed)))
    if (!is.na(last_valid) && last_valid < length(snow_processed)) {
      snow_processed[(last_valid+1):length(snow_processed)] <- snow_processed[last_valid]
    }
    
    # 3. Interpolate missing values using spline
    # Only proceed if there are still NA values to interpolate
    if (any(is.na(snow_processed))) {
      # Only interpolate if we have enough non-NA values
      if (sum(!is.na(snow_processed)) >= 3) {
        snow_interpolated <- tryCatch({
          zoo::na.spline(snow_processed)
        }, error = function(e) {
          message("  Error in spline interpolation for ", curr_combo, ": ", e$message)
          # Fall back to linear interpolation if spline fails
          approx(x = which(!is.na(snow_processed)), 
                 y = snow_processed[!is.na(snow_processed)],
                 xout = 1:length(snow_processed), 
                 rule = 2)$y
        })
      } else {
        # Not enough points for interpolation, skip this combination
        message("  Skipping ", curr_combo, " - insufficient data for interpolation")
        next
      }
    } else {
      # No NAs to interpolate
      snow_interpolated <- snow_processed
    }
    
    # 4. Apply Savitzky-Golay filter for smoothing
    # Only apply if we have enough points
    if (length(snow_interpolated) >= sgfilter_n) {
      snow_smoothed <- tryCatch({
        signal::sgolayfilt(snow_interpolated, p = sgfilter_p, n = sgfilter_n, m = 0)
      }, error = function(e) {
        message("  Error in Savitzky-Golay filter for ", curr_combo, ": ", e$message)
        # Return interpolated values if filter fails
        snow_interpolated
      })
    } else {
      # Not enough points for Savitzky-Golay filter
      snow_smoothed <- snow_interpolated
    }
    
    # 5. Create daily interpolation with expanded time series
    # Get the input date range for this location
    loc_date_range <- range(curr_data$date)
    
    # If global_date_range is provided, use it instead to cover the entire study period
    if (exists("global_date_range") && !is.null(global_date_range)) {
      expanded_range <- global_date_range
    } else {
      # Otherwise, use the date range for this location
      expanded_range <- loc_date_range
    }
    
    # Create sequence of daily dates for the expanded range
    dates_seq <- seq(expanded_range[1], expanded_range[2], by = "days")
    
    # Handle interpolation based on number of observations and expanded range
    # Check if we need extrapolation (dates outside observed range)
    needs_extrapolation <- min(dates_seq) < min(curr_data$date) || max(dates_seq) > max(curr_data$date)
    
    if (needs_extrapolation) {
      message("  Extrapolating outside observed range for ", curr_combo)
      
      # Get observed date range for this location
      loc_date_range <- range(curr_data$date)
      
      # Create new x-y pairs that include the observed range plus extrapolation points
      if (length(curr_data$date) < 4) {
        # For few points, simple linear extrapolation
        daily_values <- approx(x = as.numeric(curr_data$date), 
                               y = snow_smoothed, 
                               xout = as.numeric(dates_seq),
                               rule = 2)$y
      } else {
        # Use more robust methods for extrapolation with sufficient data
        tryCatch({
          # Try spline for smoother interpolation/extrapolation
          daily_values <- spline(x = as.numeric(curr_data$date), 
                                 y = snow_smoothed, 
                                 xout = as.numeric(dates_seq),
                                 method = "natural")$y
        }, error = function(e) {
          message("  Error in spline extrapolation for ", curr_combo, ": ", e$message)
          # Fall back to linear extrapolation
          daily_values <- approx(x = as.numeric(curr_data$date), 
                                 y = snow_smoothed, 
                                 xout = as.numeric(dates_seq),
                                 rule = 2)$y
        })
      }
    } else {
      # No extrapolation needed, just interpolate within observed range
      if (length(curr_data$date) < 4) {
        daily_values <- approx(x = as.numeric(curr_data$date), 
                               y = snow_smoothed, 
                               xout = as.numeric(dates_seq),
                               rule = 2)$y
      } else {
        # Use spline for smoother interpolation
        daily_values <- tryCatch({
          spline(x = as.numeric(curr_data$date), 
                 y = snow_smoothed, 
                 xout = as.numeric(dates_seq))$y
        }, error = function(e) {
          message("  Error in daily spline interpolation for ", curr_combo, ": ", e$message)
          # Fall back to linear interpolation
          approx(x = as.numeric(curr_data$date), 
                 y = snow_smoothed, 
                 xout = as.numeric(dates_seq),
                 rule = 2)$y
        })
      }
    }
    
    # 6. Ensure values are within valid range (0-100 for percentages)
    daily_values <- pmin(pmax(daily_values, 0), 100)
    
    # Create result dataframe with original and interpolated data
    # For original data points, create a flag to mark observed vs interpolated
    observed_dates <- curr_data$date
    
    # Create expanded result dataframe
    result_df <- data.frame(
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
      is_observed = dates_seq %in% observed_dates  # Flag for observed vs interpolated values
    )
    
    # Add year and month columns for easier analysis
    result_df <- result_df %>%
      mutate(
        year = year(date),
        month = month(date),
        day = day(date),
        doy = yday(date)  # day of year
      )
    
    # Store results
    all_results[[curr_combo]] <- result_df
  }
  
  # Combine all results
  if (length(all_results) == 0) {
    stop("No valid data could be processed for any plot-releve combination")
  }
  
  final_results <- bind_rows(all_results)
  
  message("Completed processing with ", 
          n_distinct(final_results$plot_id_releve), " plot-releve combinations and ",
          nrow(final_results), " total daily records")
  
  return(final_results)
}

# Helper function to calculate snow disappearance date
calculate_snow_disappearance <- function(snow, doy, month) {
  # Find winter peak (Jan-Mar)
  winter_idx <- which(month %in% c(1, 2, 3))
  if (length(winter_idx) == 0) return(NA)
  
  # Find peak snow in winter
  peak_idx <- winter_idx[which.max(snow[winter_idx])]
  if (length(peak_idx) == 0 || is.na(peak_idx)) return(NA)
  
  # Look for spring/summer snow disappearance
  # Only examine days between peak and day 200 (roughly mid-July)
  summer_cutoff <- which(doy >= 200)[1]
  if (is.na(summer_cutoff)) summer_cutoff <- length(snow)
  
  remaining_idx <- (peak_idx + 1):min(summer_cutoff, length(snow))
  if (length(remaining_idx) == 0) return(NA)
  
  # Find consecutive days with low snow cover
  low_snow_days <- which(snow[remaining_idx] <= 10)
  
  if (length(low_snow_days) == 0) return(NA)
  
  # Find first streak of at least 30 consecutive days with low snow
  consecutive_days <- 30
  streaks <- rle(diff(low_snow_days) == 1)
  streak_idx <- which(streaks$values & streaks$lengths >= (consecutive_days - 1))
  
  if (length(streak_idx) == 0) {
    # If no streak, find first day with low snow after peak
    snow_melt_idx <- remaining_idx[min(low_snow_days)]
  } else {
    # Find start of first qualifying streak
    start_pos <- ifelse(streak_idx[1] > 1, 
                        sum(streaks$lengths[1:(streak_idx[1]-1)]) + 1, 
                        1)
    snow_melt_idx <- remaining_idx[low_snow_days[start_pos]]
  }
  
  return(doy[snow_melt_idx])
}
# Function to calculate simplified snow metrics
calculate_snow_metrics <- function(data, stats_by_year = TRUE) {
  # Always group by plot-releve combination
  if (stats_by_year) {
    group_vars <- c("plot_id", "releve_id", "plot_id_releve", "year")
  } else {
    group_vars <- c("plot_id", "releve_id", "plot_id_releve")
  }
  
  # Calculate metrics
  snow_metrics <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      # Only the requested metrics
      snow_median = median(snow_smoothed, na.rm = TRUE),
      snow_sum = sum(snow_smoothed, na.rm = TRUE),
      days_with_snow = sum(snow_smoothed > 5, na.rm = TRUE),  # Days with >5% snow cover
      
      # Snow disappearance day
      snow_disappearance_doy = calculate_snow_disappearance(snow_smoothed, doy, month),
      
      # Time coverage
      start_date = min(date, na.rm = TRUE),
      end_date = max(date, na.rm = TRUE),
      days_covered = as.numeric(difftime(max(date, na.rm = TRUE), 
                                         min(date, na.rm = TRUE), 
                                         units = "days")) + 1,
      
      # Count observed vs total
      n_observed = sum(is_observed),
      n_total = n(),
      .groups = "drop"
    )
  
  return(snow_metrics)
}


process_temperature_data <- function(raw_temp_df, 
                                     spline_smoothing = 0.5,
                                     expand_time_series = TRUE,
                                     date_range = NULL,
                                     temp_col = "lst") {
  # Ensure library dependencies are loaded
  suppressPackageStartupMessages({
    require(zoo)
    require(lubridate)
  })
  
  # Add dplyr functions without loading the package (avoid conflicts)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }
  
  # Check inputs
  required_cols <- c("Date", "plot_id_releve", temp_col)
  missing_cols <- setdiff(required_cols, colnames(raw_temp_df))
  if (length(missing_cols) > 0) {
    stop("Input dataframe missing required columns: ", paste(missing_cols, collapse=", "))
  }
  
  # Ensure date is in proper format
  if (is.character(raw_temp_df$Date) || is.factor(raw_temp_df$Date)) {
    raw_temp_df$Date <- as.POSIXct(raw_temp_df$Date, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    if (all(is.na(raw_temp_df$Date))) {
      # Try another format
      raw_temp_df$Date <- as.POSIXct(raw_temp_df$Date, format="%Y-%m-%d", tz="UTC")
    }
  }
  
  # Create date without time component for daily aggregation
  raw_temp_df$date_only <- as.Date(raw_temp_df$Date)
  
  # Sort data
  raw_temp_df <- raw_temp_df[order(raw_temp_df$plot_id_releve, raw_temp_df$Date), ]
  
  # Define global date range for all time series
  if (expand_time_series) {
    if (!is.null(date_range)) {
      # Use user-provided date range
      global_date_range <- date_range
      message("Using user-provided date range: ", 
              format(global_date_range[1], "%Y-%m-%d"), " to ", 
              format(global_date_range[2], "%Y-%m-%d"))
    } else {
      # Use min/max dates from the entire dataset
      global_date_range <- range(raw_temp_df$date_only, na.rm = TRUE)
      message("Using global date range from data: ", 
              format(global_date_range[1], "%Y-%m-%d"), " to ", 
              format(global_date_range[2], "%Y-%m-%d"))
    }
  } else {
    global_date_range <- NULL
    message("Time series will not be expanded (using location-specific date ranges)")
  }
  
  # Get all unique plot_id_releve combinations
  plot_releve_combos <- unique(raw_temp_df$plot_id_releve)
  message("Processing ", length(plot_releve_combos), " unique plot-releve combinations...")
  
  # Initialize list to store results
  all_results <- list()
  
  # Process each plot_id_releve combination separately
  for (curr_combo in plot_releve_combos) {
    # Extract data for current combination
    curr_data <- raw_temp_df[raw_temp_df$plot_id_releve == curr_combo, ]
    
    # Skip if less than 3 valid observations (minimum needed for interpolation)
    if (nrow(curr_data) < 3 || sum(!is.na(curr_data[[temp_col]])) < 3) {
      message("Skipping ", curr_combo, " - insufficient data points")
      next
    }
    
    # Ensure data is sorted by date
    curr_data <- curr_data[order(curr_data$Date), ]
    
    message("Processing ", curr_combo, " with ", nrow(curr_data), " observations")
    
    # Consider using daily averages if multiple observations per day
    daily_data <- NULL
    try({
      daily_data <- tapply(curr_data[[temp_col]], curr_data$date_only, mean, na.rm=TRUE)
      daily_data <- data.frame(
        date = as.Date(names(daily_data)),
        temp = as.numeric(daily_data)
      )
      daily_data <- daily_data[!is.nan(daily_data$temp), ]  # Remove NaN values
    }, silent = TRUE)
    
    # If daily aggregation failed, use original data
    if (is.null(daily_data) || nrow(daily_data) < 3) {
      message("  Using original time points (daily aggregation failed)")
      daily_data <- data.frame(
        date = curr_data$date_only,
        temp = curr_data[[temp_col]]
      )
      daily_data <- daily_data[!is.na(daily_data$temp), ]  # Remove NA values
    }
    
    # Skip if still not enough valid observations after aggregation
    if (nrow(daily_data) < 3) {
      message("  Skipping ", curr_combo, " - insufficient data after aggregation")
      next
    }
    
    # Get local date range for this location
    loc_date_range <- range(daily_data$date, na.rm = TRUE)
    
    # Define the date sequence to use
    if (exists("global_date_range") && !is.null(global_date_range)) {
      expanded_range <- global_date_range
      needs_extrapolation <- TRUE
    } else {
      expanded_range <- loc_date_range
      needs_extrapolation <- FALSE
    }
    
    # Create sequence of daily dates
    dates_seq <- seq(expanded_range[1], expanded_range[2], by = "days")
    
    # Interpolate temperature for every day in the sequence
    interpolated_temps <- NULL
    
    # Flag any dates in the sequence that are in the original data
    observed_dates <- daily_data$date
    
    # Try different interpolation methods depending on data quantity/quality
    if (needs_extrapolation || length(dates_seq) > length(daily_data$date)) {
      # Need to extrapolate or interpolate missing days
      message("  Filling gaps in time series for ", curr_combo)
      
      if (nrow(daily_data) >= 5) {
        # For enough data points, use spline with smoothing
        try({
          interpolated_temps <- smooth.spline(
            x = as.numeric(daily_data$date),
            y = daily_data$temp,
            spar = spline_smoothing
          )
          
          # Predict for all dates in sequence
          interpolated_temps <- predict(
            interpolated_temps, 
            x = as.numeric(dates_seq)
          )$y
        }, silent = TRUE)
      }
      
      # Fall back to linear interpolation if spline failed
      if (is.null(interpolated_temps)) {
        message("  Using linear interpolation/extrapolation")
        interpolated_temps <- approx(
          x = as.numeric(daily_data$date),
          y = daily_data$temp,
          xout = as.numeric(dates_seq),
          rule = 2  # rule=2 means extrapolate if needed
        )$y
      }
    } else {
      # Simple case - just using existing data points
      interpolated_temps <- daily_data$temp
    }
    
    # Create result dataframe with daily temperatures
    result_df <- data.frame(
      date = dates_seq,
      temp_raw = approx(
        x = as.numeric(daily_data$date),
        y = daily_data$temp,
        xout = as.numeric(dates_seq),
        rule = 2
      )$y,
      temp_interpolated = interpolated_temps,
      plot_id_releve = curr_combo,
      is_observed = dates_seq %in% observed_dates  # Flag for observed vs interpolated values
    )
    
    # Extract lat/lon from original data
    if ("Latitude" %in% colnames(curr_data) && "Longitude" %in% colnames(curr_data)) {
      result_df$Latitude <- curr_data$Latitude[1]
      result_df$Longitude <- curr_data$Longitude[1]
    }
    
    # Add year and month columns for easier analysis
    result_df$year <- lubridate::year(result_df$date)
    result_df$month <- lubridate::month(result_df$date)
    result_df$doy <- lubridate::yday(result_df$date)  # day of year
    
    # Store results
    all_results[[curr_combo]] <- result_df
  }
  
  # Combine all results
  if (length(all_results) == 0) {
    stop("No valid data could be processed for any plot-releve combination")
  }
  
  final_results <- do.call(rbind, all_results)
  
  message("Completed processing with ", 
          length(all_results), " plot-releve combinations and ",
          nrow(final_results), " total daily records")
  
  return(final_results)
}

# Function to calculate the requested temperature metrics
calculate_temp_metrics <- function(data) {
  # Ensure required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for the metrics calculation")
  }
  
  # First calculate summer temperature for each plot_id_releve and year
  summer_temp <- data %>%
    dplyr::filter(month %in% 6:8) %>%  # June, July, August
    dplyr::group_by(plot_id_releve, year) %>%
    dplyr::summarise(
      summer_temp = mean(temp_interpolated, na.rm = TRUE),
      .groups = "drop"
    )  # Apply the function to each plot_id_releve and year
  warming_dates <- data %>%
    dplyr::group_by(plot_id_releve, year) %>%
    dplyr::summarise(
      spring_warming_doy = find_spring_warming_date(date, temp_interpolated, first(year)),
      .groups = "drop"
    )
  
  # Join the two metrics
  metrics <- dplyr::full_join(summer_temp, warming_dates, by = c("plot_id_releve", "year"))
  
  # Add Latitude/Longitude if available
  if (all(c("Latitude", "Longitude") %in% colnames(data))) {
    lat_lon <- data %>%
      dplyr::group_by(plot_id_releve) %>%
      dplyr::summarise(
        Latitude = dplyr::first(Latitude),
        Longitude = dplyr::first(Longitude),
        .groups = "drop"
      )
    
    metrics <- dplyr::left_join(metrics, lat_lon, by = "plot_id_releve")
  }
  
  return(metrics)
}
  
  # Function to find the first date temperature reaches 1°C before summer
find_spring_warming_date <- function(dates, temps, year) {
    # Get only data from March till summer to determine the growin season
    idx <- lubridate::month(dates) %in% 3:7 & lubridate::year(dates) == year
    if (sum(idx) == 0) return(NA)
    
    spring_dates <- dates[idx]
    spring_temps <- temps[idx]
    
    # Sort by date if not already sorted
    ord <- order(spring_dates)
    spring_dates <- spring_dates[ord]
    spring_temps <- spring_temps[ord]
    
    # Find first occurrence where temperature reaches 1°C
    idx_warm <- which(spring_temps >= 1)
    
    if (length(idx_warm) == 0) {
      # Temperature never reached 1°C in spring
      return(NA)
    } else {
      # Return the first date
      warming_date <- spring_dates[min(idx_warm)]
      # Return as day of year
      return(lubridate::yday(warming_date))
    }
  }
  
replace_with_closest_site_data <- function(data) {
  # Create a copy of the dataset
  filled_data <- data
  
  # Identify rows with NaN in et.annual
  na_rows <- which(is.na(data$et.annual))
  
  # Only proceed if there are NaN values to replace
  if(length(na_rows) > 0) {
    # Create spatial points for all sites
    all_sites <- data %>%
      filter(!is.na(et.annual)) %>%
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
    
    # For each row with NaN
    for(i in na_rows) {
      # Create point for the site with NaN
      target_site <- data[i, ] %>%
        st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
      
      # Calculate distances to all other sites with valid data
      distances <- st_distance(target_site, all_sites)
      
      # Find the index of the closest site
      closest_idx <- which.min(distances)
      
      # Replace NaN with value from closest site
      filled_data$et.annual[i] <- all_sites$et.annual[closest_idx]
    }
  }
  
  return(filled_data)
}

impute_environmental_data <- function(
    data,
    variables_to_impute = c("summer_temp", "fdd", "et.annual", "snow_sum"),
    m = 5,
    maxiter = 10,
    num.trees = 100,
    seed = 123,
    validation_fraction = 0.2
) {
  
  # Load required libraries
  require(tidyverse)
  require(miceRanger)
  require(sf)
  
  # Make a copy of the input data
  input_data <- data
  
  # 1. Add spatial features to enhance imputation
  # Create spatial object
  data_sf <- st_as_sf(input_data, coords = c("Longitude", "Latitude"), crs = 4326)
  
  # Transform to projected coordinate system for accurate distance calculation
  # UTM zone depends on the location - this assumes Alps region
  data_sf_projected <- st_transform(data_sf, 32632)
  
  # Calculate centroid of points
  centroid <- st_centroid(st_combine(data_sf_projected))
  
  # Calculate distance from each point to centroid
  input_data$dist_to_centroid <- st_distance(data_sf_projected, centroid) %>% as.numeric()
  
  # Calculate elevation-distance ratio (topographic position)
  input_data$elev_dist_ratio <- input_data$altitude / (input_data$dist_to_centroid + 1)
  
  # 2. Set up validation by masking a portion of the complete cases
  set.seed(seed)
  
  # Store original data for later comparison
  original_data <- input_data
  
  # Create validation set for each variable to impute
  validation_sets <- list()
  masked_indices <- list()
  
  for (var in variables_to_impute) {
    # Get indices of non-missing values
    complete_indices <- which(!is.na(input_data[[var]]))
    
    if (length(complete_indices) > 5) { # Only mask if we have enough complete cases
      # Randomly select indices to mask
      mask_size <- floor(length(complete_indices) * validation_fraction)
      mask_idx <- sample(complete_indices, mask_size)
      
      # Store original values for validation
      validation_sets[[var]] <- input_data[[var]][mask_idx]
      
      # Mask these values in the input data
      input_data[[var]][mask_idx] <- NA
      
      # Store masked indices for later
      masked_indices[[var]] <- mask_idx
    } else {
      validation_sets[[var]] <- NULL
      masked_indices[[var]] <- NULL
    }
  }
  
  # 3. Prepare data for miceRanger
  # Select variables to include in the imputation model
  numeric_cols <- names(input_data)[sapply(input_data, is.numeric)]
  
  # Remove ID column and keep only numeric variables plus our spatial additions
  mice_data <- input_data %>%
    select(all_of(numeric_cols))
  
  # 4. Run miceRanger
  set.seed(seed)
  
  mice_model <- miceRanger(
    mice_data,
    m = m,                      # Number of imputed datasets
    maxiter = maxiter,          # Number of iterations
    returnModels = TRUE,        # Save the ranger models
    verbose = FALSE,            # Don't show progress
    num.trees = num.trees,      # Trees per forest
    num.threads = parallel::detectCores() - 1  # Use parallel processing
  )
  
  # 5. Extract and combine imputed values
  completed_data <- completeData(mice_model)
  
  # Initialize results dataframe with original data
  result_data <- original_data
  
  # Add imputed values for each variable
  for (var in variables_to_impute) {
    # Extract imputed values from all datasets
    imputed_values <- sapply(completed_data, function(dataset) {
      dataset[[var]]
    })
    
    # Calculate mean imputation across all datasets
    imputed_mean <- rowMeans(imputed_values)
    
    # Create new column with imputed values
    result_data[[paste0(var, "_imputed")]] <- imputed_mean
    
    # Create final column with original values where available
    result_data[[paste0(var, "_final")]] <- ifelse(
      is.na(original_data[[var]]),
      imputed_mean,
      original_data[[var]]
    )
    
    # Add flag for imputed values
    result_data[[paste0(var, "_flag")]] <- ifelse(
      is.na(original_data[[var]]),
      "Imputed",
      "Original"
    )
  }
  
  # 6. Calculate performance metrics for validation sets
  performance <- data.frame(variable = character(),
                            r_squared = numeric(),
                            rmse = numeric(),
                            mae = numeric(),
                            stringsAsFactors = FALSE)
  
  for (var in variables_to_impute) {
    if (!is.null(validation_sets[[var]])) {
      # Get the masked indices for this variable
      mask_idx <- masked_indices[[var]]
      
      # Extract the actual values we stored earlier
      actual <- validation_sets[[var]]
      
      # Extract the predicted values for these indices
      predicted <- result_data[[paste0(var, "_imputed")]][mask_idx]
      
      # Calculate metrics
      r2 <- cor(actual, predicted)^2
      rmse <- sqrt(mean((actual - predicted)^2))
      mae <- mean(abs(actual - predicted))
      
      # Add to performance dataframe
      performance <- rbind(performance,
                           data.frame(variable = var,
                                      r_squared = r2,
                                      rmse = rmse,
                                      mae = mae))
    }
  }
  
  # 7. Return results
  return(list(
    imputed_data = result_data,
    performance = performance,
    mice_model = mice_model
  ))
}

impute_functional_traits <- function(
    data,
    variables_to_impute = c("LA", "LDMC", "N_percent", "SLA", "vegetative_height", "seed_mass"),
    m = 5,
    maxiter = 10,
    num.trees = 100,
    seed = 123,
    validation_fraction = 0.2
) {
  
  # Load required libraries
  require(tidyverse)
  require(miceRanger)
  
  # Make a copy of the input data
  input_data <- data
  
  # Store original data for later comparison
  original_data <- input_data
  
  # Create validation set for each variable to impute
  validation_sets <- list()
  masked_indices <- list()
  
  set.seed(seed)
  
  for (var in variables_to_impute) {
    # Get indices of non-missing values
    complete_indices <- which(!is.na(input_data[[var]]))
    
    if (length(complete_indices) > 5) { # Only mask if we have enough complete cases
      # Randomly select indices to mask
      mask_size <- floor(length(complete_indices) * validation_fraction)
      mask_idx <- sample(complete_indices, mask_size)
      
      # Store original values for validation
      validation_sets[[var]] <- input_data[[var]][mask_idx]
      
      # Mask these values in the input data
      input_data[[var]][mask_idx] <- NA
      
      # Store masked indices for later
      masked_indices[[var]] <- mask_idx
    } else {
      validation_sets[[var]] <- NULL
      masked_indices[[var]] <- NULL
    }
  }
  
  # Prepare data for miceRanger
  # Select variables to include in the imputation model
  numeric_cols <- names(input_data)[sapply(input_data, is.numeric)]
  
  # Remove ID column if any and keep only numeric variables
  mice_data <- input_data %>%
    select(all_of(numeric_cols))
  
  # Run miceRanger
  set.seed(seed)
  
  mice_model <- miceRanger(
    mice_data,
    m = m,                      # Number of imputed datasets
    maxiter = maxiter,          # Number of iterations
    returnModels = TRUE,        # Save the ranger models
    verbose = FALSE,            # Don't show progress
    num.trees = num.trees,      # Trees per forest
    num.threads = parallel::detectCores() - 1  # Use parallel processing
  )
  
  # Extract and combine imputed values
  completed_data <- completeData(mice_model)
  
  # Initialize results dataframe with original data
  result_data <- original_data
  
  # Add imputed values for each variable
  for (var in variables_to_impute) {
    # Extract imputed values from all datasets
    imputed_values <- sapply(completed_data, function(dataset) {
      dataset[[var]]
    })
    
    # Calculate mean imputation across all datasets
    imputed_mean <- rowMeans(imputed_values)
    
    # Create new column with imputed values
    result_data[[paste0(var, "_imputed")]] <- imputed_mean
    
    # Create final column with original values where available
    result_data[[paste0(var, "_final")]] <- ifelse(
      is.na(original_data[[var]]),
      imputed_mean,
      original_data[[var]]
    )
    
    # Add flag for imputed values
    result_data[[paste0(var, "_flag")]] <- ifelse(
      is.na(original_data[[var]]),
      "Imputed",
      "Original"
    )
  }
  
  # Evaluate the performance of the imputations
  performance <- data.frame(
    variable = character(),
    r_squared = numeric(),
    rmse = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in variables_to_impute) {
    if (!is.null(validation_sets[[var]])) {
      # Get predicted values for the masked indices
      predicted_values <- result_data[[paste0(var, "_imputed")]][masked_indices[[var]]]
      actual_values <- validation_sets[[var]]
      
      # Calculate R-squared
      ss_total <- sum((actual_values - mean(actual_values))^2)
      ss_residual <- sum((actual_values - predicted_values)^2)
      r_squared <- 1 - (ss_residual / ss_total)
      
      # Calculate RMSE
      rmse <- sqrt(mean((actual_values - predicted_values)^2))
      
      # Add to performance dataframe
      performance <- rbind(performance, data.frame(
        variable = var,
        r_squared = r_squared,
        rmse = rmse,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Return both the imputed data and performance metrics
  return(list(
    imputed_data = result_data,
    performance = performance
  ))
}

extract_mowing_events <- function(raster_pattern = "grassland-use_intensity_.*\\.tif$",
                                  coords_path = NULL,
                                  output_path = "mowing_events_at_veg_points.csv",
                                  raster_dir = ".",
                                  processed_raster_path = "processed_mowing_events.tif",
                                  force_reprocess = FALSE) {
  
  # Load required libraries
  require(tidyverse)
  require(terra)
  
  # Check if processed raster already exists and force_reprocess is FALSE
  if (file.exists(processed_raster_path) && !force_reprocess) {
    message("Using existing processed raster from: ", processed_raster_path)
    mowing_events_stack <- rast(processed_raster_path)
  } else {
    # Get all the land-use intensity rasters
    raster_files <- list.files(path = raster_dir, 
                               pattern = raster_pattern, 
                               full.names = TRUE)
    
    if (length(raster_files) == 0) {
      stop("No raster files found matching the pattern")
    }
    
    # Create a list to store the mowingEvents layers
    mowing_events_list <- list()
    
    # Process each file separately to avoid name conflicts
    for (i in seq_along(raster_files)) {
      # Load the raster
      r <- rast(raster_files[i])
      
      # Check if mowingEvents layer exists
      if ("mowingEvents" %in% names(r)) {
        # Extract just the mowingEvents layer
        mowing_layer <- r[["mowingEvents"]]
        
        # Add a unique name (using the year from the filename)
        year <- str_extract(raster_files[i], "\\d{4}")
        names(mowing_layer) <- paste0("mowingEvents_", year)
        
        # Store in list
        mowing_events_list[[i]] <- mowing_layer
      } else {
        warning(paste("mowingEvents layer not found in", raster_files[i]))
      }
    }
    
    if (length(mowing_events_list) == 0) {
      stop("No mowingEvents layers found in any of the raster files")
    }
    
    # Create a SpatRaster with all mowingEvents layers
    mowing_events_stack <- do.call(c, mowing_events_list)
    
    # Save the processed raster for future use
    writeRaster(mowing_events_stack, filename = processed_raster_path, overwrite = TRUE)
    message("Processed raster saved to: ", processed_raster_path)
  }
  
  # If coordinates are provided, extract values
  if (!is.null(coords_path)) {
    # Read the vegetation coordinates CSV file
    veg_coords <- read_csv(coords_path)
    
    # Convert to SpatVector
    veg_vector <- vect(veg_coords, geom = c("x", "y"), crs = "EPSG:4326")
    
    # Project to match the raster CRS if needed
    if (crs(veg_vector) != crs(mowing_events_stack)) {
      veg_vector <- project(veg_vector, crs(mowing_events_stack))
    }
    
    # Extract raster values at the coordinates
    extracted_values <- extract(mowing_events_stack, veg_vector)
    
    # Combine coordinates with extracted values
    result <- bind_cols(
      veg_coords,
      extracted_values %>% 
        select(-ID)  # Remove ID column from extraction if present
    )
    
    # Save the results if output path is provided
    if (!is.null(output_path)) {
      write_csv(result, output_path)
    }
    
    return(result)
  } else {
    return(mowing_events_stack)
  }
}


create_ternary <- function(data) {
  # Ensure data is properly scaled
  data <- data %>%
    mutate(
      # Scale the variables to ensure they sum to 1 for ternary plotting
      total = spa + codist + env,
      spa_scaled = spa/total,
      codist_scaled = codist/total,
      env_scaled = env/total,
      
      # Store RGB values for other plots
      R = spa_scaled,
      G = codist_scaled, 
      B = env_scaled,
      color = rgb(R, G, B)
    )
  
  # Create dense grid for background color
  n_points <- 300
  tern_grid <- crossing(
    x = seq(0, 1, length.out = n_points),
    y = seq(0, 1, length.out = n_points)
  ) %>%
    filter(x + y <= 1) %>%
    mutate(
      z = 1 - x - y,
      # Enforce strict bounds to handle floating point errors
      x = pmax(0, pmin(1, x)),
      y = pmax(0, pmin(1, y)),
      z = pmax(0, pmin(1, z)),
      color = rgb(x, y, z)
    )
  
  # Create the combined plot
  ggtern(data = tern_grid, aes(x = x, y = y, z = z)) +
    # Background color grid
    geom_point(aes(color = I(color)), size = 0.5, shape = 15) +
    # Add sites data points on top in black
    geom_point(data = data, 
               aes(x = spa_scaled, y = codist_scaled, z = env_scaled, fill = I(color)), 
               color = "black",
               shape = 21,
               size = 3,
               alpha = 0.7) +
    # Customize theme and labels
    theme_bw() +
    theme_showarrows() +
    labs(
      x = "Spatial",
      y = "Co-distribution",
      z = "Environment"
    ) +
    theme_custom() +
    theme(
      tern.axis.title = element_text(size = 12),
      tern.panel.grid.major = element_line(color = "gray90"),
      legend.position = "bottom"
    )
}

calculate_species_env_means <- function(species_data, env_data, site_id_col = "site_id") {
  # Step 1: Ensure inputs are tibbles/data frames and have the site_id column
  # Convert to tibble if not already
  species_data <- as_tibble(species_data, rownames = site_id_col)
  env_data <- as_tibble(env_data, rownames = site_id_col)
  
  # Step 2: Identify environmental variable columns (all except site_id)
  env_vars <- setdiff(colnames(env_data), site_id_col)
  
  # Step 3: Identify species columns (all except site_id)
  species_cols <- setdiff(colnames(species_data), site_id_col)
  
  # Step 4: Reshape the species data to long format
  species_long <- species_data %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "species",
      values_to = "presence"
    ) %>%
    filter(presence == 1)  # Keep only presence records
  
  # Step 5: Join with environmental data
  species_env <- species_long %>%
    left_join(env_data, by = site_id_col)
  
  # Step 6: Calculate the mean and standard deviation of each environmental variable per species
  species_means <- species_env %>%
    group_by(species) %>%
    summarise(across(
      all_of(env_vars),
      list(
        mean = ~ mean(.x, na.rm = TRUE)
      )
    ))
  
  # Return the results
  return(species_means)
}

calculate_community_traits <- function(community_data, traits_data, species_col = "species", abundance_col = NULL, trait_cols = NULL) {
  # Ensure inputs are tibbles, preserving row names as site IDs
 
    site_ids <- rownames(community_data)
    community_data <- as_tibble(community_data) %>%
      mutate(plot_id_releve = site_ids, .before = 1)
  
  # Ensure traits_data is a tibble
  traits_data <- as_tibble(traits_data)%>%
    select(-`...1`)
  
  # Identify trait columns (all columns except species column in traits_data)
  if(is.null(trait_cols)){
  trait_cols <- setdiff(colnames(traits_data), species_col)
  }
  
  # Identify species columns (all columns except site_id in community_data)
  species_cols <- setdiff(colnames(community_data), "plot_id_releve")
  
  # Reshape community data to long format
  community_long <- community_data %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "species",
      values_to = "value"
    ) %>%
    filter(value > 0)  # Keep only presence/abundance records
  
  # If abundance data is provided, use it for weighted means
  # Otherwise, use presence/absence (value = 1 for CWM calculation)
  if(is.null(abundance_col)) {
    community_long <- community_long %>%
      mutate(abundance = 1)
  } else {
    community_long <- community_long %>%
      rename(abundance = value)
  }
  
  # Join with traits data
  community_traits <- community_long %>%
    left_join(traits_data, by = c("species" = species_col))
  
  # Calculate community-weighted means and standard deviations for each trait by site
  community_weighted_means <- community_traits %>%
    group_by(plot_id_releve) %>%
    summarise(
      # Count number of species per site
      species_richness = n(),
      # Calculate weighted means for each trait
      across(
        all_of(trait_cols),
        list(
          cwm = ~ weighted.mean(.x, w = abundance, na.rm = TRUE)
        )
      )
    )
  
  return(community_weighted_means)
}

create_pca_biplot <- function(data, variables, scaling_factor = 5, point_size = 3, point_alpha = 0.6) {
  # Perform PCA
  res.pca <- PCA(data %>% select(all_of(variables)), scale.unit = TRUE)
  
  # Extract PCA coordinates
  ind_coords <- as.data.frame(get_pca_ind(res.pca)$coord)
  var_coords <- as.data.frame(get_pca_var(res.pca)$coord)
  
  # Add colors to ind_coords
  ind_coords$colors <- gsub('"', '', data$color)  # Remove quotes from hex codes
  
  # Create biplot
  p <- ggplot() +
    # Add arrows for variables
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, 
                     xend = Dim.1 * scaling_factor, 
                     yend = Dim.2 * scaling_factor),
                 arrow = arrow(length = unit(0.5, "cm")),
                 color = "black",
                 size = 0.5) +
    
    # Add variable labels
    geom_text_repel(data = var_coords,
                    aes(x = Dim.1 * scaling_factor,
                        y = Dim.2 * scaling_factor,
                        label = rownames(var_coords)),
                    color = "black",
                    max.overlaps = Inf) +
    
    # Add points with custom colors
    geom_point(data = ind_coords,
               aes(x = Dim.1, y = Dim.2, color = I(colors)),
               size = point_size,
               alpha = point_alpha) +
    
    # Add theme elements
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    # Add labels
    labs(title = "PCA - Biplot",
         x = paste0("Dim1 (", round(res.pca$eig[1,2], 1), "%)"),
         y = paste0("Dim2 (", round(res.pca$eig[2,2], 1), "%)")) +
    
    # Make it more square-shaped
    coord_fixed(ratio = 1) +
    
    # Expand the plot margins
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2))
  
  # Return list with the plot and PCA results for further analysis
  return(list(plot = p, pca = res.pca, ind_coords = ind_coords, var_coords = var_coords))
}
# ==============================================================================
# CUSTOM TERNARY PLOT FUNCTIONS FOR JSDM VARIANCE PARTITIONING
# ==============================================================================
# Adapted from sjSDM package internal functions
# Modified to color by species altitude ranges and community altitude

# Function to create abbreviated species name (first 4 letters of genus and species)
abbrev_species = function(sp_name) {
  paste0(substring(word(sp_name, 1),1,4), toupper(substring(word(sp_name, 2),1,1)), substring(word(sp_name, 2),2,4))
}


deg2rad = function(deg) {
  return(deg * pi / 180)
}

get_coords = function(X) {
  x_w = X[1] + sin(deg2rad(30)) * X[2]
  y_w = cos(deg2rad(30)) * X[2]
  return(c(x_w, y_w))
}

plot_line = function(X, Y, arrow = FALSE, ...) {
  Xc = get_coords(X)
  Yc = get_coords(Y)
  
  if(arrow) {
    arrows(Xc[1], Xc[2], Yc[1], Yc[2], ...)
  } else {
    segments(Xc[1], Xc[2], Yc[1], Yc[2], ...)
  }
}

plot_tern_base = function(data1, length = 0.12, col = "black", bg = NULL, cex = 1.0,
                          alpha = 0.7,
                          color_env = "#81caf3", color_codist = "#00bd89", color_spa = "#d00000") {
  # Draw triangle frame
  segments(0, 0, 1, 0)
  segments(0, 0.0, 0.5, 0.8660254)
  segments(1, 0, 0.5, 0.8660254)

  # Grid lines and labels
  for(i in seq(0.2, 0.8, length.out = 4)) {
    plot_line(c(1-i, 0.0, i), c(0.0, 1-i, i), col = "lightgrey")
    plot_line(c(1-i, i, 0.0), c(1-i, 0.0, i), col = "lightgrey")
    plot_line(c(0, i, 1-i), c(1-i, i, 0), col = "lightgrey")

    text(get_coords(c(1-i, 0.0, i))[1], get_coords(c(1-i, 0.0, i))[2]-0.03,
         labels = 1-i, srt = 60, xpd = NA)
    text(get_coords(c(1-i, i, 0.0))[1]+0.05, get_coords(c(1-i, i, 0.0))[2],
         labels = i, xpd = NA)
    text(get_coords(c(0, i, 1-i))[1]-0.03, get_coords(c(0, i, 1-i))[2]+0.03,
         labels = 1-i, xpd = NA, srt = -50)
  }

  # Axis labels with custom colors
  # Based on get_coords input order [spa, codist, env]:
  # Bottom-left (0,0) = env, Bottom-right (1,0) = spa, Top (0.5,0.866) = codist
  text(1, y = -0.02, pos = 1, xpd = NA, label = "Space", col = color_spa, font = 2)
  text(-0.1, y = -0.02, pos = 1, xpd = NA, label = "Environment", col = color_env, font = 2)
  text(0.5, y = 0.9, pos = 3, xpd = NA, label = "Species associations", col = color_codist, font = 2)

  # Normalize data
  data1[,1:3] = data1[,1:3]/rowSums(data1[,1:3])

  # Handle colors with alpha
  if(length(col) == 1) col = rep(col, nrow(data1))
  if(length(cex) == 1) cex = rep(cex, nrow(data1))

  # Apply alpha to colors
  if(!is.null(bg)) {
    if(length(bg) == 1) bg = rep(bg, nrow(data1))
    # Add alpha transparency to bg colors
    bg = adjustcolor(bg, alpha.f = alpha)
  } else {
    bg = adjustcolor(col, alpha.f = alpha)
  }

  # Set border color with alpha
  col = adjustcolor(col, alpha.f = alpha)

  # Plot points
  for(i in 1:nrow(data1)) {
    # Extract as [spa, codist, env] to match sjSDM coordinate system
    # where top=codist, bottom-left=env, bottom-right=spa
    coords = get_coords(unlist(data1[i, c(3, 2, 1)]))
    points(x = coords[1], y = coords[2], col = col[i], bg = bg[i], cex = cex[i], pch = 21)
  }
}

# Custom function for species ternary (colored by altitude range) - ggplot version
plot_tern_species = function(res, veg, veg.clim,
                             color_env = "#81caf3",
                             color_codist = "#00bd89",
                             color_spa = "#d00000",
                             color_palette = NULL,
                             alpha = 0.7,
                             cex = 1.5) {

  # Calculate species altitude ranges
  species_altitude_ranges = veg.clim %>%
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

  # Get species variance data
  species_data = res$internals$Species %>%
    rownames_to_column("species") %>%
    left_join(species_altitude_ranges, by = "species") %>%
    filter(!is.na(altitude_range))

  # Calculate proportions and coordinates
  species_data = species_data %>%
    mutate(
      total = env + codist + spa,
      env_prop = env / total,
      codist_prop = codist / total,
      spa_prop = spa / total
    )

  # Add x,y coordinates for plotting
  # IMPORTANT: Must match base R version which does data1[i, c(3,2,1)] on [env, codist, spa]
  # This gives [spa, codist, env] order for get_coords()
  species_data$x = NA
  species_data$y = NA
  for(i in 1:nrow(species_data)) {
    # Normalize first
    total = species_data$env[i] + species_data$codist[i] + species_data$spa[i]
    env_norm = species_data$env[i] / total
    codist_norm = species_data$codist[i] / total
    spa_norm = species_data$spa[i] / total
    # Pass as [spa, codist, env] to match sjSDM coordinate system
    coords = get_coords(c(spa_norm, codist_norm, env_norm))
    species_data$x[i] = coords[1]
    species_data$y[i] = coords[2]
  }

  # Set color palette
  if(is.null(color_palette)) {
    color_palette = colorRampPalette(c("gray90", "gray20"))
  }

  # Identify corner species for labeling
  species_data$dist_env = sqrt(species_data$x^2 + species_data$y^2)
  species_data$dist_codist = sqrt((species_data$x - 0.5)^2 + (species_data$y - 0.866)^2)
  species_data$dist_spa = sqrt((species_data$x - 1)^2 + species_data$y^2)

  # Find corner species
  corner_species = bind_rows(
    species_data %>% filter(env_prop > 0.8) %>% arrange(dist_env) %>% head(3),
    species_data %>% filter(codist_prop > 0.8) %>% arrange(dist_codist) %>% head(3),
    species_data %>% filter(spa_prop > 0.8) %>% arrange(dist_spa) %>% head(3)
  ) %>%
    distinct(species, .keep_all = TRUE)

  # Find center species (each proportion close to 1/3 = 0.33)
  # Center of ternary plot is where all three proportions are equal (≈0.33 each)
  center_species = species_data %>%
    filter(
      env_prop >= 0.25 & env_prop <= 0.40,
      codist_prop >= 0.25 & codist_prop <= 0.40,
      spa_prop >= 0.25 & spa_prop <= 0.40
    ) %>%
    mutate(dist_center = sqrt((x - 0.5)^2 + (y - 0.289)^2)) %>%  # Center is at (0.5, 0.289)
    arrange(dist_center) %>%
    head(3)

  # Combine all species to label
  labeled_species = bind_rows(corner_species, center_species) %>%
    distinct(species, .keep_all = TRUE) %>%
    mutate(label = abbrev_species(species))

  # Create ternary triangle outline with colored edges
  triangle_edges = data.frame(
    x_start = c(0, 1, 0.5),
    y_start = c(0, 0, 0.866),
    x_end = c(1, 0.5, 0),
    y_end = c(0, 0.866, 0),
    edge = c("spa", "codist", "env")
  )

  # Create grid lines and labels
  grid_lines = data.frame()
  grid_labels = data.frame()

  for(val in seq(0.2, 0.8, 0.2)) {
    # Horizontal lines (parallel to bottom, for codist axis)
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val/2, 1 - val/2),
      y = c(val * 0.866, val * 0.866),
      group = paste0("h", val),
      axis = "codist"
    ))
    # Label for codist axis (left side)
    label_pos = get_coords(c(1-val, 0.0, val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1],
      y = label_pos[2] - 0.03,
      label = as.character(1-val),
      angle = 60
    ))

    # Left diagonal lines (parallel to left edge, for spa axis)
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val, 0.5 + val/2),
      y = c(0, (1-val) * 0.866),
      group = paste0("l", val),
      axis = "spa"
    ))
    # Label for spa axis (bottom right)
    label_pos = get_coords(c(1-val, val, 0.0))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] + 0.05,
      y = label_pos[2],
      label = as.character(val),
      angle = 0
    ))

    # Right diagonal lines (parallel to right edge, for env axis)
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(1-val, 0.5 - val/2),
      y = c(0, (1-val) * 0.866),
      group = paste0("r", val),
      axis = "env"
    ))
    # Label for env axis (bottom left)
    label_pos = get_coords(c(0, val, 1-val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] - 0.03,
      y = label_pos[2] + 0.03,
      label = as.character(1-val),
      angle = -50
    ))
  }

  # Create plot
  p = ggplot() +
    # Triangle outline with colored edges
    geom_segment(data = triangle_edges,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = edge),
                 linewidth = 1.2) +
    scale_color_manual(values = c("env" = color_env, "codist" = color_codist, "spa" = color_spa),
                       guide = "none") +
    # Grid lines colored by axis they represent
    geom_line(data = grid_lines, aes(x = x, y = y, group = group, color = axis),
              linewidth = 0.5, alpha = 0.3) +
    # Grid labels
    geom_text(data = grid_labels, aes(x = x, y = y, label = label, angle = angle),
              size = 5, color = "gray30") +
    # Species points
    geom_point(data = species_data, aes(x = x, y = y, fill = mean_altitude),
               shape = 21, size = cex * 2, alpha = alpha, color = "black", stroke = 0.2) +
    # Labels for corner and center species
    geom_text_repel(data = labeled_species, aes(x = x, y = y, label = label),
                    size = 4, fontface = "italic",
                    force = 10,
                    min.segment.length = 0.01,
                    segment.color = "gray50") +
    # Color scale
    scale_fill_gradientn(colors = color_palette(100),
                         name = "Mean\naltitude (m)",
                         guide = guide_colorbar(barwidth = 0.8, barheight = 8)) +
    annotate("text", x = 0, y = 0, label = "Environment",
             color = color_env, fontface = "bold", hjust = 1.2, size = 4.5) +
    annotate("text", x = 0.5, y = 0.95, label = "Species associations",
             color = color_codist, fontface = "bold", hjust = 1, size = 4.5) +
    annotate("text", x = 1, y = 0.01, label = "Space",
             color = color_spa, fontface = "bold", hjust = -0.6, size = 4.5) +
    # Title
    ggtitle("Species") +
    # Theme
    coord_fixed(clip = "off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = -0.1, face = "bold", size = 12),
      legend.position = "right",
      plot.margin = margin(20, 10, 30, 40)
    )

  return(p)
}

# Custom function for sites/communities ternary (colored by altitude) - ggplot version
plot_tern_sites = function(res, veg.clim,
                           color_env = "#81caf3",
                           color_codist = "#00bd89",
                           color_spa = "#d00000",
                           color_palette = NULL,
                           alpha = 0.7,
                           cex = 1.5) {

  # Load required packages
  require(ggplot2)

  # Get sites variance data
  sites_data = res$internals$Sites %>%
    rownames_to_column("plot_id_releve") %>%
    left_join(veg.clim %>% select(plot_id_releve, altitude) %>% distinct(),
              by = "plot_id_releve") %>%
    filter(!is.na(altitude))

  # Add x,y coordinates for plotting
  # IMPORTANT: Must match base R version - normalize and pass as [spa, codist, env]
  sites_data$x = NA
  sites_data$y = NA
  for(i in 1:nrow(sites_data)) {
    # Normalize first
    total = sites_data$env[i] + sites_data$codist[i] + sites_data$spa[i]
    env_norm = sites_data$env[i] / total
    codist_norm = sites_data$codist[i] / total
    spa_norm = sites_data$spa[i] / total
    # Pass as [spa, codist, env] to match sjSDM coordinate system
    coords = get_coords(c(spa_norm, codist_norm, env_norm))
    sites_data$x[i] = coords[1]
    sites_data$y[i] = coords[2]
  }

  # Set color palette
  if(is.null(color_palette)) {
    color_palette = colorRampPalette(c("gray90", "gray20"))
  }

  # Create ternary triangle outline with colored edges
  triangle_edges = data.frame(
    x_start = c(0, 1, 0.5),
    y_start = c(0, 0, 0.866),
    x_end = c(1, 0.5, 0),
    y_end = c(0, 0.866, 0),
    edge = c("spa", "codist", "env")
  )

  # Create grid lines and labels
  grid_lines = data.frame()
  grid_labels = data.frame()

  for(val in seq(0.2, 0.8, 0.2)) {
    # Horizontal lines (parallel to bottom, for codist axis)
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val/2, 1 - val/2),
      y = c(val * 0.866, val * 0.866),
      group = paste0("h", val),
      axis = "codist"
    ))
    # Label for codist axis (left side)
    label_pos = get_coords(c(1-val, 0.0, val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1],
      y = label_pos[2] - 0.03,
      label = as.character(1-val),
      angle = 60
    ))

    # Left diagonal lines (parallel to left edge, for spa axis)
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(val, 0.5 + val/2),
      y = c(0, (1-val) * 0.866),
      group = paste0("l", val),
      axis = "spa"
    ))
    # Label for spa axis (bottom right)
    label_pos = get_coords(c(1-val, val, 0.0))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] + 0.05,
      y = label_pos[2],
      label = as.character(val),
      angle = 0
    ))

    # Right diagonal lines (parallel to right edge, for env axis)
    grid_lines = bind_rows(grid_lines, data.frame(
      x = c(1-val, 0.5 - val/2),
      y = c(0, (1-val) * 0.866),
      group = paste0("r", val),
      axis = "env"
    ))
    # Label for env axis (bottom left)
    label_pos = get_coords(c(0, val, 1-val))
    grid_labels = bind_rows(grid_labels, data.frame(
      x = label_pos[1] - 0.03,
      y = label_pos[2] + 0.03,
      label = as.character(1-val),
      angle = -50
    ))
  }

  # Create plot
  p = ggplot() +
    # Triangle outline with colored edges
    geom_segment(data = triangle_edges,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = edge),
                 linewidth = 1.2) +
    scale_color_manual(values = c("env" = color_env, "codist" = color_codist, "spa" = color_spa),
                       guide = "none") +
    # Grid lines colored by axis they represent
    geom_line(data = grid_lines, aes(x = x, y = y, group = group, color = axis),
              linewidth = 0.5, alpha = 0.3) +
    # Grid labels
    geom_text(data = grid_labels, aes(x = x, y = y, label = label, angle = angle),
              size = 5, color = "gray30") +
    # Sites points
    geom_point(data = sites_data, aes(x = x, y = y, fill = altitude),
               shape = 21, size = cex * 2, alpha = alpha, color = "black", stroke = 0.2) +
    # Color scale
    scale_fill_gradientn(colors = color_palette(100),
                         name = "Altitude\n(m)",
                         guide = guide_colorbar(barwidth = 0.8, barheight = 8)) +
    annotate("text", x = 0, y = 0, label = "Environment",
             color = color_env, fontface = "bold", hjust = 1.2, size = 4.5) +
    annotate("text", x = 0.5, y = 0.95, label = "Species associations",
             color = color_codist, fontface = "bold", hjust = 1, size = 4.5) +
    annotate("text", x = 1, y = 0.01, label = "Space",
             color = color_spa, fontface = "bold", hjust = -0.6, size = 4.5) +
    # Title
    ggtitle("Sites") +
    # Theme
    coord_fixed(clip = "off") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = -0.1, face = "bold", size = 12),
      legend.position = "right",
      plot.margin = margin(20, 10, 30, 40)
    )

  return(p)
}

plot.anova.custom = function(x,
                             y,
                             type = c( "McFadden", "Deviance", "Nagelkerke"),
                             fractions = c("discard", "proportional", "equal"),
                             cols = c("#81caf3","#00bd89","#d00000"),
                             alpha=0.15,
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
    if(x$spatial) {
      sapply(c("F_A", "F_B", "F_AB","F_S", "F_AS", "F_BS", "F_ABS"), function(i) which(values$models == i, arr.ind = TRUE))
    } else {
      sapply(c("F_A", "F_B", "F_AB"), function(i) which(values$models == i, arr.ind = TRUE))
    }
  values = values[select_rows,]
  col_index =
    switch (type,
            Deviance = 4,
            Nagelkerke = 5,
            McFadden = 6
    )
  graphics::plot(NULL, NULL, xlim = c(0,1), ylim =c(0,1),pty="s", axes = FALSE, xlab = "", ylab = "")
  xx = 1.1*lineSeq*cos( seq(0,2*pi, length.out=nseg))
  yy = 1.1*lineSeq*sin( seq(0,2*pi, length.out=nseg))
  graphics::polygon(xx+lineSeq,yy+(1-lineSeq), col= addA(cols[1],alpha = alpha), border = "black", lty = 1, lwd = 1)
  graphics::text(lineSeq-0.1, (1-lineSeq),labels = round(values[1,col_index],3))
  graphics::text(mean(xx+lineSeq), 0.9,labels = "Environmental", pos = 3)
  graphics::polygon(xx+1-lineSeq,yy+1-lineSeq, col= addA(cols[2],alpha = alpha), border = "black", lty = 1, lwd = 1)
  graphics::text(1-lineSeq+0.1, (1-lineSeq),labels = round(values[2,col_index],3))
  graphics::text(1-mean(xx+lineSeq), 0.9,labels = "Associations", pos = 3)
  graphics::text(0.5, (1-lineSeq),labels = round(values[3,col_index],3))
  if(x$spatial) {
    graphics::polygon(xx+0.5,yy+lineSeq, col= addA(cols[3],alpha = alpha), border = "black", lty = 1, lwd = 1)
    graphics::text(0.5, lineSeq+0.0,pos = 1,labels = round(values[4,col_index],3))
    graphics::text(0.5, 0.1,labels = "Spatial", pos = 1)
    graphics::text(0.3, 0.5,pos=1,labels   = round(values[5,col_index],3)) # AS
    graphics::text(1-0.3, 0.5,pos=1,labels = round(values[6,col_index],3)) # BS
    graphics::text(0.5, 0.5+0.05,labels    = round(values[7,col_index],3)) # ABS
  }
  out$VENN = values
  return(invisible(out))
}

addA = function(col, alpha = 0.25) apply(sapply(col, grDevices::col2rgb)/255, 2, function(x) grDevices::rgb(x[1], x[2], x[3], alpha=alpha))
