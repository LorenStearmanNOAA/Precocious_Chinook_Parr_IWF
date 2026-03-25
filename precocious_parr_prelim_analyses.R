#==============================================================================#
#   PRELIMINARY ANALYSES FOR CHINOOK SALMON PRECOCIOUS PARR FROM THE           #
#   IDAHO WILD FISH DATASET                                                    #
#==============================================================================#
#______________________________________________________________________________#
#                                                                              #

# Copyright 2025 U.S. Federal Government (in countries where
# recognized)

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software is
# distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied. See the License for the specific language governing
# permissions and limitations under the License.

# This script provides preliminary analyses of dynamics of Chinook salmon
# precocious parr across 16 streams in the Salmon River headwaters





#===============================================================================
#   ACQUISITIONS
#===============================================================================

# Libraries
    library(plotrix)
    library(sf)
    library(viridis)
    library(lubridate)
    library(dplyr)

# Acquires metadata about the sites and the sample area
    sites <- read.csv("meta/sites.csv", row.names = 1)
    sam_area <- read.csv("meta/sample_area_m2.csv", row.names = 1)
   
# Acquires the sample metadata 
  meta_raw <- read.csv("assembled/idaho_wild_fish_metadata_2026_02_06.csv")
  
# Acquires the assembled data
  data_raw <- read.csv("assembled/idaho_wild_fish_tag_data_2026_02_06.csv")
   
# Acquires the salmon river flowline shapefile
    salmon <- st_read(dsn = "spatial", layer = "Salmon_streams")

    


    
#===============================================================================
#   DATA SHAPING
#===============================================================================
    
# Merges the sample metadata with the raw data
    data_merge <- left_join(meta_raw[, c(24, 5, 7, 12)], data_raw)

# Converts the event date to an actual date
    data_merge$event_date <- parse_date_time(data_merge$event_date,
        order = c("YmdHMS", "mdYHM"))
    
# Adds a year metric
    data_merge$year <- year(data_merge$event_date)
    
# Extracts Chinook (spring/summer wild)
    data_merge <- data_merge[data_merge$srr_code %in% c("11W", "12W"), ]
    
# Extracts precocious parr
    prec_merge <- data_merge[grep("PR", data_merge$conditional_comments), ]

# Extracts regular parr
    reg_merge <- data_merge[-grep("PR", data_merge$conditional_comments), ]
    

# Removes years before 2001. This was the first year where observers were
# regularly recording precocious parr at multiple sites (n = 4 total fish
# before this year)
    prec_merge <- prec_merge[prec_merge$year > 2000, ]
    
    
    
    
#===============================================================================
#   PRECOCIOUS PARR BROAD ABUNDANCE PATTERNS
#===============================================================================

# Tallies par by site by year
    pp_y_s <- tapply(prec_merge$year, list(prec_merge$year,
      prec_merge$event_site), length)
    pp_y_s <- ifelse(is.na(pp_y_s), 0, pp_y_s)

# Calculates parr density
    pp_den <- pp_y_s / (sam_area/10000)

# Plots the density of precocious fish by site and by year
# Opens a new device
    dev.new(width = 13.33, height = 7.5)
    
# Sets the framing parameters
    par(mar = c(6, 4.5, 1, 1))
    par(mfrow = c(1, 2))
    
# Plots the first plot: precocious parr by year
    boxplot(t(pp_y_s) + 1, xaxt = "n", cex.axis = 1.15, lty = 1, cex.lab = 1.2,
        ylab = "Precocious fish / hectare")
    
# Adds the year labels
    axis(side = 1, at = 1:nrow(pp_y_s), labels = rownames(pp_y_s), las = 2,
        cex.axis = 0.8)
    
# Adds margin text
    mtext("Year", side = 1, line = 4, cex = 1.2)
    
# Plots the second plot: precocious parr by site
    boxplot(pp_y_s + 1, xaxt = "n", cex.axis = 1.15, lty = 1, ylab = "")
    
# Adds the site labels
    axis(side = 1, at = 1:ncol(pp_y_s), labels = colnames(pp_y_s),
        cex.axis = 0.8, las = 2)
    
# Adds the margin text
    mtext("Site", side = 1, line = 4, cex = 1.2)

# Plots spatial distribution of mean density
# Cuts the Salmon River to streams 3rd order and larger
    sal2 <- salmon[salmon$StreamOrde > 2, ]
    
# Opens a new device
    dev.new()
    
# Plots the streams
    plot(st_zm(st_geometry(sal2)), lwd = sal2$StreamOrde / 2,
        col = "dodgerblue2")
    
# Plots the mean density across years
    points(sites$Lon_m, sites$Lat_m, cex = rescale(apply(pp_y_s, 2, mean),
        c(1, 4)), bg = rgb(0, 0, 0, 0.4), pch = 21)


    


#===============================================================================
#   RECAPTURE PRELIM ANALYSES
#===============================================================================

# Identifies recaptures
    pp_recap <- prec_merge[grep("RE", prec_merge$conditional_comments), ]
    
# Iteratively figures out where these fish stacked up (length and condition)
# the year before
# Creates empty vectors to populate with values
    rel_prec_vec <- rel_no_dupl <- vector("numeric", nrow(pp_recap))
    rel_no_dupl <- site_init <- site_fin <- vector("character",
      length(rel_prec_vec))

# Iterative loop
    for (i in 1:length(rel_prec_vec)) {
        fish_i <- pp_recap[i, ]
        
    # Grabs the first occurrence of this individual
        first_i <- reg_merge[reg_merge$pit_tag == fish_i$pit_tag, ][1, ]
        
    # Stores the final site
        site_fin[i] <- fish_i$event_site
        
    # Next grabs the rest of that population
    # First checks to see if a prior capture exists
    # For 18//56 individuals, there isn't one. Querying PTAGIS found all 18
    # to have been tagged days prior at RSTs or fish traps
        if (!is.na(first_i$file_title)) {
            pop_i <- reg_merge[reg_merge$year == first_i$year &
                reg_merge$event_site == first_i$event_site, ]
            pop_i$len_z <- (pop_i$length - mean(pop_i$length, na.rm = TRUE)) /
                sd(pop_i$length, na.rm = TRUE)
            rel_prec_vec[[i]] <- pop_i[pop_i$pit_tag == fish_i$pit_tag, ]$len_z
            rel_no_dupl[i] <- ""
            site_init[i] <- pop_i$event_site[1]
        }
        
        else {
            rel_prec_vec[i] <- NA
            rel_no_dupl[i] <- fish_i$pit_tag
            site_init[i] <- NA
        }
        
    }
    
    
    
    
    
    
# Plots Z-score of size distribution
# Opens the plot
    dev.new(width = 5, height = 6)
    
# Sets the framing parameters
    par(mar = c(3, 5.5, 3, 0.5))
    
# Plots the boxplot
    boxplot(rel_prec_vec, boxwex = 0.7, staplewex = 0.2, lty = 1,
        xlab = "", ylab = "Z-score length relative\nto natal cohort", cex = 0,
        col = "white", cex.axis = 1.1, cex.lab = 1.15, staplewex = 0)
    
# Adds the axis label
    mtext("All recaptured\nprecocious parr (n = 38)", side = 1, line = 2,
      cex = 1.15)
    
# Adds a zero line
# Z-score of 0 = mean of the whole population
    abline(h = 0, lty = 3)
    
# Adds actual observation points
    points(rnorm(length(rel_prec_vec), 1, 0.05), rel_prec_vec, pch = 19,
        cex = 1.1, col = rgb(0, 0, 0, 0.3))


    

    
#______________________________________________________________________________#
#                                                                              #
#==============================================================================#
#   END OF SCRIPT                                                              #
#==============================================================================#