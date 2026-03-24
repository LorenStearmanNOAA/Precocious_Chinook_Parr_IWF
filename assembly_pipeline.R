#==============================================================================#
#   PIT TAGGING RAW DATA ASSEMBLY SCRIPT                                       #
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

# This script is used to assemble raw tagging files from the Idaho Wild Fish
# PIT Tag movement and survival study.





#===============================================================================
#   LIBRARIES AND DEFINITIONS                                                  
#===============================================================================

# Libraries
  library(lubridate)

# Tag file query: which files are tagging files?
    tag_files <- list.files("tagging_files", recursive = TRUE, full.names = TRUE)

# Approved colnames for metadata
    meta_colreplace <- read.csv("meta/metadata_colnames.csv")

# Acquires the site key
    site_key <- read.csv("meta/site_key.csv")
    
    
    
    
    

#===============================================================================
#   TEXT DATA ASSEMBLY
#===============================================================================
    
# Extracts the raw text files
    tag_txt <- tag_files[-grep("\\.csv", tag_files)]
    

# Stores output list objects
    meta_list <- data_list <- data_col_list <- vector("list", length(tag_txt))
    
# Iterative list to read in tag files
    for (i in 1:length(tag_txt)) {
    
    # Reads a test file
        file_i <- readLines(tag_txt[i])

    # Finds the metadata chunk
        meta_start <- grep("FILE TITLE", file_i)
        meta_stop <- grep("RELEASE RIVER KM", file_i)
        meta_foot1 <- grep("CLOSE DATE", file_i)
        meta_foot2 <- grep("CLOSE TIME", file_i)
        meta_lines <- c(meta_start:meta_stop, meta_foot1, meta_foot2)

    # Extracts the metadata chunk
        meta_chunk <- file_i[meta_lines]
        
    # Splits the metadata into columns
        meta_split <- lapply(meta_chunk, function(x) {
            out <- strsplit(x, " :")[[1]]
            if (length(out) == 1) {
                out <- c(out, NA)
            }
            out
        })
 
    # Removes spaces
        meta_split <- lapply(meta_split, function(x) gsub("(\\s+)\\s(\\s+)", "",
            x, perl = TRUE))
        
    # Grabs the metadata parameter values
        meta_vals <- lapply(meta_split, function(x) x[2])
        meta_vals <- do.call("cbind", meta_vals)
        meta_vals <- gsub("^\\s+", "", meta_vals)
        
    # Gets the column names for the metadata
        meta_names <- sapply(meta_split, function(x) x[1])
        
    # Applies column names to the metadata
        colnames(meta_vals) <- meta_names
     
    # Converts the meta values to a data frame
        meta_vals <- data.frame(meta_vals, check.names = FALSE)
        
        
    # Identifies which rows start with numbers. These should be data rows.
    # This is  done here because it's necessary for the comments in the
    # metadata to define where they should start and stop, but is used more
    # downstream with data extraction.
        number_rows <- as.numeric(unlist(lapply(file_i, function(x) {
            splits <- strsplit(x[[1]], " ")[[1]]
            splits <- splits[splits != ""]
            splits[1]
        })))

  
    # Checks to see if data exist. If not, this only saves the metadata
        if (sum(!is.na(number_rows)) != 0) {

        # Identifies the rows which contain data
            number_starts <- which.min(number_rows)
            number_rows[1:(number_starts - 1)] <- NA
            data_rows <- which(!is.na(number_rows))
            data_rows <- data_rows[data_rows > 20]
            
        # Extracts the data chunk
            data_chunk <- file_i[data_rows]
      
            
        # Extracts the record number
            rec_nums <- sapply(data_chunk, function(x) {
              as.numeric(substr(x, 1, 5))
            })
            
        # Extracts the tag code
            tag_codes <- sapply(data_chunk, function(x) {
            # Grabs the whole tag script
              raw_tag <- substr(x, 7, 20)
            # Trims whitespace
              raw_tag <- trimws(raw_tag)
            # If a mid space is present, pulls the text on the left side
            # Right side is checksum (only present on some older files)
              tag <- strsplit(raw_tag, " ")[[1]][1]
              
            })
            
        # Extracts the fork lengths
            forks <- sapply(data_chunk, function(x) {
              as.numeric(substr(x, 25, 28))
            })
            
        # Extracts the masses
            masses <- sapply(data_chunk, function(x) {
              as.numeric(substr(x, 32, 38))
            })
          
            
        # Extracts SRR
            srr <- sapply(data_chunk, function(x) {
              trimws(substr(x, 41, 43))
            })
            
        # Checks for absence of srr
            srr <- ifelse(srr == "" | nchar(srr) < 3,
              gsub(" ", "", paste0(meta_vals$SPECIES, meta_vals$RUN,
                meta_vals$'REARING TYPE')),
              srr)
            srr <- trimws(srr)
         
        # Extracts the comments
        # This section extracts text which is inferred to be comments (pos 45+),
        # and checks for/corrects pipe formatting. The assumption here is that
        # if there are no pipes and there is info, it's all conditonal comments.
            comments_raw <- sapply(data_chunk, function(x) {
              raw_vals <- substr(x, 45, nchar(x))
              trim_vals <- trimws(raw_vals)
              pipe_pos <- as.numeric(gregexpr("\\|", trim_vals)[[1]])
              pipe_count <- length(pipe_pos[!pipe_pos == -1])
              if (pipe_count == 0) {
                if (trim_vals == "" | is.null(raw_vals)) {
                  new_vals <- "||"
                } else {
                  new_vals <- paste0("|", trim_vals, "|")
                }
              }
              if (pipe_count == 1) {
                if (pipe_pos != 1) {
                  new_vals <- paste0("|", trim_vals)
                } else {
                  new_vals <- paste0(trim_vals, "|")
                }
              } else {
                if (pipe_count > 1) {
                  new_vals <- trim_vals
                }
              }
            # Removes any spaces immediately after the first pipe
              new_vals <- gsub("(?<=\\|) +", "", new_vals, perl = TRUE)
              
            # Checks to see if the last character is a pipe. If so, adds
            # an NA character
              pipe_last <- substr(new_vals, nchar(new_vals), nchar(new_vals)) ==
                "|"
              if (pipe_last == TRUE) {
                new_vals <- paste0(new_vals, "NA")
              }
              new_vals
              
          })
          

        # Splits the comments into conditional and textual comments
          comments_split <- lapply(comments_raw, function(x) {
            strsplit(x, "\\|")[[1]]
          })
          comments_split <- do.call("rbind", comments_split)
          comments_split <- ifelse(comments_split == "NA", "", comments_split)
          comments_split <- comments_split[, 2:3]
          
            
        # Creates the output data frame
          data_df <- data.frame('FILE TITLE' = meta_vals[1],
            'RECORD#' = rec_nums,
            'EVENT TYPE' = NA,
            'PIT TAG' = tag_codes,
            LENGTH = forks,
            WEIGHT = masses,
            'SRR CODE' = srr,
            'CONDITIONAL COMMENTS' = comments_split[, 1],
            'TEXT COMMENTS' = comments_split[, 2],
            check.names = FALSE)
          rownames(data_df) <- NULL
          
        # Saves the tabularized data chunk
          data_list[[i]] <- data_df
          
        }
        
        meta_list[[i]] <- meta_vals
        
        
    }

# Formats and combines the metadata
# First determines all possible columns. For simplicity's sake, the order that
# R finds these columns will be the default order for a bit.
    meta_col_all <- unique(unlist(sapply(meta_list, colnames)))
    
    meta_col_ord_mat <- sapply(meta_list, function(x) {
        col_ord <- sapply(meta_col_all, function(y) which(colnames(x) %in% y))
        col_ord <- sapply(col_ord, function(x) ifelse(is.null(x), NA, x))
        col_ord/max(col_ord, na.rm = TRUE)
    })
    meta_col_ord <- apply(meta_col_ord_mat, 1, mean, na.rm = TRUE)
    meta_col_all <- names(sort(meta_col_ord))
    
    
# Iterative column adding and reordering for each file
    txt_meta_list <- vector("list", length(meta_list))
    for (i in 1:length(txt_meta_list)) {
        meta_i <- meta_list[[i]]
        present_cols <- colnames(meta_i)
        missing_cols <- meta_col_all[!meta_col_all %in% present_cols]
        new_cols <- matrix(NA, 1, length(missing_cols))
        new_cols <- data.frame(new_cols)
        colnames(new_cols) <- missing_cols
        new_df <- cbind(meta_i, new_cols)
        new_df <- new_df[, match(meta_col_all, colnames(new_df))]
        txt_meta_list[[i]] <- new_df
    }
    
    txt_meta <- do.call("rbind", txt_meta_list)

# Combines the data 
    txt_data <- do.call("rbind", data_list)

    
    

        
#===============================================================================
#   CSV DATA ASSEMBLY
#===============================================================================
    
# Extracts the raw csv files
    tag_csv <- tag_files[grep("\\.csv", tag_files)]

# Reads the csv files
    csv_list <- lapply(tag_csv, read.csv, check.names = FALSE)
    for (i in 1:length(csv_list)) {
        csv_list[[i]]$'FILE TITLE' <- gsub(".*(?=\\/)\\/", "", tag_csv[i],
            perl = TRUE)
    }

# Determines whether multiple date/time values are present in each file
# If so, sets the date and time values to the first observation
    csv_list <- lapply(csv_list, function(x) {
        x_evdt_n <- table(x$'Event Date')
        if (length(x_evdt_n) > 1) {
          x$'Event Date' <- x$'Event Date'[1]
        }
        x_evdt_p_n <- table(x$'Event Date (PST)')
        if (length(x_evdt_p_n) > 1) {
          x$'Event Date (PST)' <- x$'Event Date (PST)'[1]
        }
        x
      
    })
    
# Combines the data files into a single object
    csv_all <- do.call("rbind", csv_list)

# Capitalizes the column names    
    colnames(csv_all) <- toupper(colnames(csv_all))

# Splits the data and metadata
# First identifies columns with no data
    na_props <- apply(csv_all, 2, function(x) sum(is.na(x))/length(x))
    nodata_cols <- names(which(na_props == 1))
    
# Next identifies data columns
    data_cols <- c(
      "FILE TITLE",
      "RECORD#",
      "EVENT TYPE",
      "PIT TAG",
      "LENGTH",
      "WEIGHT",
      "SRR CODE",
      "CONDITIONAL COMMENTS",
      "TEXT COMMENTS"
    )
    
# Last uses nodata and data columns to find metadata columns
    meta_cols <- colnames(csv_all)
    meta_cols <- meta_cols[!meta_cols %in% c(data_cols, nodata_cols)]
    meta_cols <- c(meta_cols, 'FILE TITLE')
    
# Splits the data and metadata
    csv_data <- csv_all[, colnames(csv_all) %in% data_cols]
    csv_data <- csv_data[, match(data_cols, colnames(csv_data))]
    csv_meta <- csv_all[, colnames(csv_all) %in% meta_cols]
    csv_meta <- unique(csv_meta)
    
    
    
    
    
#===============================================================================
#   METADATA ASSEMBLY AND FORMATTING
#===============================================================================

# The general philosophy here is to try to match data from multiple eras
# which may be formatted in multiple fashions over time and stored in varying
# numbers of column per era. Date/time is a great example of this.
    
# This first bit handles columns which were combined in the txt data after
# 1996
    
# Columns which can be dropped after this:
# TAG TIME
# RELEASE TIME
# CLOSE TIME
    
# TXT data
# Creates a unified tagging date/time
    new_TAG_DATE <- ifelse(!is.na(txt_meta$'TAG TIME'),
      paste(txt_meta$'TAG DATE', txt_meta$'TAG TIME'), txt_meta$'TAG DATE')
    
# Creates a unified creation date/time
# This is the same as TAG DATE
    new_CREATION_DATE <- ifelse(!is.na(txt_meta$'CREATION TIME'),
      paste(txt_meta$'CREATION DATE', txt_meta$'CREATION TIME'),
      txt_meta$'CREATION DATE')    
    
# Replaces creation date with tag date
    new_TAG_DATE <- ifelse(is.na(new_TAG_DATE), new_CREATION_DATE, new_TAG_DATE)
    
    
# Creates a unified release date/time
    new_RELEASE_DATE <- ifelse(!is.na(txt_meta$'RELEASE TIME'),
      paste(txt_meta$'RELEASE DATE', txt_meta$'RELEASE TIME'),
      txt_meta$'RELEASE DATE')
    
# Creates a unified close date/time
    new_CLOSE_DATE <- ifelse(!is.na(txt_meta$'CLOSE TIME'),
      paste(txt_meta$'CLOSE DATE', txt_meta$'CLOSE TIME'),
      txt_meta$'CLOSE DATE')   
    
# Combines release water temperature columns
    new_RELEASE_WATER_TEMP <- ifelse(is.na(txt_meta$'RELEASE TEMP'),
      txt_meta$'RELEASE WATER TEMP', txt_meta$'RELEASE TEMP')
    
# Combines release site and release location
    release_site_meta <- data.frame(OLD = c("BEAR VALLEY", "BEAR VALLEY CR.",
      "BIG CREEK", "CAPE HORN CR.", "EAST FORK SALMON", "EAST FORK SALMON R.",
      "ELK CR", "MARSH CR", "SECESH R", "SOUTH FORK SALMON R.",
      "SULFUR", "VALLEY CR"),
      NEW  = c("BEARVC", "BEARVC", "BIGC", "CAPEHC", "SALREF", "SALREF", "ELKC",
        "MARSHC", "SECESR", "SALRSF", "SULFUC", "VALEYC"))
    
    for (i in 1:nrow(txt_meta)) {
      if (txt_meta$'RELEASE LOCATION'[i] %in% release_site_meta[, 1]) {
        site_i <- which(release_site_meta[, 1] == txt_meta$'RELEASE LOCATION'[i])
        txt_meta$'RELEASE LOCATION'[i] <- release_site_meta[site_i, 2]
      }
    }
    
    new_RELEASE_SITE <- ifelse(is.na(txt_meta$'RELEASE SITE'),
      txt_meta$'RELEASE LOCATION', txt_meta$'RELEASE SITE')
    
# Combines agency and organization
    new_ORGANIZATION <- ifelse(is.na(txt_meta$'ORGANIZATION'),
      txt_meta$'AGENCY', txt_meta$'ORGANIZATION')
  
# Creates RKM MASK and RKM EXT objects
    period_pos <- sapply(gregexpr("\\.", txt_meta$'RELEASE RIVER KM'),
      tail, 1)
    RKM_MASK <- substr(txt_meta$'RELEASE RIVER KM', 1, period_pos - 1)
    RKM_EXT <- substr(txt_meta$'RELEASE RIVER KM', period_pos + 1,
      nchar(txt_meta$'RELEASE RIVER KM'))
    
    
# Updates the txt metadata
    drop_cols <- c("TAG TIME", "CREATION TIME", "RELEASE TIME", "CLOSE TIME",
      "CREATION DATE", "RELEASE LOCATION", "RELEASE WATER TEMP", "AGENCY",
      "RELEASE RIVER KM")
    txt_meta_2 <- txt_meta[, -which(colnames(txt_meta) %in% drop_cols), ]
    txt_meta_2$'TAG DATE' <- new_TAG_DATE
    txt_meta_2$'RELEASE DATE' <- new_RELEASE_DATE
    txt_meta_2$'CLOSE DATE' <- new_CLOSE_DATE
    txt_meta_2$'RELEASE SITE' <- new_RELEASE_SITE
    txt_meta_2$ORGANIZATION <- new_ORGANIZATION
    txt_meta_2$'RKM MASK' <- RKM_MASK
    txt_meta_2$'RKM EXT' <- RKM_EXT
    
    
# Fixes the column names
    for (i in 1:ncol(txt_meta_2)) {
      colnames(txt_meta_2)[i] <- meta_colreplace[meta_colreplace[, 1] ==
          colnames(txt_meta_2)[i], 2]
    }
    
# Adds empty columns to the txt data for columns present in csv but not txt
    txt_missings <- colnames(csv_meta)[!colnames(csv_meta) %in% colnames(txt_meta_2)]
    txt_empties <- matrix(NA, nrow(txt_meta_2), length(txt_missings))
    colnames(txt_empties) <- txt_missings
    txt_empties <- data.frame(txt_empties, check.names = FALSE)
    txt_meta_3 <- cbind(txt_meta_2, txt_empties)
    
# Does the same but for the csv metadata
    csv_missings <- colnames(txt_meta_2)[!colnames(txt_meta_2) %in% colnames(csv_meta)]
    csv_empties <- matrix(NA, nrow(csv_meta), length(csv_missings))
    colnames(csv_empties) <- csv_missings
    csv_empties <- data.frame(csv_empties, check.names = FALSE)
    csv_meta_2 <- cbind(csv_meta, csv_empties)
    
# Reorders the txt columns
    txt_meta_3 <- txt_meta_3[, match(colnames(csv_meta_2), colnames(txt_meta_3))]

    
    
    
    
#===============================================================================
#   DATA ASSEMBLY AND FORMATTING
#===============================================================================
    
# Combines the two metadata files
    meta_raw <- rbind(txt_meta_3, csv_meta_2)

# Stores an ordering vector to correct dates to a single format
    ords <- c("mdyHM", "mdYHMS", "mdYIMSp", "mdY")

# Finds date columns
    dt_cols <- grep("DATE", colnames(meta_raw))

# Corrects dates
    meta_raw[, dt_cols] <- lapply(meta_raw[, dt_cols], parse_date_time,
        order = ords)
    
# Corrects sites
# First identifies the site columns (minus HATCHERY.SITE)
    site_cols <- grep("SITE", colnames(meta_raw))
    hatch <- grep("HATCHERY", colnames(meta_raw))
    site_cols <- site_cols[!site_cols %in% hatch]
    
# Next matches the various input names to the output corrected names
# Note that event site and release site aren't always the same - sometimes when
# two sites are extremely close together (like a tributary confluence), fish may
# have been tagged from one and released in the other.
    meta_raw[, site_cols] <- lapply(meta_raw[, site_cols], function(x) {
        sapply(x, function(SITE) {
            site_key[site_key[, 1] == SITE, 2]
        })
    })

# Fixes temperature values
# Finds temperature columns
    temp_cols <- grep("TEMP", colnames(meta_raw))
    meta_raw[, temp_cols] <- lapply(meta_raw[, temp_cols], function(x) {
        as.numeric(gsub("[^0-9\\.\\s*]", "", x))
    })

# Fixes column names
# Names are transferred to lower case, and periods are replaced with underscores
    colnames(meta_raw) <- gsub("\\.|\\s", "_", tolower(colnames(meta_raw)))

# Merges the actual data into a single object
# Data are standardized enough to do this in a single step
    data_raw <- rbind(txt_data, csv_data)
      
# Fixes the column names
    colnames(data_raw) <- gsub("\\.|\\s", "_", tolower(colnames(data_raw)))

# Matches the file title to the metadata by removing leading directory info
#    data_raw$file_title <- gsub(".*(?=\\/)\\/", "", data_raw$file_title,
#      perl = TRUE)
    
# # Removes 3, 3.1, "PM", and "T" from conditional comments
# # Check against PTAGIS showed these were not real comments
# # T in particular seems associated with a line ending corruption in two files
#     data_raw$conditional_comments <- gsub("3|3.1|^T$|PM", "",
#       data_raw$conditional_comments)
#     
# # Splits conditional comments into unique elements
#     cond_split <- lapply(data_raw$conditional_comments, function(x) {
#         strsplit(x, " ")[[1]]
#     })
#     
# # Stores unique comments
#     unique_condcomm <- cond_split %>% unlist %>% unique %>% sort
#     
#     cond_exp <- sapply(cond_split, function(x) {
#       out <- tapply(x, factor(x, levels = unique_condcomm), length)
#       ifelse(is.na(out), 0, 1)
#     })
      
    
    

      
    
#===============================================================================
#   FILE OUTPUTS
#===============================================================================
    
    
# Stores output file names
    meta_fn <- paste0("assembled/idaho_wild_fish_metadata_", Sys.Date(),
      ".csv")
    meta_fn <- gsub("-", "_", meta_fn)
    data_fn <- paste0("assembled/idaho_wild_fish_tag_data_", Sys.Date(),
      ".csv")
    data_fn <- gsub("-", "_", data_fn)
     
# Writes output files
    write.csv(meta_raw, meta_fn, row.names = FALSE, na = "NA")
    write.csv(data_raw, data_fn, row.names = FALSE, na = "NA")

    
    
    
    
#______________________________________________________________________________#
#                                                                              #
#==============================================================================#
#   END OF SCRIPT                                                              #
#==============================================================================#
      