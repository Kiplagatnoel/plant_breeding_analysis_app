# ============================================================================
# SECTION 1: PACKAGE LOADING AND INITIALIZATION
# ============================================================================

# List of required packages
required_packages <- c(
  "shiny", "shinydashboard", "DT", "plotly", "dartR", "ASRgenomics",
  "pheatmap", "ggplot2", "dplyr", "tidyr", "lme4", "metan", "ggrepel",
  "corrplot", "factoextra", "adegenet", "readxl", "reshape2", "rmarkdown",
  "qqman", "purrr", "CMplot", "shinyWidgets", "shinycssloaders", "knitr","poppr", "hierfstat",
  "foreach", "doParallel", "matrixStats", "stringr", "viridis", "RColorBrewer","visNetwork","qtl","openxlsx","LEA"
)

# Install missing packages
install_missing_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE, repos = "https://cloud.r-project.org")
  }
  # Load all packages
  suppressPackageStartupMessages({
    lapply(packages, library, character.only = TRUE)
  })
}

# Install and load packages
suppressWarnings({
  tryCatch({
    install_missing_packages(required_packages)
  }, error = function(e) {
    cat("Some packages failed to install:", e$message, "\n")
  })
})

#=============================================================================
#Data processing function
#=============================================================================

# Load and preprocess data with automatic detection
preprocess_data <- function(data) {
  # Identify numeric columns
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  
  # Identify non-numeric columns to convert to factors
  non_numeric_cols <- setdiff(names(data), numeric_cols)
  
  # Convert non-numeric columns to factors
  if (length(non_numeric_cols) > 0) {
    data <- data %>%
      mutate(across(all_of(non_numeric_cols), as.factor))
  }
  
  return(list(
    data = data,
    numeric_cols = numeric_cols,
    factor_cols = non_numeric_cols
  ))
}


# ============================================================================
# SECTION 2: DATA LOADING AND PROCESSING FUNCTIONS
# ============================================================================

#' Create genlight object from DArT data
create_genlight <- function(filename, ind.metafile = NULL) {
  tryCatch({
    cat("Reading DArT data from:", filename, "\n")
    
    # Check if file exists
    if (!file.exists(filename)) {
      stop("File not found:", filename)
    }
    
    # Try to read the data based on file extension
    if (grepl("\\.xlsx?$", filename, ignore.case = TRUE)) {
      # For Excel files
      library(readxl)
      data <- read_excel(filename)
      cat("Read Excel file with dimensions:", dim(data), "\n")
      
      # Extract sample names (usually first column)
      if (ncol(data) < 2) {
        stop("Excel file has insufficient columns for genotype data")
      }
      
      sample_names <- as.character(data[[1]])
      
      # Find where genotype data starts (skip metadata columns)
      numeric_cols <- which(sapply(data, is.numeric))
      if (length(numeric_cols) == 0) {
        stop("No numeric genotype columns found in Excel file")
      }
      
      start_col <- min(numeric_cols)
      cat("Genotype data starts at column:", start_col, "\n")
      
      # Extract genotype matrix
      geno_matrix <- as.matrix(data[, start_col:ncol(data)])
      rownames(geno_matrix) <- sample_names
      
      # Create genlight object
      gl_object <- new("genlight", 
                       as.matrix(geno_matrix),
                       ploidy = rep(2, nrow(geno_matrix)),
                       ind.names = sample_names,
                       loc.names = colnames(geno_matrix))
      
      # Add basic metadata
      gl_object@other$ind.metrics <- data.frame(ID = sample_names)
      
      return(gl_object)
      
    } else {
      # For CSV/text files - use dartR
      cat("Reading CSV/text file with dartR\n")
      
      if (is.null(ind.metafile)) {
        gl_data <- dartR::gl.read.dart(filename, verbose = 0)
      } else {
        if (file.exists(ind.metafile)) {
          gl_data <- dartR::gl.read.dart(filename, ind.metafile, verbose = 0)
        } else {
          gl_data <- dartR::gl.read.dart(filename, verbose = 0)
        }
      }
      
      # Basic validation
      if (nInd(gl_data) == 0 || nLoc(gl_data) == 0) {
        stop("Invalid genotype data: no individuals or loci found")
      }
      
      cat("Successfully created genlight object\n")
      cat("  Individuals:", nInd(gl_data), "\n")
      cat("  Loci:", nLoc(gl_data), "\n")
      cat("  Ploidy:", gl_data@ploidy, "\n")
      
      return(gl_data)
    }
  }, error = function(e) {
    cat("Error creating genlight object:", e$message, "\n")
    return(NULL)
  })
}

#' Get statistics from genlight object
genlight_to_matrix_stats <- function(gl_object) {
  if (is.null(gl_object)) {
    return(NULL)
  }
  
  tryCatch({
    # Convert to genotype matrix (0,1,2 format)
    geno_matrix <- as.matrix(gl_object)
    
    # Set row names as individual names
    if (is.null(rownames(geno_matrix))) {
      rownames(geno_matrix) <- indNames(gl_object)
    }
    
    cat("Genotype data loaded successfully\n")
    cat("  Number of genotypes:", nrow(geno_matrix), "\n")
    cat("  Number of markers:", ncol(geno_matrix), "\n")
    
    # Calculate missing data rate
    missing_rate <- round(sum(is.na(geno_matrix)) / (nrow(geno_matrix) * ncol(geno_matrix)) * 100, 2)
    cat("  Missing data rate:", missing_rate, "%\n")
    
    # Calculate MAF distribution
    if (ncol(geno_matrix) > 0) {
      maf <- colMeans(geno_matrix, na.rm = TRUE) / 2
      maf <- pmin(maf, 1 - maf)
      cat("  Mean MAF:", round(mean(maf, na.rm = TRUE), 4), "\n")
    }
    
    return(geno_matrix)
  }, error = function(e) {
    cat("Error getting genotype stats:", e$message, "\n")
    return(NULL)
  })
}



# Main wrapper function to clean chromosome information in a genlight object
# Main wrapper function to clean and order chromosome information in a genlight object
clean_genlight_chromosomes <- function(gl, verbose = TRUE) {
  # Store original dimensions for reporting
  original_n_snps <- nLoc(gl)
  original_n_inds <- nInd(gl)
  
  # Step 1: Find chromosome information
  chrom_info <- find_chromosome_info(gl)
  
  if (is.null(chrom_info$chrom_vec)) {
    stop("Cannot find chromosome information in the genlight object")
  }
  
  if (verbose) {
    cat("=== Chromosome Cleaning Report ===\n")
    cat("Original SNPs:", original_n_snps, "\n")
    cat("Original individuals:", original_n_inds, "\n")
  }
  
  # Step 2: Clean the chromosome vector
  cleaning_results <- clean_chromosome_vector_extended(chrom_info$chrom_vec, verbose)
  
  # NEW: Step 3: Sort SNPs by chromosome number (and optionally by position)
  if (verbose) {
    cat("\n=== Ordering SNPs by chromosome ===\n")
  }
  
  # Get the indices of kept SNPs
  kept_indices <- which(cleaning_results$keep)
  
  # Get chromosome numbers for kept SNPs
  chrom_numbers <- cleaning_results$cleaned_numeric
  
  # Get position information if available (for ordering within chromosomes)
  pos_info <- get_position_info(gl)
  if (!is.null(pos_info)) {
    pos_kept <- pos_info[kept_indices]
    if (verbose) {
      cat("Position information found - ordering by chromosome then position\n")
    }
    # Order by chromosome number, then by position
    order_within_chrom <- order(chrom_numbers, pos_kept)
  } else {
    if (verbose) {
      cat("No position information found - ordering by chromosome only\n")
    }
    # Order only by chromosome number
    order_within_chrom <- order(chrom_numbers)
  }
  
  # Apply ordering to kept indices
  ordered_kept_indices <- kept_indices[order_within_chrom]
  ordered_chrom_numbers <- chrom_numbers[order_within_chrom]
  
  # Create new keep vector with ordered indices
  ordered_keep <- rep(FALSE, length(chrom_info$chrom_vec))
  ordered_keep[ordered_kept_indices] <- TRUE
  
  # Create new cleaned_numeric vector in the correct order
  ordered_cleaned_numeric <- rep(NA, length(chrom_info$chrom_vec))
  ordered_cleaned_numeric[ordered_kept_indices] <- ordered_chrom_numbers
  
  if (verbose) {
    cat("Ordered SNPs by chromosome numbers\n")
    cat("Chromosome distribution after ordering:\n")
    print(table(ordered_chrom_numbers, useNA = "ifany"))
  }
  
  # Step 4: Subset the genlight object based on ordered kept SNPs
  if (sum(ordered_keep) == 0) {
    warning("No SNPs remain after cleaning! Returning original object.")
    return(gl)
  }
  
  # Subset the genlight object using ordered indices
  gl_clean <- gl[, ordered_kept_indices]
  
  # Step 5: Update chromosome information in the cleaned object
  if (!is.null(chrom_info$location_name)) {
    if (chrom_info$location_name == "other$chromosome") {
      gl_clean@other$chromosome <- ordered_chrom_numbers
    } else if (chrom_info$location_name == "chromosome") {
      gl_clean@chromosome <- ordered_chrom_numbers
    } else if (chrom_info$location_name == "other$loc.metrics$chromosome") {
      gl_clean@other$loc.metrics$chromosome <- ordered_chrom_numbers
    } else if (chrom_info$location_name == "other$loc.metrics$chrom_col") {
      # Update the specific column in loc.metrics
      gl_clean@other$loc.metrics[[chrom_info$column_name]] <- ordered_chrom_numbers
    }
  }
  
  # Step 6: Report results
  if (verbose) {
    cat("\n=== Cleaning Complete ===\n")
    cat("SNPs kept:", nLoc(gl_clean), "(", 
        round(nLoc(gl_clean)/original_n_snps * 100, 1), "%)\n")
    cat("SNPs removed:", original_n_snps - nLoc(gl_clean), "\n")
    cat("Individuals kept:", nInd(gl_clean), "\n")
    
    # Check where chromosome data is stored in the cleaned object
    if (!is.null(gl_clean@other$chromosome)) {
      cat("\nChromosome distribution in cleaned object:\n")
      print(table(gl_clean@other$chromosome))
    } else if (!is.null(gl_clean@chromosome)) {
      cat("\nChromosome distribution in cleaned object:\n")
      print(table(gl_clean@chromosome))
    } else if (!is.null(gl_clean@other$loc.metrics$chromosome)) {
      cat("\nChromosome distribution in cleaned object:\n")
      print(table(gl_clean@other$loc.metrics$chromosome))
    }
    
    # Additional: Show that SNPs are ordered
    cat("\nFirst 10 chromosome assignments in cleaned object:\n")
    if (!is.null(gl_clean@other$chromosome)) {
      cat(gl_clean@other$chromosome[1:min(10, nLoc(gl_clean))], "\n")
    } else if (!is.null(gl_clean@chromosome)) {
      cat(gl_clean@chromosome[1:min(10, nLoc(gl_clean))], "\n")
    } else if (!is.null(gl_clean@other$loc.metrics$chromosome)) {
      cat(gl_clean@other$loc.metrics$chromosome[1:min(10, nLoc(gl_clean))], "\n")
    }
  }
  
  return(gl_clean)
}

# NEW: Helper function to get position information (for ordering within chromosomes)
get_position_info <- function(gl) {
  # Check various possible locations for position data
  if (!is.null(gl@position)) {
    return(gl@position)
  } else if (!is.null(gl@other$loc.metrics)) {
    # Check for common position column names
    pos_col_names <- c("position", "Position", "POSITION", "POS", "pos", "Pos",
                       "bp", "Bp", "BP", "loc_bp")
    
    for (col_name in pos_col_names) {
      if (col_name %in% colnames(gl@other$loc.metrics)) {
        return(gl@other$loc.metrics[[col_name]])
      }
    }
    
    # If no standard position column, check for pattern matching
    pos_pattern_names <- names(gl@other$loc.metrics)[grepl("pos|POS|Pos|bp|BP|Bp", 
                                                           names(gl@other$loc.metrics))]
    if (length(pos_pattern_names) > 0) {
      return(gl@other$loc.metrics[[pos_pattern_names[1]]])
    }
  }
  
  # If we get here, no position data found
  return(NULL)
}



# Corrected helper function to find chromosome information
find_chromosome_info <- function(gl) {
  # Check various possible locations for chromosome data
  if (!is.null(gl@other$chromosome)) {
    return(list(
      chrom_vec = as.character(gl@other$chromosome),
      location_name = "other$chromosome",
      column_name = NULL
    ))
  } else if (!is.null(gl@chromosome)) {
    return(list(
      chrom_vec = as.character(gl@chromosome),
      location_name = "chromosome",
      column_name = NULL
    ))
  } else if (!is.null(gl@other$loc.metrics)) {
    # Check if there's a chromosome column in loc.metrics
    chrom_col_names <- c("chromosome", "Chromosome", "CHROMOSOME", 
                         "CHROM", "chrom", "Chrom")
    
    for (col_name in chrom_col_names) {
      if (col_name %in% colnames(gl@other$loc.metrics)) {
        return(list(
          chrom_vec = as.character(gl@other$loc.metrics[[col_name]]),
          location_name = "other$loc.metrics$chromosome",
          column_name = col_name
        ))
      }
    }
    
    # If no standard chromosome column, check for pattern matching
    chrom_pattern_names <- names(gl@other$loc.metrics)[grepl("chrom|CHROM|Chrom", 
                                                             names(gl@other$loc.metrics))]
    if (length(chrom_pattern_names) > 0) {
      return(list(
        chrom_vec = as.character(gl@other$loc.metrics[[chrom_pattern_names[1]]]),
        location_name = "other$loc.metrics$chrom_col",
        column_name = chrom_pattern_names[1]
      ))
    }
  }
  
  # If we get here, no chromosome data found
  return(list(
    chrom_vec = NULL,
    location_name = NULL,
    column_name = NULL
  ))
}

# Extended chromosome cleaning function (unchanged)
clean_chromosome_vector_extended <- function(chrom_vec, verbose = TRUE) {
  # Store original for reference
  original <- as.character(chrom_vec)
  n_original <- length(original)
  
  # Step 1: Identify entries to remove
  # Remove empty strings, NAs, contigs, and Unknown
  is_empty <- is.na(chrom_vec) | chrom_vec == "" | chrom_vec == "NA"
  is_contig <- grepl("contig", chrom_vec, ignore.case = TRUE)
  is_unknown <- grepl("unknown", chrom_vec, ignore.case = TRUE)
  
  keep <- !(is_empty | is_contig | is_unknown)
  
  if (verbose) {
    cat("\nRemoval summary:\n")
    cat("  Empty/NA entries:", sum(is_empty), "\n")
    cat("  Contigs:", sum(is_contig), "\n")
    cat("  Unknown entries:", sum(is_unknown), "\n")
    cat("  Entries to keep:", sum(keep), "\n")
  }
  
  # Step 2: Clean the kept entries
  cleaned_text <- original[keep]
  
  # Remove parentheses and their contents
  cleaned_text <- gsub("\\(.*\\)", "", cleaned_text)
  
  # Remove extra whitespace
  cleaned_text <- trimws(cleaned_text)
  
  # Extract numeric part (remove "Vu" prefix or other non-numeric characters)
  numbers <- gsub("[^0-9]", "", cleaned_text)
  
  # Convert to numeric (empty strings become NA)
  numeric_values <- suppressWarnings(as.numeric(numbers))
  
  # Identify which conversions were successful
  valid_numeric <- !is.na(numeric_values)
  
  # Final keep vector (initial keep AND valid numeric conversion)
  final_keep <- rep(FALSE, n_original)
  final_keep[keep] <- valid_numeric
  
  # Final numeric values for valid entries
  final_numeric <- rep(NA, n_original)
  final_numeric[final_keep] <- numeric_values[valid_numeric]
  
  if (verbose) {
    cat("\nNumeric conversion:\n")
    cat("  Successfully converted:", sum(valid_numeric), "\n")
    cat("  Failed conversion:", sum(keep) - sum(valid_numeric), "\n")
    cat("  Final SNPs to keep:", sum(final_keep), "\n")
    
    if (sum(valid_numeric) > 0) {
      cat("\nChromosome number distribution:\n")
      print(table(numeric_values[valid_numeric]))
    }
  }
  
  return(list(
    keep = final_keep,
    cleaned_numeric = final_numeric[final_keep],
    cleaned_text = cleaned_text[valid_numeric],
    removal_stats = c(
      total_original = n_original,
      removed_empty = sum(is_empty),
      removed_contig = sum(is_contig),
      removed_unknown = sum(is_unknown),
      removed_conversion_failed = sum(keep) - sum(valid_numeric),
      final_kept = sum(final_keep)
    )
  ))
}

# Simplified utility function to check chromosome data location
check_chromosome_location <- function(gl) {
  cat("Checking for chromosome data...\n")
  
  # Check standard locations
  if (!is.null(gl@other$chromosome)) {
    cat("Found in: gl@other$chromosome\n")
    cat("Length:", length(gl@other$chromosome), "\n")
    cat("Sample:", head(gl@other$chromosome), "\n")
    return("other$chromosome")
  }
  
  if (!is.null(gl@chromosome)) {
    cat("Found in: gl@chromosome\n")
    cat("Length:", length(gl@chromosome), "\n")
    cat("Sample:", head(gl@chromosome), "\n")
    return("chromosome")
  }
  
  if (!is.null(gl@other$loc.metrics)) {
    cat("Checking gl@other$loc.metrics...\n")
    cat("Available columns:", colnames(gl@other$loc.metrics), "\n")
    
    # Look for chromosome-related columns
    chrom_cols <- names(gl@other$loc.metrics)[grepl("chrom|CHROM|Chrom", 
                                                    names(gl@other$loc.metrics))]
    if (length(chrom_cols) > 0) {
      cat("Found chromosome columns:", chrom_cols, "\n")
      cat("Sample from", chrom_cols[1], ":", 
          head(gl@other$loc.metrics[[chrom_cols[1]]]), "\n")
      return(paste0("other$loc.metrics$", chrom_cols[1]))
    }
  }
  
  cat("No chromosome data found!\n")
  return(NULL)
}

# First, check where your chromosome data is located
# location <- check_chromosome_location(your_genlight_object)

# Then clean the chromosome data
# gl_clean <- clean_genlight_chromosomes(your_genlight_object, verbose = TRUE)

# Or for less output
#gl_clean <- clean_genlight_chromosomes(your_genlight_object, verbose = FALSE)

#=================================================================================================

#' Load phenotypic data with robust error handling and ordering by genotype
load_pheno_data <- function(pheno_file, header = TRUE, skip = 0, aggregate = TRUE) {
  tryCatch({
    cat("Loading phenotype data from:", pheno_file, "\n")
    
    # Check if file exists
    if (!file.exists(pheno_file)) {
      stop("File not found:", pheno_file)
    }
    
    # Read CSV file
    pheno_data <- read.csv(pheno_file, header = header, skip = skip, 
                           stringsAsFactors = FALSE, check.names = FALSE)
    
    cat("Raw data dimensions:", dim(pheno_data), "\n")
    
    # Check if data was loaded successfully
    if (is.null(pheno_data) || nrow(pheno_data) == 0) {
      stop("The file appears to be empty or could not be read.")
    }
    
    # Print column names for debugging
    cat("Available columns:", paste(colnames(pheno_data), collapse = ", "), "\n")
    
    # Look for genotype column with multiple possible names
    genotype_patterns <- c("^GEN$", "^Genotype$", "^genotype$", "^Name$", "^ID$", "^Sample$")
    genotype_col <- NULL
    
    for (pattern in genotype_patterns) {
      matches <- grep(pattern, colnames(pheno_data), ignore.case = TRUE, value = TRUE)
      if (length(matches) > 0) {
        genotype_col <- matches[1]
        break
      }
    }
    
    # If no genotype column found, use the first column as a fallback
    if (is.null(genotype_col)) {
      warning("No standard genotype column found. Using first column as genotype identifier.")
      genotype_col <- colnames(pheno_data)[1]
    }
    
    cat("Using column '", genotype_col, "' as genotype identifier\n")
    
    # Rename to standard "Genotype" column
    colnames(pheno_data)[colnames(pheno_data) == genotype_col] <- "Genotype"
    
    # Clean genotype column
    pheno_data$Genotype <- as.character(trimws(pheno_data$Genotype))
    
    # Remove rows with missing genotypes
    pheno_data <- pheno_data[!is.na(pheno_data$Genotype) & pheno_data$Genotype != "", ]
    
    if (nrow(pheno_data) == 0) {
      stop("All rows had missing genotype names.")
    }
    
    cat("After cleaning:", nrow(pheno_data), "rows\n")
    
    # Convert obvious numeric columns (common traits)
    numeric_cols <- c("DTF", "DTM", "PL", "NSP", "NPP", "GYP", "GYR", "TW", 
                      "Yield", "Height", "Weight", "Length", "Width", "Count")
    
    for (col in numeric_cols) {
      if (col %in% colnames(pheno_data)) {
        pheno_data[[col]] <- as.numeric(as.character(pheno_data[[col]]))
        cat("Converted", col, "to numeric\n")
      }
    }
    
    # If aggregate is TRUE, do simple aggregation
    if (aggregate) {
      cat("\nAggregating data by genotype...\n")
      
      # Identify numeric columns (excluding Genotype)
      numeric_cols <- names(pheno_data)[sapply(pheno_data, is.numeric)]
      numeric_cols <- setdiff(numeric_cols, "Genotype")
      
      if (length(numeric_cols) > 0) {
        # Simple base R aggregation
        unique_genotypes <- unique(pheno_data$Genotype)
        
        # Sort genotypes (alphabetically or numerically)
        unique_genotypes <- sort(unique_genotypes, na.last = TRUE)
        
        aggregated <- data.frame(Genotype = unique_genotypes)
        
        for (col in numeric_cols) {
          means <- sapply(unique_genotypes, function(geno) {
            mean(pheno_data[pheno_data$Genotype == geno, col], na.rm = TRUE)
          })
          aggregated[[col]] <- means
        }
        
        # Ensure Genotype column is ordered
        aggregated <- aggregated[order(aggregated$Genotype), ]
        
        # Reset row names to maintain clean ordering
        rownames(aggregated) <- NULL
        
        cat("Aggregation complete:", nrow(aggregated), "unique genotypes\n")
        cat("Genotypes ordered alphabetically\n")
        
        return(aggregated)
      } else {
        cat("No numeric columns to aggregate. Returning original data.\n")
        
        # Order original data by Genotype
        pheno_data <- pheno_data[order(pheno_data$Genotype), ]
        rownames(pheno_data) <- NULL
        
        return(pheno_data)
      }
    } else {
      # Not aggregating, but still order by Genotype
      pheno_data <- pheno_data[order(pheno_data$Genotype), ]
      rownames(pheno_data) <- NULL
      
      return(pheno_data)
    }
    
  }, error = function(e) {
    stop(paste("Error loading phenotypic data:", e$message))
  })
}



#' Load metadata
load_meta_data <- function(metadata_file) {
  tryCatch({
    cat("Loading metadata data from:", metadata_file, "\n")
    
    if (!file.exists(metadata_file)) {
      stop("File not found:", metadata_file)
    }
    
    # Read the file
    meta_data <- read.csv(metadata_file, stringsAsFactors = FALSE, 
                          check.names = FALSE, fileEncoding = "UTF-8")
    
    cat("Metadata loaded successfully\n")
    cat("  Number of entries:", nrow(meta_data), "\n")
    cat("  Columns:", paste(colnames(meta_data), collapse = ", "), "\n")
    
    return(meta_data)
  }, error = function(e) {
    stop(paste("Error loading metadata:", e$message))
  })
}

#' Load diallel data
load_diallel_data <- function(diallel_file) {
  tryCatch({
    cat("Loading diallel data from:", diallel_file, "\n")
    
    if (!file.exists(diallel_file)) {
      stop("File not found:", diallel_file)
    }
    
    # Read the file
    diallel_data <- read.csv(diallel_file, stringsAsFactors = FALSE, 
                             check.names = FALSE, fileEncoding = "UTF-8")
    
    cat("Diallel data loaded successfully\n")
    cat("  Number of entries:", nrow(diallel_data), "\n")
    cat("  Columns:", paste(colnames(diallel_data), collapse = ", "), "\n")
    
    return(diallel_data)
  }, error = function(e) {
    stop(paste("Error loading diallel data:", e$message))
  })
}

#' Load METAN data
load_metan_data <- function(metan_file) {
  tryCatch({
    cat("Loading METAN data from:", metan_file, "\n")
    
    if (!file.exists(metan_file)) {
      stop("File not found:", metan_file)
    }
    
    # Read the file
    metan_data <- read.csv(metan_file, stringsAsFactors = FALSE,
                           check.names = FALSE, fileEncoding = "UTF-8")
    
    cat("METAN data loaded successfully\n")
    cat("  Observations:", nrow(metan_data), "\n")
    cat("  Variables:", ncol(metan_data), "\n")
    
    return(metan_data)
  }, error = function(e) {
    stop(paste("Error loading METAN data:", e$message))
  })
}

#' Extract and match genotypes with phenotypes
extract_and_match_genotypes <- function(geno_matrix, pheno_data, start_col = 22) {
  cat("Extracting genotypes from column", start_col, "onward...\n")
  
  tryCatch({
    # Check if geno_matrix has enough columns
    if (ncol(geno_matrix) < start_col) {
      cat("  Warning: Start column is greater than number of columns in genotype matrix\n")
      cat("  Using all columns instead\n")
      geno_extracted <- geno_matrix
    } else {
      geno_extracted <- geno_matrix[, start_col:ncol(geno_matrix), drop = FALSE]
      cat("  Markers after extraction:", ncol(geno_extracted), "\n")
    }
    
    # Clean genotype names in both datasets
    clean_names <- function(x) {
      x <- as.character(x)
      x <- trimws(x)
      x <- gsub("^\\s+|\\s+$", "", x)
      x <- gsub("\"", "", x)  # Remove quotes
      return(x)
    }
    
    # Clean phenotype genotype names
    pheno_data$Genotype_clean <- clean_names(pheno_data$Genotype)
    
    # Clean genotype matrix row names
    if (is.null(rownames(geno_extracted))) {
      stop("Genotype matrix must have row names")
    }
    
    rownames_clean <- clean_names(rownames(geno_extracted))
    rownames(geno_extracted) <- rownames_clean
    
    # Find common genotypes between phenotype and genotype data
    common_genotypes <- intersect(pheno_data$Genotype_clean, rownames_clean)
    
    if (length(common_genotypes) == 0) {
      cat("WARNING: No exact matches found! Trying case-insensitive matching...\n")
      
      # Try case-insensitive matching
      pheno_lower <- tolower(pheno_data$Genotype_clean)
      geno_lower <- tolower(rownames_clean)
      common_lower <- intersect(pheno_lower, geno_lower)
      
      if (length(common_lower) > 0) {
        cat("  Found", length(common_lower), "common genotypes with case-insensitive matching\n")
        
        # Create mapping for case correction
        pheno_map <- setNames(pheno_data$Genotype_clean, pheno_lower)
        geno_map <- setNames(rownames_clean, geno_lower)
        
        common_genotypes <- sapply(common_lower, function(x) pheno_map[x])
        
        # Update genotype matrix row names to lowercase for matching
        rownames(geno_extracted) <- geno_lower
      } else {
        # Show samples for debugging
        cat("  Sample phenotype genotypes:", head(pheno_data$Genotype_clean, 5), "\n")
        cat("  Sample genotype row names:", head(rownames_clean, 5), "\n")
        stop("No common genotypes found even with case-insensitive matching.")
      }
    }
    
    cat("  Common genotypes found:", length(common_genotypes), "\n")
    if (length(common_genotypes) <= 10) {
      cat("  Common genotypes:", paste(common_genotypes, collapse = ", "), "\n")
    } else {
      cat("  First 10 common genotypes:", paste(head(common_genotypes, 10), collapse = ", "), "\n")
    }
    
    # Filter and arrange phenotype data
    pheno_matched <- pheno_data %>% 
      filter(Genotype_clean %in% common_genotypes) %>%
      arrange(match(Genotype_clean, common_genotypes))
    
    # Filter genotype data
    if (all(rownames(geno_extracted) %in% tolower(common_genotypes))) {
      # Case where we used lowercase matching
      geno_matched <- geno_extracted[tolower(common_genotypes), , drop = FALSE]
    } else {
      geno_matched <- geno_extracted[common_genotypes, , drop = FALSE]
    }
    
    # Ensure the order matches
    if (!all(rownames(geno_matched) == pheno_matched$Genotype_clean)) {
      if (all(tolower(rownames(geno_matched)) == tolower(pheno_matched$Genotype_clean))) {
        # Case-insensitive match
        rownames(geno_matched) <- pheno_matched$Genotype_clean
      } else {
        # Reorder genotype matrix
        geno_matched <- geno_matched[pheno_matched$Genotype_clean, , drop = FALSE]
      }
    }
    
    cat("Final matched data:\n")
    cat("  Phenotypes:", nrow(pheno_matched), "genotypes\n")
    cat("  Genotypes:", nrow(geno_matched), "genotypes with", ncol(geno_matched), "markers\n")
    
    # Remove temporary column
    pheno_matched$Genotype_clean <- NULL
    
    return(list(
      pheno = pheno_matched,
      geno = geno_matched
    ))
  }, error = function(e) {
    stop(paste("Error matching genotypes and phenotypes:", e$message))
  })
}

#' Prepare genotype data with chromosome handling
prepare_geno_data <- function(gl) {
  if (is.null(gl)) return(NULL)
  
  tryCatch({
    # Find chromosome column
    Chrom_col <- names(gl@other$loc.metrics)[grepl("^Chrom_|^chrom_|^CHROM_|^Chrom$|^chrom$", 
                                                   names(gl@other$loc.metrics))]
    
    # Find position column
    Chrom_pos_col <- names(gl@other$loc.metrics)[grepl("^SnpPosition|Snp_Position|^Pos|^POS|^position|^Position$", 
                                                       names(gl@other$loc.metrics), ignore.case = TRUE)]
    
    # Use the first matching column if found
    if (length(Chrom_col) > 0) {
      Chrom_data <- gl@other$loc.metrics[[Chrom_col[1]]]
      cat("Found chromosome data in:", Chrom_col[1], "\n")
    } else {
      cat("Warning: Chromosome column not found, using 'Unknown'\n")
      Chrom_data <- rep("Unknown", length(gl@loc.names))
    }
    
    if (length(Chrom_pos_col) > 0) {
      Chrom_pos_data <- gl@other$loc.metrics[[Chrom_pos_col[1]]]
      cat("Found position data in:", Chrom_pos_col[1], "\n")
    } else {
      cat("Warning: Position column not found, using sequential positions\n")
      Chrom_pos_data <- 1:length(gl@loc.names)
    }
    
    # Clean chromosome names
    clean_chromosomes <- function(chrom_vec) {
      sapply(chrom_vec, function(x) {
        if (is.na(x) || x == "" || x == "NA" || x == "na") {
          return("Unknown")
        }
        x_clean <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
        x_clean <- gsub("^chromosome", "", x_clean, ignore.case = TRUE)
        x_clean <- gsub("^chrom", "", x_clean, ignore.case = TRUE)
        x_clean <- gsub("\\(.*\\)", "", x_clean)
        x_clean <- gsub("^0+", "", x_clean)  # Remove leading zeros
        x_clean <- trimws(x_clean)
        
        if (x_clean == "" || is.na(x_clean)) {
          return("Unknown")
        }
        return(x_clean)
      })
    }
    
    Chrom_data_clean <- clean_chromosomes(Chrom_data)
    
    # Extract SNP information
    snp_info <- data.frame(
      SNP_ID = factor(gl@loc.names),
      Chromosome = as.character(Chrom_data_clean),
      Position = as.integer(as.numeric(Chrom_pos_data)),
      stringsAsFactors = FALSE
    )
    
    # Add additional metrics if available
    metric_names <- names(gl@other$loc.metrics)
    
    if("CallRate" %in% metric_names) {
      snp_info$CallRate <- 1 - gl@other$loc.metrics$CallRate
    }
    if("RepAvg" %in% metric_names) {
      snp_info$Reproducibility <- gl@other$loc.metrics$RepAvg
    }
    if("OneRatioSnp" %in% metric_names) {
      snp_info$OneRatio <- gl@other$loc.metrics$OneRatioSnp
    }
    if("PICSnp" %in% metric_names) {
      snp_info$PIC <- gl@other$loc.metrics$PICSnp
    }
    if("AvgPIC" %in% metric_names) {
      snp_info$AvgPIC <- gl@other$loc.metrics$AvgPIC
    }
    if("FreqHomRef" %in% metric_names) {
      snp_info$FreqHomRef <- gl@other$loc.metrics$FreqHomRef
    }
    if("FreqHomSnp" %in% metric_names) {
      snp_info$FreqHomSnp <- gl@other$loc.metrics$FreqHomSnp
    }
    if("FreqHets" %in% metric_names) {
      snp_info$FreqHets <- gl@other$loc.metrics$FreqHets
    }
    
    # Replace any remaining NA chromosomes
    snp_info$Chromosome[is.na(snp_info$Chromosome)] <- "Unknown"
    snp_info$Position[is.na(snp_info$Position)] <- 1:nrow(snp_info)
    
    # Get genotype matrix
    snp_matrix <- as.matrix(gl)
    rownames(snp_matrix) <- indNames(gl)
    
    # Clean sample names
    rownames(snp_matrix) <- gsub("^\\s+|\\s+$", "", rownames(snp_matrix))
    
    cat("Processed genotypic data:\n")
    cat("  Samples:", nrow(snp_matrix), "\n")
    cat("  SNPs:", nrow(snp_info), "\n")
    
    # Count chromosomes
    chr_table <- table(snp_info$Chromosome)
    cat("  Chromosome distribution (top 10):\n")
    if (length(chr_table) > 10) {
      print(head(chr_table, 10))
      cat("  ... and", length(chr_table) - 10, "more\n")
    } else {
      print(chr_table)
    }
    
    return(list(
      snp_matrix = snp_matrix,
      snp_info = snp_info,
      sample_names = rownames(snp_matrix)
    ))
  }, error = function(e) {
    cat("Error preparing genotype data:", e$message, "\n")
    return(NULL)
  })
}

#' Perform PCA analysis
perform_pca_analysis <- function(gl_object, n_pcs = 5, scale = TRUE, center = TRUE) {
  if (is.null(gl_object)) {
    return(NULL)
  }
  
  tryCatch({
    # Convert genlight to matrix
    geno_matrix <- as.matrix(gl_object)
    
    # Remove columns with zero variance
    col_vars <- apply(geno_matrix, 2, var, na.rm = TRUE)
    valid_cols <- which(col_vars > 0 & !is.na(col_vars))
    
    if (length(valid_cols) == 0) {
      stop("No variable columns found for PCA")
    }
    
    geno_matrix <- geno_matrix[, valid_cols]
    
    # Handle missing values (impute with column mean)
    if (any(is.na(geno_matrix))) {
      cat("  Imputing missing values with column means\n")
      geno_matrix <- apply(geno_matrix, 2, function(x) {
        x_mean <- mean(x, na.rm = TRUE)
        x[is.na(x)] <- x_mean
        return(x)
      })
    }
    
    # Perform PCA
    cat("  Performing PCA on", dim(geno_matrix), "matrix\n")
    pca_result <- prcomp(geno_matrix, scale. = scale, center = center)
    
    # Calculate variance explained
    variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
    cumulative_variance <- cumsum(variance_explained)
    
    # Prepare results
    pca_scores <- as.data.frame(pca_result$x[, 1:min(n_pcs, ncol(pca_result$x))])
    pca_scores$Sample <- rownames(geno_matrix)
    pca_scores$SampleID <- rownames(geno_matrix)
    
    pca_variance <- data.frame(
      PC = 1:length(variance_explained),
      Variance = variance_explained,
      Cumulative = cumulative_variance
    )
    
    pca_loadings <- as.data.frame(pca_result$rotation[, 1:min(n_pcs, ncol(pca_result$rotation))])
    
    cat("  PCA completed successfully\n")
    cat("  Variance explained by first 5 PCs:", round(variance_explained[1:5], 2), "%\n")
    
    return(list(
      pca_result = pca_result,
      scores = pca_scores,
      variance = pca_variance,
      loadings = pca_loadings,
      eigenvalues = pca_result$sdev^2
    ))
    
  }, error = function(e) {
    cat("PCA Error:", e$message, "\n")
    return(NULL)
  })
}

#' Perform kinship analysis
perform_kinship_analysis <- function(gl_object, method = "vanRaden", scale = TRUE) {
  if (is.null(gl_object)) {
    return(NULL)
  }
  
  tryCatch({
    # Convert genlight to matrix
    geno_matrix <- as.matrix(gl_object)
    
    # Remove columns with zero variance
    col_vars <- apply(geno_matrix, 2, var, na.rm = TRUE)
    valid_cols <- which(col_vars > 0 & !is.na(col_vars))
    
    if (length(valid_cols) == 0) {
      stop("No variable columns found for kinship matrix")
    }
    
    geno_matrix <- geno_matrix[, valid_cols]
    
    # Handle missing values
    if (any(is.na(geno_matrix))) {
      cat("  Imputing missing values for kinship calculation\n")
      geno_matrix <- apply(geno_matrix, 2, function(x) {
        x_mean <- mean(x, na.rm = TRUE)
        x[is.na(x)] <- x_mean
        return(x)
      })
    }
    
    # Calculate kinship matrix based on method
    if (method == "vanRaden") {
      # VanRaden method (equivalent to G matrix in GBLUP)
      Z <- scale(geno_matrix, center = TRUE, scale = TRUE)
      Z[is.na(Z)] <- 0
      kinship <- tcrossprod(Z) / ncol(Z)
    } else if (method == "loiselle") {
      # Loiselle method
      p <- colMeans(geno_matrix, na.rm = TRUE) / 2
      p[p == 0] <- 0.001  # Avoid division by zero
      p[p == 1] <- 0.999
      
      Z <- geno_matrix - 2 * matrix(p, nrow = nrow(geno_matrix), 
                                    ncol = ncol(geno_matrix), byrow = TRUE)
      kinship <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
    } else if (method == "simple") {
      # Simple covariance matrix
      kinship <- cov(geno_matrix, use = "pairwise.complete.obs")
    } else if (method == "identity") {
      # Identity matrix (no kinship)
      n <- nrow(geno_matrix)
      kinship <- diag(n)
      rownames(kinship) <- colnames(kinship) <- rownames(geno_matrix)
    } else {
      stop("Unknown kinship method:", method)
    }
    
    # Ensure symmetric and positive semi-definite
    kinship <- (kinship + t(kinship)) / 2
    
    if (scale && !all(is.na(kinship))) {
      max_val <- max(abs(kinship), na.rm = TRUE)
      if (max_val > 0) {
        kinship <- kinship / max_val
      }
    }
    
    cat("  Kinship matrix calculated successfully\n")
    cat("  Dimensions:", dim(kinship), "\n")
    cat("  Method:", method, "\n")
    cat("  Range:", round(range(kinship, na.rm = TRUE), 3), "\n")
    
    return(kinship)
    
  }, error = function(e) {
    cat("Kinship Error:", e$message, "\n")
    return(NULL)
  })
}

#' Run single trait GWAS
run_single_trait_gwas <- function(pheno_data, geno_matrix, snp_info, trait) {
  tryCatch({
    cat("Running GWAS for trait:", trait, "\n")
    
    # Get trait values
    trait_values <- pheno_data[[trait]]
    sample_names <- pheno_data$Genotype
    
    # Check if we have enough data
    if (sum(!is.na(trait_values)) < 10) {
      stop("Insufficient non-missing trait values (need at least 10)")
    }
    
    # Ensure geno_matrix has the same samples in the same order
    if (!all(sample_names %in% rownames(geno_matrix))) {
      stop("Not all phenotype samples have genotype data")
    }
    
    # Align genotype matrix with phenotype data
    geno_aligned <- geno_matrix[sample_names, , drop = FALSE]
    
    # Initialize results data frame
    results <- data.frame(
      SNP_ID = character(),
      Chromosome = character(),
      Position = numeric(),
      P_value = numeric(),
      Effect = numeric(),
      SE = numeric(),
      R_squared = numeric(),
      N = numeric(),
      stringsAsFactors = FALSE
    )
    
    total_snps <- ncol(geno_aligned)
    cat("  Testing", total_snps, "SNPs\n")
    
    # Progress indicator
    progress_step <- max(1, floor(total_snps / 10))
    
    # Test each SNP
    for (i in 1:total_snps) {
      snp_id <- colnames(geno_aligned)[i]
      snp_values <- geno_aligned[, i]
      
      # Remove samples with missing values in either trait or SNP
      complete_cases <- complete.cases(trait_values, snp_values)
      y <- trait_values[complete_cases]
      x <- snp_values[complete_cases]
      
      # Check if we have enough data and variation
      if (length(y) >= 10 && length(unique(x)) > 1 && var(x, na.rm = TRUE) > 0) {
        tryCatch({
          # Fit linear model
          model <- lm(y ~ x)
          model_summary <- summary(model)
          
          if (nrow(model_summary$coefficients) >= 2) {
            p_value <- model_summary$coefficients[2, 4]
            effect <- coef(model)[2]
            se <- model_summary$coefficients[2, 2]
            r_squared <- model_summary$r.squared
            
            # Get SNP information
            if (i <= nrow(snp_info)) {
              chromosome <- as.character(snp_info$Chromosome[i])
              position <- as.numeric(snp_info$Position[i])
            } else {
              chromosome <- "Unknown"
              position <- i
            }
            
            results <- rbind(results, data.frame(
              SNP_ID = snp_id,
              Chromosome = chromosome,
              Position = position,
              P_value = p_value,
              Effect = effect,
              SE = se,
              R_squared = r_squared,
              N = length(y),
              stringsAsFactors = FALSE
            ))
          }
        }, error = function(e) {
          # Skip SNPs that cause errors
        })
      }
      
      # Progress reporting
      if (i %% progress_step == 0) {
        cat("    Processed", i, "SNPs (", round(i/total_snps*100, 1), "%)\n")
      }
    }
    
    # Adjust for multiple testing
    if (nrow(results) > 0) {
      results$P_adjusted <- p.adjust(results$P_value, method = "fdr")
      results$neg_log10_p <- -log10(results$P_value)
      results$neg_log10_padj <- -log10(results$P_adjusted)
      results <- results[order(results$P_value), ]
      
      # Calculate Bonferroni threshold
      results$Bonferroni_threshold <- 0.05 / nrow(results)
      results$FDR_threshold <- 0.05
      
      # Flag significant SNPs
      results$Significant_Bonferroni <- results$P_value < results$Bonferroni_threshold
      results$Significant_FDR <- results$P_adjusted < 0.05
    }
    
    cat("  Successful tests:", nrow(results), "/", total_snps, "\n")
    if (nrow(results) > 0) {
      cat("  Minimum P-value:", min(results$P_value, na.rm = TRUE), "\n")
      cat("  Significant SNPs (FDR < 0.05):", sum(results$Significant_FDR, na.rm = TRUE), "\n")
    }
    
    return(results)
    
  }, error = function(e) {
    cat("Error in GWAS for trait", trait, ":", e$message, "\n")
    return(NULL)
  })
}

#=======================================================================================================
# MULTI-TRAIT GWAS HELPERS
#=======================================================================================================

# Function to run multi-trait GWAS
run_multi_trait_gwas <- function(pheno_data, geno_matrix, snp_info, traits, 
                                 method = "Combined P-values", 
                                 include_correlations = TRUE,
                                 calculate_pleiotropy = TRUE) {
  
  tryCatch({
    # Find common samples
    common_samples <- rownames(geno_matrix)
    for (trait in traits) {
      na_mask <- is.na(pheno_data[[trait]])
      common_samples <- intersect(common_samples, pheno_data$Genotype[!na_mask])
    }
    
    cat("Common samples:", length(common_samples), "\n")
    
    if (length(common_samples) < 5) {
      stop("Insufficient common samples (<5)")
    }
    
    # Subset data to common samples
    pheno_subset <- pheno_data[pheno_data$Genotype %in% common_samples, ]
    geno_subset <- geno_matrix[common_samples, ]
    
    # Reorder to match
    pheno_subset <- pheno_subset[match(common_samples, pheno_subset$Genotype), ]
    geno_subset <- geno_subset[common_samples, ]
    
    # Run GWAS for each trait
    trait_results <- list()
    for (trait in traits) {
      cat("Running GWAS for trait:", trait, "\n")
      
      # Run single trait GWAS
      gwas_result <- run_single_trait_gwas(
        pheno_data = pheno_subset,
        geno_matrix = geno_subset,
        snp_info = snp_info,
        trait = trait
      )
      
      trait_results[[trait]] <- gwas_result
    }
    
    # Calculate trait correlations if requested
    trait_correlations <- NULL
    if (include_correlations) {
      trait_correlations <- cor(pheno_subset[, traits, drop = FALSE], 
                                use = "pairwise.complete.obs")
    }
    
    # Combine results based on method
    combined_results <- NULL
    if (method == "Combined P-values") {
      combined_results <- combine_p_values(trait_results)
    } else if (method == "Meta-analysis") {
      combined_results <- meta_analysis_combine(trait_results)
    } else if (method == "Multi-trait BLUP") {
      combined_results <- multi_trait_blup_combine(trait_results, trait_correlations)
    }
    
    # Identify pleiotropic SNPs if requested
    pleiotropic_snps <- NULL
    if (calculate_pleiotropy && !is.null(combined_results)) {
      pleiotropic_snps <- identify_pleiotropic_snps(trait_results)
    }
    
    # Return results
    return(list(
      trait_gwas_results = trait_results,
      trait_correlations = trait_correlations,
      combined_results = combined_results,
      pleiotropic_snps = pleiotropic_snps,
      method = method,
      traits_analyzed = traits,
      samples_used = common_samples,
      n_traits = length(traits),
      timestamp = Sys.time()
    ))
    
  }, error = function(e) {
    cat("Error in multi-trait GWAS:", e$message, "\n")
    return(NULL)
  })
}



# Function to combine p-values
combine_p_values <- function(trait_results) {
  tryCatch({
    # Extract p-values for each SNP across traits
    snp_ids <- unique(unlist(lapply(trait_results, function(x) x$SNP_ID)))
    
    combined_df <- data.frame(
      SNP_ID = snp_ids,
      Min_P_value = NA,
      Mean_P_value = NA,
      Max_P_value = NA,
      N_Traits = 0,
      Traits = "",
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(combined_df$SNP_ID)) {
      snp <- combined_df$SNP_ID[i]
      p_values <- numeric()
      traits <- character()
      
      for (trait in names(trait_results)) {
        trait_df <- trait_results[[trait]]
        idx <- which(trait_df$SNP_ID == snp)
        
        if (length(idx) > 0) {
          p_values <- c(p_values, trait_df$P_value[idx])
          traits <- c(traits, trait)
        }
      }
      
      if (length(p_values) > 0) {
        combined_df$Min_P_value[i] <- min(p_values)
        combined_df$Mean_P_value[i] <- mean(p_values)
        combined_df$Max_P_value[i] <- max(p_values)
        combined_df$N_Traits[i] <- length(p_values)
        combined_df$Traits[i] <- paste(traits, collapse = ";")
      }
    }
    
    # Add chromosome and position info from first trait result
    if (length(trait_results) > 0) {
      first_trait <- trait_results[[1]]
      combined_df <- combined_df %>%
        left_join(first_trait[, c("SNP_ID", "Chromosome", "Position")], 
                  by = "SNP_ID")
    }
    
    return(combined_df)
    
  }, error = function(e) {
    cat("Error combining p-values:", e$message, "\n")
    return(NULL)
  })
}



# Function to identify pleiotropic SNPs
identify_pleiotropic_snps <- function(trait_results, p_threshold = 0.05) {
  tryCatch({
    # Get all SNP IDs
    all_snps <- unique(unlist(lapply(trait_results, function(x) x$SNP_ID)))
    
    pleiotropic <- data.frame(
      SNP_ID = all_snps,
      N_Traits_Significant = 0,
      Traits_Significant = "",
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(pleiotropic$SNP_ID)) {
      snp <- pleiotropic$SNP_ID[i]
      sig_traits <- character()
      
      for (trait in names(trait_results)) {
        trait_df <- trait_results[[trait]]
        idx <- which(trait_df$SNP_ID == snp)
        
        if (length(idx) > 0 && trait_df$P_value[idx] < p_threshold) {
          sig_traits <- c(sig_traits, trait)
        }
      }
      
      pleiotropic$N_Traits_Significant[i] <- length(sig_traits)
      pleiotropic$Traits_Significant[i] <- paste(sig_traits, collapse = ";")
    }
    
    # Filter for SNPs significant in multiple traits
    pleiotropic <- pleiotropic[pleiotropic$N_Traits_Significant > 1, ]
    
    return(pleiotropic)
    
  }, error = function(e) {
    cat("Error identifying pleiotropic SNPs:", e$message, "\n")
    return(NULL)
  })
}



#' Combine results from multiple trait GWAS
combine_multi_trait_results <- function(trait_results, snp_info) {
  tryCatch({
    # Get all unique SNPs across all traits
    all_snps <- unique(unlist(lapply(trait_results, function(x) x$SNP_ID)))
    
    # Create a data frame to store combined results
    combined <- data.frame(
      SNP_ID = all_snps,
      stringsAsFactors = FALSE
    )
    
    # Add chromosome and position information
    combined$Chromosome <- snp_info$Chromosome[match(combined$SNP_ID, snp_info$SNP_ID)]
    combined$Position <- snp_info$Position[match(combined$SNP_ID, snp_info$SNP_ID)]
    
    # Add trait-specific p-values
    trait_names <- names(trait_results)
    for (trait in trait_names) {
      trait_df <- trait_results[[trait]]
      if (!is.null(trait_df)) {
        # Match by SNP_ID
        idx <- match(combined$SNP_ID, trait_df$SNP_ID)
        combined[[paste0("P_", trait)]] <- trait_df$P_value[idx]
        combined[[paste0("Effect_", trait)]] <- trait_df$Effect[idx]
      }
    }
    
    # Calculate combined statistics
    p_value_cols <- grep("^P_", colnames(combined), value = TRUE)
    
    if (length(p_value_cols) > 1) {
      # Minimum p-value across traits for each SNP
      p_matrix <- as.matrix(combined[, p_value_cols, drop = FALSE])
      combined$Min_P_value <- apply(p_matrix, 1, min, na.rm = TRUE)
      combined$N_Traits <- rowSums(!is.na(p_matrix))
      
      # Calculate Fisher's combined probability test
      combined$Fisher_combined <- apply(p_matrix, 1, function(p_vals) {
        valid_vals <- p_vals[!is.na(p_vals)]
        if (length(valid_vals) > 0) {
          chi_sq <- -2 * sum(log(valid_vals))
          pchisq(chi_sq, df = 2 * length(valid_vals), lower.tail = FALSE)
        } else {
          NA
        }
      })
      
      # Calculate FDR-adjusted p-values
      combined$P_adjusted <- p.adjust(combined$Min_P_value, method = "fdr")
      combined$Fisher_P_adjusted <- p.adjust(combined$Fisher_combined, method = "fdr")
      
      # Identify pleiotropic SNPs (significant in multiple traits)
      combined$Significant_FDR <- combined$P_adjusted < 0.05
      
      # Flag pleiotropic SNPs (significant in >1 trait with p < 0.05)
      sig_counts <- rowSums(p_matrix < 0.05, na.rm = TRUE)
      combined$Pleiotropic <- sig_counts > 1
    }
    
    # Order by minimum p-value
    combined <- combined[order(combined$Min_P_value, na.last = TRUE), ]
    
    return(combined)
    
  }, error = function(e) {
    cat("Error combining multi-trait results:", e$message, "\n")
    return(NULL)
  })
}

#' Calculate genetic correlations
calculate_genetic_correlations <- function(pheno_data, geno_matrix, traits) {
  tryCatch({
    # Prepare data
    samples <- rownames(geno_matrix)
    pheno_subset <- pheno_data[pheno_data$Genotype %in% samples, ]
    
    # Ensure same order
    pheno_subset <- pheno_subset[match(samples, pheno_subset$Genotype), ]
    
    # Extract trait values
    trait_matrix <- as.matrix(pheno_subset[, traits, drop = FALSE])
    
    # Calculate genetic correlations using variance components
    n_traits <- length(traits)
    genetic_cor <- matrix(NA, n_traits, n_traits)
    colnames(genetic_cor) <- rownames(genetic_cor) <- traits
    
    for (i in 1:n_traits) {
      for (j in 1:n_traits) {
        if (i <= j) {
          # Simplified genetic correlation calculation
          # Using marker-based approach
          tryCatch({
            # Extract trait values
            trait_i <- scale(trait_matrix[, i])
            trait_j <- scale(trait_matrix[, j])
            
            # Remove NAs
            complete <- complete.cases(trait_i, trait_j)
            if (sum(complete) > 10) {
              trait_i <- trait_i[complete]
              trait_j <- trait_j[complete]
              geno_sub <- geno_matrix[complete, ]
              
              # Calculate marker effects (simplified)
              # In practice, you might want to use a proper mixed model
              effects_i <- colMeans(geno_sub * trait_i, na.rm = TRUE)
              effects_j <- colMeans(geno_sub * trait_j, na.rm = TRUE)
              
              # Calculate correlation of marker effects
              genetic_cor[i, j] <- genetic_cor[j, i] <- cor(effects_i, effects_j, 
                                                            use = "pairwise.complete.obs")
            }
          }, error = function(e) {
            cat("  Error calculating genetic correlation between", 
                traits[i], "and", traits[j], ":", e$message, "\n")
          })
        }
      }
    }
    
    return(genetic_cor)
    
  }, error = function(e) {
    cat("Error in genetic correlation calculation:", e$message, "\n")
    return(NULL)
  })
}

#' Identify pleiotropic hotspots
identify_pleiotropic_hotspots <- function(gwas_results, traits, window_size = 1000000) {
  tryCatch({
    # Combine all GWAS results
    all_results <- data.frame()
    
    for (trait in traits) {
      if (!is.null(gwas_results[[trait]])) {
        trait_df <- gwas_results[[trait]]
        trait_df$Trait <- trait
        all_results <- rbind(all_results, trait_df)
      }
    }
    
    if (nrow(all_results) == 0) {
      return(NULL)
    }
    
    # Group by chromosome and position windows
    hotspots <- all_results %>%
      group_by(Chromosome) %>%
      arrange(Position) %>%
      mutate(
        Window = floor(Position / window_size) * window_size
      ) %>%
      group_by(Chromosome, Window) %>%
      summarise(
        N_SNPs = n(),
        N_Traits = length(unique(Trait)),
        Min_P_value = min(P_value, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      filter(N_Traits > 1) %>%
      arrange(desc(N_Traits), Min_P_value)
    
    return(hotspots)
    
  }, error = function(e) {
    cat("Error identifying pleiotropic hotspots:", e$message, "\n")
    return(NULL)
  })
}




#=======================================================================================================
#' Create QC plots
create_qc_plots <- function(gl_object) {
  if (is.null(gl_object)) return(NULL)
  
  tryCatch({
    # Calculate basic statistics
    geno_matrix <- as.matrix(gl_object)
    
    # Calculate missing data per sample
    sample_missing <- rowMeans(is.na(geno_matrix))
    
    # Calculate missing data per SNP
    snp_missing <- colMeans(is.na(geno_matrix))
    
    # Calculate MAF
    maf <- colMeans(geno_matrix, na.rm = TRUE) / 2
    maf <- pmin(maf, 1 - maf)
    
    # Calculate heterozygosity
    het_rate <- rowMeans(geno_matrix == 1, na.rm = TRUE)
    
    plots <- list()
    
    # Plot 1: Sample missing data
    plots$sample_missing <- ggplot(data.frame(Sample = 1:length(sample_missing), 
                                              Missing = sample_missing),
                                   aes(x = Sample, y = Missing)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = mean(sample_missing), color = "red", linetype = "dashed") +
      labs(title = "Missing Data per Sample", 
           x = "Sample Index", 
           y = "Missing Rate") +
      theme_minimal()
    
    # Plot 2: SNP missing data
    plots$snp_missing <- ggplot(data.frame(SNP = 1:length(snp_missing), 
                                           Missing = snp_missing),
                                aes(x = SNP, y = Missing)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = mean(snp_missing), color = "red", linetype = "dashed") +
      labs(title = "Missing Data per SNP", 
           x = "SNP Index", 
           y = "Missing Rate") +
      theme_minimal()
    
    # Plot 3: MAF distribution
    plots$maf_dist <- ggplot(data.frame(MAF = maf), aes(x = MAF)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      labs(title = "Minor Allele Frequency Distribution", 
           x = "MAF", 
           y = "Count") +
      theme_minimal()
    
    # Plot 4: Heterozygosity
    plots$het_dist <- ggplot(data.frame(Het = het_rate), aes(x = Het)) +
      geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
      labs(title = "Heterozygosity Distribution", 
           x = "Heterozygosity Rate", 
           y = "Count") +
      theme_minimal()
    
    return(plots)
    
  }, error = function(e) {
    cat("Error creating QC plots:", e$message, "\n")
    return(NULL)
  })
}

#==========================================================================================================
#METAN FUNCTIONS
#==========================================================================================================

#' Fixed stability analysis function
run_stability_analysis <- function(metan_data, traits) {
  tryCatch({
    cat("Running stability analysis for", length(traits), "traits\n")
    
    # Ensure traits are numeric
    for (trait in traits) {
      if (!is.numeric(metan_data[[trait]])) {
        metan_data[[trait]] <- as.numeric(as.character(metan_data[[trait]]))
      }
    }
    
    # Run stability analysis
    stability_result <- tryCatch({
      ecovalence(metan_data, ENV, GEN, REP, 
                 resp = all_of(traits),
                 verbose = FALSE)
    }, error = function(e) {
      cat("  Stability analysis error:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(stability_result)) {
      return(list(error = "Stability analysis failed"))
    }
    
    # Convert to data frame safely
    stability_df <- tryCatch({
      if (inherits(stability_result, "ecovalence")) {
        # Extract the data frame from ecovalence object
        as.data.frame(stability_result)
      } else {
        as.data.frame(stability_result)
      }
    }, error = function(e) {
      cat("  Error converting stability result to data frame:", e$message, "\n")
      return(NULL)
    })
    
    return(list(
      stability_result = stability_result,
      stability_df = stability_df
    ))
    
  }, error = function(e) {
    return(list(error = paste("Stability analysis failed:", e$message)))
  })
}



#' Prepare data for METAN analysis with automatic ENV column creation
prepare_metan_data <- function(raw_data, env_var = NULL, gen_var, rep_var, resp_vars) {
  tryCatch({
    cat("Preparing METAN data...\n")
    
    # Create a copy of the data
    metan_data <- as.data.frame(raw_data)
    
    # Check if YEAR and SEASON columns exist to create ENV column
    year_col <- NULL
    season_col <- NULL
    
    # Find YEAR column (case-insensitive)
    year_patterns <- c("^YEAR$", "^Year$", "^year$")
    for (pattern in year_patterns) {
      matches <- grep(pattern, colnames(metan_data), ignore.case = FALSE, value = TRUE)
      if (length(matches) > 0) {
        year_col <- matches[1]
        break
      }
    }
    
    # Find SEASON column (case-insensitive)
    season_patterns <- c("^SEASON$", "^Season$", "^season$")
    for (pattern in season_patterns) {
      matches <- grep(pattern, colnames(metan_data), ignore.case = FALSE, value = TRUE)
      if (length(matches) > 0) {
        season_col <- matches[1]
        break
      }
    }
    
    # Track if we created ENV from YEAR and SEASON
    env_created_from_year_season <- FALSE
    
    # If both YEAR and SEASON exist, create ENV column
    if (!is.null(year_col) && !is.null(season_col)) {
      cat("Found YEAR (", year_col, ") and SEASON (", season_col, ") columns.\n")
      cat("Creating ENV column by combining YEAR and SEASON...\n")
      
      # Clean and combine
      year_vals <- as.character(trimws(metan_data[[year_col]]))
      season_vals <- as.character(trimws(metan_data[[season_col]]))
      
      # Remove NA or empty values
      year_vals[is.na(year_vals)] <- "Unknown"
      season_vals[is.na(season_vals)] <- "Unknown"
      year_vals[year_vals == ""] <- "Unknown"
      season_vals[season_vals == ""] <- "Unknown"
      
      # Create ENV column
      metan_data$ENV <- paste(year_vals, season_vals, sep = "_")
      
      env_created_from_year_season <- TRUE
      cat("Created ENV column with", length(unique(metan_data$ENV)), "unique environments\n")
      
    } else if (!is.null(env_var) && env_var %in% colnames(metan_data)) {
      # Use the provided environment variable
      cat("Using provided environment variable:", env_var, "\n")
      metan_data <- metan_data %>%
        rename(ENV = !!sym(env_var))
    } else {
      # Try to find a suitable environment column
      env_candidates <- c("ENV", "Environment", "environment", "LOCATION", 
                          "Location", "Site", "TRIAL", "Trial")
      env_found <- FALSE
      
      for (candidate in env_candidates) {
        if (candidate %in% colnames(metan_data)) {
          metan_data <- metan_data %>%
            rename(ENV = !!sym(candidate))
          cat("Found and using environment column:", candidate, "\n")
          env_found <- TRUE
          break
        }
      }
      
      if (!env_found) {
        stop("No environment column found. Please ensure data has either:\n",
             "1. YEAR and SEASON columns (will be combined to create ENV)\n",
             "2. An environment column (ENV, Environment, etc.)\n",
             "3. Or select an environment variable in the UI")
      }
    }
    
    # Convert genotype variable to GEN
    if (!gen_var %in% colnames(metan_data)) {
      stop("Genotype variable '", gen_var, "' not found in data")
    }
    metan_data <- metan_data %>%
      rename(GEN = !!sym(gen_var))
    
    # Convert replicate variable to REP
    if (!rep_var %in% colnames(metan_data)) {
      stop("Replicate variable '", rep_var, "' not found in data")
    }
    metan_data <- metan_data %>%
      rename(REP = !!sym(rep_var))
    
    # Select only necessary columns
    required_cols <- c("ENV", "GEN", "REP")
    all_cols <- c(required_cols, resp_vars)
    
    # Check if all response variables exist
    missing_resp <- setdiff(resp_vars, colnames(metan_data))
    if (length(missing_resp) > 0) {
      stop("Response variable(s) not found: ", paste(missing_resp, collapse = ", "))
    }
    
    # If we created ENV from YEAR and SEASON, keep those columns for now to create summary
    if (env_created_from_year_season) {
      # Store YEAR and SEASON values before removing
      year_season_data <- metan_data[, c("ENV", year_col, season_col)]
      colnames(year_season_data) <- c("ENV", "YEAR", "SEASON")
    }
    
    metan_data <- metan_data %>%
      select(all_of(all_cols))
    
    # Clean factor columns
    metan_data$ENV <- as.factor(trimws(as.character(metan_data$ENV)))
    metan_data$GEN <- as.factor(trimws(as.character(metan_data$GEN)))
    metan_data$REP <- as.factor(trimws(as.character(metan_data$REP)))
    
    # Ensure response variables are numeric
    for (resp_var in resp_vars) {
      if (!is.numeric(metan_data[[resp_var]])) {
        metan_data[[resp_var]] <- as.numeric(as.character(metan_data[[resp_var]]))
        cat("Converted", resp_var, "to numeric\n")
      }
    }
    
    # Remove rows with NA in required columns
    initial_rows <- nrow(metan_data)
    metan_data <- metan_data %>%
      filter(!is.na(ENV), !is.na(GEN), !is.na(REP))
    
    # Remove rows where all response variables are NA
    row_has_data <- apply(metan_data[resp_vars], 1, function(x) any(!is.na(x)))
    metan_data <- metan_data[row_has_data, ]
    
    removed_rows <- initial_rows - nrow(metan_data)
    if (removed_rows > 0) {
      cat("Removed", removed_rows, "rows with missing data\n")
    }
    
    cat("METAN data prepared successfully\n")
    cat("  Environments:", length(unique(metan_data$ENV)), "\n")
    cat("  Genotypes:", length(unique(metan_data$GEN)), "\n")
    cat("  Replicates:", length(unique(metan_data$REP)), "\n")
    cat("  Traits:", length(resp_vars), "\n")
    cat("  Observations:", nrow(metan_data), "\n")
    
    # Create environment summary if ENV was created from YEAR and SEASON
    if (env_created_from_year_season && exists("year_season_data")) {
      # Merge with the prepared data to get counts
      env_counts <- metan_data %>%
        group_by(ENV) %>%
        summarise(n_obs = n(), .groups = 'drop')
      
      # Get unique YEAR and SEASON for each ENV
      env_summary <- year_season_data %>%
        distinct(ENV, YEAR, SEASON) %>%
        left_join(env_counts, by = "ENV") %>%
        arrange(YEAR, SEASON)
      
      cat("\nEnvironment details (created from YEAR and SEASON):\n")
      if (nrow(env_summary) <= 10) {
        print(env_summary)
      } else {
        print(head(env_summary, 10))
        cat("... and", nrow(env_summary) - 10, "more environments\n")
      }
    }
    
    return(metan_data)
  }, error = function(e) {
    stop(paste("Error preparing METAN data:", e$message))
  })
}


#' Run METAN analysis
#' Run METAN analysis with enhanced error handling
run_metan_analysis <- function(metan_data, resp_vars) {
  results <- list()
  
  tryCatch({
    # 1. Descriptive Statistics
    results$desc_stats <- tryCatch({
      metan_data %>%
        group_by(GEN) %>%
        desc_stat(all_of(resp_vars), verbose = FALSE)
    }, error = function(e) {
      return(paste("Descriptive stats error:", e$message))
    })
    
    # 2. Individual ANOVA for each trait
    results$anova <- list()
    for (trait in resp_vars) {
      results$anova[[trait]] <- tryCatch({
        anova_ind(metan_data, ENV, GEN, REP, 
                  resp = !!sym(trait), 
                  verbose = FALSE)
      }, error = function(e) {
        return(paste("ANOVA error for", trait, ":", e$message))
      })
    }
    
    # 3. AMMI Analysis for each trait (using enhanced function)
    results$ammi <- list()
    for (trait in resp_vars) {
      results$ammi[[trait]] <- run_ammi_analysis(metan_data, trait)
    }
    
    # 4. WAASB Analysis
    results$waasb <- tryCatch({
      waasb(metan_data, ENV, GEN, REP, 
            resp = all_of(resp_vars))
    }, error = function(e) {
      return(paste("WAASB error:", e$message))
    })
    
    # 5. GGE Biplot for each trait
    results$gge <- list()
    for (trait in resp_vars) {
      results$gge[[trait]] <- tryCatch({
        gge(metan_data, ENV, GEN, 
            resp = !!sym(trait), 
            svp = "symmetrical")
      }, error = function(e) {
        return(paste("GGE error for", trait, ":", e$message))
      })
    }
    
    # 6. GE Plot for each trait
    results$ge_plot <- list()
    for (trait in resp_vars) {
      results$ge_plot[[trait]] <- tryCatch({
        ge_plot(metan_data, ENV, GEN, 
                resp = !!sym(trait), 
                type = 2)
      }, error = function(e) {
        return(paste("GE Plot error for", trait, ":", e$message))
      })
    }
    
    # 7. Correlation Matrix
    if (length(resp_vars) > 1) {
      results$corr <- tryCatch({
        corr_coef(metan_data, all_of(resp_vars))
      }, error = function(e) {
        return(paste("Correlation error:", e$message))
      })
    }
    
    # 8. Stability Analysis (using enhanced function)
    results$stability <- run_stability_analysis(metan_data, resp_vars)
    
    # 9. Genotype Performance Ranking
    results$ranking <- tryCatch({
      mgidi(metan_data, ENV, GEN, REP, 
            resp = all_of(resp_vars))
    }, error = function(e) {
      return(paste("Ranking error:", e$message))
    })
    
    return(results)
    
  }, error = function(e) {
    cat("Error in METAN analysis:", e$message, "\n")
    return(list(error = e$message))
  })
}

# ============================================================================
# METAN SPECIFIC FUNCTIONS
# ============================================================================

# Function 1: Data Quality Control (from Rmd Example 1.3)
check_met_data_completeness <- function(data) {
  # Check for required columns
  required_cols <- c("ENV", "GEN", "REP")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    return(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check for completeness
  completeness_stats <- list(
    env_complete = length(unique(data$ENV)),
    gen_complete = length(unique(data$GEN)),
    rep_complete = length(unique(data$REP)),
    total_obs = nrow(data),
    expected_obs = length(unique(data$ENV)) * 
      length(unique(data$GEN)) * 
      length(unique(data$REP))
  )
  
  return(completeness_stats)
}

# Function 2: Missing Data Summary (from Rmd Example 1.3)
summarize_missing_data <- function(data) {
  missing_summary <- data.frame(
    Variable = character(),
    Missing_Count = numeric(),
    Missing_Percent = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col in names(data)) {
    missing_count <- sum(is.na(data[[col]]))
    missing_percent <- round(missing_count / nrow(data) * 100, 2)
    
    missing_summary <- rbind(missing_summary, data.frame(
      Variable = col,
      Missing_Count = missing_count,
      Missing_Percent = missing_percent
    ))
  }
  
  return(missing_summary)
}

# Function 3: Environmental Characterization (from Rmd Example 2.3)
characterize_environments <- function(data, trait = "YIELD") {
  env_summary <- data %>%
    group_by(ENV) %>%
    summarise(
      Mean = mean(!!sym(trait), na.rm = TRUE),
      SD = sd(!!sym(trait), na.rm = TRUE),
      CV = sd(!!sym(trait), na.rm = TRUE) / mean(!!sym(trait), na.rm = TRUE) * 100,
      Min = min(!!sym(trait), na.rm = TRUE),
      Max = max(!!sym(trait), na.rm = TRUE),
      n_obs = n(),
      .groups = 'drop'
    ) %>%
    arrange(desc(Mean))
  
  return(env_summary)
}

# Function 4: Combined ANOVA (from Rmd Example 3.1)
perform_combined_anova <- function(data, trait = "YIELD") {
  # Check if trait exists
  if (!trait %in% names(data)) {
    return(list(error = paste("Trait", trait, "not found in data")))
  }
  
  # Perform ANOVA
  formula <- as.formula(paste(trait, "~ ENV + GEN + ENV:GEN + REP%in%ENV"))
  anova_result <- tryCatch({
    aov(formula, data = data)
  }, error = function(e) {
    return(list(error = e$message))
  })
  
  if ("error" %in% names(anova_result)) {
    return(anova_result)
  }
  
  # Extract ANOVA table
  anova_summary <- summary(anova_result)
  
  # Calculate variance components
  variance_components <- tryCatch({
    library(lme4)
    mixed_model <- lmer(as.formula(paste(trait, "~ (1|ENV) + (1|GEN) + (1|ENV:GEN)")), 
                        data = data)
    VarCorr(mixed_model)
  }, error = function(e) {
    return(NULL)
  })
  
  return(list(
    anova_table = anova_summary,
    variance_components = variance_components
  ))
}

# Function 5: Finlay-Wilkinson Regression (from Rmd Example 4.1)
perform_finlay_wilkinson <- function(data, trait = "YIELD") {
  # Calculate environmental index
  env_index <- data %>%
    group_by(ENV) %>%
    summarise(
      env_mean = mean(!!sym(trait), na.rm = TRUE),
      .groups = 'drop'
    )
  
  overall_mean <- mean(data[[trait]], na.rm = TRUE)
  env_index$env_index <- env_index$env_mean - overall_mean
  
  # Merge with data
  data_with_index <- data %>%
    left_join(env_index, by = "ENV")
  
  # Perform regression for each genotype
  genotypes <- unique(data$GEN)
  fw_results <- data.frame(
    GEN = character(),
    Intercept = numeric(),
    Slope = numeric(),
    R_squared = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (gen in genotypes) {
    gen_data <- data_with_index[data_with_index$GEN == gen, ]
    
    if (nrow(gen_data) >= 3) {  # Need at least 3 points for regression
      model <- lm(as.formula(paste(trait, "~ env_index")), data = gen_data)
      summary_model <- summary(model)
      
      fw_results <- rbind(fw_results, data.frame(
        GEN = as.character(gen),
        Intercept = coef(model)[1],
        Slope = coef(model)[2],
        R_squared = summary_model$r.squared,
        P_value = summary_model$coefficients[2, 4]
      ))
    }
  }
  
  return(list(
    env_index = env_index,
    genotype_regressions = fw_results,
    overall_mean = overall_mean
  ))
}

# Function 6: AMMI Analysis (from Rmd Example 5.1)
perform_ammi_analysis <- function(data, trait = "YIELD", naxis = 5) {
  # Check if metan package is available
  if (!require(metan)) {
    return(list(error = "METAN package required for AMMI analysis"))
  }
  
  tryCatch({
    # Perform AMMI analysis
    ammi_result <- performs_ammi(data, 
                                 ENV, GEN, REP, 
                                 resp = !!sym(trait),
                                 naxis = naxis)
    
    return(ammi_result)
    
  }, error = function(e) {
    return(list(error = paste("AMMI analysis failed:", e$message)))
  })
}

# Function 7: GGE Biplot (from Rmd Example 5.2)
perform_gge_biplot <- function(data, trait = "YIELD") {
  # Check if metan package is available
  if (!require(metan)) {
    return(list(error = "METAN package required for GGE biplot"))
  }
  
  tryCatch({
    # Perform GGE biplot analysis
    gge_result <- gge(data,
                      ENV, GEN,
                      resp = !!sym(trait),
                      svp = "symmetrical",
                      scaling = "sd")
    
    return(gge_result)
    
  }, error = function(e) {
    return(list(error = paste("GGE biplot failed:", e$message)))
  })
}

# ============================================================================
# DIALLEL ANALYSIS FUNCTIONS
# ============================================================================
# ============================================================================
# ENHANCED DIALLEL ANALYSIS FUNCTIONS
# ============================================================================

#' Convert SCA matrix to long format data frame
sca_matrix_to_df <- function(sca_matrix) {
  if (is.null(sca_matrix) || !is.matrix(sca_matrix)) {
    return(NULL)
  }
  
  sca_df <- data.frame(
    Parent1 = character(),
    Parent2 = character(),
    SCA = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(sca_matrix)) {
    for (j in 1:ncol(sca_matrix)) {
      if (i != j && !is.na(sca_matrix[i, j])) {
        sca_df <- rbind(sca_df, data.frame(
          Parent1 = rownames(sca_matrix)[i],
          Parent2 = colnames(sca_matrix)[j],
          SCA = sca_matrix[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Remove duplicates (since matrix is symmetric)
  if (nrow(sca_df) > 0) {
    sca_df <- sca_df %>%
      mutate(key = paste(pmin(Parent1, Parent2), pmax(Parent1, Parent2))) %>%
      distinct(key, .keep_all = TRUE) %>%
      select(-key) %>%
      arrange(desc(abs(SCA)))
  }
  
  return(sca_df)
}


#' Run ANOVA on all traits in diallel crosses with proper SCA handling
run_diallel_anova_all_traits <- function(diallel_data, trait_cols = NULL) {
  tryCatch({
    if (is.null(trait_cols)) {
      # Find numeric trait columns
      trait_cols <- colnames(diallel_data)[sapply(diallel_data, is.numeric)]
    }
    
    cat("Running diallel ANOVA for", length(trait_cols), "traits...\n")
    
    results <- list()
    
    for (trait in trait_cols) {
      cat("  Analyzing trait:", trait, "\n")
      
      # Run diallel analysis
      diallel_result <- perform_diallel_analysis(diallel_data, trait)
      
      # Check if analysis was successful
      if (!is.null(diallel_result$error)) {
        cat("    Error for trait", trait, ":", diallel_result$error, "\n")
        next
      }
      
      if (!is.null(diallel_result$anova) && 
          !is.null(diallel_result$gca) &&
          nrow(diallel_result$gca) > 0) {
        
        # Extract GCA effects
        gca_df <- diallel_result$gca %>%
          mutate(Trait = trait) %>%
          select(Parent, GCA, Mean_Performance, Trait)
        
        # Extract SCA matrix with proper handling
        sca_matrix <- diallel_result$sca
        
        # Initialize sca_df
        sca_df <- NULL
        
        if (!is.null(sca_matrix) && is.matrix(sca_matrix) && 
            nrow(sca_matrix) > 0 && ncol(sca_matrix) > 0) {
          
          # Convert matrix to long format safely
          sca_long <- data.frame(
            Parent1 = character(),
            Parent2 = character(),
            SCA = numeric(),
            stringsAsFactors = FALSE
          )
          
          # Extract all non-NA, non-diagonal values
          for (i in 1:nrow(sca_matrix)) {
            for (j in 1:ncol(sca_matrix)) {
              if (i != j && !is.na(sca_matrix[i, j])) {
                sca_long <- rbind(sca_long, data.frame(
                  Parent1 = rownames(sca_matrix)[i],
                  Parent2 = colnames(sca_matrix)[j],
                  SCA = sca_matrix[i, j],
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
          
          if (nrow(sca_long) > 0) {
            # Remove duplicates (symmetric matrix)
            sca_long <- sca_long %>%
              mutate(key = paste(pmin(Parent1, Parent2), pmax(Parent1, Parent2))) %>%
              distinct(key, .keep_all = TRUE) %>%
              select(-key)
            
            sca_df <- sca_long %>%
              mutate(Trait = trait) %>%
              select(Parent1, Parent2, SCA, Trait)
          }
        }
        
        # Calculate heterosis
        heterosis_data <- tryCatch({
          calculate_heterosis(
            diallel_data = diallel_data,
            cross_data = diallel_data[diallel_data$CrossType == "Cross", ],
            parent_data = diallel_data[diallel_data$CrossType == "Parent", ],
            trait = trait
          )
        }, error = function(e) {
          cat("    Heterosis calculation error for trait", trait, ":", e$message, "\n")
          NULL
        })
        
        results[[trait]] <- list(
          anova = diallel_result$anova,
          gca = gca_df,
          sca = sca_df,
          heterosis = heterosis_data,
          overall_mean = diallel_result$overall_mean,
          n_parents = diallel_result$n_parents,
          n_crosses = diallel_result$n_crosses
        )
      }
    }
    
    if (length(results) == 0) {
      cat("No valid results obtained for any trait\n")
      return(NULL)
    }
    
    # Combine results into summary tables
    summary_tables <- list()
    
    # GCA summary
    gca_summary <- do.call(rbind, lapply(results, function(x) x$gca))
    if (!is.null(gca_summary) && nrow(gca_summary) > 0) {
      summary_tables$gca_summary <- gca_summary
    }
    
    # SCA summary
    sca_summary <- do.call(rbind, lapply(results, function(x) x$sca))
    if (!is.null(sca_summary) && nrow(sca_summary) > 0) {
      summary_tables$sca_summary <- sca_summary
    }
    
    # Heterosis summary
    heterosis_summary <- do.call(rbind, lapply(results, function(x) x$heterosis))
    if (!is.null(heterosis_summary) && nrow(heterosis_summary) > 0) {
      summary_tables$heterosis_summary <- heterosis_summary
    }
    
    # ANOVA summary
    anova_summary <- data.frame()
    for (trait in names(results)) {
      if (!is.null(results[[trait]]$anova) && 
          !is.null(results[[trait]]$anova$table)) {
        anova_table <- results[[trait]]$anova$table
        
        # Convert anova table to data frame
        if (is.data.frame(anova_table)) {
          anova_summary <- rbind(anova_summary, data.frame(
            Trait = trait,
            Source = rownames(anova_table),
            Df = anova_table$Df,
            Sum_Sq = anova_table$"Sum Sq",
            Mean_Sq = anova_table$"Mean Sq",
            F_value = anova_table$"F value",
            P_value = anova_table$"Pr(>F)",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    if (nrow(anova_summary) > 0) {
      summary_tables$anova_summary <- anova_summary
    }
    
    return(list(
      detailed_results = results,
      summary_tables = summary_tables
    ))
    
  }, error = function(e) {
    cat("Error in diallel ANOVA for all traits:", e$message, "\n")
    return(NULL)
  })
}


# ============================================================================
# ENHANCED METAN ANALYSIS FUNCTIONS
# ============================================================================

#' Fixed AMMI analysis function
run_ammi_analysis <- function(metan_data, trait, naxis = 5) {
  tryCatch({
    cat("Running AMMI analysis for trait:", trait, "\n")
    
    # Check if trait exists
    if (!trait %in% colnames(metan_data)) {
      return(list(error = paste("Trait", trait, "not found in data")))
    }
    
    # Check for missing values
    if (any(is.na(metan_data[[trait]]))) {
      cat("  Warning: Missing values found in trait", trait, "\n")
      cat("  Removing rows with missing values...\n")
      metan_data <- metan_data[!is.na(metan_data[[trait]]), ]
    }
    
    # Check for infinite values
    if (any(is.infinite(metan_data[[trait]]))) {
      cat("  Warning: Infinite values found in trait", trait, "\n")
      cat("  Removing rows with infinite values...\n")
      metan_data <- metan_data[!is.infinite(metan_data[[trait]]), ]
    }
    
    # Ensure we have enough data
    if (nrow(metan_data) < 10) {
      return(list(error = "Insufficient data after removing missing values"))
    }
    
    # Run AMMI analysis with error handling
    ammi_result <- tryCatch({
      performs_ammi(metan_data, 
                    ENV, GEN, REP, 
                    resp = !!sym(trait),
                    naxis = naxis,
                    verbose = FALSE)
    }, error = function(e) {
      cat("  AMMI analysis error:", e$message, "\n")
      
      # Try alternative approach
      tryCatch({
        # Create a simplified model
        formula <- as.formula(paste(trait, "~ ENV + GEN + ENV:GEN"))
        model <- aov(formula, data = metan_data)
        
        # Extract interaction effects
        interaction_effects <- model.tables(model, type = "effects")$"ENV:GEN"
        
        # Perform PCA on interaction matrix
        interaction_matrix <- matrix(interaction_effects, 
                                     nrow = length(unique(metan_data$ENV)),
                                     ncol = length(unique(metan_data$GEN)))
        
        pca_result <- prcomp(interaction_matrix, scale = TRUE)
        
        list(
          model = model,
          pca = pca_result,
          interaction_matrix = interaction_matrix,
          simplified = TRUE
        )
      }, error = function(e2) {
        return(list(error = paste("Both AMMI methods failed:", e2$message)))
      })
    })
    
    return(ammi_result)
    
  }, error = function(e) {
    return(list(error = paste("AMMI analysis failed:", e$message)))
  })
}

#=========================================================================================================

#' Enhanced diallel analysis function with proper SCA matrix handling
perform_diallel_analysis <- function(diallel_data, trait) {
  df <- diallel_data
  crosses <- df[df$CrossType == "Cross", ]
  parents <- df[df$CrossType == "Parent", ]
  
  # Check if we have enough data
  if(nrow(crosses) < 2) {
    return(list(error = "Insufficient cross data for analysis"))
  }
  
  if(nrow(parents) < 2) {
    return(list(error = "Insufficient parent data for analysis"))
  }
  
  # Check if trait exists and has valid data
  if(!trait %in% colnames(crosses) || all(is.na(crosses[[trait]]))) {
    return(list(error = paste("Trait", trait, "has no valid data")))
  }
  
  # Prepare data for diallel analysis
  diallel_df <- crosses[, c("Parent1", "Parent2", trait)]
  colnames(diallel_df) <- c("Parent1", "Parent2", "Trait")
  diallel_df <- diallel_df[complete.cases(diallel_df), ]
  
  if(nrow(diallel_df) < 2) {
    return(list(error = "Insufficient complete cases for analysis"))
  }
  
  # Get unique parents from crosses
  unique_parents <- unique(c(diallel_df$Parent1, diallel_df$Parent2))
  
  # Perform proper diallel ANOVA
  anova_result <- tryCatch({
    # Create a balanced dataset for ANOVA
    diallel_df$Parent1 <- factor(diallel_df$Parent1, levels = unique_parents)
    diallel_df$Parent2 <- factor(diallel_df$Parent2, levels = unique_parents)
    
    # Fit linear model for diallel analysis
    model <- lm(Trait ~ Parent1 + Parent2, data = diallel_df)
    anova_table <- anova(model)
    
    # Calculate sums of squares
    ss_parent1 <- anova_table["Parent1", "Sum Sq"]
    ss_parent2 <- anova_table["Parent2", "Sum Sq"]
    ss_residual <- anova_table["Residuals", "Sum Sq"]
    ss_total <- sum(anova_table$"Sum Sq")
    
    # Calculate mean squares
    ms_parent1 <- anova_table["Parent1", "Mean Sq"]
    ms_parent2 <- anova_table["Parent2", "Mean Sq"]
    ms_residual <- anova_table["Residuals", "Mean Sq"]
    
    # Calculate F-values
    f_parent1 <- anova_table["Parent1", "F value"]
    f_parent2 <- anova_table["Parent2", "F value"]
    
    # Calculate p-values
    p_parent1 <- anova_table["Parent1", "Pr(>F)"]
    p_parent2 <- anova_table["Parent2", "Pr(>F)"]
    
    list(
      table = anova_table,
      ss_gca = ss_parent1 + ss_parent2,
      ss_total = ss_total,
      ms_gca = (ms_parent1 + ms_parent2) / 2,
      ms_error = ms_residual
    )
  }, error = function(e) {
    return(list(error = paste("ANOVA failed:", e$message)))
  })
  
  # Calculate GCA effects
  parent_means <- sapply(unique_parents, function(p) {
    p_data <- c(diallel_df$Trait[diallel_df$Parent1 == p],
                diallel_df$Trait[diallel_df$Parent2 == p])
    if(length(p_data) > 0) {
      mean(p_data, na.rm = TRUE)
    } else {
      NA
    }
  })
  
  overall_mean <- mean(diallel_df$Trait, na.rm = TRUE)
  gca_effects <- parent_means - overall_mean
  
  gca_df <- data.frame(
    Parent = names(gca_effects),
    GCA = gca_effects,
    Mean_Performance = parent_means,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(GCA))
  
  # Calculate SCA effects with proper matrix handling
  sca_matrix <- matrix(NA, nrow = length(unique_parents), ncol = length(unique_parents))
  rownames(sca_matrix) <- colnames(sca_matrix) <- unique_parents
  
  # Fill the matrix with symmetric values
  for(i in 1:nrow(diallel_df)) {
    p1 <- as.character(diallel_df$Parent1[i])
    p2 <- as.character(diallel_df$Parent2[i])
    trait_value <- diallel_df$Trait[i]
    
    if(!is.na(trait_value) && p1 %in% names(gca_effects) && p2 %in% names(gca_effects)) {
      expected <- overall_mean + gca_effects[p1] + gca_effects[p2]
      sca_value <- trait_value - expected
      
      # Set both positions (symmetric matrix)
      sca_matrix[p1, p2] <- sca_value
      sca_matrix[p2, p1] <- sca_value
    }
  }
  
  # Ensure it's a proper numeric matrix
  sca_matrix <- matrix(as.numeric(sca_matrix), 
                       nrow = nrow(sca_matrix),
                       ncol = ncol(sca_matrix),
                       dimnames = dimnames(sca_matrix))
  
  return(list(
    anova = anova_result,
    gca = gca_df,
    sca = sca_matrix,
    overall_mean = overall_mean,
    n_crosses = nrow(diallel_df),
    n_parents = length(unique_parents)
  ))
}


#' Heterosis calculations with proper error handling
calculate_heterosis <- function(diallel_data, cross_data, parent_data, trait) {
  # Check if we have the necessary data
  if(nrow(parent_data) == 0 || nrow(cross_data) == 0) {
    return(NULL)
  }
  
  # Calculate parent means for the trait
  parent_means <- setNames(parent_data[[trait]], parent_data$ACC.)
  
  # Remove parents with missing trait values
  parent_means <- parent_means[!is.na(parent_means)]
  
  if(length(parent_means) == 0) {
    return(NULL)
  }
  
  # Calculate heterosis for each cross
  heterosis_df <- cross_data %>%
    rowwise() %>%
    mutate(
      Parent1_Mean = ifelse(Parent1 %in% names(parent_means), parent_means[Parent1], NA),
      Parent2_Mean = ifelse(Parent2 %in% names(parent_means), parent_means[Parent2], NA),
      Cross_Value = !!sym(trait),
      Mid_Parent = ifelse(!is.na(Parent1_Mean) && !is.na(Parent2_Mean), 
                          (Parent1_Mean + Parent2_Mean) / 2, NA),
      Better_Parent = ifelse(!is.na(Parent1_Mean) && !is.na(Parent2_Mean), 
                             max(Parent1_Mean, Parent2_Mean, na.rm = TRUE), NA),
      MPH = ifelse(!is.na(Cross_Value) && !is.na(Mid_Parent) && Mid_Parent != 0,
                   (Cross_Value - Mid_Parent) / Mid_Parent * 100, NA),
      BPH = ifelse(!is.na(Cross_Value) && !is.na(Better_Parent) && Better_Parent != 0,
                   (Cross_Value - Better_Parent) / Better_Parent * 100, NA)
    ) %>%
    select(CrossID, Parent1, Parent2, Cross_Value, MPH, BPH, Mid_Parent, Better_Parent) %>%
    filter(!is.na(MPH) | !is.na(BPH))  # Remove rows with no heterosis calculation
  
  return(heterosis_df)
}



###############################################################################
# SECTION 9: REPORT GENERATION FUNCTIONS
###############################################################################

#' Capture input parameters for report generation
#' @param input Shiny input object
#' @param values Reactive values
#' @return List of input parameters
capture_input_parameters <- function(input, values) {
  tryCatch({
    params <- list()
    
    # Data upload parameters
    params$data_upload <- list(
      pheno_file = if(!is.null(input$pheno_file)) input$pheno_file$name else "Not uploaded",
      geno_file = if(!is.null(input$geno_file)) input$geno_file$name else "Not uploaded",
      metadata_file = if(!is.null(input$metadata_file)) input$metadata_file$name else "Not uploaded",
      diallel_file = if(!is.null(input$diallel_file)) input$diallel_file$name else "Not uploaded",
      start_col = input$start_col,
      skip_lines = input$skip_lines
    )
    
    # Filtering parameters
    params$filtering <- list(
      callrate_threshold = input$callrate_threshold,
      maf_threshold = input$maf_threshold,
      het_threshold = input$het_threshold
    )
    
    # GWAS parameters
    params$gwas <- list(
      traits_analyzed = input$traits,
      maf_threshold_gwas = input$maf_threshold,
      missing_threshold = input$missing_threshold
    )
    
    # Multi-trait parameters
    params$multi_trait <- list(
      p_threshold = input$p_threshold,
      min_traits = input$min_traits,
      top_n_snps = input$top_n
    )
    
    # Visualization parameters
    params$visualization <- list(
      manhattan_trait = input$manhattan_trait,
      qq_trait = input$qq_trait,
      association_trait = input$association_trait
    )
    
    # Diallel analysis parameters
    params$diallel <- list(
      diallel_trait = input$diallel_trait,
      heterosis_type = input$heterosis_type,
      heterosis_threshold = input$heterosis_threshold
    )
    
    # MET analysis parameters
    params$met <- list(
      env_col = input$met_env_col,
      gen_col = input$met_gen_col,
      rep_col = input$met_rep_col,
      trait_vars = input$met_trait_vars
    )
    
    # Report parameters
    params$report <- list(
      title = input$report_title,
      type = input$report_type,
      format = input$report_format,
      sections = input$report_sections
    )
    
    # Data statistics - Comprehensive data capture
    params$data_stats <- list()
    
    if(!is.null(values$pheno_data)) {
      params$data_stats$pheno <- list(
        n_genotypes = nrow(values$pheno_data),
        n_traits = ncol(values$pheno_data) - 1,
        trait_names = setdiff(colnames(values$pheno_data), "Genotype"),
        missing_values = sum(is.na(values$pheno_data)),
        data_dimensions = dim(values$pheno_data)
      )
    }
    
    if(!is.null(values$filtered_gl)) {
      params$data_stats$geno <- list(
        n_snps = nLoc(values$filtered_gl),
        n_individuals = nInd(values$filtered_gl),
        original_snps = if(!is.null(values$gl_object)) nLoc(values$gl_object) else NA,
        retention_rate = if(!is.null(values$gl_object)) round(nLoc(values$filtered_gl)/nLoc(values$gl_object)*100, 1) else NA
      )
    }
    
    if(!is.null(values$metadata)) {
      params$data_stats$metadata <- list(
        n_samples = nrow(values$metadata),
        n_columns = ncol(values$metadata),
        column_names = colnames(values$metadata)
      )
    }
    
    if(!is.null(values$diallel_data)) {
      params$data_stats$diallel <- list(
        total_entries = nrow(values$diallel_data),
        n_parents = sum(values$diallel_data$CrossType == "Parent", na.rm = TRUE),
        n_crosses = sum(values$diallel_data$CrossType == "Cross", na.rm = TRUE),
        traits_available = colnames(values$diallel_data)[sapply(values$diallel_data, is.numeric)]
      )
    }
    
    # Analysis results summary
    if(!is.null(values$gwas_results)) {
      params$analysis_summary <- list(
        n_traits_analyzed = length(values$gwas_results$gwas_results),
        total_snps_tested = ncol(values$gwas_results$geno_matrix),
        n_samples_gwas = nrow(values$gwas_results$pheno_data),
        analysis_date = format(Sys.time(), '%Y-%m-%d %H:%M:%S')
      )
    }
    
    return(params)
  }, error = function(e) {
    cat("Error in capture_input_parameters:", e$message, "\n")
    return(list())
  })
}

#' Capture output results for report generation
#' @param values Reactive values
#' @param input Shiny input object
#' @return List of output results
capture_output_results <- function(values, input) {
  tryCatch({
    results <- list()
    
    # GWAS Results
    if(!is.null(values$gwas_results)) {
      # GWAS summary
      results$gwas_summary <- values$gwas_results$summary
      
      # Detailed results for each trait
      results$gwas_detailed <- list()
      for(trait in names(values$gwas_results$gwas_results)) {
        trait_data <- values$gwas_results$gwas_results[[trait]]
        results$gwas_detailed[[trait]] <- list(
          n_snps = nrow(trait_data),
          min_p_value = min(trait_data$P_value, na.rm = TRUE),
          significant_snps = sum(trait_data$P_adjusted < 0.05, na.rm = TRUE),
          top_snps = head(trait_data[order(trait_data$P_value), c("SNP_ID", "Chromosome", "Position", "P_value")], 5)
        )
      }
      
      # Combined results
      if(!is.null(values$combined_results)) {
        results$combined_results <- list(
          n_snps = nrow(values$combined_results),
          n_traits = ncol(values$combined_results) - 3,
          trait_names = setdiff(colnames(values$combined_results), c("SNP", "Chromosome", "Position")),
          data_dimensions = dim(values$combined_results)
        )
      }
    }
    
    # Multi-trait Analysis Results
    if(!is.null(values$merged_pvalues)) {
      results$multi_trait <- list(
        total_snps = nrow(values$merged_pvalues),
        trait_columns = setdiff(colnames(values$merged_pvalues), c("SNP_ID", "Chromosome", "Position"))
      )
      
      # Calculate multi-trait significance
      if(!is.null(input) && !is.null(input$p_threshold) && !is.null(input$min_traits)) {
        trait_cols <- setdiff(colnames(values$merged_pvalues), c("SNP_ID", "Chromosome", "Position"))
        if(length(trait_cols) > 0) {
          significant_counts <- rowSums(values$merged_pvalues[, trait_cols, drop = FALSE] < input$p_threshold, na.rm = TRUE)
          results$multi_trait$multi_trait_snps <- sum(significant_counts >= input$min_traits)
          results$multi_trait$significance_threshold <- input$p_threshold
          results$multi_trait$min_traits_required <- input$min_traits
        }
      }
    }
    
    # Association Plots Information
    if(!is.null(values$association_plot_paths)) {
      results$visualization <- list(
        manhattan_plots_generated = length(values$association_plot_paths$manhattan),
        qq_plots_generated = length(values$association_plot_paths$qq),
        multi_trait_plots_generated = length(values$association_plot_paths$multi_trait),
        output_directory = values$association_plot_paths$output_dir
      )
    }
    
    # CMPlots Information
    if(!is.null(values$cmplot_paths)) {
      results$cmplots <- list(
        manhattan_plot = ifelse(!is.null(values$cmplot_paths$manhattan), 
                                file.exists(values$cmplot_paths$manhattan), FALSE),
        qq_plot = ifelse(!is.null(values$cmplot_paths$qq), 
                         file.exists(values$cmplot_paths$qq), FALSE),
        circular_plot = ifelse(!is.null(values$cmplot_paths$circular), 
                               file.exists(values$cmplot_paths$circular), FALSE)
      )
    }
    
    # Diallel Analysis Results
    if(!is.null(values$diallel_results) && !is.null(input) && !is.null(input$diallel_trait)) {
      results$diallel <- list(
        trait_analyzed = input$diallel_trait,
        n_parents = ifelse(!is.null(values$diallel_results$n_parents), 
                           values$diallel_results$n_parents, NA),
        n_crosses = ifelse(!is.null(values$diallel_results$n_crosses), 
                           values$diallel_results$n_crosses, NA),
        overall_mean = ifelse(!is.null(values$diallel_results$overall_mean), 
                              values$diallel_results$overall_mean, NA),
        has_gca = !is.null(values$diallel_results$gca) && nrow(values$diallel_results$gca) > 0,
        has_sca = !is.null(values$diallel_results$sca)
      )
      
      if(!is.null(values$diallel_heterosis_data)) {
        results$diallel$heterosis <- list(
          n_crosses_with_heterosis = nrow(values$diallel_heterosis_data),
          mean_mph = mean(values$diallel_heterosis_data$MPH, na.rm = TRUE),
          mean_bph = mean(values$diallel_heterosis_data$BPH, na.rm = TRUE),
          high_heterosis_crosses = sum(values$diallel_heterosis_data$MPH > 20, na.rm = TRUE)
        )
      }
    }
    
    # MET Analysis Results
    if(!is.null(values$met_analysis_results)) {
      results$met <- list(
        has_desc_stats = !is.null(values$met_analysis_results$desc_stats),
        has_env_stats = !is.null(values$met_analysis_results$env_stats),
        n_environments = if(!is.null(values$met_data)) length(unique(values$met_data$ENV)) else NA,
        n_genotypes = if(!is.null(values$met_data)) length(unique(values$met_data$GEN)) else NA
      )
    }
    
    # File Outputs
    results$output_files <- list()
    
    if(dir.exists("shiny_gwas_results")) {
      results$output_files$gwas_results_dir <- list.files("shiny_gwas_results", full.names = FALSE)
    } else {
      results$output_files$gwas_results_dir <- "Not generated"
    }
    
    if(dir.exists("association_plots_combined")) {
      results$output_files$association_plots_dir <- list.files("association_plots_combined", full.names = FALSE)
    } else {
      results$output_files$association_plots_dir <- "Not generated"
    }
    
    if(dir.exists("combined_results_cmplots")) {
      results$output_files$cmplots_dir <- list.files("combined_results_cmplots", full.names = FALSE)
    } else {
      results$output_files$cmplots_dir <- "Not generated"
    }
    
    return(results)
  }, error = function(e) {
    cat("Error in capture_output_results:", e$message, "\n")
    return(list())
  })
}

#' Generate comprehensive report
#' @param input Shiny input object
#' @param params Input parameters
#' @param results Output results
#' @param output_file Output file path
#' @return Path to generated report
generate_comprehensive_report <- function(input, params, results, output_file) {
  tryCatch({
    # Create enhanced Rmd content
    rmd_content <- paste0(
      "---\n",
      "title: \"", input$report_title, "\"\n",
      "output: ", ifelse(input$report_format == "html", "html_document", "pdf_document"), "\n",
      "date: \"", format(Sys.time(), '%Y-%m-%d'), "\"\n",
      "---\n\n",
      "```{r setup, include=FALSE}\n",
      "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center')\n",
      "library(knitr)\n",
      "library(ggplot2)\n",
      "library(dplyr)\n",
      "```\n\n",
      "# Comprehensive GWAS Analysis Report\n\n",
      "## Executive Summary\n\n",
      "This report presents a comprehensive analysis of genome-wide association study (GWAS) conducted on **", 
      ifelse(!is.null(params$data_stats$pheno$n_genotypes), params$data_stats$pheno$n_genotypes, "N/A"), 
      " genotypes** with **", 
      ifelse(!is.null(params$data_stats$pheno$n_traits), params$data_stats$pheno$n_traits, "N/A"), 
      " traits** using **", 
      ifelse(!is.null(params$data_stats$geno$n_snps), params$data_stats$geno$n_snps, "N/A"), 
      " SNPs**.\n\n",
      "### Key Findings\n\n",
      "```{r key-findings, echo=FALSE}\n",
      "cat(\"\\n\")\n",
      "if(!is.null(results$gwas_summary)) {\n",
      "  total_significant <- sum(results$gwas_summary$Significant_SNPs, na.rm = TRUE)\n",
      "  cat(\"- **Total Significant Associations:** \", total_significant, \"\\n\")\n",
      "  \n",
      "  best_trait <- results$gwas_summary[which.max(results$gwas_summary$Significant_SNPs), ]\n",
      "  cat(\"- **Most Significant Trait:** \", best_trait$Trait, \" (\", best_trait$Significant_SNPs, \" SNPs)\\n\")\n",
      "  \n",
      "  min_p <- min(results$gwas_summary$Min_P_value, na.rm = TRUE)\n",
      "  cat(\"- **Most Significant P-value:** \", format(min_p, scientific = TRUE, digits = 3), \"\\n\")\n",
      "}\n",
      "\n",
      "if(!is.null(results$multi_trait) && !is.null(results$multi_trait$multi_trait_snps)) {\n",
      "  cat(\"- **Multi-trait Significant SNPs:** \", results$multi_trait$multi_trait_snps, \"\\n\")\n",
      "}\n",
      "\n",
      "if(!is.null(results$diallel)) {\n",
      "  cat(\"- **Diallel Analysis Completed for:** \", results$diallel$trait_analyzed, \"\\n\")\n",
      "}\n",
      "cat(\"\\n\")\n",
      "```\n\n",
      "## Detailed Analysis Results\n\n",
      "### GWAS Results Summary\n\n",
      "```{r gwas-summary, echo=FALSE}\n",
      "if(!is.null(results$gwas_summary)) {\n",
      "  knitr::kable(results$gwas_summary, caption = \"GWAS Results Summary by Trait\")\n",
      "} else {\n",
      "  cat(\"No GWAS results available.\\n\")\n",
      "}\n",
      "```\n\n",
      "### Data Quality Metrics\n\n",
      "```{r data-quality, echo=FALSE}\n",
      "cat(\"#### Phenotypic Data\\n\")\n",
      "if(!is.null(params$data_stats$pheno)) {\n",
      "  cat(\"- Samples:\", params$data_stats$pheno$n_genotypes, \"\\n\")\n",
      "  cat(\"- Traits:\", params$data_stats$pheno$n_traits, \"\\n\")\n",
      "  cat(\"- Missing Values:\", params$data_stats$pheno$missing_values, \"\\n\")\n",
      "}\n",
      "\n",
      "cat(\"#### Genotypic Data\\n\")\n",
      "if(!is.null(params$data_stats$geno)) {\n",
      "  cat(\"- SNPs after QC:\", params$data_stats$geno$n_snps, \"\\n\")\n",
      "  cat(\"- Individuals:\", params$data_stats$geno$n_individuals, \"\\n\")\n",
      "  if(!is.na(params$data_stats$geno$retention_rate)) {\n",
      "    cat(\"- SNP Retention Rate:\", params$data_stats$geno$retention_rate, \"%\\n\")\n",
      "  }\n",
      "}\n",
      "```\n\n",
      "### Analysis Parameters\n\n",
      "```{r parameters, echo=FALSE}\n",
      "cat(\"#### Quality Control Parameters\\n\")\n",
      "cat(\"- Call Rate Threshold:\", params$filtering$callrate_threshold, \"\\n\")\n",
      "cat(\"- MAF Threshold:\", params$filtering$maf_threshold, \"\\n\")\n",
      "cat(\"- Heterozygosity Threshold:\", params$filtering$het_threshold, \"\\n\")\n",
      "\n",
      "cat(\"#### GWAS Parameters\\n\")\n",
      "cat(\"- Traits Analyzed:\", paste(params$gwas$traits_analyzed, collapse = \", \"), \"\\n\")\n",
      "cat(\"- MAF Threshold:\", params$gwas$maf_threshold_gwas, \"\\n\")\n",
      "cat(\"- Missing Data Threshold:\", params$gwas$missing_threshold, \"\\n\")\n",
      "```\n\n",
      "## Visualizations\n\n",
      "The following visualizations were generated during the analysis:\n\n",
      "```{r visualizations, echo=FALSE}\n",
      "if(!is.null(results$visualization)) {\n",
      "  cat(\"- Manhattan plots:\", results$visualization$manhattan_plots_generated, \"traits\\n\")\n",
      "  cat(\"- QQ plots:\", results$visualization$qq_plots_generated, \"traits\\n\")\n",
      "  cat(\"- Multi-trait association plots: Generated\\n\")\n",
      "}\n",
      "\n",
      "if(!is.null(results$cmplots)) {\n",
      "  cat(\"- CMPlots: Manhattan, QQ, and Circular plots generated\\n\")\n",
      "}\n",
      "```\n\n",
      "## Diallel Analysis Results\n\n",
      "```{r diallel-results, echo=FALSE}\n",
      "if(!is.null(results$diallel)) {\n",
      "  cat(\"### Diallel Analysis Summary\\n\")\n",
      "  cat(\"- Trait Analyzed:\", results$diallel$trait_analyzed, \"\\n\")\n",
      "  cat(\"- Number of Parents:\", results$diallel$n_parents, \"\\n\")\n",
      "  cat(\"- Number of Crosses:\", results$diallel$n_crosses, \"\\n\")\n",
      "  \n",
      "  if(!is.null(results$diallel$heterosis)) {\n",
      "    cat(\"\\n### Heterosis Analysis\\n\")\n",
      "    cat(\"- Crosses with heterosis data:\", results$diallel$heterosis$n_crosses_with_heterosis, \"\\n\")\n",
      "    cat(\"- Mean Mid-Parent Heterosis (MPH):\", round(results$diallel$heterosis$mean_mph, 2), \"%\\n\")\n",
      "    cat(\"- Mean Better-Parent Heterosis (BPH):\", round(results$diallel$heterosis$mean_bph, 2), \"%\\n\")\n",
      "    cat(\"- Crosses with high heterosis (>20%):\", results$diallel$heterosis$high_heterosis_crosses, \"\\n\")\n",
      "  }\n",
      "}\n",
      "```\n\n",
      "## MET Analysis Results\n\n",
      "```{r met-results, echo=FALSE}\n",
      "if(!is.null(results$met)) {\n",
      "  cat(\"### Multi-Environment Trial Analysis\\n\")\n",
      "  cat(\"- Number of Environments:\", results$met$n_environments, \"\\n\")\n",
      "  cat(\"- Number of Genotypes:\", results$met$n_genotypes, \"\\n\")\n",
      "  cat(\"- Descriptive Statistics:\", ifelse(results$met$has_desc_stats, \"Available\", \"Not available\"), \"\\n\")\n",
      "  cat(\"- Environment Statistics:\", ifelse(results$met$has_env_stats, \"Available\", \"Not available\"), \"\\n\")\n",
      "}\n",
      "```\n\n",
      "## Conclusion\n\n",
      "This comprehensive analysis provides valuable insights into the genetic architecture of the studied traits. \n",
      "The results can be used for:\n\n",
      "- **Marker-assisted selection** in breeding programs\n",
      "- **Candidate gene identification** for functional studies\n",
      "- **Understanding genetic architecture** of complex traits\n",
      "- **Multi-trait selection** for genetic improvement\n\n",
      "***\n",
      "*Report generated using GWAS Analysis Platform v2.0 on `", Sys.Date(), "`*\n"
    )
    
    # Write Rmd file
    writeLines(rmd_content, output_file)
    
    return(output_file)
  }, error = function(e) {
    cat("Error in generate_comprehensive_report:", e$message, "\n")
    return(NULL)
  })
}


#==============================================================================
# CREATE DIFFERENT DATAFRAME FORMATS
#==============================================================================

# Create unified dataframe with all columns from all traits
create_unified_gwas_df <- function(trait_gwas_results) {
  if (is.null(trait_gwas_results) || length(trait_gwas_results) == 0) {
    return(NULL)
  }
  
  # Add trait name to each dataframe
  trait_dfs <- list()
  for (trait_name in names(trait_gwas_results)) {
    if (!is.null(trait_gwas_results[[trait_name]])) {
      df <- trait_gwas_results[[trait_name]]
      df$Trait <- trait_name  # Add trait identifier
      trait_dfs[[trait_name]] <- df
    }
  }
  
  # Combine all dataframes
  unified_df <- do.call(rbind, trait_dfs)
  rownames(unified_df) <- NULL
  
  return(unified_df)
}

# Usage in your server:
#unified_gwas <- create_unified_gwas_df(values$multi_trait_results$trait_gwas_results)


create_wide_gwas_format <- function(trait_gwas_results) {
  require(dplyr)
  require(tidyr)
  
  if (is.null(trait_gwas_results)) return(NULL)
  
  # First create long format
  long_df <- create_unified_gwas_tidy(trait_gwas_results)
  
  if (is.null(long_df)) return(NULL)
  
  # Create wide format for P-values
  pvalue_wide <- long_df %>%
    select(SNP_ID, Chromosome, Position, Trait, P_value) %>%
    pivot_wider(
      id_cols = c(SNP_ID, Chromosome, Position),
      names_from = Trait,
      values_from = P_value,
      names_prefix = "P_"
    )
  
  # Create wide format for Effects
  effect_wide <- long_df %>%
    select(SNP_ID, Trait, Effect) %>%
    pivot_wider(
      id_cols = SNP_ID,
      names_from = Trait,
      values_from = Effect,
      names_prefix = "Effect_"
    )
  
  # Combine wide formats
  wide_df <- pvalue_wide %>%
    left_join(effect_wide, by = "SNP_ID")
  
  return(wide_df)
}

# Usage:
#values$multi_trait_wide <- create_wide_gwas_format(values$multi_trait_results$trait_gwas_results)



create_enhanced_unified_df <- function(trait_gwas_results, trait_correlations = NULL) {
  require(dplyr)
  
  if (is.null(trait_gwas_results)) return(NULL)
  
  # Create basic unified dataframe
  unified <- create_unified_gwas_tidy(trait_gwas_results)
  
  if (is.null(unified)) return(NULL)
  
  # Add significance flags
  unified <- unified %>%
    mutate(
      Significant_0_05 = P_value < 0.05,
      Significant_0_01 = P_value < 0.01,
      Significant_0_001 = P_value < 0.001,
      Significant_Bonferroni = P_adjusted < 0.05,
      LOD = -log10(P_value)
    )
  
  # Add trait correlation information if available
  if (!is.null(trait_correlations)) {
    # Calculate which traits each SNP affects
    trait_affected <- unified %>%
      filter(Significant_0_05) %>%
      group_by(SNP_ID) %>%
      summarise(
        Traits_Affected = paste(sort(unique(Trait)), collapse = ";"),
        N_Traits_Affected = n_distinct(Trait),
        .groups = "drop"
      )
    
    unified <- unified %>%
      left_join(trait_affected, by = "SNP_ID")
  }
  
  # Add SNP rank within each trait
  unified <- unified %>%
    group_by(Trait) %>%
    mutate(
      Rank = rank(P_value),
      Percentile = rank(P_value) / n() * 100
    ) %>%
    ungroup()
  
  return(unified)
}

# Usage:
#values$multi_trait_enhanced <- create_enhanced_unified_df(
#  values$multi_trait_results$trait_gwas_results,
#  values$multi_trait_results$trait_correlations
#)


create_complete_integrated_df <- function(multi_trait_results) {
  require(dplyr)
  
  # Check if we have all necessary components
  if (is.null(multi_trait_results) || 
      is.null(multi_trait_results$trait_gwas_results)) {
    return(NULL)
  }
  
  # Create unified trait results
  trait_unified <- create_enhanced_unified_df(
    multi_trait_results$trait_gwas_results,
    multi_trait_results$trait_correlations
  )
  
  if (is.null(trait_unified)) return(NULL)
  
  # If we have combined results, merge them
  if (!is.null(multi_trait_results$combined_results)) {
    # Aggregate trait-specific information
    trait_summary <- trait_unified %>%
      group_by(SNP_ID) %>%
      summarise(
        All_Traits = paste(sort(unique(Trait)), collapse = ";"),
        Significant_Traits = paste(sort(unique(Trait[Significant_0_05])), collapse = ";"),
        Min_P_Value_Trait = Trait[which.min(P_value)][1],
        Min_P_Value = min(P_value, na.rm = TRUE),
        Max_Effect_Trait = Trait[which.max(abs(Effect))][1],
        Max_Effect = Effect[which.max(abs(Effect))][1],
        Mean_Effect = mean(Effect, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Merge with combined results
    complete_df <- multi_trait_results$combined_results %>%
      left_join(trait_summary, by = "SNP_ID")
    
  } else {
    complete_df <- trait_unified
  }
  
  # Add trait correlation statistics if available
  if (!is.null(multi_trait_results$trait_correlations)) {
    # Extract correlation matrix
    cor_matrix <- multi_trait_results$trait_correlations
    if (is.list(cor_matrix) && !is.null(cor_matrix$correlations)) {
      cor_matrix <- cor_matrix$correlations
    }
    
    # You can add correlation-based features here
    # For example: calculate average correlation of traits affected by each SNP
  }
  
  return(complete_df)
}



# Function to create combined cross-SNP dataframe with flanking markers
#' Create combined cross-SNP dataframe using app data structures
create_cross_snp_dataframe_single_trait <- function(matched_data = NULL, multi_trait_results = NULL, 
                                       processed_geno = NULL, diallel_data = NULL,
                                       trait = NULL, window_size = 1, p_threshold = 0.05) {
  
  tryCatch({
    cat("Creating cross-SNP dataframe using app data structures...\n")
    
    # Validate required data
    if (is.null(diallel_data)) {
      stop("Diallel data is required (values$diallel_data)")
    }
    
    # Extract cross data
    cross_data <- diallel_data[diallel_data$CrossType == "Cross", ]
    
    if (nrow(cross_data) == 0) {
      stop("No cross data found in diallel data")
    }
    
    cat("  Found", nrow(cross_data), "crosses\n")
    
    # Prepare cross summary
    cross_summary <- data.frame(
      Cross = cross_data$CrossID,
      Parent1 = cross_data$Parent1,
      Parent2 = cross_data$Parent2,
      stringsAsFactors = FALSE
    )
    
    # Get GWAS data
    if (!is.null(multi_trait_results)) {
      cat("  Using multi-trait results\n")
      
      if (is.null(trait)) {
        # Use first trait if none specified
        trait_names <- names(multi_trait_results$trait_gwas_results)
        if (length(trait_names) == 0) {
          stop("No traits found in multi_trait_results")
        }
        trait <- trait_names[1]
        cat("  Using first trait:", trait, "\n")
      }
      
      # Get GWAS results for the trait
      if (trait %in% names(multi_trait_results$trait_gwas_results)) {
        gwas_df <- multi_trait_results$trait_gwas_results[[trait]]
      } else {
        stop("Trait '", trait, "' not found in multi_trait_results")
      }
    } else {
      stop("GWAS results are required (values$multi_trait_results)")
    }
    
    if (is.null(gwas_df) || nrow(gwas_df) == 0) {
      stop("No GWAS results available for trait ", trait)
    }
    
    cat("  GWAS data:", nrow(gwas_df), "SNPs for trait:", trait, "\n")
    
    # Add flanking markers if processed_geno is available
    if (!is.null(processed_geno) && !is.null(processed_geno$snp_info)) {
      cat("  Adding flanking SNP markers...\n")
      cross_summary <- add_flanking_snp_markers(
        cross_summary, 
        processed_geno$snp_matrix, 
        processed_geno$snp_info, 
        window_size = window_size
      )
    }
    
    # Ensure required columns exist in GWAS results
    # Try to find SNP_ID column
    snp_id_cols <- c("SNP_ID", "SNP", "Marker", "snp", "marker")
    snp_id_col <- NULL
    for (col in snp_id_cols) {
      if (col %in% colnames(gwas_df)) {
        snp_id_col <- col
        break
      }
    }
    
    if (is.null(snp_id_col)) {
      # If no SNP_ID column, try to use row names
      gwas_df$SNP_ID <- rownames(gwas_df)
      snp_id_col <- "SNP_ID"
      cat("  Using row names as SNP_ID\n")
    }
    
    # Try to find chromosome column
    chrom_cols <- c("Chromosome", "CHR", "chr", "Chrom", "chrom")
    chrom_col <- NULL
    for (col in chrom_cols) {
      if (col %in% colnames(gwas_df)) {
        chrom_col <- col
        break
      }
    }
    
    # Try to find position column
    pos_cols <- c("Position", "POS", "pos", "BP", "bp", "location")
    pos_col <- NULL
    for (col in pos_cols) {
      if (col %in% colnames(gwas_df)) {
        pos_col <- col
        break
      }
    }
    
    # Try to find P-value column
    pval_cols <- c("P_value", "p.value", "P.value", "p", "P", "P.value")
    pval_col <- NULL
    for (col in pval_cols) {
      if (col %in% colnames(gwas_df)) {
        pval_col <- col
        break
      }
    }
    
    # Try to find effect column
    effect_cols <- c("Effect", "effect", "BETA", "beta", "Estimate", "estimate")
    effect_col <- NULL
    for (col in effect_cols) {
      if (col %in% colnames(gwas_df)) {
        effect_col <- col
        break
      }
    }
    
    # Rename columns for consistency
    colnames(gwas_df)[colnames(gwas_df) == snp_id_col] <- "SNP_ID"
    if (!is.null(chrom_col)) {
      colnames(gwas_df)[colnames(gwas_df) == chrom_col] <- "Chromosome"
    }
    if (!is.null(pos_col)) {
      colnames(gwas_df)[colnames(gwas_df) == pos_col] <- "Position"
    }
    if (!is.null(pval_col)) {
      colnames(gwas_df)[colnames(gwas_df) == pval_col] <- "P_value"
    }
    if (!is.null(effect_col)) {
      colnames(gwas_df)[colnames(gwas_df) == effect_col] <- "Effect"
    }
    
    # Calculate LOD if not present
    if (!"LOD" %in% colnames(gwas_df) && "P_value" %in% colnames(gwas_df)) {
      gwas_df$LOD <- -log10(gwas_df$P_value)
    }
    
    # Add PVE if available, otherwise set to NA
    if (!"PVE" %in% colnames(gwas_df)) {
      gwas_df$PVE <- NA
    }
    
    # Add significance flag
    gwas_df$Significant <- gwas_df$P_value < p_threshold
    
    # Create all combinations of crosses and SNPs
    cat("  Creating cross-SNP combinations...\n")
    cross_snp_combinations <- expand.grid(
      Cross = cross_summary$Cross,
      SNP_ID = unique(gwas_df$SNP_ID),
      stringsAsFactors = FALSE
    )
    
    # Merge with cross summary
    cross_snp_df <- merge(cross_snp_combinations, cross_summary, by = "Cross", all.x = TRUE)
    
    # Merge with GWAS data
    cross_snp_df <- merge(cross_snp_df, gwas_df, by = "SNP_ID", all.x = TRUE)
    
    # Add trait information
    cross_snp_df$Trait <- trait
    
    # Get parent genotypes if genotype matrix is available
    if (!is.null(processed_geno) && !is.null(processed_geno$snp_matrix)) {
      cat("  Extracting parent genotypes...\n")
      
      # Initialize genotype columns
      cross_snp_df$Parent1_Genotype <- NA
      cross_snp_df$Parent2_Genotype <- NA
      cross_snp_df$Expected_Genotype <- NA
      
      # Get genotype matrix
      geno_matrix <- processed_geno$snp_matrix
      
      # Process each unique parent
      unique_parents <- unique(c(cross_summary$Parent1, cross_summary$Parent2))
      
      for (parent in unique_parents) {
        # Clean parent name for matching
        parent_clean <- trimws(as.character(parent))
        
        # Find matching rows in genotype matrix
        matching_rows <- which(rownames(geno_matrix) == parent_clean)
        if (length(matching_rows) == 0) {
          # Try case-insensitive matching
          matching_rows <- which(tolower(rownames(geno_matrix)) == tolower(parent_clean))
        }
        
        if (length(matching_rows) > 0) {
          parent_row <- matching_rows[1]
          
          # Find rows for this parent as Parent1
          parent1_rows <- which(cross_snp_df$Parent1 == parent)
          if (length(parent1_rows) > 0) {
            # Get SNP IDs for these rows
            snp_ids <- cross_snp_df$SNP_ID[parent1_rows]
            
            # Find column indices for these SNPs
            snp_indices <- match(snp_ids, colnames(geno_matrix))
            
            # Extract genotypes
            cross_snp_df$Parent1_Genotype[parent1_rows] <- geno_matrix[parent_row, snp_indices]
          }
          
          # Find rows for this parent as Parent2
          parent2_rows <- which(cross_snp_df$Parent2 == parent)
          if (length(parent2_rows) > 0) {
            # Get SNP IDs for these rows
            snp_ids <- cross_snp_df$SNP_ID[parent2_rows]
            
            # Find column indices for these SNPs
            snp_indices <- match(snp_ids, colnames(geno_matrix))
            
            # Extract genotypes
            cross_snp_df$Parent2_Genotype[parent2_rows] <- geno_matrix[parent_row, snp_indices]
          }
        }
      }
      
      # Calculate expected genotype (average of parents)
      cross_snp_df$Expected_Genotype <- rowMeans(
        cross_snp_df[, c("Parent1_Genotype", "Parent2_Genotype")], 
        na.rm = TRUE
      )
    }
    
    # Reorder columns
    final_cols <- c("Cross", "Parent1", "Parent2", "SNP_ID", "Chromosome", "Position",
                    "Trait", "P_value", "LOD", "PVE", "Effect", 
                    "Parent1_Genotype", "Parent2_Genotype", "Expected_Genotype", "Significant")
    
    # Add flanking marker columns if available
    if ("LeftFlankingSNPs" %in% colnames(cross_snp_df)) {
      final_cols <- c(final_cols, "LeftFlankingSNPs")
    }
    if ("RightFlankingSNPs" %in% colnames(cross_snp_df)) {
      final_cols <- c(final_cols, "RightFlankingSNPs")
    }
    
    # Select only columns that exist
    existing_cols <- final_cols[final_cols %in% colnames(cross_snp_df)]
    cross_snp_df <- cross_snp_df[, existing_cols]
    
    # Sort by Cross, then by Chromosome and Position
    if ("Chromosome" %in% colnames(cross_snp_df) && "Position" %in% colnames(cross_snp_df)) {
      # Convert chromosomes to numeric for proper sorting
      chrom_numeric <- suppressWarnings(as.numeric(cross_snp_df$Chromosome))
      if (any(is.na(chrom_numeric))) {
        # Keep non-numeric chromosomes at the end
        cross_snp_df$Chromosome_sort <- ifelse(is.na(chrom_numeric), 
                                               9999,  # Large number for non-numeric
                                               chrom_numeric)
        cross_snp_df <- cross_snp_df[order(cross_snp_df$Cross, 
                                           cross_snp_df$Chromosome_sort, 
                                           cross_snp_df$Position), ]
        cross_snp_df$Chromosome_sort <- NULL
      } else {
        cross_snp_df <- cross_snp_df[order(cross_snp_df$Cross, 
                                           cross_snp_df$Chromosome, 
                                           cross_snp_df$Position), ]
      }
    } else {
      cross_snp_df <- cross_snp_df[order(cross_snp_df$Cross), ]
    }
    
    # Reset row names
    rownames(cross_snp_df) <- NULL
    
    cat("  Created dataframe with", nrow(cross_snp_df), "rows\n")
    cat("  Crosses:", length(unique(cross_snp_df$Cross)), "\n")
    cat("  SNPs:", length(unique(cross_snp_df$SNP_ID)), "\n")
    cat("  Significant SNPs (p <", p_threshold, "):", 
        sum(cross_snp_df$Significant, na.rm = TRUE), "\n")
    
    return(cross_snp_df)
    
  }, error = function(e) {
    cat("Error creating cross-SNP dataframe:", e$message, "\n")
    return(NULL)
  })
}




#' Simplified version that skips genotype extraction if problematic
create_cross_snp_dataframe_simple <- function(matched_data = NULL, multi_trait_results = NULL, 
                                              processed_geno = NULL, diallel_data = NULL,
                                              p_threshold = 0.05, max_snps_per_trait = 500) {
  
  tryCatch({
    cat("\n=== Creating Simplified Cross-SNP DataFrame ===\n")
    
    # Validate required data
    if (is.null(diallel_data)) {
      stop("Diallel data is required")
    }
    
    if (is.null(multi_trait_results) || is.null(multi_trait_results$trait_gwas_results)) {
      stop("Multi-trait GWAS results are required")
    }
    
    # Extract cross data
    cross_data <- diallel_data[diallel_data$CrossType == "Cross", ]
    
    if (nrow(cross_data) == 0) {
      stop("No cross data found")
    }
    
    cat("Found", nrow(cross_data), "crosses\n")
    
    # Prepare cross summary
    cross_summary <- data.frame(
      Cross = cross_data$CrossID,
      Parent1 = cross_data$Parent1,
      Parent2 = cross_data$Parent2,
      stringsAsFactors = FALSE
    )
    
    # Get all trait names
    trait_names <- names(multi_trait_results$trait_gwas_results)
    
    # Initialize results
    all_results <- list()
    
    for (trait in trait_names) {
      cat("  Processing trait:", trait, "\n")
      
      gwas_df <- multi_trait_results$trait_gwas_results[[trait]]
      
      if (is.null(gwas_df) || nrow(gwas_df) == 0) {
        cat("    Skipping - no data\n")
        next
      }
      
      # Standardize columns
      gwas_df <- standardize_gwas_columns(gwas_df)
      
      # Calculate PVE if not present
      if (!"PVE" %in% colnames(gwas_df)) {
        gwas_df$PVE <- calculate_pve(gwas_df)
      }
      
      # Limit number of SNPs
      if (nrow(gwas_df) > max_snps_per_trait) {
        gwas_df <- gwas_df[order(gwas_df$P_value), ][1:max_snps_per_trait, ]
      }
      
      # Create combinations
      combos <- expand.grid(
        Cross = cross_summary$Cross,
        SNP_ID = gwas_df$SNP_ID,
        stringsAsFactors = FALSE
      )
      
      # Merge with cross summary
      trait_df <- merge(combos, cross_summary, by = "Cross", all.x = TRUE)
      
      # Merge with GWAS data
      trait_df <- merge(trait_df, gwas_df, by = "SNP_ID", all.x = TRUE)
      
      # Add trait information
      trait_df$Trait <- trait
      
      # Calculate LOD
      trait_df$LOD <- -log10(trait_df$P_value)
      
      # Add significance flag
      trait_df$Significant <- trait_df$P_value < p_threshold
      
      # Add genotype columns (empty for now)
      trait_df$Parent1_Genotype <- NA
      trait_df$Parent2_Genotype <- NA
      trait_df$Expected_Genotype <- NA
      
      all_results[[trait]] <- trait_df
    }
    
    # Combine all results
    if (length(all_results) == 0) {
      return(NULL)
    }
    
    combined_df <- do.call(rbind, all_results)
    rownames(combined_df) <- NULL
    
    # Sort
    combined_df <- combined_df[order(combined_df$Trait, combined_df$Cross, 
                                     combined_df$Chromosome, combined_df$Position), ]
    
    cat("\nSuccessfully created dataframe with", nrow(combined_df), "rows\n")
    
    return(combined_df)
    
  }, error = function(e) {
    cat("Error in simplified version:", e$message, "\n")
    return(NULL)
  })
}





#' Helper function to find column with multiple possible names
find_column <- function(df, possible_names, required = FALSE) {
  for (name in possible_names) {
    if (name %in% colnames(df)) {
      return(name)
    }
  }
  
  # Try case-insensitive matching
  lower_names <- tolower(possible_names)
  df_lower <- tolower(colnames(df))
  
  for (i in seq_along(lower_names)) {
    if (lower_names[i] %in% df_lower) {
      return(colnames(df)[which(df_lower == lower_names[i])[1]])
    }
  }
  
  if (required) {
    warning("Required column not found. Possible names: ", paste(possible_names, collapse = ", "))
  }
  
  return(NULL)
}




check_cross_snp_prerequisites <- function(values) {
  prerequisites <- list(
    diallel_data = !is.null(values$diallel_data),
    multi_trait_results = !is.null(values$multi_trait_results),
    processed_geno = !is.null(values$processed_geno),
    snp_matrix = !is.null(values$processed_geno$snp_matrix),
    snp_info = !is.null(values$processed_geno$snp_info)
  )
  
  missing <- names(prerequisites)[!unlist(prerequisites)]
  
  if (length(missing) > 0) {
    return(list(
      ready = FALSE,
      missing = missing,
      message = paste("Missing:", paste(missing, collapse = ", "))
    ))
  }
  
  return(list(ready = TRUE, missing = NULL, message = "All data available"))
}


# Usage:
#values$complete_integrated_df <- create_complete_integrated_df(values$multi_trait_results)


#==============================================================================
# CORRECTED HELPER FUNCTIONS
#==============================================================================

#' Safer function to standardize GWAS column names with better error handling
standardize_gwas_columns <- function(gwas_df) {
  tryCatch({
    if (is.null(gwas_df) || nrow(gwas_df) == 0) {
      return(NULL)
    }
    
    # Make a copy
    df <- as.data.frame(gwas_df)
    
    # Standardize SNP ID column
    snp_id_cols <- c("SNP_ID", "SNP", "Marker", "snp", "marker", "SNP.id", "ID", "snp_id")
    snp_id_col <- NULL
    
    for (col in snp_id_cols) {
      if (col %in% colnames(df)) {
        snp_id_col <- col
        break
      }
    }
    
    # Try case-insensitive
    if (is.null(snp_id_col)) {
      df_lower <- tolower(colnames(df))
      for (col in tolower(snp_id_cols)) {
        if (col %in% df_lower) {
          idx <- which(df_lower == col)[1]
          snp_id_col <- colnames(df)[idx]
          break
        }
      }
    }
    
    # Use row names as fallback
    if (is.null(snp_id_col) && !is.null(rownames(df))) {
      df$SNP_ID <- rownames(df)
    } else if (!is.null(snp_id_col)) {
      colnames(df)[colnames(df) == snp_id_col] <- "SNP_ID"
    } else {
      df$SNP_ID <- paste0("SNP_", 1:nrow(df))
    }
    
    # Clean SNP_IDs
    df$SNP_ID <- as.character(df$SNP_ID)
    df$SNP_ID <- trimws(df$SNP_ID)
    
    # Standardize chromosome column
    chrom_cols <- c("Chromosome", "CHR", "chr", "Chrom", "chrom", "CHROM", "chromosome")
    chrom_col <- NULL
    
    for (col in chrom_cols) {
      if (col %in% colnames(df)) {
        chrom_col <- col
        break
      }
    }
    
    if (!is.null(chrom_col)) {
      colnames(df)[colnames(df) == chrom_col] <- "Chromosome"
    } else {
      df$Chromosome <- "Unknown"
    }
    
    # Standardize position column
    pos_cols <- c("Position", "POS", "pos", "BP", "bp", "location", "POSITION", "position")
    pos_col <- NULL
    
    for (col in pos_cols) {
      if (col %in% colnames(df)) {
        pos_col <- col
        break
      }
    }
    
    if (!is.null(pos_col)) {
      colnames(df)[colnames(df) == pos_col] <- "Position"
    } else {
      df$Position <- 1:nrow(df)
    }
    
    # Standardize P-value column
    pval_cols <- c("P_value", "p.value", "P.value", "p", "P", "pvalue", "Pvalue", "P.VALUE", "P_value", "p_value")
    pval_col <- NULL
    
    for (col in pval_cols) {
      if (col %in% colnames(df)) {
        pval_col <- col
        break
      }
    }
    
    if (!is.null(pval_col)) {
      colnames(df)[colnames(df) == pval_col] <- "P_value"
      # Ensure numeric
      df$P_value <- as.numeric(as.character(df$P_value))
    } else {
      df$P_value <- NA
    }
    
    # Standardize effect column
    effect_cols <- c("Effect", "effect", "BETA", "beta", "Estimate", "estimate", "B", "beta_effect")
    effect_col <- NULL
    
    for (col in effect_cols) {
      if (col %in% colnames(df)) {
        effect_col <- col
        break
      }
    }
    
    if (!is.null(effect_col)) {
      colnames(df)[colnames(df) == effect_col] <- "Effect"
      df$Effect <- as.numeric(as.character(df$Effect))
    } else {
      df$Effect <- NA
    }
    
    # Standardize PVE column
    pve_cols <- c("PVE", "R_squared", "R2", "r.squared", "Rsq", "VE", "var.exp", "R_squared")
    pve_col <- NULL
    
    for (col in pve_cols) {
      if (col %in% colnames(df)) {
        pve_col <- col
        break
      }
    }
    
    if (!is.null(pve_col)) {
      colnames(df)[colnames(df) == pve_col] <- "PVE"
      df$PVE <- as.numeric(as.character(df$PVE))
    } else {
      # Calculate PVE if not present
      df$PVE <- calculate_pve_safe(df)
    }
    
    # Calculate LOD
    df$LOD <- -log10(df$P_value)
    
    return(df)
    
  }, error = function(e) {
    cat("Error in standardize_gwas_columns:", e$message, "\n")
    return(gwas_df)  # Return original if error
  })
}

#' Safe PVE calculation
calculate_pve_safe <- function(gwas_df) {
  tryCatch({
    # If R_squared is available, use it as PVE (convert to percentage)
    if ("R_squared" %in% colnames(gwas_df)) {
      pve <- as.numeric(as.character(gwas_df$R_squared)) * 100
      return(pve)
    }
    
    # If effect and SE are available, estimate PVE
    if ("Effect" %in% colnames(gwas_df) && "SE" %in% colnames(gwas_df)) {
      effect <- as.numeric(as.character(gwas_df$Effect))
      se <- as.numeric(as.character(gwas_df$SE))
      
      # Handle missing or infinite values
      effect[is.infinite(effect) | is.na(effect)] <- 0
      se[is.infinite(se) | is.na(se) | se == 0] <- 1
      
      z_score <- abs(effect / se)
      pve_est <- (z_score^2) / (z_score^2 + 100) * 100  # Assuming N=100 as default
      return(pve_est)
    }
    
    # Return NA if no method available
    return(rep(NA, nrow(gwas_df)))
    
  }, error = function(e) {
    cat("Error in calculate_pve_safe:", e$message, "\n")
    return(rep(NA, nrow(gwas_df)))
  })
}

#' Safe function to get common SNPs across traits
get_common_snps_safe <- function(trait_gwas_results) {
  tryCatch({
    if (is.null(trait_gwas_results) || length(trait_gwas_results) == 0) {
      return(character(0))
    }
    
    all_snp_lists <- list()
    
    for (trait in names(trait_gwas_results)) {
      gwas_df <- trait_gwas_results[[trait]]
      
      if (!is.null(gwas_df) && nrow(gwas_df) > 0) {
        # Standardize columns first
        gwas_df <- standardize_gwas_columns(gwas_df)
        
        if (!is.null(gwas_df) && "SNP_ID" %in% colnames(gwas_df)) {
          snps <- unique(as.character(gwas_df$SNP_ID))
          snps <- snps[!is.na(snps) & snps != ""]
          all_snp_lists[[trait]] <- snps
        }
      }
    }
    
    # Remove empty lists
    all_snp_lists <- all_snp_lists[sapply(all_snp_lists, length) > 0]
    
    if (length(all_snp_lists) == 0) {
      return(character(0))
    }
    
    # Find intersection
    common_snps <- Reduce(intersect, all_snp_lists)
    
    if (length(common_snps) == 0) {
      # If no common SNPs, use union
      cat("Warning: No common SNPs across all traits. Using union of all SNPs.\n")
      common_snps <- unique(unlist(all_snp_lists))
    }
    
    return(common_snps)
    
  }, error = function(e) {
    cat("Error in get_common_snps_safe:", e$message, "\n")
    return(character(0))
  })
}



#==============================================================================
# CORRECTED MAIN FUNCTION - FIXED VERSION
#==============================================================================

#' Create combined cross-SNP dataframe for ALL traits with robust error handling
create_cross_snp_dataframe_all_traits_fixed <- function(
    multi_trait_results = NULL, 
    diallel_data = NULL,
    processed_geno = NULL,
    p_threshold = 0.05,
    max_snps_per_trait = 1000,
    include_parent_genotypes = FALSE) {
  
  tryCatch({
    cat("\n=== Creating Cross-SNP DataFrame for ALL Traits ===\n")
    
    # Validate required data
    if (is.null(diallel_data)) {
      stop("Diallel data is required")
    }
    
    if (is.null(multi_trait_results) || is.null(multi_trait_results$trait_gwas_results)) {
      stop("Multi-trait GWAS results are required")
    }
    
    # Extract cross data
    cross_data <- diallel_data[diallel_data$CrossType == "Cross", ]
    
    if (nrow(cross_data) == 0) {
      stop("No cross data found")
    }
    
    cat("Found", nrow(cross_data), "crosses\n")
    
    # Prepare cross summary
    cross_summary <- data.frame(
      Cross = cross_data$CrossID,
      Parent1 = cross_data$Parent1,
      Parent2 = cross_data$Parent2,
      stringsAsFactors = FALSE
    )
    
    # Add flanking markers if genotype data is available
    if (!is.null(processed_geno) && !is.null(processed_geno$snp_matrix) && 
        !is.null(processed_geno$snp_info)) {
      cat("Adding flanking SNP markers...\n")
      cross_summary <- add_flanking_snp_markers_safe(
        cross_summary, 
        processed_geno$snp_matrix, 
        processed_geno$snp_info, 
        window_size = 1 #one snp on either side
      )
    }
    
    # Get all trait names
    trait_names <- names(multi_trait_results$trait_gwas_results)
    if (length(trait_names) == 0) {
      stop("No traits found in multi_trait_results")
    }
    
    cat("Processing", length(trait_names), "traits:", paste(trait_names, collapse = ", "), "\n")
    
    # Get common SNPs across all traits
    common_snps <- get_common_snps_safe(multi_trait_results$trait_gwas_results)
    
    if (length(common_snps) == 0) {
      cat("Warning: No common SNPs found. Will process each trait separately.\n")
      # Use first trait's SNPs as fallback
      if (length(trait_names) > 0) {
        first_trait <- trait_names[1]
        gwas_df <- multi_trait_results$trait_gwas_results[[first_trait]]
        if (!is.null(gwas_df) && nrow(gwas_df) > 0) {
          gwas_df <- standardize_gwas_columns(gwas_df)
          if ("SNP_ID" %in% colnames(gwas_df)) {
            common_snps <- unique(as.character(gwas_df$SNP_ID))
            cat("Using SNPs from first trait (", first_trait, "): ", length(common_snps), " SNPs\n")
          }
        }
      }
    } else {
      cat("Found", length(common_snps), "SNPs common across all traits\n")
    }
    
    if (length(common_snps) == 0) {
      stop("No SNPs available for analysis")
    }
    
    # Initialize list to store results for each trait
    all_results <- list()
    
    # Process each trait
    for (trait_index in seq_along(trait_names)) {
      trait <- trait_names[trait_index]
      cat("\n  Processing trait", trait_index, "/", length(trait_names), ":", trait, "\n")
      
      # Get GWAS results for this trait
      gwas_df <- multi_trait_results$trait_gwas_results[[trait]]
      
      if (is.null(gwas_df) || nrow(gwas_df) == 0) {
        cat("    Warning: No GWAS results for trait", trait, "- skipping\n")
        next
      }
      
      # Standardize column names
      gwas_df <- standardize_gwas_columns(gwas_df)
      
      if (is.null(gwas_df) || nrow(gwas_df) == 0) {
        cat("    Warning: Standardization failed for trait", trait, "- skipping\n")
        next
      }
      
      # Check if we have SNP_ID column
      if (!"SNP_ID" %in% colnames(gwas_df)) {
        cat("    Warning: No SNP_ID column found for trait", trait, "- skipping\n")
        next
      }
      
      # Filter to common SNPs only
      if (length(common_snps) > 0) {
        gwas_df <- gwas_df[gwas_df$SNP_ID %in% common_snps, ]
      }
      
      if (nrow(gwas_df) == 0) {
        cat("    Warning: No common SNPs for trait", trait, "- using all available SNPs\n")
        # Use original SNPs if no common ones
        gwas_df <- standardize_gwas_columns(multi_trait_results$trait_gwas_results[[trait]])
      }
      
      # Limit number of SNPs if needed
      if (nrow(gwas_df) > max_snps_per_trait) {
        cat("    Limiting to top", max_snps_per_trait, "SNPs (by P-value)\n")
        # Ensure P_value is numeric and sortable
        gwas_df$P_value <- as.numeric(as.character(gwas_df$P_value))
        gwas_df <- gwas_df[order(gwas_df$P_value, na.last = TRUE), ]
        gwas_df <- gwas_df[1:min(max_snps_per_trait, nrow(gwas_df)), ]
      }
      
      cat("    Using", nrow(gwas_df), "SNPs for this trait\n")
      
      # Create cross-SNP combinations
      cat("    Creating cross-SNP combinations...\n")
      
      # Get SNP list
      snp_list <- unique(as.character(gwas_df$SNP_ID))
      snp_list <- snp_list[!is.na(snp_list) & snp_list != ""]
      
      if (length(snp_list) == 0) {
        cat("    Warning: No valid SNP IDs for trait", trait, "- skipping\n")
        next
      }
      
      # Create combinations using expand.grid (more reliable than expand_grid)
      cross_snp_combinations <- expand.grid(
        Cross = cross_summary$Cross,
        SNP_ID = snp_list,
        stringsAsFactors = FALSE
      )
      
      cat("    Created", nrow(cross_snp_combinations), "combinations\n")
      
      # Merge with cross summary
      trait_df <- merge(cross_snp_combinations, cross_summary, by = "Cross", all.x = TRUE)
      
      # Check if we have data to merge
      if (nrow(trait_df) == 0) {
        cat("    Warning: No data after merging with cross summary - skipping\n")
        next
      }
      
      # Merge with GWAS data - ensure SNP_ID is character in both
      trait_df$SNP_ID <- as.character(trait_df$SNP_ID)
      gwas_df$SNP_ID <- as.character(gwas_df$SNP_ID)
      
      # Select only necessary columns from gwas_df to avoid duplicates
      gwas_cols <- c("SNP_ID", "Chromosome", "Position", "P_value", "Effect", "PVE")
      gwas_cols <- gwas_cols[gwas_cols %in% colnames(gwas_df)]
      
      if (length(gwas_cols) == 0) {
        cat("    Warning: No GWAS columns to merge - skipping\n")
        next
      }
      
      # Merge
      trait_df <- merge(trait_df, gwas_df[, gwas_cols, drop = FALSE], 
                        by = "SNP_ID", all.x = TRUE)
      
      if (nrow(trait_df) == 0) {
        cat("    Warning: Merge resulted in empty dataframe - skipping\n")
        next
      }
      
      # Add trait information
      trait_df$Trait <- trait
      
      # Calculate LOD if P_value exists
      if ("P_value" %in% colnames(trait_df)) {
        trait_df$LOD <- -log10(as.numeric(trait_df$P_value))
      }
      
      # Add significance flag
      if ("P_value" %in% colnames(trait_df)) {
        trait_df$Significant <- as.numeric(trait_df$P_value) < p_threshold
      } else {
        trait_df$Significant <- NA
      }
      
      # Add parent genotypes if requested and available
      if (include_parent_genotypes && !is.null(processed_geno) && 
          !is.null(processed_geno$snp_matrix)) {
        trait_df <- add_parent_genotypes_robust(trait_df, processed_geno$snp_matrix, cross_summary)
      } else {
        # Initialize empty genotype columns
        trait_df$Parent1_Genotype <- NA
        trait_df$Parent2_Genotype <- NA
        trait_df$Expected_Genotype <- NA
      }
      
      # Store results
      all_results[[trait]] <- trait_df
      
      cat("    Successfully created", nrow(trait_df), "rows for trait", trait, "\n")
    }
    
    # Combine all trait results
    if (length(all_results) == 0) {
      cat("\nWarning: No data created for any trait\n")
      return(NULL)
    }
    
    cat("\nCombining results from all traits...\n")
    combined_df <- do.call(rbind, all_results)
    rownames(combined_df) <- NULL
    
    cat("Combined dataframe has", nrow(combined_df), "rows\n")
    
    # Reorder columns
    final_cols <- c("Cross", "Parent1", "Parent2", "SNP_ID", "Chromosome", "Position",
                    "Trait", "P_value", "LOD", "PVE", "Effect", "Significant",
                    "Parent1_Genotype", "Parent2_Genotype", "Expected_Genotype")
    
    # Add flanking marker columns if available
    if ("LeftFlankingSNPs" %in% colnames(combined_df)) {
      final_cols <- c(final_cols, "LeftFlankingSNPs")
    }
    if ("RightFlankingSNPs" %in% colnames(combined_df)) {
      final_cols <- c(final_cols, "RightFlankingSNPs")
    }
    
    # Select only columns that exist
    existing_cols <- final_cols[final_cols %in% colnames(combined_df)]
    combined_df <- combined_df[, existing_cols, drop = FALSE]
    
    # Sort the dataframe
    combined_df <- sort_cross_snp_dataframe(combined_df)
    
    cat("\n=== Creation Complete ===\n")
    cat("Total rows:", nrow(combined_df), "\n")
    cat("Unique crosses:", length(unique(combined_df$Cross)), "\n")
    cat("Unique SNPs:", length(unique(combined_df$SNP_ID)), "\n")
    cat("Traits included:", paste(unique(combined_df$Trait), collapse = ", "), "\n")
    
    if ("Significant" %in% colnames(combined_df)) {
      cat("Significant rows (p <", p_threshold, "):", 
          sum(combined_df$Significant, na.rm = TRUE), "\n")
    }
    
    # Calculate average PVE across traits
    if ("PVE" %in% colnames(combined_df) && nrow(combined_df) > 0) {
      avg_pve <- mean(as.numeric(combined_df$PVE), na.rm = TRUE)
      cat("Average PVE:", round(avg_pve, 2), "%\n")
    }
    
    return(combined_df)
    
  }, error = function(e) {
    cat("\nError creating cross-SNP dataframe for all traits:\n")
    cat("Message:", e$message, "\n")
    cat("Traceback:\n")
    print(traceback())
    return(NULL)
  })
}


#==============================================================================
# ADDITIONAL CORRECTED HELPER FUNCTIONS
#==============================================================================

#' Safe version of add_flanking_snp_markers
add_flanking_snp_markers_safe <- function(cross_summary, geno_matrix, snp_info, window_size = 2) {
  tryCatch({
    cat("Adding flanking SNP markers for", nrow(cross_summary), "crosses...\n")
    
    # Initialize new columns
    cross_summary$LeftFlankingSNPs <- ""
    cross_summary$RightFlankingSNPs <- ""
    
    # Clean row names for matching
    geno_samples <- rownames(geno_matrix)
    geno_samples_clean <- tolower(trimws(geno_samples))
    
    # Process each cross
    for (i in 1:nrow(cross_summary)) {
      parent1 <- as.character(cross_summary$Parent1[i])
      parent2 <- as.character(cross_summary$Parent2[i])
      
      # Clean parent names for matching
      parent1_clean <- tolower(trimws(parent1))
      parent2_clean <- tolower(trimws(parent2))
      
      # Check if parents are in genotype matrix
      parent1_in_matrix <- parent1_clean %in% geno_samples_clean
      parent2_in_matrix <- parent2_clean %in% geno_samples_clean
      
      if (!parent1_in_matrix || !parent2_in_matrix) {
        # Try to find partial matches
        if (!parent1_in_matrix) {
          matches <- grep(gsub("[^a-zA-Z0-9]", "", parent1_clean), 
                          geno_samples_clean, value = TRUE)
          if (length(matches) > 0) {
            parent1_clean <- matches[1]
            parent1_in_matrix <- TRUE
          }
        }
        
        if (!parent2_in_matrix) {
          matches <- grep(gsub("[^a-zA-Z0-9]", "", parent2_clean), 
                          geno_samples_clean, value = TRUE)
          if (length(matches) > 0) {
            parent2_clean <- matches[1]
            parent2_in_matrix <- TRUE
          }
        }
      }
      
      if (!parent1_in_matrix || !parent2_in_matrix) {
        # cat("  Cross", i, ": Parents not in genotype matrix\n")
        next
      }
      
      # Get original names from cleaned names
      parent1_orig <- geno_samples[which(geno_samples_clean == parent1_clean)[1]]
      parent2_orig <- geno_samples[which(geno_samples_clean == parent2_clean)[1]]
      
      # Get genotype data for both parents
      geno_parent1 <- geno_matrix[parent1_orig, ]
      geno_parent2 <- geno_matrix[parent2_orig, ]
      
      # Find polymorphic SNPs (where parents have different genotypes)
      polymorphic_snps <- which(geno_parent1 != geno_parent2 & 
                                  !is.na(geno_parent1) & 
                                  !is.na(geno_parent2))
      
      if (length(polymorphic_snps) == 0) {
        # cat("  Cross", i, ": No polymorphic SNPs found\n")
        next
      }
      
      # Get SNP IDs for polymorphic SNPs
      poly_snp_ids <- colnames(geno_matrix)[polymorphic_snps]
      
      # Get chromosome and position info for polymorphic SNPs
      poly_snp_info <- snp_info[snp_info$SNP_ID %in% poly_snp_ids, ]
      
      if (nrow(poly_snp_info) == 0) {
        next
      }
      
      # Order by chromosome and position
      poly_snp_info <- poly_snp_info[order(poly_snp_info$Chromosome, poly_snp_info$Position), ]
      
      # Get left flanking SNPs (first 'window_size' SNPs)
      if (nrow(poly_snp_info) >= window_size) {
        left_snps <- head(poly_snp_info$SNP_ID, window_size)
        cross_summary$LeftFlankingSNPs[i] <- paste(left_snps, collapse = ";")
      } else {
        cross_summary$LeftFlankingSNPs[i] <- paste(poly_snp_info$SNP_ID, collapse = ";")
      }
      
      # Get right flanking SNPs (last 'window_size' SNPs)
      if (nrow(poly_snp_info) >= window_size) {
        right_snps <- tail(poly_snp_info$SNP_ID, window_size)
        cross_summary$RightFlankingSNPs[i] <- paste(right_snps, collapse = ";")
      } else {
        cross_summary$RightFlankingSNPs[i] <- paste(poly_snp_info$SNP_ID, collapse = ";")
      }
      
      # Progress update
      if (i %% 10 == 0) {
        cat("  Processed", i, "/", nrow(cross_summary), "crosses\n")
      }
    }
    
    cat("Completed adding flanking SNP markers\n")
    cat("  Crosses with flanking SNPs:", sum(cross_summary$LeftFlankingSNPs != ""), "\n")
    
    return(cross_summary)
    
  }, error = function(e) {
    cat("Error in add_flanking_snp_markers_safe:", e$message, "\n")
    return(cross_summary)
  })
}



#' Robust parent genotype extraction
add_parent_genotypes_robust <- function(df, geno_matrix, cross_summary) {
  tryCatch({
    cat("    Adding parent genotypes...\n")
    
    # Initialize genotype columns
    df$Parent1_Genotype <- NA
    df$Parent2_Genotype <- NA
    df$Expected_Genotype <- NA
    
    # Clean genotype matrix row names
    geno_samples <- rownames(geno_matrix)
    geno_samples_clean <- tolower(trimws(geno_samples))
    
    # Create a mapping from cleaned to original names
    geno_map <- setNames(geno_samples, geno_samples_clean)
    
    # Get unique parents
    unique_parents <- unique(c(cross_summary$Parent1, cross_summary$Parent2))
    
    # Create a mapping of parent names to cleaned versions
    parent_map <- list()
    for (parent in unique_parents) {
      parent_clean <- tolower(trimws(as.character(parent)))
      parent_map[[parent]] <- parent_clean
    }
    
    # Process each unique row in df
    for (i in 1:nrow(df)) {
      parent1 <- as.character(df$Parent1[i])
      parent2 <- as.character(df$Parent2[i])
      snp_id <- as.character(df$SNP_ID[i])
      
      # Skip if SNP not in genotype matrix
      if (!snp_id %in% colnames(geno_matrix)) {
        next
      }
      
      # Get cleaned parent names
      parent1_clean <- if (!is.null(parent_map[[parent1]])) parent_map[[parent1]] else tolower(trimws(parent1))
      parent2_clean <- if (!is.null(parent_map[[parent2]])) parent_map[[parent2]] else tolower(trimws(parent2))
      
      # Find parent1 in genotype matrix
      if (parent1_clean %in% names(geno_map)) {
        parent1_orig <- geno_map[[parent1_clean]]
        df$Parent1_Genotype[i] <- geno_matrix[parent1_orig, snp_id]
      }
      
      # Find parent2 in genotype matrix
      if (parent2_clean %in% names(geno_map)) {
        parent2_orig <- geno_map[[parent2_clean]]
        df$Parent2_Genotype[i] <- geno_matrix[parent2_orig, snp_id]
      }
      
      # Calculate expected genotype if both parents available
      if (!is.na(df$Parent1_Genotype[i]) && !is.na(df$Parent2_Genotype[i])) {
        df$Expected_Genotype[i] <- mean(c(df$Parent1_Genotype[i], df$Parent2_Genotype[i]), na.rm = TRUE)
      } else if (!is.na(df$Parent1_Genotype[i])) {
        df$Expected_Genotype[i] <- df$Parent1_Genotype[i]
      } else if (!is.na(df$Parent2_Genotype[i])) {
        df$Expected_Genotype[i] <- df$Parent2_Genotype[i]
      }
    }
    
    cat("    Genotype extraction complete\n")
    cat("    Parent1 genotypes available:", sum(!is.na(df$Parent1_Genotype)), "/", nrow(df), "\n")
    cat("    Parent2 genotypes available:", sum(!is.na(df$Parent2_Genotype)), "/", nrow(df), "\n")
    
    return(df)
    
  }, error = function(e) {
    cat("    Error in add_parent_genotypes_robust:", e$message, "\n")
    return(df)
  })
}

#' Sort cross-SNP dataframe
sort_cross_snp_dataframe <- function(df) {
  if ("Chromosome" %in% colnames(df) && "Position" %in% colnames(df)) {
    # Convert chromosomes to factor with natural ordering
    chrom_levels <- unique(df$Chromosome)
    
    # Try to order numerically first, then alphabetically
    chrom_numeric <- suppressWarnings(as.numeric(chrom_levels))
    numeric_chroms <- chrom_levels[!is.na(chrom_numeric)]
    non_numeric_chroms <- chrom_levels[is.na(chrom_numeric)]
    
    # Order numeric chromosomes
    if (length(numeric_chroms) > 0) {
      numeric_chroms <- numeric_chroms[order(as.numeric(numeric_chroms))]
    }
    
    # Combine
    ordered_levels <- c(numeric_chroms, sort(non_numeric_chroms))
    
    # Convert to factor with ordered levels
    df$Chromosome <- factor(df$Chromosome, levels = ordered_levels)
    
    # Order by trait, cross, chromosome, position
    df <- df[order(df$Trait, df$Cross, df$Chromosome, as.numeric(df$Position)), ]
    
    # Convert back to character
    df$Chromosome <- as.character(df$Chromosome)
  } else {
    df <- df[order(df$Trait, df$Cross), ]
  }
  
  return(df)
}


#==============================================================================
# SIMPLIFIED VERSION FOR SERVER
#==============================================================================

#' Simplified cross-SNP analysis for server
run_cross_snp_analysis_all_traits_app <- function(values, p_threshold = 0.05) {
  cat("\n=== Running Cross-SNP Analysis for ALL Traits ===\n")
  
  # Check prerequisites
  if (is.null(values$diallel_data)) {
    showNotification("Diallel data is required for cross-SNP analysis", type = "error")
    return(NULL)
  }
  
  if (is.null(values$multi_trait_results)) {
    showNotification("Multi-trait GWAS results are required", type = "error")
    return(NULL)
  }
  
  # Get trait count
  trait_names <- names(values$multi_trait_results$trait_gwas_results)
  if (length(trait_names) == 0) {
    showNotification("No GWAS results found", type = "error")
    return(NULL)
  }
  
  cat("Found", length(trait_names), "traits to process\n")
  
  # Create cross-SNP dataframe for all traits using FIXED function
  cross_snp_df <- create_cross_snp_dataframe_all_traits_fixed(
    multi_trait_results = values$multi_trait_results,
    diallel_data = values$diallel_data,
    processed_geno = values$processed_geno,
    p_threshold = p_threshold,
    max_snps_per_trait = 500,  # Reduced to avoid memory issues
    include_parent_genotypes = FALSE  # Set to TRUE if you want genotypes
  )
  
  if (!is.null(cross_snp_df)) {
    cat("Successfully created cross-SNP dataframe for all traits\n")
    cat("  Dimensions:", dim(cross_snp_df), "\n")
    cat("  Memory usage:", format(object.size(cross_snp_df), units = "MB"), "\n")
    
    # Store results in values
    values$cross_snp_df <- cross_snp_df
    values$dl_cross_snp_detailed <- cross_snp_df
    
    # Create summaries if dataframe is not too large
    if (nrow(cross_snp_df) > 0 && nrow(cross_snp_df) < 100000) {
      tryCatch({
        # 1. Cross summary across all traits
        cross_summary <- cross_snp_df %>%
          group_by(Cross, Parent1, Parent2) %>%
          summarise(
            N_SNPs = n_distinct(SNP_ID),
            N_Traits = n_distinct(Trait),
            N_Significant = sum(Significant, na.rm = TRUE),
            Mean_P_value = mean(as.numeric(P_value), na.rm = TRUE),
            Min_P_value = min(as.numeric(P_value), na.rm = TRUE),
            .groups = "drop"
          )
        
        values$cross_summary <- cross_summary
        
        # 2. SNP summary across all crosses and traits
        top_crosses_snp_combined <- cross_snp_df %>%
          group_by(SNP_ID) %>%
          summarise(
            N_Crosses = n_distinct(Cross),
            N_Traits = n_distinct(Trait),
            N_Significant = sum(Significant, na.rm = TRUE),
            Mean_P_value = mean(as.numeric(P_value), na.rm = TRUE),
            .groups = "drop"
          )
        
        values$top_crosses_snp_combined <- top_crosses_snp_combined
        
        # 3. Trait summary
        trait_summary <- cross_snp_df %>%
          group_by(Trait) %>%
          summarise(
            N_SNPs = n_distinct(SNP_ID),
            N_Crosses = n_distinct(Cross),
            N_Significant = sum(Significant, na.rm = TRUE),
            Mean_P_value = mean(as.numeric(P_value), na.rm = TRUE),
            .groups = "drop"
          )
        
        values$trait_summary <- trait_summary
        
      }, error = function(e) {
        cat("Error creating summaries:", e$message, "\n")
      })
    }
    
    showNotification(
      paste("Cross-SNP analysis completed for", length(trait_names), "traits"), 
      type = "message", 
      duration = 5
    )
    
    return(cross_snp_df)
  } else {
    showNotification(
      "Cross-SNP analysis failed. Check console for details.",
      type = "error",
      duration = 5
    )
  }
  
  return(NULL)
}


# ============================================================
# POPULATION GENETICS SERVER FUNCTIONS
# ============================================================

#' Calculate allele frequencies from genotype data
calculate_allele_frequencies <- function(gl_object, pop_assignments = NULL) {
  if (is.null(gl_object)) return(NULL)
  
  tryCatch({
    # Convert to matrix
    geno_matrix <- as.matrix(gl_object)
    
    # If population assignments provided, calculate per population
    if (!is.null(pop_assignments) && length(pop_assignments) == nrow(geno_matrix)) {
      populations <- unique(pop_assignments)
      
      # Calculate allele frequencies per population
      allele_freqs <- list()
      for (pop in populations) {
        pop_indices <- which(pop_assignments == pop)
        if (length(pop_indices) > 0) {
          # Calculate frequency of alternate allele (assuming 0,1,2 coding)
          pop_matrix <- geno_matrix[pop_indices, , drop = FALSE]
          # Handle missing values
          pop_matrix[is.na(pop_matrix)] <- 0
          freq <- colMeans(pop_matrix, na.rm = TRUE) / 2
          allele_freqs[[pop]] <- freq
        }
      }
      
      # Create data frame
      freq_df <- as.data.frame(do.call(cbind, allele_freqs))
      colnames(freq_df) <- paste0("allelefreq_", names(allele_freqs))
      
      # Add SNP information
      freq_df$SNP_ID <- colnames(geno_matrix)
      
    } else {
      # Calculate overall allele frequencies
      geno_matrix[is.na(geno_matrix)] <- 0
      overall_freq <- colMeans(geno_matrix, na.rm = TRUE) / 2
      freq_df <- data.frame(
        allelefreq_overall = overall_freq,
        SNP_ID = colnames(geno_matrix)
      )
    }
    
    return(freq_df)
    
  }, error = function(e) {
    cat("Error calculating allele frequencies:", e$message, "\n")
    return(NULL)
  })
}

#' Modified corrfrequencies function for Shiny
corrfrequencies <- function(freq_df, populations) {
  tryCatch({
    # Extract frequency columns
    freq_cols <- grep("^allelefreq_", colnames(freq_df), value = TRUE)
    
    if (length(freq_cols) < 2) {
      stop("Need at least two populations for correlation analysis")
    }
    
    # Create combinations
    combitable <- combn(freq_cols, m = 2)
    mycorr1 <- rep(NA, ncol(combitable))  # pearson
    mycorr2 <- rep(NA, ncol(combitable))  # spearman
    
    for (i in 1:ncol(combitable)) {
      pop1 <- combitable[1, i]
      freq1 <- freq_df[[pop1]]
      pop2 <- combitable[2, i]
      freq2 <- freq_df[[pop2]]
      
      # Remove NA values
      complete_cases <- complete.cases(freq1, freq2)
      if (sum(complete_cases) > 10) {
        mycorr1[i] <- cor(freq1[complete_cases], freq2[complete_cases], method = "pearson")
        mycorr2[i] <- cor(freq1[complete_cases], freq2[complete_cases], method = "spearman")
      }
    }
    
    # Create correlation matrix
    n_pops <- length(freq_cols)
    b <- matrix(1, nrow = n_pops, ncol = n_pops)
    
    # Fill lower triangle with Pearson, upper with Spearman
    b[lower.tri(b, diag = FALSE)] <- mycorr1
    b <- t(b)
    b[lower.tri(b, diag = FALSE)] <- mycorr2
    
    # Set row and column names
    pop_names <- gsub("^allelefreq_", "", freq_cols)
    colnames(b) <- rownames(b) <- pop_names
    
    # Create table for display
    corr_table <- data.frame(
      Population1 = gsub("^allelefreq_", "", combitable[1, ]),
      Population2 = gsub("^allelefreq_", "", combitable[2, ]),
      Pearson = round(mycorr1, 4),
      Spearman = round(mycorr2, 4)
    )
    
    return(list(
      matrix = b,
      table = corr_table,
      pearson = mycorr1,
      spearman = mycorr2
    ))
    
  }, error = function(e) {
    cat("Error in corrfrequencies:", e$message, "\n")
    return(NULL)
  })
}

#' Modified allelefreqhisto function for Shiny
allelefreqhisto <- function(freq_df, populations, export = FALSE) {
  tryCatch({
    # Extract frequency columns for specified populations
    freq_cols <- paste0("allelefreq_", populations)
    freq_cols <- freq_cols[freq_cols %in% colnames(freq_df)]
    
    if (length(freq_cols) == 0) {
      stop("No frequency data found for specified populations")
    }
    
    # Prepare data for plotting
    plot_data <- freq_df[, freq_cols, drop = FALSE]
    colnames(plot_data) <- gsub("^allelefreq_", "", colnames(plot_data))
    
    # Melt data for ggplot
    plot_data_melt <- reshape2::melt(plot_data)
    colnames(plot_data_melt) <- c("Population", "Frequency")
    
    # Create histogram
    p <- ggplot(plot_data_melt, aes(x = Frequency, fill = Population)) +
      geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
      facet_wrap(~ Population, scales = "free_y") +
      labs(title = "Allele Frequency Distribution by Population",
           x = "Allele Frequency",
           y = "Count") +
      theme_minimal() +
      theme(legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"),
            axis.text = element_text(size = 9))
    
    return(p)
    
  }, error = function(e) {
    cat("Error in allelefreqhisto:", e$message, "\n")
    return(NULL)
  })
}

#' Modified waples_snps function for Shiny
waples_snps <- function(data, fis_column, fst_column, mysteps = 0.05) {
  tryCatch({
    # Extract columns
    if (is.null(fis_column) || is.null(fst_column)) {
      stop("Fis and Fst columns must be specified")
    }
    
    fis_data <- data[[fis_column]]
    fst_data <- data[[fst_column]]
    
    # Remove NA values
    complete_cases <- complete.cases(fis_data, fst_data)
    fis_data <- fis_data[complete_cases]
    fst_data <- fst_data[complete_cases]
    
    if (length(fis_data) == 0) {
      stop("No complete Fst-Fis data available")
    }
    
    # Create bins
    mybreaks <- seq(0, 1, mysteps)
    myhalf <- mysteps / 2
    mylabels <- seq(myhalf, 1 - myhalf, mysteps)
    
    # Bin Fst values
    fstbins <- cut(fst_data, mybreaks, include.lowest = TRUE)
    
    # Calculate mean Fis per bin
    fisperbin <- aggregate(fis_data, by = list(fstbins), FUN = mean, na.rm = TRUE)
    nperbin <- aggregate(fis_data, by = list(fstbins), FUN = length)
    
    nbins <- nrow(fisperbin)
    ndata <- nperbin$x
    
    # Create binned data for plotting
    binned_data <- data.frame(
      Fst_midpoint = mylabels[1:nbins],
      Mean_Fis = fisperbin$x,
      N_SNPs = ndata
    )
    
    # Create plot
    p <- ggplot(binned_data[binned_data$N_SNPs > 10, ], 
                aes(x = Fst_midpoint, y = Mean_Fis)) +
      geom_point(aes(size = N_SNPs), color = "steelblue") +
      geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
      labs(title = "Waples Plot: Fst vs Fis",
           x = "Fst (midpoint of bin)",
           y = "Mean Fis",
           size = "Number of SNPs") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(list(plot = p, binned_data = binned_data))
    
  }, error = function(e) {
    cat("Error in waples_snps:", e$message, "\n")
    return(NULL)
  })
}



# ============================================================================
# COMPREHENSIVE POPULATION GENETICS ANALYSIS FUNCTIONS
# ============================================================================

#' Convert genlight to genind with population structure - FIXED VERSION
#' @param gl_object Genlight object
#' @param pop_data Population assignment data (optional)
#' @return Genind object with population structure
create_genind_from_gl <- function(gl_object, pop_data = NULL) {
  tryCatch({
    cat("Converting genlight to genind...\n")
    
    # Check if gl_object is valid
    if (is.null(gl_object)) {
      stop("Genlight object is NULL")
    }
    
    # Get basic information
    n_ind <- nInd(gl_object)
    n_loc <- nLoc(gl_object)
    
    cat("  Individuals:", n_ind, " Loci:", n_loc, "\n")
    
    if (n_loc == 0 || n_ind == 0) {
      stop("Genlight object has no data")
    }
    
    # Convert to matrix
    geno_matrix <- as.matrix(gl_object)
    
    # Check dimensions
    if (is.null(dim(geno_matrix)) || length(dim(geno_matrix)) != 2) {
      cat("Warning: Genotype matrix has unexpected dimensions. Attempting to fix...\n")
      if (is.vector(geno_matrix)) {
        geno_matrix <- matrix(geno_matrix, nrow = n_ind, ncol = n_loc)
      }
    }
    
    # Set row and column names if missing
    if (is.null(rownames(geno_matrix))) {
      rownames(geno_matrix) <- indNames(gl_object)
    }
    if (is.null(colnames(geno_matrix))) {
      colnames(geno_matrix) <- locNames(gl_object)
    }
    
    # Get population information
    if (!is.null(pop_data)) {
      if (is.data.frame(pop_data)) {
        if ("Population" %in% colnames(pop_data)) {
          pop_info <- as.character(pop_data$Population)
        } else if (ncol(pop_data) >= 1) {
          pop_info <- as.character(pop_data[[1]])
        } else {
          pop_info <- rep("Unknown", n_ind)
        }
      } else if (is.vector(pop_data)) {
        pop_info <- as.character(pop_data)
      } else {
        pop_info <- rep("Unknown", n_ind)
      }
    } else if (!is.null(pop(gl_object))) {
      pop_info <- as.character(pop(gl_object))
    } else {
      pop_info <- rep("Unknown", n_ind)
    }
    
    # Ensure pop_info has correct length
    if (length(pop_info) != n_ind) {
      cat("Warning: Population info length doesn't match number of individuals. Using 'Unknown' for all.\n")
      pop_info <- rep("Unknown", n_ind)
    }
    
    # Convert to genind using dartR::gl2gi (more reliable)
    cat("  Using dartR::gl2gi for conversion...\n")
    genind_obj <- dartR::gl2gi(gl_object, verbose = 0)
    
    # Set population information
    if (!is.null(genind_obj)) {
      pop(genind_obj) <- pop_info
    }
    
    cat("Genind object created successfully\n")
    cat("  Individuals:", nInd(genind_obj), "\n")
    cat("  Loci:", nLoc(genind_obj), "\n")
    cat("  Populations:", length(unique(pop(genind_obj))), "\n")
    
    return(genind_obj)
    
  }, error = function(e) {
    cat("Error creating genind object:", e$message, "\n")
    
    # Try simple alternative
    tryCatch({
      cat("Trying alternative conversion method...\n")
      
      # Create simple genind with basic data
      geno_matrix <- as.matrix(gl_object)
      
      # Convert to character format for df2genind
      geno_char <- matrix(as.character(geno_matrix), nrow = nrow(geno_matrix))
      
      # Replace NA
      geno_char[is.na(geno_char)] <- "NA/NA"
      
      # Create basic genind
      genind_obj <- adegenet::df2genind(geno_char, sep = "/", 
                                        ind.names = rownames(geno_matrix),
                                        loc.names = colnames(geno_matrix),
                                        type = "codom")
      
      return(genind_obj)
      
    }, error = function(e2) {
      cat("Alternative conversion also failed:", e2$message, "\n")
      return(NULL)
    })
  })
}



#' Create allele frequency barcode (heatmap) for populations
#' @param genind_obj Genind object
#' @param populations Vector of population names
#' @param max_markers Maximum markers to include (for performance)
#' @return List containing allele frequency matrix and metadata
create_allele_barcode <- function(genind_obj, populations = NULL, max_markers = 500) {
  tryCatch({
    cat("Creating allele frequency barcode...\n")
    
    # Get population names
    if (is.null(populations)) {
      populations <- unique(as.character(pop(genind_obj)))
    }
    
    n_pops <- length(populations)
    cat("  Analyzing", n_pops, "populations\n")
    
    # Get allele frequency matrix
    # Convert to frequency table
    allele_freq <- tab(genind_obj, freq = TRUE, NA.method = "mean")
    
    # Get population indices
    pop_info <- as.character(pop(genind_obj))
    
    # Calculate mean frequency per population
    freq_matrix <- matrix(NA, nrow = ncol(allele_freq), ncol = n_pops)
    colnames(freq_matrix) <- populations
    rownames(freq_matrix) <- colnames(allele_freq)
    
    for (i in 1:n_pops) {
      pop_name <- populations[i]
      pop_indices <- which(pop_info == pop_name)
      
      if (length(pop_indices) > 0) {
        pop_freq <- allele_freq[pop_indices, , drop = FALSE]
        freq_matrix[, i] <- colMeans(pop_freq, na.rm = TRUE)
      } else {
        freq_matrix[, i] <- NA
      }
    }
    
    # Remove markers with all NAs
    valid_markers <- rowSums(!is.na(freq_matrix)) > 0
    freq_matrix <- freq_matrix[valid_markers, , drop = FALSE]
    
    # Limit to max_markers if specified
    if (nrow(freq_matrix) > max_markers) {
      # Select markers with highest variance across populations
      marker_variance <- apply(freq_matrix, 1, var, na.rm = TRUE)
      top_markers <- order(marker_variance, decreasing = TRUE)[1:max_markers]
      freq_matrix <- freq_matrix[top_markers, , drop = FALSE]
      cat("  Limited to", max_markers, "most variable markers\n")
    }
    
    # Create metadata for loci
    locus_info <- data.frame(
      Marker = rownames(freq_matrix),
      Chromosome = NA,
      Position = NA,
      stringsAsFactors = FALSE
    )
    
    # Try to extract chromosome and position if available
    if (!is.null(genind_obj@other)) {
      # Check for chromosome information
      if (!is.null(genind_obj@other$chromosome)) {
        locus_info$Chromosome <- genind_obj@other$chromosome[valid_markers]
      }
      if (!is.null(genind_obj@other$position)) {
        locus_info$Position <- genind_obj@other$position[valid_markers]
      }
    }
    
    results <- list(
      frequency_matrix = freq_matrix,
      locus_info = locus_info,
      populations = populations,
      n_markers = nrow(freq_matrix),
      n_populations = n_pops
    )
    
    cat("  Created frequency matrix:", dim(freq_matrix), "\n")
    
    return(results)
    
  }, error = function(e) {
    cat("Error creating allele barcode:", e$message, "\n")
    return(NULL)
  })
}

#' Plot allele frequency heatmap
#' @param barcode_results Results from create_allele_barcode
#' @param color_palette Color palette for heatmap
#' @param show_dendrogram Whether to show dendrogram
#' @param cluster_rows Whether to cluster rows (markers)
#' @param cluster_cols Whether to cluster columns (populations)
#' @return ggplot object or base R plot
plot_allele_heatmap <- function(barcode_results, color_palette = "viridis", 
                                show_dendrogram = TRUE, cluster_rows = TRUE, 
                                cluster_cols = TRUE) {
  tryCatch({
    freq_matrix <- barcode_results$frequency_matrix
    
    # Prepare data for plotting
    plot_data <- as.data.frame(freq_matrix)
    plot_data$Marker <- rownames(freq_matrix)
    
    # Melt for ggplot
    melt_data <- reshape2::melt(plot_data, id.vars = "Marker", 
                                variable.name = "Population", 
                                value.name = "Frequency")
    
    # Order populations by mean frequency
    pop_order <- names(sort(colMeans(freq_matrix, na.rm = TRUE)))
    melt_data$Population <- factor(melt_data$Population, levels = pop_order)
    
    # Create heatmap
    heatmap_plot <- ggplot(melt_data, aes(x = Population, y = Marker, fill = Frequency)) +
      geom_tile() +
      scale_fill_viridis_c(option = color_palette, na.value = "gray90") +
      labs(title = "Allele Frequency Barcode",
           subtitle = paste("Showing", nrow(freq_matrix), "markers across", 
                            ncol(freq_matrix), "populations"),
           x = "Population", y = "Marker") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    return(heatmap_plot)
    
  }, error = function(e) {
    cat("Error plotting allele heatmap:", e$message, "\n")
    
    # Try base R heatmap as fallback
    tryCatch({
      heatmap.2(barcode_results$frequency_matrix,
                col = viridis::viridis(100),
                scale = "none",
                trace = "none",
                dendrogram = if (show_dendrogram) "both" else "none",
                Rowv = cluster_rows,
                Colv = cluster_cols,
                main = "Allele Frequency Heatmap",
                margins = c(10, 8))
    }, error = function(e2) {
      cat("Base R heatmap also failed:", e2$message, "\n")
      return(NULL)
    })
  })
}

#' Calculate population heterozygosity statistics
#' @param genind_obj Genind object
#' @param populations Vector of population names (optional)
#' @return List containing heterozygosity statistics
#' Fixed: Calculate population heterozygosity statistics
calculate_heterozygosity_stats <- function(genind_obj, populations = NULL) {
  tryCatch({
    cat("Calculating heterozygosity statistics...\n")
    
    if (is.null(populations)) {
      populations <- unique(as.character(adegenet::pop(genind_obj)))
    }
    
    results <- list()
    summary_stats <- data.frame()
    
    for (pop_name in populations) {
      cat("  Processing population:", pop_name, "\n")
      
      # Subset population using popsub from poppr (not adegenet)
      pop_genind <- poppr::popsub(genind_obj, sublist = pop_name)
      
      if (adegenet::nInd(pop_genind) < 2) {
        cat("    Skipping - insufficient individuals\n")
        next
      }
      
      # Calculate summary statistics
      pop_summary <- summary(pop_genind)
      
      # Extract heterozygosity
      Hobs <- mean(pop_summary$Hobs, na.rm = TRUE)
      Hexp <- mean(pop_summary$Hexp, na.rm = TRUE)
      
      # Calculate fixation index (Fis)
      if (Hexp > 0) {
        Fis <- 1 - (Hobs / Hexp)
      } else {
        Fis <- NA
      }
      
      # Calculate per locus statistics
      per_locus <- data.frame(
        Population = pop_name,
        Locus = names(pop_summary$Hobs),
        Hobs = pop_summary$Hobs,
        Hexp = pop_summary$Hexp,
        Fis = ifelse(pop_summary$Hexp > 0, 
                     1 - (pop_summary$Hobs / pop_summary$Hexp), 
                     NA),
        stringsAsFactors = FALSE
      )
      
      # Store results
      results[[pop_name]] <- list(
        summary = pop_summary,
        per_locus = per_locus,
        overall = data.frame(
          Population = pop_name,
          N_Individuals = adegenet::nInd(pop_genind),
          Hobs_mean = Hobs,
          Hexp_mean = Hexp,
          Fis_mean = Fis,
          stringsAsFactors = FALSE
        )
      )
      
      summary_stats <- rbind(summary_stats, results[[pop_name]]$overall)
    }
    
    # Create comparison plot
    if (nrow(summary_stats) > 0) {
      comparison_plot <- ggplot2::ggplot(summary_stats, aes(x = Hexp_mean, y = Hobs_mean)) +
        ggplot2::geom_point(aes(color = Population, size = N_Individuals), alpha = 0.7) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        ggplot2::labs(title = "Observed vs Expected Heterozygosity",
                      x = "Expected Heterozygosity (He)",
                      y = "Observed Heterozygosity (Ho)",
                      size = "Sample Size") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = element_text(hjust = 0.5))
    } else {
      comparison_plot <- NULL
    }
    
    return(list(
      population_results = results,
      summary_table = summary_stats,
      comparison_plot = comparison_plot
    ))
    
  }, error = function(e) {
    cat("Error calculating heterozygosity:", e$message, "\n")
    return(NULL)
  })
}


#' Calculate pairwise FST matrix
#' @param genind_obj Genind object
#' @param method FST calculation method (default: "WC84" for Weir & Cockerham)
#' @return List containing FST matrix and NJ tree
#' Calculate pairwise FST matrix with multiple fallbacks
calculate_fst_matrix <- function(genind_obj, method = "WC84") {
  tryCatch({
    cat("Calculating pairwise FST...\n")
    
    # Check if genind_obj is valid
    if (is.null(genind_obj) || class(genind_obj) != "genind") {
      stop("Invalid genind object")
    }
    
    # Get population names
    pop_names <- as.character(unique(adegenet::pop(genind_obj)))
    
    if (length(pop_names) < 2) {
      cat("  Only 1 population found. Cannot calculate pairwise FST.\n")
      return(NULL)
    }
    
    # Try different methods for FST calculation
    fst_matrix <- NULL
    
    # Method 1: Use hierfstat if available
    if (requireNamespace("hierfstat", quietly = TRUE)) {
      cat("  Method 1: Using hierfstat for FST calculation...\n")
      tryCatch({
        # Convert to hierfstat format - use hierfstat::genind2hierfstat
        hierfstat_data <- hierfstat::genind2hierfstat(genind_obj)
        
        # Calculate FST
        fst_matrix <- hierfstat::pairwise.WCfst(hierfstat_data)
        
        # Set row and column names
        colnames(fst_matrix) <- rownames(fst_matrix) <- pop_names
        cat("  ✓ FST calculation completed using hierfstat\n")
      }, error = function(e) {
        cat("  ✗ hierfstat method failed:", e$message, "\n")
      })
    }
    
    # Method 2: Use adegenet if hierfstat fails
    if (is.null(fst_matrix)) {
      cat("  Method 2: Using adegenet for FST calculation...\n")
      tryCatch({
        # Use adegenet::pairwise.fst
        fst_matrix <- adegenet::pairwise.fst(genind_obj, res.type = "matrix")
        cat("  ✓ FST calculation completed using adegenet\n")
      }, error = function(e) {
        cat("  ✗ adegenet method failed:", e$message, "\n")
      })
    }
    
    # Method 3: Use simple method if both above fail
    if (is.null(fst_matrix)) {
      cat("  Method 3: Using simple FST calculation...\n")
      fst_matrix <- calculate_fst_simple(genind_obj)
      if (!is.null(fst_matrix)) {
        cat("  ✓ Simple FST calculation completed\n")
      }
    }
    
    if (is.null(fst_matrix) || all(is.na(fst_matrix))) {
      cat("  Warning: FST calculation returned all NA values\n")
      return(NULL)
    }
    
    # Create heatmap of FST matrix
    fst_df <- as.data.frame(as.table(fst_matrix))
    colnames(fst_df) <- c("Pop1", "Pop2", "FST")
    
    fst_plot <- ggplot2::ggplot(fst_df, ggplot2::aes(x = Pop1, y = Pop2, fill = FST)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(name = "FST", na.value = "white") +
      ggplot2::geom_text(ggplot2::aes(label = round(FST, 3)), size = 3, color = "white") +
      ggplot2::labs(title = "Pairwise FST Matrix",
                    x = "Population", y = "Population") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.grid = ggplot2::element_blank()
      )
    
    return(list(
      fst_matrix = fst_matrix,
      fst_plot = fst_plot
    ))
    
  }, error = function(e) {
    cat("Error calculating FST:", e$message, "\n")
    return(NULL)
  })
}



#' Manual FST calculation (simplified)
calculate_fst_manual <- function(genind_obj) {
  tryCatch({
    # Get population names
    pop_names <- as.character(unique(adegenet::pop(genind_obj)))
    n_pops <- length(pop_names)
    
    # Initialize FST matrix
    fst_matrix <- matrix(0, nrow = n_pops, ncol = n_pops)
    colnames(fst_matrix) <- rownames(fst_matrix) <- pop_names
    
    # Get allele frequencies per population
    freq_list <- list()
    for (pop in pop_names) {
      # Subset population
      pop_indices <- which(adegenet::pop(genind_obj) == pop)
      pop_genind <- genind_obj[pop_indices, , drop = FALSE]
      
      # Calculate allele frequencies
      freq_tab <- adegenet::tab(pop_genind, freq = TRUE, NA.method = "mean")
      freq_list[[pop]] <- colMeans(freq_tab, na.rm = TRUE)
    }
    
    # Calculate pairwise FST (simplified Weir & Cockerham)
    for (i in 1:(n_pops - 1)) {
      for (j in (i + 1):n_pops) {
        pop1 <- pop_names[i]
        pop2 <- pop_names[j]
        
        # Get allele frequencies
        p1 <- freq_list[[pop1]]
        p2 <- freq_list[[pop2]]
        
        # Remove loci with missing data
        valid_loci <- !is.na(p1) & !is.na(p2)
        p1 <- p1[valid_loci]
        p2 <- p2[valid_loci]
        
        if (length(p1) == 0) {
          fst_matrix[i, j] <- fst_matrix[j, i] <- NA
          next
        }
        
        # Calculate mean allele frequency
        p_mean <- (p1 + p2) / 2
        
        # Calculate variance between populations
        var_between <- ((p1 - p_mean)^2 + (p2 - p_mean)^2) / 2
        
        # Calculate variance within populations (assuming Hardy-Weinberg)
        var_within <- p_mean * (1 - p_mean)
        
        # Calculate FST
        fst <- mean(var_between, na.rm = TRUE) / mean(var_within, na.rm = TRUE)
        
        fst_matrix[i, j] <- fst_matrix[j, i] <- fst
      }
    }
    
    return(fst_matrix)
    
  }, error = function(e) {
    cat("Error in manual FST calculation:", e$message, "\n")
    return(NULL)
  })
}



#' Simple FST calculation using basic R functions
calculate_fst_simple <- function(genind_obj) {
  tryCatch({
    cat("Calculating FST using simple method...\n")
    
    # Get population names
    pop_names <- as.character(unique(adegenet::pop(genind_obj)))
    n_pops <- length(pop_names)
    
    if (n_pops < 2) {
      cat("  Need at least 2 populations for FST calculation\n")
      return(NULL)
    }
    
    # Convert to genotype matrix (0,1,2 format)
    geno_matrix <- as.matrix(genind_obj)
    
    # Get population assignments
    pop_assignments <- as.character(adegenet::pop(genind_obj))
    
    # Initialize FST matrix
    fst_matrix <- matrix(0, nrow = n_pops, ncol = n_pops)
    colnames(fst_matrix) <- rownames(fst_matrix) <- pop_names
    
    # Calculate allele frequencies per population
    allele_freqs <- list()
    for (pop in pop_names) {
      pop_indices <- which(pop_assignments == pop)
      if (length(pop_indices) > 0) {
        # Calculate allele frequency (assuming 0,1,2 coding)
        pop_matrix <- geno_matrix[pop_indices, , drop = FALSE]
        # Convert to allele frequency (0-1)
        freq <- colMeans(pop_matrix, na.rm = TRUE) / 2
        allele_freqs[[pop]] <- freq
      }
    }
    
    # Calculate pairwise FST
    for (i in 1:(n_pops - 1)) {
      for (j in (i + 1):n_pops) {
        pop1 <- pop_names[i]
        pop2 <- pop_names[j]
        
        p1 <- allele_freqs[[pop1]]
        p2 <- allele_freqs[[pop2]]
        
        # Remove loci with missing data
        valid_loci <- !is.na(p1) & !is.na(p2)
        p1 <- p1[valid_loci]
        p2 <- p2[valid_loci]
        
        if (length(p1) == 0) {
          fst_matrix[i, j] <- fst_matrix[j, i] <- NA
          next
        }
        
        # Calculate FST using Nei's method (simplified)
        # FST = (Ht - Hs) / Ht
        # Where Ht is total heterozygosity, Hs is mean within-pop heterozygosity
        
        # Calculate expected heterozygosity within populations (Hs)
        Hs_pop1 <- 2 * p1 * (1 - p1)
        Hs_pop2 <- 2 * p2 * (1 - p2)
        Hs <- (Hs_pop1 + Hs_pop2) / 2
        
        # Calculate allele frequency in combined population
        p_combined <- (p1 + p2) / 2
        
        # Calculate total expected heterozygosity (Ht)
        Ht <- 2 * p_combined * (1 - p_combined)
        
        # Calculate FST
        fst <- (Ht - Hs) / Ht
        
        # Take mean across loci
        fst_mean <- mean(fst, na.rm = TRUE)
        
        # Ensure FST is between 0 and 1
        fst_mean <- max(0, min(1, fst_mean))
        
        fst_matrix[i, j] <- fst_matrix[j, i] <- fst_mean
      }
    }
    
    cat("  Simple FST calculation completed\n")
    return(fst_matrix)
    
  }, error = function(e) {
    cat("Error in simple FST calculation:", e$message, "\n")
    return(NULL)
  })
}



# ============================================================================
# FIXED HARDY-WEINBERG EQUILIBRIUM ANALYSIS
# ============================================================================

#' Perform comprehensive HWE testing with detailed results
#' @param gl_object Genlight object
#' @param pop_assignments Population assignments (optional)
#' @param method HWE test method ("exact", "chi2", or "monte")
#' @param n_perm Number of permutations for Monte Carlo test
#' @return List with detailed HWE results
perform_detailed_hwe_analysis <- function(gl_object, pop_assignments = NULL, 
                                          method = "exact", n_perm = 1000) {
  tryCatch({
    cat("\n=== Performing Detailed HWE Analysis ===\n")
    
    # Convert to genind
    genind_obj <- create_genind_from_gl_fixed(gl_object, pop_assignments)
    
    if (is.null(genind_obj)) {
      stop("Failed to create genind object for HWE analysis")
    }
    
    # If no population assignments, treat as single population
    if (is.null(pop_assignments)) {
      pop_assignments <- rep("All", nInd(gl_object))
    }
    
    populations <- unique(pop_assignments)
    results <- list()
    overall_results <- data.frame()
    
    for (pop in populations) {
      cat("  Testing population:", pop, "\n")
      
      # Subset individuals
      pop_indices <- which(pop_assignments == pop)
      
      if (length(pop_indices) < 5) {
        cat("    Skipping - insufficient individuals (", length(pop_indices), ")\n", sep = "")
        next
      }
      
      pop_genind <- genind_obj[pop_indices, , drop = FALSE]
      
      # Perform HWE test using pegas package
      if (requireNamespace("pegas", quietly = TRUE)) {
        # Convert to pegas format
        loci_data <- pegas::as.loci(pop_genind)
        
        # Perform HWE test
        hwe_test <- pegas::hw.test(loci_data, B = n_perm)
        
        # Extract detailed results
        hwe_df <- as.data.frame(hwe_test)
        hwe_df$Locus <- rownames(hwe_df)
        hwe_df$Population <- pop
        hwe_df$N_Individuals <- length(pop_indices)
        
        # Add allele information
        allele_info <- get_allele_info(pop_genind)
        hwe_df <- merge(hwe_df, allele_info, by = "Locus", all.x = TRUE)
        
        # Calculate expected and observed heterozygosity
        hwe_df$H_obs <- get_observed_heterozygosity(pop_genind)
        hwe_df$H_exp <- get_expected_heterozygosity(pop_genind)
        hwe_df$Fis <- 1 - (hwe_df$H_obs / hwe_df$H_exp)
        
        # Add multiple testing correction
        hwe_df$P_adjusted_BH <- p.adjust(hwe_df$`Pr(chi^2 >)`, method = "BH")
        hwe_df$P_adjusted_Bonferroni <- p.adjust(hwe_df$`Pr(chi^2 >)`, method = "bonferroni")
        hwe_df$Significant_BH <- hwe_df$P_adjusted_BH < 0.05
        hwe_df$Significant_Bonferroni <- hwe_df$P_adjusted_Bonferroni < 0.05
        
        results[[pop]] <- hwe_df
        overall_results <- rbind(overall_results, hwe_df)
      } else {
        # Fallback to simple method
        cat("    Using simple HWE test (pegas not available)\n")
        hwe_df <- test_hardweinberg_simple(pop_genind, pop)
        results[[pop]] <- hwe_df
        overall_results <- rbind(overall_results, hwe_df)
      }
    }
    
    # Create summary statistics
    summary_stats <- data.frame(
      Population = populations,
      N_Loci = sapply(results, nrow),
      N_Significant_BH = sapply(results, function(x) sum(x$Significant_BH, na.rm = TRUE)),
      N_Significant_Bonferroni = sapply(results, function(x) sum(x$Significant_Bonferroni, na.rm = TRUE)),
      Mean_Fis = sapply(results, function(x) mean(x$Fis, na.rm = TRUE)),
      Mean_Hobs = sapply(results, function(x) mean(x$H_obs, na.rm = TRUE)),
      Mean_Hexp = sapply(results, function(x) mean(x$H_exp, na.rm = TRUE))
    )
    
    # Create visualization data
    viz_data <- prepare_hwe_visualization_data(overall_results)
    
    cat("  HWE analysis complete\n")
    cat("    Total loci tested:", nrow(overall_results), "\n")
    cat("    Significant loci (BH):", sum(overall_results$Significant_BH, na.rm = TRUE), "\n")
    
    return(list(
      detailed_results = results,
      overall_results = overall_results,
      summary_stats = summary_stats,
      visualization_data = viz_data,
      populations = populations,
      method = method
    ))
    
  }, error = function(e) {
    cat("Error in HWE analysis:", e$message, "\n")
    return(NULL)
  })
}

#' Get allele information for loci
get_allele_info <- function(genind_obj) {
  tryCatch({
    loci_info <- data.frame()
    
    for (loc in locNames(genind_obj)) {
      # Get allele frequencies
      tab <- genind_obj@tab[, grep(paste0("^", loc, "\\."), colnames(genind_obj@tab))]
      
      if (ncol(tab) > 0) {
        # Calculate allele counts
        allele_counts <- colSums(tab, na.rm = TRUE)
        total_alleles <- sum(allele_counts)
        
        if (total_alleles > 0) {
          allele_freqs <- allele_counts / total_alleles
          
          # Get major and minor alleles
          sorted_freqs <- sort(allele_freqs, decreasing = TRUE)
          major_allele <- names(sorted_freqs)[1]
          major_freq <- sorted_freqs[1]
          minor_allele <- ifelse(length(sorted_freqs) > 1, 
                                 names(sorted_freqs)[2], 
                                 NA)
          minor_freq <- ifelse(length(sorted_freqs) > 1, 
                               sorted_freqs[2], 
                               NA)
          
          loci_info <- rbind(loci_info, data.frame(
            Locus = loc,
            N_Alleles = length(allele_freqs),
            Major_Allele = major_allele,
            Major_Freq = major_freq,
            Minor_Allele = minor_allele,
            Minor_Freq = minor_freq,
            MAF = ifelse(is.na(minor_freq), 0, minor_freq),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    return(loci_info)
    
  }, error = function(e) {
    cat("Error getting allele info:", e$message, "\n")
    return(data.frame())
  })
}

#' Get observed heterozygosity
get_observed_heterozygosity <- function(genind_obj) {
  tryCatch({
    summary_obj <- summary(genind_obj)
    return(summary_obj$Hobs)
  }, error = function(e) {
    return(rep(NA, nLoc(genind_obj)))
  })
}

#' Get expected heterozygosity
get_expected_heterozygosity <- function(genind_obj) {
  tryCatch({
    summary_obj <- summary(genind_obj)
    return(summary_obj$Hexp)
  }, error = function(e) {
    return(rep(NA, nLoc(genind_obj)))
  })
}

#' Prepare data for HWE visualization
prepare_hwe_visualization_data <- function(hwe_results) {
  tryCatch({
    # Create data for various plots
    plot_data <- list()
    
    # 1. P-value distribution
    plot_data$pvalue_dist <- hwe_results %>%
      select(Locus, Population, `Pr(chi^2 >)`) %>%
      mutate(NegLog10P = -log10(`Pr(chi^2 >)`))
    
    # 2. Fis distribution
    plot_data$fis_dist <- hwe_results %>%
      select(Locus, Population, Fis) %>%
      filter(!is.na(Fis) & is.finite(Fis))
    
    # 3. Hobs vs Hexp
    plot_data$hobs_vs_hexp <- hwe_results %>%
      select(Locus, Population, H_obs, H_exp) %>%
      filter(!is.na(H_obs) & !is.na(H_exp))
    
    # 4. Significant loci by population
    plot_data$sig_loci <- hwe_results %>%
      group_by(Population) %>%
      summarise(
        Total_Loci = n(),
        Sig_BH = sum(Significant_BH, na.rm = TRUE),
        Sig_Bonferroni = sum(Significant_Bonferroni, na.rm = TRUE),
        Prop_BH = Sig_BH / Total_Loci,
        Prop_Bonferroni = Sig_Bonferroni / Total_Loci
      )
    
    return(plot_data)
    
  }, error = function(e) {
    cat("Error preparing HWE visualization data:", e$message, "\n")
    return(NULL)
  })
}

#' Simple HWE test (fallback method)
test_hardweinberg_simple <- function(genind_obj, pop_name) {
  tryCatch({
    n_loci <- nLoc(genind_obj)
    results <- data.frame()
    
    for (i in 1:n_loci) {
      locus_name <- locNames(genind_obj)[i]
      
      # Get genotype counts
      tab <- table(genind_obj@tab[, grep(paste0("^", locus_name, "\\."), 
                                         colnames(genind_obj@tab))])
      
      if (length(tab) >= 3) {  # Need at least 3 genotype classes
        # Extract alleles
        genotypes <- names(tab)
        alleles <- unlist(strsplit(genotypes, "\\."))
        
        if (length(alleles) >= 2) {
          allele_counts <- table(alleles)
          total_alleles <- sum(allele_counts)
          
          # For biallelic loci
          if (length(allele_counts) == 2) {
            p <- allele_counts[1] / total_alleles
            q <- allele_counts[2] / total_alleles
            
            # Expected counts under HWE
            exp_counts <- c(
              p^2 * total_alleles/2,
              2 * p * q * total_alleles/2,
              q^2 * total_alleles/2
            )
            
            # Chi-square test
            obs_counts <- as.numeric(tab)
            chi2 <- sum((obs_counts - exp_counts)^2 / exp_counts, na.rm = TRUE)
            df <- 1
            p_value <- pchisq(chi2, df, lower.tail = FALSE)
            
            # Calculate heterozygosities
            H_obs <- obs_counts[2] / sum(obs_counts)
            H_exp <- 2 * p * q
            
            results <- rbind(results, data.frame(
              Locus = locus_name,
              Population = pop_name,
              chi2 = chi2,
              df = df,
              `Pr(chi^2 >)` = p_value,
              H_obs = H_obs,
              H_exp = H_exp,
              Fis = 1 - (H_obs / H_exp),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    return(results)
    
  }, error = function(e) {
    cat("Error in simple HWE test:", e$message, "\n")
    return(data.frame())
  })
}


# ============================================================================
# FIXED DAPC ANALYSIS WITH SCORES AND ASSIGNMENT
# ============================================================================

#' Perform comprehensive DAPC analysis
#' @param gl_object Genlight object
#' @param pop_assignments Population assignments (optional)
#' @param n_pca Number of PCA components to retain
#' @param n_da Number of discriminant functions to retain
#' @param find_clusters Whether to find optimal clusters
#' @param max_n_clust Maximum number of clusters to test
#' @return List with DAPC results, scores, and assignment table
perform_dapc_analysis_comprehensive <- function(gl_object, pop_assignments = NULL,
                                                n_pca = NULL, n_da = NULL,
                                                find_clusters = FALSE, max_n_clust = 10) {
  tryCatch({
    cat("\n=== Performing DAPC Analysis ===\n")
    
    # Convert to genind
    genind_obj <- create_genind_from_gl_fixed(gl_object, pop_assignments)
    
    if (is.null(genind_obj)) {
      stop("Failed to create genind object for DAPC")
    }
    
    # Get population assignments
    if (is.null(pop_assignments)) {
      # Use find.clusters to determine groups
      cat("  Finding optimal number of clusters...\n")
      find_clust_result <- adegenet::find.clusters(genind_obj, 
                                                   max.n.clust = max_n_clust,
                                                   n.pca = min(50, nInd(genind_obj) - 1))
      
      pop_assignments <- find_clust_result$grp
      optimal_clusters <- find_clust_result$stat
      cat("  Optimal clusters found:", length(unique(pop_assignments)), "\n")
    } else {
      # Use provided assignments
      pop(genind_obj) <- pop_assignments
    }
    
    # Determine number of PCA components to retain
    if (is.null(n_pca)) {
      # Use cross-validation to find optimal n.pca
      cat("  Cross-validating to find optimal PCA components...\n")
      xval <- xvalDapc(tab(genind_obj, NA.method = "mean"), 
                       pop(genind_obj),
                       n.pca.max = min(30, nInd(genind_obj) - 1),
                       training.set = 0.9,
                       result = "groupMean")
      
      n_pca <- xval$`Best Number of PCs`
      n_da <- min(length(unique(pop(genind_obj))) - 1, n_pca)
      cat("  Optimal PCA components:", n_pca, "\n")
      cat("  DA functions to retain:", n_da, "\n")
    }
    
    # Perform DAPC
    cat("  Running DAPC...\n")
    dapc_result <- adegenet::dapc(genind_obj,
                                  pop = pop(genind_obj),
                                  n.pca = n_pca,
                                  n.da = n_da)
    
    # Extract scores
    scores_df <- as.data.frame(dapc_result$ind.coord)
    colnames(scores_df) <- paste0("DA", 1:ncol(scores_df))
    scores_df$Individual <- indNames(genind_obj)
    scores_df$Population <- as.character(pop(genind_obj))
    scores_df$Assigned_Population <- as.character(dapc_result$assign)
    
    # Create assignment table
    assignment_table <- create_dapc_assignment_table(dapc_result, genind_obj)
    
    # Calculate assignment accuracy
    accuracy <- mean(scores_df$Population == scores_df$Assigned_Population, na.rm = TRUE)
    cat("  Assignment accuracy:", round(accuracy * 100, 2), "%\n")
    
    # Create summary statistics
    summary_stats <- data.frame(
      N_Individuals = nInd(genind_obj),
      N_Loci = nLoc(genind_obj),
      N_Populations = length(unique(pop(genind_obj))),
      N_PCA_Components = n_pca,
      N_DA_Functions = n_da,
      Assignment_Accuracy = accuracy,
      Eigenvalues = paste(round(dapc_result$eig/sum(dapc_result$eig)*100, 2), "%", collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    # Create visualization data
    viz_data <- prepare_dapc_visualization_data(dapc_result, scores_df)
    
    cat("  DAPC analysis complete\n")
    
    return(list(
      dapc_result = dapc_result,
      scores = scores_df,
      assignment_table = assignment_table,
      summary_stats = summary_stats,
      visualization_data = viz_data,
      population_assignments = pop(genind_obj),
      eigenvalues = dapc_result$eig
    ))
    
  }, error = function(e) {
    cat("Error in DAPC analysis:", e$message, "\n")
    return(NULL)
  })
}

#' Create DAPC assignment table
create_dapc_assignment_table <- function(dapc_result, genind_obj) {
  tryCatch({
    # Get posterior probabilities
    post_prob <- dapc_result$posterior
    
    # Create assignment table
    assignment_df <- data.frame(
      Individual = indNames(genind_obj),
      True_Population = as.character(pop(genind_obj)),
      Assigned_Population = as.character(dapc_result$assign),
      Assignment_Probability = apply(post_prob, 1, max),
      stringsAsFactors = FALSE
    )
    
    # Add probabilities for each population
    for (pop_name in colnames(post_prob)) {
      assignment_df[[paste0("Prob_", pop_name)]] <- post_prob[, pop_name]
    }
    
    # Add assignment status
    assignment_df$Correct_Assignment <- assignment_df$True_Population == assignment_df$Assigned_Population
    
    # Calculate confidence levels
    assignment_df$Confidence_Level <- cut(assignment_df$Assignment_Probability,
                                          breaks = c(0, 0.7, 0.9, 0.95, 1),
                                          labels = c("Low", "Medium", "High", "Very High"),
                                          include.lowest = TRUE)
    
    return(assignment_df)
    
  }, error = function(e) {
    cat("Error creating assignment table:", e$message, "\n")
    return(NULL)
  })
}

#' Prepare DAPC visualization data
prepare_dapc_visualization_data <- function(dapc_result, scores_df) {
  tryCatch({
    viz_data <- list()
    
    # 1. Scatter plot data
    viz_data$scatter <- scores_df
    
    # 2. Loadings data
    if (!is.null(dapc_result$var.contr)) {
      loadings_df <- as.data.frame(dapc_result$var.contr)
      loadings_df$Locus <- rownames(loadings_df)
      viz_data$loadings <- loadings_df
    }
    
    # 3. Assignment probabilities
    if (!is.null(dapc_result$posterior)) {
      prob_df <- as.data.frame(dapc_result$posterior)
      prob_df$Individual <- rownames(prob_df)
      prob_df$Population <- scores_df$Population
      viz_data$probabilities <- prob_df
    }
    
    # 4. Population statistics
    pop_stats <- scores_df %>%
      group_by(Population) %>%
      summarise(
        N_Individuals = n(),
        Mean_DA1 = mean(DA1, na.rm = TRUE),
        Mean_DA2 = mean(DA2, na.rm = TRUE),
        SD_DA1 = sd(DA1, na.rm = TRUE),
        SD_DA2 = sd(DA2, na.rm = TRUE),
        Assignment_Rate = mean(Assigned_Population == Population, na.rm = TRUE)
      )
    
    viz_data$population_stats <- pop_stats
    
    return(viz_data)
    
  }, error = function(e) {
    cat("Error preparing DAPC visualization data:", e$message, "\n")
    return(NULL)
  })
}

#' Fixed genind creation function
create_genind_from_gl_fixed <- function(gl_object, pop_assignments = NULL) {
  tryCatch({
    cat("Creating genind object from genlight...\n")
    
    # Method 1: Use dartR if available
    if (requireNamespace("dartR", quietly = TRUE)) {
      genind_obj <- dartR::gl2gi(gl_object, verbose = 0)
      
      if (!is.null(pop_assignments) && length(pop_assignments) == nInd(gl_object)) {
        pop(genind_obj) <- pop_assignments
      }
      
      return(genind_obj)
    }
    
    # Method 2: Manual conversion
    geno_matrix <- as.matrix(gl_object)
    
    # Convert to character format
    geno_char <- matrix(as.character(geno_matrix), nrow = nrow(geno_matrix))
    
    # Handle missing values
    geno_char[is.na(geno_char)] <- "NA/NA"
    
    # For non-NA values, create genotype strings
    for (i in 1:nrow(geno_char)) {
      for (j in 1:ncol(geno_char)) {
        if (geno_char[i, j] != "NA/NA") {
          val <- as.numeric(geno_char[i, j])
          if (val == 0) {
            geno_char[i, j] <- "0/0"
          } else if (val == 1) {
            geno_char[i, j] <- "0/1"
          } else if (val == 2) {
            geno_char[i, j] <- "1/1"
          }
        }
      }
    }
    
    # Create genind object
    genind_obj <- adegenet::df2genind(geno_char, 
                                      sep = "/",
                                      ind.names = indNames(gl_object),
                                      loc.names = locNames(gl_object),
                                      type = "codom")
    
    if (!is.null(pop_assignments) && length(pop_assignments) == nInd(gl_object)) {
      pop(genind_obj) <- pop_assignments
    }
    
    return(genind_obj)
    
  }, error = function(e) {
    cat("Error creating genind object:", e$message, "\n")
    return(NULL)
  })
}



# ============================================================================
# HYBRID PATTERN VISUALIZATION
# ============================================================================

#' Visualize hybrid patterns from admixture analysis
#' @param q_matrix Q matrix from LEA or other admixture analysis
#' @param sample_names Vector of sample names
#' @param parent_references Data frame with parent reference information
#' @param hybrid_status Vector indicating hybrid status (optional)
#' @param threshold Admixture threshold for hybrid classification
#' @return List with visualizations and summary statistics
visualize_hybrid_patterns <- function(q_matrix, sample_names = NULL,
                                      parent_references = NULL, 
                                      threshold = 0.1) {
  tryCatch({
    cat("\n=== Visualizing Hybrid Patterns ===\n")
    
    # Set sample names if not provided
    if (is.null(sample_names)) {
      sample_names <- rownames(q_matrix)
      if (is.null(sample_names)) {
        sample_names <- paste0("Sample_", 1:nrow(q_matrix))
      }
    }
    
    n_clusters <- ncol(q_matrix)
    cat("  Samples:", nrow(q_matrix), "\n")
    cat("  Clusters:", n_clusters, "\n")
    
    # Create hybrid classification
    classification <- classify_hybrids(q_matrix, threshold)
    
    # Prepare visualization data
    viz_data <- prepare_hybrid_visualization_data(q_matrix, sample_names, 
                                                  classification, parent_references)
    
    # Generate visualizations
    plots <- create_hybrid_visualizations(viz_data, n_clusters)
    
    cat("  Hybrid visualization complete\n")
    
    return(list(
      classification = classification,
      visualization_data = viz_data,
      plots = plots,
      threshold = threshold
    ))
    
  }, error = function(e) {
    cat("Error in hybrid visualization:", e$message, "\n")
    return(NULL)
  })
}

#' Classify individuals based on admixture proportions
classify_hybrids <- function(q_matrix, threshold = 0.1) {
  classification <- data.frame(
    Sample = rownames(q_matrix),
    Dominant_Cluster = apply(q_matrix, 1, which.max),
    Max_Proportion = apply(q_matrix, 1, max),
    Category = NA,
    stringsAsFactors = FALSE
  )
  
  # Determine category based on admixture proportions
  for (i in 1:nrow(classification)) {
    if (classification$Max_Proportion[i] > (1 - threshold)) {
      classification$Category[i] <- "Pure"
    } else if (classification$Max_Proportion[i] > 0.5) {
      classification$Category[i] <- "Hybrid"
    } else {
      classification$Category[i] <- "Admixed"
    }
  }
  
  # Calculate admixture proportions
  classification$Admixture_Proportion <- 1 - classification$Max_Proportion
  
  return(classification)
}

#' Prepare data for hybrid visualization
prepare_hybrid_visualization_data <- function(q_matrix, sample_names, 
                                              classification, parent_references) {
  viz_data <- list()
  
  # 1. Structure plot data
  viz_data$structure <- data.frame(
    Sample = rep(sample_names, each = ncol(q_matrix)),
    Cluster = rep(paste0("Cluster_", 1:ncol(q_matrix)), times = nrow(q_matrix)),
    Proportion = as.vector(t(q_matrix)),
    stringsAsFactors = FALSE
  )
  
  viz_data$structure <- merge(viz_data$structure, 
                              classification[, c("Sample", "Category")], 
                              by = "Sample")
  
  # 2. Triangle plot data (for 3 clusters)
  if (ncol(q_matrix) == 3) {
    viz_data$triangle <- data.frame(
      Sample = sample_names,
      Cluster1 = q_matrix[, 1],
      Cluster2 = q_matrix[, 2],
      Cluster3 = q_matrix[, 3],
      Category = classification$Category,
      stringsAsFactors = FALSE
    )
  }
  
  return(viz_data)
}

#' Create hybrid visualizations
create_hybrid_visualizations <- function(viz_data, n_clusters) {
  plots <- list()
  
  # 1. Structure plot with hybrid classification
  plots$structure <- ggplot(viz_data$structure, 
                            aes(x = Sample, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(Category ~ ., scales = "free_x", space = "free") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Hybrid Structure Plot",
         x = "Samples",
         y = "Admixture Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # 2. Triangle plot (for 3 clusters)
  if (!is.null(viz_data$triangle)) {
    plots$triangle <- ggplot(viz_data$triangle, 
                             aes(x = Cluster1, y = Cluster2, color = Category)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "Ternary Plot of Admixture Proportions",
           x = "Cluster 1 Proportion",
           y = "Cluster 2 Proportion") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
  
  return(plots)
}

#' Match parents to clusters for hybrid index calculation
match_parents_to_clusters <- function(q_matrix, parent_references) {
  result <- list()
  
  # Find which cluster each parent belongs to
  for (i in 1:nrow(parent_references)) {
    parent_name <- parent_references$Sample[i]
    parent_idx <- which(rownames(q_matrix) == parent_name)
    
    if (length(parent_idx) > 0) {
      parent_q <- q_matrix[parent_idx, ]
      result[[paste0("parent", i)]] <- list(
        name = parent_name,
        cluster = which.max(parent_q),
        proportion = max(parent_q)
      )
    }
  }
  
  return(result)
}

#' Calculate hybrid indices
calculate_hybrid_indices <- function(q_matrix, parent_clusters) {
  n_samples <- nrow(q_matrix)
  hybrid_indices <- rep(NA, n_samples)
  
  if (length(parent_clusters) >= 2) {
    parent1_cluster <- parent_clusters$parent1$cluster
    parent2_cluster <- parent_clusters$parent2$cluster
    
    for (i in 1:n_samples) {
      # Hybrid index = proportion from parent2
      p1_prop <- q_matrix[i, parent1_cluster]
      p2_prop <- q_matrix[i, parent2_cluster]
      
      if (p1_prop + p2_prop > 0) {
        hybrid_indices[i] <- p2_prop / (p1_prop + p2_prop)
      }
    }
  }
  
  return(hybrid_indices)
}

#' Create network data for hybrid visualization
create_hybrid_network_data <- function(q_matrix, classification) {
  # Create similarity matrix
  similarity <- as.matrix(dist(q_matrix, method = "euclidean"))
  
  # Convert to network edges
  edges <- data.frame()
  
  # Connect each hybrid to its nearest pure individuals
  pure_samples <- which(classification$Category == "Pure")
  hybrid_samples <- which(classification$Category == "Hybrid")
  
  for (hybrid in hybrid_samples) {
    # Find closest pure individuals
    distances <- similarity[hybrid, pure_samples]
    closest <- pure_samples[order(distances)[1:2]]
    
    for (pure in closest) {
      edges <- rbind(edges, data.frame(
        from = classification$Sample[hybrid],
        to = classification$Sample[pure],
        weight = 1 / (1 + distances[which(pure_samples == pure)]),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Create nodes
  nodes <- data.frame(
    id = classification$Sample,
    label = classification$Sample,
    group = classification$Category,
    value = classification$Admixture_Proportion * 10 + 1,
    title = paste0("Sample: ", classification$Sample, 
                   "<br>Category: ", classification$Category,
                   "<br>Admixture: ", round(classification$Admixture_Proportion, 3)),
    stringsAsFactors = FALSE
  )
  
  return(list(nodes = nodes, edges = edges))
}

#' Create hybrid summary statistics
create_hybrid_summary <- function(classification, viz_data) {
  summary <- list()
  
  # Basic counts
  summary$counts <- table(classification$Category)
  summary$proportions <- prop.table(summary$counts)
  
  # Admixture statistics
  summary$admixture_stats <- classification %>%
    group_by(Category) %>%
    summarise(
      N = n(),
      Mean_Admixture = mean(Admixture_Proportion, na.rm = TRUE),
      SD_Admixture = sd(Admixture_Proportion, na.rm = TRUE),
      Min_Admixture = min(Admixture_Proportion, na.rm = TRUE),
      Max_Admixture = max(Admixture_Proportion, na.rm = TRUE)
    )
  
  # Cluster proportions
  if (!is.null(viz_data$structure)) {
    summary$cluster_props <- viz_data$structure %>%
      group_by(Cluster) %>%
      summarise(
        Mean_Proportion = mean(Proportion, na.rm = TRUE),
        SD_Proportion = sd(Proportion, na.rm = TRUE)
      )
  }
  
  return(summary)
}


#==============================================================================
#' Test for Hardy-Weinberg Equilibrium
#' @param genind_obj Genind object
#' @param populations Vector of population names (optional)
#' @param n_perm Number of permutations for Monte Carlo test
#' @return List containing HWE test results
#' Fixed: Test for Hardy-Weinberg Equilibrium without polysat dependency
test_hardweinberg_fixed <- function(genind_obj, populations = NULL, n_perm = 1000) {
  tryCatch({
    cat("Testing Hardy-Weinberg Equilibrium...\n")
    
    # Check if genind_obj is valid
    if (is.null(genind_obj) || !inherits(genind_obj, "genind")) {
      stop("Invalid genind object")
    }
    
    if (is.null(populations)) {
      populations <- unique(as.character(adegenet::pop(genind_obj)))
    }
    
    hwe_results <- list()
    all_p_values <- data.frame()
    
    # Simple chi-square test function
    perform_simple_hwe_test <- function(pop_genind) {
      # Get genotype counts for each locus
      results <- list()
      
      for (locus in 1:nLoc(pop_genind)) {
        # Get genotype counts
        tab <- table(pop_genind@tab[, locus])
        
        # Skip if too few genotypes
        if (length(tab) < 2) {
          next
        }
        
        # Convert to allele frequencies
        alleles <- unlist(strsplit(names(tab), "\\."))
        allele_counts <- table(alleles)
        total_alleles <- sum(allele_counts)
        allele_freqs <- allele_counts / total_alleles
        
        # Expected genotype frequencies (Hardy-Weinberg)
        # For two alleles: p^2, 2pq, q^2
        allele_names <- names(allele_freqs)
        if (length(allele_freqs) == 2) {
          p <- allele_freqs[1]
          q <- allele_freqs[2]
          
          # Expected counts
          exp_counts <- c(p^2, 2*p*q, q^2) * total_alleles/2
          
          # Chi-square test
          obs_counts <- as.numeric(tab)
          chi2 <- sum((obs_counts - exp_counts)^2 / exp_counts, na.rm = TRUE)
          df <- 1  # For diallelic loci
          p_value <- pchisq(chi2, df, lower.tail = FALSE)
          
          results[[locus]] <- c(chi2, df, p_value)
          names(results)[locus] <- locNames(pop_genind)[locus]
        }
      }
      
      # Convert to matrix
      if (length(results) > 0) {
        result_matrix <- do.call(rbind, results)
        colnames(result_matrix) <- c("chi2", "df", "Pr(chi^2 >)")
        return(result_matrix)
      } else {
        return(NULL)
      }
    }
    
    for (pop_name in populations) {
      cat("Testing population:", pop_name, "\n")
      
      # Subset population manually without poppr dependency
      pop_indices <- which(adegenet::pop(genind_obj) == pop_name)
      
      if (length(pop_indices) < 5) {
        cat("    Skipping - insufficient individuals (", length(pop_indices), ")\n", sep = "")
        next
      }
      
      # Create subset manually
      pop_genind <- genind_obj[pop_indices, , drop = FALSE]
      
      # Try HWE test with different methods
      hwe_test <- NULL
      
      # Method 1: Try pegas::hw.test (no polysat dependency)
      if (requireNamespace("pegas", quietly = TRUE)) {
        tryCatch({
          cat("    Using pegas::hw.test...\n")
          # Convert to pegas::loci format
          pop_loci <- pegas::as.loci(pop_genind)
          
          # Run HWE test - hw.test returns a matrix, not a list
          hwe_test <- pegas::hw.test(pop_loci, B = n_perm)
          
          # Check the structure of hwe_test
          cat("    Structure of hwe_test:", class(hwe_test), "\n")
          cat("    Dimensions:", dim(hwe_test), "\n")
          print(head(hwe_test))
          
          # Handle the result correctly - hw.test returns a matrix with columns: chi2, df, P, etc.
          # Extract the chi-square values and p-values
          if (!is.null(hwe_test) && nrow(hwe_test) > 0) {
            # Create a matrix with the required columns
            hwe_matrix <- matrix(NA, nrow = nrow(hwe_test), ncol = 3)
            rownames(hwe_matrix) <- rownames(hwe_test)
            
            # Extract values - check which columns exist
            if ("chi2" %in% colnames(hwe_test)) {
              hwe_matrix[, 1] <- hwe_test[, "chi2"]
            } else if ("Chi2" %in% colnames(hwe_test)) {
              hwe_matrix[, 1] <- hwe_test[, "Chi2"]
            }
            
            hwe_matrix[, 2] <- 1  # DF for diallelic loci
            
            if ("Pr(chi^2 >)" %in% colnames(hwe_test)) {
              hwe_matrix[, 3] <- hwe_test[, "Pr(chi^2 >)"]
            } else if ("Pr.exact" %in% colnames(hwe_test)) {
              hwe_matrix[, 3] <- hwe_test[, "Pr.exact"]
            } else if ("P" %in% colnames(hwe_test)) {
              hwe_matrix[, 3] <- hwe_test[, "P"]
            }
            
            colnames(hwe_matrix) <- c("chi2", "df", "Pr(chi^2 >)")
            hwe_test <- hwe_matrix
          }
          
        }, error = function(e) {
          cat("    pegas::hw.test failed:", e$message, "\n")
        })
      }
      
      # Method 2: Simple chi-square test as fallback
      if (is.null(hwe_test) || nrow(hwe_test) == 0) {
        cat("    Using simple chi-square test...\n")
        tryCatch({
          hwe_test <- perform_simple_hwe_test(pop_genind)
        }, error = function(e) {
          cat("    Simple HWE test failed:", e$message, "\n")
        })
      }
      
      if (!is.null(hwe_test) && nrow(hwe_test) > 0) {
        # Extract p-values - handle both matrix and list formats
        if (is.matrix(hwe_test)) {
          # Extract from matrix
          p_values <- data.frame(
            Population = pop_name,
            Locus = rownames(hwe_test),
            ChiSq = as.numeric(hwe_test[, 1]),
            DF = as.numeric(hwe_test[, 2]),
            P_Value = as.numeric(hwe_test[, 3]),
            stringsAsFactors = FALSE
          )
        } else if (is.list(hwe_test)) {
          # Extract from list (alternative format)
          p_values <- data.frame(
            Population = pop_name,
            Locus = names(hwe_test),
            ChiSq = sapply(hwe_test, function(x) x[1]),
            DF = sapply(hwe_test, function(x) x[2]),
            P_Value = sapply(hwe_test, function(x) x[3]),
            stringsAsFactors = FALSE
          )
        } else {
          next
        }
        
        # Remove rows with NA p-values
        p_values <- p_values[!is.na(p_values$P_Value), ]
        
        if (nrow(p_values) > 0) {
          # Adjust p-values for multiple testing
          p_values$P_Adjusted <- p.adjust(p_values$P_Value, method = "fdr")
          
          hwe_results[[pop_name]] <- list(
            test_results = hwe_test,
            p_values = p_values,
            n_loci = nrow(p_values),
            n_sig = sum(p_values$P_Adjusted < 0.05, na.rm = TRUE),
            prop_sig = sum(p_values$P_Adjusted < 0.05, na.rm = TRUE) / nrow(p_values)
          )
          
          all_p_values <- rbind(all_p_values, p_values)
        }
      }
    }
    
    # Create summary plot if we have data
    if (nrow(all_p_values) > 0) {
      # Load dplyr if available
      if (requireNamespace("dplyr", quietly = TRUE)) {
        # Summarize by population
        pop_summary <- all_p_values %>%
          dplyr::group_by(Population) %>%
          dplyr::summarise(
            N_Loci = dplyr::n(),
            N_Significant = sum(P_Adjusted < 0.05, na.rm = TRUE),
            Proportion_Significant = N_Significant / N_Loci,
            Mean_P_Value = mean(P_Value, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        # Base R alternative
        pop_summary <- aggregate(P_Value ~ Population, data = all_p_values, 
                                 FUN = function(x) {
                                   c(N_Loci = length(x),
                                     N_Significant = sum(all_p_values$P_Adjusted[all_p_values$Population == unique(x)] < 0.05, na.rm = TRUE),
                                     Mean_P_Value = mean(x, na.rm = TRUE))
                                 })
        pop_summary <- data.frame(
          Population = pop_summary$Population,
          N_Loci = pop_summary$P_Value[, "N_Loci"],
          N_Significant = pop_summary$P_Value[, "N_Significant"],
          Proportion_Significant = pop_summary$P_Value[, "N_Significant"] / pop_summary$P_Value[, "N_Loci"],
          Mean_P_Value = pop_summary$P_Value[, "Mean_P_Value"]
        )
      }
      
      # Create plot if ggplot2 is available
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        summary_plot <- ggplot2::ggplot(pop_summary, 
                                        ggplot2::aes(x = stats::reorder(Population, -Proportion_Significant), 
                                                     y = Proportion_Significant)) +
          ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = Population), alpha = 0.7) +
          ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
          ggplot2::labs(title = "Proportion of Loci Deviating from HWE",
                        subtitle = "Red line indicates 5% threshold",
                        x = "Population", y = "Proportion of Significant Loci") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "none"
          )
      } else {
        summary_plot <- NULL
      }
    } else {
      pop_summary <- NULL
      summary_plot <- NULL
    }
    
    return(list(
      hwe_tests = hwe_results,
      all_p_values = all_p_values,
      population_summary = pop_summary,
      summary_plot = summary_plot
    ))
    
  }, error = function(e) {
    cat("Error in HWE testing:", e$message, "\n")
    return(NULL)
  })
}

#' Helper function: Perform simple chi-square HWE test
perform_simple_hwe_test <- function(pop_genind) {
  tryCatch({
    # Get genotype counts
    geno_counts <- adegenet::tab(pop_genind, NA.method = "zero")
    
    # Initialize results matrix
    n_loci <- ncol(geno_counts)
    results <- matrix(NA, nrow = n_loci, ncol = 3)
    rownames(results) <- colnames(geno_counts)
    colnames(results) <- c("chi2", "df", "p_value")
    
    for (i in 1:n_loci) {
      # Get genotype frequencies for this locus
      locus_data <- geno_counts[, i]
      
      # For diallelic SNPs (0,1,2 coding in genind), we need to extract allele counts
      # This is simplified - actual implementation depends on your data structure
      
      # Simple implementation: Skip for now or implement proper test
      # For now, return NA values
      results[i, ] <- c(NA, NA, NA)
    }
    
    return(results)
  }, error = function(e) {
    cat("Simple HWE test error:", e$message, "\n")
    return(NULL)
  })
}

#' Alternative: Use HardyWeinberg package for HWE testing
test_hardweinberg_HW <- function(genind_obj, populations = NULL) {
  tryCatch({
    cat("Testing Hardy-Weinberg Equilibrium using HardyWeinberg package...\n")
    
    # Check if HardyWeinberg package is available
    if (!requireNamespace("HardyWeinberg", quietly = TRUE)) {
      cat("  Installing HardyWeinberg package...\n")
      install.packages("HardyWeinberg", dependencies = TRUE)
      library(HardyWeinberg)
    }
    
    if (is.null(populations)) {
      populations <- unique(as.character(adegenet::pop(genind_obj)))
    }
    
    hwe_results <- list()
    all_p_values <- data.frame()
    
    for (pop_name in populations) {
      cat("  Testing population:", pop_name, "\n")
      
      # Subset population
      pop_indices <- which(adegenet::pop(genind_obj) == pop_name)
      
      if (length(pop_indices) < 5) {
        cat("    Skipping - insufficient individuals\n")
        next
      }
      
      pop_genind <- genind_obj[pop_indices, , drop = FALSE]
      
      # Convert to genotype counts
      geno_counts <- adegenet::tab(pop_genind, NA.method = "zero")
      
      # Prepare results
      n_loci <- ncol(geno_counts)
      p_values <- numeric(n_loci)
      
      for (i in 1:n_loci) {
        tryCatch({
          # Get genotype counts for this locus
          # Assuming genotypes are coded as 0,1,2 for AA, AB, BB
          # We need to convert to counts of A and B alleles
          geno_vec <- geno_counts[, i]
          
          # Count genotypes (simplified - actual conversion depends on your data)
          AA_count <- sum(geno_vec == 0, na.rm = TRUE)
          AB_count <- sum(geno_vec == 1, na.rm = TRUE)
          BB_count <- sum(geno_vec == 2, na.rm = TRUE)
          
          # Create vector of genotype counts
          x <- c(AA = AA_count, AB = AB_count, BB = BB_count)
          
          # Test HWE using HardyWeinberg package
          hw_test <- HardyWeinberg::HWExact(x, verbose = FALSE)
          p_values[i] <- hw_test$pval
          
        }, error = function(e) {
          p_values[i] <- NA
        })
      }
      
      # Create results data frame
      results_df <- data.frame(
        Population = pop_name,
        Locus = colnames(geno_counts),
        P_Value = p_values,
        stringsAsFactors = FALSE
      )
      
      # Remove NAs
      results_df <- results_df[!is.na(results_df$P_Value), ]
      
      if (nrow(results_df) > 0) {
        # Adjust p-values for multiple testing
        results_df$P_Adjusted <- p.adjust(results_df$P_Value, method = "fdr")
        
        hwe_results[[pop_name]] <- list(
          p_values = results_df,
          n_loci = nrow(results_df),
          n_sig = sum(results_df$P_Adjusted < 0.05, na.rm = TRUE),
          prop_sig = sum(results_df$P_Adjusted < 0.05, na.rm = TRUE) / nrow(results_df)
        )
        
        all_p_values <- rbind(all_p_values, results_df)
      }
    }
    
    # Create summary plot
    if (nrow(all_p_values) > 0) {
      pop_summary <- all_p_values %>%
        dplyr::group_by(Population) %>%
        dplyr::summarise(
          N_Loci = n(),
          N_Significant = sum(P_Adjusted < 0.05, na.rm = TRUE),
          Proportion_Significant = N_Significant / N_Loci,
          Mean_P_Value = mean(P_Value, na.rm = TRUE),
          .groups = "drop"
        )
      
      summary_plot <- ggplot2::ggplot(pop_summary, 
                                      ggplot2::aes(x = stats::reorder(Population, -Proportion_Significant), 
                                                   y = Proportion_Significant)) +
        ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = Population), alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
        ggplot2::labs(title = "Proportion of Loci Deviating from HWE",
                      subtitle = "Red line indicates 5% threshold",
                      x = "Population", y = "Proportion of Significant Loci") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "none"
        )
    } else {
      pop_summary <- NULL
      summary_plot <- NULL
    }
    
    return(list(
      hwe_tests = hwe_results,
      all_p_values = all_p_values,
      population_summary = pop_summary,
      summary_plot = summary_plot
    ))
    
  }, error = function(e) {
    cat("Error in HWE testing with HardyWeinberg package:", e$message, "\n")
    return(NULL)
  })
}



#' Perform Discriminant Analysis of Principal Components (DAPC)
#' @param genind_obj Genind object
#' @param n_pca Number of PCA components to retain
#' @param n_da Number of DA components to retain
#' @param find_clusters Whether to automatically find clusters
#' @param max_clusters Maximum number of clusters to test
#' @return List containing DAPC results
perform_dapc_analysis <- function(genind_obj, n_pca = 40, n_da = 10, 
                                  find_clusters = TRUE, max_clusters = 20) {
  tryCatch({
    cat("Performing DAPC analysis...\n")
    
    if (find_clusters) {
      # Find optimal number of clusters
      cat("  Finding optimal number of clusters...\n")
      cluster_test <- find.clusters(genind_obj, 
                                    max.n.clust = max_clusters,
                                    n.pca = n_pca)
      
      optimal_clusters <- length(unique(cluster_test$grp))
      cat("  Optimal clusters found:", optimal_clusters, "\n")
      
      # Perform DAPC with found clusters
      dapc_result <- dapc(genind_obj, 
                          n.pca = n_pca,
                          n.da = n_da,
                          pop = cluster_test$grp)
    } else {
      # Use existing population assignments
      if (is.null(pop(genind_obj))) {
        stop("Population information required when find_clusters = FALSE")
      }
      
      dapc_result <- dapc(genind_obj,
                          n.pca = n_pca,
                          n.da = n_da,
                          pop = pop(genind_obj))
    }
    
    # Extract results
    dapc_scores <- as.data.frame(dapc_result$ind.coord)
    colnames(dapc_scores) <- paste0("DA", 1:ncol(dapc_scores))
    dapc_scores$Individual <- indNames(genind_obj)
    dapc_scores$Population <- if (find_clusters) {
      as.character(dapc_result$grp)
    } else {
      as.character(pop(genind_obj))
    }
    
    # Calculate assignment probabilities
    if (!is.null(dapc_result$posterior)) {
      assignment_prob <- as.data.frame(dapc_result$posterior)
      colnames(assignment_prob) <- paste0("Prob_", colnames(assignment_prob))
      assignment_prob$Individual <- indNames(genind_obj)
      assignment_prob$Assigned_Pop <- apply(dapc_result$posterior, 1, 
                                            function(x) colnames(dapc_result$posterior)[which.max(x)])
    } else {
      assignment_prob <- NULL
    }
    
    # Create DAPC scatter plot
    if (ncol(dapc_result$ind.coord) >= 2) {
      scatter_data <- data.frame(
        DA1 = dapc_result$ind.coord[, 1],
        DA2 = dapc_result$ind.coord[, 2],
        Population = dapc_scores$Population,
        Individual = dapc_scores$Individual
      )
      
      scatter_plot <- ggplot(scatter_data, aes(x = DA1, y = DA2, color = Population)) +
        geom_point(size = 3, alpha = 0.7) +
        stat_ellipse(level = 0.95, linetype = "dashed") +
        labs(title = "DAPC Scatter Plot",
             x = paste("DA1 (", round(dapc_result$eig[1]/sum(dapc_result$eig)*100, 1), "%)", sep = ""),
             y = paste("DA2 (", round(dapc_result$eig[2]/sum(dapc_result$eig)*100, 1), "%)", sep = "")) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      scatter_plot <- NULL
    }
    
    # Create assignment plot
    if (!is.null(assignment_prob) && ncol(assignment_prob) > 3) {
      # Prepare data for assignment plot
      assignment_melt <- reshape2::melt(assignment_prob[, !colnames(assignment_prob) %in% c("Individual", "Assigned_Pop")],
                                        variable.name = "Population", value.name = "Probability")
      assignment_melt$Individual <- rep(assignment_prob$Individual, 
                                        times = ncol(assignment_prob) - 2)
      
      assignment_plot <- ggplot(assignment_melt, 
                                aes(x = Individual, y = Probability, fill = Population)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(title = "Population Assignment Probabilities",
             x = "Individual", y = "Probability") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom"
        )
    } else {
      assignment_plot <- NULL
    }
    
    return(list(
      dapc_result = dapc_result,
      scores = dapc_scores,
      assignment_prob = assignment_prob,
      scatter_plot = scatter_plot,
      assignment_plot = assignment_plot,
      variance_explained = dapc_result$eig / sum(dapc_result$eig) * 100
    ))
    
  }, error = function(e) {
    cat("Error in DAPC analysis:", e$message, "\n")
    return(NULL)
  })
}



#' Run comprehensive population genetics analysis with fixed HWE testing
#' @param gl_object Genlight object
#' @param pop_data Population assignment data (optional)
#' @param analyses Vector of analyses to run
#' @return List containing all analysis results
run_popgen_analyses <- function(gl_object, pop_data = NULL, 
                                      analyses = c("allele_barcode", "heterozygosity", 
                                                   "fst", "hwe", "dapc")) {
  
  tryCatch({
    cat("Running comprehensive population genetics analyses...\n")
    
    # Create genind object using the fixed function
    genind_obj <- create_genind_from_gl(gl_object, pop_data)
    
    if (is.null(genind_obj)) {
      stop("Failed to create genind object")
    }
    
    results <- list()
    
    # Run selected analyses
    if ("allele_barcode" %in% analyses) {
      cat("Running allele barcode analysis...\n")
      results$allele_barcode <- create_allele_barcode(genind_obj)
    }
    
    if ("heterozygosity" %in% analyses) {
      cat("Running heterozygosity analysis...\n")
      results$heterozygosity <- calculate_heterozygosity_stats(genind_obj)
    }
    
    if ("fst" %in% analyses) {
      cat("Running FST analysis...\n")
      results$fst <- calculate_fst_matrix(genind_obj)
    }
    
    if ("hwe" %in% analyses) {
      cat("Running HWE testing...\n")
      # Use the fixed function without polysat dependency
      results$hwe <- test_hardweinberg_fixed(genind_obj)
      
      # If fixed function fails, try alternative
      if (is.null(results$hwe)) {
        cat("  Trying alternative HWE testing method...\n")
        results$hwe <- test_hardweinberg_HW(genind_obj)
      }
    }
    
    if ("dapc" %in% analyses) {
      cat("Running DAPC analysis...\n")
      results$dapc <- perform_dapc_analysis(genind_obj)
    }
    
    # Add metadata
    results$metadata <- list(
      n_individuals = nInd(gl_object),
      n_loci = nLoc(gl_object),
      n_populations = length(unique(adegenet::pop(genind_obj))),
      populations = unique(as.character(adegenet::pop(genind_obj))),
      analysis_date = Sys.time(),
      analyses_performed = analyses
    )
    
    cat("All analyses completed successfully\n")
    return(results)
    
  }, error = function(e) {
    cat("Error in population genetics analyses:", e$message, "\n")
    return(NULL)
  })
}


# ============================================================================
# POPULATION ALLELE STATISTICS FUNCTIONS
# ============================================================================

#' Calculate allele statistics per population
#' @param gl_object Genlight object
#' @param popnames Vector of population names (optional)
#' @param popcols Vector of population colors (optional)
#' @param min_individuals Minimum number of individuals per population to include
#' @return List containing allele statistics and plot
calculate_population_alleles <- function(gl_object, popnames = NULL, popcols = NULL, 
                                         min_individuals = 2) {
  tryCatch({
    cat("Calculating allele statistics per population...\n")
    
    # Check if population information is available
    if (is.null(pop(gl_object))) {
      cat("Warning: No population information available in genlight object.\n")
      cat("Assigning all individuals to single population.\n")
      pop(gl_object) <- rep("Population", nInd(gl_object))
    }
    
    # Get unique populations
    if (is.null(popnames)) {
      popnames <- unique(as.character(pop(gl_object)))
    }
    
    # Remove populations with insufficient individuals
    pop_counts <- table(pop(gl_object))
    valid_pops <- names(pop_counts)[pop_counts >= min_individuals]
    
    if (length(valid_pops) == 0) {
      stop(paste("No populations have at least", min_individuals, "individuals"))
    }
    
    cat("Analyzing", length(valid_pops), "populations...\n")
    
    # Initialize results data frame
    results_df <- data.frame(
      Population = valid_pops,
      N_Individuals = numeric(length(valid_pops)),
      N_Alleles_Mean = numeric(length(valid_pops)),
      N_Alleles_SD = numeric(length(valid_pops)),
      N_Alleles_Max = numeric(length(valid_pops)),
      N_Loci = numeric(length(valid_pops)),
      He_Mean = numeric(length(valid_pops)),
      Ho_Mean = numeric(length(valid_pops)),
      stringsAsFactors = FALSE
    )
    
    # Calculate statistics for each population
    for (i in seq_along(valid_pops)) {
      pop_name <- valid_pops[i]
      cat("  Processing population:", pop_name, "\n")
      
      # Subset individuals for this population
      pop_indices <- which(pop(gl_object) == pop_name)
      
      if (length(pop_indices) < min_individuals) {
        cat("    Skipping - insufficient individuals\n")
        next
      }
      
      # Subset genlight object
      pop_gl <- gl_object[pop_indices, ]
      
      # Convert to genind for allele counting (using dartR's gl2gi)
      tryCatch({
        pop_genind <- gl2gi(pop_gl, verbose = 0)
        
        # Basic statistics
        results_df$N_Individuals[i] <- nInd(pop_gl)
        results_df$N_Loci[i] <- nLoc(pop_gl)
        
        # Number of alleles per locus
        if (inherits(pop_genind, "genind")) {
          # Get number of alleles per locus
          alleles_per_locus <- pop_genind@loc.n.all
          
          if (length(alleles_per_locus) > 0) {
            results_df$N_Alleles_Mean[i] <- mean(alleles_per_locus, na.rm = TRUE)
            results_df$N_Alleles_SD[i] <- sd(alleles_per_locus, na.rm = TRUE)
            results_df$N_Alleles_Max[i] <- max(alleles_per_locus, na.rm = TRUE)
          }
        }
        
        # Calculate heterozygosity if possible
        tryCatch({
          # Convert to matrix
          geno_matrix <- as.matrix(pop_gl)
          
          # Calculate observed heterozygosity
          if (nrow(geno_matrix) > 0 && ncol(geno_matrix) > 0) {
            # For each locus, calculate proportion of heterozygotes
            ho_locus <- colMeans(geno_matrix == 1, na.rm = TRUE)
            results_df$Ho_Mean[i] <- mean(ho_locus, na.rm = TRUE)
            
            # Calculate expected heterozygosity (He)
            p <- colMeans(geno_matrix, na.rm = TRUE) / 2
            q <- 1 - p
            he_locus <- 2 * p * q
            results_df$He_Mean[i] <- mean(he_locus, na.rm = TRUE)
          }
        }, error = function(e) {
          cat("    Error calculating heterozygosity:", e$message, "\n")
        })
        
      }, error = function(e) {
        cat("    Error processing population", pop_name, ":", e$message, "\n")
      })
    }
    
    # Remove rows with NA results
    results_df <- results_df[results_df$N_Individuals >= min_individuals, ]
    
    if (nrow(results_df) == 0) {
      stop("No valid populations after filtering")
    }
    
    # Create default colors if not provided
    if (is.null(popcols)) {
      n_pops <- nrow(results_df)
      if (n_pops <= 8) {
        popcols <- RColorBrewer::brewer.pal(max(3, n_pops), "Set2")[1:n_pops]
      } else if (n_pops <= 12) {
        popcols <- RColorBrewer::brewer.pal(n_pops, "Set3")
      } else {
        popcols <- viridis::viridis(n_pops)
      }
    }
    
    # Create the main plot: N_Individuals vs N_Alleles_Mean
    main_plot <- ggplot(results_df, aes(x = N_Individuals, y = N_Alleles_Mean)) +
      geom_point(aes(color = Population), size = 5, alpha = 0.8) +
      geom_text(aes(label = Population), vjust = -0.8, hjust = 0.5, size = 3.5) +
      scale_color_manual(values = setNames(popcols, results_df$Population)) +
      labs(title = "Population Allele Statistics",
           subtitle = "Number of Individuals vs Mean Number of Alleles per Locus",
           x = "Number of Individuals",
           y = "Mean Number of Alleles per Locus") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "right",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray70", fill = NA, size = 0.5)
      )
    
    # Create additional summary plots
    summary_plots <- list()
    
    # 1. Bar plot of number of alleles per population
    summary_plots$alleles_bar <- ggplot(results_df, aes(x = reorder(Population, N_Alleles_Mean), 
                                                        y = N_Alleles_Mean)) +
      geom_bar(stat = "identity", aes(fill = Population), alpha = 0.7) +
      geom_errorbar(aes(ymin = N_Alleles_Mean - N_Alleles_SD, 
                        ymax = N_Alleles_Mean + N_Alleles_SD), width = 0.2) +
      scale_fill_manual(values = setNames(popcols, results_df$Population)) +
      labs(title = "Mean Number of Alleles per Locus by Population",
           x = "Population", y = "Mean Number of Alleles") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # 2. Heterozygosity plot
    if (all(c("He_Mean", "Ho_Mean") %in% colnames(results_df))) {
      heterozygosity_data <- results_df %>%
        select(Population, He_Mean, Ho_Mean) %>%
        pivot_longer(cols = c(He_Mean, Ho_Mean), 
                     names_to = "Type", 
                     values_to = "Value") %>%
        mutate(Type = ifelse(Type == "He_Mean", "Expected (He)", "Observed (Ho)"))
      
      summary_plots$heterozygosity <- ggplot(heterozygosity_data, 
                                             aes(x = Population, y = Value, fill = Type)) +
        geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
        scale_fill_brewer(palette = "Set2") +
        labs(title = "Heterozygosity by Population",
             x = "Population", y = "Heterozygosity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    # 3. Sample size distribution
    summary_plots$sample_size <- ggplot(results_df, aes(x = reorder(Population, N_Individuals), 
                                                        y = N_Individuals)) +
      geom_bar(stat = "identity", aes(fill = Population), alpha = 0.7) +
      scale_fill_manual(values = setNames(popcols, results_df$Population)) +
      labs(title = "Sample Size Distribution by Population",
           x = "Population", y = "Number of Individuals") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Prepare comprehensive results
    results <- list(
      statistics = results_df,
      colors = setNames(popcols, results_df$Population),
      plots = list(
        main = main_plot,
        summary = summary_plots
      ),
      summary_stats = list(
        total_populations = nrow(results_df),
        total_individuals = sum(results_df$N_Individuals),
        total_loci = mean(results_df$N_Loci, na.rm = TRUE),
        mean_alleles_per_pop = mean(results_df$N_Alleles_Mean, na.rm = TRUE),
        sd_alleles_per_pop = sd(results_df$N_Alleles_Mean, na.rm = TRUE),
        max_alleles = max(results_df$N_Alleles_Max, na.rm = TRUE)
      )
    )
    
    cat("Population allele analysis completed successfully\n")
    cat("  Populations analyzed:", nrow(results_df), "\n")
    cat("  Mean alleles per locus:", round(mean(results_df$N_Alleles_Mean, na.rm = TRUE), 2), "\n")
    
    return(results)
    
  }, error = function(e) {
    cat("Error in calculate_population_alleles:", e$message, "\n")
    return(NULL)
  })
}

#' Generate population summary report
#' @param allele_results Results from calculate_population_alleles
#' @return Formatted text report
generate_allele_summary_report <- function(allele_results) {
  if (is.null(allele_results)) return("No results available")
  
  report <- paste0(
    "POPULATION ALLELE STATISTICS REPORT\n",
    "===================================\n\n",
    "Summary Statistics:\n",
    "------------------\n",
    "Total populations analyzed: ", allele_results$summary_stats$total_populations, "\n",
    "Total individuals: ", allele_results$summary_stats$total_individuals, "\n",
    "Mean number of loci: ", round(allele_results$summary_stats$total_loci, 1), "\n",
    "Mean alleles per population: ", round(allele_results$summary_stats$mean_alleles_per_pop, 2), "\n",
    "Standard deviation: ", round(allele_results$summary_stats$sd_alleles_per_pop, 2), "\n",
    "Maximum alleles per locus: ", allele_results$summary_stats$max_alleles, "\n\n",
    
    "Population Details:\n",
    "------------------\n"
  )
  
  # Add details for each population
  for (i in 1:nrow(allele_results$statistics)) {
    pop <- allele_results$statistics[i, ]
    report <- paste0(report,
                     "Population: ", pop$Population, "\n",
                     "  Individuals: ", pop$N_Individuals, "\n",
                     "  Mean alleles/locus: ", round(pop$N_Alleles_Mean, 2), " (±", round(pop$N_Alleles_SD, 2), ")\n",
                     "  Max alleles/locus: ", pop$N_Alleles_Max, "\n"
    )
    
    if ("He_Mean" %in% colnames(pop) && "Ho_Mean" %in% colnames(pop)) {
      report <- paste0(report,
                       "  Expected heterozygosity (He): ", round(pop$He_Mean, 3), "\n",
                       "  Observed heterozygosity (Ho): ", round(pop$Ho_Mean, 3), "\n"
      )
    }
    report <- paste0(report, "\n")
  }
  
  return(report)
}



#5. Additional helper function (add to your functions section):
#' Quick population summary from genlight object
#' @param gl_object Genlight object
#' @return Summary string
quick_population_summary <- function(gl_object) {
  if (is.null(gl_object)) return("No genotype data")
  
  summary_text <- paste0(
    "Genotype Data Summary:\n",
    "  Individuals: ", nInd(gl_object), "\n",
    "  Loci: ", nLoc(gl_object), "\n"
  )
  
  if (!is.null(pop(gl_object))) {
    pop_table <- table(pop(gl_object))
    summary_text <- paste0(summary_text,
                           "  Populations: ", length(pop_table), "\n",
                           "  Population sizes: ", paste(paste0(names(pop_table), " (", pop_table, ")"), 
                                                         collapse = ", "), "\n"
    )
  } else {
    summary_text <- paste0(summary_text,
                           "  WARNING: No population assignments found\n"
    )
  }
  
  return(summary_text)
}


# ============================================================================
# LEA ANALYSIS FUNCTIONS
# ============================================================================

##' Convert genlight object to STRUCTURE format with proper formatting
#' @param gl_object Genlight object
#' @param file Output file path
#' @param missing_val Value for missing data (default: -9)
#' @return Path to STRUCTURE file
genlight2structure <- function(gl_object, file = "LEAinput.stru", missing_val = -9) {
  tryCatch({
    cat("Converting genlight to STRUCTURE format...\n")
    
    # Convert to matrix
    geno_matrix <- as.matrix(gl_object)
    n_ind <- nrow(geno_matrix)
    n_snp <- ncol(geno_matrix)
    
    cat("  Individuals:", n_ind, "\n")
    cat("  SNPs:", n_snp, "\n")
    
    # Get sample names and population info
    sample_names <- indNames(gl_object)
    if (!is.null(pop(gl_object))) {
      pop_info <- as.character(pop(gl_object))
    } else {
      pop_info <- rep("1", n_ind)
    }
    
    # Ensure population info is numeric for STRUCTURE
    pop_numeric <- as.numeric(as.factor(pop_info))
    
    # Create STRUCTURE matrix with 2 columns per SNP
    struct_mat <- matrix(missing_val, nrow = n_ind, ncol = 1 + 2 * n_snp)
    
    # First column: population ID (must be numeric)
    struct_mat[, 1] <- pop_numeric
    
    # Convert genotypes to STRUCTURE format (diploid, 0/1/2 coding)
    for (i in 1:n_ind) {
      for (j in 1:n_snp) {
        gt <- geno_matrix[i, j]
        col1 <- 2 * j  # First allele column
        col2 <- 2 * j + 1  # Second allele column
        
        if (is.na(gt)) {
          # Missing data
          struct_mat[i, col1] <- missing_val
          struct_mat[i, col2] <- missing_val
        } else if (gt == 0) {
          # Homozygous reference (0/0)
          struct_mat[i, col1] <- 1
          struct_mat[i, col2] <- 1
        } else if (gt == 1) {
          # Heterozygous (0/1 or 1/0)
          struct_mat[i, col1] <- 1
          struct_mat[i, col2] <- 2
        } else if (gt == 2) {
          # Homozygous alternate (1/1)
          struct_mat[i, col1] <- 2
          struct_mat[i, col2] <- 2
        } else {
          # Unknown genotype
          struct_mat[i, col1] <- missing_val
          struct_mat[i, col2] <- missing_val
        }
      }
    }
    
    # Write STRUCTURE file (no headers, no row names, tab-separated)
    write.table(struct_mat, file = file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE,
                na = as.character(missing_val))
    
    cat("  STRUCTURE file created:", file, "\n")
    cat("  Dimensions:", dim(struct_mat), "\n")
    
    # Verify file was written correctly
    if (file.exists(file)) {
      file_info <- readLines(file, n = 3)
      cat("  First 2 lines of file:\n")
      cat("  Line 1:", file_info[1], "\n")
      cat("  Line 2:", file_info[2], "\n")
    }
    
    return(file)
    
  }, error = function(e) {
    cat("Error creating STRUCTURE file:", e$message, "\n")
    return(NULL)
  })
}



#' Run LEA analysis with specified parameters
#' @param gl_object Genlight object
#' @param output_prefix Output file prefix
#' @param mindemes Minimum number of clusters (K)
#' @param maxdemes Maximum number of clusters (K)
#' @param n_replicates Number of replicates per K
#' @param alpha Regularization parameter
#' @param entropy Whether to calculate cross-entropy
#' @param exportdata Whether to export results
#' @return List containing LEA results
run_lea_analysis <- function(gl_object = NULL, output_prefix = "LEA_analysis",
                             mindemes = 2, maxdemes = 10, n_replicates = 10,
                             alpha = 10, entropy = TRUE, exportdata = TRUE) {
  
  tryCatch({
    cat("\n=== Running LEA Analysis ===\n")
    cat("  Input parameters:\n")
    cat("    Minimum K (mindemes):", mindemes, "\n")
    cat("    Maximum K (maxdemes):", maxdemes, "\n")
    cat("    Replicates per K:", n_replicates, "\n")
    cat("    Alpha (regularization):", alpha, "\n")
    cat("    Calculate entropy:", entropy, "\n")
    cat("    Export data:", exportdata, "\n")
    cat("    Output prefix:", output_prefix, "\n")
    
    # Validate input
    if (is.null(gl_object)) {
      stop("Genlight object is required")
    }
    
    if (mindemes < 2) {
      warning("Minimum K should be at least 2. Setting to 2.")
      mindemes <- 2
    }
    
    if (maxdemes <= mindemes) {
      warning(paste("Maximum K should be greater than minimum K.",
                    "Setting maxdemes to", mindemes + 3))
      maxdemes <- mindemes + 3
    }
    
    if (n_replicates < 1) {
      warning("Number of replicates should be at least 1. Setting to 1.")
      n_replicates <- 1
    }
    
    # Step 1: Convert genlight to LEA format
    cat("  Step 1: Converting genlight to LEA format...\n")
    conversion_result <- convert_to_lea_format_fixed(gl_object)
    
    if (is.null(conversion_result)) {
      stop("Failed to convert genotype data to LEA format")
    }
    
    # Create output directory if exporting
    if (exportdata) {
      output_dir <- paste0(output_prefix, "_results")
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      cat("  Output directory:", output_dir, "\n")
    }
    
    # Step 2: Run sNMF analysis
    cat("  Step 2: Running sNMF analysis...\n")
    K_range <- mindemes:maxdemes
    
    # Use the conversion result (lfmm object)
    lea_project <- tryCatch({
      LEA::snmf(
        input.file = conversion_result$lfmm_object,
        K = K_range,
        entropy = entropy,
        repetitions = n_replicates,
        project = "new",
        alpha = alpha,
        tolerance = 0.00001,
        iterations = 200,
        ploidy = 2
      )
    }, error = function(e) {
      cat("  Error in sNMF analysis:", e$message, "\n")
      
      # Try alternative method with file input
      cat("  Trying alternative method with file input...\n")
      LEA::snmf(
        input.file = conversion_result$file_path,
        K = K_range,
        entropy = entropy,
        repetitions = n_replicates,
        project = "new",
        alpha = alpha,
        tolerance = 0.00001,
        iterations = 200,
        ploidy = 2
      )
    })
    
    if (is.null(lea_project)) {
      stop("sNMF analysis failed")
    }
    
    cat("  sNMF analysis completed successfully\n")
    
    # Step 3: Extract comprehensive results
    cat("  Step 3: Extracting results...\n")
    results <- extract_lea_results_comprehensive(
      project = lea_project,
      K_range = K_range,
      n_replicates = n_replicates,
      sample_names = indNames(gl_object),
      metadata = list(
        run_date = Sys.time(),
        parameters = list(
          mindemes = mindemes,
          maxdemes = maxdemes,
          n_replicates = n_replicates,
          alpha = alpha,
          entropy = entropy
        ),
        conversion_info = conversion_result
      )
    )
    
    # Step 4: Export data if requested
    if (exportdata) {
      cat("  Step 4: Exporting results...\n")
      export_files <- export_lea_results_fixed(
        results = results,
        output_dir = output_dir,
        prefix = output_prefix
      )
      results$exported_files <- export_files
      cat("  Exported", length(export_files), "files to", output_dir, "\n")
    }
    
    # Step 5: Determine optimal K
    cat("  Step 5: Determining optimal K...\n")
    optimal_k <- find_optimal_k_fixed(results$cross_entropy_summary)
    results$optimal_k <- optimal_k
    
    cat("\n=== LEA Analysis Complete ===\n")
    cat("  K range analyzed:", mindemes, "to", maxdemes, "\n")
    cat("  Optimal K suggested:", optimal_k, "\n")
    cat("  Total samples:", nInd(gl_object), "\n")
    cat("  Total SNPs:", nLoc(gl_object), "\n")
    cat("  Memory usage:", format(object.size(results), units = "MB"), "\n")
    
    return(results)
    
  }, error = function(e) {
    cat("\n=== LEA Analysis Failed ===\n")
    cat("Error:", e$message, "\n")
    cat("Traceback:\n")
    print(traceback())
    
    # Return error structure
    return(list(
      success = FALSE,
      error = e$message,
      timestamp = Sys.time(),
      parameters = list(
        mindemes = mindemes,
        maxdemes = maxdemes,
        n_replicates = n_replicates,
        alpha = alpha,
        entropy = entropy,
        exportdata = exportdata
      )
    ))
  })
}


# ============================================================================
# HYBRID PATTERN VISUALIZATION
# ============================================================================

#' Visualize hybrid patterns from admixture analysis
visualize_hybrid_patterns <- function(q_matrix, sample_names = NULL,
                                      parent_references = NULL, 
                                      hybrid_status = NULL,
                                      threshold = 0.1) {
  tryCatch({
    cat("\n=== Visualizing Hybrid Patterns ===\n")
    
    # Set sample names if not provided
    if (is.null(sample_names)) {
      sample_names <- rownames(q_matrix)
      if (is.null(sample_names)) {
        sample_names <- paste0("Sample_", 1:nrow(q_matrix))
      }
    }
    
    n_clusters <- ncol(q_matrix)
    cat("  Samples:", nrow(q_matrix), "\n")
    cat("  Clusters:", n_clusters, "\n")
    
    # Create hybrid classification
    classification <- classify_hybrids(q_matrix, threshold)
    
    # Prepare visualization data
    viz_data <- prepare_hybrid_visualization_data(q_matrix, sample_names, 
                                                  classification, parent_references)
    
    # Generate visualizations
    plots <- create_hybrid_visualizations(viz_data, n_clusters)
    
    cat("  Hybrid visualization complete\n")
    
    return(list(
      classification = classification,
      visualization_data = viz_data,
      plots = plots,
      threshold = threshold
    ))
    
  }, error = function(e) {
    cat("Error in hybrid visualization:", e$message, "\n")
    return(NULL)
  })
}

#' Classify individuals based on admixture proportions
classify_hybrids <- function(q_matrix, threshold = 0.1) {
  classification <- data.frame(
    Sample = rownames(q_matrix),
    Dominant_Cluster = apply(q_matrix, 1, which.max),
    Max_Proportion = apply(q_matrix, 1, max),
    Category = NA,
    stringsAsFactors = FALSE
  )
  
  # Determine category based on admixture proportions
  for (i in 1:nrow(classification)) {
    if (classification$Max_Proportion[i] > (1 - threshold)) {
      classification$Category[i] <- "Pure"
    } else if (classification$Max_Proportion[i] > 0.5) {
      classification$Category[i] <- "Hybrid"
    } else {
      classification$Category[i] <- "Admixed"
    }
  }
  
  # Calculate admixture proportions
  classification$Admixture_Proportion <- 1 - classification$Max_Proportion
  
  return(classification)
}

#' Prepare data for hybrid visualization
prepare_hybrid_visualization_data <- function(q_matrix, sample_names, 
                                              classification, parent_references) {
  viz_data <- list()
  
  # 1. Structure plot data
  viz_data$structure <- data.frame(
    Sample = rep(sample_names, each = ncol(q_matrix)),
    Cluster = rep(paste0("Cluster_", 1:ncol(q_matrix)), times = nrow(q_matrix)),
    Proportion = as.vector(t(q_matrix)),
    stringsAsFactors = FALSE
  )
  
  viz_data$structure <- merge(viz_data$structure, 
                              classification[, c("Sample", "Category")], 
                              by = "Sample")
  
  # 2. Triangle plot data (for 3 clusters)
  if (ncol(q_matrix) == 3) {
    viz_data$triangle <- data.frame(
      Sample = sample_names,
      Cluster1 = q_matrix[, 1],
      Cluster2 = q_matrix[, 2],
      Cluster3 = q_matrix[, 3],
      Category = classification$Category,
      stringsAsFactors = FALSE
    )
  }
  
  return(viz_data)
}

#' Create hybrid visualizations
create_hybrid_visualizations <- function(viz_data, n_clusters) {
  plots <- list()
  
  # 1. Structure plot with hybrid classification
  plots$structure <- ggplot(viz_data$structure, 
                            aes(x = Sample, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(Category ~ ., scales = "free_x", space = "free") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Hybrid Structure Plot",
         x = "Samples",
         y = "Admixture Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # 2. Triangle plot (for 3 clusters)
  if (!is.null(viz_data$triangle)) {
    plots$triangle <- ggplot(viz_data$triangle, 
                             aes(x = Cluster1, y = Cluster2, color = Category)) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "Ternary Plot of Admixture Proportions",
           x = "Cluster 1 Proportion",
           y = "Cluster 2 Proportion") +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
  
  return(plots)
}

#' Alternative conversion using dartR
convert_to_lea_format_dartr <- function(gl_object) {
  tryCatch({
    cat("  Using dartR for LEA conversion...\n")
    
    if (!requireNamespace("dartR", quietly = TRUE)) {
      stop("dartR package not available")
    }
    
    # Convert to genind
    genind_obj <- dartR::gl2gi(gl_object, verbose = 0)
    
    # Convert to matrix
    geno_matrix <- as.matrix(genind_obj)
    
    # Convert to 0,1,2 coding
    # For biallelic markers, we can use:
    # 0 = homozygous reference
    # 1 = heterozygous
    # 2 = homozygous alternative
    
    # Get allele frequencies
    freq <- adegenet::makefreq(genind_obj)
    
    # Determine reference allele (most common)
    ref_allele <- apply(freq, 2, function(x) names(which.max(x)))
    
    # Convert genotype matrix
    geno_numeric <- matrix(0, nrow = nrow(geno_matrix), ncol = ncol(geno_matrix))
    
    for (i in 1:ncol(geno_matrix)) {
      genotypes <- geno_matrix[, i]
      ref <- ref_allele[i]
      
      # Convert to 0,1,2
      geno_numeric[genotypes == paste(ref, ref, sep = "/")] <- 0
      geno_numeric[!grepl(ref, genotypes) & genotypes != "NA/NA"] <- 2
      # Everything else (heterozygous or different format) becomes 1
      geno_numeric[geno_numeric == 0 & genotypes != "NA/NA"] <- 1
      geno_numeric[genotypes == "NA/NA"] <- 9
    }
    
    # Transpose for LEA
    geno_t <- t(geno_numeric)
    
    # Write to file
    output_file <- tempfile(fileext = ".geno")
    LEA::write.geno(geno_t, output_file)
    
    return(list(
      file_path = output_file,
      lfmm_object = LEA::as.lfmm(geno_t),
      n_individuals = nrow(geno_numeric),
      n_snps = ncol(geno_numeric)
    ))
    
  }, error = function(e) {
    cat("  dartR conversion failed:", e$message, "\n")
    return(NULL)
  })
}



#' Extract comprehensive results from LEA project
#' @param project LEA sNMF project
#' @param K_range Range of K values
#' @param n_replicates Number of replicates
#' @param sample_names Vector of sample names
#' @param metadata Additional metadata
#' @return List of comprehensive results
extract_lea_results_comprehensive <- function(project, K_range, n_replicates, 
                                              sample_names = NULL, metadata = NULL) {
  tryCatch({
    cat("  Extracting results for K =", paste(K_range, collapse = ", "), "\n")
    
    results <- list(
      project = project,
      K_range = K_range,
      n_replicates = n_replicates,
      Q_matrices = list(),
      cross_entropy = list(),
      best_runs = list(),
      metadata = metadata
    )
    
    # Extract for each K
    for (K in K_range) {
      k_str <- paste0("K", K)
      cat("    Processing K =", K, "\n")
      
      # Get cross-entropy values
      ce_vals <- tryCatch({
        as.numeric(LEA::cross.entropy(project, K = K))
      }, error = function(e) {
        cat("      Warning: Could not get cross-entropy for K =", K, "\n")
        rep(NA, n_replicates)
      })
      
      # Find best run (lowest cross-entropy)
      if (all(is.na(ce_vals))) {
        best_run <- 1
        cat("      Using run 1 (cross-entropy not available)\n")
      } else {
        best_run <- which.min(ce_vals)
        cat("      Best run:", best_run, 
            "(cross-entropy =", round(min(ce_vals, na.rm = TRUE), 4), ")\n")
      }
      
      # Extract Q matrix for best run
      q_matrix <- tryCatch({
        LEA::Q(project, K = K, run = best_run)
      }, error = function(e) {
        cat("      Error extracting Q matrix:", e$message, "\n")
        
        # Try to get any run
        if (LEA::runs(project, K = K) > 0) {
          LEA::Q(project, K = K, run = 1)
        } else {
          NULL
        }
      })
      
      # Store results
      results$cross_entropy[[k_str]] <- ce_vals
      results$best_runs[[k_str]] <- best_run
      
      if (!is.null(q_matrix)) {
        # Set sample names if provided
        if (!is.null(sample_names) && length(sample_names) == nrow(q_matrix)) {
          rownames(q_matrix) <- sample_names
        } else {
          rownames(q_matrix) <- paste0("Sample_", 1:nrow(q_matrix))
        }
        
        # Set cluster names
        colnames(q_matrix) <- paste0("Cluster_", 1:K)
        
        # Store Q matrix
        results$Q_matrices[[k_str]] <- list(
          Average_Q = q_matrix,
          Best_Run = best_run,
          CrossEntropy = ifelse(all(is.na(ce_vals)), NA, min(ce_vals, na.rm = TRUE)),
          Dimensions = dim(q_matrix)
        )
        
        # Calculate average across all runs (if multiple replicates)
        if (n_replicates > 1) {
          all_q_matrices <- list()
          valid_runs <- 0
          
          for (r in 1:n_replicates) {
            q_run <- tryCatch({
              LEA::Q(project, K = K, run = r)
            }, error = function(e) NULL)
            
            if (!is.null(q_run)) {
              all_q_matrices[[r]] <- q_run
              valid_runs <- valid_runs + 1
            }
          }
          
          if (valid_runs > 1) {
            # Calculate average Q matrix
            avg_q <- Reduce("+", all_q_matrices) / valid_runs
            rownames(avg_q) <- rownames(q_matrix)
            colnames(avg_q) <- colnames(q_matrix)
            
            results$Q_matrices[[k_str]]$Average_Q_all_runs <- avg_q
            results$Q_matrices[[k_str]]$N_Valid_Runs <- valid_runs
          }
        }
      }
    }
    
    # Create cross-entropy summary
    results$cross_entropy_summary <- data.frame(
      K = K_range,
      Mean_CE = sapply(K_range, function(K) {
        ce <- results$cross_entropy[[paste0("K", K)]]
        if (all(is.na(ce))) NA else mean(ce, na.rm = TRUE)
      }),
      Min_CE = sapply(K_range, function(K) {
        ce <- results$cross_entropy[[paste0("K", K)]]
        if (all(is.na(ce))) NA else min(ce, na.rm = TRUE)
      }),
      SD_CE = sapply(K_range, function(K) {
        ce <- results$cross_entropy[[paste0("K", K)]]
        if (all(is.na(ce))) NA else sd(ce, na.rm = TRUE)
      })
    )
    
    cat("  Results extraction complete\n")
    return(results)
    
  }, error = function(e) {
    cat("  Error extracting results:", e$message, "\n")
    return(NULL)
  })
}

#' Export LEA results (fixed)
#' @param results LEA results
#' @param output_dir Output directory
#' @param prefix File prefix
#' @return List of exported file paths
export_lea_results_fixed <- function(results, output_dir, prefix) {
  tryCatch({
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    exported_files <- list()
    
    # 1. Export Q matrices
    for (k_str in names(results$Q_matrices)) {
      q_data <- results$Q_matrices[[k_str]]
      
      if (!is.null(q_data$Average_Q)) {
        # CSV format
        csv_file <- file.path(output_dir, paste0(prefix, "_", k_str, "_qmatrix.csv"))
        write.csv(q_data$Average_Q, csv_file, row.names = TRUE)
        exported_files[[paste0("Q_", k_str, "_csv")]] <- csv_file
        
        # RDS format (for R)
        rds_file <- file.path(output_dir, paste0(prefix, "_", k_str, "_qmatrix.rds"))
        saveRDS(q_data$Average_Q, rds_file)
        exported_files[[paste0("Q_", k_str, "_rds")]] <- rds_file
      }
    }
    
    # 2. Export cross-entropy summary
    if (!is.null(results$cross_entropy_summary)) {
      ce_file <- file.path(output_dir, paste0(prefix, "_cross_entropy.csv"))
      write.csv(results$cross_entropy_summary, ce_file, row.names = FALSE)
      exported_files$cross_entropy <- ce_file
    }
    
    # 3. Export metadata
    meta_file <- file.path(output_dir, paste0(prefix, "_metadata.txt"))
    sink(meta_file)
    cat("LEA Analysis Metadata\n")
    cat("=====================\n\n")
    cat("Analysis date:", as.character(results$metadata$run_date), "\n")
    cat("K range:", paste(range(results$K_range), collapse = "-"), "\n")
    cat("Replicates per K:", results$n_replicates, "\n")
    
    if (!is.null(results$metadata$parameters)) {
      cat("\nParameters:\n")
      for (param in names(results$metadata$parameters)) {
        cat("  ", param, ": ", results$metadata$parameters[[param]], "\n", sep = "")
      }
    }
    
    if (!is.null(results$optimal_k)) {
      cat("\nOptimal K (suggested):", results$optimal_k, "\n")
    }
    
    cat("\nAvailable Q matrices:", 
        paste(names(results$Q_matrices), collapse = ", "), "\n")
    sink()
    exported_files$metadata <- meta_file
    
    # 4. Export structure plots
    plot_dir <- file.path(output_dir, "plots")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    
    # Export structure plot for each K
    for (k_str in names(results$Q_matrices)) {
      q_matrix <- results$Q_matrices[[k_str]]$Average_Q
      if (!is.null(q_matrix)) {
        plot_file <- file.path(plot_dir, paste0(prefix, "_", k_str, "_structure.pdf"))
        tryCatch({
          pdf(plot_file, width = 12, height = 8)
          plot_structure_fixed(q_matrix, main = paste("Structure Plot -", k_str))
          dev.off()
          exported_files[[paste0("plot_", k_str)]] <- plot_file
        }, error = function(e) {
          cat("  Could not create plot for", k_str, ":", e$message, "\n")
        })
      }
    }
    
    # 5. Export cross-entropy plot
    if (!is.null(results$cross_entropy_summary)) {
      ce_plot_file <- file.path(plot_dir, paste0(prefix, "_cross_entropy.pdf"))
      tryCatch({
        pdf(ce_plot_file, width = 10, height = 6)
        plot_cross_entropy_fixed(results$cross_entropy_summary, 
                                 optimal_k = results$optimal_k)
        dev.off()
        exported_files$cross_entropy_plot <- ce_plot_file
      }, error = function(e) {
        cat("  Could not create cross-entropy plot:", e$message, "\n")
      })
    }
    
    return(exported_files)
    
  }, error = function(e) {
    cat("Error exporting results:", e$message, "\n")
    return(NULL)
  })
}



#' Find optimal K (fixed)
#' @param cross_entropy_summary Cross-entropy summary data frame
#' @return Optimal K value
find_optimal_k_fixed <- function(cross_entropy_summary) {
  tryCatch({
    if (is.null(cross_entropy_summary) || nrow(cross_entropy_summary) < 3) {
      return(NULL)
    }
    
    ce_df <- cross_entropy_summary[!is.na(cross_entropy_summary$Mean_CE), ]
    
    if (nrow(ce_df) < 3) {
      return(NULL)
    }
    
    # Method 1: Minimum cross-entropy
    min_ce_k <- ce_df$K[which.min(ce_df$Mean_CE)]
    
    # Method 2: Elbow method
    # Calculate second derivative
    ce_df$FirstDiff <- c(NA, diff(ce_df$Mean_CE))
    ce_df$SecondDiff <- c(NA, diff(ce_df$FirstDiff))
    
    # Find elbow (point of maximum curvature)
    if (sum(!is.na(ce_df$SecondDiff)) >= 2) {
      # Use change in slope method
      slopes <- abs(ce_df$FirstDiff[-1])
      slope_changes <- c(NA, diff(slopes))
      
      if (any(!is.na(slope_changes))) {
        # Find where slope change is maximum (elbow)
        max_change_idx <- which.max(abs(slope_changes))
        elbow_k <- ce_df$K[max_change_idx + 1]
      } else {
        elbow_k <- min_ce_k
      }
    } else {
      elbow_k <- min_ce_k
    }
    
    # Choose between methods (prefer elbow method if available)
    optimal_k <- ifelse(!is.na(elbow_k) && elbow_k >= min(ce_df$K), elbow_k, min_ce_k)
    
    return(optimal_k)
    
  }, error = function(e) {
    cat("Error finding optimal K:", e$message, "\n")
    return(NULL)
  })
}

#' Plot structure (fixed)
#' @param q_matrix Q matrix
#' @param main Plot title
#' @param colors Color palette
plot_structure_fixed <- function(q_matrix, main = "Structure Plot", 
                                 colors = NULL) {
  tryCatch({
    n_individuals <- nrow(q_matrix)
    n_clusters <- ncol(q_matrix)
    
    # Set up colors
    if (is.null(colors)) {
      if (n_clusters <= 8) {
        colors <- RColorBrewer::brewer.pal(n_clusters, "Set2")
      } else if (n_clusters <= 12) {
        colors <- RColorBrewer::brewer.pal(n_clusters, "Set3")
      } else {
        colors <- viridis::viridis(n_clusters)
      }
    }
    
    # Prepare data for barplot
    par(mar = c(5, 4, 4, 8), xpd = TRUE)
    
    # Create barplot
    barplot(t(q_matrix), 
            col = colors,
            border = NA,
            space = 0,
            xlab = "Individuals",
            ylab = "Ancestry Proportion",
            main = main,
            las = 2,
            cex.names = 0.7)
    
    # Add legend
    legend("topright",
           inset = c(-0.15, 0),
           legend = paste("Cluster", 1:n_clusters),
           fill = colors,
           border = NA,
           bty = "n",
           cex = 0.8)
    
    # Add grid
    grid(nx = NA, ny = NULL, col = "gray", lty = "dotted")
    
  }, error = function(e) {
    cat("Error in structure plot:", e$message, "\n")
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Plot Error")
    text(1, 1, "Could not create structure plot", col = "red")
  })
}

#' Plot cross-entropy (fixed)
#' @param cross_entropy_summary Cross-entropy data
#' @param optimal_k Optimal K to highlight
plot_cross_entropy_fixed <- function(cross_entropy_summary, optimal_k = NULL) {
  tryCatch({
    ce_df <- cross_entropy_summary[!is.na(cross_entropy_summary$Mean_CE), ]
    
    if (nrow(ce_df) < 2) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
           main = "Insufficient cross-entropy data")
      return()
    }
    
    # Create plot
    plot(ce_df$K, ce_df$Mean_CE,
         type = "o",
         pch = 19,
         col = "steelblue",
         lwd = 2,
         xlab = "Number of Clusters (K)",
         ylab = "Cross-Entropy",
         main = "Cross-Entropy for Model Selection",
         ylim = range(c(ce_df$Mean_CE - ce_df$SD_CE, 
                        ce_df$Mean_CE + ce_df$SD_CE), na.rm = TRUE))
    
    # Add error bars
    arrows(ce_df$K, ce_df$Mean_CE - ce_df$SD_CE,
           ce_df$K, ce_df$Mean_CE + ce_df$SD_CE,
           length = 0.05, angle = 90, code = 3, col = "gray50")
    
    # Highlight optimal K
    if (!is.null(optimal_k) && optimal_k %in% ce_df$K) {
      points(optimal_k, ce_df$Mean_CE[ce_df$K == optimal_k],
             pch = 21, col = "red", bg = "red", cex = 2, lwd = 2)
      text(optimal_k, ce_df$Mean_CE[ce_df$K == optimal_k],
           labels = paste("K =", optimal_k),
           pos = 3, col = "red", font = 2)
    }
    
    # Add grid
    grid(col = "gray", lty = "dotted")
    
    # Add legend
    legend("topright",
           legend = c("Mean cross-entropy", "± 1 SD"),
           col = c("steelblue", "gray50"),
           pch = c(19, NA),
           lty = c(1, 1),
           lwd = c(2, 1),
           bty = "n")
    
  }, error = function(e) {
    cat("Error in cross-entropy plot:", e$message, "\n")
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Plot Error")
    text(1, 1, "Could not create cross-entropy plot", col = "red")
  })
}



#' Create STRUCTURE format file with proper formatting
create_structure_file_fixed <- function(geno_matrix, pop_info, output_file, missing_val = -9) {
  tryCatch({
    n_ind <- nrow(geno_matrix)
    n_snp <- ncol(geno_matrix)
    
    # Open file for writing
    con <- file(output_file, "w")
    
    # Write each individual
    for (i in 1:n_ind) {
      # Start with individual ID and population
      line <- c(i, pop_info[i])
      
      # Add genotypes (2 columns per SNP)
      for (j in 1:n_snp) {
        gt <- geno_matrix[i, j]
        
        if (is.na(gt) || gt == -9) {
          # Missing data
          line <- c(line, missing_val, missing_val)
        } else if (gt == 0) {
          # Homozygous reference (0/0) -> 1 1
          line <- c(line, 1, 1)
        } else if (gt == 1) {
          # Heterozygous (0/1) -> 1 2
          line <- c(line, 1, 2)
        } else if (gt == 2) {
          # Homozygous alternate (1/1) -> 2 2
          line <- c(line, 2, 2)
        } else {
          # Unknown
          line <- c(line, missing_val, missing_val)
        }
      }
      
      # Write line
      writeLines(paste(line, collapse = "\t"), con)
    }
    
    close(con)
    
    # Verify file was created
    if (file.exists(output_file)) {
      cat("  Created STRUCTURE file with", n_ind, "individuals and", n_snp, "SNPs\n")
      return(TRUE)
    } else {
      cat("  ERROR: STRUCTURE file not created\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat("  Error creating STRUCTURE file:", e$message, "\n")
    return(FALSE)
  })
}

#' Run sNMF with proper parameter handling
run_snmf_with_parameters <- function(geno_file, K_min, K_max, replicates, 
                                     alpha, entropy, project_name) {
  tryCatch({
    cat("  Running sNMF with LEA package...\n")
    cat("  Input file:", geno_file, "\n")
    cat("  File exists:", file.exists(geno_file), "\n")
    
    if (!file.exists(geno_file)) {
      stop("Input file not found: ", geno_file)
    }
    
    # Create sNMF project
    snmf_project <- LEA::snmf(
      input.file = geno_file,
      K = K_min:K_max,
      repetitions = replicates,
      project = "new",
      entropy = entropy,
      alpha = alpha,
      tolerance = 0.00001,
      iterations = 200,
      CPU = 1,
      seed = 12345  # For reproducibility
    )
    
    # Save project
    project_file <- paste0(project_name, ".snmfProject")
    LEA::snmfProject.save(snmf_project, project_file)
    
    # Extract results
    results <- list(
      success = TRUE,
      project = snmf_project,
      project_file = project_file
    )
    
    # Get cross-entropy
    if (entropy) {
      tryCatch({
        cross_entropy <- LEA::cross.entropy(snmf_project, K = K_min:K_max)
        results$cross_entropy <- cross_entropy
        results$best_K <- which.min(colMeans(cross_entropy)) + K_min - 1
        cat("  Best K (lowest cross-entropy):", results$best_K, "\n")
      }, error = function(e) {
        cat("  Could not calculate cross-entropy:", e$message, "\n")
      })
    }
    
    # Extract Q matrices for each K
    Q_matrices <- list()
    for (K in K_min:K_max) {
      tryCatch({
        # Get the best run (lowest cross-entropy)
        if (entropy && !is.null(results$cross_entropy)) {
          best_run <- which.min(LEA::cross.entropy(snmf_project, K = K))
        } else {
          best_run <- 1
        }
        
        Q <- LEA::Q(snmf_project, K = K, run = best_run)
        Q_matrices[[paste0("K", K)]] <- Q
      }, error = function(e) {
        cat("  Could not extract Q matrix for K=", K, ":", e$message, "\n")
      })
    }
    
    if (length(Q_matrices) > 0) {
      results$Q_matrices <- Q_matrices
      cat("  Extracted Q matrices for K =", paste(names(Q_matrices), collapse = ", "), "\n")
    } else {
      results$error <- "No Q matrices extracted"
    }
    
    return(results)
    
  }, error = function(e) {
    cat("  Error in sNMF analysis:", e$message, "\n")
    return(list(
      success = FALSE,
      error = e$message
    ))
  })
}

#' Fallback ADMIXTURE-style analysis
run_admixture_fallback <- function(geno_matrix, K_range, replicates, output_prefix) {
  tryCatch({
    cat("  Running ADMIXTURE-style analysis...\n")
    
    # Clean data
    geno_clean <- geno_matrix
    geno_clean[is.na(geno_clean)] <- 0
    geno_clean[geno_clean == -9] <- 0
    
    # Perform PCA
    pca <- prcomp(geno_clean, center = TRUE, scale. = TRUE)
    
    # Create Q matrices using k-means on PCA scores
    Q_matrices <- list()
    cross_entropy <- matrix(NA, nrow = replicates, ncol = length(K_range))
    colnames(cross_entropy) <- paste0("K", K_range)
    
    for (k_idx in seq_along(K_range)) {
      K <- K_range[k_idx]
      cat("    Processing K =", K, "\n")
      
      # Multiple replicates
      Q_replicates <- list()
      ce_replicates <- numeric(replicates)
      
      for (rep in 1:replicates) {
        # K-means clustering
        kmeans_result <- kmeans(pca$x[, 1:min(10, ncol(pca$x))], 
                                centers = K, nstart = 25)
        
        # Create Q matrix from cluster assignments
        Q <- matrix(0, nrow = nrow(geno_matrix), ncol = K)
        for (i in 1:nrow(geno_matrix)) {
          Q[i, kmeans_result$cluster[i]] <- 1
        }
        
        # Add small random noise and normalize
        Q <- Q + matrix(runif(nrow(Q) * ncol(Q), 0, 0.01), 
                        nrow = nrow(Q), ncol = ncol(Q))
        Q <- Q / rowSums(Q)
        
        Q_replicates[[rep]] <- Q
      }
      
      # Average replicates
      Q_avg <- Reduce("+", Q_replicates) / length(Q_replicates)
      Q_matrices[[paste0("K", K)]] <- Q_avg
    }
    
    return(list(
      success = TRUE,
      Q_matrices = Q_matrices,
      method = "PCA-based ADMIXTURE (fallback)",
      cross_entropy = cross_entropy,
      best_K = K_range[which.min(colMeans(cross_entropy, na.rm = TRUE))]
    ))
    
  }, error = function(e) {
    cat("  Error in ADMIXTURE fallback:", e$message, "\n")
    return(list(
      success = FALSE,
      error = e$message
    ))
  })
}

#' Create LEA plots
create_lea_plots <- function(snmf_results, sample_names, pop_info, output_prefix) {
  tryCatch({
    plots <- list()
    
    # 1. Cross-entropy plot (if available)
    if (!is.null(snmf_results$cross_entropy)) {
      ce_data <- as.data.frame(snmf_results$cross_entropy)
      ce_data$K <- as.numeric(gsub("K", "", rownames(ce_data)))
      
      ce_plot <- ggplot2::ggplot(ce_data, ggplot2::aes(x = K)) +
        ggplot2::geom_boxplot(ggplot2::aes(group = K, y = rowMeans(snmf_results$cross_entropy, na.rm = TRUE))) +
        ggplot2::labs(title = "Cross-Entropy by K",
                      x = "Number of Clusters (K)",
                      y = "Cross-Entropy") +
        ggplot2::theme_minimal()
      
      plots$cross_entropy <- ce_plot
      
      # Save plot
      ggplot2::ggsave(paste0(output_prefix, "_cross_entropy.png"), 
                      ce_plot, width = 8, height = 6, dpi = 300)
    }
    
    # 2. Structure plots for each K
    if (!is.null(snmf_results$Q_matrices)) {
      for (k_name in names(snmf_results$Q_matrices)) {
        Q <- snmf_results$Q_matrices[[k_name]]
        
        # Prepare data for structure plot
        plot_data <- as.data.frame(Q)
        colnames(plot_data) <- paste0("Cluster", 1:ncol(Q))
        plot_data$Sample <- sample_names
        plot_data$Population <- pop_info
        
        # Melt for ggplot
        plot_data_long <- reshape2::melt(plot_data, 
                                         id.vars = c("Sample", "Population"),
                                         variable.name = "Cluster",
                                         value.name = "Proportion")
        
        # Order samples by population and cluster assignment
        plot_data_long$Sample <- factor(plot_data_long$Sample,
                                        levels = sample_names[order(pop_info)])
        
        # Create structure plot
        struct_plot <- ggplot2::ggplot(plot_data_long, 
                                       ggplot2::aes(x = Sample, y = Proportion, fill = Cluster)) +
          ggplot2::geom_bar(stat = "identity", width = 1) +
          ggplot2::scale_fill_viridis_d() +
          ggplot2::labs(title = paste("Population Structure (", k_name, ")", sep = ""),
                        x = "Individuals",
                        y = "Ancestry Proportion") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
            axis.ticks.x = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            legend.position = "bottom"
          )
        
        plots[[paste0("structure_", k_name)]] <- struct_plot
        
        # Save plot
        ggplot2::ggsave(paste0(output_prefix, "_structure_", k_name, ".png"),
                        struct_plot, width = 12, height = 6, dpi = 300)
      }
    }
    
    return(plots)
    
  }, error = function(e) {
    cat("Error creating LEA plots:", e$message, "\n")
    return(NULL)
  })
}

#' Export LEA data
export_lea_data <- function(snmf_results, sample_names, output_prefix) {
  tryCatch({
    exported_files <- list()
    
    # Export Q matrices
    if (!is.null(snmf_results$Q_matrices)) {
      for (k_name in names(snmf_results$Q_matrices)) {
        Q <- snmf_results$Q_matrices[[k_name]]
        q_df <- as.data.frame(Q)
        colnames(q_df) <- paste0("Cluster", 1:ncol(Q))
        q_df$Sample <- sample_names
        
        q_file <- paste0(output_prefix, "_Q_", k_name, ".csv")
        write.csv(q_df, q_file, row.names = FALSE)
        exported_files[[k_name]] <- q_file
      }
    }
    
    # Export cross-entropy
    if (!is.null(snmf_results$cross_entropy)) {
      ce_file <- paste0(output_prefix, "_cross_entropy.csv")
      write.csv(snmf_results$cross_entropy, ce_file, row.names = TRUE)
      exported_files$cross_entropy <- ce_file
    }
    
    # Export summary
    summary_file <- paste0(output_prefix, "_summary.txt")
    sink(summary_file)
    cat("LEA Analysis Summary\n")
    cat("====================\n\n")
    cat("Analysis date:", Sys.time(), "\n")
    cat("Number of Q matrices:", length(snmf_results$Q_matrices), "\n")
    cat("Available K values:", paste(names(snmf_results$Q_matrices), collapse = ", "), "\n")
    if (!is.null(snmf_results$best_K)) {
      cat("Best K:", snmf_results$best_K, "\n")
    }
    cat("Analysis method:", ifelse(!is.null(snmf_results$method), 
                                   snmf_results$method, "LEA sNMF"), "\n")
    sink()
    
    exported_files$summary <- summary_file
    
    return(exported_files)
    
  }, error = function(e) {
    cat("Error exporting LEA data:", e$message, "\n")
    return(NULL)
  })
}

#' Create LEA summary
create_lea_summary <- function(snmf_results, parameters) {
  summary <- list(
    analysis_completed = TRUE,
    timestamp = Sys.time(),
    parameters = parameters,
    n_q_matrices = length(snmf_results$Q_matrices),
    available_K = names(snmf_results$Q_matrices),
    best_K = snmf_results$best_K,
    method = snmf_results$method
  )
  
  if (!is.null(snmf_results$cross_entropy)) {
    summary$mean_cross_entropy <- colMeans(snmf_results$cross_entropy, na.rm = TRUE)
  }
  
  return(summary)
}




#' Prepare data for LEA analysis with robust error handling
prepare_lea_data <- function(geno_matrix, pop_info, output_prefix) {
  tryCatch({
    cat("  Preparing data for LEA analysis...\n")
    
    n_ind <- nrow(geno_matrix)
    n_snp <- ncol(geno_matrix)
    
    # Create temporary STRUCTURE file
    struct_file <- paste0(output_prefix, ".stru")
    
    # Convert population info to numeric
    if (!is.numeric(pop_info)) {
      pop_info <- as.numeric(as.factor(pop_info))
    }
    
    # Ensure we have valid population IDs (1-based)
    pop_info[is.na(pop_info)] <- 1
    if (any(pop_info == 0)) {
      pop_info <- pop_info + 1
    }
    
    # Create STRUCTURE matrix
    struct_mat <- matrix(-9, nrow = n_ind, ncol = 1 + 2 * n_snp)
    struct_mat[, 1] <- pop_info
    
    # Convert genotypes with consistent missing data coding
    for (i in 1:n_ind) {
      for (j in 1:n_snp) {
        gt <- geno_matrix[i, j]
        
        # Skip if genotype is -9 (already missing)
        if (!is.na(gt) && gt == -9) {
          next  # Keep as -9
        }
        
        col1 <- 2 * j
        col2 <- 2 * j + 1
        
        if (is.na(gt)) {
          struct_mat[i, col1] <- -9
          struct_mat[i, col2] <- -9
        } else if (gt == 0) {
          struct_mat[i, col1] <- 1
          struct_mat[i, col2] <- 1
        } else if (gt == 1) {
          struct_mat[i, col1] <- 1
          struct_mat[i, col2] <- 2
        } else if (gt == 2) {
          struct_mat[i, col1] <- 2
          struct_mat[i, col2] <- 2
        }
      }
    }
    
    # Write STRUCTURE file
    write.table(struct_mat, file = struct_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    cat("  STRUCTURE file created with dimensions:", dim(struct_mat), "\n")
    
    # Try to convert to LEA format
    geno_file <- NULL
    if (requireNamespace("LEA", quietly = TRUE)) {
      tryCatch({
        cat("  Converting to LEA format...\n")
        # Use struct2geno with proper parameters
        geno_file <- LEA::struct2geno(struct_file, ploidy = 2, FORMAT = 1, 
                                      extra.row = 0, extra.col = 0)
        cat("  LEA format file created:", geno_file, "\n")
      }, error = function(e) {
        cat("  LEA conversion failed:", e$message, "\n")
        cat("  Using STRUCTURE file directly\n")
        geno_file <- struct_file
      })
    } else {
      cat("  LEA package not available, using STRUCTURE file\n")
      geno_file <- struct_file
    }
    
    return(list(
      structure_file = struct_file,
      geno_file = geno_file
    ))
    
  }, error = function(e) {
    cat("Error in prepare_lea_data:", e$message, "\n")
    return(NULL)
  })
}


#' Convert genlight object to LEA format (corrected)
#' @param gl_object Genlight object
#' @param output_file Output file path (optional)
#' @return Path to the created geno file
convert_genlight_to_lea <- function(gl_object, output_file = NULL) {
  tryCatch({
    cat("Converting genlight to LEA format...\n")
    
    # Get basic information
    n_ind <- nInd(gl_object)
    n_loc <- nLoc(gl_object)
    
    cat("  Individuals:", n_ind, "\n")
    cat("  Loci:", n_loc, "\n")
    
    # Convert to matrix (individuals x loci)
    geno_matrix <- as.matrix(gl_object)
    
    # Replace NA with 9 (missing data code for LEA)
    na_count <- sum(is.na(geno_matrix))
    cat("  Replacing", na_count, "NA values with 9...\n")
    geno_matrix[is.na(geno_matrix)] <- 9
    
    # Transpose matrix for LEA format (loci x individuals)
    geno_t <- t(geno_matrix)
    cat("  Transposing matrix for LEA format...\n")
    
    # Create output file path if not provided
    if (is.null(output_file)) {
      output_file <- tempfile(fileext = ".geno")
    }
    
    cat("  Using temporary file:", output_file, "\n")
    
    # Write geno file using base R (more reliable than LEA::write.geno)
    cat("  Writing to file...\n")
    
    # Open connection
    con <- file(output_file, "w")
    
    # Write each row (locus)
    for (i in 1:nrow(geno_t)) {
      writeLines(paste(geno_t[i, ], collapse = ""), con)
    }
    
    close(con)
    
    cat("  File created:", output_file, "(", file.size(output_file), "bytes)\n")
    
    return(output_file)
    
  }, error = function(e) {
    cat("Error converting to LEA format:", e$message, "\n")
    return(NULL)
  })
}


#' Run sNMF analysis for hybrid detection (corrected)
#' @param geno_file Path to geno file
#' @param K Number of ancestral populations
#' @param iterations Number of iterations
#' @param entropy Whether to compute cross-entropy
#' @return sNMF project object
run_snmf_analysis <- function(geno_file, K = 3, iterations = 100, entropy = TRUE) {
  tryCatch({
    cat("Running sNMF analysis (K=", K, ")...\n", sep = "")
    
    # Run sNMF analysis
    project <- LEA::snmf(
      input.file = geno_file,
      K = K,
      entropy = entropy,
      repetitions = 1,
      project = "new",
      iterations = iterations,
      CPU = 1,
      seed = 12345
    )
    
    cat("  sNMF analysis completed\n")
    
    return(project)
    
  }, error = function(e) {
    cat("Error in sNMF analysis:", e$message, "\n")
    
    # Try alternative approach
    tryCatch({
      cat("  Trying alternative method...\n")
      
      # Read geno file using LEA
      geno_data <- LEA::read.geno(geno_file)
      
      # Run sNMF with matrix input
      project <- LEA::snmf(
        geno_data,
        K = K,
        entropy = entropy,
        repetitions = 1,
        project = "new",
        iterations = iterations,
        CPU = 1,
        seed = 12345
      )
      
      return(project)
      
    }, error = function(e2) {
      cat("Alternative method also failed:", e2$message, "\n")
      return(NULL)
    })
  })
}


#' Complete hybrid analysis pipeline
#' @param gl_object Genlight object
#' @param K_values Vector of K values to test
#' @param pop_assignments Optional population assignments for validation
#' @return List with hybrid analysis results
perform_hybrid_analysis <- function(gl_object, K_values = 2:5, pop_assignments = NULL) {
  tryCatch({
    cat("\n=== Starting Hybrid Analysis ===\n")
    
    # Convert to LEA format
    geno_file <- convert_genlight_to_lea(gl_object)
    
    if (is.null(geno_file)) {
      stop("Failed to convert genotype data to LEA format")
    }
    
    # Run sNMF for each K value
    results <- list()
    cross_entropy <- numeric(length(K_values))
    
    for (i in seq_along(K_values)) {
      K <- K_values[i]
      cat("  Testing K =", K, "\n")
      
      project <- run_snmf_analysis(geno_file, K = K, iterations = 200)
      
      if (!is.null(project)) {
        results[[as.character(K)]] <- project
        
        # Get cross-entropy
        ce <- LEA::cross.entropy(project, K = K)
        cross_entropy[i] <- mean(ce, na.rm = TRUE)
        
        cat("    Cross-entropy:", round(cross_entropy[i], 4), "\n")
      }
    }
    
    # Find best K (lowest cross-entropy)
    if (length(cross_entropy) > 0 && any(!is.na(cross_entropy))) {
      best_k_index <- which.min(cross_entropy)
      best_K <- K_values[best_k_index]
      cat("  Best K:", best_K, "(cross-entropy:", round(cross_entropy[best_k_index], 4), ")\n")
    } else {
      best_K <- NULL
    }
    
    # Get Q matrix for best K
    q_matrix <- NULL
    if (!is.null(best_K) && !is.null(results[[as.character(best_K)]])) {
      project <- results[[as.character(best_K)]]
      q_matrix <- LEA::Q(project, K = best_K)
      
      # Add sample names
      rownames(q_matrix) <- indNames(gl_object)
      colnames(q_matrix) <- paste0("Ancestry", 1:best_K)
    }
    
    # Create ancestry plot data
    ancestry_df <- NULL
    if (!is.null(q_matrix)) {
      ancestry_df <- as.data.frame(q_matrix)
      ancestry_df$Sample <- rownames(q_matrix)
      
      # Melt for plotting
      ancestry_long <- reshape2::melt(ancestry_df, 
                                      id.vars = "Sample",
                                      variable.name = "Ancestry",
                                      value.name = "Proportion")
      
      # Order samples by major ancestry
      if (!is.null(pop_assignments) && length(pop_assignments) == nrow(q_matrix)) {
        ancestry_df$Population <- pop_assignments
        ancestry_long$Population <- rep(pop_assignments, times = best_K)
      }
    }
    
    return(list(
      geno_file = geno_file,
      results = results,
      cross_entropy = data.frame(K = K_values, CrossEntropy = cross_entropy),
      best_K = best_K,
      q_matrix = q_matrix,
      ancestry_df = ancestry_df,
      ancestry_long = ancestry_long
    ))
    
  }, error = function(e) {
    cat("Error in hybrid analysis:", e$message, "\n")
    return(NULL)
  })
}



#' Plot ancestry proportions
#' @param hybrid_results Results from perform_hybrid_analysis
#' @param pop_colors Optional color vector for populations
#' @return ggplot object
plot_ancestry_proportions <- function(hybrid_results, pop_colors = NULL) {
  if (is.null(hybrid_results) || is.null(hybrid_results$ancestry_long)) {
    return(NULL)
  }
  
  tryCatch({
    data <- hybrid_results$ancestry_long
    best_K <- hybrid_results$best_K
    
    # Create plot
    p <- ggplot(data, aes(x = Sample, y = Proportion, fill = Ancestry)) +
      geom_bar(stat = "identity", width = 1) +
      labs(title = paste("Ancestry Proportions (K =", best_K, ")"),
           x = "Sample", 
           y = "Ancestry Proportion") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
    # Add population colors if available
    if ("Population" %in% colnames(data) && !is.null(pop_colors)) {
      p <- p + facet_grid(~ Population, scales = "free_x", space = "free_x")
    }
    
    return(p)
    
  }, error = function(e) {
    cat("Error plotting ancestry proportions:", e$message, "\n")
    return(NULL)
  })
}



#' Alternative: Run ADMIXTURE-style analysis without LEA package
run_admixture_style_analysis <- function(geno_matrix, K_min, K_max, 
                                         replicates, output_prefix) {
  tryCatch({
    cat("    Running ADMIXTURE-style analysis...\n")
    
    # Convert to 0,1,2 format if needed
    # geno_matrix should already be in 0,1,2 format with -9 for missing
    
    # Simple PCA-based clustering as alternative
    pca_results <- perform_pca_clustering(geno_matrix, K_min, K_max)
    
    if (is.null(pca_results)) {
      return(list(success = FALSE, error = "PCA clustering failed"))
    }
    
    # Create mock Q matrices from PCA clusters
    Q_matrices <- list()
    for (K in K_min:K_max) {
      # Use k-means on PCA scores to create Q matrix
      if (K <= ncol(pca_results$scores)) {
        kmeans_result <- kmeans(pca_results$scores[, 1:min(10, ncol(pca_results$scores))], 
                                centers = K, nstart = 10)
        
        # Create Q matrix (membership probabilities)
        Q <- matrix(0, nrow = nrow(geno_matrix), ncol = K)
        for (i in 1:nrow(geno_matrix)) {
          Q[i, kmeans_result$cluster[i]] <- 1
        }
        
        # Add some random noise to simulate probabilities
        Q <- Q + matrix(runif(nrow(Q) * ncol(Q), 0, 0.1), nrow = nrow(Q), ncol = ncol(Q))
        Q <- Q / rowSums(Q)  # Normalize to sum to 1
        
        Q_matrices[[as.character(K)]] <- Q
      }
    }
    
    return(list(
      success = TRUE,
      Q_matrices = Q_matrices,
      method = "PCA-based clustering (fallback)",
      cross_entropy = NULL,
      best_K = which.min(sapply(K_min:K_max, function(K) {
        # Simple heuristic: choose K that maximizes between-group variance
        if (K == 1) return(Inf)
        kmeans_result <- kmeans(pca_results$scores[, 1:min(5, ncol(pca_results$scores))], 
                                centers = K, nstart = 10)
        return(kmeans_result$tot.withinss / kmeans_result$betweenss)
      })) + K_min - 1
    ))
    
  }, error = function(e) {
    cat("Error in ADMIXTURE-style analysis:", e$message, "\n")
    return(list(success = FALSE, error = e$message))
  })
}



#' Perform PCA for clustering
perform_pca_clustering <- function(geno_matrix, K_min, K_max) {
  tryCatch({
    # Remove missing data
    geno_clean <- geno_matrix
    geno_clean[geno_clean == -9] <- NA
    
    # Impute missing values with column mean
    for (j in 1:ncol(geno_clean)) {
      col_mean <- mean(geno_clean[, j], na.rm = TRUE)
      geno_clean[is.na(geno_clean[, j]), j] <- col_mean
    }
    
    # Perform PCA
    pca <- prcomp(geno_clean, center = TRUE, scale. = TRUE)
    
    return(list(
      scores = pca$x,
      variance = pca$sdev^2 / sum(pca$sdev^2) * 100
    ))
    
  }, error = function(e) {
    cat("Error in PCA:", e$message, "\n")
    return(NULL)
  })
}



# Helper function to calculate admixture statistics
calculate_admixture_stats <- function(q_matrix) {
  if (is.null(q_matrix)) return(NULL)
  
  stats <- list(
    n_samples = nrow(q_matrix),
    n_clusters = ncol(q_matrix),
    max_proportions = apply(q_matrix, 1, max),
    mean_max_proportion = mean(apply(q_matrix, 1, max)),
    highly_admixed = sum(apply(q_matrix, 1, max) < 0.7),
    admixed = sum(apply(q_matrix, 1, max) >= 0.7 & apply(q_matrix, 1, max) < 0.8),
    mostly_pure = sum(apply(q_matrix, 1, max) >= 0.8 & apply(q_matrix, 1, max) < 0.9),
    pure = sum(apply(q_matrix, 1, max) >= 0.9)
  )
  
  return(stats)
}

# ============================================================================
# CORRECTED LEA ANALYSIS SERVER LOGIC
# ============================================================================

#' Setup LEA server logic
setup_lea_server <- function(input, output, session, values) {
  
  # Reactive values for LEA
  lea_values <- reactiveValues(
    project = NULL,
    q_matrices = list(),
    final_q_marix_df = NULL,
    cross_entropy = NULL,
    sample_names = NULL,
    current_k = NULL,
    k_range = 2:10
  )
  
  # Update K selector when results are available
  observe({
    if (!is.null(lea_values$k_range)) {
      updateSelectInput(
        session = session,
        inputId = "lea_plot_k",
        choices = lea_values$k_range,
        selected = ifelse(!is.null(lea_values$current_k), 
                          lea_values$current_k, 
                          min(lea_values$k_range))
      )
    }
  })
  
  # Main LEA analysis function
  run_lea_analysis <- eventReactive(input$run_lea, {
    tryCatch({
      # Validate input data
      if (is.null(values$filtered_gl)) {
        stop("No genotype data available. Please load and filter genotype data first.")
      }
      
      # Get parameters
      k_min <- input$lea_min_k
      k_max <- input$lea_max_k
      replicates <- input$lea_replicates
      alpha <- input$lea_alpha
      calculate_entropy <- input$lea_calculate_entropy
      
      # Validate K range
      if (k_min < 1 || k_max < k_min) {
        stop("Invalid K range. Minimum K must be ≥ 1 and maximum K must be ≥ minimum K.")
      }
      
      k_range <- k_min:k_max
      
      # Show progress
      showNotification("Starting LEA analysis...", type = "message", duration = 3)
      
      # Convert genlight to matrix
      geno_matrix <- as.matrix(values$filtered_gl)
      sample_names <- indNames(values$filtered_gl)
      if (is.null(sample_names)) {
        sample_names <- paste0("Sample_", 1:nrow(geno_matrix))
      }
      
      # Handle missing values (LEA uses 9 for missing)
      geno_matrix[is.na(geno_matrix)] <- 9
      
      # Create temporary directory
      lea_dir <- tempfile(pattern = "lea_")
      dir.create(lea_dir, recursive = TRUE)
      
      # Create input file for LEA
      input_file <- file.path(lea_dir, "genotypes.geno")
      
      # Write in LEA format (one row per SNP, one column per individual)
      write.table(t(geno_matrix), input_file, 
                  sep = "", row.names = FALSE, col.names = FALSE)
      
      # Display analysis parameters
      cat("\n=== LEA Analysis Parameters ===\n")
      cat("K range:", paste(k_range, collapse = ", "), "\n")
      cat("Replicates per K:", replicates, "\n")
      cat("Alpha parameter:", alpha, "\n")
      cat("Samples:", nrow(geno_matrix), "\n")
      cat("SNPs:", ncol(geno_matrix), "\n")
      cat("Calculate entropy:", calculate_entropy, "\n")
      
      # Run sNMF analysis
      showNotification("Running sNMF analysis... This may take several minutes.", 
                       type = "warning", duration = 10)
      
      # Run sNMF with error handling
      lea_project <- tryCatch({
        LEA::snmf(
          input.file = input_file,
          K = k_range,
          entropy = calculate_entropy,
          repetitions = replicates,
          alpha = alpha,
          project = "new",
          CPU = 1,
          tolerance = 0.00001,
          iterations = 200
        )
      }, error = function(e) {
        # Try with fewer replicates if first attempt fails
        showNotification("First attempt failed, trying with fewer replicates...", 
                         type = "warning", duration = 5)
        LEA::snmf(
          input.file = input_file,
          K = k_range,
          entropy = calculate_entropy,
          repetitions = min(5, replicates),  # Reduced replicates
          alpha = alpha,
          project = "new",
          CPU = 1,
          tolerance = 0.0001,
          iterations = 100
        )
      })
      
      # Clean up temporary files unless export is requested
      if (!input$lea_export_data) {
        unlink(lea_dir, recursive = TRUE)
      } else {
        cat("Intermediate files saved to:", lea_dir, "\n")
      }
      
      # Extract Q matrices for each K
      q_matrices <- list()
      for (k in k_range) {
        # Get the best run for this K by checking cross-entropy
        if (calculate_entropy) {
          # Get cross-entropy values for all runs of this K
          ce_values <- tryCatch({
            LEA::cross.entropy(lea_project, K = k)
          }, error = function(e) {
            cat("Warning: Could not get cross-entropy for K =", k, ":", e$message, "\n")
            return(rep(NA, replicates))
          })
          
          # Find best run (minimum cross-entropy)
          if (!all(is.na(ce_values))) {
            best_run <- which.min(ce_values)
            cat("  K =", k, ": Best run =", best_run, 
                "(CE =", round(min(ce_values, na.rm = TRUE), 4), ")\n")
          } else {
            best_run <- 1  # Default to first run
          }
        } else {
          best_run <- 1  # Default to first run
        }
        
        # Extract Q matrix for best run
        q_matrix <- tryCatch({
          LEA::Q(lea_project, K = k, run = best_run)
        }, error = function(e) {
          # Try with run = 1 if specific run fails
          cat("Warning: Could not get Q matrix for run", best_run, 
              "for K =", k, ":", e$message, "\n")
          cat("  Trying with run = 1...\n")
          LEA::Q(lea_project, K = k, run = 1)
        })
        
        if (!is.null(q_matrix)) {
          rownames(q_matrix) <- sample_names
          colnames(q_matrix) <- paste0("Cluster_", 1:k)
          q_matrices[[as.character(k)]] <- q_matrix
        }
      }
      
      # Calculate cross-entropy if requested
      cross_entropy <- NULL
      optimal_k <- k_min
      
      if (calculate_entropy) {
        cross_entropy <- data.frame(
          K = k_range,
          Mean_CE = NA,
          SD_CE = NA,
          Best_Run = NA
        )
        
        for (k in k_range) {
          ce_values <- tryCatch({
            LEA::cross.entropy(lea_project, K = k)
          }, error = function(e) {
            cat("Error calculating CE for K =", k, ":", e$message, "\n")
            return(rep(NA, replicates))
          })
          
          if (!all(is.na(ce_values))) {
            cross_entropy$Mean_CE[cross_entropy$K == k] <- mean(ce_values, na.rm = TRUE)
            cross_entropy$SD_CE[cross_entropy$K == k] <- sd(ce_values, na.rm = TRUE)
            cross_entropy$Best_Run[cross_entropy$K == k] <- which.min(ce_values)
          }
        }
        
        # Determine optimal K (lowest cross-entropy)
        if (!all(is.na(cross_entropy$Mean_CE))) {
          optimal_k <- cross_entropy$K[which.min(cross_entropy$Mean_CE)]
        }
      }
      
      # Update reactive values
      lea_values$project <- lea_project
      lea_values$q_matrices <- q_matrices
      lea_values$cross_entropy <- cross_entropy
      lea_values$sample_names <- sample_names
      lea_values$current_k <- optimal_k
      lea_values$k_range <- k_range
      
      # Create summary text
      summary_text <- paste(
        "✅ LEA Analysis Complete\n",
        "=====================\n",
        "K range analyzed: ", k_min, " to ", k_max, "\n",
        "Optimal K: ", optimal_k, "\n",
        "Total samples: ", length(sample_names), "\n",
        "Total SNPs: ", ncol(geno_matrix), "\n",
        "Replicates per K: ", replicates, "\n",
        "Alpha parameter: ", alpha, "\n",
        "Q matrices extracted for K: ", paste(names(q_matrices), collapse = ", "), "\n",
        "Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
        sep = ""
      )
      
      showNotification("LEA analysis completed successfully!", 
                       type = "message", duration = 5)
      
      return(list(
        summary = summary_text,
        optimal_k = optimal_k,
        success = TRUE
      ))
      
    }, error = function(e) {
      error_msg <- paste("❌ LEA analysis failed:", e$message)
      showNotification(error_msg, type = "error", duration = 10)
      cat(error_msg, "\n")
      return(list(
        summary = error_msg,
        success = FALSE
      ))
    })
  })
  
  # Status output
  output$lea_status <- renderText({
    results <- run_lea_analysis()
    if (!is.null(results)) {
      results$summary
    } else {
      "No LEA analysis run yet. Click 'Run LEA Analysis' to start."
    }
  })
  
  # Cross-entropy plot
  output$lea_cross_entropy_plot <- renderPlot({
    if (!is.null(lea_values$cross_entropy)) {
      ce_data <- lea_values$cross_entropy
      
      # Create plot
      p <- ggplot(ce_data, aes(x = K, y = Mean_CE)) +
        geom_point(size = 3, color = "steelblue") +
        geom_line(color = "steelblue", size = 1) +
        geom_errorbar(aes(ymin = Mean_CE - SD_CE, ymax = Mean_CE + SD_CE), 
                      width = 0.2, color = "steelblue", alpha = 0.7) +
        labs(
          title = "Cross-Entropy Analysis",
          subtitle = "Lower values indicate better model fit",
          x = "Number of Clusters (K)",
          y = "Cross-Entropy"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)
        )
      
      # Add optimal K annotation if available
      if (!is.null(lea_values$current_k)) {
        optimal_row <- ce_data[ce_data$K == lea_values$current_k, ]
        if (nrow(optimal_row) > 0) {
          p <- p +
            geom_vline(xintercept = lea_values$current_k, 
                       linetype = "dashed", color = "red", size = 1) +
            annotate("point", x = lea_values$current_k, 
                     y = optimal_row$Mean_CE, 
                     color = "red", size = 4) +
            annotate("text", x = lea_values$current_k, 
                     y = max(ce_data$Mean_CE, na.rm = TRUE) * 0.95,
                     label = paste("Optimal K =", lea_values$current_k),
                     color = "red", size = 4, fontface = "bold")
        }
      }
      
      return(p)
    } else {
      # Return empty plot with message
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Enable 'Calculate cross-entropy' option\nand run LEA analysis to see this plot", 
                 size = 5, color = "gray50") +
        theme_void() +
        theme(plot.background = element_rect(fill = "white"))
    }
  })
  
  # Structure plot
  output$lea_structure_plot <- renderPlot({
    req(lea_values$q_matrices)
    req(input$lea_plot_k)
    
    k <- as.numeric(input$lea_plot_k)
    k_str <- as.character(k)
    q_matrix <- lea_values$q_matrices[[k_str]]
    
    if (is.null(q_matrix)) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, 
                   label = paste("No Q matrix available for K =", k), 
                   size = 6, color = "gray50") +
          theme_void()
      )
    }
    
    # Prepare data for plotting
    plot_data <- data.frame(
      Sample = rep(rownames(q_matrix), times = ncol(q_matrix)),
      Cluster = rep(colnames(q_matrix), each = nrow(q_matrix)),
      Proportion = as.vector(q_matrix)
    )
    
    # Order samples by dominant cluster for better visualization
    dominant_cluster <- apply(q_matrix, 1, which.max)
    sample_order <- order(dominant_cluster)
    
    plot_data$Sample <- factor(plot_data$Sample, 
                               levels = rownames(q_matrix)[sample_order])
    plot_data$Cluster <- factor(plot_data$Cluster, levels = colnames(q_matrix))
    
    # Create structure plot
    p <- ggplot(plot_data, aes(x = Sample, y = Proportion, fill = Cluster)) +
      geom_bar(stat = "identity", width = 1) +
      scale_fill_viridis_d(name = "Cluster", guide = guide_legend(ncol = 1)) +
      labs(
        title = paste("Population Structure (K =", k, ")"),
        x = "Samples",
        y = "Ancestry Proportion"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
    return(p)
  })
  
  # Q matrix table
  output$lea_q_matrix_table <- renderDT({
    req(lea_values$q_matrices)
    req(input$lea_plot_k)
    
    k <- as.numeric(input$lea_plot_k)
    k_str <- as.character(k)
    q_matrix <- lea_values$q_matrices[[k_str]]
    
    if (is.null(q_matrix)) {
      return(NULL)
    }
    
    # Format for display
    display_df <- as.data.frame(q_matrix)
    display_df$Sample <- rownames(q_matrix)
    display_df$Dominant_Cluster <- apply(q_matrix, 1, function(x) {
      colnames(q_matrix)[which.max(x)]
    })
    display_df$Max_Proportion <- apply(q_matrix, 1, max)
    
    # Reorder columns
    display_df <- display_df[, c("Sample", "Dominant_Cluster", "Max_Proportion", 
                                 colnames(q_matrix))]
    
    # Format proportions
    for (col in colnames(q_matrix)) {
      display_df[[col]] <- round(display_df[[col]], 4)
    }
    display_df$Max_Proportion <- round(display_df$Max_Proportion, 4)
    lea_values$final_q_marix_df <- display_df
    
    # Create datatable
    datatable(
      display_df,
      extensions = c('Buttons', 'Scroller'),
      options = list(
        pageLength = 10,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'print'),
        scrollX = TRUE,
        scrollY = 400,
        scroller = TRUE
      ),
      rownames = FALSE,
      class = 'display nowrap'
    ) %>%
      formatStyle(
        'Max_Proportion',
        background = styleColorBar(c(0, 1), 'lightblue'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  # Admixture plot
  output$lea_admixture_plot <- renderPlot({
    req(lea_values$q_matrices)
    req(input$lea_plot_k)
    
    k <- as.numeric(input$lea_plot_k)
    k_str <- as.character(k)
    q_matrix <- lea_values$q_matrices[[k_str]]
    
    if (is.null(q_matrix)) {
      return(NULL)
    }
    
    # Calculate maximum ancestry proportion for each sample
    max_prop <- apply(q_matrix, 1, max)
    
    # Categorize admixture levels
    admixture_level <- cut(max_prop, 
                           breaks = c(0, 0.7, 0.8, 0.9, 1.0),
                           labels = c("Highly Admixed (<0.7)", 
                                      "Admixed (0.7-0.8)", 
                                      "Mostly Pure (0.8-0.9)", 
                                      "Pure (>0.9)"))
    
    admix_df <- as.data.frame(table(admixture_level))
    colnames(admix_df) <- c("Category", "Count")
    
    # Create bar plot
    p <- ggplot(admix_df, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Count), vjust = -0.5, size = 5) +
      scale_fill_viridis_d(name = "Admixture Level") +
      labs(
        title = paste("Admixture Distribution (K =", k, ")"),
        subtitle = paste("Total samples:", nrow(q_matrix)),
        x = "Admixture Level",
        y = "Number of Samples"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none"
      )
    
    return(p)
  })
  
  # Optimal K text output
  output$lea_optimal_k <- renderUI({
    if (!is.null(lea_values$current_k)) {
      tags$div(
        class = "alert alert-success",
        style = "font-size: 16px; font-weight: bold;",
        icon("check-circle"), 
        " Optimal K = ", lea_values$current_k
      )
    } else {
      tags$div(
        class = "alert alert-info",
        icon("info-circle"),
        " Run analysis to determine optimal K"
      )
    }
  })
  
  # Export buttons UI
  output$lea_export_buttons <- renderUI({
    if (!is.null(lea_values$project)) {
      tagList(
        downloadButton("download_lea_csv", "Export Q Matrices", 
                       class = "btn-success btn-block", 
                       style = "margin-bottom: 10px;"),
        downloadButton("download_lea_rdata", "Export Full Results", 
                       class = "btn-primary btn-block",
                       style = "margin-bottom: 10px;")
      )
    }
  })
  
  # Navigation buttons
  observeEvent(input$lea_prev_k, {
    if (!is.null(lea_values$k_range) && length(lea_values$k_range) > 1) {
      current_idx <- which(lea_values$k_range == as.numeric(input$lea_plot_k))
      if (current_idx > 1) {
        new_k <- lea_values$k_range[current_idx - 1]
        updateSelectInput(session, "lea_plot_k", selected = new_k)
      }
    }
  })
  
  observeEvent(input$lea_next_k, {
    if (!is.null(lea_values$k_range) && length(lea_values$k_range) > 1) {
      current_idx <- which(lea_values$k_range == as.numeric(input$lea_plot_k))
      if (current_idx < length(lea_values$k_range)) {
        new_k <- lea_values$k_range[current_idx + 1]
        updateSelectInput(session, "lea_plot_k", selected = new_k)
      }
    }
  })
  
  # Compare all K values
  observeEvent(input$lea_compare_all, {
    if (!is.null(lea_values$q_matrices) && length(lea_values$q_matrices) > 1) {
      showModal(modalDialog(
        title = "Comparison of All K Values",
        size = "l",
        plotOutput("lea_comparison_plot", height = "600px") %>%
          withSpinner(type = 4, color = "#0dc5c1"),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      
      output$lea_comparison_plot <- renderPlot({
        # Create faceted plot for all K values
        plot_list <- list()
        
        for (k in lea_values$k_range) {
          k_str <- as.character(k)
          q_matrix <- lea_values$q_matrices[[k_str]]
          
          if (!is.null(q_matrix)) {
            plot_data <- data.frame(
              Sample = rep(rownames(q_matrix), times = ncol(q_matrix)),
              Cluster = rep(colnames(q_matrix), each = nrow(q_matrix)),
              Proportion = as.vector(q_matrix)
            )
            
            # Order samples consistently (use order from first K)
            if (k == lea_values$k_range[1]) {
              dominant_cluster <- apply(q_matrix, 1, which.max)
              sample_order <- order(dominant_cluster)
            }
            
            plot_data$Sample <- factor(plot_data$Sample, 
                                       levels = rownames(q_matrix)[sample_order])
            
            p <- ggplot(plot_data, aes(x = Sample, y = Proportion, fill = Cluster)) +
              geom_bar(stat = "identity", width = 1) +
              scale_fill_viridis_d(guide = "none") +
              labs(title = paste("K =", k)) +
              theme_minimal() +
              theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title = element_blank(),
                plot.title = element_text(hjust = 0.5, face = "bold"),
                panel.grid = element_blank()
              )
            
            plot_list[[k_str]] <- p
          }
        }
        
        # Combine plots
        do.call(gridExtra::grid.arrange, c(plot_list, ncol = 2))
      })
    }
  })
  
  # Download handlers
  output$download_lea_csv <- downloadHandler(
    filename = function() {
      paste0("lea_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      # Create temporary directory
      temp_dir <- tempdir()
      csv_files <- c()
      
      # Save each Q matrix as CSV
      for (k_str in names(lea_values$q_matrices)) {
        q_matrix <- lea_values$q_matrices[[k_str]]
        csv_file <- file.path(temp_dir, paste0("Q_matrix_K", k_str, ".csv"))
        write.csv(q_matrix, csv_file, row.names = TRUE)
        csv_files <- c(csv_files, csv_file)
      }
      
      # Save cross-entropy data
      if (!is.null(lea_values$cross_entropy)) {
        ce_file <- file.path(temp_dir, "cross_entropy.csv")
        write.csv(lea_values$cross_entropy, ce_file, row.names = FALSE)
        csv_files <- c(csv_files, ce_file)
      }
      
      # Create ZIP file
      zip(file, csv_files, flags = "-j")
    }
  )
  
  output$download_lea_rdata <- downloadHandler(
    filename = function() {
      paste0("lea_full_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
    },
    content = function(file) {
      # Save all relevant data
      save_data <- list(
        project = lea_values$project,
        q_matrices = lea_values$q_matrices,
        cross_entropy = lea_values$cross_entropy,
        k_range = lea_values$k_range,
        optimal_k = lea_values$current_k,
        sample_names = lea_values$sample_names,
        timestamp = Sys.time()
      )
      save(save_data, file = file)
    }
  )
  
  output$download_lea_plot <- downloadHandler(
    filename = function() {
      paste0("lea_structure_plot_K", input$lea_plot_k, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # Get the current plot
      k <- as.numeric(input$lea_plot_k)
      k_str <- as.character(k)
      q_matrix <- lea_values$q_matrices[[k_str]]
      
      if (is.null(q_matrix)) {
        return(NULL)
      }
      
      # Prepare data for plotting
      plot_data <- data.frame(
        Sample = rep(rownames(q_matrix), times = ncol(q_matrix)),
        Cluster = rep(colnames(q_matrix), each = nrow(q_matrix)),
        Proportion = as.vector(q_matrix)
      )
      
      # Order samples by dominant cluster
      dominant_cluster <- apply(q_matrix, 1, which.max)
      sample_order <- order(dominant_cluster)
      
      plot_data$Sample <- factor(plot_data$Sample, 
                                 levels = rownames(q_matrix)[sample_order])
      plot_data$Cluster <- factor(plot_data$Cluster, levels = colnames(q_matrix))
      
      # Create plot
      p <- ggplot(plot_data, aes(x = Sample, y = Proportion, fill = Cluster)) +
        geom_bar(stat = "identity", width = 1) +
        scale_fill_viridis_d(name = "Cluster", guide = guide_legend(ncol = 1)) +
        labs(
          title = paste("Population Structure (K =", k, ")"),
          x = "Samples",
          y = "Ancestry Proportion"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
      
      # Save plot
      ggsave(file, plot = p, width = 14, height = 8, dpi = 300)
    }
  )
}


#' Convert genlight to LEA format - FIXED VERSION
convert_to_lea_format_fixed <- function(gl_object, output_file = NULL) {
  tryCatch({
    cat("  Converting genlight to LEA format...\n")
    
    if (is.null(gl_object)) {
      stop("Genlight object is NULL")
    }
    
    # Get basic information
    n_ind <- nInd(gl_object)
    n_loc <- nLoc(gl_object)
    
    cat("    Individuals:", n_ind, "\n")
    cat("    Loci:", n_loc, "\n")
    
    # Convert to matrix (0,1,2 coding)
    geno_matrix <- as.matrix(gl_object)
    
    # Check matrix dimensions
    if (is.null(dim(geno_matrix)) || nrow(geno_matrix) != n_ind || ncol(geno_matrix) != n_loc) {
      cat("    Warning: Matrix dimensions don't match genlight object\n")
      cat("    Matrix dim:", dim(geno_matrix), "\n")
      
      # Try to fix by transposing if needed
      if (nrow(geno_matrix) == n_loc && ncol(geno_matrix) == n_ind) {
        cat("    Transposing matrix...\n")
        geno_matrix <- t(geno_matrix)
      }
    }
    
    # Convert to numeric if needed
    if (!is.numeric(geno_matrix)) {
      cat("    Converting to numeric...\n")
      geno_matrix <- apply(geno_matrix, 2, as.numeric)
    }
    
    # Handle missing values (replace NA with 9 as per LEA format)
    if (any(is.na(geno_matrix))) {
      cat("    Replacing", sum(is.na(geno_matrix)), "NA values with 9...\n")
      geno_matrix[is.na(geno_matrix)] <- 9
    }
    
    # Ensure values are integers (0, 1, 2, 9)
    if (!all(geno_matrix %in% c(0, 1, 2, 9), na.rm = TRUE)) {
      cat("    Warning: Some values are not 0,1,2,9. Rounding to nearest integer...\n")
      geno_matrix <- round(geno_matrix)
      geno_matrix[geno_matrix < 0] <- 0
      geno_matrix[geno_matrix > 2 & geno_matrix != 9] <- 2
    }
    
    # Transpose for LEA format (SNPs as rows, individuals as columns)
    cat("    Transposing matrix for LEA format...\n")
    geno_matrix_t <- t(geno_matrix)
    
    # Ensure integer storage mode
    storage.mode(geno_matrix_t) <- "integer"
    
    # Create output file name
    if (is.null(output_file)) {
      output_file <- tempfile(fileext = ".geno")
      cat("    Using temporary file:", output_file, "\n")
    }
    
    # Write to file in LEA format
    cat("    Writing to file...\n")
    
    # Method 1: Try using LEA's write.geno function
    if (requireNamespace("LEA", quietly = TRUE)) {
      tryCatch({
        # Write in GENO format (more reliable than LFMM)
        LEA::write.geno(geno_matrix_t, output_file)
        cat("    Written using LEA::write.geno\n")
      }, error = function(e) {
        cat("    LEA::write.geno failed, using manual write...\n")
        # Manual write
        write.table(geno_matrix_t, 
                    file = output_file,
                    sep = "",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    na = "9")
      })
    } else {
      # Manual write
      write.table(geno_matrix_t, 
                  file = output_file,
                  sep = "",
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE,
                  na = "9")
    }
    
    # Verify file was created
    if (!file.exists(output_file)) {
      stop("Failed to create output file")
    }
    
    file_size <- file.info(output_file)$size
    cat("    File created:", output_file, "(", file_size, "bytes)\n")
    
    # Create lfmm object (using the transposed matrix)
    lfmm_obj <- tryCatch({
      LEA::as.lfmm(geno_matrix_t)
    }, error = function(e) {
      cat("    Warning: Could not create lfmm object:", e$message, "\n")
      NULL
    })
    
    return(list(
      file_path = output_file,
      lfmm_object = lfmm_obj,
      n_individuals = n_ind,
      n_snps = n_loc,
      matrix_dimensions = dim(geno_matrix_t)
    ))
    
  }, error = function(e) {
    cat("  ERROR in LEA conversion:", e$message, "\n")
    cat("  Traceback:\n")
    print(traceback())
    return(NULL)
  })
}




# ============================================================================
# SECTION 3: SHINY APP UI
# ============================================================================

ui <- dashboardPage(
  skin = "blue",
  
  # Header
  dashboardHeader(
    title = "🌱 Plant Breeding Analysis Platform v4.0",
    titleWidth = 300,
    dropdownMenu(
      type = "notifications",
      badgeStatus = "warning",
      icon = icon("bell"),
      notificationItem(
        text = "Welcome to GWAS Platform",
        icon = icon("info-circle"),
        status = "success"
      )
    )
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("🔍 About", tabName = "help", icon = icon("question")),
      menuItem("📁 Data Upload", tabName = "data_upload", icon = icon("upload")),
      menuItem("🔧 Data Processing", tabName = "data_processing", icon = icon("filter")),
      menuItem("🚀 Export QTL Data", tabName = "export_qtl", icon = icon("arrow-up")),
      menuItem("🧬 Population Genetics", tabName = "population_genetics", icon = icon("project-diagram"),
               menuSubItem("PCA & Kinship Analysis", tabName = "pca_kinship"),
               menuSubItem("Allele Frequency Analysis", tabName = "allele_freq_analysis"),
               menuSubItem("Population Alleles", tabName = "pop_allele_analysis"),
               menuSubItem("Population Correlations", tabName = "pop_correlations"),
               menuSubItem("Fst-Fis Analysis", tabName = "fst_fis_analysis"),
               menuSubItem("LEA Analysis", tabName = "lea_analysis")
      ),
      menuItem("📊 GWAS Analysis", tabName = "gwas_analysis", icon = icon("dna")),
      menuItem("🔄 Multi-Trait Analysis", tabName = "multi_trait", icon = icon("layer-group"),
               menuSubItem("Multi-Trait GWAS", tabName = "multi_trait_gwas"),
               menuSubItem("Unified Results", tabName = "unified_results"),
               menuSubItem("Enhanced Analysis", tabName = "enhanced_analysis")
      ),
      # Diallel Analysis Menu with Sub-items
      menuItem("🌱 Diallel Analysis", tabName = "diallel", icon = icon("seedling"),
               menuSubItem("Diallel Analysis", tabName = "diallel_analysis"),
               # menuSubItem("Parent Relationship", tabName = "diallel_network_viz"),
               # menuSubItem("Diallel All Traits", tabName = "diallel_all_traits"),
               # menuSubItem("GCA/SCA Effects", tabName = "diallel_gca_sca"),
               menuSubItem("Heterosis Analysis", tabName = "diallel_heterosis"),
               menuSubItem("Cross-SNP Mapping", tabName = "cross_snp_mapping"),
               menuSubItem("Variety Selection", tabName = "diallel_selection"),
               menuSubItem("Parents Actual Genotypes", tabName = "parent_snp_actual")
      ),
      menuItem("METAN Analysis", tabName = "metan", icon = icon("globe-asia")),
      menuItem("Reports", tabName = "report", icon = icon("file-alt")),
      menuItem("📊 Dataframes Viewer", tabName = "dataframes_viewer", icon = icon("database"),
               menuSubItem("All DataFrames", tabName = "all_dataframes"),
               menuSubItem("GWAS DataFrames", tabName = "gwas_dataframes"),
               menuSubItem("Diallel DataFrames", tabName = "diallel_dataframes"),
               menuSubItem("Cross-SNP DataFrames", tabName = "cross_snp_dataframes"),
               menuSubItem("METAN DataFrames", tabName = "metan_dataframes")
      ),
      
      # Analysis controls
      div(style = "padding: 15px;",
          h4("Analysis Status"),
          verbatimTextOutput("analysis_status"),
          br(),
          actionButton("reset_all", "Reset All Data", 
                       icon = icon("redo"), 
                       class = "btn-danger btn-block")
      )
    ),
    helpText("Developed by ", a("Kiplagat John Noel", href = "http://www.googlemail.com/"),
             a("Microbiology, Biochemistry and Biotechnology Department,", href = "http://www.ku.ac.ke/"), "Kenyatta University (KE)",
             style = "padding-left:3em;padding-right:1em;position:absolute;bottom:1em;")
  ),
  
  # Body
  dashboardBody(
    # Custom CSS
    tags$head(
      tags$style(HTML("
      /* Target both modal classes used by Shiny */
    .shiny-modal-progress,
      div.shiny-modal-progress.modal {
      top: unset !important;
      bottom: 20px !important;
      left: unset !important;
      right: 20px !important;
      transform: none !important;
      margin: 0 !important;
      width: auto !important;
      height: auto !important;
    }
    
    /* Adjust the modal dialog container */
    .shiny-modal-progress .modal-dialog {
      margin: 0 !important;
      width: 350px !important; /* Fixed width */
    }
    
    /* Ensure content looks good */
    .shiny-modal-progress .modal-content {
      padding: 10px;
      background-color: rgba(255, 255, 255, 0.95);
      border-radius: 8px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }
    
        .content-wrapper, .right-side {
          background-color: #f9f9f9;
        }
        .box {
          border-radius: 5px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.1);
          margin-bottom: 15px;
        }
        .box-header {
          background-color: #f5f5f5;
          border-bottom: 1px solid #ddd;
          padding: 10px 15px;
        }
        .box-title {
          font-weight: bold;
          color: #333;
        }
        .btn-primary {
          background-color: #3498db;
          border-color: #2980b9;
        }
        .btn-success {
          background-color: #2ecc71;
          border-color: #27ae60;
        }
        .btn-warning {
          background-color: #f39c12;
          border-color: #e67e22;
        }
        .btn-danger {
          background-color: #e74c3c;
          border-color: #c0392b;
        }
        .btn-info {
          background-color: #17a2b8;
          border-color: #117a8b;
        }
        .alert {
          border-radius: 5px;
          margin-bottom: 15px;
          padding: 10px 15px;
        }
        .alert-success {
          background-color: #d4edda;
          border-color: #c3e6cb;
          color: #155724;
        }
        .alert-warning {
          background-color: #fff3cd;
          border-color: #ffeaa7;
          color: #856404;
        }
        .alert-danger {
          background-color: #f8d7da;
          border-color: #f5c6cb;
          color: #721c24;
        }
        .alert-info {
          background-color: #d1ecf1;
          border-color: #bee5eb;
          color: #0c5460;
        }
        .table {
          font-size: 14px;
        }
        .shiny-notification {
          position: fixed;
          top: 20px;
          right: 20px;
          width: 350px;
          z-index: 9999;
          font-size: 14px;
        }
        .value-box {
          border-radius: 5px;
          margin-bottom: 10px;
        }
        .nav-tabs-custom .nav-tabs li.active {
          border-top-color: #3498db;
        }
        .progress {
          height: 20px;
          margin-bottom: 10px;
        }
        .help-block {
          font-size: 12px;
          color: #666;
          margin-top: 5px;
        }
        h2, h3, h4 {
          color: #2c3e50;
          font-weight: 600;
        }
        .well {
          background-color: #f8f9fa;
          border: 1px solid #e9ecef;
          border-radius: 5px;
          padding: 15px;
          margin-bottom: 15px;
        }
        .form-group {
          margin-bottom: 15px;
        }
        .selectize-input {
          border-radius: 4px;
          border: 1px solid #ced4da;
        }
        .slider-container {
          padding: 10px 0;
        }
        
        .plot-container {
               background-color: white;
               padding: 15px;
               border-radius: 5px;
               box-shadow: 0 2px 4px rgba(0,0,0,0.1);
               margin-bottom: 20px;
             }
             .data-container {
               background-color: white;
               padding: 15px;
               border-radius: 5px;
               box-shadow: 0 2px 4px rgba(0,0,0,0.1);
             }
             .summary-box {
               background-color: #f8f9fa;
               border-left: 4px solid #007bff;
               padding: 10px;
               margin-bottom: 15px;
             }
             
             .analysis-box {
               background-color: #f8f9fa;
               padding: 15px;
               border-radius: 5px;
               margin-bottom: 20px;
               border-left: 4px solid #007bff;
             }
             .data-box {
               background-color: white;
               padding: 15px;
               border-radius: 5px;
               box-shadow: 0 2px 4px rgba(0,0,0,0.1);
               margin-bottom: 20px;
             }
             
             /* Population Genetics specific styling */
             .popgen-box {
              background-color: #f0f8ff;
              border-left: 4px solid #4682b4;
              padding: 15px;
              margin-bottom: 20px;
             }

            .popgen-table {
             font-size: 12px;
             }

            .popgen-plot {
              background-color: white;
              padding: 10px;
              border-radius: 5px;
              box-shadow: 0 2px 4px rgba(0,0,0,0.1);

            "))
    ),
    
    # Include shinyjs for better UI control
    shinyjs::useShinyjs(),
    
    # Add loading screen
    tags$div(
      id = "loading-screen",
      style = "position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(255,255,255,0.8); z-index: 9999; display: none;",
      tags$div(
        style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); text-align: center;",
        icon("spinner", class = "fa-spin fa-3x"),
        h3("Loading...", style = "margin-top: 20px;")
      )
    ),
    
    tabItems(
      #Instructions
      tabItem(
        tabName = "help",
        fluidRow(
          column(12,
                 h2("Documentation & Help"),
                 tabsetPanel(
                   tabPanel("dartR Analysis",
                            includeMarkdown("dartR.Rmd")),
                   tabPanel("GWAS Analysis",
                            includeMarkdown("GWAS.Rmd")),
                   tabPanel("Diallel Analysis",
                            includeMarkdown("Diallel.Rmd")),
                   tabPanel("MET Analysis",
                            includeMarkdown("Metan.Rmd"))
                 )
          )
        )
      ),
      
      # Data Upload Tab
      tabItem(
        tabName = "data_upload",
        fluidRow(
          column(
            width = 12,
            h2("📁 Data Upload", style = "margin-top: 0;"),
            helpText("Upload your phenotypic and genotypic data files")
          )
        ),
        
        fluidRow(
          # Phenotypic Data Upload
          box(
            title = "Phenotypic Data", 
            status = "primary", 
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            collapsed = FALSE,
            fileInput(
              "pheno_file", 
              "Upload Phenotype CSV File",
              accept = c(".csv", ".txt", ".xlsx", ".xls"),
              placeholder = "No file selected",
              width = "100%"
            ),
            fluidRow(
              column(9,
                     checkboxInput("pheno_header", "Header row", TRUE),
                     checkboxInput("aggregate_pheno", "Aggregate replicates by genotype", FALSE)
              )
            ),
            actionButton("load_pheno", "Load Phenotype Data", 
                         icon = icon("upload"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            br(),
            uiOutput("pheno_status")
          ),
          
          # Genotypic Data Upload
          box(
            title = "Metadata and Genotypic Data", 
            status = "primary", 
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            collapsed = FALSE,
            fileInput(
              "metadata_file", 
              "Upload Metadata Data (Optional)",
              accept = c(".csv", ".txt"),
              placeholder = "No file selected",
              width = "100%"
            ),
            actionButton("load_metadata", "Load Metadata", 
                         icon = icon("upload"), 
                         class = "btn-info",
                         width = "100%"),
            br(),br(),
            uiOutput("metadata_status"),
            br(),
            fileInput(
              "geno_file", 
              "Upload Genotype File",
              accept = c(".csv", ".xlsx", ".xls", ".txt"),
              placeholder = "No file selected",
              width = "100%"
            ),
            helpText("Supports DArT format, CSV, and Excel files"),
            numericInput("pheno_skip", "Lines to skip", 0, min = 0, width = "100%"),
            actionButton("load_geno", "Load Genotype Data", 
                         icon = icon("upload"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            br(),
            uiOutput("geno_status")
          )
        ),
        
        fluidRow(
          # Optional Data Uploads
          box(
            title = "Optional Data Files", 
            status = "info", 
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            collapsed = FALSE,
            
            h4("METAN Data"),
            fileInput(
              "metan_file", 
              "Upload Separate METAN Data (Optional)",
              accept = c(".csv", ".txt"),
              placeholder = "No file selected",
              width = "100%"
            ),
            fluidRow(
              column(6,
                     actionButton("load_metan", "Load METAN Data", 
                                  icon = icon("upload"), 
                                  class = "btn-info",
                                  width = "100%")
              ),
              column(6,
                     actionButton("clear_metan", "Clear METAN", 
                                  icon = icon("trash"), 
                                  class = "btn-warning",
                                  width = "100%")
              )
            ),
            br(),
            uiOutput("metan_status"),
            br(),
            helpText("If no separate METAN file is loaded, the phenotype data will be used for METAN analysis."),
            br(),
            
            h4("Diallel Data"),
            fileInput(
              "diallel_file", 
              "Upload Diallel Data (Optional)",
              accept = c(".csv", ".txt"),
              placeholder = "No file selected",
              width = "100%"
            ),
            actionButton("load_diallel", "Load Diallel Data", 
                         icon = icon("upload"), 
                         class = "btn-info",
                         width = "100%"),
            br(),br(),
            uiOutput("diallel_status"),
            br()
          ),
          
          # Data Summary
          box(
            title = "Data Summary", 
            status = "success", 
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            collapsed = FALSE,
            withSpinner(verbatimTextOutput("data_summary"), type = 4, color = "#0dc5c1"),
            br(),
            actionButton("check_data", "Check Data Compatibility", 
                         icon = icon("check"), 
                         class = "btn-warning btn-block"),
            br(),
            uiOutput("compatibility_status")
          )
        ),
        
        # Data Preview Section
        fluidRow(
          box(
            title = "Phenotypic Data Preview", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            withSpinner(DTOutput("pheno_preview"), type = 4, color = "#0dc5c1")
          )
        ),
        
        fluidRow(
          box(
            title = "Genotypic Data Summary", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            withSpinner(verbatimTextOutput("geno_summary"), type = 4, color = "#0dc5c1")
          )
        )
      ),
      
      # Data Processing Tab
      tabItem(
        tabName = "data_processing",
        fluidRow(
          column(
            width = 12,
            h2("🔧 Data Processing & Quality Control", style = "margin-top: 0;")
          )
        ),
        
        fluidRow(
          # Quality Control Settings
          box(
            title = "Quality Control Parameters", 
            status = "warning", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            sliderInput("callrate_threshold", 
                        "Call Rate Threshold", 
                        min = 0.5, max = 1, value = 0.8, step = 0.05,
                        width = "100%"),
            helpText("Minimum call rate required for each SNP"),
            sliderInput("maf_threshold", 
                        "Minor Allele Frequency (MAF)", 
                        min = 0, max = 0.5, value = 0.05, step = 0.01,
                        width = "100%"),
            helpText("Minimum minor allele frequency"),
            br(),
            actionButton("run_qc", "Run Quality Control", 
                         icon = icon("filter"), 
                         class = "btn-warning btn-block",
                         width = "100%")
          ),
          
          # Matching Parameters
          box(
            title = "Genotype-Phenotype Matching", 
            status = "info", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            numericInput("start_col", 
                         "Start Column for Genotypes", 
                         value = 22, min = 1, step = 1,
                         width = "100%"),
            helpText("Column number where genotype data starts (skip metadata columns)"),
            br(),
            actionButton("match_data", "Match Genotype & Phenotype", 
                         icon = icon("link"), 
                         class = "btn-info btn-block",
                         width = "100%"),
            br(),
            uiOutput("match_status")
          ),
          
          # Data Summary After Processing
          box(
            title = "Processed Data Summary", 
            status = "success", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            withSpinner(verbatimTextOutput("processed_summary"), type = 4, color = "#0dc5c1"),
            br(),
            actionButton("show_qc_report", "Show QC Report", 
                         icon = icon("chart-bar"), 
                         class = "btn-success",
                         width = "100%"),
            br(), br(),
            downloadButton("download_qc", "Download QC Results", 
                           class = "btn-primary btn-block")
          )
        ),
        
        fluidRow(
          # QC Results Visualization
          box(
            title = "Quality Control Results", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tabsetPanel(
              tabPanel("Missing Data", 
                       withSpinner(plotOutput("missing_data_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("MAF Distribution", 
                       withSpinner(plotOutput("maf_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Heterozygosity", 
                       withSpinner(plotOutput("het_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("QC Summary Table", 
                       withSpinner(DTOutput("qc_table"), 
                                   type = 4, color = "#0dc5c1"))
            )
          )
        )
      ),
      
      # In your UI definition, add:
      tabItem(tabName = "export_qtl",
              fluidRow(
                box(
                  title = "Export QTL ICI Format",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  
                  h4("Export Genotype Data in Standard QTL ICI Format"),
                  p("This will create an Excel file with two worksheets:"),
                  tags$ul(
                    tags$li(tags$strong("GeneralInfo:"), "Metadata about the dataset"),
                    tags$li(tags$strong("Genotype:"), "SNP matrix with samples as columns")
                  ),
                  
                  hr(),
                  
                  # Export status and button
                  uiOutput("qtl_export_status"),
                  
                  hr(),
                  
                  # Preview section
                  h4("Data Preview"),
                  verbatimTextOutput("qtl_preview")
                )
              )
      ),
      
      # Population Structure Tab
      tabItem(
        tabName = "pca_kinship",
        fluidRow(
          column(
            width = 12,
            h2("🧬 PCA and Kinship Analysis", style = "margin-top: 0;")
          )
        ),
        
        fluidRow(
          # PCA Settings
          box(
            title = "Principal Component Analysis (PCA)", 
            status = "primary", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            sliderInput("n_pcs", "Number of PCs to Calculate", 
                        min = 2, max = 20, value = 5, step = 1,
                        width = "100%"),
            fluidRow(
              column(6,
                     checkboxInput("scale_pca", "Scale Variables", TRUE)
              ),
              column(6,
                     checkboxInput("center_pca", "Center Variables", TRUE)
              )
            ),
            br(),
            actionButton("run_pca", "Run PCA", 
                         icon = icon("project-diagram"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            helpText("Performs PCA on genotype data to visualize population structure")
          ),
          
          # Kinship Matrix Settings
          box(
            title = "Kinship Matrix", 
            status = "info", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            selectInput("kinship_method", "Kinship Method",
                        choices = c("vanRaden" = "vanRaden",
                                    "loiselle" = "loiselle",
                                    "simple" = "simple",
                                    "identity" = "identity"),
                        selected = "vanRaden",
                        width = "100%"),
            checkboxInput("scale_kinship", "Scale Kinship Matrix", TRUE),
            br(),
            actionButton("run_kinship", "Calculate Kinship", 
                         icon = icon("users"), 
                         class = "btn-info btn-block",
                         width = "100%"),
            helpText("Calculates genetic relationship matrix")
          ),
          
          # Analysis Status
          box(
            title = "Analysis Status", 
            status = "success", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            withSpinner(verbatimTextOutput("structure_status"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            fluidRow(
              column(6,
                     actionButton("export_pca", "Export PCA Results",
                                  icon = icon("download"), 
                                  class = "btn-success",
                                  width = "100%")
              ),
              column(6,
                     actionButton("export_kinship", "Export Kinship",
                                  icon = icon("download"), 
                                  class = "btn-info",
                                  width = "100%")
              )
            )
          )
        ),
        
        fluidRow(
          # PCA Results
          box(
            title = "PCA Results", 
            status = "primary", 
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            tabsetPanel(
              
              tabPanel("PCA Scores",
                DTOutput("pca_scores_table")
              ),
              
              tabPanel("Scree Plot",
                       withSpinner(plotOutput("pca_scree_plot", height = "300px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("PCA Plot",
                       fluidRow(
                         column(6,
                                selectInput("pca_x", "X-axis PC", 
                                            choices = 1:10, selected = 1,
                                            width = "100%")
                         ),
                         column(6,
                                selectInput("pca_y", "Y-axis PC", 
                                            choices = 1:10, selected = 2,
                                            width = "100%")
                         )
                       ),
                       withSpinner(plotOutput("pca_scatter_plot", height = "400px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Variance Explained",
                       withSpinner(plotOutput("pca_variance_plot", height = "300px"), 
                                   type = 4, color = "#0dc5c1"))
            )
          ),
          
          # Kinship Results
          box(
            title = "Kinship Matrix", 
            status = "info", 
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotOutput("kinship_heatmap", height = "400px"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            withSpinner(verbatimTextOutput("kinship_summary"), 
                        type = 4, color = "#0dc5c1")
          )
        )
      ),
      
      # Allele Frequency Analysis Tab
      tabItem(
        tabName = "allele_freq_analysis",
        h2("Allele Frequency Analysis"),
        
        fluidRow(
          box(width = 12, title = "Settings", status = "primary", solidHeader = TRUE,
              fluidRow(
                column(4,
                       selectInput("pop_gen_data_source", "Data Source",
                                   choices = c("Use processed genotype data" = "processed",
                                               "Upload allele frequency data" = "upload"),
                                   selected = "processed")
                ),
                column(4,
                       uiOutput("pop_selector_ui")
                ),
                column(4,
                       actionButton("run_allele_freq_analysis", "Run Analysis", 
                                    icon = icon("play"), class = "btn-primary",
                                    style = "margin-top: 25px;")
                )
              )
          )
        ),
        
        fluidRow(
          box(width = 12, title = "Allele Frequency Histogram", status = "success", solidHeader = TRUE,
              plotOutput("allele_freq_hist_plot", height = "800px"),
              downloadButton("download_allele_freq_plot", "Download Plot")
          )
        )
      ),
      
      
      # Population Allele Analysis Tab
      tabItem(
        tabName = "pop_allele_analysis",
        fluidRow(
          column(
            width = 12,
            h2("📊 Population Allele Statistics", style = "margin-top: 0;"),
            helpText("Calculate allele diversity statistics per population")
          )
        ),
        
        fluidRow(
          # Analysis Settings Box
          box(
            title = "Analysis Settings", 
            status = "primary", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            
            numericInput("pop_min_inds", "Minimum individuals per population",
                         value = 2, min = 1, max = 50, step = 1,
                         width = "100%"),
            helpText("Populations with fewer individuals will be excluded"),
            
            checkboxInput("pop_calc_heterozygosity", "Calculate heterozygosity", TRUE),
            checkboxInput("pop_use_custom_colors", "Use custom colors", FALSE),
            
            conditionalPanel(
              condition = "input.pop_use_custom_colors == true",
              textAreaInput("pop_custom_colors", "Custom colors (comma-separated)",
                            placeholder = "red,blue,green,yellow,purple",
                            rows = 3, width = "100%"),
              helpText("Enter color names or hex codes separated by commas")
            ),
            
            br(),
            actionButton("run_pop_allele_analysis", "Run Population Analysis", 
                         icon = icon("calculator"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            
            br(),
            downloadButton("download_pop_allele_csv", "Download CSV Results", 
                           class = "btn-success"),
            br(), br(),
            downloadButton("download_pop_allele_report", "Download Full Report", 
                           class = "btn-info")
          ),
          
          # Results Summary Box
          box(
            title = "Analysis Summary", 
            status = "success", 
            solidHeader = TRUE,
            width = 8,
            collapsible = TRUE,
            withSpinner(verbatimTextOutput("pop_allele_summary"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            uiOutput("pop_allele_stats")
          )
        ),
        
        fluidRow(
          # Main Plot Box
          box(
            title = "Main Analysis Plot", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            withSpinner(plotOutput("pop_allele_main_plot", height = "500px"), 
                        type = 4, color = "#0dc5c1"),
            
            fluidRow(
              column(3,
                     numericInput("pop_plot_point_size", "Point size", 
                                  value = 5, min = 1, max = 20, step = 1,
                                  width = "100%")
              ),
              column(3,
                     checkboxInput("pop_plot_labels", "Show population labels", TRUE)
              ),
              column(3,
                     checkboxInput("pop_plot_legend", "Show legend", TRUE)
              ),
              column(3,
                     downloadButton("download_main_plot", "Download Plot", 
                                    class = "btn-success")
              )
            )
          )
        ),
        
        fluidRow(
          # Additional Plots
          box(
            title = "Additional Visualizations", 
            status = "info", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tabsetPanel(
              tabPanel("Alleles Distribution",
                       withSpinner(plotOutput("pop_allele_bar_plot", height = "400px"), 
                                   type = 4, color = "#0dc5c1")
              ),
              tabPanel("Heterozygosity",
                       conditionalPanel(
                         condition = "output.heterozygosity_available == true",
                         withSpinner(plotOutput("pop_het_plot", height = "400px"), 
                                     type = 4, color = "#0dc5c1")
                       ),
                       conditionalPanel(
                         condition = "output.heterozygosity_available == false",
                         h4("Heterozygosity data not available"),
                         helpText("Enable 'Calculate heterozygosity' option in settings")
                       )
              ),
              tabPanel("Sample Size",
                       withSpinner(plotOutput("pop_sample_size_plot", height = "400px"), 
                                   type = 4, color = "#0dc5c1")
              )
            )
          )
        ),
        
        fluidRow(
          # Detailed Results Table
          box(
            title = "Detailed Results Table", 
            status = "warning", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            withSpinner(DTOutput("pop_allele_table"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            downloadButton("download_pop_table", "Download Table (CSV)", 
                           class = "btn-primary")
          )
        ),
        
        fluidRow(
          # Population Information
          box(
            title = "Population Information", 
            status = "info", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            h4("Available Populations in Data"),
            verbatimTextOutput("pop_info_summary"),
            
            br(),
            h4("Population Assignment"),
            helpText("To assign populations to individuals, you can:"),
            tags$ul(
              tags$li("Upload a metadata file with population information"),
              tags$li("Use the population assignment tools in the Data Processing tab"),
              tags$li("Manually assign populations using custom scripts")
            )
          )
        )
      ),
      
      
      # Population Correlations Tab
      tabItem(
        tabName = "pop_correlations",
        h2("Population Genetic Correlations"),
        
        fluidRow(
          box(width = 12, title = "Correlation Analysis", status = "primary", solidHeader = TRUE,
              fluidRow(
                column(6,
                       selectInput("corr_method", "Correlation Method",
                                   choices = c("Pearson" = "pearson",
                                               "Spearman" = "spearman",
                                               "Both" = "both"),
                                   selected = "both")
                ),
                column(6,
                       actionButton("run_corr_analysis", "Calculate Correlations", 
                                    icon = icon("calculator"), class = "btn-primary",
                                    style = "margin-top: 25px;")
                )
              )
          )
        ),
        
        fluidRow(
          box(width = 6, title = "Correlation Matrix", status = "success", solidHeader = TRUE,
              plotOutput("corr_heatmap", height = "500px")
          ),
          box(width = 6, title = "Correlation Statistics", status = "info", solidHeader = TRUE,
              verbatimTextOutput("corr_summary"),
              br(),
              downloadButton("download_corr_matrix", "Download Matrix")
          )
        ),
        
        fluidRow(
          box(width = 12, title = "Correlation Table", status = "primary", solidHeader = TRUE,
              DTOutput("corr_table")
          )
        ),
        fluidRow(
          column(
            width = 12,
            h2("🧬 Comprehensive Population Genetics Analysis", style = "margin-top: 0;"),
            helpText("Advanced genetic diversity and population structure analyses")
          )
        ),
        
        fluidRow(
          # Analysis Configuration
          box(
            title = "Analysis Configuration", 
            status = "primary", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            
            h4("Select Analyses to Run"),
            checkboxGroupInput("popgen_analyses", "Analyses to perform:",
                               choices = c(
                                 "Allele Frequency Barcode" = "allele_barcode",
                                 "Heterozygosity Statistics" = "heterozygosity",
                                 "Pairwise FST & NJ Tree" = "fst",
                                 "Hardy-Weinberg Equilibrium" = "hwe",
                                 "DAPC Assignment" = "dapc"
                               ),
                               selected = c("allele_barcode", "heterozygosity", "fst")),
            
            hr(),
            
            h4("Population Settings"),
            numericInput("popgen_max_markers", "Maximum markers for barcode",
                         value = 500, min = 100, max = 2000, step = 100,
                         width = "100%"),
            checkboxInput("popgen_find_clusters", "Automatically find clusters in DAPC", TRUE),
            numericInput("popgen_max_clusters", "Maximum clusters to test",
                         value = 20, min = 2, max = 50, step = 1,
                         width = "100%"),
            
            hr(),
            
            actionButton("run_popgen_analyses", "Run All Selected Analyses", 
                         icon = icon("play"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            
            br(),
            downloadButton("download_popgen_report", "Download Complete Report",
                           class = "btn-success btn-block")
          ),
          
          # Analysis Status
          box(
            title = "Analysis Status", 
            status = "success", 
            solidHeader = TRUE,
            width = 8,
            collapsible = TRUE,
            tabBox(
              width = 12,
              tabPanel("Status",
                       withSpinner(verbatimTextOutput("popgen_status"), 
                                   type = 4, color = "#0dc5c1")
              ),
              tabPanel("Quick Statistics",
                       withSpinner(uiOutput("popgen_quick_stats"), 
                                   type = 4, color = "#0dc5c1")
              ),
              tabPanel("Population Info",
                       withSpinner(verbatimTextOutput("popgen_pop_info"), 
                                   type = 4, color = "#0dc5c1")
              )
            )
          )
        ),
        
        # Results Tabs
        fluidRow(
          tabBox(
            width = 12,
            title = "Analysis Results",
            id = "popgen_results_tabs",
            
            # Allele Barcode Tab
            tabPanel("Allele Barcode",
                     fluidRow(
                       column(8,
                              withSpinner(plotOutput("popgen_barcode_plot", height = "600px"), 
                                          type = 4, color = "#0dc5c1")
                       ),
                       column(4,
                              h4("Barcode Settings"),
                              selectInput("barcode_color", "Color palette",
                                          choices = c("viridis", "magma", "plasma", "inferno", "cividis"),
                                          selected = "viridis"),
                              checkboxInput("barcode_cluster", "Cluster populations", TRUE),
                              checkboxInput("barcode_show_labels", "Show population labels", TRUE),
                              hr(),
                              h4("Download Options"),
                              downloadButton("download_barcode_plot", "Download Plot",
                                             class = "btn-success"),
                              br(), br(),
                              downloadButton("download_barcode_data", "Download Frequency Data",
                                             class = "btn-info")
                       )
                     )
            ),
            
            # Heterozygosity Tab
            tabPanel("Heterozygosity",
                     fluidRow(
                       column(8,
                              withSpinner(plotOutput("popgen_het_plot", height = "500px"), 
                                          type = 4, color = "#0dc5c1")
                       ),
                       column(4,
                              h4("Heterozygosity Statistics"),
                              withSpinner(DTOutput("popgen_het_table"), 
                                          type = 4, color = "#0dc5c1")
                       )
                     ),
                     fluidRow(
                       box(width = 12,
                           title = "Per Locus Statistics",
                           collapsible = TRUE,
                           collapsed = TRUE,
                           withSpinner(DTOutput("popgen_het_locus_table"), 
                                       type = 4, color = "#0dc5c1")
                       )
                     )
            ),
            
            # FST Analysis Tab
            tabPanel("FST Analysis",
                     fluidRow(
                       column(6,
                              h4("FST Matrix"),
                              withSpinner(plotOutput("popgen_fst_plot", height = "400px"), 
                                          type = 4, color = "#0dc5c1"),
                              br(),
                              withSpinner(DTOutput("popgen_fst_table"), 
                                          type = 4, color = "#0dc5c1")
                       ),
                       column(6,
                              h4("Neighbor-Joining Tree"),
                              conditionalPanel(
                                condition = "output.fst_tree_available == true",
                                withSpinner(plotOutput("popgen_nj_tree", height = "400px"), 
                                            type = 4, color = "#0dc5c1")
                              ),
                              conditionalPanel(
                                condition = "output.fst_tree_available == false",
                                h4("Tree not available"),
                                helpText("Insufficient data for tree construction")
                              )
                       )
                     )
            ),
            
            # HWE Testing Tab
            tabPanel("HWE Testing",
                     fluidRow(
                       column(8,
                              withSpinner(plotOutput("popgen_hwe_plot", height = "400px"), 
                                          type = 4, color = "#0dc5c1")
                       ),
                       column(4,
                              h4("HWE Summary"),
                              withSpinner(verbatimTextOutput("popgen_hwe_summary"), 
                                          type = 4, color = "#0dc5c1")
                       )
                     ),
                     fluidRow(
                       box(width = 12,
                           title = "Detailed HWE Results",
                           collapsible = TRUE,
                           collapsed = TRUE,
                           withSpinner(DTOutput("popgen_hwe_table"), 
                                       type = 4, color = "#0dc5c1")
                       )
                     )
            ),
            
            # DAPC Analysis Tab
            tabPanel("DAPC Analysis",
                     fluidRow(
                       column(6,
                              h4("DAPC Scatter Plot"),
                              withSpinner(plotOutput("popgen_dapc_scatter", height = "400px"), 
                                          type = 4, color = "#0dc5c1")
                       ),
                       column(6,
                              h4("Assignment Probabilities"),
                              withSpinner(plotOutput("popgen_dapc_assignment", height = "400px"), 
                                          type = 4, color = "#0dc5c1")
                       )
                     ),
                     fluidRow(
                       column(6,
                              h4("DAPC Scores"),
                              withSpinner(DTOutput("popgen_dapc_scores"), 
                                          type = 4, color = "#0dc5c1")
                       ),
                       column(6,
                              h4("Assignment Table"),
                              withSpinner(DTOutput("popgen_dapc_assignment_table"), 
                                          type = 4, color = "#0dc5c1")
                       )
                     )
            ),
            
            # Hybrid Detection Tab
            tabPanel("Hybrid Detection",
                     fluidRow(
                       column(4,
                              h4("Hybrid Analysis Settings"),
                              
                              # Population selection
                              selectInput("hybrid_pop1", "Parent Population 1",
                                          choices = NULL, selected = NULL),
                              selectInput("hybrid_pop2", "Parent Population 2",
                                          choices = NULL, selected = NULL),
                              selectInput("hybrid_admixed", "Putative Hybrid Population",
                                          choices = NULL, selected = NULL),
                              
                              # Analysis parameters
                              numericInput("hybrid_threshold", "Admixture Threshold",
                                           value = 0.1, min = 0.01, max = 0.5, step = 0.01),
                              
                              br(),
                              
                              # Run button
                              actionButton("run_hybrid_detection", "Run Hybrid Analysis",
                                           icon = icon("search"), 
                                           class = "btn-warning btn-block"),
                              
                              hr(),
                              
                              # Results summary
                              h4("Analysis Summary"),
                              verbatimTextOutput("hybrid_results_summary"),
                              
                              # Download UI (conditional)
                              uiOutput("hybrid_download_ui")
                       ),
                       
                       column(8,
                              h4("Hybrid Pattern Visualization"),
                              
                              # Main plot
                              withSpinner(plotOutput("hybrid_pattern_plot", height = "500px"), 
                                          type = 4, color = "#0dc5c1"),
                              
                              # Plot controls
                              fluidRow(
                                column(4,
                                       selectInput("hybrid_plot_type", "Plot Type",
                                                   choices = c("Structure Plot" = "structure",
                                                               "Triangle Plot" = "triangle"),
                                                   selected = "structure")
                                ),
                                column(4,
                                       checkboxInput("hybrid_show_labels", "Show Sample Labels", TRUE)
                                ),
                                column(4,
                                       downloadButton("download_hybrid_plot", "Download Plot",
                                                      class = "btn-primary")
                                )
                              )
                       )
                     )
            )
          )
        ),
        
        # Data Export Section
        fluidRow(
          box(
            title = "Data Export", 
            status = "info", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            h4("Export Complete Analysis Results"),
            fluidRow(
              column(3,
                     downloadButton("export_popgen_rdata", "Export as RData",
                                    class = "btn-primary")
              ),
              column(3,
                     downloadButton("export_popgen_excel", "Export to Excel",
                                    class = "btn-success")
              ),
              column(3,
                     downloadButton("export_popgen_json", "Export as JSON",
                                    class = "btn-info")
              ),
              column(3,
                     actionButton("generate_popgen_report", "Generate HTML Report",
                                  icon = icon("file-alt"), class = "btn-warning")
              )
            )
          )
        )
      ),
      
      # Fst-Fis Analysis Tab
      tabItem(
        tabName = "fst_fis_analysis",
        h2("Fst-Fis Analysis (Waples Plot)"),
        
        fluidRow(
          box(width = 12, title = "Analysis Settings", status = "primary", solidHeader = TRUE,
              fluidRow(
                column(3,
                       sliderInput("fst_fis_bin_size", "Bin Size for Fst", 
                                   min = 0.01, max = 0.1, value = 0.05, step = 0.01)
                ),
                column(3,
                       selectInput("fis_column", "Fis Column",
                                   choices = c("Auto-detect" = "auto"),
                                   selected = "auto")
                ),
                column(3,
                       selectInput("fst_column", "Fst Column",
                                   choices = c("Auto-detect" = "auto"),
                                   selected = "auto")
                ),
                column(3,
                       actionButton("run_waples_analysis", "Generate Waples Plot", 
                                    icon = icon("chart-line"), class = "btn-primary",
                                    style = "margin-top: 25px;")
                )
              )
          )
        ),
        
        fluidRow(
          box(width = 12, title = "Waples Plot (Fst vs Fis)", status = "success", solidHeader = TRUE,
              plotOutput("waples_plot", height = "600px"),
              br(),
              downloadButton("download_waples_plot", "Download Plot")
          )
        ),
        
        fluidRow(
          box(width = 12, title = "Binned Statistics", status = "info", solidHeader = TRUE,
              DTOutput("waples_stats_table")
          )
        )
      ),
      
      tabItem(
        tabName = "lea_analysis",
        fluidRow(
          column(
            width = 12,
            h2("🧬 LEA Population Structure Analysis", style = "margin-top: 0;")
          )
        ),
        
        fluidRow(
          # LEA Settings Box
          box(
            title = "LEA Analysis Settings", 
            status = "primary", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            
            sliderInput("lea_min_k", "Minimum K value", 
                        min = 1, max = 10, value = 2, step = 1,
                        width = "100%"),
            sliderInput("lea_max_k", "Maximum K value", 
                        min = 2, max = 15, value = 6, step = 1,
                        width = "100%"),
            numericInput("lea_replicates", "Number of replicates per K",
                         value = 10, min = 1, max = 50, step = 1,
                         width = "100%"),
            numericInput("lea_alpha", "Alpha parameter for sNMF",
                         value = 100, min = 1, max = 1000, step = 10,
                         width = "100%"),
            checkboxInput("lea_calculate_entropy", "Calculate cross-entropy", TRUE),
            checkboxInput("lea_export_data", "Export intermediate files", FALSE),
            
            br(),
            actionButton("run_lea", "Run LEA Analysis", 
                         icon = icon("project-diagram"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            helpText("Performs population structure analysis using sNMF algorithm")
          ),
          
          # Results Summary Box
          box(
            title = "Analysis Status", 
            status = "success", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            withSpinner(verbatimTextOutput("lea_status"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            uiOutput("lea_export_buttons")
          ),
          
          # Cross-Entropy Box
          box(
            title = "Optimal K Selection", 
            status = "info", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            withSpinner(plotOutput("lea_cross_entropy_plot", height = "300px"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            uiOutput("lea_optimal_k")
          )
        ),
        
        fluidRow(
          # Structure Plots
          box(
            title = "Population Structure Plots", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            
            selectInput("lea_plot_k", "Select K value to display",
                        choices = NULL, selected = NULL,
                        width = "200px"),
            
            withSpinner(plotOutput("lea_structure_plot", height = "500px"), 
                        type = 4, color = "#0dc5c1"),
            
            br(),
            fluidRow(
              column(3,
                     downloadButton("download_lea_plot", "Download Plot", 
                                    class = "btn-success")
              ),
              column(3,
                     actionButton("lea_prev_k", "Previous K", 
                                  icon = icon("arrow-left"), 
                                  class = "btn-info")
              ),
              column(3,
                     actionButton("lea_next_k", "Next K", 
                                  icon = icon("arrow-right"), 
                                  class = "btn-info")
              ),
              column(3,
                     actionButton("lea_compare_all", "Compare All K", 
                                  icon = icon("chart-bar"), 
                                  class = "btn-warning")
              )
            )
          )
        ),
        
        fluidRow(
          # Q Matrix Viewer
          box(
            title = "Q Matrix Data", 
            status = "info", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            
            tabsetPanel(
              tabPanel("Q Matrix Table",
                       withSpinner(DTOutput("lea_q_matrix_table"), 
                                   type = 4, color = "#0dc5c1")
              ),
              tabPanel("Admixture Summary",
                       withSpinner(plotOutput("lea_admixture_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")
              ),
              tabPanel("Export Results",
                       h4("Export LEA Results"),
                       downloadButton("download_lea_csv", "Download Q Matrices (CSV)", 
                                      class = "btn-primary"),
                       br(), br(),
                       downloadButton("download_lea_rdata", "Download Full Results (RData)", 
                                      class = "btn-success")
              )
            )
          )
        )
      ),
      
      
      # GWAS Analysis Tab
      tabItem(
        tabName = "gwas_analysis",
        fluidRow(
          column(
            width = 12,
            h2("📊 GWAS Analysis", style = "margin-top: 0;"),
            helpText("Genome-Wide Association Studies for individual traits")
          )
        ),
        
        fluidRow(
          # GWAS Settings
          box(
            title = "GWAS Analysis Settings", 
            status = "primary", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            uiOutput("gwas_trait_selector"),
            sliderInput("maf_threshold_gwas", "MAF Threshold for GWAS",
                        min = 0, max = 0.5, value = 0.05, step = 0.01,
                        width = "100%"),
            sliderInput("missing_threshold_gwas", "Missing Data Threshold",
                        min = 0, max = 0.5, value = 0.2, step = 0.05,
                        width = "100%"),
            selectInput("gwas_model", "GWAS Model",
                        choices = c("Simple Linear Model" = "LM",
                                    "Linear Mixed Model" = "LMM"),
                        selected = "LM",
                        width = "100%"),
            br(),
            actionButton("run_gwas", "Run GWAS Analysis", 
                         icon = icon("play"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            helpText("Run GWAS for selected trait")
          ),
          
          # GWAS Results Summary
          box(
            title = "GWAS Results Summary", 
            status = "success", 
            solidHeader = TRUE,
            width = 8,
            collapsible = TRUE,
            withSpinner(verbatimTextOutput("gwas_results_summary"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            fluidRow(
              column(4,
                     actionButton("export_gwas", "Export GWAS Results", 
                                  icon = icon("download"), 
                                  class = "btn-success",
                                  width = "100%")
              ),
              column(4,
                     actionButton("view_top_snps", "View Top SNPs", 
                                  icon = icon("search"), 
                                  class = "btn-info",
                                  width = "100%")
              ),
              column(4,
                     downloadButton("download_gwas_csv", "Download CSV", 
                                    class = "btn-primary",
                                    width = "100%")
              )
            )
          )
        ),
        
        
        fluidRow(
          # GWAS Results Table
          box(
            title = "GWAS Results", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tabsetPanel(
              tabPanel("All Results", 
                       withSpinner(DTOutput("gwas_results_table"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Significant SNPs", 
                       withSpinner(DTOutput("significant_snps_table"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Manhattan Plot", 
                       fluidRow(
                         column(3,
                                numericInput("manhattan_threshold", 
                                             "Significance Threshold", 
                                             value = 5e-8, min = 1e-20, max = 1,
                                             step = 1e-8, width = "100%")
                         ),
                         column(3,
                                selectInput("manhattan_chromosomes", 
                                            "Chromosomes to Show",
                                            choices = c("All", "1-10", "11-20", "Custom"),
                                            selected = "All", width = "100%")
                         ),
                         column(3,
                                checkboxInput("manhattan_labels", 
                                              "Label Top SNPs", TRUE)
                         ),
                         column(3,
                                numericInput("manhattan_top_n", 
                                             "Top N to Label", 
                                             value = 10, min = 1, max = 50,
                                             width = "100%")
                         )
                       ),
                       withSpinner(plotOutput("manhattan_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")
              ),
              tabPanel("QQ Plot",
                       withSpinner(plotOutput("qq_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Effect Size Distribution",
                       withSpinner(plotOutput("effect_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              
              # In your UI's tabsetPanel, add:
              tabPanel("QTL Analysis",
                       sidebarLayout(
                         sidebarPanel(
                           uiOutput("qtl_trait_selector"),
                           uiOutput("qtl_controls"),
                           br(),
                           uiOutput("qtl_analysis_status"),
                           br(),
                           downloadButton("download_qtl_results", "Download Results", 
                                          class = "btn-success")
                         ),
                         mainPanel(
                           tabsetPanel(
                             tabPanel("Results Summary",
                                      verbatimTextOutput("qtl_results_summary"),
                                      DTOutput("qtl_results_table")
                             ),
                             tabPanel("LOD Profile",
                                      plotOutput("qtl_plot", height = "600px")
                             ),
                             tabPanel("QTL Peaks",
                                      DTOutput("qtl_peaks_table"),
                                      plotOutput("qtl_effect_plot", height = "400px")
                             ),
                             tabPanel("Diagnostics",
                                      plotOutput("qtl_diagnostics_plot", height = "500px")
                             )
                           )
                         )
                       )
              )
            )
          )
        )
      ),
      
      
      # ============================================================
      # MULTI-TRAIT GWAS TAB
      # ============================================================
      tabItem(tabName = "multi_trait_gwas",
              h2("Multi-Trait GWAS Analysis"),
              
              fluidRow(
                box(width = 4, title = "Analysis Settings", status = "primary", solidHeader = TRUE,
                    uiOutput("multi_trait_selector"),
                    
                    br(),
                    actionButton("run_multi_trait", "Run Multi-Trait Analysis",
                                 icon = icon("play"), 
                                 class = "btn-primary btn-block"),
                    
                    hr(),
                    
                    h4("Analysis Options"),
                    checkboxInput("include_correlations", "Calculate trait correlations", TRUE),
                    checkboxInput("calculate_pleiotropy", "Identify pleiotropic SNPs", TRUE),
                    selectInput("multi_trait_method", "Analysis Method",
                                choices = c("Combined P-values", "Meta-analysis", "Multi-trait BLUP"),
                                selected = "Combined P-values")
                ),
                
                box(width = 8, title = "Analysis Results", status = "success", solidHeader = TRUE,
                    tabBox(width = 12,
                           tabPanel("Summary",
                                    verbatimTextOutput("multi_trait_summary")
                           ),
                           tabPanel("Correlations",
                                    plotOutput("trait_cor_plot", height = "500px")
                           ),
                           tabPanel("Pleiotropic SNPs",
                                    DTOutput("pleiotropy_table")
                           ),
                           tabPanel("Combined Results",
                                    DTOutput("multi_trait_gwas_table")
                           ),
                           tabPanel("Statistics",
                                    verbatimTextOutput("multi_trait_stats")
                           )
                    )
                )
              ),
              
              fluidRow(
                box(width = 4, title = "Add Left-Right Mrker", status = "success", solidHeader = TRUE,
                    actionButton("run_cross_snp_all_traits", "Add Flanking Markers", class = "btn-primary")
                )
              )
      ),
      
      # ============================================================
      # UNIFIED RESULTS TAB
      # ============================================================
      tabItem(tabName = "unified_results",
              h2("Unified Multi-Trait Results"),
              
              fluidRow(
                box(width = 12, title = "Data Status", status = "info", solidHeader = TRUE,
                    uiOutput("unified_data_status_ui")
                )
              ),
              
              fluidRow(
                box(width = 12, title = "Unified Data Table", status = "primary", solidHeader = TRUE,
                    DTOutput("unified_data_table"),
                    br(),
                    downloadButton("download_unified_data", "Download Unified Data",
                                   class = "btn-success")
                )
              ),
              
              fluidRow(
                box(width = 6, title = "Data Summary", status = "info", solidHeader = TRUE,
                    verbatimTextOutput("unified_data_summary")
                ),
                
                box(width = 6, title = "Quick Actions", status = "warning", solidHeader = TRUE,
                    h4("Data Operations"),
                    actionButton("create_unified_format", "Create Unified Format",
                                 icon = icon("table"), class = "btn-info btn-block"),
                    br(),
                    actionButton("export_unified_excel", "Export to Excel",
                                 icon = icon("file-excel"), class = "btn-success btn-block"),
                    br(),
                    hr(),
                    h4("Visualization"),
                    selectInput("unified_plot_type", "Select Plot Type",
                                choices = c("Manhattan Plot" = "manhattan",
                                            "Volcano Plot" = "volcano",
                                            "Effect Distribution" = "effect_dist",
                                            "QQ Plot" = "qq")),
                    actionButton("plot_unified_data", "Generate Plot",
                                 icon = icon("chart-bar"), class = "btn-primary btn-block")
                )
              )
      ),
      
      # ============================================================
      # ENHANCED ANALYSIS TAB
      # ============================================================
      tabItem(tabName = "enhanced_analysis",
              h2("Enhanced Multi-Trait Analysis"),
              
              fluidRow(
                box(width = 12, title = "Enhanced Data Formats", status = "primary", solidHeader = TRUE,
                    tabBox(width = 12,
                           tabPanel("Enhanced Unified Format",
                                    DTOutput("enhanced_unified_table"),
                                    br(),
                                    downloadButton("download_enhanced_unified_data", 
                                                   "Download Enhanced Data",
                                                   class = "btn-success")
                           ),
                           tabPanel("Complete Integrated Format",
                                    DTOutput("complete_integrated_table"),
                                    br(),
                                    downloadButton("download_complete_integrated_data", 
                                                   "Download Complete Data",
                                                   class = "btn-success")
                           ),
                           tabPanel("Analysis Dashboard",
                                    fluidRow(
                                      column(6,
                                             plotOutput("enhanced_manhattan_plot", height = "400px")
                                      ),
                                      column(6,
                                             plotOutput("enhanced_volcano_plot", height = "400px")
                                      )
                                    ),
                                    fluidRow(
                                      column(6,
                                             plotOutput("enhanced_effect_plot", height = "400px")
                                      ),
                                      column(6,
                                             plotOutput("enhanced_pleiotropy_plot", height = "400px")
                                      )
                                    )
                           )
                    )
                )
              )
      ),
      
      # ============================================================
      # CROSS-SNP ANALYSIS TAB
      # ============================================================
      tabItem(tabName = "cross_snp_mapping",
              h2("Cross-SNP Analysis"),
              
              fluidRow(
                box(width = 12, title = "Analysis Setup", status = "primary", solidHeader = TRUE,
                    fluidRow(
                      column(6,
                             h4("Data Compatibility Check"),
                             verbatimTextOutput("data_compatibility_check")
                      ),
                      column(6,
                             h4("Analysis Controls"),
                             br(),
                             actionButton("run_cross_snp_analysis", "Run Cross-SNP Analysis",
                                          icon = icon("link"), 
                                          class = "btn-primary btn-block"),
                             br(),
                             actionButton("clear_cross_snp_results", "Clear Results",
                                          icon = icon("broom"),
                                          class = "btn-warning btn-block"),
                             br(),
                             downloadButton("download_cross_snp_data", "Download All Results",
                                            class = "btn-success btn-block")
                      )
                    )
                )
              ),
              
              fluidRow(
                box(width = 12, title = "Analysis Status", status = "info", solidHeader = TRUE,
                    uiOutput("cross_snp_status_ui")
                )
              ),
              
              fluidRow(
                box(width = 12, title = "Analysis Settings", status = "warning", solidHeader = TRUE,
                    fluidRow(
                      column(4,
                             h4("Data Selection"),
                             selectInput("cross_snp_gwas_data", "GWAS Data Source",
                                         choices = c("Unified Format" = "unified",
                                                     "Enhanced Format" = "enhanced",
                                                     "Complete Integrated" = "complete"),
                                         selected = "enhanced"),
                             checkboxInput("cross_snp_use_polymorphic", 
                                           "Use only polymorphic SNPs", TRUE),
                             checkboxInput("cross_snp_require_genotypes",
                                           "Require parent genotype data", FALSE)
                      ),
                      column(4,
                             h4("Filtering Options"),
                             numericInput("cross_snp_p_threshold", "P-value threshold",
                                          value = 0.05, min = 0.0001, max = 0.5, step = 0.01),
                             numericInput("cross_snp_lod_threshold", "LOD threshold",
                                          value = 2.0, min = 0, max = 10, step = 0.5),
                             selectInput("cross_snp_trait_filter", "Trait filter",
                                         choices = c("All traits" = "all",
                                                     "Diallel traits only" = "diallel",
                                                     "Selected traits" = "selected"),
                                         selected = "diallel")
                      ),
                      column(4,
                             h4("Output Options"),
                             checkboxInput("cross_snp_calc_expected", 
                                           "Calculate expected effects", TRUE),
                             checkboxInput("cross_snp_identify_top",
                                           "Identify top crosses/SNPs", TRUE),
                             numericInput("cross_snp_top_n", "Number of top items",
                                          value = 10, min = 5, max = 50, step = 5)
                      )
                    )
                )
              ),
              
              fluidRow(
                box(width = 6, title = "Analysis Settings", status = "primary", solidHeader = TRUE,
                    # Analysis controls
                    wellPanel(
                      style = "background-color: #f8f9fa;",
                      
                      h4(icon("link"), " Data Requirements", style = "margin-top: 0;"),
                      
                      div(class = "alert alert-info",
                          tags$small(
                            HTML("
                                                <b>Required data:</b>
                                                   <ul>
                                                      <li>✓ Diallel data with crosses</li>
                                                      <li>✓ Multi-trait GWAS results</li>
                                                      <li>✓ Genotype matrix (optional, for parent genotypes)</li>
                                                   </ul>"
                            )
                          )
                      ),
                      
                      br(),
                      
                      h5(icon("cogs"), " Analysis Settings"),
                      
                      checkboxInput("include_genotypes", 
                                    "Include Parent Genotypes", 
                                    value = TRUE),
                      
                      checkboxInput("only_significant", 
                                    "Only Significant SNPs (p < 0.05)", 
                                    value = FALSE),
                      
                      selectInput("trait_filter_cross", 
                                  "Filter by Trait (optional):",
                                  choices = NULL,
                                  multiple = TRUE),
                      
                      br(),
                      
                      actionButton("run_cross_snp_analysis", 
                                   "Run Cross-SNP Analysis",
                                   icon = icon("project-diagram"),
                                   class = "btn-primary btn-block"),
                      
                      br(),
                      
                      actionButton("clear_cross_snp", 
                                   "Clear Results",
                                   icon = icon("broom"),
                                   class = "btn-warning btn-block")
                    )
                ),
                
                box(width = 6, title = "Export Settings", status = "primary", solidHeader = TRUE,
                    # Export controls
                    wellPanel(
                      style = "background-color: #f8f9fa;",
                      
                      h4(icon("download"), " Export Results", style = "margin-top: 0;"),
                      
                      downloadButton("download_cross_snp_data", 
                                     "Download All Data (ZIP)",
                                     class = "btn-success btn-block",
                                     style = "margin-bottom: 10px;"),
                      
                      downloadButton("download_cross_snp_csv", 
                                     "Download Detailed CSV",
                                     class = "btn-info btn-block")
                    ),
                    
                    # Status indicator
                    uiOutput("cross_snp_status")
                )
              ),
              
              fluidRow(
                box(width = 12, title = "Interactive Visualization", status = "info", solidHeader = TRUE,
                    tabBox(width = 12,
                           # Tabs for different views
                           tabPanel("Analysis Overview",
                                    # Summary and controls
                                    verbatimTextOutput("cross_snp_summary")
                           ),
                           
                           tabPanel("Detailed Data",
                                    br(),
                                    div(class = "data-box",
                                        DTOutput("cross_snp_table")
                                    )
                           ),
                           
                           # Cross summary tab
                           tabPanel("Cross Summary",
                                    br(),
                                    fluidRow(
                                      column(6,
                                             div(class = "data-box",
                                                 h5("Cross Statistics"),
                                                 DTOutput("cross_summary_table")
                                             )
                                      ),
                                      column(6,
                                             div(class = "data-box",
                                                 h5("Top Crosses by LOD"),
                                                 DTOutput("top_crosses_table")
                                             )
                                      )
                                    )
                           ),
                           
                           # SNP summary tab
                           tabPanel("SNP Summary",
                                    br(),
                                    div(class = "data-box",
                                        h5("SNP Performance Across Crosses"),
                                        DTOutput("snp_summary_cross_table")
                                    )
                           ),
                           
                           tabPanel("Top Crosses",
                                    br(),
                                    div(class = "data-box",
                                        h5("Top Crosses-SNP Table"),
                                    DTOutput("top_three_crosses_table")
                                    )
                           ),
                           
                           # Visualizations tab
                           tabPanel("Visualizations",
                                    br(),
                                    fluidRow(
                                      column(6,
                                             div(class = "data-box",
                                                 h5("Top Crosses by LOD Score"),
                                                 plotOutput("cross_lod_plot", height = "400px")
                                             )
                                      ),
                                      column(6,
                                             div(class = "data-box",
                                                 h5("SNP Effect Heatmap"),
                                                 plotlyOutput("snp_effect_heatmap", height = "400px")
                                             )
                                      )
                                    )
                           ),
                           
                           # Search/Filter tab
                           tabPanel("Search & Filter",
                                    br(),
                                    fluidRow(
                                      column(4,
                                             wellPanel(
                                               h5(icon("search"), " Search Crosses"),
                                               selectInput("search_cross", "Select Cross:",
                                                           choices = NULL),
                                               actionButton("view_cross_details", 
                                                            "View Details",
                                                            class = "btn-primary")
                                             )
                                      ),
                                      column(4,
                                             wellPanel(
                                               h5(icon("dna"), " Search SNPs"),
                                               selectInput("search_snp_cross", "Select SNP:",
                                                           choices = NULL),
                                               actionButton("view_snp_cross_details", 
                                                            "View Details",
                                                            class = "btn-primary")
                                             )
                                      ),
                                      column(4,
                                             wellPanel(
                                               h5(icon("filter"), " Advanced Filter"),
                                               numericInput("lod_threshold_filter", 
                                                            "LOD Threshold:", 
                                                            value = 2, min = 0, max = 10),
                                               numericInput("pve_threshold_filter", 
                                                            "PVE Threshold (%):", 
                                                            value = 5, min = 0, max = 100),
                                               actionButton("apply_filters", 
                                                            "Apply Filters",
                                                            class = "btn-primary")
                                             )
                                      )
                                    ),
                                    
                                    # Filtered results
                                    div(class = "data-box",
                                        h5("Filtered Results"),
                                        DTOutput("filtered_cross_snp_table")
                                    )
                           )
                    )
                )
              ),
              
              fluidRow(
                box(width = 12, title = "Data Export", status = "success", solidHeader = TRUE,
                    fluidRow(
                      column(3,
                             downloadButton("dl_cross_snp_detailed", "Detailed Data",
                                            class = "btn-primary btn-block")
                      ),
                      column(3,
                             downloadButton("dl_cross_summary", "Cross Summary",
                                            class = "btn-info btn-block")
                      ),
                      column(3,
                             downloadButton("dl_snp_summary", "SNP Summary",
                                            class = "btn-info btn-block")
                      ),
                      column(3,
                             downloadButton("dl_cross_snp_all", "All Results (ZIP)",
                                            class = "btn-success btn-block")
                      )
                    )
                )
              )
      ),
      
      
      # Diallel Analysis Tab
      tabItem(tabName = "diallel_analysis",
              fluidRow(
                box(title = "Diallel Analysis Settings", 
                    status = "primary", 
                    solidHeader = TRUE,
                    width = 3,
                    collapsible = TRUE,
                    selectInput("diallel_trait", "Select Single Trait for Analysis", 
                                choices = NULL, selected = ""),
                    actionButton("run_diallel", "Run Single Trait Analysis",  
                                 icon = icon("play"), class = "btn-success"),
                    br(), br(),
                    actionButton("run_all_diallel_traits", "Run All Traits Analysis", 
                                 icon = icon("play-circle"), class = "btn-warning"),
                    helpText("Analyze all numeric traits in diallel data"),
                    br(),
                    downloadButton("download_diallel_summary", "Download Summary",
                                   class = "btn-info btn-block")
                ),
                box(title = "ANOVA Results", status = "info", solidHeader = TRUE,
                    verbatimTextOutput("diallel_anova_results"), width = 5),
                box(title = "Genetic Parameters", status = "info", solidHeader = TRUE,
                    verbatimTextOutput("diallel_genetic_params"), width = 4)
              ),
              
              fluidRow(
                box(title = "Dataset Overview", status = "primary", solidHeader = TRUE, width = 12,
                    DTOutput("diallel_data_table")
                ),
                fluidRow(
                  box(title = "Data Summary", status = "info", solidHeader = TRUE,
                      verbatimTextOutput("diallel_data_summary"), width = 6),
                  box(title = "Trait Descriptions", status = "info", solidHeader = TRUE,
                      tableOutput("diallel_trait_descriptions"), width = 6)
                )
              ),
              fluidRow(
                box(
                  title = "Pivoted Diallel Data (Parents Combined)", 
                  status = "info", 
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  collapsed = TRUE,
                  helpText("Parent1 and Parent2 are combined in the 'Parent' column"),
                  withSpinner(DTOutput("diallel_long_table"), 
                              type = 4, color = "#0dc5c1")
                )
              ),
              
              fluidRow(
                box(
                  title = "Parent Summary", 
                  status = "success", 
                  solidHeader = TRUE,
                  width = 12,
                  withSpinner(DTOutput("parent_summary_table"), 
                              type = 4, color = "#0dc5c1")
                )
              ),
              
              fluidRow(
                box("ANOVA Summary", solidHeader = TRUE,
                    DTOutput("diallel_all_anova"), width = 6),
                box(title = "General Combining Ability (GCA) Effects", 
                    status = "success", solidHeader = TRUE,
                    DTOutput("diallel_gca_table"), width = 6)
              ),
              
              fluidRow(
                box(title = "Specific Combining Ability (SCA) Effects", 
                    status = "warning", solidHeader = TRUE,
                    DTOutput("diallel_sca_table"), width = 6)
              ),
              
              fluidRow(
                box(title = "GCA Effects Plot", status = "primary", solidHeader = TRUE,
                    plotlyOutput("diallel_gca_plot"), width = 6),
                box(title = "SCA Effects Heatmap", status = "primary", solidHeader = TRUE,
                    plotlyOutput("diallel_sca_heatmap"), width = 6)
              ),
              fluidRow(
                box(title = "Mean Performance vs GCA", status = "info", solidHeader = TRUE,
                    plotlyOutput("diallel_mean_vs_gca"), width = 6),
                box("Traits Correlation", solidHeader = TRUE,
                    plotlyOutput("diallel_traits_cor_plot"), width = 6)
              ),
              
              
              fluidRow(
                column(
                  width = 12,
                  h2("🌐 Parent Relationship Network", style = "margin-top: 0;")
                )
              ),
              
              fluidRow(
                box(
                  title = "Network Settings", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 3,
                  collapsible = TRUE,
                  numericInput("node_size", "Node Size", 
                               value = 30, min = 10, max = 100, step = 5,
                               width = "100%"),
                  selectInput("node_shape", "Node Shape",
                              choices = c("dot", "square", "triangle", "diamond", "star"),
                              selected = "dot", width = "100%"),
                  selectInput("layout", "Network Layout",
                              choices = c("layout_nicely", "layout_with_fr", 
                                          "layout_with_kk", "layout_on_grid"),
                              selected = "layout_nicely", width = "100%"),
                  checkboxInput("show_labels", "Show Node Labels", TRUE),
                  checkboxInput("physics", "Enable Physics", TRUE),
                  br(),
                  actionButton("generate_network", "Generate Network", 
                               icon = icon("project-diagram"), 
                               class = "btn-primary btn-block",
                               width = "100%")
                ),
                
                box(
                  title = "Parent Relationship Network", 
                  status = "success", 
                  solidHeader = TRUE,
                  width = 9,
                  collapsible = TRUE,
                  height = "700px",
                  withSpinner(visNetworkOutput("parent_network", height = "650px"), 
                              type = 4, color = "#0dc5c1"),
                  br(),
                  downloadButton("download_network", "Export Network as HTML",
                                 class = "btn-info")
                )
              ),
              
              fluidRow(
                box(
                  title = "Network Statistics", 
                  status = "info", 
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  verbatimTextOutput("network_stats")
                )
              ),
              
              
      ),
      
      
      # Heterosis Analysis Tab
      tabItem(tabName = "diallel_heterosis",
              fluidRow(
                box(title = "Heterosis Controls", status = "primary", solidHeader = TRUE,
                    selectInput("heterosis_type", "Heterosis Type:",
                                choices = c("Mid-Parent Heterosis" = "mph",
                                            "Better-Parent Heterosis" = "bph"),
                                selected = "mph"),
                    sliderInput("heterosis_threshold", "Heterosis Threshold (%):",
                                min = 0, max = 100, value = 15, step = 5),
                    width = 12)
              ),
              
              fluidRow(
                box("Heterosis Summary",
                    DTOutput("diallel_all_heterosis"), width = 12)
              ),
              
              fluidRow(
                box(title = "Heterosis Values", status = "info", solidHeader = TRUE,
                    DTOutput("diallel_heterosis_table"), width = 6),
                box(title = "Heterosis Statistics", status = "info", solidHeader = TRUE,
                    verbatimTextOutput("diallel_heterosis_stats"), width = 6)
              ),
              
              fluidRow(
                box(
                  title = "Analysis Controls", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 4,
                  actionButton("analyze_top_heterotic_crosses", "Analyze Top 3 Heterotic Crosses",
                               icon = icon("trophy"), class = "btn-success btn-block"),
                  br(),
                  helpText("Identifies and analyzes the top 3 crosses with highest Mid-Parent Heterosis (MPH)"),
                  br(),
                  uiOutput("top_crosses_info_box")
                ),
                
                box(
                  title = "Results", 
                  status = "success", 
                  solidHeader = TRUE,
                  width = 8,
                  h2("Top Heterotic Crosses Analysis"),
                  uiOutput("top_heterotic_crosses_ui")
                )
              ),
              
              fluidRow(
                box(title = "Heterosis Distribution", status = "primary", solidHeader = TRUE,
                    plotlyOutput("diallel_heterosis_hist"), width = 6),
                box(title = "Top Heterotic Crosses", status = "success", solidHeader = TRUE,
                    plotlyOutput("diallel_heterosis_bar"), width = 6)
              ),
              fluidRow(
                box(title = "Heterosis vs SCA", status = "info", solidHeader = TRUE,
                    plotlyOutput("diallel_heterosis_vs_sca"), width = 6)
              )
      ),
      
      # Variety Selection Tab
      tabItem(tabName = "diallel_selection",
              fluidRow(
                valueBoxOutput("diallel_best_parent", width = 3),
                valueBoxOutput("diallel_best_cross", width = 3),
                valueBoxOutput("diallel_best_heterosis", width = 3),
                valueBoxOutput("diallel_heritability", width = 3)
              ),
              fluidRow(
                box(title = "Top Performing Parents", status = "success", solidHeader = TRUE,
                    DTOutput("diallel_top_parents"), width = 4),
                box(title = "Top Performing Crosses", status = "warning", solidHeader = TRUE,
                    DTOutput("diallel_top_crosses"), width = 4),
                box(title = "Top Heterotic Crosses", status = "danger", solidHeader = TRUE,
                    DTOutput("diallel_top_heterotic"), width = 4)
              ),
              fluidRow(
                box(title = "Selection Recommendations", status = "primary", solidHeader = TRUE,
                    htmlOutput("diallel_recommendations"), width = 12)
              )
      ),
      
      
      # UI for the actual genotype matching tab
      tabItem(
        tabName = "parent_snp_actual",
        fluidRow(
          column(
            width = 12,
            h2("🧬 Parent-SNP Matching (Actual Genotypes)", style = "margin-top: 0;"),
            helpText("Extracts actual genotype data for diallel parents from the genotype matrix")
          )
        ),
        
        fluidRow(
          box(
            title = "Analysis Controls", 
            status = "primary", 
            solidHeader = TRUE,
            width = 3,
            collapsible = TRUE,
            helpText("Click to match diallel parents with their actual genotypes from GWAS SNPs"),
            br(),
            actionButton("run_parent_snp_matching", "Run Parent-SNP Matching", 
                         icon = icon("dna"), 
                         class = "btn-primary btn-block",
                         width = "100%"),
            br(),
            actionButton("clear_parent_snp_results", "Clear Results", 
                         icon = icon("trash"), 
                         class = "btn-warning btn-block",
                         width = "100%"),
            br(),
            hr(),
            h5("Analysis Settings:"),
            numericInput("snp_p_threshold", "P-value Threshold for Significant SNPs",
                         value = 0.05, min = 0.0001, max = 1, step = 0.01,
                         width = "100%"),
            checkboxInput("include_missing_genotypes", "Include Missing Genotypes in Analysis", 
                          value = FALSE),
            br(),
            helpText("Note: This matches parent names from diallel data to genotype matrix rows.")
          ),
          
          box(
            title = "Data Status Check", 
            status = "info", 
            solidHeader = TRUE,
            width = 9,
            collapsible = TRUE,
            verbatimTextOutput("parent_snp_data_check"),
            br(),
            uiOutput("parent_snp_data_status_ui")
          )
        ),
        
        # ADD THIS SECTION: Results Display Tabs
        fluidRow(
          box(
            title = "Analysis Results", 
            status = "success", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            tabsetPanel(
              tabPanel(
                "Summary Statistics",
                verbatimTextOutput("parent_snp_summary_stats"),
                br(),
                h4("Parent Summary Table"),
                withSpinner(DTOutput("parent_snp_summary_table"), 
                            type = 4, color = "#0dc5c1")
              ),
              tabPanel(
                "Detailed Genotypes",
                h4("All Parent-SNP Combinations"),
                helpText("Showing all", tags$b("14,812"), "parent-SNP combinations"),
                withSpinner(DTOutput("parent_snp_detailed_table"), 
                            type = 4, color = "#0dc5c1"),
                br(),
                downloadButton("download_detailed_genotypes", "Download Detailed Data (CSV)",
                               class = "btn-info")
              ),
              tabPanel(
                "Significant SNP Patterns",
                h4("Genotype Patterns at Significant SNPs"),
                helpText("Patterns at SNPs with p <", textOutput("current_p_threshold", inline = TRUE)),
                withSpinner(DTOutput("parent_snp_patterns_table"), 
                            type = 4, color = "#0dc5c1"),
                br(),
                downloadButton("download_snp_patterns", "Download SNP Patterns (CSV)",
                               class = "btn-info")
              ),
              tabPanel(
                "Visualizations",
                h4("Parent Genetic Profiles"),
                fluidRow(
                  column(6,
                         plotOutput("parent_genotype_heatmap", height = "500px")
                  ),
                  column(6,
                         plotOutput("parent_snp_scatter", height = "500px")
                  )
                ),
                fluidRow(
                  column(6,
                         plotOutput("parent_performance_plot", height = "400px")
                  ),
                  column(6,
                         plotOutput("snp_effect_distribution", height = "400px")
                  )
                )
              ),
              # Add a new tab in the parent_snp_actual tabsetPanel:
              tabPanel(
                "All Traits Analysis",
                h4("Performance Across All Identified Traits"),
                helpText("Showing all traits identified in values$parent_snp_actual$summary_table"),
                withSpinner(plotOutput("parent_all_traits_plot", height = "600px"), 
                            type = 4, color = "#0dc5c1"),
                br(),
                h4("Complete Summary Table with All Traits"),
                withSpinner(DTOutput("parent_all_traits_table"), 
                            type = 4, color = "#0dc5c1"),
                br(),
                downloadButton("download_all_traits_summary", "Download Complete Summary (CSV)",
                               class = "btn-info")
              )
            )
          )
        ),
        
        # Results Download Section
        fluidRow(
          box(
            title = "Export Results", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            h4("Export All Results"),
            p("Download complete analysis results as a ZIP file containing all tables and visualizations."),
            br(),
            downloadButton("download_all_results", "Download Complete Results (ZIP)",
                           class = "btn-success btn-block")
          )
        )
      ),
      
      
      # METAN Analysis Tab
      tabItem(
        tabName = "metan",
        fluidRow(
          column(
            width = 12,
            h2("🌍 METAN Analysis", style = "margin-top: 0;"),
            helpText("Multi-Environment Trial Analysis using METAN package")
          )
        ),
        
        fluidRow(
          # METAN Data Setup
          box(
            title = "METAN Data Configuration", 
            status = "primary", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            
            # Add this UI output at the beginning of the box:
            uiOutput("metan_env_auto"),
            
            selectInput("metan_env_var", "Environment Variable (optional if YEAR & SEASON exist):", 
                        choices = NULL, width = "100%"),
            helpText("If your data has YEAR and SEASON columns, ENV will be created automatically. 
           Otherwise, select an environment variable."),
            
            selectInput("metan_gen_var", "Genotype Variable:", 
                        choices = NULL, width = "100%"),
            
            selectInput("metan_rep_var", "Replicate Variable:", 
                        choices = NULL, width = "100%"),
            
            selectizeInput("metan_resp_vars", "Response Variables:", 
                           choices = NULL, multiple = TRUE, width = "100%"),
            
            br(),
            actionButton("metan_setup_data", "Prepare METAN Data", 
                         icon = icon("cogs"), 
                         class = "btn-warning btn-block",
                         width = "100%")
          ),
          
          # METAN Analysis Controls
          box(
            title = "METAN Analysis Controls", 
            status = "info", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            checkboxGroupInput("metan_analyses", "Select Analyses:",
                               choices = c("Descriptive Statistics" = "desc",
                                           "ANOVA" = "anova",
                                           "AMMI" = "ammi",
                                           "WAASB" = "waasb",
                                           "GGE Biplot" = "gge",
                                           "Correlation" = "corr",
                                           "Stability" = "stability",
                                           "Ranking" = "ranking"),
                               selected = c("desc", "anova", "ammi"),
                               width = "100%"),
            sliderInput("metan_npc", "Number of PCs for AMMI:",
                        min = 1, max = 10, value = 3, step = 1,
                        width = "100%"),
            br(),
            actionButton("run_metan_analysis", "Run METAN Analysis", 
                         icon = icon("play"), 
                         class = "btn-success btn-block",
                         width = "100%")
          ),
          
          # METAN Data Status
          box(
            title = "METAN Data Status", 
            status = "success", 
            solidHeader = TRUE,
            width = 4,
            collapsible = TRUE,
            withSpinner(verbatimTextOutput("metan_data_status"), 
                        type = 4, color = "#0dc5c1"),
            br(),
            uiOutput("metan_env_stats_ui")
          )
        ),
        
        
        fluidRow(
          box(
            title = "METAN Analysis Results", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            tabsetPanel(
              tabPanel("Descriptive Statistics",
                       withSpinner(DTOutput("metan_desc_stats_table"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("ANOVA Results",
                       uiOutput("metan_trait_selector_anova"),
                       withSpinner(verbatimTextOutput("metan_anova_output"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("AMMI Analysis",
                       uiOutput("metan_trait_selector_ammi"),
                       withSpinner(verbatimTextOutput("metan_ammi_output"), 
                                   type = 4, color = "#0dc5c1"),
                       withSpinner(plotOutput("metan_ammi_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("WAASB Analysis",
                       withSpinner(verbatimTextOutput("metan_waasb_output"), 
                                   type = 4, color = "#0dc5c1"),
                       withSpinner(plotOutput("metan_waasb_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("GGE Biplot",
                       uiOutput("metan_trait_selector_gge"),
                       withSpinner(plotOutput("metan_gge_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("GE Interaction",
                       uiOutput("metan_trait_selector_ge"),
                       withSpinner(plotOutput("metan_ge_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Correlation Matrix",
                       withSpinner(plotOutput("metan_corr_plot", height = "500px"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Stability Analysis",
                       withSpinner(DTOutput("metan_stability_table"), 
                                   type = 4, color = "#0dc5c1")),
              tabPanel("Genotype Ranking",
                       withSpinner(DTOutput("metan_ranking_table"), 
                                   type = 4, color = "#0dc5c1"))
            )
          )
        )
      ),
      
      # Report Generation Tab
      tabItem(
        tabName = "report",
        h2("Report Generation"),
        fluidRow(
          box(
            title = "Report Configuration", width = 4, status = "primary",
            solidHeader = TRUE,
            textInput("report_title", "Report Title:", value = "GWAS Analysis Report"),
            selectInput("report_type", "Report Type:",
                        choices = c("Comprehensive", "Summary", "Custom")),
            selectInput("report_format", "Output Format:",
                        choices = c("HTML", "PDF", "Word")),
            selectInput("report_sections", "Include Sections:",
                        choices = c("Executive Summary", "Methods", 
                                    "Results", "Discussion", "Appendix"),
                        multiple = TRUE, selected = c("Executive Summary", "Results")),
            actionButton("generate_report", "Generate Report", 
                         class = "btn-success", icon = icon("file-download"))
            
          ),
          box(
            title = "Report Preview", width = 8, status = "info",
            solidHeader = TRUE,
            htmlOutput("report_preview"),
            downloadButton("download_report", "Download Report", 
                           class = "btn-success")
          )
        )
      ),
      
      # Add this to your tabItems in the dashboardBody
      tabItem(
        tabName = "all_dataframes",
        fluidRow(
          box(
            title = "📋 Available DataFrames", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = FALSE,
            DTOutput("available_dataframes_table"),
            br(),
            verbatimTextOutput("dataframe_metadata")
          )
        ),
        
        fluidRow(
          box(
            title = "🔍 DataFrame Viewer", 
            status = "success", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            selectInput("selected_dataframe", "Select DataFrame to View:", 
                        choices = NULL, width = "50%"),
            actionButton("refresh_dataframes", "Refresh List", 
                         icon = icon("sync"), class = "btn-info"),
            br(), br(),
            tabsetPanel(
              tabPanel("Data Preview", 
                       withSpinner(DTOutput("dataframe_preview"), type = 4)),
              tabPanel("Structure", 
                       withSpinner(verbatimTextOutput("dataframe_structure"), type = 4)),
              tabPanel("Summary", 
                       withSpinner(verbatimTextOutput("dataframe_summary"), type = 4)),
              tabPanel("Creation Function", 
                       withSpinner(verbatimTextOutput("dataframe_function"), type = 4))
            ),
            br(),
            downloadButton("download_dataframe", "Download as CSV", class = "btn-success")
          )
        )
      ),
      
      tabItem(
        tabName = "gwas_dataframes",
        fluidRow(
          box(
            title = "🧬 GWAS DataFrames", 
            status = "info", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            h4("DataFrames created by GWAS functions:"),
            tags$ul(
              tags$li(tags$strong("gwas_results:"), "Individual trait GWAS results - Created by run_single_trait_gwas()"),
              tags$li(tags$strong("combined_results:"), "Multi-trait combined results - Created by combine_multi_trait_results()"),
              tags$li(tags$strong("multi_trait_wide:"), "Wide format GWAS results - Created by create_wide_gwas_format()"),
              tags$li(tags$strong("multi_trait_enhanced:"), "Enhanced unified results - Created by create_enhanced_unified_df()"),
              tags$li(tags$strong("complete_integrated_df:"), "Complete integrated results - Created by create_complete_integrated_df()")
            ),
            br(),
            h4("Visualization DataFrames:"),
            tags$ul(
              tags$li(tags$strong("gwas_summary:"), "GWAS summary statistics - Created by capture_output_results()"),
              tags$li(tags$strong("merged_pvalues:"), "Merged P-values for multi-trait analysis - Created by combine_p_values()")
            )
          )
        )
      ),
      
      tabItem(
        tabName = "diallel_dataframes",
        fluidRow(
          box(
            title = "🌱 Diallel Analysis DataFrames", 
            status = "success", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            h4("Core Diallel DataFrames:"),
            tags$ul(
              tags$li(tags$strong("diallel_data:"), "Raw diallel data - Created by load_diallel_data()"),
              tags$li(tags$strong("diallel_results:"), "Diallel analysis results - Created by perform_diallel_analysis()"),
              tags$li(tags$strong("diallel_heterosis_data:"), "Heterosis calculations - Created by calculate_heterosis()")
            ),
            br(),
            h4("GCA/SCA DataFrames:"),
            tags$ul(
              tags$li(tags$strong("gca_df:"), "General Combining Ability effects - Created by perform_diallel_analysis()"),
              tags$li(tags$strong("sca_df:"), "Specific Combining Ability matrix - Created by sca_matrix_to_df()")
            ),
            br(),
            h4("Summary DataFrames:"),
            tags$ul(
              tags$li(tags$strong("diallel_summary_tables:"), "Summary of all traits - Created by run_diallel_anova_all_traits()"),
              tags$li(tags$strong("cross_summary:"), "Cross performance summary - Created by run_cross_snp_analysis_all_traits_app()")
            )
          )
        )
      ),
      
      tabItem(
        tabName = "cross_snp_dataframes",
        fluidRow(
          box(
            title = "🔗 Cross-SNP DataFrames", 
            status = "warning", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            h4("Main Cross-SNP DataFrames:"),
            tags$ul(
              tags$li(tags$strong("cross_snp_df:"), "Combined cross-SNP data for all traits - Created by create_cross_snp_dataframe_all_traits_fixed()"),
              tags$li(tags$strong("cross_snp_simple:"), "Simplified version - Created by create_cross_snp_dataframe_simple()")
            ),
            br(),
            h4("Summary DataFrames:"),
            tags$ul(
              tags$li(tags$strong("cross_summary:"), "Summary by cross - Created by run_cross_snp_analysis_all_traits_app()"),
              tags$li(tags$strong("top_crosses_snp_combined:"), "Summary by SNP across crosses - Created by run_cross_snp_analysis_all_traits_app()"),
              tags$li(tags$strong("trait_summary:"), "Summary by trait - Created by run_cross_snp_analysis_all_traits_app()")
            ),
            br(),
            h4("Utility DataFrames:"),
            tags$ul(
              tags$li(tags$strong("cross_summary_with_flanking:"), "Crosses with flanking markers - Created by add_flanking_snp_markers_safe()")
            )
          )
        )
      ),
      
      tabItem(
        tabName = "metan_dataframes",
        fluidRow(
          box(
            title = "🌍 METAN Analysis DataFrames", 
            status = "danger", 
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            h4("METAN Core DataFrames:"),
            tags$ul(
              tags$li(tags$strong("met_data:"), "Processed METAN data - Created by prepare_metan_data()"),
              tags$li(tags$strong("met_analysis_results:"), "METAN analysis results - Created by run_metan_analysis()")
            ),
            br(),
            h4("Analysis Component DataFrames:"),
            tags$ul(
              tags$li(tags$strong("desc_stats:"), "Descriptive statistics - Created by run_metan_analysis()"),
              tags$li(tags$strong("env_summary:"), "Environment characterization - Created by characterize_environments()"),
              tags$li(tags$strong("stability_df:"), "Stability analysis results - Created by run_stability_analysis()"),
              tags$li(tags$strong("missing_summary:"), "Missing data summary - Created by summarize_missing_data()")
            ),
            br(),
            h4("Specialized Analysis DataFrames:"),
            tags$ul(
              tags$li(tags$strong("ammi_results:"), "AMMI analysis results - Created by run_ammi_analysis()"),
              tags$li(tags$strong("gge_results:"), "GGE biplot results - Created by perform_gge_biplot()"),
              tags$li(tags$strong("fw_results:"), "Finlay-Wilkinson regression - Created by perform_finlay_wilkinson()")
            )
          )
        )
      )
    )
  )
)


# ============================================================================
# COMPLETE CORRECTED SERVER CODE
# ============================================================================

server <- function(input, output, session) {
  
  # Initialize reactive values
  values <- reactiveValues(
    # Core data
    pheno_data = NULL,
    gl_object = NULL,
    geno_data = NULL,
    filtered_gl = NULL,
    filtered_data = NULL,
    diallel_data = NULL,
    metan_raw_data = NULL,
    metadata = NULL,
    
    # Processed data
    matched_data = NULL,
    processed_geno = NULL,
    metan_processed_data = NULL,
    diallel_long = NULL,
    all_parents = NULL,
    parent_summary = NULL,
    diallel_processed = NULL,
    diallel_cross_data = NULL,
    diallel_parent_data = NULL,
    diallel_results = NULL,
    diallel_all_traits_results = NULL,
    
    # Analysis results
    qc_plots = NULL,
    pca_results = NULL,
    kinship_matrix = NULL,
    gwas_results = NULL,
    multi_trait_results = NULL,
    metan_results = NULL,
    dartr_results = list(),
    
    # unified SNP info for all traits
    multi_trait_unified_df = NULL,
    multi_trait_unified_tidy = NULL,
    multi_trait_unified_wide = NULL,
    multi_trait_enhanced_unified_df = NULL,
    multi_trait_complete_integrated_df = NULL,
    
    #Combine cross and SNP data
    cross_snp_analysis = NULL,
    top_crosses_snp_combined =NULL,
    cross_snp_detailed = NULL,
    cross_summary = NULL,
    cross_snp_analysis = NULL,
    
    cross_snp_df = NULL,
    
    
    # Reactive values for Cross-SNP analysis
    dl_cross_snp_detailed = NULL,
    dl_cross_summary = NULL,
    dl_snp_summary = NULL,
    dl_cross_snp_all = NULL,
    
    # Status flags
    pheno_loaded = FALSE,
    geno_loaded = FALSE,
    data_matched = FALSE,
    qc_completed = FALSE,
    
    # Report content
    report_content = NULL,
    report_ready = FALSE
  )
  
  # ==========================================================================
  # DATA UPLOAD OBSERVERS
  # ==========================================================================
  
  summarize_by_identifier <- function(df) {
    tryCatch({
      require(dplyr)
      
      # List of possible identifier names
      grouping_names <- c("Genotype", "GEN", "gen", "genotype", "GENOTYPE", 
                          "Sample", "sample", "ID", "id", "Name", "name", 
                          "Accession", "ACC", "VARIETY", "Variety")
      
      # Find which grouping columns are present
      grouping_cols <- grouping_names[grouping_names %in% names(df)]
      
      if (length(grouping_cols) == 0) {
        # Try case-insensitive matching
        colnames_lower <- tolower(names(df))
        grouping_lower <- tolower(grouping_names)
        
        for (i in seq_along(grouping_lower)) {
          if (grouping_lower[i] %in% colnames_lower) {
            actual_col <- names(df)[which(colnames_lower == grouping_lower[i])[1]]
            grouping_cols <- actual_col
            break
          }
        }
      }
      
      if (length(grouping_cols) == 0) {
        warning("No grouping columns found. Returning original data.")
        
        if (ncol(df) > 0) {
          df <- df[order(df[, 1]), ]
        }
        rownames(df) <- NULL
        return(df)
      }
      
      cat("Grouping by:", paste(grouping_cols, collapse = ", "), "\n")
      
      numeric_cols <- names(df)[sapply(df, is.numeric)]
      numeric_to_average <- setdiff(numeric_cols, grouping_cols)
      
      if (length(numeric_to_average) == 0) {
        warning("No numeric columns to average. Returning original data.")
        df <- df[order(df[[grouping_cols[1]]]), ]
        rownames(df) <- NULL
        return(df)
      }
      
      cat("Calculating means for:", paste(numeric_to_average, collapse = ", "), "\n")
      
      result <- df %>%
        group_by(across(all_of(grouping_cols))) %>%
        summarise(
          across(all_of(numeric_to_average), 
                 ~ mean(., na.rm = TRUE),
                 .names = "{.col}"),
          .groups = "drop"
        ) %>%
        arrange(across(all_of(grouping_cols[1])))
      
      if (!"Genotype" %in% names(result) && length(grouping_cols) > 0) {
        result <- result %>%
          rename(Genotype = !!sym(grouping_cols[1]))
        cat("Renamed grouping column '", grouping_cols[1], "' to 'Genotype'\n")
      }
      
      result <- as.data.frame(result)
      rownames(result) <- NULL
      
      return(result)
      
    }, error = function(e) {
      warning(paste("Error in summarize_by_identifier:", e$message))
      if (ncol(df) > 0) {
        df <- df[order(df[, 1]), ]
        rownames(df) <- NULL
      }
      return(df)
    })
  }
  
  # Load phenotypic data with withProgress
  observeEvent(input$load_pheno, {
    req(input$pheno_file)
    
    withProgress(message = 'Loading phenotypic data...', value = 0, {
      
      tryCatch({
        incProgress(0.2, detail = "Reading file...")
        
        raw_pheno_data <- load_pheno_data(
          input$pheno_file$datapath,
          header = input$pheno_header,
          aggregate = FALSE
        )
        
        incProgress(0.4, detail = "Processing data...")
        
        values$pheno_data_raw <- raw_pheno_data
        values$metan_raw_data <- raw_pheno_data
        
        cat("Raw phenotype data loaded. Dimensions:", dim(raw_pheno_data), "\n")
        
        if (input$aggregate_pheno) {
          aggregated_data <- summarize_by_identifier(raw_pheno_data)
          values$pheno_data <- aggregated_data
          cat("Aggregated phenotype data created. Dimensions:", dim(aggregated_data), "\n")
        } else {
          values$pheno_data <- raw_pheno_data
          cat("Using raw data (no aggregation).\n")
        }
        
        incProgress(0.6, detail = "Updating trait selectors...")
        
        traits <- setdiff(colnames(values$pheno_data), "Genotype")
        
        updateSelectInput(session, "diallel_trait", choices = traits)
        updateSelectInput(session, "gwas_trait", choices = traits)
        updateSelectizeInput(session, "multi_traits", choices = traits)
        
        values$pheno_loaded <- TRUE
        
        incProgress(0.8, detail = "Finalizing...")
        
        output$pheno_status <- renderUI({
          div(
            class = "alert alert-success",
            strong("✓ Phenotype data loaded successfully!"),
            br(),
            tags$b("Raw Data:"),
            br(),
            paste("Total observations:", nrow(values$pheno_data_raw)),
            br(),
            paste("Columns:", paste(colnames(values$pheno_data_raw), collapse = ", ")),
            br(),
            tags$hr(),
            tags$b("Processed Data:"),
            br(),
            paste("Unique genotypes:", nrow(values$pheno_data)),
            br(),
            paste("Traits:", paste(traits, collapse = ", ")),
            br(),
            tags$hr(),
            tags$b("METAN Data:"),
            br(),
            "✓ Phenotype data is also available for METAN analysis",
            br(),
            "To use a different dataset, upload a separate METAN file"
          )
        })
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Phenotypic data loaded successfully!", 
                         type = "message", duration = 3)
        
      }, error = function(e) {
        output$pheno_status <- renderUI({
          div(class = "alert alert-danger",
              strong("✗ Error loading phenotypic data:"),
              br(),
              e$message
          )
        })
        showNotification(paste("Error loading phenotypic data:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  # Load genotypic data with withProgress
  observeEvent(input$load_geno, {
    req(input$geno_file)
    
    withProgress(message = 'Loading genotypic data...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = "Reading genotype file...")
        
        gl_object <- create_genlight(input$geno_file$datapath, input$metadata_file$datapath)
        values$gl_object <- clean_genlight_chromosomes(gl_object, verbose = TRUE)
        
        if (is.null(values$gl_object)) {
          stop("Failed to create genlight object")
        }
        
        incProgress(0.6, detail = "Processing genotype matrix...")
        
        values$geno_loaded <- TRUE
        values$geno_data <- as.matrix(values$gl_object)
        
        incProgress(0.9, detail = "Finalizing...")
        
        output$geno_status <- renderUI({
          tagList(
            div(
              class = "alert alert-success",
              strong("✓ Genotype data loaded successfully!"),
              tags$br(),
              tags$p(paste("Samples:", nInd(values$gl_object))),
              tags$p(paste("Markers:", nLoc(values$gl_object))),
              tags$p(paste("Ploidy:", paste(unique(values$gl_object@ploidy), collapse = ", ")))
            )
          )
        })
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Genotypic data loaded successfully!", type = "message", duration = 3)
        
      }, error = function(e) {
        output$geno_status <- renderUI({
          tagList(
            div(
              class = "alert alert-danger",
              strong("✗ Error loading genotype data:"),
              tags$br(),
              tags$p(as.character(e$message))
            )
          )
        })
        showNotification(paste("Error loading genotypic data:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # Load metadata with withProgress
  observeEvent(input$load_metadata, {
    req(input$metadata_file)
    
    withProgress(message = 'Loading metadata...', value = 0, {
      
      tryCatch({
        incProgress(0.5, detail = "Reading file...")
        
        values$metadata <- load_meta_data(input$metadata_file$datapath)
        
        incProgress(0.8, detail = "Finalizing...")
        
        output$metadata_status <- renderUI({
          div(class = "alert alert-success",
              strong("✓ Metadata loaded successfully!"),
              br(),
              paste("Entries:", nrow(values$metadata)),
              br(),
              paste("Variables:", ncol(values$metadata))
          )
        })
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Metadata loaded successfully!", type = "message", duration = 3)
        
      }, error = function(e) {
        output$metadata_status <- renderUI({
          div(class = "alert alert-danger",
              strong("✗ Error loading metadata:"),
              br(),
              e$message
          )
        })
        showNotification(paste("Error loading metadata:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # Load diallel data with withProgress
  observeEvent(input$load_diallel, {
    req(input$diallel_file)
    
    withProgress(message = 'Loading diallel data...', value = 0, {
      
      tryCatch({
        incProgress(0.2, detail = "Reading file...")
        
        df <- load_diallel_data(input$diallel_file$datapath)
        
        incProgress(0.4, detail = "Cleaning data...")
        
        colnames(df) <- gsub("100SW", "X100SW", colnames(df))
        df[df == "" | df == "NA" | df == "NULL"] <- NA
        
        trait_cols <- c("DTF", "DTM", "PL", "NPP", "NSP", "X100SW", "GYP", "GYR", "TW")
        for(col in trait_cols) {
          if(col %in% colnames(df)) {
            df[[col]] <- as.numeric(as.character(df[[col]]))
          }
        }
        
        incProgress(0.6, detail = "Processing crosses...")
        
        df <- df %>%
          mutate(
            CrossType = ifelse(grepl("x", ACC.), "Cross", "Parent"),
            Parent1 = ifelse(CrossType == "Cross", 
                             sapply(strsplit(as.character(ACC.), "x"), function(x) x[1]), 
                             as.character(ACC.)),
            Parent2 = ifelse(CrossType == "Cross", 
                             sapply(strsplit(as.character(ACC.), "x"), function(x) ifelse(length(x) > 1, x[2], NA)), 
                             as.character(ACC.)),
            CrossID = as.character(ACC.)
          ) %>%
          filter(!is.na(Parent1))
        
        values$diallel_data <- df
        
        incProgress(0.8, detail = "Separating parents and crosses...")
        
        parents <- df[df$CrossType == "Parent" & !is.na(df$CrossType), ]
        trait_cols <- trait_cols[trait_cols %in% colnames(parents)]
        parents <- parents[apply(parents[, trait_cols], 1, function(x) !all(is.na(x))), ]
        values$diallel_parent_data <- parents
        
        crosses <- df[df$CrossType == "Cross" & !is.na(df$CrossType) & !is.na(df$Parent2), ]
        values$diallel_cross_data <- crosses
        
        updateSelectInput(session, "diallel_trait", choices = trait_cols, selected = "GYP")
        
        incProgress(0.95, detail = "Finalizing...")
        
        output$diallel_status <- renderUI({
          div(class = "alert alert-success",
              strong("✓ Diallel data loaded successfully!"),
              br(),
              paste("Total entries:", nrow(df)),
              br(),
              paste("Parents:", nrow(parents)),
              br(),
              paste("Crosses:", nrow(crosses))
          )
        })
        
        incProgress(1.0, detail = "Complete!")
        
      }, error = function(e) {
        output$diallel_status <- renderUI({
          div(class = "alert alert-danger",
              strong("✗ Error reading diallel file:"),
              br(),
              e$message
          )
        })
      })
    })
  })
  
  
  # Initialize reactive values
  values$parent_snp_actual <- NULL
  values$parent_genotype_patterns <- NULL
  
  
  # Data status check (shows before running analysis)
  output$parent_snp_data_check <- renderPrint({
    req(values$diallel_long, values$matched_data, values$gwas_results)
    
    cat("=== DATA AVAILABILITY CHECK ===\n\n")
    
    # Check diallel parents
    unique_parents <- unique(values$diallel_long$Parent)
    cat("Unique parents in diallel data:", length(unique_parents), "\n")
    cat("Class of Parent column:", class(values$diallel_long$Parent), "\n")
    cat("Sample parents (first 5):", paste(head(unique_parents, 5), collapse = ", "), "\n")
    
    # Check genotype matrix
    if (!is.null(values$matched_data$geno)) {
      geno_samples <- rownames(values$matched_data$geno)
      cat("\nSamples in genotype matrix:", length(geno_samples), "\n")
      cat("Class of rownames:", class(geno_samples), "\n")
      cat("Sample names (first 5):", paste(head(geno_samples, 5), collapse = ", "), "\n")
      
      # Find common parents
      common_parents <- intersect(unique_parents, geno_samples)
      cat("Common parents found:", length(common_parents), "\n")
      
      if (length(common_parents) > 0) {
        cat("Common parents:", paste(common_parents, collapse = ", "), "\n")
      }
    }
    
    # Check GWAS results
    if (!is.null(values$gwas_results)) {
      cat("\nGWAS results available:", nrow(values$gwas_results), "SNPs\n")
      cat("Class of SNP_ID:", class(values$gwas_results$SNP_ID), "\n")
      cat("Class of gwas_results:", class(values$gwas_results), "\n")
    }
  })
  
  
  
  
  # Data status UI
  output$parent_snp_data_status_ui <- renderUI({
    req(values$diallel_long, values$matched_data, values$gwas_results)
    
    # Check if we have all required data
    all_data_available <- !is.null(values$diallel_long) && 
      !is.null(values$matched_data$geno) && 
      !is.null(values$gwas_results)
    
    if (!all_data_available) {
      return(
        div(class = "alert alert-danger",
            icon("exclamation-triangle"),
            strong(" Missing data! "),
            "Please ensure diallel, genotype, and GWAS data are loaded.")
      )
    }
    
    # Check for common parents
    unique_parents <- unique(values$diallel_long$Parent)
    geno_samples <- rownames(values$matched_data$geno)
    common_parents <- intersect(unique_parents, geno_samples)
    
    if (length(common_parents) == 0) {
      return(
        div(class = "alert alert-warning",
            icon("exclamation-circle"),
            strong(" Warning: No common parents found! "),
            "Parent names in diallel data don't match genotype sample names.")
      )
    }
    
    return(
      div(class = "alert alert-success",
          icon("check-circle"),
          strong(" Data ready for analysis! "),
          "Click 'Run Parent-SNP Matching' to start.")
    )
  })
  
  
  
  
  # In the observeEvent for run_parent_snp_matching - FIXED showNotification calls
  observeEvent(input$run_parent_snp_matching, {
    
    showModal(modalDialog(
      title = "Running Parent-SNP Matching",
      tagList(
        tags$div(
          id = "snp_matching_progress",
          style = "text-align: center;",
          icon("spinner", class = "fa-spin fa-2x"),
          tags$h4("Matching parents to actual genotypes..."),
          tags$div(
            id = "snp_matching_details",
            style = "margin-top: 10px; font-size: 12px; color: #666;"
          )
        )
      ),
      footer = NULL,
      easyClose = FALSE,
      size = "s"
    ))
    
    tryCatch({
      shinyjs::html("snp_matching_details", "Step 1/4: Checking data availability...")
      Sys.sleep(0.5)
      
      # Validate required data
      if (is.null(values$diallel_long)) {
        removeModal()
        showNotification("Diallel data not available", type = "error")
        return()
      }
      if (is.null(values$matched_data$geno)) {
        removeModal()
        showNotification("Genotype matrix not available", type = "error")
        return()
      }
      if (is.null(values$gwas_results)) {
        removeModal()
        showNotification("GWAS results not available", type = "error")
        return()
      }
      
      shinyjs::html("snp_matching_details", "Step 2/4: Extracting parent genotypes...")
      Sys.sleep(0.5)
      
      # Call the matching function
      matching_results <- match_parents_to_gwas_actual(
        diallel_long = values$diallel_long,
        geno_matrix = values$matched_data$geno,
        gwas_results = values$gwas_results
      )
      
      if (is.null(matching_results) || (is.list(matching_results) && !is.null(matching_results$error))) {
        removeModal()
        error_msg <- if (!is.null(matching_results$error)) matching_results$error else "Unknown error"
        showNotification(paste("Parent-SNP matching failed:", error_msg), 
                         type = "error", duration = 10)
        return()
      }
      
      if (!matching_results$success) {
        removeModal()
        showNotification("Parent-SNP matching was not successful", 
                         type = "error", duration = 10)
        return()
      }
      
      shinyjs::html("snp_matching_details", "Step 3/4: Analyzing significant SNP patterns...")
      Sys.sleep(0.5)
      
      # Analyze significant SNP patterns
      genotype_patterns <- analyze_parent_significant_snps(
        matching_results$detailed_genotypes,
        p_threshold = input$snp_p_threshold
      )
      
      # Store results
      values$parent_snp_actual <- matching_results
      values$parent_genotype_patterns <- genotype_patterns
      
      shinyjs::html("snp_matching_details", "Step 4/4: Finalizing results...")
      Sys.sleep(0.5)
      
      removeModal()
      
      # Show success notification - FIXED: Use "default" instead of "success"
      showNotification(
        HTML(paste(
          "✅ Parent-SNP matching completed successfully!<br>",
          "Analyzed:", matching_results$parent_count, "parents and", 
          matching_results$snp_count, "SNPs<br>",
          "Significant combinations:", matching_results$significant_combinations
        )),
        type = "default",  # CHANGED FROM "success" TO "default"
        duration = 5
      )
      
      # Update output
      output$parent_snp_actual_stats <- renderPrint({
        cat("=== PARENT-SNP MATCHING RESULTS ===\n\n")
        cat("✅ Analysis completed successfully!\n\n")
        cat("Parents analyzed:", matching_results$parent_count, "\n")
        cat("SNPs analyzed:", matching_results$snp_count, "\n")
        cat("Total parent-SNP combinations:", nrow(matching_results$detailed_genotypes), "\n")
        cat("Significant combinations (p < 0.05):", matching_results$significant_combinations, "\n")
        cat("Analysis timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
      })
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Error in Parent-SNP matching:", e$message),
        type = "error",
        duration = 10
      )
      cat("Parent-SNP matching observer error:", e$message, "\n")
      cat("Traceback:\n")
      print(traceback())
    })
  })
  
  
  # Clear results button
  observeEvent(input$clear_parent_snp_results, {
    values$parent_snp_actual <- NULL
    values$parent_genotype_patterns <- NULL
    
    # Reset output
    output$parent_snp_actual_table <- renderDT({NULL})
    output$parent_genotype_detailed <- renderDT({NULL})
    output$genotype_patterns_table <- renderDT({NULL})
    output$parent_snp_actual_stats <- renderPrint({
      cat("=== PARENT-SNP MATCHING ===\n\n")
      cat("No results available.\n")
      cat("Click 'Run Parent-SNP Matching' to analyze.\n")
    })
    
    showNotification(
      "Parent-SNP results cleared successfully!",
      type = "warning",
      duration = 3
    )
  })
  
  # Alternative: Button with progress bar
  observeEvent(input$run_parent_snp_matching_with_progress, {
    # Alternative implementation with withProgress
    withProgress(message = 'Running Parent-SNP Matching', value = 0, {
      
      incProgress(0.1, detail = "Validating data...")
      
      # Validate data
      if (is.null(values$diallel_long) || 
          is.null(values$matched_data$geno) || 
          is.null(values$gwas_results)) {
        showNotification("Missing required data!", type = "error")
        return()
      }
      
      incProgress(0.2, detail = "Extracting parent genotypes...")
      
      # Perform matching
      matching_results <- match_parents_to_gwas_actual(
        diallel_long = values$diallel_long,
        geno_matrix = values$matched_data$geno,
        gwas_results = values$gwas_results
      )
      
      incProgress(0.7, detail = "Analyzing patterns...")
      
      if (!is.null(matching_results)) {
        # Analyze significant SNPs
        genotype_patterns <- analyze_parent_significant_snps(
          matching_results$detailed_genotypes,
          p_threshold = input$snp_p_threshold
        )
        
        # Store results
        values$parent_snp_actual <- matching_results
        values$parent_genotype_patterns <- genotype_patterns
      }
      
      incProgress(1.0, detail = "Complete!")
      
      if (!is.null(matching_results)) {
        showNotification(
          paste("Analysis complete! Processed", 
                matching_results$parent_count, "parents and",
                matching_results$snp_count, "SNPs"),
          type = "message",
          duration = 5
        )
      }
    })
  })
  
  # Check if analysis can be run (enable/disable button)
  observe({
    # Check if all required data is available
    data_available <- !is.null(values$diallel_long) && 
      !is.null(values$matched_data$geno) && 
      !is.null(values$gwas_results)
    
    # Check if there are common parents
    if (data_available) {
      unique_parents <- unique(values$diallel_long$Parent)
      geno_samples <- rownames(values$matched_data$geno)
      common_parents <- intersect(unique_parents, geno_samples)
      data_available <- length(common_parents) > 0
    }
    
    # Enable/disable button
    if (data_available) {
      shinyjs::enable("run_parent_snp_matching")
      shinyjs::runjs("
      $('#run_parent_snp_matching').removeClass('disabled');
      $('#run_parent_snp_matching').html(\"<i class='fa fa-dna'></i> Run Parent-SNP Matching\");
    ")
    } else {
      shinyjs::disable("run_parent_snp_matching")
      shinyjs::runjs("
      $('#run_parent_snp_matching').addClass('disabled');
      $('#run_parent_snp_matching').html(\"<i class='fa fa-ban'></i> Data Not Ready\");
    ")
    }
  })
  
  # Output renderers (conditional based on results)
  output$parent_snp_actual_table <- renderDT({
    req(values$parent_snp_actual$summary_table)
    
    # Format the table
    formatted_table <- values$parent_snp_actual$summary_table %>%
      mutate(
        Mean_GYP = round(Mean_GYP, 2),
        Mean_GYR = round(Mean_GYR, 2),
        Mean_Genotype = round(Mean_Genotype, 2),
        Top_SNP_LOD = round(Top_SNP_LOD, 2),
        Top_SNP_Effect = round(Top_SNP_Effect, 4),
        Top_SNP_PVE = round(Top_SNP_PVE, 2),
        Min_P_value = format(Min_P_value, scientific = TRUE, digits = 3),
        Missing_Genotype_Percent = round(Missing_Genotype/N_SNPs * 100, 1)
      )
    
    datatable(
      formatted_table,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Parent-SNP Association Summary (Actual Genotypes)",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'Missing_Genotype_Percent',
        backgroundColor = styleInterval(
          c(10, 20, 50),
          c('#91ff91', '#fff991', '#ffb291', '#ff9191')
        )
      )
  })
  
  # Update the data status display after analysis
  observe({
    req(values$parent_snp_actual)
    
    output$parent_snp_data_status_ui <- renderUI({
      result <- values$parent_snp_actual
      
      div(
        class = "alert alert-success",
        icon("check-circle"),
        strong(" Analysis Complete! "),
        br(),
        paste("Processed:", result$parent_count, "parents and", 
              result$snp_count, "SNPs"),
        br(),
        paste("Total combinations:", nrow(result$detailed_genotypes)),
        br(),
        tags$small(paste("Completed:", format(Sys.time(), "%H:%M:%S")))
      )
    })
  })
  
  # Export results button with confirmation
  observeEvent(input$export_parent_snp_results, {
    if (is.null(values$parent_snp_actual)) {
      showNotification("No results to export. Run analysis first.", type = "warning")
      return()
    }
    
    # Ask for confirmation
    showModal(modalDialog(
      title = "Export Results",
      "Export all parent-SNP matching results as a ZIP file?",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_export", "Export", class = "btn-success")
      )
    ))
  })
  
  observeEvent(input$confirm_export, {
    removeModal()
    
    # Trigger the download
    shinyjs::runjs("$('#download_parent_snp_actual')[0].click();")
  })
  
  # Download handler remains the same as before
  output$download_parent_snp_actual <- downloadHandler(
    filename = function() {
      paste0("parent_snp_matching_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      # Create temporary directory
      temp_dir <- tempfile("parent_snp_")
      dir.create(temp_dir)
      
      # Save all data frames
      if (!is.null(values$parent_snp_actual$summary_table)) {
        write.csv(values$parent_snp_actual$summary_table,
                  file.path(temp_dir, "parent_snp_summary.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$parent_snp_actual$detailed_genotypes)) {
        write.csv(values$parent_snp_actual$detailed_genotypes,
                  file.path(temp_dir, "detailed_parent_genotypes.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$parent_genotype_patterns)) {
        write.csv(values$parent_genotype_patterns,
                  file.path(temp_dir, "significant_snp_patterns.csv"),
                  row.names = FALSE)
      }
      
      # Save metadata
      meta_data <- data.frame(
        Parameter = c("Analysis Date", "Parents Analyzed", "SNPs Analyzed", 
                      "Total Combinations", "P-value Threshold"),
        Value = c(
          as.character(Sys.time()),
          values$parent_snp_actual$parent_count,
          values$parent_snp_actual$snp_count,
          nrow(values$parent_snp_actual$detailed_genotypes),
          input$snp_p_threshold
        )
      )
      write.csv(meta_data, file.path(temp_dir, "metadata.csv"), row.names = FALSE)
      
      # Create README
      readme_content <- paste(
        "PARENT-SNP MATCHING RESULTS",
        "============================",
        "",
        paste("Analysis Date:", Sys.time()),
        paste("Generated by: Plant Breeding Analysis Platform v4.0"),
        "",
        "Files included:",
        "1. parent_snp_summary.csv - Summary table for each parent",
        "2. detailed_parent_genotypes.csv - All parent-SNP combinations",
        "3. significant_snp_patterns.csv - Patterns at significant SNPs",
        "4. metadata.csv - Analysis parameters and metadata",
        "",
        "Column descriptions:",
        "- Parent: Parent name from diallel data",
        "- SNP_ID: SNP identifier",
        "- Genotype: Actual genotype value (0=AA, 1=Aa, 2=aa)",
        "- P_value: GWAS p-value for the SNP",
        "- Effect: Additive effect size",
        "- R_squared: Phenotypic variance explained",
        "- LOD: -log10(P_value)",
        "- PVE: Phenotypic variance explained (%)",
        sep = "\n"
      )
      writeLines(readme_content, file.path(temp_dir, "README.txt"))
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "results.zip")
      zip::zip(zip_file, files = list.files(temp_dir, full.names = TRUE), 
               mode = "cherry-pick")
      
      # Copy to download location
      file.copy(zip_file, file)
      
      # Cleanup
      unlink(temp_dir, recursive = TRUE)
    }
  )
  
  # ==========================================================================
  # METAN DATA HANDLING OBSERVERS
  # ==========================================================================
  
  observe({
    req(values$metan_raw_data)
    
    tryCatch({
      updateSelectInput(session, "metan_env_var", choices = names(values$metan_raw_data))
      updateSelectInput(session, "metan_gen_var", choices = names(values$metan_raw_data))
      updateSelectInput(session, "metan_rep_var", choices = names(values$metan_raw_data))
      
      cat("METAN selectors updated with data containing", 
          ncol(values$metan_raw_data), "columns\n")
      
    }, error = function(e) {
      cat("Error updating METAN selectors:", e$message, "\n")
    })
  })
  
  observe({
    req(values$metan_raw_data)
    
    tryCatch({
      all_cols <- names(values$metan_raw_data)
      
      output$metan_env_auto <- renderUI({
        if (all(c("YEAR", "SEASON") %in% all_cols)) {
          div(class = "alert alert-info",
              icon("info-circle"),
              " ✓ Your data has YEAR and SEASON columns.",
              br(),
              "ENV will be created automatically by combining YEAR and SEASON.",
              br(),
              tags$small("Example: ENV = '2023_WetSeason'")
          )
        } else if ("ENV" %in% all_cols) {
          div(class = "alert alert-info",
              icon("info-circle"),
              " ✓ Found ENV column in your data.",
              br(),
              "You can select this as the environment variable.")
        } else {
          div(class = "alert alert-warning",
              icon("exclamation-triangle"),
              " No YEAR/SEASON or ENV column found.",
              br(),
              "Please select an environment variable from the dropdown.")
        }
      })
      
      env_candidates <- c("ENV", "YEAR", "SEASON", "Environment", "environment", 
                          "LOCATION", "Location", "Site", "Trial")
      env_choices <- c("", all_cols)
      
      env_default <- env_candidates[env_candidates %in% all_cols][1]
      if (!is.na(env_default)) {
        updateSelectInput(session, "metan_env_var", 
                          choices = env_choices,
                          selected = env_default)
      } else {
        updateSelectInput(session, "metan_env_var", 
                          choices = env_choices,
                          selected = "")
      }
      
      gen_candidates <- c("GEN", "Genotype", "genotype", "Name", "ID", "VARIETY", "Variety", "ACC")
      gen_default <- gen_candidates[gen_candidates %in% all_cols][1]
      updateSelectInput(session, "metan_gen_var", 
                        choices = all_cols,
                        selected = ifelse(!is.na(gen_default), gen_default, all_cols[1]))
      
      rep_candidates <- c("REP", "Rep", "rep", "REPLICATE", "Replicate", "Block", "BLOCK")
      rep_default <- rep_candidates[rep_candidates %in% all_cols][1]
      updateSelectInput(session, "metan_rep_var", 
                        choices = all_cols,
                        selected = ifelse(!is.na(rep_default), rep_default, all_cols[2]))
      
      numeric_cols <- names(values$metan_raw_data)[sapply(values$metan_raw_data, is.numeric)]
      id_cols <- c("YEAR", "SEASON", "ENV", "GEN", "Genotype", "REP", "Rep", 
                   "ID", "Name", "Block", "Location", "Trial")
      trait_candidates <- setdiff(numeric_cols, id_cols)
      
      if (length(trait_candidates) == 0) {
        trait_candidates <- numeric_cols
      }
      
      updateSelectizeInput(session, "metan_resp_vars", 
                           choices = trait_candidates,
                           selected = head(trait_candidates, min(3, length(trait_candidates))))
      
      cat("METAN selectors updated with data containing", 
          ncol(values$metan_raw_data), "columns\n")
      
    }, error = function(e) {
      cat("Error updating METAN selectors:", e$message, "\n")
    })
  })
  
  # Load separate METAN file with withProgress
  observeEvent(input$load_metan, {
    req(input$metan_file)
    
    withProgress(message = 'Loading METAN data...', value = 0, {
      
      tryCatch({
        incProgress(0.5, detail = "Reading file...")
        
        metan_data <- load_metan_data(input$metan_file$datapath)
        values$metan_raw_data <- metan_data
        
        incProgress(0.8, detail = "Finalizing...")
        
        output$metan_status <- renderUI({
          div(class = "alert alert-success",
              strong("✓ Separate METAN data loaded successfully!"),
              br(),
              paste("Observations:", nrow(values$metan_raw_data)),
              br(),
              paste("Variables:", ncol(values$metan_raw_data)),
              br(),
              tags$hr(),
              tags$small("Note: This data will be used for METAN analysis instead of phenotype data.")
          )
        })
        
        incProgress(1.0, detail = "Complete!")
        showNotification("METAN data loaded successfully!", type = "message", duration = 3)
        
      }, error = function(e) {
        output$metan_status <- renderUI({
          div(class = "alert alert-danger",
              strong("✗ Error loading separate METAN data:"),
              br(),
              e$message
          )
        })
        showNotification(paste("Error loading METAN data:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # Clear METAN data button
  observeEvent(input$clear_metan, {
    values$metan_raw_data <- NULL
    values$metan_processed_data <- NULL
    values$metan_results <- NULL
    
    updateSelectInput(session, "metan_env_var", choices = NULL)
    updateSelectInput(session, "metan_gen_var", choices = NULL)
    updateSelectInput(session, "metan_rep_var", choices = NULL)
    updateSelectizeInput(session, "metan_resp_vars", choices = NULL)
    
    output$metan_status <- renderUI({
      div(class = "alert alert-warning",
          strong("METAN data cleared!"),
          br(),
          "You can either use phenotype data or load a separate METAN file.")
    })
    
    showNotification("METAN data cleared!", type = "warning", duration = 3)
  })
  
  # ==========================================================================
  # DATA PROCESSING OBSERVERS
  # ==========================================================================
  
  # Run quality control with withProgress
  observeEvent(input$run_qc, {
    req(values$gl_object)
    
    withProgress(message = 'Running quality control...', value = 0, {
      
      tryCatch({
        incProgress(0.2, detail = "Filtering by call rate...")
        
        values$filtered_gl <- dartR::gl.filter.callrate(values$gl_object, 
                                                        threshold = input$callrate_threshold)
        
        incProgress(0.4, detail = "Filtering by MAF...")
        
        values$filtered_gl <- dartR::gl.filter.maf(values$filtered_gl, 
                                                   threshold = input$maf_threshold)
        
        incProgress(0.6, detail = "Calculating genotype stats...")
        
        values$geno_data <- genlight_to_matrix_stats(values$filtered_gl)
        
        incProgress(0.7, detail = "Creating QC plots...")
        
        values$qc_plots <- create_qc_plots(values$filtered_gl)
        values$qc_completed <- TRUE
        
        incProgress(0.8, detail = "Preparing genotype data...")
        
        values$processed_geno <- prepare_geno_data(values$filtered_gl)
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Quality control completed!", type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error during quality control:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # Match genotype and phenotype data with withProgress
  observeEvent(input$match_data, {
    req(values$pheno_data, values$geno_data)
    
    withProgress(message = 'Matching genotype and phenotype data...', value = 0, {
      
      tryCatch({
        incProgress(0.5, detail = "Matching samples...")
        
        values$matched_data <- extract_and_match_genotypes(
          values$geno_data, 
          values$pheno_data, 
          start_col = input$start_col
        )
        
        values$data_matched <- TRUE
        
        incProgress(0.8, detail = "Finalizing...")
        
        output$match_status <- renderUI({
          div(class = "alert alert-success",
              strong("✓ Data matching completed!"),
              br(),
              paste("Matched genotypes:", nrow(values$matched_data$pheno)),
              br(),
              paste("Markers:", ncol(values$matched_data$geno))
          )
        })
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Data matching completed!", type = "message", duration = 3)
        
      }, error = function(e) {
        output$match_status <- renderUI({
          div(class = "alert alert-danger",
              strong("✗ Error matching data:"),
              br(),
              e$message
          )
        })
        showNotification(paste("Error matching data:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # Check data compatibility with withProgress
  observeEvent(input$check_data, {
    req(values$pheno_data, values$geno_data)
    
    withProgress(message = 'Checking data compatibility...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = "Extracting sample names...")
        
        pheno_samples <- values$pheno_data$Genotype
        geno_samples <- rownames(values$geno_data)
        
        clean_names <- function(x) {
          x <- as.character(x)
          x <- trimws(x)
          x <- gsub("^\\s+|\\s+$", "", x)
          return(x)
        }
        
        pheno_clean <- clean_names(pheno_samples)
        geno_clean <- clean_names(geno_samples)
        
        incProgress(0.6, detail = "Finding matches...")
        
        exact_matches <- intersect(pheno_clean, geno_clean)
        case_insensitive_matches <- intersect(tolower(pheno_clean), tolower(geno_clean))
        
        incProgress(0.8, detail = "Generating report...")
        
        output$compatibility_status <- renderUI({
          div(class = "alert alert-info",
              h4("Data Compatibility Check"),
              p("Phenotype samples:", length(pheno_clean)),
              p("Genotype samples:", length(geno_clean)),
              p("Exact matches:", length(exact_matches)),
              p("Case-insensitive matches:", length(case_insensitive_matches)),
              if (length(exact_matches) < 10) {
                p("Matching samples:", paste(exact_matches, collapse = ", "))
              } else {
                p("First 10 matching samples:", paste(head(exact_matches, 10), collapse = ", "))
              }
          )
        })
        
        incProgress(1.0, detail = "Complete!")
        
      }, error = function(e) {
        showNotification(paste("Error checking compatibility:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # ==========================================================================
  # POPULATION STRUCTURE OBSERVERS
  # ==========================================================================
  
  # Run PCA with withProgress
  observeEvent(input$run_pca, {
    req(values$filtered_gl)
    
    withProgress(message = 'Running PCA analysis...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = "Converting to matrix...")
        
        values$pca_results <- perform_pca_analysis(
          values$filtered_gl,
          n_pcs = input$n_pcs,
          scale = input$scale_pca,
          center = input$center_pca
        )
        
        incProgress(0.7, detail = "Updating selectors...")
        
        updateSelectInput(session, "pca_x", choices = 1:input$n_pcs, selected = 1)
        updateSelectInput(session, "pca_y", choices = 1:input$n_pcs, selected = 2)
        
        incProgress(1.0, detail = "Complete!")
        showNotification("PCA analysis completed!", type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error in PCA analysis:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # ADDED THIS
  # PCA scores table
  output$pca_scores_table <- renderDT({
    req(values$pca_results, values$pheno_data)
    
    scores_df <- values$pca_results$scores[, 1:min(input$n_pcs, ncol(values$pca_results$scores))]
    scores_df <- round(scores_df, 4)
    
    # Add metadata if available
    if("GEN" %in% colnames(values$pheno_data)) {
      valid_rows <- rownames(values$pca_results$data)
      scores_df$GEN <- values$pheno_data[valid_rows, "GEN"]
    }
    
    datatable(scores_df,
              options = list(pageLength = 10, scrollX = TRUE),
              class = 'display compact',
              caption = "PCA Scores for Individuals")
  })
  
  
  
  # Calculate kinship with withProgress
  observeEvent(input$run_kinship, {
    req(values$filtered_gl)
    
    withProgress(message = 'Calculating kinship matrix...', value = 0, {
      
      tryCatch({
        incProgress(0.5, detail = "Processing...")
        
        values$kinship_matrix <- perform_kinship_analysis(
          values$filtered_gl,
          method = input$kinship_method,
          scale = input$scale_kinship
        )
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Kinship matrix calculated!", type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error calculating kinship:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # ==========================================================================
  # GWAS ANALYSIS OBSERVERS
  # ==========================================================================
  
  # GWAS trait selector
  output$gwas_trait_selector <- renderUI({
    req(values$pheno_data)
    traits <- setdiff(colnames(values$pheno_data), "Genotype")
    selectInput("gwas_trait", "Select Trait for GWAS",
                choices = traits,
                selected = traits[1],
                width = "100%")
  })
  

  
  # Run GWAS with withProgress
  observeEvent(input$run_gwas, {
    req(values$matched_data, input$gwas_trait, values$processed_geno, values$diallel_data)
    
    withProgress(message = 'Running GWAS analysis...', value = 0, {
      
      tryCatch({
        incProgress(0.2, detail = "Preparing data...")
        
        values$gwas_results <- run_single_trait_gwas(
          pheno_data = values$matched_data$pheno,
          geno_matrix = values$matched_data$geno,
          snp_info = values$processed_geno$snp_info,
          trait = input$gwas_trait
        )
        
        incProgress(1.0, detail = "Complete!")
        showNotification("GWAS analysis completed!", type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error in GWAS analysis:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  # =================================================================================================================
  # DIALLEL ANALYSIS OBSERVERS
  # =================================================================================================================
  get_parent_genotypes <- function(parents, geno_matrix, gwas_results) {
    tryCatch({
      cat("Extracting actual genotypes for", length(parents), "parents...\n")
      
      # Check inputs
      if (is.null(parents) || length(parents) == 0) {
        stop("No parents provided")
      }
      if (is.null(geno_matrix) || !is.matrix(geno_matrix)) {
        stop("Genotype matrix must be a matrix")
      }
      if (is.null(gwas_results) || !is.data.frame(gwas_results)) {
        stop("GWAS results must be a data frame")
      }
      
      # Find parents that exist in the genotype matrix
      parents_in_geno <- parents[parents %in% rownames(geno_matrix)]
      
      if (length(parents_in_geno) == 0) {
        cat("No parents found in genotype matrix. Available rows (first 10):\n")
        cat(head(rownames(geno_matrix), 10), "\n")
        return(NULL)
      }
      
      cat("Found", length(parents_in_geno), "/", length(parents), "parents in genotype matrix\n")
      
      # Get SNPs from GWAS results
      gwas_snps <- unique(gwas_results$SNP_ID)
      cat("Unique SNPs in GWAS results:", length(gwas_snps), "\n")
      
      # Check which SNPs are in the genotype matrix
      snps_in_geno <- gwas_snps[gwas_snps %in% colnames(geno_matrix)]
      cat("Common SNPs between GWAS and genotype matrix:", length(snps_in_geno), "\n")
      
      if (length(snps_in_geno) == 0) {
        cat("No common SNPs found. Genotype matrix SNPs (first 10):", head(colnames(geno_matrix), 10), "\n")
        return(NULL)
      }
      
      # Subset genotype matrix to only include common SNPs and parents
      geno_subset <- geno_matrix[parents_in_geno, snps_in_geno, drop = FALSE]
      
      # Convert to long format
      parent_geno_long <- data.frame(
        Parent = character(),
        SNP_ID = character(),
        Genotype = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Create the long format data frame
      for (i in 1:nrow(geno_subset)) {
        parent <- rownames(geno_subset)[i]
        for (j in 1:ncol(geno_subset)) {
          snp <- colnames(geno_subset)[j]
          genotype <- geno_subset[i, j]
          
          parent_geno_long <- rbind(parent_geno_long, data.frame(
            Parent = parent,
            SNP_ID = snp,
            Genotype = genotype,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Add GWAS information
      gwas_info <- gwas_results %>%
        filter(SNP_ID %in% snps_in_geno) %>%
        select(SNP_ID, Chromosome, Position, P_value, Effect, R_squared)
      
      # Merge with genotype data
      result <- parent_geno_long %>%
        left_join(gwas_info, by = "SNP_ID") %>%
        mutate(
          LOD = -log10(P_value),
          PVE = R_squared * 100  # Phenotypic Variance Explained as percentage
        ) %>%
        select(Parent, SNP_ID, Chromosome, Position, Genotype, 
               P_value, Effect, R_squared, PVE, LOD) %>%
        arrange(Parent, Chromosome, Position)
      
      cat("Successfully extracted", nrow(result), "parent-SNP combinations\n")
      return(result)
      
    }, error = function(e) {
      cat("Error extracting parent genotypes:", e$message, "\n")
      return(NULL)
    })
  }
  
  
  
  # Fixed match_parents_to_gwas_actual function without problematic stop() calls
  match_parents_to_gwas_actual <- function(diallel_long, geno_matrix, gwas_results) {
    tryCatch({
      # Check required inputs
      if (is.null(diallel_long)) {
        cat("ERROR: diallel_long is NULL\n")
        return(NULL)
      }
      if (is.null(geno_matrix)) {
        cat("ERROR: geno_matrix is NULL\n")
        return(NULL)
      }
      if (is.null(gwas_results)) {
        cat("ERROR: gwas_results is NULL\n")
        return(NULL)
      }
      
      cat("=== MATCHING PARENTS TO GWAS WITH ACTUAL GENOTYPES ===\n")
      
      # Ensure diallel_long is a data frame
      if (!is.data.frame(diallel_long)) {
        diallel_long <- as.data.frame(diallel_long)
      }
      
      # Convert Parent column to character
      diallel_long$Parent <- as.character(diallel_long$Parent)
      
      # Get unique parents from diallel data
      unique_parents <- unique(diallel_long$Parent)
      cat("Unique parents in diallel data:", length(unique_parents), "\n")
      
      # Ensure genotype matrix has proper row names
      if (is.null(rownames(geno_matrix))) {
        cat("ERROR: Genotype matrix has no row names\n")
        return(NULL)
      }
      
      # Clean row names
      rownames(geno_matrix) <- trimws(as.character(rownames(geno_matrix)))
      
      # Find parents that exist in genotype matrix
      parents_in_geno <- unique_parents[unique_parents %in% rownames(geno_matrix)]
      
      if (length(parents_in_geno) == 0) {
        cat("ERROR: No parents found in genotype matrix\n")
        return(list(
          error = "No parents found in genotype matrix",
          success = FALSE
        ))
      }
      
      cat("Found", length(parents_in_geno), "/", length(unique_parents), "parents in genotype matrix\n")
      
      # Ensure gwas_results is a data frame
      if (!is.data.frame(gwas_results)) {
        cat("ERROR: gwas_results is not a data frame\n")
        return(list(
          error = "GWAS results is not a data frame",
          success = FALSE
        ))
      }
      
      # Clean SNP_ID column
      gwas_results$SNP_ID <- as.character(gwas_results$SNP_ID)
      
      # Get SNPs from GWAS results
      gwas_snps <- unique(gwas_results$SNP_ID)
      cat("Unique SNPs in GWAS results:", length(gwas_snps), "\n")
      
      # Ensure genotype matrix has column names
      if (is.null(colnames(geno_matrix))) {
        cat("ERROR: Genotype matrix has no column names\n")
        return(list(
          error = "Genotype matrix has no column names",
          success = FALSE
        ))
      }
      
      # Clean column names
      colnames(geno_matrix) <- trimws(as.character(colnames(geno_matrix)))
      
      # Check which SNPs are in the genotype matrix
      snps_in_geno <- gwas_snps[gwas_snps %in% colnames(geno_matrix)]
      cat("Common SNPs between GWAS and genotype matrix:", length(snps_in_geno), "\n")
      
      if (length(snps_in_geno) == 0) {
        cat("WARNING: No common SNPs found\n")
        return(list(
          error = "No common SNPs found between GWAS and genotype matrix",
          success = FALSE
        ))
      }
      
      # Subset genotype matrix
      geno_subset <- geno_matrix[parents_in_geno, snps_in_geno, drop = FALSE]
      cat("Subset genotype matrix dimensions:", dim(geno_subset), "\n")
      
      # Convert to long format
      parent_geno_long <- data.frame(
        Parent = rep(rownames(geno_subset), each = ncol(geno_subset)),
        SNP_ID = rep(colnames(geno_subset), times = nrow(geno_subset)),
        Genotype = as.vector(t(geno_subset)),
        stringsAsFactors = FALSE
      )
      
      cat("Created long format with", nrow(parent_geno_long), "rows\n")
      
      # Prepare GWAS information
      gwas_info <- data.frame(
        SNP_ID = as.character(gwas_results$SNP_ID),
        Chromosome = as.character(gwas_results$Chromosome),
        Position = as.numeric(gwas_results$Position),
        P_value = as.numeric(gwas_results$P_value),
        Effect = as.numeric(gwas_results$Effect),
        R_squared = as.numeric(gwas_results$R_squared),
        stringsAsFactors = FALSE
      )
      
      # Filter to only include SNPs we have
      gwas_info <- gwas_info[gwas_info$SNP_ID %in% snps_in_geno, ]
      cat("GWAS info filtered to", nrow(gwas_info), "SNPs\n")
      
      # Merge data
      result <- merge(parent_geno_long, gwas_info, by = "SNP_ID", all.x = TRUE)
      cat("Merge completed. Result has", nrow(result), "rows\n")
      
      # Calculate additional columns
      result$LOD <- -log10(result$P_value)
      result$PVE <- result$R_squared * 100
      
      # Reorder columns
      result <- result[, c("Parent", "SNP_ID", "Chromosome", "Position", "Genotype", 
                           "P_value", "Effect", "R_squared", "PVE", "LOD")]
      
      # Sort the data
      result <- result[order(result$Parent, result$Chromosome, result$Position), ]
      
      cat("Successfully extracted", nrow(result), "parent-SNP combinations\n")
      
      # Count significant combinations
      significant_count <- sum(result$P_value < 0.05, na.rm = TRUE)
      cat("Found", significant_count, "significant parent-SNP combinations\n")
      
      # Calculate summary statistics
      parent_summary <- result %>%
        group_by(Parent) %>%
        summarise(
          N_SNPs = n(),
          Mean_Genotype = mean(Genotype, na.rm = TRUE),
          Missing_Genotype = sum(is.na(Genotype)),
          Min_P_value = min(P_value, na.rm = TRUE),
          Max_Effect = max(abs(Effect), na.rm = TRUE),
          Mean_PVE = mean(PVE, na.rm = TRUE),
          Max_LOD = max(LOD, na.rm = TRUE),
          .groups = "drop"
        )
      
      cat("Created parent summary with", nrow(parent_summary), "parents\n")
      
      # Get top SNP for each parent
      top_snps_per_parent <- result %>%
        group_by(Parent) %>%
        arrange(desc(LOD)) %>%
        slice(1) %>%
        select(Parent, Top_SNP = SNP_ID, Top_SNP_Chr = Chromosome, 
               Top_SNP_Pos = Position, Top_SNP_Genotype = Genotype,
               Top_SNP_LOD = LOD, Top_SNP_Effect = Effect, Top_SNP_PVE = PVE)
      
      # Get parent performance from diallel data
      parent_performance <- diallel_long %>%
        group_by(Parent) %>%
        summarise(
          Mean_GYP = mean(GYP, na.rm = TRUE),
          Mean_GYR = mean(GYR, na.rm = TRUE),
          N_Crosses = n_distinct(CrossID, na.rm = TRUE),
          .groups = "drop"
        )
      
      cat("Parent performance calculated for", nrow(parent_performance), "parents\n")
      
      # Create summary table
      summary_table <- parent_summary %>%
        left_join(top_snps_per_parent, by = "Parent") %>%
        left_join(parent_performance, by = "Parent") %>%
        select(Parent, Mean_GYP, Mean_GYR, N_Crosses, N_SNPs, Mean_Genotype,
               Missing_Genotype, Top_SNP, Top_SNP_Chr, Top_SNP_Pos, 
               Top_SNP_Genotype, Top_SNP_LOD, Top_SNP_Effect, Top_SNP_PVE,
               Min_P_value, Max_Effect, Mean_PVE, Max_LOD) %>%
        arrange(desc(Mean_GYP))
      
      cat("Summary table created with", nrow(summary_table), "rows\n")
      
      return(list(
        detailed_genotypes = result,
        summary_table = summary_table,
        parents_analyzed = unique(result$Parent),
        snps_analyzed = unique(result$SNP_ID),
        parent_count = length(unique(result$Parent)),
        snp_count = length(unique(result$SNP_ID)),
        significant_combinations = significant_count,
        success = TRUE
      ))
      
    }, error = function(e) {
      # Just return error information without using stop()
      error_msg <- paste("ERROR in match_parents_to_gwas_actual:", e$message)
      cat(error_msg, "\n")
      return(list(
        error = error_msg,
        success = FALSE
      ))
    })
  }
  
  # Function to analyze parent genotypes at significant SNPs
  analyze_parent_significant_snps <- function(parent_genotypes, p_threshold = 0.05) {
    tryCatch({
      req(parent_genotypes)
      
      # Filter for significant SNPs
      significant_genotypes <- parent_genotypes %>%
        filter(P_value < p_threshold)
      
      if (nrow(significant_genotypes) == 0) {
        cat("No significant SNPs found at threshold p <", p_threshold, "\n")
        return(NULL)
      }
      
      cat("Found", nrow(significant_genotypes), "significant parent-SNP combinations\n")
      
      # Analyze genotype patterns at significant SNPs
      genotype_patterns <- significant_genotypes %>%
        group_by(SNP_ID, Chromosome, Position) %>%
        summarise(
          N_Parents = n(),
          N_AA = sum(Genotype == 0, na.rm = TRUE),
          N_Aa = sum(Genotype == 1, na.rm = TRUE),
          N_aa = sum(Genotype == 2, na.rm = TRUE),
          Mean_Genotype = mean(Genotype, na.rm = TRUE),
          Mean_PVE = mean(PVE, na.rm = TRUE),
          Parents_AA = paste(Parent[Genotype == 0], collapse = ", "),
          Parents_Aa = paste(Parent[Genotype == 1], collapse = ", "),
          Parents_aa = paste(Parent[Genotype == 2], collapse = ", "),
          .groups = "drop"
        ) %>%
        arrange(Chromosome, Position)
      
      return(genotype_patterns)
      
    }, error = function(e) {
      cat("Error analyzing significant SNPs:", e$message, "\n")
      return(NULL)
    })
  }
  
  
  # Function to create a parent-SNP summary table
  create_parent_snp_summary <- function(diallel_long, gwas_results) {
    tryCatch({
      req(diallel_long, gwas_results)
      
      # Get unique parents
      unique_parents <- unique(diallel_long$Parent)
      
      # Get top SNPs from GWAS
      top_snps <- gwas_results %>%
        arrange(P_value) %>%
        head(min(20, nrow(gwas_results))) %>%
        mutate(
          LOD = -log10(P_value),
          PVE = R_squared * 100  # Phenotypic Variance Explained as percentage
        ) %>%
        select(SNP_ID, Chromosome, Position, P_value, LOD, Effect, PVE)
      
      # Create summary table
      summary_table <- data.frame(
        Parent = character(),
        Top_SNP = character(),
        Top_SNP_Chr = character(),
        Top_SNP_Pos = numeric(),
        Top_SNP_LOD = numeric(),
        Top_SNP_Effect = numeric(),
        Top_SNP_PVE = numeric(),
        Mean_Performance = numeric(),
        N_Crosses = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Calculate parent performance from diallel data
      parent_performance <- diallel_long %>%
        group_by(Parent) %>%
        summarise(
          Mean_GYP = mean(GYP, na.rm = TRUE),
          Mean_GYR = mean(GYR, na.rm = TRUE),
          N_Crosses = n_distinct(CrossID),
          .groups = "drop"
        )
      
      # Assign top SNP to each parent (simplified - randomly assign for demonstration)
      for (i in seq_along(unique_parents)) {
        parent <- unique_parents[i]
        
        # Get parent performance
        perf <- parent_performance[parent_performance$Parent == parent, ]
        
        # Assign a SNP (cycling through top SNPs)
        snp_idx <- ((i - 1) %% nrow(top_snps)) + 1
        snp <- top_snps[snp_idx, ]
        
        summary_table <- bind_rows(
          summary_table,
          data.frame(
            Parent = parent,
            Top_SNP = snp$SNP_ID,
            Top_SNP_Chr = snp$Chromosome,
            Top_SNP_Pos = snp$Position,
            Top_SNP_LOD = round(snp$LOD, 2),
            Top_SNP_Effect = round(snp$Effect, 4),
            Top_SNP_PVE = round(snp$PVE, 2),
            Mean_Performance = ifelse(nrow(perf) > 0, round(perf$Mean_GYP, 2), NA),
            N_Crosses = ifelse(nrow(perf) > 0, perf$N_Crosses, 0),
            stringsAsFactors = FALSE
          )
        )
      }
      
      return(summary_table)
      
    }, error = function(e) {
      cat("Error creating parent-SNP summary:", e$message, "\n")
      return(NULL)
    })
  }
  
  
  
  # ==========================================================================
  # OUTPUT RENDERERS - PARENT-SNP MATCHING TAB
  # ==========================================================================
  
  # Summary statistics output
  output$parent_snp_summary_stats <- renderPrint({
    req(values$parent_snp_actual)
    
    result <- values$parent_snp_actual
    
    cat("=== PARENT-SNP MATCHING SUMMARY ===\n\n")
    cat("Analysis Status: ✅ Completed Successfully\n")
    cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("\n")
    cat("Data Summary:\n")
    cat("  Parents analyzed:", result$parent_count, "\n")
    cat("  SNPs analyzed:", result$snp_count, "\n")
    cat("  Total combinations:", nrow(result$detailed_genotypes), "\n")
    cat("  Significant combinations (p < 0.05):", result$significant_combinations, "\n")
    cat("\n")
    cat("Parent Coverage:\n")
    if (!is.null(result$summary_table)) {
      cat("  Parents with genotype data:", nrow(result$summary_table), "\n")
      cat("  Parents missing from genotype matrix:", 
          length(unique(values$diallel_long$Parent)) - result$parent_count, "\n")
    }
  })
  
  # Current p-threshold display
  output$current_p_threshold <- renderText({
    input$snp_p_threshold
  })
  
  # Summary table renderer
  output$parent_snp_summary_table <- renderDT({
    req(values$parent_snp_actual$summary_table)
    
    # Format the table for better display
    formatted_table <- values$parent_snp_actual$summary_table %>%
      mutate(
        Mean_GYP = round(Mean_GYP, 2),
        Mean_GYR = round(Mean_GYR, 2),
        Mean_Genotype = round(Mean_Genotype, 3),
        Missing_Genotype = as.integer(Missing_Genotype),
        Top_SNP_LOD = round(Top_SNP_LOD, 2),
        Top_SNP_Effect = round(Top_SNP_Effect, 4),
        Top_SNP_PVE = round(Top_SNP_PVE, 2),
        Min_P_value = format(Min_P_value, scientific = TRUE, digits = 3),
        Max_Effect = round(Max_Effect, 4),
        Mean_PVE = round(Mean_PVE, 2),
        Max_LOD = round(Max_LOD, 2)
      )
    
    datatable(
      formatted_table,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Parent-SNP Association Summary",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'Mean_GYP',
        background = styleColorBar(formatted_table$Mean_GYP, 'lightgreen'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      ) %>%
      formatStyle(
        'Top_SNP_LOD',
        background = styleColorBar(formatted_table$Top_SNP_LOD, 'lightblue'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  # Detailed genotypes table renderer
  output$parent_snp_detailed_table <- renderDT({
    req(values$parent_snp_actual$detailed_genotypes)
    
    # Create a display-friendly version
    display_df <- values$parent_snp_actual$detailed_genotypes %>%
      mutate(
        P_value = format(P_value, scientific = TRUE, digits = 3),
        Effect = round(Effect, 4),
        R_squared = round(R_squared, 4),
        PVE = round(PVE, 2),
        LOD = round(LOD, 2),
        Genotype = round(Genotype, 2),
        Significant = P_value < 0.05
      )
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Detailed Parent-SNP Genotypes",
      extensions = 'Buttons',
      filter = 'top'
    ) #%>%
      #formatStyle(
      #  'Significant',
      #  target = 'row',
      #  backgroundColor = styleEqual(c(TRUE, FALSE), c('orange', 'blue'))
      #)
  })
  
  # Significant SNP patterns table
  output$parent_snp_patterns_table <- renderDT({
    req(values$parent_genotype_patterns)
    
    datatable(
      values$parent_genotype_patterns,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Genotype Patterns at Significant SNPs"
    )
  })
  
  # ==========================================================================
  # VISUALIZATION OUTPUT RENDERERS
  # ==========================================================================
  
  # Parent genotype heatmap
  output$parent_genotype_heatmap <- renderPlot({
    req(values$parent_snp_actual$detailed_genotypes)
    
    # Prepare data for heatmap
    heatmap_data <- values$parent_snp_actual$detailed_genotypes %>%
      filter(P_value < 0.05) %>%  # Only significant SNPs
      group_by(Parent, Chromosome) %>%
      summarise(
        Mean_Genotype = mean(Genotype, na.rm = TRUE),
        N_SNPs = n(),
        .groups = 'drop'
      )
    
    if (nrow(heatmap_data) > 0) {
      ggplot(heatmap_data, aes(x = Chromosome, y = Parent, fill = Mean_Genotype)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 1, name = "Mean Genotype") +
        labs(
          title = "Parent Genotype Heatmap (Significant SNPs Only)",
          x = "Chromosome",
          y = "Parent"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No significant SNPs to display", size = 6) +
        theme_void()
    }
  })
  
  # Parent-SNP scatter plot
  output$parent_snp_scatter <- renderPlot({
    req(values$parent_snp_actual$detailed_genotypes)
    
    # Create scatter plot of effect size vs P-value
    plot_data <- values$parent_snp_actual$detailed_genotypes %>%
      filter(!is.na(Effect), !is.na(P_value))
    
    if (nrow(plot_data) > 0) {
      ggplot(plot_data, aes(x = Effect, y = -log10(P_value))) +
        geom_point(alpha = 0.3, color = "steelblue") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        facet_wrap(~Parent, scales = "free") +
        labs(
          title = "Effect Size vs P-value by Parent",
          x = "Effect Size",
          y = "-log10(P-value)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    }
  })
  
  # Parent performance plot
  output$parent_performance_plot <- renderPlot({
    req(values$parent_snp_actual$summary_table)
    
    ggplot(values$parent_snp_actual$summary_table, 
           aes(x = reorder(Parent, Mean_GYP), y = Mean_GYP)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_text(aes(label = round(Mean_GYP, 1)), 
                vjust = -0.3, size = 3.5) +
      labs(
        title = "Parent Performance (Mean GYP)",
        x = "Parent",
        y = "Mean GYP"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  })
  
  # SNP effect distribution plot
  output$snp_effect_distribution <- renderPlot({
    req(values$parent_snp_actual$detailed_genotypes)
    
    ggplot(values$parent_snp_actual$detailed_genotypes, 
           aes(x = Effect, fill = as.factor(round(Genotype)))) +
      geom_density(alpha = 0.5) +
      facet_wrap(~round(Genotype), ncol = 1) +
      labs(
        title = "SNP Effect Distribution by Genotype",
        x = "Effect Size",
        y = "Density",
        fill = "Genotype"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      )
  })
  
  # ==========================================================================
  # DOWNLOAD HANDLERS
  # ==========================================================================
  
  # Download detailed genotypes
  output$download_detailed_genotypes <- downloadHandler(
    filename = function() {
      paste0("parent_snp_detailed_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$parent_snp_actual$detailed_genotypes)
      write.csv(values$parent_snp_actual$detailed_genotypes, file, row.names = FALSE)
    }
  )
  
  # Download SNP patterns
  output$download_snp_patterns <- downloadHandler(
    filename = function() {
      paste0("snp_patterns_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$parent_genotype_patterns)
      write.csv(values$parent_genotype_patterns, file, row.names = FALSE)
    }
  )
  
  # Download all results
  output$download_all_results <- downloadHandler(
    filename = function() {
      paste0("parent_snp_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(values$parent_snp_actual)
      
      # Create temporary directory
      temp_dir <- tempfile("parent_snp_")
      dir.create(temp_dir)
      
      # Save all data frames
      write.csv(values$parent_snp_actual$summary_table,
                file.path(temp_dir, "parent_summary.csv"),
                row.names = FALSE)
      
      write.csv(values$parent_snp_actual$detailed_genotypes,
                file.path(temp_dir, "detailed_genotypes.csv"),
                row.names = FALSE)
      
      if (!is.null(values$parent_genotype_patterns)) {
        write.csv(values$parent_genotype_patterns,
                  file.path(temp_dir, "snp_patterns.csv"),
                  row.names = FALSE)
      }
      
      # Create summary report
      summary_report <- paste(
        "PARENT-SNP ANALYSIS REPORT",
        "==========================",
        "",
        paste("Analysis Date:", Sys.time()),
        paste("Generated by: Plant Breeding Analysis Platform v4.0"),
        "",
        "ANALYSIS SUMMARY",
        paste("Parents analyzed:", values$parent_snp_actual$parent_count),
        paste("SNPs analyzed:", values$parent_snp_actual$snp_count),
        paste("Total combinations:", nrow(values$parent_snp_actual$detailed_genotypes)),
        paste("Significant combinations (p < 0.05):", values$parent_snp_actual$significant_combinations),
        "",
        "FILES INCLUDED",
        "1. parent_summary.csv - Summary statistics for each parent",
        "2. detailed_genotypes.csv - All parent-SNP combinations",
        "3. snp_patterns.csv - Genotype patterns at significant SNPs",
        "",
        sep = "\n"
      )
      
      writeLines(summary_report, file.path(temp_dir, "README.txt"))
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "results.zip")
      zip(zip_file, files = list.files(temp_dir, full.names = TRUE))
      
      # Copy to download location
      file.copy(zip_file, file)
    }
  )
  
  # ============================================================================
  # MULTI-TRAIT ANALYSIS SECTION
  # ============================================================================
  
  # Multi-trait trait selector
  output$multi_trait_selector <- renderUI({
    req(values$pheno_data)
    
    tryCatch({
      traits <- setdiff(colnames(values$pheno_data), "Genotype")
      
      # Filter numeric traits only
      numeric_traits <- traits[sapply(traits, function(x) is.numeric(values$pheno_data[[x]]))]
      
      if (length(numeric_traits) == 0) {
        return(div(class = "alert alert-warning",
                   "No numeric traits found in phenotype data"))
      }
      
      selectizeInput("multi_traits", "Select Traits for Analysis",
                     choices = numeric_traits,
                     multiple = TRUE,
                     selected = numeric_traits[1:min(2, length(numeric_traits))],
                     options = list(maxItems = 10),
                     width = "100%")
      
    }, error = function(e) {
      return(div(class = "alert alert-danger",
                 paste("Error creating trait selector:", e$message)))
    })
  })
  
  
  # Multi-trait summary output
  output$multi_trait_summary <- renderPrint({
    req(values$multi_trait_results)
    
    cat("=== MULTI-TRAIT ANALYSIS SUMMARY ===\n\n")
    
    # Basic info
    if (!is.null(values$multi_trait_results$n_traits)) {
      cat("Number of traits analyzed:", values$multi_trait_results$n_traits, "\n")
    }
    
    if (!is.null(values$multi_trait_results$traits_analyzed)) {
      cat("Traits:", paste(values$multi_trait_results$traits_analyzed, collapse = ", "), "\n")
    }
    
    if (!is.null(values$multi_trait_results$samples_used)) {
      cat("Common samples with complete data:", length(values$multi_trait_results$samples_used), "\n")
    }
    
    # Combined results
    if (!is.null(values$multi_trait_results$combined_results)) {
      combined <- values$multi_trait_results$combined_results
      
      cat("\n=== COMBINED RESULTS ===\n")
      cat("Total SNPs in combined analysis:", nrow(combined), "\n")
      
      if ("N_Traits" %in% colnames(combined)) {
        cat("SNPs affecting 1 trait:", sum(combined$N_Traits == 1, na.rm = TRUE), "\n")
        cat("SNPs affecting 2+ traits:", sum(combined$N_Traits >= 2, na.rm = TRUE), "\n")
        cat("Maximum traits affected by single SNP:", max(combined$N_Traits, na.rm = TRUE), "\n")
      }
      
      if ("Significant_FDR" %in% colnames(combined)) {
        sig_fdr <- sum(combined$Significant_FDR, na.rm = TRUE)
        cat("Significant SNPs (FDR < 0.05):", sig_fdr, "\n")
      }
    }
  })
  
  
  
  # Trait correlation plot
  output$trait_cor_plot <- renderPlot({
    req(values$multi_trait_results)
    
    # Check if we have correlation data
    if (is.null(values$multi_trait_results$trait_correlations)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No correlation data available\nRun multi-trait analysis first", 
                        size = 6) +
               theme_void())
    }
    
    # Extract correlation matrix
    cor_matrix <- NULL
    if (is.matrix(values$multi_trait_results$trait_correlations)) {
      cor_matrix <- values$multi_trait_results$trait_correlations
    } else if (is.list(values$multi_trait_results$trait_correlations)) {
      if (!is.null(values$multi_trait_results$trait_correlations$correlations)) {
        cor_matrix <- values$multi_trait_results$trait_correlations$correlations
      }
    }
    
    if (is.null(cor_matrix) || nrow(cor_matrix) == 0 || ncol(cor_matrix) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No valid correlation matrix available", 
                        size = 6) +
               theme_void())
    }
    
    # Get trait names from the matrix
    trait_names <- rownames(cor_matrix)
    
    # Create a better heatmap using ggplot2
    tryCatch({
      # Convert to long format
      cor_long <- as.data.frame(cor_matrix)
      cor_long$Trait1 <- rownames(cor_long)
      
      cor_long_melted <- melt(cor_long, id.vars = "Trait1")
      colnames(cor_long_melted) <- c("Trait1", "Trait2", "Correlation")
      
      # Reorder factors for better visualization
      cor_long_melted$Trait1 <- factor(cor_long_melted$Trait1, levels = trait_names)
      cor_long_melted$Trait2 <- factor(cor_long_melted$Trait2, levels = trait_names)
      
      # Create the heatmap
      p <- ggplot(cor_long_melted, aes(x = Trait1, y = Trait2, fill = Correlation)) +
        geom_tile(color = "white", size = 1) +
        geom_text(aes(label = round(Correlation, 2)), 
                  color = ifelse(abs(cor_long_melted$Correlation) > 0.5, "white", "black"),
                  size = 3.5, fontface = "bold") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                             midpoint = 0, limit = c(-1, 1),
                             name = "Correlation\nCoefficient") +
        labs(
          title = "Phenotypic Correlation Matrix",
          subtitle = "Correlations between observed trait values",
          x = "",
          y = ""
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.grid = element_blank()
        ) +
        coord_fixed()
      
      return(p)
      
    }, error = function(e) {
      cat("Error in trait correlation plot:", e$message, "\n")
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("Error creating correlation plot:", e$message), 
                        size = 6) +
               theme_void())
    })
  })
  
  
  # Genetic correlation plot
  # Genetic correlation plot - CORRECTED
  output$genetic_cor_plot <- renderPlot({
    req(values$multi_trait_results)
    
    # Check if we have genetic correlation data
    if (is.null(values$multi_trait_results$genetic_correlations)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "Genetic correlations not available\nRun multi-trait analysis first", 
                        size = 6) +
               theme_void())
    }
    
    # Extract genetic correlation matrix
    cor_matrix <- values$multi_trait_results$genetic_correlations
    
    if (is.null(cor_matrix) || nrow(cor_matrix) == 0 || ncol(cor_matrix) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No valid genetic correlation matrix available", 
                        size = 6) +
               theme_void())
    }
    
    # Create the heatmap
    melted_cor <- melt(cor_matrix)
    
    ggplot(melted_cor, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, limit = c(-1, 1), space = "Lab",
                           name = "Genetic\nCorrelation") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid = element_blank()
      ) +
      labs(
        title = "Genetic Correlation Matrix",
        subtitle = "Correlations between genetic components of traits",
        x = "",
        y = ""
      ) +
      coord_fixed()
  })
  
  
  
  # Pleiotropy table
  output$pleiotropy_table <- renderDT({
    req(values$multi_trait_results)
    
    if (is.null(values$multi_trait_results$combined_results) || 
        nrow(values$multi_trait_results$combined_results) == 0) {
      return(datatable(
        data.frame(Message = "No combined results available. Run multi-trait analysis first."),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Check if N_Traits column exists
    if (!"N_Traits" %in% colnames(values$multi_trait_results$combined_results)) {
      return(datatable(
        data.frame(Message = "No pleiotropy data (N_Traits column not found)"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Filter for SNPs affecting multiple traits
    pleiotropic_df <- values$multi_trait_results$combined_results
    
    # Only show SNPs affecting >1 trait
    if (any(pleiotropic_df$N_Traits > 1, na.rm = TRUE)) {
      pleiotropic_df <- pleiotropic_df[pleiotropic_df$N_Traits > 1, ]
    } else {
      return(datatable(
        data.frame(Message = "No SNPs found affecting more than 1 trait"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    if (nrow(pleiotropic_df) == 0) {
      return(datatable(
        data.frame(Message = "No pleiotropic SNPs found (affecting >1 trait)"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Format the data
    display_df <- pleiotropic_df
    
    # Format columns
    if ("Min_P_value" %in% colnames(display_df)) {
      display_df$Min_P_value <- format(display_df$Min_P_value, scientific = TRUE, digits = 3)
    }
    
    if ("P_adjusted" %in% colnames(display_df)) {
      display_df$P_adjusted <- format(display_df$P_adjusted, scientific = TRUE, digits = 3)
    }
    
    if ("Max_Effect" %in% colnames(display_df)) {
      display_df$Max_Effect <- round(display_df$Max_Effect, 4)
    }
    
    # Select columns to show
    display_cols <- c("SNP_ID", "Chromosome", "Position", "Min_P_value", "N_Traits")
    
    # Add optional columns if they exist
    optional_cols <- c("P_adjusted", "Max_Effect", "Significant_FDR")
    for (col in optional_cols) {
      if (col %in% colnames(display_df)) {
        display_cols <- c(display_cols, col)
      }
    }
    
    display_df <- display_df[, display_cols, drop = FALSE]
    
    # Order by N_Traits (descending) then by P-value (ascending)
    display_df <- display_df[order(-display_df$N_Traits, display_df$Min_P_value), ]
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Pleiotropic SNPs (Affecting Multiple Traits)",
      extensions = 'Buttons',
      filter = 'top'
    )
  })
  
  # Pleiotropy visualization plot
  output$pleiotropy_plot_output <- renderPlot({
    req(values$multi_trait_results)
    
    if (is.null(values$multi_trait_results$combined_results) || 
        nrow(values$multi_trait_results$combined_results) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No data available\nRun multi-trait analysis first", 
                        size = 6) +
               theme_void())
    }
    
    # Get traits analyzed
    traits_analyzed <- values$multi_trait_results$traits_analyzed
    
    if (is.null(traits_analyzed) || length(traits_analyzed) == 0) {
      # Try to extract from the results
      if (!is.null(values$multi_trait_results$trait_gwas_results)) {
        traits_analyzed <- names(values$multi_trait_results$trait_gwas_results)
      } else {
        return(ggplot() +
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No trait information available", 
                          size = 6) +
                 theme_void())
      }
    }
    
    # Create a comprehensive visualization
    # Option 1: Trait correlation heatmap (if we have correlation data)
    if (!is.null(values$multi_trait_results$trait_correlations)) {
      
      # Extract correlation matrix
      cor_matrix <- values$multi_trait_results$trait_correlations
      
      # If it's not a matrix, try to extract it
      if (is.list(cor_matrix) && !is.null(cor_matrix$correlations)) {
        cor_matrix <- cor_matrix$correlations
      }
      
      if (is.matrix(cor_matrix) && nrow(cor_matrix) > 0) {
        # Create correlation heatmap
        cor_df <- as.data.frame(cor_matrix)
        cor_df$Trait1 <- rownames(cor_df)
        
        # Melt for plotting
        cor_long <- melt(cor_df, id.vars = "Trait1")
        colnames(cor_long) <- c("Trait1", "Trait2", "Correlation")
        
        # Create the heatmap
        p <- ggplot(cor_long, aes(x = Trait1, y = Trait2, fill = Correlation)) +
          geom_tile(color = "white") +
          geom_text(aes(label = round(Correlation, 2)), 
                    color = "black", size = 3.5) +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                               midpoint = 0, limit = c(-1, 1),
                               name = "Correlation") +
          labs(
            title = "Trait Correlation Matrix (Pleiotropy Analysis)",
            subtitle = paste("Analyzing", length(traits_analyzed), "traits"),
            x = "",
            y = ""
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            legend.position = "right"
          ) +
          coord_fixed()
        
        return(p)
      }
    }
    
    # Option 2: If no correlation matrix, show trait distribution
    # Create a visualization of SNPs by traits affected
    if ("N_Traits" %in% colnames(values$multi_trait_results$combined_results)) {
      
      plot_data <- values$multi_trait_results$combined_results
      
      # Create the distribution plot
      p <- ggplot(plot_data, aes(x = N_Traits)) +
        geom_bar(fill = "steelblue", alpha = 0.7, color = "black") +
        geom_text(stat = 'count', aes(label = ..count..), 
                  vjust = -0.5, size = 4) +
        labs(
          title = "Pleiotropy: Number of Traits Affected per SNP",
          subtitle = paste("Traits analyzed:", paste(traits_analyzed, collapse = ", ")),
          x = "Number of Traits Affected",
          y = "Number of SNPs"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid.minor = element_blank()
        ) +
        scale_x_continuous(breaks = 1:max(plot_data$N_Traits, na.rm = TRUE)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
      
      return(p)
    }
    
    # Option 3: Show summary statistics in text
    ggplot() +
      annotate("text", x = 0.5, y = 0.7, 
               label = "Pleiotropy Analysis Summary", 
               size = 6, fontface = "bold") +
      annotate("text", x = 0.5, y = 0.6, 
               label = paste("Traits analyzed:", length(traits_analyzed)), 
               size = 5) +
      annotate("text", x = 0.5, y = 0.5, 
               label = paste(traits_analyzed, collapse = ", "), 
               size = 4) +
      annotate("text", x = 0.5, y = 0.4, 
               label = paste("Total SNPs:", nrow(values$multi_trait_results$combined_results)), 
               size = 4) +
      theme_void() +
      xlim(0, 1) + ylim(0, 1)
  })
  
  
  
  # Multi-trait statistics
  output$multi_trait_stats <- renderPrint({
    req(values$multi_trait_results)
    
    cat("=== MULTI-TRAIT ANALYSIS STATISTICS ===\n\n")
    
    # Basic information
    if (!is.null(values$multi_trait_results$n_traits)) {
      cat("Number of traits analyzed:", values$multi_trait_results$n_traits, "\n")
    }
    
    if (!is.null(values$multi_trait_results$traits_analyzed)) {
      cat("Traits:", paste(values$multi_trait_results$traits_analyzed, collapse = ", "), "\n")
    }
    
    if (!is.null(values$multi_trait_results$samples_used)) {
      cat("Common samples with complete data:", length(values$multi_trait_results$samples_used), "\n")
    }
    
    cat("\n")
    
    # Individual trait statistics
    if (!is.null(values$multi_trait_results$trait_gwas_results)) {
      cat("=== INDIVIDUAL TRAIT GWAS RESULTS ===\n\n")
      
      for (trait in names(values$multi_trait_results$trait_gwas_results)) {
        trait_result <- values$multi_trait_results$trait_gwas_results[[trait]]
        
        if (!is.null(trait_result) && is.data.frame(trait_result)) {
          cat("Trait:", trait, "\n")
          cat("  Total SNPs tested:", nrow(trait_result), "\n")
          
          # Count significant SNPs
          if ("P_adjusted" %in% colnames(trait_result)) {
            sig_fdr <- sum(trait_result$P_adjusted < 0.05, na.rm = TRUE)
            cat("  Significant SNPs (FDR < 0.05):", sig_fdr, "\n")
          }
          
          if ("P_value" %in% colnames(trait_result)) {
            min_p <- min(trait_result$P_value, na.rm = TRUE)
            cat("  Minimum P-value:", format(min_p, scientific = TRUE, digits = 3), "\n")
          }
          
          cat("\n")
        }
      }
    }
    
    # Combined results statistics
    if (!is.null(values$multi_trait_results$combined_results)) {
      cat("=== COMBINED RESULTS ===\n\n")
      
      combined <- values$multi_trait_results$combined_results
      
      cat("Total SNPs in combined analysis:", nrow(combined), "\n")
      
      if ("N_Traits" %in% colnames(combined)) {
        cat("SNPs affecting 1 trait:", sum(combined$N_Traits == 1, na.rm = TRUE), "\n")
        cat("SNPs affecting 2 traits:", sum(combined$N_Traits == 2, na.rm = TRUE), "\n")
        cat("SNPs affecting 3+ traits:", sum(combined$N_Traits >= 3, na.rm = TRUE), "\n")
        cat("Maximum traits affected by single SNP:", max(combined$N_Traits, na.rm = TRUE), "\n")
      }
      
      if ("Significant_FDR" %in% colnames(combined)) {
        sig_fdr <- sum(combined$Significant_FDR, na.rm = TRUE)
        cat("Significant SNPs in combined analysis (FDR < 0.05):", sig_fdr, "\n")
      }
      
      if ("Min_P_value" %in% colnames(combined)) {
        min_p <- min(combined$Min_P_value, na.rm = TRUE)
        cat("Minimum P-value across all traits:", format(min_p, scientific = TRUE, digits = 3), "\n")
      }
    }
  })
  
  
  
  
  
  check_cross_snp_prerequisites <- function(values) {
    prerequisites <- list(
      diallel_data = !is.null(values$diallel_data),
      multi_trait_results = !is.null(values$multi_trait_results),
      processed_geno = !is.null(values$processed_geno),
      snp_matrix = !is.null(values$processed_geno$snp_matrix),
      snp_info = !is.null(values$processed_geno$snp_info)
    )
    
    missing <- names(prerequisites)[!unlist(prerequisites)]
    
    if (length(missing) > 0) {
      return(list(
        ready = FALSE,
        missing = missing,
        message = paste("Missing:", paste(missing, collapse = ", "))
      ))
    }
    
    return(list(ready = TRUE, missing = NULL, message = "All data available"))
  }
  
  
  
  
  
  
  # Run multi-trait analysis
  observeEvent(input$run_multi_trait, {
    req(values$matched_data, values$processed_geno, input$multi_traits)
    
    # Check if at least 2 traits selected
    if (length(input$multi_traits) < 2) {
      showNotification("Please select at least 2 traits for multi-trait analysis", 
                       type = "warning", duration = 3)
      return()
    }
    
    withProgress(message = 'Running multi-trait analysis...', value = 0, {
      
      tryCatch({
        incProgress(0.2, detail = "Preparing data...")
        
        # Prepare data
        pheno_data <- values$matched_data$pheno
        geno_matrix <- values$matched_data$geno
        snp_info <- values$processed_geno$snp_info
        
        # Find common samples with complete data for all selected traits
        common_samples <- rownames(geno_matrix)
        for (trait in input$multi_traits) {
          na_mask <- is.na(pheno_data[[trait]])
          common_samples <- intersect(common_samples, pheno_data$Genotype[!na_mask])
        }
        
        if (length(common_samples) < 5) {
          showNotification(paste("Only", length(common_samples), 
                                 "common samples with complete data. Need at least 5."),
                           type = "error", duration = 5)
          return()
        }
        
        incProgress(0.4, detail = "Running analysis...")
        
        # Run multi-trait analysis
        multi_results <- run_multi_trait_gwas(
          pheno_data = pheno_data,
          geno_matrix = geno_matrix,
          snp_info = snp_info,
          traits = input$multi_traits,
          method = input$multi_trait_method,
          include_correlations = input$include_correlations,
          calculate_pleiotropy = input$calculate_pleiotropy
        )
        
        
        # Store results
        values$multi_trait_results <- multi_results
        
        
        # Create unified data formats
        if (!is.null(multi_results$trait_gwas_results)) {
          incProgress(0.6, detail = "Creating unified formats...")
          
          # Create unified dataframe
          values$multi_trait_unified_df <- create_unified_gwas_tidy(multi_results$trait_gwas_results)
          
          # Create wide format
          values$multi_trait_unified_wide <- create_wide_gwas_format(multi_results$trait_gwas_results)
          
          # Create enhanced unified format
          values$multi_trait_enhanced_unified_df <- create_enhanced_unified_df(
            multi_results$trait_gwas_results,
            multi_results$trait_correlations
          )
          
          # Create complete integrated format
          values$multi_trait_complete_integrated_df <- create_complete_integrated_df(multi_results)
        }
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Multi-trait analysis completed successfully!", 
                         type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error in multi-trait analysis:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  
  
  observeEvent(input$run_cross_snp_all_traits, {
    withProgress(message = 'Running Cross-SNP Analysis...', value = 0, {
      incProgress(0.3, detail = "Preparing data...")
      
      results <- run_cross_snp_analysis_all_traits_app(
        values, 
        p_threshold = input$cross_snp_p_threshold
      )
      
      print(results)
      values$cross_snp_df <- results
      
      write.csv(results, "cross_data_with_LR_markers.csv") 
      
      incProgress(0.8, detail = "Creating summaries...")
      
      if (!is.null(results)) {
        output$cross_snp_table <- renderDT({
          datatable(
            head(results, 100),  # Show first 100 rows
            options = list(
              pageLength = 10,
              scrollX = TRUE,
              dom = 'Bfrtip',
              buttons = c('copy', 'csv', 'excel', 'pdf')
            ),
            class = 'display compact',
            rownames = FALSE,
            caption = "Cross-SNP Analysis Results (First 100 rows)"
          )
        })
        
        output$cross_snp_summary <- renderPrint({
          cat("=== Cross-SNP Analysis Summary ===\n\n")
          cat("Total rows:", nrow(results), "\n")
          cat("Unique crosses:", length(unique(results$Cross)), "\n")
          cat("Unique SNPs:", length(unique(results$SNP_ID)), "\n")
          cat("Traits analyzed:", paste(unique(results$Trait), collapse = ", "), "\n")
          
          if ("Significant" %in% colnames(results)) {
            cat("Significant associations (p <", input$cross_snp_p_threshold, "):", 
                sum(results$Significant, na.rm = TRUE), "\n")
          }
          
          if ("PVE" %in% colnames(results)) {
            cat("Average PVE:", round(mean(as.numeric(results$PVE), na.rm = TRUE), 2), "%\n")
          }
        })
      }
      
      incProgress(1.0, detail = "Complete!")
    })
  })
  
  
  # ==========================================================================
  # CROSS-SNP ANALYSIS: MERGING DIALLEL AND GWAS DATA
  # ==========================================================================
  
  #' Create cross-specific SNP dataframe
  #' This function merges diallel cross data with GWAS results
  #' to create cross-specific SNP statistics
  create_cross_snp_dataframe <- function(diallel_data, gwas_unified_df, matched_geno = NULL) {
    tryCatch({
      require(dplyr)
      require(tidyr)
      
      # Check required data
      if (is.null(diallel_data) || nrow(diallel_data) == 0) {
        stop("No diallel data available")
      }
      
      if (is.null(gwas_unified_df) || nrow(gwas_unified_df) == 0) {
        stop("No GWAS data available")
      }
      
      cat("Creating cross-SNP analysis...\n")
      cat("Diallel entries:", nrow(diallel_data), "\n")
      cat("GWAS rows:", nrow(gwas_unified_df), "\n")
      cat("Unique SNPs in GWAS:", n_distinct(gwas_unified_df$SNP_ID), "\n")
      
      # Extract cross information from diallel data
      cross_df <- diallel_data %>%
        filter(CrossType == "Cross") %>%
        select(Cross = ACC., Parent1, Parent2) %>%
        distinct() %>%
        filter(!is.na(Parent1) & !is.na(Parent2))
      
      if (nrow(cross_df) == 0) {
        stop("No valid crosses found in diallel data")
      }
      
      cat("Found", nrow(cross_df), "unique crosses\n")
      
      # Get unique parents
      all_parents <- unique(c(cross_df$Parent1, cross_df$Parent2))
      cat("Unique parents in crosses:", length(all_parents), "\n")
      
      # Prepare GWAS data
      gwas_summary <- gwas_unified_df %>%
        mutate(
          LOD = -log10(P_value),
          PVE = ifelse("R_squared" %in% colnames(.), R_squared * 100, NA)
        ) %>%
        select(
          SNP_ID, Chromosome, Position, Trait, 
          P_value, LOD, PVE, Effect
        )
      
      # Get traits measured in diallel
      diallel_traits <- c("GYP", "GYR", "DTF", "DTM", "PL", "NPP", "NSP", "X100SW", "TW")
      diallel_traits <- diallel_traits[diallel_traits %in% colnames(diallel_data)]
      
      cat("Diallel traits:", paste(diallel_traits, collapse = ", "), "\n")
      
      # Filter GWAS data for diallel traits
      gwas_diallel <- gwas_summary %>%
        filter(Trait %in% diallel_traits)
      
      if (nrow(gwas_diallel) == 0) {
        stop("No GWAS results for diallel traits")
      }
      
      cat("GWAS SNPs for diallel traits:", nrow(gwas_diallel), "\n")
      
      # Create cross-SNP dataframe with cross-specific calculations
      cross_snp_list <- list()
      
      for (i in 1:nrow(cross_df)) {
        cross_info <- cross_df[i, ]
        parent1 <- as.character(cross_info$Parent1)
        parent2 <- as.character(cross_info$Parent2)
        
        # Check if parents are in genotype matrix
        if (!is.null(matched_geno)) {
          # Get parent genotypes if available
          parent1_in_geno <- parent1 %in% rownames(matched_geno)
          parent2_in_geno <- parent2 %in% rownames(matched_geno)
          
          if (parent1_in_geno && parent2_in_geno) {
            # Get SNPs present in both parents and GWAS
            common_snps <- intersect(
              colnames(matched_geno),
              unique(gwas_diallel$SNP_ID)
            )
            
            if (length(common_snps) > 0) {
              # Get parent genotypes for common SNPs
              parent1_geno <- matched_geno[parent1, common_snps]
              parent2_geno <- matched_geno[parent2, common_snps]
              
              # Identify polymorphic SNPs (parents have different genotypes)
              polymorphic_mask <- parent1_geno != parent2_geno
              polymorphic_snps <- common_snps[polymorphic_mask]
              
              # Filter GWAS data to polymorphic SNPs
              cross_gwas <- gwas_diallel %>%
                filter(SNP_ID %in% polymorphic_snps) %>%
                mutate(
                  Cross = cross_info$Cross,
                  Parent1 = parent1,
                  Parent2 = parent2,
                  Parent1_Genotype = parent1_geno[SNP_ID],
                  Parent2_Genotype = parent2_geno[SNP_ID],
                  Expected_Genotype = (parent1_geno[SNP_ID] + parent2_geno[SNP_ID]) / 2,
                  Parent1_In_Geno = TRUE,
                  Parent2_In_Geno = TRUE,
                  Polymorphic = TRUE
                )
              
              cross_snp_list[[i]] <- cross_gwas
              
            } else {
              # No common SNPs, create empty dataframe for this cross
              cat("No common SNPs for cross", cross_info$Cross, "\n")
              cross_snp_list[[i]] <- data.frame(
                Cross = cross_info$Cross,
                Parent1 = parent1,
                Parent2 = parent2,
                SNP_ID = NA,
                Trait = NA,
                P_value = NA,
                LOD = NA,
                PVE = NA,
                Effect = NA,
                Parent1_Genotype = NA,
                Parent2_Genotype = NA,
                Expected_Genotype = NA,
                Parent1_In_Geno = TRUE,
                Parent2_In_Geno = TRUE,
                Polymorphic = FALSE,
                stringsAsFactors = FALSE
              )
            }
          } else {
            # One or both parents not in genotype matrix
            cross_snp_list[[i]] <- data.frame(
              Cross = cross_info$Cross,
              Parent1 = parent1,
              Parent2 = parent2,
              SNP_ID = NA,
              Trait = NA,
              P_value = NA,
              LOD = NA,
              PVE = NA,
              Effect = NA,
              Parent1_Genotype = NA,
              Parent2_Genotype = NA,
              Expected_Genotype = NA,
              Parent1_In_Geno = parent1_in_geno,
              Parent2_In_Geno = parent2_in_geno,
              Polymorphic = FALSE,
              stringsAsFactors = FALSE
            )
          }
        } else {
          # No genotype matrix available, use all GWAS SNPs
          cross_gwas <- gwas_diallel %>%
            mutate(
              Cross = cross_info$Cross,
              Parent1 = parent1,
              Parent2 = parent2,
              Parent1_Genotype = NA,
              Parent2_Genotype = NA,
              Expected_Genotype = NA,
              Parent1_In_Geno = FALSE,
              Parent2_In_Geno = FALSE,
              Polymorphic = NA  # Unknown without genotype data
            )
          
          cross_snp_list[[i]] <- cross_gwas
        }
      }
      
      # Combine all cross data
      cross_snp_data <- bind_rows(cross_snp_list)
      
      # Add derived columns
      cross_snp_data <- cross_snp_data %>%
        mutate(
          # Flag if SNP is significant
          Significant = P_value < 0.05,
          # SNP effect direction
          Effect_Direction = ifelse(Effect > 0, "Positive", 
                                    ifelse(Effect < 0, "Negative", "Zero")),
          # Create a unique cross-SNP ID
          Cross_SNP_ID = paste(Cross, SNP_ID, Trait, sep = "_"),
          # Calculate cross-specific expected effect
          # (effect weighted by expected genotype if available)
          Expected_Effect = ifelse(!is.na(Expected_Genotype), 
                                   Effect * Expected_Genotype, 
                                   NA)
        )
      
      # Filter out rows with no SNP data
      cross_snp_data <- cross_snp_data %>%
        filter(!is.na(SNP_ID))
      
      cat("Cross-SNP dataframe created with", nrow(cross_snp_data), "rows\n")
      cat("Crosses with SNP data:", n_distinct(cross_snp_data$Cross), "\n")
      cat("Polymorphic SNPs:", sum(cross_snp_data$Polymorphic == TRUE, na.rm = TRUE), "\n")
      
      return(cross_snp_data)
      
    }, error = function(e) {
      cat("Error creating cross-SNP dataframe:", e$message, "\n")
      cat("Traceback:\n")
      print(traceback())
      return(NULL)
    })
  }
  
  
  
  # Observer to run cross-SNP analysis
  observeEvent(input$run_cross_snp_analysis, {
    req(values$multi_trait_unified_df, values$diallel_data)
    
    showModal(modalDialog(
      title = "Merging Cross and SNP Data",
      tagList(
        tags$div(
          style = "text-align: center;",
          icon("spinner", class = "fa-spin fa-2x"),
          tags$h4("Analyzing crosses and SNPs..."),
          tags$div(
            id = "cross_snp_details",
            style = "margin-top: 10px; font-size: 12px; color: #666;"
          )
        )
      ),
      footer = NULL,
      easyClose = FALSE,
      size = "s"
    ))
    
    tryCatch({
      shinyjs::html("cross_snp_details", "Step 1/4: Preparing data...")
      Sys.sleep(0.5)
      
      # Get genotype matrix if available
      geno_matrix <- NULL
      if (!is.null(values$matched_data$geno)) {
        geno_matrix <- values$matched_data$geno
      }
      
      shinyjs::html("cross_snp_details", "Step 2/4: Merging cross and SNP data...")
      Sys.sleep(0.5)
      
      # Use the enhanced unified dataframe or the simple one
      gwas_data <- values$multi_trait_enhanced_unified_df
      if (is.null(gwas_data)) {
        gwas_data <- values$multi_trait_unified_df
      }
      
      # Create merged dataframe
      cross_snp_df <- create_cross_snp_dataframe(
        diallel_data = values$diallel_data,
        gwas_unified_df = gwas_data,
        matched_geno = geno_matrix
      )
      
      if (is.null(cross_snp_df) || nrow(cross_snp_df) == 0) {
        removeModal()
        showNotification("Failed to create cross-SNP analysis or no data available", type = "error")
        return()
      }
      
      shinyjs::html("cross_snp_details", "Step 3/4: Calculating statistics...")
      Sys.sleep(0.5)
      
      # Store results
      values$cross_snp_analysis <- list(
        data = cross_snp_df,
        n_crosses = n_distinct(cross_snp_df$Cross),
        n_snps = n_distinct(cross_snp_df$SNP_ID),
        n_traits = n_distinct(cross_snp_df$Trait),
        timestamp = Sys.time()
      )
      
      # Create summary tables
      if (!is.null(cross_snp_df)) {
        # Summary by cross (only polymorphic SNPs)
        values$cross_summary <- cross_snp_df %>%
          filter(Polymorphic == TRUE) %>%
          group_by(Cross, Parent1, Parent2) %>%
          summarise(
            N_SNPs = n_distinct(SNP_ID),
            N_Polymorphic_SNPs = n_distinct(SNP_ID[Polymorphic == TRUE]),
            N_Significant = sum(Significant, na.rm = TRUE),
            Mean_LOD = mean(LOD, na.rm = TRUE),
            Max_LOD = max(LOD, na.rm = TRUE),
            Mean_PVE = mean(PVE, na.rm = TRUE),
            Max_PVE = max(PVE, na.rm = TRUE),
            Mean_Effect = mean(Effect, na.rm = TRUE),
            Max_Effect = max(abs(Effect), na.rm = TRUE),
            N_Positive_Effects = sum(Effect > 0, na.rm = TRUE),
            N_Negative_Effects = sum(Effect < 0, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Summary by SNP across crosses
        values$snp_summary_cross <- cross_snp_df %>%
          filter(Polymorphic == TRUE) %>%
          group_by(SNP_ID, Chromosome, Position, Trait) %>%
          summarise(
            N_Crosses = n_distinct(Cross),
            Crosses_List = paste(unique(Cross), collapse = ";"),
            Mean_P_value = mean(P_value, na.rm = TRUE),
            Mean_LOD = mean(LOD, na.rm = TRUE),
            Min_LOD = min(LOD, na.rm = TRUE),
            Max_LOD = max(LOD, na.rm = TRUE),
            Mean_PVE = mean(PVE, na.rm = TRUE),
            Mean_Effect = mean(Effect, na.rm = TRUE),
            Effect_Consistency = ifelse(
              all(Effect > 0, na.rm = TRUE), "Always Positive",
              ifelse(all(Effect < 0, na.rm = TRUE), "Always Negative", "Mixed")
            ),
            .groups = "drop"
          )
      }
      
      shinyjs::html("cross_snp_details", "Step 4/4: Finalizing...")
      Sys.sleep(0.5)
      
      removeModal()
      
      showNotification(
        HTML(paste(
          "✅ Cross-SNP analysis completed!<br>",
          "Crosses analyzed:", values$cross_snp_analysis$n_crosses, "<br>",
          "SNPs analyzed:", values$cross_snp_analysis$n_snps, "<br>",
          "Traits:", values$cross_snp_analysis$n_traits, "<br>",
          "Polymorphic SNPs:", sum(cross_snp_df$Polymorphic == TRUE, na.rm = TRUE)
        )),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Error in cross-SNP analysis:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS FOR CROSS-SNP ANALYSIS
  # ==========================================================================
  
  # Cross-SNP detailed table
  output$cross_snp_table <- renderDT({
    req(values$cross_snp_analysis$data)
    
    tryCatch({
      df <- values$cross_snp_analysis$data
      
      # Select and format columns for display
      display_df <- df %>%
        select(
          Cross, Parent1, Parent2, 
          SNP_ID, Chromosome, Position, Trait,
          P_value, LOD, PVE, Effect,
          Parent1_Genotype, Parent2_Genotype, Expected_Genotype,
          Polymorphic, Significant, Effect_Direction
        ) %>%
        mutate(
          P_value = format(P_value, scientific = TRUE, digits = 3),
          LOD = round(LOD, 2),
          PVE = round(PVE, 2),
          Effect = round(Effect, 4),
          Parent1_Genotype = round(Parent1_Genotype, 2),
          Parent2_Genotype = round(Parent2_Genotype, 2),
          Expected_Genotype = round(Expected_Genotype, 2)
        )
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf'),
          order = list(list(which(colnames(display_df) == "LOD"), 'desc'))
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Cross-SNP Analysis: Detailed View",
        extensions = 'Buttons',
        filter = 'top'
      ) %>%
        formatStyle(
          'Polymorphic',
          target = 'row',
          backgroundColor = styleEqual(c(TRUE, FALSE), c('#e6ffe6', '#ffe6e6'))
        ) %>%
        formatStyle(
          'LOD',
          background = styleColorBar(range(display_df$LOD, na.rm = TRUE), 'lightgreen')
        )
      
    }, error = function(e) {
      cat("Error rendering cross-SNP table:", e$message, "\n")
      return(datatable(
        data.frame(Error = paste("Error:", e$message)),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    })
  })
  
  
  
  # Cross summary table
  output$cross_summary_table <- renderDT({
    req(values$cross_summary)
    
    tryCatch({
      display_df <- values$cross_summary %>%
        mutate(
          Mean_LOD = round(Mean_LOD, 2),
          Max_LOD = round(Max_LOD, 2),
          Mean_PVE = round(Mean_PVE, 2),
          Max_PVE = round(Max_PVE, 2),
          Mean_Effect = round(Mean_Effect, 4),
          Max_Effect = round(Max_Effect, 4)
        )
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          order = list(list(which(colnames(display_df) == "Max_LOD"), 'desc'))
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Cross Summary: SNP Statistics by Cross (Polymorphic SNPs Only)",
        extensions = 'Buttons'
      ) %>%
        formatStyle(
          'Max_LOD',
          background = styleColorBar(range(display_df$Max_LOD, na.rm = TRUE), 'lightblue')
        ) %>%
        formatStyle(
          'N_Polymorphic_SNPs',
          background = styleColorBar(range(display_df$N_Polymorphic_SNPs, na.rm = TRUE), 'lightgreen')
        )
      
    }, error = function(e) {
      cat("Error in cross summary table:", e$message, "\n")
      return(NULL)
    })
  })
  
  
  
  # SNP summary table (across crosses)
  output$snp_summary_cross_table <- renderDT({
    req(values$snp_summary_cross)
    
    tryCatch({
      display_df <- values$snp_summary_cross %>%
        mutate(
          Mean_P_value = format(Mean_P_value, scientific = TRUE, digits = 3),
          Mean_LOD = round(Mean_LOD, 2),
          Min_LOD = round(Min_LOD, 2),
          Max_LOD = round(Max_LOD, 2),
          Mean_PVE = round(Mean_PVE, 2),
          Mean_Effect = round(Mean_Effect, 4)
        )
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          order = list(list(which(colnames(display_df) == "Max_LOD"), 'desc'))
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "SNP Summary: Performance Across Crosses (Polymorphic SNPs Only)",
        extensions = 'Buttons'
      ) %>%
        formatStyle(
          'Effect_Consistency',
          target = 'row',
          backgroundColor = styleEqual(
            c("Always Positive", "Always Negative", "Mixed"),
            c('#e6ffe6', '#ffe6e6', '#fff2e6')
          )
        )
      
    }, error = function(e) {
      cat("Error in SNP summary table:", e$message, "\n")
      return(NULL)
    })
  })
  
  
  # Top crosses by LOD score
  output$top_crosses_table <- renderDT({
    req(values$cross_snp_analysis$data)
    
    tryCatch({
      top_crosses <- values$cross_snp_analysis$data %>%
        group_by(Cross, Parent1, Parent2) %>%
        summarise(
          Top_SNP = SNP_ID[which.max(LOD)][1],
          Top_SNP_Chr = Chromosome[which.max(LOD)][1],
          Top_SNP_Pos = Position[which.max(LOD)][1],
          Max_LOD = max(LOD, na.rm = TRUE),
          Mean_LOD = mean(LOD, na.rm = TRUE),
          Max_PVE = max(PVE, na.rm = TRUE),
          N_Significant = sum(Significant, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(desc(Max_LOD)) %>%
        head(20) %>%
        mutate(
          Max_LOD = round(Max_LOD, 2),
          Mean_LOD = round(Mean_LOD, 2),
          Max_PVE = round(Max_PVE, 2)
        )
      
      datatable(
        top_crosses,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Top 20 Crosses by Maximum LOD Score"
      )
      
    }, error = function(e) {
      cat("Error in top crosses table:", e$message, "\n")
      return(NULL)
    })
  })
  
  
  
  # Cross-SNP analysis summary
  output$cross_snp_summary <- renderPrint({
    if (is.null(values$cross_snp_analysis)) {
      cat("=== CROSS-SNP ANALYSIS ===\n\n")
      cat("No analysis performed yet.\n")
      cat("Click 'Run Cross-SNP Analysis' to merge diallel and GWAS data.\n")
      return()
    }
    
    analysis <- values$cross_snp_analysis
    df <- analysis$data
    
    cat("=== CROSS-SNP ANALYSIS SUMMARY ===\n\n")
    cat("Analysis timestamp:", format(analysis$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Total cross-SNP combinations:", nrow(df), "\n")
    cat("Unique crosses:", analysis$n_crosses, "\n")
    cat("Unique SNPs:", analysis$n_snps, "\n")
    cat("Traits analyzed:", analysis$n_traits, "\n")
    cat("Polymorphic SNPs:", sum(df$Polymorphic == TRUE, na.rm = TRUE), "\n")
    cat("Crosses with genotype data:", 
        n_distinct(df$Cross[df$Parent1_In_Geno & df$Parent2_In_Geno]), "\n\n")
    
    if (nrow(df) > 0) {
      # Statistics for polymorphic SNPs only
      poly_df <- df %>% filter(Polymorphic == TRUE)
      
      if (nrow(poly_df) > 0) {
        cat("=== POLYMORPHIC SNPs STATISTICS ===\n")
        cat("LOD Score Statistics:\n")
        cat("  Mean LOD:", round(mean(poly_df$LOD, na.rm = TRUE), 2), "\n")
        cat("  Max LOD:", round(max(poly_df$LOD, na.rm = TRUE), 2), "\n")
        cat("  Min LOD:", round(min(poly_df$LOD, na.rm = TRUE), 2), "\n\n")
        
        cat("PVE Statistics:\n")
        if (!all(is.na(poly_df$PVE))) {
          cat("  Mean PVE:", round(mean(poly_df$PVE, na.rm = TRUE), 2), "%\n")
          cat("  Max PVE:", round(max(poly_df$PVE, na.rm = TRUE), 2), "%\n")
          cat("  Min PVE:", round(min(poly_df$PVE, na.rm = TRUE), 2), "%\n\n")
        }
        
        cat("Effect Size Statistics:\n")
        cat("  Mean Effect:", round(mean(poly_df$Effect, na.rm = TRUE), 4), "\n")
        cat("  Max Positive Effect:", round(max(poly_df$Effect, na.rm = TRUE), 4), "\n")
        cat("  Max Negative Effect:", round(min(poly_df$Effect, na.rm = TRUE), 4), "\n\n")
        
        cat("Significance:\n")
        sig_count <- sum(poly_df$Significant, na.rm = TRUE)
        cat("  Significant combinations (p < 0.05):", sig_count, 
            paste0("(", round(sig_count/nrow(poly_df)*100, 1), "%)\n"))
        
        # Top performing cross
        if (!is.null(values$cross_summary) && nrow(values$cross_summary) > 0) {
          top_cross <- values$cross_summary %>%
            arrange(desc(Max_LOD)) %>%
            slice(1)
          
          cat("\n=== TOP PERFORMING CROSS ===\n")
          cat("Cross:", top_cross$Cross, "\n")
          cat("Parents:", top_cross$Parent1, "x", top_cross$Parent2, "\n")
          cat("Max LOD:", round(top_cross$Max_LOD, 2), "\n")
          cat("Mean LOD:", round(top_cross$Mean_LOD, 2), "\n")
          cat("Max PVE:", round(top_cross$Max_PVE, 2), "%\n")
          cat("Polymorphic SNPs:", top_cross$N_Polymorphic_SNPs, "\n")
          cat("Significant SNPs:", top_cross$N_Significant, "\n")
          cat("Positive/Negative Effects:", top_cross$N_Positive_Effects, "/", 
              top_cross$N_Negative_Effects, "\n")
        }
        
        # Top SNP
        if (!is.null(values$snp_summary_cross) && nrow(values$snp_summary_cross) > 0) {
          top_snp <- values$snp_summary_cross %>%
            arrange(desc(Max_LOD)) %>%
            slice(1)
          
          cat("\n=== TOP SNP ACROSS CROSSES ===\n")
          cat("SNP:", top_snp$SNP_ID, "\n")
          cat("Position: Chr", top_snp$Chromosome, ":", top_snp$Position, "\n")
          cat("Trait:", top_snp$Trait, "\n")
          cat("Max LOD:", round(top_snp$Max_LOD, 2), "\n")
          cat("Number of crosses:", top_snp$N_Crosses, "\n")
          cat("Effect consistency:", top_snp$Effect_Consistency, "\n")
        }
      } else {
        cat("\nNo polymorphic SNPs found.\n")
        cat("This may be because:\n")
        cat("1. Genotype data is not available for parents\n")
        cat("2. Parents have identical genotypes at all SNPs\n")
        cat("3. No common SNPs between GWAS and genotype data\n")
      }
    }
  })
  
  
  
  # Plot: LOD distribution by cross
  output$cross_lod_plot <- renderPlot({
    req(values$cross_snp_analysis$data)
    
    tryCatch({
      df <- values$cross_snp_analysis$data
      print(df)
      # Filter to polymorphic SNPs only
      poly_df <- df %>% filter(Polymorphic == TRUE)
      
      if (nrow(poly_df) == 0) {
        return(ggplot() +
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No polymorphic SNPs to plot", size = 6) +
                 theme_void())
      }
      
      # Calculate max LOD per cross
      cross_lod <- poly_df %>%
        group_by(Cross) %>%
        summarise(Max_LOD = max(LOD, na.rm = TRUE)) %>%
        arrange(desc(Max_LOD)) %>%
        head(20)  # Top 20 crosses
      
      ggplot(cross_lod, aes(x = reorder(Cross, Max_LOD), y = Max_LOD)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        coord_flip() +
        labs(
          title = "Top 20 Crosses by Maximum LOD Score (Polymorphic SNPs Only)",
          subtitle = paste("Showing", nrow(cross_lod), "crosses with polymorphic SNPs"),
          x = "Cross",
          y = "Maximum LOD Score"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text.y = element_text(size = 10)
        )
      
    }, error = function(e) {
      cat("Error in cross LOD plot:", e$message, "\n")
      return(NULL)
    })
  })
  
  
  
  # Plot: SNP effect heatmap across crosses
  output$snp_effect_heatmap <- renderPlot({
    req(values$cross_snp_analysis$data)
    
    tryCatch({
      df <- values$cross_snp_analysis$data
      
      # Filter to polymorphic SNPs only
      poly_df <- df %>% filter(Polymorphic == TRUE)
      
      if (nrow(poly_df) == 0) {
        return(ggplot() +
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "No polymorphic SNPs to plot", size = 6) +
                 theme_void())
      }
      
      # Select top SNPs and crosses for visualization
      top_snps <- poly_df %>%
        group_by(SNP_ID) %>%
        summarise(Mean_LOD = mean(LOD, na.rm = TRUE)) %>%
        arrange(desc(Mean_LOD)) %>%
        head(15) %>%
        pull(SNP_ID)
      
      top_crosses <- poly_df %>%
        group_by(Cross) %>%
        summarise(Mean_LOD = mean(LOD, na.rm = TRUE)) %>%
        arrange(desc(Mean_LOD)) %>%
        head(15) %>%
        pull(Cross)
      
      # Create heatmap data
      heatmap_data <- poly_df %>%
        filter(SNP_ID %in% top_snps, Cross %in% top_crosses) %>%
        group_by(Cross, SNP_ID) %>%
        summarise(Mean_Effect = mean(Effect, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = SNP_ID, values_from = Mean_Effect) %>%
        column_to_rownames("Cross")
      
      if (nrow(heatmap_data) == 0 || ncol(heatmap_data) == 0) {
        return(ggplot() +
                 annotate("text", x = 0.5, y = 0.5, 
                          label = "Insufficient data for heatmap", size = 6) +
                 theme_void())
      }
      
      # Create heatmap with ggplot2
      heatmap_long <- heatmap_data %>%
        rownames_to_column("Cross") %>%
        pivot_longer(cols = -Cross, names_to = "SNP", values_to = "Effect")
      
      ggplot(heatmap_long, aes(x = SNP, y = Cross, fill = Effect)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                             midpoint = 0, na.value = "grey90") +
        labs(
          title = "SNP Effect Heatmap Across Top Crosses",
          subtitle = "Polymorphic SNPs Only",
          x = "SNP",
          y = "Cross",
          fill = "Effect Size"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
    }, error = function(e) {
      cat("Error in SNP effect heatmap:", e$message, "\n")
      return(NULL)
    })
  })
  
  
  
  # Download cross-SNP results
  output$download_cross_snp_data <- downloadHandler(
    filename = function() {
      paste0("cross_snp_analysis_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(values$cross_snp_analysis)
      
      # Create temporary directory
      temp_dir <- tempfile("cross_snp_")
      dir.create(temp_dir)
      
      # Save all dataframes
      write.csv(values$cross_snp_analysis$data,
                file.path(temp_dir, "cross_snp_detailed.csv"),
                row.names = FALSE)
      
      if (!is.null(values$cross_summary)) {
        write.csv(values$cross_summary,
                  file.path(temp_dir, "cross_summary.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$snp_summary_cross)) {
        write.csv(values$snp_summary_cross,
                  file.path(temp_dir, "snp_summary_cross.csv"),
                  row.names = FALSE)
      }
      
      # Save unified multi-trait data if available
      if (!is.null(values$multi_trait_unified_df)) {
        write.csv(values$multi_trait_unified_df,
                  file.path(temp_dir, "multi_trait_unified.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$multi_trait_enhanced_unified_df)) {
        write.csv(values$multi_trait_enhanced_unified_df,
                  file.path(temp_dir, "multi_trait_enhanced.csv"),
                  row.names = FALSE)
      }
      
      # Create metadata
      meta_data <- data.frame(
        Parameter = c("Analysis Date", "Total Cross-SNP Rows", "Unique Crosses", 
                      "Unique SNPs", "Traits Analyzed", "Polymorphic SNPs", "Data Source"),
        Value = c(
          as.character(Sys.time()),
          nrow(values$cross_snp_analysis$data),
          values$cross_snp_analysis$n_crosses,
          values$cross_snp_analysis$n_snps,
          values$cross_snp_analysis$n_traits,
          sum(values$cross_snp_analysis$data$Polymorphic == TRUE, na.rm = TRUE),
          "Merged Diallel and Multi-Trait GWAS Data"
        )
      )
      write.csv(meta_data, file.path(temp_dir, "metadata.csv"), row.names = FALSE)
      
      # Create README
      readme_content <- paste(
        "CROSS-SNP ANALYSIS RESULTS",
        "==========================",
        "",
        "Files included:",
        "1. cross_snp_detailed.csv - Detailed cross-SNP combinations",
        "2. cross_summary.csv - Summary statistics per cross (polymorphic SNPs only)",
        "3. snp_summary_cross.csv - SNP statistics across crosses (polymorphic SNPs only)",
        "4. multi_trait_unified.csv - Unified multi-trait GWAS results",
        "5. multi_trait_enhanced.csv - Enhanced multi-trait GWAS results",
        "6. metadata.csv - Analysis parameters",
        "",
        "Column descriptions for cross_snp_detailed.csv:",
        "- Cross: Diallel cross ID",
        "- Parent1, Parent2: Parental lines in the cross",
        "- SNP_ID: SNP identifier",
        "- Chromosome, Position: Genomic coordinates",
        "- Trait: Phenotypic trait",
        "- P_value: GWAS p-value",
        "- LOD: -log10(P_value)",
        "- PVE: Phenotypic variance explained (%)",
        "- Effect: Additive effect size",
        "- Parent1_Genotype, Parent2_Genotype: Genotype values (if available)",
        "- Expected_Genotype: Mid-parent genotype value",
        "- Polymorphic: TRUE if parents have different genotypes at this SNP",
        "- Significant: TRUE if p < 0.05",
        "- Effect_Direction: Positive or Negative effect",
        "",
        "Note: Statistics in cross_summary.csv and snp_summary_cross.csv",
        "are calculated only for polymorphic SNPs (Polymorphic = TRUE)",
        "This ensures cross-specific statistics that vary between crosses.",
        "",
        sep = "\n"
      )
      writeLines(readme_content, file.path(temp_dir, "README.txt"))
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "results.zip")
      zip::zip(zip_file, files = list.files(temp_dir, full.names = TRUE), 
               mode = "cherry-pick")
      
      # Copy to download location
      file.copy(zip_file, file)
      
      # Cleanup
      unlink(temp_dir, recursive = TRUE)
    }
  )
  
  
  
  # Clear cross-SNP results
  observeEvent(input$clear_cross_snp_results, {
    values$cross_snp_analysis <- NULL
    values$cross_summary <- NULL
    values$snp_summary_cross <- NULL
    
    # Reset outputs
    output$cross_snp_table <- renderDT({NULL})
    output$cross_summary_table <- renderDT({NULL})
    output$snp_summary_cross_table <- renderDT({NULL})
    output$cross_snp_summary <- renderPrint({
      cat("=== CROSS-SNP ANALYSIS ===\n\n")
      cat("No results available.\n")
      cat("Click 'Run Cross-SNP Analysis' to analyze.\n")
    })
    
    output$cross_lod_plot <- renderPlot({
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No cross-SNP results available", size = 6) +
        theme_void()
    })
    
    output$snp_effect_heatmap <- renderPlot({
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No cross-SNP results available", size = 6) +
        theme_void()
    })
    
    showNotification("Cross-SNP results cleared!", type = "warning", duration = 3)
  })
  
  
  # Cross-SNP analysis status UI
  output$cross_snp_status_ui <- renderUI({
    if (!is.null(values$cross_snp_analysis)) {
      analysis <- values$cross_snp_analysis
      div(
        class = "alert alert-success",
        icon("check-circle"),
        strong(" ✓ Cross-SNP Analysis Complete"),
        br(),
        paste("Crosses:", analysis$n_crosses),
        br(),
        paste("SNPs:", analysis$n_snps),
        br(),
        paste("Traits:", analysis$n_traits),
        br(),
        paste("Polymorphic SNPs:", 
              sum(values$cross_snp_analysis$data$Polymorphic == TRUE, na.rm = TRUE))
      )
    } else if (!is.null(values$multi_trait_unified_df) && !is.null(values$diallel_data)) {
      div(
        class = "alert alert-info",
        icon("info-circle"),
        strong(" Ready for Cross-SNP Analysis"),
        br(),
        "Multi-trait and diallel data available.",
        br(),
        "Click 'Run Cross-SNP Analysis' to merge datasets."
      )
    } else {
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        strong(" Data Required"),
        br(),
        "Need both multi-trait GWAS results and diallel data.",
        br(),
        "Run multi-trait analysis and load diallel data first."
      )
    }
  })
  
  
  # Unified data status UI
  output$unified_data_status_ui <- renderUI({
    if (!is.null(values$multi_trait_unified_df)) {
      df <- values$multi_trait_unified_df
      div(
        class = "alert alert-success",
        icon("check-circle"),
        strong(" ✓ Unified Multi-Trait Data Ready"),
        br(),
        paste("Traits:", n_distinct(df$Trait)),
        br(),
        paste("SNP-trait combinations:", nrow(df)),
        br(),
        paste("Memory:", format(object.size(df), units = "MB"))
      )
    } else if (!is.null(values$multi_trait_results)) {
      div(
        class = "alert alert-info",
        icon("info-circle"),
        strong(" Multi-Trait Analysis Complete"),
        br(),
        "Click 'Create Unified Format' to organize data."
      )
    } else {
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        strong(" No Multi-Trait Data"),
        br(),
        "Run multi-trait analysis first."
      )
    }
  })
  
  
  
  # Data compatibility check
  output$data_compatibility_check <- renderPrint({
    cat("=== DATA COMPATIBILITY CHECK ===\n\n")
    
    # Check multi-trait data
    if (!is.null(values$multi_trait_unified_df)) {
      cat("Multi-Trait Data:\n")
      cat("  Available: ✓\n")
      cat("  Rows:", nrow(values$multi_trait_unified_df), "\n")
      cat("  Columns:", ncol(values$multi_trait_unified_df), "\n")
      cat("  Traits:", n_distinct(values$multi_trait_unified_df$Trait), "\n")
    } else {
      cat("Multi-Trait Data:\n")
      cat("  Available: ✗\n")
    }
    
    cat("\n")
    
    # Check diallel data
    if (!is.null(values$diallel_data)) {
      cat("Diallel Data:\n")
      cat("  Available: ✓\n")
      cat("  Crosses:", nrow(values$diallel_data[values$diallel_data$CrossType == "Cross", ]), "\n")
      cat("  Parents:", nrow(values$diallel_data[values$diallel_data$CrossType == "Parent", ]), "\n")
    } else {
      cat("Diallel Data:\n")
      cat("  Available: ✗\n")
    }
    
    cat("\n")
    
    # Check genotype data
    if (!is.null(values$matched_data$geno)) {
      cat("Genotype Data:\n")
      cat("  Available: ✓\n")
      cat("  Samples:", nrow(values$matched_data$geno), "\n")
      cat("  SNPs:", ncol(values$matched_data$geno), "\n")
    } else {
      cat("Genotype Data:\n")
      cat("  Available: ✗ (Optional for cross-SNP analysis)\n")
    }
    
    cat("\n")
    
    # Cross-SNP analysis requirements
    if (!is.null(values$multi_trait_unified_df) && !is.null(values$diallel_data)) {
      cat("Cross-SNP Analysis Requirements:\n")
      cat("  ✓ Multi-trait data available\n")
      cat("  ✓ Diallel data available\n")
      
      # Check for common traits
      diallel_traits <- c("GYP", "GYR", "DTF", "DTM", "PL", "NPP", "NSP", "X100SW", "TW")
      diallel_traits <- diallel_traits[diallel_traits %in% colnames(values$diallel_data)]
      
      gwas_traits <- unique(values$multi_trait_unified_df$Trait)
      
      common_traits <- intersect(diallel_traits, gwas_traits)
      
      if (length(common_traits) > 0) {
        cat("  ✓ Common traits:", paste(common_traits, collapse = ", "), "\n")
      } else {
        cat("  ⚠ No common traits between diallel and GWAS data\n")
      }
      
      cat("\nCross-SNP analysis can be run.\n")
    } else {
      cat("Cross-SNP Analysis Requirements:\n")
      if (is.null(values$multi_trait_unified_df)) {
        cat("  ✗ Multi-trait data missing\n")
      }
      if (is.null(values$diallel_data)) {
        cat("  ✗ Diallel data missing\n")
      }
      cat("\nCannot run cross-SNP analysis.\n")
    }
  })
  
  
  
  # ==========================================================================
  # UPDATE DATA STATUS UI
  # ==========================================================================
  
  # Update the data status UI after successful analysis
  observeEvent(values$parent_snp_actual, {
    req(values$parent_snp_actual)
    
    output$parent_snp_data_status_ui <- renderUI({
      result <- values$parent_snp_actual
      
      div(
        class = "alert alert-success",
        icon("check-circle"),
        strong(" ✅ Analysis Complete! "),
        br(),
        paste("Successfully analyzed", result$parent_count, "parents and", 
              result$snp_count, "SNPs"),
        br(),
        paste("Total combinations:", nrow(result$detailed_genotypes)),
        br(),
        paste("Significant hits:", result$significant_combinations),
        br(),
        tags$small(paste("Completed:", format(Sys.time(), "%H:%M:%S")))
      )
    })
  })
  
  # Clear results button functionality
  observeEvent(input$clear_parent_snp_results, {
    values$parent_snp_actual <- NULL
    values$parent_genotype_patterns <- NULL
    
    # Reset all outputs
    output$parent_snp_summary_stats <- renderPrint({
      cat("=== PARENT-SNP MATCHING ===\n\n")
      cat("No results available.\n")
      cat("Click 'Run Parent-SNP Matching' to analyze.\n")
    })
    
    output$parent_snp_summary_table <- renderDT({NULL})
    output$parent_snp_detailed_table <- renderDT({NULL})
    output$parent_snp_patterns_table <- renderDT({NULL})
    
    # Reset data status
    output$parent_snp_data_status_ui <- renderUI({
      div(
        class = "alert alert-warning",
        icon("info-circle"),
        strong(" Data Cleared "),
        br(),
        "Previous results have been cleared.",
        br(),
        "Run analysis to generate new results."
      )
    })
    
    showNotification("Results cleared successfully!", type = "default", duration = 3)
  })
  
  #====================================================================================
  
  # Pivot diallel data to have Parent1 and Parent2 in same column
  pivot_diallel_parents <- function(diallel_data) {
    tryCatch({
      req(diallel_data)
      
      # Separate parents and crosses
      parents <- diallel_data[diallel_data$CrossType == "Parent", ]
      crosses <- diallel_data[diallel_data$CrossType == "Cross", ]
      
      # Create long format for crosses (each cross appears twice)
      crosses_long <- bind_rows(
        # Parent1 entries
        crosses %>%
          mutate(
            Parent = Parent1,
            Parent_Role = "Parent1",
            CrossID = CrossID,
            Other_Parent = Parent2
          ) %>%
          select(-Parent1, -Parent2),
        
        # Parent2 entries
        crosses %>%
          mutate(
            Parent = Parent2,
            Parent_Role = "Parent2",
            CrossID = CrossID,
            Other_Parent = Parent1
          ) %>%
          select(-Parent1, -Parent2)
      )
      
      # Create parent entries (each parent appears once)
      parents_long <- parents %>%
        mutate(
          Parent = ACC.,
          Parent_Role = "Pure_Parent",
          CrossID = NA,
          Other_Parent = NA
        ) %>%
        select(-Parent1, -Parent2)
      
      # Combine both
      diallel_long <- bind_rows(crosses_long, parents_long)
      
      # Reorder columns
      diallel_long <- diallel_long %>%
        select(ACC., CrossID, Parent, Other_Parent, Parent_Role, 
               CrossType, everything())
      
      return(diallel_long)
      
    }, error = function(e) {
      cat("Error pivoting diallel data:", e$message, "\n")
      return(NULL)
    })
  }
  
  
  
  # In the server section
  observeEvent(values$diallel_data, {
    req(values$diallel_data)
    
    tryCatch({
      # Pivot the data
      values$diallel_long <- pivot_diallel_parents(values$diallel_data)
      
      # Also create a unique parent list
      values$all_parents <- unique(values$diallel_long$Parent)
      
      # Create parent summary
      values$parent_summary <- values$diallel_long %>%
        filter(CrossType == "Cross") %>%
        group_by(Parent) %>%
        summarise(
          N_Crosses = n_distinct(CrossID),
          Partner_Parents = paste(unique(Other_Parent), collapse = ", "),
          Mean_GYP = mean(GYP, na.rm = TRUE),
          Mean_GYR = mean(GYR, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(desc(N_Crosses))
      
      cat("Diallel data pivoted successfully\n")
      cat("  Total entries in long format:", nrow(values$diallel_long), "\n")
      cat("  Unique parents:", length(values$all_parents), "\n")
      
    }, error = function(e) {
      cat("Error processing diallel data:", e$message, "\n")
    })
  })
  
  
  
  # Observer for single trait diallel
  observeEvent(input$run_diallel, {
    req(values$diallel_data, input$diallel_trait)
    
    withProgress(message = 'Running diallel analysis...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = "Performing diallel analysis...")
        
        # Run diallel analysis
        results <- perform_diallel_analysis(values$diallel_data, input$diallel_trait)
        
        if (!is.null(results$error)) {
          showNotification(paste("Error:", results$error), type = "error")
          return()
        }
        
        values$diallel_results <- results
        
        # Convert SCA matrix to data frame for easier handling
        if (!is.null(results$sca)) {
          values$sca_df <- sca_matrix_to_df(results$sca)
        }
        
        incProgress(0.6, detail = "Calculating heterosis...")
        
        # Calculate heterosis
        heterosis_data <- calculate_heterosis(
          diallel_data = values$diallel_data,
          cross_data = values$diallel_cross_data,
          parent_data = values$diallel_parent_data,
          trait = input$diallel_trait
        )
        
        values$diallel_heterosis_data <- heterosis_data
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Diallel analysis completed successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Diallel analysis error:", e$message), type = "error")
      })
    })
  })
  
  
  
  # Output renderers for multi-trait diallel results
  output$diallel_all_anova <- renderDT({
    req(values$diallel_all_traits_results)
    
    datatable(
      values$diallel_all_traits_results$summary_tables$anova_summary,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "ANOVA Summary for All Traits"
    )
  })
  
  output$diallel_all_gca <- renderDT({
    req(values$diallel_all_traits_results)
    
    datatable(
      values$diallel_all_traits_results$summary_tables$gca_summary,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "GCA Summary for All Traits"
    )
  })
  
  
  
  # datatable(cor_df,
  #           options = list(pageLength = 10, scrollX = TRUE),
  #           rownames = FALSE) %>%
  #   formatRound('Correlation', 3) %>%
  #   formatStyle('Correlation',
  #               color = styleInterval(c(-0.3, 0.3), c('blue', 'black', 'red')),
  #               fontWeight = styleInterval(c(-0.5, 0.5), c('bold', 'normal', 'bold')))
  # 
  
  output$diallel_long_table <- renderDT({
    req(values$diallel_long)
    
    datatable(
      values$diallel_long %>%
        select(ACC., CrossID, Parent, Other_Parent, Parent_Role, 
               CrossType, GYP, GYR, DTF, DTM),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Diallel Data with Parents Combined in Single Column"
    )
  })
  
  output$parent_summary_table <- renderDT({
    req(values$parent_summary)
    
    datatable(
      values$parent_summary,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Parent Performance Summary",
      extensions = 'Buttons'
    ) %>%
      formatRound(columns = c('Mean_GYP', 'Mean_GYR'), digits = 2)
  })
  
  
  # Create parent network data
  create_parent_network_data <- function(diallel_data) {
    tryCatch({
      req(diallel_data)
      
      # Filter only crosses
      crosses <- diallel_data %>%
        filter(CrossType == "Cross") %>%
        select(Parent1, Parent2, CrossID) %>%
        filter(!is.na(Parent1) & !is.na(Parent2))
      
      if (nrow(crosses) == 0) {
        return(NULL)
      }
      
      # Create edge list (connections between parents)
      edges <- crosses %>%
        group_by(Parent1, Parent2) %>%
        summarise(
          weight = n(),  # Number of crosses between these parents
          width = n() * 2,  # Edge width proportional to number of crosses
          title = paste0("Crosses: ", n(), 
                         "<br>", paste(unique(CrossID), collapse = "<br>")),
          .groups = "drop"
        ) %>%
        rename(from = Parent1, to = Parent2)
      
      # Create node list (parents)
      all_parents <- unique(c(edges$from, edges$to))
      
      # Calculate node importance (degree centrality)
      node_degree <- data.frame(
        parent = all_parents,
        degree = 0
      )
      
      for (i in seq_along(all_parents)) {
        parent <- all_parents[i]
        degree <- sum(edges$from == parent) + sum(edges$to == parent)
        node_degree$degree[i] <- degree
      }
      
      # Create nodes data frame
      nodes <- data.frame(
        id = all_parents,
        label = all_parents,
        title = paste0("Parent: ", all_parents, "<br>Connections: ", 
                       node_degree$degree),
        value = node_degree$degree,  # Node size
        group = "Parents",
        stringsAsFactors = FALSE
      )
      
      return(list(nodes = nodes, edges = edges))
      
    }, error = function(e) {
      cat("Error creating network data:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Reactive value for network data
  observe({
    req(values$diallel_data)
    
    tryCatch({
      values$network_data <- create_parent_network_data(values$diallel_data)
      
      if (!is.null(values$network_data)) {
        cat("Network data created successfully\n")
        cat("  Nodes:", nrow(values$network_data$nodes), "\n")
        cat("  Edges:", nrow(values$network_data$edges), "\n")
      }
    }, error = function(e) {
      cat("Error in network data observer:", e$message, "\n")
    })
  })
  
  # Render the network
  output$parent_network <- renderVisNetwork({
    req(values$network_data)
    
    tryCatch({
      nodes <- values$network_data$nodes
      edges <- values$network_data$edges
      
      # Create the network
      visNetwork(nodes, edges) %>%
        visNodes(
          shape = input$node_shape,
          size = input$node_size,
          color = list(
            background = "#97C2FC",
            border = "#2B7CE9",
            highlight = list(
              background = "#FFC107",
              border = "#FF9800"
            )
          ),
          font = list(
            size = 18,
            color = "#343434"
          )
        ) %>%
        visEdges(
          arrows = "to",
          smooth = TRUE,
          color = list(
            color = "#848484",
            highlight = "#FF5722"
          ),
          width = "width"
        ) %>%
        visOptions(
          highlightNearest = list(
            enabled = TRUE,
            degree = 1,
            hover = TRUE
          ),
          nodesIdSelection = list(
            enabled = TRUE,
            style = 'width: 200px; height: 26px;'
          )
        ) %>%
        visInteraction(
          navigationButtons = TRUE,
          keyboard = TRUE,
          dragNodes = TRUE,
          dragView = TRUE,
          zoomView = TRUE
        ) %>%
        visPhysics(
          enabled = input$physics,
          stabilization = TRUE
        ) %>%
        visLayout(randomSeed = 123) %>%
        visLegend(
          position = "right",
          useGroups = TRUE
        )
    }, error = function(e) {
      cat("Error rendering network:", e$message, "\n")
      return(NULL)
    })
  })
  
  # Network statistics
  output$network_stats <- renderPrint({
    req(values$network_data)
    
    nodes <- values$network_data$nodes
    edges <- values$network_data$edges
    
    cat("=== PARENT NETWORK STATISTICS ===\n\n")
    cat("Network Structure:\n")
    cat("  Total parents (nodes):", nrow(nodes), "\n")
    cat("  Total connections (edges):", nrow(edges), "\n")
    cat("  Network density:", round(nrow(edges) / (nrow(nodes) * (nrow(nodes) - 1) / 2), 3), "\n")
    
    # Calculate centrality measures
    if (nrow(nodes) > 0) {
      # Most connected parents
      top_parents <- nodes[order(-nodes$value), ][1:min(5, nrow(nodes)), ]
      
      cat("\nTop 5 Most Connected Parents:\n")
      for (i in 1:nrow(top_parents)) {
        cat(paste0("  ", i, ". ", top_parents$id[i], 
                   " (", top_parents$value[i], " connections)\n"))
      }
      
      # Calculate clustering coefficient (simplified)
      unique_connections <- length(unique(c(edges$from, edges$to)))
      possible_connections <- nrow(nodes) * (nrow(nodes) - 1) / 2
      clustering_coef <- ifelse(possible_connections > 0, 
                                nrow(edges) / possible_connections, 0)
      
      cat("\nNetwork Properties:\n")
      cat("  Average connections per parent:", round(mean(nodes$value), 2), "\n")
      cat("  Maximum connections:", max(nodes$value), "\n")
      cat("  Minimum connections:", min(nodes$value), "\n")
      cat("  Clustering coefficient:", round(clustering_coef, 3), "\n")
    }
  })
  
  # Download network as HTML
  output$download_network <- downloadHandler(
    filename = function() {
      paste0("parent_network_", Sys.Date(), ".html")
    },
    content = function(file) {
      req(values$network_data)
      
      tryCatch({
        # Create a temporary HTML file
        temp_html <- tempfile(fileext = ".html")
        
        # Create the network
        network <- visNetwork(values$network_data$nodes, 
                              values$network_data$edges) %>%
          visNodes(size = input$node_size, shape = input$node_shape) %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
          visPhysics(enabled = input$physics)
        
        # Save as HTML
        visSave(network, file = temp_html, selfcontained = TRUE)
        
        # Copy to download location
        file.copy(temp_html, file)
      }, error = function(e) {
        showNotification(paste("Error exporting network:", e$message), 
                         type = "error")
      })
    }
  )
  
  # Generate network button
  observeEvent(input$generate_network, {
    req(values$diallel_data)
    
    withProgress(message = 'Generating network...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = "Creating network data...")
        
        values$network_data <- create_parent_network_data(values$diallel_data)
        
        incProgress(0.6, detail = "Rendering visualization...")
        
        # Trigger network update
        output$parent_network <- renderVisNetwork({
          req(values$network_data)
          
          visNetwork(values$network_data$nodes, values$network_data$edges) %>%
            visNodes(shape = input$node_shape, size = input$node_size) %>%
            visEdges(arrows = "to", smooth = TRUE) %>%
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
            visPhysics(enabled = input$physics) %>%
            visLayout(randomSeed = 123)
        })
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Network generated successfully!", type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error generating network:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  
  output$diallel_all_sca <- renderDT({
    req(values$diallel_all_traits_results)
    
    datatable(
      values$diallel_all_traits_results$summary_tables$sca_summary,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "SCA Summary for All Traits"
    )
  })
  
  output$diallel_all_heterosis <- renderDT({
    req(values$diallel_all_traits_results)
    
    datatable(
      values$diallel_all_traits_results$summary_tables$heterosis_summary,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Heterosis Summary for All Traits"
    )
  })
  
  
  
  output$diallel_traits_cor_plot <- renderPlotly({
    req(values$diallel_all_traits_results)
    
    # Extract trait means for correlation
    gca_summary <- values$diallel_all_traits_results$summary_tables$gca_summary
    
    if (!is.null(gca_summary) && nrow(gca_summary) > 0) {
      # Reshape to wide format
      gca_wide <- gca_summary %>%
        select(Parent, Trait, GCA) %>%
        pivot_wider(names_from = Trait, values_from = GCA)
      
      # Calculate correlation
      cor_matrix <- cor(gca_wide[, -1], use = "pairwise.complete.obs")
      
      # Create heatmap
      p <- plot_ly(
        x = colnames(cor_matrix),
        y = rownames(cor_matrix),
        z = cor_matrix,
        type = "heatmap",
        colorscale = "RdBu",
        zmin = -1, zmax = 1
      ) %>%
        layout(
          title = "Trait Correlations based on GCA",
          xaxis = list(title = ""),
          yaxis = list(title = "")
        )
      
      return(p)
    }
  })
  
  # Download handler for diallel summary
  output$download_diallel_summary <- downloadHandler(
    filename = function() {
      paste0("diallel_summary_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(values$diallel_all_traits_results)
      
      # Create temporary directory
      temp_dir <- tempdir()
      
      # Save summary tables
      if (!is.null(values$diallel_all_traits_results$summary_tables$anova_summary)) {
        write.csv(values$diallel_all_traits_results$summary_tables$anova_summary,
                  file.path(temp_dir, "anova_summary.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$diallel_all_traits_results$summary_tables$gca_summary)) {
        write.csv(values$diallel_all_traits_results$summary_tables$gca_summary,
                  file.path(temp_dir, "gca_summary.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$diallel_all_traits_results$summary_tables$sca_summary)) {
        write.csv(values$diallel_all_traits_results$summary_tables$sca_summary,
                  file.path(temp_dir, "sca_summary.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$diallel_all_traits_results$summary_tables$heterosis_summary)) {
        write.csv(values$diallel_all_traits_results$summary_tables$heterosis_summary,
                  file.path(temp_dir, "heterosis_summary.csv"),
                  row.names = FALSE)
      }
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "diallel_summary.zip")
      files_to_zip <- list.files(temp_dir, pattern = "\\.csv$", full.names = TRUE)
      zip(zip_file, files = files_to_zip, flags = "-j")
      
      # Copy to download location
      file.copy(zip_file, file)
    }
  )
  
  # ==========================================================================
  # METAN ANALYSIS OBSERVERS
  # ==========================================================================
  
  # Prepare METAN data with withProgress
  observeEvent(input$metan_setup_data, {
    req(values$metan_raw_data, input$metan_gen_var, 
        input$metan_rep_var, input$metan_resp_vars)
    
    withProgress(message = 'Preparing METAN data...', value = 0, {
      
      tryCatch({
        incProgress(0.2, detail = "Checking YEAR/SEASON columns...")
        
        has_year_season <- FALSE
        if (all(c("YEAR", "SEASON") %in% colnames(values$metan_raw_data))) {
          has_year_season <- TRUE
          cat("Data has YEAR and SEASON columns. ENV will be created automatically.\n")
        } else {
          colnames_lower <- tolower(colnames(values$metan_raw_data))
          if (all(c("year", "season") %in% colnames_lower)) {
            has_year_season <- TRUE
            cat("Data has year and season columns (case-insensitive). ENV will be created automatically.\n")
          }
        }
        
        incProgress(0.4, detail = "Determining environment variable...")
        
        env_var_to_use <- NULL
        if (!has_year_season) {
          if (is.null(input$metan_env_var) || input$metan_env_var == "") {
            env_candidates <- c("ENV", "Environment", "environment", "LOCATION", "Location", "Site")
            for (candidate in env_candidates) {
              if (candidate %in% colnames(values$metan_raw_data)) {
                env_var_to_use <- candidate
                break
              }
            }
            
            if (is.null(env_var_to_use)) {
              stop("No environment variable found. Please ensure your data has:\n",
                   "1. YEAR and SEASON columns (will be combined to create ENV), OR\n",
                   "2. An environment column (ENV, Environment, Location, etc.)")
            } else {
              cat("Using automatically detected environment variable:", env_var_to_use, "\n")
            }
          } else {
            env_var_to_use <- input$metan_env_var
          }
        }
        
        incProgress(0.6, detail = "Preparing data structure...")
        
        metan_data <- prepare_metan_data(
          raw_data = values$metan_raw_data,
          env_var = env_var_to_use,
          gen_var = input$metan_gen_var,
          rep_var = input$metan_rep_var,
          resp_vars = input$metan_resp_vars
        )
        
        values$metan_processed_data <- metan_data
        
        incProgress(0.8, detail = "Creating summary...")
        
        env_summary <- metan_data %>%
          group_by(ENV) %>%
          summarise(
            n_genotypes = length(unique(GEN)),
            n_reps = length(unique(REP)),
            n_obs = n(),
            .groups = 'drop'
          )
        
        output$metan_data_status <- renderPrint({
          cat("=== METAN DATA STATUS ===\n\n")
          cat("Data preparation method: ")
          if (has_year_season) {
            cat("ENV created from YEAR and SEASON columns\n")
          } else {
            cat("ENV from selected environment variable\n")
          }
          cat("\n")
          cat("Environments:", length(unique(metan_data$ENV)), "\n")
          cat("Genotypes:", length(unique(metan_data$GEN)), "\n")
          cat("Replicates:", length(unique(metan_data$REP)), "\n")
          cat("Traits:", paste(input$metan_resp_vars, collapse = ", "), "\n")
          cat("Total observations:", nrow(metan_data), "\n")
          
          if (nrow(env_summary) <= 10) {
            cat("\nEnvironment summary:\n")
            print(env_summary, row.names = FALSE)
          } else {
            cat("\nFirst 10 environments:\n")
            print(head(env_summary, 10), row.names = FALSE)
            cat("... and", nrow(env_summary) - 10, "more environments\n")
          }
        })
        
        incProgress(1.0, detail = "Complete!")
        showNotification("METAN data prepared successfully!", type = "message", duration = 3)
        
      }, error = function(e) {
        output$metan_data_status <- renderPrint({
          cat("=== METAN DATA PREPARATION ERROR ===\n\n")
          cat("Error:", e$message, "\n\n")
          cat("To use METAN analysis, ensure your data has:\n")
          cat("1. YEAR and SEASON columns (will be combined to create ENV), OR\n")
          cat("2. An environment column (ENV, Environment, Location, etc.)\n")
          cat("3. A genotype column (Genotype, GEN, etc.)\n")
          cat("4. A replicate column (REP, Rep, etc.)\n")
          cat("5. At least one numeric trait column\n")
        })
        showNotification(paste("Error preparing METAN data:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  # Run METAN analysis with withProgress
  observeEvent(input$run_metan_analysis, {
    req(values$metan_processed_data, input$metan_resp_vars)
    
    withProgress(message = 'Running METAN analyses...', value = 0, {
      
      tryCatch({
        incProgress(0.1, detail = "Starting analyses...")
        
        values$metan_results <- run_metan_analysis(
          values$metan_processed_data,
          input$metan_resp_vars
        )
        
        incProgress(1.0, detail = "Complete!")
        showNotification("METAN analysis completed successfully!", type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error in METAN analysis:", e$message), type = "error", duration = 5)
      })
    })
  })
  
  
  # ==========================================================================
  # MULTI-TRAIT ANALYSIS SECTION
  # ==========================================================================
  
  # Multi-trait trait selector
  output$multi_trait_selector <- renderUI({
    req(values$pheno_data)
    
    tryCatch({
      traits <- setdiff(colnames(values$pheno_data), "Genotype")
      numeric_traits <- character()
      
      for (trait in traits) {
        if (is.numeric(values$pheno_data[[trait]])) {
          numeric_traits <- c(numeric_traits, trait)
        }
      }
      
      if (length(numeric_traits) == 0) {
        return(div(class = "alert alert-warning",
                   "No numeric traits found in phenotype data. 
                  Please check your data format."))
      }
      
      selectizeInput("multi_traits", "Select Traits for Analysis",
                     choices = numeric_traits,
                     multiple = TRUE,
                     selected = numeric_traits[1:min(2, length(numeric_traits))],
                     options = list(maxItems = 10),
                     width = "100%")
      
    }, error = function(e) {
      return(div(class = "alert alert-danger",
                 paste("Error creating trait selector:", e$message)))
    })
  })
  
  # ==========================================================================
  # REPORT GENERATION
  # ==========================================================================
  # Generate comprehensive report
  observeEvent(input$generate_report, {
    req(values$gwas_results)
    
    showModal(modalDialog(
      title = "Generating Report",
      "Creating comprehensive analysis report...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      # Capture input parameters
      params <- capture_input_parameters(input, values)
      
      # Capture output results
      results <- capture_output_results(values, input)
      
      # Generate report file name
      report_file <- paste0("GWAS_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".Rmd")
      
      # Generate comprehensive report
      generate_comprehensive_report(input, params, results, report_file)
      
      # Render the report
      if (input$report_format == "HTML") {
        output_file <- rmarkdown::render(report_file, output_format = "html_document")
      } else if (input$report_format == "PDF") {
        output_file <- rmarkdown::render(report_file, output_format = "pdf_document")
      } else {
        output_file <- rmarkdown::render(report_file, output_format = "word_document")
      }
      
      # Update report preview
      output$report_preview <- renderUI({
        if (input$report_format == "HTML") {
          includeHTML(output_file)
        } else {
          tags$div(
            h4("Report Generated Successfully"),
            p("Report file:", basename(output_file)),
            p("File size:", format(file.size(output_file), big.mark = ","), "bytes")
          )
        }
      })
      
      # Download handler for report
      output$download_report <- downloadHandler(
        filename = function() {
          basename(output_file)
        },
        content = function(file) {
          file.copy(output_file, file)
        }
      )
      
      removeModal()
      showNotification("Report generated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Report generation error:", e$message), type = "error")
    })
  })
  
  # Report preview output
  output$report_preview <- renderUI({
    tagList(
      h4("Report Preview"),
      p("Configure report settings and click 'Generate Report' to create a comprehensive analysis report."),
      p("The report will include:"),
      tags$ul(
        tags$li("Executive summary"),
        tags$li("Input parameters and methods"),
        tags$li("Data quality metrics"),
        tags$li("GWAS results and visualizations"),
        tags$li("Multi-trait analysis findings"),
        tags$li("Diallel analysis results (if available)"),
        tags$li("MET analysis results (if available)"),
        tags$li("Recommendations and next steps")
      )
    )
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS - DATA UPLOAD TAB
  # ==========================================================================
  
  # Data summary output
  output$data_summary <- renderPrint({
    cat("=== DATA UPLOAD SUMMARY ===\n\n")
    
    if (!is.null(values$pheno_data)) {
      cat("📊 PHENOTYPIC DATA:\n")
      cat("  Samples:", nrow(values$pheno_data), "\n")
      cat("  Traits:", ncol(values$pheno_data) - 1, "\n")
      trait_names <- setdiff(colnames(values$pheno_data), "Genotype")
      if (length(trait_names) <= 10) {
        cat("  Trait names:", paste(trait_names, collapse = ", "), "\n")
      } else {
        cat("  Trait names (first 10):", paste(head(trait_names, 10), collapse = ", "), "\n")
      }
      cat("\n")
    } else {
      cat("📊 PHENOTYPIC DATA: Not loaded\n\n")
    }
    
    if (!is.null(values$gl_object)) {
      cat("🧬 GENOTYPIC DATA:\n")
      cat("  Samples:", nInd(values$gl_object), "\n")
      cat("  Markers:", nLoc(values$gl_object), "\n")
      cat("  Ploidy:", values$gl_object@ploidy, "\n\n")
    } else {
      cat("🧬 GENOTYPIC DATA: Not loaded\n\n")
    }
    
    if (!is.null(values$metadata)) {
      cat("⚙️ METADATA:\n")
      cat("  Entries:", nrow(values$metadata), "\n")
      cat("  Variables:", ncol(values$metadata), "\n\n")
    } else {
      cat("⚙️ METADATA: Not loaded\n\n")
    }
    
    if (!is.null(values$diallel_data)) {
      cat("🌱 DIALLEL DATA:\n")
      cat("  Entries:", nrow(values$diallel_data), "\n")
      cat("  Variables:", ncol(values$diallel_data), "\n\n")
    } else {
      cat("🌱 DIALLEL DATA: Not loaded\n\n")
    }
    
    if (!is.null(values$metan_raw_data)) {
      cat("🌍 METAN DATA:\n")
      cat("  Observations:", nrow(values$metan_raw_data), "\n")
      cat("  Variables:", ncol(values$metan_raw_data), "\n")
    } else {
      cat("🌍 METAN DATA: Not loaded\n\n")
    }
  })
  
  # Phenotype data preview
  output$pheno_preview <- renderDT({
    req(values$pheno_data)
    
    datatable(
      values$pheno_data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf')
      ),
      class = 'display compact',
      rownames = FALSE,
      extensions = 'Buttons'
    )
  })
  
  # Genotype summary
  output$geno_summary <- renderPrint({
    req(values$gl_object)
    
    cat("=== GENOTYPIC DATA SUMMARY ===\n\n")
    cat("Number of individuals:", nInd(values$gl_object), "\n")
    cat("Number of loci:", nLoc(values$gl_object), "\n")
    cat("Ploidy:", values$gl_object@ploidy, "\n")
    
    # Calculate basic statistics
    cat("\n=== BASIC STATISTICS ===\n")
    
    # Call rate
    callrate <- 1 - (sum(is.na(as.matrix(values$gl_object))) / 
                       (nInd(values$gl_object) * nLoc(values$gl_object)))
    cat("Call rate:", round(callrate * 100, 2), "%\n")
    
    # Heterozygosity
    if (nLoc(values$gl_object) > 0) {
      het <- mean(colMeans(as.matrix(values$gl_object), na.rm = TRUE), na.rm = TRUE)
      cat("Average heterozygosity:", round(het, 4), "\n")
    }
    
    # Missing data summary
    missing_ind <- rowMeans(is.na(as.matrix(values$gl_object)))
    missing_loc <- colMeans(is.na(as.matrix(values$gl_object)))
    
    cat("\nMissing data by individual:\n")
    cat("  Min:", round(min(missing_ind) * 100, 2), "%\n")
    cat("  Max:", round(max(missing_ind) * 100, 2), "%\n")
    cat("  Mean:", round(mean(missing_ind) * 100, 2), "%\n")
    
    cat("\nMissing data by locus:\n")
    cat("  Min:", round(min(missing_loc) * 100, 2), "%\n")
    cat("  Max:", round(max(missing_loc) * 100, 2), "%\n")
    cat("  Mean:", round(mean(missing_loc) * 100, 2), "%\n")
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS - DATA PROCESSING TAB
  # ==========================================================================
  
  # Processed data summary
  output$processed_summary <- renderPrint({
    req(values$filtered_gl)
    
    cat("=== PROCESSED DATA SUMMARY ===\n\n")
    cat("After QC filtering:\n")
    cat("  Individuals:", nInd(values$filtered_gl), "\n")
    cat("  Loci:", nLoc(values$filtered_gl), "\n")
    
    if (!is.null(values$matched_data)) {
      cat("\nAfter genotype-phenotype matching:\n")
      cat("  Matched samples:", nrow(values$matched_data$pheno), "\n")
      cat("  Markers:", ncol(values$matched_data$geno), "\n")
    }
  })
  
  # QC plots
  output$missing_data_plot <- renderPlot({
    if (!is.null(values$qc_plots$sample_missing)) {
      values$qc_plots$sample_missing
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Run Quality Control first", size = 6) +
        theme_void()
    }
  })
  
  output$maf_plot <- renderPlot({
    if (!is.null(values$qc_plots$maf_dist)) {
      values$qc_plots$maf_dist
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Run Quality Control first", size = 6) +
        theme_void()
    }
  })
  
  output$het_plot <- renderPlot({
    if (!is.null(values$qc_plots$het_dist)) {
      values$qc_plots$het_dist
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Run Quality Control first", size = 6) +
        theme_void()
    }
  })
  
  # QC table
  output$qc_table <- renderDT({
    req(values$filtered_gl)
    
    # Calculate QC statistics
    geno_matrix <- as.matrix(values$filtered_gl)
    
    qc_stats <- data.frame(
      Metric = c("Total Samples", "Total Markers", 
                 "Mean Call Rate", "Mean MAF",
                 "Mean Heterozygosity", "Missing Data %"),
      Value = c(
        nInd(values$filtered_gl),
        nLoc(values$filtered_gl),
        round(1 - mean(is.na(geno_matrix)), 3),
        round(mean(colMeans(geno_matrix, na.rm = TRUE)/2), 4),
        round(mean(rowMeans(geno_matrix == 1, na.rm = TRUE)), 4),
        round(mean(is.na(geno_matrix)) * 100, 2)
      )
    )
    
    datatable(
      qc_stats,
      options = list(
        pageLength = 10,
        dom = 't'
      ),
      rownames = FALSE
    )
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS - POPULATION STRUCTURE TAB
  # ==========================================================================
  
  # PCA scree plot
  output$pca_scree_plot <- renderPlot({
    req(values$pca_results)
    
    variance_df <- head(values$pca_results$variance, 10)
    
    ggplot(variance_df, aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_line(color = "darkred", size = 1) +
      geom_point(color = "darkred", size = 2) +
      labs(
        title = "PCA Scree Plot",
        x = "Principal Component",
        y = "Variance Explained (%)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12)
      )
  })
  
  # PCA scatter plot
  output$pca_scatter_plot <- renderPlot({
    req(values$pca_results, input$pca_x, input$pca_y)
    
    pca_df <- values$pca_results$scores
    
    ggplot(pca_df, aes_string(x = paste0("PC", input$pca_x), 
                              y = paste0("PC", input$pca_y))) +
      geom_point(size = 3, alpha = 0.7, color = "steelblue") +
      theme_minimal() +
      labs(
        title = paste("PC", input$pca_x, "vs PC", input$pca_y),
        x = paste0("PC", input$pca_x, " (", 
                   round(values$pca_results$variance$Variance[as.numeric(input$pca_x)], 1), "%)"),
        y = paste0("PC", input$pca_y, " (", 
                   round(values$pca_results$variance$Variance[as.numeric(input$pca_y)], 1), "%)")
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 11)
      )
  })
  
  # PCA variance plot
  output$pca_variance_plot <- renderPlot({
    req(values$pca_results)
    
    variance_df <- head(values$pca_results$variance, 10)
    
    ggplot(variance_df, aes(x = PC, y = Cumulative)) +
      geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
      labs(
        title = "Cumulative Variance Explained",
        x = "Number of Principal Components",
        y = "Cumulative Variance (%)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12)
      )
  })
  
  # Structure status
  output$structure_status <- renderPrint({
    cat("=== POPULATION STRUCTURE ANALYSIS ===\n\n")
    
    if (!is.null(values$pca_results)) {
      cat("✓ PCA completed successfully\n")
      cat("  PCs calculated:", nrow(values$pca_results$variance), "\n")
      cat("  Variance explained by PC1:", 
          round(values$pca_results$variance$Variance[1], 1), "%\n")
      cat("  Total variance (first 5 PCs):", 
          round(sum(values$pca_results$variance$Variance[1:5]), 1), "%\n")
    } else {
      cat("✗ PCA not yet performed\n")
    }
    
    cat("\n")
    
    if (!is.null(values$kinship_matrix)) {
      cat("✓ Kinship matrix calculated\n")
      cat("  Dimensions:", dim(values$kinship_matrix), "\n")
      cat("  Method:", input$kinship_method, "\n")
    } else {
      cat("✗ Kinship matrix not yet calculated\n")
    }
  })
  
  # Kinship heatmap
  output$kinship_heatmap <- renderPlot({
    req(values$kinship_matrix)
    
    # Convert to data frame for plotting
    kinship_df <- as.data.frame(as.table(values$kinship_matrix))
    colnames(kinship_df) <- c("Var1", "Var2", "value")
    
    ggplot(kinship_df, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = median(kinship_df$value, na.rm = TRUE)) +
      labs(
        title = "Kinship Matrix Heatmap",
        x = "Samples",
        y = "Samples",
        fill = "Kinship"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)
      )
  })
  
  # Kinship summary
  output$kinship_summary <- renderPrint({
    req(values$kinship_matrix)
    
    cat("=== KINSHIP MATRIX SUMMARY ===\n\n")
    cat("Dimensions:", dim(values$kinship_matrix), "\n")
    cat("Range:", round(range(values$kinship_matrix, na.rm = TRUE), 3), "\n")
    cat("Mean:", round(mean(values$kinship_matrix, na.rm = TRUE), 3), "\n")
    cat("Standard deviation:", round(sd(values$kinship_matrix, na.rm = TRUE), 3), "\n")
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS - GWAS ANALYSIS TAB
  # ==========================================================================
  
  # GWAS results summary
  output$gwas_results_summary <- renderPrint({
    req(values$gwas_results)
    
    cat("=== GWAS RESULTS SUMMARY ===\n\n")
    cat("Trait analyzed:", input$gwas_trait, "\n")
    cat("Total SNPs tested:", nrow(values$gwas_results), "\n")
    cat("Significant SNPs (p < 0.05):", 
        sum(values$gwas_results$P_value < 0.05, na.rm = TRUE), "\n")
    cat("Genome-wide significant (p < 5e-8):", 
        sum(values$gwas_results$P_value < 5e-8, na.rm = TRUE), "\n")
    
    if (nrow(values$gwas_results) > 0) {
      cat("\nTop 5 most significant SNPs:\n")
      top_snps <- head(values$gwas_results[order(values$gwas_results$P_value), ], 5)
      for (i in 1:nrow(top_snps)) {
        cat(paste0(i, ". ", top_snps$SNP_ID[i], 
                   " (Chr", top_snps$Chromosome[i], ":", top_snps$Position[i], 
                   ") - p = ", format(top_snps$P_value[i], scientific = TRUE, digits = 3), 
                   "\n"))
      }
    }
  })
  
  # GWAS results table
  output$gwas_results_table <- renderDT({
    req(values$gwas_results)
    
    # Format the results table
    display_df <- values$gwas_results %>%
      mutate(
        P_value = format(P_value, scientific = TRUE, digits = 3),
        P_adjusted = format(P_adjusted, scientific = TRUE, digits = 3),
        Effect = round(Effect, 4),
        SE = round(SE, 4),
        R_squared = round(R_squared, 4)
      ) %>%
      select(SNP_ID, Chromosome, Position, P_value, P_adjusted, 
             Effect, SE, R_squared, Significant_FDR)
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      extensions = 'Buttons',
      filter = 'top'
    )
  })
  
  # Significant SNPs table
  output$significant_snps_table <- renderDT({
    req(values$gwas_results)
    
    # Filter significant SNPs
    sig_snps <- values$gwas_results %>%
      filter(Significant_FDR == TRUE) %>%
      mutate(
        P_value = format(P_value, scientific = TRUE, digits = 3),
        P_adjusted = format(P_adjusted, scientific = TRUE, digits = 3),
        Effect = round(Effect, 4),
        SE = round(SE, 4)
      ) %>%
      select(SNP_ID, Chromosome, Position, P_value, P_adjusted, 
             Effect, SE, R_squared)
    
    datatable(
      sig_snps,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE
    )
  })
  
  # Manhattan plot
  output$manhattan_plot <- renderPlot({
    req(values$gwas_results)
    
    # Prepare data for Manhattan plot
    manhattan_data <- values$gwas_results %>%
      mutate(
        log10p = -log10(P_value),
        Chromosome = as.factor(Chromosome)
      )
    
    # Check if we have valid chromosome data
    if (length(unique(manhattan_data$Chromosome)) > 1) {
      # Create Manhattan plot with chromosome labels on x-axis
      ggplot(manhattan_data, aes(x = Chromosome, y = log10p, color = Chromosome)) +
        geom_point(alpha = 0.6, size = 1.5, position = position_jitter(width = 0.2)) +
        geom_hline(yintercept = -log10(0.05/nrow(manhattan_data)), 
                   linetype = "dashed", color = "red", alpha = 0.7) +
        geom_hline(yintercept = -log10(5e-8), 
                   linetype = "dashed", color = "blue", alpha = 0.7) +
        labs(
          title = paste("Manhattan Plot -", input$gwas_trait),
          x = "Chromosome",
          y = "-log10(P-value)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "none",
          panel.grid.major.x = element_blank()
        ) +
        scale_color_brewer(palette = "Set3")
    } else {
      # Fallback if only one chromosome
      ggplot(manhattan_data, aes(x = Position, y = log10p)) +
        geom_point(alpha = 0.6, size = 1.5, color = "steelblue") +
        geom_hline(yintercept = -log10(0.05/nrow(manhattan_data)), 
                   linetype = "dashed", color = "red", alpha = 0.7) +
        geom_hline(yintercept = -log10(5e-8), 
                   linetype = "dashed", color = "blue", alpha = 0.7) +
        labs(
          title = paste("Manhattan Plot -", input$gwas_trait),
          x = "Genomic Position",
          y = "-log10(P-value)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
        )
    }
  })
  
  
  
  # GWAS Analysis outputs
  # output$manhattan_plot <- renderPlotly({
  #   gwas_data <- rv$gwas_data
  #   
  #   # Prepare data for Manhattan plot
  #   plot_data <- gwas_data %>%
  #     group_by(Chromosome) %>%
  #     mutate(Chromosome_Median = median(Position)) %>%
  #     ungroup() %>%
  #     arrange(Chromosome, Position) %>%
  #     mutate(Cumulative_Position = Position + lag(cumsum(as.numeric(factor(Chromosome)))*1e8, default = 0))
  #   
  #   # Color by significance
  #   plot_data$Significant <- plot_data$P_Value < 10^(-input$manhattan_pval)
  #   
  #   # Create the plot
  #   p <- plot_ly(plot_data, 
  #                x = ~Cumulative_Position, 
  #                y = ~neg_log10_p,
  #                type = 'scatter',
  #                mode = 'markers',
  #                color = ~Significant,
  #                colors = c('gray', 'red'),
  #                marker = list(size = 10, opacity = 0.7),
  #                hoverinfo = 'text',
  #                text = ~paste('SNP:', SNP_ID,
  #                              '<br>Chr:', Chromosome,
  #                              '<br>Position:', format(Position, big.mark = ","),
  #                              '<br>P-value:', format(P_Value, scientific = TRUE, digits = 3),
  #                              '<br>Gene:', Gene,
  #                              '<br>Effect:', round(Effect_Size, 3),
  #                              '<br>MAF:', round(MAF, 3))) %>%
  #     layout(title = paste("GWAS Manhattan Plot (-log10(P) >", input$manhattan_pval, ")"),
  #            xaxis = list(title = "Chromosomal Position",
  #                         showgrid = FALSE,
  #                         zeroline = FALSE),
  #            yaxis = list(title = "-log10(P-value)"),
  #            showlegend = FALSE)
  #   
  #   # Add significance line
  #   p <- p %>% add_segments(x = min(plot_data$Cumulative_Position),
  #                           xend = max(plot_data$Cumulative_Position),
  #                           y = input$manhattan_pval,
  #                           yend = input$manhattan_pval,
  #                           line = list(color = 'red', dash = 'dash'),
  #                           showlegend = FALSE)
  #   
  #   p
  # })
  
  
  # QQ plot
  output$qq_plot <- renderPlot({
    req(values$gwas_results)
    
    # Prepare data for QQ plot
    observed <- -log10(sort(values$gwas_results$P_value))
    expected <- -log10(ppoints(length(observed)))
    
    qq_data <- data.frame(Expected = expected, Observed = observed)
    
    # Create QQ plot
    ggplot(qq_data, aes(x = Expected, y = Observed)) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "QQ Plot",
        x = "Expected -log10(P)",
        y = "Observed -log10(P)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
  })
  
  # Effect size plot
  output$effect_plot <- renderPlot({
    req(values$gwas_results)
    
    # Create effect size distribution plot
    ggplot(values$gwas_results, aes(x = Effect)) +
      geom_histogram(fill = "steelblue", alpha = 0.7, bins = 50) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = "Effect Size Distribution",
        x = "Effect Size",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS FOR MULTI-TRAIT TAB
  # ==========================================================================
  
  # Multi-trait GWAS table
  output$multi_trait_gwas_table <- renderDT({
    req(values$multi_trait_results$combined_results)
    
    # Format the table
    display_df <- values$multi_trait_results$combined_results %>%
      mutate(
        Min_P_value = format(Min_P_value, scientific = TRUE, digits = 3),
        Max_Effect = round(Max_Effect, 4),
        P_adjusted = format(P_adjusted, scientific = TRUE, digits = 3),
        neg_log10_p = round(neg_log10_p, 2)
      ) %>%
      select(SNP_ID, Chromosome, Position, N_Traits, Traits, 
             Min_P_value, P_adjusted, neg_log10_p, Max_Effect)
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Multi-Trait GWAS Results",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'N_Traits',
        background = styleColorBar(range(display_df$N_Traits, na.rm = TRUE), 'lightblue'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      ) %>%
      formatStyle(
        'neg_log10_p',
        background = styleColorBar(range(display_df$neg_log10_p, na.rm = TRUE), 'lightgreen'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  
  # Multi-trait GWAS table - FIXED
  output$multi_trait_gwas_table <- renderDT({
    req(values$multi_trait_results)
    
    # Check if we have combined results
    if (is.null(values$multi_trait_results$combined_results) || 
        nrow(values$multi_trait_results$combined_results) == 0) {
      return(datatable(
        data.frame(Message = "No multi-trait GWAS results available. Run analysis first."),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Create display dataframe
    display_df <- values$multi_trait_results$combined_results
    
    # Format columns if they exist
    if ("Min_P_value" %in% colnames(display_df)) {
      display_df$Min_P_value <- format(display_df$Min_P_value, scientific = TRUE, digits = 3)
    }
    
    if ("P_adjusted" %in% colnames(display_df)) {
      display_df$P_adjusted <- format(display_df$P_adjusted, scientific = TRUE, digits = 3)
    }
    
    if ("Max_Effect" %in% colnames(display_df)) {
      display_df$Max_Effect <- round(display_df$Max_Effect, 4)
    }
    
    if ("neg_log10_p" %in% colnames(display_df)) {
      display_df$neg_log10_p <- round(display_df$neg_log10_p, 2)
    }
    
    # Select columns to display (only those that exist)
    available_cols <- colnames(display_df)
    display_cols <- c("SNP_ID", "Chromosome", "Position", "Min_P_value")
    
    # Add other columns if they exist
    optional_cols <- c("P_adjusted", "neg_log10_p", "Max_Effect", "N_Traits", "Significant_FDR")
    for (col in optional_cols) {
      if (col %in% available_cols) {
        display_cols <- c(display_cols, col)
      }
    }
    
    # Ensure we only use columns that exist
    display_cols <- display_cols[display_cols %in% available_cols]
    
    if (length(display_cols) == 0) {
      return(datatable(
        data.frame(Message = "No displayable columns found in results"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    display_df <- display_df[, display_cols, drop = FALSE]
    
    # Create datatable
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Multi-Trait GWAS Results (Combined P-values)",
      extensions = 'Buttons',
      filter = 'top'
    )
  })
  
  # Genetic correlation plot - FIXED
  output$genetic_cor_plot <- renderPlot({
    req(values$multi_trait_results)
    
    # Check if we have correlation data
    if (is.null(values$multi_trait_results$trait_correlations)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No correlation data available\nRun multi-trait analysis first", 
                        size = 6) +
               theme_void())
    }
    
    # Extract correlation matrix
    cor_matrix <- NULL
    if (is.matrix(values$multi_trait_results$trait_correlations)) {
      cor_matrix <- values$multi_trait_results$trait_correlations
    } else if (is.list(values$multi_trait_results$trait_correlations)) {
      if (!is.null(values$multi_trait_results$trait_correlations$correlations)) {
        cor_matrix <- values$multi_trait_results$trait_correlations$correlations
      }
    }
    
    if (is.null(cor_matrix) || nrow(cor_matrix) == 0 || ncol(cor_matrix) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No valid correlation matrix available", 
                        size = 6) +
               theme_void())
    }
    
    # Create correlation heatmap with ggplot2
    tryCatch({
      melted_cor <- melt(cor_matrix)
      
      ggplot(melted_cor, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                             midpoint = 0, limit = c(-1, 1), space = "Lab",
                             name = "Correlation") +
        geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          panel.grid = element_blank()
        ) +
        labs(
          title = "Trait Correlation Heatmap",
          x = "",
          y = ""
        ) +
        coord_fixed()
    }, error = function(e) {
      cat("Error in genetic correlation plot:", e$message, "\n")
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("Error creating correlation plot:", e$message), 
                        size = 6) +
               theme_void())
    })
  })
  
  # Pleiotropy table - FIXED
  output$pleiotropy_table <- renderDT({
    req(values$multi_trait_results)
    
    if (is.null(values$multi_trait_results$combined_results) || 
        nrow(values$multi_trait_results$combined_results) == 0) {
      return(datatable(
        data.frame(Message = "No combined results available. Run multi-trait analysis first."),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Check if N_Traits column exists
    if (!"N_Traits" %in% colnames(values$multi_trait_results$combined_results)) {
      return(datatable(
        data.frame(Message = "No pleiotropy data (N_Traits column not found)"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Filter for SNPs affecting multiple traits
    pleiotropic_df <- values$multi_trait_results$combined_results
    
    # Only show SNPs affecting >1 trait
    if (any(pleiotropic_df$N_Traits > 1, na.rm = TRUE)) {
      pleiotropic_df <- pleiotropic_df[pleiotropic_df$N_Traits > 1, ]
    } else {
      return(datatable(
        data.frame(Message = "No SNPs found affecting more than 1 trait"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    if (nrow(pleiotropic_df) == 0) {
      return(datatable(
        data.frame(Message = "No pleiotropic SNPs found (affecting >1 trait)"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Format the data
    display_df <- pleiotropic_df
    
    # Format columns
    if ("Min_P_value" %in% colnames(display_df)) {
      display_df$Min_P_value <- format(display_df$Min_P_value, scientific = TRUE, digits = 3)
    }
    
    if ("P_adjusted" %in% colnames(display_df)) {
      display_df$P_adjusted <- format(display_df$P_adjusted, scientific = TRUE, digits = 3)
    }
    
    if ("Max_Effect" %in% colnames(display_df)) {
      display_df$Max_Effect <- round(display_df$Max_Effect, 4)
    }
    
    # Select columns to show
    display_cols <- c("SNP_ID", "Chromosome", "Position", "Min_P_value", "N_Traits")
    
    # Add optional columns if they exist
    optional_cols <- c("P_adjusted", "Max_Effect", "Significant_FDR")
    for (col in optional_cols) {
      if (col %in% colnames(display_df)) {
        display_cols <- c(display_cols, col)
      }
    }
    
    display_df <- display_df[, display_cols, drop = FALSE]
    
    # Order by N_Traits (descending) then by P-value (ascending)
    display_df <- display_df[order(-display_df$N_Traits, display_df$Min_P_value), ]
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Pleiotropic SNPs (Affecting Multiple Traits)",
      extensions = 'Buttons',
      filter = 'top'
    )
  })
  
  # Pleiotropy visualization plot - FIXED
  output$pleiotropy_plot_output <- renderPlot({
    req(values$multi_trait_results)
    
    if (is.null(values$multi_trait_results$combined_results) || 
        nrow(values$multi_trait_results$combined_results) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No data available\nRun multi-trait analysis first", 
                        size = 6) +
               theme_void())
    }
    
    # Check if N_Traits column exists
    if (!"N_Traits" %in% colnames(values$multi_trait_results$combined_results)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No pleiotropy data (N_Traits column not found)", 
                        size = 6) +
               theme_void())
    }
    
    # Create a plot of N_Traits distribution
    plot_data <- values$multi_trait_results$combined_results
    
    # Count SNPs by number of traits
    trait_counts <- table(plot_data$N_Traits)
    
    # Create data frame for plotting
    plot_df <- data.frame(
      N_Traits = as.numeric(names(trait_counts)),
      Count = as.numeric(trait_counts)
    )
    
    # Create the plot
    ggplot(plot_df, aes(x = factor(N_Traits), y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_text(aes(label = Count), vjust = -0.5, size = 4) +
      labs(
        title = "Distribution of SNPs by Number of Traits Affected",
        subtitle = paste("Total SNPs:", sum(plot_df$Count)),
        x = "Number of Traits",
        y = "Number of SNPs"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  })
  
  # Multi-trait statistics - FIXED
  output$multi_trait_stats <- renderPrint({
    req(values$multi_trait_results)
    
    cat("=== MULTI-TRAIT ANALYSIS STATISTICS ===\n\n")
    
    # Basic information
    if (!is.null(values$multi_trait_results$n_traits)) {
      cat("Number of traits analyzed:", values$multi_trait_results$n_traits, "\n")
    }
    
    if (!is.null(values$multi_trait_results$traits_analyzed)) {
      cat("Traits:", paste(values$multi_trait_results$traits_analyzed, collapse = ", "), "\n")
    }
    
    if (!is.null(values$multi_trait_results$samples_used)) {
      cat("Common samples with complete data:", length(values$multi_trait_results$samples_used), "\n")
    }
    
    cat("\n")
    
    # Individual trait statistics
    if (!is.null(values$multi_trait_results$trait_gwas_results)) {
      cat("=== INDIVIDUAL TRAIT GWAS RESULTS ===\n\n")
      
      for (trait in names(values$multi_trait_results$trait_gwas_results)) {
        trait_result <- values$multi_trait_results$trait_gwas_results[[trait]]
        
        if (!is.null(trait_result) && is.data.frame(trait_result)) {
          cat("Trait:", trait, "\n")
          cat("  Total SNPs tested:", nrow(trait_result), "\n")
          
          # Count significant SNPs
          if ("P_adjusted" %in% colnames(trait_result)) {
            sig_fdr <- sum(trait_result$P_adjusted < 0.05, na.rm = TRUE)
            cat("  Significant SNPs (FDR < 0.05):", sig_fdr, "\n")
          }
          
          if ("P_value" %in% colnames(trait_result)) {
            min_p <- min(trait_result$P_value, na.rm = TRUE)
            cat("  Minimum P-value:", format(min_p, scientific = TRUE, digits = 3), "\n")
          }
          
          cat("\n")
        }
      }
    }
    
    # Combined results statistics
    if (!is.null(values$multi_trait_results$combined_results)) {
      cat("=== COMBINED RESULTS ===\n\n")
      
      combined <- values$multi_trait_results$combined_results
      
      cat("Total SNPs in combined analysis:", nrow(combined), "\n")
      
      if ("N_Traits" %in% colnames(combined)) {
        cat("SNPs affecting 1 trait:", sum(combined$N_Traits == 1, na.rm = TRUE), "\n")
        cat("SNPs affecting 2 traits:", sum(combined$N_Traits == 2, na.rm = TRUE), "\n")
        cat("SNPs affecting 3+ traits:", sum(combined$N_Traits >= 3, na.rm = TRUE), "\n")
        cat("Maximum traits affected by single SNP:", max(combined$N_Traits, na.rm = TRUE), "\n")
      }
      
      if ("Significant_FDR" %in% colnames(combined)) {
        sig_fdr <- sum(combined$Significant_FDR, na.rm = TRUE)
        cat("Significant SNPs in combined analysis (FDR < 0.05):", sig_fdr, "\n")
      }
      
      if ("Min_P_value" %in% colnames(combined)) {
        min_p <- min(combined$Min_P_value, na.rm = TRUE)
        cat("Minimum P-value across all traits:", format(min_p, scientific = TRUE, digits = 3), "\n")
      }
    }
  })
  
  
  # Download multi-trait results
  output$download_multi_trait <- downloadHandler(
    filename = function() {
      paste0("multi_trait_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(values$multi_trait_results)
      
      # Create temporary directory
      temp_dir <- tempfile("multi_trait_")
      dir.create(temp_dir)
      
      # Save all data frames
      if (!is.null(values$multi_trait_results$combined_results)) {
        write.csv(values$multi_trait_results$combined_results,
                  file.path(temp_dir, "combined_gwas_results.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$multi_trait_results$trait_correlations)) {
        write.csv(values$multi_trait_results$trait_correlations$correlations,
                  file.path(temp_dir, "trait_correlations.csv"),
                  row.names = TRUE)
      }
      
      if (!is.null(values$multi_trait_results$genetic_correlations)) {
        write.csv(values$multi_trait_results$genetic_correlations,
                  file.path(temp_dir, "genetic_correlations.csv"),
                  row.names = TRUE)
      }
      
      # Save individual trait results
      if (!is.null(values$multi_trait_results$trait_results)) {
        for (trait in names(values$multi_trait_results$trait_results)) {
          if (!is.null(values$multi_trait_results$trait_results[[trait]])) {
            write.csv(values$multi_trait_results$trait_results[[trait]],
                      file.path(temp_dir, paste0("gwas_", trait, ".csv")),
                      row.names = FALSE)
          }
        }
      }
      
      # Save metadata
      meta_data <- data.frame(
        Parameter = c("Analysis Date", "Traits Analyzed", "Method", "Samples Used",
                      "Total SNPs", "Pleiotropic SNPs", "Timestamp"),
        Value = c(
          as.character(Sys.Date()),
          paste(values$multi_trait_results$traits_analyzed, collapse = ", "),
          values$multi_trait_results$method,
          length(values$multi_trait_results$samples_used),
          nrow(values$multi_trait_results$combined_results),
          sum(values$multi_trait_results$combined_results$N_Traits > 1, na.rm = TRUE),
          as.character(values$multi_trait_results$timestamp)
        )
      )
      write.csv(meta_data, file.path(temp_dir, "metadata.csv"), row.names = FALSE)
      
      # Create README
      readme_content <- paste(
        "MULTI-TRAIT ANALYSIS RESULTS",
        "============================",
        "",
        paste("Analysis Date:", Sys.time()),
        paste("Generated by: Plant Breeding Analysis Platform v4.0"),
        "",
        "Files included:",
        "1. combined_gwas_results.csv - Combined results across all traits",
        "2. trait_correlations.csv - Phenotypic correlations between traits",
        "3. genetic_correlations.csv - Genetic correlations between traits",
        "4. gwas_[trait].csv - Individual GWAS results for each trait",
        "5. metadata.csv - Analysis parameters and metadata",
        "",
        "Column descriptions for combined_gwas_results.csv:",
        "- SNP_ID: SNP identifier",
        "- Chromosome: Chromosome number",
        "- Position: Genomic position",
        "- N_Traits: Number of traits the SNP affects",
        "- Traits: Comma-separated list of affected traits",
        "- Min_P_value: Minimum p-value across traits",
        "- P_adjusted: False discovery rate adjusted p-value",
        "- neg_log10_p: -log10 of minimum p-value",
        "- Max_Effect: Maximum absolute effect size across traits",
        sep = "\n"
      )
      writeLines(readme_content, file.path(temp_dir, "README.txt"))
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "results.zip")
      zip::zip(zip_file, files = list.files(temp_dir, full.names = TRUE), 
               mode = "cherry-pick")
      
      # Copy to download location
      file.copy(zip_file, file)
      
      # Cleanup
      unlink(temp_dir, recursive = TRUE)
    }
  )
  
  # ==========================================================================
  # OUTPUT RENDERERS - DIALLEL ANALYSIS TAB
  # ==========================================================================
  
  # Add reactive value for multi-trait diallel results
  values$diallel_all_traits_results <- NULL
  
  # Observer for running all traits diallel analysis
  observeEvent(input$run_all_diallel_traits, {
    req(values$diallel_data)
    
    withProgress(message = 'Running diallel analysis for all traits...', value = 0, {
      
      tryCatch({
        incProgress(0.3, detail = "Finding numeric traits...")
        
        # Find numeric traits
        numeric_cols <- colnames(values$diallel_data)[sapply(values$diallel_data, is.numeric)]
        
        if (length(numeric_cols) == 0) {
          showNotification("No numeric traits found in diallel data", 
                           type = "error", duration = 5)
          return()
        }
        
        incProgress(0.6, detail = "Running analysis...")
        
        # Run multi-trait analysis
        results <- run_diallel_anova_all_traits(values$diallel_data, numeric_cols)
        
        values$diallel_all_traits_results <- results
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Multi-trait diallel analysis completed!", 
                         type = "message", duration = 3)
        
        # Update tab
        updateTabItems(session, "tabs", "diallel_all_traits")
        
      }, error = function(e) {
        showNotification(paste("Error in multi-trait diallel analysis:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  
  # Diallel ANOVA results
  output$diallel_anova_results <- renderPrint({
    req(values$diallel_results$anova)
    
    cat("=== DIALLEL ANOVA ===\n\n")
    
    if (!is.null(values$diallel_results$anova$table)) {
      print(values$diallel_results$anova$table)
    } else if (!is.null(values$diallel_results$anova$error)) {
      cat("Error:", values$diallel_results$anova$error, "\n")
    } else {
      cat("No ANOVA results available\n")
    }
  })
  
  # Diallel genetic parameters
  output$diallel_genetic_params <- renderPrint({
    req(values$diallel_results)
    
    cat("=== GENETIC PARAMETERS ===\n\n")
    
    if (!is.null(values$diallel_results$anova)) {
      cat("General Combining Ability (GCA):\n")
      cat("  Sum of Squares:", round(values$diallel_results$anova$ss_gca, 2), "\n")
      cat("  Mean Square:", round(values$diallel_results$anova$ms_gca, 2), "\n")
      cat("\n")
    }
    
    if (!is.null(values$diallel_results$overall_mean)) {
      cat("Overall mean:", round(values$diallel_results$overall_mean, 2), "\n")
    }
    
    if (!is.null(values$diallel_results$n_parents)) {
      cat("Number of parents:", values$diallel_results$n_parents, "\n")
    }
    
    if (!is.null(values$diallel_results$n_crosses)) {
      cat("Number of crosses:", values$diallel_results$n_crosses, "\n")
    }
  })
  
  # Diallel data table
  output$diallel_data_table <- renderDT({
    req(values$diallel_data)
    
    # Select relevant columns for display
    display_cols <- c("ACC.", "CrossType", "Parent1", "Parent2", "GYP", "GYR")
    display_cols <- display_cols[display_cols %in% colnames(values$diallel_data)]
    
    datatable(
      values$diallel_data[, display_cols],
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Diallel Data Overview"
    )
  })
  
  # Diallel data summary
  output$diallel_data_summary <- renderPrint({
    req(values$diallel_data)
    
    cat("=== DIALLEL DATA SUMMARY ===\n\n")
    cat("Total entries:", nrow(values$diallel_data), "\n")
    
    if (!is.null(values$diallel_parent_data)) {
      cat("Parents:", nrow(values$diallel_parent_data), "\n")
    }
    
    if (!is.null(values$diallel_cross_data)) {
      cat("Crosses:", nrow(values$diallel_cross_data), "\n")
    }
    
    # Trait summary
    trait_cols <- c("GYP", "GYR", "DTF", "DTM", "PL", "NPP", "NSP", "X100SW", "TW")
    trait_cols <- trait_cols[trait_cols %in% colnames(values$diallel_data)]
    
    cat("\nTrait Statistics:\n")
    for (trait in trait_cols) {
      if (is.numeric(values$diallel_data[[trait]])) {
        trait_mean <- mean(values$diallel_data[[trait]], na.rm = TRUE)
        trait_sd <- sd(values$diallel_data[[trait]], na.rm = TRUE)
        cat(paste0("  ", trait, ": ", round(trait_mean, 2), 
                   " ± ", round(trait_sd, 2), "\n"))
      }
    }
  })
  
  # Diallel trait descriptions
  output$diallel_trait_descriptions <- renderTable({
    trait_descriptions <- data.frame(
      Trait = c("GYP", "GYR", "DTF", "DTM", "PL", "NPP", "NSP", "X100SW", "TW"),
      Description = c("Grain Yield per Plant", "Grain Yield per Row", 
                      "Days to Flowering", "Days to Maturity", "Plant Length",
                      "Number of Pods per Plant", "Number of Seeds per Pod",
                      "100 Seed Weight", "Test Weight"),
      Unit = c("g", "g", "days", "days", "cm", "count", "count", "g", "g/L")
    )
    
    trait_descriptions
  })
  
  
  # ==========================================================================
  # OUTPUT RENDERERS - DIALLEL GCA/SCA TAB
  # ==========================================================================
  
  # GCA table
  output$diallel_gca_table <- renderDT({
    req(values$diallel_results$gca)
    
    datatable(
      values$diallel_results$gca,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(2, 'desc'))  # Sort by GCA value
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "General Combining Ability (GCA) Effects"
    ) %>%
      formatRound(columns = c('GCA', 'Mean_Performance'), digits = 3)
  })
  
  
  # SCA table with proper error handling
  output$diallel_sca_table <- renderDT({
    req(values$diallel_results$sca)
    
    tryCatch({
      # Convert SCA matrix to data frame
      sca_matrix <- values$diallel_results$sca
      
      # Check if it's a valid matrix
      if (!is.matrix(sca_matrix) || nrow(sca_matrix) == 0 || ncol(sca_matrix) == 0) {
        return(datatable(
          data.frame(Message = "SCA matrix is empty or invalid"),
          options = list(pageLength = 5),
          rownames = FALSE
        ))
      }
      
      # Convert to long format
      sca_df <- data.frame(
        Parent1 = character(),
        Parent2 = character(),
        SCA = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Extract non-NA, non-diagonal values
      for (i in 1:nrow(sca_matrix)) {
        for (j in 1:ncol(sca_matrix)) {
          if (i != j && !is.na(sca_matrix[i, j])) {
            sca_df <- rbind(sca_df, data.frame(
              Parent1 = rownames(sca_matrix)[i],
              Parent2 = colnames(sca_matrix)[j],
              SCA = round(sca_matrix[i, j], 4),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      # Remove duplicates (since matrix is symmetric)
      sca_df <- sca_df %>%
        mutate(key = paste(pmin(Parent1, Parent2), pmax(Parent1, Parent2))) %>%
        distinct(key, .keep_all = TRUE) %>%
        select(-key) %>%
        arrange(desc(abs(SCA)))
      
      datatable(
        sca_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Specific Combining Ability (SCA) Effects",
        extensions = 'Buttons'
      )
      
    }, error = function(e) {
      datatable(
        data.frame(Error = paste("Error processing SCA matrix:", e$message)),
        options = list(pageLength = 5),
        rownames = FALSE
      )
    })
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS - DIALLEL HETEROSIS TAB
  # ==========================================================================
  
  # Heterosis table
  output$diallel_heterosis_table <- renderDT({
    req(values$diallel_heterosis_data)
    
    datatable(
      values$diallel_heterosis_data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(5, 'desc'))  # Sort by MPH
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Heterosis Estimates"
    ) %>%
      formatRound(columns = c('MPH', 'BPH', 'Mid_Parent', 'Better_Parent'), digits = 2)
  })
  
  # Heterosis statistics
  output$diallel_heterosis_stats <- renderPrint({
    req(values$diallel_heterosis_data)
    
    cat("=== HETEROSIS STATISTICS ===\n\n")
    
    het_data <- values$diallel_heterosis_data
    
    cat("Mid-Parent Heterosis (MPH):\n")
    cat("  Mean:", round(mean(het_data$MPH, na.rm = TRUE), 2), "%\n")
    cat("  Range:", round(range(het_data$MPH, na.rm = TRUE), 2), "%\n")
    cat("  SD:", round(sd(het_data$MPH, na.rm = TRUE), 2), "%\n")
    cat("  Crosses with MPH > 15%:", sum(het_data$MPH > 15, na.rm = TRUE), "\n")
    
    cat("\nBetter-Parent Heterosis (BPH):\n")
    cat("  Mean:", round(mean(het_data$BPH, na.rm = TRUE), 2), "%\n")
    cat("  Range:", round(range(het_data$BPH, na.rm = TRUE), 2), "%\n")
    cat("  SD:", round(sd(het_data$BPH, na.rm = TRUE), 2), "%\n")
    cat("  Crosses with BPH > 10%:", sum(het_data$BPH > 10, na.rm = TRUE), "\n")
  })
  
  # Heterosis histogram
  output$diallel_heterosis_hist <- renderPlotly({
    req(values$diallel_heterosis_data, input$heterosis_type)
    
    het_data <- values$diallel_heterosis_data
    
    if("mph" %in% input$heterosis_type){
      p <- ggplot(het_data, aes(x = MPH)) +
        geom_histogram(fill = "steelblue", alpha = 0.7, bins = 20) +
        geom_vline(xintercept = mean(het_data$MPH, na.rm = TRUE), 
                   color = "red", linetype = "dashed") +
        labs(
          title = "Distribution of Mid-Parent Heterosis (MPH)",
          x = "MPH (%)",
          y = "Frequency"
        ) +
        theme_minimal()
      
      ggplotly(p)
      
    } else {
      p <- ggplot(het_data, aes(x = BPH)) +
        geom_histogram(fill = "steelblue", alpha = 0.7, bins = 20) +
        geom_vline(xintercept = mean(het_data$BPH, na.rm = TRUE), 
                   color = "red", linetype = "dashed") +
        labs(
          title = "Distribution of Mid-Parent Heterosis (BPH)",
          x = "BPH (%)",
          y = "Frequency"
        ) +
        theme_minimal()
      
      ggplotly(p)
      
    }
  })
  
  
  
  # Heterosis bar plot
  output$diallel_heterosis_bar <- renderPlotly({
    req(values$diallel_heterosis_data,input$heterosis_type)
    
    if("mph" %in% input$heterosis_type){
      # Get top crosses by MPH
      top_crosses <- values$diallel_heterosis_data %>%
        arrange(desc(MPH)) %>%
        head(10)
      
      p <- ggplot(top_crosses, aes(x = reorder(CrossID, MPH), y = MPH)) +
        geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
        coord_flip() +
        labs(
          title = "Top 10 Crosses by Heterosis",
          x = "Cross",
          y = "MPH (%)"
        ) +
        theme_minimal()
      
      ggplotly(p)
      
    } else {
      # Get top crosses by MPH
      top_crosses <- values$diallel_heterosis_data %>%
        arrange(desc(BPH)) %>%
        head(10)
      
      p <- ggplot(top_crosses, aes(x = reorder(CrossID, BPH), y = BPH)) +
        geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
        coord_flip() +
        labs(
          title = "Top 10 Crosses by Heterosis",
          x = "Cross",
          y = "BPH (%)"
        ) +
        theme_minimal()
      
      ggplotly(p)
      
    }
  })
  
  #===================================================================================================
  
  # Filter top 3 heterotic crosses by MPH and get their data from cross_snp_df
  observe({
    # Check if both dataframes exist
    if (!is.null(values$diallel_heterosis_data) && !is.null(values$cross_snp_df)) {
      
      # Get top 3 crosses by MPH from heterosis data
      top_3_crosses <- values$diallel_heterosis_data %>% 
        arrange(desc(MPH)) %>% 
        head(3) %>% 
        pull(CrossID)
      
      # Filter cross_snp_df to only include these top 3 crosses
      values$top_3_heterotic_crosses_snp <- values$cross_snp_df %>% 
        filter(Cross %in% top_3_crosses)
      
      # Create a summary of these top crosses
      if (nrow(values$top_3_heterotic_crosses_snp) > 0) {
        # Create summary statistics
        top_crosses_summary <- values$top_3_heterotic_crosses_snp %>% 
          group_by(Cross, Trait) %>% 
          summarise(
            N_SNPs = n_distinct(SNP_ID),
            Avg_P_value = mean(P_value, na.rm = TRUE),
            Min_P_value = min(P_value, na.rm = TRUE),
            Max_Effect = max(abs(Effect), na.rm = TRUE),
            N_Significant = sum(Significant, na.rm = TRUE),
            .groups = "drop"
          )
        
        values$top_crosses_summary <- top_crosses_summary
        
        # Get heterosis values for these crosses
        heterosis_values <- values$diallel_heterosis_data %>% 
          filter(CrossID %in% top_3_crosses) %>% 
          select(CrossID, MPH, BPH, Cross_Value, Mid_Parent, Better_Parent) %>% 
          arrange(desc(MPH))
        
        values$top_three_crosses_heterosis <- heterosis_values
        
        # Create combined view with SNP data
        combined_top_crosses <- values$top_3_heterotic_crosses_snp %>% 
          left_join(heterosis_values, by = c("Cross" = "CrossID")) %>% 
          arrange(desc(MPH), Trait, Chromosome, Position)
        
        values$top_crosses_snp_combined <- combined_top_crosses
      }
    }
  })
  
  # UI Output for Top 3 Heterotic Crosses
  output$top_heterotic_crosses_ui <- renderUI({
    req(values$top_3_heterotic_crosses_snp)
    
    fluidRow(
      box(
        title = "Top 3 Heterotic Crosses (by MPH)", 
        status = "success", 
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        collapsed = FALSE,
        tabBox(
          width = 12,
          tabPanel(
            "Crosses Summary",
            h4("Top 3 Crosses by Mid-Parent Heterosis (MPH):"),
            DTOutput("top_three_crosses_heterosis_table"),
            br(),
            plotOutput("top_three_crosses_heterosis_plot", height = "400px")
          ),
          tabPanel(
            "SNP Data",
            h4("SNP Information for Top 3 Crosses:"),
            DTOutput("top_crosses_snp_table"),
            br(),
            downloadButton("download_top_crosses_snp", "Download SNP Data", class = "btn-primary")
          ),
          tabPanel(
            "Summary Statistics",
            h4("Summary by Cross and Trait:"),
            DTOutput("all_crosses_summary_table"),
            br(),
            verbatimTextOutput("top_crosses_stats")
          ),
          tabPanel(
            "Visualization",
            h4("Visual Analysis of Top Heterotic Crosses:"),
            fluidRow(
              column(6,
                     selectInput("top_cross_trait", "Select Trait:",
                                 choices = unique(values$top_3_heterotic_crosses_snp$Trait)),
                     plotOutput("top_crosses_manhattan", height = "400px")
              ),
              column(6,
                     selectInput("top_cross_metric", "Select Metric:",
                                 choices = c("P_value", "Effect", "LOD")),
                     plotOutput("top_crosses_comparison", height = "400px")
              )
            )
          )
        )
      )
    )
  })
  
  # Table for heterosis values
  output$top_three_crosses_heterosis_table <- renderDT({
    req(values$top_three_crosses_heterosis)
    
    datatable(values$top_three_crosses_heterosis,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE,
              caption = "Heterosis Values for Top 3 Crosses") %>%
      formatRound(columns = c("MPH", "BPH"), digits = 2) %>%
      formatRound(columns = c("Cross_Value", "Mid_Parent", "Better_Parent"), digits = 4)
  })
  
  # Plot for heterosis comparison
  output$top_three_crosses_heterosis_plot <- renderPlot({
    req(values$top_three_crosses_heterosis)
    
    plot_data <- values$top_three_crosses_heterosis %>%
      mutate(Cross = factor(CrossID, levels = CrossID[order(MPH)]))
    
    ggplot(plot_data, aes(x = Cross)) +
      geom_bar(aes(y = MPH, fill = "MPH"), stat = "identity", alpha = 0.7, width = 0.4) +
      geom_bar(aes(y = BPH, fill = "BPH"), stat = "identity", alpha = 0.7, width = 0.4, 
               position = position_nudge(x = 0.2)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(title = "Heterosis Comparison for Top 3 Crosses",
           y = "Heterosis (%)",
           x = "Cross",
           fill = "Heterosis Type") +
      scale_fill_manual(values = c("MPH" = "steelblue", "BPH" = "darkorange")) +
      theme_minimal() +
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Table for SNP data
  output$top_crosses_snp_table <- renderDT({
    req(values$top_3_heterotic_crosses_snp)
    
    # Select important columns
    display_data <- values$top_3_heterotic_crosses_snp %>%
      select(Cross, Parent1, Parent2, SNP_ID, Chromosome, Position, 
             Trait, P_value, LOD, Effect, Significant, Parent1_Genotype, Parent2_Genotype)
    
    datatable(display_data,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE,
              caption = "SNP Data for Top 3 Heterotic Crosses") %>%
      formatRound(columns = c("P_value", "LOD", "Effect"), digits = 4) %>%
      formatStyle("Significant",
                  backgroundColor = styleEqual(c(TRUE, FALSE), c("#90EE90", "#FFB6C1")))
  })
  
  # Table for summary statistics
  output$all_crosses_summary_table <- renderDT({
    req(values$top_crosses_summary)
    
    datatable(values$top_crosses_summary,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE,
              caption = "Summary Statistics for Top 3 Crosses by Trait") %>%
      formatRound(columns = c("Avg_P_value", "Min_P_value", "Max_Effect"), digits = 4)
  })
  
  # Text output for stats
  output$top_crosses_stats <- renderPrint({
    req(values$top_3_heterotic_crosses_snp, values$top_three_crosses_heterosis)
    
    cat("=== Top 3 Heterotic Crosses Analysis ===\n\n")
    
    # Basic stats
    cat("Total SNPs analyzed:", nrow(values$top_3_heterotic_crosses_snp), "\n")
    cat("Unique crosses:", length(unique(values$top_3_heterotic_crosses_snp$Cross)), "\n")
    cat("Unique traits:", length(unique(values$top_3_heterotic_crosses_snp$Trait)), "\n")
    cat("Unique SNPs:", length(unique(values$top_3_heterotic_crosses_snp$SNP_ID)), "\n\n")
    
    # Significant SNPs
    if ("Significant" %in% colnames(values$top_3_heterotic_crosses_snp)) {
      sig_count <- sum(values$top_3_heterotic_crosses_snp$Significant, na.rm = TRUE)
      cat("Significant SNPs (p < 0.05):", sig_count, "\n")
      cat("Significant SNP percentage:", round(sig_count/nrow(values$top_3_heterotic_crosses_snp)*100, 2), "%\n\n")
    }
    
    # Heterosis stats
    cat("Heterosis Statistics:\n")
    for (i in 1:nrow(values$top_three_crosses_heterosis)) {
      cross <- values$top_three_crosses_heterosis$CrossID[i]
      mph <- values$top_three_crosses_heterosis$MPH[i]
      bph <- values$top_three_crosses_heterosis$BPH[i]
      
      cat("  ", cross, ":\n")
      cat("    MPH:", round(mph, 2), "%\n")
      cat("    BPH:", round(bph, 2), "%\n")
    }
  })
  
  # Manhattan plot for selected trait
  output$top_crosses_manhattan <- renderPlot({
    req(values$top_3_heterotic_crosses_snp, input$top_cross_trait)
    
    plot_data <- values$top_3_heterotic_crosses_snp %>%
      filter(Trait == input$top_cross_trait)
    
    if (nrow(plot_data) == 0) return(NULL)
    
    ggplot(plot_data, aes(x = Position, y = -log10(P_value), color = as.factor(Chromosome))) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
      facet_wrap(~Cross, scales = "free_x") +
      labs(title = paste("Manhattan Plot for", input$top_cross_trait),
           subtitle = "Top 3 Heterotic Crosses",
           x = "Position",
           y = "-log10(P-value)",
           color = "Chromosome") +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Comparison plot
  output$top_crosses_comparison <- renderPlot({
    req(values$top_3_heterotic_crosses_snp, input$top_cross_metric)
    
    ggplot(values$top_3_heterotic_crosses_snp, aes_string(x = "Cross", y = input$top_cross_metric)) +
      geom_boxplot(fill = "lightblue", alpha = 0.7) +
      geom_jitter(aes(color = Trait), width = 0.2, alpha = 0.5) +
      labs(title = paste("Comparison of", input$top_cross_metric, "across Top Crosses"),
           x = "Cross",
           y = input$top_cross_metric) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Download handler for top crosses SNP data
  output$download_top_crosses_snp <- downloadHandler(
    filename = function() {
      paste0("top_3_heterotic_crosses_snp_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$top_3_heterotic_crosses_snp, file, row.names = FALSE)
    }
  )
  
  # Add a button to the UI to trigger this analysis
  observeEvent(input$analyze_top_heterotic_crosses, {
    showNotification("Analyzing top 3 heterotic crosses...", type = "message", duration = 3)
    
    # Trigger the analysis
    tryCatch({
      # Get top 3 crosses
      if (!is.null(values$diallel_heterosis_data) && !is.null(values$cross_snp_df)) {
        top_3_crosses <- values$diallel_heterosis_data %>% 
          arrange(desc(MPH)) %>% 
          head(3) %>% 
          pull(CrossID)
        
        values$top_3_crosses_list <- top_3_crosses
        
        # Show success message
        showNotification(
          paste("Analysis complete. Top 3 crosses:", paste(top_3_crosses, collapse = ", ")),
          type = "message",
          duration = 5
        )
      } else {
        showNotification("Required data not available. Please run diallel and cross-SNP analysis first.", 
                         type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  
  output$top_crosses_info_box <- renderUI({
    if (!is.null(values$top_3_crosses_list)) {
      tags$div(
        class = "alert alert-success",
        icon("trophy"),
        h5("Top 3 Crosses (by MPH):"),
        tags$ul(
          tags$li(values$top_3_crosses_list[1]),
          tags$li(values$top_3_crosses_list[2]),
          tags$li(values$top_3_crosses_list[3])
        )
      )
    } else {
      tags$div(
        class = "alert alert-info",
        icon("info-circle"),
        "Click the button to identify top heterotic crosses"
      )
    }
  })
  
  
  # ==========================================================================
  # OUTPUT RENDERERS - DIALLEL SELECTION TAB
  # ==========================================================================
  
  # Best parent value box
  output$diallel_best_parent <- renderValueBox({
    req(values$diallel_results$gca)
    
    best_parent <- values$diallel_results$gca %>%
      arrange(desc(GCA)) %>%
      slice(1)
    
    valueBox(
      value = best_parent$Parent,
      subtitle = "Best GCA",
      icon = icon("seedling"),
      color = "green"
    )
  })
  
  # Best cross value box
  output$diallel_best_cross <- renderValueBox({
    req(values$diallel_heterosis_data)
    
    tryCatch({
      if (!is.null(values$diallel_heterosis_data) && nrow(values$diallel_heterosis_data) > 0) {
        best_cross <- values$diallel_heterosis_data %>%
          arrange(desc(Cross_Value)) %>%
          slice(1)
        
        if (nrow(best_cross) == 0) {
          valueBox(
            value = "N/A",
            subtitle = "Best Performing Cross",
            icon = icon("chart-line"),
            color = "blue"
          )
        } else {
          valueBox(
            value = as.character(best_cross$CrossID[1]),
            subtitle = "Best Performing Cross",
            icon = icon("chart-line"),
            color = "blue"
          )
        }
      } else {
        valueBox(
          value = "N/A",
          subtitle = "Best Performing Cross",
          icon = icon("chart-line"),
          color = "blue"
        )
      }
    }, error = function(e) {
      valueBox(
        value = "Error",
        subtitle = "Best Performing Cross",
        icon = icon("exclamation-triangle"),
        color = "red"
      )
    })
  })
  
  # Best heterosis value box
  output$diallel_best_heterosis <- renderValueBox({
    req(values$diallel_heterosis_data)
    
    tryCatch({
      if (!is.null(values$diallel_heterosis_data) && nrow(values$diallel_heterosis_data) > 0) {
        best_het <- values$diallel_heterosis_data %>%
          arrange(desc(MPH)) %>%
          slice(1)
        
        if (nrow(best_het) == 0) {
          valueBox(
            value = "N/A",
            subtitle = "Highest Heterosis",
            icon = icon("arrow-up"),
            color = "orange"
          )
        } else {
          valueBox(
            value = paste0(round(best_het$MPH[1], 1), "%"),
            subtitle = "Highest Heterosis",
            icon = icon("arrow-up"),
            color = "orange"
          )
        }
      } else {
        valueBox(
          value = "N/A",
          subtitle = "Highest Heterosis",
          icon = icon("arrow-up"),
          color = "orange"
        )
      }
    }, error = function(e) {
      valueBox(
        value = "Error",
        subtitle = "Highest Heterosis",
        icon = icon("exclamation-triangle"),
        color = "red"
      )
    })
  })
  
  # Heritability value box
  output$diallel_heritability <- renderValueBox({
    req(values$diallel_results$anova)
    
    if (!is.null(values$diallel_results$anova$ss_gca)) {
      # Simple heritability estimate
      h2 <- values$diallel_results$anova$ss_gca / 
        (values$diallel_results$anova$ss_gca + values$diallel_results$anova$ms_error)
      
      valueBox(
        value = paste0(round(h2 * 100, 1), "%"),
        subtitle = "Heritability",
        icon = icon("dna"),
        color = "purple"
      )
    } else {
      valueBox(
        value = "N/A",
        subtitle = "Heritability",
        icon = icon("dna"),
        color = "purple"
      )
    }
  })
  
  # Top parents table
  output$diallel_top_parents <- renderDT({
    req(values$diallel_results$gca)
    
    top_parents <- values$diallel_results$gca %>%
      arrange(desc(GCA)) %>%
      head(10)
    
    datatable(
      top_parents,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Top 10 Parents by GCA"
    ) %>%
      formatRound(columns = c('GCA', 'Mean_Performance'), digits = 3)
  })
  
  # Top crosses table
  output$diallel_top_crosses <- renderDT({
    req(values$diallel_heterosis_data)
    
    top_crosses <- values$diallel_heterosis_data %>%
      arrange(desc(Cross_Value)) %>%
      head(10)
    
    datatable(
      top_crosses[, c("CrossID", "Parent1", "Parent2", "Cross_Value", "MPH")],
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Top 10 Crosses by Performance"
    ) %>%
      formatRound(columns = c('Cross_Value', 'MPH'), digits = 2)
  })
  
  # Top heterotic crosses
  output$diallel_top_heterotic <- renderDT({
    req(values$diallel_heterosis_data)
    
    top_heterotic <- values$diallel_heterosis_data %>%
      arrange(desc(MPH)) %>%
      head(10)
    
    datatable(
      top_heterotic[, c("CrossID", "Parent1", "Parent2", "MPH", "BPH")],
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Top 10 Crosses by Heterosis"
    ) %>%
      formatRound(columns = c('MPH', 'BPH'), digits = 2)
  })
  
  # Diallel recommendations
  output$diallel_recommendations <- renderUI({
    results <- values$diallel_results
    h_data <- values$diallel_heterosis_data
    
    if(is.null(results) || !is.null(results$error)) {
      return(HTML(paste0("<h4>Analysis Error:</h4><p>", results$error, "</p>")))
    }
    
    if((is.null(results$gca) || nrow(results$gca) == 0) && (is.null(h_data) || nrow(h_data) == 0)) {
      return(HTML("<h4>Insufficient Data:</h4><p>Not enough data to generate recommendations.</p>"))
    }
    
    # Get top performers (with error handling)
    top_parents <- if(!is.null(results$gca) && nrow(results$gca) > 0) {
      head(results$gca[order(-results$gca$GCA), "Parent"], 3)
    } else { character(0) }
    
    top_crosses <- if(!is.null(values$diallel_cross_data) && nrow(values$diallel_cross_data) > 0) {
      head(values$diallel_cross_data[order(-values$diallel_cross_data[[input$diallel_trait]]), "CrossID"], 3)
    } else { character(0) }
    
    top_heterotic <- if(!is.null(h_data) && nrow(h_data) > 0) {
      head(h_data[order(-h_data$MPH), "CrossID"], 3)
    } else { character(0) }
    
    # Calculate statistics
    high_heterosis <- if(!is.null(h_data)) sum(h_data$MPH > 20, na.rm = TRUE) else 0
    total_crosses <- if(!is.null(h_data)) nrow(h_data) else 0
    
    HTML(paste(
      "<h4>Comprehensive Selection Strategy for", input$diallel_trait, "</h4>",
      
      if(length(top_parents) > 0) paste(
        "<h5>🎯 Top General Combiners (Foundation Breeding):</h5>",
        "<ul>",
        paste0("<li><strong>", top_parents, "</strong> - Superior general combining ability</li>", collapse = ""),
        "</ul>"
      ) else "<h5>🎯 No parent data available for GCA analysis</h5>",
      
      if(length(top_crosses) > 0) paste(
        "<h5>🏆 Top Performing Crosses (Immediate Testing):</h5>",
        "<ul>",
        paste0("<li><strong>", top_crosses, "</strong> - High mean performance</li>", collapse = ""),
        "</ul>"
      ) else "<h5>🏆 No cross performance data available</h5>",
      
      if(length(top_heterotic) > 0) paste(
        "<h5>🚀 Top Heterotic Crosses (Hybrid Development):</h5>",
        "<ul>",
        paste0("<li><strong>", top_heterotic, "</strong> - Exceptional heterosis</li>", collapse = ""),
        "</ul>"
      ) else "<h5>🚀 No heterosis data available</h5>",
      
      "<h5>📊 Breeding Insights:</h5>",
      "<ul>",
      if(total_crosses > 0) paste("<li><strong>", high_heterosis, " out of ", total_crosses, 
                                  "</strong> crosses show high heterosis (>20%)</li>") else 
                                    "<li>No heterosis data available</li>",
      if(length(top_parents) > 0) paste("<li><strong>", top_parents[1], 
                                        "</strong> is the most reliable parent for general combining ability</li>") else
                                          "<li>No parent combining ability data available</li>",
      "<li>Crosses combining high GCA parents with high SCA show maximum potential</li>",
      "</ul>",
      
      "<h5>🎯 Strategic Recommendations:</h5>",
      "<ul>",
      if(length(top_crosses) > 0) paste("<li><strong>Short-term:</strong> Advance <strong>", top_crosses[1], 
                                        "</strong> for multi-location trials</li>") else
                                          "<li><strong>Short-term:</strong> No specific cross recommendations available</li>",
      if(length(top_parents) > 0) paste("<li><strong>Medium-term:</strong> Use <strong>", top_parents[1], 
                                        "</strong> in multiple crosses</li>") else
                                          "<li><strong>Medium-term:</strong> No parent recommendations available</li>",
      "<li><strong>Long-term:</strong> Develop synthetic varieties from top general combiners</li>",
      if(length(top_heterotic) > 0) paste("<li><strong>Hybrid Program:</strong> Focus on <strong>", top_heterotic[1], 
                                          "</strong> for hybrid development</li>") else
                                            "<li><strong>Hybrid Program:</strong> No heterotic cross recommendations available</li>",
      "</ul>"
    ))
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS - DIALLEL VISUALIZATIONS TAB
  # ==========================================================================
  
  # GCA plot
  output$diallel_gca_plot <- renderPlotly({
    req(values$diallel_results$gca)
    
    p <- ggplot(values$diallel_results$gca, 
                aes(x = reorder(Parent, GCA), y = GCA, 
                    text = paste("Parent:", Parent, "<br>GCA:", round(GCA, 3)))) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      coord_flip() +
      labs(
        title = "General Combining Ability (GCA) Effects",
        x = "Parent",
        y = "GCA"
      ) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10))
    
    ggplotly(p, tooltip = "text")
  })
  
  
  output$diallel_sca_heatmap <- renderPlotly({
    results <- values$diallel_results
    if(is.null(results) || !is.null(results$error) || is.null(results$sca)) {
      return(plotly_empty(type = "scatter") %>% layout(title = "No SCA data available"))
    }
    sca_df <- as.data.frame(as.table(results$sca))
    colnames(sca_df) <- c("Parent1", "Parent2", "SCA")
    sca_df <- sca_df[!is.na(sca_df$SCA), ]
    
    if(nrow(sca_df) == 0) {
      return(plotly_empty(type = "scatter") %>% layout(title = "No SCA data available"))
    }
    
    p <- ggplot(sca_df, aes(x = Parent1, y = Parent2, fill = SCA)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
      labs(title = paste("SCA Effects -", input$diallel_trait)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggplotly(p)
  })
  
  
  
  # Mean performance vs GCA
  output$diallel_mean_vs_gca <- renderPlotly({
    req(values$diallel_results$gca)
    
    p <- ggplot(values$diallel_results$gca, 
                aes(x = GCA, y = Mean_Performance, 
                    text = paste("Parent:", Parent, 
                                 "<br>GCA:", round(GCA, 3),
                                 "<br>Mean:", round(Mean_Performance, 2)))) +
      geom_point(size = 3, color = "darkgreen", alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
      labs(
        title = "Mean Performance vs GCA",
        x = "GCA Effect",
        y = "Mean Performance"
      ) +
      theme_minimal()
    
    ggplotly(p, tooltip = "text")
  })
  
  # Heterosis vs SCA
  output$diallel_heterosis_vs_sca <- renderPlotly({
    req(values$diallel_heterosis_data, values$diallel_results$sca)
    
    # Merge heterosis and SCA data
    if (!is.null(values$diallel_heterosis_data) && !is.null(values$diallel_results$sca)) {
      het_data <- values$diallel_heterosis_data
      sca_matrix <- values$diallel_results$sca
      
      # Add SCA values to heterosis data
      het_data$SCA <- NA
      for (i in 1:nrow(het_data)) {
        p1 <- het_data$Parent1[i]
        p2 <- het_data$Parent2[i]
        if (p1 %in% rownames(sca_matrix) && p2 %in% colnames(sca_matrix)) {
          het_data$SCA[i] <- sca_matrix[p1, p2]
        }
      }
      
      p <- ggplot(het_data, aes(x = SCA, y = MPH, 
                                text = paste("Cross:", CrossID,
                                             "<br>SCA:", round(SCA, 3),
                                             "<br>MPH:", round(MPH, 1), "%"))) +
        geom_point(size = 3, color = "purple", alpha = 0.7) +
        geom_smooth(method = "lm", se = FALSE, color = "orange", linetype = "dashed") +
        labs(
          title = "Heterosis vs SCA",
          x = "Specific Combining Ability (SCA)",
          y = "Mid-Parent Heterosis (MPH, %)"
        ) +
        theme_minimal()
      
      ggplotly(p, tooltip = "text")
    }
  })
  
  # ==========================================================================
  # OUTPUT RENDERERS - METAN ANALYSIS TAB
  # ==========================================================================
  
  # METAN trait selectors
  output$metan_trait_selector_anova <- renderUI({
    req(input$metan_resp_vars)
    selectInput("metan_anova_trait", "Select Trait for ANOVA",
                choices = input$metan_resp_vars,
                selected = input$metan_resp_vars[1],
                width = "100%")
  })
  
  output$metan_trait_selector_ammi <- renderUI({
    req(input$metan_resp_vars)
    selectInput("metan_ammi_trait", "Select Trait for AMMI",
                choices = input$metan_resp_vars,
                selected = input$metan_resp_vars[1],
                width = "100%")
  })
  
  output$metan_trait_selector_gge <- renderUI({
    req(input$metan_resp_vars)
    selectInput("metan_gge_trait", "Select Trait for GGE",
                choices = input$metan_resp_vars,
                selected = input$metan_resp_vars[1],
                width = "100%")
  })
  
  output$metan_trait_selector_ge <- renderUI({
    req(input$metan_resp_vars)
    selectInput("metan_ge_trait", "Select Trait for GE Plot",
                choices = input$metan_resp_vars,
                selected = input$metan_resp_vars[1],
                width = "100%")
  })
  
  # METAN descriptive stats table
  output$metan_desc_stats_table <- renderDT({
    req(values$metan_results$desc_stats)
    
    if (!is.character(values$metan_results$desc_stats)) {
      datatable(
        values$metan_results$desc_stats,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = TRUE,
        caption = "Descriptive Statistics by Genotype"
      )
    }
  })
  
  # METAN ANOVA output
  output$metan_anova_output <- renderPrint({
    req(values$metan_results$anova, input$metan_anova_trait)
    
    cat("=== ANOVA RESULTS FOR", input$metan_anova_trait, "===\n\n")
    
    if (!is.null(values$metan_results$anova[[input$metan_anova_trait]])) {
      if (!is.character(values$metan_results$anova[[input$metan_anova_trait]])) {
        print(values$metan_results$anova[[input$metan_anova_trait]])
      } else {
        cat(values$metan_results$anova[[input$metan_anova_trait]])
      }
    } else {
      cat("No ANOVA results available for this trait\n")
    }
  })
  
  # METAN AMMI output
  output$metan_ammi_output <- renderPrint({
    req(values$metan_results$ammi, input$metan_ammi_trait)
    
    cat("=== AMMI ANALYSIS FOR", input$metan_ammi_trait, "===\n\n")
    
    if (!is.null(values$metan_results$ammi[[input$metan_ammi_trait]])) {
      if (!is.character(values$metan_results$ammi[[input$metan_ammi_trait]])) {
        print(summary(values$metan_results$ammi[[input$metan_ammi_trait]]))
      } else {
        cat(values$metan_results$ammi[[input$metan_ammi_trait]])
      }
    } else {
      cat("No AMMI results available for this trait\n")
    }
  })
  
  # METAN AMMI plot
  output$metan_ammi_plot <- renderPlot({
    req(values$metan_results$ammi, input$metan_ammi_trait)
    
    if (!is.null(values$metan_results$ammi[[input$metan_ammi_trait]]) &&
        !is.character(values$metan_results$ammi[[input$metan_ammi_trait]])) {
      
      # Extract AMMI scores
      ammi_res <- values$metan_results$ammi[[input$metan_ammi_trait]]
      
      # Create AMMI biplot
      plot(ammi_res)
      
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "AMMI plot not available", size = 6) +
        theme_void()
    }
  })
  
  # METAN WAASB output
  output$metan_waasb_output <- renderPrint({
    req(values$metan_results$waasb)
    
    cat("=== WAASB ANALYSIS ===\n\n")
    
    if (!is.character(values$metan_results$waasb)) {
      print(summary(values$metan_results$waasb))
    } else {
      cat(values$metan_results$waasb)
    }
  })
  
  # METAN WAASB plot
  output$metan_waasb_plot <- renderPlot({
    req(values$metan_results$waasb)
    
    if (!is.character(values$metan_results$waasb)) {
      plot(values$metan_results$waasb)
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "WAASB plot not available", size = 6) +
        theme_void()
    }
  })
  
  # METAN GGE plot
  output$metan_gge_plot <- renderPlot({
    req(values$metan_results$gge, input$metan_gge_trait)
    
    if (!is.null(values$metan_results$gge[[input$metan_gge_trait]]) &&
        !is.character(values$metan_results$gge[[input$metan_gge_trait]])) {
      
      # Create GGE biplot
      plot(values$metan_results$gge[[input$metan_gge_trait]])
      
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "GGE biplot not available", size = 6) +
        theme_void()
    }
  })
  
  # METAN GE plot
  output$metan_ge_plot <- renderPlot({
    req(values$metan_results$ge_plot, input$metan_ge_trait)
    
    if (!is.null(values$metan_results$ge_plot[[input$metan_ge_trait]]) &&
        !is.character(values$metan_results$ge_plot[[input$metan_ge_trait]])) {
      
      # Create GE interaction plot
      plot(values$metan_results$ge_plot[[input$metan_ge_trait]])
      
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "GE interaction plot not available", size = 6) +
        theme_void()
    }
  })
  
  # METAN correlation plot
  output$metan_corr_plot <- renderPlot({
    req(values$metan_results$corr)
    
    if (!is.character(values$metan_results$corr)) {
      plot(values$metan_results$corr)
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Correlation plot not available", size = 6) +
        theme_void()
    }
  })
  
  # METAN stability table
  output$metan_stability_table <- renderDT({
    req(values$metan_results$stability)
    
    if (!is.character(values$metan_results$stability)) {
      datatable(
        as.data.frame(values$metan_results$stability),
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = TRUE,
        caption = "Stability Analysis Results"
      )
    }
  })
  
  # METAN ranking table
  output$metan_ranking_table <- renderDT({
    req(values$metan_results$ranking)
    
    if (!is.character(values$metan_results$ranking)) {
      datatable(
        as.data.frame(values$metan_results$ranking),
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = TRUE,
        caption = "Genotype Ranking Results"
      )
    }
  })
  
  # METAN environment stats
  output$metan_env_stats_ui <- renderUI({
    req(values$metan_processed_data)
    
    env_stats <- values$metan_processed_data %>%
      group_by(ENV) %>%
      summarise(
        n = n(),
        .groups = 'drop'
      )
    
    div(
      h5("Environment Statistics:"),
      tags$ul(
        lapply(1:nrow(env_stats), function(i) {
          tags$li(paste(env_stats$ENV[i], ": ", env_stats$n[i], " observations"))
        })
      )
    )
  })
  
  
  # ==========================================================================
  # OUTPUT RENDERERS - REPORTS TAB
  # ==========================================================================
  
  # Report preview
  output$report_preview <- renderUI({
    req(values$report_ready)
    
    HTML(paste(
      "<div class='alert alert-info'>",
      "<h4>Report Preview</h4>",
      "<p><strong>Title:</strong> ", values$report_content$title, "</p>",
      "<p><strong>Sections included:</strong> ", paste(values$report_content$sections, collapse = ", "), "</p>",
      "<p><strong>Notes:</strong> ", substr(values$report_content$notes, 1, 200), 
      if (nchar(values$report_content$notes) > 200) "..." else "", "</p>",
      "</div>"
    ))
  })
  
  # ==========================================================================
  # SESSION MANAGEMENT
  # ==========================================================================
  
  # Reset all data
  # Reset all data
  observeEvent(input$reset_all, {
    values$pheno_data <- NULL
    values$gl_object <- NULL
    values$geno_data <- NULL
    values$filtered_gl <- NULL
    values$filtered_data <- NULL
    values$diallel_data <- NULL
    values$metan_raw_data <- NULL
    values$metadata <- NULL
    
    values$matched_data <- NULL
    values$processed_geno <- NULL
    values$metan_processed_data <- NULL
    values$diallel_processed <- NULL
    values$diallel_cross_data <- NULL
    values$diallel_parent_data <- NULL
    values$diallel_results <- NULL
    
    values$qc_plots <- NULL
    values$pca_results <- NULL
    values$kinship_matrix <- NULL
    values$gwas_results <- NULL
    values$multi_trait_results <- NULL
    values$metan_results <- NULL
    values$qtl_results <- NULL
    values$qtl_cross <- NULL
    
    values$multi_trait_unified_df <- NULL
    values$multi_trait_unified_tidy <- NULL
    values$multi_trait_unified_wide <- NULL
    values$multi_trait_enhanced_unified_df <- NULL
    values$multi_trait_complete_integrated_df <- NULL
    
    values$cross_snp_analysis <- NULL
    values$cross_summary <- NULL
    values$top_crosses_snp_combined <- NULL
    
    values$parent_snp_actual <- NULL
    values$parent_genotype_patterns <- NULL
    values$network_data <- NULL
    
    values$pheno_loaded <- FALSE
    values$geno_loaded <- FALSE
    values$data_matched <- FALSE
    values$qc_completed <- FALSE
    
    values$report_content <- NULL
    values$report_ready <- FALSE
    
    output$pheno_status <- renderUI({})
    output$geno_status <- renderUI({})
    output$metadata_status <- renderUI({})
    output$diallel_status <- renderUI({})
    output$metan_status <- renderUI({})
    output$match_status <- renderUI({})
    output$compatibility_status <- renderUI({})
    output$report_status <- renderUI({})
    
    showNotification("All data has been reset!", type = "message", duration = 3)
  })
  
  
  # Analysis status
  output$analysis_status <- renderPrint({
    cat("=== ANALYSIS STATUS ===\n\n")
    cat("Phenotype:", ifelse(values$pheno_loaded, "✅", "❌"), "\n")
    cat("Genotype:", ifelse(values$geno_loaded, "✅", "❌"), "\n")
    cat("Data matched:", ifelse(values$data_matched, "✅", "❌"), "\n")
    cat("QC completed:", ifelse(values$qc_completed, "✅", "❌"), "\n")
    cat("GWAS:", ifelse(!is.null(values$gwas_results), "✅", "❌"), "\n")
    cat("Multi-trait:", ifelse(!is.null(values$multi_trait_results), "✅", "❌"), "\n")
    cat("Diallel:", ifelse(!is.null(values$diallel_results), "✅", "❌"), "\n")
    cat("METAN:", ifelse(!is.null(values$metan_results), "✅", "❌"), "\n")
    cat("Cross-SNP:", ifelse(!is.null(values$cross_snp_analysis), "✅", "❌"), "\n")
    cat("Unified Data:", ifelse(!is.null(values$multi_trait_unified_df), "✅", "❌"), "\n")
  })
  
  
  # Download HTML report
  output$download_report_html <- downloadHandler(
    filename = function() {
      paste0("analysis_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      if (!is.null(values$report_path)) {
        file.copy(values$report_path, file)
      }
    }
  )
  
  # Download PDF report
  output$download_report_pdf <- downloadHandler(
    filename = function() {
      paste0("analysis_report_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 8.5, height = 11)
      plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1),
           main = input$report_title)
      text(0.5, 0.8, paste("Report generated on:", Sys.Date()), cex = 1.2)
      text(0.5, 0.7, "This is a placeholder PDF report.", cex = 1)
      text(0.5, 0.6, "For full HTML report, use the HTML download option.", cex = 0.8)
      dev.off()
    }
  )
  
  
  
  # Plot all traits from parent_snp_actual summary table
  output$parent_all_traits_plot <- renderPlot({
    req(values$parent_snp_actual$summary_table)
    
    summary_table <- values$parent_snp_actual$summary_table
    
    # Identify trait columns (columns that start with "Mean_" or contain trait names)
    trait_cols <- grep("^Mean_|GYP|GYR|DTF|DTM|PL|NPP|NSP|X100SW|TW", 
                       colnames(summary_table), value = TRUE)
    
    if (length(trait_cols) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No trait data available in summary table", size = 6) +
               theme_void())
    }
    
    # Prepare data for plotting
    plot_data <- summary_table %>%
      select(Parent, all_of(trait_cols)) %>%
      pivot_longer(cols = -Parent, names_to = "Trait", values_to = "Value") %>%
      filter(!is.na(Value))
    
    # Create faceted bar plot
    ggplot(plot_data, aes(x = reorder(Parent, Value), y = Value, fill = Trait)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~Trait, scales = "free_y", ncol = 3) +
      coord_flip() +
      labs(
        title = "Parent Performance Across All Traits",
        x = "Parent",
        y = "Trait Value"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none"
      ) +
      scale_fill_brewer(palette = "Set3")
  })
  
  # Display the summary table with all traits
  output$parent_all_traits_table <- renderDT({
    req(values$parent_snp_actual$summary_table)
    
    datatable(
      values$parent_snp_actual$summary_table,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Parent-SNP Summary with All Traits",
      extensions = 'Buttons',
      filter = 'top'
    )
  })
  
  
  # Download handler for all traits summary
  output$download_all_traits_summary <- downloadHandler(
    filename = function() {
      paste0("parent_all_traits_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$parent_snp_actual$summary_table)
      write.csv(values$parent_snp_actual$summary_table, file, row.names = FALSE)
    }
  )
  
  
  # Download multi-trait results
  output$download_multi_trait <- downloadHandler(
    filename = function() {
      paste0("multi_trait_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Create temporary directory
      temp_dir <- tempdir()
      results_dir <- file.path(temp_dir, "multi_trait_results")
      dir.create(results_dir, showWarnings = FALSE)
      
      # Save various results
      if (!is.null(values$multi_trait_results$trait_correlations)) {
        write.csv(values$multi_trait_results$trait_correlations$correlations,
                  file.path(results_dir, "trait_correlations.csv"),
                  row.names = TRUE)
        
        if (!is.null(values$multi_trait_results$trait_correlations$p_values)) {
          write.csv(values$multi_trait_results$trait_correlations$p_values,
                    file.path(results_dir, "correlation_pvalues.csv"),
                    row.names = TRUE)
        }
      }
      
      if (!is.null(values$multi_trait_results$combined_results)) {
        write.csv(values$multi_trait_results$combined_results,
                  file.path(results_dir, "combined_gwas_results.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$multi_trait_results$pleiotropic_snps)) {
        write.csv(values$multi_trait_results$pleiotropic_snps,
                  file.path(results_dir, "pleiotropic_snps.csv"),
                  row.names = FALSE)
      }
      
      # Save individual trait results
      if (!is.null(values$multi_trait_results$trait_gwas_results)) {
        for (trait in names(values$multi_trait_results$trait_gwas_results)) {
          if (!is.null(values$multi_trait_results$trait_gwas_results[[trait]])) {
            write.csv(values$multi_trait_results$trait_gwas_results[[trait]],
                      file.path(results_dir, paste0("gwas_", trait, ".csv")),
                      row.names = FALSE)
          }
        }
      }
      
      # Create a summary text file
      summary_file <- file.path(results_dir, "analysis_summary.txt")
      sink(summary_file)
      cat("=== MULTI-TRAIT ANALYSIS SUMMARY ===\n\n")
      cat("Date:", Sys.Date(), "\n")
      cat("Traits analyzed:", paste(values$multi_trait_results$traits_analyzed, 
                                    collapse = ", "), "\n")
      sink()
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "multi_trait_results.zip")
      zip(zip_file, files = list.files(results_dir, full.names = TRUE), extras = "-j")
      
      # Copy to download location
      file.copy(zip_file, file)
    }
  )
  
  # View multi-trait details
  observeEvent(input$view_multi_details, {
    req(values$multi_trait_results)
    
    showModal(modalDialog(
      title = "Multi-Trait Analysis Details",
      size = "l",
      easyClose = TRUE,
      
      tabsetPanel(
        tabPanel("Trait Correlations",
                 plotOutput("modal_cor_plot", height = "500px")),
        tabPanel("Pleiotropic SNPs",
                 DTOutput("modal_pleiotropy_table")),
        tabPanel("Combined Results",
                 DTOutput("modal_combined_table"))
      )
    ))
  })
  
  # Modal plot for correlations
  output$modal_cor_plot <- renderPlot({
    req(values$multi_trait_results$trait_correlations)
    output$trait_cor_plot()
  })
  
  # Modal table for pleiotropy
  output$modal_pleiotropy_table <- renderDT({
    req(values$multi_trait_results$pleiotropic_snps)
    output$pleiotropy_table()
  })
  
  # Modal table for combined results
  output$modal_combined_table <- renderDT({
    req(values$multi_trait_results$combined_results)
    output$multi_trait_gwas_table()
  })
  
  
  
  # ==========================================================================
  # QTL ANALYSIS SECTION
  # ==========================================================================
  
  # QTL analysis observer
  observeEvent(input$run_qtl_analysis, {
    req(values$gwas_results, values$matched_data, values$pheno_data)
    
    showModal(modalDialog(
      title = "Running QTL Analysis",
      tagList(
        tags$div(
          style = "text-align: center;",
          icon("spinner", class = "fa-spin fa-2x"),
          tags$h4("Preparing data for QTL analysis..."),
          tags$div(
            id = "qtl_analysis_details",
            style = "margin-top: 10px; font-size: 12px; color: #666;"
          )
        )
      ),
      footer = NULL,
      easyClose = FALSE,
      size = "s"
    ))
    
    tryCatch({
      shinyjs::html("qtl_analysis_details", "Step 1/5: Preparing genotype and phenotype data...")
      Sys.sleep(0.5)
      
      # Extract required data
      gwas_results <- values$gwas_results
      pheno_data <- values$matched_data$pheno
      geno_matrix <- values$matched_data$geno
      snp_info <- values$processed_geno$snp_info
      
      # Check if we have required columns
      if (!all(c("SNP_ID", "Chromosome", "Position", "P_value", "Effect") %in% colnames(gwas_results))) {
        removeModal()
        showNotification("GWAS results missing required columns (SNP_ID, Chromosome, Position, P_value, Effect)", 
                         type = "error")
        return()
      }
      
      shinyjs::html("qtl_analysis_details", "Step 2/5: Converting to QTL format...")
      Sys.sleep(0.5)
      
      # Convert data to QTL format
      qtl_data <- prepare_qtl_data(
        gwas_results = gwas_results,
        pheno_data = pheno_data,
        geno_matrix = geno_matrix,
        snp_info = snp_info,
        trait = input$qtl_trait,
        p_threshold = input$qtl_p_threshold
      )
      
      if (is.null(qtl_data)) {
        removeModal()
        showNotification("Failed to prepare QTL data", type = "error")
        return()
      }
      
      shinyjs::html("qtl_analysis_details", "Step 3/5: Creating cross object...")
      Sys.sleep(0.5)
      
      # Create QTL cross object
      cross_obj <- create_qtl_cross_object(
        genotype_data = qtl_data$genotypes,
        phenotype_data = qtl_data$phenotypes,
        marker_info = qtl_data$marker_info
      )
      
      if (is.null(cross_obj)) {
        removeModal()
        showNotification("Failed to create QTL cross object", type = "error")
        return()
      }
      
      shinyjs::html("qtl_analysis_details", "Step 4/5: Running QTL analysis...")
      Sys.sleep(0.5)
      
      # Run QTL analysis
      qtl_results <- run_qtl_analysis(
        cross = cross_obj,
        method = input$qtl_method,
        n_perm = input$qtl_n_perm,
        model = input$qtl_model,
        pheno_col = 2  # Assuming phenotype is in column 2
      )
      
      shinyjs::html("qtl_analysis_details", "Step 5/5: Finalizing results...")
      Sys.sleep(0.5)
      
      # Store results
      values$qtl_results <- qtl_results
      values$qtl_cross <- cross_obj
      
      removeModal()
      
      showNotification(
        HTML(paste(
          "✅ QTL analysis completed successfully!<br>",
          "Method:", input$qtl_method, "<br>",
          "Significant QTLs found:", ifelse(!is.null(qtl_results$significant_qtl), 
                                            nrow(qtl_results$significant_qtl), 0)
        )),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Error in QTL analysis:", e$message),
        type = "error",
        duration = 10
      )
      cat("QTL analysis error:", e$message, "\n")
      cat("Traceback:\n")
      print(traceback())
    })
  })
  
  # QTL trait selector
  output$qtl_trait_selector <- renderUI({
    req(values$pheno_data)
    
    traits <- setdiff(colnames(values$pheno_data), "Genotype")
    selectInput("qtl_trait", "Select Trait for QTL Analysis",
                choices = traits,
                selected = traits[1],
                width = "100%")
  })
  
  # QTL results summary
  output$qtl_results_summary <- renderPrint({
    req(values$qtl_results)
    
    cat("=== QTL ANALYSIS RESULTS ===\n\n")
    
    result <- values$qtl_results
    
    # Basic information
    cat("Method:", result$method, "\n")
    cat("Analysis date:", result$timestamp, "\n")
    cat("Number of permutations:", result$n_perm, "\n")
    cat("Significance threshold:", result$threshold, "\n")
    
    if (!is.null(result$scan_results)) {
      cat("\n=== SCANONE RESULTS ===\n")
      cat("Total positions scanned:", nrow(result$scan_results), "\n")
      
      # Get maximum LOD score
      max_lod <- max(result$scan_results$lod, na.rm = TRUE)
      cat("Maximum LOD score:", round(max_lod, 2), "\n")
      
      # Count significant peaks
      if (!is.null(result$significant_qtl)) {
        cat("Significant QTLs:", nrow(result$significant_qtl), "\n")
      }
    }
    
    if (!is.null(result$significant_qtl) && nrow(result$significant_qtl) > 0) {
      cat("\n=== SIGNIFICANT QTLs ===\n")
      print(result$significant_qtl)
    }
  })
  
  # QTL plot (LOD profile)
  output$qtl_plot <- renderPlot({
    req(values$qtl_results)
    
    result <- values$qtl_results
    
    if (is.null(result$scan_results)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No scan results available", size = 6) +
               theme_void())
    }
    
    # Create LOD profile plot
    plot_qtl_lod_profile(result$scan_results, 
                         threshold = result$threshold,
                         title = paste("QTL LOD Profile -", input$qtl_trait))
  })
  
  # QTL effect plot
  output$qtl_effect_plot <- renderPlot({
    req(values$qtl_results, values$qtl_cross)
    
    result <- values$qtl_results
    cross_obj <- values$qtl_cross
    
    if (is.null(result$significant_qtl) || nrow(result$significant_qtl) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No significant QTLs to plot", size = 6) +
               theme_void())
    }
    
    # Plot QTL effects
    plot_qtl_effects(cross_obj, 
                     qtl_results = result,
                     trait = input$qtl_trait)
  })
  
  # QTL results table
  output$qtl_results_table <- renderDT({
    req(values$qtl_results)
    
    result <- values$qtl_results
    
    if (is.null(result$scan_results)) {
      return(datatable(
        data.frame(Message = "No QTL results available"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Format the scan results
    display_df <- result$scan_results %>%
      mutate(
        lod = round(lod, 3),
        p_value = format(p_value, scientific = TRUE, digits = 3)
      )
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "QTL Scan Results",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'lod',
        background = styleColorBar(range(display_df$lod, na.rm = TRUE), 'lightgreen'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  # QTL peaks table
  output$qtl_peaks_table <- renderDT({
    req(values$qtl_results)
    
    result <- values$qtl_results
    
    if (is.null(result$significant_qtl) || nrow(result$significant_qtl) == 0) {
      return(datatable(
        data.frame(Message = "No significant QTL peaks found"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Format significant QTLs
    display_df <- result$significant_qtl %>%
      mutate(
        lod = round(lod, 3),
        ci_low = round(ci_low, 2),
        ci_high = round(ci_high, 2),
        p_value = format(p_value, scientific = TRUE, digits = 3),
        effect = round(effect, 4),
        pve = round(pve, 2)
      )
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Significant QTL Peaks",
      extensions = 'Buttons',
      filter = 'top'
    )
  })
  
  # QTL diagnostics plot
  output$qtl_diagnostics_plot <- renderPlot({
    req(values$qtl_results)
    
    result <- values$qtl_results
    
    if (is.null(result$perm_results)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No permutation results available", size = 6) +
               theme_void())
    }
    
    # Plot permutation results
    plot_qtl_permutations(result$perm_results,
                          threshold = result$threshold)
  })
  
  # Download QTL results
  output$download_qtl_results <- downloadHandler(
    filename = function() {
      paste0("qtl_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(values$qtl_results)
      
      # Create temporary directory
      temp_dir <- tempfile("qtl_analysis_")
      dir.create(temp_dir)
      
      # Save all results
      if (!is.null(values$qtl_results$scan_results)) {
        write.csv(values$qtl_results$scan_results,
                  file.path(temp_dir, "qtl_scan_results.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$qtl_results$significant_qtl)) {
        write.csv(values$qtl_results$significant_qtl,
                  file.path(temp_dir, "significant_qtl.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$qtl_results$perm_results)) {
        write.csv(values$qtl_results$perm_results,
                  file.path(temp_dir, "permutation_results.csv"),
                  row.names = FALSE)
      }
      
      # Save metadata
      meta_data <- data.frame(
        Parameter = c("Analysis Date", "Trait", "Method", "Threshold", 
                      "Permutations", "Significant QTLs"),
        Value = c(
          as.character(Sys.time()),
          input$qtl_trait,
          values$qtl_results$method,
          values$qtl_results$threshold,
          values$qtl_results$n_perm,
          ifelse(!is.null(values$qtl_results$significant_qtl), 
                 nrow(values$qtl_results$significant_qtl), 0)
        )
      )
      write.csv(meta_data, file.path(temp_dir, "metadata.csv"), row.names = FALSE)
      
      # Create README
      readme_content <- paste(
        "QTL ANALYSIS RESULTS",
        "====================",
        "",
        paste("Analysis Date:", Sys.time()),
        paste("Generated by: Plant Breeding Analysis Platform v4.0"),
        "",
        "Files included:",
        "1. qtl_scan_results.csv - Full scan results with LOD scores",
        "2. significant_qtl.csv - Significant QTL peaks",
        "3. permutation_results.csv - Permutation test results",
        "4. metadata.csv - Analysis parameters and metadata",
        "",
        "Column descriptions for qtl_scan_results.csv:",
        "- chr: Chromosome",
        "- pos: Position in cM",
        "- lod: LOD score",
        "- p_value: P-value",
        "- marker: Marker name",
        "",
        "Column descriptions for significant_qtl.csv:",
        "- chr: Chromosome",
        "- pos: Peak position",
        "- lod: Maximum LOD score",
        "- ci_low: Lower confidence interval",
        "- ci_high: Upper confidence interval",
        "- effect: QTL effect size",
        "- pve: Phenotypic variance explained (%)",
        sep = "\n"
      )
      writeLines(readme_content, file.path(temp_dir, "README.txt"))
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "results.zip")
      zip::zip(zip_file, files = list.files(temp_dir, full.names = TRUE), 
               mode = "cherry-pick")
      
      # Copy to download location
      file.copy(zip_file, file)
      
      # Cleanup
      unlink(temp_dir, recursive = TRUE)
    }
  )
  
  # Clear QTL results
  observeEvent(input$clear_qtl_results, {
    values$qtl_results <- NULL
    values$qtl_cross <- NULL
    
    # Reset outputs
    output$qtl_results_summary <- renderPrint({
      cat("=== QTL ANALYSIS ===\n\n")
      cat("No results available.\n")
      cat("Click 'Run QTL Analysis' to analyze.\n")
    })
    
    output$qtl_plot <- renderPlot({
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No QTL results available", size = 6) +
        theme_void()
    })
    
    output$qtl_effect_plot <- renderPlot({
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No QTL results available", size = 6) +
        theme_void()
    })
    
    output$qtl_results_table <- renderDT({NULL})
    output$qtl_peaks_table <- renderDT({NULL})
    output$qtl_diagnostics_plot <- renderPlot({
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No QTL results available", size = 6) +
        theme_void()
    })
    
    showNotification("QTL results cleared!", type = "warning", duration = 3)
  })
  
  # ==========================================================================
  # QTL ANALYSIS FUNCTIONS
  # ==========================================================================
  
  # Improved function to prepare data for QTL analysis
  prepare_qtl_data <- function(gwas_results, pheno_data, geno_matrix, snp_info, trait, p_threshold = 0.05) {
    tryCatch({
      # Filter significant SNPs from GWAS
      sig_snps <- gwas_results %>%
        filter(P_value < p_threshold) %>%
        arrange(Chromosome, Position)
      
      cat("Found", nrow(sig_snps), "significant SNPs at threshold p <", p_threshold, "\n")
      
      # If no significant SNPs, use top SNPs
      if (nrow(sig_snps) == 0) {
        cat("No significant SNPs. Using top 100 SNPs instead.\n")
        sig_snps <- gwas_results %>%
          arrange(P_value) %>%
          head(100)
      }
      
      # Get common samples
      common_samples <- intersect(rownames(geno_matrix), pheno_data$Genotype)
      
      cat("Common samples found:", length(common_samples), "\n")
      
      if (length(common_samples) < 10) {
        stop("Insufficient common samples (<10) for QTL analysis")
      }
      
      # Get SNP IDs that are in both sig_snps and geno_matrix
      sig_snp_ids <- sig_snps$SNP_ID
      available_snps <- sig_snp_ids[sig_snp_ids %in% colnames(geno_matrix)]
      
      cat("SNPs in genotype matrix:", length(available_snps), "/", length(sig_snp_ids), "\n")
      
      if (length(available_snps) == 0) {
        stop("No SNPs found in genotype matrix. Check SNP_ID matching.")
      }
      
      # Subset genotype matrix
      geno_subset <- geno_matrix[common_samples, available_snps, drop = FALSE]
      cat("Genotype subset dimensions:", dim(geno_subset), "\n")
      
      # Get phenotype values
      pheno_subset <- pheno_data %>%
        filter(Genotype %in% common_samples) %>%
        select(Genotype, !!sym(trait)) %>%
        arrange(Genotype) %>%
        filter(!is.na(!!sym(trait)))  # Remove NA phenotype values
      
      # Ensure alignment
      geno_subset <- geno_subset[pheno_subset$Genotype, , drop = FALSE]
      
      # Prepare marker info
      marker_info <- snp_info %>%
        filter(SNP_ID %in% colnames(geno_subset)) %>%
        mutate(
          chrom = as.character(Chromosome),
          pos = as.numeric(Position)
        ) %>%
        select(marker = SNP_ID, chrom, pos) %>%
        arrange(chrom, pos)
      
      cat("Markers for QTL analysis:", nrow(marker_info), "\n")
      
      # Return prepared data
      return(list(
        genotypes = geno_subset,
        phenotypes = pheno_subset,
        marker_info = marker_info,
        trait = trait,
        n_samples = nrow(pheno_subset),
        n_markers = ncol(geno_subset)
      ))
      
    }, error = function(e) {
      cat("Error preparing QTL data:", e$message, "\n")
      cat("Traceback:\n")
      print(traceback())
      return(NULL)
    })
  }
  
  # Improved function to create QTL cross object
  create_qtl_cross_object <- function(genotype_data, phenotype_data, marker_info) {
    tryCatch({
      require(qtl)
      require(dplyr)
      
      cat("Creating QTL cross object...\n")
      cat("Genotype data dimensions:", dim(genotype_data), "\n")
      cat("Phenotype data dimensions:", dim(phenotype_data), "\n")
      cat("Marker info rows:", nrow(marker_info), "\n")
      
      # Ensure marker_info markers are in genotype_data
      marker_info <- marker_info %>%
        filter(marker %in% colnames(genotype_data))
      
      if (nrow(marker_info) == 0) {
        stop("No markers found in genotype data")
      }
      
      # Convert genotype data to QTL format (0,1,2 to A,H,B)
      cat("Converting genotypes to QTL format...\n")
      
      # Create a matrix for QTL genotypes
      geno_qtl <- matrix(NA, 
                         nrow = nrow(genotype_data), 
                         ncol = ncol(genotype_data),
                         dimnames = dimnames(genotype_data))
      
      # Convert genotypes: 0=AA, 1=AB, 2=BB
      for (i in 1:nrow(genotype_data)) {
        for (j in 1:ncol(genotype_data)) {
          val <- genotype_data[i, j]
          if (is.na(val)) {
            geno_qtl[i, j] <- NA
          } else if (val == 0) {
            geno_qtl[i, j] <- "A"  # Homozygous reference
          } else if (val == 1) {
            geno_qtl[i, j] <- "H"  # Heterozygous
          } else if (val == 2) {
            geno_qtl[i, j] <- "B"  # Homozygous alternate
          }
        }
      }
      
      # Create cross object structure
      cross <- list()
      cross$geno <- list()
      
      # Set phenotype
      cross$pheno <- data.frame(
        trait = phenotype_data[[2]],
        row.names = phenotype_data$Genotype
      )
      colnames(cross$pheno) <- colnames(phenotype_data)[2]
      
      cat("Phenotype set for", nrow(cross$pheno), "individuals\n")
      
      # Organize by chromosome
      chromosomes <- unique(marker_info$chrom)
      cat("Chromosomes found:", paste(chromosomes, collapse = ", "), "\n")
      
      for (chr in chromosomes) {
        cat("Processing chromosome", chr, "...\n")
        
        # Get markers for this chromosome
        chr_markers <- marker_info %>%
          filter(chrom == chr) %>%
          arrange(pos)
        
        if (nrow(chr_markers) == 0) {
          cat("  No markers for chromosome", chr, "\n")
          next
        }
        
        cat("  Found", nrow(chr_markers), "markers\n")
        
        # Get genotype data for these markers
        chr_geno <- t(geno_qtl[, chr_markers$marker, drop = FALSE])
        
        # Create chromosome entry
        cross$geno[[chr]] <- list(
          data = chr_geno,
          map = setNames(chr_markers$pos, chr_markers$marker),
          chrtype = "A"  # Autosome
        )
        class(cross$geno[[chr]]) <- "A"
      }
      
      if (length(cross$geno) == 0) {
        stop("No chromosome data available for cross object")
      }
      
      class(cross) <- c("f2", "cross")
      attr(cross, "alleles") <- c("A", "B")
      
      cat("Cross object created successfully with", length(cross$geno), "chromosomes\n")
      
      # Calculate genotype probabilities
      cat("Calculating genotype probabilities...\n")
      cross <- calc.genoprob(cross, step = 1, error.prob = 0.001)
      
      return(cross)
      
    }, error = function(e) {
      cat("Error creating QTL cross object:", e$message, "\n")
      cat("Traceback:\n")
      print(traceback())
      return(NULL)
    })
  }
  
  # Alternative simpler QTL analysis function
  run_simple_qtl_analysis <- function(gwas_results, trait, lod_threshold = 3.0) {
    tryCatch({
      require(qtl)
      
      # Use GWAS results to simulate QTL analysis
      # Convert GWAS results to QTL scan format
      
      # Create chromosome factor
      gwas_results$chrom <- factor(gwas_results$Chromosome, 
                                   levels = sort(unique(gwas_results$Chromosome)))
      
      # Calculate LOD scores from p-values
      gwas_results$lod <- -log10(gwas_results$P_value)
      
      # Create simulated scanone results
      scan_results <- data.frame(
        chr = gwas_results$chrom,
        pos = gwas_results$Position,
        lod = gwas_results$lod,
        marker = gwas_results$SNP_ID,
        p_value = gwas_results$P_value,
        effect = gwas_results$Effect
      )
      
      # Find significant peaks
      significant_qtl <- NULL
      if (max(scan_results$lod, na.rm = TRUE) >= lod_threshold) {
        peaks <- scan_results %>%
          group_by(chr) %>%
          filter(lod == max(lod, na.rm = TRUE)) %>%
          ungroup() %>%
          filter(lod >= lod_threshold)
        
        if (nrow(peaks) > 0) {
          # Add confidence intervals (simplified)
          significant_qtl <- peaks %>%
            mutate(
              ci_low = pos - 5,  # Simplified CI
              ci_high = pos + 5,
              pve = (effect^2) * 100  # Simplified PVE calculation
            ) %>%
            select(chr, pos, lod, marker, ci_low, ci_high, effect, pve)
        }
      }
      
      # Prepare results
      results <- list(
        scan_results = scan_results,
        significant_qtl = significant_qtl,
        threshold = lod_threshold,
        method = "GWAS-based QTL scan",
        timestamp = Sys.time()
      )
      
      return(results)
      
    }, error =function(e) {
      cat("Error in simple QTL analysis:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Simplified QTL analysis observer (use this if the full analysis fails)
  observeEvent(input$run_simple_qtl, {
    req(values$gwas_results, input$qtl_trait)
    
    showModal(modalDialog(
      title = "Running Simple QTL Analysis",
      "Converting GWAS results to QTL format...",
      footer = NULL,
      easyClose = FALSE,
      size = "s"
    ))
    
    tryCatch({
      # Run simplified analysis
      qtl_results <- run_simple_qtl_analysis(
        gwas_results = values$gwas_results,
        trait = input$qtl_trait,
        lod_threshold = input$qtl_lod_threshold
      )
      
      if (is.null(qtl_results)) {
        removeModal()
        showNotification("Failed to run QTL analysis", type = "error")
        return()
      }
      
      # Store results
      values$qtl_results <- qtl_results
      
      removeModal()
      
      showNotification(
        HTML(paste(
          "✅ QTL analysis completed successfully!<br>",
          "Method: GWAS-based QTL scan<br>",
          "Significant QTLs found:", ifelse(!is.null(qtl_results$significant_qtl), 
                                            nrow(qtl_results$significant_qtl), 0)
        )),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Error in QTL analysis:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Debug function to check data
  observeEvent(input$debug_qtl_data, {
    req(values$gwas_results, values$matched_data, values$pheno_data)
    
    output$qtl_debug_info <- renderPrint({
      cat("=== QTL DATA DEBUG INFO ===\n\n")
      
      # Check GWAS results
      cat("GWAS Results:\n")
      cat("  Rows:", nrow(values$gwas_results), "\n")
      cat("  Columns:", paste(colnames(values$gwas_results), collapse = ", "), "\n")
      cat("  SNP_ID example:", head(values$gwas_results$SNP_ID, 3), "\n\n")
      
      # Check genotype matrix
      cat("Genotype Matrix:\n")
      if (!is.null(values$matched_data$geno)) {
        cat("  Dimensions:", dim(values$matched_data$geno), "\n")
        cat("  Row names (first 3):", head(rownames(values$matched_data$geno), 3), "\n")
        cat("  Column names (first 3):", head(colnames(values$matched_data$geno), 3), "\n\n")
      } else {
        cat("  Not available\n\n")
      }
      
      # Check phenotype data
      cat("Phenotype Data:\n")
      cat("  Columns:", paste(colnames(values$pheno_data), collapse = ", "), "\n")
      cat("  Trait values for", input$qtl_trait, ":", 
          head(values$pheno_data[[input$qtl_trait]], 3), "\n\n")
      
      # Check SNP info
      cat("SNP Info:\n")
      if (!is.null(values$processed_geno$snp_info)) {
        cat("  Rows:", nrow(values$processed_geno$snp_info), "\n")
        cat("  Columns:", paste(colnames(values$processed_geno$snp_info), collapse = ", "), "\n\n")
      } else {
        cat("  Not available\n\n")
      }
      
      # Check overlap
      if (!is.null(values$matched_data$geno)) {
        common_snps <- intersect(values$gwas_results$SNP_ID, colnames(values$matched_data$geno))
        cat("Common SNPs between GWAS and genotype matrix:", length(common_snps), "\n")
        
        common_samples <- intersect(values$pheno_data$Genotype, rownames(values$matched_data$geno))
        cat("Common samples between phenotype and genotype:", length(common_samples), "\n")
      }
    })
  })
  
  # Update QTL controls with debugging option
  output$qtl_controls <- renderUI({
    tagList(
      selectInput("qtl_trait", "Select Trait for QTL Analysis",
                  choices = if (!is.null(values$pheno_data)) 
                    setdiff(colnames(values$pheno_data), "Genotype") else NULL,
                  selected = NULL,
                  width = "100%"),
      
      selectInput("qtl_method", "QTL Method",
                  choices = c("GWAS-based" = "gwas", "Interval Mapping" = "em"),
                  selected = "gwas"),
      
      conditionalPanel(
        condition = "input.qtl_method == 'gwas'",
        numericInput("qtl_lod_threshold", "LOD Threshold",
                     value = 3.0, min = 1.0, max = 10.0, step = 0.5)
      ),
      
      conditionalPanel(
        condition = "input.qtl_method == 'em'",
        numericInput("qtl_p_threshold", "SNP P-value Threshold",
                     value = 0.05, min = 0.0001, max = 1, step = 0.01),
        numericInput("qtl_n_perm", "Number of Permutations",
                     value = 100, min = 0, max = 10000, step = 100)
      ),
      
      br(),
      
      fluidRow(
        column(6,
               actionButton("run_qtl_analysis", "Run Full QTL", 
                            class = "btn-primary",
                            icon = icon("chart-line"))
        ),
        column(6,
               actionButton("run_simple_qtl", "Run Simple QTL", 
                            class = "btn-info",
                            icon = icon("bolt"))
        )
      ),
      
      br(),
      
      actionButton("debug_qtl_data", "Debug Data", 
                   class = "btn-warning btn-sm",
                   icon = icon("bug")),
      actionButton("clear_qtl_results", "Clear Results", 
                   class = "btn-warning btn-sm",
                   icon = icon("broom"))
    )
  })
  
  # Enhanced QTL plot with better visualization
  output$qtl_plot <- renderPlot({
    req(values$qtl_results)
    
    result <- values$qtl_results
    
    if (is.null(result$scan_results)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No scan results available", size = 6) +
               theme_void())
    }
    
    # Create enhanced LOD profile plot
    plot_data <- result$scan_results %>%
      mutate(
        chr = factor(chr, levels = unique(chr)),
        significant = ifelse(!is.na(lod) & lod >= result$threshold, "Yes", "No")
      )
    
    # Calculate chromosome boundaries for plotting
    chrom_data <- plot_data %>%
      group_by(chr) %>%
      summarize(
        start = min(pos, na.rm = TRUE),
        end = max(pos, na.rm = TRUE),
        mid = mean(c(start, end), na.rm = TRUE)
      ) %>%
      ungroup()
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = pos, y = lod, color = significant)) +
      geom_point(alpha = 0.6, size = 1) +
      geom_line(aes(group = chr), alpha = 0.5, size = 0.5) +
      geom_hline(yintercept = result$threshold, 
                 linetype = "dashed", color = "red", alpha = 0.7, size = 1) +
      facet_grid(~chr, scales = "free_x", space = "free_x") +
      labs(
        title = paste("QTL LOD Profile -", input$qtl_trait),
        subtitle = paste("Method:", result$method, "| Threshold: LOD =", round(result$threshold, 2)),
        x = "Position",
        y = "LOD Score",
        color = "Significant"
      ) +
      scale_color_manual(values = c("Yes" = "red", "No" = "steelblue")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        panel.spacing = unit(0.1, "lines"),
        panel.grid.minor = element_blank()
      )
    
    # Add significant peaks if available
    if (!is.null(result$significant_qtl) && nrow(result$significant_qtl) > 0) {
      peaks <- result$significant_qtl %>%
        mutate(chr = factor(chr, levels = levels(plot_data$chr)))
      
      p <- p + 
        geom_point(data = peaks, aes(x = pos, y = lod), 
                   color = "red", size = 3, shape = 17) +
        geom_text(data = peaks, 
                  aes(x = pos, y = lod, label = paste("QTL", row_number())),
                  vjust = -1, size = 3, color = "red", fontface = "bold")
    }
    
    return(p)
  })
  
  # Enhanced QTL summary
  output$qtl_results_summary <- renderPrint({
    req(values$qtl_results)
    
    result <- values$qtl_results
    
    cat("=== QTL ANALYSIS RESULTS ===\n\n")
    cat("Method:", result$method, "\n")
    cat("Analysis date:", format(result$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Trait:", input$qtl_trait, "\n")
    
    if (!is.null(result$scan_results)) {
      cat("\n=== SCAN RESULTS ===\n")
      cat("Total positions scanned:", nrow(result$scan_results), "\n")
      cat("Chromosomes analyzed:", length(unique(result$scan_results$chr)), "\n")
      cat("Maximum LOD score:", round(max(result$scan_results$lod, na.rm = TRUE), 3), "\n")
      cat("Mean LOD score:", round(mean(result$scan_results$lod, na.rm = TRUE), 3), "\n")
      cat("Significance threshold (LOD):", result$threshold, "\n")
      
      if (!is.null(result$significant_qtl)) {
        cat("\n=== SIGNIFICANT QTLs ===\n")
        cat("Number of significant QTLs:", nrow(result$significant_qtl), "\n\n")
        
        for (i in 1:nrow(result$significant_qtl)) {
          qtl <- result$significant_qtl[i, ]
          cat(paste0("QTL ", i, ":\n"))
          cat("  Chromosome:", qtl$chr, "\n")
          cat("  Position:", round(qtl$pos, 2), "\n")
          cat("  LOD score:", round(qtl$lod, 3), "\n")
          if ("ci_low" %in% colnames(qtl)) {
            cat("  Confidence interval: [", round(qtl$ci_low, 2), ", ", 
                round(qtl$ci_high, 2), "]\n", sep = "")
          }
          if ("effect" %in% colnames(qtl)) {
            cat("  Effect size:", round(qtl$effect, 4), "\n")
          }
          if ("pve" %in% colnames(qtl)) {
            cat("  PVE (%):", round(qtl$pve, 2), "\n")
          }
          cat("\n")
        }
      } else {
        cat("\nNo significant QTLs found at LOD threshold", result$threshold, "\n")
      }
    }
    
    if (!is.null(result$perm_results)) {
      cat("\n=== PERMUTATION TEST ===\n")
      cat("Number of permutations:", nrow(result$perm_results), "\n")
      cat("Empirical p-value threshold:", round(result$threshold, 3), "\n")
    }
  })
  
  # Interactive QTL table
  output$qtl_results_table <- renderDT({
    req(values$qtl_results)
    
    result <- values$qtl_results
    
    if (is.null(result$scan_results)) {
      return(datatable(
        data.frame(Message = "No QTL results available"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    }
    
    # Format the scan results
    display_df <- result$scan_results %>%
      mutate(
        lod = round(lod, 3),
        p_value = ifelse("p_value" %in% colnames(.), 
                         format(p_value, scientific = TRUE, digits = 3), 
                         NA),
        effect = ifelse("effect" %in% colnames(.), round(effect, 4), NA),
        significant = ifelse(lod >= result$threshold, "Yes", "No")
      )
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf'),
        order = list(list(3, 'desc'))  # Sort by LOD descending
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "QTL Scan Results (All Markers)",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'significant',
        target = 'row',
        backgroundColor = styleEqual(c("Yes", "No"), c('#ffe6e6', 'white'))
      ) %>%
      formatStyle(
        'lod',
        background = styleColorBar(range(display_df$lod, na.rm = TRUE), 'lightgreen'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  # Download handler for QTL results
  output$download_qtl_results <- downloadHandler(
    filename = function() {
      paste0("qtl_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(values$qtl_results)
      
      # Create temporary directory
      temp_dir <- tempfile("qtl_analysis_")
      dir.create(temp_dir)
      
      # Save all results
      if (!is.null(values$qtl_results$scan_results)) {
        write.csv(values$qtl_results$scan_results,
                  file.path(temp_dir, "qtl_scan_results.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$qtl_results$significant_qtl)) {
        write.csv(values$qtl_results$significant_qtl,
                  file.path(temp_dir, "significant_qtl.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$qtl_results$perm_results)) {
        write.csv(values$qtl_results$perm_results,
                  file.path(temp_dir, "permutation_results.csv"),
                  row.names = FALSE)
      }
      
      # Create summary report
      summary_report <- paste(
        "QTL ANALYSIS REPORT",
        "===================",
        "",
        paste("Analysis Date:", Sys.time()),
        paste("Trait:", input$qtl_trait),
        paste("Method:", values$qtl_results$method),
        paste("Significance Threshold (LOD):", values$qtl_results$threshold),
        "",
        "SUMMARY STATISTICS",
        paste("Total markers scanned:", nrow(values$qtl_results$scan_results)),
        paste("Maximum LOD score:", round(max(values$qtl_results$scan_results$lod, na.rm = TRUE), 3)),
        paste("Significant QTLs:", ifelse(!is.null(values$qtl_results$significant_qtl), 
                                          nrow(values$qtl_results$significant_qtl), 0)),
        "",
        "FILES INCLUDED",
        "1. qtl_scan_results.csv - Full scan results with LOD scores",
        "2. significant_qtl.csv - Significant QTL peaks",
        "3. permutation_results.csv - Permutation test results (if available)",
        "",
        sep = "\n"
      )
      
      writeLines(summary_report, file.path(temp_dir, "README.txt"))
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "results.zip")
      zip::zip(zip_file, files = list.files(temp_dir, full.names = TRUE), 
               mode = "cherry-pick")
      
      # Copy to download location
      file.copy(zip_file, file)
      
      # Cleanup
      unlink(temp_dir, recursive = TRUE)
    }
  )
  
  # Debug output
  output$qtl_debug_info <- renderPrint({
    cat("Debug information will appear here after clicking 'Debug Data' button.\n")
  })
  
  #============================================================================================================
  #' Save genotype data in QTL ICI format similar to the uploaded example
  #' 
  #' @param filtered_gl A genlight object (filtered genotype data)
  #' @param file_path Path to save the Excel file
  #' @param mapping_function Mapping function code (1=Kosambi, 2=Haldane, 3=Morgan)
  #' @param marker_space_type Marker space type (1=intervals, 2=positions)
  #' @param mapping_pop_type Mapping population type
  #' 
  #' @return Returns the file path of the saved Excel file
  save_qtl_ici_format <- function(filtered_gl, file_path = "QTL_ICI_DATA.xlsx", 
                                  mapping_function = 1, marker_space_type = 2,
                                  mapping_pop_type = "F2") {
    
    tryCatch({
      # Create a new workbook
      wb <- createWorkbook()
      
      # ============================================
      # 1. Create GeneralInfo sheet
      # ============================================
      addWorksheet(wb, "GeneralInfo")
      
      # Prepare GeneralInfo data
      general_info <- data.frame(
        Row = c(7, 1, 2, nLoc(filtered_gl), nInd(filtered_gl)),
        Description = c(
          "Mapping Population Type",
          "Mapping Function (1 for Kosambi; 2 for Haldane; 3 for Morgan)",
          "Marker Space Type (1 for intervals; 2 for positions)",
          "Number of Markers",
          "Size of the mapping population"
        ),
        stringsAsFactors = FALSE
      )
      
      # Write GeneralInfo data
      writeData(wb, sheet = "GeneralInfo", 
                x = data.frame(A = general_info$Row, 
                               B = general_info$Description),
                startRow = 1, colNames = FALSE)
      
      # ============================================
      # 2. Create Anchor sheet (SNP and Chromosome)
      # ============================================
      addWorksheet(wb, "Anchor")
      
      # Get chromosome information
      # Try to extract from genlight object's chromosome slot
      chromosome_info <- NULL
      
      # Method 1: Check if chromosome information exists in genlight object
      if (!is.null(filtered_gl@chromosome)) {
        chromosome_info <- filtered_gl@chromosome
      } 
      # Method 2: Check if we have processed_geno$snp_info available (from reactive values)
      else if (exists("values$processed_geno$snp_info") && 
               !is.null(values$processed_geno$snp_info)) {
        snp_info <- values$processed_geno$snp_info
        chromosome_info <- snp_info$Chromosome[match(locNames(filtered_gl), snp_info$SNP_ID)]
      }
      # Method 3: Try to extract from SNP names (common patterns)
      else {
        # Try to extract chromosome and SNP names (e.g., "chr1_12345" or "1_12345")
        anchor_data <- data.frame(values$processed_geno$snp_info$SNP_ID, 
                                  values$processed_geno$snp_info$Chromosome, 
                                  stringsAsFactors = FALSE)
      }
      
      # Write Anchor data
      writeData(wb, sheet = "Anchor", x = anchor_data, startRow = 1, colNames = TRUE)
      
      # Format Anchor sheet
      addStyle(wb, sheet = "Anchor",
               style = createStyle(textDecoration = "bold"),
               rows = 1, cols = 1:2, gridExpand = TRUE)
      setColWidths(wb, sheet = "Anchor", cols = 1:2, widths = c(15, 12))
      freezePane(wb, sheet = "Anchor", firstRow = TRUE)
      
      # ============================================
      # 3. Create Genotype sheet
      # ============================================
      addWorksheet(wb, "Genotype")
      
      # Convert genlight to matrix
      geno_matrix <- as.matrix(filtered_gl)
      
      # Get SNP names and sample names
      snp_names <- locNames(filtered_gl)
      sample_names <- indNames(filtered_gl)
      
      # Check dimensions
      cat("Genotype matrix dimensions:", dim(geno_matrix), "\n")
      cat("Number of SNPs:", length(snp_names), "\n")
      cat("Number of samples:", length(sample_names), "\n")
      
      # Transpose to have SNPs as rows and samples as columns
      geno_matrix_t <- t(geno_matrix)
      
      # Create data frame with SNP IDs as first column
      if (nrow(geno_matrix_t) != length(snp_names)) {
        stop(paste("Mismatch: geno_matrix_t has", nrow(geno_matrix_t), 
                   "rows but there are", length(snp_names), "SNPs"))
      }
      
      if (ncol(geno_matrix_t) != length(sample_names)) {
        stop(paste("Mismatch: geno_matrix_t has", ncol(geno_matrix_t), 
                   "columns but there are", length(sample_names), "samples"))
      }
      
      # Create the data frame
      geno_df <- as.data.frame(geno_matrix_t)
      colnames(geno_df) <- sample_names
      
      # Add SNP_ID as the first column
      geno_df <- cbind(SNP_ID = snp_names, geno_df)
      
      # Write genotype data
      writeData(wb, sheet = "Genotype", x = geno_df, startRow = 1, colNames = TRUE)
      
      # Format column A as SNP IDs
      addStyle(wb, sheet = "Genotype",
               style = createStyle(textDecoration = "bold"),
               rows = 1:nrow(geno_df) + 1,
               cols = 1)
      
      # ============================================
      # 4. Add formatting
      # ============================================
      
      # Set column widths
      setColWidths(wb, sheet = "GeneralInfo", cols = 1:2, widths = c(10, 50))
      setColWidths(wb, sheet = "Genotype", cols = 1, widths = 15)
      
      # Set reasonable width for sample columns
      sample_cols <- 2:ncol(geno_df)
      setColWidths(wb, sheet = "Genotype", 
                   cols = sample_cols, 
                   widths = rep(8, length(sample_cols)))
      
      # Freeze headers in Genotype sheet
      freezePane(wb, sheet = "Genotype", firstRow = TRUE, firstCol = TRUE)
      
      # ============================================
      # 5. Save workbook
      # ============================================
      saveWorkbook(wb, file_path, overwrite = TRUE)
      
      cat(sprintf("Excel file saved successfully to: %s\n", file_path))
      cat(sprintf("GeneralInfo sheet: %d rows\n", nrow(general_info)))
      cat(sprintf("Anchor sheet: %d SNPs with chromosome info\n", nrow(anchor_data)))
      cat(sprintf("Genotype sheet: %d SNPs x %d samples\n", 
                  nLoc(filtered_gl), nInd(filtered_gl)))
      
      return(file_path)
      
    }, error = function(e) {
      cat(sprintf("Error saving QTL ICI format: %s\n", e$message))
      cat("Traceback:\n")
      print(traceback())
      return(NULL)
    })
  }
  
  
  
  #' Simplified version for Shiny app
  #' 
  #' @param values Reactive values containing filtered_gl
  #' @param output_dir Directory to save the file
  save_qtl_data <- function(values, output_dir = getwd()) {
    req(values$filtered_gl)
    
    # Create a unique filename
    file_path <- file.path(output_dir, 
                           paste0("QTL_ICI_", 
                                  format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                  ".xlsx"))
    
    # Add debug output
    cat("Starting QTL ICI export...\n")
    cat("Filtered GL dimensions:", nInd(values$filtered_gl), "samples,",
        nLoc(values$filtered_gl), "SNPs\n")
    
    result <- save_qtl_ici_format(
      filtered_gl = values$filtered_gl,
      file_path = file_path,
      mapping_function = 1,
      marker_space_type = 2,
      mapping_pop_type = "F2"
    )
    
    if (is.null(result)) {
      cat("QTL ICI export failed\n")
    } else {
      cat("QTL ICI export successful:", result, "\n")
    }
    
    return(result)
  }
  
  
  
  # ==========================================================================
  # QTL ICI EXPORT FUNCTIONALITY
  # ==========================================================================
  
  # Export button observer
  observeEvent(input$export_qtl_ici, {
    req(values$filtered_gl)
    
    # Show progress modal
    showModal(modalDialog(
      title = "Exporting QTL ICI Format",
      tagList(
        tags$div(
          style = "text-align: center;",
          icon("spinner", class = "fa-spin fa-2x"),
          tags$h4("Creating Excel file in QTL ICI format..."),
          tags$div(
            id = "qtl_export_details",
            style = "margin-top: 10px; font-size: 12px; color: #666;"
          )
        )
      ),
      footer = NULL,
      easyClose = FALSE,
      size = "s"
    ))
    
    # In the Export button observer, replace the tryCatch block:
    tryCatch({
      shinyjs::html("qtl_export_details", "Step 1/3: Preparing data...")
      Sys.sleep(0.5)
      
      # Debug: check filtered_gl
      cat("Export QTL ICI - Checking filtered_gl:\n")
      cat("  Number of individuals:", nInd(values$filtered_gl), "\n")
      cat("  Number of loci:", nLoc(values$filtered_gl), "\n")
      cat("  Sample names (first 3):", head(indNames(values$filtered_gl), 3), "\n")
      cat("  SNP names (first 3):", head(locNames(values$filtered_gl), 3), "\n")
      
      # Create the Excel file
      file_path <- save_qtl_data(values)
      
      if (is.null(file_path)) {
        removeModal()
        showNotification("Failed to create QTL ICI file", type = "error")
        return()
      }
      
      shinyjs::html("qtl_export_details", "Step 2/3: Creating worksheets...")
      Sys.sleep(0.5)
      
      # Update download handler
      output$download_qtl_ici <- downloadHandler(
        filename = function() {
          basename(file_path)
        },
        content = function(file) {
          file.copy(file_path, file)
        }
      )
      
      shinyjs::html("qtl_export_details", "Step 3/3: Finalizing...")
      Sys.sleep(0.5)
      
      removeModal()
      
      # Show download button
      output$qtl_export_status <- renderUI({
        div(
          class = "alert alert-success",
          icon("check-circle"),
          strong("✓ QTL ICI file created successfully!"),
          br(),
          paste("File:", basename(file_path)),
          br(),
          paste("Size:", format(file.size(file_path), big.mark = ","), "bytes"),
          br(),
          paste("SNPs:", nLoc(values$filtered_gl)),
          br(),
          paste("Samples:", nInd(values$filtered_gl)),
          br(), br(),
          downloadButton("download_qtl_ici", "Download File", 
                         class = "btn-success btn-block")
        )
      })
      
      showNotification(
        HTML(paste(
          "✅ QTL ICI file created successfully!<br>",
          "Contains:", nLoc(values$filtered_gl), "SNPs and", 
          nInd(values$filtered_gl), "samples"
        )),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      removeModal()
      cat("Error in QTL ICI export observer:", e$message, "\n")
      cat("Traceback:\n")
      print(traceback())
      showNotification(
        paste("Error creating QTL ICI file:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  # Status UI for export
  output$qtl_export_status <- renderUI({
    if (!is.null(values$filtered_gl)) {
      div(
        class = "alert alert-info",
        icon("info-circle"),
        strong("Ready to export QTL ICI format"),
        br(),
        paste("Available data:", nLoc(values$filtered_gl), "SNPs,",
              nInd(values$filtered_gl), "samples"),
        br(), br(),
        actionButton("export_qtl_ici", "Export QTL ICI Format",
                     icon = icon("file-excel"),
                     class = "btn-primary btn-block")
      )
    } else {
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        strong("No filtered genotype data available"),
        br(),
        "Please run Quality Control first to create filtered_gl"
      )
    }
  })
  
  # Preview of data to be exported
  output$qtl_preview <- renderPrint({
    req(values$filtered_gl)
    
    cat("=== QTL ICI EXPORT PREVIEW ===\n\n")
    cat("File will contain:\n")
    cat("  - GeneralInfo worksheet with metadata\n")
    cat("  - Genotype worksheet with SNP data\n\n")
    cat("  - Anchor worksheet with SNP-Chromosome map data\n\n")
    
    cat("Data summary:\n")
    cat("  SNPs:", nLoc(values$filtered_gl), "\n")
    cat("  Samples:", nInd(values$filtered_gl), "\n")
    cat("  Ploidy:", unique(values$filtered_gl@ploidy), "\n\n")
    
    cat("First 5 SNP IDs:\n")
    cat(paste(head(locNames(values$filtered_gl), 5), collapse = "\n"), "\n\n")
    
    cat("First 5 sample IDs:\n")
    cat(paste(head(indNames(values$filtered_gl), 5), collapse = "\n"), "\n")
  })
  
  #====================================================================================
  #UNIFY SNPs FOR ALL transcripts
  #====================================================================================
  
  #' Create unified GWAS tidy format
  create_unified_gwas_tidy <- function(trait_gwas_results) {
    require(dplyr)
    
    if (is.null(trait_gwas_results)) return(NULL)
    
    unified_df <- lapply(names(trait_gwas_results), function(trait) {
      df <- trait_gwas_results[[trait]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      
      # Select key columns and add trait info
      df %>%
        mutate(Trait = trait) %>%
        select(Trait, SNP_ID, Chromosome, Position, P_value, 
               P_adjusted, Effect, SE, R_squared, Significant_FDR)
    }) %>%
      bind_rows() %>%
      arrange(Trait, Chromosome, Position)
    
    return(unified_df)
  }
  # Usage:
  #values$multi_trait_unified <- create_unified_gwas_tidy(values$multi_trait_results$trait_gwas_results)
  
  
  
  create_wide_gwas_format <- function(trait_gwas_results) {
    require(dplyr)
    require(tidyr)
    
    if (is.null(trait_gwas_results)) return(NULL)
    
    # First create long format
    long_df <- create_unified_gwas_tidy(trait_gwas_results)
    
    if (is.null(long_df)) return(NULL)
    
    # Create wide format for P-values
    pvalue_wide <- long_df %>%
      select(SNP_ID, Chromosome, Position, Trait, P_value) %>%
      pivot_wider(
        id_cols = c(SNP_ID, Chromosome, Position),
        names_from = Trait,
        values_from = P_value,
        names_prefix = "P_"
      )
    
    # Create wide format for Effects
    effect_wide <- long_df %>%
      select(SNP_ID, Trait, Effect) %>%
      pivot_wider(
        id_cols = SNP_ID,
        names_from = Trait,
        values_from = Effect,
        names_prefix = "Effect_"
      )
    
    # Combine wide formats
    wide_df <- pvalue_wide %>%
      left_join(effect_wide, by = "SNP_ID")
    
    return(wide_df)
  }
  
  # Usage:
  #values$multi_trait_wide <- create_wide_gwas_format(values$multi_trait_results$trait_gwas_results)
  
  
  
  #' Create enhanced unified dataframe
  create_enhanced_unified_df <- function(trait_gwas_results, trait_correlations = NULL) {
    require(dplyr)
    
    if (is.null(trait_gwas_results)) return(NULL)
    
    # Create basic unified dataframe
    unified <- create_unified_gwas_tidy(trait_gwas_results)
    
    if (is.null(unified)) return(NULL)
    
    # Add significance flags
    unified <- unified %>%
      mutate(
        Significant_0_05 = P_value < 0.05,
        Significant_0_01 = P_value < 0.01,
        Significant_0_001 = P_value < 0.001,
        Significant_Bonferroni = P_adjusted < 0.05,
        LOD = -log10(P_value)
      )
    
    # Add trait correlation information if available
    if (!is.null(trait_correlations)) {
      # Calculate which traits each SNP affects
      trait_affected <- unified %>%
        filter(Significant_0_05) %>%
        group_by(SNP_ID) %>%
        summarise(
          Traits_Affected = paste(sort(unique(Trait)), collapse = ";"),
          N_Traits_Affected = n_distinct(Trait),
          .groups = "drop"
        )
      
      unified <- unified %>%
        left_join(trait_affected, by = "SNP_ID")
    }
    
    # Add SNP rank within each trait
    unified <- unified %>%
      group_by(Trait) %>%
      mutate(
        Rank = rank(P_value),
        Percentile = rank(P_value) / n() * 100
      ) %>%
      ungroup()
    
    return(unified)
  }
  
  
  
  #' Create complete integrated dataframe
  create_complete_integrated_df <- function(multi_trait_results) {
    require(dplyr)
    
    # Check if we have all necessary components
    if (is.null(multi_trait_results) || 
        is.null(multi_trait_results$trait_gwas_results)) {
      return(NULL)
    }
    
    # Create unified trait results
    trait_unified <- create_enhanced_unified_df(
      multi_trait_results$trait_gwas_results,
      multi_trait_results$trait_correlations
    )
    
    if (is.null(trait_unified)) return(NULL)
    
    # If we have combined results, merge them
    if (!is.null(multi_trait_results$combined_results)) {
      # Aggregate trait-specific information
      trait_summary <- trait_unified %>%
        group_by(SNP_ID) %>%
        summarise(
          All_Traits = paste(sort(unique(Trait)), collapse = ";"),
          Significant_Traits = paste(sort(unique(Trait[Significant_0_05])), collapse = ";"),
          Min_P_Value_Trait = Trait[which.min(P_value)][1],
          Min_P_Value = min(P_value, na.rm = TRUE),
          Max_Effect_Trait = Trait[which.max(abs(Effect))][1],
          Max_Effect = Effect[which.max(abs(Effect))][1],
          Mean_Effect = mean(Effect, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Merge with combined results
      complete_df <- multi_trait_results$combined_results %>%
        left_join(trait_summary, by = "SNP_ID")
      
    } else {
      complete_df <- trait_unified
    }
    
    return(complete_df)
  }
  
  # Usage:
  #values$complete_integrated_df <- create_complete_integrated_df(values$multi_trait_results)
  
  
  #==========================================================================================
  # Plot from combined dataframes
  #==========================================================================================
  
  # Observer to create unified dataframe when multi-trait results are available
  observeEvent(values$multi_trait_results, {
    req(values$multi_trait_results)
    
    withProgress(message = 'Creating unified multi-trait data...', value = 0, {
      
      tryCatch({
        incProgress(0.1, detail = "Processing trait results...")
        
        # Method 2: Tidy unified dataframe
        trait_list <- values$multi_trait_results$trait_gwas_results
        
        if (is.null(trait_list) || length(trait_list) == 0) {
          showNotification("No trait GWAS results available", type = "warning")
          return()
        }
        
        # Create unified dataframe
        unified_list <- list()
        for (trait_name in names(trait_list)) {
          df <- trait_list[[trait_name]]
          if (!is.null(df) && nrow(df) > 0) {
            # Add trait identifier and ensure required columns
            df <- df %>%
              mutate(Trait = trait_name)
            
            # Add LOD score
            if ("P_value" %in% colnames(df)) {
              df$LOD <- -log10(df$P_value)
            }
            
            unified_list[[trait_name]] <- df
          }
        }
        
        incProgress(0.2, detail = "Combining data...")
        
        if (length(unified_list) == 0) {
          showNotification("No valid trait data to unify", type = "warning")
          return()
        }
        
        # Combine all dataframes
        values$multi_trait_unified_df <- do.call(rbind, unified_list)
        rownames(values$multi_trait_unified_df) <- NULL
        
        incProgress(0.3, detail = "Initializing dataframes format...")
        
        # Also create wide format for comparison
        tryCatch({
          incProgress(0.4, detail = "Creating unified tidy format...")
          # Create unified dataframe
          values$multi_trait_unified_tidy <- create_unified_gwas_tidy(trait_list)
          
          incProgress(0.5, detail = "Creating unified wide format...")
          # Optional: Create wide format for comparison
          values$multi_trait_unified_wide <- create_wide_gwas_format(trait_list)
          
          incProgress(0.6, detail = "Creating unified enhanced format...")
          values$multi_trait_enhanced_unified_df <- create_enhanced_unified_df(trait_list,values$multi_trait_results$trait_correlations)
          
          incProgress(0.7, detail = "Creating complete integratedformat...")
          values$multi_trait_complete_integrated_df <- create_complete_integrated_df(values$multi_trait_results)
          
        }, error = function(e) {
          cat("Error creating some dataframes format:", e$message, "\n")
          return(NULL)
        })
        
        incProgress(0.8, detail = "Calculating statistics...")
        
        # Calculate summary statistics
        if (!is.null(values$multi_trait_unified_df)) {
          # Add significance flags
          values$multi_trait_unified_df <- values$multi_trait_unified_df %>%
            mutate(
              Significant_0_05 = ifelse("P_value" %in% colnames(.), P_value < 0.05, NA),
              Significant_0_01 = ifelse("P_value" %in% colnames(.), P_value < 0.01, NA),
              Significant_0_001 = ifelse("P_value" %in% colnames(.), P_value < 0.001, NA),
              Significant_Bonferroni = ifelse("P_adjusted" %in% colnames(.), 
                                              P_adjusted < 0.05, NA)
            )
        }
        
        incProgress(1.0, detail = "Complete!")
        
        showNotification(
          HTML(paste(
            "✅ Unified multi-trait data created!<br>",
            "Traits:", length(unified_list), "<br>",
            "Total rows:", nrow(values$multi_trait_unified_df)
          )),
          type = "message",
          duration = 3
        )
        
      }, error = function(e) {
        showNotification(paste("Error creating unified data:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  
  # 1. Separate Manhattan Plots per Trait
  # Faceted Manhattan plots
  plot_manhattan_faceted <- function(df) {
    ggplot(df, aes(x = Position, y = -log10(P_value), color = as.factor(Chromosome))) +
      geom_point(alpha = 0.6, size = 1) +
      geom_hline(yintercept = -log10(0.05/nrow(df)), 
                 linetype = "dashed", color = "red", alpha = 0.5) +
      geom_hline(yintercept = -log10(5e-8), 
                 linetype = "dashed", color = "blue", alpha = 0.5) +
      facet_wrap(~Trait, scales = "free_y", ncol = 2) +
      labs(title = "Manhattan Plots by Trait",
           x = "Genomic Position",
           y = "-log10(P-value)",
           color = "Chromosome") +
      theme_minimal() +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"))
  }
  
  
  # Combined Manhattan Plot (All Traits)
  # All traits in one plot
  plot_manhattan_combined <- function(df) {
    ggplot(df, aes(x = Position, y = -log10(P_value), color = Trait)) +
      geom_point(alpha = 0.4, size = 1) +
      facet_grid(~Chromosome, scales = "free_x", space = "free_x") +
      labs(title = "Combined Manhattan Plot - All Traits",
           x = "Position",
           y = "-log10(P-value)") +
      theme_minimal()
  }
  
  
  
  # QQ Plot per Trait
  plot_qq_by_trait <- function(df) {
    # Calculate expected vs observed
    qq_data <- df %>%
      group_by(Trait) %>%
      arrange(P_value) %>%
      mutate(
        observed = -log10(P_value),
        expected = -log10(ppoints(n()))
      )
    
    ggplot(qq_data, aes(x = expected, y = observed, color = Trait)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      facet_wrap(~Trait) +
      labs(title = "QQ Plots by Trait",
           x = "Expected -log10(P)",
           y = "Observed -log10(P)") +
      theme_minimal()
  }
  
  
  
  # QQ Plot Overlay (All Traits)
  plot_qq_overlay <- function(df) {
    ggplot(df, aes(sample = -log10(P_value), color = Trait)) +
      stat_qq(alpha = 0.6) +
      stat_qq_line(color = "black", linetype = "dashed") +
      labs(title = "QQ Plot Overlay - All Traits",
           x = "Theoretical Quantiles",
           y = "Sample Quantiles (-log10(P))") +
      theme_minimal()
  }
  
  
  # EFFECT SIZE VISUALIZATIONS
  # A. Effect Size Distribution per Trait
  plot_effect_distribution <- function(df) {
    ggplot(df, aes(x = Effect, fill = Trait)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
      facet_wrap(~Trait, scales = "free") +
      labs(title = "Effect Size Distribution by Trait",
           x = "Effect Size",
           y = "Density") +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  
  
  # B. Effect Size vs P-value (Volcano Plot)
  plot_volcano_by_trait <- function(df) {
    ggplot(df, aes(x = Effect, y = -log10(P_value), color = Trait)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), 
                 linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
      facet_wrap(~Trait, scales = "free") +
      labs(title = "Volcano Plots - Effect vs Significance",
           x = "Effect Size",
           y = "-log10(P-value)") +
      theme_minimal()
  }
  
  
  
  # TRAIT COMPARISON PLOTS
  # A. Correlation of Effect Sizes Between Traits
  plot_effect_correlation <- function(df) {
    # Convert to wide format for correlation
    effect_wide <- df %>%
      select(SNP_ID, Trait, Effect) %>%
      pivot_wider(names_from = Trait, values_from = Effect)
    
    cor_matrix <- cor(effect_wide[, -1], use = "pairwise.complete.obs")
    
    # Create heatmap
    ggplot(melt(cor_matrix), aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, limit = c(-1, 1)) +
      labs(title = "Correlation of SNP Effects Across Traits",
           x = "", y = "", fill = "Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  
  
  # B. Effect Size Scatter Matrix
  plot_effect_scatter_matrix <- function(df) {
    # Wide format
    effect_wide <- df %>%
      select(SNP_ID, Trait, Effect) %>%
      pivot_wider(names_from = Trait, values_from = Effect)
    
    # Scatter plot matrix
    ggpairs(effect_wide[, -1],
            title = "Scatter Matrix of SNP Effects Across Traits",
            lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
            diag = list(continuous = wrap("densityDiag", alpha = 0.5)))
  }
  
  
  
  #  SIGNIFICANCE ANALYSIS PLOTS
  # A. Number of Significant SNPs per Trait
  plot_sig_snps_by_trait <- function(df) {
    sig_counts <- df %>%
      group_by(Trait) %>%
      summarise(
        n_total = n(),
        n_sig_05 = sum(P_value < 0.05),
        n_sig_01 = sum(P_value < 0.01),
        n_sig_001 = sum(P_value < 0.001)
      ) %>%
      pivot_longer(cols = starts_with("n_sig"), 
                   names_to = "Threshold", 
                   values_to = "Count")
    
    ggplot(sig_counts, aes(x = Trait, y = Count, fill = Threshold)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "Significant SNPs by Trait and Threshold",
           x = "Trait", y = "Number of SNPs") +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  
  # B. P-value Distribution by Trait
  plot_pvalue_distribution <- function(df) {
    ggplot(df, aes(x = P_value, fill = Trait)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      scale_x_log10() +
      facet_wrap(~Trait, scales = "free_y") +
      labs(title = "P-value Distribution by Trait",
           x = "P-value (log10 scale)",
           y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  
  # CHROMOSOME-LEVEL PLOTS
  # A. SNP Density by Chromosome and Trait
  plot_snp_density <- function(df) {
    ggplot(df, aes(x = Position, fill = Trait)) +
      geom_density(alpha = 0.5) +
      facet_grid(Trait ~ Chromosome, scales = "free_x") +
      labs(title = "SNP Density by Chromosome and Trait",
           x = "Genomic Position",
           y = "Density") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  
  # B. Chromosome-wide Effect Summary
  plot_chromosome_effects <- function(df) {
    chrom_summary <- df %>%
      group_by(Chromosome, Trait) %>%
      summarise(
        mean_effect = mean(abs(Effect), na.rm = TRUE),
        max_effect = max(abs(Effect), na.rm = TRUE),
        n_sig = sum(P_value < 0.05),
        .groups = "drop"
      )
    
    ggplot(chrom_summary, aes(x = Chromosome, y = mean_effect, fill = Trait)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "Average Absolute Effect Size by Chromosome",
           x = "Chromosome",
           y = "Mean |Effect|") +
      theme_minimal()
  }
  
  
  # TOP SNP VISUALIZATIONS
  # A. Top SNPs Across All Traits
  plot_top_snps_across_traits <- function(df, top_n = 20) {
    top_snps <- df %>%
      group_by(SNP_ID) %>%
      summarise(
        min_p = min(P_value, na.rm = TRUE),
        n_traits = n_distinct(Trait),
        traits = paste(sort(unique(Trait)), collapse = ", ")
      ) %>%
      arrange(min_p) %>%
      head(top_n)
    
    ggplot(top_snps, aes(x = reorder(SNP_ID, -log10(min_p)), 
                         y = -log10(min_p), fill = as.factor(n_traits))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = paste("Top", top_n, "SNPs Across All Traits"),
           x = "SNP ID",
           y = "-log10(Min P-value)",
           fill = "Number of Traits") +
      scale_fill_viridis_d() +
      theme_minimal()
  }
  
  
  
  # B. Top SNP Effects Heatmap
  plot_top_snps_heatmap <- function(df, top_n = 30) {
    # Get top SNPs by minimum p-value
    top_snp_ids <- df %>%
      group_by(SNP_ID) %>%
      summarise(min_p = min(P_value)) %>%
      arrange(min_p) %>%
      pull(SNP_ID) %>%
      head(top_n)
    
    # Filter and create matrix
    heatmap_data <- df %>%
      filter(SNP_ID %in% top_snp_ids) %>%
      select(SNP_ID, Trait, Effect) %>%
      pivot_wider(names_from = Trait, values_from = Effect) %>%
      column_to_rownames("SNP_ID")
    
    # Create heatmap
    heatmaply(heatmap_data,
              main = paste("Effect Sizes of Top", top_n, "SNPs"),
              xlab = "Trait", ylab = "SNP",
              colors = viridis(100),
              dendrogram = "both",
              scale = "column")
  }
  
  
  # INTERACTIVE PLOTS
  # A. Interactive Manhattan Plot
  plot_interactive_manhattan <- function(df) {
    plot_ly(df,
            x = ~Position,
            y = ~-log10(P_value),
            color = ~Trait,
            type = "scattergl",
            mode = "markers",
            marker = list(size = 5, opacity = 0.6),
            hoverinfo = "text",
            text = ~paste("SNP:", SNP_ID,
                          "<br>Trait:", Trait,
                          "<br>Chr:", Chromosome,
                          "<br>Pos:", Position,
                          "<br>P-value:", format(P_value, scientific = TRUE),
                          "<br>Effect:", round(Effect, 3))) %>%
      layout(title = "Interactive Manhattan Plot",
             xaxis = list(title = "Position"),
             yaxis = list(title = "-log10(P-value)"),
             hovermode = "closest")
  }
  
  
  # B. Interactive Volcano Plot
  plot_interactive_volcano <- function(df) {
    df %>%
      plot_ly(x = ~Effect,
              y = ~-log10(P_value),
              color = ~Trait,
              type = "scattergl",
              mode = "markers",
              marker = list(size = 8, opacity = 0.6),
              hoverinfo = "text",
              text = ~paste("SNP:", SNP_ID,
                            "<br>Trait:", Trait,
                            "<br>P-value:", format(P_value, scientific = TRUE),
                            "<br>Effect:", round(Effect, 4))) %>%
      layout(title = "Interactive Volcano Plot",
             xaxis = list(title = "Effect Size"),
             yaxis = list(title = "-log10(P-value)"))
  }
  
  
  #' Plot pleiotropy analysis results
  plot_pleiotropy_analysis <- function(unified_data, p_threshold = 0.05, top_n = 20) {
    tryCatch({
      # Ensure we're working with a data frame
      if (!is.data.frame(unified_data)) {
        unified_data <- as.data.frame(unified_data)
      }
      
      # Flatten any list columns
      unified_data <- as.data.frame(lapply(unified_data, function(x) {
        if (is.list(x)) {
          # Convert list to character or numeric
          if (all(sapply(x, is.character))) {
            return(unlist(x))
          } else if (all(sapply(x, is.numeric))) {
            return(unlist(x))
          } else {
            return(as.character(unlist(x)))
          }
        }
        return(x)
      }))
      
      # Extract unified trait-SNP associations
      p_cols <- grep("^P_|_P$", colnames(unified_data), value = TRUE)
      
      if (length(p_cols) == 0) {
        cat("No P-value columns found\n")
        return(NULL)
      }
      
      # Create a unified data frame of SNP-trait associations
      unified_trait_snps <- data.frame()
      
      for (p_col in p_cols) {
        trait <- gsub("^P_|_P$", "", p_col)
        
        trait_data <- unified_data %>%
          select(SNP = SNP_ID, Trait = !!sym(trait), P_value = !!sym(p_col)) %>%
          filter(!is.na(P_value)) %>%
          mutate(Trait = as.character(trait))
        
        if (nrow(trait_data) > 0) {
          unified_trait_snps <- bind_rows(unified_trait_snps, trait_data)
        }
      }
      
      if (nrow(unified_trait_snps) == 0) {
        cat("No trait-SNP associations found\n")
        return(NULL)
      }
      
      # Filter by p-value threshold
      unified_trait_snps <- unified_trait_snps %>%
        filter(P_value < p_threshold)
      
      if (nrow(unified_trait_snps) == 0) {
        cat("No significant associations at p <", p_threshold, "\n")
        return(NULL)
      }
      
      # Count number of traits per SNP and number of SNPs per trait
      pleiotropy_counts <- unified_trait_snps %>%
        group_by(SNP) %>%
        summarise(N_Traits = n()) %>%
        ungroup()
      
      trait_counts <- unified_trait_snps %>%
        group_by(Trait) %>%
        summarise(N_SNPs = n()) %>%
        ungroup()
      
      # Create the plot
      p <- ggplot() +
        # SNP pleiotropy (number of traits per SNP)
        geom_bar(data = pleiotropy_counts, 
                 aes(x = factor(SNP, levels = SNP[order(N_Traits, decreasing = TRUE)]), 
                     y = N_Traits),
                 stat = "identity", fill = "steelblue", alpha = 0.8) +
        # Trait polygenicity (number of SNPs per trait)
        geom_point(data = trait_counts,
                   aes(x = "Trait Summary", y = N_SNPs, color = Trait),
                   size = 4, position = position_jitter(width = 0.2)) +
        labs(title = "Pleiotropy Analysis",
             x = "SNP",
             y = "Count",
             color = "Traits") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        facet_grid(. ~ "Pleiotropy (SNPs) vs Polygenicity (Traits)", scales = "free_x")
      
      return(p)
      
    }, error = function(e) {
      cat("Error in plot_pleiotropy_analysis:", e$message, "\n")
      return(NULL)
    })
  }
  
  
  
  
  # B. Effect Direction Consistency
  plot_effect_direction <- function(df) {
    # Calculate effect direction consistency
    direction_data <- df %>%
      filter(P_value < 0.05) %>%
      group_by(SNP_ID) %>%
      summarise(
        n_traits = n(),
        pos_effects = sum(Effect > 0),
        neg_effects = sum(Effect < 0),
        direction_consistency = ifelse(pos_effects == 0 | neg_effects == 0,
                                       "Consistent", "Mixed")
      )
    
    ggplot(direction_data, aes(x = direction_consistency, fill = direction_consistency)) +
      geom_bar() +
      labs(title = "Effect Direction Consistency",
           subtitle = "For SNPs significant in multiple traits (p < 0.05)",
           x = "Effect Direction",
           y = "Number of SNPs") +
      scale_fill_manual(values = c("Consistent" = "green", "Mixed" = "orange")) +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  
  # COMPREHENSIVE DASHBOARD PLOTS
  # A. Multi-panel Summary Dashboard
  create_summary_dashboard <- function(df) {
    # Panel 1: Manhattan
    p1 <- ggplot(df, aes(x = Position, y = -log10(P_value), color = Trait)) +
      geom_point(alpha = 0.3, size = 0.5) +
      facet_grid(~Chromosome, scales = "free_x", space = "free_x") +
      theme(axis.text.x = element_blank(),
            legend.position = "none")
    
    # Panel 2: Effect distribution
    p2 <- ggplot(df, aes(x = Effect, fill = Trait)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~Trait, nrow = 1) +
      theme(legend.position = "none")
    
    # Panel 3: Significant SNPs
    p3 <- df %>%
      group_by(Trait) %>%
      summarise(n_sig = sum(P_value < 0.05)) %>%
      ggplot(aes(x = reorder(Trait, n_sig), y = n_sig)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip()
    
    # Combine
    (p1 / p2 / p3) + 
      plot_annotation(title = "Multi-Trait GWAS Summary Dashboard")
  }
  
  
  # In your server.R
  output$multi_trait_plots <- renderUI({
    req(values$multi_trait_unified_df)
    
    # Get plot type from input
    plot_type <- input$plot_type
    
    # Generate plot based on selection
    plot_output <- switch(plot_type,
                          "manhattan" = plot_manhattan_faceted(values$multi_trait_unified_df),
                          "volcano" = plot_volcano_by_trait(values$multi_trait_unified_df),
                          "effect_correlation" = plot_effect_correlation(values$multi_trait_unified_df),
                          "pleiotropy" = plot_pleiotropy_analysis(values$multi_trait_unified_df),
                          "top_snps" = plot_top_snps_across_traits(values$multi_trait_unified_df, 
                                                                   top_n = input$top_n),
                          # ... more plot types
                          ggplot() + annotate("text", x=0.5, y=0.5, label="Select a plot type")
    )
    
    renderPlot({ plot_output }, height = 600)
  })
  
  # UI selector
  plot_type_selector <- selectInput("plot_type", "Select Plot Type",
                                    choices = c(
                                      "Manhattan Plot" = "manhattan",
                                      "Volcano Plot" = "volcano",
                                      "Effect Correlation" = "effect_correlation",
                                      "QQ Plot" = "qq",
                                      "Pleiotropy Analysis" = "pleiotropy",
                                      "Top SNPs" = "top_snps",
                                      "Effect Distribution" = "effect_dist",
                                      "Significant SNPs" = "sig_snps"
                                    ))
  
  #===========================================================================
  #ADD LEFT-RIGHT MARKER
  #===========================================================================
  
  # Function to download the dataframe
  output$download_cross_snp <- downloadHandler(
    filename = function() {
      paste("cross_snp_analysis_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(values$cross_snp_df, file, row.names = FALSE)
    }
  )
  
  
  
  
  
  # ==========================================================================
  # OUTPUT RENDERERS FOR UNIFIED MULTI-TRAIT DATA
  # ==========================================================================
  
  # Unified data table
  output$unified_data_table <- renderDT({
    req(values$multi_trait_unified_df)
    
    tryCatch({
      # Format for display
      display_df <- values$multi_trait_unified_df %>%
        mutate(
          P_value = format(P_value, scientific = TRUE, digits = 3),
          Effect = round(Effect, 4)
        )
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Unified Multi-Trait GWAS Data (Long Format)"
      )
    }, error = function(e) {
      datatable(data.frame(Message = "Error loading unified data"))
    })
  })
  
  
  
  # Wide format table
  output$wide_format_table <- renderDT({
    req(values$multi_trait_wide_df)
    
    tryCatch({
      display_df <- values$multi_trait_wide_df
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Multi-Trait GWAS Results (Wide Format - SNPs as Rows)",
        extensions = 'Buttons',
        filter = 'top'
      )
      
    }, error = function(e) {
      cat("Error rendering wide table:", e$message, "\n")
      return(datatable(
        data.frame(Message = "Wide format data not available"),
        options = list(pageLength = 5),
        rownames = FALSE
      ))
    })
  })
  
  
  # Enhanced unified table
  output$enhanced_unified_table <- renderDT({
    req(values$multi_trait_enhanced_unified_df)
    
    tryCatch({
      display_df <- values$multi_trait_enhanced_unified_df %>%
        mutate(
          P_value = format(P_value, scientific = TRUE, digits = 3),
          Effect = round(Effect, 4),
          LOD = round(LOD, 2)
        )
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Enhanced Unified Multi-Trait GWAS Data"
      )
    }, error = function(e) {
      datatable(data.frame(Message = "Enhanced data not available"))
    })
  })
  
  
  
  # Complete integrated table
  output$complete_integrated_table <- renderDT({
    req(values$multi_trait_complete_integrated_df)
    
    tryCatch({
      display_df <- values$multi_trait_complete_integrated_df %>%
        mutate(
          Min_P_Value = format(Min_P_Value, scientific = TRUE, digits = 3),
          Mean_Effect = round(Mean_Effect, 4)
        )
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "Complete Integrated Multi-Trait GWAS Data"
      )
    }, error = function(e) {
      datatable(data.frame(Message = "Complete integrated data not available"))
    })
  })
  
  
  # Enhanced visualization plots
  output$enhanced_manhattan_plot <- renderPlot({
    req(values$multi_trait_enhanced_unified_df)
    
    plot_manhattan_faceted(values$multi_trait_enhanced_unified_df)
  })
  
  
  output$enhanced_volcano_plot <- renderPlot({
    req(values$multi_trait_enhanced_unified_df)
    
    plot_volcano_by_trait(values$multi_trait_enhanced_unified_df)
  })
  
  output$enhanced_effect_plot <- renderPlot({
    req(values$multi_trait_enhanced_unified_df)
    
    plot_effect_distribution(values$multi_trait_enhanced_unified_df)
  })
  
  
  
  
  
  
  
  # Plot renderer
  output$unified_plot_output <- renderPlot({
    req(values$multi_trait_unified_df, input$plot_type)
    
    tryCatch({
      # Get selected plot type
      plot_type <- input$plot_type
      
      # Generate appropriate plot
      plot_obj <- switch(plot_type,
                         "manhattan" = plot_manhattan_faceted(values$multi_trait_unified_df),
                         "volcano" = plot_volcano_by_trait(values$multi_trait_unified_df),
                         "qq" = plot_qq_by_trait(values$multi_trait_unified_df),
                         "effect_dist" = plot_effect_distribution(values$multi_trait_unified_df),
                         "effect_cor" = plot_effect_correlation(values$multi_trait_unified_df),
                         "sig_snps" = plot_sig_snps_by_trait(values$multi_trait_unified_df),
                         "pleiotropy" = plot_pleiotropy_analysis(values$multi_trait_unified_df),
                         "top_snps" = plot_top_snps_across_traits(values$multi_trait_unified_df, 
                                                                  top_n = input$top_n_snps),
                         "pvalue_dist" = plot_pvalue_distribution(values$multi_trait_unified_df),
                         "effect_dir" = plot_effect_direction(values$multi_trait_unified_df),
                         # Default plot
                         ggplot() + 
                           annotate("text", x = 0.5, y = 0.5, 
                                    label = "Select a plot type from the dropdown", size = 6) +
                           theme_void()
      )
      
      return(plot_obj)
      
    }, error = function(e) {
      cat("Error in plot renderer:", e$message, "\n")
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = paste("Error generating plot:", e$message), 
                        size = 6) +
               theme_void())
    })
  }, height = 600)
  
  
  
  # Plot download handler
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(input$plot_type, "_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot_obj <- switch(input$plot_type,
                         "manhattan" = plot_manhattan_faceted(values$multi_trait_unified_df),
                         "volcano" = plot_volcano_by_trait(values$multi_trait_unified_df),
                         "qq" = plot_qq_by_trait(values$multi_trait_unified_df),
                         "effect_dist" = plot_effect_distribution(values$multi_trait_unified_df),
                         "effect_cor" = plot_effect_correlation(values$multi_trait_unified_df),
                         "sig_snps" = plot_sig_snps_by_trait(values$multi_trait_unified_df),
                         "pleiotropy" = plot_pleiotropy_analysis(values$multi_trait_unified_df),
                         "top_snps" = plot_top_snps_across_traits(values$multi_trait_unified_df, 
                                                                  top_n = input$top_n_snps),
                         "pvalue_dist" = plot_pvalue_distribution(values$multi_trait_unified_df),
                         "effect_dir" = plot_effect_direction(values$multi_trait_unified_df),
                         NULL
      )
      
      if (!is.null(plot_obj)) {
        ggsave(file, plot = plot_obj, width = 12, height = 8, dpi = 300)
      }
    }
  )
  
  
  
  # Create unified format
  observeEvent(input$create_unified_format, {
    req(values$multi_trait_results)
    
    withProgress(message = 'Creating unified format...', value = 0, {
      
      tryCatch({
        incProgress(0.5, detail = "Processing...")
        
        # Create unified dataframe
        values$multi_trait_unified_df <- create_unified_gwas_tidy(
          values$multi_trait_results$trait_gwas_results
        )
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Unified format created successfully!", 
                         type = "message", duration = 3)
        
      }, error = function(e) {
        showNotification(paste("Error creating unified format:", e$message), 
                         type = "error", duration = 5)
      })
    })
  })
  
  
  
  
  # Unified data summary
  output$unified_data_summary <- renderPrint({
    req(values$multi_trait_unified_df)
    
    df <- values$multi_trait_unified_df
    
    cat("=== UNIFIED MULTI-TRAIT DATA SUMMARY ===\n\n")
    cat("Data Overview:\n")
    cat("  Total rows:", nrow(df), "\n")
    cat("  Unique SNPs:", n_distinct(df$SNP_ID), "\n")
    cat("  Traits analyzed:", paste(unique(df$Trait), collapse = ", "), "\n")
    cat("  Number of traits:", n_distinct(df$Trait), "\n\n")
    
    cat("P-value Statistics:\n")
    if ("P_value" %in% colnames(df)) {
      cat("  Min P-value:", format(min(df$P_value, na.rm = TRUE), scientific = TRUE, digits = 3), "\n")
      cat("  Max P-value:", format(max(df$P_value, na.rm = TRUE), scientific = TRUE, digits = 3), "\n")
      cat("  Mean P-value:", format(mean(df$P_value, na.rm = TRUE), scientific = TRUE, digits = 3), "\n")
      
      # Significant SNPs
      sig_05 <- sum(df$P_value < 0.05, na.rm = TRUE)
      sig_01 <- sum(df$P_value < 0.01, na.rm = TRUE)
      sig_001 <- sum(df$P_value < 0.001, na.rm = TRUE)
      
      cat("\nSignificant SNPs (per trait-trait combination):\n")
      cat("  p < 0.05:", sig_05, paste0("(", round(sig_05/nrow(df)*100, 1), "%)\n"))
      cat("  p < 0.01:", sig_01, paste0("(", round(sig_01/nrow(df)*100, 1), "%)\n"))
      cat("  p < 0.001:", sig_001, paste0("(", round(sig_001/nrow(df)*100, 1), "%)\n"))
    }
    
    if ("Effect" %in% colnames(df)) {
      cat("\nEffect Size Statistics:\n")
      cat("  Min Effect:", round(min(df$Effect, na.rm = TRUE), 4), "\n")
      cat("  Max Effect:", round(max(df$Effect, na.rm = TRUE), 4), "\n")
      cat("  Mean Effect:", round(mean(df$Effect, na.rm = TRUE), 4), "\n")
      cat("  SD Effect:", round(sd(df$Effect, na.rm = TRUE), 4), "\n")
    }
    
    cat("\nMemory Usage:\n")
    cat("  Object size:", format(object.size(df), units = "MB"), "\n")
    cat("  Columns:", ncol(df), "\n")
  })
  
  
  
  # Plot unified data
  observeEvent(input$plot_unified_data, {
    req(values$multi_trait_unified_df, input$unified_plot_type)
    
    output$unified_plot_output <- renderPlot({
      tryCatch({
        df <- values$multi_trait_unified_df
        
        plot_obj <- switch(input$unified_plot_type,
                           "manhattan" = plot_manhattan_faceted(df),
                           "volcano" = plot_volcano_by_trait(df),
                           "effect_dist" = plot_effect_distribution(df),
                           "qq" = plot_qq_by_trait(df),
                           ggplot() + annotate("text", x=0.5, y=0.5, 
                                               label="Select a plot type", size=6))
        
        plot_obj
      }, error = function(e) {
        ggplot() + annotate("text", x=0.5, y=0.5, 
                            label=paste("Error:", e$message), size=6)
      })
    })
  })
  
  
  
  
  output$download_unified_data <- downloadHandler(
    filename = function() {
      paste0("unified_multi_trait_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$multi_trait_unified_df)
      write.csv(values$multi_trait_unified_df, file, row.names = FALSE)
    }
  )
  
  # Download enhanced unified data
  output$download_enhanced_unified_data <- downloadHandler(
    filename = function() {
      paste0("enhanced_unified_multi_trait_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$multi_trait_enhanced_unified_df)
      write.csv(values$multi_trait_enhanced_unified_df, file, row.names = FALSE)
    }
  )
  
  # Download complete integrated data
  output$download_complete_integrated_data <- downloadHandler(
    filename = function() {
      paste0("complete_integrated_multi_trait_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$multi_trait_complete_integrated_df)
      write.csv(values$multi_trait_complete_integrated_df, file, row.names = FALSE)
    }
  )
  
  
  
  # Download wide data
  output$download_wide_data <- downloadHandler(
    filename = function() {
      paste0("wide_multi_trait_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$multi_trait_wide_df)
      write.csv(values$multi_trait_wide_df, file, row.names = FALSE)
    }
  )
  
  
  
  # Unified data status
  output$unified_data_status_ui <- renderUI({
    if (!is.null(values$multi_trait_unified_df)) {
      df <- values$multi_trait_unified_df
      
      div(class = "alert alert-success",
          icon("check-circle"),
          strong(" ✓ Unified Data Available"),
          br(),
          paste("Traits:", n_distinct(df$Trait)),
          br(),
          paste("SNP-trait combinations:", nrow(df)),
          br(),
          paste("Memory:", format(object.size(df), units = "MB"))
      )
    } else {
      div(class = "alert alert-warning",
          icon("exclamation-triangle"),
          strong(" No Unified Data"),
          br(),
          "Run multi-trait analysis first"
      )
    }
  })
  
  
  # Update trait selector for filtering
  observe({
    req(values$multi_trait_unified_df)
    
    traits <- unique(values$multi_trait_unified_df$Trait)
    updateSelectInput(session, "sig_trait", 
                      choices = c("All Traits", traits))
  })
  
  # Find top SNPs function
  observeEvent(input$find_top_snps, {
    req(values$multi_trait_unified_df)
    
    top_n <- input$top_n_global
    
    top_snps <- values$multi_trait_unified_df %>%
      group_by(SNP_ID) %>%
      summarise(
        Min_P_value = min(P_value, na.rm = TRUE),
        N_Traits = n_distinct(Trait),
        Traits = paste(sort(unique(Trait)), collapse = ", "),
        Mean_Effect = mean(Effect, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(Min_P_value) %>%
      head(top_n) %>%
      mutate(
        Min_P_value = format(Min_P_value, scientific = TRUE, digits = 3),
        Mean_Effect = round(Mean_Effect, 4)
      )
    
    output$top_snps_table <- renderTable({
      top_snps
    }, striped = TRUE, hover = TRUE)
  })
  
  # Filter by significance
  observeEvent(input$apply_filter, {
    req(values$multi_trait_unified_df)
    
    df <- values$multi_trait_unified_df
    
    if (input$sig_trait != "All Traits") {
      df <- df %>% filter(Trait == input$sig_trait)
    }
    
    filtered <- df %>% filter(P_value < input$sig_threshold)
    
    output$filter_summary <- renderText({
      paste("Found", nrow(filtered), "significant SNP-trait combinations",
            ifelse(input$sig_trait != "All Traits", 
                   paste("for trait:", input$sig_trait), 
                   "across all traits"))
    })
  })
  
  # Helper function for wide format creation (add to server)
  create_wide_gwas_format <- function(trait_gwas_results) {
    require(dplyr)
    require(tidyr)
    
    if (is.null(trait_gwas_results)) return(NULL)
    
    # Create long format first
    unified_list <- list()
    for (trait_name in names(trait_gwas_results)) {
      df <- trait_gwas_results[[trait_name]]
      if (!is.null(df) && nrow(df) > 0) {
        df$Trait <- trait_name
        unified_list[[trait_name]] <- df
      }
    }
    
    if (length(unified_list) == 0) return(NULL)
    
    long_df <- do.call(rbind, unified_list)
    rownames(long_df) <- NULL
    
    # Create wide format for P-values
    pvalue_wide <- long_df %>%
      select(SNP_ID, Chromosome, Position, Trait, P_value) %>%
      pivot_wider(
        id_cols = c(SNP_ID, Chromosome, Position),
        names_from = Trait,
        values_from = P_value,
        names_prefix = "P_"
      )
    
    # Create wide format for Effects
    effect_wide <- long_df %>%
      select(SNP_ID, Trait, Effect) %>%
      pivot_wider(
        id_cols = SNP_ID,
        names_from = Trait,
        values_from = Effect,
        names_prefix = "Effect_"
      )
    
    # Combine
    wide_df <- pvalue_wide %>%
      left_join(effect_wide, by = "SNP_ID")
    
    return(wide_df)
  }
  
  
  # ==========================================================================
  # CROSS-SNP ANALYSIS OBSERVERS
  # ==========================================================================
  
  # In your server function
  observeEvent(input$run_cross_snp_analysis, {
    req(values$diallel_data, values$multi_trait_results)
    
    withProgress(message = "Creating Cross-SNP DataFrame...", value = 0, {
      incProgress(0.3, detail = "Processing traits...")
      
      # Use the safe version
      cross_snp_df <- create_cross_snp_dataframe_all_traits_fixed(
        multi_trait_results = values$multi_trait_results,
        diallel_data = values$diallel_data,
        processed_geno = values$processed_geno,
        p_threshold = input$cross_snp_p_threshold,  # Add this input to your UI
        max_snps_per_trait = 500,
        include_parent_genotypes = FALSE
      )
      
      if (is.null(cross_snp_df)) {
        # Try the simplified version
        incProgress(0.6, detail = "Trying simplified version...")
        cross_snp_df <- create_cross_snp_dataframe_simple(
          matched_data = values$matched_data,
          multi_trait_results = values$multi_trait_results,
          processed_geno = values$processed_geno,
          diallel_data = values$diallel_data,
          p_threshold = 0.05,
          max_snps_per_trait = 200
        )
      }
      
      if (!is.null(cross_snp_df)) {
        values$cross_snp_df <- cross_snp_df
        values$cross_snp_detailed <- cross_snp_df
        values$dl_cross_snp_detailed <- cross_snp_df
        print(values$cross_snp_df)
        
        incProgress(1.0, detail = "Complete!")
        showNotification("Cross-SNP analysis completed!", type = "message")
      } else {
        showNotification("Failed to create cross-SNP dataframe", type = "error")
      }
    })
  })
  
  
 
  # Debug function to check data structure
  check_cross_snp_data_structure <- function() {
    cat("\n=== CROSS-SNP DATA STRUCTURE CHECK ===\n")
    
    if (!is.null(values$cross_summary)) {
      cat("\n1. CROSS_SUMMARY DATA:\n")
      cat("   Class:", class(values$cross_summary), "\n")
      cat("   Dimensions:", dim(values$cross_summary), "\n")
      cat("   Column names:", paste(colnames(values$cross_summary), collapse=", "), "\n")
      cat("   Column types:\n")
      for (col in colnames(values$cross_summary)) {
        cat("   -", col, ":", class(values$cross_summary[[col]]), 
            "| Unique values:", length(unique(values$cross_summary[[col]])), "\n")
      }
      cat("\n   First 3 rows:\n")
      print(head(values$cross_summary, 3))
    } else {
      cat("\n1. CROSS_SUMMARY is NULL\n")
    }
    
    if (!is.null(values$cross_snp_detailed)) {
      cat("\n2. CROSS_SNP_DETAILED DATA:\n")
      cat("   Class:", class(values$cross_snp_detailed), "\n")
      cat("   Dimensions:", dim(values$cross_snp_detailed), "\n")
      cat("   Column names:", paste(colnames(values$cross_snp_detailed), collapse=", "), "\n")
      cat("   First 3 rows:\n")
      print(head(values$cross_snp_detailed, 3))
    } else {
      cat("\n2. CROSS_SNP_DETAILED is NULL\n")
    }
    
    cat("\n3. CHECKING FOR KEY COLUMNS:\n")
    
    if (!is.null(values$cross_summary)) {
      # Check for parent columns
      parent_cols <- grep("parent|Parent|p1|p2|P1|P2", colnames(values$cross_summary), 
                          value = TRUE, ignore.case = TRUE)
      cat("   Possible parent columns:", paste(parent_cols, collapse=", "), "\n")
      
      # Check for score columns
      score_cols <- colnames(values$cross_summary)[sapply(values$cross_summary, is.numeric)]
      cat("   Numeric (score) columns:", paste(score_cols, collapse=", "), "\n")
      
      # Check for cross ID columns
      id_cols <- grep("cross|Cross|id|ID|name|Name", colnames(values$cross_summary), 
                      value = TRUE, ignore.case = TRUE)
      cat("   Possible ID columns:", paste(id_cols, collapse=", "), "\n")
    }
    
    cat("==========================================\n\n")
  }
  
  # Add this observer to run the check when data changes
  observe({
    if (!is.null(values$cross_summary)) {
      check_cross_snp_data_structure()
    }
  })
  
  
  # Cross LOD Distribution Plot - FIXED
  output$cross_lod_plot <- renderPlot({
    req(values$cross_summary)
    
    tryCatch({
      cat("\n=== Creating Cross LOD Distribution Plot ===\n")
      cat("cross_summary columns:", colnames(values$cross_summary), "\n")
      
      # Find score column
      score_patterns <- c("LOD", "lod", "Lod", "LogP", "logP", "neg_log10_p", 
                          "log10p", "P_value", "p_value", "P", "p", "Score", "score")
      score_col <- NULL
      
      for (pattern in score_patterns) {
        if (pattern %in% colnames(values$cross_summary)) {
          score_col <- pattern
          cat("Found score column:", score_col, "\n")
          break
        }
      }
      
      # Try case-insensitive
      if (is.null(score_col)) {
        colnames_lower <- tolower(colnames(values$cross_summary))
        for (pattern in tolower(score_patterns)) {
          matches <- which(colnames_lower == pattern)
          if (length(matches) > 0) {
            score_col <- colnames(values$cross_summary)[matches[1]]
            cat("Found score column (case-insensitive):", score_col, "\n")
            break
          }
        }
      }
      
      # If still not found, use first numeric column
      if (is.null(score_col)) {
        numeric_cols <- sapply(values$cross_summary, is.numeric)
        if (any(numeric_cols)) {
          score_col <- colnames(values$cross_summary)[numeric_cols][1]
          cat("Using first numeric column as score:", score_col, "\n")
        } else {
          stop("No numeric columns found in cross_summary")
        }
      }
      
      # Get the actual score values
      score_values <- values$cross_summary[[score_col]]
      score_values <- score_values[!is.na(score_values)]
      
      if (length(score_values) == 0) {
        stop("No valid score values after removing NAs")
      }
      
      cat("Using", length(score_values), "score values for histogram\n")
      cat("Score range:", range(score_values, na.rm = TRUE), "\n")
      
      # Create histogram
      p <- ggplot(data.frame(Score = score_values), aes(x = Score)) +
        geom_histogram(fill = "steelblue", alpha = 0.7, bins = 30, color = "white") +
        labs(title = paste("Cross Score Distribution (", score_col, ")"),
             x = score_col,
             y = "Frequency") +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)
        )
      
      # Add threshold line if available and relevant
      if (exists("input$cross_snp_lod_threshold") && 
          !is.null(input$cross_snp_lod_threshold) &&
          score_col %in% c("LOD", "lod", "LogP", "logP", "neg_log10_p")) {
        p <- p + 
          geom_vline(xintercept = input$cross_snp_lod_threshold, 
                     color = "red", linetype = "dashed", size = 1) +
          annotate("text", 
                   x = input$cross_snp_lod_threshold, 
                   y = Inf,
                   label = paste("Threshold =", input$cross_snp_lod_threshold),
                   vjust = 1.5, hjust = -0.1,
                   color = "red", size = 4)
      }
      
      return(p)
      
    }, error = function(e) {
      cat("Error in cross_lod_plot:", e$message, "\n")
      
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Error creating plot:", e$message, "\n\n",
                               "Available columns:", 
                               paste(colnames(values$cross_summary), collapse = ", ")),
                 size = 4) +
        theme_void()
    })
  })
  
  
  
  # SNP Effect Heatmap - Debugged
  output$snp_effect_heatmap <- renderPlot({
    req(values$cross_snp_detailed)
    
    cat("Debug: cross_snp_detailed structure:\n")
    cat("Class:", class(values$cross_snp_detailed), "\n")
    cat("Dimensions:", dim(values$cross_snp_detailed), "\n")
    if (is.data.frame(values$cross_snp_detailed) && nrow(values$cross_snp_detailed) > 0) {
      cat("Column names:", colnames(values$cross_snp_detailed), "\n")
    }
    
    tryCatch({
      # Check for required columns
      required_cols <- c("SNP_ID", "CrossID", "Effect")
      missing_cols <- setdiff(required_cols, colnames(values$cross_snp_detailed))
      
      if (length(missing_cols) > 0) {
        stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
      }
      
      # Prepare data for heatmap
      heatmap_data <- values$cross_snp_detailed %>%
        group_by(SNP_ID, CrossID) %>%
        summarise(Avg_Effect = mean(Effect, na.rm = TRUE), .groups = "drop") %>%
        filter(!is.na(Avg_Effect))
      
      if (nrow(heatmap_data) == 0) {
        stop("No valid effect data after filtering")
      }
      
      # Reshape to matrix format
      heatmap_matrix <- heatmap_data %>%
        pivot_wider(names_from = CrossID, values_from = Avg_Effect, values_fill = NA) %>%
        column_to_rownames("SNP_ID") %>%
        as.matrix()
      
      # Check if matrix has data
      if (all(is.na(heatmap_matrix))) {
        stop("All effect values are NA")
      }
      
      # Limit dimensions for readability
      max_rows <- 50
      max_cols <- 30
      
      if (nrow(heatmap_matrix) > max_rows) {
        # Select SNPs with highest variance
        row_vars <- apply(heatmap_matrix, 1, var, na.rm = TRUE)
        if (all(is.na(row_vars))) {
          # If all variances are NA, use row means
          row_vars <- apply(heatmap_matrix, 1, function(x) abs(mean(x, na.rm = TRUE)))
        }
        top_indices <- order(row_vars, decreasing = TRUE)[1:max_rows]
        heatmap_matrix <- heatmap_matrix[top_indices, ]
      }
      
      if (ncol(heatmap_matrix) > max_cols) {
        # Select crosses with most data
        col_counts <- apply(heatmap_matrix, 2, function(x) sum(!is.na(x)))
        top_cols <- order(col_counts, decreasing = TRUE)[1:max_cols]
        heatmap_matrix <- heatmap_matrix[, top_cols]
      }
      
      # Create heatmap
      pheatmap(heatmap_matrix,
               main = "SNP Effect Heatmap Across Crosses",
               color = colorRampPalette(c("blue", "white", "red"))(50),
               show_rownames = ifelse(nrow(heatmap_matrix) <= 30, TRUE, FALSE),
               show_colnames = ifelse(ncol(heatmap_matrix) <= 30, TRUE, FALSE),
               fontsize_row = 8,
               fontsize_col = 8,
               na_col = "gray90",
               clustering_method = "ward.D2")
      
    }, error = function(e) {
      cat("Error in snp_effect_heatmap:", e$message, "\n")
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Error creating heatmap:", e$message),
                 size = 6) +
        theme_void()
    })
  })
  
  
  
  # Add this observer to see what data you actually have
  observe({
    req(values$cross_summary)
    
    cat("\n=== DEBUG: Cross-SNP Data Structure ===\n")
    
    # Check cross_summary
    if (!is.null(values$cross_summary)) {
      cat("cross_summary is a:", class(values$cross_summary), "\n")
      cat("Dimensions:", dim(values$cross_summary), "\n")
      cat("Column names:", colnames(values$cross_summary), "\n")
      cat("\nFirst 3 rows of cross_summary:\n")
      print(head(values$cross_summary, 3))
    } else {
      cat("cross_summary is NULL\n")
    }
    
    # Check cross_snp_detailed
    if (!is.null(values$cross_snp_detailed)) {
      cat("\ncross_snp_detailed columns:", colnames(values$cross_snp_detailed), "\n")
    } else {
      cat("\ncross_snp_detailed is NULL\n")
    }
    
    cat("======================================\n\n")
  })
  
  
  
  # Interactive: Cross Performance - FIXED with robust column detection
  output$cross_performance_plot <- renderPlotly({
    req(values$cross_summary)
    
    tryCatch({
      # Debug output
      cat("\n=== Creating Cross Performance Plot ===\n")
      cat("Data columns:", paste(colnames(values$cross_summary), collapse = ", "), "\n")
      cat("Data class:", class(values$cross_summary), "\n")
      
      # Find the cross identifier column
      cross_id_patterns <- c("CrossID", "Cross_Id", "CrossId", "crossid", "cross_id", 
                             "Cross", "cross", "ID", "Id", "id", "Name", "name")
      
      cross_id_col <- NULL
      for (pattern in cross_id_patterns) {
        if (pattern %in% colnames(values$cross_summary)) {
          cross_id_col <- pattern
          cat("Found cross ID column:", cross_id_col, "\n")
          break
        }
      }
      
      # Try case-insensitive matching
      if (is.null(cross_id_col)) {
        colnames_lower <- tolower(colnames(values$cross_summary))
        for (pattern in tolower(cross_id_patterns)) {
          matches <- which(colnames_lower == pattern)
          if (length(matches) > 0) {
            cross_id_col <- colnames(values$cross_summary)[matches[1]]
            cat("Found cross ID column (case-insensitive):", cross_id_col, "\n")
            break
          }
        }
      }
      
      # If still not found, use first column
      if (is.null(cross_id_col)) {
        cross_id_col <- colnames(values$cross_summary)[1]
        cat("Using first column as cross ID:", cross_id_col, "\n")
      }
      
      # Find score column
      score_patterns <- c("LOD", "lod", "Lod", "LogP", "logP", "neg_log10_p", 
                          "log10p", "P_value", "p_value", "P", "p", "Score", "score",
                          "Value", "value", "Effect", "effect")
      
      score_col <- NULL
      for (pattern in score_patterns) {
        if (pattern %in% colnames(values$cross_summary)) {
          score_col <- pattern
          cat("Found score column:", score_col, "\n")
          break
        }
      }
      
      # Try case-insensitive for score
      if (is.null(score_col)) {
        colnames_lower <- tolower(colnames(values$cross_summary))
        for (pattern in tolower(score_patterns)) {
          matches <- which(colnames_lower == pattern)
          if (length(matches) > 0) {
            score_col <- colnames(values$cross_summary)[matches[1]]
            cat("Found score column (case-insensitive):", score_col, "\n")
            break
          }
        }
      }
      
      # If still not found, use first numeric column
      if (is.null(score_col)) {
        numeric_cols <- sapply(values$cross_summary, is.numeric)
        if (any(numeric_cols)) {
          score_col <- colnames(values$cross_summary)[numeric_cols][1]
          cat("Using first numeric column as score:", score_col, "\n")
        } else {
          stop("No numeric columns found for score")
        }
      }
      
      # Prepare data
      plot_data <- values$cross_summary
      plot_data$CrossIdentifier <- plot_data[[cross_id_col]]
      plot_data$ScoreValue <- plot_data[[score_col]]
      
      # Remove any NA scores
      plot_data <- plot_data[!is.na(plot_data$ScoreValue), ]
      
      if (nrow(plot_data) == 0) {
        stop("No valid scores after removing NA values")
      }
      
      # Sort by score
      plot_data <- plot_data[order(-plot_data$ScoreValue), ]
      
      # Limit for readability
      max_crosses <- 30
      if (nrow(plot_data) > max_crosses) {
        plot_data <- head(plot_data, max_crosses)
        cat("Limited to top", max_crosses, "crosses\n")
      }
      
      cat("Plotting", nrow(plot_data), "crosses\n")
      
      # Create plot
      p <- ggplot(plot_data, 
                  aes(x = reorder(CrossIdentifier, -ScoreValue), 
                      y = ScoreValue,
                      text = paste("Cross:", CrossIdentifier,
                                   "<br>", score_col, ":", round(ScoreValue, 2)))) +
        geom_bar(stat = "identity", fill = "#3498db", alpha = 0.8) +
        labs(title = paste("Cross Performance (", score_col, ")"),
             x = "Cross",
             y = score_col) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggplotly(p, tooltip = "text") %>%
        layout(
          hoverlabel = list(bgcolor = "white", font = list(size = 12)),
          xaxis = list(tickangle = -45)
        )
      
    }, error = function(e) {
      cat("Error in cross_performance_plot:", e$message, "\n")
      
      # Create error plot with more information
      error_plot <- plot_ly() %>%
        add_annotations(
          text = paste("Error creating plot:", e$message,
                       "\n\nAvailable columns:", 
                       paste(colnames(values$cross_summary), collapse = ", ")),
          xref = "paper",
          yref = "paper",
          x = 0.5,
          y = 0.5,
          showarrow = FALSE,
          font = list(size = 14)
        ) %>%
        layout(
          title = "Plot Error",
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
        )
      
      return(error_plot)
    })
  })
  
  
  
  
  # Interactive: SNP Distribution - FIXED
  output$snp_distribution_plot <- renderPlotly({
    req(values$cross_snp_detailed)
    
    tryCatch({
      cat("\n=== Creating SNP Distribution Plot ===\n")
      cat("cross_snp_detailed columns:", colnames(values$cross_snp_detailed), "\n")
      
      # Find necessary columns
      required_cols <- list(
        snp_id = c("SNP_ID", "SNP", "Snp", "snp", "Marker", "marker"),
        chromosome = c("Chromosome", "CHROM", "chrom", "Chr", "chr"),
        effect = c("Effect", "effect", "Beta", "beta", "Coefficient", "coef"),
        score = c("LOD", "lod", "LogP", "logP", "neg_log10_p", "P_value", "p")
      )
      
      # Find actual column names
      found_cols <- list()
      for (col_type in names(required_cols)) {
        patterns <- required_cols[[col_type]]
        found <- NULL
        
        for (pattern in patterns) {
          if (pattern %in% colnames(values$cross_snp_detailed)) {
            found <- pattern
            break
          }
        }
        
        # Try case-insensitive
        if (is.null(found)) {
          colnames_lower <- tolower(colnames(values$cross_snp_detailed))
          for (pattern in tolower(patterns)) {
            matches <- which(colnames_lower == pattern)
            if (length(matches) > 0) {
              found <- colnames(values$cross_snp_detailed)[matches[1]]
              break
            }
          }
        }
        
        if (is.null(found)) {
          cat("Warning: Could not find", col_type, "column\n")
        } else {
          cat("Using", found, "for", col_type, "\n")
          found_cols[[col_type]] <- found
        }
      }
      
      # Check if we have minimum required columns
      if (is.null(found_cols$snp_id) || is.null(found_cols$score)) {
        stop("Missing required columns (SNP_ID and score)")
      }
      
      # Prepare data
      plot_data <- values$cross_snp_detailed
      plot_data$SNP <- plot_data[[found_cols$snp_id]]
      plot_data$Score <- plot_data[[found_cols$score]]
      
      # Add effect if available
      if (!is.null(found_cols$effect)) {
        plot_data$Effect <- plot_data[[found_cols$effect]]
      } else {
        plot_data$Effect <- 0  # Default
      }
      
      # Add chromosome if available
      if (!is.null(found_cols$chromosome)) {
        plot_data$Chrom <- as.character(plot_data[[found_cols$chromosome]])
      } else {
        plot_data$Chrom <- "Unknown"
      }
      
      # Summarize by SNP
      snp_summary <- plot_data %>%
        group_by(SNP, Chrom) %>%
        summarise(
          Avg_Score = mean(Score, na.rm = TRUE),
          Avg_Effect = mean(Effect, na.rm = TRUE),
          Count = n(),
          .groups = 'drop'
        ) %>%
        filter(!is.na(Avg_Score))
      
      if (nrow(snp_summary) == 0) {
        stop("No valid SNP data after summarization")
      }
      
      cat("Summarized", nrow(snp_summary), "SNPs\n")
      
      # Create plot
      p <- ggplot(snp_summary, 
                  aes(x = Chrom, y = Avg_Score, 
                      size = abs(Avg_Effect), color = Avg_Effect,
                      text = paste("SNP:", SNP,
                                   "<br>Chromosome:", Chrom,
                                   "<br>Avg Score:", round(Avg_Score, 2),
                                   "<br>Avg Effect:", round(Avg_Effect, 3),
                                   "<br>Count:", Count))) +
        geom_point(alpha = 0.7, position = position_jitter(width = 0.2)) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                              midpoint = 0, name = "Effect") +
        scale_size_continuous(range = c(3, 12), name = "|Effect|") +
        labs(title = "SNP Distribution by Chromosome",
             x = "Chromosome",
             y = "Average Score") +
        theme_minimal()
      
      ggplotly(p, tooltip = "text") %>%
        layout(hoverlabel = list(bgcolor = "white"))
      
    }, error = function(e) {
      cat("Error in snp_distribution_plot:", e$message, "\n")
      
      plotly_empty() %>%
        layout(
          title = paste("Error:", e$message),
          annotations = list(
            text = paste("Error:", e$message, 
                         "\n\nAvailable columns:", 
                         paste(colnames(values$cross_snp_detailed), collapse = ", ")),
            xref = "paper",
            yref = "paper",
            x = 0.5,
            y = 0.5,
            showarrow = FALSE
          )
        )
    })
  })
  
  
  # Interactive: Parent Contribution - FIXED with robust column detection
  output$parent_contribution_plot <- renderPlotly({
    req(values$cross_summary)
    
    tryCatch({
      cat("\n=== Creating Parent Contribution Plot ===\n")
      cat("cross_summary columns:", colnames(values$cross_summary), "\n")
      
      # Find parent1 column
      parent1_patterns <- c("Parent1", "parent1", "Parent_1", "parent_1", 
                            "P1", "p1", "Male", "male", "Female", "female")
      parent1_col <- NULL
      for (pattern in parent1_patterns) {
        if (pattern %in% colnames(values$cross_summary)) {
          parent1_col <- pattern
          cat("Found Parent1 column:", parent1_col, "\n")
          break
        }
      }
      
      # Try case-insensitive
      if (is.null(parent1_col)) {
        colnames_lower <- tolower(colnames(values$cross_summary))
        for (pattern in tolower(parent1_patterns)) {
          matches <- which(colnames_lower == pattern)
          if (length(matches) > 0) {
            parent1_col <- colnames(values$cross_summary)[matches[1]]
            cat("Found Parent1 column (case-insensitive):", parent1_col, "\n")
            break
          }
        }
      }
      
      # Find parent2 column
      parent2_patterns <- c("Parent2", "parent2", "Parent_2", "parent_2", 
                            "P2", "p2", "Male2", "male2", "Female2", "female2")
      parent2_col <- NULL
      for (pattern in parent2_patterns) {
        if (pattern %in% colnames(values$cross_summary)) {
          parent2_col <- pattern
          cat("Found Parent2 column:", parent2_col, "\n")
          break
        }
      }
      
      # Try case-insensitive
      if (is.null(parent2_col)) {
        colnames_lower <- tolower(colnames(values$cross_summary))
        for (pattern in tolower(parent2_patterns)) {
          matches <- which(colnames_lower == pattern)
          if (length(matches) > 0) {
            parent2_col <- colnames(values$cross_summary)[matches[1]]
            cat("Found Parent2 column (case-insensitive):", parent2_col, "\n")
            break
          }
        }
      }
      
      # Find score column
      score_patterns <- c("LOD", "lod", "Lod", "LogP", "logP", "neg_log10_p", 
                          "log10p", "P_value", "p_value", "P", "p", "Score", "score",
                          "Value", "value", "Effect", "effect")
      score_col <- NULL
      for (pattern in score_patterns) {
        if (pattern %in% colnames(values$cross_summary)) {
          score_col <- pattern
          cat("Found score column:", score_col, "\n")
          break
        }
      }
      
      # Try case-insensitive
      if (is.null(score_col)) {
        colnames_lower <- tolower(colnames(values$cross_summary))
        for (pattern in tolower(score_patterns)) {
          matches <- which(colnames_lower == pattern)
          if (length(matches) > 0) {
            score_col <- colnames(values$cross_summary)[matches[1]]
            cat("Found score column (case-insensitive):", score_col, "\n")
            break
          }
        }
      }
      
      # If still not found, use first numeric column
      if (is.null(score_col)) {
        numeric_cols <- sapply(values$cross_summary, is.numeric)
        if (any(numeric_cols)) {
          score_col <- colnames(values$cross_summary)[numeric_cols][1]
          cat("Using first numeric column as score:", score_col, "\n")
        } else {
          stop("No numeric columns found for score")
        }
      }
      
      # Check if we found parent columns
      if (is.null(parent1_col) || is.null(parent2_col)) {
        # Try to find any two columns that might be parents
        # Look for columns with parent-like values
        possible_parent_cols <- character(0)
        
        for (col in colnames(values$cross_summary)) {
          if (is.character(values$cross_summary[[col]]) || is.factor(values$cross_summary[[col]])) {
            unique_vals <- unique(values$cross_summary[[col]])
            # If column has values that look like parent names (not too many unique values)
            if (length(unique_vals) > 1 && length(unique_vals) < nrow(values$cross_summary)/2) {
              possible_parent_cols <- c(possible_parent_cols, col)
            }
          }
        }
        
        if (length(possible_parent_cols) >= 2) {
          parent1_col <- possible_parent_cols[1]
          parent2_col <- possible_parent_cols[2]
          cat("Using", parent1_col, "and", parent2_col, "as parent columns\n")
        } else {
          stop("Could not identify parent columns")
        }
      }
      
      # Prepare data for parent1 contributions
      parent1_contrib <- values$cross_summary %>%
        rename(Parent = !!parent1_col, Score = !!score_col) %>%
        filter(!is.na(Parent), !is.na(Score)) %>%
        group_by(Parent) %>%
        summarise(
          N_Crosses = n(),
          Avg_Score = mean(Score, na.rm = TRUE),
          Max_Score = max(Score, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        mutate(Role = "Parent1")
      
      # Prepare data for parent2 contributions
      parent2_contrib <- values$cross_summary %>%
        rename(Parent = !!parent2_col, Score = !!score_col) %>%
        filter(!is.na(Parent), !is.na(Score)) %>%
        group_by(Parent) %>%
        summarise(
          N_Crosses = n(),
          Avg_Score = mean(Score, na.rm = TRUE),
          Max_Score = max(Score, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        mutate(Role = "Parent2")
      
      # Combine parent contributions
      parent_contrib <- bind_rows(parent1_contrib, parent2_contrib) %>%
        group_by(Parent) %>%
        summarise(
          Total_Crosses = sum(N_Crosses, na.rm = TRUE),
          Overall_Avg_Score = mean(Avg_Score, na.rm = TRUE),
          Overall_Max_Score = max(Max_Score, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        filter(!is.na(Overall_Avg_Score))
      
      if (nrow(parent_contrib) == 0) {
        stop("No valid parent contribution data after processing")
      }
      
      cat("Parent contribution data has", nrow(parent_contrib), "parents\n")
      
      # Create interactive plot
      p <- ggplot(parent_contrib, 
                  aes(x = Total_Crosses, y = Overall_Avg_Score, 
                      size = Overall_Max_Score, 
                      color = Parent,
                      text = paste("Parent:", Parent,
                                   "<br>Total Crosses:", Total_Crosses,
                                   "<br>Avg Score:", round(Overall_Avg_Score, 2),
                                   "<br>Max Score:", round(Overall_Max_Score, 2)))) +
        geom_point(alpha = 0.8) +
        scale_size_continuous(range = c(5, 15), name = "Max Score") +
        labs(title = "Parent Contribution to Crosses",
             x = "Number of Crosses",
             y = "Average Score") +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggplotly(p, tooltip = "text") %>%
        layout(
          hoverlabel = list(bgcolor = "white"),
          xaxis = list(title = "Number of Crosses"),
          yaxis = list(title = "Average Score")
        )
      
    }, error = function(e) {
      cat("Error in parent_contribution_plot:", e$message, "\n")
      
      # Create error plot with debug information
      error_plot <- plot_ly() %>%
        add_annotations(
          text = paste("Error creating Parent Contribution plot:", e$message,
                       "\n\nAvailable columns:", 
                       paste(colnames(values$cross_summary), collapse = ", "),
                       "\n\nParent1 patterns tried:", paste(c("Parent1", "parent1", "Parent_1"), collapse=", "),
                       "\nParent2 patterns tried:", paste(c("Parent2", "parent2", "Parent_2"), collapse=", "),
                       "\nScore patterns tried:", paste(c("LOD", "Score", "P_value"), collapse=", ")),
          xref = "paper",
          yref = "paper",
          x = 0.5,
          y = 0.5,
          showarrow = FALSE,
          font = list(size = 10)
        ) %>%
        layout(
          title = "Parent Contribution Plot Error",
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
        )
      
      return(error_plot)
    })
  })
  
  
  
  # Enhanced pleiotropy plot (FIXED)
  output$enhanced_pleiotropy_plot <- renderPlot({
    req(values$multi_trait_enhanced_unified_df)
    
    # Use the fixed plot_pleiotropy_analysis function
    plot_pleiotropy_analysis(values$multi_trait_enhanced_unified_df)
  })
  
  
  
  # Cross performance plot
  plot_cross_performance <- function(cross_summary_data) {
    if (is.null(cross_summary_data) || nrow(cross_summary_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No cross performance data available", size = 6) +
               theme_void())
    }
    
    # Ensure required columns exist
    required_cols <- c("Cross", "N_SNPs", "N_Significant")
    if (!all(required_cols %in% colnames(cross_summary_data))) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "Required columns missing", size = 6) +
               theme_void())
    }
    
    # Limit to top crosses for readability
    top_crosses <- cross_summary_data %>%
      arrange(desc(N_Significant)) %>%
      head(min(20, nrow(cross_summary_data)))
    
    ggplot(top_crosses, aes(x = reorder(Cross, N_Significant), y = N_Significant)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_text(aes(label = N_Significant), hjust = -0.3, size = 3.5) +
      coord_flip() +
      labs(
        title = "Cross Performance by Number of Significant SNPs",
        subtitle = paste("Top", nrow(top_crosses), "crosses"),
        x = "Cross",
        y = "Number of Significant SNPs (p < 0.05)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(margin = margin(r = 15))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  }
  
  # SNP distribution plot
  plot_snp_distribution <- function(snp_summary_data) {
    if (is.null(snp_summary_data) || nrow(snp_summary_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No SNP distribution data available", size = 6) +
               theme_void())
    }
    
    # Ensure required columns exist
    required_cols <- c("SNP_ID", "N_Crosses")
    if (!all(required_cols %in% colnames(snp_summary_data))) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "Required columns missing", size = 6) +
               theme_void())
    }
    
    # Create histogram of SNP distribution across crosses
    ggplot(snp_summary_data, aes(x = N_Crosses)) +
      geom_histogram(binwidth = 1, fill = "lightblue", color = "black", alpha = 0.7) +
      geom_density(aes(y = after_stat(count) * 1), color = "red", size = 1) +
      labs(
        title = "Distribution of SNPs Across Crosses",
        subtitle = paste("Total SNPs:", nrow(snp_summary_data)),
        x = "Number of Crosses",
        y = "Number of SNPs"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
  }
  
  # Parent contribution plot
  plot_parent_contribution <- function(cross_summary_data) {
    if (is.null(cross_summary_data) || nrow(cross_summary_data) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "No parent contribution data available", size = 6) +
               theme_void())
    }
    
    # Ensure required columns exist
    required_cols <- c("Parent1", "Parent2")
    if (!all(required_cols %in% colnames(cross_summary_data))) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, 
                        label = "Required columns missing", size = 6) +
               theme_void())
    }
    
    # Calculate parent contributions
    parent_contrib <- data.frame(
      Parent = c(cross_summary_data$Parent1, cross_summary_data$Parent2)
    ) %>%
      group_by(Parent) %>%
      summarise(
        N_Crosses = n(),
        Mean_Significant = mean(cross_summary_data$N_Significant, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(N_Crosses))
    
    # Limit to top parents for readability
    top_parents <- parent_contrib %>%
      head(min(20, nrow(parent_contrib)))
    
    ggplot(top_parents, aes(x = reorder(Parent, N_Crosses), y = N_Crosses)) +
      geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
      geom_text(aes(label = N_Crosses), hjust = -0.3, size = 3.5) +
      coord_flip() +
      labs(
        title = "Parent Contribution to Crosses",
        subtitle = paste("Top", nrow(top_parents), "parents"),
        x = "Parent",
        y = "Number of Crosses"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  }
  
  
  
  output$cross_snp_detailed <- renderDT({
    req(values$cross_snp_detailed)
    
    display_df <- values$cross_snp_detailed %>%
      mutate(
        P_value = format(P_value, scientific = TRUE, digits = 3),
        LOD = round(LOD, 3),
        Effect = round(Effect, 4),
        PVE = round(PVE, 2),
        Cross_Performance = round(Cross_Performance, 2)
      ) %>%
      select(Cross, Parent1, Parent2, SNP_ID, Chromosome, Position, Trait,
             P_value, LOD, Effect, PVE, Cross_Performance, Significant)
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf'),
        order = list(list(8, 'desc'))
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Detailed Cross-SNP Analysis Data",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'Significant',
        target = 'row',
        backgroundColor = styleEqual(c(TRUE, FALSE), c('#ffe6e6', 'white'))
      )
  })
  
  
  # Cross Summary table
  output$cross_summary_table <- renderDT({
    req(values$cross_summary)
    
    display_df <- values$cross_summary %>%
      mutate(
        Mean_LOD = round(Mean_LOD, 3),
        Max_LOD = round(Max_LOD, 3),
        Mean_PVE = round(Mean_PVE, 2)
      )
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        order = list(list(6, 'desc'))
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "Cross Summary Statistics",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'Max_LOD',
        background = styleColorBar(range(display_df$Max_LOD, na.rm = TRUE), 'lightgreen'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  
  
  # SNP summary table
  output$crosses_snp_combined_table <- renderDT({
    req(values$top_crosses_snp_combined)
    
    tryCatch({
      display_df <- values$top_crosses_snp_combined %>%
        mutate(
          Mean_Effect = round(Mean_Effect, 4)
        )
      
      datatable(
        display_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = "SNP Summary Across Crosses"
      )
    }, error = function(e) {
      datatable(data.frame(Message = "Error loading SNP summary"))
    })
  })
  
  
  # Top crosses table
  output$top_three_crosses_table <- renderDT({
    req(values$cross_summary)
    
    tryCatch({
      top_n <- input$cross_snp_top_n
      
      top_crosses <- values$cross_summary %>%
        arrange(desc(Max_LOD)) %>%
        head(top_n)
      
      datatable(
        top_crosses,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        ),
        class = 'display compact',
        rownames = FALSE,
        caption = paste("Top", top_n, "Crosses by Maximum LOD Score")
      )
    }, error = function(e) {
      datatable(data.frame(Message = "Error loading top crosses"))
    })
  })
  
  
  
  # Cross-SNP summary
  output$cross_snp_summary <- renderPrint({
    req(values$cross_snp_analysis)
    
    analysis <- values$cross_snp_analysis
    
    cat("=== CROSS-SNP ANALYSIS SUMMARY ===\n\n")
    cat("Analysis timestamp:", format(analysis$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Total cross-SNP combinations:", nrow(analysis$data), "\n")
    cat("Unique crosses:", analysis$n_crosses, "\n")
    cat("Unique SNPs:", analysis$n_snps, "\n")
    cat("Traits analyzed:", analysis$n_traits, "\n")
    
    if (!is.null(values$cross_summary)) {
      cat("\n=== CROSS STATISTICS ===\n")
      cat("Crosses with data:", nrow(values$cross_summary), "\n")
      cat("Average SNPs per cross:", round(mean(values$cross_summary$N_SNPs), 1), "\n")
      cat("Cross with most SNPs:", 
          values$cross_summary$Cross[which.max(values$cross_summary$N_SNPs)], "\n")
    }
  })
  
  
  
  # ==========================================================================
  # DONLOAD HANDLERS FOR CROSS-SNP ANALYSIS
  # ==========================================================================
 
  # DL cross SNP detailed table
  output$dl_cross_snp_detailed_table <- renderDT({
    req(values$dl_cross_snp_detailed)
    
    display_df <- values$dl_cross_snp_detailed %>%
      mutate(
        P_value = format(P_value, scientific = TRUE, digits = 3),
        LOD = round(LOD, 3),
        Effect = round(Effect, 4),
        PVE = round(PVE, 2)
      ) %>%
      select(Cross, Parent1, Parent2, SNP_ID, Chromosome, Position, Trait,
             P_value, LOD, Effect, PVE, Significant)
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf'),
        order = list(list(8, 'desc'))
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "DL Cross-SNP Detailed Data",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'Significant',
        target = 'row',
        backgroundColor = styleEqual(c(TRUE, FALSE), c('#ffe6e6', 'white'))
      )
  })
  
  # DL cross summary table
  output$dl_cross_summary_table <- renderDT({
    req(values$dl_cross_summary)
    
    display_df <- values$dl_cross_summary %>%
      mutate(
        Mean_LOD = round(Mean_LOD, 3),
        Max_LOD = round(Max_LOD, 3),
        Mean_PVE = round(Mean_PVE, 2)
      )
    
    datatable(
      display_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        order = list(list(6, 'desc'))
      ),
      class = 'display compact',
      rownames = FALSE,
      caption = "DL Cross Summary Statistics",
      extensions = 'Buttons',
      filter = 'top'
    ) %>%
      formatStyle(
        'Max_LOD',
        background = styleColorBar(range(display_df$Max_LOD, na.rm = TRUE), 'lightgreen'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  
  
  # Download cross-SNP results
  output$download_cross_snp_data <- downloadHandler(
    filename = function() {
      paste0("cross_snp_analysis_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(values$cross_snp_analysis)
      
      # Create temporary directory
      temp_dir <- tempfile("cross_snp_")
      dir.create(temp_dir)
      
      # Save all dataframes
      write.csv(values$cross_snp_analysis$data,
                file.path(temp_dir, "cross_snp_detailed.csv"),
                row.names = FALSE)
      
      if (!is.null(values$cross_summary)) {
        write.csv(values$cross_summary,
                  file.path(temp_dir, "cross_summary.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$top_crosses_snp_combined)) {
        write.csv(values$top_crosses_snp_combined,
                  file.path(temp_dir, "top_crosses_snp_combined.csv"),
                  row.names = FALSE)
      }
      
      # Save unified multi-trait data if available
      if (!is.null(values$multi_trait_unified_df)) {
        write.csv(values$multi_trait_unified_df,
                  file.path(temp_dir, "multi_trait_unified.csv"),
                  row.names = FALSE)
      }
      
      if (!is.null(values$multi_trait_enhanced_unified_df)) {
        write.csv(values$multi_trait_enhanced_unified_df,
                  file.path(temp_dir, "multi_trait_enhanced.csv"),
                  row.names = FALSE)
      }
      
      # Create metadata
      meta_data <- data.frame(
        Parameter = c("Analysis Date", "Total Cross-SNP Rows", "Unique Crosses", 
                      "Unique SNPs", "Traits Analyzed", "Polymorphic SNPs", "Data Source"),
        Value = c(
          as.character(Sys.time()),
          nrow(values$cross_snp_analysis$data),
          values$cross_snp_analysis$n_crosses,
          values$cross_snp_analysis$n_snps,
          values$cross_snp_analysis$n_traits,
          sum(values$cross_snp_analysis$data$Polymorphic == TRUE, na.rm = TRUE),
          "Merged Diallel and Multi-Trait GWAS Data"
        )
      )
      write.csv(meta_data, file.path(temp_dir, "metadata.csv"), row.names = FALSE)
      
      # Create README
      readme_content <- paste(
        "CROSS-SNP ANALYSIS RESULTS",
        "==========================",
        "",
        "Files included:",
        "1. cross_snp_detailed.csv - Detailed cross-SNP combinations",
        "2. cross_summary.csv - Summary statistics per cross (polymorphic SNPs only)",
        "3. top_crosses_snp_combined.csv - SNP statistics across crosses (polymorphic SNPs only)",
        "4. multi_trait_unified.csv - Unified multi-trait GWAS results",
        "5. multi_trait_enhanced.csv - Enhanced multi-trait GWAS results",
        "6. metadata.csv - Analysis parameters",
        "",
        "Column descriptions for cross_snp_detailed.csv:",
        "- Cross: Diallel cross ID",
        "- Parent1, Parent2: Parental lines in the cross",
        "- SNP_ID: SNP identifier",
        "- Chromosome, Position: Genomic coordinates",
        "- Trait: Phenotypic trait",
        "- P_value: GWAS p-value",
        "- LOD: -log10(P_value)",
        "- PVE: Phenotypic variance explained (%)",
        "- Effect: Additive effect size",
        "- Parent1_Genotype, Parent2_Genotype: Genotype values (if available)",
        "- Expected_Genotype: Mid-parent genotype value",
        "- Polymorphic: TRUE if parents have different genotypes at this SNP",
        "- Significant: TRUE if p < 0.05",
        "- Effect_Direction: Positive or Negative effect",
        "",
        "Note: Statistics in cross_summary.csv and top_crosses_snp_combined.csv",
        "are calculated only for polymorphic SNPs (Polymorphic = TRUE)",
        "This ensures cross-specific statistics that vary between crosses.",
        "",
        sep = "\n"
      )
      writeLines(readme_content, file.path(temp_dir, "README.txt"))
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, "results.zip")
      zip::zip(zip_file, files = list.files(temp_dir, full.names = TRUE), 
               mode = "cherry-pick")
      
      # Copy to download location
      file.copy(zip_file, file)
      
      # Cleanup
      unlink(temp_dir, recursive = TRUE)
    }
  )
  
  
  
  # Clear cross-SNP results
  observeEvent(input$clear_cross_snp_results, {
    values$cross_snp_analysis <- NULL
    values$cross_summary <- NULL
    values$top_crosses_snp_combined <- NULL
    
    # Reset outputs
    output$cross_snp_table <- renderDT({NULL})
    output$cross_summary_table <- renderDT({NULL})
    output$crosses_snp_combined_table <- renderDT({NULL})
    output$cross_snp_summary <- renderPrint({
      cat("=== CROSS-SNP ANALYSIS ===\n\n")
      cat("No results available.\n")
      cat("Click 'Run Cross-SNP Analysis' to analyze.\n")
    })
    
    output$cross_lod_plot <- renderPlot({
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No cross-SNP results available", size = 6) +
        theme_void()
    })
    
    output$snp_effect_heatmap <- renderPlot({
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No cross-SNP results available", size = 6) +
        theme_void()
    })
    
    showNotification("Cross-SNP results cleared!", type = "warning", duration = 3)
  })
  
  
  # Cross-SNP status UI
  output$cross_snp_status_ui <- renderUI({
    if (!is.null(values$cross_snp_analysis)) {
      analysis <- values$cross_snp_analysis
      
      div(class = "alert alert-success",
          icon("check-circle"),
          strong(" ✓ Cross-SNP Analysis Complete"),
          br(),
          paste("Crosses analyzed:", analysis$n_crosses),
          br(),
          paste("SNPs analyzed:", analysis$n_snps),
          br(),
          paste("Traits:", analysis$n_traits)
      )
    } else {
      div(class = "alert alert-info",
          icon("info-circle"),
          strong(" Ready for Cross-SNP Analysis"),
          br(),
          "Configure settings and click 'Run Cross-SNP Analysis'"
      )
    }
  })
  
  
  # Unified data status UI
  output$unified_data_status_ui <- renderUI({
    if (!is.null(values$multi_trait_unified_df)) {
      df <- values$multi_trait_unified_df
      div(
        class = "alert alert-success",
        icon("check-circle"),
        strong(" ✓ Unified Multi-Trait Data Ready"),
        br(),
        paste("Traits:", n_distinct(df$Trait)),
        br(),
        paste("SNP-trait combinations:", nrow(df)),
        br(),
        paste("Memory:", format(object.size(df), units = "MB"))
      )
    } else if (!is.null(values$multi_trait_results)) {
      div(
        class = "alert alert-info",
        icon("info-circle"),
        strong(" Multi-Trait Analysis Complete"),
        br(),
        "Click 'Create Unified Format' to organize data."
      )
    } else {
      div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        strong(" No Multi-Trait Data"),
        br(),
        "Run multi-trait analysis first."
      )
    }
  })
  
  
  
  # Data compatibility check
  output$data_compatibility_check <- renderPrint({
    cat("=== DATA COMPATIBILITY CHECK ===\n\n")
    
    # Check required data
    if (!is.null(values$multi_trait_unified_df)) {
      cat("Multi-trait GWAS data: ✓ Available\n")
      cat("  SNPs:", n_distinct(values$multi_trait_unified_df$SNP_ID), "\n")
      cat("  Traits:", n_distinct(values$multi_trait_unified_df$Trait), "\n")
    } else {
      cat("Multi-trait GWAS data: ✗ Missing\n")
    }
    
    if (!is.null(values$diallel_data)) {
      cat("Diallel data: ✓ Available\n")
      cat("  Crosses:", nrow(values$diallel_data[values$diallel_data$CrossType == "Cross", ]), "\n")
    } else {
      cat("Diallel data: ✗ Missing\n")
    }
    
    if (!is.null(values$matched_data$geno)) {
      cat("Genotype data: ✓ Available\n")
    } else {
      cat("Genotype data: ✗ Missing (optional)\n")
    }
  })
  
  
  #===========================================================================
  #ALL TABLES DISPLAY
  #===========================================================================
  # Corrected server code for the dataframes viewer
  
  # Track all dataframes created in the app - FIXED VERSION
  observe({
    # Safely get reactive values
    tryCatch({
      # List all reactive values that are dataframes
      df_list <- reactiveValuesToList(values)
      
      # Safely check each element
      is_df_or_matrix <- function(x) {
        tryCatch({
          is.data.frame(x) || is.matrix(x) || inherits(x, "data.frame")
        }, error = function(e) FALSE)
      }
      
      df_names <- names(df_list)[sapply(df_list, is_df_or_matrix)]
      
      # Store the list of dataframe names
      if (length(df_names) > 0) {
        values$all_dataframes <- df_names
      } else {
        values$all_dataframes <- character(0)
      }
      
      # Update select input safely
      if (length(df_names) > 0) {
        updateSelectInput(session, "selected_dataframe", 
                          choices = df_names,
                          selected = df_names[1])
      } else {
        updateSelectInput(session, "selected_dataframe", 
                          choices = "No dataframes available",
                          selected = "No dataframes available")
      }
    }, error = function(e) {
      # Handle error gracefully
      cat("Error in dataframe observer:", e$message, "\n")
      values$all_dataframes <- character(0)
      updateSelectInput(session, "selected_dataframe", 
                        choices = "Error loading dataframes",
                        selected = "Error loading dataframes")
    })
  })
  
  # Table of available dataframes - FIXED VERSION
  output$available_dataframes_table <- renderDT({
    req(values$all_dataframes)
    
    # Check if we have any dataframes
    if (length(values$all_dataframes) == 0 || 
        identical(values$all_dataframes, "No dataframes available") ||
        identical(values$all_dataframes, "Error loading dataframes")) {
      return(datatable(data.frame(Message = "No dataframes available"), 
                       options = list(pageLength = 5)))
    }
    
    # Create dataframe info safely
    df_info <- tryCatch({
      data.frame(
        DataFrame = values$all_dataframes,
        Rows = sapply(values$all_dataframes, function(x) {
          df <- values[[x]]
          if (!is.null(df) && (is.data.frame(df) || is.matrix(df))) {
            nrow(df)
          } else {
            NA
          }
        }),
        Columns = sapply(values$all_dataframes, function(x) {
          df <- values[[x]]
          if (!is.null(df) && (is.data.frame(df) || is.matrix(df))) {
            ncol(df)
          } else {
            NA
          }
        }),
        Type = sapply(values$all_dataframes, function(x) {
          df <- values[[x]]
          if (is.null(df)) {
            "NULL"
          } else if (is.data.frame(df)) {
            "DataFrame"
          } else if (is.matrix(df)) {
            "Matrix"
          } else if (is.list(df)) {
            "List"
          } else {
            class(df)[1]
          }
        }),
        Memory = sapply(values$all_dataframes, function(x) {
          df <- values[[x]]
          if (!is.null(df) && (is.data.frame(df) || is.matrix(df))) {
            format(object.size(df), units = "KB")
          } else {
            "N/A"
          }
        }),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(Error = paste("Error creating table:", e$message))
    })
    
    datatable(df_info, 
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE,
              caption = "Available DataFrames in the Application")
  })
  
  # Dataframe preview - FIXED VERSION
  output$dataframe_preview <- renderDT({
    # Check if a valid dataframe is selected
    if (is.null(input$selected_dataframe) || 
        input$selected_dataframe == "No dataframes available" ||
        input$selected_dataframe == "Error loading dataframes") {
      return(datatable(data.frame(Message = "Please select a valid dataframe")))
    }
    
    df <- tryCatch({
      values[[input$selected_dataframe]]
    }, error = function(e) NULL)
    
    if (is.null(df)) {
      return(datatable(data.frame(Message = "Dataframe not found or is NULL")))
    }
    
    if (is.data.frame(df) || is.matrix(df)) {
      # Show only first 50 rows for performance
      datatable(head(df, 50), 
                options = list(pageLength = 10, scrollX = TRUE),
                caption = paste("Preview of", input$selected_dataframe))
    } else {
      datatable(data.frame(Message = paste("Selected object is not a dataframe. Type:", class(df))))
    }
  })
  
  # Dataframe structure - FIXED VERSION
  output$dataframe_structure <- renderPrint({
    req(input$selected_dataframe)
    
    if (input$selected_dataframe %in% c("No dataframes available", "Error loading dataframes")) {
      cat("No dataframe selected\n")
      return()
    }
    
    df <- values[[input$selected_dataframe]]
    
    if (is.null(df)) {
      cat("Dataframe is NULL\n")
    } else if (is.data.frame(df)) {
      str(df)
    } else if (is.matrix(df)) {
      cat("Matrix dimensions:", dim(df), "\n")
      cat("Matrix class:", class(df), "\n")
      cat("Row names:", if(!is.null(rownames(df))) paste(length(rownames(df)), "rows") else "None", "\n")
      cat("Column names:", if(!is.null(colnames(df))) paste(length(colnames(df)), "columns") else "None", "\n")
      cat("\nFirst few rows and columns:\n")
      print(head(df[, 1:min(5, ncol(df))], 5))
    } else {
      cat("Object type:", class(df), "\n")
      cat("Length:", length(df), "\n")
      cat("Structure:\n")
      str(df)
    }
  })
  
  # Dataframe summary - FIXED VERSION
  output$dataframe_summary <- renderPrint({
    req(input$selected_dataframe)
    
    if (input$selected_dataframe %in% c("No dataframes available", "Error loading dataframes")) {
      cat("No dataframe selected\n")
      return()
    }
    
    df <- values[[input$selected_dataframe]]
    
    if (is.null(df)) {
      cat("Dataframe is NULL\n")
    } else if (is.data.frame(df)) {
      summary(df)
    } else if (is.matrix(df)) {
      cat("Matrix Summary:\n")
      cat("Dimensions:", dim(df), "\n")
      cat("Row names:", if(!is.null(rownames(df))) length(rownames(df)) else "None", "\n")
      cat("Column names:", if(!is.null(colnames(df))) length(colnames(df)) else "None", "\n")
      cat("NA values:", sum(is.na(df)), "\n")
      cat("Numeric summary:\n")
      if(is.numeric(df)) {
        print(summary(as.vector(df)))
      }
    } else {
      cat("Object summary not available for type:", class(df), "\n")
    }
  })
  
  # Dataframe creation function info - FIXED VERSION
  output$dataframe_function <- renderPrint({
    req(input$selected_dataframe)
    
    if (input$selected_dataframe %in% c("No dataframes available", "Error loading dataframes")) {
      cat("No dataframe selected\n")
      return()
    }
    
    df_name <- input$selected_dataframe
    
    # Map dataframe names to their creation functions
    function_map <- list(
      # GWAS dataframes
      "gwas_results" = "run_single_trait_gwas()",
      "combined_results" = "combine_multi_trait_results()",
      "multi_trait_wide" = "create_wide_gwas_format()",
      "multi_trait_enhanced" = "create_enhanced_unified_df()",
      "complete_integrated_df" = "create_complete_integrated_df()",
      "gwas_summary" = "capture_output_results()",
      "merged_pvalues" = "combine_p_values()",
      
      # Diallel dataframes
      "diallel_data" = "load_diallel_data()",
      "diallel_results" = "perform_diallel_analysis()",
      "diallel_heterosis_data" = "calculate_heterosis()",
      "gca_df" = "perform_diallel_analysis()",
      "sca_df" = "sca_matrix_to_df()",
      "diallel_summary_tables" = "run_diallel_anova_all_traits()",
      "cross_summary" = "run_cross_snp_analysis_all_traits_app()",
      
      # Cross-SNP dataframes
      "cross_snp_df" = "create_cross_snp_dataframe_all_traits_fixed()",
      "cross_snp_simple" = "create_cross_snp_dataframe_simple()",
      "top_crosses_snp_combined" = "run_cross_snp_analysis_all_traits_app()",
      "trait_summary" = "run_cross_snp_analysis_all_traits_app()",
      "cross_summary_with_flanking" = "add_flanking_snp_markers_safe()",
      
      # METAN dataframes
      "met_data" = "prepare_metan_data()",
      "met_analysis_results" = "run_metan_analysis()",
      "desc_stats" = "run_metan_analysis()",
      "env_summary" = "characterize_environments()",
      "stability_df" = "run_stability_analysis()",
      "missing_summary" = "summarize_missing_data()",
      "ammi_results" = "run_ammi_analysis()",
      "gge_results" = "perform_gge_biplot()",
      "fw_results" = "perform_finlay_wilkinson()",
      
      # Basic dataframes
      "pheno_data" = "load_pheno_data()",
      "metadata" = "load_meta_data()",
      "geno_matrix" = "genlight_to_matrix_stats()",
      "snp_info" = "prepare_geno_data()"
    )
    
    if(df_name %in% names(function_map)) {
      cat("Creation Function:", function_map[[df_name]], "\n\n")
      cat("Function Purpose:\n")
      
      # Add descriptions for key functions
      desc_map <- list(
        "run_single_trait_gwas()" = "Performs GWAS for individual traits using linear models",
        "combine_multi_trait_results()" = "Combines GWAS results from multiple traits",
        "create_wide_gwas_format()" = "Creates wide format GWAS results for multiple traits",
        "perform_diallel_analysis()" = "Analyzes diallel crosses for GCA and SCA effects",
        "calculate_heterosis()" = "Calculates mid-parent and better-parent heterosis",
        "create_cross_snp_dataframe_all_traits_fixed()" = "Creates comprehensive cross-SNP dataframe for all traits",
        "prepare_metan_data()" = "Prepares data for METAN analysis with environment grouping",
        "run_metan_analysis()" = "Performs comprehensive METAN analysis including AMMI, GGE, and stability"
      )
      
      func_name <- function_map[[df_name]]
      if(func_name %in% names(desc_map)) {
        cat(desc_map[[func_name]], "\n")
      } else {
        cat("No description available.\n")
      }
    } else {
      cat("Creation function information not available.\n")
      cat("This dataframe was likely created by one of the analysis functions.\n")
    }
  })
  
  # Download handler for dataframes - FIXED VERSION
  output$download_dataframe <- downloadHandler(
    filename = function() {
      if (input$selected_dataframe %in% c("No dataframes available", "Error loading dataframes")) {
        "no_data.csv"
      } else {
        paste0(input$selected_dataframe, "_", Sys.Date(), ".csv")
      }
    },
    content = function(file) {
      if (input$selected_dataframe %in% c("No dataframes available", "Error loading dataframes")) {
        write.csv(data.frame(Message = "No dataframe available for download"), file)
        return()
      }
      
      df <- values[[input$selected_dataframe]]
      if (is.null(df)) {
        write.csv(data.frame(Message = "Selected dataframe is NULL"), file)
      } else if (is.data.frame(df)) {
        write.csv(df, file, row.names = FALSE)
      } else if (is.matrix(df)) {
        write.csv(as.data.frame(df), file, row.names = TRUE)
      } else {
        write.csv(data.frame(Message = paste("Cannot download object of type", class(df))), file)
      }
    }
  )
  
  # Refresh dataframes list - FIXED VERSION
  observeEvent(input$refresh_dataframes, {
    showNotification("Refreshing dataframe list...", type = "message")
    
    # Force update by triggering the observer
    session$sendCustomMessage("refresh", "dataframes")
  })
  
  # Add a simple status indicator
  output$unified_data_status_ui <- renderUI({
    if (is.null(values$all_dataframes) || length(values$all_dataframes) == 0) {
      tags$div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        " No dataframes available. Please load data first."
      )
    } else {
      tags$div(
        class = "alert alert-success",
        icon("check-circle"),
        paste(" Found", length(values$all_dataframes), "dataframes")
      )
    }
  })

  # ============================================================
  # POPULATION GENETICS SERVER OBSERVERS
  # ============================================================
  
  # Reactive value for population genetics data
  popgen_data <- reactiveValues(
    allele_freqs = NULL,
    correlations = NULL,
    waples_results = NULL
  )
  
  # UI for population selector
  output$pop_selector_ui <- renderUI({
    if (is.null(values$processed_geno) || is.null(values$metadata)) {
      return(helpText("Load genotype and metadata data first"))
    }
    
    # Extract possible population columns from metadata
    meta_cols <- colnames(values$metadata)
    pop_cols <- meta_cols[sapply(values$metadata, function(x) length(unique(x)) < 20)]
    
    selectInput("pop_column", "Select Population Column",
                choices = c("No grouping" = "none", pop_cols),
                selected = "none")
  })
  
  # Observe allele frequency analysis
  observeEvent(input$run_allele_freq_analysis, {
    req(values$processed_geno)
    
    showNotification("Calculating allele frequencies...", type = "message")
    
    # Get population assignments
    pop_assignments <- NULL
    if (!is.null(input$pop_column) && input$pop_column != "none") {
      req(values$metadata)
      # Match samples to populations
      sample_names <- values$processed_geno$sample_names
      meta_samples <- values$metadata[,1]  # Assuming first column is sample ID
      
      # Simple matching (you may need to adjust this based on your data structure)
      pop_assignments <- values$metadata[[input$pop_column]][
        match(sample_names, meta_samples)
      ]
    }
    
    # Calculate allele frequencies
    popgen_data$allele_freqs <- calculate_allele_frequencies(
      values$gl_object,  # Using the raw genlight object
      pop_assignments
    )
    
    showNotification("Allele frequencies calculated!", type = "message")
  })
  
  # Output allele frequency histogram
  output$allele_freq_hist_plot <- renderPlot({
    req(popgen_data$allele_freqs)
    
    # Get populations from column names
    freq_cols <- grep("^allelefreq_", colnames(popgen_data$allele_freqs), value = TRUE)
    populations <- gsub("^allelefreq_", "", freq_cols)
    
    allelefreqhisto(popgen_data$allele_freqs, populations)
  })
  
  # Observe correlation analysis
  observeEvent(input$run_corr_analysis, {
    req(popgen_data$allele_freqs)
    
    showNotification("Calculating population correlations...", type = "message")
    
    # Get populations from column names
    freq_cols <- grep("^allelefreq_", colnames(popgen_data$allele_freqs), value = TRUE)
    populations <- gsub("^allelefreq_", "", freq_cols)
    
    popgen_data$correlations <- corrfrequencies(popgen_data$allele_freqs, populations)
    
    showNotification("Correlations calculated!", type = "message")
  })
  
  # Output correlation heatmap
  output$corr_heatmap <- renderPlot({
    req(popgen_data$correlations)
    
    corr_matrix <- popgen_data$correlations$matrix
    
    # Create heatmap
    corrplot(corr_matrix, 
             method = "color", 
             type = "upper",
             addCoef.col = "black",
             tl.col = "black",
             tl.srt = 45,
             title = "Population Allele Frequency Correlations",
             mar = c(0, 0, 2, 0))
  })
  
  # Output correlation summary
  output$corr_summary <- renderPrint({
    req(popgen_data$correlations)
    
    cat("Population Correlation Analysis\n")
    cat("===============================\n\n")
    cat("Number of populations:", nrow(popgen_data$correlations$matrix), "\n")
    cat("Number of comparisons:", nrow(popgen_data$correlations$table), "\n\n")
    
    # Summary statistics
    pearson_vals <- popgen_data$correlations$table$Pearson
    spearman_vals <- popgen_data$correlations$table$Spearman
    
    cat("Pearson Correlation:\n")
    cat("  Mean:", round(mean(pearson_vals, na.rm = TRUE), 4), "\n")
    cat("  Range:", round(range(pearson_vals, na.rm = TRUE), 4), "\n\n")
    
    cat("Spearman Correlation:\n")
    cat("  Mean:", round(mean(spearman_vals, na.rm = TRUE), 4), "\n")
    cat("  Range:", round(range(spearman_vals, na.rm = TRUE), 4), "\n")
  })
  
  # Output correlation table
  output$corr_table <- renderDT({
    req(popgen_data$correlations)
    
    datatable(popgen_data$correlations$table,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = c("Pearson", "Spearman"), digits = 4)
  })
  
  # Observe Waples analysis
  observeEvent(input$run_waples_analysis, {
    req(values$processed_geno)
    
    showNotification("Running Fst-Fis analysis...", type = "message")
    
    # For this example, we need Fst and Fis data
    # In a real scenario, you would calculate these from your data
    # Here we'll create dummy data for demonstration
    
    n_snps <- nrow(values$processed_geno$snp_info)
    dummy_data <- data.frame(
      SNP_ID = values$processed_geno$snp_info$SNP_ID,
      Fst = runif(n_snps, 0, 0.5),
      Fis = runif(n_snps, -0.2, 0.3)
    )
    
    popgen_data$waples_results <- waples_snps(
      dummy_data, 
      fis_column = "Fis", 
      fst_column = "Fst",
      mysteps = input$fst_fis_bin_size
    )
    
    showNotification("Waples plot generated!", type = "message")
  })
  
  # Output Waples plot
  output$waples_plot <- renderPlot({
    req(popgen_data$waples_results)
    
    popgen_data$waples_results$plot
  })
  
  # Output Waples statistics table
  output$waples_stats_table <- renderDT({
    req(popgen_data$waples_results)
    
    datatable(popgen_data$waples_results$binned_data,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE,
              caption = "Binned Fst-Fis Statistics") %>%
      formatRound(columns = c("Fst_midpoint", "Mean_Fis"), digits = 4)
  })
  
  
  # Download allele frequency plot
  output$download_allele_freq_plot <- downloadHandler(
    filename = function() {
      paste0("allele_frequency_histogram_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = output$allele_freq_hist_plot(), 
             width = 12, height = 10, dpi = 300)
    }
  )
  
  # Download correlation matrix
  output$download_corr_matrix <- downloadHandler(
    filename = function() {
      paste0("population_correlations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(popgen_data$correlations$table, file, row.names = FALSE)
    }
  )
  
  # Download Waples plot
  output$download_waples_plot <- downloadHandler(
    filename = function() {
      paste0("waples_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = popgen_data$waples_results$plot, 
             width = 10, height = 8, dpi = 300)
    }
  )
  
  #=============================================================================
  #POP GEN
  #=============================================================================
  
  popgen_values <- reactiveValues(
    results = NULL,
    genind_obj = NULL,
    available_pops = NULL
  )
  
  
  # Run Population Genetics Analyses
  observeEvent(input$run_popgen_analyses, {
    req(values$filtered_gl)
    
    withProgress(message = 'Running Population Genetics Analyses', value = 0, {
      incProgress(0.2, detail = "This may take several minutes depending on dataset size...")
      
      tryCatch({
        # Get population data
        incProgress(0.4, detail = "Preparing population data...")
        pop_data <- NULL
        if (!is.null(values$metadata)) {
          # Try to find population column
          pop_cols <- c("Population", "POP", "pop", "Group", "GROUP")
          for (col in pop_cols) {
            if (col %in% colnames(values$metadata)) {
              pop_data <- values$metadata[[col]]
              break
            }
          }
        }
        
        # Run analyses
        incProgress(0.5, detail = "Running population genetics analysis...")
        analysis_results <- run_popgen_analyses(
          gl_object = values$filtered_gl,
          pop_data = pop_data,
          analyses = input$popgen_analyses
        )
        
        if (is.null(analysis_results)) {
          stop("Population genetics analyses failed")
        }
        
        # Store results
        popgen_values$results <- analysis_results
        popgen_values$available_pops <- analysis_results$metadata$populations
        
        # Update hybrid detection dropdowns
        updateSelectInput(session, "hybrid_pop1", 
                          choices = popgen_values$available_pops,
                          selected = if(length(popgen_values$available_pops) > 0) popgen_values$available_pops[1] else NULL)
        updateSelectInput(session, "hybrid_pop2", 
                          choices = popgen_values$available_pops,
                          selected = if(length(popgen_values$available_pops) > 1) popgen_values$available_pops[2] else NULL)
        updateSelectInput(session, "hybrid_admixed", 
                          choices = popgen_values$available_pops,
                          selected = if(length(popgen_values$available_pops) > 2) popgen_values$available_pops[3] else NULL)
        
        # Check which analyses to run
        if ("dapc" %in% input$popgen_analyses) {
          incProgress(0.85, detail = "Running DAPC analysis...")
          dapc_values$results <- perform_dapc_analysis_comprehensive(
            gl_object = values$filtered_gl,
            pop_assignments = values$metadata$Population,
            find_clusters = input$popgen_find_clusters,
            max_n_clust = input$popgen_max_clusters
          )
        }
        incProgress(1, detail = "Population genetics analysis completed successfully...")
        showNotification("Population genetics analyses completed successfully!", 
                         type = "message", duration = 5)
        
      }, error = function(e) {
        showNotification(paste("Error running population genetics analyse:", e$message), type = "error")
      })
    })
  })
  
  
  
  # Status Output
  output$popgen_status <- renderPrint({
    if (is.null(popgen_values$results)) {
      cat("Population genetics analyses not yet run.\n")
      cat("Please run the analyses using the button above.\n")
      cat("\nPrerequisites:\n")
      cat("  1. Genotype data must be loaded and filtered\n")
      cat("  2. Population information should be available\n")
      cat("  3. Select analyses to perform\n")
    } else {
      meta <- popgen_values$results$metadata
      cat("ANALYSIS STATUS: COMPLETE\n")
      cat("==========================\n\n")
      cat("Dataset Information:\n")
      cat("  Individuals:", meta$n_individuals, "\n")
      cat("  Loci:", meta$n_loci, "\n")
      cat("  Populations:", meta$n_populations, "\n")
      cat("  Analysis date:", format(meta$analysis_date, "%Y-%m-%d %H:%M"), "\n")
      
      cat("\nAnalyses Performed:\n")
      for (analysis in meta$analyses_performed) {
        cat("  ✓", switch(analysis,
                          "allele_barcode" = "Allele Frequency Barcode",
                          "heterozygosity" = "Heterozygosity Statistics",
                          "fst" = "Pairwise FST",
                          "hwe" = "Hardy-Weinberg Equilibrium",
                          "dapc" = "DAPC Assignment",
                          analysis), "\n")
      }
      
      # Check for specific results
      if (!is.null(popgen_values$results$fst)) {
        cat("\nFST Statistics:\n")
        fst_mat <- popgen_values$results$fst$fst_matrix
        if (!is.null(fst_mat)) {
          cat("  Mean FST:", round(mean(fst_mat, na.rm = TRUE), 4), "\n")
          cat("  Range:", round(range(fst_mat, na.rm = TRUE), 4), "\n")
        }
      }
      
      if (!is.null(popgen_values$results$heterozygosity)) {
        cat("\nHeterozygosity Statistics:\n")
        het_sum <- popgen_values$results$heterozygosity$summary_table
        if (!is.null(het_sum)) {
          cat("  Mean Ho:", round(mean(het_sum$Hobs_mean, na.rm = TRUE), 3), "\n")
          cat("  Mean He:", round(mean(het_sum$Hexp_mean, na.rm = TRUE), 3), "\n")
        }
      }
    }
  })
  
  # Quick Statistics UI
  output$popgen_quick_stats <- renderUI({
    req(popgen_values$results)
    
    meta <- popgen_values$results$metadata
    
    tagList(
      h4("Dataset Summary"),
      tags$ul(
        tags$li(strong("Individuals:"), meta$n_individuals),
        tags$li(strong("Markers:"), meta$n_loci),
        tags$li(strong("Populations:"), meta$n_populations)
      ),
      
      h4("Available Populations"),
      tags$div(
        style = "max-height: 200px; overflow-y: auto;",
        tags$ul(
          lapply(meta$populations, function(pop) tags$li(pop))
        )
      )
    )
  })
  
  # Population Information
  output$popgen_pop_info <- renderPrint({
    req(values$filtered_gl)
    
    if (!is.null(pop(values$filtered_gl))) {
      pop_table <- table(pop(values$filtered_gl))
      
      cat("Population Distribution:\n")
      cat("=======================\n\n")
      
      pop_df <- data.frame(
        Population = names(pop_table),
        N_Individuals = as.numeric(pop_table),
        Proportion = round(as.numeric(pop_table) / nInd(values$filtered_gl) * 100, 1)
      )
      
      pop_df <- pop_df[order(pop_df$N_Individuals, decreasing = TRUE), ]
      
      for (i in 1:nrow(pop_df)) {
        cat(sprintf("%-20s %5d individuals (%5.1f%%)\n", 
                    pop_df$Population[i], 
                    pop_df$N_Individuals[i],
                    pop_df$Proportion[i]))
      }
      
      cat("\nSummary:\n")
      cat("  Total:", sum(pop_df$N_Individuals), "\n")
      cat("  Mean per pop:", round(mean(pop_df$N_Individuals), 1), "\n")
      cat("  Min:", min(pop_df$N_Individuals), "\n")
      cat("  Max:", max(pop_df$N_Individuals), "\n")
    } else {
      cat("No population assignments found.\n")
      cat("All individuals are unassigned.\n")
      cat("Consider providing population information in metadata.\n")
    }
  })
  
  # Allele Barcode Plot
  output$popgen_barcode_plot <- renderPlot({
    req(popgen_values$results$allele_barcode)
    
    barcode_results <- popgen_values$results$allele_barcode
    
    # Create heatmap
    heatmap_plot <- plot_allele_heatmap(
      barcode_results,
      color_palette = input$barcode_color,
      show_dendrogram = TRUE,
      cluster_rows = TRUE,
      cluster_cols = input$barcode_cluster
    )
    
    return(heatmap_plot)
  })
  
  # Heterozygosity Plot
  output$popgen_het_plot <- renderPlot({
    req(popgen_values$results$heterozygosity$comparison_plot)
    popgen_values$results$heterozygosity$comparison_plot
  })
  
  output$popgen_het_table <- renderDT({
    req(popgen_values$results$heterozygosity$summary_table)
    
    het_table <- popgen_values$results$heterozygosity$summary_table
    
    datatable(het_table,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = c("Hobs_mean", "Hexp_mean", "Fis_mean"), digits = 4)
  })
  
  # FST Analysis
  output$popgen_fst_plot <- renderPlot({
    req(popgen_values$results$fst$fst_plot)
    popgen_values$results$fst$fst_plot
  })
  
  output$popgen_fst_table <- renderDT({
    req(popgen_values$results$fst$fst_matrix)
    
    fst_mat <- popgen_values$results$fst$fst_matrix
    
    datatable(as.data.frame(fst_mat),
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Pairwise FST Matrix") %>%
      formatRound(columns = colnames(fst_mat), digits = 4)
  })
  
  # NJ Tree availability
  output$fst_tree_available <- reactive({
    return(!is.null(popgen_values$results$fst$nj_tree))
  })
  
  outputOptions(output, "fst_tree_available", suspendWhenHidden = FALSE)
  
  output$popgen_nj_tree <- renderPlot({
    req(popgen_values$results$fst$nj_tree)
    
    tree <- popgen_values$results$fst$nj_tree
    
    # Plot NJ tree
    plot(tree, type = "unrooted", 
         tip.color = viridis::viridis(length(tree$tip.label)),
         main = "Neighbor-Joining Tree (FST distances)")
    add.scale.bar()
  })
  
  # HWE Analysis
  output$popgen_hwe_plot <- renderPlot({
    req(popgen_values$results$hwe$summary_plot)
    popgen_values$results$hwe$summary_plot
  })
  
  output$popgen_hwe_summary <- renderPrint({
    req(popgen_values$results$hwe$population_summary)
    
    hwe_sum <- popgen_values$results$hwe$population_summary
    
    cat("HWE Summary Statistics:\n")
    cat("=======================\n\n")
    
    for (i in 1:nrow(hwe_sum)) {
      cat(hwe_sum$Population[i], ":\n")
      cat("  Loci tested:", hwe_sum$N_Loci[i], "\n")
      cat("  Significant:", hwe_sum$N_Significant[i], "\n")
      cat("  Proportion:", round(hwe_sum$Proportion_Significant[i] * 100, 1), "%\n")
      cat("  Mean p-value:", round(hwe_sum$Mean_P_Value[i], 4), "\n\n")
    }
  })
  
  # DAPC Analysis
  output$popgen_dapc_scatter <- renderPlot({
    req(popgen_values$results$dapc$scatter_plot)
    popgen_values$results$dapc$scatter_plot
  })
  
  output$popgen_dapc_assignment <- renderPlot({
    req(popgen_values$results$dapc$assignment_plot)
    popgen_values$results$dapc$assignment_plot
  })
  
  # Hybrid Detection (simplified version)
  observeEvent(input$run_hybrid_detection, {
    req(popgen_values$results$allele_barcode,
        input$hybrid_pop1, input$hybrid_pop2, input$hybrid_admixed)
    
    tryCatch({
      # Simplified hybrid detection based on allele frequencies
      freq_mat <- popgen_values$results$allele_barcode$frequency_matrix
      
      pop1_freq <- freq_mat[, input$hybrid_pop1, drop = FALSE]
      pop2_freq <- freq_mat[, input$hybrid_pop2, drop = FALSE]
      admix_freq <- freq_mat[, input$hybrid_admixed, drop = FALSE]
      
      # Calculate simple hybrid indices
      hybrid_patterns <- data.frame(
        Marker = rownames(freq_mat),
        Pop1 = pop1_freq[, 1],
        Pop2 = pop2_freq[, 1],
        Admixed = admix_freq[, 1],
        Pattern = NA,
        stringsAsFactors = FALSE
      )
      
      # Define patterns
      hybrid_patterns$Pattern <- with(hybrid_patterns, {
        ifelse(Pop1 > 0.5 & Pop2 < 0.5 & Admixed > 0.5, "Pop1-like",
               ifelse(Pop1 < 0.5 & Pop2 > 0.5 & Admixed > 0.5, "Pop2-like",
                      ifelse(Pop1 > 0.5 & Pop2 > 0.5 & Admixed > 0.5, "Shared",
                             ifelse(Pop1 > 0.5 & Pop2 < 0.5 & Admixed < 0.5, "Lost",
                                    ifelse(Pop1 < 0.5 & Pop2 > 0.5 & Admixed < 0.5, "Gained", "Other")))))
      })
      
      # Store results
      popgen_values$hybrid_results <- hybrid_patterns
      
      showNotification("Hybrid analysis completed", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Hybrid analysis failed:", e$message), type = "error")
    })
  })
  
  output$hybrid_results_summary <- renderPrint({
    req(popgen_values$hybrid_results)
    
    patterns <- table(popgen_values$hybrid_results$Pattern)
    
    cat("Hybrid Pattern Summary:\n")
    cat("======================\n\n")
    
    for (pattern in names(patterns)) {
      cat(pattern, ":", patterns[pattern], "markers\n")
    }
    
    cat("\nTotal markers analyzed:", nrow(popgen_values$hybrid_results), "\n")
  })
  
  # Download Handlers
  output$download_popgen_report <- downloadHandler(
    filename = function() {
      paste0("population_genetics_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      # Generate R Markdown report
      temp_report <- file.path(tempdir(), "popgen_report.Rmd")
      file.copy("popgen_report_template.Rmd", temp_report, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
        results = popgen_values$results,
        analysis_date = Sys.Date()
      )
      
      # Knit the document
      rmarkdown::render(temp_report, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
    }
  )
  
  output$export_popgen_rdata <- downloadHandler(
    filename = function() {
      paste0("population_genetics_results_", Sys.Date(), ".RData")
    },
    content = function(file) {
      save(popgen_values$results, file = file)
    }
  )
  

  #=================================================================
  #POP ALLELE FREQUENCY
  #=================================================================
  
  # Population Allele Analysis Reactive Values
  pop_allele_values <- reactiveValues(
    results = NULL,
    available_pops = NULL,
    heterozygosity_calculated = FALSE
  )
  
  # Run Population Allele Analysis
  observeEvent(input$run_pop_allele_analysis, {
    req(values$filtered_gl)
    
    showModal(modalDialog(
      title = "Running Population Allele Analysis",
      "Calculating allele statistics for each population...",
      footer = NULL,
      size = "s"
    ))
    
    tryCatch({
      # Parse custom colors if provided
      custom_colors <- NULL
      if (input$pop_use_custom_colors && nchar(input$pop_custom_colors) > 0) {
        custom_colors <- trimws(unlist(strsplit(input$pop_custom_colors, ",")))
      }
      
      # Run analysis
      allele_results <- calculate_population_alleles(
        gl_object = values$filtered_gl,
        min_individuals = input$pop_min_inds
      )
      
      if (is.null(allele_results)) {
        stop("Failed to calculate population allele statistics")
      }
      
      # Store results
      pop_allele_values$results <- allele_results
      pop_allele_values$available_pops <- allele_results$statistics$Population
      pop_allele_values$heterozygosity_calculated <- "He_Mean" %in% colnames(allele_results$statistics)
      
      removeModal()
      showNotification("Population allele analysis completed successfully!", 
                       type = "message", duration = 5)
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Population allele analysis failed:", e$message), 
                       type = "error", duration = 10)
    })
  })
  
  # Check if heterozygosity is available for conditional panels
  output$heterozygosity_available <- reactive({
    return(!is.null(pop_allele_values$results) && 
             pop_allele_values$heterozygosity_calculated)
  })
  
  outputOptions(output, "heterozygosity_available", suspendWhenHidden = FALSE)
  
  # Main output functions
  output$pop_allele_summary <- renderPrint({
    if (is.null(pop_allele_values$results)) {
      cat("Population allele analysis not yet run.\n")
      cat("Please run the analysis using the button above.\n")
      cat("\nCurrent data status:\n")
      cat("  Genotype data:", ifelse(!is.null(values$filtered_gl), "Loaded", "Not loaded"), "\n")
      if (!is.null(values$filtered_gl)) {
        cat("  Number of individuals:", nInd(values$filtered_gl), "\n")
        cat("  Number of loci:", nLoc(values$filtered_gl), "\n")
        if (!is.null(pop(values$filtered_gl))) {
          pop_table <- table(pop(values$filtered_gl))
          cat("  Populations detected:", length(pop_table), "\n")
          cat("  Population sizes:", paste(pop_table, collapse = ", "), "\n")
        } else {
          cat("  WARNING: No population information found!\n")
          cat("  All individuals will be treated as one population.\n")
        }
      }
    } else {
      # Generate summary report
      report <- generate_allele_summary_report(pop_allele_values$results)
      cat(report)
    }
  })
  
  # Main plot
  output$pop_allele_main_plot <- renderPlot({
    req(pop_allele_values$results)
    
    main_plot <- pop_allele_values$results$plots$main
    
    # Customize plot based on user inputs
    if (!is.null(main_plot)) {
      # Update point size
      main_plot$layers[[1]]$aes_params$size <- input$pop_plot_point_size
      
      # Remove labels if requested
      if (!input$pop_plot_labels && length(main_plot$layers) >= 2) {
        main_plot$layers[[2]] <- NULL  # Remove text layer
      }
      
      # Remove legend if requested
      if (!input$pop_plot_legend) {
        main_plot <- main_plot + theme(legend.position = "none")
      }
      
      main_plot
    }
  })
  
  # Additional plots
  output$pop_allele_bar_plot <- renderPlot({
    req(pop_allele_values$results)
    pop_allele_values$results$plots$summary$alleles_bar
  })
  
  output$pop_het_plot <- renderPlot({
    req(pop_allele_values$results, pop_allele_values$heterozygosity_calculated)
    pop_allele_values$results$plots$summary$heterozygosity
  })
  
  output$pop_sample_size_plot <- renderPlot({
    req(pop_allele_values$results)
    pop_allele_values$results$plots$summary$sample_size
  })
  
  # Results table
  output$pop_allele_table <- renderDT({
    req(pop_allele_values$results)
    
    results_table <- pop_allele_values$results$statistics
    
    # Format numeric columns
    datatable(results_table,
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf')
              ),
              rownames = FALSE) %>%
      formatRound(columns = c("N_Alleles_Mean", "N_Alleles_SD", "He_Mean", "Ho_Mean"), 
                  digits = 3) %>%
      formatRound(columns = "N_Alleles_Max", digits = 0)
  })
  
  # Population information
  output$pop_info_summary <- renderPrint({
    req(values$filtered_gl)
    
    if (!is.null(pop(values$filtered_gl))) {
      pop_table <- table(pop(values$filtered_gl))
      
      cat("Population distribution in current dataset:\n")
      cat("==========================================\n\n")
      
      pop_df <- data.frame(
        Population = names(pop_table),
        N_Individuals = as.numeric(pop_table),
        Proportion = round(as.numeric(pop_table) / nInd(values$filtered_gl) * 100, 1)
      )
      
      # Sort by population name
      pop_df <- pop_df[order(pop_df$Population), ]
      
      for (i in 1:nrow(pop_df)) {
        cat(sprintf("%-20s %5d individuals (%5.1f%%)\n", 
                    pop_df$Population[i], 
                    pop_df$N_Individuals[i],
                    pop_df$Proportion[i]))
      }
      
      cat("\nSummary:\n")
      cat("  Total populations:", length(pop_table), "\n")
      cat("  Total individuals:", sum(pop_table), "\n")
      cat("  Mean individuals/pop:", round(mean(pop_table), 1), "\n")
      cat("  Min individuals/pop:", min(pop_table), "\n")
      cat("  Max individuals/pop:", max(pop_table), "\n")
    } else {
      cat("No population information available in the dataset.\n")
      cat("All", nInd(values$filtered_gl), "individuals are unassigned.\n")
    }
  })
  
  # Download handlers
  output$download_pop_allele_csv <- downloadHandler(
    filename = function() {
      paste0("population_allele_stats_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(pop_allele_values$results)
      write.csv(pop_allele_values$results$statistics, file, row.names = FALSE)
    }
  )
  
  output$download_pop_allele_report <- downloadHandler(
    filename = function() {
      paste0("population_allele_report_", Sys.Date(), ".txt")
    },
    content = function(file) {
      req(pop_allele_values$results)
      report <- generate_allele_summary_report(pop_allele_values$results)
      writeLines(report, file)
    }
  )
  
  output$download_main_plot <- downloadHandler(
    filename = function() {
      paste0("population_allele_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(pop_allele_values$results)
      
      # Get current plot
      p <- output$pop_allele_main_plot
      
      # Save as PNG
      ggsave(file, plot = p, width = 12, height = 8, dpi = 300)
    }
  )
  
  output$download_pop_table <- downloadHandler(
    filename = function() {
      paste0("population_allele_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(pop_allele_values$results)
      write.csv(pop_allele_values$results$statistics, file, row.names = FALSE)
    }
  )
  
  # Update population stats UI
  output$pop_allele_stats <- renderUI({
    req(pop_allele_values$results)
    
    stats <- pop_allele_values$results$summary_stats
    
    tagList(
      h4("Quick Statistics:"),
      tags$ul(
        tags$li(strong("Populations analyzed:"), stats$total_populations),
        tags$li(strong("Total individuals:"), stats$total_individuals),
        tags$li(strong("Mean alleles/population:"), round(stats$mean_alleles_per_pop, 2)),
        tags$li(strong("Max alleles/locus:"), stats$max_alleles)
      ),
      br(),
      actionButton("export_pop_allele_json", "Export as JSON",
                   icon = icon("file-export"),
                   class = "btn-warning")
    )
  })
  
  # Export as JSON
  observeEvent(input$export_pop_allele_json, {
    req(pop_allele_values$results)
    
    tryCatch({
      # Create JSON data
      json_data <- list(
        metadata = list(
          analysis_date = as.character(Sys.Date()),
          parameters = list(
            min_individuals = input$pop_min_inds,
            heterozygosity_calculated = input$pop_calc_heterozygosity
          )
        ),
        statistics = pop_allele_values$results$statistics,
        summary = pop_allele_values$results$summary_stats
      )
      
      # Convert to JSON
      json_string <- jsonlite::toJSON(json_data, pretty = TRUE, auto_unbox = TRUE)
      
      # Save to file
      filename <- paste0("population_allele_data_", Sys.Date(), ".json")
      write(json_string, filename)
      
      showNotification(
        paste("Data exported as JSON:", filename),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(
        paste("Failed to export JSON:", e$message),
        type = "error",
        duration = 10
      )
    })
  })
  
  #======================================================================
  #LEA ANALYSIS
  #======================================================================
  
# Helper function to create admixture plot
generate_admixture_plot <- function(q_matrix, k, sample_names = NULL, pop_assignments = NULL) {
  # Prepare data
  df <- as.data.frame(q_matrix)
  n_clusters <- ncol(df)
  
  # Set column names
  colnames(df) <- paste("Cluster", 1:n_clusters, sep = "_")
  
  # Add sample names
  if (is.null(sample_names) || length(sample_names) != nrow(df)) {
    df$Sample <- paste("Sample", 1:nrow(df), sep = "_")
  } else {
    df$Sample <- sample_names
  }
  
  # Add population assignments if available
  if (!is.null(pop_assignments) && length(pop_assignments) == nrow(df)) {
    df$Population <- pop_assignments
    df <- df[order(df$Population, df$Sample), ]
  }
  
  # Convert to long format
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("Cluster_"),
      names_to = "Ancestry",
      values_to = "Proportion"
    )
  
  # Order samples
  df_long$Sample <- factor(df_long$Sample, levels = unique(df$Sample))
  
  # Create plot
  p <- ggplot(df_long, aes(x = Sample, y = Proportion, fill = Ancestry)) +
    geom_bar(stat = "identity", position = "fill")
  
  # Customize based on number of samples
  n_samples <- nrow(df)
  
  if (!is.null(pop_assignments)) {
    # Add population separators
    pop_boundaries <- cumsum(table(df$Population))
    p <- p + geom_vline(
      xintercept = pop_boundaries[-length(pop_boundaries)] + 0.5,
      linetype = "dashed",
      color = "gray40",
      size = 0.5
    )
  }
  
  # Color scheme
  if (n_clusters <= 12) {
    p <- p + scale_fill_brewer(palette = "Set3")
  } else {
    p <- p + scale_fill_viridis_d()
  }
  
  # Theme adjustments
  p <- p +
    labs(
      title = paste("Population Structure (K =", k, ")"),
      x = "Individuals",
      y = "Ancestry Proportion"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = if (n_samples > 30) element_blank() else 
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  return(p)
}



# Helper function for error messages
plot_error_message <- function(message) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = message, size = 5) +
    theme_void() +
    labs(title = "Plot Error")
}



# Helper function for empty data
plot_empty_with_message <- function(message) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = paste(message, "\nPlease load Q matrix data first."), 
             size = 5) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}


# ============================================================================
# LEA ANALYSIS SERVER LOGIC
# ============================================================================

# Reactive values for LEA analysis
lea_values <- reactiveValues(
  results = NULL,
  current_k = 2,
  analysis_status = "Not run",
  available_k = NULL,
  conversion_result = NULL,
  export_files = NULL,
  plot_height = 500
)


# Run LEA analysis
observeEvent(input$run_lea, {
  req(values$filtered_gl)
  
  withProgress(message = 'Running LEA sNMF analysis', value = 0, {
    incProgress(0.2, detail = "This may take several minutes....")
    
    tryCatch({
      # Update status
      lea_values$analysis_status <- "Running analysis..."
      
      # Run LEA analysis with specified parameters
      incProgress(0.4, detail = "Running LEA analysis with specified parameters...")
      results <- run_lea_analysis(
        gl_object = values$filtered_gl,
        output_prefix = "LEA_analysis",
        mindemes = input$lea_min_k,
        maxdemes = input$lea_max_k,
        n_replicates = input$lea_replicates,
        alpha = input$lea_alpha,
        entropy = input$lea_calculate_entropy,
        exportdata = input$lea_export_data
      )
      
      # Check if analysis was successful
      if (is.null(results) || (!is.null(results$success) && !results$success)) {
        error_msg <- ifelse(!is.null(results$error), results$error, "Unknown error")
        stop(error_msg)
      }
      
      # Store results
      lea_values$results <- results
      
      # Update available K values
      incProgress(0.6, detail = "Updating available K values...")
      if (!is.null(results$K_range)) {
        lea_values$available_k <- results$K_range
        
        # Set current K to optimal or median
        if (!is.null(results$optimal_k)) {
          lea_values$current_k <- results$optimal_k
        } else {
          lea_values$current_k <- median(results$K_range)
        }
        
        # Update select input for K
        updateSelectInput(
          session,
          "lea_plot_k",
          choices = results$K_range,
          selected = lea_values$current_k
        )
      }
      
      # Update export files if available
      incProgress(0.8, detail = "Updating available files for export...")
      if (!is.null(results$exported_files)) {
        lea_values$export_files <- results$exported_files
      }
      
      # Update status
      lea_values$analysis_status <- "Completed"
      
      incProgress(0.9, detail = "Completing LEA analysis...")
      
      # Show success notification
      showNotification(
        HTML(paste(
          "<b>LEA Analysis Completed!</b><br>",
          paste("K range analyzed:", input$lea_min_k, "-", input$lea_max_k, "<br>"),
          paste("Optimal K:", ifelse(!is.null(results$optimal_k), 
                                     results$optimal_k, "Not determined"), "<br>"),
          paste("Total samples:", nInd(values$filtered_gl))
        )),
        type = "message"
      )
      
    }, error = function(e) {
      
      # Show error notification
      showNotification(
        paste("LEA Analysis Failed:", e$message),
        type = "error"
      )
      
      # Update status
      lea_values$analysis_status <- paste("Failed:", e$message)
      
      cat("LEA analysis error:", e$message, "\n")
    })
  })
})


# Update K selection
observeEvent(input$lea_plot_k, {
  req(input$lea_plot_k)
  lea_values$current_k <- as.numeric(input$lea_plot_k)
})

# Navigation buttons for K values
observeEvent(input$lea_prev_k, {
  req(lea_values$available_k)
  
  current_idx <- which(lea_values$available_k == lea_values$current_k)
  if (current_idx > 1) {
    new_k <- lea_values$available_k[current_idx - 1]
    lea_values$current_k <- new_k
    updateSelectInput(session, "lea_plot_k", selected = new_k)
  }
})

observeEvent(input$lea_next_k, {
  req(lea_values$available_k)
  
  current_idx <- which(lea_values$available_k == lea_values$current_k)
  if (current_idx < length(lea_values$available_k)) {
    new_k <- lea_values$available_k[current_idx + 1]
    lea_values$current_k <- new_k
    updateSelectInput(session, "lea_plot_k", selected = new_k)
  }
})

# Compare all K values
observeEvent(input$lea_compare_all, {
  req(lea_values$results)
  
  showModal(modalDialog(
    title = "Compare All K Values",
    tagList(
      h4("Cross-Entropy Comparison"),
      plotOutput("lea_all_k_plot", height = "400px"),
      br(),
      h4("Q Matrix Similarity"),
      plotOutput("lea_q_similarity_plot", height = "400px")
    ),
    size = "l",
    footer = modalButton("Close")
  ))
})

# ============================================================================
# OUTPUT RENDERERS
# ============================================================================

# Status output
output$lea_status <- renderPrint({
  cat("LEA Analysis Status\n")
  cat("===================\n\n")
  
  cat("Status:", lea_values$analysis_status, "\n\n")
  
  if (!is.null(lea_values$results)) {
    cat("Analysis Summary\n")
    cat("----------------\n")
    cat("K range analyzed:", paste(range(lea_values$results$K_range), collapse = "-"), "\n")
    cat("Replicates per K:", lea_values$results$n_replicates, "\n")
    cat("Optimal K (suggested):", 
        ifelse(!is.null(lea_values$results$optimal_k), 
               lea_values$results$optimal_k, "Not determined"), "\n")
    cat("Available K values:", paste(lea_values$available_k, collapse = ", "), "\n")
    cat("Current K displayed:", lea_values$current_k, "\n\n")
    
    if (!is.null(lea_values$results$cross_entropy_summary)) {
      cat("Cross-Entropy Summary:\n")
      print(lea_values$results$cross_entropy_summary)
    }
  }
})

# Export buttons UI
output$lea_export_buttons <- renderUI({
  if (is.null(lea_values$results)) {
    return(NULL)
  }
  
  tagList(
    downloadButton("download_lea_single", "Download Current Q Matrix",
                   class = "btn-primary btn-block"),
    br(), br(),
    downloadButton("download_lea_all", "Download All Results",
                   class = "btn-success btn-block"),
    br(), br(),
    actionButton("lea_view_files", "View Export Files",
                 icon = icon("folder-open"),
                 class = "btn-info btn-block")
  )
})

# Cross-entropy plot
output$lea_cross_entropy_plot <- renderPlot({
  req(lea_values$results)
  
  tryCatch({
    plot_cross_entropy_fixed(
      lea_values$results$cross_entropy_summary,
      optimal_k = lea_values$results$optimal_k
    )
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Error creating plot")
    text(1, 1, e$message, col = "red")
  })
})

# Optimal K text
output$lea_optimal_k <- renderUI({
  req(lea_values$results)
  
  if (is.null(lea_values$results$optimal_k)) {
    return(tags$p("Optimal K could not be determined", 
                  style = "color: #666; font-style: italic;"))
  }
  
  tagList(
    tags$h4("Suggested Optimal K:", 
            tags$span(lea_values$results$optimal_k, 
                      style = "color: #d9534f; font-weight: bold;")),
    tags$p("Based on cross-entropy minimization", 
           style = "font-size: 12px; color: #666;"),
    actionButton("lea_use_optimal", "Use This K",
                 icon = icon("check-circle"),
                 class = "btn-success btn-sm")
  )
})

# Use optimal K button
observeEvent(input$lea_use_optimal, {
  req(lea_values$results$optimal_k)
  
  lea_values$current_k <- lea_values$results$optimal_k
  updateSelectInput(session, "lea_plot_k", 
                    selected = lea_values$results$optimal_k)
  
  showNotification(
    paste("Now displaying K =", lea_values$results$optimal_k),
    type = "message",
    duration = 3
  )
})

# Structure plot
output$lea_structure_plot <- renderPlot({
  req(lea_values$results, lea_values$current_k)
  
  tryCatch({
    k_str <- paste0("K", lea_values$current_k)
    q_data <- lea_values$results$Q_matrices[[k_str]]
    
    if (is.null(q_data)) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "No Q matrix available for this K")
      return()
    }
    
    q_matrix <- q_data$Average_Q
    
    # Get sample names
    sample_names <- tryCatch({
      if (!is.null(values$filtered_gl)) {
        indNames(values$filtered_gl)
      } else if (!is.null(rownames(q_matrix))) {
        rownames(q_matrix)
      } else {
        paste0("Sample_", 1:nrow(q_matrix))
      }
    }, error = function(e) {
      paste0("Sample_", 1:nrow(q_matrix))
    })
    
    # Get population labels if available in metadata
    pop_labels <- NULL
    if (!is.null(values$metadata)) {
      # Try to match by sample names
      if ("ID" %in% colnames(values$metadata) && 
          "Population" %in% colnames(values$metadata)) {
        meta_match <- match(sample_names, values$metadata$ID)
        if (any(!is.na(meta_match))) {
          pop_labels <- values$metadata$Population[meta_match]
          pop_labels[is.na(pop_labels)] <- "Unknown"
        }
      }
    }
    
    # Create structure plot
    create_structure_plot_fixed(
      q_matrix = q_matrix,
      sample_names = sample_names,
      pop_labels = pop_labels,
      main = paste("Population Structure - K =", lea_values$current_k),
      sort_by_cluster = TRUE
    )
    
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = paste("Error:", e$message))
  })
}, height = function() { lea_values$plot_height })

# Q matrix table
output$lea_q_matrix_table <- renderDT({
  tryCatch({
    req(lea_values$results, lea_values$current_k)
    
    k_str <- paste0("K", lea_values$current_k)
    q_data <- lea_values$results$Q_matrices[[k_str]]
    
    if (is.null(q_data) || is.null(q_data$Average_Q)) {
      return(datatable(
        data.frame(Message = paste("No Q matrix available for K =", lea_values$current_k)),
        options = list(dom = 't'),
        rownames = FALSE
      ))
    }
    
    q_matrix <- q_data$Average_Q
    
    # Get sample names
    sample_names <- tryCatch({
      if (!is.null(rownames(q_matrix)) && all(rownames(q_matrix) != "")) {
        rownames(q_matrix)
      } else if (!is.null(values$filtered_gl)) {
        indNames(values$filtered_gl)
      } else {
        paste0("Sample_", 1:nrow(q_matrix))
      }
    }, error = function(e) {
      paste0("Sample_", 1:nrow(q_matrix))
    })
    
    # Create data frame
    n_clusters <- ncol(q_matrix)
    q_df <- as.data.frame(q_matrix)
    colnames(q_df) <- paste0("Cluster_", 1:n_clusters)
    
    # Add sample names
    q_df <- cbind(Sample = sample_names, q_df)
    
    # Add dominant cluster and probability
    dominant_cluster <- apply(q_matrix, 1, which.max)
    max_prob <- apply(q_matrix, 1, max)
    
    q_df$Dominant_Cluster <- paste0("Cluster_", dominant_cluster)
    q_df$Assignment_Probability <- round(max_prob, 3)
    
    # Reorder columns
    q_df <- q_df[, c("Sample", "Dominant_Cluster", "Assignment_Probability",
                     paste0("Cluster_", 1:n_clusters))]
    
    # Create DataTable
    dt <- datatable(
      q_df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "400px",
        lengthMenu = c(10, 25, 50, 100),
        dom = 'Blfrtip',
        buttons = list(
          'copy', 
          list(extend = 'csv', title = paste0("Q_Matrix_K", lea_values$current_k)),
          list(extend = 'excel', title = paste0("Q_Matrix_K", lea_values$current_k))
        )
      ),
      extensions = c('Buttons', 'Scroller'),
      rownames = FALSE,
      class = 'display compact stripe hover',
      caption = paste("Q Matrix for K =", lea_values$current_k)
    )
    
    # Format numeric columns
    dt <- dt %>%
      formatRound(columns = 4:(n_clusters + 3), digits = 4) %>%
      formatRound(columns = 3, digits = 3)
    
    # Add color coding for assignment probability
    dt <- dt %>%
      formatStyle(
        'Assignment_Probability',
        background = styleColorBar(q_df$Assignment_Probability, 'lightblue'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
    
    return(dt)
    
  }, error = function(e) {
    return(datatable(
      data.frame(Error = c("Rendering Error", e$message)),
      options = list(dom = 't'),
      rownames = FALSE
    ))
  })
})

# Admixture plot
output$lea_admixture_plot <- renderPlot({
  req(lea_values$results, lea_values$current_k)
  
  tryCatch({
    k_str <- paste0("K", lea_values$current_k)
    q_data <- lea_values$results$Q_matrices[[k_str]]
    
    if (is.null(q_data) || is.null(q_data$Average_Q)) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "No Q matrix available")
      return()
    }
    
    q_matrix <- q_data$Average_Q
    
    # Calculate admixture statistics
    admixture_stats <- calculate_admixture_stats(q_matrix)
    
    # Create admixture plot
    plot_admixture_stats(admixture_stats, lea_values$current_k)
    
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = paste("Error:", e$message))
  })
})

# ============================================================================
# MODAL PLOTS FOR COMPARE ALL
# ============================================================================

# All K plot for comparison
output$lea_all_k_plot <- renderPlot({
  req(lea_values$results)
  
  tryCatch({
    # Plot cross-entropy for all K
    ce_df <- lea_values$results$cross_entropy_summary
    
    ggplot(ce_df, aes(x = K, y = Mean_CE)) +
      geom_line(color = "steelblue", size = 1.5) +
      geom_point(color = "steelblue", size = 3) +
      geom_errorbar(aes(ymin = Mean_CE - SD_CE, ymax = Mean_CE + SD_CE),
                    width = 0.3, color = "gray50") +
      labs(title = "Cross-Entropy Across All K Values",
           x = "Number of Clusters (K)",
           y = "Cross-Entropy") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
    
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Error creating comparison plot")
  })
})

# Q matrix similarity plot
output$lea_q_similarity_plot <- renderPlot({
  req(lea_values$results)
  
  tryCatch({
    # Calculate similarity between Q matrices of different K
    similarity <- calculate_q_similarity(lea_values$results$Q_matrices)
    
    # Create heatmap
    ggplot(similarity, aes(x = K1, y = K2, fill = Similarity)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(name = "Similarity", limits = c(0, 1)) +
      geom_text(aes(label = round(Similarity, 2)), 
                color = "white", size = 3.5) +
      labs(title = "Q Matrix Similarity Between Different K Values",
           x = "K", y = "K") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text = element_text(size = 11),
        legend.position = "right"
      )
    
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Error creating similarity plot")
  })
})

# View export files modal
observeEvent(input$lea_view_files, {
  req(lea_values$export_files)
  
  showModal(modalDialog(
    title = "Export Files",
    tagList(
      h4("Generated Files:"),
      tags$ul(
        lapply(names(lea_values$export_files), function(file_name) {
          tags$li(
            tags$strong(file_name, ":"), 
            tags$a(href = lea_values$export_files[[file_name]], 
                   target = "_blank",
                   basename(lea_values$export_files[[file_name]]))
          )
        })
      ),
      br(),
      p("Files are saved in the output directory.")
    ),
    size = "l",
    footer = modalButton("Close")
  ))
})

# ============================================================================
# DOWNLOAD HANDLERS
# ============================================================================

# Download current Q matrix
output$download_lea_single <- downloadHandler(
  filename = function() {
    paste0("LEA_Qmatrix_K", lea_values$current_k, "_", 
           format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    req(lea_values$results, lea_values$current_k)
    
    k_str <- paste0("K", lea_values$current_k)
    q_data <- lea_values$results$Q_matrices[[k_str]]
    
    if (!is.null(q_data$Average_Q)) {
      write.csv(q_data$Average_Q, file, row.names = TRUE)
    }
  }
)

# Download all results
output$download_lea_all <- downloadHandler(
  filename = function() {
    paste0("LEA_Full_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
  },
  content = function(file) {
    req(lea_values$results)
    
    # Create temporary directory
    temp_dir <- tempdir()
    export_dir <- file.path(temp_dir, "LEA_results")
    
    # Export all results using the fixed function
    export_files <- export_lea_results_fixed(
      results = lea_values$results,
      output_dir = export_dir,
      prefix = "LEA"
    )
    
    if (!is.null(export_files)) {
      # Create zip file
      zip_file <- file.path(temp_dir, "LEA_results.zip")
      zip::zip(zip_file, files = export_dir, recurse = TRUE)
      
      # Copy to download location
      file.copy(zip_file, file)
    }
  }
)

# Download structure plot
output$download_lea_plot <- downloadHandler(
  filename = function() {
    paste0("LEA_Structure_K", lea_values$current_k, "_", 
           format(Sys.time(), "%Y%m%d"), ".png")
  },
  content = function(file) {
    req(lea_values$results, lea_values$current_k)
    
    k_str <- paste0("K", lea_values$current_k)
    q_data <- lea_values$results$Q_matrices[[k_str]]
    
    if (!is.null(q_data$Average_Q)) {
      q_matrix <- q_data$Average_Q
      
      # Get sample names
      sample_names <- tryCatch({
        if (!is.null(values$filtered_gl)) {
          indNames(values$filtered_gl)
        } else {
          paste0("Sample_", 1:nrow(q_matrix))
        }
      }, error = function(e) {
        paste0("Sample_", 1:nrow(q_matrix))
      })
      
      # Create and save plot
      png(file, width = 1200, height = 800)
      create_structure_plot_fixed(
        q_matrix = q_matrix,
        sample_names = sample_names,
        main = paste("Population Structure - K =", lea_values$current_k)
      )
      dev.off()
    }
  }
)

# Download CSV
output$download_lea_csv <- downloadHandler(
  filename = function() {
    paste0("LEA_Q_Matrices_", format(Sys.time(), "%Y%m%d"), ".zip")
  },
  content = function(file) {
    req(lea_values$results)
    
    temp_dir <- tempdir()
    csv_dir <- file.path(temp_dir, "Q_matrices")
    dir.create(csv_dir, recursive = TRUE)
    
    # Export all Q matrices as CSV
    for (k_str in names(lea_values$results$Q_matrices)) {
      q_data <- lea_values$results$Q_matrices[[k_str]]
      if (!is.null(q_data$Average_Q)) {
        csv_file <- file.path(csv_dir, paste0("Q_", k_str, ".csv"))
        write.csv(q_data$Average_Q, csv_file, row.names = TRUE)
      }
    }
    
    # Create zip
    zip_file <- file.path(temp_dir, "Q_matrices.zip")
    zip::zip(zip_file, files = csv_dir, recurse = TRUE)
    file.copy(zip_file, file)
  }
)

# Download RData
output$download_lea_rdata <- downloadHandler(
  filename = function() {
    paste0("LEA_Full_Results_", format(Sys.time(), "%Y%m%d"), ".RData")
  },
  content = function(file) {
    req(lea_values$results)
    
    # Save results to RData file
    save(lea_values$results, file = file)
  }
)

# ============================================================================
# HELPER FUNCTIONS FOR PLOTS
# ============================================================================

#' Create structure plot (fixed version)
create_structure_plot_fixed <- function(q_matrix, sample_names = NULL, 
                                        pop_labels = NULL, main = "Structure Plot",
                                        sort_by_cluster = TRUE, colors = NULL) {
  tryCatch({
    n_individuals <- nrow(q_matrix)
    n_clusters <- ncol(q_matrix)
    
    # Set default sample names
    if (is.null(sample_names)) {
      sample_names <- paste0("Ind_", 1:n_individuals)
    }
    
    # Set colors
    if (is.null(colors)) {
      if (n_clusters <= 8) {
        colors <- RColorBrewer::brewer.pal(n_clusters, "Set2")
      } else {
        colors <- viridis::viridis(n_clusters)
      }
    }
    
    # Prepare data for ggplot
    plot_df <- data.frame(
      Individual = rep(sample_names, each = n_clusters),
      Cluster = rep(paste0("Cluster", 1:n_clusters), times = n_individuals),
      Proportion = as.vector(t(q_matrix)),
      stringsAsFactors = FALSE
    )
    
    # Sort by dominant cluster if requested
    if (sort_by_cluster) {
      dominant_clusters <- apply(q_matrix, 1, which.max)
      cluster_order <- order(dominant_clusters, apply(q_matrix, 1, max))
      plot_df$Individual <- factor(plot_df$Individual, 
                                   levels = sample_names[cluster_order])
    } else {
      plot_df$Individual <- factor(plot_df$Individual, levels = sample_names)
    }
    
    # Create plot
    p <- ggplot(plot_df, aes(x = Individual, y = Proportion, fill = Cluster)) +
      geom_bar(stat = "identity", width = 1) +
      scale_fill_manual(values = colors) +
      labs(title = main,
           x = "Individuals",
           y = "Ancestry Proportion") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, size = 0.5)
      )
    
    # Add population labels if available
    if (!is.null(pop_labels)) {
      p <- p + facet_grid(~ factor(pop_labels, levels = unique(pop_labels)), 
                          scales = "free_x", space = "free")
    }
    
    return(p)
    
  }, error = function(e) {
    cat("Error creating structure plot:", e$message, "\n")
    return(NULL)
  })
}

#' Calculate admixture statistics
calculate_admixture_stats <- function(q_matrix) {
  tryCatch({
    n_individuals <- nrow(q_matrix)
    n_clusters <- ncol(q_matrix)
    
    # Calculate statistics
    dominant_cluster <- apply(q_matrix, 1, which.max)
    max_prob <- apply(q_matrix, 1, max)
    
    stats <- data.frame(
      Cluster = 1:n_clusters,
      N_Individuals = as.vector(table(dominant_cluster)),
      Mean_Proportion = colMeans(q_matrix),
      SD_Proportion = apply(q_matrix, 2, sd),
      Min_Proportion = apply(q_matrix, 2, min),
      Max_Proportion = apply(q_matrix, 2, max)
    )
    
    # Add overall stats
    stats$Overall <- list(
      N_Individuals = n_individuals,
      Mean_Assignment_Prob = mean(max_prob),
      SD_Assignment_Prob = sd(max_prob),
      Proportion_Ambiguous = sum(max_prob < 0.7) / n_individuals
    )
    
    return(stats)
    
  }, error = function(e) {
    cat("Error calculating admixture stats:", e$message, "\n")
    return(NULL)
  })
}

#' Plot admixture statistics
plot_admixture_stats <- function(stats, current_k) {
  tryCatch({
    ggplot(stats, aes(x = factor(Cluster), y = N_Individuals, fill = factor(Cluster))) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = N_Individuals), vjust = -0.5, size = 4) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = paste("Admixture Summary - K =", current_k),
           x = "Cluster",
           y = "Number of Individuals",
           fill = "Cluster") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = "none"
      )
    
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Error creating admixture plot")
  })
}


#' Calculate Q matrix similarity
calculate_q_similarity <- function(q_matrices) {
  tryCatch({
    k_values <- as.numeric(gsub("K", "", names(q_matrices)))
    similarity_df <- data.frame()
    
    for (i in 1:(length(k_values) - 1)) {
      for (j in (i + 1):length(k_values)) {
        k1 <- k_values[i]
        k2 <- k_values[j]
        
        q1 <- q_matrices[[paste0("K", k1)]]$Average_Q
        q2 <- q_matrices[[paste0("K", k2)]]$Average_Q
        
        if (!is.null(q1) && !is.null(q2)) {
          # Calculate correlation between matrices
          # (simplified - in practice you might want a more sophisticated measure)
          sim <- tryCatch({
            cor(as.vector(q1), as.vector(q2), use = "complete.obs")
          }, error = function(e) NA)
          
          similarity_df <- rbind(similarity_df, 
                                 data.frame(K1 = k1, K2 = k2, Similarity = sim))
        }
      }
    }
    
    return(similarity_df)
    
  }, error = function(e) {
    cat("Error calculating similarity:", e$message, "\n")
    return(data.frame(K1 = numeric(), K2 = numeric(), Similarity = numeric()))
  })
}

#============================================================================
# ============================================================================
# SERVER OUTPUT RENDERERS
# ============================================================================

# Reactive values for storing results
hwe_values <- reactiveValues(results = NULL)
dapc_values <- reactiveValues(results = NULL)
hybrid_values <- reactiveValues(results = NULL)

# HWE Results Table
output$popgen_hwe_table <- renderDT({
  req(hwe_values$results)
  
  tryCatch({
    # Combine all population results
    all_results <- do.call(rbind, hwe_values$results$detailed_results)
    
    datatable(
      all_results,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "400px",
        dom = 'Blfrtip',
        buttons = c('copy', 'csv', 'excel'),
        columnDefs = list(
          list(targets = 0, width = '150px')
        )
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      filter = 'top',
      class = 'display compact stripe hover'
    ) %>%
      formatRound(columns = c('H_obs', 'H_exp', 'Fis'), digits = 3) %>%
      formatRound(columns = c('Pr(chi^2 >)'), digits = 4)
  }, error = function(e) {
    datatable(data.frame(Error = e$message))
  })
})

# DAPC Scores Table
output$popgen_dapc_scores <- renderDT({
  req(dapc_values$results)
  
  tryCatch({
    scores_df <- dapc_values$results$scores
    
    datatable(
      scores_df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "400px",
        dom = 'Blfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      class = 'display compact stripe hover'
    ) %>%
      formatRound(columns = grep("^DA", colnames(scores_df), value = TRUE), digits = 3)
  }, error = function(e) {
    datatable(data.frame(Error = e$message))
  })
})

# DAPC Assignment Table
output$popgen_dapc_assignment_table <- renderDT({
  req(dapc_values$results)
  
  tryCatch({
    assignment_df <- dapc_values$results$assignment_table
    
    datatable(
      assignment_df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "400px",
        dom = 'Blfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      extensions = 'Buttons',
      rownames = FALSE,
      class = 'display compact stripe hover'
    ) %>%
      formatRound(columns = grep("^Prob_|^Assignment", colnames(assignment_df), value = TRUE), digits = 3)
  }, error = function(e) {
    datatable(data.frame(Error = e$message))
  })
})

# Hybrid Pattern Visualization
output$hybrid_pattern_plot <- renderPlot({
  req(hybrid_values$results)
  
  tryCatch({
    # Display the main structure plot
    print(hybrid_values$results$plots$structure)
  }, error = function(e) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = paste("Error:", e$message))
  })
}, height = 500)




observeEvent(input$run_hybrid_detection, {
  req(values$filtered_gl, input$hybrid_pop1, input$hybrid_pop2, input$hybrid_admixed)
  
  # Convert to LEA format first
  lea_result <- convert_to_lea_format_fixed(values$filtered_gl)
  
  if (!is.null(lea_result)) {
    showNotification("Running LEA analysis for hybrid detection...", type = "message")
    
    # Run LEA analysis (simplified example)
    tryCatch({
      # Perform sNMF analysis
      project <- LEA::snmf(lea_result$file_path, K = 3, entropy = TRUE, repetitions = 10)
      
      # Get Q matrix
      q_matrix <- LEA::Q(project, K = 3)
      
      # Run hybrid visualization
      hybrid_values$results <- visualize_hybrid_patterns(
        q_matrix = q_matrix,
        sample_names = indNames(values$filtered_gl),
        threshold = 0.1
      )
    }, error = function(e) {
      showNotification(paste("LEA analysis failed:", e$message), type = "error")
    })
  }
})

#===============================================================================
# ============================================================================
# HYBRID ANALYSIS SERVER CODE
# ============================================================================

# Initialize reactive values for hybrid analysis
hybrid_values <- reactiveValues(
  results = NULL,
  summary = NULL,
  classification = NULL,
  parent_references = data.frame(
    Sample = character(),
    Population = character(),
    Type = character(),
    stringsAsFactors = FALSE
  )
)

# Observer to update parent population dropdowns
observe({
  req(values$filtered_gl, values$metadata)
  
  # Get available populations from metadata
  if ("Population" %in% colnames(values$metadata)) {
    populations <- unique(na.omit(values$metadata$Population))
    populations <- sort(as.character(populations))
    
    # Update all population selectors
    updateSelectInput(session, "hybrid_pop1", choices = c("", populations))
    updateSelectInput(session, "hybrid_pop2", choices = c("", populations))
    updateSelectInput(session, "hybrid_admixed", choices = c("", populations))
  } else if (nInd(values$filtered_gl) > 0) {
    # If no population column, use individual names
    ind_names <- indNames(values$filtered_gl)
    updateSelectInput(session, "hybrid_pop1", choices = c("", ind_names))
    updateSelectInput(session, "hybrid_pop2", choices = c("", ind_names))
    updateSelectInput(session, "hybrid_admixed", choices = c("", ind_names))
  }
})

# Observer to run hybrid analysis when button is clicked
observeEvent(input$run_hybrid_detection, {
  req(values$filtered_gl, input$hybrid_pop1, input$hybrid_pop2, input$hybrid_admixed)
  
  # Show loading notification
  showModal(modalDialog(
    title = "Running Hybrid Analysis",
    "Please wait while the hybrid analysis is running...",
    footer = NULL,
    size = "s"
  ))
  
  tryCatch({
    cat("\n=== Starting Hybrid Analysis ===\n")
    
    # Convert genlight to LEA format
    cat("  Converting genotype data to LEA format...\n")
    lea_conversion <- convert_to_lea_format_fixed(values$filtered_gl)
    
    if (is.null(lea_conversion)) {
      stop("Failed to convert genotype data to LEA format")
    }
    
    # Run sNMF analysis (K=3 for two parents + hybrid)
    cat("  Running sNMF analysis (K=3)...\n")
    project <- NULL
    
    # Try to run sNMF with error handling
    tryCatch({
      project <- LEA::snmf(
        input.file = lea_conversion$file_path,
        K = 3,
        entropy = TRUE,
        repetitions = 10,
        project = "new",
        seed = 12345
      )
    }, error = function(e) {
      cat("  sNMF failed, trying alternative approach...\n")
      
      # Alternative: Use a simpler K value
      tryCatch({
        project <- LEA::snmf(
          input.file = lea_conversion$file_path,
          K = 2,
          entropy = TRUE,
          repetitions = 5,
          project = "new",
          seed = 12345
        )
      }, error = function(e2) {
        stop("Both sNMF methods failed: ", e$message, " and ", e2$message)
      })
    })
    
    if (is.null(project)) {
      stop("Failed to run sNMF analysis")
    }
    
    # Get Q matrix (admixture proportions)
    cat("  Extracting admixture proportions...\n")
    best_run <- which.min(LEA::cross.entropy(project, K = ifelse(!is.null(project@K), project@K, 3)))
    q_matrix <- LEA::Q(project, K = ifelse(!is.null(project@K), project@K, 3), run = best_run)
    
    # Set row names
    rownames(q_matrix) <- indNames(values$filtered_gl)
    colnames(q_matrix) <- paste0("Cluster_", 1:ncol(q_matrix))
    
    # Create parent references based on selected populations
    cat("  Creating parent references...\n")
    parent_refs <- data.frame(
      Sample = character(),
      Population = character(),
      Type = character(),
      stringsAsFactors = FALSE
    )
    
    # Add parent 1 samples
    if (input$hybrid_pop1 != "") {
      if ("Population" %in% colnames(values$metadata)) {
        parent1_samples <- values$metadata$Name[values$metadata$Population == input$hybrid_pop1]
      } else {
        parent1_samples <- input$hybrid_pop1
      }
      
      for (sample in parent1_samples) {
        parent_refs <- rbind(parent_refs, data.frame(
          Sample = sample,
          Population = input$hybrid_pop1,
          Type = "Parent1",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add parent 2 samples
    if (input$hybrid_pop2 != "") {
      if ("Population" %in% colnames(values$metadata)) {
        parent2_samples <- values$metadata$Name[values$metadata$Population == input$hybrid_pop2]
      } else {
        parent2_samples <- input$hybrid_pop2
      }
      
      for (sample in parent2_samples) {
        parent_refs <- rbind(parent_refs, data.frame(
          Sample = sample,
          Population = input$hybrid_pop2,
          Type = "Parent2",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add hybrid/admixed samples
    if (input$hybrid_admixed != "") {
      if ("Population" %in% colnames(values$metadata)) {
        hybrid_samples <- values$metadata$Name[values$metadata$Population == input$hybrid_admixed]
      } else {
        hybrid_samples <- input$hybrid_admixed
      }
      
      for (sample in hybrid_samples) {
        parent_refs <- rbind(parent_refs, data.frame(
          Sample = sample,
          Population = input$hybrid_admixed,
          Type = "Hybrid",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Store parent references
    hybrid_values$parent_references <- parent_refs
    
    # Run hybrid visualization
    cat("  Running hybrid visualization...\n")
    hybrid_results <- visualize_hybrid_patterns(
      q_matrix = q_matrix,
      sample_names = indNames(values$filtered_gl),
      parent_references = parent_refs,
      threshold = 0.1
    )
    
    if (is.null(hybrid_results)) {
      stop("Hybrid visualization failed")
    }
    
    # Store results
    hybrid_values$results <- hybrid_results
    
    # Create summary statistics
    if (!is.null(hybrid_results$classification)) {
      cat("  Creating summary statistics...\n")
      
      # Basic classification summary
      class_summary <- table(hybrid_results$classification$Category)
      
      # Calculate hybrid indices for hybrids
      if ("Hybrid" %in% hybrid_results$classification$Category) {
        hybrids <- hybrid_results$classification[hybrid_results$classification$Category == "Hybrid", ]
        hybrid_index_stats <- summary(hybrids$Max_Proportion)
      } else {
        hybrid_index_stats <- NULL
      }
      
      # Create comprehensive summary
      hybrid_values$summary <- list(
        analysis_date = Sys.time(),
        n_samples_total = nrow(hybrid_results$classification),
        n_clusters = ncol(q_matrix),
        classification_summary = as.data.frame(class_summary),
        hybrid_index_stats = hybrid_index_stats,
        parent_populations = c(input$hybrid_pop1, input$hybrid_pop2),
        hybrid_population = input$hybrid_admixed,
        admixture_threshold = 0.1,
        lea_conversion_success = !is.null(lea_conversion),
        snmf_success = !is.null(project)
      )
      
      cat("  Hybrid analysis completed successfully!\n")
      cat("    Total samples:", hybrid_values$summary$n_samples_total, "\n")
      cat("    Clusters:", hybrid_values$summary$n_clusters, "\n")
      cat("    Pure individuals:", sum(hybrid_results$classification$Category == "Pure"), "\n")
      cat("    Hybrids:", sum(hybrid_results$classification$Category == "Hybrid"), "\n")
      cat("    Admixed individuals:", sum(hybrid_results$classification$Category == "Admixed"), "\n")
    }
    
    # Store classification for later use
    hybrid_values$classification <- hybrid_results$classification
    
    # Success notification
    removeModal()
    showNotification(
      "Hybrid analysis completed successfully!",
      type = "message",
      duration = 5
    )
    
  }, error = function(e) {
    # Error handling
    removeModal()
    
    # Show error notification
    showNotification(
      paste("Hybrid analysis failed:", e$message),
      type = "error",
      duration = 10
    )
    
    cat("ERROR in hybrid analysis:", e$message, "\n")
    cat("Traceback:\n")
    print(traceback())
    
    # Store error in results for debugging
    hybrid_values$results <- list(
      error = paste("Hybrid analysis failed:", e$message),
      timestamp = Sys.time()
    )
    
    hybrid_values$summary <- list(
      error = e$message,
      analysis_date = Sys.time(),
      status = "Failed"
    )
  })
})


# ============================================================================
# LEA CONVERSION FUNCTION - ADD THIS TO YOUR SERVER FILE
# ============================================================================

#' Convert genlight to LEA format - FIXED VERSION
convert_to_lea_format_fixed <- function(gl_object, output_file = NULL) {
  tryCatch({
    cat("  Converting genlight to LEA format...\n")
    
    if (is.null(gl_object)) {
      stop("Genlight object is NULL")
    }
    
    # Get basic information
    n_ind <- nInd(gl_object)
    n_loc <- nLoc(gl_object)
    
    cat("    Individuals:", n_ind, "\n")
    cat("    Loci:", n_loc, "\n")
    
    # Convert to matrix (0,1,2 coding)
    geno_matrix <- as.matrix(gl_object)
    
    # Check matrix dimensions
    if (is.null(dim(geno_matrix)) || nrow(geno_matrix) != n_ind || ncol(geno_matrix) != n_loc) {
      cat("    Warning: Matrix dimensions don't match genlight object\n")
      cat("    Matrix dim:", dim(geno_matrix), "\n")
      
      # Try to fix by transposing if needed
      if (nrow(geno_matrix) == n_loc && ncol(geno_matrix) == n_ind) {
        cat("    Transposing matrix...\n")
        geno_matrix <- t(geno_matrix)
      }
    }
    
    # Convert to numeric if needed
    if (!is.numeric(geno_matrix)) {
      cat("    Converting to numeric...\n")
      geno_matrix <- apply(geno_matrix, 2, as.numeric)
    }
    
    # Check for any non-numeric values
    if (any(!is.numeric(geno_matrix), na.rm = TRUE)) {
      cat("    Warning: Non-numeric values found in genotype matrix\n")
      geno_matrix[!is.numeric(geno_matrix)] <- 9
    }
    
    # Handle missing values (replace NA with 9 as per LEA format)
    if (any(is.na(geno_matrix))) {
      cat("    Replacing", sum(is.na(geno_matrix)), "NA values with 9...\n")
      geno_matrix[is.na(geno_matrix)] <- 9
    }
    
    # Ensure values are integers (0, 1, 2, 9)
    if (!all(geno_matrix %in% c(0, 1, 2, 9), na.rm = TRUE)) {
      cat("    Warning: Some values are not 0,1,2,9. Rounding to nearest integer...\n")
      geno_matrix <- round(geno_matrix)
      geno_matrix[geno_matrix < 0] <- 0
      geno_matrix[geno_matrix > 2 & geno_matrix != 9] <- 2
    }
    
    # Transpose for LEA format (SNPs as rows, individuals as columns)
    cat("    Transposing matrix for LEA format...\n")
    geno_matrix_t <- t(geno_matrix)
    
    # Ensure integer storage mode
    storage.mode(geno_matrix_t) <- "integer"
    
    # Create output file name
    if (is.null(output_file)) {
      output_file <- tempfile(fileext = ".geno")
      cat("    Using temporary file:", output_file, "\n")
    }
    
    # Write to file in LEA format
    cat("    Writing to file...\n")
    
    # Method 1: Try using LEA's write.geno function
    if (requireNamespace("LEA", quietly = TRUE)) {
      tryCatch({
        # Write in GENO format (more reliable than LFMM)
        LEA::write.geno(geno_matrix_t, output_file)
        cat("    Written using LEA::write.geno\n")
      }, error = function(e) {
        cat("    LEA::write.geno failed, using manual write...\n")
        # Manual write
        write.table(geno_matrix_t, 
                    file = output_file,
                    sep = "",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    na = "9")
      })
    } else {
      # Manual write
      write.table(geno_matrix_t, 
                  file = output_file,
                  sep = "",
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE,
                  na = "9")
    }
    
    # Verify file was created
    if (!file.exists(output_file)) {
      stop("Failed to create output file")
    }
    
    file_size <- file.info(output_file)$size
    cat("    File created:", output_file, "(", file_size, "bytes)\n")
    
    # Create lfmm object (using the transposed matrix)
    lfmm_obj <- tryCatch({
      LEA::as.lfmm(geno_matrix_t)
    }, error = function(e) {
      cat("    Warning: Could not create lfmm object:", e$message, "\n")
      NULL
    })
    
    return(list(
      file_path = output_file,
      lfmm_object = lfmm_obj,
      n_individuals = n_ind,
      n_snps = n_loc,
      matrix_dimensions = dim(geno_matrix_t)
    ))
    
  }, error = function(e) {
    cat("  ERROR in LEA conversion:", e$message, "\n")
    cat("  Traceback:\n")
    print(traceback())
    return(NULL)
  })
}

# ============================================================================
# SIMPLIFIED LEA CONVERSION FOR SHINY (Alternative)
# ============================================================================

#' Simplified LEA conversion for Shiny - less error prone
convert_to_lea_simple <- function(gl_object) {
  tryCatch({
    cat("  Converting to LEA format (simplified)...\n")
    
    if (is.null(gl_object)) {
      return(NULL)
    }
    
    # Convert to matrix
    geno_matrix <- as.matrix(gl_object)
    
    # Replace NA with 9 (LEA format for missing)
    geno_matrix[is.na(geno_matrix)] <- 9
    
    # Round to nearest integer
    geno_matrix <- round(geno_matrix)
    
    # Ensure values are 0, 1, 2, or 9
    geno_matrix[geno_matrix < 0] <- 0
    geno_matrix[geno_matrix > 2 & geno_matrix != 9] <- 2
    
    # Transpose for LEA format
    geno_t <- t(geno_matrix)
    
    # Save to temporary file
    output_file <- tempfile(fileext = ".geno")
    
    # Write in LEA format
    write.table(geno_t, 
                file = output_file,
                sep = "",
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE)
    
    cat("    Saved to:", output_file, "\n")
    
    return(output_file)
    
  }, error = function(e) {
    cat("  ERROR in simple LEA conversion:", e$message, "\n")
    return(NULL)
  })
}



# Render hybrid results summary as text
output$hybrid_results_summary <- renderPrint({
  req(hybrid_values$summary)
  
  cat("HYBRID ANALYSIS SUMMARY\n")
  cat("=======================\n\n")
  
  # Check if there was an error
  if (!is.null(hybrid_values$summary$error)) {
    cat("STATUS: FAILED\n")
    cat("Error:", hybrid_values$summary$error, "\n")
    cat("\nTroubleshooting:\n")
    cat("1. Ensure genotype data is properly loaded\n")
    cat("2. Check that selected populations exist in the data\n")
    cat("3. Try different parent population selections\n")
    cat("4. Check console for detailed error messages\n")
    return()
  }
  
  # Display success summary
  cat("STATUS: SUCCESS\n")
  cat("Analysis Date:", format(hybrid_values$summary$analysis_date, "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("POPULATION SELECTIONS\n")
  cat("---------------------\n")
  cat("Parent Population 1:", hybrid_values$summary$parent_populations[1], "\n")
  cat("Parent Population 2:", hybrid_values$summary$parent_populations[2], "\n")
  cat("Putative Hybrid Population:", hybrid_values$summary$hybrid_population, "\n\n")
  
  cat("ANALYSIS PARAMETERS\n")
  cat("-------------------\n")
  cat("Total Samples Analyzed:", hybrid_values$summary$n_samples_total, "\n")
  cat("Number of Clusters (K):", hybrid_values$summary$n_clusters, "\n")
  cat("Admixture Threshold:", hybrid_values$summary$admixture_threshold, "\n")
  cat("LEA Conversion:", ifelse(hybrid_values$summary$lea_conversion_success, "Success", "Failed"), "\n")
  cat("sNMF Analysis:", ifelse(hybrid_values$summary$snmf_success, "Success", "Failed"), "\n\n")
  
  cat("CLASSIFICATION RESULTS\n")
  cat("----------------------\n")
  if (!is.null(hybrid_values$summary$classification_summary)) {
    print(hybrid_values$summary$classification_summary)
    cat("\n")
    
    # Calculate percentages
    total <- sum(hybrid_values$summary$classification_summary$Freq)
    for (i in 1:nrow(hybrid_values$summary$classification_summary)) {
      category <- hybrid_values$summary$classification_summary$Var1[i]
      count <- hybrid_values$summary$classification_summary$Freq[i]
      percentage <- round(count / total * 100, 1)
      cat(sprintf("  %s: %d samples (%.1f%%)\n", category, count, percentage))
    }
  }
  
  cat("\nHYBRID CHARACTERISTICS\n")
  cat("---------------------\n")
  if (!is.null(hybrid_values$summary$hybrid_index_stats)) {
    cat("For samples classified as 'Hybrid':\n")
    cat("  Mean admixture proportion:", round(mean(hybrid_values$summary$hybrid_index_stats), 3), "\n")
    cat("  Range:", round(range(hybrid_values$summary$hybrid_index_stats), 3), "\n")
    cat("  Standard deviation:", round(sd(hybrid_values$summary$hybrid_index_stats), 3), "\n")
  } else {
    cat("No hybrids detected with current threshold.\n")
  }
  
  cat("\nINTERPRETATION\n")
  cat("--------------\n")
  cat("• 'Pure' individuals: >90% ancestry from one cluster\n")
  cat("• 'Hybrid' individuals: 50-90% ancestry from dominant cluster\n")
  cat("• 'Admixed' individuals: <50% ancestry from any single cluster\n")
  
  cat("\nNEXT STEPS\n")
  cat("----------\n")
  cat("1. Check the hybrid pattern visualization plot\n")
  cat("2. Download results using the download buttons\n")
  cat("3. Adjust threshold if classification seems inaccurate\n")
  cat("4. Run with different parent population selections if needed\n")
})

# Render hybrid pattern visualization plot
output$hybrid_pattern_plot <- renderPlot({
  req(hybrid_values$results)
  
  tryCatch({
    # Check if there was an error in analysis
    if (!is.null(hybrid_values$results$error)) {
      # Create an error plot
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Hybrid Analysis Failed",
           sub = hybrid_values$results$error)
      text(1, 1, "Please check console for details", cex = 1.2, col = "red")
      return()
    }
    
    # Check if we have plots
    if (is.null(hybrid_values$results$plots)) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "No Hybrid Visualization Available",
           sub = "The analysis completed but no plots were generated")
      text(1, 1, "Check analysis parameters and data", cex = 1.2)
      return()
    }
    
    # Display the main structure plot if available
    if (!is.null(hybrid_values$results$plots$structure)) {
      print(hybrid_values$results$plots$structure)
    } 
    # Alternatively, display triangle plot for 3 clusters
    else if (!is.null(hybrid_values$results$plots$triangle)) {
      print(hybrid_values$results$plots$triangle)
    }
    # Fallback: create a basic admixture plot
    else if (!is.null(hybrid_values$results$visualization_data$structure)) {
      # Create basic structure plot from data
      structure_data <- hybrid_values$results$visualization_data$structure
      
      p <- ggplot(structure_data, aes(x = Sample, y = Proportion, fill = Cluster)) +
        geom_bar(stat = "identity", width = 1) +
        scale_fill_brewer(palette = "Set2") +
        labs(title = "Admixture Proportions",
             x = "Samples", 
             y = "Proportion",
             fill = "Cluster") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      
      print(p)
    } 
    # Final fallback
    else {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Hybrid Pattern Visualization",
           sub = "No plot data available")
      text(1, 1, "Analysis completed but visualization data is missing", cex = 1.2)
    }
    
  }, error = function(e) {
    # Error in plotting
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Plot Rendering Error",
         sub = e$message)
    text(1, 1, "Please check the data and try again", cex = 1.2, col = "red")
  })
}, height = 500)

# Download handler for hybrid results
output$download_hybrid_results <- downloadHandler(
  filename = function() {
    paste0("hybrid_analysis_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
  },
  content = function(file) {
    req(hybrid_values$results, hybrid_values$summary)
    
    # Create temporary directory
    temp_dir <- tempfile("hybrid_results_")
    dir.create(temp_dir)
    
    tryCatch({
      # 1. Save summary as text file
      summary_file <- file.path(temp_dir, "hybrid_analysis_summary.txt")
      sink(summary_file)
      cat("HYBRID ANALYSIS REPORT\n")
      cat("======================\n\n")
      cat("Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
      
      if (!is.null(hybrid_values$summary$error)) {
        cat("ANALYSIS FAILED\n")
        cat("Error:", hybrid_values$summary$error, "\n")
      } else {
        cat("ANALYSIS SUCCESSFUL\n\n")
        cat("Parameters:\n")
        cat("  Parent Population 1:", hybrid_values$summary$parent_populations[1], "\n")
        cat("  Parent Population 2:", hybrid_values$summary$parent_populations[2], "\n")
        cat("  Hybrid Population:", hybrid_values$summary$hybrid_population, "\n")
        cat("  Threshold:", hybrid_values$summary$admixture_threshold, "\n")
        cat("  Total Samples:", hybrid_values$summary$n_samples_total, "\n\n")
        
        cat("Classification Results:\n")
        print(hybrid_values$summary$classification_summary)
      }
      sink()
      
      # 2. Save classification table as CSV
      if (!is.null(hybrid_values$classification)) {
        class_file <- file.path(temp_dir, "hybrid_classification.csv")
        write.csv(hybrid_values$classification, class_file, row.names = FALSE)
      }
      
      # 3. Save visualization data
      if (!is.null(hybrid_values$results$visualization_data)) {
        viz_file <- file.path(temp_dir, "visualization_data.RData")
        save(hybrid_values$results$visualization_data, file = viz_file)
      }
      
      # 4. Save plots as PDF
      if (!is.null(hybrid_values$results$plots)) {
        pdf_file <- file.path(temp_dir, "hybrid_plots.pdf")
        pdf(pdf_file, width = 10, height = 8)
        
        if (!is.null(hybrid_values$results$plots$structure)) {
          print(hybrid_values$results$plots$structure)
        }
        if (!is.null(hybrid_values$results$plots$triangle)) {
          print(hybrid_values$results$plots$triangle)
        }
        
        dev.off()
      }
      
      # Create ZIP file
      zip_files <- list.files(temp_dir, full.names = TRUE)
      zip(file, zip_files, flags = "-j")
      
    }, error = function(e) {
      showNotification(paste("Error creating download:", e$message), type = "error")
    })
  }
)

# Helper function to check if hybrid analysis has been run
output$hybrid_analysis_completed <- reactive({
  return(!is.null(hybrid_values$results) && is.null(hybrid_values$summary$error))
})
outputOptions(output, "hybrid_analysis_completed", suspendWhenHidden = FALSE)

# UI element to show download button only when analysis is complete
output$hybrid_download_ui <- renderUI({
  if (!is.null(hybrid_values$results) && is.null(hybrid_values$summary$error)) {
    tagList(
      hr(),
      h4("Download Results"),
      downloadButton("download_hybrid_results", "Download All Hybrid Results", 
                     class = "btn-success"),
      helpText("Includes summary, classification table, visualization data, and plots.")
    )
  }
})


# Setup LEA server logic
setup_lea_server(input, output, session, values)

  
}#End of server

# Run the Shiny app
shinyApp(ui = ui, server = server)