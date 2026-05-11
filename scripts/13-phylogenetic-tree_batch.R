# =============================================================================
# Script:      phylogenetic-tree_batch.R
# Description: Automated batch processing of phylogenetic trees for multiple
#              gene families. Iterates over all genes within a hormone/project
#              category, renders phylogenetic-tree.Rmd for each gene in both
#              cladogram and phylogram modes, and builds a consolidated
#              presence/absence table across all gene families.
#
# Usage:       Rscript phylogenetic-tree_batch.R
#              (or source from RStudio after editing USER CONFIGURATION)
#
# Input:       phylo/{HORMONE}/07_iqtree/{GENE}/{GENE}.nwk
#              data/metadata.tsv
#              scripts/12-phylogenetic-tree.Rmd  — visualization template
#
# Output:      phylo/{HORMONE}/08_figures/{GENE}/  — SVG and PNG figures
#              phylo/{HORMONE}/09_presence-absence/ — per-gene tables
#              phylo/{HORMONE}/merged-genes.tsv     — consolidated table
#
# Dependencies (R): rmarkdown, here, stringr
#                   (all packages required by phylogenetic-tree.Rmd)
#
# Author:      David Alberto García Estrada
#              ORCID: 0009-0007-1169-5329
# Last updated: 2026-01-27
# =============================================================================


#############################################
### RENV ENVIRONMENT ACTIVATION
#############################################
# Activate renv for reproducible package management
# This ensures all collaborators use the same package versions

if (file.exists("renv/activate.R")) {
  cat("🔧 Activating renv environment...\n")
  source("renv/activate.R")
  cat("✓ renv environment activated\n\n")
} else {
  warning("⚠️  renv/activate.R not found. Proceeding without renv.")
  cat("   Recommendation: Initialize renv with renv::init()\n\n")
}

#############################################
### PACKAGE DEPENDENCIES CHECK
#############################################
# Verify that required packages are installed

if (!require("rmarkdown", quietly = TRUE)) {
  stop("❌ ERROR: rmarkdown package not found.\n",
       "   Solution: Install with install.packages('rmarkdown') or renv::install('rmarkdown')")
}

if (!require("here", quietly = TRUE)) {
  stop("❌ ERROR: here package not found.\n",
       "   Solution: Install with install.packages('here') or renv::install('here')")
}

library(rmarkdown)
library(here)

cat("✓ Required packages loaded: rmarkdown, here\n\n")

#############################################
### USER CONFIGURATION
#############################################
# EDIT THESE PARAMETERS FOR YOUR ANALYSIS

# Hormone/project category to process
# Examples: "ABA", "CK", "Auxin"...
HORMONE <- "ABA"  # <-- CHANGE THIS

# R Markdown template filename
# This should be the improved phylogenetic tree template
RMD_TEMPLATE <- "scripts/12-phylogenetic-tree.Rmd"

# Visualization modes to generate
# "none"         = Cladogram (topology-focused, tips aligned)
# "proportional" = Phylogram (branch lengths represent evolutionary distance)
MODES <- c("none", "proportional")

# Optional: Override automatic gene detection
# Leave as NULL to auto-detect all .nwk files
# Or specify: GENES_MANUAL <- c("IPT", "CYP450", "NCED")
GENES_MANUAL <- "RiHHK6"

#############################################
### AUTOMATIC DETECTION OF .NWK FILES
#############################################
# Scans the IQ-TREE output directory for gene folders containing .nwk files

base_dir <- here::here()
iqtree_dir <- file.path(base_dir, "phylogenetic", HORMONE, "07_iqtree")

cat("📁 Scanning for phylogenetic trees...\n")
cat("   Base directory:", base_dir, "\n")
cat("   IQ-TREE directory:", iqtree_dir, "\n\n")

# Verify that IQ-TREE directory exists
if (!dir.exists(iqtree_dir)) {
  stop("❌ ERROR: IQ-TREE directory not found: ", iqtree_dir, "\n",
       "   Check that HORMONE = '", HORMONE, "' is correct\n",
       "   Expected structure: phylo/{HORMONE}/07_iqtree/")
}

# Determine which genes to process
if (!is.null(GENES_MANUAL)) {
  # Manual mode: use specified genes
  genes_to_process <- GENES_MANUAL
  cat("📋 Manual gene selection mode\n")
  cat("   Processing specified genes:", paste(genes_to_process, collapse = ", "), "\n\n")
  
  # Verify that .nwk files exist for all specified genes
  for (gene in genes_to_process) {
    nwk_file <- file.path(iqtree_dir, gene, paste0(gene, ".nwk"))
    if (!file.exists(nwk_file)) {
      warning("⚠️  .nwk file not found for gene: ", gene, "\n   Expected: ", nwk_file)
    }
  }
  
} else {
  # Automatic mode: detect all genes with .nwk files
  cat("🔍 Auto-detection mode\n")
  
  gene_dirs <- list.dirs(iqtree_dir, full.names = FALSE, recursive = FALSE)
  genes_to_process <- character()
  
  for (gene_dir in gene_dirs) {
    gene_path <- file.path(iqtree_dir, gene_dir)
    nwk_file <- file.path(gene_path, paste0(gene_dir, ".nwk"))
    
    if (file.exists(nwk_file)) {
      genes_to_process <- c(genes_to_process, gene_dir)
    }
  }
  
  cat("   Found", length(genes_to_process), "gene(s) with .nwk files\n\n")
}

# Verify that we have genes to process
if (length(genes_to_process) == 0) {
  stop("❌ ERROR: No .nwk files found in ", iqtree_dir, "\n",
       "   Each gene should have structure: {gene}/{gene}.nwk")
}

#############################################
### PROCESSING CONFIGURATION SUMMARY
#############################################

cat("╔════════════════════════════════════════╗\n")
cat("║   BATCH PROCESSING CONFIGURATION       ║\n")
cat("╚════════════════════════════════════════╝\n\n")

cat("📁 Hormone/Category:", HORMONE, "\n")
cat("🧬 Genes to process:", length(genes_to_process), "\n")
cat("🎨 Visualization modes:", paste(MODES, collapse = ", "), "\n")
cat("📊 Total trees to generate:", length(genes_to_process) * length(MODES), "\n")
cat("📋 R Markdown template:", RMD_TEMPLATE, "\n\n")

cat("📊 Presence/absence tables: ENABLED\n")
cat("   Individual tables: phylo/", HORMONE, "/09_presence-absence/{gene}.tsv\n", sep = "")
cat("   Consolidated table: phylo/", HORMONE, "/merged-genes.tsv\n\n", sep = "")

cat("🧬 Genes to process:\n")
for (i in seq_along(genes_to_process)) {
  cat("   ", sprintf("%2d", i), ". ", genes_to_process[i], "\n", sep = "")
}

cat("\n")

# User confirmation
cat("⚠️  This will generate", length(genes_to_process) * length(MODES), "phylogenetic trees\n")
cat("   Processing time estimate: ~", ceiling(length(genes_to_process) * length(MODES) * 0.5), 
    " minutes\n\n", sep = "")

response <- readline(prompt = "▶ Press [Enter] to continue or [Ctrl+C] to cancel... ")

#############################################
### VERIFY TEMPLATE EXISTS
#############################################

if (!file.exists(RMD_TEMPLATE)) {
  stop("❌ ERROR: R Markdown template not found: ", RMD_TEMPLATE, "\n",
       "   Make sure the template file is in the current directory\n",
       "   Expected: ", file.path(base_dir, RMD_TEMPLATE))
}

cat("\n✓ R Markdown template found\n\n")


#############################################
### BATCH PROCESSING LOOP
#############################################

cat("╔════════════════════════════════════════╗\n")
cat("║   STARTING BATCH PROCESSING            ║\n")
cat("╚════════════════════════════════════════╝\n\n")

# Initialize counters
total_processed <- 0
total_errors <- 0
error_log <- character()
start_time <- Sys.time()

# Process each gene in each mode
for (gene in genes_to_process) {
  for (mode in MODES) {
    
    cat("\n═══════════════════════════════════════\n")
    cat("📍 Processing:", gene, "|", mode, "\n")
    cat("═══════════════════════════════════════\n")
    
    # Create new environment for rendering
    # This prevents variable conflicts between iterations
    render_env <- new.env()
    render_env$params <- list(
      hormone = HORMONE,
      file = gene,
      mode = mode
    )
    
    # Render R Markdown with parameters
    result <- tryCatch({
      rmarkdown::render(
        input = RMD_TEMPLATE,
        envir = render_env,
        quiet = FALSE  # Show rendering progress
      )
      total_processed <- total_processed + 1
      cat("✅ SUCCESS:", gene, "|", mode, "\n")
      TRUE
      
    }, error = function(e) {
      total_errors <- total_errors + 1
      error_msg <- paste0(gene, " | ", mode, ": ", e$message)
      error_log <<- c(error_log, error_msg)
      cat("❌ ERROR:", e$message, "\n")
      FALSE
    })
    
    # Brief pause between iterations to prevent system overload
    Sys.sleep(0.5)
  }
}

# Calculate elapsed time
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

#############################################
### GENERATE CONSOLIDATED PRESENCE/ABSENCE TABLE
#############################################

cat("\n\n╔═══════════════════════════════════════╗\n")
cat("║   CREATING CONSOLIDATED TABLE          ║\n")
cat("╚═══════════════════════════════════════╝\n\n")

# Define paths
presence_absence_dir <- file.path(base_dir, "phylogenetic", HORMONE, "09_presence-absence")
merged_table_file <- file.path(base_dir, "phylogenetic", HORMONE, "merged-genes.tsv")
metadata_file <- file.path(base_dir, "data", "Metadata.tsv")

# Check if required files/directories exist
if (!file.exists(metadata_file)) {
  cat("⚠️  WARNING: Metadata file not found:", metadata_file, "\n")
  cat("   Skipping consolidated table generation\n\n")
  
} else if (!dir.exists(presence_absence_dir)) {
  cat("⚠️  WARNING: Presence/absence directory not found:", presence_absence_dir, "\n")
  cat("   Skipping consolidated table generation\n\n")
  
} else {
  
  # Load metadata
  cat("📖 Loading metadata...\n")
  metadata <- tryCatch({
    read.table(metadata_file, header = TRUE, sep = "\t", 
               stringsAsFactors = FALSE, quote = "", 
               fill = TRUE, comment.char = "")
  }, error = function(e) {
    stop("❌ ERROR reading metadata: ", e$message)
  })
  
  metadata$Organism_clean <- stringr::str_trim(metadata$Organism)
  cat("   ✓ Loaded", nrow(metadata), "organisms from metadata\n\n")
  
  # Find all individual gene tables
  gene_files <- list.files(presence_absence_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  if (length(gene_files) == 0) {
    cat("⚠️  No individual gene tables found in:", presence_absence_dir, "\n")
    cat("   Skipping consolidated table generation\n\n")
    
  } else {
    cat("🔍 Found", length(gene_files), "individual gene table(s)\n")
    cat("🔄 Merging tables...\n\n")
    
    # Initialize merged table with metadata columns
    merged <- data.frame(
      Organism = metadata$Organism_clean,
      TaxIDs = metadata$TaxIDs,
      Kingdom = if("Kingdom" %in% colnames(metadata)) metadata$Kingdom else NA_character_,
      Phylum = if("Phylum" %in% colnames(metadata)) metadata$Phylum else NA_character_,
      Early = if("Early" %in% colnames(metadata)) metadata$Early else NA_character_,
      stringsAsFactors = FALSE
    )
    
    # Track organisms found in trees but not in metadata
    all_organisms <- metadata$Organism_clean
    new_organisms_list <- list()
    
    # Process each individual gene table
    for (i in seq_along(gene_files)) {
      gene_file <- gene_files[i]
      gene_name <- tools::file_path_sans_ext(basename(gene_file))
      
      cat("   Processing:", gene_name, "...")
      
      # Read gene table with proper column types
      gene_table <- tryCatch({
        read.table(
          gene_file,
          header = TRUE,
          sep = "\t",
          stringsAsFactors = FALSE,
          quote = "",
          comment.char = "",
          colClasses = c("character", "logical", rep("character", 3))
        )
      }, error = function(e) {
        cat(" ❌ ERROR\n")
        warning("   Could not read file: ", e$message)
        return(NULL)
      })
      
      if (is.null(gene_table)) next
      
      # Identify organisms in tree but not in metadata
      new_orgs <- setdiff(gene_table$Organism, all_organisms)
      if (length(new_orgs) > 0) {
        new_organisms_list[[gene_name]] <- new_orgs
      }
      
      # Add gene column to merged table
      merged[[gene_name]] <- FALSE  # Initialize with FALSE (absent)
      
      # Transfer presence values from gene table to merged table
      for (j in seq_len(nrow(gene_table))) {
        org <- gene_table$Organism[j]
        presence <- gene_table[[gene_name]][j]
        
        # Convert to logical if stored as character
        if (is.character(presence)) {
          presence <- as.logical(presence)
        }
        
        # Find organism in merged table and update
        idx <- which(merged$Organism == org)
        if (length(idx) > 0) {
          merged[[gene_name]][idx[1]] <- presence
        }
      }
      
      cat(" ✓\n")
    }
    
    # Handle organisms found in trees but not in metadata
    all_new_organisms <- unique(unlist(new_organisms_list))
    
    if (length(all_new_organisms) > 0) {
      cat("\n   📝 Adding", length(all_new_organisms), 
          "organism(s) found in trees but not in metadata:\n")
      
      for (org in all_new_organisms) {
        cat("      •", org, "\n")
      }
      
      # Create rows for new organisms
      new_rows <- data.frame(
        Organism = all_new_organisms,
        TaxIDs = NA_character_,
        Uniprot = NA_character_,
        Protein = NA_character_,
        Kingdom = "Unknown",
        Phylum = "Unknown",
        Early = NA_character_,
        stringsAsFactors = FALSE
      )
      
      # Initialize all gene columns to FALSE
      gene_cols_current <- setdiff(
        colnames(merged), 
        c("Organism", "TaxIDs", "Kingdom", "Phylum", "Early")
      )
      
      for (col in gene_cols_current) {
        new_rows[[col]] <- FALSE
      }
      
      # Re-read gene tables to populate presence values for new organisms
      for (gene_file in gene_files) {
        gene_name <- tools::file_path_sans_ext(basename(gene_file))
        gene_table <- read.table(
          gene_file, 
          header = TRUE, 
          sep = "\t", 
          stringsAsFactors = FALSE, 
          quote = "", 
          comment.char = "",
          colClasses = c("character", "logical", rep("character", 3))
        )
        
        for (k in seq_len(nrow(new_rows))) {
          org <- new_rows$Organism[k]
          gene_idx <- which(gene_table$Organism == org)
          
          if (length(gene_idx) > 0) {
            presence <- gene_table[[gene_name]][gene_idx[1]]
            if (is.character(presence)) {
              presence <- as.logical(presence)
            }
            new_rows[[gene_name]][k] <- presence
          }
        }
      }
      
      # Append new organisms to merged table
      merged <- rbind(merged, new_rows)
    }
    
    # Calculate N_Genes (number of genes each organism possesses)
    metadata_cols <- c("Organism", "TaxIDs", "Uniprot", "Protein", "Kingdom", "Phylum", "Early")
    gene_cols <- setdiff(colnames(merged), metadata_cols)
    
    # Ensure all gene columns are logical type
    for (col in gene_cols) {
      if (is.character(merged[[col]])) {
        merged[[col]] <- as.logical(merged[[col]])
      }
    }
    
    # Calculate gene count per organism
    merged$N_Genes <- rowSums(merged[, gene_cols, drop = FALSE], na.rm = TRUE)
    
    # Sort by gene count (descending) then by organism name
    merged <- merged[order(-merged$N_Genes, merged$Organism), ]
    
    # Export consolidated table
    write.table(
      merged, 
      file = merged_table_file, 
      sep = "\t", 
      row.names = FALSE, 
      quote = FALSE, 
      na = ""
    )
    
    #############################################
    ### CONSOLIDATED TABLE SUMMARY REPORT
    #############################################
    
    cat("\n╔═══════════════════════════════════════╗\n")
    cat("║   CONSOLIDATED TABLE SUMMARY          ║\n")
    cat("╚═══════════════════════════════════════╝\n\n")
    
    cat("📊 Table dimensions:\n")
    cat("   Total organisms:", nrow(merged), "\n")
    cat("   Total genes:", length(gene_cols), "\n")
    cat("   Mean genes per organism:", round(mean(merged$N_Genes), 2), "\n")
    cat("   Median genes per organism:", median(merged$N_Genes), "\n\n")
    
    cat("🧬 Gene coverage (organisms with gene present):\n")
    for (gene_col in gene_cols) {
      n_present <- sum(merged[[gene_col]], na.rm = TRUE)
      percentage <- round(n_present / nrow(merged) * 100, 1)
      cat("   ", sprintf("%-15s", gene_col), ": ", 
          sprintf("%3d", n_present), "/", nrow(merged), 
          " (", sprintf("%5.1f", percentage), "%)\n", sep = "")
    }
    
    cat("\n🌍 Organism distribution:\n")
    if ("Kingdom" %in% colnames(merged)) {
      kingdom_table <- table(merged$Kingdom)
      for (kingdom in names(sort(kingdom_table, decreasing = TRUE))) {
        cat("   ", sprintf("%-15s", kingdom), ": ", kingdom_table[kingdom], "\n", sep = "")
      }
    }
    
    cat("\n📁 Output files:\n")
    cat("   Location:", merged_table_file, "\n")
    cat("   Dimensions:", nrow(merged), "organisms ×", 
        length(gene_cols), "genes\n\n")
    
    cat("✅ Consolidated table created successfully!\n\n")
  }
}

#############################################
### FINAL BATCH PROCESSING REPORT
#############################################

cat("\n╔════════════════════════════════════════╗\n")
cat("║   BATCH PROCESSING COMPLETED           ║\n")
cat("╚════════════════════════════════════════╝\n\n")

cat("⏱️  Processing time:", round(as.numeric(elapsed_time), 2), "minutes\n\n")

cat("📊 Summary:\n")
cat("   Total attempts:", length(genes_to_process) * length(MODES), "\n")
cat("   ✅ Successful:", total_processed, "\n")
cat("   ❌ Errors:", total_errors, "\n\n")

# Display error log if there were failures
if (total_errors > 0) {
  cat("⚠️  Error log:\n")
  for (i in seq_along(error_log)) {
    cat("   ", i, ". ", error_log[i], "\n", sep = "")
  }
  cat("\n")
}

# Success message
if (total_errors == 0) {
  cat("🎉 All trees generated successfully!\n\n")
  
  cat("📁 Output locations:\n")
  cat("   Figures: phylo/", HORMONE, "/08_figures/{gene}/\n", sep = "")
  cat("   Individual tables: phylo/", HORMONE, "/09_presence-absence/\n", sep = "")
  cat("   Consolidated table: phylo/", HORMONE, "/merged-genes.tsv\n\n", sep = "")
  
} else {
  cat("⚠️  Processing completed with", total_errors, "error(s)\n")
  cat("   Review error log above for details\n\n")
}

cat("✓ Batch processing complete!\n")