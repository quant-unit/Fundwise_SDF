# =============================================================================
# CSV Directory Comparison Script for Research Reproducibility
# =============================================================================
# This script compares CSV files between two directories to verify that
# research results are reproducible. Files are matched by content (excluding
# date/datetime columns), not by filename.
#
# Author: Auto-generated
# Date: 2026-01-11
# =============================================================================

library(digest) # For content hashing

# =============================================================================
# Configuration
# =============================================================================

# Columns to exclude from comparison (case-insensitive matching)
DATE_COLUMN_PATTERNS <- c("datetime", "date", "timestamp", "time", "created", "modified")

# Tolerance for floating-point comparisons
NUMERIC_TOLERANCE <- 1e-3

# =============================================================================
# Helper Functions
# =============================================================================

#' Find all CSV files recursively in a directory
#'
#' @param dir_path Path to the directory to search
#' @return Character vector of full file paths to CSV files
find_csv_files <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        stop(paste("Directory does not exist:", dir_path))
    }

    files <- list.files(
        path = dir_path,
        pattern = "\\.csv$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
    )

    return(files)
}

#' Identify date/datetime columns to exclude from comparison
#'
#' @param df Data frame to analyze
#' @return Character vector of column names to exclude
identify_date_columns <- function(df) {
    col_names <- tolower(names(df))
    date_cols <- c()

    for (i in seq_along(col_names)) {
        # Check if column name matches date patterns
        for (pattern in DATE_COLUMN_PATTERNS) {
            if (grepl(pattern, col_names[i], ignore.case = TRUE)) {
                date_cols <- c(date_cols, names(df)[i])
                break
            }
        }
    }

    return(date_cols)
}

#' Remove date columns from a data frame
#'
#' @param df Data frame to process
#' @return Data frame with date columns removed
remove_date_columns <- function(df) {
    date_cols <- identify_date_columns(df)

    if (length(date_cols) > 0) {
        df <- df[, !(names(df) %in% date_cols), drop = FALSE]
    }

    return(df)
}

#' Sort a data frame by all columns for consistent comparison
#'
#' @param df Data frame to sort
#' @return Data frame with rows sorted
sort_dataframe_rows <- function(df) {
    if (nrow(df) == 0) {
        return(df)
    }

    # Convert all columns to character for consistent sorting
    df_char <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

    # Create a sort key by concatenating all columns
    sort_key <- do.call(paste, c(df_char, sep = "|"))

    # Sort and return original data frame in new order
    df[order(sort_key), , drop = FALSE]
}

#' Generate a content fingerprint for a data frame
#'
#' @param df Data frame to fingerprint
#' @return Character string hash representing the content
generate_fingerprint <- function(df) {
    # Sort columns alphabetically for consistent comparison
    df <- df[, order(names(df)), drop = FALSE]

    # Sort rows for order-independent comparison
    df <- sort_dataframe_rows(df)

    # Convert to character representation
    content_str <- paste(
        paste(names(df), collapse = "|"),
        nrow(df),
        digest(df, algo = "md5"),
        sep = "||"
    )

    return(content_str)
}

#' Compare two data frames for equality (within tolerance)
#'
#' @param df1 First data frame
#' @param df2 Second data frame
#' @param tolerance Numeric tolerance for floating-point comparisons
#' @return List with match status and details
compare_dataframes <- function(df1, df2, tolerance = NUMERIC_TOLERANCE) {
    result <- list(
        match = FALSE,
        reason = NULL,
        differences = NULL
    )

    # Check dimensions
    if (nrow(df1) != nrow(df2)) {
        result$reason <- sprintf("Row count differs: %d vs %d", nrow(df1), nrow(df2))
        return(result)
    }

    if (ncol(df1) != ncol(df2)) {
        result$reason <- sprintf("Column count differs: %d vs %d", ncol(df1), ncol(df2))
        return(result)
    }

    # Check column names (sorted)
    cols1 <- sort(names(df1))
    cols2 <- sort(names(df2))

    if (!identical(cols1, cols2)) {
        missing_in_2 <- setdiff(cols1, cols2)
        missing_in_1 <- setdiff(cols2, cols1)
        result$reason <- sprintf(
            "Column names differ. In file1 but not file2: %s. In file2 but not file1: %s",
            paste(missing_in_2, collapse = ", "),
            paste(missing_in_1, collapse = ", ")
        )
        return(result)
    }

    # Reorder columns to match
    df1 <- df1[, cols1, drop = FALSE]
    df2 <- df2[, cols2, drop = FALSE]

    # Sort rows for order-independent comparison
    df1 <- sort_dataframe_rows(df1)
    df2 <- sort_dataframe_rows(df2)

    # Compare values column by column
    differences <- list()

    for (col in cols1) {
        v1 <- df1[[col]]
        v2 <- df2[[col]]

        if (is.numeric(v1) && is.numeric(v2)) {
            # Numeric comparison with tolerance
            if (!isTRUE(all.equal(v1, v2, tolerance = tolerance, check.attributes = FALSE))) {
                diff_idx <- which(abs(v1 - v2) > tolerance | (is.na(v1) != is.na(v2)))
                if (length(diff_idx) > 0) {
                    differences[[col]] <- list(
                        rows = head(diff_idx, 5),
                        sample_values = data.frame(
                            row = head(diff_idx, 5),
                            value1 = v1[head(diff_idx, 5)],
                            value2 = v2[head(diff_idx, 5)]
                        )
                    )
                }
            }
        } else {
            # Character/factor comparison
            if (!identical(as.character(v1), as.character(v2))) {
                diff_idx <- which(as.character(v1) != as.character(v2))
                if (length(diff_idx) > 0) {
                    differences[[col]] <- list(
                        rows = head(diff_idx, 5),
                        sample_values = data.frame(
                            row = head(diff_idx, 5),
                            value1 = as.character(v1[head(diff_idx, 5)]),
                            value2 = as.character(v2[head(diff_idx, 5)])
                        )
                    )
                }
            }
        }
    }

    if (length(differences) > 0) {
        result$reason <- sprintf(
            "Values differ in %d column(s): %s",
            length(differences),
            paste(names(differences), collapse = ", ")
        )
        result$differences <- differences
        return(result)
    }

    result$match <- TRUE
    result$reason <- "Content matches"
    return(result)
}

#' Load and prepare a CSV file for comparison
#'
#' @param file_path Path to CSV file
#' @return List with data frame (excluding date columns) and metadata
load_csv_for_comparison <- function(file_path) {
    tryCatch(
        {
            df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
            date_cols <- identify_date_columns(df)
            df_clean <- remove_date_columns(df)

            list(
                success = TRUE,
                data = df_clean,
                original_cols = names(df),
                excluded_cols = date_cols,
                nrows = nrow(df),
                ncols_original = ncol(df),
                ncols_clean = ncol(df_clean),
                fingerprint = generate_fingerprint(df_clean)
            )
        },
        error = function(e) {
            list(
                success = FALSE,
                error = e$message
            )
        }
    )
}

# =============================================================================
# Main Comparison Function
# =============================================================================

#' Compare CSV files between two directories
#'
#' @param dir1 Path to first directory
#' @param dir2 Path to second directory
#' @param verbose Print progress messages
#' @return List with comparison results
compare_csv_directories <- function(dir1, dir2, verbose = TRUE) {
    if (verbose) cat("=== CSV Directory Comparison ===\n\n")

    # Find all CSV files
    if (verbose) cat("Searching for CSV files...\n")
    files1 <- find_csv_files(dir1)
    files2 <- find_csv_files(dir2)

    if (verbose) {
        cat(sprintf("  Directory 1: %s\n", dir1))
        cat(sprintf("    Found: %d CSV files\n", length(files1)))
        cat(sprintf("  Directory 2: %s\n", dir2))
        cat(sprintf("    Found: %d CSV files\n\n", length(files2)))
    }

    if (length(files1) == 0 || length(files2) == 0) {
        warning("One or both directories contain no CSV files")
        return(NULL)
    }

    # Load and prepare all files from both directories
    if (verbose) cat("Loading and preparing files from Directory 1...\n")
    data1 <- lapply(files1, function(f) {
        result <- load_csv_for_comparison(f)
        result$path <- f
        result$name <- basename(f)
        result
    })
    names(data1) <- files1

    if (verbose) cat("Loading and preparing files from Directory 2...\n")
    data2 <- lapply(files2, function(f) {
        result <- load_csv_for_comparison(f)
        result$path <- f
        result$name <- basename(f)
        result
    })
    names(data2) <- files2

    # Check for load errors
    errors1 <- Filter(function(x) !x$success, data1)
    errors2 <- Filter(function(x) !x$success, data2)

    if (length(errors1) > 0) {
        cat(sprintf("\nWarning: %d files failed to load from Directory 1\n", length(errors1)))
    }
    if (length(errors2) > 0) {
        cat(sprintf("\nWarning: %d files failed to load from Directory 2\n", length(errors2)))
    }

    # Filter to successfully loaded files
    data1 <- Filter(function(x) x$success, data1)
    data2 <- Filter(function(x) x$success, data2)

    # Group files by fingerprint (quick matching)
    if (verbose) cat("\nGrouping files by content fingerprint...\n")

    fingerprints1 <- sapply(data1, function(x) x$fingerprint)
    fingerprints2 <- sapply(data2, function(x) x$fingerprint)

    # Find matches
    matches <- list()
    unmatched_from_1 <- c()
    matched_from_2 <- c()

    if (verbose) cat("Matching files...\n\n")

    for (i in seq_along(data1)) {
        fp1 <- fingerprints1[i]
        file1 <- names(data1)[i]

        # Find potential matches by fingerprint
        potential_matches <- which(fingerprints2 == fp1)

        if (length(potential_matches) > 0) {
            # Verify with detailed comparison
            for (j in potential_matches) {
                if (!(names(data2)[j] %in% matched_from_2)) {
                    comparison <- compare_dataframes(data1[[i]]$data, data2[[j]]$data)

                    if (comparison$match) {
                        matches[[length(matches) + 1]] <- list(
                            file1 = file1,
                            file2 = names(data2)[j],
                            file1_name = data1[[i]]$name,
                            file2_name = data2[[j]]$name,
                            rows = data1[[i]]$nrows,
                            cols = data1[[i]]$ncols_clean,
                            excluded_cols = data1[[i]]$excluded_cols
                        )
                        matched_from_2 <- c(matched_from_2, names(data2)[j])
                        break
                    }
                }
            }
        }

        # If no match found, add to unmatched
        if (!any(sapply(matches, function(m) m$file1 == file1))) {
            unmatched_from_1 <- c(unmatched_from_1, file1)
        }
    }

    unmatched_from_2 <- setdiff(names(data2), matched_from_2)

    # Generate report
    results <- list(
        dir1 = dir1,
        dir2 = dir2,
        total_files_dir1 = length(files1),
        total_files_dir2 = length(files2),
        successfully_loaded_dir1 = length(data1),
        successfully_loaded_dir2 = length(data2),
        matched_count = length(matches),
        unmatched_from_dir1 = length(unmatched_from_1),
        unmatched_from_dir2 = length(unmatched_from_2),
        matches = matches,
        unmatched_files_dir1 = unmatched_from_1,
        unmatched_files_dir2 = unmatched_from_2
    )

    # Print summary
    if (verbose) {
        cat("=== COMPARISON RESULTS ===\n\n")
        cat(sprintf("Directory 1: %s\n", dir1))
        cat(sprintf("Directory 2: %s\n\n", dir2))

        cat("Summary:\n")
        cat(sprintf(
            "  Files in Dir1: %d (loaded: %d)\n",
            results$total_files_dir1, results$successfully_loaded_dir1
        ))
        cat(sprintf(
            "  Files in Dir2: %d (loaded: %d)\n",
            results$total_files_dir2, results$successfully_loaded_dir2
        ))
        cat(sprintf("  Matched pairs: %d\n", results$matched_count))
        cat(sprintf("  Unmatched from Dir1: %d\n", results$unmatched_from_dir1))
        cat(sprintf("  Unmatched from Dir2: %d\n\n", results$unmatched_from_dir2))

        if (results$matched_count > 0) {
            match_rate <- results$matched_count / min(
                results$successfully_loaded_dir1,
                results$successfully_loaded_dir2
            ) * 100
            cat(sprintf("Match rate: %.1f%%\n\n", match_rate))

            if (match_rate == 100) {
                cat("✓ FULL MATCH: All files from the smaller directory have matching content!\n")
                cat("  Research results appear to be REPRODUCIBLE.\n\n")
            } else if (match_rate >= 90) {
                cat("⚠ HIGH MATCH: Most files match, some discrepancies found.\n\n")
            } else {
                cat("✗ PARTIAL MATCH: Significant differences found between directories.\n\n")
            }
        }

        # Show some matches
        if (length(matches) > 0 && verbose) {
            cat("Sample matched pairs (first 5):\n")
            for (i in seq_len(min(5, length(matches)))) {
                m <- matches[[i]]
                cat(sprintf(
                    "  %d. %s <-> %s (%d rows, %d cols)\n",
                    i, m$file1_name, m$file2_name, m$rows, m$cols
                ))
            }
            cat("\n")
        }

        # Show unmatched files
        if (length(unmatched_from_1) > 0 && length(unmatched_from_1) <= 10) {
            cat("Unmatched files from Dir1:\n")
            for (f in unmatched_from_1) {
                cat(sprintf("  - %s\n", basename(f)))
            }
            cat("\n")
        } else if (length(unmatched_from_1) > 10) {
            cat(sprintf(
                "Unmatched files from Dir1: %d files (too many to list)\n\n",
                length(unmatched_from_1)
            ))
        }

        if (length(unmatched_from_2) > 0 && length(unmatched_from_2) <= 10) {
            cat("Unmatched files from Dir2:\n")
            for (f in unmatched_from_2) {
                cat(sprintf("  - %s\n", basename(f)))
            }
            cat("\n")
        } else if (length(unmatched_from_2) > 10) {
            cat(sprintf(
                "Unmatched files from Dir2: %d files (too many to list)\n\n",
                length(unmatched_from_2)
            ))
        }
    }

    return(invisible(results))
}

# =============================================================================
# Single File Comparison Function
# =============================================================================

#' Compare two specific CSV files
#'
#' @param file1 Path to first CSV file
#' @param file2 Path to second CSV file
#' @param verbose Print detailed comparison
#' @return Comparison result
compare_two_files <- function(file1, file2, verbose = TRUE) {
    if (verbose) cat(sprintf("Comparing:\n  File 1: %s\n  File 2: %s\n\n", file1, file2))

    data1 <- load_csv_for_comparison(file1)
    data2 <- load_csv_for_comparison(file2)

    if (!data1$success) {
        stop(sprintf("Failed to load file 1: %s", data1$error))
    }
    if (!data2$success) {
        stop(sprintf("Failed to load file 2: %s", data2$error))
    }

    if (verbose) {
        cat(sprintf(
            "File 1: %d rows, %d cols (excluded: %s)\n",
            data1$nrows, data1$ncols_clean,
            if (length(data1$excluded_cols) > 0) paste(data1$excluded_cols, collapse = ", ") else "none"
        ))
        cat(sprintf(
            "File 2: %d rows, %d cols (excluded: %s)\n\n",
            data2$nrows, data2$ncols_clean,
            if (length(data2$excluded_cols) > 0) paste(data2$excluded_cols, collapse = ", ") else "none"
        ))
    }

    result <- compare_dataframes(data1$data, data2$data)

    if (verbose) {
        if (result$match) {
            cat("✓ FILES MATCH (excluding date columns)\n")
        } else {
            cat(sprintf("✗ FILES DO NOT MATCH\n  Reason: %s\n", result$reason))

            if (!is.null(result$differences)) {
                cat("\nSample differences:\n")
                for (col in names(result$differences)) {
                    cat(sprintf("  Column '%s':\n", col))
                    print(result$differences[[col]]$sample_values)
                }
            }
        }
    }

    return(invisible(result))
}

# =============================================================================
# Example Usage / Test
# =============================================================================

# Uncomment the following to run a test with the example directories:

test_comparison_results <- function() {
    # Set your base path
    base_path <- "simulation/data_results_sim"

    dir1 <- file.path(base_path, "data_out_2020")
    dir1 <- file.path(base_path, "data_out_2025-2")
    dir2 <- file.path(base_path, "data_out_2025")

    results <- compare_csv_directories(dir1, dir2)

    return(results)
}

# Run the test:
results <- test_comparison_results()

test_comparison_prepared <- function() {
    # Set your base path
    base_path <- "simulation"

    dir1 <- file.path(base_path, "data_prepared_sim")
    dir2 <- file.path(base_path, "data_prepared_sim-test2026")

    results <- compare_csv_directories(dir1, dir2)

    return(results)
}

# Run the test:
results <- test_comparison_prepared()

cat("CSV Directory Comparison Script loaded successfully.\n")
cat("Use compare_csv_directories(dir1, dir2) to compare two directories.\n")
cat("Use compare_two_files(file1, file2) to compare two specific files.\n")
