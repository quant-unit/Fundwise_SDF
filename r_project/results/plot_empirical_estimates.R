# Plot Empirical Estimation Results
# ============================================================================
# Visualize estimated coefficients from empirical estimation across horizons,
# comparing EW/FW weighting and asymptotic/cross-validation inference.
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

#' Plot Empirical Estimation Results Across Horizons
#'
#' Creates a publication-quality multi-panel chart showing estimated coefficients
#' from empirical estimation results. Each column represents a different factor
#' model, each subplot shows 4 lines (EW/FW x asymptotic/cross-validation)
#' over estimation horizons.
#'
#' @param data_dir Path to the data_out_*** folder containing cache subdirectories
#' @param fund_type Character, fund type to plot: "PE", "VC", "BO", "RE", etc.
#' @param factors Character vector of factor names for the columns (e.g.
#'   c("MKT", "EG", "IA", "ME", "ROE")). "MKT" produces a single-factor model
#'   column (one plot row); other factors produce two-factor model columns
#'   (two plot rows: MKT + second factor).
#' @param export_svg Logical, whether to export as SVG file
#' @param export_png Logical, whether to export as PNG file
#' @param export_pdf Logical, whether to export as PDF file
#' @param export_csv Logical, whether to export combined plot data as CSV
#' @param export_latex Logical, whether to export a LaTeX summary table
#' @param output_file Path for output (extension added automatically)
#' @param width Plot width in inches (default: 14)
#' @param height Plot height in inches (default: 6)
#' @param png_dpi PNG resolution (default: 300)
#' @param y.max.mkt Optional numeric, maximum y-axis value for MKT row
#' @param y.min.mkt Optional numeric, minimum y-axis value for MKT row
#' @param y.lim.second Optional named list of y-axis limits for second factor row.
#'   Each element is a numeric vector c(min, max) keyed by factor name.
#'   Use "default" for factors not explicitly listed.
#'   Example: list(Alpha = c(-0.01, 0.01), default = c(-5, 5))
#' @param x.max Optional numeric, maximum x-axis value (in months)
#' @param x.min Optional numeric, minimum x-axis value (in months)
#'
#' @return A ggplot object (invisibly if exported)
#'
#' @examples
#' plot_empirical_estimates(
#'     data_dir = "data_out_2026-emp-max-vin-2019",
#'     fund_type = "PE",
#'     factors = c("MKT", "EG", "IA", "ME", "ROE")
#' )
#'
plot_empirical_estimates <- function(
    data_dir,
    fund_type = "PE",
    factors = c("MKT", "EG", "IA", "ME", "ROE"),
    export_svg = FALSE,
    export_png = FALSE,
    export_pdf = FALSE,
    export_csv = FALSE,
    export_latex = FALSE,
    output_file = "empirical_estimates_plot",
    width = 14,
    height = 6,
    png_dpi = 300,
    y.max.mkt = NULL,
    y.min.mkt = NULL,
    y.lim.second = NULL,
    x.max = NULL,
    x.min = NULL,
    v.lines = NULL,
    h.lines = NULL,
    v.colors = c("black", "black"),
    h.colors = c("black", "black"),
    main.linewidth = 0.8,
    abline.linewidth = 0.7,
    cex = 1.0) {
    # -------------------------------------------------------------------------
    # Define file paths for the 4 CSV sources
    # -------------------------------------------------------------------------
    csv_sources <- list(
        list(
            path = file.path(data_dir, "cache_q_factors_preqin_EW_VYP", "0_asymptotic_inference_summary.csv"),
            label = "EW Asymptotic",
            weighting = "EW",
            method = "Asymptotic"
        ),
        list(
            path = file.path(data_dir, "cache_q_factors_preqin_EW_VYP", "0_cross_validation_summary.csv"),
            label = "EW Cross-Validation",
            weighting = "EW",
            method = "Cross-Validation"
        ),
        list(
            path = file.path(data_dir, "cache_q_factors_preqin_FW_VYP", "0_asymptotic_inference_summary.csv"),
            label = "FW Asymptotic",
            weighting = "FW",
            method = "Asymptotic"
        ),
        list(
            path = file.path(data_dir, "cache_q_factors_preqin_FW_VYP", "0_cross_validation_summary.csv"),
            label = "FW Cross-Validation",
            weighting = "FW",
            method = "Cross-Validation"
        )
    )

    # -------------------------------------------------------------------------
    # Read and combine data from all 4 sources
    # -------------------------------------------------------------------------
    all_data_full <- lapply(csv_sources, function(src) {
        if (!file.exists(src$path)) {
            warning(paste("File not found:", src$path))
            return(NULL)
        }
        df <- read.csv(src$path, stringsAsFactors = FALSE)
        df$source <- src$label
        df$weighting <- src$weighting
        df$method <- src$method
        df
    })

    all_data_full <- bind_rows(Filter(Negate(is.null), all_data_full))

    # Subset for plotting (only the columns needed)
    all_data <- lapply(csv_sources, function(src) {
        if (!file.exists(src$path)) {
            return(NULL)
        }
        df <- read.csv(src$path, stringsAsFactors = FALSE)
        keep_cols <- c("Type", "max.month", "MKT", "Factor", "Coef")
        df <- df[, intersect(keep_cols, names(df)), drop = FALSE]
        df$source <- src$label
        df$weighting <- src$weighting
        df$method <- src$method
        df
    })

    all_data <- do.call(rbind, Filter(Negate(is.null), all_data))

    if (is.null(all_data) || nrow(all_data) == 0) {
        stop("No data could be loaded from the specified directory.")
    }

    # Filter by fund type
    plot_data <- all_data %>%
        filter(Type == fund_type)

    if (nrow(plot_data) == 0) {
        stop(paste(
            "No data found for fund type:", fund_type,
            "\nAvailable types:", paste(unique(all_data$Type), collapse = ", ")
        ))
    }

    # Convert max.month to numeric and compute horizon in years
    plot_data$max.month <- as.numeric(plot_data$max.month)
    plot_data$horizon_years <- plot_data$max.month / 12

    # Apply x-axis limits (filter data by horizon range)
    if (!is.null(x.min)) {
        plot_data <- plot_data %>% filter(max.month >= x.min)
    }
    if (!is.null(x.max)) {
        plot_data <- plot_data %>% filter(max.month <= x.max)
    }

    # Order sources for consistent legend
    source_levels <- c(
        "EW Asymptotic", "EW Cross-Validation",
        "FW Asymptotic", "FW Cross-Validation"
    )
    plot_data$source <- factor(plot_data$source, levels = source_levels)

    # -------------------------------------------------------------------------
    # Define publication-quality theme
    # -------------------------------------------------------------------------
    theme_publication <- theme_minimal(base_size = 11 * cex, base_family = "serif") +
        theme(
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
            panel.background = element_rect(fill = "white"),
            strip.text = element_text(face = "bold", size = 10 * cex, margin = margin(b = 5, t = 5)),
            strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.5),
            axis.title = element_text(face = "bold", size = 10 * cex),
            axis.text = element_text(size = 9 * cex, color = "grey20"),
            axis.ticks = element_line(color = "grey40", linewidth = 0.3),
            axis.line = element_blank(),
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 10 * cex),
            legend.text = element_text(size = 9 * cex),
            legend.key.size = unit(0.8 * cex, "cm"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.margin = margin(t = 5),
            plot.title = element_text(
                face = "bold", size = 13 * cex, hjust = 0.5,
                margin = margin(b = 10)
            ),
            plot.subtitle = element_text(
                size = 10 * cex, hjust = 0.5, color = "grey30",
                margin = margin(b = 15)
            ),
            plot.margin = margin(5, 5, 5, 5)
        )

    # -------------------------------------------------------------------------
    # Color and linetype palette
    # -------------------------------------------------------------------------
    # EW = green tones, FW = purple tones
    # Asymptotic = solid, Cross-Validation = dashed
    colors_source <- c(
        "EW Asymptotic"        = "#048072",
        "EW Cross-Validation"  = "#65E6D7",
        "FW Asymptotic"        = "#AB23CC",
        "FW Cross-Validation"  = "#E3A7F2"
    )

    linetypes_source <- c(
        "EW Asymptotic"        = "solid",
        "EW Cross-Validation"  = "dashed",
        "FW Asymptotic"        = "solid",
        "FW Cross-Validation"  = "dashed"
    )

    shapes_source <- c(
        "EW Asymptotic"        = 16,
        "EW Cross-Validation"  = 1,
        "FW Asymptotic"        = 17,
        "FW Cross-Validation"  = 2
    )

    # -------------------------------------------------------------------------
    # Determine which factors are single-factor vs two-factor
    # -------------------------------------------------------------------------
    has_second_factor <- any(factors != "MKT")
    single_factor_only <- all(factors == "MKT")

    # All factors define the columns (including MKT if present)
    # Both rows iterate over the same vector for alignment
    column_factors <- factors

    # -------------------------------------------------------------------------
    # Build MKT row: one subplot per column factor
    # -------------------------------------------------------------------------
    mkt_plots <- list()

    for (i in seq_along(column_factors)) {
        cf <- column_factors[i]

        if (cf == "MKT") {
            # Single-factor model: use rows where Factor == "MKT"
            mkt_df <- plot_data %>%
                filter(Factor == "MKT") %>%
                select(horizon_years, MKT, source) %>%
                rename(beta_MKT = MKT)
        } else {
            # Two-factor model: use rows where Factor == cf, take MKT column
            mkt_df <- plot_data %>%
                filter(Factor == cf) %>%
                select(horizon_years, MKT, source) %>%
                rename(beta_MKT = MKT)
        }

        p <- ggplot(mkt_df, aes(
            x = horizon_years, y = beta_MKT,
            color = source, linetype = source, shape = source
        )) +
            geom_line(linewidth = main.linewidth, alpha = 0.9) +
            geom_point(size = 2.5 * cex, fill = "white", stroke = 0.8) +
            scale_color_manual(values = colors_source, name = "Estimation Method") +
            scale_linetype_manual(values = linetypes_source, name = "Estimation Method") +
            scale_shape_manual(values = shapes_source, name = "Estimation Method") +
            scale_x_continuous(
                breaks = seq(0, max(mkt_df$horizon_years, na.rm = TRUE), by = 5),
                expand = c(0.02, 0)
            ) +
            labs(
                title = if (cf == "MKT") "MKT (Single Factor)" else paste0("MKT (", cf, " Model)"),
                x = if (has_second_factor && !single_factor_only) NULL else "Horizon (Years)",
                y = if (i == 1) expression(beta[MKT]) else NULL
            ) +
            {
                if (!is.null(y.max.mkt) || !is.null(y.min.mkt)) {
                    coord_cartesian(ylim = c(
                        if (is.null(y.min.mkt)) NA else y.min.mkt,
                        if (is.null(y.max.mkt)) NA else y.max.mkt
                    ))
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(v.lines)) {
                    lapply(seq_along(v.lines), function(idx) {
                        geom_vline(xintercept = v.lines[idx], color = v.colors[min(idx, length(v.colors))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines)) {
                    lapply(seq_along(h.lines), function(idx) {
                        geom_hline(yintercept = h.lines[idx], color = h.colors[min(idx, length(h.colors))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            theme_publication +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 10 * cex, hjust = 0.5, margin = margin(b = 5)),
                axis.title.y = if (i > 1) element_blank() else element_text(face = "bold", size = 10 * cex)
            )

        mkt_plots[[i]] <- p
    }

    # -------------------------------------------------------------------------
    # Build Second Factor row (only if we have non-MKT factors)
    # -------------------------------------------------------------------------
    second_plots <- list()
    first_coef_plot <- TRUE # Track first real coef plot for y-axis label

    if (has_second_factor && !single_factor_only) {
        for (i in seq_along(column_factors)) {
            cf <- column_factors[i]

            if (cf == "MKT") {
                # MKT column => blank spacer in row 2
                blank_plot <- ggplot() +
                    theme_void() +
                    theme(plot.margin = margin(5, 5, 5, 5))
                second_plots[[i]] <- blank_plot
            } else {
                coef_df <- plot_data %>%
                    filter(Factor == cf) %>%
                    select(horizon_years, Coef, source)

                p <- ggplot(coef_df, aes(
                    x = horizon_years, y = Coef,
                    color = source, linetype = source, shape = source
                )) +
                    geom_line(linewidth = main.linewidth, alpha = 0.9) +
                    geom_point(size = 2.5 * cex, fill = "white", stroke = 0.8) +
                    scale_color_manual(values = colors_source, name = "Estimation Method") +
                    scale_linetype_manual(values = linetypes_source, name = "Estimation Method") +
                    scale_shape_manual(values = shapes_source, name = "Estimation Method") +
                    scale_x_continuous(
                        breaks = seq(0, max(coef_df$horizon_years, na.rm = TRUE), by = 5),
                        expand = c(0.02, 0)
                    ) +
                    labs(
                        title = cf,
                        x = "Horizon (Years)",
                        y = if (first_coef_plot) expression(beta["Second"]) else NULL
                    ) +
                    {
                        # Per-factor y-limits: look up cf in y.lim.second, fall back to "default"
                        if (!is.null(y.lim.second)) {
                            ylim_vec <- if (!is.null(y.lim.second[[cf]])) {
                                y.lim.second[[cf]]
                            } else if (!is.null(y.lim.second[["default"]])) {
                                y.lim.second[["default"]]
                            } else {
                                NULL
                            }
                            if (!is.null(ylim_vec)) {
                                coord_cartesian(ylim = ylim_vec)
                            } else {
                                NULL
                            }
                        } else {
                            NULL
                        }
                    } +
                    {
                        if (!is.null(v.lines)) {
                            lapply(seq_along(v.lines), function(idx) {
                                geom_vline(xintercept = v.lines[idx], color = v.colors[min(idx, length(v.colors))], linetype = "dotted", linewidth = abline.linewidth)
                            })
                        } else {
                            NULL
                        }
                    } +
                    {
                        if (!is.null(h.lines)) {
                            lapply(seq_along(h.lines), function(idx) {
                                geom_hline(yintercept = h.lines[idx], color = h.colors[min(idx, length(h.colors))], linetype = "dotted", linewidth = abline.linewidth)
                            })
                        } else {
                            NULL
                        }
                    } +
                    theme_publication +
                    theme(
                        legend.position = "none",
                        plot.title = element_text(size = 10 * cex, hjust = 0.5, margin = margin(b = 5)),
                        axis.title.y = if (!first_coef_plot) element_blank() else element_text(face = "bold", size = 10 * cex)
                    )

                second_plots[[i]] <- p
                first_coef_plot <- FALSE
            }
        }
    }

    # -------------------------------------------------------------------------
    # Assemble combined plot using patchwork
    # -------------------------------------------------------------------------

    # Create row 1 (MKT)
    if (length(mkt_plots) == 1) {
        row1 <- mkt_plots[[1]]
    } else {
        row1 <- Reduce(`|`, mkt_plots)
    }

    if (length(second_plots) > 0) {
        # Create Panel A label as a dedicated text element
        panel_a_label <- wrap_elements(
            panel = textGrob(
                expression(bold("Panel A: Market Factor (" * beta[MKT] * ")")),
                x = 0, y = 0, hjust = 0, vjust = 0,
                gp = gpar(fontsize = 12 * cex, fontfamily = "serif", fontface = "bold")
            )
        )

        # Create row 2 (second factor)
        if (length(second_plots) == 1) {
            row2 <- second_plots[[1]]
        } else {
            row2 <- Reduce(`|`, second_plots)
        }

        # Create Panel B label as a dedicated text element
        panel_b_label <- wrap_elements(
            panel = textGrob(
                "Panel B: Second Factor",
                x = 0, y = 0, hjust = 0, vjust = 0,
                gp = gpar(fontsize = 12 * cex, fontfamily = "serif", fontface = "bold")
            )
        )

        # Stack: Panel A label, row1, Panel B label, row2
        combined_plot <- panel_a_label / row1 / panel_b_label / row2 +
            plot_layout(heights = c(0.1, 1, 0.1, 1))
    } else {
        # Single-factor only — no panel labels needed
        combined_plot <- row1
    }

    # -------------------------------------------------------------------------
    # Add shared legend at bottom
    # -------------------------------------------------------------------------
    # Collect legends and unify at bottom via patchwork
    combined_plot <- combined_plot +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")

    # -------------------------------------------------------------------------
    # Export or display
    # -------------------------------------------------------------------------
    any_export <- export_svg || export_png || export_pdf || export_csv || export_latex

    if (any_export) {
        base_file <- sub("\\.(svg|png|pdf)$", "", output_file, ignore.case = TRUE)

        output_dir_path <- dirname(base_file)
        if (!dir.exists(output_dir_path) && output_dir_path != "." && output_dir_path != "") {
            dir.create(output_dir_path, recursive = TRUE)
        }

        if (export_svg) {
            svg_file <- paste0(base_file, ".svg")
            tryCatch(
                {
                    ggsave(
                        filename = svg_file,
                        plot = combined_plot,
                        width = width,
                        height = height,
                        device = "svg"
                    )
                    message(paste("SVG exported to:", svg_file))
                },
                error = function(e) {
                    warning("SVG export failed. Falling back to PDF.")
                    pdf_file <- paste0(base_file, ".pdf")
                    ggsave(
                        filename = pdf_file,
                        plot = combined_plot,
                        width = width,
                        height = height,
                        device = "pdf"
                    )
                    message(paste("PDF exported to:", pdf_file))
                }
            )
        }

        if (export_png) {
            png_file <- paste0(base_file, ".png")
            ggsave(
                filename = png_file,
                plot = combined_plot,
                width = width,
                height = height,
                device = "png",
                dpi = png_dpi,
                bg = "white"
            )
            message(paste("PNG exported to:", png_file))
        }

        if (export_pdf) {
            pdf_file <- paste0(base_file, ".pdf")
            ggsave(
                filename = pdf_file,
                plot = combined_plot,
                width = width,
                height = height,
                device = "pdf"
            )
            message(paste("PDF exported to:", pdf_file))
        }

        # ---------------------------------------------------------------------
        # CSV Export: full data with all columns
        # ---------------------------------------------------------------------
        if (export_csv) {
            csv_export <- all_data_full %>% filter(Type == fund_type)
            csv_file <- paste0(base_file, "_data.csv")
            write.csv(csv_export, file = csv_file, row.names = FALSE)
            message(paste("Data CSV exported to:", csv_file))
        }

        # ---------------------------------------------------------------------
        # LaTeX Table Export: compact booktabs table per factor
        # ---------------------------------------------------------------------
        if (export_latex) {
            .export_latex_tables(
                all_data_full = all_data_full,
                fund_type = fund_type,
                factors = factors,
                base_file = base_file
            )
        }

        invisible(combined_plot)
    } else {
        print(combined_plot)
        invisible(combined_plot)
    }
}


# =============================================================================
# Helper: LaTeX table export
# =============================================================================

.export_latex_tables <- function(all_data_full, fund_type, factors, base_file) {
    # Filter to fund type
    ft_data <- all_data_full %>% filter(Type == fund_type)

    if (nrow(ft_data) == 0) {
        warning(paste("No data for fund type", fund_type, "- skipping LaTeX export."))
        return(invisible(NULL))
    }

    # Available horizons (sorted)
    horizons <- sort(unique(as.numeric(ft_data$max.month)))

    # Significance stars based on |t|
    sig_stars <- function(t_val) {
        t_val <- abs(t_val)
        ifelse(is.na(t_val), "",
            ifelse(t_val >= 3.29, "$^{***}$",
                ifelse(t_val >= 2.58, "$^{**}$",
                    ifelse(t_val >= 1.96, "$^{*}$", "")
                )
            )
        )
    }

    # Format numeric with rounding; return "--" for NA
    fmt <- function(x, digits = 3) {
        ifelse(is.na(x), "--", formatC(round(x, digits), format = "f", digits = digits))
    }
    fmt2 <- function(x) fmt(x, 2)
    fmt3 <- function(x) fmt(x, 3)

    # Accumulate all table lines across factors into one file
    all_lines <- c("% Auto-generated by plot_empirical_estimates.R")

    # For each factor, build a table
    for (fct in factors) {
        # Determine if this is single-factor (MKT) or two-factor
        is_single <- (fct == "MKT")

        # -- Gather data for each weighting x method combination --
        get_row_data <- function(w, m, factor_name) {
            sub <- ft_data %>%
                filter(
                    .data$weighting == w,
                    .data$method == m,
                    .data$Factor == factor_name
                ) %>%
                mutate(max.month = as.numeric(max.month))
            sub
        }

        ew_asym <- get_row_data("EW", "Asymptotic", if (is_single) "MKT" else fct)
        ew_cv <- get_row_data("EW", "Cross-Validation", if (is_single) "MKT" else fct)
        fw_asym <- get_row_data("FW", "Asymptotic", if (is_single) "MKT" else fct)
        fw_cv <- get_row_data("FW", "Cross-Validation", if (is_single) "MKT" else fct)

        # Build table lines
        lines <- character()

        # LaTeX preamble
        if (is_single) {
            caption_text <- paste0("Empirical Estimates: ", fund_type, " --- Single-Factor (MKT) Model")
        } else {
            caption_text <- paste0("Empirical Estimates: ", fund_type, " --- MKT + ", fct, " Model")
        }

        lines <- c(
            lines,
            "",
            "\\begin{table}[htbp]",
            "\\centering",
            "\\footnotesize",
            paste0("\\caption{", caption_text, "}"),
            paste0("\\label{tab:empirical_", tolower(fund_type), "_", tolower(fct), "}"),
            "\\begin{tabular}{l *{5}{r} *{5}{r}}",
            "\\toprule",
            " & \\multicolumn{5}{c}{Equal-Weighted (EW)} & \\multicolumn{5}{c}{Value-Weighted (FW)} \\\\",
            "\\cmidrule(lr){2-6} \\cmidrule(lr){7-11}",
            paste0(
                "Horizon & $\\hat{\\beta}$ & SE$_A$ & SE$_{CV}$ & $t_A$ & $t_{CV}$",
                " & $\\hat{\\beta}$ & SE$_A$ & SE$_{CV}$ & $t_A$ & $t_{CV}$ \\\\"
            ),
            "\\midrule"
        )

        # ----- Panel A: MKT -----
        if (!is_single) {
            lines <- c(
                lines,
                "\\multicolumn{11}{l}{\\textit{Panel A: Market Factor ($\\beta_{MKT}$)}} \\\\",
                "\\addlinespace[2pt]"
            )
        }

        for (h in horizons) {
            horizon_label <- paste0(round(h / 12, 1), " yr")

            # EW
            ew_a <- ew_asym %>% filter(max.month == h)
            ew_c <- ew_cv %>% filter(max.month == h)
            # FW
            fw_a <- fw_asym %>% filter(max.month == h)
            fw_c <- fw_cv %>% filter(max.month == h)

            # MKT values
            ew_beta <- if (nrow(ew_a) > 0) ew_a$MKT[1] else NA
            ew_se_a <- if (nrow(ew_a) > 0 && "SE.MKT" %in% names(ew_a)) ew_a$SE.MKT[1] else NA
            ew_se_cv <- if (nrow(ew_c) > 0 && "SE.MKT" %in% names(ew_c)) ew_c$SE.MKT[1] else NA
            ew_t_a <- if (nrow(ew_a) > 0 && "t.MKT" %in% names(ew_a)) ew_a$t.MKT[1] else NA
            ew_t_cv <- if (nrow(ew_c) > 0 && "t.MKT" %in% names(ew_c)) ew_c$t.MKT[1] else NA

            fw_beta <- if (nrow(fw_a) > 0) fw_a$MKT[1] else NA
            fw_se_a <- if (nrow(fw_a) > 0 && "SE.MKT" %in% names(fw_a)) fw_a$SE.MKT[1] else NA
            fw_se_cv <- if (nrow(fw_c) > 0 && "SE.MKT" %in% names(fw_c)) fw_c$SE.MKT[1] else NA
            fw_t_a <- if (nrow(fw_a) > 0 && "t.MKT" %in% names(fw_a)) fw_a$t.MKT[1] else NA
            fw_t_cv <- if (nrow(fw_c) > 0 && "t.MKT" %in% names(fw_c)) fw_c$t.MKT[1] else NA

            row <- paste0(
                horizon_label,
                " & ", fmt3(ew_beta), sig_stars(ew_t_a),
                " & ", fmt3(ew_se_a),
                " & ", fmt3(ew_se_cv),
                " & ", fmt2(ew_t_a),
                " & ", fmt2(ew_t_cv),
                " & ", fmt3(fw_beta), sig_stars(fw_t_a),
                " & ", fmt3(fw_se_a),
                " & ", fmt3(fw_se_cv),
                " & ", fmt2(fw_t_a),
                " & ", fmt2(fw_t_cv),
                " \\\\"
            )
            lines <- c(lines, row)
        }

        # ----- Panel B: Second Factor (if two-factor) -----
        if (!is_single) {
            lines <- c(
                lines,
                "\\addlinespace[4pt]",
                paste0("\\multicolumn{11}{l}{\\textit{Panel B: ", fct, " ($\\beta_{", fct, "}$)}} \\\\"),
                "\\addlinespace[2pt]"
            )

            for (h in horizons) {
                horizon_label <- paste0(round(h / 12, 1), " yr")

                ew_a <- ew_asym %>% filter(max.month == h)
                ew_c <- ew_cv %>% filter(max.month == h)
                fw_a <- fw_asym %>% filter(max.month == h)
                fw_c <- fw_cv %>% filter(max.month == h)

                # Second factor (Coef) values
                ew_beta <- if (nrow(ew_a) > 0) ew_a$Coef[1] else NA
                ew_se_a <- if (nrow(ew_a) > 0 && "SE.Coef" %in% names(ew_a)) ew_a$SE.Coef[1] else NA
                ew_se_cv <- if (nrow(ew_c) > 0 && "SE.Coef" %in% names(ew_c)) ew_c$SE.Coef[1] else NA
                ew_t_a <- if (nrow(ew_a) > 0 && "t.Coef" %in% names(ew_a)) ew_a$t.Coef[1] else NA
                ew_t_cv <- if (nrow(ew_c) > 0 && "t.Coef" %in% names(ew_c)) ew_c$t.Coef[1] else NA

                fw_beta <- if (nrow(fw_a) > 0) fw_a$Coef[1] else NA
                fw_se_a <- if (nrow(fw_a) > 0 && "SE.Coef" %in% names(fw_a)) fw_a$SE.Coef[1] else NA
                fw_se_cv <- if (nrow(fw_c) > 0 && "SE.Coef" %in% names(fw_c)) fw_c$SE.Coef[1] else NA
                fw_t_a <- if (nrow(fw_a) > 0 && "t.Coef" %in% names(fw_a)) fw_a$t.Coef[1] else NA
                fw_t_cv <- if (nrow(fw_c) > 0 && "t.Coef" %in% names(fw_c)) fw_c$t.Coef[1] else NA

                row <- paste0(
                    horizon_label,
                    " & ", fmt3(ew_beta), sig_stars(ew_t_a),
                    " & ", fmt3(ew_se_a),
                    " & ", fmt3(ew_se_cv),
                    " & ", fmt2(ew_t_a),
                    " & ", fmt2(ew_t_cv),
                    " & ", fmt3(fw_beta), sig_stars(fw_t_a),
                    " & ", fmt3(fw_se_a),
                    " & ", fmt3(fw_se_cv),
                    " & ", fmt2(fw_t_a),
                    " & ", fmt2(fw_t_cv),
                    " \\\\"
                )
                lines <- c(lines, row)
            }
        }

        # Table footer
        lines <- c(
            lines,
            "\\bottomrule",
            "\\end{tabular}",
            "\\par\\vspace{2pt}",
            "{\\scriptsize\\textit{Notes:} $\\hat{\\beta}$ = parameter estimate; SE$_A$ = asymptotic standard error; SE$_{CV}$ = cross-validation standard error; $t_A$/$t_{CV}$ = corresponding $t$-ratios. Significance: $^{*}$\\,$p<0.05$, $^{**}$\\,$p<0.01$, $^{***}$\\,$p<0.005$.}",
            "\\end{table}"
        )

        # Append this factor's table to the combined output
        all_lines <- c(all_lines, lines)
    }

    # Write all tables to a single file
    tex_file <- paste0(base_file, "_tables.tex")
    writeLines(all_lines, con = tex_file)
    message(paste("LaTeX tables exported to:", tex_file))

    invisible(NULL)
}


# =============================================================================
# Plot Max Vintage–Year Cutoff Analysis
# =============================================================================

#' Plot Empirical Estimates Across Max Vintage–Year Cutoffs
#'
#' Creates a multi-panel chart comparing estimated coefficients across different
#' maximum vintage-year cutoffs (e.g., 2011–2021). For each factor column, two
#' side-by-side panels (EW and FW) show 11 lines—one per cutoff—over horizons.
#' Uses only asymptotic inference results.
#'
#' @param data_dir  Path to the data_out folder containing the max-vintage
#'   subdirectories (e.g., cache_q_factors_preqin_EW_VYP_max_vin_2011, …).
#' @param fund_type Character, fund type to plot (default "PE").
#' @param factors   Character vector of factor names for the columns. "MKT"
#'   produces a single-factor column (one row); other factors produce two-factor
#'   columns (MKT row + second factor row).
#' @param vintages  Integer vector of vintage cutoff years (default 2011:2021).
#' @param export_svg, export_png, export_pdf  Logical, whether to export.
#' @param output_file Base path for output (extension added automatically).
#' @param width, height Plot dimensions in inches.
#' @param png_dpi   PNG resolution (default 300).
#' @param y.max.mkt, y.min.mkt Optional y-axis limits for the MKT row.
#' @param y.lim.second Optional named list of y-axis limits for the second
#'   factor row. Keyed by factor name; use "default" as fallback.
#' @param x.max, x.min Optional x-axis limits (in months).
#'
#' @return A ggplot object (invisibly if exported).
plot_max_vintage_cutoff <- function(
    data_dir,
    fund_type = "PE",
    factors = c("MKT", "Alpha", "EG", "IA", "ME", "ROE"),
    vintages = 2011:2021,
    nc_tag = "",
    export_svg = FALSE,
    export_png = FALSE,
    export_pdf = FALSE,
    export_csv = FALSE,
    export_latex = FALSE,
    output_file = "max_vintage_cutoff_plot",
    width = 14,
    height = 6,
    png_dpi = 300,
    y.max.mkt = NULL,
    y.min.mkt = NULL,
    y.lim.second = NULL,
    x.max = NULL,
    x.min = NULL,
    v.lines = NULL,
    h.lines = NULL,
    v.colors = c("black", "black"),
    h.colors = c("black", "black"),
    main.linewidth = 0.7,
    abline.linewidth = 0.7,
    cex = 1) {
    # -------------------------------------------------------------------------
    # 1. Load data from all vintage × weighting combinations
    # -------------------------------------------------------------------------
    all_rows <- list()

    for (v in vintages) {
        for (w in c("EW", "FW")) {
            folder <- file.path(
                data_dir,
                paste0("cache_q_factors_preqin_", w, "_VYP", nc_tag, "_max_vin_", v)
            )
            csv_path <- file.path(folder, "0_asymptotic_inference_summary.csv")

            if (!file.exists(csv_path)) {
                warning(paste("File not found (skipped):", csv_path))
                next
            }

            df <- read.csv(csv_path, stringsAsFactors = FALSE)
            df$vintage <- v
            df$weighting <- w
            all_rows[[length(all_rows) + 1L]] <- df
        }
    }

    all_data <- bind_rows(all_rows)

    if (nrow(all_data) == 0) {
        stop("No data could be loaded from the specified directory.")
    }

    # -------------------------------------------------------------------------
    # 2. Filter and prepare plot data
    # -------------------------------------------------------------------------
    plot_data <- all_data %>%
        filter(Type == fund_type)

    if (nrow(plot_data) == 0) {
        stop(paste(
            "No data found for fund type:", fund_type,
            "\nAvailable types:", paste(unique(all_data$Type), collapse = ", ")
        ))
    }

    plot_data$max.month <- as.numeric(plot_data$max.month)
    plot_data$horizon_years <- plot_data$max.month / 12
    plot_data$vintage <- factor(plot_data$vintage)

    # Apply x-axis limits
    if (!is.null(x.min)) plot_data <- plot_data %>% filter(max.month >= x.min)
    if (!is.null(x.max)) plot_data <- plot_data %>% filter(max.month <= x.max)

    # -------------------------------------------------------------------------
    # 3. Publication-quality theme
    # -------------------------------------------------------------------------
    theme_publication <- theme_minimal(base_size = 11, base_family = "serif") +
        theme(
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
            panel.background = element_rect(fill = "white"),
            strip.text = element_text(
                face = "bold", size = 10,
                margin = margin(b = 5, t = 5)
            ),
            strip.background = element_rect(
                fill = "grey95", color = "grey40",
                linewidth = 0.5
            ),
            axis.title = element_text(face = "bold", size = 10),
            axis.text = element_text(size = 9, color = "grey20"),
            axis.ticks = element_line(color = "grey40", linewidth = 0.3),
            axis.line = element_blank(),
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9),
            legend.key.size = unit(0.8, "cm"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.margin = margin(t = 5),
            plot.title = element_text(
                face = "bold", size = 10,
                hjust = 0.5, margin = margin(b = 5)
            ),
            plot.margin = margin(5, 5, 5, 5)
        )

    # -------------------------------------------------------------------------
    # 4. Color palette for vintage lines (viridis turbo, 11 levels)
    # -------------------------------------------------------------------------
    n_vintages <- length(levels(plot_data$vintage))

    # -------------------------------------------------------------------------
    # 5. Helper: build a single panel (one factor-coefficient × one weighting)
    # -------------------------------------------------------------------------
    .build_panel <- function(df, y_col, title, show_x, show_y_label,
                             y_label, ylim_vec = NULL) {
        p <- ggplot(df, aes(
            x     = horizon_years,
            y     = .data[[y_col]],
            color = vintage,
            group = vintage
        )) +
            geom_line(linewidth = main.linewidth, alpha = 0.85) +
            geom_point(size = 1.5 * cex, alpha = 0.85) +
            scale_color_viridis_d(
                option = "turbo",
                name   = "Max Vintage Year",
                guide  = guide_legend(nrow = 1)
            ) +
            scale_x_continuous(
                breaks = seq(0, max(df$horizon_years, na.rm = TRUE), by = 5),
                expand = c(0.02, 0)
            ) +
            labs(
                title = title,
                x     = if (show_x) "Horizon (Years)" else NULL,
                y     = if (show_y_label) y_label else NULL
            ) +
            {
                if (!is.null(v.lines)) {
                    lapply(seq_along(v.lines), function(idx) {
                        geom_vline(xintercept = v.lines[idx], color = v.colors[min(idx, length(v.colors))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines)) {
                    lapply(seq_along(h.lines), function(idx) {
                        geom_hline(yintercept = h.lines[idx], color = h.colors[min(idx, length(h.colors))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            theme_publication +
            theme(
                legend.position = "none",
                axis.title.y = if (!show_y_label) {
                    element_blank()
                } else {
                    element_text(face = "bold", size = 10 * cex)
                }
            )

        if (!is.null(ylim_vec)) {
            p <- p + coord_cartesian(ylim = ylim_vec)
        }
        p
    }

    # -------------------------------------------------------------------------
    # 6. Determine layout
    # -------------------------------------------------------------------------
    has_second_factor <- any(factors != "MKT")
    single_factor_only <- all(factors == "MKT")
    column_factors <- factors

    # -------------------------------------------------------------------------
    # 7. Build MKT row: 2 panels per factor column (EW + FW side-by-side)
    # -------------------------------------------------------------------------
    mkt_plots <- list()

    for (i in seq_along(column_factors)) {
        cf <- column_factors[i]

        for (w in c("EW", "FW")) {
            if (cf == "MKT") {
                mkt_df <- plot_data %>%
                    filter(Factor == "MKT", weighting == w)
            } else {
                mkt_df <- plot_data %>%
                    filter(Factor == cf, weighting == w)
            }

            title_text <- paste0(w, ": MKT")

            ylim_mkt <- NULL
            if (!is.null(y.max.mkt) || !is.null(y.min.mkt)) {
                ylim_mkt <- c(
                    if (is.null(y.min.mkt)) NA else y.min.mkt,
                    if (is.null(y.max.mkt)) NA else y.max.mkt
                )
            }

            show_x <- !has_second_factor || single_factor_only
            # y-label only for the first panel (EW of first column)
            show_y <- (i == 1 && w == "EW")

            p <- .build_panel(
                df            = mkt_df,
                y_col         = "MKT",
                title         = title_text,
                show_x        = show_x,
                show_y_label  = show_y,
                y_label       = expression(beta[MKT]),
                ylim_vec      = ylim_mkt
            )

            mkt_plots[[length(mkt_plots) + 1L]] <- p
        }
    }

    # -------------------------------------------------------------------------
    # 8. Build Second Factor row (only for non-MKT columns)
    # -------------------------------------------------------------------------
    second_plots <- list()
    first_coef_plot <- TRUE

    if (has_second_factor && !single_factor_only) {
        for (i in seq_along(column_factors)) {
            cf <- column_factors[i]

            if (cf == "MKT") {
                # Blank spacer for each of the 2 sub-panels (EW + FW)
                for (w in c("EW", "FW")) {
                    second_plots[[length(second_plots) + 1L]] <-
                        ggplot() +
                        theme_void() +
                        theme(plot.margin = margin(5, 5, 5, 5))
                }
            } else {
                for (w in c("EW", "FW")) {
                    coef_df <- plot_data %>%
                        filter(Factor == cf, weighting == w)

                    title_text <- paste0(w, ": ", cf)

                    # Per-factor y-limits
                    ylim_vec <- NULL
                    if (!is.null(y.lim.second)) {
                        ylim_vec <- if (!is.null(y.lim.second[[cf]])) {
                            y.lim.second[[cf]]
                        } else if (!is.null(y.lim.second[["default"]])) {
                            y.lim.second[["default"]]
                        }
                    }

                    show_y <- (first_coef_plot && w == "EW")

                    p <- .build_panel(
                        df            = coef_df,
                        y_col         = "Coef",
                        title         = title_text,
                        show_x        = TRUE,
                        show_y_label  = show_y,
                        y_label       = expression(beta["Second"]),
                        ylim_vec      = ylim_vec
                    )

                    second_plots[[length(second_plots) + 1L]] <- p
                }
                first_coef_plot <- FALSE
            }
        }
    }

    # -------------------------------------------------------------------------
    # 9. Assemble with patchwork
    # -------------------------------------------------------------------------
    n_cols <- length(column_factors) * 2 # 2 panels (EW + FW) per factor column

    # Row 1: MKT
    row1 <- Reduce(`|`, mkt_plots)

    if (length(second_plots) > 0) {
        # Panel labels
        panel_a_label <- wrap_elements(
            panel = textGrob(
                expression(bold("Panel A: Market Factor (" * beta[MKT] * ")")),
                x = 0, y = 0, hjust = 0, vjust = 0,
                gp = gpar(fontsize = 12 * cex, fontfamily = "serif", fontface = "bold")
            )
        )
        panel_b_label <- wrap_elements(
            panel = textGrob(
                "Panel B: Second Factor",
                x = 0, y = 0, hjust = 0, vjust = 0,
                gp = gpar(fontsize = 12 * cex, fontfamily = "serif", fontface = "bold")
            )
        )

        row2 <- Reduce(`|`, second_plots)

        combined_plot <- panel_a_label / row1 / panel_b_label / row2 +
            plot_layout(heights = c(0.1, 1, 0.1, 1))
    } else {
        combined_plot <- row1
    }

    # Shared legend at bottom
    combined_plot <- combined_plot +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")

    # -------------------------------------------------------------------------
    # 10. Export or display
    # -------------------------------------------------------------------------
    any_export <- export_svg || export_png || export_pdf || export_csv || export_latex

    if (any_export) {
        base_file <- sub("\\.(svg|png|pdf)$", "", output_file, ignore.case = TRUE)

        output_dir_path <- dirname(base_file)
        if (!dir.exists(output_dir_path) &&
            output_dir_path != "." && output_dir_path != "") {
            dir.create(output_dir_path, recursive = TRUE)
        }

        if (export_svg) {
            svg_file <- paste0(base_file, ".svg")
            tryCatch(
                {
                    ggsave(
                        filename = svg_file, plot = combined_plot,
                        width = width, height = height, device = "svg"
                    )
                    message(paste("SVG exported to:", svg_file))
                },
                error = function(e) {
                    warning("SVG export failed. Falling back to PDF.")
                    pdf_file <- paste0(base_file, ".pdf")
                    ggsave(
                        filename = pdf_file, plot = combined_plot,
                        width = width, height = height, device = "pdf"
                    )
                    message(paste("PDF exported to:", pdf_file))
                }
            )
        }

        if (export_png) {
            png_file <- paste0(base_file, ".png")
            ggsave(
                filename = png_file, plot = combined_plot,
                width = width, height = height, device = "png",
                dpi = png_dpi, bg = "white"
            )
            message(paste("PNG exported to:", png_file))
        }

        if (export_pdf) {
            pdf_file <- paste0(base_file, ".pdf")
            ggsave(
                filename = pdf_file, plot = combined_plot,
                width = width, height = height, device = "pdf"
            )
            message(paste("PDF exported to:", pdf_file))
        }

        # CSV Export: full data with all columns
        if (export_csv) {
            csv_export <- all_data %>% filter(Type == fund_type)
            csv_file <- paste0(base_file, "_data.csv")
            write.csv(csv_export, file = csv_file, row.names = FALSE)
            message(paste("Data CSV exported to:", csv_file))
        }

        # LaTeX Export: compact vintage-cutoff tables
        if (export_latex) {
            .export_vintage_latex_tables(
                all_data  = all_data,
                fund_type = fund_type,
                factors   = factors,
                vintages  = vintages,
                base_file = base_file
            )
        }

        invisible(combined_plot)
    } else {
        print(combined_plot)
        invisible(combined_plot)
    }
}


# =============================================================================
# Helper: LaTeX table export for vintage-cutoff analysis
# =============================================================================

.export_vintage_latex_tables <- function(all_data, fund_type, factors, vintages,
                                         base_file) {
    ft_data <- all_data %>%
        filter(Type == fund_type) %>%
        mutate(max.month = as.numeric(max.month))

    if (nrow(ft_data) == 0) {
        warning(paste(
            "No data for fund type", fund_type,
            "- skipping LaTeX export."
        ))
        return(invisible(NULL))
    }

    horizons <- sort(unique(ft_data$max.month))
    vintages <- sort(vintages)
    n_vin <- length(vintages)

    # Format helpers
    fmt <- function(x, digits = 3) {
        ifelse(is.na(x), "--",
            formatC(round(x, digits), format = "f", digits = digits)
        )
    }

    all_lines <- c("% Auto-generated by plot_max_vintage_cutoff()")

    for (fct in factors) {
        is_single <- (fct == "MKT")

        # Caption
        if (is_single) {
            cap <- paste0(
                "Vintage-Year Cutoff Sensitivity: ", fund_type,
                " --- Single-Factor (MKT) Model"
            )
        } else {
            cap <- paste0(
                "Vintage-Year Cutoff Sensitivity: ", fund_type,
                " --- MKT + ", fct, " Model"
            )
        }

        lines <- c(
            "",
            "\\begin{table}[htbp]",
            "\\centering",
            "\\scriptsize",
            paste0("\\caption{", cap, "}"),
            paste0(
                "\\label{tab:vin_", tolower(fund_type), "_",
                tolower(fct), "}"
            ),
            paste0(
                "\\begin{tabular}{l", paste(rep("r", n_vin), collapse = ""),
                "}"
            ),
            "\\toprule"
        )

        # Column headers: vintage years
        header <- paste0(
            "Horizon & ",
            paste(vintages, collapse = " & "), " \\\\"
        )
        lines <- c(lines, header, "\\midrule")

        # ---- EW panel ----
        lines <- c(
            lines,
            "\\multicolumn{1}{l}{\\textit{Equal-Weighted (EW)}} \\\\"
        )

        # Panel A: MKT (β_MKT)
        if (!is_single) {
            lines <- c(
                lines,
                "\\addlinespace[2pt]",
                paste0(
                    "\\multicolumn{", n_vin + 1,
                    "}{l}{\\textit{Panel A: $\\beta_{MKT}$}} \\\\"
                )
            )
        }

        for (h in horizons) {
            h_label <- paste0(round(h / 12, 1), " yr")
            vals <- sapply(vintages, function(v) {
                row <- ft_data %>%
                    filter(
                        max.month == h, weighting == "EW",
                        vintage == v,
                        Factor == if (is_single) "MKT" else fct
                    )
                if (nrow(row) > 0) fmt(row$MKT[1]) else "--"
            })
            lines <- c(lines, paste0(
                h_label, " & ",
                paste(vals, collapse = " & "),
                " \\\\"
            ))
        }

        # Panel B: Second factor (β_Coef) — EW
        if (!is_single) {
            lines <- c(
                lines,
                "\\addlinespace[4pt]",
                paste0(
                    "\\multicolumn{", n_vin + 1,
                    "}{l}{\\textit{Panel B: $\\beta_{", fct, "}$}} \\\\"
                ),
                "\\addlinespace[2pt]"
            )

            for (h in horizons) {
                h_label <- paste0(round(h / 12, 1), " yr")
                vals <- sapply(vintages, function(v) {
                    row <- ft_data %>%
                        filter(
                            max.month == h, weighting == "EW",
                            vintage == v, Factor == fct
                        )
                    if (nrow(row) > 0) fmt(row$Coef[1]) else "--"
                })
                lines <- c(lines, paste0(
                    h_label, " & ",
                    paste(vals, collapse = " & "),
                    " \\\\"
                ))
            }
        }

        # ---- FW panel ----
        lines <- c(
            lines,
            "\\addlinespace[6pt]",
            "\\midrule",
            "\\multicolumn{1}{l}{\\textit{Value-Weighted (FW)}} \\\\"
        )

        if (!is_single) {
            lines <- c(
                lines,
                "\\addlinespace[2pt]",
                paste0(
                    "\\multicolumn{", n_vin + 1,
                    "}{l}{\\textit{Panel A: $\\beta_{MKT}$}} \\\\"
                )
            )
        }

        for (h in horizons) {
            h_label <- paste0(round(h / 12, 1), " yr")
            vals <- sapply(vintages, function(v) {
                row <- ft_data %>%
                    filter(
                        max.month == h, weighting == "FW",
                        vintage == v,
                        Factor == if (is_single) "MKT" else fct
                    )
                if (nrow(row) > 0) fmt(row$MKT[1]) else "--"
            })
            lines <- c(lines, paste0(
                h_label, " & ",
                paste(vals, collapse = " & "),
                " \\\\"
            ))
        }

        if (!is_single) {
            lines <- c(
                lines,
                "\\addlinespace[4pt]",
                paste0(
                    "\\multicolumn{", n_vin + 1,
                    "}{l}{\\textit{Panel B: $\\beta_{", fct, "}$}} \\\\"
                ),
                "\\addlinespace[2pt]"
            )

            for (h in horizons) {
                h_label <- paste0(round(h / 12, 1), " yr")
                vals <- sapply(vintages, function(v) {
                    row <- ft_data %>%
                        filter(
                            max.month == h, weighting == "FW",
                            vintage == v, Factor == fct
                        )
                    if (nrow(row) > 0) fmt(row$Coef[1]) else "--"
                })
                lines <- c(lines, paste0(
                    h_label, " & ",
                    paste(vals, collapse = " & "),
                    " \\\\"
                ))
            }
        }

        # Footer
        lines <- c(
            lines,
            "\\bottomrule",
            "\\end{tabular}",
            "\\par\\vspace{2pt}",
            paste0(
                "{\\scriptsize\\textit{Notes:} Each cell shows the ",
                "estimated $\\hat{\\beta}$ (asymptotic) for the respective ",
                "maximum vintage-year cutoff. Horizons are in years.}"
            ),
            "\\end{table}"
        )

        all_lines <- c(all_lines, lines)
    }

    tex_file <- paste0(base_file, "_tables.tex")
    writeLines(all_lines, con = tex_file)
    message(paste("LaTeX vintage-cutoff tables exported to:", tex_file))

    invisible(NULL)
}

# =============================================================================
# Plot Combined MKT Max Vintage–Year Cutoff Analysis
# =============================================================================

#' Plot Empirical Estimates Across Max Vintage-Year Cutoffs (Combined NC Tags for MKT)
#'
#' @param data_dir Path to the data_out folder
#' @param fund_type Fund type to plot
#' @param vintages Vintage years
#' @param nc_tags Vector of tags to combine
#' @param nc_labels Headers for the different tags
plot_max_vintage_cutoff_combined_mkt <- function(
    data_dir,
    fund_type = "PE",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_svg = FALSE,
    export_png = FALSE,
    export_pdf = FALSE,
    export_csv = FALSE,
    output_file = "max_vintage_cutoff_combined_plot",
    width = 14,
    height = 4,
    png_dpi = 300,
    y.max.mkt = NULL,
    y.min.mkt = NULL,
    x.max = NULL,
    x.min = NULL,
    v.lines = NULL,
    h.lines = NULL,
    v.colors = c("black", "black"),
    h.colors = c("black", "black"),
    main.linewidth = 0.7,
    abline.linewidth = 0.7,
    cex = 1.0) {
    all_rows <- list()
    for (v in vintages) {
        for (i in seq_along(nc_tags)) {
            nt <- nc_tags[i]
            nl <- nc_labels[i]
            for (w in c("EW", "FW")) {
                folder <- file.path(data_dir, paste0("cache_q_factors_preqin_", w, "_VYP", nt, "_max_vin_", v))
                csv_path <- file.path(folder, "0_asymptotic_inference_summary.csv")
                if (!file.exists(csv_path)) {
                    warning(paste("File not found (skipped):", csv_path))
                    next
                }
                df <- read.csv(csv_path, stringsAsFactors = FALSE)
                df$vintage <- v
                df$weighting <- w
                df$nc_tag <- nt
                df$nc_label <- nl
                all_rows[[length(all_rows) + 1L]] <- df
            }
        }
    }
    all_data <- bind_rows(all_rows)
    if (nrow(all_data) == 0) stop("No data could be loaded from the specified directory.")

    plot_data <- all_data %>% filter(Type == fund_type)
    if (nrow(plot_data) == 0) stop("No data found for fund type.")

    plot_data$max.month <- as.numeric(plot_data$max.month)
    plot_data$horizon_years <- plot_data$max.month / 12
    plot_data$vintage <- factor(plot_data$vintage)

    if (!is.null(x.min)) plot_data <- plot_data %>% filter(max.month >= x.min)
    if (!is.null(x.max)) plot_data <- plot_data %>% filter(max.month <= x.max)

    theme_publication <- theme_minimal(base_size = 11 * cex, base_family = "serif") +
        theme(
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
            panel.background = element_rect(fill = "white"),
            strip.text = element_text(face = "bold", size = 10 * cex, margin = margin(b = 5, t = 5)),
            strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.5),
            axis.title = element_text(face = "bold", size = 10 * cex),
            axis.text = element_text(size = 9 * cex, color = "grey20"),
            axis.ticks = element_line(color = "grey40", linewidth = 0.3),
            axis.line = element_blank(),
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 10 * cex),
            legend.text = element_text(size = 9 * cex),
            legend.key.size = unit(0.8 * cex, "cm"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.margin = margin(t = 5),
            plot.title = element_text(face = "bold", size = 10 * cex, hjust = 0.5, margin = margin(b = 5)),
            plot.margin = margin(5, 5, 5, 5)
        )

    .build_panel <- function(df, y_col, title, show_y_label, y_label, ylim_vec = NULL) {
        p <- ggplot(df, aes(x = horizon_years, y = .data[[y_col]], color = vintage, group = vintage)) +
            geom_line(linewidth = main.linewidth, alpha = 0.85) +
            geom_point(size = 1.5 * cex, alpha = 0.85) +
            scale_color_viridis_d(option = "turbo", name = "Max Vintage Year", guide = guide_legend(nrow = 1)) +
            scale_x_continuous(breaks = seq(0, max(df$horizon_years, na.rm = TRUE), by = 5), expand = c(0.02, 0)) +
            labs(title = title, x = "Horizon (Years)", y = if (show_y_label) y_label else NULL) +
            {
                if (!is.null(v.lines)) {
                    lapply(seq_along(v.lines), function(idx) {
                        geom_vline(xintercept = v.lines[idx], color = v.colors[min(idx, length(v.colors))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines)) {
                    lapply(seq_along(h.lines), function(idx) {
                        geom_hline(yintercept = h.lines[idx], color = h.colors[min(idx, length(h.colors))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            theme_publication +
            theme(
                legend.position = "none",
                axis.title.y = if (!show_y_label) element_blank() else element_text(face = "bold", size = 10 * cex)
            )
        if (!is.null(ylim_vec)) {
            p <- p + coord_cartesian(ylim = ylim_vec)
        }
        p
    }

    mkt_plots <- list()
    for (i in seq_along(nc_tags)) {
        nt <- nc_tags[i]
        nl <- nc_labels[i]
        for (w in c("EW", "FW")) {
            mkt_df <- plot_data %>%
                filter(Factor == "MKT", weighting == w, nc_tag == nt)

            title_text <- paste0(nl, " (", w, ")")
            ylim_mkt <- NULL
            if (!is.null(y.max.mkt) || !is.null(y.min.mkt)) {
                ylim_mkt <- c(if (is.null(y.min.mkt)) NA else y.min.mkt, if (is.null(y.max.mkt)) NA else y.max.mkt)
            }
            show_y <- (i == 1 && w == "EW")

            p <- .build_panel(
                df = mkt_df, y_col = "MKT", title = title_text,
                show_y_label = show_y, y_label = expression(beta[MKT]), ylim_vec = ylim_mkt
            )
            mkt_plots[[length(mkt_plots) + 1L]] <- p
        }
    }

    combined_plot <- Reduce(`|`, mkt_plots) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")

    any_export <- export_svg || export_png || export_pdf || export_csv
    if (any_export) {
        base_file <- sub("\\.(svg|png|pdf)$", "", output_file, ignore.case = TRUE)
        output_dir_path <- dirname(base_file)
        if (!dir.exists(output_dir_path) && output_dir_path != "." && output_dir_path != "") {
            dir.create(output_dir_path, recursive = TRUE)
        }
        if (export_svg) {
            svg_file <- paste0(base_file, ".svg")
            tryCatch(
                {
                    ggsave(filename = svg_file, plot = combined_plot, width = width, height = height, device = "svg")
                    message(paste("SVG exported to:", svg_file))
                },
                error = function(e) {
                    warning("SVG export failed. Falling back to PDF.")
                    pdf_file <- paste0(base_file, ".pdf")
                    ggsave(filename = pdf_file, plot = combined_plot, width = width, height = height, device = "pdf")
                    message(paste("PDF exported to:", pdf_file))
                }
            )
        }
        if (export_png) {
            png_file <- paste0(base_file, ".png")
            ggsave(filename = png_file, plot = combined_plot, width = width, height = height, device = "png", dpi = png_dpi, bg = "white")
            message(paste("PNG exported to:", png_file))
        }
        if (export_pdf) {
            pdf_file <- paste0(base_file, ".pdf")
            ggsave(filename = pdf_file, plot = combined_plot, width = width, height = height, device = "pdf")
            message(paste("PDF exported to:", pdf_file))
        }
        if (export_csv) {
            csv_export <- all_data %>% filter(Type == fund_type)
            csv_file <- paste0(base_file, "_data.csv")
            write.csv(csv_export, file = csv_file, row.names = FALSE)
            message(paste("Data CSV exported to:", csv_file))
        }
        invisible(combined_plot)
    } else {
        print(combined_plot)
        invisible(combined_plot)
    }
}

# ============================================================================
# Example Usage (uncomment to run)
# ============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

file.folder <- "results/data_out_2026-emp-max-vin-2019"
file.folder <- "data_out_2026_02_26"
out.folder <- "figures2"

# # Two-factor model with all factors
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("Alpha", "EG", "IA", "ME", "ROE"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    y.max.mkt = 1.5, y.min.mkt = -0.25,
    y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/empirical_PE_all_factors")
)

# Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 4,
    y.max.mkt = 1.5, y.min.mkt = 0.5,
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/empirical_PE_MKT")
)

# Max vintage-year cutoff analysis (all factors)
plot_max_vintage_cutoff(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("Alpha", "EG", "IA", "ME", "ROE"),
    vintages = 2011:2021,
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    y.max.mkt = 2.5, y.min.mkt = 0,
    y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/max_vintage_PE_all_factors")
)

if (FALSE) {
    # Max vintage-year cutoff analysis (MKT only)
    plot_max_vintage_cutoff(
        data_dir = file.folder,
        fund_type = "PE",
        factors = c("MKT"),
        vintages = 2011:2021,
        export_pdf = TRUE,
        export_svg = TRUE,
        export_csv = TRUE,
        export_latex = TRUE,
        height = 3,
        y.max.mkt = 2.5, y.min.mkt = 0,
        output_file = paste0(out.folder, "/max_vintage_PE_MKT")
    )

    # Max vintage-year cutoff analysis (MKT only, NC50)
    plot_max_vintage_cutoff(
        data_dir = file.folder,
        fund_type = "PE",
        factors = c("MKT"),
        vintages = 2011:2021,
        nc_tag = "_NC50",
        export_pdf = TRUE,
        export_svg = TRUE,
        export_csv = TRUE,
        export_latex = TRUE,
        height = 3,
        y.max.mkt = 2.5, y.min.mkt = 0,
        output_file = paste0(out.folder, "/max_vintage_PE_NC50_MKT")
    )
} else {
    # Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
    plot_max_vintage_cutoff_combined_mkt(
        data_dir = file.folder,
        fund_type = "PE",
        vintages = 2011:2021,
        nc_tags = c("", "_NC50"),
        nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
        export_pdf = TRUE,
        export_svg = TRUE,
        export_csv = TRUE,
        width = 14,
        height = 4,
        y.max.mkt = 2.5, y.min.mkt = 0,
        output_file = paste0(out.folder, "/max_vintage_PE_MKT_combined_NC")
    )
}


# Max vintage-year cutoff analysis (all factors)
plot_max_vintage_cutoff(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("Alpha", "EG", "IA", "ME", "ROE"),
    vintages = 2011:2021,
    nc_tag = "_NC50",
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    y.max.mkt = 2.5, y.min.mkt = 0,
    y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/max_vintage_PE_NC50_all_factors")
)



### PRESENTATION

out.folder <- "figures3"


# Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 6,
    width = 8,
    y.max.mkt = 1.5, y.min.mkt = -0.1,
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    v.lines = c(0, 15),
    h.lines = c(0, 0.75),
    v.colors = c("blue", "orange"),
    h.colors = c("gray", "red"),
    main.linewidth = 2,
    abline.linewidth = 1.5,
    cex = 1,
    output_file = paste0(out.folder, "/empirical_PE_MKT_alines")
)

# Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    y.max.mkt = 1.5, y.min.mkt = 0.5,
    main.linewidth = 2,
    cex = 1.5,
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/empirical_PE_MKT_larger")
)
