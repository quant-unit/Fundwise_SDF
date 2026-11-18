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
source(here::here("helper/result_cache.R"))

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
    data_source = "preqin",
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
    v.lines.mkt = NULL,
    h.lines.mkt = NULL,
    v.colors.mkt = c("black", "black"),
    h.colors.mkt = c("black", "black"),
    v.lines.second = NULL,
    h.lines.second = NULL,
    v.colors.second = c("black", "black"),
    h.colors.second = c("black", "black"),
    main.linewidth = 0.8,
    abline.linewidth = 0.7,
    cex = 1.0,
    weighting_filter = c("EW", "FW"),
    exclude_cv = FALSE,
    region = NULL,
    error_function = "L2_Lasso") {
    # -------------------------------------------------------------------------
    # Define file paths for the 4 CSV sources
    # -------------------------------------------------------------------------
    region_tag <- if (!is.null(region) && region != "") paste0("_", region) else ""
    ef_tag <- if (!is.null(error_function) && error_function != "L2_Lasso") paste0("_", error_function) else ""

    csv_sources <- list(
        list(
            path = file.path(data_dir, paste0("cache_q_factors_", data_source, "_EW_VYP", region_tag, ef_tag), "0_asymptotic_inference_summary.csv"),
            label = "EW Asymptotic",
            weighting = "EW",
            method = "Asymptotic"
        ),
        list(
            path = file.path(data_dir, paste0("cache_q_factors_", data_source, "_EW_VYP", region_tag, ef_tag), "0_cross_validation_summary.csv"),
            label = "EW Cross-Validation",
            weighting = "EW",
            method = "Cross-Validation"
        ),
        list(
            path = file.path(data_dir, paste0("cache_q_factors_", data_source, "_FW_VYP", region_tag, ef_tag), "0_asymptotic_inference_summary.csv"),
            label = "FW Asymptotic",
            weighting = "FW",
            method = "Asymptotic"
        ),
        list(
            path = file.path(data_dir, paste0("cache_q_factors_", data_source, "_FW_VYP", region_tag, ef_tag), "0_cross_validation_summary.csv"),
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
        df <- read.result.summary(src$path)
        if (nrow(df) == 0) return(NULL)
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
        df <- read.result.summary(src$path)
        if (nrow(df) == 0) return(NULL)
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

    # Optionally filter by weighting
    if (!is.null(weighting_filter)) {
        plot_data <- plot_data %>%
            filter(weighting %in% weighting_filter)
    }

    # Optionally exclude cross-validation lines
    if (exclude_cv) {
        plot_data <- plot_data %>%
            filter(method != "Cross-Validation")
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
            legend.title = element_text(face = "bold", size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 6 else 9) * cex),
            legend.text = element_text(size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 5 else 8) * cex),
            legend.key.size = unit((if (!is.null(weighting_filter) && length(weighting_filter) == 1) 0.5 else 0.8) * cex, "cm"),
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
                breaks = seq(0, max(c(0, mkt_df$horizon_years), na.rm = TRUE), by = 5),
                expand = c(0.02, 0)
            ) +
            labs(
                title = if (cf == "MKT") paste0(fund_type, " - MKT (Single Factor)") else paste0(fund_type, " - MKT (", cf, " Model)"),
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
                if (!is.null(v.lines.mkt)) {
                    lapply(seq_along(v.lines.mkt), function(idx) {
                        geom_vline(xintercept = v.lines.mkt[idx], color = v.colors.mkt[min(idx, length(v.colors.mkt))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines.mkt)) {
                    lapply(seq_along(h.lines.mkt), function(idx) {
                        geom_hline(yintercept = h.lines.mkt[idx], color = h.colors.mkt[min(idx, length(h.colors.mkt))], linetype = "dotted", linewidth = abline.linewidth)
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

    # Determine y-axis label for second factor
    non_mkt_factors <- column_factors[column_factors != "MKT"]
    if (length(non_mkt_factors) == 1) {
        if (tolower(non_mkt_factors[1]) == "alpha") {
            sec_ylab <- expression(alpha)
        } else {
            sec_ylab <- bquote(beta[.(non_mkt_factors[1])])
        }
    } else {
        sec_ylab <- expression(beta["Second"])
    }

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
                        breaks = seq(0, max(c(0, coef_df$horizon_years), na.rm = TRUE), by = 5),
                        expand = c(0.02, 0)
                    ) +
                    labs(
                        title = paste0(fund_type, " - ", cf),
                        x = "Horizon (Years)",
                        y = if (first_coef_plot) sec_ylab else NULL
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
                        if (!is.null(v.lines.second)) {
                            lapply(seq_along(v.lines.second), function(idx) {
                                geom_vline(xintercept = v.lines.second[idx], color = v.colors.second[min(idx, length(v.colors.second))], linetype = "dotted", linewidth = abline.linewidth)
                            })
                        } else {
                            NULL
                        }
                    } +
                    {
                        if (!is.null(h.lines.second)) {
                            lapply(seq_along(h.lines.second), function(idx) {
                                geom_hline(yintercept = h.lines.second[idx], color = h.colors.second[min(idx, length(h.colors.second))], linetype = "dotted", linewidth = abline.linewidth)
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

# B9: provenance comment lines for every generated LaTeX table: generation
# time and commit, and the commit(s), scenario(s) and data folder(s) of the
# estimation results the table is built from (columns written by
# estim_model_optimized.R since B9; "unknown" for older cached results).
.table_provenance <- function(dat) {
    source(here::here("helper/result_reporting.R"), local = TRUE)
    result.table.provenance(dat)
}

# Inference robustness table (decisions of 17 November 2026): next to the main
# table, SE_A (Bartlett, S = ceil(1.3 sqrt V)), the leave-one-vintage-out
# jackknife SE and the effective number of vintages carrying the coefficient;
# S, V and b in the notes. SEs are suppressed for estimates on a bound of the
# admissible domain (as in the main table); SE_JK carries a double dagger when
# a jackknife refit hit a bound of the domain or left the basin of the
# full-sample estimate. Not produced when the results have no SEs (L1) or
# predate these columns.
.inference_robustness_table <- function(ft_data, fund_type, fct, is_single, horizons, fmt2, fmt3) {
    if (!all(c("hac.S", "hac.V", "hac.b") %in% names(ft_data))) return(character())
    source(here::here("helper/result_reporting.R"), local = TRUE) # boundary.mask()
    se.col <- if (is_single) "SE.MKT" else "SE.Coef"
    jk.col <- paste0(se.col, ".jk")
    # N_eff of the reported coefficient (influence-based since r5); older
    # results only have neff.units, which belongs to beta_MKT
    neff.col <- if (is_single) intersect(c("neff.MKT", "neff.units"), names(ft_data))[1] else intersect("neff.Coef", names(ft_data))[1]
    fac <- if (is_single) "MKT" else fct
    rows <- ft_data[ft_data$method == "Asymptotic" & ft_data$Factor == fac, , drop = FALSE]
    if (!(se.col %in% names(rows)) || all(is.na(rows[[se.col]]))) return(character())
    row.of <- function(w, h) ft_data[ft_data$weighting == w & ft_data$method == "Asymptotic" & ft_data$Factor == fac & as.numeric(ft_data$max.month) == h, , drop = FALSE]
    val <- function(d, col) if (nrow(d) == 0 || is.na(col) || !(col %in% names(d))) NA else d[[col]][1]
    bound <- function(d) nrow(d) > 0 && boundary.mask(d)[1]
    jk.flag <- function(d) {
        n <- sum(c(val(d, "jk.n.bound"), val(d, "jk.n.switch")), na.rm = TRUE)
        if (n > 0) "$^{\\ddagger}$" else ""
    }
    cells <- function(w, h) {
        d <- row.of(w, h)
        if (bound(d)) return(c("--", "--", fmt2(val(d, neff.col))))
        jk <- val(d, jk.col)
        c(fmt3(val(d, se.col)), if (is.na(jk)) "--" else paste0(fmt3(jk), jk.flag(d)), fmt2(val(d, neff.col)))
    }
    par.label <- if (is_single) "$\\hat{\\beta}$" else if (tolower(fct) == "alpha") "$\\hat{\\alpha}$" else paste0("$\\hat{\\beta}_{", fct, "}$")
    lines <- c("", "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
        paste0("\\caption{Inference robustness: ", fund_type, " --- ", if (is_single) "Single-Factor (MKT)" else paste0("MKT + ", fct), " Model, SEs of ", par.label, "}"),
        paste0("\\label{tab:inference_", tolower(fund_type), "_", tolower(fct), "}"),
        "\\begin{tabular}{l *{3}{r} *{3}{r}}", "\\toprule",
        " & \\multicolumn{3}{c}{Equal-Weighted (EW)} & \\multicolumn{3}{c}{Value-Weighted (FW)} \\\\",
        "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
        "Horizon & SE$_A$ & SE$_{JK}$ & $N_{\\mathrm{eff}}$ & SE$_A$ & SE$_{JK}$ & $N_{\\mathrm{eff}}$ \\\\", "\\midrule")
    for (h in horizons) {
        lines <- c(lines, paste0(round(h / 12, 1), " yr & ", paste(c(cells("EW", h), cells("FW", h)), collapse = " & "), " \\\\"))
    }
    sv <- unique(na.omit(rows[, c("hac.S", "hac.V", "hac.b"), drop = FALSE]))
    sv.text <- if (nrow(sv)) paste(sprintf("$S=%g$, $V=%g$, $b=%.3f$", sv$hac.S, sv$hac.V, sv$hac.b), collapse = "; ") else "--"
    c(lines, "\\bottomrule", "\\end{tabular}", "\\par\\vspace{2pt}",
        paste0("{\\scriptsize\\textit{Notes:} SE$_A$ = Bartlett (Newey--West) sandwich over vintage years, truncation $S=\\lceil 1.3\\sqrt{V}\\rceil$; SE$_{JK}$ = leave-one-vintage-out jackknife SE (assumes independent vintages; every refit on the same admissible domain); $^{\\ddagger}$\\,at least one jackknife refit lies on a bound of the admissible domain or in a different local minimum than the full-sample estimate. ",
            "$N_{\\mathrm{eff}}=(\\sum_v c_v^2)^2/\\sum_v c_v^4$, $c_v$ the influence contribution of vintage $v$ to ", par.label, " (the score of $v$ premultiplied by the inverse Hessian): the effective number of vintages carrying the estimate, a concentration measure, not degrees of freedom. ",
            "``--'' for estimates on a bound of the admissible domain. ", sv.text, ". With few effective vintages all standard errors are approximate.}"),
        "\\end{table}")
}

.export_latex_tables <- function(all_data_full, fund_type, factors, base_file) {
    source(here::here("helper/result_reporting.R"), local = TRUE)
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

    # B7: SR_CV = full-sample estimate / SD_CV (fold dispersion), as defined in
    # the text. The former column showed t.MKT of the CV summary, i.e. the
    # fold MEAN / SD_CV. Fold means are reported in 0_cross_validation_summary.csv.
    stability.ratio <- function(est, sd.cv) {
        ifelse(is.na(est) | is.na(sd.cv) | sd.cv == 0, NA, est / sd.cv)
    }

    # Stars from two-sided fixed-b p-values (analyze_result.R: p.MKT, p.Coef;
    # decision 2026-11-17); results without them fall back to normal |t| cutoffs
    sig_stars_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "$^{***}$", ifelse(p < 0.01, "$^{**}$", ifelse(p < 0.05, "$^{*}$", ""))))
    pick_p <- function(d, col) if (nrow(d) > 0 && col %in% names(d)) d[[col]][1] else NA

    # Estimates on a parameter bound (at.bound, B9): dagger instead of
    # significance stars, asymptotic SE and t not reported (review fix)
    on.bound <- function(d) nrow(d) > 0 && boundary.mask(d)[1]
    star.or.bound <- function(d, t_val, p_val = NA) {
        if (on.bound(d)) return("$^{\\dagger}$")
        if (!is.na(p_val)) sig_stars_p(p_val) else sig_stars(t_val)
    }
    na.if.bound <- function(d, cell) if (on.bound(d)) "--" else cell

    # L1 tables (output file *_L1): own labels; no inference (decision of
    # 17 November 2026: no valid inference established for the L1 criterion)
    is.l1 <- grepl("_L1$", basename(base_file))

    # Accumulate all table lines across factors into one file
    all_lines <- c("% Auto-generated by plot_empirical_estimates.R", .table_provenance(ft_data))

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
            paste0("\\label{tab:empirical_", tolower(fund_type), "_", tolower(fct), if (is.l1) "_l1", "}"),
            "\\begin{tabular}{l *{5}{r} *{5}{r}}",
            "\\toprule",
            " & \\multicolumn{5}{c}{Equal-Weighted (EW)} & \\multicolumn{5}{c}{Value-Weighted (FW)} \\\\",
            "\\cmidrule(lr){2-6} \\cmidrule(lr){7-11}",
            paste0(
                "Horizon & $\\hat{\\beta}$ & SE$_A$ & SD$_{CV}$ & $t_A$ & $\\mathrm{SR}_{CV}$",
                " & $\\hat{\\beta}$ & SE$_A$ & SD$_{CV}$ & $t_A$ & $\\mathrm{SR}_{CV}$ \\\\"
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
                " & ", fmt3(ew_beta), star.or.bound(ew_a, ew_t_a, pick_p(ew_a, "p.MKT")),
                " & ", na.if.bound(ew_a, fmt3(ew_se_a)),
                " & ", fmt3(ew_se_cv),
                " & ", na.if.bound(ew_a, fmt2(ew_t_a)),
                " & ", fmt2(stability.ratio(ew_beta, ew_se_cv)),
                " & ", fmt3(fw_beta), star.or.bound(fw_a, fw_t_a, pick_p(fw_a, "p.MKT")),
                " & ", na.if.bound(fw_a, fmt3(fw_se_a)),
                " & ", fmt3(fw_se_cv),
                " & ", na.if.bound(fw_a, fmt2(fw_t_a)),
                " & ", fmt2(stability.ratio(fw_beta, fw_se_cv)),
                " \\\\"
            )
            lines <- c(lines, row)
        }

        # ----- Panel B: Second Factor (if two-factor) -----
        if (!is_single) {
            fct_label <- if (tolower(fct) == "alpha") "$\\alpha$" else paste0("$\\beta_{", fct, "}$")
            
            lines <- c(
                lines,
                "\\addlinespace[4pt]",
                paste0("\\multicolumn{11}{l}{\\textit{Panel B: ", fct, " (", fct_label, ")}} \\\\"),
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
                    " & ", fmt3(ew_beta), star.or.bound(ew_a, ew_t_a, pick_p(ew_a, "p.Coef")),
                    " & ", na.if.bound(ew_a, fmt3(ew_se_a)),
                    " & ", fmt3(ew_se_cv),
                    " & ", na.if.bound(ew_a, fmt2(ew_t_a)),
                    " & ", fmt2(stability.ratio(ew_beta, ew_se_cv)),
                    " & ", fmt3(fw_beta), star.or.bound(fw_a, fw_t_a, pick_p(fw_a, "p.Coef")),
                    " & ", na.if.bound(fw_a, fmt3(fw_se_a)),
                    " & ", fmt3(fw_se_cv),
                    " & ", na.if.bound(fw_a, fmt2(fw_t_a)),
                    " & ", fmt2(stability.ratio(fw_beta, fw_se_cv)),
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
            paste0("{\\scriptsize\\textit{Notes:} $\\hat{\\beta}$ = parameter estimate; SE$_A$ = asymptotic standard error; SD$_{CV}$ = cross-validation fold dispersion (standard deviation across hold-out folds); $t_A$ = asymptotic $t$-ratio; $\\mathrm{SR}_{CV}$ = signed stability ratio (estimate / fold dispersion), so negative values reflect negative coefficients and $|\\mathrm{SR}_{CV}|$ measures the estimate relative to fold dispersion. $^{\\dagger}$\\,estimate on a bound of the admissible domain ($\\beta_{MKT}\\in[-1,4]$, $\\alpha\\in[-2\\%,2\\%]$ per month, monthly gross benchmark return $1+r_f+\\theta^\\top F_t\\ge 0.08$ over the estimation window); SE$_A$ and $t_A$ not valid and not reported. ",
            if (is.l1) "L1 criterion: no valid inference is established for this estimator (it sits at a vintage whose pricing error is zero), so SE$_A$, $t_A$ and significance markers are not reported; SD$_{CV}$ measures sensitivity to omitting blocks of vintages, not sampling uncertainty.}" else
            "SE$_A$: Bartlett (Newey--West) sandwich over vintage years with truncation $S=\\lceil 1.3\\sqrt{V}\\rceil$, $V$ the number of vintages (Lazarus, Lewis, Stock and Watson, 2018). Significance markers from two-sided Kiefer--Vogelsang fixed-$b$ $p$-values of $t_A$ at $b=S/V$: $^{*}$\\,$p<0.05$, $^{**}$\\,$p<0.01$, $^{***}$\\,$p<0.001$.}"),
            "\\end{table}"
        )


        # Append this factor's table to the combined output
        all_lines <- c(all_lines, lines, .inference_robustness_table(ft_data, fund_type, fct, is_single, horizons, fmt2, fmt3))
    }

    # Write all tables to a single file
    tex_file <- paste0(base_file, "_tables.tex")
    writeLines(all_lines, con = tex_file)
    message(paste("LaTeX tables exported to:", tex_file))

    invisible(NULL)
}


# =============================================================================
# Cross-fund pair-product (U_{v,H}) tables, Gate G2
# =============================================================================
# Decision of 18 November 2026 (revision plan, log entry "G2: the U_{v,H}
# estimator"): the U_{v,H} minimizer over the admissible domain, equal vintage
# weights (q_v = 1) as the main specification and fund-count weights
# (q_v = n_v) as robustness, for EW and FW fund weighting alike; estimates on a
# bound get a dagger and no SE or marker (D2). U results are full-sample fits
# only (no CV folds, no jackknife), so each cell is estimate and SE_A.
#
# Writes <output_file>_tables.tex (one table per model: MKT, or MKT plus one
# second factor) and <output_file>_data.csv (every printed cell).
export_u_pair_tables <- function(data_dir, fund_type, output_file,
                                 factors = c("MKT", "Alpha", "EG", "IA", "ME", "ROE"),
                                 region = "North America", data_source = "preqin") {
    source(here::here("helper/result_reporting.R"), local = TRUE)
    # Table columns, left to right: fund weighting (EW, FW) x vintage weights q_v
    cols <- data.frame(weighting = c("EW", "EW", "FW", "FW"), q = c("equal", "n.funds", "equal", "n.funds"))
    folder.of <- function(w, q) file.path(data_dir, paste0("cache_q_factors_", data_source, "_", w, "_", region, "_U_pair",
        if (q == "n.funds") "_qn.funds" else ""))
    dat <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) {
        f <- file.path(folder.of(cols$weighting[j], cols$q[j]), "0_asymptotic_inference_summary.csv")
        if (!file.exists(f)) stop("U summary not found: ", f)
        s <- read.result.summary(f)
        s <- s[s$Type == fund_type, , drop = FALSE]
        if (nrow(s) == 0) stop("No U results for ", fund_type, " in ", f)
        s$weighting <- cols$weighting[j]; s$q <- cols$q[j]
        s
    }))
    dat$max.month <- as.numeric(dat$max.month)
    horizons <- sort(unique(dat$max.month))

    fmt3 <- function(x) if (is.na(x)) "--" else formatC(round(x, 3), format = "f", digits = 3)
    stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "$^{***}$" else if (p < 0.01) "$^{**}$" else if (p < 0.05) "$^{*}$" else ""
    # one printed cell pair (estimate with marker, SE_A) and its CSV record
    cell <- function(model, par, w, q, h) {
        d <- dat[dat$Factor == model & dat$weighting == w & dat$q == q & dat$max.month == h, , drop = FALSE]
        if (nrow(d) == 0) return(list(tex = c("--", "--"), rec = NULL))
        if (nrow(d) > 1) stop("Duplicate U rows: ", fund_type, " ", model, " ", w, " ", q, " ", h)
        est <- d[[if (par == "MKT") "MKT" else "Coef"]]
        se <- d[[if (par == "MKT") "SE.MKT" else "SE.Coef"]]
        p <- d[[if (par == "MKT") "p.MKT" else "p.Coef"]]
        bound <- boundary.mask(d)
        list(tex = c(paste0(fmt3(est), if (bound) "$^{\\dagger}$" else stars(p)), if (bound) "--" else fmt3(se)),
             rec = data.frame(Type = fund_type, model = model, parameter = par, weighting = w, vintage.weights = q,
                              max.month = h, estimate = est, SE.A = if (bound) NA else se, p = if (bound) NA else p,
                              at.bound = bound, active.constraints = if ("active.constraints" %in% names(d)) d$active.constraints else NA,
                              hac.V = d$hac.V, hac.S = d$hac.S, hac.b = d$hac.b, git.commit = d$git.commit))
    }
    par.label <- function(par) if (par == "MKT") "$\\beta_{MKT}$" else if (tolower(par) == "alpha") "$\\alpha$" else paste0("$\\beta_{", par, "}$")

    lines <- c("% Auto-generated by plot_empirical_estimates.R (export_u_pair_tables)", .table_provenance(dat))
    recs <- list()
    for (model in factors) {
        pars <- if (model == "MKT") "MKT" else c("MKT", model)
        lines <- c(lines, "", "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
            paste0("\\caption{Cross-fund pair-product estimates ($U_{v,H}$): ", fund_type, " --- ",
                if (model == "MKT") "Single-Factor (MKT) Model" else paste0("MKT + ", model, " Model"), "}"),
            paste0("\\label{tab:upair_", tolower(fund_type), "_", tolower(model), "}"),
            "\\begin{tabular}{l *{4}{r} *{4}{r}}", "\\toprule",
            " & \\multicolumn{4}{c}{Equal-weighted funds (EW)} & \\multicolumn{4}{c}{Value-weighted funds (FW)} \\\\",
            "\\cmidrule(lr){2-5} \\cmidrule(lr){6-9}",
            " & \\multicolumn{2}{c}{$q_v=1$} & \\multicolumn{2}{c}{$q_v=n_v$} & \\multicolumn{2}{c}{$q_v=1$} & \\multicolumn{2}{c}{$q_v=n_v$} \\\\",
            "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}",
            "Horizon & Est. & SE$_A$ & Est. & SE$_A$ & Est. & SE$_A$ & Est. & SE$_A$ \\\\", "\\midrule")
        for (k in seq_along(pars)) {
            par <- pars[k]
            if (length(pars) > 1) {
                lines <- c(lines, if (k > 1) "\\addlinespace[4pt]",
                    paste0("\\multicolumn{9}{l}{\\textit{Panel ", LETTERS[k], ": ", if (par == "MKT") "Market Factor" else par,
                        " (", par.label(par), ")}} \\\\"),
                    "\\addlinespace[2pt]")
            }
            for (h in horizons) {
                cs <- lapply(seq_len(nrow(cols)), function(j) cell(model, par, cols$weighting[j], cols$q[j], h))
                recs <- c(recs, lapply(cs, `[[`, "rec"))
                lines <- c(lines, paste0(round(h / 12, 1), " yr & ", paste(unlist(lapply(cs, `[[`, "tex")), collapse = " & "), " \\\\"))
            }
        }
        sv <- unique(na.omit(dat[dat$Factor == model, c("hac.S", "hac.V", "hac.b"), drop = FALSE]))
        sv.text <- if (nrow(sv)) paste(sprintf("$S=%g$, $V=%g$, $b=%.3f$", sv$hac.S, sv$hac.V, sv$hac.b), collapse = "; ") else "--"
        lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\par\\vspace{2pt}",
            paste0("{\\scriptsize\\textit{Notes:} Estimates minimize $V^{-1}\\sum_v w_vU_{v,H}(\\theta)$ over the admissible domain; $U_{v,H}$ is the within-vintage cross-fund pair-product criterion (the squared pooled pricing error without each fund's own squared error), unbiased for the squared shared vintage pricing error for any $n_v\\ge2$; vintages with fewer than two funds are dropped. ",
                "EW/FW: fund weighting within a vintage (equal, or by fund size). $q_v$: vintage weights, $w_v=Vq_v/\\sum_uq_u$; $q_v=1$ is the main specification, $q_v=n_v$ (the vintage's fund count) a robustness check. ",
                "SE$_A$: Bartlett (Newey--West) sandwich over vintage years with truncation $S=\\lceil 1.3\\sqrt{V}\\rceil$; markers from two-sided Kiefer--Vogelsang fixed-$b$ $p$-values at $b=S/V$: $^{*}$\\,$p<0.05$, $^{**}$\\,$p<0.01$, $^{***}$\\,$p<0.001$. ",
                "$^{\\dagger}$\\,estimate on a bound of the admissible domain ($\\beta_{MKT}\\in[-1,4]$, $\\alpha\\in[-2\\%,2\\%]$ per month, monthly gross benchmark return $1+r_f+\\theta^\\top F_t\\ge 0.08$ over the estimation window); SE$_A$ and marker not reported. ",
                "Full-sample fits only (no cross-validation or jackknife for this criterion). ", sv.text, ".}"),
            "\\end{table}")
    }
    writeLines(lines, con = paste0(output_file, "_tables.tex"))
    write.csv(do.call(rbind, Filter(Negate(is.null), recs)), paste0(output_file, "_data.csv"), row.names = FALSE)
    message("U tables exported to: ", output_file, "_tables.tex")
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
    data_source = "preqin",
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
    v.lines.mkt = NULL,
    h.lines.mkt = NULL,
    v.colors.mkt = c("black", "black"),
    h.colors.mkt = c("black", "black"),
    v.lines.second = NULL,
    h.lines.second = NULL,
    v.colors.second = c("black", "black"),
    h.colors.second = c("black", "black"),
    main.linewidth = 0.7,
    abline.linewidth = 0.7,
    cex = 1,
    weighting_filter = c("EW", "FW"),
    region = NULL) {
    # -------------------------------------------------------------------------
    # 1. Load data from all vintage × weighting combinations
    # -------------------------------------------------------------------------
    region_tag <- if (!is.null(region) && region != "") paste0("_", region) else ""
    all_rows <- list()

    for (v in vintages) {
        for (w in c("EW", "FW")) {
            folder <- file.path(
                data_dir,
                paste0("cache_q_factors_", data_source, "_", w, "_VYP", region_tag, nc_tag, "_max_vin_", v)
            )
            csv_path <- file.path(folder, "0_asymptotic_inference_summary.csv")

            if (!file.exists(csv_path)) {
                warning(paste("File not found (skipped):", csv_path))
                next
            }

            df <- read.result.summary(csv_path)
            if (nrow(df) == 0) next
            df$vintage <- v
            df$weighting <- w
            df$method <- "Asymptotic"
            all_rows[[length(all_rows) + 1L]] <- df
        }
    }

    all_data <- bind_rows(all_rows)

    if (nrow(all_data) == 0) {
        stop("No data could be loaded from the specified directory.")
    }

    if (!is.null(weighting_filter)) {
        all_data <- all_data %>% filter(weighting %in% weighting_filter)
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
            legend.title = element_text(face = "bold", size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 6 else 9)),
            legend.text = element_text(size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 5 else 8)),
            legend.key.size = unit((if (!is.null(weighting_filter) && length(weighting_filter) == 1) 0.5 else 0.8), "cm"),
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
                             y_label, ylim_vec = NULL,
                             v.lines_vec = NULL, v.colors_vec = NULL,
                             h.lines_vec = NULL, h.colors_vec = NULL) {
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
                breaks = seq(0, max(c(0, df$horizon_years), na.rm = TRUE), by = 5),
                expand = c(0.02, 0)
            ) +
            labs(
                title = title,
                x     = if (show_x) "Horizon (Years)" else NULL,
                y     = if (show_y_label) y_label else NULL
            ) +
            {
                if (!is.null(v.lines_vec)) {
                    lapply(seq_along(v.lines_vec), function(idx) {
                        geom_vline(xintercept = v.lines_vec[idx], color = v.colors_vec[min(idx, length(v.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines_vec)) {
                    lapply(seq_along(h.lines_vec), function(idx) {
                        geom_hline(yintercept = h.lines_vec[idx], color = h.colors_vec[min(idx, length(h.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
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

        for (w in weighting_filter) {
            if (cf == "MKT") {
                mkt_df <- plot_data %>%
                    filter(Factor == "MKT", weighting == w)
            } else {
                mkt_df <- plot_data %>%
                    filter(Factor == cf, weighting == w)
            }

            title_text <- paste0(fund_type, " - ", w, ": MKT")

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
                ylim_vec      = ylim_mkt,
                v.lines_vec   = v.lines.mkt,
                v.colors_vec  = v.colors.mkt,
                h.lines_vec   = h.lines.mkt,
                h.colors_vec  = h.colors.mkt
            )

            mkt_plots[[length(mkt_plots) + 1L]] <- p
        }
    }

    # -------------------------------------------------------------------------
    # 8. Build Second Factor row (only for non-MKT columns)
    # -------------------------------------------------------------------------
    second_plots <- list()
    first_coef_plot <- TRUE

    # Determine y-axis label for second factor dynamically
    non_mkt_factors <- column_factors[column_factors != "MKT"]
    if (length(non_mkt_factors) == 1) {
        if (tolower(non_mkt_factors[1]) == "alpha") {
            sec_ylab <- expression(alpha)
        } else {
            sec_ylab <- bquote(beta[.(non_mkt_factors[1])])
        }
    } else {
        sec_ylab <- expression(beta["Second"])
    }

    if (has_second_factor && !single_factor_only) {
        for (i in seq_along(column_factors)) {
            cf <- column_factors[i]

            if (cf == "MKT") {
                # Blank spacer for each sub-panel based on weighting
                for (w in weighting_filter) {
                    second_plots[[length(second_plots) + 1L]] <-
                        ggplot() +
                        theme_void() +
                        theme(plot.margin = margin(5, 5, 5, 5))
                }
            } else {
                for (w in weighting_filter) {
                    coef_df <- plot_data %>%
                        filter(Factor == cf, weighting == w)

                    title_text <- paste0(fund_type, " - ", w, ": ", cf)

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
                        y_label       = sec_ylab,
                        ylim_vec      = ylim_vec,

                        v.lines_vec   = v.lines.second,
                        v.colors_vec  = v.colors.second,
                        h.lines_vec   = h.lines.second,
                        h.colors_vec  = h.colors.second
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
    source(here::here("helper/result_reporting.R"), local = TRUE)
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

    estimate.cell <- function(row, column) {
        if (nrow(row) == 0) return("--")
        paste0(fmt(row[[column]][1]), boundary.marker(row)[1])
    }

    all_lines <- c("% Auto-generated by plot_max_vintage_cutoff()", .table_provenance(all_data))

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
                estimate.cell(row, "MKT")
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
                    estimate.cell(row, "Coef")
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
                estimate.cell(row, "MKT")
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
                    estimate.cell(row, "Coef")
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
                "maximum vintage-year cutoff. Horizons are in years. ",
                "$^{\\dagger}$ denotes an estimate on a parameter bound; asymptotic inference is not valid there.}"
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
    data_source = "preqin",
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
    v.lines.mkt = NULL,
    h.lines.mkt = NULL,
    v.colors.mkt = c("black", "black"),
    h.colors.mkt = c("black", "black"),
    main.linewidth = 0.7,
    abline.linewidth = 0.7,
    cex = 1.0,
    weighting_filter = c("EW", "FW"),
    region = NULL) {
    region_tag <- if (!is.null(region) && region != "") paste0("_", region) else ""
    all_rows <- list()
    for (v in vintages) {
        for (i in seq_along(nc_tags)) {
            nt <- nc_tags[i]
            nl <- nc_labels[i]
            for (w in c("EW", "FW")) {
                folder <- file.path(data_dir, paste0("cache_q_factors_", data_source, "_", w, "_VYP", region_tag, nt, "_max_vin_", v))
                csv_path <- file.path(folder, "0_asymptotic_inference_summary.csv")
                if (!file.exists(csv_path)) {
                    warning(paste("File not found (skipped):", csv_path))
                    next
                }
                df <- read.result.summary(csv_path)
                if (nrow(df) == 0) next
                df$vintage <- v
                df$weighting <- w
                df$nc_tag <- nt
                df$nc_label <- nl
                df$method <- "Asymptotic"
                all_rows[[length(all_rows) + 1L]] <- df
            }
        }
    }
    all_data <- bind_rows(all_rows)
    if (nrow(all_data) == 0) stop("No data could be loaded from the specified directory.")

    if (!is.null(weighting_filter)) {
        all_data <- all_data %>% filter(weighting %in% weighting_filter)
    }

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
            legend.title = element_text(face = "bold", size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 6 else 9) * cex),
            legend.text = element_text(size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 5 else 8) * cex),
            legend.key.size = unit((if (!is.null(weighting_filter) && length(weighting_filter) == 1) 0.5 else 0.8) * cex, "cm"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.margin = margin(t = 5),
            plot.title = element_text(face = "bold", size = 10 * cex, hjust = 0.5, margin = margin(b = 5)),
            plot.margin = margin(5, 5, 5, 5)
        )

    .build_panel <- function(df, y_col, title, show_y_label, y_label, ylim_vec = NULL,
                             v.lines_vec = NULL, v.colors_vec = NULL,
                             h.lines_vec = NULL, h.colors_vec = NULL) {
        p <- ggplot(df, aes(x = horizon_years, y = .data[[y_col]], color = vintage, group = vintage)) +
            geom_line(linewidth = main.linewidth, alpha = 0.85) +
            geom_point(size = 1.5 * cex, alpha = 0.85) +
            scale_color_viridis_d(option = "turbo", name = "Max Vintage Year", guide = guide_legend(nrow = 1)) +
            scale_x_continuous(breaks = seq(0, max(c(0, df$horizon_years), na.rm = TRUE), by = 5), expand = c(0.02, 0)) +
            labs(title = title, x = "Horizon (Years)", y = if (show_y_label) y_label else NULL) +
            {
                if (!is.null(v.lines_vec)) {
                    lapply(seq_along(v.lines_vec), function(idx) {
                        geom_vline(xintercept = v.lines_vec[idx], color = v.colors_vec[min(idx, length(v.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines_vec)) {
                    lapply(seq_along(h.lines_vec), function(idx) {
                        geom_hline(yintercept = h.lines_vec[idx], color = h.colors_vec[min(idx, length(h.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
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
        for (w in weighting_filter) {
            mkt_df <- plot_data %>%
                filter(Factor == "MKT", weighting == w, nc_tag == nt)

            title_text <- paste0(fund_type, ": ", nl, " (", w, ")")
            ylim_mkt <- NULL
            if (!is.null(y.max.mkt) || !is.null(y.min.mkt)) {
                ylim_mkt <- c(if (is.null(y.min.mkt)) NA else y.min.mkt, if (is.null(y.max.mkt)) NA else y.max.mkt)
            }
            show_y <- (i == 1 && w == "EW")

            p <- .build_panel(
                df = mkt_df, y_col = "MKT", title = title_text,
                show_y_label = show_y, y_label = expression(beta[MKT]), ylim_vec = ylim_mkt,
                v.lines_vec = v.lines.mkt, v.colors_vec = v.colors.mkt,
                h.lines_vec = h.lines.mkt, h.colors_vec = h.colors.mkt
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

# =============================================================================
# Plot Combined Alpha + MKT Max Vintage–Year Cutoff Analysis (Two-Panel)
# =============================================================================

#' Plot Two-Factor (Alpha + MKT) Estimates Across Max Vintage-Year Cutoffs
#'
#' Creates a two-row, multi-column chart comparing estimated coefficients from
#' the two-factor Alpha + MKT model across different maximum vintage-year
#' cutoffs (e.g., 2011–2021) and NAV treatment specifications.
#' Top row:  β_MKT   (Panel A)
#' Bottom row: α      (Panel B)
#' Each column corresponds to one nc_tag × weighting combination.
#'
#' @param data_dir   Path to the data_out folder
#' @param fund_type  Fund type to plot
#' @param vintages   Vintage years
#' @param nc_tags    Vector of NC tags to combine
#' @param nc_labels  Human-readable headers for the different tags
#' @param y.max.mkt, y.min.mkt  Optional y-axis limits for the MKT (top) row
#' @param y.max.alpha, y.min.alpha  Optional y-axis limits for the Alpha (bottom) row
plot_max_vintage_cutoff_combined_alpha_mkt <- function(
    data_dir,
    data_source = "preqin",
    fund_type = "PE",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_svg = FALSE,
    export_png = FALSE,
    export_pdf = FALSE,
    export_csv = FALSE,
    output_file = "max_vintage_cutoff_combined_alpha_mkt_plot",
    width = 14,
    height = 7,
    png_dpi = 300,
    y.max.mkt = NULL,
    y.min.mkt = NULL,
    y.max.alpha = NULL,
    y.min.alpha = NULL,
    x.max = NULL,
    x.min = NULL,
    v.lines.mkt = NULL,
    h.lines.mkt = NULL,
    v.colors.mkt = c("black", "black"),
    h.colors.mkt = c("black", "black"),
    v.lines.alpha = NULL,
    h.lines.alpha = NULL,
    v.colors.alpha = c("black", "black"),
    h.colors.alpha = c("black", "black"),
    main.linewidth = 0.7,
    abline.linewidth = 0.7,
    cex = 1.0,
    weighting_filter = c("EW", "FW"),
    region = NULL) {

    # -------------------------------------------------------------------------
    # 1. Load data
    # -------------------------------------------------------------------------
    region_tag <- if (!is.null(region) && region != "") paste0("_", region) else ""
    all_rows <- list()
    for (v in vintages) {
        for (i in seq_along(nc_tags)) {
            nt <- nc_tags[i]
            nl <- nc_labels[i]
            for (w in c("EW", "FW")) {
                folder <- file.path(data_dir, paste0(
                    "cache_q_factors_", data_source, "_", w,
                    "_VYP", region_tag, nt, "_max_vin_", v
                ))
                csv_path <- file.path(folder, "0_asymptotic_inference_summary.csv")
                if (!file.exists(csv_path)) {
                    warning(paste("File not found (skipped):", csv_path))
                    next
                }
                df <- read.result.summary(csv_path)
                if (nrow(df) == 0) next
                df$vintage <- v
                df$weighting <- w
                df$nc_tag <- nt
                df$nc_label <- nl
                df$method <- "Asymptotic"
                all_rows[[length(all_rows) + 1L]] <- df
            }
        }
    }
    all_data <- bind_rows(all_rows)
    if (nrow(all_data) == 0) stop("No data could be loaded from the specified directory.")

    if (!is.null(weighting_filter)) {
        all_data <- all_data %>% filter(weighting %in% weighting_filter)
    }

    plot_data <- all_data %>% filter(Type == fund_type, Factor == "Alpha")
    if (nrow(plot_data) == 0) stop("No Alpha-factor data found for fund type.")

    plot_data$max.month <- as.numeric(plot_data$max.month)
    plot_data$horizon_years <- plot_data$max.month / 12
    plot_data$vintage <- factor(plot_data$vintage)

    if (!is.null(x.min)) plot_data <- plot_data %>% filter(max.month >= x.min)
    if (!is.null(x.max)) plot_data <- plot_data %>% filter(max.month <= x.max)

    # -------------------------------------------------------------------------
    # 2. Publication theme
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
            legend.title = element_text(face = "bold", size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 6 else 9) * cex),
            legend.text = element_text(size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 5 else 8) * cex),
            legend.key.size = unit((if (!is.null(weighting_filter) && length(weighting_filter) == 1) 0.5 else 0.8) * cex, "cm"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.margin = margin(t = 5),
            plot.title = element_text(face = "bold", size = 10 * cex, hjust = 0.5, margin = margin(b = 5)),
            plot.margin = margin(5, 5, 5, 5)
        )

    # -------------------------------------------------------------------------
    # 3. Panel builder
    # -------------------------------------------------------------------------
    .build_panel <- function(df, y_col, title, show_x, show_y_label, y_label,
                             ylim_vec = NULL,
                             v.lines_vec = NULL, v.colors_vec = NULL,
                             h.lines_vec = NULL, h.colors_vec = NULL) {
        p <- ggplot(df, aes(x = horizon_years, y = .data[[y_col]], color = vintage, group = vintage)) +
            geom_line(linewidth = main.linewidth, alpha = 0.85) +
            geom_point(size = 1.5 * cex, alpha = 0.85) +
            scale_color_viridis_d(option = "turbo", name = "Max Vintage Year", guide = guide_legend(nrow = 1)) +
            scale_x_continuous(breaks = seq(0, max(c(0, df$horizon_years), na.rm = TRUE), by = 5), expand = c(0.02, 0)) +
            labs(title = title, x = if (show_x) "Horizon (Years)" else NULL, y = if (show_y_label) y_label else NULL) +
            {
                if (!is.null(v.lines_vec)) {
                    lapply(seq_along(v.lines_vec), function(idx) {
                        geom_vline(xintercept = v.lines_vec[idx], color = v.colors_vec[min(idx, length(v.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines_vec)) {
                    lapply(seq_along(h.lines_vec), function(idx) {
                        geom_hline(yintercept = h.lines_vec[idx], color = h.colors_vec[min(idx, length(h.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
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

    # -------------------------------------------------------------------------
    # 4. Build MKT row (Panel A — top) and Alpha row (Panel B — bottom)
    # -------------------------------------------------------------------------
    mkt_plots <- list()
    alpha_plots <- list()

    ylim_mkt <- NULL
    if (!is.null(y.max.mkt) || !is.null(y.min.mkt)) {
        ylim_mkt <- c(if (is.null(y.min.mkt)) NA else y.min.mkt, if (is.null(y.max.mkt)) NA else y.max.mkt)
    }
    ylim_alpha <- NULL
    if (!is.null(y.max.alpha) || !is.null(y.min.alpha)) {
        ylim_alpha <- c(if (is.null(y.min.alpha)) NA else y.min.alpha, if (is.null(y.max.alpha)) NA else y.max.alpha)
    }

    for (i in seq_along(nc_tags)) {
        nt <- nc_tags[i]
        nl <- nc_labels[i]
        for (w in weighting_filter) {
            sub_df <- plot_data %>% filter(weighting == w, nc_tag == nt)
            title_text <- paste0(fund_type, ": ", nl, " (", w, ")")

            # -- MKT panel (top) --
            show_y_mkt <- (i == 1 && w == weighting_filter[1])
            p_mkt <- .build_panel(
                df = sub_df, y_col = "MKT", title = title_text,
                show_x = FALSE, show_y_label = show_y_mkt,
                y_label = expression(beta[MKT]), ylim_vec = ylim_mkt,
                v.lines_vec = v.lines.mkt, v.colors_vec = v.colors.mkt,
                h.lines_vec = h.lines.mkt, h.colors_vec = h.colors.mkt
            )
            mkt_plots[[length(mkt_plots) + 1L]] <- p_mkt

            # -- Alpha panel (bottom) --
            show_y_alpha <- (i == 1 && w == weighting_filter[1])
            p_alpha <- .build_panel(
                df = sub_df, y_col = "Coef", title = NULL,
                show_x = TRUE, show_y_label = show_y_alpha,
                y_label = expression(alpha), ylim_vec = ylim_alpha,
                v.lines_vec = v.lines.alpha, v.colors_vec = v.colors.alpha,
                h.lines_vec = h.lines.alpha, h.colors_vec = h.colors.alpha
            )
            alpha_plots[[length(alpha_plots) + 1L]] <- p_alpha
        }
    }

    # -------------------------------------------------------------------------
    # 5. Assemble with patchwork (Panel A label / MKT row / Panel B label / Alpha row)
    # -------------------------------------------------------------------------
    panel_a_label <- wrap_elements(
        panel = textGrob(
            expression(bold("Panel A: Market Factor (" * beta[MKT] * ")")),
            x = 0, y = 0, hjust = 0, vjust = 0,
            gp = gpar(fontsize = 12 * cex, fontfamily = "serif", fontface = "bold")
        )
    )
    panel_b_label <- wrap_elements(
        panel = textGrob(
            expression(bold("Panel B: Alpha (" * alpha * ")")),
            x = 0, y = 0, hjust = 0, vjust = 0,
            gp = gpar(fontsize = 12 * cex, fontfamily = "serif", fontface = "bold")
        )
    )

    row1 <- Reduce(`|`, mkt_plots)
    row2 <- Reduce(`|`, alpha_plots)

    combined_plot <- panel_a_label / row1 / panel_b_label / row2 +
        plot_layout(heights = c(0.08, 1, 0.08, 1), guides = "collect") &
        theme(legend.position = "bottom")

    # -------------------------------------------------------------------------
    # 6. Export or display
    # -------------------------------------------------------------------------
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
            csv_export <- all_data %>% filter(Type == fund_type, Factor == "Alpha")
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
# Plot NAV Discounts Combined Market Analysis
# ============================================================================

#' Plot Empirical Estimates Across Different NAV Discounts
#'
#' Creates a multi-panel chart comparing estimated coefficients across different
#' NAV discount specifications. For each weighting filter, a panel shows lines—
#' one per NAV discount—over horizons.
plot_nav_discounts_combined_mkt <- function(
    data_dir,
    data_source = "preqin",
    fund_type = "PE",
    vintages = 2010,
    nc_tags = c("_NC0", "_NC25", "_NC50", "_NC75"),
    nc_labels = c("0% NAV as Cashflow", "25% NAV-discount sample", "50% NAV-discount sample", "75% NAV-discount sample"),
    export_svg = FALSE,
    export_png = FALSE,
    export_pdf = FALSE,
    export_csv = FALSE,
    output_file = "nav_discount_combined_plot",
    width = 14,
    height = 4,
    png_dpi = 300,
    y.max.mkt = NULL,
    y.min.mkt = NULL,
    x.max = NULL,
    x.min = NULL,
    v.lines.mkt = NULL,
    h.lines.mkt = NULL,
    v.colors.mkt = c("black", "black"),
    h.colors.mkt = c("black", "black"),
    main.linewidth = 0.7,
    abline.linewidth = 0.7,
    cex = 1.0,
    weighting_filter = c("EW", "FW"),
    region = NULL) {
    region_tag <- if (!is.null(region) && region != "") paste0("_", region) else ""
    all_rows <- list()
    for (i in seq_along(nc_tags)) {
        nt <- nc_tags[i]
        nl <- nc_labels[i]
        for (w in c("EW", "FW")) {
            folder <- file.path(data_dir, paste0("cache_q_factors_", data_source, "_", w, "_VYP", region_tag, nt))
            if (!is.null(vintages)) {
                # If vintages is specified, it just adds the parameter to the plot object, but doesn't change path
                v <- vintages[1]
            } else {
                v <- NA
            }
            csv_path <- file.path(folder, "0_asymptotic_inference_summary.csv")
            if (!file.exists(csv_path)) {
                warning(paste("File not found (skipped):", csv_path))
                next
            }
            df <- read.result.summary(csv_path)
            if (nrow(df) == 0) next
            df$vintage <- v
            df$weighting <- w
            df$nc_tag <- nt
            df$nc_label <- nl
            df$method <- "Asymptotic"
            all_rows[[length(all_rows) + 1L]] <- df
        }
    }
    all_data <- bind_rows(all_rows)
    if (nrow(all_data) == 0) stop("No data could be loaded from the specified directory.")

    if (!is.null(weighting_filter)) {
        all_data <- all_data %>% filter(weighting %in% weighting_filter)
    }

    plot_data <- all_data %>% filter(Type == fund_type)
    if (nrow(plot_data) == 0) stop("No data found for fund type.")

    plot_data$max.month <- as.numeric(plot_data$max.month)
    plot_data$horizon_years <- plot_data$max.month / 12
    # Ensure factor order follows the input vector for nc_labels
    plot_data$nc_label <- factor(plot_data$nc_label, levels = unique(nc_labels))

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
            legend.title = element_text(face = "bold", size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 6 else 7) * cex),
            legend.text = element_text(size = (if (!is.null(weighting_filter) && length(weighting_filter) == 1) 5 else 6) * cex),
            legend.key.size = unit((if (!is.null(weighting_filter) && length(weighting_filter) == 1) 0.5 else 0.5) * cex, "cm"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.margin = margin(t = 5),
            plot.title = element_text(face = "bold", size = 10 * cex, hjust = 0.5, margin = margin(b = 5)),
            plot.margin = margin(5, 5, 5, 5)
        )

    .build_panel <- function(df, y_col, title, show_y_label, y_label, ylim_vec = NULL,
                             v.lines_vec = NULL, v.colors_vec = NULL,
                             h.lines_vec = NULL, h.colors_vec = NULL) {
        p <- ggplot(df, aes(x = horizon_years, y = .data[[y_col]], color = nc_label, group = nc_label)) +
            geom_line(linewidth = main.linewidth, alpha = 0.85) +
            geom_point(size = 1.5 * cex, alpha = 0.85) +
            scale_color_viridis_d(option = "plasma", name = "NAV Treatment", guide = guide_legend(nrow = 1)) +
            scale_x_continuous(breaks = seq(0, max(c(0, df$horizon_years), na.rm = TRUE), by = 5), expand = c(0.02, 0)) +
            labs(title = title, x = "Horizon (Years)", y = if (show_y_label) y_label else NULL) +
            {
                if (!is.null(v.lines_vec)) {
                    lapply(seq_along(v.lines_vec), function(idx) {
                        geom_vline(xintercept = v.lines_vec[idx], color = v.colors_vec[min(idx, length(v.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
                    })
                } else {
                    NULL
                }
            } +
            {
                if (!is.null(h.lines_vec)) {
                    lapply(seq_along(h.lines_vec), function(idx) {
                        geom_hline(yintercept = h.lines_vec[idx], color = h.colors_vec[min(idx, length(h.colors_vec))], linetype = "dotted", linewidth = abline.linewidth)
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
    for (i in seq_along(weighting_filter)) {
        w <- weighting_filter[i]

        mkt_df <- plot_data %>%
            filter(Factor == "MKT", weighting == w)

        title_text <- paste0(fund_type, ": NAV Discounts (", w, ")")
        ylim_mkt <- NULL
        if (!is.null(y.max.mkt) || !is.null(y.min.mkt)) {
            ylim_mkt <- c(if (is.null(y.min.mkt)) NA else y.min.mkt, if (is.null(y.max.mkt)) NA else y.max.mkt)
        }
        show_y <- (i == 1)

        p <- .build_panel(
            df = mkt_df, y_col = "MKT", title = title_text,
            show_y_label = show_y, y_label = expression(beta[MKT]), ylim_vec = ylim_mkt,
            v.lines_vec = v.lines.mkt, v.colors_vec = v.colors.mkt,
            h.lines_vec = h.lines.mkt, h.colors_vec = h.colors.mkt
        )
        mkt_plots[[length(mkt_plots) + 1L]] <- p
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

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
getwd()

file.folder <- "results/data_out_2026-emp-max-vin-2019"
file.folder <- "data_out_2026_02_26"
out.folder <- "figures"

# # Two-factor models with all factors
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("Alpha", "EG", "IA", "ME", "ROE"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    y.max.mkt = 1.5, y.min.mkt = -0.25,
    y.lim.second = list(Alpha = c(-0.01, 0.01)),
    main.linewidth = 1,
    output_file = paste0(out.folder, "/empirical_PE_all_factors")
)

# # Two-factor model with only Alpha
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "PE",
    factors = c("Alpha"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    y.max.mkt = 1.5, y.min.mkt = -0.25,
    y.lim.second = list(Alpha = c(-0.005, 0.01)),
    main.linewidth = 1.5,
    output_file = paste0(out.folder, "/empirical_PE_Alpha")
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
    main.linewidth = 1.5,
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
    main.linewidth = 0.8,
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
        main.linewidth = 0.8,
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
        main.linewidth = 0.8,
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
        main.linewidth = 0.8,
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
    main.linewidth = 0.8,
    output_file = paste0(out.folder, "/max_vintage_PE_NC50_all_factors")
)


### Paper: Buyout (BO) and Venture Capital (VC)

file.folder <- "data_out_2026_02_26"
out.folder <- "figures_bovc_paper"

# BO: Two-factor model with only Alpha
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("Alpha"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 5,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.647399009343063),
    # h.lines.second = c(1.0194438416743136^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"),
    output_file = paste0(out.folder, "/empirical_BO_Alpha")
)

# BO: Two-factor model with q^5 factors
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("EG", "IA", "ME", "ROE"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 14,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.647399009343063),
    # h.lines.second = c(1.0194438416743136^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"),
    output_file = paste0(out.folder, "/empirical_BO_Multi")
)

# BO: Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 4,
    width = 6,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.769855673660812),
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"),
    # exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/empirical_BO_MKT")
)

# BO Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_mkt(
    data_dir = file.folder,
    fund_type = "BO",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 4,
    y.max.mkt = 2.5, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    output_file = paste0(out.folder, "/max_vintage_BO_MKT_combined_NC")
)

# VC: Two-factor model with only Alpha
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("Alpha"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 5,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.97551979450791),
    # h.lines.second = c((1 - 0.00390550172544275)^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    # exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    output_file = paste0(out.folder, "/empirical_VC_Alpha")
)

# BO: Two-factor model with q^5 factors
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("EG", "IA", "ME", "ROE"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 14,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.647399009343063),
    # h.lines.second = c(1.0194438416743136^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"),
    output_file = paste0(out.folder, "/empirical_VC_Multi")
)

# VC: Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 4,
    width = 6,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.925920685301403),
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"),
    # exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/empirical_VC_MKT")
)

# VC Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_mkt(
    data_dir = file.folder,
    fund_type = "VC",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 4,
    y.max.mkt = 2.5, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    output_file = paste0(out.folder, "/max_vintage_VC_MKT_combined_NC")
)



### Presentation: Buyout (BO) and Venture Capital (VC)

file.folder <- "data_out_2026_02_26"
out.folder <- "figures_bovc"

# BO: Two-factor model with only Alpha
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("Alpha"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    width = 7,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1.5,
    cex = 1,
    h.lines.mkt = c(0.647399009343063),
    h.lines.second = c(1.0194438416743136^(1 / 3) - 1), # alpha quarterly
    abline.linewidth = 1.5,
    h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    output_file = paste0(out.folder, "/empirical_BO_Alpha")
)

# BO: Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 4, width = 7,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    main.linewidth = 1.5,
    cex = 1,
    h.lines.mkt = c(0.769855673660812),
    abline.linewidth = 1.5,
    h.colors.mkt = c("orange"),
    exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/empirical_BO_MKT")
)

# BO Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_mkt(
    data_dir = file.folder,
    fund_type = "BO",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 4,
    y.max.mkt = 2.5, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    output_file = paste0(out.folder, "/max_vintage_BO_MKT_combined_NC")
)

# VC: Two-factor model with only Alpha
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("Alpha"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    width = 7,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1.5,
    cex = 1,
    h.lines.mkt = c(0.97551979450791),
    h.lines.second = c((1 - 0.00390550172544275)^(1 / 3) - 1), # alpha quarterly
    abline.linewidth = 1.5,
    h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    output_file = paste0(out.folder, "/empirical_VC_Alpha")
)

# VC: Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 4, width = 7,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    main.linewidth = 1.5,
    cex = 1,
    h.lines.mkt = c(0.925920685301403),
    abline.linewidth = 1.5,
    h.colors.mkt = c("orange"),
    exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = paste0(out.folder, "/empirical_VC_MKT")
)

# VC Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_mkt(
    data_dir = file.folder,
    fund_type = "VC",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 4,
    y.max.mkt = 2.5, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    output_file = paste0(out.folder, "/max_vintage_VC_MKT_combined_NC")
)


### PRESENTATION

out.folder <- "figures_presentation"


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
    v.lines.mkt = c(0, 15),
    h.lines.mkt = c(0, 0.8),
    v.colors.mkt = c("blue", "orange"),
    h.colors.mkt = c("gray", "red"),
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

### Pitchbook

file.folder <- "data_out_2026_02_27_pitchbook"
out.folder <- "figures_pb"

# Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_mkt(
    data_dir = file.folder,
    data_source = "pitchbook",
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


### Preqin - North America
### Paper: Buyout (BO) and Venture Capital (VC)

file.folder <- "data_out_2026_03_27"
out.folder <- "figures_bovc_us" # North America Results

# BO: Two-factor model with only Alpha
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("Alpha"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 5,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.647399009343063),
    # h.lines.second = c(1.0194438416743136^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"),
    region = "North America",
    output_file = paste0(out.folder, "/empirical_BO_US_Alpha")
)

# BO: Two-factor model with q^5 factors
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("EG", "IA", "ME", "ROE"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 14,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.647399009343063),
    # h.lines.second = c(1.0194438416743136^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"),
    region = "North America",
    output_file = paste0(out.folder, "/empirical_BO_US_Multi")
)

# BO: Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "BO",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 3.5,
    width = 6,
    y.max.mkt = 2.2, y.min.mkt = 0,
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.769855673660812),
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"),
    # exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    region = "North America",
    output_file = paste0(out.folder, "/empirical_BO_US_MKT")
)

# BO Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_mkt(
    data_dir = file.folder,
    fund_type = "BO",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 3,
    y.max.mkt = 3, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    region = "North America",
    output_file = paste0(out.folder, "/max_vintage_BO_US_MKT_combined_NC")
)

# BO Combined Alpha+MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_alpha_mkt(
    data_dir = file.folder,
    fund_type = "BO",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 6,
    y.max.mkt = 3, y.min.mkt = 0,
    y.max.alpha = 0.02, y.min.alpha = -0.02,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    region = "North America",
    output_file = paste0(out.folder, "/max_vintage_BO_US_Alpha_MKT_combined_NC")
)

# VC: Two-factor model with only Alpha
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("Alpha"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 5,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.97551979450791),
    # h.lines.second = c((1 - 0.00390550172544275)^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    # exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    region = "North America",
    output_file = paste0(out.folder, "/empirical_VC_US_Alpha")
)

# VC: Two-factor model with q^5 factors
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("EG", "IA", "ME", "ROE"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 5,
    width = 14,
    y.max.mkt = 2.5, y.min.mkt = -0.5,
    y.lim.second = list(Alpha = c(-0.02, 0.02)),
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.647399009343063),
    # h.lines.second = c(1.0194438416743136^(1 / 3) - 1), # alpha quarterly
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"), h.colors.second = c("orange"),
    exclude_cv = FALSE,
    weighting_filter = c("FW", "EW"),
    region = "North America",
    output_file = paste0(out.folder, "/empirical_VC_US_Multi")
)

# VC: Single-factor model with MKT only
plot_empirical_estimates(
    data_dir = file.folder,
    fund_type = "VC",
    factors = c("MKT"),
    export_pdf = TRUE,
    export_csv = TRUE,
    export_latex = TRUE,
    height = 3.5,
    width = 6,
    y.max.mkt = 2.2, y.min.mkt = 0,
    main.linewidth = 1,
    cex = 0.8,
    # h.lines.mkt = c(0.925920685301403),
    # abline.linewidth = 1.5,
    # h.colors.mkt = c("orange"),
    # exclude_cv = TRUE,
    weighting_filter = c("FW", "EW"),
    # y.lim.second = list(Alpha = c(-0.01, 0.01)),
    region = "North America",
    output_file = paste0(out.folder, "/empirical_VC_US_MKT")
)

# VC Combined MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_mkt(
    data_dir = file.folder,
    fund_type = "VC",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 3,
    y.max.mkt = 3, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    region = "North America",
    output_file = paste0(out.folder, "/max_vintage_VC_US_MKT_combined_NC")
)

# VC Combined Alpha+MKT Max vintage-year cutoff analysis (100% NAV and 50% NAV-discount)
plot_max_vintage_cutoff_combined_alpha_mkt(
    data_dir = file.folder,
    fund_type = "VC",
    vintages = 2011:2021,
    nc_tags = c("", "_NC50"),
    nc_labels = c("100% NAV as Cashflow", "50% NAV-discount sample"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 7,
    height = 6,
    y.max.mkt = 3, y.min.mkt = 0,
    y.max.alpha = 0.02, y.min.alpha = -0.02,
    main.linewidth = 0.8,
    weighting_filter = c("FW"),
    region = "North America",
    output_file = paste0(out.folder, "/max_vintage_VC_US_Alpha_MKT_combined_NC")
)

# Analyze Different NAV Discounts

file.folder <- "data_out_2026_03_27"
out.folder <- "figures_bovc_us" # North America Results

plot_nav_discounts_combined_mkt(
    data_dir = file.folder,
    fund_type = "VC",
    vintages = 2010,
    nc_tags = c("", "_NC75", "_NC50", "_NC25", "_NC0"),
    nc_labels = c("NAV as Cashflow", "25% NAV Disc.", "50% NAV Disc.", "75% NAV Disc.", "Full Write-Off"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 6,
    height = 3,
    y.max.mkt = 2.1, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW", "EW"),
    region = "North America",
    output_file = paste0(out.folder, "/nav_disc_VC_US_MKT_combined_NC")
)

plot_nav_discounts_combined_mkt(
    data_dir = file.folder,
    fund_type = "BO",
    vintages = 2010,
    nc_tags = c("", "_NC75", "_NC50", "_NC25", "_NC0"),
    nc_labels = c("NAV as Cashflow", "25% NAV Disc.", "50% NAV Disc.", "75% NAV Disc.", "Full Write-Off"),
    export_pdf = TRUE,
    export_svg = TRUE,
    export_csv = TRUE,
    width = 6,
    height = 3,
    y.max.mkt = 2.1, y.min.mkt = 0,
    main.linewidth = 0.8,
    weighting_filter = c("FW", "EW"),
    region = "North America",
    output_file = paste0(out.folder, "/nav_disc_BO_US_MKT_combined_NC")
)

# L1 Loss Function 
# BO: Single-factor model with MKT only (L1)
plot_empirical_estimates(
  data_dir = file.folder,
  fund_type = "BO",
  factors = c("MKT"),
  export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
  height = 3, width = 6,
  y.max.mkt = 2.5, y.min.mkt = -0.5,
  main.linewidth = 1, cex = 0.8,
  weighting_filter = c("FW", "EW"),
  region = "North America",
  error_function = "L1_Ridge",
  output_file = paste0(out.folder, "/empirical_BO_US_MKT_L1")
)

# BO: Two-factor model with Alpha (L1)
plot_empirical_estimates(
  data_dir = file.folder,
  fund_type = "BO",
  factors = c("Alpha"),
  export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
  height = 4, width = 5,
  y.max.mkt = 2.5, y.min.mkt = -0.5,
  y.lim.second = list(Alpha = c(-0.02, 0.02)),
  main.linewidth = 1, cex = 0.8,
  exclude_cv = FALSE,
  weighting_filter = c("FW", "EW"),
  region = "North America",
  error_function = "L1_Ridge",
  output_file = paste0(out.folder, "/empirical_BO_US_Alpha_L1")
)

# VC: Single-factor model with MKT only (L1)
plot_empirical_estimates(
  data_dir = file.folder,
  fund_type = "VC",
  factors = c("MKT"),
  export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
  height = 3, width = 6,
  y.max.mkt = 2.5, y.min.mkt = -0.5,
  main.linewidth = 1, cex = 0.8,
  weighting_filter = c("FW", "EW"),
  region = "North America",
  error_function = "L1_Ridge",
  output_file = paste0(out.folder, "/empirical_VC_US_MKT_L1")
)

# VC: Two-factor model with Alpha (L1)
plot_empirical_estimates(
  data_dir = file.folder,
  fund_type = "VC",
  factors = c("Alpha"),
  export_pdf = TRUE, export_csv = TRUE, export_latex = TRUE,
  height = 4, width = 5,
  y.max.mkt = 2.5, y.min.mkt = -0.5,
  y.lim.second = list(Alpha = c(-0.02, 0.02)),
  main.linewidth = 1, cex = 0.8,
  weighting_filter = c("FW", "EW"),
  region = "North America",
  error_function = "L1_Ridge",
  output_file = paste0(out.folder, "/empirical_VC_US_Alpha_L1")
)

# Cross-fund pair-product (U_{v,H}) tables, both vintage weightings (Gate G2)
export_u_pair_tables(data_dir = file.folder, fund_type = "BO", region = "North America",
  output_file = paste0(out.folder, "/upair_BO_US"))
export_u_pair_tables(data_dir = file.folder, fund_type = "VC", region = "North America",
  output_file = paste0(out.folder, "/upair_VC_US"))
