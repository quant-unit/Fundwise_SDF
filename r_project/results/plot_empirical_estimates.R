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
    output_file = "empirical_estimates_plot",
    width = 14,
    height = 6,
    png_dpi = 300,
    y.max.mkt = NULL,
    y.min.mkt = NULL,
    y.lim.second = NULL,
    x.max = NULL,
    x.min = NULL) {
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
    all_data <- lapply(csv_sources, function(src) {
        if (!file.exists(src$path)) {
            warning(paste("File not found:", src$path))
            return(NULL)
        }
        df <- read.csv(src$path, stringsAsFactors = FALSE)
        # Select only the columns needed for plotting (CSVs have different extras)
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
    theme_publication <- theme_minimal(base_size = 11, base_family = "serif") +
        theme(
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
            panel.background = element_rect(fill = "white"),
            strip.text = element_text(face = "bold", size = 10, margin = margin(b = 5, t = 5)),
            strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.5),
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
                face = "bold", size = 13, hjust = 0.5,
                margin = margin(b = 10)
            ),
            plot.subtitle = element_text(
                size = 10, hjust = 0.5, color = "grey30",
                margin = margin(b = 15)
            ),
            plot.margin = margin(15, 15, 10, 15)
        )

    # -------------------------------------------------------------------------
    # Color and linetype palette
    # -------------------------------------------------------------------------
    # EW = blue tones, FW = orange tones
    # Asymptotic = solid, Cross-Validation = dashed
    colors_source <- c(
        "EW Asymptotic"        = "#0072B2",
        "EW Cross-Validation"  = "#56B4E9",
        "FW Asymptotic"        = "#D55E00",
        "FW Cross-Validation"  = "#E69F00"
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
            geom_line(linewidth = 0.8, alpha = 0.9) +
            geom_point(size = 2.5, fill = "white", stroke = 0.8) +
            scale_color_manual(values = colors_source, name = "Source") +
            scale_linetype_manual(values = linetypes_source, name = "Source") +
            scale_shape_manual(values = shapes_source, name = "Source") +
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
            theme_publication +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 5)),
                axis.title.y = if (i > 1) element_blank() else element_text(face = "bold", size = 10)
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
                    theme(plot.margin = margin(15, 15, 10, 15))
                second_plots[[i]] <- blank_plot
            } else {
                coef_df <- plot_data %>%
                    filter(Factor == cf) %>%
                    select(horizon_years, Coef, source)

                p <- ggplot(coef_df, aes(
                    x = horizon_years, y = Coef,
                    color = source, linetype = source, shape = source
                )) +
                    geom_line(linewidth = 0.8, alpha = 0.9) +
                    geom_point(size = 2.5, fill = "white", stroke = 0.8) +
                    scale_color_manual(values = colors_source, name = "Source") +
                    scale_linetype_manual(values = linetypes_source, name = "Source") +
                    scale_shape_manual(values = shapes_source, name = "Source") +
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
                    theme_publication +
                    theme(
                        legend.position = "none",
                        plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 5)),
                        axis.title.y = if (!first_coef_plot) element_blank() else element_text(face = "bold", size = 10)
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

    # Create Panel A label as a dedicated text element
    panel_a_label <- wrap_elements(
        panel = textGrob(
            expression(bold("Panel A: Market Factor (" * beta[MKT] * ")")),
            x = 0, y = 0, hjust = 0, vjust = 0,
            gp = gpar(fontsize = 12, fontfamily = "serif", fontface = "bold")
        )
    )

    if (length(second_plots) > 0) {
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
                gp = gpar(fontsize = 12, fontfamily = "serif", fontface = "bold")
            )
        )

        # Stack: Panel A label, row1, Panel B label, row2
        combined_plot <- panel_a_label / row1 / panel_b_label / row2 +
            plot_layout(heights = c(0.1, 1, 0.1, 1))
    } else {
        # Single-factor only
        combined_plot <- panel_a_label / row1 +
            plot_layout(heights = c(0.1, 1))
    }

    # -------------------------------------------------------------------------
    # Add shared legend at bottom
    # -------------------------------------------------------------------------
    # Collect legends and unify at bottom via patchwork
    combined_plot <- combined_plot +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom") &
        scale_color_manual(values = colors_source, name = "Estimation Method") &
        scale_linetype_manual(values = linetypes_source, name = "Estimation Method") &
        scale_shape_manual(values = shapes_source, name = "Estimation Method")

    # -------------------------------------------------------------------------
    # Export or display
    # -------------------------------------------------------------------------
    if (export_svg || export_png || export_pdf) {
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

        invisible(combined_plot)
    } else {
        print(combined_plot)
        invisible(combined_plot)
    }
}


# ============================================================================
# Example Usage (uncomment to run)
# ============================================================================

# # Two-factor model with all factors
plot_empirical_estimates(
    data_dir = "results/data_out_2026-emp-max-vin-2019",
    fund_type = "PE",
    factors = c("Alpha", "EG", "IA", "ME", "ROE"),
    export_pdf = TRUE,
    y.max.mkt = 2.5, y.min.mkt = 1,
    y.lim.second = list(Alpha = c(-0.01, 0.01)),
    output_file = "results/figures/empirical_PE_all_factors"
)
