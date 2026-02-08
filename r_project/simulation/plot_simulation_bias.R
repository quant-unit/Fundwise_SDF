# Plot simulation bias
# ============================================================================
# Visualize simulation bias results across scenarios and horizons
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Plot Simulation Bias Across Scenarios and Horizons
#'
#' Creates a publication-quality multi-panel chart showing estimation bias
#' for MKT and second factor across different time horizons.
#'
#' @param bias_file Path to the bias_by_scenario_horizon.csv file
#' @param scenarios Character vector of scenario IDs to plot (4-5 recommended)
#' @param export_svg Logical, whether to export as SVG file
#' @param output_file Path for SVG output (only used if export_svg = TRUE)
#' @param width Plot width in inches (default: 14)
#' @param height Plot height in inches (default: 8)
#' @param show_relative_bias Logical, whether to show relative bias (%) instead of absolute
#' @param y.max.mkt Optional numeric, maximum y-axis value for MKT panel (default: NULL, auto-scale)
#' @param y.max.second Optional numeric, maximum y-axis value for second factor panel (default: NULL, auto-scale)
#'
#' @return A ggplot object (invisibly if exported)
#'
#' @examples
#' plot_simulation_bias(
#'     bias_file = "data_out_2026_new/bias_analysis/2026-02-07_194037_bias_by_scenario_horizon.csv",
#'     scenarios = c("base_case_vyp", "base_case_cross_sectional", "base_case_ME")
#' )
#'
plot_simulation_bias <- function(
    bias_file,
    scenarios = c(
        "base_case_vyp", "base_case_cross_sectional",
        "base_case_ME", "base_case_ROE", "base_case_IA"
    ),
    export_svg = FALSE,
    export_png = FALSE,
    export_pdf = FALSE,
    output_file = "simulation_bias_plot",
    width = 14,
    height = 8,
    png_dpi = 300,
    show_relative_bias = FALSE,
    y.max.mkt = NULL,
    y.max.second = NULL) {
    # -------------------------------------------------------------------------
    # Read and prepare data
    # -------------------------------------------------------------------------
    bias_data <- read.csv(bias_file, stringsAsFactors = FALSE)

    # Filter to selected scenarios
    plot_data <- bias_data %>%
        filter(scenario_id %in% scenarios)

    # Check which scenarios were found
    found_scenarios <- unique(plot_data$scenario_id)
    missing_scenarios <- setdiff(scenarios, found_scenarios)

    if (length(missing_scenarios) > 0) {
        warning(paste("Scenarios not found in data:", paste(missing_scenarios, collapse = ", ")))
    }

    if (nrow(plot_data) == 0) {
        stop("No data found for the specified scenarios.")
    }

    # Order scenarios as specified
    plot_data$scenario_id <- factor(plot_data$scenario_id, levels = scenarios)

    # Create cleaner scenario labels
    scenario_labels <- gsub("_", " ", scenarios)
    scenario_labels <- tools::toTitleCase(scenario_labels)
    names(scenario_labels) <- scenarios

    # -------------------------------------------------------------------------
    # Define publication-quality theme
    # -------------------------------------------------------------------------
    theme_publication <- theme_minimal(base_size = 11, base_family = "serif") +
        theme(
            # Panel appearance
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
            panel.background = element_rect(fill = "white"),

            # Strip (facet labels)
            strip.text = element_text(face = "bold", size = 10, margin = margin(b = 5, t = 5)),
            strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.5),

            # Axis appearance
            axis.title = element_text(face = "bold", size = 10),
            axis.text = element_text(size = 9, color = "grey20"),
            axis.ticks = element_line(color = "grey40", linewidth = 0.3),
            axis.line = element_blank(),

            # Legend
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9),
            legend.key.size = unit(0.8, "cm"),
            legend.background = element_rect(fill = "white", color = NA),
            legend.margin = margin(t = 5),

            # Overall plot
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
    # Color palette (colorblind-friendly)
    # -------------------------------------------------------------------------
    # Using a professional palette suitable for journals
    colors_est <- "#0072B2" # Blue - estimated values
    colors_true <- "#D55E00" # Orange/Red - true values

    # -------------------------------------------------------------------------
    # Create MKT bias plot (top row)
    # -------------------------------------------------------------------------

    mkt_data <- plot_data %>%
        select(scenario_id, max_month, est_beta_MKT, true_beta_MKT, first_factor_name) %>%
        pivot_longer(
            cols = c(est_beta_MKT, true_beta_MKT),
            names_to = "type",
            values_to = "beta_MKT"
        ) %>%
        mutate(
            type = ifelse(type == "est_beta_MKT", "Estimated", "True (DGP)")
        )

    # Convert horizon to years for better readability
    mkt_data$horizon_years <- mkt_data$max_month / 12

    p_mkt <- ggplot(mkt_data, aes(
        x = horizon_years, y = beta_MKT,
        color = type, shape = type
    )) +
        geom_line(linewidth = 0.8, alpha = 0.9) +
        geom_point(size = 2.5, fill = "white", stroke = 0.8) +
        facet_wrap(~scenario_id,
            nrow = 1,
            labeller = labeller(scenario_id = scenario_labels)
        ) +
        scale_color_manual(
            values = c("Estimated" = colors_est, "True (DGP)" = colors_true),
            name = "Value Type"
        ) +
        scale_shape_manual(
            values = c("Estimated" = 16, "True (DGP)" = 17),
            name = "Value Type"
        ) +
        scale_x_continuous(
            breaks = seq(0, max(mkt_data$horizon_years, na.rm = TRUE), by = 5),
            expand = c(0.02, 0)
        ) +
        labs(
            title = expression(bold("Panel A: Market Factor (" * beta[MKT] * ")")),
            x = "Horizon (Years)",
            y = expression(beta[MKT])
        ) +
        {
            if (!is.null(y.max.mkt)) coord_cartesian(ylim = c(NA, y.max.mkt)) else NULL
        } +
        theme_publication +
        theme(
            legend.position = "none", # Will use combined legend at bottom
            plot.title = element_text(size = 11, hjust = 0, margin = margin(b = 8))
        )

    # -------------------------------------------------------------------------
    # Create Second Factor bias plot (bottom row)
    # -------------------------------------------------------------------------

    # Check if any scenarios have a second factor
    has_second_factor <- any(!is.na(plot_data$est_second))

    if (has_second_factor) {
        # Keep ALL scenarios but with NA values for those without second factor
        # This ensures proper alignment with the top row
        second_data <- plot_data %>%
            select(scenario_id, max_month, est_second, true_second, second_factor_name) %>%
            pivot_longer(
                cols = c(est_second, true_second),
                names_to = "type",
                values_to = "beta_second"
            ) %>%
            mutate(
                type = ifelse(type == "est_second", "Estimated", "True (DGP)")
            )

        second_data$horizon_years <- second_data$max_month / 12

        # Create dynamic y-axis label based on factor names
        second_factor_names <- unique(na.omit(plot_data$second_factor_name))
        if (length(second_factor_names) == 1) {
            y_label <- bquote(beta[.(second_factor_names[1])])
        } else {
            y_label <- expression(beta["Second Factor"])
        }

        p_second <- ggplot(second_data, aes(
            x = horizon_years, y = beta_second,
            color = type, shape = type
        )) +
            geom_line(linewidth = 0.8, alpha = 0.9, na.rm = TRUE) +
            geom_point(size = 2.5, fill = "white", stroke = 0.8, na.rm = TRUE) +
            facet_wrap(~scenario_id,
                nrow = 1,
                drop = FALSE, # Keep all factor levels even if empty
                labeller = labeller(scenario_id = scenario_labels)
            ) +
            scale_color_manual(
                values = c("Estimated" = colors_est, "True (DGP)" = colors_true),
                name = "Value Type"
            ) +
            scale_shape_manual(
                values = c("Estimated" = 16, "True (DGP)" = 17),
                name = "Value Type"
            ) +
            scale_x_continuous(
                breaks = seq(0, max(second_data$horizon_years, na.rm = TRUE), by = 5),
                expand = c(0.02, 0)
            ) +
            labs(
                title = "Panel B: Second Factor",
                x = "Horizon (Years)",
                y = y_label
            ) +
            {
                if (!is.null(y.max.second)) coord_cartesian(ylim = c(NA, y.max.second)) else NULL
            } +
            theme_publication +
            theme(
                legend.position = "bottom",
                plot.title = element_text(size = 11, hjust = 0, margin = margin(b = 8))
            )

        # Combine panels using patchwork
        combined_plot <- p_mkt / p_second +
            plot_layout(heights = c(1, 1), guides = "collect") +
            plot_annotation(
                title = "Simulation Bias Analysis: Estimated vs. True Factor Loadings",
                subtitle = "Comparison across scenarios and estimation horizons",
                theme = theme(
                    plot.title = element_text(
                        face = "bold", size = 14, hjust = 0.5,
                        family = "serif"
                    ),
                    plot.subtitle = element_text(
                        size = 11, hjust = 0.5, color = "grey30",
                        family = "serif", margin = margin(b = 5)
                    )
                )
            ) &
            theme(legend.position = "bottom")
    } else {
        # Only MKT panel if no second factors
        combined_plot <- p_mkt +
            labs(title = "Simulation Bias Analysis: Estimated vs. True Market Factor Loading") +
            theme(
                legend.position = "bottom",
                plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
            )
    }

    # -------------------------------------------------------------------------
    # Export or display
    # -------------------------------------------------------------------------

    if (export_svg || export_png || export_pdf) {
        # Strip any existing extension from output_file
        base_file <- sub("\\.(svg|png|pdf)$", "", output_file, ignore.case = TRUE)

        # Ensure directory exists
        output_dir <- dirname(base_file)
        if (!dir.exists(output_dir) && output_dir != "." && output_dir != "") {
            dir.create(output_dir, recursive = TRUE)
        }

        if (export_svg) {
            # Export to SVG (requires svglite package)
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
                    # Fallback to PDF if svglite is not available
                    warning("SVG export failed (svglite not installed?). Falling back to PDF.")
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
            # Export to PNG
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
            # Export to PDF (best for LaTeX, handles fonts properly)
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
        # Display in RStudio
        print(combined_plot)
        invisible(combined_plot)
    }
}


#' Quick preview function for examining bias by scenario
#'
#' @param bias_file Path to the bias CSV file
#' @param scenario_id Single scenario ID to examine in detail
#'
plot_single_scenario_bias <- function(bias_file, scenario_id) {
    bias_data <- read.csv(bias_file, stringsAsFactors = FALSE) %>%
        filter(scenario_id == !!scenario_id)

    if (nrow(bias_data) == 0) {
        stop(paste("Scenario not found:", scenario_id))
    }

    # Create a simple detailed view
    has_second <- any(!is.na(bias_data$est_second))

    # Prepare data
    plot_df <- bias_data %>%
        mutate(horizon_years = max_month / 12) %>%
        select(
            horizon_years, est_beta_MKT, true_beta_MKT,
            est_second, true_second, bias_MKT, bias_second
        )

    # MKT plot
    p1 <- ggplot(plot_df, aes(x = horizon_years)) +
        geom_ribbon(
            aes(
                ymin = pmin(est_beta_MKT, true_beta_MKT),
                ymax = pmax(est_beta_MKT, true_beta_MKT)
            ),
            fill = "grey80", alpha = 0.5
        ) +
        geom_line(aes(y = est_beta_MKT, color = "Estimated"), linewidth = 1) +
        geom_line(aes(y = true_beta_MKT, color = "True"), linewidth = 1, linetype = "dashed") +
        geom_point(aes(y = est_beta_MKT, color = "Estimated"), size = 3) +
        geom_point(aes(y = true_beta_MKT, color = "True"), size = 3, shape = 17) +
        scale_color_manual(values = c("Estimated" = "#0072B2", "True" = "#D55E00")) +
        labs(
            title = paste("Market Factor Bias:", scenario_id),
            subtitle = "Shaded area indicates bias magnitude",
            x = "Horizon (Years)",
            y = expression(beta[MKT]),
            color = "Type"
        ) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")

    if (has_second) {
        p2 <- ggplot(plot_df %>% filter(!is.na(est_second)), aes(x = horizon_years)) +
            geom_ribbon(
                aes(
                    ymin = pmin(est_second, true_second),
                    ymax = pmax(est_second, true_second)
                ),
                fill = "grey80", alpha = 0.5
            ) +
            geom_line(aes(y = est_second, color = "Estimated"), linewidth = 1) +
            geom_line(aes(y = true_second, color = "True"), linewidth = 1, linetype = "dashed") +
            geom_point(aes(y = est_second, color = "Estimated"), size = 3) +
            geom_point(aes(y = true_second, color = "True"), size = 3, shape = 17) +
            scale_color_manual(values = c("Estimated" = "#0072B2", "True" = "#D55E00")) +
            labs(
                title = "Second Factor Bias",
                x = "Horizon (Years)",
                y = expression(beta["Second"]),
                color = "Type"
            ) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "bottom")

        combined <- p1 / p2 + plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
        print(combined)
        invisible(combined)
    } else {
        print(p1)
        invisible(p1)
    }
}

# ============================================================================
# Example Usage (uncomment to run)
# ============================================================================
#
# # Basic usage with default scenarios
plot_simulation_bias(
    scenarios = c("base_case_vyp", "base_case_cross_sectional", "base_case_positive_alpha", "base_case_negative_alpha"),
    bias_file = "simulation/data_out_2026_new/bias_analysis/2026-02-08_184350_bias_by_scenario_horizon.csv",
    export_pdf = TRUE,
    output_file = "figures/bias_comparison1",
    y.max.mkt = 1.25,
    y.max.second = 0.005
)
plot_simulation_bias(
    scenarios = c("big_n_v_40funds", "big_v_10funds_1967", "big_v_20funds_1967", "small_v_1986_1995", "small_v_1996_2005"),
    bias_file = "simulation/data_out_2026_new/bias_analysis/2026-02-08_184350_bias_by_scenario_horizon.csv",
    export_pdf = TRUE,
    height = 4,
    output_file = "figures/bias_comparison2"
)
plot_simulation_bias(
    scenarios = c("base_case_ME", "base_case_IA", "base_case_ROE", "base_case_EG"),
    bias_file = "simulation/data_out_2026_new/bias_analysis/2026-02-08_184350_bias_by_scenario_horizon.csv",
    export_pdf = TRUE,
    output_file = "figures/bias_comparison3"
)
