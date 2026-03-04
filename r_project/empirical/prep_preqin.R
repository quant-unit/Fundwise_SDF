## ----setup, message=FALSE, warning=FALSE--------------------------------------
## prepare preqin
# prep ----
if (sys.nframe() == 0L) rm(list = ls())

library(readxl)
library(here)

data.prepared.folder <- here("empirical/data_prepared_2026_b")


# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
if (!dir.exists(data.prepared.folder)) dir.create(data.prepared.folder)


## ----load_data----------------------------------------------------------------
path <- "empirical/data_in/Preqin_Cashflow_export-26_Feb_20e003d5aa-5a18-4c12-aba8-7586a7435ac9.xlsx"
path <- "empirical/data_in/Preqin_Cashflow_export-14_Apr_22272f54e5-79fe-43c0-8bc0-9206381de0b7.xlsx"
path <- here(path)
sheet <- "Preqin_Export"
df.xl <- data.frame(readxl::read_excel(path = path, sheet = sheet, guess_max = 100000))
if ("GEOGRAPHIC.FOCUS" %in% colnames(df.xl)) df.xl$PRIMARY.GEOGRAPHIC.FOCUS <- df.xl$GEOGRAPHIC.FOCUS
table(df.xl$PRIMARY.GEOGRAPHIC.FOCUS)
table(df.xl$REGION)
colnames(df.xl)

# Exclude too recent vintages
preqin.cutoff.vintage <- 2021
df.xl <- df.xl[df.xl$VINTAGE...INCEPTION.YEAR <= preqin.cutoff.vintage, ]

# df <- df.xl
table(df.xl$ASSET.CLASS[!duplicated(df.xl$FUND.ID)])
table(df.xl$STRATEGY[!duplicated(df.xl$FUND.ID)])
table(df.xl$VINTAGE...INCEPTION.YEAR[!duplicated(df.xl$FUND.ID)])


## ----strategy_summary---------------------------------------------------------
make.strategy.summary <- function() {
    l <- list()
    for (asset.class in levels(as.factor(df.xl$ASSET.CLASS))) {
        # print(asset.class)
        # print(table(df.xl$STRATEGY[df.xl$ASSET.CLASS == asset.class]))
        df <- data.frame(table(df.xl$STRATEGY[(df.xl$ASSET.CLASS == asset.class) & (!duplicated(df.xl$FUND.ID))]))
        colnames(df) <- c("Strategy", "Fund Count")
        df["Asset Class"] <- asset.class
        df <- df[, c(3, 1, 2)]
        l[[asset.class]] <- df
    }
    df <- data.frame(do.call(rbind, l))

    print(
        xtable::xtable(df,
            caption = "Summary of asset classes in Preqin dataset.",
            digits = 3, label = "tab:summary_preqin_strategy"
        ),
        include.rownames = FALSE
        # , format.args = list(big.mark = ",")
    )
    return(df)
}
make.strategy.summary()


## ----vintage_summary----------------------------------------------------------
make.vintage.summary <- function(asset.class = "Private Equity", region = NA) {
    if (asset.class == "Buyout") {
        df <- df.xl[df.xl$STRATEGY == asset.class, ]
    } else {
        df <- df.xl[df.xl$ASSET.CLASS == asset.class, ]
    }

    if (!is.na(region)) {
        df <- df[df$REGION == region, ]
        if (region == "North America") {
            df <- df[!is.na(df$FUND.CURRENCY) & df$FUND.CURRENCY == "USD", ]
        }
    }

    df <- df[!duplicated(df$FUND.ID), ]

    # Create the frequency table
    tbl <- table(df$VINTAGE...INCEPTION.YEAR)

    # Convert table to the desired string format: "Year (Count)"
    # names(tbl) gets the years, and tbl itself contains the counts
    output_string <- paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", ")

    asset.class.region <- asset.class
    if (!is.na(region)) {
        asset.class.region <- paste0(asset.class, " (", region, ")")
    }

    msg <- paste0(
        "The ", asset.class.region, " sample contains ", length(unique(df$FUND.ID)),
        " distinct funds spreading over ", length(tbl), " vintage years."
    )
    # cat(paste(strwrap(msg, width = 80), collapse = "\n"), "\n\n")
    print(msg)

    msg <- paste0("The vintage year distribution is as follows: ", output_string, ".")
    # cat(paste(strwrap(msg, width = 80), collapse = "\n"), "\n")
    print(msg)

    # Region split
    if (!is.na(region)) {
        tbl <- table(df$PRIMARY.GEOGRAPHIC.FOCUS)
    } else {
        tbl <- table(df$REGION)
    }
    output_string <- paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", ")
    msg <- paste0("The region distribution is as follows: ", output_string, ".")
    # cat(paste(strwrap(msg, width = 80), collapse = "\n"), "\n")
    print(msg)

    # Strategy split (sub-strategy of the asst class)
    tbl <- table(df$STRATEGY)
    output_string <- paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", ")
    msg <- paste0("The strategy distribution is as follows: ", output_string, ".")
    # cat(paste(strwrap(msg, width = 80), collapse = "\n"), "\n")
    print(msg)

    # Currency split
    tbl <- table(df$FUND.CURRENCY)
    output_string <- paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", ")
    msg <- paste0("The fund currency distribution is as follows: ", output_string, ".")
    print(msg)
}
make.vintage.summary("Private Equity")
make.vintage.summary("Buyout")
make.vintage.summary("Venture Capital")
make.vintage.summary("Buyout", "North America")
make.vintage.summary("Venture Capital", "North America")


## ----core_methodology---------------------------------------------------------
make.preqin.df <- function(
    df = df.xl,
    fund.size.weighting = TRUE,
    acs.filter = "Fund of Funds",
    out.name = "FOF",
    region.filter = NA,
    vin.year.pfs = FALSE,
    final.nav.discount = 1) {
    col.to.keep <- c(
        "FUND.ID", "FUND.SIZE..USD.MN.", "VINTAGE...INCEPTION.YEAR",
        "TRANSACTION.DATE", "TRANSACTION.AMOUNT", "NET.CASHFLOW",
        "CUMULATIVE.CONTRIBUTION", "CUMULATIVE.DISTRIBUTION"
    )

    df <- df[df$TRANSACTION.TYPE == "Value", ]

    if (length(acs.filter) > 0) {
        test <- df[df$ASSET.CLASS %in% acs.filter, ]
        if (nrow(test) > 0) df <- test

        test <- df[df$STRATEGY %in% acs.filter, ]
        if (nrow(test) > 0) df <- test
    }

    if (!is.na(region.filter)) {
        df <- df[df$REGION %in% region.filter, ]
        if (region.filter == "North America") {
            print(paste("nrow1", nrow(df)))
            df <- df[!is.na(df$FUND.CURRENCY) & df$FUND.CURRENCY == "USD", ]
            print(paste("nrow2", nrow(df)))
        }
    }

    df <- df[, col.to.keep]

    for (col in c("FUND.ID")) {
        df[, col] <- as.factor(df[, col])
    }
    df$FUND.SIZE..USD.MN. <- ifelse(is.na(df$FUND.SIZE..USD.MN.),
        mean(df$FUND.SIZE..USD.MN., na.rm = TRUE),
        df$FUND.SIZE..USD.MN.
    )
    df <- df[order(df$FUND.ID, df$TRANSACTION.DATE), ]

    list.df <- list()
    for (fund.id in levels(df$FUND.ID)) {
        df.ss <- df[df$FUND.ID == fund.id, ]
        df.ss <- df.ss[order(df.ss$TRANSACTION.DATE), ]
        df.ss$CF <- c(df.ss$NET.CASHFLOW[1], diff(df.ss$NET.CASHFLOW))
        df.ss$CUMULATIVE.CONTRIBUTION <- c(df.ss$CUMULATIVE.CONTRIBUTION[1], diff(df.ss$CUMULATIVE.CONTRIBUTION))
        df.ss$CUMULATIVE.DISTRIBUTION <- c(df.ss$CUMULATIVE.DISTRIBUTION[1], diff(df.ss$CUMULATIVE.DISTRIBUTION))
        # for last date: CF = CF + NAV, i.e., regard final NAV as final distribution
        final.nav <- df.ss$TRANSACTION.AMOUNT[nrow(df.ss)] * final.nav.discount
        df.ss$CF[nrow(df.ss)] <- df.ss$CF[nrow(df.ss)] + final.nav

        year.diff <- min(as.numeric(format(df.ss$TRANSACTION.DATE, "%Y"))) - as.numeric(df.ss$VINTAGE...INCEPTION.YEAR[1])

        if (is.na(final.nav)) print(paste("Missing Final NAV:", fund.id))

        if ((abs(year.diff) < 2) & (!is.na(final.nav))) {
            list.df[[fund.id]] <- df.ss
        }
    }
    df <- data.frame(do.call(rbind, list.df))

    if (fund.size.weighting) {
        df$CF <- df$CF * df$FUND.SIZE..USD.MN.
        df$TRANSACTION.AMOUNT <- df$TRANSACTION.AMOUNT * df$FUND.SIZE..USD.MN. # aka NAV
        df$CUMULATIVE.CONTRIBUTION <- df$CUMULATIVE.CONTRIBUTION * df$FUND.SIZE..USD.MN.
        df$CUMULATIVE.DISTRIBUTION <- df$CUMULATIVE.DISTRIBUTION * df$FUND.SIZE..USD.MN.
    } else {
        # df$CF <- df$CF * mean(df$FUND.SIZE..USD.MN.) # why?
        df$CF <- df$CF
    }

    ten.mio <- 10 * 1000 * 1000
    df$CF <- df$CF / ten.mio
    df$TRANSACTION.AMOUNT <- df$TRANSACTION.AMOUNT / ten.mio # aka NAV
    df$CUMULATIVE.CONTRIBUTION <- df$CUMULATIVE.CONTRIBUTION / ten.mio
    df$CUMULATIVE.DISTRIBUTION <- df$CUMULATIVE.DISTRIBUTION / ten.mio
    df <- df[, c(
        "FUND.ID", "TRANSACTION.DATE", "CF",
        "TRANSACTION.AMOUNT", "CUMULATIVE.CONTRIBUTION", "CUMULATIVE.DISTRIBUTION",
        "VINTAGE...INCEPTION.YEAR"
    )]
    colnames(df) <- c("Fund.ID", "Date", "CF", "NAV", "CON", "DIS", "Vintage")
    df$type <- out.name

    # Vintage Year Porftolio Formation
    agg.vin.year.pf <- function(df) {
        Vintage <- df$Vintage[1]
        type <- df$type[1]
        df.out <- aggregate(CF ~ Date, data = df, mean)
        df.nav <- aggregate(NAV ~ Date, data = df, mean)
        df.out <- merge(df.out, df.nav, by = "Date")
        df.con <- aggregate(CON ~ Date, data = df, mean)
        df.out <- merge(df.out, df.con, by = "Date")
        df.dis <- aggregate(DIS ~ Date, data = df, mean)
        df.out <- merge(df.out, df.dis, by = "Date")
        df.out$type <- type
        df.out$Vintage <- Vintage
        df.out$Fund.ID <- paste0(Vintage, "_", type)
        df.out <- df.out[order(df.out$Date), ]
        return(df.out)
    }
    if (vin.year.pfs) {
        dfx <- split(df, df$Vintage)
        df <- data.frame(do.call(rbind, lapply(dfx, agg.vin.year.pf)))
    }

    return(df)
}
df <- make.preqin.df(acs.filter = "Private Equity", region.filter = "North America")
length(df$Fund.ID[!duplicated(df$Fund.ID)])


## ----execution----------------------------------------------------------------
# run ----

acs <- list(
    INF = "Infrastructure",
    NATRES = "Natural Resources",
    PD = "Private Debt",
    PE = "Private Equity",
    RE = "Real Estate",
    VC = "Venture Capital",
    MEZZ = "Mezzanine",
    DD = c("Special Situations", "Distressed Debt"),
    BO = "Buyout",
    GroBO = c("Growth", "Buyout"),
    FOF = "Fund of Funds",
    SEC = "Secondaries"
)


make.preqin.csv <- function(fund.size.weighting, vin.year.pfs, region.filter, final.nav.discount = 1) {
    l <- list()
    for (i in names(acs)) {
        print(acs[[i]])
        df0 <- make.preqin.df(
            fund.size.weighting = fund.size.weighting,
            acs.filter = acs[[i]],
            region.filter = region.filter,
            vin.year.pfs = vin.year.pfs,
            out.name = i,
            final.nav.discount = final.nav.discount
        )
        l[[i]] <- df0
    }
    df.out <- data.frame(do.call(rbind, l))

    tag <- ifelse(fund.size.weighting, "FW", "EW")
    tag <- ifelse(vin.year.pfs, paste0(tag, "_VYP"), tag)
    tag <- ifelse(is.na(region.filter), tag, paste0(tag, "_", region.filter))
    tag <- ifelse(final.nav.discount != 1, paste0(tag, "_NC", final.nav.discount * 100), tag)
    file <- paste0(data.prepared.folder, "/preqin_cashflows_2022_", tag, "_NAV.csv")
    file <- paste0(data.prepared.folder, "/preqin_cashflows_", tag, "_NAV.csv")
    write.csv(df.out, file, row.names = FALSE)
    invisible(df.out)
}
# Execute for various configurations
df.out <- make.preqin.csv(TRUE, TRUE, "North America")
df.out <- make.preqin.csv(FALSE, TRUE, "North America")
df.out <- make.preqin.csv(TRUE, TRUE, "North America", 0.5) # 50% final NAV discount
df.out <- make.preqin.csv(FALSE, TRUE, "North America", 0.5) # 50% final NAV discount

df.out <- make.preqin.csv(TRUE, TRUE, NA)
df.out <- make.preqin.csv(FALSE, TRUE, NA)
df.out.fw <- make.preqin.csv(TRUE, FALSE, NA)
df.out.ew <- make.preqin.csv(FALSE, FALSE, NA)
df.out <- make.preqin.csv(TRUE, TRUE, NA, 0.5) # 50% final NAV discount
df.out <- make.preqin.csv(FALSE, TRUE, NA, 0.5) # 50% final NAV discount

df.out <- make.preqin.csv(TRUE, TRUE, "Europe")
df.out <- make.preqin.csv(FALSE, TRUE, "Europe")



summary(df.out)


## ----plot_cum_cf, fig.width=10, fig.height=5, warning=FALSE, message=FALSE----
library(ggplot2)
library(dplyr)

# Use all Private Equity vintages available in the dataset

# Helper function to process the cashflow data into percentage of commitment proxy
process_plot_data <- function(df, weight_label, fund_type) {
    df.subset <- df[df$type == fund_type, ]

    df.subset$Date <- as.Date(df.subset$Date)
    df.subset$Year <- as.numeric(format(df.subset$Date, "%Y"))
    df.subset$Age <- df.subset$Year - as.numeric(as.character(df.subset$Vintage))

    # Filter up to year 20 (t+20)
    df.subset <- df.subset[df.subset$Age >= 0 & df.subset$Age <= 20, ]

    agg <- aggregate(cbind(CF, CON) ~ Vintage + Age, data = df.subset, sum, na.rm = TRUE)
    agg <- agg[order(agg$Vintage, agg$Age), ]

    # Calculate cumulative net cashflows and cumulative contributions
    agg$CumCF <- unlist(tapply(agg$CF, agg$Vintage, cumsum))
    agg$CumCON <- unlist(tapply(agg$CON, agg$Vintage, cumsum))

    # Drawdown proxy: maximum drawn capital for the vintage (maximum absolute Cumulative CON)
    max_con <- aggregate(CumCON ~ Vintage, data = agg, function(x) max(abs(x)))
    colnames(max_con)[2] <- "MaxDrawn"
    agg <- merge(agg, max_con, by = "Vintage")

    # Scale to percent of proxy commitment
    agg$CumCF_Pct <- agg$CumCF / agg$MaxDrawn
    agg$Weighting <- weight_label
    return(agg)
}

for (f_type in c("PE", "BO", "VC")) {
    agg.fw <- process_plot_data(df.out.fw, "Fund-size weighted", f_type)
    agg.ew <- process_plot_data(df.out.ew, "Equal weighted", f_type)

    # Combine for side-by-side plotting
    agg.all <- rbind(agg.fw, agg.ew)
    agg.all$Vintage <- as.factor(agg.all$Vintage)
    agg.all$Weighting <- factor(agg.all$Weighting, levels = c("Fund-size weighted", "Equal weighted"))

    # Create faceted plot
    p_jcurve <- ggplot(agg.all, aes(x = Age, y = CumCF_Pct, color = Vintage)) +
        geom_line(linewidth = 1) +
        facet_wrap(~Weighting) +
        scale_y_continuous(labels = scales::percent) +
        scale_x_continuous(breaks = 0:20, labels = c("Inception", paste0("t+", 1:20))) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            legend.position = "right",
            strip.text = element_text(size = 11, face = "bold"),
            panel.grid.minor = element_blank()
        ) +
        scale_color_viridis_d(option = "plasma", end = 0.9) +
        labs(
            title = paste0("Cumulative net cash flows since inception (all available ", f_type, " vintages in Preqin dataset)"),
            y = "",
            x = "",
            color = "Vintage"
        )

    print(p_jcurve)

    # Export the plot as PDF
    ggsave(paste0("cumcashflows", f_type, ".pdf"), p_jcurve, width = 10, height = 5)
}

