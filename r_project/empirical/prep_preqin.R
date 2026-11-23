## ----setup, message=FALSE, warning=FALSE--------------------------------------
## prepare preqin
# prep ----
if (sys.nframe() == 0L) rm(list = ls())

library(readxl)
library(here)

source(here("helper/pool_vintage.R"))

# 2026_c: constant-composition VYPs and commitment-normalized FW pools
# (revision plan B4/B6). 2026_b holds the earlier reporter-mean pools.
data.prepared.folder <- here("empirical/data_prepared_2026_c")

# Sample rule (main sample): a fund has a COMPLETE HISTORY if its last Preqin
# record lies more than min.history.years after January 1 of its vintage year,
# OR if it is still reporting at the data end (last record within
# live.window.months of the export's last NAV date). A young fund that is still
# reporting is censored by the data end, not truncated; its terminal NAV is
# treated like that of every live fund (NAV haircut exercise).
# In the April 2022 export, 218 of 1,033 North America BO/VC funds (V <= 2010)
# have histories that stop within 5 years, mostly a vendor artifact (the Feb
# 2020 export has later records for these funds). Booking their last NAV as a
# liquidating distribution is not defensible, so the main sample excludes
# them. For V <= 2010 the live clause adds no fund (all live funds are older
# than 5 years); it keeps recent vintages usable for the vintage-cutoff check.
# Files in data.prepared.folder hold the main sample; files with all funds
# (flagged by the Complete column) go to data.allfunds.folder and serve as a
# robustness sample. A funnel table records every exclusion.
min.history.years <- 5
live.window.months <- 12 # Preqin NAVs lag 1-3 quarters; most live funds' latest NAV is 2021-06 or 2021-09 vs. data end 2022-03
data.allfunds.folder <- paste0(data.prepared.folder, "_allfunds")
if (!dir.exists(data.allfunds.folder)) dir.create(data.allfunds.folder)


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

# Data end = last NAV (Value) date in the export; funds with a last record on
# or after live.cutoff count as still reporting (sample rule above).
data.end <- as.Date(max(df.xl$TRANSACTION.DATE[df.xl$TRANSACTION.TYPE %in% "Value"], na.rm = TRUE))
live.cutoff <- seq(data.end + 1, by = paste0("-", live.window.months, " months"), length.out = 2)[2] - 1
print(paste("Data end:", data.end, "| live if last record >=", live.cutoff))

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
# Impute missing fund sizes (revision fix, Nov 2026):
#   - "vintage.median" (baseline): observed same-vintage median across funds;
#     if a vintage has no observed size (e.g. singleton 1982 VC), expand
#     symmetrically by calendar distance d = 1, 2, ... until observed sizes
#     exist and take their median (fallback 1 if a strategy has none).
#   - "pooled.mean" (sensitivity / legacy): pooled cross-vintage mean over funds.
#   - "drop" (sensitivity): exclude funds with missing sizes.
# Always records Size.Imputed (logical) before imputation.
impute.fund.sizes <- function(df, size.imputation = c("vintage.median", "pooled.mean", "drop")) {
    size.imputation <- match.arg(size.imputation)
    df$Size.Imputed <- is.na(df$FUND.SIZE..USD.MN.)
    if (size.imputation == "drop") {
        return(df[!df$Size.Imputed, , drop = FALSE])
    }
    u <- df[!duplicated(df$FUND.ID), c("FUND.ID", "VINTAGE...INCEPTION.YEAR", "FUND.SIZE..USD.MN."), drop = FALSE]
    if (size.imputation == "pooled.mean") {
        fill <- mean(u$FUND.SIZE..USD.MN., na.rm = TRUE)
        if (!is.finite(fill)) fill <- 1
        df$FUND.SIZE..USD.MN.[df$Size.Imputed] <- fill
        return(df)
    }
    # vintage.median with expanding nearby-vintage fallback
    obs <- u[!is.na(u$FUND.SIZE..USD.MN.), , drop = FALSE]
    miss.v <- unique(df$VINTAGE...INCEPTION.YEAR[df$Size.Imputed])
    for (v in miss.v) {
        vals <- obs$FUND.SIZE..USD.MN.[obs$VINTAGE...INCEPTION.YEAR == v]
        if (!length(vals) && nrow(obs) > 0L) {
            dists <- abs(obs$VINTAGE...INCEPTION.YEAR - v)
            vals <- obs$FUND.SIZE..USD.MN.[dists == min(dists)]
        }
        fill <- if (length(vals)) stats::median(vals) else 1
        idx <- df$Size.Imputed & df$VINTAGE...INCEPTION.YEAR == v
        df$FUND.SIZE..USD.MN.[idx] <- fill
    }
    df
}

make.preqin.df <- function(
    df = df.xl,
    fund.size.weighting = TRUE,
    acs.filter = "Fund of Funds",
    out.name = "FOF",
    region.filter = NA,
    vin.year.pfs = FALSE,
    final.nav.discount = 1,
    size.imputation = "vintage.median") {
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
    df <- impute.fund.sizes(df, size.imputation = size.imputation)
    df <- df[order(df$FUND.ID, df$TRANSACTION.DATE), ]

    list.df <- list()
    list.status <- list()
    for (df.ss in split(df, df$FUND.ID, drop = TRUE)) {
        fund.id <- as.character(df.ss$FUND.ID[1])
        df.ss <- df.ss[order(df.ss$TRANSACTION.DATE), ]
        df.ss$CF <- c(df.ss$NET.CASHFLOW[1], diff(df.ss$NET.CASHFLOW))
        df.ss$CUMULATIVE.CONTRIBUTION <- c(df.ss$CUMULATIVE.CONTRIBUTION[1], diff(df.ss$CUMULATIVE.CONTRIBUTION))
        df.ss$CUMULATIVE.DISTRIBUTION <- c(df.ss$CUMULATIVE.DISTRIBUTION[1], diff(df.ss$CUMULATIVE.DISTRIBUTION))
        # for last date: CF = CF + NAV, i.e., regard final NAV as final distribution
        final.nav <- df.ss$TRANSACTION.AMOUNT[nrow(df.ss)] * final.nav.discount
        df.ss$CF[nrow(df.ss)] <- df.ss$CF[nrow(df.ss)] + final.nav

        year.diff <- min(as.numeric(format(df.ss$TRANSACTION.DATE, "%Y"))) - as.numeric(df.ss$VINTAGE...INCEPTION.YEAR[1])

        if (is.na(final.nav)) print(paste("Missing Final NAV:", fund.id))

        vintage <- as.numeric(df.ss$VINTAGE...INCEPTION.YEAR[1])
        last.date <- as.Date(max(df.ss$TRANSACTION.DATE))
        age.last <- as.numeric(last.date - as.Date(paste0(vintage, "-01-01"))) / 365.25
        df.ss$Last.Date <- last.date
        df.ss$Age.Last <- age.last
        complete <- (age.last > min.history.years) || (last.date >= live.cutoff)
        df.ss$Complete <- complete

        keep <- (abs(year.diff) < 2) & (!is.na(final.nav))
        if (keep) list.df[[fund.id]] <- df.ss
        list.status[[fund.id]] <- data.frame(
            Fund.ID = fund.id, Vintage = vintage,
            status = if (is.na(final.nav)) {
                "missing.nav"
            } else if (abs(year.diff) >= 2) {
                "late.start"
            } else if (complete) {
                "complete"
            } else {
                "truncated"
            },
            stringsAsFactors = FALSE
        )
    }
    df <- data.frame(do.call(rbind, list.df))
    df.status <- data.frame(do.call(rbind, list.status))

    # Preqin reports cash flows per USD 10mn of commitment, so dividing by
    # 10mn expresses EW flows per unit of commitment. FW fund-level flows are
    # additionally multiplied by fund size (USD mn), i.e. they are in USD mn;
    # the FW vintage pools are normalized by total vintage commitment in
    # pool.vintage() (B6), so both EW and FW pools are per unit of commitment.
    if (fund.size.weighting) {
        df$CF <- df$CF * df$FUND.SIZE..USD.MN.
        df$TRANSACTION.AMOUNT <- df$TRANSACTION.AMOUNT * df$FUND.SIZE..USD.MN. # aka NAV
        df$CUMULATIVE.CONTRIBUTION <- df$CUMULATIVE.CONTRIBUTION * df$FUND.SIZE..USD.MN.
        df$CUMULATIVE.DISTRIBUTION <- df$CUMULATIVE.DISTRIBUTION * df$FUND.SIZE..USD.MN.
    }

    ten.mio <- 10 * 1000 * 1000
    df$CF <- df$CF / ten.mio
    df$TRANSACTION.AMOUNT <- df$TRANSACTION.AMOUNT / ten.mio # aka NAV
    df$CUMULATIVE.CONTRIBUTION <- df$CUMULATIVE.CONTRIBUTION / ten.mio
    df$CUMULATIVE.DISTRIBUTION <- df$CUMULATIVE.DISTRIBUTION / ten.mio
    df <- df[, c(
        "FUND.ID", "TRANSACTION.DATE", "CF",
        "TRANSACTION.AMOUNT", "CUMULATIVE.CONTRIBUTION", "CUMULATIVE.DISTRIBUTION",
        "VINTAGE...INCEPTION.YEAR", "FUND.SIZE..USD.MN.",
        "Last.Date", "Age.Last", "Complete", "Size.Imputed"
    )]
    colnames(df) <- c(
        "Fund.ID", "Date", "CF", "NAV", "CON", "DIS", "Vintage", "Size",
        "Last.Date", "Age.Last", "Complete", "Size.Imputed"
    )
    df$Fund.ID <- as.character(df$Fund.ID)
    df$Date <- as.Date(df$Date)
    df$type <- out.name
    df <- df[, c(
        "Fund.ID", "Date", "CF", "NAV", "CON", "DIS", "Vintage", "type", "Size",
        "Last.Date", "Age.Last", "Complete", "Size.Imputed"
    )]

    # Vintage year portfolio formation: constant composition (B4, B6)
    if (vin.year.pfs) df <- pool.all.vintages(df, fund.size.weighting)

    df.status$type <- out.name
    attr(df, "status") <- df.status
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


# Builds the fund-level panel once and writes, from that single panel,
#   - the main sample (complete-history funds): fund-level file (if
#     write.fund.level) and constant-composition VYP file in data.prepared.folder;
#   - the all-funds robustness sample in data.allfunds.folder: VYP file always,
#     fund-level file and (for EW) the sample funnel table if write.fund.level.
# Fund-level and VYP files of a sample are always mutually consistent (B4 test).
make.preqin.csv <- function(fund.size.weighting, region.filter, final.nav.discount = 1,
                            write.fund.level = TRUE, size.imputation = "vintage.median",
                            out.main = data.prepared.folder, out.all = data.allfunds.folder) {
    if (!dir.exists(out.main)) dir.create(out.main, recursive = TRUE)
    if (!is.null(out.all) && !dir.exists(out.all)) dir.create(out.all, recursive = TRUE)
    l <- list()
    l.status <- list()
    for (i in names(acs)) {
        print(acs[[i]])
        l[[i]] <- make.preqin.df(
            fund.size.weighting = fund.size.weighting,
            acs.filter = acs[[i]],
            region.filter = region.filter,
            vin.year.pfs = FALSE,
            out.name = i,
            final.nav.discount = final.nav.discount,
            size.imputation = size.imputation
        )
        l.status[[i]] <- attr(l[[i]], "status")
    }
    df.all <- data.frame(do.call(rbind, l))
    rownames(df.all) <- NULL

    pool <- function(df.fund) {
        df.vyp <- do.call(rbind, lapply(
            split(df.fund, df.fund$type, drop = TRUE),
            pool.all.vintages,
            fund.size.weighting = fund.size.weighting
        ))
        df.vyp <- df.vyp[order(match(df.vyp$type, names(acs)), df.vyp$Vintage, df.vyp$Date), ]
        rownames(df.vyp) <- NULL
        df.vyp
    }
    region.tag <- ifelse(is.na(region.filter), "", paste0("_", region.filter))
    make.file <- function(folder, vin.year.pfs) {
        tag <- ifelse(fund.size.weighting, "FW", "EW")
        tag <- ifelse(vin.year.pfs, paste0(tag, "_VYP"), tag)
        tag <- paste0(tag, region.tag)
        tag <- ifelse(final.nav.discount != 1, paste0(tag, "_NC", final.nav.discount * 100), tag)
        paste0(folder, "/preqin_cashflows_", tag, "_NAV.csv")
    }

    # main sample: complete-history funds
    df.main <- df.all[df.all$Complete, ]
    if (write.fund.level) write.csv(df.main, make.file(out.main, FALSE), row.names = FALSE)
    write.csv(pool(df.main), make.file(out.main, TRUE), row.names = FALSE)

    # robustness sample: all funds (truncated histories flagged Complete = FALSE).
    # VYP files are written for every configuration (incl. NAV haircuts), so
    # every scenario also runs on the all-funds folder.
    if (!is.null(out.all)) {
        write.csv(pool(df.all), make.file(out.all, TRUE), row.names = FALSE)
        if (write.fund.level) {
            write.csv(df.all, make.file(out.all, FALSE), row.names = FALSE)
            if (!fund.size.weighting && size.imputation == "vintage.median") {
                st <- do.call(rbind, l.status)
                st$status <- factor(st$status, levels = c("complete", "truncated", "late.start", "missing.nav"))
                tab <- as.data.frame.matrix(table(paste(st$type, st$Vintage), st$status))
                funnel <- data.frame(
                    type = sub(" .*", "", rownames(tab)), Vintage = as.integer(sub(".* ", "", rownames(tab))),
                    N.raw = rowSums(tab), N.late.start = tab$late.start, N.missing.nav = tab$missing.nav,
                    N.truncated = tab$truncated, N.complete = tab$complete
                )
                funnel <- funnel[order(match(funnel$type, names(acs)), funnel$Vintage), ]
                write.csv(funnel, paste0(out.main, "/preqin_sample_funnel", region.tag, ".csv"), row.names = FALSE)
            }
        }
    }
    invisible(NULL)
}

# Execute for various configurations (fund-level files only without NAV haircut)
# Baseline FW imputation: within-vintage median with nearby-vintage fallback.
for (fsw in c(TRUE, FALSE)) {
    for (nc in c(1, 0, 0.25, 0.5, 0.75)) {
        make.preqin.csv(fsw, "North America", nc, write.fund.level = (nc == 1))
    }
    for (nc in c(1, 0.5)) {
        make.preqin.csv(fsw, NA, nc, write.fund.level = (nc == 1))
    }
    make.preqin.csv(fsw, "Europe")
}

# FW fund-size imputation sensitivities (North America, main sample, nc = 1):
# pooled cross-vintage mean (legacy rule) and dropping missing-size funds.
make.preqin.csv(TRUE, "North America", 1, write.fund.level = TRUE, size.imputation = "pooled.mean",
                out.main = paste0(data.prepared.folder, "_fwmean"), out.all = NULL)
make.preqin.csv(TRUE, "North America", 1, write.fund.level = TRUE, size.imputation = "drop",
                out.main = paste0(data.prepared.folder, "_fwdrop"), out.all = NULL)

# Checked-in aggregate diagnostics of fund-size missingness & imputed vintage shares
# (North America complete-history sample; aggregate counts/shares only, no fund IDs)
local({
    f.med <- read.csv(file.path(data.prepared.folder, "preqin_cashflows_FW_North America_NAV.csv"), stringsAsFactors = FALSE)
    f.mean <- read.csv(file.path(paste0(data.prepared.folder, "_fwmean"), "preqin_cashflows_FW_North America_NAV.csv"), stringsAsFactors = FALSE)
    u.med <- f.med[!duplicated(paste(f.med$type, f.med$Fund.ID)), c("type", "Vintage", "Fund.ID", "Size", "Size.Imputed")]
    u.mean <- f.mean[!duplicated(paste(f.mean$type, f.mean$Fund.ID)), c("type", "Vintage", "Fund.ID", "Size")]
    u <- merge(u.med, u.mean, by = c("type", "Vintage", "Fund.ID"), suffixes = c(".med", ".mean"))
    rows <- lapply(split(u, paste(u$type, u$Vintage)), function(x) {
        tot.med <- sum(x$Size.med); tot.mean <- sum(x$Size.mean)
        data.frame(
            type = x$type[1], Vintage = x$Vintage[1],
            N.complete = nrow(x), N.size.missing = sum(x$Size.Imputed),
            Observed.Median = if (any(!x$Size.Imputed)) stats::median(x$Size.med[!x$Size.Imputed]) else NA_real_,
            Imputed.Median.Size = if (any(x$Size.Imputed)) x$Size.med[x$Size.Imputed][1] else NA_real_,
            Imputed.Mean.Size = if (any(x$Size.Imputed)) x$Size.mean[x$Size.Imputed][1] else NA_real_,
            Share.Imputed.Median = if (tot.med > 0) sum(x$Size.med[x$Size.Imputed]) / tot.med else NA_real_,
            Share.Imputed.Mean = if (tot.mean > 0) sum(x$Size.mean[x$Size.Imputed]) / tot.mean else NA_real_,
            stringsAsFactors = FALSE
        )
    })
    diag <- do.call(rbind, rows)
    diag <- diag[order(match(diag$type, names(acs)), diag$Vintage), ]
    rownames(diag) <- NULL
    write.csv(diag, file.path(data.prepared.folder, "preqin_size_imputation_North America.csv"), row.names = FALSE)
})
# all funds (for the cumulative cash-flow plots below)
df.out.fw <- read.csv(paste0(data.allfunds.folder, "/preqin_cashflows_FW_NAV.csv"))
df.out.ew <- read.csv(paste0(data.allfunds.folder, "/preqin_cashflows_EW_NAV.csv"))

summary(df.out.ew)


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

