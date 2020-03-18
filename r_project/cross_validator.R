
df <- read.csv("data_prepared/preqin_cashflows_FW.csv")

vintage.blocks <- list()
for(year in seq(1988, 2015, by=3)) {
  validate <- year:(year + 2)
  block <- (year - 3):(year - 1)
  block <- c(block, (year + 3):(year + 5))
  
  estimate <- 1980:2020
  estimate <- estimate[!(estimate %in% validate)]
  estimate <- estimate[!(estimate %in% block)]
  
  vintage.blocks[[paste(year)]] <- list(validate=validate, block=block, estimate=estimate)
}
for(year in names(vintage.blocks)) {
  for(key in c("estimate", "validate")) {
    vintage.blocks[[year]][[paste0("df.", key)]] <- df[df$Vintage %in% vintage.blocks[[year]][[key]], ]
  }
}

