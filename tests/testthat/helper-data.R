lungData <- function() {
    lung <- survival::lung
    data.frame(
        time = as.numeric(lung$time),
        event = factor(ifelse(lung$status == 2, "dead", "censored")),
        status01 = as.integer(lung$status == 2),
        status12 = as.integer(lung$status),
        ecog = factor(lung$ph.ecog, levels = 0:2, labels = c("e0", "e1", "e2")),
        sex = factor(lung$sex, levels = 1:2, labels = c("male", "female"))
    )
}

# lung data where level e9 of ecog exists but has no complete cases
lungDataWithEmptyLevel <- function() {
    data <- lungData()
    data$ecog <- factor(data$ecog, levels = c("e0", "e9", "e1", "e2"))
    data$ecog[1:3] <- "e9"
    data$time[1:3] <- NA
    data
}

# runs surv() with defaults suited to the lung data; do.call() keeps
# argument values from being read as column names
runSurv <- function(...) {
    args <- list(
        data = lungData(),
        elapsed = "time",
        event = "event",
        eventLevel = "dead",
        groups = NULL,
        tests = NULL,
        sc = FALSE
    )
    overrides <- list(...)
    args[names(overrides)] <- overrides
    do.call(deathwatch::surv, args)
}

# the complete cases the analysis works with, coded 0/1
completeCases <- function(data, group) {
    tmpDat <- na.omit(data.frame(
        times = data$time,
        status = as.integer(data$event == "dead"),
        group = data[[group]]
    ))
    tmpDat$group <- droplevels(tmpDat$group)
    tmpDat
}

coinTest <- function(data, type) {
    coin::logrank_test(
        survival::Surv(times, status) ~ group,
        data = data,
        type = type
    )
}

envStatsP <- function(data, x, y) {
    x <- data[data$group == x, ]
    y <- data[data$group == y, ]
    EnvStats::twoSampleLinearRankTestCensored(
        test = "logrank",
        x = x$times,
        x.censored = x$status == 0,
        y = y$times,
        y.censored = y$status == 0,
        censoring.side = "right"
    )$p.value
}
