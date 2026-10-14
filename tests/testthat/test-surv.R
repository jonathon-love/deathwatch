test_that("pairwise p-values are Bonferroni-corrected", {
    data <- lungData()

    r <- deathwatch::surv(
        data = data,
        elapsed = "time",
        event = "event",
        eventLevel = "dead",
        groups = "ecog",
        tests = "logrank",
        testspw = TRUE,
        sc = FALSE)

    tmpDat <- na.omit(data.frame(
        times = data$time,
        status = as.integer(data$event == "dead"),
        group = data$ecog))
    x <- tmpDat[tmpDat$group == "e0", ]
    y <- tmpDat[tmpDat$group == "e2", ]
    raw <- EnvStats::twoSampleLinearRankTestCensored(
        test = "logrank",
        x = x$times,
        x.censored = x$status == 0,
        y = y$times,
        y.censored = y$status == 0,
        censoring.side = "right")$p.value

    table <- r$testspw
    p <- table$getCell(rowKey = c("e0", "e2"), "ppw[logrank]")$value
    expect_equal(p, min(raw * 3, 1))
})
