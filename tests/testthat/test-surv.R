test_that("overall tests are matched to the test by name, not position", {
    data <- lungData()

    r <- deathwatch::surv(
        data = data,
        elapsed = "time",
        event = "event",
        eventLevel = "dead",
        groups = "sex",
        tests = "gehan",
        sc = FALSE)

    tmpDat <- na.omit(data.frame(
        times = data$time,
        status = as.integer(data$event == "dead"),
        group = data$sex))
    expected <- coin::statistic(coin::logrank_test(
        survival::Surv(times, status) ~ group,
        data = tmpDat,
        type = "Gehan-Breslow"))^2

    table <- r$tests
    expect_equal(table$getCell(rowNo = 1, "chisqr[gehan]")$value, unname(expected[1]))
    expect_equal(table$getCell(rowNo = 1, "df[gehan]")$value, 1)
    expect_equal(
        table$getCell(rowNo = 1, "pvalue[gehan]")$value,
        pchisq(unname(expected[1]), 1, lower.tail = FALSE))
})

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
