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

test_that("Gehan and Tarone-Ware agree with coin for more than two groups", {
    data <- lungData()

    r <- deathwatch::surv(
        data = data,
        elapsed = "time",
        event = "event",
        eventLevel = "dead",
        groups = "ecog",
        tests = c("gehan", "taroneware"),
        sc = FALSE)

    tmpDat <- na.omit(data.frame(
        times = data$time,
        status = as.integer(data$event == "dead"),
        group = data$ecog))

    table <- r$tests
    for (test in c("gehan", "taroneware")) {
        type <- c(gehan = "Gehan-Breslow", taroneware = "Tarone-Ware")[[test]]
        fit <- coin::logrank_test(
            survival::Surv(times, status) ~ group,
            data = tmpDat,
            type = type)

        expect_equal(table$getCell(rowNo = 1, paste0("chisqr[", test, "]"))$value, unname(coin::statistic(fit)))
        expect_equal(table$getCell(rowNo = 1, paste0("df[", test, "]"))$value, 2)
        expect_equal(table$getCell(rowNo = 1, paste0("pvalue[", test, "]"))$value, as.numeric(coin::pvalue(fit)))
    }
})

test_that("numeric 1/2 event coding gives the same results as 0/1 coding", {
    data <- lungData()

    args <- list(
        data = data,
        elapsed = "time",
        eventLevel = NULL,
        groups = "ecog",
        tests = c("logrank", "gehan"),
        testspw = TRUE,
        sc = FALSE)

    r01 <- do.call(deathwatch::surv, c(args, event = "status01"))
    r12 <- do.call(deathwatch::surv, c(args, event = "status12"))

    expect_equal(r12$summary$asDF, r01$summary$asDF)
    expect_equal(r12$medianestimates$asDF, r01$medianestimates$asDF)
    expect_equal(r12$tests$asDF, r01$tests$asDF)
    expect_equal(r12$testspw$asDF, r01$testspw$asDF)
})

test_that("numeric event variables not coded 0/1 or 1/2 are rejected", {
    data <- lungData()
    data$status012 <- data$status01
    data$status012[1:5] <- 2

    expect_error(
        deathwatch::surv(
            data = data,
            elapsed = "time",
            event = "status012",
            eventLevel = NULL,
            groups = NULL,
            tests = NULL,
            sc = FALSE),
        "coded 0 \\(censored\\) and 1 \\(event\\)")
})

test_that("a grouping variable with a single level is rejected", {
    data <- lungData()
    data$one <- factor(rep("a", nrow(data)))

    expect_error(
        deathwatch::surv(
            data = data,
            elapsed = "time",
            event = "event",
            eventLevel = "dead",
            groups = "one",
            tests = "logrank",
            sc = FALSE),
        "at least two levels")
})
