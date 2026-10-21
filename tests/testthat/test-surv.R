test_that("a Gehan-only selection reports the Gehan chi-squared", {
    # GIVEN the lung data grouped by sex
    data <- lungData()

    # WHEN only the Gehan test is requested
    results <- runSurv(data = data, groups = "sex", tests = "gehan")

    # THEN the Gehan row holds coin's Gehan statistic, not the log-rank one
    fit <- coinTest(completeCases(data, "sex"), "Gehan-Breslow")
    expected <- unname(coin::statistic(fit)^2)
    expect_equal(
        results$tests$getCell(rowNo = 1, "chisqr[gehan]")$value,
        expected[1]
    )
})

test_that("the Gehan chi-squared matches coin for more than two groups", {
    # GIVEN the lung data grouped by its three ECOG levels
    data <- lungData()

    # WHEN the Gehan test is requested
    results <- runSurv(data = data, groups = "ecog", tests = "gehan")

    # THEN the chi-squared is coin's k-sample statistic, not its square
    fit <- coinTest(completeCases(data, "ecog"), "Gehan-Breslow")
    expect_equal(
        results$tests$getCell(rowNo = 1, "chisqr[gehan]")$value,
        unname(coin::statistic(fit))
    )
})

test_that("the Tarone-Ware chi-squared matches coin for more than two groups", {
    # GIVEN the lung data grouped by its three ECOG levels
    data <- lungData()

    # WHEN the Tarone-Ware test is requested
    results <- runSurv(data = data, groups = "ecog", tests = "taroneware")

    # THEN the chi-squared is coin's k-sample statistic, not its square
    fit <- coinTest(completeCases(data, "ecog"), "Tarone-Ware")
    expect_equal(
        results$tests$getCell(rowNo = 1, "chisqr[taroneware]")$value,
        unname(coin::statistic(fit))
    )
})

test_that("the Gehan p-value matches coin for more than two groups", {
    # GIVEN the lung data grouped by its three ECOG levels
    data <- lungData()

    # WHEN the Gehan test is requested
    results <- runSurv(data = data, groups = "ecog", tests = "gehan")

    # THEN the p-value equals coin's
    fit <- coinTest(completeCases(data, "ecog"), "Gehan-Breslow")
    expect_equal(
        results$tests$getCell(rowNo = 1, "pvalue[gehan]")$value,
        as.numeric(coin::pvalue(fit))
    )
})

test_that("pairwise p-values are Bonferroni-corrected", {
    # GIVEN the lung data grouped by its three ECOG levels
    data <- lungData()

    # WHEN pairwise log-rank comparisons are requested
    results <- runSurv(
        data = data,
        groups = "ecog",
        tests = "logrank",
        testspw = TRUE
    )

    # THEN each p-value is the raw p-value times the three comparisons
    raw <- envStatsP(completeCases(data, "ecog"), "e0", "e2")
    expect_equal(
        results$testspw$getCell(rowKey = c("e0", "e2"), "ppw[logrank]")$value,
        min(raw * 3, 1)
    )
})

test_that("numeric 1/2 event coding gives the same pairwise results as 0/1", {
    # GIVEN the same events coded as 0/1 and as 1/2
    args <- list(
        eventLevel = NULL,
        groups = "ecog",
        tests = "logrank",
        testspw = TRUE
    )

    # WHEN the analysis is run with each coding
    results01 <- do.call(runSurv, c(args, event = "status01"))
    results12 <- do.call(runSurv, c(args, event = "status12"))

    # THEN the pairwise comparisons are identical
    expect_equal(results12$testspw$asDF, results01$testspw$asDF)
})

test_that("numeric event variables not coded 0/1 or 1/2 are rejected", {
    # GIVEN an event variable that also contains the value 2 alongside 0/1
    data <- lungData()
    data$status012 <- data$status01
    data$status012[1:5] <- 2

    # WHEN the analysis is run on it
    run <- function() {
        runSurv(data = data, event = "status012", eventLevel = NULL)
    }

    # THEN it is rejected with an explanation of the coding
    expect_error(run(), "coded 0 \\(censored\\) and 1 \\(event\\)")
})

test_that("a factor event variable without an event level is rejected", {
    # GIVEN the lung data with a factor event variable

    # WHEN the analysis is run without choosing the event level
    run <- function() runSurv(eventLevel = NULL)

    # THEN it asks for the event level
    expect_error(run(), "choose the event level")
})

test_that("a grouping variable with a single level is rejected", {
    # GIVEN a grouping variable with only one level
    data <- lungData()
    data$one <- factor(rep("a", nrow(data)))

    # WHEN the analysis is run with it
    run <- function() runSurv(data = data, groups = "one", tests = "logrank")

    # THEN it is rejected
    expect_error(run(), "at least two levels")
})

test_that("negative times are rejected", {
    # GIVEN a few negative survival times
    data <- lungData()
    data$time[1:3] <- -5

    # WHEN the analysis is run
    run <- function() runSurv(data = data)

    # THEN it is rejected
    expect_error(run(), "cannot be negative")
})

test_that("summary rows are filled by group when a level has no data", {
    # GIVEN an ECOG level, listed second, with no complete cases
    data <- lungDataWithEmptyLevel()

    # WHEN the analysis is run grouped by ECOG
    results <- runSurv(data = data, groups = "ecog")

    # THEN each level's row holds that level's own count
    tmpDat <- completeCases(data, "ecog")
    counts <- sapply(c("e0", "e1", "e2"), function(level) {
        unname(results$summary$getCell(rowKey = level, "n")$value)
    })
    expect_equal(counts, c(table(tmpDat$group)[c("e0", "e1", "e2")]))
})

test_that("median rows are filled by group when a level has no data", {
    # GIVEN an ECOG level, listed second, with no complete cases
    data <- lungDataWithEmptyLevel()

    # WHEN the analysis is run grouped by ECOG
    results <- runSurv(data = data, groups = "ecog")

    # THEN each level's row holds that level's own median
    fit <- survival::survfit(
        survival::Surv(times, status) ~ group,
        data = completeCases(data, "ecog")
    )
    expected <- summary(fit)$table[, "median"]
    medians <- sapply(c("e0", "e1", "e2"), function(level) {
        results$medianestimates$getCell(rowKey = level, "median")$value
    })
    expect_equal(unname(medians), unname(expected))
})

test_that("a level with no data is left empty", {
    # GIVEN an ECOG level with no complete cases
    data <- lungDataWithEmptyLevel()

    # WHEN the analysis is run grouped by ECOG
    results <- runSurv(data = data, groups = "ecog")

    # THEN that level's summary row has no values
    expect_null(results$summary$getCell(rowKey = "e9", "n")$value)
})

test_that("pairwise comparisons skip a level with no data", {
    # GIVEN an ECOG level with no complete cases
    data <- lungDataWithEmptyLevel()

    # WHEN pairwise comparisons are requested
    results <- runSurv(
        data = data,
        groups = "ecog",
        tests = "logrank",
        testspw = TRUE
    )

    # THEN comparisons involving that level are not computed
    expect_null(
        results$testspw$getCell(rowKey = c("e0", "e9"), "ppw[logrank]")$value
    )
})

test_that("the Bonferroni count only includes comparisons that are made", {
    # GIVEN four ECOG levels, one of them without complete cases
    data <- lungDataWithEmptyLevel()

    # WHEN pairwise comparisons are requested
    results <- runSurv(
        data = data,
        groups = "ecog",
        tests = "logrank",
        testspw = TRUE
    )

    # THEN p-values are corrected for 3 comparisons, not 6
    raw <- envStatsP(completeCases(data, "ecog"), "e0", "e2")
    expect_equal(
        results$testspw$getCell(rowKey = c("e0", "e2"), "ppw[logrank]")$value,
        min(raw * 3, 1)
    )
})

test_that("pairwise comparisons explain that a test is needed", {
    # GIVEN the lung data grouped by ECOG

    # WHEN pairwise comparisons are requested without any test
    results <- runSurv(groups = "ecog", testspw = TRUE)

    # THEN the table notes that a test must be selected
    expect_match(
        results$testspw$notes$notests$note,
        "Select at least one test"
    )
})

test_that("no pairwise rows are built when pairwise comparisons are off", {
    # GIVEN the lung data grouped by ECOG

    # WHEN the analysis is run without pairwise comparisons
    results <- runSurv(groups = "ecog", tests = "logrank", testspw = FALSE)

    # THEN the pairwise table has no rows
    expect_length(results$testspw$rowKeys, 0)
})

test_that("a pairwise row is built for each pair of groups", {
    # GIVEN the lung data grouped by its three ECOG levels

    # WHEN pairwise comparisons are requested
    results <- runSurv(groups = "ecog", tests = "logrank", testspw = TRUE)

    # THEN there is one row per pair
    expect_length(results$testspw$rowKeys, 3)
})

test_that("plots do not store the data in their state", {
    # GIVEN the lung data grouped by sex

    # WHEN both plots are requested
    results <- runSurv(groups = "sex", sc = TRUE, chf = TRUE)

    # THEN neither plot carries state
    expect_null(results$sc$state)
    expect_null(results$chf$state)
})

test_that("plots render from the dataset", {
    # GIVEN an analysis with both plots and the risk table requested
    results <- runSurv(groups = "sex", sc = TRUE, chf = TRUE, risk = TRUE)
    paths <- c(tempfile(fileext = ".png"), tempfile(fileext = ".png"))

    # WHEN the plots are saved
    results$sc$saveAs(paths[1])
    results$chf$saveAs(paths[2])

    # THEN both image files are written
    expect_true(all(file.exists(paths)))
})
