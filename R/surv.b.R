

#' @importFrom R6 R6Class
#' @importFrom jmvcore toNumeric
#' @importFrom jmvcore .
survClass <- R6::R6Class(
    "survClass",
    inherit = survBase,
    private = list(
        .init = function() {

            groups <- private$.groups()


            ####################################
            ### Add levels dinamically to first
            ### column in Table 1
            ### Events Summary
            ####################################
            summary <- self$results$summary
            for (group in groups)
                summary$addRow(rowKey=group, list(group=group))

            ####################################
            ### Add levels dinamically to first
            ### column in Table 3
            ### Median Estimates
            ####################################
            summary <- self$results$medianestimates
            for (group in groups)
                summary$addRow(rowKey=group, list(group=group))

            tests <- self$results$tests
            if ( ! is.null(self$options$groups)) {
                tests$addRow(rowKey=1)
            }
            else {
                tests$setVisible(FALSE)
            }


            ####################################
            ### Add levels dinamically to first
            ### column in Table with pairwise
            ### comparisons
            ####################################
            testspw <- self$results$testspw
            if ( ! is.null(self$options$groups)) {
                if (length(groups) >= 2) {
                    comparisons <- combn(groups, 2)
                    for (i in seq_len(ncol(comparisons))) {
                        key = comparisons[,i]
                        testspw$addRow(rowKey=key)
                    }
                }
            }
            else {
                testspw$setVisible(FALSE)
            }


            ################################
            ### Change the size of the
            ### survival plot if
            ### risk table option is ticked
            ################################
            if (self$options$risk)
                self$results$sc$setSize(600, 600)

        },
        .groups = function() {
            if (is.null(self$options$groups))
                group <- NULL
            else
                group <- self$data[[self$options$groups]]

            groups <- levels(group)
            if (length(groups) == 0)
                groups = ''

            groups
        },
        .buildData = function() {

            eventVarName <- self$options$event
            elapsedVarName <- self$options$elapsed

            if (nrow(self$data) == 0)
                jmvcore::reject(.('Data contains no (complete) rows'))

            if ( ! is.null(self$options$groups) && length(private$.groups()) < 2)
                jmvcore::reject(.('The grouping variable needs at least two levels'))


            ###########################
            ### Creates the data frame
            ###########################
            if (is.numeric(self$data[[eventVarName]])) {
                #---- If censor indicator is continuous
                status <- self$data[[eventVarName]]
                # follow survival::Surv(): 0/1 = censored/event, 1/2 = censored/event
                if ( ! all(status %in% c(0, 1, NA))) {
                    if ( ! all(status %in% c(1, 2, NA)))
                        jmvcore::reject(.('A numeric event variable must be coded 0 (censored) and 1 (event), or 1 (censored) and 2 (event). Otherwise, use a factor and choose the event level.'))
                    status <- status - 1
                }
                tmpDat <- data.frame(times = jmvcore::toNumeric(self$data[[elapsedVarName]]),
                                     status = status)
            } else {
                #---- If censor indicator is categorical
                if (is.null(self$options$eventLevel))
                    jmvcore::reject(.('Please choose the event level'))
                tmpDat <- data.frame(times = jmvcore::toNumeric(self$data[[elapsedVarName]]),
                                     status = ifelse(self$data[[eventVarName]] == self$options$eventLevel, 1, 0))
            }

            ##########################
            ### Add grouping variable
            ##########################
            group <- 1
            if ( ! is.null(self$options$groups))
                group <- as.factor(self$data[[self$options$groups]])
            if ( is.null(self$options$groups))
                group <- 1
            tmpDat$group <- group

            ############################
            ### use only complete cases
            ############################
            tmpDat <- na.omit(tmpDat)

            if (any(tmpDat$times < 0))
                jmvcore::reject(.('Time elapsed cannot be negative'))


            ############################
            ### Check if the data.frame
            ### contains data at all
            ############################
            if (nrow(tmpDat) == 0)
                jmvcore::reject(.('Data contains no (complete) rows'))

            ############################
            ### Drop levels of the
            ### grouping variable with
            ### no complete cases
            ############################
            if ( ! is.null(self$options$groups)) {
                tmpDat$group <- droplevels(tmpDat$group)
                if (nlevels(tmpDat$group) < 2)
                    jmvcore::reject(.('The grouping variable needs at least two levels with complete cases'))
            }

            tmpDat
        },
        .run = function() {

            eventVarName <- self$options$event
            elapsedVarName <- self$options$elapsed

            if (is.null(eventVarName) || is.null(elapsedVarName))
                return()

            tmpDat <- private$.buildData()

            empty <- setdiff(private$.groups(), levels(tmpDat$group))
            if ( ! is.null(self$options$groups) && length(empty) > 0) {
                note <- jmvcore::format(.('No complete cases for: {levels}'), levels=paste(empty, collapse=', '))
                self$results$summary$setNote('empty', note)
                self$results$medianestimates$setNote('empty', note)
            }

            private$.populateSummary(tmpDat)
            private$.populateMedians(tmpDat)

            if ( ! is.null(self$options$groups)) {
                private$.populateTests(tmpDat)
                private$.populatePairwise(tmpDat)
            }

        },
        .populateSummary = function(tmpDat) {

            ##################################
            ### Populate first table with the
            ### event summary
            ##################################
            table1 <- self$results$summary
            if (is.null(self$options$groups)) {
                fit <- survival::survfit(survival::Surv(times, status) ~ group, data=tmpDat)
                table1$setRow(rowNo=1, values=list(
                        n=fit$n,
                        censored = fit$n - sum(fit$n.event),
                        obs = sum(fit$n.event)))
            }
            else {
                fit <- survival::survdiff(survival::Surv(times, status) ~ group, data=tmpDat)
                for (i in seq_along(levels(tmpDat$group))) {
                    table1$setRow(rowKey=levels(tmpDat$group)[i], values=list(
                        n=fit$n[i],
                        censored = fit$n[i] - fit$obs[i],
                        obs=fit$obs[i],
                        exp=fit$exp[i]
                    ))
                }
            }

        },
        .populateMedians = function(tmpDat) {

            #######################################
            ### Populate the table with the median
            ### summaries
            #######################################
            s <- survival::survfit(survival::Surv(times, status) ~ group, data=tmpDat)
            if (is.null(self$options$groups))
                temp <- as.data.frame(rbind(summary(s)$table,summary(s)$table))[1,]
            else
                temp <- as.data.frame(summary(s)$table)
            table3 <- self$results$medianestimates
            for (i in seq_len(nrow(temp))) {
                values <- list(
                    median=temp$median[i],
                    cilb=temp$`0.95LCL`[i],
                    ciub=temp$`0.95UCL`[i])
                if (is.null(self$options$groups))
                    table3$setRow(rowNo=1, values=values)
                else
                    table3$setRow(rowKey=levels(tmpDat$group)[i], values=values)
            }

        },
        .populateTests = function(tmpDat) {

            #################################
            ### Populates the table with the
            ### tests results
            #################################
            tt <- self$results$tests
            if (tt$isFilled() || length(self$options$tests) == 0)
                return()

            df <- nlevels(tmpDat$group) - 1

            for (test in self$options$tests) {

                if (tt$isFilled(rowKey=1, col=paste0('chisqr[', test, ']')))
                    next()

                private$.checkpoint()

                chisqr <- switch(test,
                    logrank = survival::survdiff(survival::Surv(times, status) ~ group, data=tmpDat)$chisq,
                    petopeto = survival::survdiff(survival::Surv(times, status) ~ group, data=tmpDat, rho=1)$chisq,
                    gehan = private$.coinChisqr(tmpDat, "Gehan-Breslow"),
                    taroneware = private$.coinChisqr(tmpDat, "Tarone-Ware"))

                row <- list()
                row[[paste0('chisqr[', test, ']')]] <- chisqr
                row[[paste0('df[', test, ']')]] <- df
                row[[paste0('pvalue[', test, ']')]] <- pchisq(chisqr, df, lower.tail=FALSE)

                tt$setRow(rowKey=1, values=row)
            }

        },
        .populatePairwise = function(tmpDat) {

            #---------------------------------
            #-- Populates the table with the
            #-- tests results with multiple
            #-- comparisons.
            #-- Only if the number of levels
            #-- is >2
            #---------------------------------
            tt <- self$results$testspw
            groups <- levels(tmpDat$group)
            ngroups <- length(groups)

            if (ngroups <= 2) {
                tt$setNote("note0", jmvcore::format(.('Grouping variable, {var}, must have 3 or more levels to perform this analysis.'), var=self$options$groups))
                return()
            }
            if ( ! self$options$testspw)
                return()

            if (length(self$options$tests) == 0) {
                tt$setNote('notests', .('Select at least one test to perform pairwise comparisons.'))
                return()
            }

            if (tt$isFilled())
                return()

            nComparisons <- ngroups * (ngroups - 1) / 2

            groupsData <- list()
            for (group in groups) {
                ss <- (tmpDat$group == group)
                groupsData[[group]] <- subset(tmpDat, subset=ss, select=c('times', 'status'))
            }

            for (pair in tt$rowKeys) {

                if ( ! all(pair %in% groups))
                    next()

                x <- groupsData[[pair[1]]]
                y <- groupsData[[pair[2]]]

                for (test in self$options$tests) {

                    if (tt$isFilled(rowKey=pair, col=paste0('nupw[', test, ']')))
                        next()

                    private$.checkpoint()

                    testArg <- test
                    if (test == 'taroneware') {
                        testArg <- 'tarone-ware'
                    } else if (test == 'petopeto') {
                        testArg <- 'peto-peto'
                    }

                    result <- EnvStats::twoSampleLinearRankTestCensored(
                        test = testArg,
                        x = x$times,
                        x.censored = ! x$status,
                        y = y$times,
                        y.censored = ! y$status,
                        censoring.side = 'right')

                    row <- list()

                    row[[paste0('zpw[', test, ']')]] <- result$statistic['z']
                    row[[paste0('nupw[', test, ']')]] <- result$statistic['nu']
                    row[[paste0('nusepw[', test, ']')]] <- sqrt(result$statistic['var.nu'])
                    row[[paste0('ppw[', test, ']')]] <- min(result$p.value * nComparisons, 1)

                    tt$setRow(rowKey=pair, values=row)
                }
            }

            tt$setNote("note", .("P-values are Bonferroni-corrected."))

        },
        .coinChisqr = function(data, type) {
            fit <- coin::logrank_test(survival::Surv(times, status) ~ group, data=data, type=type)
            stat <- coin::statistic(fit)
            # two-sample tests return a standardised z, k-sample tests a chi-squared
            if (inherits(fit, "ScalarIndependenceTest"))
                stat <- stat^2
            unname(stat[1])
        },
        .plot=function(image, theme, ggtheme, ...) {

            eventVarName <- self$options$event
            elapsedVarName <- self$options$elapsed

            if (is.null(eventVarName) || is.null(elapsedVarName))
                return()

            tmpDat <- private$.buildData()

            s <- survival::survfit(survival::Surv(times, status) ~ group, data=tmpDat)
            if (!is.null(names(s$strata))) names(s$strata) <- gsub("group=", "", names(s$strata))

            ci <- self$options$ci
            cens <- self$options$cens
            risk <- self$options$risk
            legend <- "top"
            legendtitle <- self$options$groups
            if (is.null(self$options$groups)) {
                legend <-"none"
                legendtitle <- NULL
                }
            xlab <- .("Time")
            if (!self$options$units=="") xlab <- jmvcore::format(.("Time ({units})"), units=self$options$units)
            survmedianline <-"none"
            if (self$options$median) survmedianline <- "hv"


            if (identical(image, self$results$sc)) {
                func <- NULL
                ylab <- .("Survival Probability")
            } else {
                if (!self$options$logscale) {
                    func <- "cumhaz"
                    ylab <- .("Cumulative Hazard")
                    ci <- FALSE
                    cens <- FALSE
                    risk <- FALSE
                    survmedianline <-"none"
                } else {
                    func <- function(x) log(-log(x))
                    ylab <- .("Cumulative Hazard (log scale)")
                    ci <- FALSE
                    cens <- FALSE
                    risk <- FALSE
                    survmedianline <-"none"
                    }
                }

            p <- survminer::ggsurvplot(s,
                                       data = tmpDat,
                                       xlab=xlab,
                                       ylab = ylab,
                                       conf.int = ci,
                                       legend = legend,
                                       legend.title=legendtitle,
                                       censor = cens,
                                       fun = func,
                                       censor.shape = "|",
                                       censor.size = 6,
                                       risk.table = risk,
                                       surv.median.line = survmedianline,
                                       ggtheme = ggtheme,
                                       tables.height = length(unique(tmpDat$group))*0.02+0.21
                                       )
            return(p)
        })
)
