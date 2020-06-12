

#' @importFrom R6 R6Class
#' @importFrom jmvcore toNumeric
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
            if (length(groups) >= 2) {
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
            if (length(groups) >= 2) {
                comparisons <- combn(groups, 2)
                for (i in seq_len(ncol(comparisons))) {
                    key = comparisons[,i]
                    testspw$addRow(rowKey=key)
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
        .run = function() {
            
            eventVarName <- self$options$event
            elapsedVarName <- self$options$elapsed

            if (is.null(eventVarName) || is.null(elapsedVarName))
                return()
            
            if (nrow(self$data) == 0)
                stop('Data contains no (complete) rows')
            
            
            ###########################
            ### Creates the data frame
            ###########################
            #---- If censor indicator is continuous
            if (class(self$data[[eventVarName]])=="numeric" | class(self$data[[eventVarName]])=="integer") {
                tmpDat <- data.frame(times = as.numeric(self$data[[elapsedVarName]]),
                                     status = self$data[[eventVarName]])
            }
            #---- If censor indicator is categorical
            if ( ! (class(self$data[[eventVarName]])=="numeric" | class(self$data[[eventVarName]])=="integer") ) {
                tmpDat <- data.frame(times = as.numeric(self$data[[elapsedVarName]]),
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
            
            
            ############################
            ### Check if the data.frame
            ### contains data at all
            ############################
            if (nrow(tmpDat) == 0)
                stop('Data contains no (complete) rows')
            

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
            if (!is.null(self$options$groups)) {
                fit <- survival::survdiff(survival::Surv(times, status) ~ group, data=tmpDat)
                for (i in 1:length(unique(tmpDat$group)) ) {
                    table1$setRow(rowNo=i, values=list(
                        n=fit$n[i],
                        censored = fit$n[i] - fit$obs[i],
                        obs=fit$obs[i],
                        exp=fit$exp[i]
                    ))
                    }
                }            


            #######################################
            ### Populate the table with the median
            ### summaries
            #######################################
            s <- survival::survfit(survival::Surv(times, status) ~ group, data=tmpDat)
            if (is.null(self$options$groups)) {
                temp <- as.data.frame(rbind(summary(s)$table,summary(s)$table))[1,]
            } 
            if (!is.null(self$options$groups)) temp <- as.data.frame(summary(s)$table)
            table3 <- self$results$medianestimates
            for (i in 1:length(unique(tmpDat$group)) ) {
                table3$setRow(rowNo=i, values=list(
                    median=temp$median[i],
                    cilb=temp$`0.95LCL`[i],
                    ciub=temp$`0.95UCL`[i]
                ))
            }
            
            
            #################################
            #### Set the state for the plots
            #################################
            self$results$sc$setState(tmpDat)
            self$results$chf$setState(tmpDat)
            
            
            #---------------------------------
            #-- Populates the table with the  
            #-- tests results
            #-- Only if there is a grouping
            #-- variable
            #---------------------------------
            if (!is.null(self$options$groups)) {
                #######################################
                ### Calculates data frame with results
                ### of the tests
                #######################################
                df <- length(unique(tmpDat$group)) - 1
                # -------   Logrank test
                fit <- survival::survdiff(survival::Surv(times, status) ~ group, data=tmpDat)
                dat1 <- data.frame(var = "Log-Rank",chisqr = fit$chisq)
                # -------   Peto-Peto
                fit <- survival::survdiff(survival::Surv(times, status) ~ group, data=tmpDat,rho=1)
                dat2 <- data.frame(var = "Peto-Peto",chisqr = fit$chisq)
                # -------  Gehan
                fit <- coin::logrank_test(survival::Surv(times, status) ~ group, data=tmpDat,type="Gehan-Breslow")
                dat3 <- data.frame(var = "Gehan",chisqr = coin::statistic(fit)^2)
                # -------  Tarone-Ware
                fit <- coin::logrank_test(survival::Surv(times, status) ~ group, data=tmpDat,type="Tarone-Ware")
                dat4 <- data.frame(var = "Tarone-Ware",chisqr = coin::statistic(fit)^2)
                # ------- bind the methods
                DF <- rbind(dat1,dat2,dat3,dat4)
                # ------- add the df and the p-values
                DF$df <- df
                DF$pvalue <- pchisq(DF$chisqr, df, lower.tail = FALSE)
                
                
                #################################
                ### Populates the table with the  
                ### tests results
                #################################
                tt <- self$results$tests
                if (tt$isNotFilled() && length(self$options$tests) > 0) {
                    
                    for (i in seq_along(self$options$tests)) {
                        
                        test <- self$options$tests[i]
                        
                        if (tt$isFilled(rowKey=1, col=paste0('chisqr[', test, ']')))
                            next()
                        
                        private$.checkpoint()
                        row <- list()
                        row[[paste0('chisqr[', test, ']')]] <- DF$chisqr[i]
                        row[[paste0('df[', test, ']')]] <- DF$df[i]
                        row[[paste0('pvalue[', test, ']')]] <- DF$pvalue[i]
                        
                        tt$setRow(rowKey=1, values=row)
                        
                    }

                }

            }
            
            #---------------------------------
            #-- Populates the table with the  
            #-- tests results with multiple
            #-- comparisons.
            #-- Only if there is a grouping
            #-- variable and the number of
            #-- levels is >2
            #---------------------------------
            if (!is.null(self$options$groups) ) {
                tt <- self$results$testspw
                if (length(unique(tmpDat$group))<=2) {
                    tt$setNote("note0",paste("Grouping variable, ",self$options$groups,", must have 3 or more levels to perform this analysis.",sep="") )
                    return()
                }
                if (tt$isNotFilled() && self$options$testspw) {
                    
                    groups <- private$.groups()
                    ngroups <- length(groups)
                    
                    if (length(groups) > 2) {
                        groupsData <- list()
                        for (group in groups) {
                            ss <- (tmpDat$group == group)
                            groupsData[[group]] <- subset(tmpDat, subset=ss, select=c('times', 'status'))
                        }
                    }
                    
                    for (pair in tt$rowKeys) {
                        
                        x <- groupsData[[pair[1]]]
                        y <- groupsData[[pair[2]]]
                        
                        for (i in seq_along(self$options$tests)) {
                            test <- self$options$tests[i]
                            
                            if (tt$isFilled(rowKey=pair, col=paste0('nupw[', test, ']')))
                                next()
                            
                            private$.checkpoint()
                            
                            result <- EnvStats::twoSampleLinearRankTestCensored(
                                test = test,
                                x = x$times,
                                x.censored = ! x$status,
                                y = y$times,
                                y.censored = ! y$status,
                                censoring.side = 'right')
                            
                            row <- list()
                            
                            row[[paste0('zpw[', test, ']')]] <- result$statistic['z']
                            row[[paste0('nupw[', test, ']')]] <- result$statistic['nu']
                            row[[paste0('nusepw[', test, ']')]] <- sqrt(result$statistic['var.nu'])
                            row[[paste0('ppw[', test, ']')]] <- result$p.value * min(ngroups*(ngroups - 1)/2,1)
                            
                            tt$setRow(rowKey=pair, values=row)
                        }
                       
                        tt$setNote("note","P-values are Bonferroni-corrected.") 
                    }
                }
            }

        },
        .plot=function(image, theme, ggtheme, ...) {
            
            eventVarName <- self$options$event
            elapsedVarName <- self$options$elapsed
            
            if (is.null(eventVarName) || is.null(elapsedVarName))
                return()
            
            tmpDat <- image$state 
            
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
            xlab <- "Time"
            if (!self$options$units=="") xlab <- paste("Time (",self$options$units,")",sep="")
            survmedianline <-"none"
            if (self$options$median) survmedianline <- "hv"
            
            
            if (identical(image, self$results$sc)) {
                func <- NULL
                ylab <- "Survival Probabilities"
            } else {
                if (!self$options$logscale) {
                    func <- "cumhaz"
                    ylab <- "Cumulative Hazards"
                    ci <- FALSE
                    cens <- FALSE
                    risk <- FALSE
                    survmedianline <-"none"
                } else {
                    func <- function(x) log(-log(x))
                    ylab <- "Cumulative Hazards (log scale)"
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
