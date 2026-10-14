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
