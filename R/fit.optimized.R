fit.optimized  <- function(survival = survival, n.ROTS = n.ROTS){
    if (survival == TRUE) {
        fit <- testStatSurvivalOptimized(lapply(split(seq_len(
            length(groups)
        ), groups), function(x)
            data[, x]), groups, event)
    } else if (n.ROTS == FALSE) {
        fit <- testStatistic_with_covariates_Fit(
            x = lapply(split(seq_len(
                length(groups)
            ), groups), function(x)
                data[, x]),
            group.name = group.name,
            meta.info = meta.info,
            formula.str = formula.str,
            trend = trend,
            robust = robust
        )
    } else {
        fit <- testStatOptimized(paired, lapply(split(seq_len(
            length(groups)
        ), groups), function(x)
            data[, x]))
    }


    return(fit)
}
