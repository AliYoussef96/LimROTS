bootstrap.loading <- function() {
    D <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    S <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pD <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))
    pS <- matrix(nrow = nrow(as.matrix(data)), ncol = nrow(samples))

    if (is.null(cluster)) {
        cluster <- makeCluster(2)
        registerDoParallel(cluster)
        message("No cluster found; only two cores will be used!")
    } else {
        registerDoParallel(cluster)
    }

    clusterSetRNGStream(cluster, iseed = seed.cl)
    clusterExport(
        cluster,
        varlist = c(
            "pb",
            "samples",
            "pSamples",
            "D",
            "data",
            "S",
            "pD",
            "pS",
            "time",
            "formula.str",
            "group.name",
            "groups",
            "event",
            "meta.info",
            "a1",
            "a2",
            "trend",
            "robust",
            "n.ROTS",
            "survival"
        ),
        envir = environment()
    )


    results_list <- foreach(
        i = seq_len(nrow(samples)),
        .combine = "c",
        .packages = c("utils", "stringr", "stats" , "limma"),
        .export = c(
            "testStatSurvivalOptimized" ,
            "testStatistic_with_covariates" ,
            "testStatOptimized",
            "testStatistic_with_covariates_permutating"
        )
    ) %dorng% {
        samples.R <- split(samples[i, ], groups)
        d_result <- s_result <- pd_result <- ps_result <- NULL
        if (is.null(a1) | is.null(a2)) {
            if (survival == TRUE) {
                fit <- testStatSurvivalOptimized(lapply(samples.R, function(x)
                    data[, x]),
                    groups,
                    event)
            } else if (n.ROTS == FALSE) {
                fit <- testStatistic_with_covariates(
                    x = lapply(samples.R, function(x)
                        data[, x]),
                    group.name = group.name,
                    meta.info = meta.info,
                    formula.str = formula.str,
                    trend =
                        trend,
                    robust = robust
                )
            } else {
                fit <- testStatOptimized(paired, lapply(samples.R, function(x)
                    data[, x]))
            }
            d_result <- fit$d
            s_result <- fit$s
            df1 <- data.frame(d_result = d_result, s_result = s_result)
        }
        if (survival == TRUE) {
            pSamples.R <- split(pSamples[i, ], groups)
            pFit <- testStatSurvivalOptimized(lapply(pSamples.R, function(x)
                data[, x]), groups, event)
        } else if (n.ROTS == FALSE) {
            pFit <- testStatistic_with_covariates_permutating(
                x = lapply(split(seq_len(
                    length(groups)
                ), groups), function(x)
                    data[, x]),
                group.name = group.name,
                meta.info = meta.info,
                formula.str = formula.str,
                trend =
                    trend,
                robust = robust ,
                permutating.group = permutating.group
            )
        } else {
            pSamples.R <- split(pSamples[i, ], groups)

            pFit <- testStatOptimized(paired, lapply(pSamples.R, function(x)
                data[, x]))
        }
        pd_result <- pFit$d
        ps_result <- pFit$s
        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
        list(ds = df1, pdps = df2)
    }

    stopCluster(cluster)


    names(results_list) <- paste0(names(results_list), seq(1, length(names(results_list))))
    j <- 0
    q <- 0
    for (i in seq_along(results_list)) {
        if (grepl("ds", names(results_list)[i], fixed = TRUE)) {
            j <- j + 1
            D[, j] <- results_list[[names(results_list)[i]]]$d_result
            S[, j] <- results_list[[names(results_list)[i]]]$s_result
        } else {
            q <- q + 1
            pD[, q] <- results_list[[names(results_list)[i]]]$pd_result
            pS[, q] <- results_list[[names(results_list)[i]]]$ps_result
        }

    }

    return(list(D = D, S = S, pD = pD , pS))
}
