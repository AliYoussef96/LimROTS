Boot_parallel <- function(cluster, seed.cl , samples, pSamples, D , data,
                          S, pD, pS, formula.str, group.name, groups,
                          meta.info, a1, a2, trend, robust,
                          permutating.group) {
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
        varlist = c("samples", "pSamples", "D", "data", "S", "pD", "pS",
                        "formula.str", "group.name", "groups", "meta.info",
                            "a1", "a2", "trend", "robust" , "permutating.group"),
        envir = environment()
                    )
    results_list <- foreach(
        i = seq_len(nrow(samples)),
        .combine = "c",
        .packages = c("utils", "stringr", "stats", "limma"),
        .export = c("testStatistic_with_covariates",
                        "testStatistic_with_covariates_permutating")
    ) %dorng% {
        samples.R <- split(samples[i, ], groups)
        # Initialize placeholders for results
        d_result <- s_result <- pd_result <- ps_result <- NULL
        # Compute D and S if conditions are met
        if (is.null(a1) | is.null(a2)) {
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
        }
        d_result <- fit$d
        s_result <- fit$s
        df1 <- data.frame(d_result = d_result, s_result = s_result)

        # Compute pD and pS
        pFit <- testStatistic_with_covariates_permutating(
            x = lapply(split(seq_len(
                length(groups)
            ), groups), function(x)
                data[, x]),
            group.name = group.name,
            meta.info = meta.info,
            formula.str = formula.str,
            trend = trend,
            robust = robust,
            permutating.group = permutating.group
        )
        pd_result <- pFit$d
        ps_result <- pFit$s
        df2 <- data.frame(pd_result = pd_result, ps_result = ps_result)
        # Return results for this iteration as a data frame
        list(ds = df1, pdps = df2)
    }
    stopCluster(cluster)
    return(results_list)
}
