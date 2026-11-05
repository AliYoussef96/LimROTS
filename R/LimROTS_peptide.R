LimROTS_peptide <- function(    
    assay.type = NULL,
    niter = 1000,
    K = NULL,
    meta.info,
    peptide_metadata,
    BPPARAM  = NULL,
    group,
    formula.str,
    pairwise_contrasts,
    strata){

    ### Perform sanity checks on the input parameters and data.
    SanityChecK.list <- SanityChecK(
        x,
        assay.type = assay.type,
        niter = niter,
        K = K,
        meta.info = meta.info,
        group = group,
        model = "lmer", 
        formula.str = formula.str,
        pairwise_contrasts = pairwise_contrasts,
        correlation_block = NULL,
        strata = strata,
        survival = FALSE,
        competing_risks = FALSE
    )


    meta.info <- SanityChecK.list[["meta.info"]]
    data <- SanityChecK.list[["data"]]
    group <- SanityChecK.list[["group"]]
    K <- SanityChecK.list[["K"]]

    ### Extra Perform sanity checks on the peptide_metadata

    formula_variables <- as.formula(formula.str)
    formula_variables <- all.vars(formula_variables)

    if (length(formula_variables) < 2) {
        stop("The formula must contain a response and a group variable.")
    }

    if (formula_variables[2] != group) {
        stop(
            paste0(
                "The first variable after '~' must be the group variable: ",
                group
            )
        )
    }

    formula.str <- paste0(
        "y ~" ,formula.str,
        " + peptide_id + (1 | sample_id)"
    )

    model_formula <- as.formula(formula.str)
    
    message("Using the following formula for the model: ", formula.str)

    #============================================================
    # Check the intensity matrix
    #============================================================

    if (!is.matrix(data)) {
        stop("data must be a matrix.")
    }

    if (!is.numeric(data)) {
        stop("data must contain numeric intensities.")
    }

    if (is.null(rownames(data))) {
        stop("data must have peptide IDs as row names.")
    }

    if (is.null(colnames(data))) {
        stop("data must have sample IDs as column names.")
    }

    if (anyDuplicated(rownames(data))) {
        stop("Peptide IDs in data must be unique.")
    }

    if (anyDuplicated(colnames(data))) {
        stop("Sample IDs in data must be unique.")
    }
    #============================================================
    # Check peptide metadata
    #============================================================

    if (
        !all(
            c("peptide_id", "protein_id") %in%
            colnames(peptide_metadata)
        )
    ) {
        stop("peptide_metadata must contain peptide_id and protein_id.")
    }

    if (anyNA(peptide_metadata[["peptide_id"]])) {
        stop("peptide_metadata contains missing peptide IDs.")
    }

    if (anyNA(peptide_metadata[["protein_id"]])) {
        stop("peptide_metadata contains missing protein IDs.")
    }

    if (anyDuplicated(peptide_metadata[["peptide_id"]])) {
        stop("peptide IDs in peptide_metadata must be unique.")
    }

    if (
        !setequal(
            rownames(data),
            peptide_metadata[["peptide_id"]]
        )
    ) {
        stop("Peptide IDs do not align between data and peptide_metadata.")
    }


    #============================================================
    # Check sample metadata
    #============================================================
    if (anyNA(meta.info)) {
        stop("meta.info contains missing values.")
    }

    if (anyDuplicated(meta.info[["sample_id"]])) {
        stop("sample IDs in meta.info must be unique.")
    }

    if (
        !identical()(
            colnames(data),
            meta.info[["sample_id"]]
        )
    ) {
        stop("Sample IDs do not align between data and meta.info.")
    }


    #============================================================
    # Reorder metadata to exactly match the matrix
    #============================================================

    peptide_metadata <- peptide_metadata[
        match(
            rownames(data),
            peptide_metadata[["peptide_id"]]
        ),
        ,
        drop = FALSE
    ]

    if(!identical(rownames(data), peptide_metadata[["peptide_id"]])) {
        stop("Peptide IDs do not align between data and peptide_metadata after
         reordering.")
    }


    ### bootstrap and permutation samples
    list_boot_perm <- .generate_bootstrap_permutations(
        niter_2 = 2 * niter,
        meta.info = meta.info,
        group = group,
        correlation_block = NULL,
        strata = strata
    )


    ### parallelized bootstrapping and permutation fit
    results_list <- parallel_peptide(
        BPPARAM  = BPPARAM ,
        samples = list_boot_perm[["boot_Samples"]], 
        data = data,
        model_formula = model_formula,
        group = group,
        meta.info = meta.info,
        pSamples = list_boot_perm[["perm_Samples"]],
        pairwise_contrasts = pairwise_contrasts
    )


    }