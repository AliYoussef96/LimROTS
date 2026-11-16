# ==============================================================================
# limrots_peptide.jl
#
# Julia helper for LimROTS_peptide(backend = "julia").
#
# WHAT THIS FILE DOES
#   For ONE replicate (the original data, one bootstrap sample, or one
#   permutation) it fits one linear mixed model per protein:
#
#       y ~ 1 + x1 + x2 + ... + (1 | peptide_id) + (1 | sample_id)     (REML)
#
#   and returns, for every protein, the group contrast estimate and its
#   standard error:
#
#       estimate = L' * beta_hat              (L = contrast vector from R)
#       se       = sqrt( L' * Vcov(beta) * L )
#
#   This is exactly what emmeans(..., lmer.df = "asymptotic") + contrast() does
#   in the R code, without the (slow) emmeans machinery.
#
# WHAT THIS FILE DOES NOT DO
#   * No random numbers: bootstrap / permutation sampling is done in R.
#   * No formula parsing: R builds the design matrix X (S x k) and the
#     contrast vector L (length k + 1, intercept first), so the contrast
#     direction and the factor coding are decided in R, exactly as before.
#   * No absolute values: R applies abs() for bootstrap / permutation.
# ==============================================================================

using MixedModels      # fits the mixed models
using DataFrames       # small per-protein tables
using LinearAlgebra    # BLAS threading control

# Every model is tiny, so it is faster to run MANY models at the same time
# (one per Julia thread) than to let each model use several BLAS threads.
BLAS.set_num_threads(1)


# Helper: JuliaConnectoR turns R vectors of length 1 into plain scalars.
# This makes sure we always work with a Vector.
_as_vector(x) = x isa AbstractArray ? vec(x) : [x]


"""
    limrots_fit_proteins(Y, protein_index, X, L, protein_names)

Arguments (all come from R):
  * Y             : peptides x samples intensity matrix (no missing values)
  * protein_index : for each row of Y, the number (1..n_proteins) of its protein
  * X             : samples x k matrix of the fixed-effect columns,
                    WITHOUT the intercept column (rows = columns of Y)
  * L             : contrast vector of length k + 1 (intercept first)
  * protein_names : names of the proteins, in the order of protein_index 1..n

Returns a NamedTuple with
  * protein  : the protein names, in the order of the results
  * estimate : signed contrast estimate per protein
  * se       : standard error per protein
"""
function limrots_fit_proteins(Y, protein_index, X, L, protein_names)

    # ---- 1. Convert whatever R sent into plain Julia types -----------------
    Y             = Matrix{Float64}(Y)
    X             = Matrix{Float64}(X)
    protein_index = Vector{Int}(_as_vector(protein_index))
    L             = Vector{Float64}(_as_vector(L))
    protein_names = Vector{String}(_as_vector(protein_names))

    n_peptides, n_samples = size(Y)
    n_proteins = length(protein_names)
    k = size(X, 2)                         # number of non-intercept columns

    # ---- 2. Basic safety checks (R checks too, this is a second guard) -----
    length(protein_index) == n_peptides ||
        throw(ArgumentError("protein_index must have one entry per row of Y"))
    size(X, 1) == n_samples ||
        throw(ArgumentError("X must have one row per column of Y"))
    length(L) == k + 1 ||
        throw(ArgumentError("L must have length ncol(X) + 1 (intercept first)"))
    maximum(protein_index) == n_proteins ||
        throw(ArgumentError("protein_index does not match protein_names"))

    # ---- 3. Build the model formula once (it is the same for all proteins) -
    # The fixed columns were renamed x1, x2, ... by R.
    x_names = ["x$j" for j in 1:k]
    rhs     = join(vcat("1", x_names), " + ")
    formula = eval(Meta.parse(
        "@formula(y ~ $rhs + (1 | peptide_id) + (1 | sample_id))"))

    # peptide_id / sample_id are only labels, so tell MixedModels to treat
    # them as grouping factors (no contrasts needed).
    grouping = Dict(:peptide_id => Grouping(), :sample_id => Grouping())

    # ---- 4. Which rows of Y belong to which protein ------------------------
    rows_of_protein = [findall(==(p), protein_index) for p in 1:n_proteins]

    # ---- 5. Result containers (one slot per protein, so threads never clash)
    estimate = fill(NaN, n_proteins)
    se       = fill(NaN, n_proteins)
    errors   = fill("", n_proteins)        # error message per protein, if any

    # ---- 6. Fit all proteins, one protein per thread at a time -------------
    Threads.@threads for p in 1:n_proteins
        try
            rows = rows_of_protein[p]      # the peptides of this protein
            m    = length(rows)            # number of peptides

            # Long format table: one row per (peptide, sample) pair.
            # vec(Y[rows, :]) stacks the samples one after another, so
            #   peptide_id = 1..m, 1..m, ...   (repeat)
            #   sample_id  = 1,1,..,1, 2,2,..  (inner repeat)
            df = DataFrame(
                y          = vec(Y[rows, :]),
                peptide_id = repeat(1:m, n_samples),
                sample_id  = repeat(1:n_samples, inner = m),
            )
            # Sample-level covariates are repeated for every peptide.
            for j in 1:k
                df[!, Symbol(x_names[j])] = repeat(X[:, j], inner = m)
            end

            # REML fit with the same random structure as the R model.
            model = fit(MixedModel, formula, df;
                        REML = true, progress = false, contrasts = grouping)

            beta = coef(model)             # fixed effects (intercept first)
            V    = vcov(model)             # covariance matrix of the fixed effects

            estimate[p] = dot(L, beta)
            se[p]       = sqrt(dot(L, V * L))
        catch err
            errors[p] = sprint(showerror, err)   # remember why it failed
        end
    end

    # ---- 7. Abort (like the R code does) if any fit failed -----------------
    failed = [p for p in 1:n_proteins
              if !isempty(errors[p]) || !isfinite(estimate[p]) || !isfinite(se[p])]
    if !isempty(failed)
        shown = join(protein_names[failed[1:min(5, end)]], ", ")
        msg   = isempty(errors[failed[1]]) ? "non-finite estimate or SE" :
                                             errors[failed[1]]
        error("Julia fit failed for $(length(failed)) protein(s) " *
              "(first: $shown). First error: $msg")
    end

    # ---- 8. Return in the same order as protein_names ----------------------
    return (protein = protein_names, estimate = estimate, se = se)
end
