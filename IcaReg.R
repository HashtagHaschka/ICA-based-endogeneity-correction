## -----------------------------------------------------------------------------
##  ICA-based endogeneity correction in R
##  https://github.com/HashtagHaschka/ICA-based-endogeneity-correction
##
##  Copyright (C) 2026 Rouven E. Haschka
##  ORCID: https://orcid.org/0000-0002-2916-9745
##
##  If this code contributes to work you publish, please cite the software
##
##    Haschka, R. E. (2026). ICA-based endogeneity correction in R.
##    https://github.com/HashtagHaschka/ICA-based-endogeneity-correction
##
##  and the paper this estimator implements
##
##    Haschka, R. E. and F. Dost (2025). ICA at the cocktail party: Casting
##      instrument-free omitted variable bias correction as a blind source
##      separation problem. Available at SSRN: https://ssrn.com/abstract=5361801
##
##  This program is free software: you can redistribute it and/or modify it
##  under the terms of the GNU General Public License as published by the Free
##  Software Foundation, either version 3 of the License, or (at your option)
##  any later version.
##
##  This program is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
##  more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program.  If not, see <https://www.gnu.org/licenses/>.
## -----------------------------------------------------------------------------


## =============================================================================
##  IcaReg -- instrument-free correction of omitted variable bias
## =============================================================================
##
##  METHOD
##
##  The response and the endogenous regressors are decomposed by Independent
##  Component Analysis into statistically independent sources. Under the
##  maintained assumptions -- the omitted confounder is normally distributed and
##  enters the endogenous regressor additively, and the regressor also carries
##  non-normal exogenous variation -- exactly one of the recovered sources is
##  Gaussian. That source estimates the confounder and is used as a control
##  function in the structural equation, which restores consistency without
##  external instruments.
##
##  Exogenous regressors are partialled out by OLS before the decomposition.
##  Standard errors are delete-one jackknife standard errors: every replicate
##  re-runs the complete pipeline, so the additional sampling variability from
##  the first-stage projection and from the estimated control function is
##  accounted for automatically.
##
##
##  USAGE
##
##    IcaReg(formula, data, method = "jade", CF = FALSE,
##           select = c("kurtosis", "ks"), kurt_tol = 1,
##           parallel = FALSE, ncores = NULL)
##
##
##  ARGUMENTS
##
##  formula   depvar ~ endog_var1 + endog_var2 + ... | exog_var1 + exog_var2 + ...
##            - depvar is the dependent variable, endog_var1 etc. are the
##              continuous endogenous variables, exog_var1 etc. are exogenous
##            - the exogenous block may be omitted entirely
##            - transformations are allowed on both sides, e.g.
##              log(sales) ~ log(price) | promotion
##            - dummy variables via as.factor(exog_var1)
##            - accepts -1 to remove the intercept; place it in the FIRST part
##            - a variable must not appear in both blocks
##
##  data      a data.frame
##
##  method    ICA algorithm, passed to ica::ica(). One of "jade" (default),
##            "fast", "imax". "jade" separates sources using fourth-order
##            cumulants and matches the identification argument of the paper.
##
##  CF        FALSE (default) reports the structural coefficients only.
##            TRUE additionally reports the coefficient on the control function,
##            which is a Hausman-type test statistic for endogeneity: a
##            coefficient indistinguishable from zero indicates no detectable
##            endogeneity. All other coefficients are identical either way.
##
##  select    rule for identifying the Gaussian source.
##            "kurtosis" (default) selects the source with the smallest absolute
##            excess kurtosis. This is the quantity identification turns on and
##            the quantity JADE optimises.
##            "ks" selects the smallest Kolmogorov-Smirnov distance to normality.
##            The two rules agree whenever identification is clear-cut; comparing
##            them is a useful robustness check.
##
##  kurt_tol  threshold for the diagnostic warning on the selected source.
##            A warning is issued if its absolute excess kurtosis exceeds this
##            value, i.e. if no source looks Gaussian at all. Default 1.
##
##  parallel  compute the jackknife replicates in parallel. Recommended for
##            n greater than a few thousand.
##
##  ncores    number of cores; defaults to detectCores() - 1.
##
##
##  VALUE
##
##  A list with two elements:
##    [[1]]  matrix of coefficients and delete-one jackknife standard errors
##    [[2]]  the estimated control function, one value per observation
##
##
##  DIAGNOSTICS
##
##  The excess kurtosis of every recovered source is reported. Identification
##  requires exactly one source close to zero while the others are clearly away
##  from it. Under normality the sampling standard deviation of the sample excess
##  kurtosis is about sqrt(24 / n), which gives a scale for judging the readout.
##  A reading such as
##
##      Excess kurtosis of the ICA sources: +2.103  +53.669  -0.101
##
##  at n = 250 is clear-cut: one source sits within sampling noise of zero and
##  the nearest competitor is several standard deviations away. Two sources close
##  to zero mean the Gaussian source is not pinned down and the estimates should
##  not be trusted.
##
##  Note that the reported value for the selected source is a minimum over all
##  sources and is therefore biased towards zero, so the warning is conservative.
##
##
##  REQUIREMENTS AND LIMITATIONS
##
##  - endogenous regressors must be continuous; ties are reported because
##    discreteness destabilises the higher-order moments the method relies on
##  - identification is weak when the exogenous variation in the endogenous
##    regressor is close to normal, which shows up as large standard errors
##  - the jackknife assumes independent observations. With clustered or panel
##    data the standard errors will be too small
##
## =============================================================================


## ---- packages ---------------------------------------------------------------
pacman::p_load(
  nlme,
  ica,
  pbapply,
  parallel
)


## ---- internal helpers -------------------------------------------------------

## Deparse an expression into a single string.
.IcaReg_dep1 <- function(x) paste(deparse(x), collapse = " ")

## Rebuild a language object from its deparsed form.
## Language objects must not travel through do.call(): with quote = FALSE, the
## default, do.call evaluates every element of the argument list in the calling
## frame, so a bare symbol would be looked up as a variable. The argument list
## therefore carries strings only, and expressions are reconstructed here.
.IcaReg_lang <- function(txt) parse(text = txt, keep.source = FALSE)[[1]]

## Excess kurtosis, i.e. the standardised fourth cumulant. Zero under normality.
.IcaReg_kurt <- function(x) {
  x  <- x - mean(x)
  m2 <- mean(x * x)
  if (!is.finite(m2) || m2 <= 0) return(NA_real_)
  mean(x^4) / (m2 * m2) - 3
}

## Selection rule "kurtosis": the Gaussian source is the one whose fourth
## cumulant vanishes in population.
.IcaReg_select_kurtosis <- function(S) {
  k <- apply(S, 2, .IcaReg_kurt)
  k[!is.finite(k)] <- Inf
  which.min(abs(k))
}

## Selection rule "ks": smallest Kolmogorov-Smirnov distance to normality. This
## is a different notion of non-Gaussianity from the one JADE exploits, since a
## source can be far from normal in distribution and still have zero excess
## kurtosis.
.IcaReg_select_ks <- function(S) {
  ks <- apply(S, 2, function(x) {
    suppressWarnings(
      stats::ks.test(x, "pnorm", mean = mean(x), sd = stats::sd(x))$statistic
    )
  })
  which.min(ks)
}


## Run the complete pipeline once: partial out the exogenous block, decompose,
## select the Gaussian source, fit the structural equation.
##
## select_fn(S) returns the index of the column of the source matrix to be used
## as the control function. In the full sample this is the rule chosen via
## `select`; in a jackknife replicate it is correlation-matching to the
## full-sample control function, because ICA identifies sources only up to
## permutation and sign and replicates must therefore be aligned.
.IcaReg_fit <- function(df, lhs_label, P_rhs, X_rhs, has_exog,
                        has_intercept, CF, method, select_fn,
                        env = parent.frame()) {

  n <- nrow(df)

  ## response
  lhs <- .IcaReg_lang(lhs_label)
  y   <- eval(lhs, envir = df, enclos = env)
  if (!is.numeric(y) || length(y) != n)
    stop("The left-hand side '", lhs_label,
         "' does not evaluate to a numeric vector of length nrow(data).")

  ## design matrices
  Pmm  <- stats::model.matrix(stats::as.formula(paste("~", P_rhs)), data = df)
  Pcol <- setdiff(colnames(Pmm), "(Intercept)")
  Pmat <- Pmm[, Pcol, drop = FALSE]

  if (has_exog) {
    Xmm  <- stats::model.matrix(stats::as.formula(paste("~", X_rhs)), data = df)
    Xcol <- setdiff(colnames(Xmm), "(Intercept)")
    Xmat <- Xmm[, Xcol, drop = FALSE]
  } else {
    Xcol <- character(0)
    Xmat <- Pmat[, 0L, drop = FALSE]
  }

  clash <- intersect(Pcol, Xcol)
  if (length(clash))
    stop("The following regressors appear in BOTH blocks of the formula: ",
         paste(clash, collapse = ", "),
         ". After partialling X out of P they would be identically zero.")

  ## model.matrix() drops rows whose terms evaluate to NA/NaN. Such rows are
  ## removed up front, so this guard should never fire; if it ever does, fail
  ## loudly rather than misalign y, P and X against each other.
  if (nrow(Pmat) != n || nrow(Xmat) != n)
    stop("Row-count mismatch after evaluating the model terms: ", n,
         " observations, ", nrow(Pmat), " rows in the endogenous block, ",
         nrow(Xmat), " in the exogenous block. A transformed term probably ",
         "evaluated to NA/NaN.")

  ## First stage. [1, X] is always partialled out, so the recovered sources are
  ## exactly orthogonal to the constant and to the exogenous block; whether the
  ## second stage carries an intercept is governed by has_intercept alone.
  W   <- cbind(`(Intercept)` = rep(1, n), Xmat)
  qrW <- qr(W)
  yr  <- as.numeric(qr.resid(qrW, y))
  Pr  <- qr.resid(qrW, Pmat)
  Pr  <- matrix(Pr, nrow = n, ncol = length(Pcol),
                dimnames = list(NULL, Pcol))

  sdP  <- apply(Pr, 2, stats::sd)
  tolP <- sqrt(.Machine$double.eps) * max(c(1, sdP[is.finite(sdP)]))
  if (any(!is.finite(sdP)) || any(sdP < tolP))
    stop("After partialling out the exogenous block these endogenous ",
         "regressors have (numerically) zero variance: ",
         paste(Pcol[!is.finite(sdP) | sdP < tolP], collapse = ", "))

  Z <- cbind(yr, Pr)
  colnames(Z) <- c(lhs_label, Pcol)
  if (qr(Z)$rank < ncol(Z))
    stop("The ICA input [y_res, P_res] is rank deficient.")

  ## decomposition and selection of the Gaussian source
  ic <- ica::ica(Z, nc = ncol(Pmat) + 1L, method = method)
  S  <- ic$S
  k  <- select_fn(S)
  if (length(k) != 1L || is.na(k))
    stop("Control-function selection failed to return a single component index.")
  cf <- as.numeric(S[, k])

  ## Second stage.
  if (CF) {
    ## y ~ P + X + control_func
    Pdes  <- Pmat
    extra <- matrix(cf, ncol = 1L, dimnames = list(NULL, "control_func"))
  } else {
    ## Residualise P on the control function without an intercept, P~ = M_cf P.
    ## Since cf is centred this reproduces the Frisch-Waugh-Lovell slopes of the
    ## control-function regression and leaves colMeans(P~) = colMeans(P), so the
    ## reported intercept is the structural one.
    ss <- sum(cf * cf)
    if (!is.finite(ss) || ss <= 0)
      stop("Selected ICA component is degenerate (zero variance).")
    Pdes <- Pmat - outer(cf, as.numeric(crossprod(cf, Pmat)) / ss)
    dimnames(Pdes) <- dimnames(Pmat)
    extra <- Pmat[, 0L, drop = FALSE]
  }

  D <- cbind(Pdes, Xmat, extra)
  if (has_intercept) D <- cbind(`(Intercept)` = rep(1, n), D)

  if (anyDuplicated(colnames(D)))
    stop("Duplicated column names in the second-stage design matrix: ",
         paste(unique(colnames(D)[duplicated(colnames(D))]), collapse = ", "))

  fit <- stats::lm.fit(x = D, y = y)
  b   <- fit$coefficients
  if (is.null(names(b))) names(b) <- colnames(D)

  list(coef = b, cf = cf, k = k, Pmat = Pmat, ica_obj = ic)
}


## One jackknife replicate. A replicate that fails for any reason returns an
## all-NA vector of the correct length rather than aborting the whole run, and
## results are aligned by coefficient name rather than by position.
.IcaReg_jack <- function(i, df, point_cf, coef_names, fit_args) {

  out <- stats::setNames(rep(NA_real_, length(coef_names)), coef_names)

  res <- tryCatch({
    df_i   <- df[-i, , drop = FALSE]
    target <- point_cf[-i]
    sel <- function(S) {
      r <- suppressWarnings(as.numeric(stats::cor(S, target)))
      r[!is.finite(r)] <- 0
      which.max(abs(r))
    }
    do.call(.IcaReg_fit, c(list(df = df_i, select_fn = sel), fit_args))$coef
  }, error = function(e) NULL)

  if (is.null(res) || is.null(names(res))) return(out)
  nm <- intersect(names(res), coef_names)
  out[nm] <- as.numeric(res[nm])
  out
}


## ---- main function ----------------------------------------------------------

IcaReg <- function(formula, data, method = "jade", CF = FALSE,
                   select = c("kurtosis", "ks"), kurt_tol = 1,
                   parallel = FALSE, ncores = NULL, nboots = NULL) {

  ## `nboots` is accepted so that scripts written for earlier versions keep
  ## running. It has no effect: the delete-one jackknife uses exactly n
  ## replicates and has no tuning parameter.

  select  <- match.arg(select)
  sel_fun <- switch(select,
                    kurtosis = .IcaReg_select_kurtosis,
                    ks       = .IcaReg_select_ks)

  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  data <- as.data.frame(data)

  ## ---- parse the formula --------------------------------------------------
  if (length(formula) != 3L)
    stop("'formula' must be two-sided, e.g. y ~ p1 + p2 | x1 + x2.")

  lhs       <- formula[[2]]
  lhs_label <- .IcaReg_dep1(lhs)

  f1            <- nlme::splitFormula(formula, sep = "|")
  has_exog      <- length(f1) > 1
  has_intercept <- attr(stats::terms(f1[[1]]), "intercept") == 1

  if (has_exog && attr(stats::terms(f1[[2]]), "intercept") == 0)
    stop("Suppress the intercept in the FIRST part of the formula, ",
         "e.g. 'y ~ p1 + p2 - 1 | x1 + x2'. A '-1' in the exogenous block is ",
         "ambiguous and is not supported.")

  P_rhs <- .IcaReg_dep1(f1[[1]][[2]])
  X_rhs <- if (has_exog) .IcaReg_dep1(f1[[2]][[2]]) else NULL

  Y_vars <- all.vars(lhs)
  P_vars <- all.vars(f1[[1]])
  X_vars <- if (has_exog) all.vars(f1[[2]]) else character(0)

  ## ---- validity checks -----------------------------------------------------
  missing_vars <- setdiff(c(Y_vars, P_vars, X_vars), names(data))
  if (length(missing_vars) > 0)
    stop("The following variables are missing in the data: ",
         paste(missing_vars, collapse = ", "))

  numeric_P <- vapply(data[P_vars], is.numeric, logical(1))
  if (!all(numeric_P))
    stop("Only continuous variables can be endogenous. The following are not numeric: ",
         paste(P_vars[!numeric_P], collapse = ", "))

  both <- intersect(P_vars, X_vars)
  if (length(both) > 0)
    stop("The following variables appear in BOTH blocks of the formula: ",
         paste(both, collapse = ", "),
         ". Partialling X out of P would make them identically zero, which ",
         "makes the ICA input singular.")

  all_chk <- unique(c(Y_vars, P_vars, X_vars))
  constant_var <- vapply(data[all_chk], function(x) length(unique(x)) == 1L,
                         logical(1))
  if (any(constant_var))
    stop("The following variables are constant: ",
         paste(all_chk[constant_var], collapse = ", "))

  ## ---- analysis data frame -------------------------------------------------
  df <- data[, unique(c(Y_vars, P_vars, X_vars)), drop = FALSE]
  df <- droplevels(stats::na.omit(df))
  n  <- nrow(df)
  if (n < 5L)
    stop("Only ", n, " complete observations remain after na.omit().")

  full_rhs <- if (has_exog) paste(P_rhs, "+", X_rhs) else P_rhs

  ## na.omit() above only sees the raw columns. A term such as log(z) evaluates
  ## to NaN for non-positive z without any NA being present in z itself, and
  ## model.matrix() would then drop those rows from one block but not from
  ## another. Resolve this once, on the full model frame including the response,
  ## so that every design matrix is built from an identical row set.
  mf <- stats::model.frame(
          stats::as.formula(paste(lhs_label, "~", full_rhs)),
          data = df, na.action = stats::na.omit)
  drop_idx <- attr(mf, "na.action")
  if (!is.null(drop_idx)) {
    warning(length(drop_idx), " row(s) dropped: a transformed term evaluated ",
            "to NA/NaN (e.g. log of a non-positive value).", call. = FALSE)
    df <- droplevels(df[-as.integer(drop_idx), , drop = FALSE])
    n  <- nrow(df)
    if (n < 5L)
      stop("Only ", n, " observations remain after removing rows whose ",
           "transformed terms are NA/NaN.")
    mf <- stats::model.frame(
            stats::as.formula(paste(lhs_label, "~", full_rhs)),
            data = df, na.action = stats::na.omit)
  }

  ## Inf survives na.omit() and would silently poison the cumulants.
  num_col <- vapply(mf, is.numeric, logical(1))
  if (any(vapply(mf[num_col], function(z) any(!is.finite(z)), logical(1))))
    stop("Non-finite values (Inf) after evaluating the model terms; ",
         "check for log(0) or division by zero.")

  dm <- stats::model.matrix(stats::as.formula(paste("~", full_rhs)), data = df)
  if (!has_intercept)
    dm <- dm[, setdiff(colnames(dm), "(Intercept)"), drop = FALSE]
  dm_rank <- qr(dm)$rank
  if (dm_rank < ncol(dm))
    stop("Design matrix is rank deficient (rank ", dm_rank, " < ",
         ncol(dm), " columns).")

  ## The argument list must contain no language objects; see .IcaReg_lang above.
  fit_args <- list(lhs_label = lhs_label,
                   P_rhs = P_rhs, X_rhs = X_rhs,
                   has_exog = has_exog, has_intercept = has_intercept,
                   CF = CF, method = method, env = env)

  ## ---- point estimate ------------------------------------------------------
  point <- do.call(.IcaReg_fit,
                   c(list(df = df, select_fn = sel_fun), fit_args))
  Estimates    <- point$coef
  control_func <- point$cf
  coef_names   <- names(Estimates)

  ## ---- jackknife standard errors -------------------------------------------
  message("Estimation done. Calculating jackknife standard errors (", n,
          " replicates)")

  cl <- NULL
  if (parallel) {
    if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
    if (.Platform$OS.type == "windows") {
      ## Windows has no forking: build a PSOCK cluster and export what the
      ## workers need.
      cl <- parallel::makeCluster(ncores)
      parallel::clusterEvalQ(cl, suppressMessages(library(ica)))
      parallel::clusterExport(cl,
                              c(".IcaReg_fit", ".IcaReg_jack", ".IcaReg_lang"),
                              envir = environment())
      on.exit(parallel::stopCluster(cl), add = TRUE)
    } else {
      cl <- ncores   # integer => forking, copy-on-write
    }
  }

  reps <- pbapply::pbsapply(
    seq_len(n),
    function(i) .IcaReg_jack(i, df, point_cf = control_func,
                             coef_names = coef_names, fit_args = fit_args),
    cl = cl
  )

  if (is.null(dim(reps)))
    reps <- matrix(reps, nrow = length(coef_names),
                   dimnames = list(coef_names, NULL))
  reps <- reps[coef_names, , drop = FALSE]

  m_ok     <- rowSums(is.finite(reps))
  n_failed <- sum(colSums(is.finite(reps)) == 0L)

  if (n_failed > 0)
    warning(n_failed, " of ", n, " jackknife replicates failed entirely and ",
            "were dropped from the variance.", call. = FALSE)
  if (any(m_ok < n & m_ok > 0))
    warning("Some coefficients are based on fewer than ", n, " replicates ",
            "(min = ", min(m_ok), "); a factor level may drop under deletion.",
            call. = FALSE)

  theta_bar <- ifelse(m_ok > 0L, rowSums(reps, na.rm = TRUE) / pmax(m_ok, 1L),
                      NA_real_)
  V_jk      <- ifelse(m_ok > 1L,
                      (m_ok - 1L) / pmax(m_ok, 1L) *
                        rowSums((reps - theta_bar)^2, na.rm = TRUE),
                      NA_real_)
  ses <- sqrt(V_jk)

  Estimates1 <- cbind(Estimates, ses)
  colnames(Estimates1) <- c("Estimate", "Std. Error")
  rownames(Estimates1) <- coef_names

  ## ---- identification diagnostics ------------------------------------------
  kurt_all <- apply(point$ica_obj$S, 2, .IcaReg_kurt)
  kurt_cf  <- kurt_all[point$k]

  message("Excess kurtosis of the ICA sources: ",
          paste(sprintf("%+.3f", kurt_all), collapse = "  "),
          "  [control function: source ", point$k, "]")

  if (is.finite(kurt_cf) && abs(kurt_cf) > kurt_tol)
    warning("The selected control function has excess kurtosis ",
            sprintf("%+.3f", kurt_cf), " (threshold ", kurt_tol,
            "). Identification requires one source with a vanishing fourth ",
            "cumulant; no source here looks Gaussian.", call. = FALSE)

  tie_frac <- apply(point$Pmat, 2, function(x) 1 - length(unique(x)) / length(x))
  if (any(tie_frac > 0)) {
    bad <- tie_frac[tie_frac > 0]
    warning("Endogenous regressors contain ties (repeated values): ",
            paste0(names(bad), " ", sprintf("%.1f%%", 100 * bad),
                   collapse = ", "),
            ". ICA identification assumes continuously distributed sources.",
            call. = FALSE)
  }

  return(list(Estimates1, control_func))
}
