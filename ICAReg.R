# =============================================================================
#  ICA-based endogeneity correction with delete-one JACKKNIFE standard errors
# =============================================================================

# Load required packages
pacman::p_load(
  nlme,
  dplyr,
  ica,
  Matrix,
  pbapply,
  parallel
)

# -----------------------------------------------------------------------------
#  Helper 1: pick the Gaussian source in the FULL sample (most normal column),
#            via the Kolmogorov-Smirnov distance to normality
# -----------------------------------------------------------------------------
select_gaussian_ks <- function(S) {
  ks <- apply(S, 2, function(x) {
    suppressWarnings(stats::ks.test(x, "pnorm", mean = mean(x), sd = sd(x))$statistic)
  })
  which.min(ks)
}

# -----------------------------------------------------------------------------
#  Helper 2: the core estimator. Runs the WHOLE pipeline once (first-stage OLS
#            partialling-out of exogenous regressors if present, ICA, control-
#            function selection, final regression) and returns the coefficient
#            vector together with the selected control function.
#
#  `select_fn(S)` returns the index of the column of the ICA source matrix S that
#  is to be used as the control function.  In the full sample this is the KS rule;
#  in a leave-one-out replicate it is correlation-matching to the full-sample source.
# -----------------------------------------------------------------------------
ica_fit_core <- function(df, dep, P_rhs, X_rhs, has_exog, has_intercept,
                         CF, method, select_fn) {
  
  drop_int <- if (!has_intercept) " - 1" else ""
  
  if (!has_exog) {
    
    ## ---- No exogenous regressors -------------------------------------------
    # ICA input = [Y , design matrix of the endogenous part without intercept]
    Pmm  <- model.matrix(as.formula(paste("~", P_rhs)), data = df)
    Pcol <- setdiff(colnames(Pmm), "(Intercept)")
    Pmat <- Pmm[, Pcol, drop = FALSE]
    Yv   <- df[[dep]]
    
    ic <- ica::ica(cbind(Yv, Pmat), nc = ncol(Pmat) + 1L, method = method)
    S  <- ic$S
    cf <- S[, select_fn(S)]
    
    if (CF) {
      # Control-function regression:  Y ~ P + control_func
      reg <- as.data.frame(Pmat)
      reg[[dep]]        <- Yv
      reg$control_func  <- cf
      rhs  <- paste(c(Pcol, "control_func"), collapse = " + ")
      fit  <- lm(as.formula(paste(dep, "~", rhs, drop_int)), data = reg)
    } else {
      # Residualise each endogenous regressor on control_func, then Y ~ P_tilde
      Pres <- Pmat
      for (j in seq_len(ncol(Pmat))) Pres[, j] <- residuals(lm(Pmat[, j] ~ cf))
      reg <- as.data.frame(Pres)
      reg[[dep]] <- Yv
      rhs  <- paste(Pcol, collapse = " + ")
      fit  <- lm(as.formula(paste(dep, "~", rhs, drop_int)), data = reg)
    }
    
  } else {
    
    ## ---- With exogenous regressors -----------------------------------------
    # First stage: partial X out of Y and of each endogenous regressor (OLS).
    Xmm  <- model.matrix(as.formula(paste("~", X_rhs)), data = df)
    qrX  <- qr(Xmm)
    Yres <- qr.resid(qrX, df[[dep]])
    
    Pmm  <- model.matrix(as.formula(paste("~", P_rhs)), data = df)
    Pcol <- setdiff(colnames(Pmm), "(Intercept)")
    Pmat <- Pmm[, Pcol, drop = FALSE]
    Pres <- qr.resid(qrX, Pmat)
    
    # ICA on the first-stage residuals = [Y_res , P_res]
    ic <- ica::ica(cbind(Yres, Pres), nc = ncol(Pmat) + 1L, method = method)
    S  <- ic$S
    cf <- S[, select_fn(S)]
    
    if (CF) {
      # Y ~ P + X + control_func   (X enters as in lm(as.Formula(Y ~ P | X)))
      reg <- df
      reg$control_func <- cf
      rhs <- paste(P_rhs, "+", X_rhs, "+ control_func")
      fit <- lm(as.formula(paste(dep, "~", rhs, drop_int)), data = reg)
    } else {
      # Residualise raw endogenous regressors on control_func, then Y ~ P_tilde + X
      Pcf <- Pmat
      for (j in seq_len(ncol(Pmat))) Pcf[, j] <- residuals(lm(Pmat[, j] ~ cf))
      reg <- df
      reg[, Pcol] <- Pcf
      rhs <- paste(P_rhs, "+", X_rhs)
      fit <- lm(as.formula(paste(dep, "~", rhs, drop_int)), data = reg)
    }
  }
  
  list(coef = coef(fit), cf = cf)
}

# -----------------------------------------------------------------------------
#  Helper 3: one jackknife replicate. Deletes observation i, re-runs the whole
#            pipeline, and aligns the ICA component by correlation-matching the
#            replicate sources to the full-sample control function (minus row i).
# -----------------------------------------------------------------------------
jack_replicate <- function(i, df, point_cf, fit_args) {
  df_i   <- df[-i, , drop = FALSE]
  target <- point_cf[-i]
  sel    <- function(S) which.max(abs(as.numeric(cor(S, target))))
  do.call(ica_fit_core, c(list(df = df_i, select_fn = sel), fit_args))$coef
}

# -----------------------------------------------------------------------------
#  Main function
# -----------------------------------------------------------------------------
ica_reg <- function(formula, data, method = "jade", CF = FALSE,
                    parallel = FALSE, ncores = NULL, nboots = NULL) {
  # `nboots` is accepted for backward compatibility but ignored: the delete-one
  # jackknife uses exactly n replicates and has no tuning parameter.
  
  ## ---- parse the formula --------------------------------------------------
  f1            <- nlme::splitFormula(formula, sep = "|")
  has_intercept <- attr(terms(f1[[1]]), "intercept") == 1
  has_exog      <- length(f1) > 1
  dep           <- all.vars(formula)[1]
  P_rhs         <- paste(deparse(f1[[1]][[2]]), collapse = " ")
  X_rhs         <- if (has_exog) paste(deparse(f1[[2]][[2]]), collapse = " ") else NULL
  
  P_vars <- all.vars(f1[[1]])
  X_vars <- if (has_exog) all.vars(f1[[2]]) else character(0)
  
  ## ---- validity checks (as in the original) -------------------------------
  missing_vars <- setdiff(c(P_vars, X_vars), names(data))
  if (length(missing_vars) > 0)
    stop(paste("The following variables are missing in the data:",
               paste(missing_vars, collapse = ", ")))
  
  numeric_P  <- sapply(data[P_vars], is.numeric)
  if (!all(numeric_P))
    stop("Only continuous variables can be endogenous. The following are not numeric: ",
         paste(P_vars[!numeric_P], collapse = ", "))
  
  all_chk      <- c(P_vars, X_vars)
  constant_var <- sapply(data[all_chk], function(x) length(unique(x)) == 1)
  if (any(constant_var))
    stop("The following variables are constant: ",
         paste(all_chk[constant_var], collapse = ", "))
  
  ## ---- assemble analysis data frame ---------------------------------------
  df <- data %>%
    dplyr::select(dplyr::all_of(c(dep, P_vars, X_vars))) %>%
    na.omit() %>%
    as.data.frame()
  
  ## ---- full-sample design rank check --------------------------------------
  full_formula <- as.formula(paste(dep, "~",
                                   if (has_exog) paste(P_rhs, "+", X_rhs) else P_rhs))
  dm   <- model.matrix(full_formula, data = df)
  if (Matrix::rankMatrix(crossprod(dm))[1] != ncol(dm))
    stop("Design matrix is rank deficient")
  
  n        <- nrow(df)
  fit_args <- list(dep = dep, P_rhs = P_rhs, X_rhs = X_rhs,
                   has_exog = has_exog, has_intercept = has_intercept,
                   CF = CF, method = method)
  
  ## ---- point estimate (KS selection of the Gaussian source) ---------------
  point <- do.call(ica_fit_core,
                   c(list(df = df, select_fn = select_gaussian_ks), fit_args))
  Estimates    <- point$coef
  control_func <- point$cf
  coef_names   <- names(Estimates)
  
  ## ---- jackknife standard errors ------------------------------------------
  message("Estimation done. Calculating jackknife standard errors (", n, " replicates)")
  
  cl <- NULL
  if (parallel) {
    if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
    if (.Platform$OS.type == "windows") {
      # Windows has no forking: build a PSOCK cluster and export the workers' needs.
      cl <- parallel::makeCluster(ncores)
      parallel::clusterEvalQ(cl, suppressMessages(library(ica)))
      parallel::clusterExport(cl, c("ica_fit_core", "jack_replicate"),
                              envir = environment())
      on.exit(parallel::stopCluster(cl), add = TRUE)
    } else {
      cl <- ncores   # integer => forking via mclapply inside pbsapply (copy-on-write)
    }
  }
  
  reps <- pbapply::pbsapply(
    seq_len(n),
    function(i) jack_replicate(i, df, point_cf = control_func, fit_args = fit_args),
    cl = cl
  )
  
  ## ---- assemble jackknife variance ----------------------------------------
  if (is.list(reps)) {
    # A replicate returned a different coefficient set (e.g. a factor level
    # vanished under deletion). Align by name, fill gaps with NA, and warn.
    reps <- sapply(reps, function(v) {
      out <- setNames(rep(NA_real_, length(coef_names)), coef_names)
      out[names(v)] <- v
      out
    })
    warning("Some jackknife replicates produced a different set of coefficients ",
            "(a factor level may have dropped under deletion); aligned by name.",
            call. = FALSE)
  }
  if (is.null(dim(reps))) reps <- matrix(reps, nrow = 1,
                                         dimnames = list(coef_names, NULL))
  reps      <- reps[coef_names, , drop = FALSE]
  theta_bar <- rowMeans(reps, na.rm = TRUE)
  V_jk      <- (n - 1) / n * rowSums((reps - theta_bar)^2, na.rm = TRUE)
  ses       <- sqrt(V_jk)
  
  Estimates1 <- cbind(Estimates, ses)
  colnames(Estimates1) <- c("Estimate", "Std. Error")
  
  ## ---- identification diagnostics (as in the original) --------------------
  ks_test1 <- suppressWarnings(stats::ks.test(control_func, "pnorm"))
  if (ks_test1$p.value < .1)
    warning("Joint component may not be normally distributed: Kolmogorov-Smirnov p = ",
            ks_test1$p.value, call. = FALSE)
  if (any(duplicated(control_func)))
    warning("Endogenous regressors contain ties (repeated values)", call. = FALSE)
  
  return(list(Estimates1, control_func))
}


### Simulated data -----------------------------------------------------------

N <- 100
testdata <- as.data.frame(matrix(NA, nrow = N, ncol = 0))

xi1 <- rnorm(N)
xi2 <- rnorm(N)

testdata$z1 <- rlnorm(N) + xi1
testdata$z2 <- rlnorm(N) + xi2

testdata$x1 <- rnorm(N)
testdata$x2 <- rnorm(N)

testdata$e <- xi1 + xi2
testdata$Y <- 2 + testdata$z1 + testdata$z2 + testdata$x1 + testdata$x2 + testdata$e

mod1 <- ica_reg(formula = Y ~ z1 + z2, data = testdata)

mod1[[1]]        # estimates with jackknife standard errors
hist(mod1[[2]])  # control function


### 1. Real world data example -----------------------------------------------

library(AER)

data("CPS1988")
data1 <- CPS1988
data1$lwage <- log(data1$wage)
data1$experience_sq <- data1$experience^2

mod1 <- lm(lwage ~ education + experience + experience_sq + ethnicity + smsa + region + parttime,
           data1)
# n ~ 28k => 28k re-estimations; run in parallel
mod2 <- ica_reg(lwage ~ education | experience + experience_sq + ethnicity + smsa + region + parttime,
                data = data1, parallel = TRUE)
mod2[[1]]
hist(mod2[[2]])


### 2. Real world data example -----------------------------------------------

library(ISLR)

data("Carseats")
dat1 <- Carseats

dat1$lsales <- log(dat1$Sales)
dat1$lprice <- log(dat1$Price)
dat1$lcompprice <- log(dat1$CompPrice)
dat1 <- subset(dat1, is.finite(log(Sales)))

mod1 <- lm(lsales ~ lprice + lcompprice + ShelveLoc + Income + Advertising + Population + Age + Education + Urban + US,
           dat1)
mod2 <- ica_reg(lsales ~ lprice + lcompprice | ShelveLoc + Income + Advertising + Population + Age + Education + Urban + US,
                data = dat1)
mod2[[1]]
hist(mod2[[2]])


### 3. Real world data example -----------------------------------------------

library(mlbench)

data("BostonHousing")
dat1 <- BostonHousing

mod1 <- lm(medv ~ crim + zn + indus + chas + nox + rm + age + dis + rad + tax +
             ptratio + b + lstat, dat1)
mod2 <- ica_reg(medv ~ crim + zn + indus | chas + nox + rm + age + dis + rad + tax +
                  ptratio + b + lstat, data = dat1)
mod2[[1]]
hist(mod2[[2]])


### 4. Real world example ----------------------------------------------------

library(bayesm)
# https://www.rdocumentation.org/packages/bayesm/versions/3.1-6/topics/orangeJuice

data("orangeJuice")
dat1 <- orangeJuice[[1]]

# store with highest sales
dat1 %>%
  mutate(sales = exp(logmove)) %>%
  group_by(store) %>%
  summarise(total_sales = sum(sales, na.rm = TRUE)) %>%
  arrange(desc(total_sales)) %>%
  slice(1)

dat111 <- dat1 %>% filter(store == 111) %>% filter(brand == 4) %>%
  mutate(across(starts_with("price"), log)) # Tropicana

mod0 <- lm(logmove ~ price4 + price1 + price2 + price3 + price4 + price5 +
             price6 + price7 + price8 + price9 + price10 + price11 + deal + feat,
           dat111)
mod1 <- ica_reg(formula = logmove ~ price4 | price1 + price2 + price3 +
                  price4 + price5 + price6 + price7 + price8 + price9 +
                  price10 + price11 + deal + feat, data = dat111)
mod1[[1]]
hist(mod1[[2]])


# most profitable store
dat1 %>%
  mutate(profit = profit) %>%
  group_by(store) %>%
  summarise(total_profit = sum(profit, na.rm = TRUE)) %>%
  arrange(desc(total_profit)) %>%
  slice(1)

dat124 <- dat1 %>% filter(store == 124) %>% filter(brand == 4) %>%
  mutate(across(starts_with("price"), log)) # Tropicana

mod0 <- lm(logmove ~ price4 + price1 + price2 + price3 + price4 + price5 +
             price6 + price7 + price8 + price9 + price10 + price11 + deal + feat,
           dat124)
mod1 <- ica_reg(formula = logmove ~ price4 | price1 + price2 + price3 +
                  price4 + price5 + price6 + price7 + price8 + price9 +
                  price10 + price11 + deal + feat, data = dat124)
mod1[[1]]
hist(mod1[[2]])


# all stores (large n => parallel recommended)
dat_all <- dat1 %>% filter(brand == 4) %>% mutate(across(starts_with("price"),
                                                         log)) # Tropicana
mod0 <- lm(logmove ~ price4 + price1 + price2 + price3 + price4 + price5 +
             price6 + price7 + price8 + price9 + price10 + price11 + deal + feat,
           dat_all)
mod1 <- ica_reg(formula = logmove ~ price4 | price1 + price2 + price3 +
                  price4 + price5 + price6 + price7 + price8 + price9 +
                  price10 + price11 + deal + feat, data = dat_all, parallel = TRUE)
mod1[[1]]
hist(mod1[[2]])

mod2 <- ica_reg(formula = logmove ~ price4 + price1 + price2 + price3 +
                  price4 + price5 + price6 + price7 + price8 + price9 +
                  price10 + price11 | deal + feat, data = dat_all, parallel = TRUE)
mod2[[1]]
hist(mod2[[2]])

mod3 <- ica_reg(formula = logmove ~ price4 + price1 + price2 + price3 +
                  price4 + price5 + price6 + price7 + price8 + price9 +
                  price10 + price11 | deal + feat + as.factor(week), data = dat_all,
                parallel = TRUE)
mod3[[1]]
hist(mod3[[2]])
