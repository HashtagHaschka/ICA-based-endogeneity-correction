## -----------------------------------------------------------------------------
##  ICA-based endogeneity correction in R -- examples
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
##
##  HOW TO RUN
##
##  Open this file in R or RStudio and run it. The estimator is downloaded
##  automatically from GitHub, all datasets come from freely available packages.
##
##  Standard errors are delete-one jackknife standard errors, so the estimator
##  refits the model once per observation. Runtime therefore grows with the
##  sample size. Examples 1 to 4 finish in roughly a minute in total. Examples 5
##  and 6 use samples of about 28,000 and 9,600 observations and take
##  considerably longer; set run_slow to TRUE below to include them.
##
## -----------------------------------------------------------------------------

run_slow <- FALSE


## ---- load the estimator and the packages used by the examples ---------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

source(paste0("https://raw.githubusercontent.com/HashtagHaschka/",
              "ICA-based-endogeneity-correction/main/IcaReg.R"))

pacman::p_load(dplyr, AER, ISLR, mlbench, bayesm)


## =============================================================================
##  1. Simulated data -- does the correction recover the truth?
## =============================================================================
##
##  Two endogenous regressors z1 and z2, each the sum of a non-normal exogenous
##  part (lognormal) and a normal confounder xi. The confounder also enters the
##  structural equation, so OLS is biased. True coefficients: intercept 2,
##  z1 = 1, z2 = 1, x1 = 1, x2 = 1.

N <- 100

testdata <- as.data.frame(matrix(NA, nrow = N, ncol = 0))
xi1 <- rnorm(N)
xi2 <- rnorm(N)
testdata$z1 <- rlnorm(N) + xi1
testdata$z2 <- rlnorm(N) + xi2
testdata$x1 <- rnorm(N)
testdata$x2 <- rnorm(N)
testdata$e  <- xi1 + xi2
testdata$Y  <- 2 + testdata$z1 + testdata$z2 +
                   testdata$x1 + testdata$x2 + testdata$e

## OLS: biased upward on z1 and z2
coef(lm(Y ~ z1 + z2 + x1 + x2, data = testdata))

## ICA correction
mod1 <- IcaReg(formula = Y ~ z1 + z2 | x1 + x2, data = testdata)
mod1[[1]]

## the recovered confounder
hist(mod1[[2]], main = "Estimated control function", xlab = "")

## CF = TRUE additionally reports the coefficient on the control function.
## It is a Hausman-type test for endogeneity; all other coefficients are
## unchanged.
mod2 <- IcaReg(formula = Y ~ z1 + z2 | x1 + x2, data = testdata, CF = TRUE)
mod2[[1]]

## The two selection rules should agree whenever identification is clear-cut.
## Comparing them is a cheap robustness check.
mod3 <- IcaReg(formula = Y ~ z1 + z2 | x1 + x2, data = testdata,
               select = "ks")
cbind(kurtosis = mod1[[1]][, "Estimate"],
      ks       = mod3[[1]][rownames(mod1[[1]]), "Estimate"])


## =============================================================================
##  2. Carseats (ISLR) -- child car seat sales across 400 stores
## =============================================================================
##
##  Price and competitor price are treated as endogenous: unobserved local
##  demand conditions plausibly drive both prices and sales.

data("Carseats")
dat1 <- Carseats
dat1$lsales     <- log(dat1$Sales)
dat1$lprice     <- log(dat1$Price)
dat1$lcompprice <- log(dat1$CompPrice)
dat1 <- subset(dat1, is.finite(lsales))

mod0 <- lm(lsales ~ lprice + lcompprice + ShelveLoc + Income + Advertising +
             Population + Age + Education + Urban + US, data = dat1)
summary(mod0)

mod1 <- IcaReg(formula = lsales ~ lprice + lcompprice | ShelveLoc + Income +
                 Advertising + Population + Age + Education + Urban + US,
               data = dat1)
mod1[[1]]
hist(mod1[[2]], main = "Estimated control function", xlab = "")


## =============================================================================
##  3. Boston housing (mlbench) -- median home values across 506 tracts
## =============================================================================
##
##  Note on the data rather than the method: zn is zero for about three quarters
##  of the tracts and indus is recorded at town level, so both carry many
##  repeated values. Expect the ties warning. The example is included because it
##  is widely known, not because it is a good illustration of a continuous
##  endogenous regressor.

data("BostonHousing")
dat1 <- BostonHousing

mod0 <- lm(medv ~ crim + zn + indus + chas + nox + rm + age + dis + rad +
             tax + ptratio + b + lstat, data = dat1)
summary(mod0)

mod1 <- IcaReg(formula = medv ~ crim + zn + indus | chas + nox + rm + age +
                 dis + rad + tax + ptratio + b + lstat, data = dat1)
mod1[[1]]
hist(mod1[[2]], main = "Estimated control function", xlab = "")


## =============================================================================
##  4. Orange juice (bayesm) -- weekly store-level scanner data
## =============================================================================
##
##  Log sales of Tropicana on its own price and the prices of the competing
##  brands. Prices respond to unobserved local demand shocks, so OLS is suspect.
##  A single store gives roughly 120 weekly observations.
##
##  Note that a variable must not appear in both blocks of the formula: the
##  endogenous price is not repeated among the exogenous controls.

data("orangeJuice")
dat1 <- orangeJuice[[1]]

## store with the highest sales
dat1 %>%
  mutate(sales = exp(logmove)) %>%
  group_by(store) %>%
  summarise(total_sales = sum(sales, na.rm = TRUE)) %>%
  arrange(desc(total_sales)) %>%
  slice(1)

dat111 <- dat1 %>%
  filter(store == 111, brand == 4) %>%          # brand 4 = Tropicana
  mutate(across(starts_with("price"), log))

mod0 <- lm(logmove ~ price4 + price1 + price2 + price3 + price5 + price6 +
             price7 + price8 + price9 + price10 + price11 + deal + feat,
           data = dat111)
summary(mod0)

mod1 <- IcaReg(formula = logmove ~ price4 | price1 + price2 + price3 + price5 +
                 price6 + price7 + price8 + price9 + price10 + price11 +
                 deal + feat, data = dat111)
mod1[[1]]
hist(mod1[[2]], main = "Estimated control function", xlab = "")

## most profitable store
dat1 %>%
  group_by(store) %>%
  summarise(total_profit = sum(profit, na.rm = TRUE)) %>%
  arrange(desc(total_profit)) %>%
  slice(1)

dat124 <- dat1 %>%
  filter(store == 124, brand == 4) %>%
  mutate(across(starts_with("price"), log))

mod0 <- lm(logmove ~ price4 + price1 + price2 + price3 + price5 + price6 +
             price7 + price8 + price9 + price10 + price11 + deal + feat,
           data = dat124)
summary(mod0)

mod1 <- IcaReg(formula = logmove ~ price4 | price1 + price2 + price3 + price5 +
                 price6 + price7 + price8 + price9 + price10 + price11 +
                 deal + feat, data = dat124)
mod1[[1]]
hist(mod1[[2]], main = "Estimated control function", xlab = "")


## =============================================================================
##  The remaining examples use large samples. Set run_slow <- TRUE at the top
##  of this file to run them.
## =============================================================================

if (isTRUE(run_slow)) {

  ## ===========================================================================
  ##  5. CPS 1988 (AER) -- returns to schooling, about 28,000 observations
  ## ===========================================================================
  ##
  ##  Education is endogenous in the classic sense: unobserved ability raises
  ##  both schooling and wages. Roughly 28,000 jackknife replicates, so run this
  ##  in parallel and expect it to take a while.

  data("CPS1988")
  data1 <- CPS1988
  data1$lwage         <- log(data1$wage)
  data1$experience_sq <- data1$experience^2

  mod0 <- lm(lwage ~ education + experience + experience_sq + ethnicity +
               smsa + region + parttime, data = data1)
  summary(mod0)

  mod1 <- IcaReg(formula = lwage ~ education | experience + experience_sq +
                   ethnicity + smsa + region + parttime,
                 data = data1, parallel = TRUE)
  mod1[[1]]
  hist(mod1[[2]], main = "Estimated control function", xlab = "")


  ## ===========================================================================
  ##  6. Orange juice, all stores pooled
  ## ===========================================================================
  ##
  ##  Pooling stores and weeks gives about 9,600 observations. Note that the
  ##  jackknife assumes independent observations; with repeated observations per
  ##  store the standard errors will understate the true sampling variability.

  dat_all <- dat1 %>%
    filter(brand == 4) %>%
    mutate(across(starts_with("price"), log))

  mod0 <- lm(logmove ~ price4 + price1 + price2 + price3 + price5 + price6 +
               price7 + price8 + price9 + price10 + price11 + deal + feat,
             data = dat_all)
  summary(mod0)

  ## one endogenous price
  mod1 <- IcaReg(formula = logmove ~ price4 | price1 + price2 + price3 +
                   price5 + price6 + price7 + price8 + price9 + price10 +
                   price11 + deal + feat, data = dat_all, parallel = TRUE)
  mod1[[1]]

  ## all prices endogenous
  mod2 <- IcaReg(formula = logmove ~ price4 + price1 + price2 + price3 +
                   price5 + price6 + price7 + price8 + price9 + price10 +
                   price11 | deal + feat, data = dat_all, parallel = TRUE)
  mod2[[1]]

  ## all prices endogenous, with week fixed effects
  mod3 <- IcaReg(formula = logmove ~ price4 + price1 + price2 + price3 +
                   price5 + price6 + price7 + price8 + price9 + price10 +
                   price11 | deal + feat + as.factor(week),
                 data = dat_all, parallel = TRUE)
  mod3[[1]]
  hist(mod3[[2]], main = "Estimated control function", xlab = "")

}
