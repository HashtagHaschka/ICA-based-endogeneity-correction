The function implements the ICA-based endogeneity correction by Dost & Haschka (2025). The required arguments should be specified as follows:
- formula should be depvar ~ endog_var1 + endog_var2 + ... | exog_var1 + exog_var2 + ...
- depvar is the dependent variable, endog_var1, etc., are the continuous endogenous variables, exog_var1, etc., are the exogenous variables
- formula accepts -1 to remove intercept
- dummy variables can be modelled using as.factor(exog_var1), etc.
- data should be a data.frame
- CF = TRUE indicates if the control function approach should be used, or the control function should be residualised out before
- Example: ica_reg(formula = Y ~ endog1 + endog2 + endog3 | exog1 + exog2 + as.factor(exog3), data = data1)

REFERENCES
- Dost, F. and Haschka, R. E. (2025). XYC. SSRN Working Paper. 
