### News for package ARCensReg

<font color='grey'>**ARCensReg 3.0.0 (2022-06-25)**</font>

* Incorporate the `ARtCensReg` function to fit a univariate censored linear regression model with autoregressive errors considering Student-t innovations.

* Generic functions `print`, `summary`, and `plot` have methods for objects given as an output of `ARCensReg` and `ARtCensReg` functions.

* Function `residuals` was incorporated to compute the conditional and quantile residuals for objects inheriting from class ARpCRM or ARtpCRM, given as an output of `ARCensReg` and `ARtCensReg` function, respectively.

* Function `plot` is also valid for objects returned by `residuals`. This procedure returns the following plots for the quantile residuals: a dispersion plot (residual vs. time), an autocorrelation plot, a histogram, and a Quantile-Quantile (Q-Q) plot.

* Argument `pit` was substituted by `phi` in the function `rARCens`. See the documentation.

* Function `rARCens` also simulates datasets with Student-t innovations.

* Some modifications in the arguments of the `InfDiag` function were made. The generic function `plot` also works for objects given as an output of function `InfDiag`. 


<font color='grey'>**ARCensReg 2.1 (2016-09-11)**</font>

*


<font color='grey'>**ARCensReg 2.0 (2016-09-09)**</font>

* Incorporate the `rARCens` function to simulate censored autoregressive data.

* Incorporate the `InfDiag` function to perform influence diagnostic by a local influence approach with three possible perturbation schemes: response perturbation, scale matrix perturbation, or explanatory variable perturbation.

* The Phosphorus concentration data was added.


<font color='grey'>**ARCensReg 1.0 (2016-05-16)**</font>

* Initial release.
