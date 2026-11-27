
# Smith POT MLE reproduction


**Theoretical specification, implementation and saved numerical results**

2021 Pacific Northwest heatwave | Existing regional ERA5 data | Prepared 19 November 2026

## Scope of this document

This report explains the actual reproduction in ERA5_Smith_PNW_POT_MLE.R and helpers/smith_pnw_pot.R, section by section. It describes the saved MLE run from 17 November 2026. No model is refitted in preparing this document. The separate Bayesian implementation and the earlier fixed-threshold POT analysis are outside its numerical scope.

The target is Smith's Section 0.6.2: a seasonal nonstationary point-process (PP) model fitted to declustered daily temperature extremes from 1950-2020. The existing regional ERA5 series and smoothed GMST are used. The 2021 data supply an evaluation temperature and prediction season, and are excluded from fitting.

## Reading guide

| Stage | Statistical purpose | Code sections |

| --- | --- | --- |

| 1. Inputs and seasonal design | Define observations, exposure and covariates | 1-2 |

| 2. Quantile threshold | Choose a changing high-temperature boundary | 3 |

| 3. Declustering | Reduce persistence to distinct event peaks | 4 |

| 4. PP model and MLE | Estimate event rate and tail severity jointly | 5 |

| 5. Count assessment | Compare historical counts with fitted expectations | 6 |

| 6. Attribution and endpoints | Compare fixed climate states in 2021 | 7 |

| 7. Numerical products | Interpret saved estimates and their limitations | 8 |

**Main result.** The fitted seasonal endpoints are 36.0205 C in the factual climate and 33.6877 C in the counterfactual climate. Both lie below the observed 2021 regional seasonal maximum, 39.5174 C. Consequently the plug-in MLE assigns zero probability in both scenarios at that observed level, and RR/FAR there are undefined. These are model outputs, not physical temperature limits.


# 1. Observations, season and exposure


## What enters the analysis

smith_read_daily reads the one-dimensional tmax series from ERA5_reference_daily_tmax.nc, checks Celsius units and a regular daily date axis, and creates date, year and a 365-day seasonal coordinate. Days after February in leap years are shifted back by one, so the same calendar date has the same seasonal phase.

smith_prepare keeps 1950-2020 and retains a year only when all selected seasonal days and its GMST value are available. This saved run contains 71 complete seasons, each with 121 days: 8,591 daily exposure records. A season without peaks would still contribute all its exposure days to the likelihood.

**Actual window:** window_mode="paper_doy" uses days 135-255 inclusive on the 365-day calendar, corresponding to May 15-September 12. The script also contains May 15-September 15 calendar settings, but they are inactive in this run. The paper describes those calendar dates while also specifying days 135-255; the reconstruction follows the numbered-day limits.

## Covariates

\[
\phi_d=2\pi(d-1)/365,\quad c_{1d}=\sin\phi_d,\quad c_{2d}=\cos\phi_d
\]

\[
c_{3d}=\sin(2\phi_d),\quad c_{4d}=\cos(2\phi_d),\quad z_y=g_y-g_0
\]

Here d is the no-leap day of year; g_y is the existing smoothed annual GMST; and g_0 is the mean GMST over retained fitting years (approximately 0.2857717 C). Centering changes the interpretation of the intercept, not the fitted climate response. The threshold trend is a separate covariate, t_y=(y-1950)/10.

The input GMST file is consumed as prepared; this script does not estimate the smoothing itself. The prediction states are g_F=0.9322916667 C and g_C=-0.2677083333 C, with a prescribed difference of 1.2 C.

## Event versus observation

A daily row represents observation time, not an independent heatwave. A declustered peak represents an event in the PP event term. A seasonal probability later means at least one such event above a specified temperature during this 121-day window. The implementation is therefore neither an annual-maxima fit nor a fit to all raw exceedance days as independent events.


# 2. Estimating the changing threshold


smith_threshold fits a conditional 0.90 quantile by quantile regression using all 8,591 selected seasonal temperatures. There is no initial restriction to temperatures above 30 C. The value 0.90 is a probability level, whereas u_y,d is a fitted temperature in C that varies by date and year.

\[
u_{y,d}=a_0+a_Tt_y+a_1c_{1d}+a_2c_{2d}+a_3c_{3d}+a_4c_{4d}
\]

\[
\widehat{a}=\arg\min_a\sum_{y,d}\rho_{0.9}(X_{y,d}-u_{y,d}(a))
\]

\[
\rho_\tau(r)=r\,[\tau-\mathbf{1}(r<0)]
\]

The asymmetric check loss penalizes positive residuals more strongly than negative residuals at tau=0.90. quantreg::rq uses method="fn". The fitted surface approximates the conditional 90th percentile jointly across seasons; it is not a separately estimated empirical percentile for each calendar day.

| Threshold term | Saved estimate |

| --- | --- |

| (Intercept) | 22.5954 |

| trend_decades | 0.289688 |

| sin1 | 0.081947 |

| cos1 | -2.1969 |

| sin2 | 3.68741 |

| cos2 | 1.59013 |

The fitted thresholds range from 20.30429 to 30.47474 C. A daily exceedance is defined strictly by X_y,d > u_y,d; there are 858 such days. The fraction is close to 10%, but the method does not require exactly 10% exceedances at every date or in every year.

**Why 30 C appears:** it is a historical count-check level, not the threshold fitted above. Smith selected tau after trying alternatives and assessing high-temperature counts. The exact quantile-regression formula was not published; this script explicitly reconstructs it using a linear calendar-year trend and two harmonic pairs. Threshold coefficients are estimated first and then held fixed during PP optimization. Their uncertainty is not propagated in these MLE products.


# 3. Declustering persistent hot days


Successive exceedance days often belong to the same heat episode. Treating them as independent events would inflate counts and distort event-rate inference. smith_decluster applies texmex's Ferro-Segers intervals procedure to the residual series R_j=X_j-u_j, using zero as its exceedance threshold.

## Intervals estimator and automatic gap

Let j_1,...,j_N be exceedance positions on the concatenated seasonal observation axis, and T_i=j_(i+1)-j_i their interexceedance intervals. For the branch used when the largest interval exceeds two, texmex computes:

\[
\widehat{\theta}=\min\left\{1,\frac{2[\sum_{i=1}^{N-1}(T_i-1)]^2}{(N-1)\sum_{i=1}^{N-1}(T_i-1)(T_i-2)}\right\}
\]

The extremal index theta summarizes clustering; heuristically its reciprocal is a limiting mean cluster size at high thresholds. It is not the probability of a heatwave. The observed value is 0.2822494. texmex chooses a gap cutoff using ordered intervals and the estimated number of clusters, with tie handling. In this run the cutoff r is 4: a new automatic cluster starts when T_i > 4. This describes index gaps on the concatenated seasonal series, not four winter-aware calendar days.

## Selecting the event peak

The helper then forces a boundary between different years, giving 244 retained clusters rather than 240 automatic clusters. For each resulting cluster it selects the largest actual temperature, not the largest residual X-u. The selected row supplies its own date, GMST covariates and threshold to the event likelihood.

| Saved quantity | Value |

| --- | --- |

| Exceedance days before declustering | 858 |

| Automatic clusters before year splitting | 240 |

| Retained event peaks | 244 |

| Mean exceedance days per retained event | 858 / 244 = 3.5164 |

The peak table retains dates, start/end dates, threshold, excess and exceedance_days for inspection. Only the peak row and its temperature enter the PP event term; exceedance_days and duration are not intensity weights. This reproduction models peak temperature severity, not cumulative heat.

**Implementation assumption:** winter gaps are omitted when estimating the extremal index on concatenated seasons. Year-boundary splitting is imposed afterward. This is an explicit reconstruction convention; neither the paper nor the saved fit establishes that all residual dependence has been removed.


# 4. The seasonal PP model


smith_fit_pp fits ten parameters jointly. Define a vector of six location coefficients, three scale coefficients and one shape. GMST enters location; the seasonal cycle enters both location and scale. The PP scale uses an identity link, with positivity required on fitted days.

\[
\mu_{y,d}=b_0+b_Gz_y+b_1c_{1d}+b_2c_{2d}+b_3c_{3d}+b_4c_{4d}
\]

\[
\sigma_d=s_0+s_1c_{1d}+s_2c_{2d}>0,\qquad \xi=\mathrm{constant}
\]

## Tail measure and intensity

\[
h_{y,d}(x)=1+\xi\frac{x-\mu_{y,d}}{\sigma_d},\qquad A_{y,d}(x)=h_{y,d}(x)^{-1/\xi}
\]

\[
\lambda_{y,d}(x)=-\frac{\partial A_{y,d}(x)}{\partial x}=\frac{1}{\sigma_d}h_{y,d}(x)^{-1/\xi-1}
\]

These expressions apply inside the support h(x)>0. A(x) is an annualized expected count above x per unit model time, and lambda(x) is its temperature density. Neither is a normalized PDF for daily temperatures. With daily exposure Delta=1/365, the expected event count on a day is approximately Delta A(x).

At xi=0 the continuous limits are A(x)=exp[-(x-mu)/sigma] and lambda(x)=A(x)/sigma. For xi<0, A(x) is exactly zero at and above the finite upper endpoint mu-sigma/xi. This support rule drives the zero probabilities at high evaluation temperatures.

## Where the GPD appears

\[
\Pr(X>u+e\mid X>u,y,d)=\frac{A_{y,d}(u+e)}{A_{y,d}(u)}=\left(1+\xi\frac{e}{\sigma_u(y,d)}\right)^{-1/\xi}
\]

\[
\sigma_u(y,d)=\sigma_d+\xi[u_{y,d}-\mu_{y,d}]
\]

\[
f_E(e\mid y,d)=\frac{1}{\sigma_u(y,d)}\left(1+\xi\frac{e}{\sigma_u(y,d)}\right)^{-1/\xi-1}
\]

The conditional excess density is therefore GPD with scale sigma_u and shape xi. This is a consequence of the fitted PP model, not a second GPD fit after pp.fit. Even though PP scale has no GMST term, the conditional GPD scale can vary with climate through mu and the threshold. Keeping xi constant is a model restriction that pools sparse tail information, not evidence that the physical tail shape cannot change.

The excess density applies for e>0 with 1+xi e/sigma_u>0. For xi<0 this means e<-sigma_u/xi; for xi=0 its limit is the exponential density exp(-e/sigma_u)/sigma_u.


# 5. Constructing and maximizing the likelihood


Let J be all 8,591 exposure days and P the 244 retained peak rows. A Poisson-process likelihood combines the probability of the observed event count with the event temperatures. Up to constants independent of the parameters:

\[
L(\vartheta)\ \propto\ \exp\left[-\frac{1}{365}\sum_{j\in J}A_j(u_j)\right]\prod_{i\in P}\lambda_i(x_i)
\]

\[
\mathrm{NLL}(\vartheta)=\frac{1}{365}\sum_{j\in J}A_j(u_j)+\sum_{i\in P}\left[\log\sigma_i+\left(1+\frac{1}{\xi}\right)\log h_i(x_i)\right]
\]

The first term is the expected number of threshold events over all observed seasons. It carries the information from quiet days and years. The second term rewards parameter values assigning suitable intensity to the selected event temperatures. Parameters jointly control event frequency and peak severity; mu is not estimated by a separate regression of observed peaks.

Equivalently, on day j the threshold count rate is A_j(u_j), while the conditional excess density is lambda_j(x)/A_j(u_j). Multiplying count and severity contributions gives the same PP likelihood. This explains why fitting only the conditional GPD would lose rate information needed for seasonal attribution.

## How pp.fit receives these data

The helper constructs xdat equal to u_j for every nonpeak row and equal to the original temperature for peak rows. Because pp.fit selects strict exceedances, only the 244 peaks enter the event term; all 8,591 rows remain in the exposure sum. Replacing nonpeaks in xdat does not claim their observed temperatures equalled the threshold. It is an implementation device for removing their event contributions.

The covariate matrix contains centered GMST and four harmonics. mul=1:5 selects all five for location; sigl=2:3 selects annual harmonics for scale; shl=NULL keeps constant shape. npy=365 sets the time normalization. The extremal index is not multiplied into the fitted rates again after declustering.

## Optimization and checks

The routine uses Nelder-Mead, maxit=50,000 and reltol=1e-10, with initial shape values -0.20, -0.05 and 0.10. Location starts combine a least-squares representation of the fitted threshold with a rate-based intercept shift. The scale start is max(SD of peak temperatures, 0.5), with zero seasonal scale coefficients. It selects the converged run with the lowest NLL and independently reevaluates the likelihood with smith_nll, including its xi=0 limit and support checks.

Numerical convergence does not prove a global optimum: the three saved starts converge to different objective values. These products contain no bootstrap intervals or Bayesian uncertainty. The threshold and declustering decisions are treated as fixed in the PP MLE.


# 6. Historical count assessment


smith_count_check integrates the fitted daily tail measures over all 71 historical seasons, retaining each year's own GMST. At each level x it compares the predicted count with the number of retained cluster peaks exceeding x:

\[
\widehat{E}_{\rm hist}(x)=\frac{1}{365}\sum_{j\in J}A_j(x),\qquad O(x)=\sum_{i\in P}\mathbf{1}(x_i>x)
\]

\[
C(x)=\frac{[O(x)-\widehat{E}_{\rm hist}(x)]^2}{\widehat{E}_{\rm hist}(x)}
\]

C(x) is a descriptive Pearson component, reported only for positive finite expectations. The same events contribute to multiple nested levels, parameters are estimated from the data, and expected tail counts can be very small. The saved components therefore do not have an automatically valid independent chi-square testing interpretation.

| x (C) | Observed | Expected | C(x) | Days with u>x |

| --- | --- | --- | --- | --- |

| 30 | 79 | 78.7081 | 0.00108219 | 310 |

| 30.5 | 66 | 57.3778 | 1.29565 | 0 |

| 31 | 47 | 40.2384 | 1.13619 | 0 |

| 31.5 | 27 | 26.9024 | 0.000354029 | 0 |

| 32 | 18 | 16.9423 | 0.0660356 | 0 |

| 32.5 | 8 | 9.88453 | 0.359294 | 0 |

| 33 | 3 | 5.22302 | 0.946159 | 0 |

| 33.5 | 1 | 2.42849 | 0.840266 | 0 |

| 34 | 0 | 0.961457 | 0.961457 | 0 |

| 35 | 0 | 0.0631818 | 0.0631818 | 0 |

| 36 | 0 | 5.25104e-07 | 5.25104e-07 | 0 |

At 30 C, 310 historical days have thresholds above the evaluation level. On those days the model is extrapolated below its fitted tail threshold; the output flags this. At 30.5 C and higher, all evaluation levels lie above the fitted historical thresholds.

The agreement at 30 C is close (79 observed, 78.708 expected). At 30.5 and 31 C the components exceed one, whereas Smith reports components below one at temperatures of 30 C or greater. This aspect is not exactly numerically reproduced. At 34 C there are zero observed historical peaks but an expected count near 0.96; the absence of events at a high level does not require an expectation of zero.

**GOF status:** these are count diagnostics, not calibrated GOF p-values. This POT MLE product does not contain a KS test or a fitted-parameter bootstrap calibration. Parameter convergence and count agreement cannot establish accurate extrapolation to the 2021 extreme.


# 7. Attribution and upper endpoints


For prediction, the same 121 seasonal dates are assigned either factual 2021 GMST or GMST reduced by 1.2 C. The coefficients, harmonics, scale function and shape remain fixed between scenarios. Since only location depends directly on GMST, mu_F(d)-mu_C(d)=1.2 b_G, approximately 2.3328 C.

\[
\Lambda_s(x)=\frac{1}{365}\sum_{d=135}^{255}A_{s,d}(x),\qquad s\in\{F,C\}
\]

\[
p_s(x)=\Pr(N_s(x)\geq1)=1-\exp[-\Lambda_s(x)]
\]

Lambda is the expected number of cluster events above x in the season. p is the probability of at least one. They are close for rare events but are distinct quantities. The saved return period 1/p is a waiting time in repeated seasons at a fixed climate state, not a forecast under future warming.

\[
RR(x)=\frac{p_F(x)}{p_C(x)},\quad FAR(x)=1-\frac{p_C(x)}{p_F(x)},\quad RR^{-1}(x)=\frac{p_C(x)}{p_F(x)}
\]

smith_attribution computes probabilities and ratios in log space to reduce numerical underflow. If p_C=0 while p_F>0, it reports RR=Inf, FAR=1 and inverse_RR=0. If both probabilities are zero, all ratios and FAR are NA because 0/0 is undefined. The MLE does not use the separate posterior-summary convention that sometimes replaces a zero-over-zero inverse ratio by zero.

## Finite endpoint calculation

\[
U_s(d)=\mu_s(d)-\frac{\sigma_d}{\xi},\qquad U_s^{\rm season}=\max_{d=135,\ldots,255}U_s(d)\quad(\xi<0)
\]

The maximum daily endpoint is the seasonal upper support limit on the implementation's daily grid. Above it, every daily tail measure is zero and so is the seasonal probability. The bound is estimated from the fitted distribution; it is neither a confidence limit nor a physical upper bound on future heat.

The evaluation vector contains 34-39 C and the observed 2021 seasonal maximum, 39.517395 C on 29 June. These evaluation temperatures do not define the threshold used in fitting. No 2021 observations are used to force the fitted endpoint upward.

The counterfactual is a prescribed 1.2 C GMST perturbation interpreted through the fitted climate response. Its attribution meaning depends on that scenario and on the adequacy of the GMST-based model; it is not independently identified from this likelihood alone.


# 8. Saved MLE estimates and comparison


| PP parameter | Saved estimate |

| --- | --- |

| mu_intercept | 33.299 |

| mu_gmst_centered | 1.944 |

| mu_sin1 | 2.85199 |

| mu_cos1 | 4.47533 |

| mu_sin2 | 4.00342 |

| mu_cos2 | 2.29544 |

| sigma_intercept | 2.07042 |

| sigma_sin1 | 0.607845 |

| sigma_cos1 | 1.08128 |

| xi | -0.386048 |

The intercept is defined at centered GMST zero and zero-valued seasonal basis functions. It is not a directly observed summer mean. The coefficient b_G=1.9440 gives the fitted location response per 1 C increase in GMST. sigma coefficients have temperature units; xi is dimensionless. The shape estimate is negative, implying finite fitted upper endpoints.

| Initial shape | Convergence | Negative log likelihood |

| --- | --- | --- |

| -0.2 | 0 | 58.4328 |

| -0.05 | 0 | 57.8042 |

| 0.1 | 0 | 60.4281 |

All three starts report convergence code zero; the selected run starts at -0.05 and has NLL 57.80424775. The objective contains the package's time-normalization convention and parameter-independent constants are omitted. Its absolute value is not a standalone GOF score or directly comparable with an annual GEV likelihood on different observations.

| Seasonal endpoint (C) | Reproduction MLE | Smith reported MLE | Difference |

| --- | --- | --- | --- |

| Factual | 36.02046 | 36.4 | -0.37954 |

| Counterfactual | 33.68766 | 34.0 | -0.31234 |

The qualitative endpoint issue is reproduced, but exact endpoint values are not. Smith's Table 0.2 is based on Bayesian posterior calculations; its nonzero posterior mean probabilities must not be used as direct targets for this plug-in MLE table. Posterior parameter uncertainty can admit support above x even when the MLE does not.


# 9. Attribution products and interpretation


| x (C) | Lambda F | p F | p C | RR | FAR |

| --- | --- | --- | --- | --- | --- |

| 34 | 0.0796396 | 0.0765509 | 0 | Inf | 1 |

| 35 | 0.00987036 | 0.00982181 | 0 | Inf | 1 |

| 36 | 5.74548e-08 | 5.74548e-08 | 0 | Inf | 1 |

| 37 | 0 | 0 | 0 | NA | NA |

| 38 | 0 | 0 | 0 | NA | NA |

| 39 | 0 | 0 | 0 | NA | NA |

| 39.5174* | 0 | 0 | 0 | NA | NA |

* Observed regional seasonal maximum on 29 June 2021. Counterfactual expected event counts are zero at every displayed level. NA denotes an undefined ratio; Inf denotes a model-implied infinite quantity, not a finite estimate with high precision.

| x (C) | Factual return period | Counterfactual return period |

| --- | --- | --- |

| 34 | 13.0632 | Inf |

| 35 | 101.814 | Inf |

| 36 | 1.7405e+07 | Inf |

| 37 | Inf | Inf |

| 38 | Inf | Inf |

| 39 | Inf | Inf |

| 39.5174* | Inf | Inf |

At 34 C, Lambda_F=0.0796396 but p_F=0.0765509; their difference is visible. At 36 C, the level is only about 0.0205 C below the fitted factual endpoint and its probability is approximately 5.74548e-8. The corresponding waiting time is approximately 17.4 million fixed-climate seasons. Such an extrapolated number is extremely sensitive to the endpoint and is not a validated forecast.

At 37 C and above, both climates have zero MLE probability. At the actual event level this produces neither a finite RR nor evidence for an infinite RR: the ratio is undefined. The zero arises from the fitted support rather than a finite simulation failing to observe the event.

## What the result establishes

The reconstruction produces a coherent saved PP likelihood fit, declustered event sample, count checks and climate-state comparisons. It reproduces the qualitative finding that the pre-2021 fitted endpoint falls below the 2021 heatwave. It does not establish exact numerical reproduction, calibrated goodness of fit, or precise event probabilities at the observed extreme. No uncertainty intervals are provided by this MLE run.


# 10. Reproducibility and remaining assumptions


## Code-to-theory map

| Helper | Role |

| --- | --- |

| smith_read_daily / smith_prepare | Read inputs, select complete seasons and retain exposure. |

| smith_features / smith_threshold | Construct covariates and estimate the 0.90 quantile surface. |

| smith_decluster | Estimate extremal index, split clusters and select temperature peaks. |

| smith_parameters / smith_nll | Map coefficients to parameters and evaluate the PP objective. |

| smith_fit_pp | Optimize pp.fit from multiple starts and check its likelihood. |

| smith_count_check | Compare historical observed and predicted cluster counts. |

| smith_expected_count / smith_attribution | Integrate scenario counts and derive probabilities, RR and FAR. |

| smith_daily_parameters | Produce daily parameters and endpoint summaries. |

The main script is Code/smith_reproduce/ERA5_Smith_PNW_POT_MLE.R. Its only analysis helper is Code/smith_reproduce/helpers/smith_pnw_pot.R. Saved numerical inputs to this report are in Results/ERA5_Smith_PNW_POT_MLE; the RDS snapshot also preserves settings, fitting data, coverage, threshold model, clusters, fitted model and scenario data.

## Interpretation limits that affect reproduction

The local ERA5 series and smoothed GMST are available inputs, but exact identity with Smith's historical KNMI-derived series and GMST processing has not been established. A spatial-domain difference cannot be inferred from quoted outer bounds alone, since they may describe cell edges. Daily temperature definitions, regional averaging and time-window conventions remain possible sources of differences, not demonstrated explanations.

The paper does not specify every reconstruction detail: the exact threshold design matrix, winter-gap treatment and all computational conventions are not available. The daily integration, harmonic design and year splitting described here are the actual local choices. The likelihood further assumes a suitable Poisson approximation for retained peaks conditional on covariates. Declustering is a practical step toward this assumption, not its proof.

## References and provenance

[1] Richard L. Smith. Detection and Attribution of Extreme Weather Events. Revised 18 April 2025. Section 0.6.2, printed pp. 24-26 (PDF pp. 26-28). Local source: C:/Users/siyan/Downloads/SmithDetectionAttributionPreprint.pdf.

[2] ismev::pp.fit documentation and local implementation. The reproduction uses its PP parameterization, npy=365 and identity parameter links. Documentation: https://search.r-project.org/CRAN/refmans/ismev/html/pp.fit.html

[3] texmex::extremalIndex and declust.extremalIndex, locally inspected version 2.4.9. The intervals estimator and gap rule in Section 3 follow this implementation, associated with Ferro-Segers intervals declustering.

[4] Saved local products: threshold_coefficients.csv; extremal_index.csv; declustered_peaks.csv; pp_parameters.csv; optimizer_audit.csv; historical_count_check.csv; attribution_mle.csv; prediction_daily_parameters.csv; prediction_endpoints.csv; smith_pnw_pot_mle.rds. Results were read for this report; no analysis output was overwritten.
