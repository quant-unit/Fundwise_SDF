Referee: 1

Comments to the Author
Referee report “Identification and Stability of SDF Estimators for Pooled,
Non-traded Cash Flows”, submitted to the Journal of Financial Econometrics

Summary
This paper develops a semiparametric Least-Mean-Distance framework for estimating stochastic discount factors from pooled, non-traded cash flows, generalizing the estimator of Driessen, Lin and Phalippou (2012). The central theoretical contribution is to consider an estimator that averages pricing errors over discounting dates beyond fund inception. This could stabilize finite-sample estimation but shifts the population target through a timing-related covariance mechanism. Applied to Preqin data on Buyout and Venture Capital funds, single-factor market models yield stable betas (roughly 0.8–1.2 for Buyout, 1.0–1.9 for Venture Capital), but adding an alpha intercept or second factor produces unstable, sign-reversing estimates, indicative of weak identification.

Discussion
Various econometric approaches to assess the risk and return of private equity funds (and more generally non-traded assets) have been proposed. This study focuses on two of these methods: the methods of Driessen, Lin and Phalippou (DLP, 2012) and Korteweg and Nagel (KN, 2016). Using a simulation study the authors analyze whether discounting cash flows to a set of dates, instead of a single date, helps the small-sample performance of the estimator. The drawback of such multiple-date discounting is that the estimator may be asymptotically inconsistent.
Given the different existing approaches and the lack of consensus in the literature, it is very welcome that the authors study these methods in further detail. However, as discussed below in more detail, the current version of the paper does not clearly integrate the different approaches.
I would advise the authors to: (i) describe the similarities and differences of the approaches mathematically more clearly, (ii) derive the conditions under which using multiple discounting dates gives consistent estimates, and (iii) extend the simulation study to compare the different approaches. I will elaborate on these issues below.

1: Foundation of the different approaches
The KN-2016 approach is fundamentally different from the DLP-2012 approach. KN start with the SDF pricing equation: P = E[MX] with P price, M the SDF and X the payoff. Their benchmark analysis focuses on a pricing kernel that is a linear function of the market return: M = a + bR, with R the market return. Their GMM approach gives estimates of the parameters of the SDF, when fitted to private equity payoff data. One can test if the pricing errors of this model are zero. This approach does not give direct estimates of the alpha and beta of the non-traded asset, but these can be obtained with some additional assumptions.
The DLP-2012 approach does not use the SDF approach. Instead, it makes assumptions about the return-generating process and then derives an estimation procedure to estimate the parameters of this process, using cross-sectional asymptotics. In case of a linear one-factor return model, this gives estimates of the alpha and beta of the non-traded asset. This approach thus explicitly incorporates a “pricing error” alpha, which is not part of the SDF estimation approach.
This paper seems to mix these approaches, and this is confusing and unnecessary. The theorems in section 1 start with the SDF approach (equations 1 and 2). However, the authors then compare to DLP-2012 in section 1.4.1, but as discussed above this is not an SDF-based approach. More importantly, in their empirical analysis the authors define the SDF using equation (19), but this equation uses the return-generating process of DLP as an SDF. Specifically, with inception date tau=0, the “SDF” here is the inverse of (1+alpha+rf+beta*R).
In other words, calculating an NPV with realized market returns (or factor returns), which is also what the well-known PME method does, is not automatically the same as using the SDF approach.
It would therefore be good if the authors can discuss the differences of the approaches, properly set up their approach and relate their horizon-based extension clearly to these approaches.

2: Compounding-horizon wedge
Part (3) of Lemma 1 introduces the compounding-horizon wedge: if one evaluates the pricing error at a date beyond the inception date, the pricing error may not have a zero expectation because of a covariance term. Because of this, the estimation approach that the authors propose may be asymptotically inconsistent (as acknowledged by the authors). It would be good if the authors work out the conditions under which consistency is obtained. My guess is that if returns follow an i.i.d. linear factor process with idiosyncratic shocks that are independent from the factors the estimator is consistent (when the cross-section of funds per time period tends to infinity), also when using multiple discounting dates. This is important because it is not very satisfactory that the new approach proposed by the authors can lead to inconsistent estimates.

3: Simulation study
It would be useful if the simulation study would also compare the different approaches in more detail. Ideally, some comparison of the KN-2016 and DLP-2012 approaches would be useful, ideally combined with the horizon extension.

4: Empirical analysis
A comparison with existing empirical results is missing. It would be useful to have this.

5: Writing
The paper is sometimes hard to follow.
