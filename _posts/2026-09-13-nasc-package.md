---
layout: post
title: Network-Aware Synthetic Control (NASC)
date: 2026-09-13
description: R package for bayesian Network-Aware Synthetic Control estimation when SUTVA is violated as treatment spillovers travel through networks. It corrects donor-pool contamination and recovers implied spillovers onto untreated units by including an treatment exposure penalty in the synthetic weight estimation and a bias correction factor.
categories: packages
package: true
thumbnail: assets/img/publication_preview/Tankrabatt3.png
_styles: |
  .nasc-figure-bottom-right {
    float: right;
    max-width: 32%;
    margin: 1rem 0 0.75rem 1.25rem;
  }

  .nasc-top-figure {
    justify-content: center;
  }

  .nasc-top-figure figure {
    max-width: 72%;
    margin: 0 auto;
  }

  .post-meta,
  .post-tags {
    display: none;
  }

  .post {
    display: grid;
    grid-template-columns: 9rem minmax(0, 1fr);
    column-gap: 3rem;
    align-items: start;
    width: calc(100% + 2rem);
  }

  .post-header {
    display: flex;
    flex-direction: column;
  }

  .post-title {
    color: var(--global-text-color-light);
    font-size: 0.85rem;
    font-weight: 500;
    letter-spacing: 0.08em;
    line-height: 1.4;
    text-transform: uppercase;
    margin: 2rem 0 0;
    padding-top: 0.2rem;
  }

  .post-content {
    width: calc(100% + 5rem);
  }

  .nasc-top-figure img,
  .nasc-figure-bottom-right img {
    mix-blend-mode: multiply;
    box-shadow: none !important;
  }

  html[data-theme="dark"] .nasc-top-figure img,
  html[data-theme="dark"] .nasc-figure-bottom-right img {
    filter: invert(1);
    mix-blend-mode: screen;
  }

  .nasc-figure-bottom-right figure {
    margin: 0;
  }

  .nasc-repository-card {
    clear: both;
    margin-top: 2rem;
  }

  .post-content h4 {
    color: var(--global-text-color-light);
    font-size: 0.78rem;
    font-weight: 500;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    margin: 2rem 0 0.75rem;
  }

  .post-content mjx-container[display="true"] {
    overflow-x: auto;
    overflow-y: hidden;
    padding: 0.35rem 0;
  }

  @media (max-width: 767px) {
    .post {
      display: block;
      width: 100%;
    }

    .post-header {
      margin-bottom: 1rem;
    }

    .post-content {
      width: 100%;
    }

    .nasc-figure-bottom-right {
      float: none;
      max-width: 100%;
      margin: 1rem 0;
    }
  }
---

<div class="row nasc-top-figure">
  <div class="col-12 mt-3">
    {% include figure.liquid path="assets/img/small world networks.png" title="nasc multiple solutions" class="img-fluid rounded z-depth-1" %}
  </div>
</div>

<br>

The synthetic control method recovers credible counterfactuals under SUTVA, an
assumption rarely plausible in spatial settings such as transport investments,
place-based policies, or regional taxation. Effects propagate through trade, commuting,
and competition networks to precisely those donor units from which the counterfactual
is constructed, so that a convex combination of contaminated donors inherits part of
the policy's impact and yields a biased treatment effect estimate.

The Network-Aware Synthetic Control (NASC) estimator embeds the donor pool in a spatial
autoregressive outcome model on a known network and addresses interference through two
complementary channels. A contamination penalty discounts each donor in proportion to
its network exposure, sidestepping the necessity of restricting the donor pool ex ante,
while a bias correction accounts for the remaining contamination that the selected
donors inherit. The estimator requires neither a clean donor nor an exposure pattern
specified beyond the network itself, and its bias vanishes asymptotically under perfect
pre-treatment balance.

#### Exposure along the network

Interference is assumed to travel across a known undirected graph with row-standardized
weight matrix $$W$$. The common factor model underlying the conventional synthetic control
is embedded in a spatial Durbin process, so that potential outcomes are pre-multiplied by
the spatial multiplier $$(I_{J+1}-\rho W)^{-1}$$,

$$
y_t(0) = (I_{J+1} - \rho W)^{-1}\big(Z\beta_t + WZ\theta_t + M\lambda_t + \nu_t\iota_{J+1} + \varepsilon_t\big),
$$

where the autoregressive parameter $$\rho \in (-1,1)$$ proxies the strength of spillovers.
Partitioning $$W$$ into the links $$w_{J1}$$ running from the treated unit to the donors and
the donor-to-donor block $$W_J$$ yields the *contamination vector*

$$
s \;\equiv\; \rho\,(I_J - \rho W_J)^{-1} w_{J1},
\qquad
\delta_t \;=\; s\,\tau_{1t},
$$

so that donor $$i$$ inherits the share $$s_i$$ of the treatment effect. Expanding the
inverse as a Neumann series $$\sum_{k\ge 0}(\rho W_J)^k$$ decomposes $$s_i$$ into successive
orders of network exposure: it accumulates over paths of increasing length connecting
donor $$i$$ to the treated unit, each discounted by a power of $$\rho$$, so exposure decays
with network distance.

This makes the bias of the conventional estimator explicit. For a weight vector
$$\gamma^{sc}$$ on the simplex, the synthetic control estimate decomposes as

$$
\tau^{sc}_{1t} \;=\; \kappa^{sc}_{1t} \;+\; \big(1 - \gamma^{sc\prime} s\big)\,\tau_{1t},
$$

where $$\kappa^{sc}_{1t}$$ is the usual post-treatment prediction error. Interference enters
as a multiplicative distortion through the scalar $$\gamma^{sc\prime}s$$, additively separable
from $$\kappa^{sc}_{1t}$$ and, unlike it, not vanishing as the pre-treatment window grows. For
positive spillovers the distortion is strictly attenuating; for negative spillovers its sign
depends on the network and the realized weights.

#### Regularization and correction

NASC attacks that scalar from both sides. First, the synthetic control objective is augmented
with a penalty on the weighted absolute contamination, discounting each donor in proportion
to its exposure:

$$
\hat\gamma^{\,nasc}(\lambda) \;=\; \arg\min_{\gamma \in \Delta^{J-1}}\;
\lVert X_1 - X_0\gamma \rVert_V^2 \;+\; \lambda\,\gamma'\lvert s\rvert .
$$

As $$\lambda \to 0$$ the objective collapses to the conventional synthetic control; as
$$\lambda$$ grows, mass shifts from heavily exposed donors towards weakly contaminated ones.
Donor selection is thereby endogenized rather than restricted ex ante: a contaminated donor
still contributes whenever its value for the pre-treatment fit outweighs the contamination it
introduces. Any $$\lambda > 0$$ additionally resolves the familiar non-uniqueness of synthetic
weights by selecting the least contaminated among fit-equivalent solutions, and $$\lambda$$
itself is chosen by cross-validation on the pre-treatment period.

<div class="nasc-figure-bottom-right">
  {% include figure.liquid path="assets/img/nasc multiple solutions.png" title="nasc multiple solutions" class="img-fluid rounded z-depth-1" %}
</div>

Second, a bias correction rescales the contamination the selected donors still carry. With
$$\tau^{cr}_{1t} = y_{1t} - \gamma^{nasc\prime}y^{d}_{t}$$ denoting the penalized synthetic
control estimate, the NASC estimators of the direct and indirect effects are

$$
\hat\tau^{\,nasc}_{1t} \;=\; \frac{\tau^{cr}_{1t}}{1-\gamma^{nasc\prime}s},
\qquad
\hat\delta^{\,nasc}_{t} \;=\; \frac{s}{1-\gamma^{nasc\prime}s}\,\tau^{cr}_{1t}.
$$

Because regularization has already shrunk $$\gamma'\lvert s\rvert$$, the correction factor sits
closer to unity and amplifies the remaining prediction error less than correcting unpenalized
weights would. The second expression delivers the implied spillover onto every untreated unit
as a by-product. Under perfect pre-treatment balance the estimation error reduces to
$$(1-\gamma^{nasc\prime}s)^{-1}\kappa^{nasc}_{1t}$$, which is bounded by a quantity that vanishes
as the pre-treatment period expands.

#### Modular Bayesian inference

Since $$s$$ is a function of the estimated $$\rho$$, the correction is itself a random variable.
Inference places an exponentially tilted Dirichlet prior on the weights, which encodes the
contamination penalty directly on the simplex,

$$
\gamma \mid \lambda, \rho \;\sim\; \mathrm{ETDir}\big(\alpha, \lambda\lvert s\rvert\big)
\;\propto\; \mathrm{Dir}(\gamma;\alpha)\,\exp\!\big(-\lambda\,\gamma'\lvert s\rvert\big),
\qquad \gamma \in \Delta^{J-1},
$$

and whose maximum-a-posteriori estimate under a uniform base measure reproduces the frequentist
NASC weights exactly. Estimation then proceeds along a cut posterior,

$$
p_{\mathrm{cut}}(\zeta \mid Y, Z, W) \;=\; p(\zeta_2 \mid X_1, W, \zeta_1)\;p(\zeta_1 \mid Y, Z, W),
$$

with $$\zeta_1=\{\rho,\beta,\sigma_1\}$$ estimated on pre-treatment data in an initial module and
the weights $$\zeta_2=\{\gamma,\lambda,\sigma_2\}$$ in a second module conditional on it.
Information flows strictly forward: uncertainty about $$\rho$$ propagates into the weights and
the bias correction, while the synthetic control matching cannot feed back and misspecify the
underlying network dependence. Monte Carlo simulations across regular and small-world
topologies show that the proposed estimator outperforms the conventional and the Bayesian
synthetic control across various interference structures.

The method is implemented in the R package nasc, which extends the estimator to a
spatial Durbin random-effects panel specification. NASC is open-source software. Visit
the repository for installation instructions, examples, and the latest development
version.


<div class="nasc-repository-card repositories d-flex flex-wrap flex-md-row flex-column justify-content-center">
  {% include repository/repo.liquid repository="fvwaldow/nasc" %}
</div>
