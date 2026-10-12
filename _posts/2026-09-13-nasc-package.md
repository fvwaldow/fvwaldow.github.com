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

  .nasc-explorer {
    clear: both;
    margin: 1.5rem 0 2rem;
    border: 1px solid var(--global-divider-color);
    border-radius: 6px;
    padding: 1rem 1.1rem 1.2rem;
  }

  .nasc-controls {
    display: flex;
    flex-direction: column;
    gap: 0.6rem;
    padding-bottom: 0.9rem;
    border-bottom: 1px solid var(--global-divider-color);
  }

  .nasc-row {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 0.45rem 0.9rem;
  }

  .nasc-label {
    color: var(--global-text-color-light);
    font-size: 0.65rem;
    font-weight: 500;
    letter-spacing: 0.08em;
    text-transform: uppercase;
  }

  .nasc-explorer button {
    background: transparent;
    border: 1px solid var(--global-divider-color);
    border-radius: 3px;
    color: var(--global-text-color);
    cursor: pointer;
    font-size: 0.72rem;
    letter-spacing: 0.02em;
    padding: 0.22rem 0.6rem;
  }

  .nasc-explorer button:hover {
    border-color: var(--global-text-color-light);
  }

  .nasc-explorer button.is-active {
    background: var(--global-text-color);
    border-color: var(--global-text-color);
    color: var(--global-bg-color);
  }

  .nasc-explorer button.nasc-ghost {
    margin-left: auto;
    color: var(--global-text-color-light);
  }

  .nasc-sliders label,
  .nasc-inline-slider label {
    display: flex;
    align-items: center;
    gap: 0.45rem;
    margin: 0;
  }

  .nasc-sliders input[type="range"],
  .nasc-inline-slider input[type="range"] {
    width: 108px;
    accent-color: var(--global-theme-color);
  }

  .nasc-sliders output,
  .nasc-inline-slider output {
    color: var(--global-text-color);
    font-size: 0.72rem;
    font-variant-numeric: tabular-nums;
    min-width: 2.2rem;
  }

  .nasc-inline-slider {
    display: flex;
    justify-content: center;
    margin-top: 0.4rem;
  }

  .nasc-inline-slider input[type="range"] {
    width: 150px;
  }

  .nasc-grid {
    display: grid;
    grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
    gap: 0.4rem 1.4rem;
    margin-top: 0.9rem;
  }

  .nasc-cell {
    margin: 0;
    min-width: 0;
  }

  .nasc-cell canvas {
    display: block;
    touch-action: pan-y;
    width: 100%;
  }

  .nasc-cell figcaption {
    color: var(--global-text-color-light);
    font-size: 0.68rem;
    line-height: 1.45;
    margin-top: 0.3rem;
  }

  .nasc-key {
    display: inline-block;
    height: 2px;
    margin: 0 0.3rem 0.18rem 0.55rem;
    vertical-align: middle;
    width: 14px;
  }

  .nasc-key:first-child {
    margin-left: 0;
  }

  .nasc-key-obs {
    background: var(--global-text-color);
  }
  .nasc-key-truth {
    background: var(--global-text-color-light);
  }
  .nasc-key-sc {
    background: #e07a3f;
  }
  .nasc-key-nasc {
    background: #3f7fb5;
  }

  .nasc-readout {
    align-content: start;
    display: grid;
    gap: 0.9rem;
    padding-top: 0.4rem;
  }

  .nasc-table {
    border-collapse: collapse;
    font-size: 0.73rem;
    width: 100%;
  }

  .nasc-table th,
  .nasc-table td {
    border-bottom: 1px solid var(--global-divider-color);
    padding: 0.22rem 0;
  }

  .nasc-table th {
    color: var(--global-text-color-light);
    font-weight: 400;
    text-align: left;
  }

  .nasc-table td {
    font-variant-numeric: tabular-nums;
    text-align: right;
  }

  .nasc-fallback {
    color: var(--global-text-color-light);
    font-size: 0.8rem;
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

    .nasc-grid {
      grid-template-columns: minmax(0, 1fr);
      gap: 1.1rem;
    }

    .nasc-explorer button.nasc-ghost {
      margin-left: 0;
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
weight matrix $$W$$. The common factor model underlying the conventional synthetic  control (see Abadie et al. 2010)
is embedded in a spatial Durbin model,

$$
y_t(0) = (I_{J+1} - \rho W)^{-1}\big(Z\beta_t + WZ\theta_t + M\lambda_t + \nu_t\iota_{J+1} + \varepsilon_t\big),
$$

where the autoregressive parameter $$\rho \in (-1,1)$$ proxies the strength of spillovers.
Partitioning $$W$$ into the links $$w_{J1}$$ running from the treated unit to the donors and
the donor-to-donor block $$W_J$$ yields the contamination vector

$$
s \;\equiv\; \rho\,(I_J - \rho W_J)^{-1} w_{J1},
\qquad
\delta_t \;=\; s\,\tau_{1t},
$$

so that donor $$i$$ inherits the share $$s_i$$ of the treatment effect. Expanding the
inverse as a Neumann series $$\sum_{k\ge 0}(\rho W_J)^k$$ decomposes $$s_i$$ into successive
orders of network exposure, accumulating over paths of increasing length connecting
donor $$i$$ to the treated unit, each discounted by a power of $$\rho$$, so exposure decays
with network distance.

This makes the bias of the conventional estimator explicit. For a weight vector
$$\gamma^{sc}$$ on the simplex, the synthetic control estimate decomposes as

$$
\tau^{sc}_{1t} \;=\; \kappa^{sc}_{1t} \;+\; \big(1 - \gamma^{sc\prime} s\big)\,\tau_{1t},
$$

where $$\kappa^{sc}_{1t}$$ is the usual post-treatment prediction error. Interference enters
as a multiplicative bias through the scalar $$\gamma^{sc\prime}s$$, additively separable
from $$\kappa^{sc}_{1t}$$ and, unlike it, not vanishing as the pre-treatment period grows. For
positive spillovers the bias is strictly attenuating; for negative spillovers its sign
depends on the network and the realized weights.

#### Regularization and correction

NASC accounts for the interference bias via two channels. First, the synthetic control objective is augmented
with a penalty on the weighted absolute contamination, discounting each donor in proportion
to its exposure:

$$
\hat\gamma^{\,nasc}(\lambda) \;=\; \arg\min_{\gamma \in \Delta^{J-1}}\;
\lVert X_1 - X_0\gamma \rVert_V^2 \;+\; \lambda\,\gamma'\lvert s\rvert .
$$

As $$\lambda \to 0$$ the objective collapses to the conventional synthetic control; as
$$\lambda$$ grows, mass shifts from heavily exposed donors towards weakly contaminated ones.
Donor selection is thereby endogenized rather than restricted ex ante. A contaminated donor
still contributes whenever its value for the pre-treatment fit outweighs the contamination it
introduces. Any $$\lambda > 0$$ additionally resolves the familiar non-uniqueness of synthetic
weights by selecting the least contaminated among fit-equivalent solutions, and $$\lambda$$
itself is chosen by cross-validation on the pre-treatment period.

Second, a bias correction rescales the contamination the selected donors still inherit. With
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

and whose maximum-a-posteriori estimate under $$\alpha = 1$$ reproduces the frequentist
NASC weights exactly. To preserve uncertainty quantification estimation then proceeds along a cut posterior,

$$
p_{\mathrm{cut}}(\zeta \mid Y, Z, W) \;=\; p(\zeta_2 \mid X_1, W, \zeta_1)\;p(\zeta_1 \mid Y, Z, W),
$$

with $$\zeta_1=\{\rho,\beta,\sigma_1\}$$ estimated on pre-treatment data in an initial module and
the weights $$\zeta_2=\{\gamma,\lambda,\sigma_2\}$$ in a second module conditional on the first.
Information flows strictly forward as uncertainty about $$\rho$$ propagates into the weights and
the bias correction, while the synthetic control matching cannot feed back and misspecify the
underlying network dependence. Monte Carlo simulations across regular and small-world
topologies show that the proposed estimator outperforms the conventional and the Bayesian
synthetic control across various interference structures.

#### Hands-On

The hands-on example below calculates the NASC estimator on a four-unit network with a single treated unit. In this case the weight simplex $$\Delta^{2}$$ reduces to a two-dimensional triangle. Toggle edges to change the network,
and the contamination vector, the weights, the counterfactual and the tilted prior all
update together.

Donor untreated potential outcomes solve the spatial autoregression

$$
y^{d}_{t}(0) = \rho\,w_{J1}\,y_{1t}(0) + \rho\,W_{J}\,y^{d}_{t}(0) + \mu + Z_{t}'\beta + \varepsilon_{t},
$$

with $$\beta = (1, 0.5)'$$, unit effects $$\mu_i \sim N(0,1)$$, standard normal covariates
and $$\varepsilon_{it} \sim N(0, 0.5^2)$$. The treated unit's own untreated
path is the convex combination $$\gamma^{*\prime}(\mu + Z_t'\beta)$$ with
$$\gamma^{*} = (0.2, 0.5, 0.3)'$$, which places the largest true weight on the donor
adjacent to the treated unit. Donors that contribute substantially to the 
pre-treatment fit are also the most contaminated. Treatment occurs at $$t = 25$$
with $$\tau_t \sim N(5,1)$$, and donor $$i$$ inherits $$s_i\tau_t$$.

<div id="nasc-explorer" class="nasc-explorer">
  <div class="nasc-controls">
    <div class="nasc-row">
      <span class="nasc-label">Network</span>
      <button type="button" data-nasc-preset="dense">Dense</button>
      <button type="button" data-nasc-preset="asymmetric" class="is-active">Asymmetric</button>
      <button type="button" data-nasc-preset="partial">Partial</button>
      <button type="button" data-nasc-preset="none">No interference</button>
      <button type="button" data-nasc="reseed" class="nasc-ghost">New sample</button>
    </div>
    <div class="nasc-row nasc-sliders">
      <label>
        <span class="nasc-label">Spillover ρ</span>
        <input type="range" min="-90" max="90" step="5" value="60" data-nasc="rho" />
        <output data-nasc="rhoVal">0.60</output>
      </label>
      <label>
        <span class="nasc-label">Penalty λ</span>
        <input type="range" min="0" max="500" step="5" value="25" data-nasc="lambda" />
        <output data-nasc="lambdaVal">25</output>
      </label>
    </div>
  </div>
  <div class="nasc-grid">
    <figure class="nasc-cell">
      <canvas data-nasc="net" role="img" aria-label="Editable four-unit network with contamination shading"></canvas>
      <figcaption>Click an edge to connect or disconnect two units. Node shading is absolute contamination.</figcaption>
    </figure>
    <figure class="nasc-cell">
      <canvas data-nasc="chart" role="img" aria-label="Synthetic control outcome trajectories before and after treatment"></canvas>
      <figcaption>
        <span class="nasc-key nasc-key-obs"></span>treated
        <span class="nasc-key nasc-key-truth"></span>true counterfactual
        <span class="nasc-key nasc-key-sc"></span>SCM
        <span class="nasc-key nasc-key-nasc"></span>NASC
      </figcaption>
    </figure>
  </div>
  <div class="nasc-grid">
    <figure class="nasc-cell">
      <canvas data-nasc="simplex" role="img" aria-label="Exponentially tilted Dirichlet density over the donor weight simplex"></canvas>
      <div class="nasc-inline-slider">
        <label>
          <span class="nasc-label">Prior α</span>
          <input type="range" min="10" max="30" step="1" value="10" data-nasc="alpha" />
          <output data-nasc="alphaVal">1.0</output>
        </label>
      </div>
      <figcaption>Tilted Dirichlet density over the weight simplex. Drag to rotate; the dashed path traces the weights as λ grows. At α = 1 the base measure is uniform, so the surface shows the tilt alone.</figcaption>
    </figure>
    <div class="nasc-cell nasc-readout" data-nasc="readout"></div>
  </div>
</div>

<noscript>
  <p class="nasc-fallback">
    The interactive explorer needs JavaScript. The static figures above show the same
    objects for the thesis example networks.
  </p>
</noscript>

<script src="{{ '/assets/js/nasc-explorer.js' | relative_url | bust_file_cache }}"></script>


The method is implemented in the R package nasc, which extends the estimator to a
spatial Durbin random-effects panel specification. NASC is open-source software. Visit
the repository for installation instructions, examples, and the latest development
version.

<div class="nasc-repository-card repositories d-flex flex-wrap flex-md-row flex-column justify-content-center">
  {% include repository/repo.liquid repository="fvwaldow/nasc" %}
</div>
