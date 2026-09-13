---
layout: post
title: Network-Aware Synthetic Control (NASC)
date: 2026-09-13
description: NASC is an R package for Network-Aware Synthetic Control estimation when SUTVA is violated as treatment spillovers travel through networks. It corrects donor-pool contamination and recovers implied spillovers onto untreated units.
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

  @media (max-width: 767px) {
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

<div class="nasc-figure-bottom-right">
  {% include figure.liquid path="assets/img/nasc multiple solutions.png" title="nasc multiple solutions" class="img-fluid rounded z-depth-1" %}
</div>

Inference is modular Bayesian, based on a cut posterior that preserves uncertainty about
the estimated bias correction, with an exponentially tilted Dirichlet prior encoding the
contamination penalty on the simplex. As a by-product, the estimator recovers the
implied spillovers onto untreated units. Monte Carlo simulations across regular and
small-world topologies show that the proposed estimator outperforms the conventional and
the Bayesian synthetic control across various interference structures.

The method is implemented in the R package nasc, which extends the estimator to a
spatial Durbin random-effects panel specification. NASC is open-source software. Visit
the repository for installation instructions, examples, and the latest development
version.


<div class="nasc-repository-card repositories d-flex flex-wrap flex-md-row flex-column justify-content-center">
  {% include repository/repo.liquid repository="fvwaldow/nasc" %}
</div>
