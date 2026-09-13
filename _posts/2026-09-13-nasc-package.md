---
layout: post
title: NASC
date: 2026-09-13
description: NASC is an R package for Network-Aware Synthetic Control estimation when spillovers travel through networks. It corrects donor-pool contamination and recovers implied spillovers onto untreated units.
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
    max-width: 82%;
    margin: 0 auto;
  }

  .nasc-top-figure img,
  .nasc-figure-bottom-right img {
    mix-blend-mode: multiply;
    box-shadow: none !important;
  }

  .nasc-figure-bottom-right figure {
    margin: 0;
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

Synthetic control methods assume that a policy leaves the untreated comparison units unaffected. In spatial settings — transport investments, place-based policies, regional taxes — this assumption often fails, because effects travel through trade, commuting, or competition networks to exactly those donor units used to construct the counterfactual. I develop the Network-Aware Synthetic Control (NASC) estimator, which models this propagation as a spatial autoregressive process on a known network and corrects the treatment effect for contamination of the donor pool.

<div class="nasc-figure-bottom-right">
  {% include figure.liquid path="assets/img/nasc multiple solutions.png" title="nasc multiple solutions" class="img-fluid rounded z-depth-1" %}
</div>

Inference is Bayesian, based on a modular (cut) posterior with an exponentially tilted Dirichlet prior on the weights. As a by-product, the estimator recovers the implied spillovers onto untreated units. Monte Carlo simulations on network data-generating processes compare NASC with conventional and Bayesian synthetic control and show that its performance depends critically on the quality of the pre-treatment fit. The method is implemented in the R package nasc.

NASC is open-source software. Visit the repository for installation instructions, examples, and the latest development version.

<a href="https://github.com/fvwaldow/nasc" class="btn btn-sm z-depth-0" role="button" target="_blank" rel="noopener noreferrer">View NASC on GitHub</a>
