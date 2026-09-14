---
layout: page
permalink: /repositories/
title: Repositories
description: Open-source code and research software by Frederik von Waldow, including statistical packages for econometrics and causal inference.
nav: true
nav_order: 5
_styles: |
  .post-header {
    display: none;
  }

  /* cards are requested with bg_color=00000000, so they inherit the page background
     in both themes — no blend-mode or inversion needed */
---

<div class="research-section">
  <h2>GitHub Repositories</h2>
  <div class="publications repositories d-flex flex-wrap flex-md-row flex-column justify-content-center">
    {% for repository in site.data.repositories.github_repos %}
      {% include repository/repo.liquid repository=repository %}
    {% endfor %}
  </div>
</div>

