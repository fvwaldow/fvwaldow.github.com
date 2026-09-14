---
layout: page
permalink: /repositories/
title: Repositories
description:
nav: true
nav_order: 5
_styles: |
  .post-header {
    display: none;
  }

  .repo img {
    mix-blend-mode: multiply;
  }

  html[data-theme="dark"] .repo img {
    filter: invert(1);
    mix-blend-mode: screen;
  }
---

<div class="research-section">
  <h2>GitHub Repositories</h2>
  <div class="publications repositories d-flex flex-wrap flex-md-row flex-column justify-content-center">
    {% for repository in site.data.repositories.github_repos %}
      {% include repository/repo.liquid repository=repository %}
    {% endfor %}
  </div>
</div>

