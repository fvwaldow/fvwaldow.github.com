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
---

## GitHub Repositories

<div class="repositories d-flex flex-wrap flex-md-row flex-column justify-content-center">
  {% for repository in site.data.repositories.github_repos %}
    {% include repository/repo.liquid repository=repository %}
  {% endfor %}
</div>

