---
layout: default
title: Packages
description: Statistical software packages by Frederik von Waldow, including NASC for network-aware synthetic control estimation under interference.
permalink: /packages/
nav: true
nav_order: 4
---

<style>
  .package-entry {
    display: grid;
    grid-template-columns: minmax(0, 1fr) minmax(140px, 22%);
    gap: 1.5rem;
    align-items: center;
  }

  .package-entry figure {
    max-width: 180px;
    margin: 0;
  }

  .package-entry > a {
    justify-self: end;
    width: 100%;
  }

  .package-entry img {
    mix-blend-mode: multiply;
    box-shadow: none !important;
  }

  html[data-theme="dark"] .package-entry img {
    filter: invert(1);
    mix-blend-mode: screen;
  }

  @media (max-width: 767px) {
    .package-entry {
      display: block;
    }

    .package-entry figure {
      margin-top: 1rem;
    }
  }
</style>

<div class="post">
  {% assign package_posts = site.posts | where_exp: "post", "post.package == true" | sort: "date" | reverse %}
  {% for post in package_posts %}
    <div class="research-section">
      <h2><a href="{{ post.url | relative_url }}">{{ post.title }}</a></h2>
      <div class="publications package-entry">
        <div>
          <p>{{ post.description }}</p>
        </div>
        <a href="{{ post.url | relative_url }}" aria-label="Open {{ post.title }} package page">
          {% include figure.liquid path="assets/img/nasc plot.png" title="NASC estimator" class="img-fluid rounded z-depth-1" %}
        </a>
      </div>
    </div>
  {% endfor %}
</div>
