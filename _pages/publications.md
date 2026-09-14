---
layout: page
permalink: /publications/
title: Research
description: Publications and working papers by Frederik von Waldow on spatial competition, fuel taxes, transport policy, and causal inference.
nav: true
nav_order: 2
_styles: |
  .post-header {
    display: none;
  }
  #bibsearch,
  .bibsearch-form-input {
    display: none;
  }
---

<!-- _pages/publications.md -->

<div class="research-section">
  <h2>Publications</h2>
  <div class="publications">
  {% bibliography --query @article %}
  </div>
</div>

<div class="research-section">
  <h2>Working Papers</h2>
  <div class="publications">
  {% bibliography --query @techreport %}
  </div>
</div>
