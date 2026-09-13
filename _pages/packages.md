---
layout: default
title: Packages
permalink: /packages/
nav: true
nav_order: 4
---

<div class="post">
  <ul class="post-list">
    {% assign package_posts = site.posts | where_exp: "post", "post.package == true" | sort: "date" | reverse %}
    {% for post in package_posts %}
      <li>
        <h3><a class="post-title" href="{{ post.url | relative_url }}">{{ post.title }}</a></h3>
        <p>{{ post.description }}</p>
        <p class="post-meta">{{ post.date | date: '%B %d, %Y' }}</p>
      </li>
    {% endfor %}
  </ul>
</div>
