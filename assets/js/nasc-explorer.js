/**
 * NASC explorer — an interactive companion to the Network-Aware Synthetic Control post.
 *
 * Three linked panels, all driven by the same state:
 *   1. a 4-node network (1 treated unit, 3 donors) whose edges the reader toggles;
 *   2. the synthetic control trajectories implied by that network and by rho;
 *   3. the exponentially tilted Dirichlet density over the simplex of donor weights.
 *
 * Self-contained: no plotting library, no external dependency.
 */
(function () {
  "use strict";

  var ROOT_ID = "nasc-explorer";
  var N_UNITS = 4; // unit 0 is treated, units 1..3 are donors
  var J = N_UNITS - 1;

  /* ------------------------------------------------------------------ *
   * Linear algebra (tiny dense matrices)
   * ------------------------------------------------------------------ */

  function zeros(n, m) {
    var out = [];
    for (var i = 0; i < n; i++) {
      var row = [];
      for (var j = 0; j < m; j++) row.push(0);
      out.push(row);
    }
    return out;
  }

  /** Gauss-Jordan inverse with partial pivoting. Returns null if singular. */
  function matInverse(A) {
    var n = A.length;
    var M = [];
    for (var i = 0; i < n; i++) {
      M.push(A[i].slice());
      for (var j = 0; j < n; j++) M[i].push(i === j ? 1 : 0);
    }
    for (var col = 0; col < n; col++) {
      var piv = col;
      for (var r = col + 1; r < n; r++) {
        if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
      }
      if (Math.abs(M[piv][col]) < 1e-12) return null;
      var tmp = M[col];
      M[col] = M[piv];
      M[piv] = tmp;
      var d = M[col][col];
      for (var c = 0; c < 2 * n; c++) M[col][c] /= d;
      for (var r2 = 0; r2 < n; r2++) {
        if (r2 === col) continue;
        var f = M[r2][col];
        if (f === 0) continue;
        for (var c2 = 0; c2 < 2 * n; c2++) M[r2][c2] -= f * M[col][c2];
      }
    }
    var inv = [];
    for (var k = 0; k < n; k++) inv.push(M[k].slice(n));
    return inv;
  }

  function matVec(A, v) {
    var out = [];
    for (var i = 0; i < A.length; i++) {
      var s = 0;
      for (var j = 0; j < v.length; j++) s += A[i][j] * v[j];
      out.push(s);
    }
    return out;
  }

  function dot(a, b) {
    var s = 0;
    for (var i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  /* ------------------------------------------------------------------ *
   * Contamination vector  s = rho (I - rho W_J)^{-1} w_J1
   * ------------------------------------------------------------------ */

  /** Row-standardize an adjacency matrix; all-zero rows stay all-zero. */
  function rowStandardize(adj) {
    var n = adj.length;
    var W = zeros(n, n);
    for (var i = 0; i < n; i++) {
      var sum = 0;
      for (var j = 0; j < n; j++) sum += adj[i][j];
      if (sum <= 0) continue;
      for (var k = 0; k < n; k++) W[i][k] = adj[i][k] / sum;
    }
    return W;
  }

  /**
   * Splits the row-standardized W into the treated row/column blocks and
   * returns the spatial multiplier of the donor block together with the
   * contamination vector s = rho (I_J - rho W_J)^{-1} w_J1 of Equation (6).
   */
  function spatialParts(adj, rho) {
    var W = rowStandardize(adj);
    var WJ = zeros(J, J);
    var wJ1 = [];
    for (var i = 0; i < J; i++) {
      wJ1.push(W[i + 1][0]);
      for (var j = 0; j < J; j++) WJ[i][j] = W[i + 1][j + 1];
    }
    var M = zeros(J, J);
    for (var a = 0; a < J; a++) {
      for (var b = 0; b < J; b++) M[a][b] = (a === b ? 1 : 0) - rho * WJ[a][b];
    }
    var Minv = matInverse(M);
    if (!Minv) {
      Minv = zeros(J, J);
      for (var d = 0; d < J; d++) Minv[d][d] = 1;
    }
    var s = matVec(Minv, wJ1).map(function (x) {
      return rho * x;
    });
    return { wJ1: wJ1, WJ: WJ, Minv: Minv, s: s };
  }

  function contamination(adj, rho) {
    return spatialParts(adj, rho).s;
  }

  /* ------------------------------------------------------------------ *
   * Synthetic panel data (seeded, so sliders do not reshuffle the data)
   * ------------------------------------------------------------------ */

  function mulberry32(a) {
    return function () {
      a |= 0;
      a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function normal(rnd) {
    var u = 1 - rnd();
    var v = rnd();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  var T0 = 24; // pre-treatment periods
  var T1 = 12; // post-treatment periods
  var T = T0 + T1;

  /* Monte Carlo design of Section 8.2, scaled down to three donors. */
  var BETA = [1, 0.5]; // covariate coefficients
  var SIGMA = 0.5; // sd of the idiosyncratic errors
  var GAMMA_TRUE = [0.2, 0.5, 0.3]; // true donor weights on donors 2, 3, 4

  /**
   * Structural shocks, drawn once per seed and held fixed while the sliders
   * move: unit fixed effects mu_i ~ N(0,1), two standard normal covariates
   * with beta = (1, 0.5)', errors eps_it ~ N(0, 0.5^2), and a treatment effect
   * tau_t ~ N(5, 1) switched on at T0. The treated unit's untreated potential
   * outcome is the convex combination of the donors' systematic components,
   * which makes GAMMA_TRUE the target the synthetic control is chasing.
   *
   * GAMMA_TRUE puts its largest mass on donor 3, the donor adjacent to the
   * treated unit under the default network. As in the thesis, the donors that
   * help most with pre-treatment fit are also the most contaminated ones.
   */
  function buildShocks(seed) {
    var rnd = mulberry32(seed);
    var systematic = []; // mu_i + Z_it' beta
    for (var i = 0; i < J; i++) {
      var mu = normal(rnd);
      var row = [];
      for (var t = 0; t < T; t++) {
        row.push(mu + BETA[0] * normal(rnd) + BETA[1] * normal(rnd));
      }
      systematic.push(row);
    }
    var eps = [];
    for (var i2 = 0; i2 < J; i2++) {
      var erow = [];
      for (var t2 = 0; t2 < T; t2++) erow.push(SIGMA * normal(rnd));
      eps.push(erow);
    }
    var tau = [];
    var y1Untreated = [];
    for (var t3 = 0; t3 < T; t3++) {
      tau.push(t3 < T0 ? 0 : 5 + normal(rnd));
      var c = 0;
      for (var k = 0; k < J; k++) c += GAMMA_TRUE[k] * systematic[k][t3];
      y1Untreated.push(c);
    }
    return {
      gammaTrue: GAMMA_TRUE,
      systematic: systematic,
      eps: eps,
      tau: tau,
      y1Untreated: y1Untreated,
    };
  }

  /**
   * Realized panel for one network and one value of rho. Donor untreated
   * potential outcomes solve the spatial autoregression
   *
   *   y^d_t(0) = rho w_J1 y_1t(0) + rho W_J y^d_t(0) + mu + Z_t' beta + eps_t,
   *
   * so donors are spatially filtered and already co-move with the treated unit
   * before any treatment. Observed outcomes then add the treatment itself:
   * tau_t for the treated unit and s_i tau_t for donor i, which is exactly the
   * indirect effect delta_t = s tau_1t of Equation (7).
   */
  function buildPanel(shocks, adj, rho) {
    var sp = spatialParts(adj, rho);
    var s = sp.s;

    var y0 = [];
    for (var i = 0; i < J; i++) y0.push([]);
    for (var t = 0; t < T; t++) {
      var rhs = [];
      for (var i2 = 0; i2 < J; i2++) {
        rhs.push(rho * sp.wJ1[i2] * shocks.y1Untreated[t] + shocks.systematic[i2][t] + shocks.eps[i2][t]);
      }
      var col = matVec(sp.Minv, rhs);
      for (var i3 = 0; i3 < J; i3++) y0[i3].push(col[i3]);
    }

    var yd = y0.map(function (path, i4) {
      return path.map(function (v, t2) {
        return v + s[i4] * shocks.tau[t2];
      });
    });
    var y1 = shocks.y1Untreated.map(function (v, t3) {
      return v + shocks.tau[t3];
    });

    return {
      s: s,
      y0: y0,
      yd: yd,
      y1: y1,
      y1Untreated: shocks.y1Untreated,
      tau: shocks.tau,
      gammaTrue: shocks.gammaTrue,
    };
  }

  /* ------------------------------------------------------------------ *
   * Weight estimation on the simplex
   * ------------------------------------------------------------------ */

  /** Euclidean projection onto the probability simplex. */
  function projectSimplex(v) {
    var n = v.length;
    var u = v.slice().sort(function (a, b) {
      return b - a;
    });
    var css = 0;
    var rho = 0;
    var theta = 0;
    for (var i = 0; i < n; i++) {
      css += u[i];
      var t = (css - 1) / (i + 1);
      if (u[i] - t > 0) {
        rho = i + 1;
        theta = t;
      }
    }
    return v.map(function (x) {
      return Math.max(x - theta, 0);
    });
  }

  /**
   * argmin over the simplex of  ||X1 - X0 g||^2 + lambda * g'|s|
   * by projected gradient descent. X0 is T0-by-J, X1 is length T0.
   */
  function solveWeights(X0, X1, absS, lambda) {
    var n = X0.length;
    var G = zeros(J, J); // X0' X0
    var h = []; // X0' X1
    for (var a = 0; a < J; a++) {
      var ha = 0;
      for (var t = 0; t < n; t++) ha += X0[t][a] * X1[t];
      h.push(ha);
      for (var b = 0; b < J; b++) {
        var g = 0;
        for (var t2 = 0; t2 < n; t2++) g += X0[t2][a] * X0[t2][b];
        G[a][b] = g;
      }
    }
    // Step size from the top eigenvalue of G, found by power iteration. Using the
    // trace instead would be valid but several times more conservative.
    var pv = [];
    for (var z = 0; z < J; z++) pv.push(1);
    var top = 0;
    for (var it0 = 0; it0 < 60; it0++) {
      var nv = matVec(G, pv);
      var norm = Math.sqrt(dot(nv, nv));
      if (norm < 1e-12) break;
      top = norm;
      pv = nv.map(function (x) {
        return x / norm;
      });
    }
    var step = 1 / (2 * Math.max(top, 1e-6));

    var gamma = [];
    for (var i = 0; i < J; i++) gamma.push(1 / J);
    for (var it = 0; it < 4000; it++) {
      var grad = [];
      for (var p = 0; p < J; p++) {
        var gv = 0;
        for (var q = 0; q < J; q++) gv += G[p][q] * gamma[q];
        grad.push(2 * (gv - h[p]) + lambda * absS[p]);
      }
      var next = [];
      for (var m = 0; m < J; m++) next.push(gamma[m] - step * grad[m]);
      gamma = projectSimplex(next);
    }
    return gamma;
  }

  /** Everything the three panels need, for one setting of the controls. */
  function computeModel(shocks, adj, rho, lambda) {
    var panel = buildPanel(shocks, adj, rho);
    var s = panel.s;
    var absS = s.map(Math.abs);
    var obs = { yd: panel.yd, y1: panel.y1 };

    var X0 = [];
    var X1 = [];
    for (var t = 0; t < T0; t++) {
      X0.push([obs.yd[0][t], obs.yd[1][t], obs.yd[2][t]]);
      X1.push(obs.y1[t]);
    }

    var gSc = solveWeights(X0, X1, absS, 0);
    var gNasc = solveWeights(X0, X1, absS, lambda);

    var scFactor = 1 - dot(gSc, s);
    var nascFactor = 1 - dot(gNasc, s);

    var scPath = []; // conventional synthetic control outcome
    var nascPath = []; // bias-corrected counterfactual under penalized weights
    var attTrue = 0;
    var attSc = 0;
    var attBc = 0;
    var attCr = 0;
    var attNasc = 0;

    for (var t2 = 0; t2 < T; t2++) {
      var ySc = 0;
      var yCr = 0;
      for (var i = 0; i < J; i++) {
        ySc += gSc[i] * obs.yd[i][t2];
        yCr += gNasc[i] * obs.yd[i][t2];
      }
      scPath.push(ySc);
      var tauSc = obs.y1[t2] - ySc;
      var tauCr = obs.y1[t2] - yCr;
      var tauNasc = nascFactor !== 0 ? tauCr / nascFactor : tauCr;
      // Before treatment there is nothing to correct, so the penalized synthetic
      // control is just the weighted donor average.
      nascPath.push(t2 < T0 ? yCr : obs.y1[t2] - tauNasc);
      if (t2 >= T0) {
        attTrue += panel.tau[t2];
        attSc += tauSc;
        attBc += scFactor !== 0 ? tauSc / scFactor : tauSc;
        attCr += tauCr;
        attNasc += tauNasc;
      }
    }

    return {
      panel: panel,
      s: s,
      absS: absS,
      obs: obs,
      gammaTrue: panel.gammaTrue,
      gSc: gSc,
      gNasc: gNasc,
      scPath: scPath,
      nascPath: nascPath,
      scFactor: scFactor,
      nascFactor: nascFactor,
      tilt: lambda / (2 * SIGMA * SIGMA * T0), // prior-scale tilt, Section 9.1
      contamSc: dot(gSc, absS),
      contamNasc: dot(gNasc, absS),
      att: {
        truth: attTrue / T1,
        sc: attSc / T1,
        bc: attBc / T1,
        cr: attCr / T1,
        nasc: attNasc / T1,
      },
    };
  }

  /* ------------------------------------------------------------------ *
   * Drawing helpers
   * ------------------------------------------------------------------ */

  var SERIES = {
    sc: "#e07a3f", // orange, as in the thesis figures
    nasc: "#3f7fb5", // blue, as in the thesis figures
  };

  function theme() {
    var cs = getComputedStyle(document.documentElement);
    var dark = document.documentElement.getAttribute("data-theme") === "dark";
    function v(name, fallback) {
      var out = (cs.getPropertyValue(name) || "").trim();
      return out || fallback;
    }
    return {
      dark: dark,
      text: v("--global-text-color", dark ? "#e8e8e8" : "#1c1c1d"),
      muted: v("--global-text-color-light", "#828282"),
      divider: v("--global-divider-color", dark ? "#424246" : "rgba(0,0,0,0.1)"),
      theme: v("--global-theme-color", "#b509ac"),
      bg: v("--global-bg-color", dark ? "#1c1c1d" : "#faf9f6"),
    };
  }

  /** Size a canvas for the device pixel ratio and return its 2d context. */
  function prepare(canvas, cssHeight) {
    var ratio = window.devicePixelRatio || 1;
    var width = canvas.parentNode.clientWidth;
    if (!width) width = 320;
    canvas.style.width = "100%";
    canvas.style.height = cssHeight + "px";
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(cssHeight * ratio);
    var ctx = canvas.getContext("2d");
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.clearRect(0, 0, width, cssHeight);
    return { ctx: ctx, w: width, h: cssHeight };
  }

  function lerp(a, b, u) {
    return a + (b - a) * u;
  }

  /** Yellow to dark red ramp, echoing the thesis ternary plots. */
  function heatColor(u) {
    var stops = [
      [255, 247, 188],
      [254, 153, 41],
      [153, 52, 4],
    ];
    u = Math.max(0, Math.min(1, u));
    var scaled = u * (stops.length - 1);
    var i = Math.min(Math.floor(scaled), stops.length - 2);
    var f = scaled - i;
    var c0 = stops[i];
    var c1 = stops[i + 1];
    return "rgb(" + Math.round(lerp(c0[0], c1[0], f)) + "," + Math.round(lerp(c0[1], c1[1], f)) + "," + Math.round(lerp(c0[2], c1[2], f)) + ")";
  }

  /* ------------------------------------------------------------------ *
   * Panel 1 — the network
   * ------------------------------------------------------------------ */

  var NODE_POS = [
    { x: 0.22, y: 0.24, label: "1" },
    { x: 0.78, y: 0.24, label: "2" },
    { x: 0.22, y: 0.78, label: "3" },
    { x: 0.78, y: 0.78, label: "4" },
  ];

  var EDGES = [];
  for (var ei = 0; ei < N_UNITS; ei++) {
    for (var ej = ei + 1; ej < N_UNITS; ej++) EDGES.push([ei, ej]);
  }

  function drawNetwork(canvas, state, model) {
    var p = prepare(canvas, 230);
    var ctx = p.ctx;
    var th = theme();
    var pad = 30;
    var size = Math.min(p.w, p.h) - 2 * pad;
    var ox = (p.w - size) / 2;
    var oy = (p.h - size) / 2;

    function pos(i) {
      return { x: ox + NODE_POS[i].x * size, y: oy + NODE_POS[i].y * size };
    }

    var maxAbs = Math.max.apply(null, model.absS.concat([1e-6]));

    // Edges are the click targets themselves: present ones solid, absent ones
    // dashed affordances, and whichever one is under the cursor picked out in
    // the accent colour.
    EDGES.forEach(function (e, k) {
      var a = pos(e[0]);
      var b = pos(e[1]);
      var on = state.adj[e[0]][e[1]] === 1;
      var hot = state.hoverEdge === k;
      ctx.beginPath();
      ctx.moveTo(a.x, a.y);
      ctx.lineTo(b.x, b.y);
      ctx.lineWidth = hot ? 3 : on ? 1.6 : 1;
      ctx.setLineDash(on ? [] : [3, 4]);
      ctx.strokeStyle = hot ? th.theme : on ? th.text : th.divider;
      ctx.stroke();
      ctx.setLineDash([]);
    });

    // nodes
    for (var i = 0; i < N_UNITS; i++) {
      var c = pos(i);
      ctx.beginPath();
      ctx.arc(c.x, c.y, 17, 0, 2 * Math.PI);
      if (i === 0) {
        ctx.fillStyle = th.text;
      } else {
        var intensity = model.absS[i - 1] / maxAbs;
        ctx.fillStyle = "rgba(209,73,91," + (0.12 + 0.8 * intensity).toFixed(3) + ")";
      }
      ctx.fill();
      ctx.lineWidth = 1.2;
      ctx.strokeStyle = th.divider;
      ctx.stroke();

      ctx.fillStyle = i === 0 ? th.bg : th.text;
      ctx.font = "600 12px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(NODE_POS[i].label, c.x, c.y);

      ctx.fillStyle = th.muted;
      ctx.font = "10px -apple-system, BlinkMacSystemFont, sans-serif";
      if (i === 0) {
        ctx.fillText("treated", c.x, c.y - 28);
      } else {
        ctx.fillText("s = " + model.s[i - 1].toFixed(3), c.x, c.y + (NODE_POS[i].y > 0.5 ? 29 : -28));
      }
    }
    canvas._hit = { pos: pos };
  }

  var EDGE_TOL = 10; // px from the line that still counts as a click on it
  var NODE_GUARD = 19; // px around a node centre (radius 17) where edges tie
  var TIE_MARGIN = 3; // px by which the nearest edge must beat the runner-up

  /** Distance from (px, py) to the segment ab. */
  function segDist(px, py, a, b) {
    var dx = b.x - a.x;
    var dy = b.y - a.y;
    var len2 = dx * dx + dy * dy;
    var t = len2 > 0 ? ((px - a.x) * dx + (py - a.y) * dy) / len2 : 0;
    t = Math.max(0, Math.min(1, t));
    return Math.hypot(px - (a.x + t * dx), py - (a.y + t * dy));
  }

  /**
   * Index into EDGES of the edge under the pointer, or -1. Two neighbourhoods
   * are genuinely ambiguous and are excluded rather than resolved arbitrarily:
   * a node, where every incident edge is under the cursor, and the point where
   * the two diagonals cross. Everywhere else the nearest edge wins.
   */
  function nearestEdge(canvas, x, y) {
    if (!canvas._hit) return -1;
    for (var n = 0; n < N_UNITS; n++) {
      var c = canvas._hit.pos(n);
      if (Math.hypot(x - c.x, y - c.y) < NODE_GUARD) return -1;
    }
    var best = -1;
    var bestD = Infinity;
    var runnerUp = Infinity;
    for (var k = 0; k < EDGES.length; k++) {
      var a = canvas._hit.pos(EDGES[k][0]);
      var b = canvas._hit.pos(EDGES[k][1]);
      var d = segDist(x, y, a, b);
      if (d < bestD) {
        runnerUp = bestD;
        bestD = d;
        best = k;
      } else if (d < runnerUp) {
        runnerUp = d;
      }
    }
    if (bestD > EDGE_TOL) return -1;
    return runnerUp - bestD < TIE_MARGIN ? -1 : best;
  }

  /** Canvas coordinates of a pointer or touch event. */
  function canvasPoint(canvas, evt) {
    var rect = canvas.getBoundingClientRect();
    return { x: evt.clientX - rect.left, y: evt.clientY - rect.top };
  }

  function networkHit(canvas, state, evt) {
    var pt = canvasPoint(canvas, evt);
    var k = nearestEdge(canvas, pt.x, pt.y);
    if (k < 0) return false;
    var e = EDGES[k];
    var cur = state.adj[e[0]][e[1]];
    state.adj[e[0]][e[1]] = cur ? 0 : 1;
    state.adj[e[1]][e[0]] = cur ? 0 : 1;
    return true;
  }

  /* ------------------------------------------------------------------ *
   * Panel 2 — synthetic control trajectories
   * ------------------------------------------------------------------ */

  function drawChart(canvas, panel, model) {
    var p = prepare(canvas, 260);
    var ctx = p.ctx;
    var th = theme();
    var padL = 38;
    var padR = 10;
    var padT = 14;
    var padB = 30;
    var plotW = p.w - padL - padR;
    var plotH = p.h - padT - padB;

    var series = [
      { data: model.obs.y1, color: th.text, width: 1.8 },
      { data: panel.y1Untreated, color: th.muted, width: 1.4, dash: [4, 3] },
      { data: model.scPath, color: SERIES.sc, width: 1.6 },
      { data: model.nascPath, color: SERIES.nasc, width: 1.6 },
    ];

    var lo = Infinity;
    var hi = -Infinity;
    series.forEach(function (s) {
      s.data.forEach(function (v) {
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      });
    });
    var span = hi - lo || 1;
    lo -= span * 0.12;
    hi += span * 0.12;

    function X(t) {
      return padL + (t / (T - 1)) * plotW;
    }
    function Y(v) {
      return padT + (1 - (v - lo) / (hi - lo)) * plotH;
    }

    // horizontal gridlines
    ctx.strokeStyle = th.divider;
    ctx.lineWidth = 1;
    ctx.font = "10px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = th.muted;
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (var g = 0; g <= 4; g++) {
      var val = lo + ((hi - lo) * g) / 4;
      var yy = Y(val);
      ctx.beginPath();
      ctx.moveTo(padL, yy);
      ctx.lineTo(padL + plotW, yy);
      ctx.stroke();
      ctx.fillText(val.toFixed(1), padL - 6, yy);
    }

    // treatment date
    ctx.save();
    ctx.setLineDash([4, 4]);
    ctx.strokeStyle = th.muted;
    ctx.beginPath();
    ctx.moveTo(X(T0 - 0.5), padT);
    ctx.lineTo(X(T0 - 0.5), padT + plotH);
    ctx.stroke();
    ctx.restore();
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillStyle = th.muted;
    ctx.fillText("treatment", X(T0 - 0.5), padT + plotH + 6);
    ctx.textAlign = "left";
    ctx.fillText("t = 0", padL, padT + plotH + 6);
    ctx.textAlign = "right";
    ctx.fillText("t = " + (T - 1), padL + plotW, padT + plotH + 6);

    series.forEach(function (s) {
      ctx.beginPath();
      ctx.lineWidth = s.width;
      ctx.strokeStyle = s.color;
      ctx.setLineDash(s.dash || []);
      s.data.forEach(function (v, t) {
        if (t === 0) ctx.moveTo(X(t), Y(v));
        else ctx.lineTo(X(t), Y(v));
      });
      ctx.stroke();
      ctx.setLineDash([]);
    });
  }

  /* ------------------------------------------------------------------ *
   * Panel 3 — the exponentially tilted Dirichlet over the simplex
   * ------------------------------------------------------------------ */

  var TRI = [
    { x: -0.5, y: -Math.sqrt(3) / 6 },
    { x: 0.5, y: -Math.sqrt(3) / 6 },
    { x: 0.0, y: Math.sqrt(3) / 3 },
  ];

  function baryToPlane(g) {
    return {
      x: g[0] * TRI[0].x + g[1] * TRI[1].x + g[2] * TRI[2].x,
      y: g[0] * TRI[0].y + g[1] * TRI[1].y + g[2] * TRI[2].y,
    };
  }

  /** Unnormalized ETDir density: prod g_i^(alpha-1) * exp(-lambda g'|s|). */
  function etdirDensity(g, alpha, lambda, absS) {
    var v = Math.exp(-lambda * (g[0] * absS[0] + g[1] * absS[1] + g[2] * absS[2]));
    for (var i = 0; i < 3; i++) {
      var gi = Math.max(g[i], 1e-9);
      v *= Math.pow(gi, alpha - 1);
    }
    return v;
  }

  function drawSimplex(canvas, state, model) {
    var p = prepare(canvas, 300);
    var ctx = p.ctx;
    var th = theme();
    var alpha = state.alpha;
    // The penalty lambda multiplies a sum of squares over T0 periods, whereas
    // the prior tilt multiplies an O(1) log density. Section 9.1 of the thesis
    // bridges the two with a sigma^-2 (T0 + sum v_m) factor; here the matching
    // per-period tilt is lambda / (2 sigma^2 T0), which keeps the surface in
    // the regime drawn in Figure 2.
    var lambda = model.tilt;
    var absS = model.absS;

    var N = 72; // barycentric mesh resolution
    var dens = [];
    var maxD = 0;
    for (var i = 0; i <= N; i++) {
      dens.push([]);
      for (var j = 0; j <= N - i; j++) {
        var g = [i / N, j / N, 1 - i / N - j / N];
        var d = etdirDensity(g, alpha, lambda, absS);
        if (isFinite(d) && d > maxD) maxD = d;
        dens[i].push(isFinite(d) ? d : 0);
      }
    }
    if (maxD <= 0) maxD = 1;

    var az = state.az;
    var el = state.el;
    var ca = Math.cos(az);
    var sa = Math.sin(az);
    var ce = Math.cos(el);
    var se = Math.sin(el);
    var scale = Math.min(p.w, p.h) * 0.72;
    var cx = p.w / 2;
    var cy = p.h / 2 + 24;
    var height = 0.46;

    function project(x, y, z) {
      var x1 = x * ca - y * sa;
      var y1 = x * sa + y * ca;
      return {
        x: cx + x1 * scale,
        y: cy - (y1 * se + z * ce) * scale,
        depth: -y1 * ce + z * se,
      };
    }

    // Project every mesh node once; each interior node is shared by six facets.
    var verts = [];
    for (var vi = 0; vi <= N; vi++) {
      verts.push([]);
      for (var vj = 0; vj <= N - vi; vj++) {
        var a = vi / N;
        var b = vj / N;
        var pl = baryToPlane([a, b, 1 - a - b]);
        var u = dens[vi][vj] / maxD;
        var pr = project(pl.x, pl.y, u * height);
        pr.u = u;
        verts[vi].push(pr);
      }
    }

    // build and depth-sort the facets
    var facets = [];
    for (var a2 = 0; a2 < N; a2++) {
      for (var b2 = 0; b2 < N - a2; b2++) {
        var v1 = verts[a2][b2];
        var v2 = verts[a2 + 1][b2];
        var v3 = verts[a2][b2 + 1];
        facets.push({
          pts: [v1, v2, v3],
          depth: (v1.depth + v2.depth + v3.depth) / 3,
          u: (v1.u + v2.u + v3.u) / 3,
        });
        if (b2 < N - a2 - 1) {
          var v4 = verts[a2 + 1][b2 + 1];
          facets.push({
            pts: [v2, v4, v3],
            depth: (v2.depth + v4.depth + v3.depth) / 3,
            u: (v2.u + v4.u + v3.u) / 3,
          });
        }
      }
    }
    facets.sort(function (f1, f2) {
      return f1.depth - f2.depth;
    });

    // Constant for every facet; the hairline stroke closes antialiasing seams.
    ctx.lineWidth = 0.4;
    ctx.lineJoin = "round";
    facets.forEach(function (f) {
      ctx.beginPath();
      ctx.moveTo(f.pts[0].x, f.pts[0].y);
      ctx.lineTo(f.pts[1].x, f.pts[1].y);
      ctx.lineTo(f.pts[2].x, f.pts[2].y);
      ctx.closePath();
      ctx.fillStyle = heatColor(f.u);
      ctx.fill();
      ctx.strokeStyle = ctx.fillStyle;
      ctx.stroke();
    });

    // base triangle
    ctx.beginPath();
    for (var c = 0; c < 3; c++) {
      var pr = project(TRI[c].x, TRI[c].y, 0);
      if (c === 0) ctx.moveTo(pr.x, pr.y);
      else ctx.lineTo(pr.x, pr.y);
    }
    ctx.closePath();
    ctx.lineWidth = 1.1;
    ctx.strokeStyle = th.muted;
    ctx.stroke();

    // trajectory of the penalized weights as lambda grows
    if (state.trajectory && state.trajectory.length > 1) {
      ctx.beginPath();
      state.trajectory.forEach(function (g, idx) {
        var pl2 = baryToPlane(g);
        var pr2 = project(pl2.x, pl2.y, 0);
        if (idx === 0) ctx.moveTo(pr2.x, pr2.y);
        else ctx.lineTo(pr2.x, pr2.y);
      });
      ctx.lineWidth = 1.6;
      ctx.setLineDash([3, 3]);
      ctx.strokeStyle = SERIES.nasc;
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // current weight vector: stem from the base up to the density surface
    var pl3 = baryToPlane(model.gNasc);
    var dCur = etdirDensity(model.gNasc, alpha, lambda, absS) / maxD;
    var base = project(pl3.x, pl3.y, 0);
    var top = project(pl3.x, pl3.y, Math.min(dCur, 1) * height);
    ctx.beginPath();
    ctx.moveTo(base.x, base.y);
    ctx.lineTo(top.x, top.y);
    ctx.lineWidth = 1.4;
    ctx.strokeStyle = SERIES.nasc;
    ctx.stroke();
    [base, top].forEach(function (pt, idx) {
      ctx.beginPath();
      ctx.arc(pt.x, pt.y, idx === 0 ? 2.6 : 4.2, 0, 2 * Math.PI);
      ctx.fillStyle = SERIES.nasc;
      ctx.fill();
    });

    // vertex labels
    ctx.font = "10px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = th.muted;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ["donor 2", "donor 3", "donor 4"].forEach(function (lab, idx) {
      var pr3 = project(TRI[idx].x * 1.17, TRI[idx].y * 1.17, 0);
      ctx.fillText(lab, pr3.x, pr3.y);
    });
  }

  /* ------------------------------------------------------------------ *
   * Controller
   * ------------------------------------------------------------------ */

  var PRESETS = {
    dense: [
      [0, 1],
      [0, 2],
      [0, 3],
      [1, 2],
      [1, 3],
      [2, 3],
    ],
    asymmetric: [
      [0, 2],
      [0, 3],
      [1, 3],
      [2, 3],
    ],
    partial: [
      [0, 2],
      [0, 3],
      [2, 3],
    ],
    none: [[2, 3]],
  };

  function adjFromEdges(edges) {
    var adj = zeros(N_UNITS, N_UNITS);
    edges.forEach(function (e) {
      adj[e[0]][e[1]] = 1;
      adj[e[1]][e[0]] = 1;
    });
    return adj;
  }

  function fmt(x, d) {
    if (!isFinite(x)) return "—";
    return x.toFixed(d === undefined ? 3 : d);
  }

  function init() {
    var root = document.getElementById(ROOT_ID);
    if (!root) return;

    var state = {
      adj: adjFromEdges(PRESETS.asymmetric),
      rho: 0.6,
      lambda: 25,
      alpha: 1.0,
      seed: 7,
      az: -0.6,
      el: 0.5,
      trajectory: [],
      hoverEdge: -1,
    };

    var shocks = buildShocks(state.seed);
    var el = {};
    ["net", "chart", "simplex", "rho", "lambda", "alpha", "rhoVal", "lambdaVal", "alphaVal", "readout", "reseed"].forEach(function (id) {
      el[id] = root.querySelector('[data-nasc="' + id + '"]');
    });

    var pending = false;
    function render() {
      if (pending) return;
      pending = true;
      window.requestAnimationFrame(function () {
        pending = false;
        draw();
      });
    }

    function draw() {
      var model = computeModel(shocks, state.adj, state.rho, state.lambda);

      // lambda path, for the dashed trajectory on the simplex floor
      var X0 = [];
      for (var t = 0; t < T0; t++) {
        X0.push([model.obs.yd[0][t], model.obs.yd[1][t], model.obs.yd[2][t]]);
      }
      var X1 = model.obs.y1.slice(0, T0);
      var traj = [];
      for (var k = 0; k <= 16; k++) {
        var lam = Math.pow(k / 16, 1.6) * 500;
        traj.push(solveWeights(X0, X1, model.absS, lam));
      }
      state.trajectory = traj;

      if (el.net) drawNetwork(el.net, state, model);
      if (el.chart) drawChart(el.chart, model.panel, model);
      if (el.simplex) drawSimplex(el.simplex, state, model);

      if (el.rhoVal) el.rhoVal.textContent = fmt(state.rho, 2);
      if (el.lambdaVal) el.lambdaVal.textContent = fmt(state.lambda, 0);
      if (el.alphaVal) el.alphaVal.textContent = fmt(state.alpha, 1);

      if (el.readout) {
        var rows = [
          ["True ATT", fmt(model.att.truth)],
          ["SCM", fmt(model.att.sc)],
          ["BC", fmt(model.att.bc)],
          ["CR", fmt(model.att.cr)],
          ["NASC", fmt(model.att.nasc)],
        ];
        var html =
          '<table class="nasc-table"><tbody>' +
          rows
            .map(function (r) {
              return "<tr><th>" + r[0] + "</th><td>" + r[1] + "</td></tr>";
            })
            .join("") +
          "</tbody></table>" +
          '<table class="nasc-table"><tbody>' +
          [
            ["γ<sup>sc</sup>′|s|", fmt(model.contamSc)],
            ["γ<sup>nasc</sup>′|s|", fmt(model.contamNasc)],
            ["(1−γ<sup>nasc</sup>′s)<sup>−1</sup>", fmt(1 / model.nascFactor)],
            [
              "γ<sup>nasc</sup>",
              model.gNasc
                .map(function (g) {
                  return fmt(g, 2);
                })
                .join(" · "),
            ],
            [
              "γ<sup>*</sup> (true)",
              model.gammaTrue
                .map(function (g) {
                  return fmt(g, 2);
                })
                .join(" · "),
            ],
          ]
            .map(function (r) {
              return "<tr><th>" + r[0] + "</th><td>" + r[1] + "</td></tr>";
            })
            .join("") +
          "</tbody></table>";
        el.readout.innerHTML = html;
      }
    }

    // sliders
    function bindSlider(node, key, scale) {
      if (!node) return;
      node.addEventListener("input", function () {
        state[key] = parseFloat(node.value) * (scale || 1);
        render();
      });
    }
    bindSlider(el.rho, "rho", 0.01);
    bindSlider(el.lambda, "lambda", 1);
    bindSlider(el.alpha, "alpha", 0.1);

    // presets
    Array.prototype.forEach.call(root.querySelectorAll("[data-nasc-preset]"), function (btn) {
      btn.addEventListener("click", function () {
        var key = btn.getAttribute("data-nasc-preset");
        if (!PRESETS[key]) return;
        state.adj = adjFromEdges(PRESETS[key]);
        Array.prototype.forEach.call(root.querySelectorAll("[data-nasc-preset]"), function (b) {
          b.classList.toggle("is-active", b === btn);
        });
        render();
      });
    });

    // network editing
    if (el.net) {
      el.net.addEventListener("click", function (evt) {
        if (networkHit(el.net, state, evt)) {
          Array.prototype.forEach.call(root.querySelectorAll("[data-nasc-preset]"), function (b) {
            b.classList.remove("is-active");
          });
          render();
        }
      });

      // With no handles to aim at, the hover highlight is the only affordance.
      el.net.addEventListener("mousemove", function (evt) {
        var pt = canvasPoint(el.net, evt);
        var k = nearestEdge(el.net, pt.x, pt.y);
        if (k === state.hoverEdge) return;
        state.hoverEdge = k;
        el.net.style.cursor = k >= 0 ? "pointer" : "default";
        render();
      });

      el.net.addEventListener("mouseleave", function () {
        if (state.hoverEdge === -1) return;
        state.hoverEdge = -1;
        el.net.style.cursor = "default";
        render();
      });
    }

    // rotate the simplex by dragging
    if (el.simplex) {
      var dragging = false;
      var lastX = 0;
      var lastY = 0;
      var start = function (x, y) {
        dragging = true;
        lastX = x;
        lastY = y;
      };
      var move = function (x, y) {
        if (!dragging) return;
        state.az += (x - lastX) * 0.01;
        state.el = Math.max(0.08, Math.min(1.35, state.el + (y - lastY) * 0.006));
        lastX = x;
        lastY = y;
        render();
      };
      var stop = function () {
        dragging = false;
      };
      el.simplex.addEventListener("mousedown", function (e) {
        start(e.clientX, e.clientY);
      });
      window.addEventListener("mousemove", function (e) {
        move(e.clientX, e.clientY);
      });
      window.addEventListener("mouseup", stop);
      el.simplex.addEventListener(
        "touchstart",
        function (e) {
          start(e.touches[0].clientX, e.touches[0].clientY);
        },
        { passive: true }
      );
      el.simplex.addEventListener(
        "touchmove",
        function (e) {
          if (!dragging) return;
          e.preventDefault();
          move(e.touches[0].clientX, e.touches[0].clientY);
        },
        { passive: false }
      );
      el.simplex.addEventListener("touchend", stop);
    }

    if (el.reseed) {
      el.reseed.addEventListener("click", function () {
        state.seed = Math.floor(Math.random() * 100000);
        shocks = buildShocks(state.seed);
        render();
      });
    }

    window.addEventListener("resize", render);

    // redraw when the site theme flips
    new MutationObserver(render).observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["data-theme"],
    });

    root.classList.add("is-ready");
    render();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
