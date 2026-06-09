/* Neurovascular Molecular Signatures in Neurodegeneration -- client logic.
   Consumes window.NV (see build_data.py). No framework, no build step. */
(function () {
  "use strict";
  var NV = window.NV;
  var D = NV.diseases, DO = NV.diseaseOrder, L = NV.lineages, deg = NV.deg;
  var colorByCode = {}; D.forEach(function (d) { colorByCode[d.code] = d.color; });
  var dispByCode = {}; D.forEach(function (d) { dispByCode[d.code] = d; });
  var linLabel = L.map(function (x) { return x.label; });

  // ---- gene index: gene -> [rows] -------------------------------------------
  var geneIndex = {};
  for (var i = 0; i < deg.length; i++) {
    var g = deg[i][0];
    (geneIndex[g] || (geneIndex[g] = [])).push(deg[i]);
  }
  // which lineages exist per disease (for dropdowns)
  var linByDisease = {}; DO.forEach(function (c) { linByDisease[c] = {}; });
  deg.forEach(function (rw) { linByDisease[DO[rw[1]]][rw[2]] = true; });
  var labelByCode = {}; L.forEach(function (l) { labelByCode[l.code] = l.label; });

  // ---- expression matrix (absolute view) + full gene list -------------------
  var EXPR = window.NV_EXPR || null;
  var exprIdx = {};
  var allGenes = NV.genes.slice();
  if (EXPR) {
    DO.forEach(function (c) { exprIdx[c] = {}; var ed = EXPR[c]; if (ed) ed.genes.forEach(function (g, i) { exprIdx[c][g] = i; }); });
    var seenG = {}; allGenes.forEach(function (g) { seenG[g] = 1; });
    (EXPR.allGenes || []).forEach(function (g) { if (!seenG[g]) { seenG[g] = 1; allGenes.push(g); } });
    allGenes.sort();
  }

  var $ = function (id) { return document.getElementById(id); };
  function el(tag, cls, html) { var e = document.createElement(tag); if (cls) e.className = cls; if (html != null) e.innerHTML = html; return e; }
  function fmtFC(v) { return (v >= 0 ? "+" : "") + v.toFixed(2); }
  function fmtQ(q) {
    if (q == null) return "n/s";
    if (q === 0) return "<1e-300";
    if (q < 0.001) return q.toExponential(1);
    return q.toFixed(3);
  }

  // ---- color scale (diverging: green = lower, orange = higher in disease) ----
  function lerp(a, b, t) { return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]; }
  function rampColor(stops, t) {
    for (var k = 1; k < stops.length; k++) {
      if (t <= stops[k][0]) {
        var p = (t - stops[k - 1][0]) / (stops[k][0] - stops[k - 1][0]);
        var c = lerp(stops[k - 1][1], stops[k][1], p);
        return c;
      }
    }
    return stops[stops.length - 1][1];
  }
  var POS = [[0, [245, 245, 245]], [0.4, [253, 186, 110]], [0.75, [217, 119, 6]], [1, [140, 77, 0]]];
  var NEG = [[0, [245, 245, 245]], [0.4, [120, 210, 172]], [0.75, [5, 150, 105]], [1, [3, 92, 62]]];
  function fcColor(fc) {
    if (fc == null) return [245, 245, 245];
    var m = Math.min(Math.abs(fc), 3) / 3;
    return rampColor(fc >= 0 ? POS : NEG, m);
  }
  function rgb(c) { return "rgb(" + Math.round(c[0]) + "," + Math.round(c[1]) + "," + Math.round(c[2]) + ")"; }
  function textOn(c) { return (0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]) < 140 ? "#fff" : "#333"; }

  // ===========================================================================
  // header
  // ===========================================================================
  function renderHeader() {
    var b = $("badges");
    D.forEach(function (d) {
      var s = el("span", "badge");
      s.style.background = "#f4f4f6";
      s.innerHTML = '<span class="dot" style="background:' + d.color + '"></span>' + d.label + ' · ' + d.region;
      b.appendChild(s);
    });
    var stats = [
      [3, "Diseases"], [L.length, "Cell types"],
      [NV.meta.nGenes.toLocaleString(), "Genes"],
      [NV.meta.nDeg.toLocaleString(), "DE tests"], ["MAST", "DE method"]
    ];
    var sc = $("stats");
    stats.forEach(function (s) {
      var c = el("div", "stat-card");
      c.appendChild(el("div", "num", String(s[0])));
      c.appendChild(el("div", "lbl", s[1]));
      sc.appendChild(c);
    });
    $("footer").innerHTML =
      "Neurovascular Molecular Signatures in Neurodegeneration · Rust Lab, University of Southern California · " +
      "Single-cell / single-nucleus RNA-seq, harmonised AD, HD and FTD-GRN human brain vasculature · " +
      "<b>Unpublished data — for collaborator review.</b>";
  }

  // ===========================================================================
  // gene search
  // ===========================================================================
  var input = $("geneInput"), ac = $("ac"), acItems = [], acSel = -1;
  function searchGenes(q) {
    q = q.toUpperCase();
    var starts = [], incl = [];
    for (var i = 0; i < allGenes.length && starts.length < 12; i++) {
      var g = allGenes[i];
      if (g.toUpperCase().indexOf(q) === 0) starts.push(g);
    }
    if (starts.length < 12) {
      for (var j = 0; j < allGenes.length && incl.length < 12 - starts.length; j++) {
        var gg = allGenes[j];
        if (gg.toUpperCase().indexOf(q) > 0) incl.push(gg);
      }
    }
    return starts.concat(incl).slice(0, 12);
  }
  function showAC() {
    var q = input.value.trim();
    if (!q) { ac.style.display = "none"; return; }
    acItems = searchGenes(q); acSel = -1;
    if (!acItems.length) { ac.style.display = "none"; return; }
    ac.innerHTML = "";
    acItems.forEach(function (g, idx) {
      var d = el("div");
      var tag = geneIndex[g] ? (geneIndex[g].length + ' DE tests') : 'expression';
      d.innerHTML = '<span class="g">' + g + '</span> <span class="muted" style="font-size:12px">' + tag + '</span>';
      d.onclick = function () { pick(g); };
      ac.appendChild(d);
    });
    ac.style.display = "block";
  }
  function pick(g) { input.value = g; ac.style.display = "none"; renderGene(g); }
  input.addEventListener("input", showAC);
  input.addEventListener("keydown", function (e) {
    var items = ac.children;
    if (e.key === "ArrowDown") { acSel = Math.min(acSel + 1, items.length - 1); e.preventDefault(); }
    else if (e.key === "ArrowUp") { acSel = Math.max(acSel - 1, 0); e.preventDefault(); }
    else if (e.key === "Enter") { if (acSel >= 0 && acItems[acSel]) pick(acItems[acSel]); else if (acItems[0]) pick(acItems[0]); return; }
    else return;
    for (var k = 0; k < items.length; k++) items[k].className = (k === acSel ? "sel" : "");
  });
  document.addEventListener("click", function (e) { if (!ac.contains(e.target) && e.target !== input) ac.style.display = "none"; });

  var curGene = null, searchView = "diff";
  function renderGene(gene) { curGene = gene; renderView(); }
  function renderView() { if (!curGene) return; (searchView === "abs" ? renderAbs : renderDiff)(curGene); }
  function setView(v) {
    searchView = v;
    $("search-view").querySelectorAll("button").forEach(function (b) { b.classList.toggle("active", b.dataset.v === v); });
    renderView();
  }
  function initSearchView() {
    $("search-view").querySelectorAll("button").forEach(function (btn) { btn.onclick = function () { setView(btn.dataset.v); }; });
  }

  function renderDiff(gene) {
    var rows = geneIndex[gene];
    var box = $("result");
    box.style.display = "block";
    if (!rows || !rows.length) {
      box.innerHTML = '<div class="rh"><h3>' + gene + '</h3><span class="meta">not a significant DEG in any cell type</span></div>' +
        '<div class="empty">No significant disease-vs-control change for <b>' + gene + '</b>. Switch to <b class="goabs" style="cursor:pointer;color:#1565c0">Absolute expression</b> to see its levels per cell type.</div>';
      var ga = box.querySelector(".goabs"); if (ga) ga.onclick = function () { setView("abs"); };
      return;
    }
    // lineages present (in canonical order), and lookup by (disease,lineage)
    var present = {}, byKey = {};
    rows.forEach(function (rw) { present[rw[2]] = true; byKey[rw[1] + "_" + rw[2]] = rw; });
    var linIdxs = Object.keys(present).map(Number).sort(function (a, b) { return a - b; });

    var html = '<div class="rh"><h3>' + gene + '</h3>' +
      '<span class="meta">significant change in ' + rows.length + ' cell-type · disease combination' + (rows.length > 1 ? 's' : '') +
      ' &nbsp;·&nbsp; <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + encodeURIComponent(gene) + '" target="_blank" rel="noopener">GeneCards ↗</a></span></div>';

    // heatmap
    html += '<div class="hm"><table class="heat"><thead><tr><th class="row">Cell type</th>';
    DO.forEach(function (c) { html += '<th><span style="color:' + colorByCode[c] + '">●</span> ' + dispByCode[c].short + '</th>'; });
    html += '</tr></thead><tbody>';
    linIdxs.forEach(function (li) {
      html += '<tr><td class="lbl">' + linLabel[li] + ' <span class="grp">' + L[li].group + '</span></td>';
      DO.forEach(function (c, di) {
        var rw = byKey[di + "_" + li];
        if (!rw) { html += '<td class="cell na">·</td>'; return; }
        var col = fcColor(rw[3]); var sig = (rw[6] != null && rw[6] <= 0.05);
        var tip = dispByCode[c].short + " · " + linLabel[li] + "\nlog2FC " + fmtFC(rw[3]) +
          "\npct " + (rw[4] == null ? "?" : rw[4]) + " vs " + (rw[5] == null ? "?" : rw[5]) + "\nq = " + fmtQ(rw[6]);
        html += '<td class="cell' + (sig ? ' sig' : '') + '" style="background:' + rgb(col) + ';color:' + textOn(col) + '" title="' + tip + '">' + fmtFC(rw[3]) + '</td>';
      });
      html += '</tr>';
    });
    html += '</tbody></table></div>';

    // legend
    html += '<div class="hm-legend">' +
      '<div class="scalebar"><span>Lower in disease</span><span class="bar"></span><span>Higher in disease</span></div>' +
      '<div class="lg-item"><span class="lg-sw" style="box-shadow:inset 0 0 0 2px #1a1a1a;background:#fde9d0"></span>FDR ≤ 0.05</div>' +
      '<div class="lg-item"><span class="lg-sw" style="background:#f5f5f5"></span>not significant / absent</div></div>';

    // detail table
    html += '<table class="dt"><thead><tr><th>Disease</th><th>Cell type</th><th>log2FC</th><th>pct (dis)</th><th>pct (ctrl)</th><th>FDR</th></tr></thead><tbody>';
    rows.slice().sort(function (a, b) { return a[1] - b[1] || a[2] - b[2]; }).forEach(function (rw) {
      html += '<tr><td><span style="color:' + colorByCode[DO[rw[1]]] + '">●</span> ' + dispByCode[DO[rw[1]]].short + '</td>' +
        '<td>' + linLabel[rw[2]] + '</td>' +
        '<td class="' + (rw[3] >= 0 ? 'fc-pos' : 'fc-neg') + '">' + fmtFC(rw[3]) + '</td>' +
        '<td>' + (rw[4] == null ? '—' : rw[4]) + '</td><td>' + (rw[5] == null ? '—' : rw[5]) + '</td>' +
        '<td>' + fmtQ(rw[6]) + '</td></tr>';
    });
    html += '</tbody></table>';
    box.innerHTML = html;
  }

  // ---- absolute expression view ---------------------------------------------
  var BLUES = [[0, [247, 247, 247]], [0.25, [198, 219, 239]], [0.6, [66, 146, 198]], [1, [8, 48, 107]]];
  function exprColor(v, max) { if (v == null) return [245, 245, 245]; return rampColor(BLUES, max > 0 ? Math.min(v / max, 1) : 0); }
  function exprVal(d, gene, linCode, cond) {
    var ed = EXPR && EXPR[d]; if (!ed) return null;
    var gi = exprIdx[d][gene]; if (gi == null) return null;
    var j = ed.lin.indexOf(linCode); if (j < 0) return null;
    return (cond === "C" ? ed.avgC : ed.avgD)[gi * ed.lin.length + j] / 100;
  }
  function renderAbs(gene) {
    var box = $("result"); box.style.display = "block";
    var head = '<div class="rh"><h3>' + gene + '</h3>';
    if (!EXPR) { box.innerHTML = head + '</div><div class="empty">Expression data not loaded.</div>'; return; }
    var present = {};
    DO.forEach(function (d) { var ed = EXPR[d]; if (ed && exprIdx[d][gene] != null) ed.lin.forEach(function (Lc) { present[Lc] = true; }); });
    var linCodes = L.map(function (l) { return l.code; }).filter(function (c) { return present[c]; });
    if (!linCodes.length) {
      box.innerHTML = head + '<span class="meta">below detection floor</span></div>' +
        '<div class="empty">No expression record for <b>' + gene + '</b> in these datasets (very low or undetected).</div>';
      return;
    }
    var gmax = 0;
    linCodes.forEach(function (Lc) { DO.forEach(function (d) { ["C", "D"].forEach(function (cc) { var v = exprVal(d, gene, Lc, cc); if (v != null && v > gmax) gmax = v; }); }); });
    var html = head + '<span class="meta">mean log-normalised expression · colour scaled to this gene\'s max (' + gmax.toFixed(2) + ') · <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + encodeURIComponent(gene) + '" target="_blank" rel="noopener">GeneCards ↗</a></span></div>';
    html += '<div class="hm"><table class="heat"><thead><tr><th class="row">Cell type</th>';
    DO.forEach(function (d) { html += '<th colspan="2" style="border-bottom:2px solid ' + colorByCode[d] + '"><span style="color:' + colorByCode[d] + '">●</span> ' + dispByCode[d].short + '</th>'; });
    html += '</tr><tr><th class="row"></th>';
    DO.forEach(function () { html += '<th style="font-weight:500;color:#999;font-size:10.5px">Ctrl</th><th style="font-weight:500;color:#999;font-size:10.5px">Dis</th>'; });
    html += '</tr></thead><tbody>';
    linCodes.forEach(function (Lc) {
      html += '<tr><td class="lbl">' + labelByCode[Lc] + '</td>';
      DO.forEach(function (d) {
        ["C", "D"].forEach(function (cc) {
          var v = exprVal(d, gene, Lc, cc);
          if (v == null) { html += '<td class="cell na">·</td>'; return; }
          var col = exprColor(v, gmax);
          html += '<td class="cell" style="background:' + rgb(col) + ';color:' + textOn(col) +
            '" title="' + dispByCode[d].short + ' · ' + labelByCode[Lc] + ' · ' + (cc === "C" ? "Control" : "Disease") + '\nmean expr ' + v.toFixed(2) + '">' + v.toFixed(1) + '</td>';
        });
      });
      html += '</tr>';
    });
    html += '</tbody></table></div>';
    html += '<div class="hm-legend"><div class="scalebar"><span>0</span><span class="bar" style="background:linear-gradient(90deg,#f7f7f7,#c6dbef,#4292c6,#08306b)"></span><span>' + gmax.toFixed(1) + '</span></div>' +
      '<div class="lg-item">Mean log-normalised expression, control vs disease per cell type</div></div>';
    box.innerHTML = html;
  }

  // ===========================================================================
  // cell-type DEGs
  // ===========================================================================
  function fillDiseaseSelect(sel) { sel.innerHTML = ""; D.forEach(function (d) { sel.appendChild(new Option(d.label, d.code)); }); }
  function fillLineageSelect(sel, code, def) {
    sel.innerHTML = "";
    L.forEach(function (l, idx) {
      if (linByDisease[code][idx]) { var o = new Option(l.label, idx); sel.appendChild(o); }
    });
    if (def != null) { for (var i = 0; i < sel.options.length; i++) if (+sel.options[i].value === def) sel.value = def; }
  }
  var depDir = "up";
  function renderDeps() {
    var code = $("dep-disease").value, li = +$("dep-lineage").value;
    var list = deg.filter(function (rw) {
      return DO[rw[1]] === code && rw[2] === li && rw[6] != null && rw[6] <= 0.05 && (depDir === "up" ? rw[3] > 0 : rw[3] < 0);
    }).sort(function (a, b) { return depDir === "up" ? b[3] - a[3] : a[3] - b[3]; }).slice(0, 60);
    var box = $("dep-list"); box.innerHTML = "";
    if (!list.length) { box.appendChild(el("div", "empty", "No significant genes for this selection.")); return; }
    list.forEach(function (rw, idx) {
      var it = el("div", "dep-item");
      it.innerHTML = '<div class="dep-rank">' + (idx + 1) + '</div><div class="dep-gene">' + rw[0] + '</div>' +
        '<div class="dep-desc">pct ' + (rw[4] == null ? '?' : rw[4]) + ' vs ' + (rw[5] == null ? '?' : rw[5]) + '</div>' +
        '<div class="' + (rw[3] >= 0 ? 'fc-pos' : 'fc-neg') + '">' + fmtFC(rw[3]) + '</div>' +
        '<div class="dep-q">q=' + fmtQ(rw[6]) + '</div>';
      it.onclick = function () { goGene(rw[0]); };
      box.appendChild(it);
    });
  }
  function initDeps() {
    fillDiseaseSelect($("dep-disease"));
    fillLineageSelect($("dep-lineage"), "AD", null);
    $("dep-disease").onchange = function () { fillLineageSelect($("dep-lineage"), this.value, null); renderDeps(); };
    $("dep-lineage").onchange = renderDeps;
    $("dep-dir").querySelectorAll("button").forEach(function (btn) {
      btn.onclick = function () {
        $("dep-dir").querySelectorAll("button").forEach(function (b) { b.classList.remove("active"); });
        btn.classList.add("active"); depDir = btn.dataset.dir; renderDeps();
      };
    });
    renderDeps();
  }

  // ===========================================================================
  // cross-disease overlap (Venn + concordant genes)
  // ===========================================================================
  function vennData(li) {
    var sets = { AD: {}, HD: {}, FTD: {} };
    deg.forEach(function (rw) {
      if (rw[2] !== li) return;
      if (rw[6] == null || rw[6] > 0.05 || Math.abs(rw[3]) < 0.5) return;
      sets[DO[rw[1]]][rw[0]] = true;
    });
    var AD = sets.AD, HD = sets.HD, FTD = sets.FTD;
    function has(s, g) { return !!s[g]; }
    var keys = {}; [AD, HD, FTD].forEach(function (s) { for (var g in s) keys[g] = true; });
    var c = { AD_only: 0, HD_only: 0, FTD_only: 0, AD_HD: 0, AD_FTD: 0, HD_FTD: 0, all3: 0 };
    for (var g in keys) {
      var a = has(AD, g), h = has(HD, g), f = has(FTD, g);
      if (a && h && f) c.all3++;
      else if (a && h) c.AD_HD++;
      else if (a && f) c.AD_FTD++;
      else if (h && f) c.HD_FTD++;
      else if (a) c.AD_only++;
      else if (h) c.HD_only++;
      else if (f) c.FTD_only++;
    }
    c.totAD = Object.keys(AD).length; c.totHD = Object.keys(HD).length; c.totFTD = Object.keys(FTD).length;
    return c;
  }
  function drawVenn(c) {
    var col = { AD: colorByCode.AD, HD: colorByCode.HD, FTD: colorByCode.FTD };
    function circ(cx, cy, fill) { return '<circle cx="' + cx + '" cy="' + cy + '" r="66" fill="' + fill + '" fill-opacity="0.20" stroke="' + fill + '" stroke-width="1.5"/>'; }
    function txt(x, y, s, fill, size, w) { return '<text x="' + x + '" y="' + y + '" text-anchor="middle" font-size="' + (size || 13) + '" font-weight="' + (w || 400) + '" fill="' + (fill || '#333') + '">' + s + '</text>'; }
    var s = '<svg viewBox="0 0 280 250" width="280" height="250">';
    s += circ(100, 96, col.AD) + circ(180, 96, col.HD) + circ(140, 158, col.FTD);
    s += txt(70, 40, 'AD', col.AD, 13, 700) + txt(210, 40, 'HD', col.HD, 13, 700) + txt(140, 238, 'FTD', col.FTD, 13, 700);
    s += txt(74, 92, c.AD_only, '#333', 14, 700) + txt(206, 92, c.HD_only, '#333', 14, 700) + txt(140, 196, c.FTD_only, '#333', 14, 700);
    s += txt(140, 74, c.AD_HD, '#555') + txt(102, 150, c.AD_FTD, '#555') + txt(178, 150, c.HD_FTD, '#555');
    s += txt(140, 122, c.all3, '#000', 15, 700);
    s += '</svg>';
    s += '<div class="pw-sub" style="margin-top:6px">Significant DEGs (FDR ≤ 0.05, |log2FC| ≥ 0.5): AD ' + c.totAD + ' · HD ' + c.totHD + ' · FTD ' + c.totFTD + '</div>';
    return s;
  }
  function renderOverlap() {
    var li = +$("ov-lineage").value;
    $("venn").innerHTML = drawVenn(vennData(li));
    var entry = null;
    NV.common.forEach(function (cc) { if (cc.lin === li) entry = cc; });
    var box = $("shared-genes");
    if (!entry || (entry.nUp + entry.nDown) === 0) {
      box.innerHTML = '<div class="empty">No genes are dysregulated in the same direction across all three diseases in this cell type.</div>';
      return;
    }
    function chips(arr, cls) { return arr.map(function (g) { return '<span class="chip ' + cls + '" data-g="' + g + '">' + g + '</span>'; }).join(""); }
    box.innerHTML =
      '<div class="pw-sub"><b style="color:#9a5b00">▲ ' + entry.nUp + ' higher</b> in disease</div><div class="chips">' + chips(entry.up, "up") + '</div>' +
      '<div class="pw-sub" style="margin-top:12px"><b style="color:#0a6b48">▼ ' + entry.nDown + ' lower</b> in disease</div><div class="chips">' + chips(entry.down, "down") + '</div>';
    box.querySelectorAll(".chip").forEach(function (ch) { ch.onclick = function () { goGene(ch.dataset.g); }; });
  }
  function initOverlap() {
    var sel = $("ov-lineage"); sel.innerHTML = "";
    // only lineages that have a concordant set, most-convergent first
    NV.common.forEach(function (cc) { if (cc.nUp + cc.nDown > 0) sel.appendChild(new Option(L[cc.lin].label + " (" + (cc.nUp + cc.nDown) + " shared)", cc.lin)); });
    sel.onchange = renderOverlap;
    renderOverlap();
  }

  // ===========================================================================
  // pericytes
  // ===========================================================================
  function initPericytes() {
    var p = NV.pericyte;
    // proportion chart
    var labels = DO.map(function (c) { return dispByCode[c].short; });
    var ctrl = DO.map(function (c) { return p.subtype[c] ? p.subtype[c].M_ctrl : null; });
    var dis = DO.map(function (c) { return p.subtype[c] ? p.subtype[c].M_dis : null; });
    new Chart($("peri-prop").getContext("2d"), {
      type: "bar",
      data: {
        labels: labels,
        datasets: [
          { label: "Control", data: ctrl, backgroundColor: "#b8c2cc" },
          { label: "Disease", data: dis, backgroundColor: "#9333ea" }
        ]
      },
      options: {
        responsive: true, maintainAspectRatio: false,
        plugins: { legend: { position: "bottom" }, tooltip: { callbacks: { label: function (c) { return c.dataset.label + ": " + c.raw + "% M-peri"; } } } },
        scales: { y: { beginAtZero: true, max: 60, title: { display: true, text: "% M-pericytes" } } }
      }
    });
    var o = p.overlap;
    $("peri-overlap").innerHTML = "Pericyte DEGs per disease: <b>AD " + (o.AD || 0) + "</b> · <b>HD " + (o.HD || 0) +
      "</b> · <b>FTD " + (o.FTD || 0) + "</b>. <b>" + (o["AD∩FTD∩HD"] || 0) + "</b> shared by all three.";
    function chipset(arr) { return arr.map(function (g) { return '<span class="chip" data-g="' + g + '">' + g + '</span>'; }).join(""); }
    $("peri-convergent").innerHTML = chipset(p.convergent);
    $("peri-mk-m").innerHTML = chipset(p.markers.M);
    $("peri-mk-t").innerHTML = chipset(p.markers.T);
    document.querySelectorAll("#pericytes .chip").forEach(function (ch) { ch.onclick = function () { goGene(ch.dataset.g); }; });
  }

  // ===========================================================================
  // endothelial zonation
  // ===========================================================================
  var endoDir = "up";
  function renderEndo() {
    var code = $("endo-disease").value, seg = +$("endo-seg").value;
    var di = DO.indexOf(code);
    var list = NV.endo.filter(function (rw) {
      return rw[1] === di && rw[2] === seg && rw[6] != null && rw[6] <= 0.05 && (endoDir === "up" ? rw[3] > 0 : rw[3] < 0);
    }).sort(function (a, b) { return endoDir === "up" ? b[3] - a[3] : a[3] - b[3]; }).slice(0, 60);
    var box = $("endo-list"); box.innerHTML = "";
    if (!list.length) { box.appendChild(el("div", "empty", "No significant genes for this selection.")); return; }
    list.forEach(function (rw, idx) {
      var it = el("div", "dep-item");
      it.innerHTML = '<div class="dep-rank">' + (idx + 1) + '</div><div class="dep-gene">' + rw[0] + '</div>' +
        '<div class="dep-desc">pct ' + (rw[4] == null ? '?' : rw[4]) + ' vs ' + (rw[5] == null ? '?' : rw[5]) + '</div>' +
        '<div class="' + (rw[3] >= 0 ? 'fc-pos' : 'fc-neg') + '">' + fmtFC(rw[3]) + '</div>' +
        '<div class="dep-q">q=' + fmtQ(rw[6]) + '</div>';
      it.onclick = function () { goGene(rw[0]); };
      box.appendChild(it);
    });
  }
  function initEndo() {
    fillDiseaseSelect($("endo-disease"));
    var sel = $("endo-seg"); sel.innerHTML = "";
    NV.segs.forEach(function (s, idx) { sel.appendChild(new Option(s, idx)); });
    // default to Capillary if present
    var capIdx = NV.segs.indexOf("Capillary"); if (capIdx >= 0) sel.value = capIdx;
    $("endo-disease").onchange = renderEndo; sel.onchange = renderEndo;
    $("endo-dir").querySelectorAll("button").forEach(function (btn) {
      btn.onclick = function () {
        $("endo-dir").querySelectorAll("button").forEach(function (b) { b.classList.remove("active"); });
        btn.classList.add("active"); endoDir = btn.dataset.dir; renderEndo();
      };
    });
    renderEndo();
  }

  // ===========================================================================
  // pathways
  // ===========================================================================
  var pwOnt = "BP";
  function renderPathways() {
    var code = $("pw-disease").value, scope = $("pw-scope").value;
    var terms = (((NV.pathways[code] || {})[scope] || {})[pwOnt] || []).slice();
    terms.sort(function (a, b) { return (b.fold || 0) - (a.fold || 0); });
    var box = $("pw-list"); box.innerHTML = "";
    if (!terms.length) { box.appendChild(el("div", "empty", "No enriched terms for this selection.")); return; }
    var maxFold = Math.max.apply(null, terms.map(function (t) { return t.fold || 0; }));
    terms.forEach(function (t) {
      var row = el("div", "pw-row");
      var w = Math.max(2, (t.fold / maxFold) * 100);
      row.innerHTML = '<div class="pw-label">' + t.desc + ' <span class="pw-sub">(' + t.id + ')</span></div>' +
        '<div class="pw-bar-wrap"><div class="pw-bar" style="width:' + w + '%" title="fold enrichment ' + t.fold + ' · q=' + fmtQ(t.q) + ' · ' + t.count + ' genes"></div>' +
        '<span class="pw-sub">' + t.fold + '× · ' + t.count + ' genes</span></div>';
      box.appendChild(row);
    });
  }
  function initPathways() {
    fillDiseaseSelect($("pw-disease"));
    var sc = $("pw-scope"); sc.innerHTML = ""; sc.appendChild(new Option("All cell types", "All")); sc.appendChild(new Option("Pericytes", "Pericyte"));
    $("pw-disease").onchange = renderPathways; sc.onchange = renderPathways;
    $("pw-ont").querySelectorAll("button").forEach(function (btn) {
      btn.onclick = function () {
        $("pw-ont").querySelectorAll("button").forEach(function (b) { b.classList.remove("active"); });
        btn.classList.add("active"); pwOnt = btn.dataset.ont; renderPathways();
      };
    });
    renderPathways();
  }

  // ===========================================================================
  // signaling (curated LR)
  // ===========================================================================
  function chgCell(v) {
    if (v > 0) return '<span class="chg up" title="enhanced in disease">▲</span>';
    if (v < 0) return '<span class="chg down" title="reduced in disease">▼</span>';
    return '<span class="chg zero" title="not prominently changed">–</span>';
  }
  function initLR() {
    var body = $("lr-body");
    NV.lr.forEach(function (r) {
      var tr = el("tr");
      tr.innerHTML = '<td><b>' + r.pathway + '</b></td><td>' + r.sender + '</td><td class="arrow">→</td><td>' + r.receiver + '</td>' +
        '<td>' + chgCell(r.chg.AD) + '</td><td>' + chgCell(r.chg.HD) + '</td><td>' + chgCell(r.chg.FTD) + '</td>' +
        '<td><span class="tag ' + (r.shared ? 'tag-shared">shared' : 'tag-spec">disease-specific') + '</span></td>' +
        '<td class="pw-sub">' + r.note + '</td>';
      body.appendChild(tr);
    });
  }

  // ===========================================================================
  // cell atlas (UMAP scatter)
  // ===========================================================================
  var LINPAL = {
    Endothelial: "#d6336c", Pericyte: "#7048e8", vSMC: "#f76707", M_Fibro: "#0ca678",
    P_Fibro: "#66a80f", P_Mac: "#c2255c", Microglia: "#1c7ed6", Astro: "#f59f00",
    Oligo: "#495057", OPC: "#adb5bd", Neuron_exc: "#9c36b5", Neuron_inh: "#3bc9db"
  };
  var atlasBy = "lineage";
  function renderAtlas() {
    var code = $("atlas-disease").value;
    var pts = window.NV_UMAP.filter(function (p) { return p[3] === code; });
    var cv = $("atlas-canvas"), ctx = cv.getContext("2d");
    ctx.clearRect(0, 0, cv.width, cv.height);
    if (!pts.length) return;
    var xs = pts.map(function (p) { return p[0]; }), ys = pts.map(function (p) { return p[1]; });
    var minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs);
    var minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);
    var pad = 22, W = cv.width, H = cv.height;
    function sx(x) { return pad + (x - minX) / (maxX - minX) * (W - 2 * pad); }
    function sy(y) { return H - pad - (y - minY) / (maxY - minY) * (H - 2 * pad); }
    ctx.globalAlpha = 0.55;
    pts.forEach(function (p) {
      ctx.fillStyle = atlasBy === "lineage" ? (LINPAL[p[2]] || "#999") : (p[4] === 1 ? colorByCode[code] : "#c2c8ce");
      ctx.beginPath(); ctx.arc(sx(p[0]), sy(p[1]), 1.8, 0, 6.2832); ctx.fill();
    });
    ctx.globalAlpha = 1;
    var leg = $("atlas-legend");
    if (atlasBy === "lineage") {
      $("atlas-legend-title").textContent = "Cell type";
      var present = {}; pts.forEach(function (p) { present[p[2]] = (present[p[2]] || 0) + 1; });
      leg.innerHTML = L.filter(function (l) { return present[l.code]; }).map(function (l) {
        return '<div class="lg-item" style="margin:4px 0"><span class="lg-sw" style="background:' + (LINPAL[l.code] || "#999") + '"></span>' + l.label + ' <span class="muted">' + present[l.code] + '</span></div>';
      }).join("");
    } else {
      $("atlas-legend-title").textContent = "Condition";
      leg.innerHTML = '<div class="lg-item" style="margin:4px 0"><span class="lg-sw" style="background:#c2c8ce"></span>Control</div>' +
        '<div class="lg-item" style="margin:4px 0"><span class="lg-sw" style="background:' + colorByCode[code] + '"></span>Disease (' + dispByCode[code].short + ')</div>';
    }
  }
  function initAtlas() {
    var sec = document.getElementById("atlas");
    if (!window.NV_UMAP || !window.NV_UMAP.length) { if (sec) sec.style.display = "none"; return; }
    fillDiseaseSelect($("atlas-disease"));
    $("atlas-disease").onchange = renderAtlas;
    $("atlas-color").querySelectorAll("button").forEach(function (btn) {
      btn.onclick = function () {
        $("atlas-color").querySelectorAll("button").forEach(function (b) { b.classList.remove("active"); });
        btn.classList.add("active"); atlasBy = btn.dataset.by; renderAtlas();
      };
    });
    // force the intended default (colour by cell type)
    atlasBy = "lineage";
    $("atlas-color").querySelectorAll("button").forEach(function (b) { b.classList.toggle("active", b.dataset.by === "lineage"); });
    renderAtlas();
  }

  // ===========================================================================
  // navigation helpers
  // ===========================================================================
  function goGene(g) {
    input.value = g; renderGene(g);
    document.getElementById("search").scrollIntoView({ behavior: "smooth" });
  }
  function initNav() {
    var links = Array.prototype.slice.call(document.querySelectorAll("#nav a"));
    links.forEach(function (a) {
      a.onclick = function (e) { e.preventDefault(); document.querySelector(a.getAttribute("href")).scrollIntoView({ behavior: "smooth" }); };
    });
    var secs = links.map(function (a) { return document.querySelector(a.getAttribute("href")); });
    window.addEventListener("scroll", function () {
      var y = window.scrollY + 90, cur = 0;
      secs.forEach(function (s, i) { if (s.offsetTop <= y) cur = i; });
      links.forEach(function (a, i) { a.classList.toggle("active", i === cur); });
    });
  }
  document.querySelectorAll(".hint b").forEach(function (b) { b.onclick = function () { goGene(b.textContent); }; });

  // ---- boot -----------------------------------------------------------------
  renderHeader();
  initDeps();
  initOverlap();
  initPericytes();
  initEndo();
  initPathways();
  initLR();
  initAtlas();
  initSearchView();
  initNav();
  renderGene("ADAMTS9"); // a convergent, vascular headline gene
})();
