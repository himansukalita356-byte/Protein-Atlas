/* Protein Atlas - the analysis workspace shared by protein pages and Sequence Lab.
   Every number shown here is read from ONE analysis object produced by the server-side engine
   (atlas_engine.analyze_sequence); nothing is recalculated in the browser. */
(function (global) {
  "use strict";
  const Atlas = global.Atlas;
  const { el, clear, fmt, api, states, charts, table, COLORS } = Atlas;
  const ok = (section) => section && section.status === "success";

  // ---------------------------------------------------------------- small builders
  function tile(label, value, unit, hint) {
    return el("div", { class: "stat-tile", title: hint || null }, [
      el("div", { class: "stat-label" }, label),
      el("div", { class: "stat-value" }, [typeof value === "number" ? fmt.text(value) : value, unit ? el("small", {}, " " + unit) : null]),
    ]);
  }
  function card(title, note, children) {
    return el("div", { class: "card" }, [el("h3", { class: "panel-title" }, title), note ? el("p", { class: "panel-note" }, note) : null, children]);
  }
  function kv(rows) {
    return el("div", { class: "kv" }, rows.map(([k, v]) => el("div", { class: "kv-row" }, [el("span", {}, k), el("b", {}, typeof v === "number" ? fmt.text(v) : v)])));
  }
  function guard(section, render, label) {
    if (!section) return states.noDataBlock(label || "Data unavailable");
    if (section.status === "error") return states.stateBlock("\u26a0", "Analysis unavailable \u2014 calculation error", section.message || null);
    if (section.status === "no_data") return states.noDataBlock(section.message || "Data unavailable");
    if (section.status === "invalid_input") return states.stateBlock("\u270e", "Check your input", section.message || null);
    return render(section);
  }
  const residueLabel = (r) => r.name + " (" + r.three + " \u00b7 " + r.code + ")";
  const seqOf = (A) => A.analysis.sequence.residues;

  // ---------------------------------------------------------------- sequence viewer / motif map
  function rangesIn(hits, from, to) {
    return hits.filter((h) => h.end >= from + 1 && h.start <= to).map((h) => [Math.max(h.start - 1, from), Math.min(h.end, to)]).sort((a, b) => a[0] - b[0]);
  }
  function spaced(seq, from, to) {           // blocks of 10 keyed to absolute position
    let out = "";
    for (let i = from; i < to; i++) { if (i > from && i % 10 === 0) out += " "; out += seq[i]; }
    return out;
  }
  function sequenceViewer(sequence, options) {
    options = options || {};
    const line = options.line || 60, hits = options.hits || [], hits2 = options.hits2 || [];
    const wrap = el("div", { class: "sequence-viewer", role: "region", "aria-label": "Protein sequence" });
    const rows = Math.ceil(sequence.length / line);
    let drawn = 0;
    function draw(count) {
      const end = Math.min(rows, drawn + count);
      for (let r = drawn; r < end; r++) {
        const from = r * line, to = Math.min(sequence.length, from + line);
        const cuts = rangesIn(hits, from, to).map((x) => ({ x, cls: "hl" })).concat(rangesIn(hits2, from, to).map((x) => ({ x, cls: "hl2" }))).sort((a, b) => a.x[0] - b.x[0]);
        const letters = el("span", { class: "letters" });
        let cursor = from;
        cuts.forEach((c) => {
          const [a, b] = [Math.max(c.x[0], cursor), c.x[1]];
          if (b <= a) return;
          if (a > cursor) letters.appendChild(document.createTextNode(spaced(sequence, cursor, a)));
          letters.appendChild(el("span", { class: c.cls }, spaced(sequence, a, b)));
          cursor = b;
        });
        if (cursor < to) letters.appendChild(document.createTextNode(spaced(sequence, cursor, to)));
        wrap.appendChild(el("div", { class: "seq-row" }, [el("span", { class: "pos" }, String(from + 1)), letters]));
      }
      drawn = end;
    }
    draw(250);
    if (rows > 250) {
      const more = el("button", { class: "btn btn-ghost", type: "button", onclick: () => { draw(400); if (drawn >= rows) more.remove(); } }, "Show more of the sequence");
      return el("div", {}, [wrap, more]);
    }
    return wrap;
  }
  Atlas.sequenceViewer = sequenceViewer;

  function motifTrack(length, rows) {
    const ticks = [0, 0.25, 0.5, 0.75, 1].map((f) => Math.max(1, Math.round(f * length)));
    const lanes = rows.map((row) => el("div", { class: "motif-row" }, [
      el("div", { class: "motif-label", title: row.title || row.label }, row.label),
      el("div", { class: "motif-lane" }, row.hits.map((h) => el("div", {
        class: "motif-seg", style: "left:" + ((h.start - 1) / length * 100).toFixed(3) + "%;width:" + Math.max((h.end - h.start + 1) / length * 100, 0.6).toFixed(3) + "%;background:" + (row.color || COLORS.teal),
        title: row.label + " \u00b7 " + h.start + "\u2013" + h.end + (h.sequence ? " \u00b7 " + h.sequence : ""),
      }))),
    ]));
    return el("div", { class: "motif-map", role: "img", "aria-label": "Map of matches along the sequence" }, [
      ...lanes,
      el("div", { class: "motif-axis" }, ticks.map((t) => el("span", {}, String(t)))),
    ]);
  }
  Atlas.motifTrack = motifTrack;
  const hitsText = (hits) => hits.slice(0, 6).map((h) => (h.start === h.end ? h.start : h.start + "\u2013" + h.end)).join(", ") + (hits.length > 6 ? " \u2026 (+" + (hits.length - 6) + ")" : "");

  // ---------------------------------------------------------------- pattern search (used twice)
  const PATTERN_EXAMPLES = [["N-glycosylation", "N-{P}-[ST]-{P}"], ["Phosphorylation (PKC)", "[ST]-x-[RK]"], ["Cys-Cys", "C-x(2,4)-C"], ["Poly-Gly", "GGG"], ["Single residue", "G"]];
  function patternSearchPanel(getSequence, options) {
    options = options || {};
    const input = el("input", { type: "text", placeholder: "Literal fragment (GGG) or PROSITE pattern (N-{P}-[ST]-{P})", "aria-label": "Pattern", value: options.initial || "", autocomplete: "off", spellcheck: "false" });
    const out = el("div", { class: "mt" });
    const busy = { n: 0 };
    async function run() {
      const pattern = input.value.trim();
      clear(out);
      if (!pattern) { out.appendChild(states.stateBlock("\u270e", "Enter a pattern", "Type a residue fragment such as G, or a PROSITE pattern such as N-{P}-[ST]-{P}.")); return; }
      out.appendChild(states.loadingBlock("Searching\u2026"));
      const ticket = ++busy.n;
      let result;
      try { result = await api("/pattern", { body: { sequence: getSequence(), pattern } }); }
      catch (error) { if (ticket === busy.n) { clear(out); out.appendChild(states.errorBlock(error, run)); } return; }
      if (ticket !== busy.n) return;
      clear(out);
      const n = result.count;
      out.appendChild(el("div", { class: "flex flex-wrap", style: "margin-bottom:12px" }, [
        el("span", { class: "badge " + (n ? "badge-ok" : "") }, n ? fmt.int(n) + " match" + (n === 1 ? "" : "es") : "Search completed \u2014 no matches"),
        el("span", { class: "muted" }, (result.kind === "prosite" ? "PROSITE-style pattern" : "Literal fragment") + " \u00b7 " + fmt.int(result.sequence_length) + " residues searched"),
      ]));
      if (n) {
        out.appendChild(table([
          { label: "Match", align: "num", render: (r) => String(r.match) },
          { label: "Start", align: "num", render: (r) => String(r.start) },
          { label: "End", align: "num", render: (r) => String(r.end) },
          { label: "Sequence", align: "mono", render: (r) => r.sequence },
        ], result.matches.slice(0, 500), { id: "pattern-results" }));
        if (n > 500) out.appendChild(el("p", { class: "muted" }, "Showing the first 500 of " + fmt.int(n) + " matches."));
        out.appendChild(el("h4", { class: "mt" }, "Matches in the sequence"));
        out.appendChild(motifTrack(result.sequence_length, [{ label: pattern.length > 22 ? pattern.slice(0, 20) + "\u2026" : pattern, title: pattern, hits: result.matches }]));
        out.appendChild(sequenceViewer(getSequence().replace(/\s+/g, "").toUpperCase(), { hits: result.matches }));
      }
    }
    const examples = el("div", { class: "flex flex-wrap", style: "margin:10px 0 0" }, PATTERN_EXAMPLES.map(([label, pat]) =>
      el("button", { class: "btn btn-ghost", type: "button", style: "font-size:12px;padding:5px 12px", onclick: () => { input.value = pat; run(); } }, label)));
    input.addEventListener("keydown", (e) => { if (e.key === "Enter") run(); });
    const root = el("div", {}, [
      el("div", { class: "flex" }, [input, el("button", { class: "btn btn-primary", type: "button", onclick: run }, "Search")]),
      examples, out,
    ]);
    root.run = run; root.input = input;
    return root;
  }
  Atlas.patternSearchPanel = patternSearchPanel;

  // ---------------------------------------------------------------- PROSITE panel
  function prositePanel(A, options) {
    options = options || {};
    const root = el("div", {});
    let exclude = false, ticket = 0, alive = true;
    Atlas.onCleanup(() => { alive = false; });
    async function scan() {
      const mine = ++ticket;
      clear(root); root.appendChild(states.loadingBlock("Scanning the sequence against the PROSITE pattern database\u2026"));
      let result;
      try { result = await A.scanFn(exclude); }
      catch (error) { if (mine === ticket && alive) { clear(root); root.appendChild(states.errorBlock(error, scan)); } return; }
      if (mine !== ticket || !alive) return;
      draw(result);
    }
    function draw(r) {
      clear(root);
      if (r.status === "loading" || r.status === "idle") {
        root.appendChild(states.loadingBlock("The PROSITE database is downloading for the first time (about 25 MB). No scan has been run yet \u2014 this page will scan automatically when it is ready."));
        A.poll = setTimeout(() => { if (alive && document.body.contains(root)) scan(); }, 4000);
        Atlas.onCleanup(() => clearTimeout(A.poll));
        return;
      }
      if (r.status === "database_unavailable") {
        root.appendChild(states.stateBlock("\u26a0", "PROSITE database unavailable", (r.message || "") + " No scan was performed \u2014 this is not the same as \u201cno motifs found\u201d.",
          el("button", { class: "btn", type: "button", onclick: async () => { try { await api("/prosite/retry", { method: "POST", body: {} }); } catch (e) { /* shown by scan */ } scan(); } }, "Retry download")));
        return;
      }
      if (r.status === "failed") { root.appendChild(states.stateBlock("\u26a0", "The PROSITE scan failed", r.message || null, el("button", { class: "btn", type: "button", onclick: scan }, "Try again"))); return; }
      if (r.status !== "completed") { root.appendChild(states.noDataBlock(r.message || "Data unavailable")); return; }
      const db = r.database || {};
      root.appendChild(el("div", { class: "grid grid-4" }, [
        tile("Patterns tested", fmt.int(r.patterns_tested)), tile("Patterns matched", fmt.int(r.patterns_matched)),
        tile("Total occurrences", fmt.int(r.total_hits)), tile("PROSITE release", fmt.text(db.release)),
      ]));
      root.appendChild(el("p", { class: "panel-note", style: "margin-top:10px" }, [
        el("span", { class: "badge badge-ok" }, "Scan completed"), " " + (r.message || ""),
        " Source: " + fmt.text(db.source) + " \u00b7 " + fmt.int(db.patterns_usable) + " scannable patterns" + (db.patterns_unsupported ? " (" + fmt.int(db.patterns_unsupported) + " use syntax this scanner does not support)" : "") + ".",
      ]));
      root.appendChild(el("label", { class: "flex", style: "margin:10px 0" }, [
        el("input", { type: "checkbox", checked: exclude, style: "width:auto", onchange: (e) => { exclude = e.target.checked; scan(); } }),
        "Hide frequently-matching patterns (PROSITE \u201cskip\u201d flag \u2014 usually not informative)"]));
      if (!r.matches.length) {
        root.appendChild(states.stateBlock("\u2713", "No PROSITE pattern matches this sequence",
          "Scan completed: " + fmt.int(r.patterns_tested) + " patterns were tested against " + fmt.int(r.sequence_length) + " residues and none matched."));
        return;
      }
      const rows = r.matches;
      root.appendChild(table([
        { label: "Accession", render: (m) => el("a", { href: m.url, target: "_blank", rel: "noopener" }, m.accession) },
        { label: "Pattern name", render: (m) => m.id },
        { label: "Description", render: (m) => m.description },
        { label: "PROSITE pattern", align: "mono", render: (m) => m.pattern },
        { label: "Hits", align: "num", render: (m) => String(m.hits.length) },
        { label: "Start \u2013 End", render: (m) => hitsText(m.hits) },
        { label: "Matched sequence", align: "mono", render: (m) => m.hits.slice(0, 3).map((h) => h.sequence).join(", ") },
      ], rows.slice(0, options.compact ? 12 : 300), { id: "prosite-results" }));
      if (options.compact && rows.length > 12) root.appendChild(el("p", { class: "muted" }, "Showing 12 of " + rows.length + " matching patterns. The Function Scan tab lists all of them."));
      root.appendChild(el("h4", { class: "mt" }, "Where the motifs sit"));
      root.appendChild(motifTrack(r.sequence_length, rows.slice(0, 24).map((m, i) => ({
        label: m.id, title: m.accession + " \u00b7 " + m.description, hits: m.hits, color: [COLORS.teal, COLORS.amber, COLORS.violet, COLORS.coral, COLORS.blue][i % 5] }))));
      A.lastScan = r;
    }
    root.scan = scan;
    scan();
    return root;
  }

  // ---------------------------------------------------------------- tabs
  function tabOverview(panel, A) {
    const an = A.analysis, comp = an.composition, phys = an.physicochemical;
    panel.appendChild(el("p", { class: "panel-note" }, [el("span", { class: "badge badge-ok" }, "Validated"), " " + fmt.text(an.sequence.validation.message)]));
    const mw = ok(phys) && ok(phys.molecular_weight) ? phys.molecular_weight : null;
    panel.appendChild(el("div", { class: "grid grid-4" }, [
      tile("Length", fmt.int(an.sequence.length), "aa"),
      tile("Molecular weight", mw ? fmt.n(mw.value_kda, 2) : fmt.UNAVAILABLE || Atlas.UNAVAILABLE, mw ? "kDa" : ""),
      tile("Isoelectric point", ok(phys) && ok(phys.isoelectric_point) ? fmt.n(phys.isoelectric_point.value, 2) : Atlas.UNAVAILABLE, ""),
      tile("Net charge at pH 7", ok(phys) ? fmt.signed(phys.net_charge_ph7, 2) : Atlas.UNAVAILABLE, "e"),
      tile("GRAVY", ok(phys) ? fmt.n(phys.gravy, 3) : Atlas.UNAVAILABLE, ""),
      tile("Aromaticity", ok(phys) ? fmt.pct(phys.aromaticity_percent, 2) : Atlas.UNAVAILABLE, ""),
      tile("Aliphatic index", ok(phys) ? fmt.n(phys.aliphatic_index, 1) : Atlas.UNAVAILABLE, ""),
      tile("Residue types present", ok(comp) ? fmt.int(comp.distinct_standard_residues) + " / 20" : Atlas.UNAVAILABLE, ""),
    ]));
    if (A.protein && (A.protein.function || (A.protein.keywords && A.protein.keywords.length))) {
      panel.appendChild(el("div", { class: "mt" }, card("Annotated function (UniProt)", null, [
        A.protein.function ? el("p", {}, A.protein.function) : null,
        A.protein.keywords && A.protein.keywords.length ? el("div", { class: "flex flex-wrap" }, A.protein.keywords.map((k) => el("span", { class: "badge" }, k))) : null])));
    }
    panel.appendChild(el("div", { class: "mt" }, guard(comp, (c) => {
      const abundant = (title, block, isLeast) => card(title, null, block ? [
        ...block.residues.map((r) => kv([["Residue", residueLabel(r)]])),
        kv([["Count", fmt.int(block.count)], ["Share of sequence", fmt.pct(block.percent, 2)]]),
        block.residues.length > 1 ? el("p", { class: "panel-note" }, block.residues.length + " residues are tied at this count.") : null,
        isLeast && block.only_one_type_present ? el("p", { class: "panel-note" }, "Only one residue type occurs in this sequence, so it is both the most and the least abundant.") : null,
        isLeast && c.absent.length ? el("p", { class: "panel-note" }, c.absent.length + " of the 20 standard residues do not occur: " + c.absent.join(" ") + ".") : null,
      ] : states.noDataBlock("No standard residues are present."));
      return el("div", { class: "grid grid-2" }, [abundant("Most abundant residue", c.most_abundant, false), abundant("Least abundant residue (among those present)", c.least_abundant, true)]);
    })));
    panel.appendChild(el("div", { class: "mt grid grid-2" }, [
      card("Chemical family breakdown", null, guard(comp, (c) => el("div", {}, [
        table([
          { label: "Family", render: (f) => Atlas.classPill(f.family) }, { label: "Residues", align: "mono", render: (f) => f.residues },
          { label: "Count", align: "num", render: (f) => fmt.int(f.count) }, { label: "Share", align: "num", render: (f) => fmt.pct(f.percent, 2) },
        ], c.families, { id: "family-table" }),
        el("p", { class: "panel-note", style: "margin-top:10px" }, fmt.text(c.classification_note)),
      ]))),
      card("Family composition", "Hover a slice for its exact share.", (function () {
        const holder = el("div", {});
        if (ok(comp)) charts.doughnut(holder, { labels: comp.families.map((f) => f.family), datasets: [{ data: comp.families.map((f) => f.count), backgroundColor: comp.families.map((f) => Atlas.CLASS_COLOR[f.family] || "#74897f"), borderWidth: 0 }], tall: true, ariaLabel: "Family composition" });
        else holder.appendChild(states.noDataBlock("Data unavailable"));
        return holder;
      })()),
    ]));
    panel.appendChild(el("div", { class: "mt" }, card("Amino acid composition", "Percentage of the sequence made up by each residue.", (function () {
      const holder = el("div", {});
      if (ok(comp)) compositionChart(holder, comp, "percent");
      else holder.appendChild(states.noDataBlock("Data unavailable"));
      return holder;
    })())));
  }
  function compositionChart(holder, comp, field) {
    const rows = comp.rows;
    return charts.bar(holder, {
      labels: rows.map((r) => r.code), legend: false, tall: true,
      datasets: [{ label: field === "count" ? "Count" : "Percent of sequence", data: rows.map((r) => r[field]), backgroundColor: rows.map((r) => Atlas.CLASS_COLOR[r.family] || COLORS.dim), borderWidth: 0 }],
      xTitle: "Residue", yTitle: field === "count" ? "Residues" : "% of sequence", ariaLabel: "Amino acid composition",
      tooltip: { title: (items) => { const r = rows[items[0].dataIndex]; return r ? r.name + " (" + r.three + ")" : ""; },
        label: (item) => { const r = rows[item.dataIndex]; return r ? fmt.int(r.count) + " residues \u00b7 " + fmt.pct(r.percent, 2) : ""; } },
    });
  }

  function tabPhysico(panel, A) {
    const phys = A.analysis.physicochemical;
    panel.appendChild(guard(phys, (p) => {
      const root = el("div", {});
      root.appendChild(guard(p.molecular_weight, (mw) => card("Molecular weight", mw.convention, el("div", { class: "grid grid-4" }, [
        tile("Weight", fmt.n(mw.value_da, 2), "Da"), tile("In kilodaltons", fmt.n(mw.value_kda, 3), "kDa"),
        tile("Residues counted", fmt.int(mw.residues_counted)), tile("Approximate-mass residues", fmt.int(mw.approximated_residues), "", "X, B, Z, J use average masses"),
      ]))));
      root.appendChild(el("div", { class: "mt" }, guard(p.mass_breakdown, (mb) => card("Molecular-weight breakdown by residue", "Residue mass = free amino acid mass \u2212 one water. Contribution = count \u00d7 residue mass.", el("div", {}, [
        table([
          { label: "Residue", render: (r) => r.name + " (" + r.code + ")" }, { label: "Count", align: "num", render: (r) => fmt.int(r.count) },
          { label: "Residue mass (Da)", align: "num", render: (r) => fmt.n(r.residue_mass_da, 2) }, { label: "Total mass (Da)", align: "num", render: (r) => fmt.n(r.mass_da, 2) },
          { label: "Share", align: "num", render: (r) => fmt.pct(r.percent, 2) },
        ], mb.rows, { id: "mass-table" }),
        el("div", { class: "reconcile" }, [
          kv([["Sum of residue contributions", fmt.n(mb.residue_sum_da, 2) + " Da"], ["+ one water (termini)", fmt.n(mb.water_da, 3) + " Da"], ["= Total", fmt.n(mb.total_da, 2) + " Da"]]),
          el("span", { class: "badge " + (mb.reconciles ? "badge-ok" : "badge-warn") }, mb.reconciles ? "Reconciles with the molecular weight above" : "Difference: " + fmt.n(mb.difference_da, 4) + " Da"),
        ]),
        (function () { const h = el("div", { class: "mt" }); charts.bar(h, { labels: mb.rows.map((r) => r.code), legend: false, datasets: [{ label: "Mass contribution (Da)", data: mb.rows.map((r) => r.mass_da), backgroundColor: COLORS.amber, borderWidth: 0 }], xTitle: "Residue", yTitle: "Mass contribution (Da)", ariaLabel: "Mass contribution per residue",
          tooltip: { label: (item) => { const r = mb.rows[item.dataIndex]; return r ? r.count + " \u00d7 " + fmt.n(r.residue_mass_da, 2) + " = " + fmt.n(r.mass_da, 2) + " Da" : ""; } } }); return h; })(),
      ])))));
      root.appendChild(el("div", { class: "mt grid grid-3" }, [
        tile("Isoelectric point", ok(p.isoelectric_point) ? fmt.n(p.isoelectric_point.value, 2) : Atlas.UNAVAILABLE, "pH", p.isoelectric_point && p.isoelectric_point.method),
        tile("Net charge at pH 7", fmt.signed(p.net_charge_ph7, 2), "e"), tile("Aromaticity", fmt.pct(p.aromaticity_percent, 2)),
        tile("Aliphatic index", fmt.n(p.aliphatic_index, 1)), tile("GRAVY", fmt.n(p.gravy, 3)),
        tile("Extinction coeff. (reduced)", fmt.int(p.extinction_coefficient && p.extinction_coefficient.reduced), "M\u207b\u00b9cm\u207b\u00b9"),
      ]));
      root.appendChild(el("div", { class: "mt grid grid-2" }, [
        card("Net charge vs pH", "Henderson\u2013Hasselbalch titration using the " + fmt.text(p.isoelectric_point && p.isoelectric_point.pka_set) + " pKa set. The dashed line marks the isoelectric point.", (function () {
          const h = el("div", {});
          const curve = p.charge_curve || [];
          const lo = Math.min.apply(null, curve.map((c) => c.charge)), hi = Math.max.apply(null, curve.map((c) => c.charge));
          const datasets = [{ label: "Net charge", data: curve.map((c) => ({ x: c.ph, y: c.charge })), borderColor: COLORS.teal, backgroundColor: COLORS.teal, pointRadius: 0, borderWidth: 2, tension: 0.2 },
            { label: "Zero charge", data: [{ x: 0, y: 0 }, { x: 14, y: 0 }], borderColor: "rgba(255,255,255,0.3)", borderDash: [4, 4], pointRadius: 0, borderWidth: 1 }];
          if (ok(p.isoelectric_point) && Number.isFinite(lo) && Number.isFinite(hi)) datasets.push({ label: "pI " + fmt.n(p.isoelectric_point.value, 2), data: [{ x: p.isoelectric_point.value, y: lo }, { x: p.isoelectric_point.value, y: hi }], borderColor: COLORS.amber, borderDash: [6, 4], pointRadius: 0, borderWidth: 1.5 });
          charts.line(h, { datasets, xTitle: "pH", yTitle: "Net charge (e)", xRange: [0, 14], tall: true, ariaLabel: "Net charge versus pH" });
          return h;
        })()),
        card("Charged residues", null, (function () {
          const cc = p.charge_composition || {};
          return el("div", {}, [
            table([{ label: "Group" }, { label: "Residues", align: "mono" }, { label: "Count", align: "num" }, { label: "Share", align: "num" }].map((c, i) => Object.assign(c, { render: (r) => r[i] })), [
              ["Acidic", cc.acidic_residues, fmt.int(cc.acidic_count), fmt.pct(cc.acidic_percent, 2)], ["Basic", cc.basic_residues, fmt.int(cc.basic_count), fmt.pct(cc.basic_percent, 2)],
              ["Total charged", "", fmt.int(cc.charged_count), fmt.pct(cc.charged_percent, 2)], ["Uncharged", "", fmt.int(cc.uncharged_count), fmt.pct(cc.uncharged_percent, 2)]], { id: "charge-table" }),
            el("p", { class: "panel-note", style: "margin-top:10px" }, "Extinction coefficient at 280 nm: " + fmt.int(p.extinction_coefficient && p.extinction_coefficient.reduced) + " (reduced), " + fmt.int(p.extinction_coefficient && p.extinction_coefficient.oxidised) + " (" + fmt.int(p.extinction_coefficient && p.extinction_coefficient.cystines) + " cystine pair(s) assumed)."),
          ]);
        })()),
      ]));
      return root;
    }));
  }

  function repeatsCard(A) {
    const feat = A.analysis.features;
    return card("Repeats and low-complexity runs", null, guard(feat, (f) => guard(f.repeats, (r) => el("div", {}, [
      el("p", { class: "panel-note" }, fmt.int(r.total_windows) + " overlapping " + r.size + "-mers, " + fmt.int(r.distinct_windows) + " distinct."),
      table([{ label: r.size + "-mer", align: "mono", render: (x) => x.unit }, { label: "Occurrences", align: "num", render: (x) => fmt.int(x.occurrences) }], r.top_kmers, { id: "kmer-table" }),
      r.longest_run ? el("p", { class: "panel-note", style: "margin-top:10px" }, "Longest single-residue run: " + r.longest_run.length + " \u00d7 " + r.longest_run.residue + " at positions " + r.longest_run.start + "\u2013" + r.longest_run.end + ".") : null,
    ]))));
  }
  function tabFeatures(panel, A) {
    panel.appendChild(card("Custom pattern search", "Search this sequence for a literal fragment (e.g. GGG) or a PROSITE-style pattern. Overlapping matches are all reported.", patternSearchPanel(() => seqOf(A))));
    panel.appendChild(el("div", { class: "mt" }, card("Known motifs (PROSITE)", "Every pattern in the PROSITE database is tested against this sequence.", prositePanel(A, { compact: true }))));
    panel.appendChild(el("div", { class: "mt" }, repeatsCard(A)));
  }

  function tabAminoAcids(panel, A) {
    const comp = A.analysis.composition;
    panel.appendChild(guard(comp, (c) => {
      const root = el("div", {});
      const chartHolder = el("div", {});
      let field = "percent";
      const seg = el("div", { class: "seg-group", style: "width:max-content;margin-bottom:12px" }, [["percent", "Percent"], ["count", "Count"]].map(([v, l]) =>
        el("button", { class: "seg" + (v === field ? " active" : ""), type: "button", onclick: (e) => { field = v; Array.from(seg.children).forEach((b) => b.classList.toggle("active", b === e.currentTarget)); clear(chartHolder); compositionChart(chartHolder, c, field); } }, l)));
      root.appendChild(card("Composition of the 20 standard amino acids", "Coloured by side-chain class. Zoom with the wheel over the chart; drag to pan.", [seg, chartHolder]));
      compositionChart(chartHolder, c, field);
      const max = Math.max.apply(null, c.rows.map((r) => r.count)) || 1;
      root.appendChild(el("div", { class: "mt" }, card("Residue table", "Count and percentage for every residue, computed from this sequence.", table([
        { label: "Residue", render: (r) => r.name }, { label: "Code", align: "mono", render: (r) => r.three + " \u00b7 " + r.code },
        { label: "Class", render: (r) => Atlas.classPill(r.family) }, { label: "Count", align: "num", render: (r) => fmt.int(r.count) },
        { label: "Percentage", align: "num", render: (r) => fmt.pct(r.percent, 2) },
        { label: "", render: (r) => el("div", { class: "bar-cell" }, el("div", { class: "bar-fill", style: "width:" + (r.count / max * 100).toFixed(1) + "%;background:" + (Atlas.CLASS_COLOR[r.family] || COLORS.dim) })) },
      ], c.rows, { id: "aa-table" }))));
      // enrichment against the complete dataset
      const enrich = el("div", {}, states.loadingBlock("Comparing with the dataset statistics\u2026"));
      api("/dataset/stats?scope=human_reviewed").then((st) => {
        clear(enrich);
        if (st.state !== "ready") { enrich.appendChild(states.noDataBlock(st.state === "building" ? "The dataset statistics are still being built (first run only) \u2014 reopen this tab in a minute." : "Dataset statistics are unavailable: " + fmt.text(st.error))); return; }
        const ds = {}; st.stats.amino_acids.forEach((r) => { ds[r.code] = r.percent; });
        const rows = c.rows.filter((r) => ds[r.code] !== undefined).map((r) => ({ code: r.code, name: r.name, mine: r.percent, ref: ds[r.code], diff: r.percent - ds[r.code], fold: ds[r.code] > 0 ? r.percent / ds[r.code] : null }));
        const holder = el("div", {});
        charts.bar(holder, { labels: rows.map((r) => r.code), legend: false, tall: true, xTitle: "Residue", yTitle: "This sequence minus dataset (percentage points)", ariaLabel: "Enrichment versus dataset",
          datasets: [{ label: "Difference", data: rows.map((r) => r.diff), backgroundColor: rows.map((r) => (r.diff >= 0 ? COLORS.teal : COLORS.coral)), borderWidth: 0 }], beginAtZero: false,
          tooltip: { label: (item) => { const r = rows[item.dataIndex]; return r ? fmt.pct(r.mine, 2) + " here vs " + fmt.pct(r.ref, 2) + " in the dataset (" + (r.fold === null ? "n/a" : fmt.n(r.fold, 2) + "\u00d7") + ")" : ""; } } });
        enrich.appendChild(holder);
        enrich.appendChild(el("div", { class: "mt" }, table([{ label: "Residue", render: (r) => r.name + " (" + r.code + ")" }, { label: "This sequence", align: "num", render: (r) => fmt.pct(r.mine, 2) },
          { label: "Dataset", align: "num", render: (r) => fmt.pct(r.ref, 3) }, { label: "Enrichment", align: "num", render: (r) => (r.fold === null ? Atlas.UNAVAILABLE : fmt.n(r.fold, 2) + "\u00d7") }], rows, { id: "enrichment-table" })));
        enrich.appendChild(el("p", { class: "panel-note", style: "margin-top:10px" }, "Dataset: " + fmt.text(st.label) + " (" + fmt.int(st.stats.total_proteins) + " proteins). Enrichment = share in this sequence \u00f7 share in the dataset."));
      }).catch((error) => { clear(enrich); enrich.appendChild(states.errorBlock(error, null)); });
      root.appendChild(el("div", { class: "mt" }, card("Compared with the dataset", "Which residues this protein has more (or less) of than the average protein in the reference dataset.", enrich)));
      // residue position finder
      const present = c.rows.filter((r) => r.count > 0);
      const positions = (A.analysis.features && A.analysis.features.residue_positions) || {};
      const out = el("div", { class: "mt" });
      const select = el("select", { "aria-label": "Residue" }, present.map((r) => el("option", { value: r.code }, r.name + " (" + r.code + ") \u2014 " + r.count)));
      function show() {
        clear(out);
        const code = select.value, list = positions[code] || [];
        const hits = list.map((p) => ({ start: p, end: p }));
        out.appendChild(el("p", { class: "panel-note" }, fmt.int(list.length) + " occurrence(s) of " + code + (list.length ? ": " + list.slice(0, 80).join(", ") + (list.length > 80 ? " \u2026" : "") : "")));
        if (list.length) { out.appendChild(motifTrack(A.analysis.sequence.length, [{ label: code, hits, color: Atlas.CLASS_COLOR[(c.rows.find((r) => r.code === code) || {}).family] || COLORS.teal }])); out.appendChild(sequenceViewer(seqOf(A), { hits })); }
      }
      select.addEventListener("change", show);
      root.appendChild(el("div", { class: "mt" }, card("Residue position finder", "Choose a residue to see where it sits in the sequence.", [el("div", { style: "max-width:340px" }, select), out])));
      if (present.length) show();
      return root;
    }));
  }

  function tabHydropathy(panel, A) {
    const holder = el("div", {});
    panel.appendChild(holder);
    function draw() {
      clear(holder);
      const hyd = A.analysis.hydropathy;
      holder.appendChild(guard(hyd, (h) => {
        const root = el("div", {});
        const input = el("input", { type: "text", value: String(h.window), style: "width:90px", "aria-label": "Window size" });
        const apply = async () => {
          const w = parseInt(input.value, 10);
          if (!Number.isFinite(w) || w < 3 || w > 101) { Atlas.toast("Window size must be a whole number from 3 to 101.", "error"); return; }
          if (w === h.window) return;
          try { await A.setWindow(w); draw(); } catch (error) { Atlas.toast(error.message, "error"); }
        };
        input.addEventListener("keydown", (e) => { if (e.key === "Enter") apply(); });
        const segs = h.membrane_segments || [];
        root.appendChild(el("div", { class: "grid grid-4" }, [
          tile("GRAVY score", fmt.n(h.gravy, 3), "", "Grand average of hydropathy (Kyte\u2013Doolittle)"), tile("Window size", fmt.int(h.window), "aa"),
          tile("Membrane-like segments", h.membrane_status === "success" ? fmt.int(segs.length) : Atlas.UNAVAILABLE),
          tile("Scale", fmt.text(h.scale)),
        ]));
        const chartCard = card("Kyte\u2013Doolittle hydropathy profile", "Mean hydropathy in a sliding window (positive = hydrophobic). Shaded bands are membrane-like segments.", [
          el("div", { class: "flex", style: "margin-bottom:12px" }, [el("label", { style: "margin:0" }, "Window (3\u2013101)"), input, el("button", { class: "btn", type: "button", onclick: apply }, "Apply")]),
        ]);
        const box = el("div", {}); chartCard.appendChild(box);
        if (h.profile && h.profile.positions.length) {
          const pts = h.profile.positions.map((x, i) => ({ x, y: h.profile.values[i] }));
          const last = h.profile.positions[h.profile.positions.length - 1];
          charts.line(box, {
            datasets: [{ label: "Hydropathy (window " + h.window + ")", data: pts, borderColor: COLORS.teal, backgroundColor: COLORS.teal, pointRadius: 0, borderWidth: 1.6, tension: 0 },
              { label: "Membrane threshold " + h.membrane_threshold, data: [{ x: 1, y: h.membrane_threshold }, { x: A.analysis.sequence.length, y: h.membrane_threshold }], borderColor: COLORS.coral, borderDash: [6, 4], pointRadius: 0, borderWidth: 1 },
              { label: "Zero", data: [{ x: 1, y: 0 }, { x: A.analysis.sequence.length, y: 0 }], borderColor: "rgba(255,255,255,0.25)", borderDash: [3, 3], pointRadius: 0, borderWidth: 1 }],
            xRange: [1, Math.max(A.analysis.sequence.length, last)], xTitle: "Residue position (window centre)", yTitle: "Mean hydropathy", tall: true, ariaLabel: "Hydropathy profile",
            bands: segs.map((s) => ({ from: s.start, to: s.end, color: "rgba(226,114,91,0.16)" })),
          });
        } else if (h.per_residue && h.per_residue.length) {
          box.appendChild(el("p", { class: "panel-note" }, fmt.text(h.message)));
          charts.bar(box, { labels: seqOf(A).split(""), legend: false, tall: true, datasets: [{ label: "Hydropathy", data: h.per_residue, backgroundColor: h.per_residue.map((v) => (v >= 0 ? COLORS.teal : COLORS.blue)), borderWidth: 0 }], xTitle: "Residue", yTitle: "Kyte\u2013Doolittle value", ariaLabel: "Per-residue hydropathy" });
        } else box.appendChild(states.noDataBlock(h.message || "Data unavailable"));
        root.appendChild(el("div", { class: "mt" }, chartCard));
        root.appendChild(el("div", { class: "mt" }, card("Membrane-like segments", "Windows of " + fmt.int(h.membrane_window) + " residues averaging above " + h.membrane_threshold + " (a simple heuristic, not a trained predictor).",
          h.membrane_status !== "success" ? states.noDataBlock(h.membrane_message || "Data unavailable")
            : !segs.length ? el("p", { class: "panel-note" }, [el("span", { class: "badge badge-ok" }, "Analysis completed"), " No membrane-like segment was found."])
              : table([{ label: "#", align: "num", render: (s, i) => String(i + 1) }, { label: "Start", align: "num", render: (s) => String(s.start) }, { label: "End", align: "num", render: (s) => String(s.end) },
                { label: "Length", align: "num", render: (s) => fmt.int(s.length) }, { label: "Mean hydropathy", align: "num", render: (s) => fmt.n(s.mean, 3) }], segs, { id: "membrane-table" }))));
        return root;
      }));
    }
    draw();
  }

  function tabFunction(panel, A) {
    panel.appendChild(card("Protein function analysis (PROSITE scan)", "Function is inferred only from patterns that genuinely match. The scan status below tells you exactly what was run.", prositePanel(A, {})));
    if (A.protein && (A.protein.function || (A.protein.keywords || []).length)) {
      panel.appendChild(el("div", { class: "mt" }, card("Annotated function (UniProt)", "Curated annotation, independent of the pattern scan.", [
        A.protein.function ? el("p", {}, A.protein.function) : null,
        (A.protein.keywords || []).length ? el("div", { class: "flex flex-wrap" }, A.protein.keywords.map((k) => el("span", { class: "badge" }, k))) : null])));
    }
  }

  function fastaOf(A) {
    const seq = seqOf(A), lines = [];
    for (let i = 0; i < seq.length; i += 60) lines.push(seq.slice(i, i + 60));
    const header = A.protein ? (A.protein.header || (A.protein.id + " " + (A.protein.name || ""))) : (A.label || "sequence");
    return ">" + header + "\n" + lines.join("\n") + "\n";
  }
  function tabSequence(panel, A) {
    panel.appendChild(card("Full sequence", "Numbered every 60 residues, grouped in blocks of 10.", [
      el("div", { class: "flex flex-wrap", style: "margin-bottom:12px" }, [
        el("button", { class: "btn", type: "button", onclick: () => Atlas.copyText(seqOf(A)) }, "Copy sequence"),
        el("button", { class: "btn", type: "button", onclick: () => Atlas.copyText(fastaOf(A)) }, "Copy as FASTA"),
        el("button", { class: "btn", type: "button", onclick: () => Atlas.download(A.fileBase + ".fasta", fastaOf(A)) }, "Download FASTA")]),
      sequenceViewer(seqOf(A), {}),
    ]));
  }

  // ---------------------------------------------------------------- export
  const csvCell = (v) => { const s = v === null || v === undefined ? "" : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
  const csv = (header, rows) => [header.map(csvCell).join(",")].concat(rows.map((r) => r.map(csvCell).join(","))).join("\n") + "\n";
  Atlas.csv = csv;
  const exporters = {
    fasta: { label: "FASTA sequence", ext: "fasta", make: (A) => fastaOf(A) },
    composition: { label: "Amino-acid composition (CSV)", ext: "composition.csv", make: (A) => csv(["Code", "Amino acid", "Family", "Count", "Percent"], A.analysis.composition.rows.map((r) => [r.code, r.name, r.family, r.count, fmt.plain(r.percent, 4)])) },
    mass: { label: "Molecular-weight breakdown (CSV)", ext: "mass.csv", make: (A) => csv(["Code", "Amino acid", "Count", "Residue mass (Da)", "Total mass (Da)", "Percent"], A.analysis.physicochemical.mass_breakdown.rows.map((r) => [r.code, r.name, r.count, fmt.plain(r.residue_mass_da, 3), fmt.plain(r.mass_da, 3), fmt.plain(r.percent, 3)])) },
    hydropathy: { label: "Hydropathy profile (CSV)", ext: "hydropathy.csv", make: (A) => { const h = A.analysis.hydropathy, seq = seqOf(A); const win = {}; if (h.profile) h.profile.positions.forEach((p, i) => { win[p] = h.profile.values[i]; });
      return csv(["Position", "Residue", "Residue hydropathy", "Window mean (window " + h.window + ")"], seq.split("").map((c, i) => [i + 1, c, h.per_residue ? h.per_residue[i] : "", win[i + 1] !== undefined ? fmt.plain(win[i + 1], 4) : ""])); } },
    prosite: { label: "PROSITE matches (CSV)", ext: "prosite.csv", needsScan: true, make: (A, scan) => csv(["Accession", "Name", "Description", "Pattern", "Start", "End", "Sequence"], (scan.matches || []).reduce((acc, m) => acc.concat(m.hits.map((h) => [m.accession, m.id, m.description, m.pattern, h.start, h.end, h.sequence])), [])) },
    json: { label: "Complete JSON report", ext: "report.json", needsScan: true, make: (A, scan) => JSON.stringify({ protein: A.protein || null, analysis: A.analysis, prosite: scan || null }, null, 2) },
    text: { label: "Text report", ext: "report.txt", make: (A) => textReport(A) },
  };
  function textReport(A) {
    const an = A.analysis, comp = an.composition, phys = an.physicochemical, hyd = an.hydropathy, L = [];
    L.push("PROTEIN ATLAS REPORT", "====================", A.protein ? A.protein.id + " - " + (A.protein.name || "") : (A.label || "Sequence"), "");
    L.push("Length: " + fmt.int(an.sequence.length) + " aa", "Validation: " + fmt.text(an.sequence.validation.message), "");
    if (ok(phys)) {
      L.push("PHYSICOCHEMICAL", "Molecular weight: " + fmt.da(phys.molecular_weight.value_da), "Isoelectric point: " + fmt.n(phys.isoelectric_point.value, 2),
        "Net charge at pH 7: " + fmt.signed(phys.net_charge_ph7, 2), "GRAVY: " + fmt.n(phys.gravy, 3), "Aromaticity: " + fmt.pct(phys.aromaticity_percent, 2), "");
    }
    if (ok(comp)) {
      L.push("MOST ABUNDANT: " + comp.most_abundant.residues.map((r) => r.name).join(", ") + " (" + comp.most_abundant.count + ")",
        "LEAST ABUNDANT: " + comp.least_abundant.residues.map((r) => r.name).join(", ") + " (" + comp.least_abundant.count + ")", "", "COMPOSITION");
      comp.rows.forEach((r) => L.push("  " + r.code + "  " + r.name.padEnd(14) + String(r.count).padStart(6) + "  " + fmt.pct(r.percent, 2)));
    }
    if (ok(hyd)) L.push("", "HYDROPATHY (window " + hyd.window + "): GRAVY " + fmt.n(hyd.gravy, 3) + ", membrane-like segments: " + (hyd.membrane_segments || []).length);
    return L.join("\n") + "\n";
  }
  Atlas.exporters = exporters;
  function tabExport(panel, A) {
    panel.appendChild(card("Export analysis", "Files are generated in your browser from the analysis already on screen \u2014 nothing is re-computed differently.", el("div", { class: "grid grid-2" },
      Object.keys(exporters).map((key) => {
        const ex = exporters[key];
        return el("div", { class: "export-item" }, [el("div", {}, [el("b", {}, ex.label), el("div", { class: "muted" }, A.fileBase + "." + ex.ext)]),
          el("button", { class: "btn", type: "button", onclick: async () => {
            try {
              let scan = null;
              if (ex.needsScan) { try { scan = A.lastScan || await A.scanFn(false); } catch (e) { scan = { status: "unavailable", message: e.message }; } }
              Atlas.download(A.fileBase + "." + ex.ext, ex.make(A, scan), ex.ext.endsWith("json") ? "application/json" : "text/plain");
              Atlas.toast("Saved " + A.fileBase + "." + ex.ext);
            } catch (error) { Atlas.toast("That export is unavailable for this sequence.", "error"); }
          } }, "Download")]);
      }))));
  }

  // ---------------------------------------------------------------- the tabbed workspace
  const TABS = [
    ["overview", "Overview", tabOverview], ["physico", "Physicochemical", tabPhysico], ["features", "Sequence Features", tabFeatures],
    ["aminoacids", "Amino Acids", tabAminoAcids], ["hydropathy", "Hydropathy", tabHydropathy], ["function", "Function Scan", tabFunction],
    ["sequence", "Full Sequence", tabSequence], ["export", "Export", tabExport],
  ];
  Atlas.TAB_IDS = TABS.map((t) => t[0]);
  /* A = { analysis, protein?, label, fileBase, scanFn(exclude), setWindow(w) } */
  function analysisTabs(A, options) {
    options = options || {};
    const tabs = el("div", { class: "tabs", role: "tablist" });
    const panels = el("div", { class: "tab-panels" });
    const built = {};
    TABS.forEach(([id, label]) => {
      tabs.appendChild(el("button", { class: "tab-btn", type: "button", role: "tab", "data-tab": id, onclick: () => { if (options.onTab) options.onTab(id); else activate(id); } }, label));
      panels.appendChild(el("div", { class: "tab-panel", id: "tab-" + id, role: "tabpanel" }));
    });
    function activate(id) {
      if (!TABS.some((t) => t[0] === id)) id = "overview";
      Atlas.$$(".tab-btn", tabs).forEach((b) => b.classList.toggle("active", b.getAttribute("data-tab") === id));
      Atlas.$$(".tab-panel", panels).forEach((p) => p.classList.toggle("active", p.id === "tab-" + id));
      if (!built[id]) {
        built[id] = true;
        const panel = Atlas.$("#tab-" + id, panels);
        try { TABS.find((t) => t[0] === id)[2](panel, A); }
        catch (error) { Atlas.diagnostics.errors.push("tab " + id + ": " + (error && error.stack || error)); clear(panel); panel.appendChild(states.stateBlock("\u26a0", "Analysis unavailable \u2014 calculation error", "This section could not be drawn.")); }
        Atlas.scrub(panel);
      }
      A.activeTab = id;
    }
    A.invalidate = () => { Object.keys(built).forEach((k) => { built[k] = false; clear(Atlas.$("#tab-" + k, panels)); }); activate(A.activeTab || "overview"); };
    const root = el("div", {}, [tabs, panels]);
    root.activate = activate; root.built = built;
    activate(options.tab || "overview");
    return root;
  }
  Atlas.analysisTabs = analysisTabs;
  Atlas.exportCard = (A) => { const p = el("div", {}); tabExport(p, A); return p; };
  Atlas.prositePanel = prositePanel;
  Atlas.helpers = { tile, card, kv, guard };
})(typeof window !== "undefined" ? window : globalThis);
