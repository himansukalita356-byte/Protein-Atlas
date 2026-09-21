/* Protein Atlas - pages (routes): Discover, Protein, Sequence Lab, Pattern Lab, Compare, Reference, Dataset, Export. */
(function (global) {
  "use strict";
  const Atlas = global.Atlas;
  const { el, $, clear, fmt, api, qs, states, charts, table, Router, Session, COLORS } = Atlas;
  const { tile, card, kv, guard } = Atlas.helpers;
  Atlas.pages = { protein: null };
  Atlas.cache = { reference: null };
  const view = () => $("#view");
  const head = (label, title, sub) => [el("div", { class: "section-label" }, label), el("h1", { class: "page-title" }, title), sub ? el("p", { class: "page-sub" }, sub) : null];
  function setPage(nodes) { const v = view(); clear(v); (Array.isArray(nodes) ? nodes : [nodes]).forEach((n) => n && v.appendChild(n)); return v; }
  const intParam = (v, d) => { const n = parseInt(v, 10); return Number.isFinite(n) && n > 0 ? n : d; };

  // ================================================================= Discover
  const QUICK = [["Human \u00b7 reviewed", "reviewed:true AND organism_id:9606"], ["TP53", "gene:TP53"], ["Insulin", "insulin AND reviewed:true"],
    ["Kinases", "keyword:Kinase AND reviewed:true AND organism_id:9606"], ["Membrane proteins", "keyword:Membrane AND reviewed:true AND organism_id:9606"], ["Mouse \u00b7 reviewed", "reviewed:true AND organism_id:10090"]];
  function pager(q, data) {
    const go = (page) => Router.go("/", { q, page, size: data.page_size });
    const p = data.page, n = data.total_pages;
    const jump = el("input", { type: "text", value: String(p), style: "width:70px;text-align:center", "aria-label": "Go to page number", inputmode: "numeric" });
    const doJump = () => { const t = parseInt(jump.value, 10); if (Number.isFinite(t) && t >= 1 && t <= n) go(t); else Atlas.toast("Enter a page from 1 to " + fmt.int(n) + ".", "error"); };
    jump.addEventListener("keydown", (e) => { if (e.key === "Enter") doJump(); });
    const size = el("select", { style: "width:auto", "aria-label": "Results per page", onchange: (e) => Router.go("/", { q, page: 1, size: e.target.value }) },
      [12, 24, 48, 96].map((s) => el("option", { value: String(s), selected: s === data.page_size }, s + " per page")));
    return el("nav", { class: "pager", "aria-label": "Pagination" }, [
      el("button", { class: "btn", type: "button", disabled: p <= 1, onclick: () => go(1) }, "\u00ab First"),
      el("button", { class: "btn", type: "button", disabled: p <= 1, onclick: () => go(p - 1) }, "\u2039 Previous"),
      el("span", {}, "Page"), jump, el("span", {}, "of " + fmt.int(n)),
      el("button", { class: "btn btn-ghost", type: "button", onclick: doJump }, "Go"),
      el("button", { class: "btn", type: "button", disabled: p >= n, onclick: () => go(p + 1) }, "Next \u203a"),
      el("button", { class: "btn", type: "button", disabled: p >= n, onclick: () => go(n) }, "Last \u00bb"), size]);
  }
  function proteinCard(p) {
    return el("a", { class: "protein-card", href: "#/protein/" + encodeURIComponent(p.id) }, [
      el("div", { class: "pc-top" }, [el("span", { class: "pc-id" }, p.id), el("span", { class: "pc-len" }, p.length ? fmt.int(p.length) + " aa" : "")]),
      el("div", { class: "pc-name" }, fmt.text(p.name)),
      el("div", { class: "pc-meta" }, [p.gene ? el("b", {}, p.gene) : null, p.gene ? " \u00b7 " : "", fmt.text(p.organism), p.reviewed === true ? " \u00b7 reviewed" : p.reviewed === false ? " \u00b7 unreviewed" : ""])]);
  }
  async function viewDiscover(ctx) {
    const q = ctx.query.q || "", page = intParam(ctx.query.page, 1), size = intParam(ctx.query.size, 24);
    const input = $("#global-search-input"); if (input) input.value = q;
    const results = el("div", {}, el("div", { class: "results-grid" }, Array.from({ length: 6 }, () => states.skeleton(96))));
    setPage([...head("DISCOVER", "Explore the protein universe", "Search UniProt live by accession, gene or keyword, then open any protein for full analysis, 3D structure and export."),
      el("div", { class: "flex flex-wrap", style: "margin-bottom:22px" }, QUICK.map(([label, query]) => el("a", { class: "btn btn-ghost", href: Router.href("/", { q: query }) }, label))), results]);
    let data;
    try { data = await api("/proteins" + qs({ query: q, page, page_size: size })); }
    catch (error) {
      if (!ctx.isCurrent()) return;
      clear(results);
      if (error.kind === "page_out_of_range" || (error.extra && error.extra.kind === "page_out_of_range")) {
        const last = error.extra.total_pages || 1;
        results.appendChild(states.stateBlock("\u21c4", "That page does not exist", error.message, el("a", { class: "btn", href: Router.href("/", { q, page: last, size }) }, "Go to the last page (" + fmt.int(last) + ")")));
      } else results.appendChild(states.errorBlock(error, () => Router.render()));
      if (error.kind === "page_unreachable") results.appendChild(el("a", { class: "btn", href: Router.href("/", { q, page: 1, size }) }, "Back to page 1"));
      return;
    }
    if (!ctx.isCurrent()) return;
    clear(results);
    const first = (data.page - 1) * data.page_size + 1, last = Math.min(data.total_items, first + data.items.length - 1);
    if (!data.items.length) { results.appendChild(states.stateBlock("\ud83d\udd0d", "No proteins matched", "Try a gene name (TP53), an accession (P04637) or a keyword (kinase).")); return; }
    results.appendChild(el("p", { class: "panel-note" }, [el("b", {}, fmt.int(first) + "\u2013" + fmt.int(last)), " of ", el("b", {}, fmt.int(data.total_items)), " proteins" + (q ? " for \u201c" + q + "\u201d" : " (default: reviewed human proteins)") + " \u00b7 ordered by " + fmt.text(data.ordering)]));
    results.appendChild(el("div", { class: "results-grid" }, data.items.map(proteinCard)));
    results.appendChild(pager(q, data));
  }

  // ================================================================= Protein page
  async function buildProteinA(id, win) {
    const [record, bundle] = await Promise.all([api("/proteins/" + encodeURIComponent(id)), api("/proteins/" + encodeURIComponent(id) + "/analysis" + qs({ hydropathy_window: win }))]);
    const A = { analysis: bundle.analysis, protein: record, label: record.id, fileBase: record.id, scanCache: {} };
    A.scanFn = async (exclude) => {
      const key = exclude ? "x" : "a";
      if (A.scanCache[key]) return A.scanCache[key];
      const result = await api("/proteins/" + encodeURIComponent(record.id) + "/prosite" + qs({ exclude_frequent: exclude ? "true" : "" }));
      if (result.status === "completed") A.scanCache[key] = result;
      return result;
    };
    A.setWindow = async (w) => { const res = await api("/proteins/" + encodeURIComponent(record.id) + "/analysis" + qs({ hydropathy_window: w })); A.analysis = res.analysis; A.window = w; };
    A.window = win;
    Session.add(record); Session.setAnalysis(record.id, A.analysis, record);
    return A;
  }
  Atlas.buildProteinA = buildProteinA;

  async function viewProtein(ctx) {
    const id = decodeURIComponent(ctx.params.id);
    const win = intParam(ctx.query.win, 9);
    setPage([el("a", { class: "back-link", href: "#/", onclick: (e) => { e.preventDefault(); Router.back("/"); } }, "\u2190 Back"), states.loadingBlock("Loading " + id + " from UniProt and running the analysis\u2026")]);
    let A;
    try { A = await buildProteinA(id, win); }
    catch (error) {
      if (!ctx.isCurrent()) return;
      setPage([el("a", { class: "back-link", href: "#/", onclick: (e) => { e.preventDefault(); Router.back("/"); } }, "\u2190 Back"), states.errorBlock(error, () => Router.render())]);
      return;
    }
    if (!ctx.isCurrent()) return;
    const p = A.protein;
    const path = "/protein/" + encodeURIComponent(p.id);
    const state = { id: p.id, A, openSrc: null };
    const withWin = (extra) => Object.assign({}, A.window && A.window !== 9 ? { win: A.window } : {}, extra);
    const tabsNode = Atlas.analysisTabs(A, { tab: ctx.query.tab || "overview", onTab: (tab) => Router.go(path, withWin({ tab })) });
    const origSet = A.setWindow;
    A.setWindow = async (w) => { await origSet(w); Router.go(path, withWin({ tab: "hydropathy" }), { replace: true }); };
    const teaser = el("button", { class: "structure-teaser", type: "button", onclick: () => Router.go(path, withWin({ tab: A.activeTab, view: "3d" })), "aria-label": "Open interactive 3D structure" }, [
      el("div", { class: "glyph" }, "\ud83e\uddec"), el("div", { class: "cta" }, "Open 3D structure"), el("small", {}, "PDB \u00b7 AlphaFold \u00b7 SWISS-MODEL")]);
    const root = el("div", { "data-page": "protein" }, [
      el("a", { class: "back-link", href: "#/", onclick: (e) => { e.preventDefault(); Router.back("/"); } }, "\u2190 Back"),
      el("div", { class: "protein-hero" }, [
        el("div", {}, [el("div", { class: "protein-hero-id" }, p.id + (p.entry_name ? " \u00b7 " + p.entry_name : "")), el("h1", { class: "protein-hero-name" }, fmt.text(p.name)),
          el("div", { class: "protein-hero-meta" }, [el("span", {}, ["Gene ", el("b", {}, fmt.text(p.gene))]), el("span", {}, [el("b", {}, fmt.text(p.organism))]), el("span", {}, [el("b", {}, fmt.int(p.length)), " aa"]),
            p.reviewed === true ? el("span", { class: "badge badge-ok" }, "Reviewed (Swiss-Prot)") : p.reviewed === false ? el("span", { class: "badge" }, "Unreviewed (TrEMBL)") : null])]),
        el("div", { class: "protein-hero-actions" }, [teaser, el("div", { class: "flex flex-wrap", style: "max-width:260px" }, [
          el("a", { class: "btn", href: Router.href("/compare", { a: p.id }) }, "Compare\u2026"), el("a", { class: "btn", href: Router.href("/patterns", { protein: p.id }) }, "Pattern search"),
          el("a", { class: "btn btn-ghost", href: "https://www.uniprot.org/uniprotkb/" + encodeURIComponent(p.id), target: "_blank", rel: "noopener" }, "UniProt \u2197")])]),
      ]),
      tabsNode,
    ]);
    state.root = root;
    state.apply = (query) => {
      const tab = query.tab || "overview";
      if (tab !== A.activeTab) tabsNode.activate(tab);
      if (query.view === "3d") {
        if (!Atlas.structure.isOpen() || state.openSrc !== (query.src || "")) {
          state.openSrc = query.src || "";
          Atlas.structure.openProtein(p, query.src, (value) => Router.go(path, withWin({ tab: A.activeTab, view: "3d", src: value }), { replace: true }));
        }
      } else if (Atlas.structure.isOpen()) { Atlas.structure.close(); state.openSrc = null; }
    };
    Atlas.pages.protein = state;
    Atlas.onCleanup(() => { Atlas.pages.protein = null; Atlas.structure.close(); });
    setPage(root);
    state.apply(ctx.query);
  }
  const softProtein = (ctx) => {
    const s = Atlas.pages.protein;
    if (!s || s.id !== decodeURIComponent(ctx.params.id) || !document.body.contains(s.root)) return false;
    s.apply(ctx.query); return true;
  };

  // ================================================================= Sequence Lab
  const LAB_EXAMPLES = [["Insulin B-chain", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"], ["MOTS-c peptide", "MRWQEMGYIFYPRKLR"], ["Poly-glycine", "GGGGG"],
    ["Ubiquitin", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"]];
  Atlas.lab = { sequence: "", analysis: null };
  async function viewLab(ctx) {
    const area = el("textarea", { placeholder: "Paste a protein sequence (plain letters or FASTA) \u2026", "aria-label": "Protein sequence", spellcheck: "false" }); area.value = Atlas.lab.sequence || "";
    const out = el("div", { class: "mt" });
    const run = async () => {
      clear(out);
      const sequence = area.value;
      if (!sequence.trim()) { out.appendChild(states.stateBlock("\u270e", "Enter a sequence", "Paste single-letter amino-acid codes (or a FASTA record).")); return; }
      out.appendChild(states.loadingBlock("Analysing\u2026"));
      let analysis;
      try { analysis = await api("/analyze", { body: { sequence, hydropathy_window: 9 } }); }
      catch (error) { clear(out); out.appendChild(states.errorBlock(error, run)); return; }
      Atlas.lab = { sequence, analysis };
      clear(out);
      const A = { analysis, protein: null, label: "Pasted sequence", fileBase: "sequence", scanCache: {} };
      A.scanFn = (exclude) => api("/prosite/scan", { body: { sequence, exclude_frequent: !!exclude } });
      A.setWindow = async (w) => { A.analysis = await api("/analyze", { body: { sequence, hydropathy_window: w } }); };
      out.appendChild(el("p", { class: "panel-note" }, [el("span", { class: "badge badge-ok" }, "Analysed"), " Same engine as the protein pages: " + fmt.int(analysis.sequence.length) + " residues."]));
      out.appendChild(Atlas.analysisTabs(A, { tab: ctx.query.tab || "overview" }));
      Atlas.scrub(out);
    };
    setPage([...head("SEQUENCE LAB", "Analyse any sequence", "No accession needed. Paste a sequence and get the complete analysis \u2014 identical to what a UniProt protein page shows."),
      card("Sequence", null, [area, el("div", { class: "flex flex-wrap", style: "margin-top:12px" }, [
        el("button", { class: "btn btn-primary", type: "button", onclick: run }, "Analyse"),
        ...LAB_EXAMPLES.map(([l, s]) => el("button", { class: "btn btn-ghost", type: "button", onclick: () => { area.value = s; run(); } }, l))])]), out]);
    if (Atlas.lab.sequence && Atlas.lab.analysis) run();
  }

  // ================================================================= Pattern Lab
  async function viewPatterns(ctx) {
    const holder = el("div", { class: "mt" });
    const acc = el("input", { type: "text", placeholder: "UniProt accession, e.g. P04637", value: ctx.query.protein || "", "aria-label": "UniProt accession" });
    const area = el("textarea", { placeholder: "\u2026or paste a sequence here", "aria-label": "Sequence", spellcheck: "false", style: "min-height:80px" });
    async function load() {
      clear(holder);
      let sequence = area.value, label = "pasted sequence";
      if (!sequence.trim()) {
        const id = acc.value.trim();
        if (!id) { holder.appendChild(states.stateBlock("\u270e", "Choose a sequence", "Enter a UniProt accession or paste a sequence.")); return; }
        holder.appendChild(states.loadingBlock("Loading " + id + "\u2026"));
        try { const rec = await api("/proteins/" + encodeURIComponent(id)); sequence = rec.sequence; label = rec.id + " \u2014 " + fmt.text(rec.name); Session.add(rec); }
        catch (error) { clear(holder); holder.appendChild(states.errorBlock(error, load)); return; }
        clear(holder);
      }
      holder.appendChild(el("p", { class: "panel-note" }, "Searching " + label + " (" + fmt.int(sequence.replace(/\s+/g, "").length) + " residues)."));
      holder.appendChild(card("Pattern", "Literal fragment or PROSITE syntax: x = any residue, [ST] = S or T, {P} = anything but P, x(2,4) = 2 to 4 residues, < and > anchor the termini.", Atlas.patternSearchPanel(() => sequence)));
    }
    setPage([...head("PATTERN LAB", "Search for motifs", "Find where a residue fragment or PROSITE-style pattern occurs in a protein."),
      card("Sequence source", null, [el("div", { class: "flex" }, [acc, el("button", { class: "btn btn-primary", type: "button", onclick: () => { area.value = ""; load(); } }, "Load protein")]), el("div", { class: "mt" }, area),
        el("button", { class: "btn", type: "button", onclick: load }, "Use pasted sequence")]), holder]);
    if (ctx.query.protein) load();
  }

  // ================================================================= Compare
  async function viewCompare(ctx) {
    const a = el("input", { type: "text", value: ctx.query.a || "", placeholder: "e.g. P04637", "aria-label": "First accession" });
    const b = el("input", { type: "text", value: ctx.query.b || "", placeholder: "e.g. P01308", "aria-label": "Second accession" });
    const k = el("input", { type: "text", value: ctx.query.k || "5", style: "width:90px", "aria-label": "k-mer size", inputmode: "numeric" });
    const sa = el("textarea", { placeholder: "First sequence (optional \u2014 overrides the accession)", spellcheck: "false", style: "min-height:70px" });
    const sb = el("textarea", { placeholder: "Second sequence (optional \u2014 overrides the accession)", spellcheck: "false", style: "min-height:70px" });
    const out = el("div", { class: "mt" });
    async function run(fromRoute) {
      const kk = parseInt(k.value, 10);
      if (!fromRoute && !sa.value.trim() && !sb.value.trim()) { Router.go("/compare", { a: a.value.trim(), b: b.value.trim(), k: k.value.trim() }); return; }
      clear(out); out.appendChild(states.loadingBlock("Comparing\u2026"));
      let res;
      try { res = await api("/compare", { body: { first: a.value.trim(), second: b.value.trim(), first_sequence: sa.value || null, second_sequence: sb.value || null, kmer_size: Number.isFinite(kk) ? kk : 0 } }); }
      catch (error) { clear(out); out.appendChild(states.errorBlock(error, () => run(fromRoute))); return; }
      clear(out); drawCompare(out, res); Atlas.scrub(out);
    }
    setPage([...head("COMPARE", "Compare two proteins", "Physicochemistry, composition and shared k-mers side by side \u2014 computed by the same engine as every other page."),
      card("Choose two proteins", null, [el("div", { class: "grid grid-2" }, [el("div", {}, [el("label", {}, "First protein (UniProt accession)"), a, el("div", { class: "mt", style: "margin-top:10px" }, sa)]),
        el("div", {}, [el("label", {}, "Second protein (UniProt accession)"), b, el("div", { style: "margin-top:10px" }, sb)])]),
        el("div", { class: "flex flex-wrap", style: "margin-top:14px" }, [el("label", { style: "margin:0" }, "k-mer size (1\u201320, at most the shorter length)"), k,
          el("button", { class: "btn btn-primary", type: "button", onclick: () => run(false) }, "Compare")])]), out]);
    if (ctx.query.a && ctx.query.b) run(true);
  }
  function drawCompare(out, r) {
    const s = r.similarity;
    out.appendChild(el("div", { class: "grid grid-2" }, [r.first, r.second].map((p, i) => card(i ? "Second" : "First", null, kv([["Protein", fmt.text(p.name)], ["Accession", fmt.text(p.id)], ["Length", fmt.int(p.length) + " aa"]])))));
    out.appendChild(el("div", { class: "mt grid grid-4" }, [
      tile("Word similarity (k = " + s.k + ")", fmt.pct(s.word_similarity_percent, 2), "", "Jaccard index of distinct k-mers"), tile("Shared " + s.k + "-mers", fmt.int(s.shared_kmers), "of " + fmt.int(s.union_kmers)),
      tile("Distinct " + s.k + "-mers", fmt.int(s.distinct_kmers[0]) + " / " + fmt.int(s.distinct_kmers[1])), tile("Composition distance", fmt.n(s.composition_distance, 3)),
    ]));
    out.appendChild(el("p", { class: "panel-note", style: "margin-top:12px" }, [el("span", { class: "badge " + (s.shared_kmers ? "badge-ok" : "") }, s.shared_kmers ? "Shared words found" : "Computed result: no shared words"), " " + fmt.text(s.interpretation)]));
    out.appendChild(el("div", { class: "mt" }, card("Physicochemical properties", null, table([
      { label: "Property", render: (p) => p.label + (p.unit ? " (" + p.unit + ")" : "") }, { label: "First", align: "num", render: (p) => fmt.n(p.first, p.key === "length" ? 0 : 3) },
      { label: "Second", align: "num", render: (p) => fmt.n(p.second, p.key === "length" ? 0 : 3) }, { label: "Difference", align: "num", render: (p) => fmt.signed(fmt.isNum(p.first) && fmt.isNum(p.second) ? p.first - p.second : null, p.key === "length" ? 0 : 3) },
    ], r.properties, { id: "compare-table" }))));
    const rel = r.properties.filter((p) => fmt.isNum(p.first) && fmt.isNum(p.second) && Math.max(Math.abs(p.first), Math.abs(p.second)) > 0);
    const relHolder = el("div", {});
    charts.bar(relHolder, { labels: rel.map((p) => p.label), rotate: true, tall: true, zoom: false, beginAtZero: true, yTitle: "% of the larger absolute value", ariaLabel: "Relative property comparison",
      datasets: [{ label: "First", data: rel.map((p) => Math.abs(p.first) / Math.max(Math.abs(p.first), Math.abs(p.second)) * 100), backgroundColor: COLORS.teal, borderWidth: 0 },
        { label: "Second", data: rel.map((p) => Math.abs(p.second) / Math.max(Math.abs(p.first), Math.abs(p.second)) * 100), backgroundColor: COLORS.amber, borderWidth: 0 }],
      tooltip: { label: (item) => { const p = rel[item.dataIndex]; return p ? item.dataset.label + ": " + fmt.n(item.datasetIndex ? p.second : p.first, 3) + " " + p.unit : ""; } } });
    out.appendChild(el("div", { class: "mt" }, card("Relative comparison", "Each property scaled so the larger absolute value is 100%.", relHolder)));
    const compHolder = el("div", {});
    charts.bar(compHolder, { labels: r.composition.map((c) => c.code), tall: true, yTitle: "% of sequence", xTitle: "Residue", ariaLabel: "Composition comparison",
      datasets: [{ label: "First", data: r.composition.map((c) => c.first_percent), backgroundColor: COLORS.teal, borderWidth: 0 }, { label: "Second", data: r.composition.map((c) => c.second_percent), backgroundColor: COLORS.amber, borderWidth: 0 }] });
    out.appendChild(el("div", { class: "mt" }, card("Amino-acid composition", "Percentage of each residue in the two sequences.", compHolder)));
    const profHolder = el("div", {});
    charts.line(profHolder, { xTitle: "k-mer size (k)", yTitle: "Word similarity (%)", tall: false, ariaLabel: "Similarity versus k",
      datasets: [{ label: "Jaccard similarity", data: s.profile.map((p) => ({ x: p.k, y: p.similarity_percent })), borderColor: COLORS.teal, backgroundColor: COLORS.teal, pointRadius: 4, borderWidth: 2 }],
      tooltip: { label: (item) => { const p = s.profile[item.dataIndex]; return p ? "k = " + p.k + ": " + fmt.pct(p.similarity_percent, 2) + " (" + fmt.int(p.shared) + " shared)" : ""; } } });
    out.appendChild(el("div", { class: "mt" }, card("Similarity versus word length", "Short words are shared by almost any two proteins; long shared words indicate real relatedness. " + fmt.text(s.method), profHolder)));
    const dots = r.dotplot;
    const dotHolder = el("div", {});
    if (dots.points.length) charts.scatter(dotHolder, { xTitle: "Position in second sequence", yTitle: "Position in first sequence", tall: true, legend: false, xRange: [1, Math.max(dots.length_second, 2)], ariaLabel: "Dot plot",
      datasets: [{ label: "Shared " + dots.k + "-mers", data: dots.points, backgroundColor: COLORS.teal, pointRadius: 2 }] });
    else dotHolder.appendChild(states.stateBlock("\u2205", "No shared " + dots.k + "-mers to plot", "The dot plot marks positions where both sequences contain the same " + dots.k + "-residue word."));
    out.appendChild(el("div", { class: "mt" }, card("Dot plot", "Each dot is a word shared by both sequences (word length " + dots.k + (dots.k !== dots.requested_k ? ", raised from " + dots.requested_k + " to keep the plot readable" : "") + ")." + (dots.truncated ? " Showing the first " + fmt.int(dots.points.length) + " of " + fmt.int(dots.pairs_total) + " matches." : ""), dotHolder)));
  }

  // ================================================================= Reference (amino acids)
  async function referenceRows() {
    if (!Atlas.cache.reference) Atlas.cache.reference = await api("/reference/amino-acids");
    return Atlas.cache.reference;
  }
  function residueDetail(r) {
    const list = (arr) => (arr && arr.length ? arr.join(", ") : "none");
    return el("div", { class: "grid grid-3" }, [
      card("Identity", null, kv([["Name", r.name], ["One / three letter", r.code + " \u00b7 " + r.three], ["Side-chain class", Atlas.classPill(r.side_chain_class)], ["Residue family", r.ring_type], ["Polarity", r.polarity], ["Stereochemistry", r.stereo], ["Aromatic", r.aromatic ? "Yes" : "No"]])),
      card("Chemistry", null, kv([["Free amino acid formula", r.formula], ["Residue formula (in a protein)", r.residue_formula], ["Free mass (average)", fmt.n(r.free_mass_da, 2) + " Da"], ["Residue mass (average)", fmt.n(r.residue_mass_da, 3) + " Da"],
        ["Free mass (monoisotopic)", fmt.n(r.monoisotopic_free_mass_da, 4) + " Da"], ["Residue mass (monoisotopic)", fmt.n(r.monoisotopic_residue_mass_da, 4) + " Da"], ["Hydropathy (Kyte\u2013Doolittle)", fmt.n(r.hydropathy, 1)]])),
      card("Ionisation and bonding", null, kv([["Charge at pH 7", r.charge_text], ["\u03b1-carboxyl pKa", fmt.n(r.pk1, 2)], ["\u03b1-amino pKa", fmt.n(r.pk2, 2)],
        ["Side-chain pKa (free amino acid)", r.pkr === null || r.pkr === undefined ? "none (not ionisable)" : fmt.n(r.pkr, 2) + " \u00b7 " + fmt.text(r.pkr_group)], ["Side-chain pKa (EMBOSS, used for protein pI)", r.side_chain_pka === null || r.side_chain_pka === undefined ? "none" : fmt.n(r.side_chain_pka, 1)],
        ["Isoelectric point (free amino acid)", fmt.n(r.isoelectric_point_free, 2)], ["Side-chain H-bond role", r.hbond_role], ["Donor atoms", list(r.hbond_donor_atoms)], ["Acceptor atoms", list(r.hbond_acceptor_atoms)]])),
    ]);
  }
  async function viewReference(ctx) {
    const holder = el("div", {});
    Atlas.onCleanup(() => Atlas.structure.close());
    setPage([...head("REFERENCE", "The 20 amino acids", "Every building block of a protein: chemistry, ionisation, hydrogen-bonding and an interactive 3D model with atoms and bonds."), holder]);
    holder.appendChild(states.loadingBlock("Loading\u2026"));
    let rows;
    try { rows = await referenceRows(); } catch (error) { clear(holder); holder.appendChild(states.errorBlock(error, () => Router.render())); return; }
    if (!ctx.isCurrent()) return;
    clear(holder);
    const sel = (ctx.query.aa || "").toUpperCase();
    const current = rows.find((r) => r.code === sel) || null;
    holder.appendChild(el("div", { class: "aa-grid" }, rows.map((r) => el("a", { class: "aa-tile", href: Router.href("/reference", { aa: r.code }), style: "color:" + (Atlas.CLASS_COLOR[r.side_chain_class] || COLORS.dim) + ";" + (current && current.code === r.code ? "border-color:currentColor;background:rgba(87,217,196,0.08)" : ""), title: r.name }, [
      el("div", { class: "code" }, r.code), el("div", { class: "abbr" }, r.three)]))));
    const detail = el("div", { class: "aa-detail" }); holder.appendChild(detail);
    if (current) {
      detail.appendChild(el("h2", { class: "panel-title", style: "font-size:26px;margin-top:18px" }, current.name + " (" + current.three + " \u00b7 " + current.code + ")"));
      detail.appendChild(el("p", { class: "panel-note" }, fmt.text(current.note)));
      detail.appendChild(residueDetail(current));
      const viewerCard = el("div", { class: "card mt" }); detail.appendChild(viewerCard);
      viewerCard.appendChild(el("h3", { class: "panel-title" }, "Interactive 3D model"));
      viewerCard.appendChild(el("p", { class: "panel-note" }, "Drag to rotate \u00b7 scroll to zoom \u00b7 right-drag to pan \u00b7 hover an atom for details. Bonds are drawn from the molecular graph; double bonds are shown as such."));
      const host = el("div", { class: "mini-structure", style: "height:400px" }); viewerCard.appendChild(host);
      const controls = el("div", { class: "flex flex-wrap", style: "margin-top:10px" }); viewerCard.appendChild(controls);
      const info = el("p", { class: "panel-note", style: "margin-top:8px" }, "Loading 3D model\u2026"); viewerCard.appendChild(info);
      let viewer = null, model = null;
      Atlas.onCleanup(() => { if (viewer) viewer.destroy(); });
      try {
        model = await api("/reference/amino-acids/" + encodeURIComponent(current.code));
        if (!ctx.isCurrent()) return;
        if (typeof global.$3Dmol === "undefined") throw new Atlas.AtlasError("database_unavailable", "The 3D viewer library (3Dmol.js) could not be downloaded. Check your internet connection or ad-blocker and refresh with Ctrl+F5.");
        viewer = new Atlas.MolViewer(host, { small: true });
        viewer.onHover = (atom) => { info.textContent = atom ? atom.elem + " atom #" + atom.serial + " \u00b7 bonded to " + (atom.bonds || []).length + " atoms" : baseInfo; };
        const atoms = viewer.load(model.structure.data, model.structure.format);
        const baseInfo = current.name + ": " + atoms + " atoms (" + current.heavy_atoms + " heavy) and " + current.bonds_total + " bonds. " + fmt.text(model.structure.note);
        info.textContent = baseInfo;
        const styleSeg = (items, key, cur) => { const g = el("div", { class: "seg-group" }); items.forEach(([v, l]) => { const bt = el("button", { class: "seg" + (v === cur ? " active" : ""), type: "button", onclick: () => { Array.from(g.children).forEach((c) => c.classList.toggle("active", c === bt)); if (key === "style") viewer.setStyle(v); else viewer.setColor(v); } }, l); g.appendChild(bt); }); return g; };
        controls.appendChild(styleSeg([["ballstick", "Ball & stick"], ["stick", "Stick"], ["sphere", "Spacefill"], ["line", "Line"]], "style", "ballstick"));
        controls.appendChild(styleSeg([["element", "Element"], ["part", "Backbone / side chain"]], "color", "element"));
        controls.appendChild(el("button", { class: "btn btn-ghost", type: "button", onclick: (e) => { const on = viewer.toggle("showLabels"); e.currentTarget.classList.toggle("active", on); } }, "Atom labels"));
        controls.appendChild(el("button", { class: "btn btn-ghost", type: "button", onclick: () => viewer.toggle("spinning") }, "\u27f2 Spin"));
        controls.appendChild(el("button", { class: "btn btn-ghost", type: "button", onclick: () => viewer.resetCamera() }, "Reset"));
        controls.appendChild(el("button", { class: "btn", type: "button", onclick: () => Router.go("/reference", { aa: current.code, view: "3d" }) }, "Expand full screen"));
        controls.appendChild(el("button", { class: "btn btn-ghost", type: "button", onclick: () => Atlas.copyText(current.smiles) }, "Copy SMILES"));
        controls.appendChild(el("a", { class: "btn btn-ghost", href: current.pubchem_url, target: "_blank", rel: "noopener" }, "PubChem \u2197"));
      } catch (error) { clear(host); host.appendChild(states.errorBlock(error, () => Router.render())); info.textContent = ""; }
      if (ctx.query.view === "3d" && model) Atlas.structure.openMolecule({ title: current.name + " (" + current.code + ")", subtitle: current.formula + " \u00b7 " + fmt.text(model.structure.source), format: model.structure.format, data: model.structure.data });
      detail.appendChild(el("p", { class: "panel-note mt" }, "SMILES: " + current.smiles));
    } else detail.appendChild(el("p", { class: "panel-note mt" }, "Click an amino acid above to see its full chemistry and an interactive 3D model."));
    holder.appendChild(el("div", { class: "mt" }, card("All 20 side by side", null, table([
      { label: "Code", align: "mono", render: (r) => el("a", { href: Router.href("/reference", { aa: r.code }) }, r.code) }, { label: "Name", render: (r) => r.name }, { label: "Class", render: (r) => Atlas.classPill(r.side_chain_class) },
      { label: "Formula", align: "mono", render: (r) => r.formula }, { label: "Residue mass", align: "num", render: (r) => fmt.n(r.residue_mass_da, 2) }, { label: "Charge pH 7", render: (r) => r.charge_text.split(" (")[0] },
      { label: "pKa side chain", align: "num", render: (r) => (r.pkr === null || r.pkr === undefined ? "\u2014" : fmt.n(r.pkr, 2)) }, { label: "pI (free)", align: "num", render: (r) => fmt.n(r.isoelectric_point_free, 2) },
      { label: "Hydropathy", align: "num", render: (r) => fmt.n(r.hydropathy, 1) }, { label: "H-bond role", render: (r) => r.hbond_role }], rows, { id: "reference-table" }))));
  }
  // ================================================================= Dataset statistics
  async function viewDataset(ctx) {
    const scope = ctx.query.scope || "human_reviewed";
    const holder = el("div", {}); const selectBox = el("div", {});
    setPage([...head("DATASET", "Dataset statistics", "Complete statistics for an entire UniProt dataset \u2014 every sequence is streamed once, analysed with the same engine, and the result is cached."), selectBox, holder]);
    let alive = true, timer = null;
    Atlas.onCleanup(() => { alive = false; clearTimeout(timer); });
    let scopes = {};
    try { scopes = (await api("/dataset")).scopes.reduce((m, s) => { m[s.key] = s.label; return m; }, {}); } catch (e) { scopes = { [scope]: scope }; }
    selectBox.appendChild(el("div", { class: "flex flex-wrap", style: "margin-bottom:18px" }, [el("label", { style: "margin:0" }, "Dataset"),
      el("select", { style: "width:auto", "aria-label": "Dataset scope", onchange: (e) => Router.go("/dataset", { scope: e.target.value }) }, Object.keys(scopes).map((k) => el("option", { value: k, selected: k === scope }, scopes[k])))]));
    const universe = el("div", {}); holder.appendChild(universe);
    api("/dataset/universe").then((u) => {
      if (!alive) return;
      const t = [["UniProtKB (all entries)", u.uniprotkb_total], ["Reviewed (Swiss-Prot)", u.reviewed_total], ["Human (all)", u.human_total], ["Human reviewed", u.human_reviewed_total]].filter((x) => fmt.isNum(x[1]));
      if (t.length) universe.appendChild(el("div", { class: "card", style: "margin-bottom:18px" }, [el("h3", { class: "panel-title" }, "UniProt right now"), el("div", { class: "grid grid-4" }, t.map(([l, v]) => tile(l, fmt.int(v), "entries")))]));
    }).catch(() => { /* optional context only */ });
    const body = el("div", {}); holder.appendChild(body);
    async function load(refresh) {
      let st;
      try { st = await api("/dataset/stats" + qs({ scope, refresh: refresh ? "true" : "" })); }
      catch (error) { if (alive) { clear(body); body.appendChild(states.errorBlock(error, () => load(false))); } return; }
      if (!alive) return;
      clear(body);
      if (st.state === "building") {
        const pr = st.progress || {};
        body.appendChild(el("div", { class: "card" }, [el("h3", { class: "panel-title" }, "Building the complete statistics\u2026"),
          el("p", { class: "panel-note" }, "One-off job: the whole dataset (" + fmt.text(st.label) + ") is streamed from UniProt and analysed. " + fmt.text(pr.phase)),
          el("div", { class: "progress" }, el("div", { class: "progress-fill", style: "width:" + (fmt.isNum(pr.percent) ? Math.min(100, pr.percent).toFixed(1) : 3) + "%" })),
          el("p", { class: "muted" }, fmt.int(pr.processed) + (fmt.isNum(pr.total) ? " of " + fmt.int(pr.total) : "") + " sequences processed" + (fmt.isNum(pr.percent) ? " (" + fmt.n(pr.percent, 1) + "%)" : ""))]));
        timer = setTimeout(() => load(false), 2500); return;
      }
      if (st.state === "failed") { body.appendChild(states.stateBlock("\u26a0", "Dataset statistics unavailable", st.error, el("button", { class: "btn", type: "button", onclick: () => load(true) }, "Rebuild"))); return; }
      drawStats(body, st, () => load(true));
    }
    await load(false);
  }
  function statsTable(d, unit, digits) {
    const rows = [["Mean", d.mean], ["Median", d.median], ["Minimum", d.min], ["Maximum", d.max], ["Standard deviation", d.std], ["5th percentile", d.p05], ["25th percentile", d.p25], ["75th percentile", d.p75], ["95th percentile", d.p95]];
    return table([{ label: "Statistic" }, { label: "Value" + (unit ? " (" + unit + ")" : ""), align: "num" }].map((c, i) => Object.assign(c, { render: (r) => (i ? fmt.n(r[1], digits) : r[0]) })), rows);
  }
  function datasetCsv(st) {
    const s = st.stats, rows = [["dataset", fmt.text(st.label), "", ""], ["query", fmt.text(st.query), "", ""], ["proteins", s.total_proteins, "", ""], ["residues", s.total_residues, "", ""]];
    [["length", "aa"], ["molecular_weight", "Da"], ["isoelectric_point", "pH"], ["gravy", ""], ["aromaticity", "%"]].forEach(([key, unit]) => {
      ["mean", "median", "min", "max", "std", "p05", "p25", "p75", "p95"].forEach((stat) => rows.push([key, stat, fmt.plain(s[key][stat], 4), unit]));
    });
    rows.push(["", "", "", ""], ["amino_acid", "count", "percent", "name"]);
    s.amino_acids.forEach((r) => rows.push([r.code, r.count, fmt.plain(r.percent, 4), r.name]));
    return Atlas.csv(["section", "item", "value", "unit"], rows);
  }
  function drawStats(body, st, rebuild) {
    const s = st.stats;
    body.appendChild(el("div", { class: "card" }, [el("div", { class: "flex flex-wrap between" }, [el("div", {}, [el("h3", { class: "panel-title" }, fmt.text(st.label)), el("p", { class: "panel-note" }, "Query: " + fmt.text(st.query) + " \u00b7 built " + fmt.text(st.built_at) + " \u00b7 complete for this dataset: all " + fmt.int(s.total_proteins) + " sequences were analysed" + (fmt.isNum(s.reported_total) ? " (UniProt reports " + fmt.int(s.reported_total) + ")" : "") + ".")]),
      el("div", { class: "flex" }, [el("button", { class: "btn", type: "button", onclick: () => Atlas.download("dataset-" + fmt.text(st.scope) + "-summary.csv", datasetCsv(st), "text/csv") }, "Download summary CSV"),
        el("button", { class: "btn btn-ghost", type: "button", onclick: rebuild }, "Rebuild statistics")])]),
      el("div", { class: "grid grid-4" }, [tile("Proteins", fmt.int(s.total_proteins)), tile("Total residues", fmt.int(s.total_residues)), tile("Mean length", fmt.n(s.length.mean, 1), "aa"), tile("Median length", fmt.n(s.length.median, 0), "aa"),
        tile("Shortest", fmt.int(s.length.min), "aa"), tile("Longest", fmt.int(s.length.max), "aa"), tile("Mean molecular weight", fmt.n(s.molecular_weight.mean / 1000, 2), "kDa"), tile("Median pI", fmt.n(s.isoelectric_point.median, 2))])]));
    const dist = (title, key, unit, digits, note, color) => {
      const c = s[key], h = el("div", {});
      charts.histogram(h, c.histogram, { xTitle: title + (unit ? " (" + unit + ")" : ""), color, ariaLabel: title + " distribution" });
      return el("div", { class: "mt grid grid-2" }, [card(title + " distribution", note, h), card(title + " statistics", null, statsTable(c, unit, digits))]);
    };
    body.appendChild(dist("Length", "length", "aa", 1, "Number of proteins in each length bin (last bin = longer).", COLORS.teal));
    body.appendChild(dist("Molecular weight", "molecular_weight", "Da", 1, fmt.text(s.molecular_weight.note), COLORS.amber));
    body.appendChild(dist("Isoelectric point", "isoelectric_point", "pH", 2, "Bimodal: acidic and basic proteins.", COLORS.violet));
    body.appendChild(dist("GRAVY", "gravy", "", 3, "Negative = hydrophilic, positive = hydrophobic.", COLORS.blue));
    body.appendChild(dist("Aromaticity", "aromaticity", "%", 2, "Share of F, W and Y residues.", COLORS.coral));
    const aaHolder = el("div", {}), famHolder = el("div", {});
    charts.bar(aaHolder, { labels: s.amino_acids.map((r) => r.code), legend: false, tall: true, yTitle: "% of all residues", xTitle: "Residue", ariaLabel: "Amino-acid frequencies",
      datasets: [{ label: "Frequency", data: s.amino_acids.map((r) => r.percent), backgroundColor: s.amino_acids.map((r) => Atlas.CLASS_COLOR[r.family] || COLORS.dim), borderWidth: 0 }],
      tooltip: { title: (it) => { const r = s.amino_acids[it[0].dataIndex]; return r ? r.name : ""; }, label: (it) => { const r = s.amino_acids[it.dataIndex]; return r ? fmt.int(r.count) + " residues \u00b7 " + fmt.pct(r.percent, 3) : ""; } } });
    charts.doughnut(famHolder, { labels: s.families.map((f) => f.family), datasets: [{ data: s.families.map((f) => f.count), backgroundColor: s.families.map((f) => Atlas.CLASS_COLOR[f.family]), borderWidth: 0 }], tall: true, ariaLabel: "Family composition" });
    body.appendChild(el("div", { class: "mt grid grid-2" }, [card("Amino-acid frequencies", "Across every residue of every sequence in the dataset.", aaHolder), card("Chemical families", null, famHolder)]));
    body.appendChild(el("div", { class: "mt" }, card("Amino-acid table", null, table([{ label: "Residue", render: (r) => r.name }, { label: "Code", align: "mono", render: (r) => r.code }, { label: "Count", align: "num", render: (r) => fmt.int(r.count) }, { label: "Percentage", align: "num", render: (r) => fmt.pct(r.percent, 3) }], s.amino_acids, { id: "dataset-aa-table" }))));
    const v = s.validation;
    body.appendChild(el("div", { class: "mt" }, card("Sequence validation", "Every sequence was checked before it was counted.", el("div", { class: "grid grid-4" }, [tile("Sequences checked", fmt.int(v.proteins_checked)), tile("With non-standard codes", fmt.int(v.proteins_with_nonstandard_codes)),
      tile("With unexpected characters", fmt.int(v.proteins_with_unexpected_characters)), tile("Standard residues", fmt.pct(v.standard_residue_share_percent, 3))]))));
  }

  // ================================================================= Export
  async function viewExport(ctx) {
    const acc = el("input", { type: "text", placeholder: "UniProt accession", value: ctx.query.protein || "", "aria-label": "UniProt accession" });
    const out = el("div", { class: "mt" });
    async function load() {
      clear(out); const id = acc.value.trim();
      if (!id) { out.appendChild(states.stateBlock("\u270e", "Choose a protein", "Enter an accession or pick one you already opened.")); return; }
      out.appendChild(states.loadingBlock("Preparing " + id + "\u2026"));
      try { const A = await buildProteinA(id, 9); clear(out); out.appendChild(Atlas.exportCard(A)); Atlas.scrub(out); }
      catch (error) { clear(out); out.appendChild(states.errorBlock(error, load)); }
    }
    const opened = Session.list();
    setPage([...head("EXPORT", "Export analysis", "Download FASTA, CSV tables, a JSON report or a text report for any protein."),
      card("Protein", null, [el("div", { class: "flex" }, [acc, el("button", { class: "btn btn-primary", type: "button", onclick: load }, "Prepare files")]),
        opened.length ? el("div", { class: "flex flex-wrap", style: "margin-top:12px" }, [el("span", { class: "muted" }, "Opened this session:"), ...opened.map((p) => el("button", { class: "btn btn-ghost", type: "button", onclick: () => { acc.value = p.id; load(); } }, p.id))]) : null]), out]);
    if (ctx.query.protein) load();
  }

  function viewNotFound() { setPage([...head("404", "Page not found", "That address does not exist in Protein Atlas."), el("a", { class: "btn btn-primary", href: "#/" }, "Back to Discover")]); }

  Router.add("discover", /^\/$/, viewDiscover);
  Router.add("protein", /^\/protein\/(?<id>[^/]+)\/?$/, viewProtein, softProtein);
  Router.add("lab", /^\/lab$/, viewLab);
  Router.add("patterns", /^\/patterns$/, viewPatterns);
  Router.add("compare", /^\/compare$/, viewCompare);
  Router.add("reference", /^\/reference$/, async (ctx) => { Atlas.pages.refAa = (ctx.query.aa || "").toUpperCase(); return viewReference(ctx); });
  Router.add("dataset", /^\/dataset$/, viewDataset);
  Router.add("export", /^\/export$/, viewExport);
  Router.add("notfound", /^.*$/, viewNotFound);
})(typeof window !== "undefined" ? window : globalThis);
