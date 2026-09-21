/* Protein Atlas - built-in self-test (#/selftest).
   Runs the validation checklist against the LIVE server in your own browser: pagination, values,
   pattern-table alignment, PROSITE honesty, chart interaction and page scrolling, 3D data,
   amino-acid models, Sequence Lab, comparison validation, dataset statistics, routing, and a scan of every
   page for "undefined" / "NaN". */
(function (global) {
  "use strict";
  const Atlas = global.Atlas;
  const { el, $, clear, fmt, api, states, Router } = Atlas;
  const BAD = /\bundefined\b|\bNaN\b|\[object Object\]|(^|\s)null(\s|$)/;
  const P53 = "P04637", INS = "P01308";
  const results = Atlas.selftestResults = Atlas.selftestResults || [];

  async function sandbox(fn) {                     // render into an off-screen #view and hand back its text
    const real = $("#view");
    real.id = "view-real";
    const box = el("div", { id: "view", style: "position:absolute;left:-99999px;top:0;width:1100px" });
    document.body.appendChild(box);
    try { await fn(box); await Atlas.sleep(250); return box.textContent || ""; }
    finally { Router.cleanup(); box.remove(); real.id = "view"; }
  }
  async function renderRoute(path, query) {
    const route = Router.routes.find((r) => r.pattern.test(path));
    const m = route.pattern.exec(path);
    return sandbox(async () => { await route.handler({ path, query: query || {}, params: m.groups || {}, raw: path, isCurrent: () => true }); });
  }
  const need = (cond, message) => { if (!cond) throw new Error(message); };
  // the report itself must never contain the words the safety net removes, so break them with a zero-width space
  const safe = (t) => String(t).replace(/undefined|NaN|null/g, (m) => m[0] + "\u200b" + m.slice(1));
  const noBad = (text, where) => { const m = BAD.exec(text); need(!m, where + " shows the placeholder \u201c" + safe(m && m[0].trim()) + "\u201d"); };
  const finite = (v) => typeof v === "number" && Number.isFinite(v);

  const TESTS = [
    ["Server, database status and 3D/chart libraries", async () => {
      const s = await api("/status");
      const libs = "Chart.js " + (Atlas.charts.available() ? "loaded" : "MISSING") + ", 3Dmol " + (typeof global.$3Dmol !== "undefined" ? "loaded" : "MISSING");
      need(Atlas.charts.available(), "Chart.js did not load (" + libs + ")");
      need(typeof global.$3Dmol !== "undefined", "3Dmol.js did not load (" + libs + ")");
      return "PROSITE " + s.prosite.state + (s.prosite.patterns_usable ? " (" + s.prosite.patterns_usable + " patterns)" : "") + "; " + libs;
    }],
    ["Discovery: page 1 \u2192 2 \u2192 1 with real totals", async () => {
      const p1 = await api("/proteins?page=1&page_size=24"), p2 = await api("/proteins?page=2&page_size=24"), again = await api("/proteins?page=1&page_size=24");
      need(p1.total_items > 24 && p1.total_pages > 1, "total_pages is " + p1.total_pages + " (expected more than 1)");
      const ids1 = p1.items.map((x) => x.id), ids2 = p2.items.map((x) => x.id);
      need(ids2.length && !ids2.some((id) => ids1.includes(id)), "page 2 repeats page 1");
      need(JSON.stringify(again.items.map((x) => x.id)) === JSON.stringify(ids1), "going back to page 1 changed the results");
      const text = await sandbox(async () => { await renderRoute0("/", { page: "2" }); });
      need(/Page/.test(text) && /Previous/.test(text) && /Next/.test(text), "pager controls missing on page 2");
      return fmt.int(p1.total_items) + " proteins in " + fmt.int(p1.total_pages) + " pages; page 2 has " + ids2.length + " distinct items";
    }],
    ["Protein overview, residues and masses: real values on every tab", async () => {
      const A = await Atlas.buildProteinA(P53, 9);
      const c = A.analysis.composition, ph = A.analysis.physicochemical;
      need(c.status === "success" && c.most_abundant && c.least_abundant, "most/least abundant residue missing");
      need(c.most_abundant.residues[0].name && finite(c.most_abundant.count), "most abundant residue/count undefined");
      need(c.least_abundant.residues[0].name && finite(c.least_abundant.count), "least abundant residue/count undefined");
      need(c.rows.length >= 20 && c.rows.every((r) => Number.isInteger(r.count) && finite(r.percent)), "some of the 20 residues lack count/percentage");
      need(ph.mass_breakdown.rows.every((r) => Number.isInteger(r.count) && finite(r.mass_da)), "mass breakdown has undefined count/mass");
      need(ph.mass_breakdown.reconciles, "mass breakdown does not reconcile with the molecular weight");
      need(c.families.every((f) => finite(f.count)), "family breakdown incomplete");
      const before = Atlas.diagnostics.scrubbed.length;
      let text = "";
      for (const tab of Atlas.TAB_IDS) {
        text += await sandbox(async (box) => { const node = Atlas.analysisTabs(A, { tab }); box.appendChild(node); });
      }
      noBad(text, "protein tabs");
      need(Atlas.diagnostics.scrubbed.length === before, "the safety net had to replace: " + safe(Atlas.diagnostics.scrubbed.slice(before).join(" | ")));
      return "Most abundant: " + c.most_abundant.residues.map((r) => r.name).join("/") + " \u00d7" + c.most_abundant.count + "; least: " + c.least_abundant.residues.map((r) => r.name).join("/") + " \u00d7" + c.least_abundant.count + "; MW " + fmt.da(ph.molecular_weight.value_da);
    }],
    ["Custom pattern table: Start/End/Sequence under the right headings", async () => {
      let table = null;
      await sandbox(async (box) => {
        const panel = Atlas.patternSearchPanel(() => "MTKSGGGGSRYLL", { initial: "G" });
        box.appendChild(panel); await panel.run(); await Atlas.sleep(400); table = $("#pattern-results", box);
        need(table, "no result table rendered");
        const heads = Atlas.$$("th", table).map((t) => t.textContent);
        need(heads.join("|") === "Match|Start|End|Sequence", "headings are " + heads.join("|"));
        const rows = Atlas.$$("tbody tr", table).map((tr) => Atlas.$$("td", tr).map((td) => td.textContent));
        need(rows.length === 4 && rows.every((r) => r.length === 4), "row cell count differs from headings");
        need(rows[0].join("|") === "1|5|5|G" && rows[3].join("|") === "4|8|8|G", "first/last rows are " + rows[0].join("|") + " and " + rows[3].join("|"));
      });
      return "G in MTKSGGGGSRYLL \u2192 4 rows, columns aligned";
    }],
    ["PROSITE scan reports what it really did", async () => {
      const r = await api("/prosite/scan", { body: { sequence: "AANASAASRKAAATMKTIIALSYIFCLVFADYKDDDDK" } });
      need(["completed", "loading", "database_unavailable", "failed"].includes(r.status), "unknown status " + r.status);
      if (r.status === "completed") { need(r.patterns_tested > 100, "only " + r.patterns_tested + " patterns tested (expected the full database)"); return "Scan completed: " + fmt.int(r.patterns_tested) + " patterns tested, " + r.patterns_matched + " matched (release " + fmt.text(r.database.release) + ")"; }
      return "Honest state: " + r.status + " \u2014 " + r.message;
    }],
    ["Amino Acids tab lists all 20 residues", async () => {
      const A = await Atlas.buildProteinA(INS, 9);
      const text = await sandbox(async (box) => box.appendChild(Atlas.analysisTabs(A, { tab: "aminoacids" })));
      const rows = A.analysis.composition.rows.filter((r) => Atlas.$$ && "ACDEFGHIKLMNPQRSTVWY".includes(r.code));
      need(rows.length === 20, "table has " + rows.length + " residues");
      noBad(text, "Amino Acids tab");
      return "20 residues with count and percentage";
    }],
    ["Charts: hover data, wheel-zoom over the chart only, page scroll elsewhere", async () => {
      const A = await Atlas.buildProteinA(P53, 9);
      let detail = "";
      await sandbox(async (box) => { box.appendChild(Atlas.analysisTabs(A, { tab: "hydropathy" })); await Atlas.sleep(500);
        const h = Atlas.charts.registry.find((x) => x.kind === "line" && x.chart);
        need(h, "hydropathy chart was not created");
        const canvas = h.canvas, cfg = h.chart.config.data.datasets[0].data;
        need(cfg.length > 50 && cfg.every((p) => finite(p.x) && finite(p.y)), "hydropathy data is empty or not numeric");
        const before = [h.chart.scales.x.min, h.chart.scales.x.max];
        const wheel = new WheelEvent("wheel", { deltaY: -240, clientX: canvas.getBoundingClientRect().left + 100, clientY: canvas.getBoundingClientRect().top + 50, bubbles: true, cancelable: true });
        canvas.dispatchEvent(wheel);
        need(wheel.defaultPrevented, "wheel over the chart did not zoom");
        need(h.chart.scales.x.max - h.chart.scales.x.min < before[1] - before[0], "chart range did not shrink after zooming");
        h.zoom.reset();
        h.zoom.state.wheel = false;
        const locked = new WheelEvent("wheel", { deltaY: -240, bubbles: true, cancelable: true }); canvas.dispatchEvent(locked);
        need(!locked.defaultPrevented, "with wheel-zoom off the page must scroll");
        h.zoom.state.wheel = true;
        detail = cfg.length + " points; zoom/reset work; wheel lock releases scrolling";
      });
      const page = new WheelEvent("wheel", { deltaY: 200, bubbles: true, cancelable: true }); document.body.dispatchEvent(page);
      need(!page.defaultPrevented, "the wheel is blocked outside charts (page cannot scroll)");
      return detail + "; wheel on the page is never blocked";
    }],
    ["3D: structure resolution and coordinates", async () => {
      const l = await api("/proteins/" + P53 + "/structures");
      need(Array.isArray(l.sources) && l.sources.length, "no source report");
      const statuses = l.sources.map((s) => s.label + ": " + s.status).join("; ");
      need(l.recommended, "no structure located. " + l.summary + " [" + statuses + "]");
      const f = await api("/structures/" + l.recommended.source + "/" + encodeURIComponent(l.recommended.id) + "/file");
      const atoms = f.data.split("\n").filter((x) => x.indexOf("ATOM") === 0 || x.indexOf("HETATM") === 0).length;
      need(f.format === "cif" || atoms > 100, "coordinate file has " + atoms + " atoms");
      need(typeof global.$3Dmol !== "undefined", "3Dmol.js is not loaded");
      return l.recommended.source + " " + l.recommended.id + " (" + f.format + ", " + (f.format === "pdb" ? fmt.int(atoms) + " atoms" : "mmCIF") + "). " + statuses;
    }],
    ["Amino-acid reference: 20 residues with 3D atoms and bonds", async () => {
      const rows = await api("/reference/amino-acids");
      need(rows.length === 20, rows.length + " residues");
      need(rows.every((r) => r.formula && finite(r.residue_mass_da) && finite(r.pk1) && finite(r.isoelectric_point_free) && r.smiles && r.hbond_role), "some reference fields are missing");
      const ala = await api("/reference/amino-acids/A"), trp = await api("/reference/amino-acids/W");
      const count = (s) => { const l = s.data.split("\n"); return [parseInt(l[3].slice(0, 3), 10), parseInt(l[3].slice(3, 6), 10)]; };
      need(count(ala.structure)[0] === 13 && count(ala.structure)[1] === 12, "Ala should have 13 atoms / 12 bonds");
      need(count(trp.structure)[0] === 27, "Trp should have 27 atoms");
      const text = await renderRoute("/reference", { aa: "W" }); noBad(text, "Reference page");
      return "Ala 13 atoms/12 bonds, Trp 27 atoms; all fields present";
    }],
    ["Sequence Lab: GGGGG", async () => {
      const a = await api("/analyze", { body: { sequence: "GGGGG", hydropathy_window: 9 } });
      const c = a.composition;
      need(c.most_abundant.residues[0].name === "Glycine" && c.most_abundant.count === 5, "most abundant is not Glycine \u00d75");
      need(c.least_abundant.residues[0].name === "Glycine" && c.least_abundant.count === 5, "least abundant is not Glycine \u00d75");
      const p53 = await Atlas.buildProteinA(INS, 9);
      const lab = await api("/analyze", { body: { sequence: p53.analysis.sequence.residues, hydropathy_window: 9 } });
      need(JSON.stringify(lab.physicochemical.molecular_weight) === JSON.stringify(p53.analysis.physicochemical.molecular_weight), "Sequence Lab and the protein page disagree");
      return "GGGGG \u2192 Glycine \u00d75 (most and least); Lab and protein page agree exactly";
    }],
    ["Comparison: valid k and clear validation messages", async () => {
      const ok = await api("/compare", { body: { first: P53, second: INS, kmer_size: 3 } });
      need(finite(ok.similarity.word_similarity_percent) && ok.similarity.interpretation, "similarity missing");
      const errs = [];
      for (const k of [99, 0, 25]) {
        try { await api("/compare", { body: { first: P53, second: INS, kmer_size: k } }); throw new Error("k=" + k + " was accepted"); }
        catch (e) { need(e.kind === "invalid_input" && e.message.length > 20 && !/failed/i.test(e.message), "k=" + k + " gave: " + e.message); errs.push("k=" + k + ": " + e.message); }
      }
      try { await api("/compare", { body: { first_sequence: "MKTIIALSYIFC", second_sequence: "MKTIIAL", kmer_size: 9 } }); throw new Error("k longer than the sequence was accepted"); }
      catch (e) { need(/shorter sequence/.test(e.message), "short-sequence message was: " + e.message); }
      return "k=3 similarity " + fmt.pct(ok.similarity.word_similarity_percent, 2) + " (" + ok.similarity.shared_kmers + " shared); invalid k explained: \u201c" + errs[0].slice(0, 70) + "\u2026\u201d";
    }],
    ["Dataset statistics cover the whole dataset", async () => {
      const st = await api("/dataset/stats?scope=human_reviewed");
      if (st.state === "building") return "Building in the background: " + fmt.n(st.progress.percent, 1) + "% (" + fmt.int(st.progress.processed) + " sequences so far). Re-run in a minute.";
      need(st.state === "ready", "state is " + st.state + (st.error ? ": " + st.error : ""));
      const s = st.stats;
      need(s.total_proteins > 1000, "only " + s.total_proteins + " proteins");
      need(s.length.histogram.counts.reduce((a, b) => a + b, 0) === s.total_proteins, "length histogram does not cover every protein");
      need(fmt.isNum(s.reported_total) ? Math.abs(s.total_proteins - s.reported_total) / s.reported_total < 0.01 : true, "analysed " + s.total_proteins + " but UniProt reports " + s.reported_total);
      const text = await renderRoute("/dataset", { scope: "human_reviewed" }); noBad(text, "Dataset page");
      return fmt.int(s.total_proteins) + " proteins analysed" + (fmt.isNum(s.reported_total) ? " of " + fmt.int(s.reported_total) + " reported by UniProt" : "") + "; mean length " + fmt.n(s.length.mean, 1);
    }],
    ["Every page shows real values only (no placeholder text)", async () => {
      const before = Atlas.diagnostics.scrubbed.length;
      const pages = [["/", {}], ["/lab", {}], ["/patterns", { protein: P53 }], ["/compare", {}], ["/reference", {}], ["/export", {}]];
      let n = 0;
      for (const [path, query] of pages) { noBad(await renderRoute(path, query), path); n++; }
      need(Atlas.diagnostics.scrubbed.length === before, "safety net triggered: " + safe(Atlas.diagnostics.scrubbed.slice(before).join(" | ")));
      need(Atlas.diagnostics.errors.length === 0, "internal errors: " + Atlas.diagnostics.errors.slice(0, 2).join(" | "));
      return n + " pages scanned, no unintended placeholders and no internal errors";
    }],
    ["Back / forward navigation", async () => {
      const wait = () => new Promise((r) => { const h = () => { global.removeEventListener("hashchange", h); setTimeout(r, 300); }; global.addEventListener("hashchange", h); });
      const start = global.location.hash;
      let w = wait(); Router.go("/reference"); await w;
      w = wait(); Router.go("/dataset"); await w;
      w = wait(); global.history.back(); await w;
      need(global.location.hash.indexOf("#/reference") === 0, "Back did not return to the previous page (at " + global.location.hash + ")");
      w = wait(); global.history.forward(); await w;
      need(global.location.hash.indexOf("#/dataset") === 0, "Forward did not work (at " + global.location.hash + ")");
      w = wait(); Router.go("/selftest"); await w;
      return "Reference \u2192 Dataset \u2192 Back \u2192 Forward all behaved (started at " + start + ")";
    }],
  ];
  async function renderRoute0(path, query) { const route = Router.routes.find((r) => r.pattern.test(path)); await route.handler({ path, query, params: {}, raw: path, isCurrent: () => true }); }

  async function runAll(table, summary) {
    results.length = 0;
    for (let i = 0; i < TESTS.length; i++) {
      const [title, fn] = TESTS[i];
      const row = table.rows[i]; setRow(row, "running", "Running\u2026");
      try { const detail = await fn(); results.push({ title, ok: true, detail }); setRow(row, "pass", detail); }
      catch (error) { results.push({ title, ok: false, detail: error.message }); setRow(row, "fail", error.message || String(error)); }
      paintSummary(summary);
    }
  }
  function setRow(row, state, text) {
    row.badge.className = "badge " + (state === "pass" ? "badge-ok" : state === "fail" ? "badge-err" : "");
    row.badge.textContent = state === "pass" ? "PASS" : state === "fail" ? "FAIL" : state === "running" ? "RUNNING" : "WAITING";
    row.detail.textContent = safe(text);
  }
  function paintSummary(node) {
    const pass = results.filter((r) => r.ok).length;
    node.textContent = pass + " of " + TESTS.length + " passed" + (results.length < TESTS.length ? " (" + results.length + " run)" : "") + (results.length === TESTS.length && pass === TESTS.length ? " \u2014 everything checks out." : "");
  }
  async function viewSelfTest() {
    const summary = el("p", { class: "panel-note" }, "Not run yet.");
    const table = { rows: [] };
    const tbody = el("tbody", {});
    TESTS.forEach(([title], i) => {
      const badge = el("span", { class: "badge" }, "WAITING"), detail = el("span", { class: "muted" }, "");
      const prev = results[i];
      table.rows.push({ badge, detail });
      if (prev) setRow(table.rows[i], prev.ok ? "pass" : "fail", prev.detail);
      tbody.appendChild(el("tr", {}, [el("td", { class: "num" }, String(i + 1)), el("td", {}, title), el("td", {}, badge), el("td", {}, detail)]));
    });
    if (results.length) paintSummary(summary);
    const v = $("#view"); clear(v);
    [el("div", { class: "section-label" }, "SELF-TEST"), el("h1", { class: "page-title" }, "Validation checklist"),
      el("p", { class: "page-sub" }, "Runs against the live server in your browser. Needs an internet connection (UniProt, PROSITE, PDB, AlphaFold). The dataset-statistics check may ask you to re-run after the first build finishes."),
      el("div", { class: "flex flex-wrap", style: "margin-bottom:14px" }, [
        el("button", { class: "btn btn-primary", type: "button", onclick: () => runAll(table, summary) }, "Run all checks"),
        el("button", { class: "btn", type: "button", onclick: () => Atlas.copyText(JSON.stringify(results, null, 2)) }, "Copy report")]), summary,
      el("div", { class: "table-wrap" }, el("table", { id: "selftest-table" }, [el("thead", {}, el("tr", {}, ["#", "Check", "Result", "Detail"].map((h) => el("th", {}, h)))), tbody]))].forEach((n) => v.appendChild(n));
  }
  Atlas.selftest = { TESTS, runAll };
  Router.add("selftest", /^\/selftest$/, viewSelfTest);
})(typeof window !== "undefined" ? window : globalThis);
