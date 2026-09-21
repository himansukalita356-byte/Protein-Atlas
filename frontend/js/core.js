/* Protein Atlas - core: DOM helper, safe formatting, API client, router, session, state blocks.
   Everything the interface prints goes through these helpers so that a missing value can only ever
   appear as an honest message ("Data unavailable"), never as "undefined" or "NaN". */
(function (global) {
  "use strict";
  const Atlas = global.Atlas = global.Atlas || {};
  const UNAVAILABLE = "Data unavailable";
  Atlas.UNAVAILABLE = UNAVAILABLE;
  Atlas.diagnostics = { scrubbed: [], errors: [] };

  // ------------------------------------------------------------------ DOM
  const $ = (sel, root) => (root || document).querySelector(sel);
  const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));

  function appendChildren(parent, child) {
    if (child === null || child === undefined || child === false || child === true) return;
    if (Array.isArray(child)) { child.forEach((c) => appendChildren(parent, c)); return; }
    if (typeof child === "number") child = Number.isFinite(child) ? String(child) : UNAVAILABLE;
    if (typeof child === "string") { parent.appendChild(document.createTextNode(child)); return; }
    parent.appendChild(child);
  }
  function el(tag, attrs, children) {
    const node = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach((key) => {
        const value = attrs[key];
        if (value === null || value === undefined || value === false) return;
        if (key === "class") node.className = value;
        else if (key === "style") node.setAttribute("style", value);
        else if (key.indexOf("on") === 0 && typeof value === "function") node.addEventListener(key.slice(2), value);
        else if (key === "value" || key === "checked" || key === "disabled" || key === "selected") node[key] = value;
        else node.setAttribute(key, value === true ? "" : String(value));
      });
    }
    appendChildren(node, children);
    return node;
  }
  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); return node; }
  Object.assign(Atlas, { $, $$, el, clear });

  // ------------------------------------------------------------------ formatting
  const isNum = (v) => typeof v === "number" && Number.isFinite(v);
  const fmt = {
    isNum,
    n(v, d) { d = d === undefined ? 2 : d; return isNum(v) ? v.toLocaleString("en-US", { minimumFractionDigits: d, maximumFractionDigits: d }) : UNAVAILABLE; },
    int(v) { return isNum(v) ? Math.round(v).toLocaleString("en-US") : UNAVAILABLE; },
    pct(v, d) { return isNum(v) ? fmt.n(v, d === undefined ? 2 : d) + "%" : UNAVAILABLE; },
    signed(v, d) { return isNum(v) ? (v > 0 ? "+" : "") + fmt.n(v, d === undefined ? 2 : d) : UNAVAILABLE; },
    da(v) { return isNum(v) ? fmt.n(v, 2) + " Da" : UNAVAILABLE; },
    kda(v) { return isNum(v) ? fmt.n(v / 1000, 2) + " kDa" : UNAVAILABLE; },
    text(v) { return v === null || v === undefined || v === "" || (typeof v === "number" && !isNum(v)) ? UNAVAILABLE : String(v); },
    plain(v, d) { return isNum(v) ? String(Math.round(v * Math.pow(10, d)) / Math.pow(10, d)) : ""; },  // for CSV
  };
  Atlas.fmt = fmt;

  // ------------------------------------------------------------------ errors & API
  class AtlasError extends Error {
    constructor(kind, message, status, extra) {
      super(message);
      this.kind = kind; this.status = status || 0; this.extra = extra || {};
    }
  }
  Atlas.AtlasError = AtlasError;
  const KIND_BY_STATUS = { 404: "not_found", 409: "page_unreachable", 422: "invalid_input", 500: "calculation_error", 502: "database_unavailable", 503: "database_unavailable" };
  const FALLBACK_MESSAGE = {
    not_found: "That item could not be found.",
    invalid_input: "That input could not be analysed.",
    database_unavailable: "An external database is unavailable right now.",
    calculation_error: "Analysis unavailable \u2014 calculation error",
    page_unreachable: "That page cannot be reached directly.",
  };
  async function api(path, options) {
    options = options || {};
    const url = path.indexOf("/api/") === 0 ? path : "/api/v1" + path;
    const init = { method: options.method || "GET", headers: {} };
    if (options.body !== undefined) { init.method = options.method || "POST"; init.headers["Content-Type"] = "application/json"; init.body = JSON.stringify(options.body); }
    let response;
    try { response = await fetch(url, init); }
    catch (e) { throw new AtlasError("network", "Cannot reach the Atlas server. Make sure the black server window is still open, then try again.", 0); }
    let payload = null;
    try { const text = await response.text(); payload = text ? JSON.parse(text) : null; } catch (e) { payload = null; }
    if (!response.ok) {
      const kind = (payload && payload.kind) || KIND_BY_STATUS[response.status] || "calculation_error";
      const detail = payload && typeof payload.detail === "string" ? payload.detail : FALLBACK_MESSAGE[kind] || FALLBACK_MESSAGE.calculation_error;
      throw new AtlasError(kind, detail, response.status, payload || {});
    }
    return payload;
  }
  Atlas.api = api;
  const qs = (obj) => {
    const p = new URLSearchParams();
    Object.keys(obj).forEach((k) => { if (obj[k] !== undefined && obj[k] !== null && obj[k] !== "") p.set(k, obj[k]); });
    const s = p.toString();
    return s ? "?" + s : "";
  };
  Atlas.qs = qs;

  // ------------------------------------------------------------------ state blocks
  function stateBlock(glyph, title, body, extra) {
    return el("div", { class: "state-block", role: "status" }, [
      el("div", { class: "glyph" }, glyph), el("h3", {}, title), body ? el("p", {}, body) : null, extra || null]);
  }
  const ERROR_TITLES = {
    invalid_input: ["\u270e", "Check your input"], not_found: ["\ud83e\uddec", "Not found"],
    database_unavailable: ["\u26a0", "Database unavailable"], calculation_error: ["\u26a0", "Analysis unavailable \u2014 calculation error"],
    network: ["\ud83d\udd0c", "Server unreachable"], page_unreachable: ["\u21c4", "That page is too far to jump to"],
    page_out_of_range: ["\u21c4", "Page does not exist"],
  };
  function errorBlock(error, retry) {
    const kind = error && error.kind ? error.kind : "calculation_error";
    const meta = ERROR_TITLES[kind] || ERROR_TITLES.calculation_error;
    const message = error && error.kind ? error.message : FALLBACK_MESSAGE.calculation_error;
    if (!(error && error.kind)) Atlas.diagnostics.errors.push(String(error && error.stack || error));
    const actions = retry ? el("button", { class: "btn", onclick: retry }, "Try again") : null;
    return stateBlock(meta[0], meta[1], message, actions);
  }
  const loadingBlock = (text) => el("div", { class: "state-block", role: "status" }, [
    el("div", { class: "spinner" }), el("p", {}, text || "Loading\u2026")]);
  const noDataBlock = (text) => stateBlock("\u2205", "No data", text || UNAVAILABLE);
  const skeleton = (h) => el("div", { class: "skeleton", style: "height:" + (h || 200) + "px;border-radius:14px;" });
  Atlas.states = { stateBlock, errorBlock, loadingBlock, noDataBlock, skeleton };

  // ------------------------------------------------------------------ misc helpers
  function toast(message, kind) {
    const stack = $("#toast-stack");
    if (!stack) return;
    const node = el("div", { class: "toast " + (kind === "error" ? "error" : "") }, String(message));
    stack.appendChild(node);
    setTimeout(() => node.remove(), 4200);
  }
  function download(filename, text, mime) {
    const blob = new Blob([text], { type: mime || "text/plain" });
    const url = URL.createObjectURL(blob);
    const link = el("a", { href: url, download: filename });
    document.body.appendChild(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  async function copyText(text) {
    try { await navigator.clipboard.writeText(text); toast("Copied to clipboard."); }
    catch (e) {
      const area = el("textarea", { style: "position:fixed;opacity:0;" }); area.value = text;
      document.body.appendChild(area); area.select();
      try { document.execCommand("copy"); toast("Copied to clipboard."); } catch (e2) { toast("Copy failed \u2014 select the text and press Ctrl+C.", "error"); }
      area.remove();
    }
  }
  const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
  function debounce(fn, ms) { let t = null; return function () { const a = arguments; clearTimeout(t); t = setTimeout(() => fn.apply(null, a), ms); }; }
  const CLASS_COLOR = { Nonpolar: "#d7cd8a", Aromatic: "#9c8cf2", Polar: "#57d9c4", Acidic: "#e2725b", Basic: "#5b9ee2", "Non-standard": "#74897f" };
  const classPill = (cls) => el("span", { class: "pill pill-" + String(cls || "").toLowerCase().replace(/[^a-z]/g, "") }, fmt.text(cls));
  Object.assign(Atlas, { toast, download, copyText, sleep, debounce, CLASS_COLOR, classPill });

  // last line of defence: never let "undefined"/"NaN" reach the screen
  function scrub(root) {
    if (!root) return 0;
    let fixed = 0;
    const walk = (node) => {
      if (node.nodeType === 3) {
        const text = node.nodeValue || "";
        if (/\bundefined\b|\bNaN\b|^null$|\[object Object\]/.test(text)) {
          Atlas.diagnostics.scrubbed.push(text.slice(0, 80));
          node.nodeValue = text.replace(/\bundefined\b|\bNaN\b|\[object Object\]/g, UNAVAILABLE).replace(/^null$/, UNAVAILABLE);
          fixed++;
        }
        return;
      }
      const kids = node.childNodes || [];
      for (let i = 0; i < kids.length; i++) walk(kids[i]);
    };
    walk(root);
    return fixed;
  }
  Atlas.scrub = scrub;

  // ------------------------------------------------------------------ tables
  /* columns: [{key, label, align: 'num'|'mono'|'', render(row, i) -> node|string}] */
  function table(columns, rows, options) {
    options = options || {};
    const head = el("tr", {}, columns.map((c) => el("th", { class: c.align === "num" ? "num-h" : "" , style: c.align === "num" ? "text-align:right" : null }, c.label)));
    const body = rows.map((row, i) => el("tr", {}, columns.map((c) => {
      const raw = c.render ? c.render(row, i) : row[c.key];
      const content = typeof raw === "number" ? (Number.isFinite(raw) ? String(raw) : UNAVAILABLE) : (raw === null || raw === undefined ? UNAVAILABLE : raw);
      return el("td", { class: c.align === "num" ? "num" : c.align === "mono" ? "mono" : "" }, content);
    })));
    const t = el("table", options.id ? { id: options.id } : {}, [el("thead", {}, head), el("tbody", {}, body)]);
    return el("div", { class: "table-wrap" }, t);
  }
  Atlas.table = table;

  // ------------------------------------------------------------------ session
  const Session = {
    items: new Map(),
    add(record) {
      if (!record || !record.id) return;
      const prev = this.items.get(record.id) || {};
      this.items.set(record.id, Object.assign({}, prev, { id: record.id, name: record.name, gene: record.gene, organism: record.organism, length: record.length }));
      this.save(); this.paint();
    },
    get(id) { return this.items.get(id); },
    setAnalysis(id, analysis, record) {
      const prev = this.items.get(id) || { id };
      this.items.set(id, Object.assign({}, prev, { analysis }, record ? { record } : {}));
    },
    list() { return Array.from(this.items.values()); },
    save() {
      try { sessionStorage.setItem("atlas.session", JSON.stringify(this.list().map((x) => ({ id: x.id, name: x.name, gene: x.gene, organism: x.organism, length: x.length })))); } catch (e) { /* storage may be blocked */ }
    },
    load() {
      try { (JSON.parse(sessionStorage.getItem("atlas.session") || "[]")).forEach((x) => { if (x && x.id) this.items.set(x.id, x); }); } catch (e) { /* ignore */ }
    },
    paint() {
      const count = $("#session-count"); if (count) count.textContent = String(this.items.size);
      const list = $("#session-list"); if (!list) return;
      clear(list);
      if (!this.items.size) list.appendChild(el("p", { class: "muted" }, "Nothing opened yet."));
      this.list().forEach((p) => list.appendChild(el("a", { class: "session-item", href: "#/protein/" + encodeURIComponent(p.id) }, [
        el("b", {}, p.id), el("span", {}, fmt.text(p.gene || p.name))])));
    },
  };
  Atlas.Session = Session;

  // ------------------------------------------------------------------ router
  const routes = [];
  const cleanups = [];
  const scrollMemory = {};
  const Router = Atlas.Router = { routes, current: null, pushed: false, token: 0 };
  Atlas.onCleanup = (fn) => cleanups.push(fn);
  Router.cleanup = () => runCleanups();
  function runCleanups() { while (cleanups.length) { try { cleanups.pop()(); } catch (e) { Atlas.diagnostics.errors.push(String(e)); } } }

  Router.add = (name, pattern, handler, soft) => routes.push({ name, pattern, handler, soft });
  Router.parse = () => {
    const raw = (global.location.hash || "").replace(/^#/, "") || "/";
    const cut = raw.indexOf("?");
    const path = cut === -1 ? raw : raw.slice(0, cut);
    const query = {};
    new URLSearchParams(cut === -1 ? "" : raw.slice(cut + 1)).forEach((v, k) => { query[k] = v; });
    return { raw, path: path || "/", query };
  };
  Router.href = (path, query) => "#" + path + qs(query || {});
  Router.go = (path, query, options) => {
    options = options || {};
    const target = Router.href(path, query);
    if (target === (global.location.hash || "#/")) { return Router.render(); }
    if (options.replace) { global.history.replaceState(null, "", target); return Router.render(); }
    if (query && query.view && !(Router.current && Router.current.query && Router.current.query.view)) Atlas.viewerPushed = true;
    Router.pushed = true;
    global.location.hash = target;
    return Promise.resolve();
  };
  Router.back = (fallbackPath) => {
    if (global.history.length > 1) global.history.back(); else Router.go(fallbackPath || "/");
  };
  function setActiveNav(name) {
    $$(".topnav a").forEach((a) => a.classList.toggle("active", a.getAttribute("data-route") === name));
  }
  Router.render = async function () {
    const token = ++Router.token;
    const location = Router.parse();
    const wasPush = Router.pushed; Router.pushed = false;
    let matched = null, groups = {};
    for (const r of routes) {
      if (r.name === "notfound") continue;
      const m = r.pattern.exec(location.path);
      if (m) { matched = r; groups = m.groups || {}; break; }
    }
    const ctx = { path: location.path, query: location.query, params: groups, raw: location.raw, token, isCurrent: () => token === Router.token };
    Router.current = ctx;
    if (!matched) matched = routes.find((r) => r.name === "notfound");
    setActiveNav(matched.name);
    if (matched.soft && matched.soft(ctx)) return;              // e.g. only the tab changed on the same protein page
    runCleanups();
    const view = $("#view");
    try { await matched.handler(ctx); }
    catch (error) {
      if (ctx.isCurrent()) { clear(view); view.appendChild(errorBlock(error, () => Router.render())); }
    }
    if (ctx.isCurrent()) {
      scrub(view);
      if (!wasPush && scrollMemory[location.raw] !== undefined) global.scrollTo(0, scrollMemory[location.raw]);
      else if (wasPush || !scrollMemory[location.raw]) global.scrollTo(0, 0);
    }
  };
  Router.start = () => {
    global.addEventListener("hashchange", () => Router.render());
    global.addEventListener("scroll", debounce(() => { scrollMemory[Router.parse().raw] = global.scrollY || 0; }, 80), { passive: true });
    return Router.render();
  };
})(typeof window !== "undefined" ? window : globalThis);
