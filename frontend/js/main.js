/* Protein Atlas - start-up. */
(function (global) {
  "use strict";
  const Atlas = global.Atlas;
  const { $, el, fmt, api, Router, Session } = Atlas;

  function closeViewer() {
    const cur = Router.current;
    if (cur && cur.query && cur.query.view) {
      if (Atlas.viewerPushed) { Atlas.viewerPushed = false; global.history.back(); }
      else { const q = Object.assign({}, cur.query); delete q.view; delete q.src; Router.go(cur.path, q, { replace: true }); }
    } else Atlas.structure.close();
  }

  async function pollDatabaseStatus() {
    const node = $("#db-status"); if (!node) return;
    try {
      const s = await api("/prosite/status");
      const text = s.state === "ready" ? "PROSITE ready \u00b7 " + fmt.int(s.patterns_usable) + " patterns" + (s.release ? " (release " + s.release + ")" : "")
        : s.state === "failed" ? "PROSITE database unavailable" : "PROSITE database downloading\u2026";
      node.textContent = text;
      if (s.state === "loading" || s.state === "idle") setTimeout(pollDatabaseStatus, 5000);
    } catch (e) { node.textContent = ""; }
  }

  function boot() {
    Session.load(); Session.paint();
    const form = $("#global-search");
    if (form) form.addEventListener("submit", (e) => { e.preventDefault(); Router.go("/", { q: $("#global-search-input").value.trim() }); });
    const drawer = $("#session-drawer");
    $("#session-toggle").addEventListener("click", () => drawer.classList.toggle("hidden"));
    $("#session-close").addEventListener("click", () => drawer.classList.add("hidden"));
    document.addEventListener("click", (e) => {
      const a = e.target && e.target.closest ? e.target.closest("a[href^='#/']") : null;
      if (a) { Router.pushed = true; drawer.classList.add("hidden"); }
    });
    Atlas.structure.init(closeViewer);
    Router.start();
    pollDatabaseStatus();
  }
  Atlas.boot = boot;
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})(typeof window !== "undefined" ? window : globalThis);
