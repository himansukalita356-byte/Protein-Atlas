/* Protein Atlas - interactive charts (Chart.js) with wheel-zoom, drag-pan, reset and tooltips.
   Zoom/pan is implemented here (no plugin needed).  The wheel is captured ONLY while the pointer is
   over a chart that has wheel-zoom switched on, so normal page scrolling works everywhere else. */
(function (global) {
  "use strict";
  const Atlas = global.Atlas;
  const { el, fmt } = Atlas;
  const C = { teal: "#57d9c4", amber: "#e8a33d", coral: "#e2725b", violet: "#9c8cf2", blue: "#5b9ee2", sand: "#d7cd8a", dim: "#a9bdb4", grid: "rgba(255,255,255,0.06)" };
  Atlas.COLORS = C;
  const available = () => typeof global.Chart !== "undefined" && !!global.Chart;
  let defaultsSet = false;

  function setDefaults() {
    if (defaultsSet || !available() || !global.Chart.defaults) return;
    try {
      global.Chart.defaults.color = C.dim;
      global.Chart.defaults.borderColor = C.grid;
      if (global.Chart.defaults.font) global.Chart.defaults.font.family = "'Inter', system-ui, sans-serif";
    } catch (e) { /* cosmetic only */ }
    defaultsSet = true;
  }

  const bandsPlugin = {
    id: "atlasBands",
    beforeDatasetsDraw(chart) {
      const opts = chart.options && chart.options.plugins && chart.options.plugins.atlasBands;
      const bands = (opts && opts.bands) || [];
      if (!bands.length || !chart.scales || !chart.scales.x) return;
      const { ctx, chartArea } = chart; const x = chart.scales.x;
      bands.forEach((b) => {
        const a = x.getPixelForValue(b.from), z = x.getPixelForValue(b.to);
        const left = Math.max(chartArea.left, Math.min(a, z)), right = Math.min(chartArea.right, Math.max(a, z));
        if (right > left) { ctx.save(); ctx.fillStyle = b.color; ctx.fillRect(left, chartArea.top, right - left, chartArea.bottom - chartArea.top); ctx.restore(); }
      });
    },
  };

  function extent(datasets) {
    let lo = Infinity, hi = -Infinity;
    datasets.forEach((d) => (d.data || []).forEach((p) => { const v = typeof p === "object" && p !== null ? p.x : NaN; if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } }));
    return Number.isFinite(lo) ? [lo, hi] : null;
  }

  function buildConfig(kind, spec) {
    const category = kind === "bar" || kind === "doughnut";
    const datasets = spec.datasets || [];
    const options = {
      responsive: true, maintainAspectRatio: false, animation: spec.animate === false ? false : { duration: 250 },
      interaction: { mode: kind === "scatter" ? "nearest" : "index", intersect: false },
      plugins: {
        legend: { display: spec.legend !== false && (datasets.length > 1 || kind === "doughnut"), labels: { color: C.dim, boxWidth: 12 } },
        tooltip: { callbacks: spec.tooltip || {} },
        atlasBands: { bands: spec.bands || [] },
      },
    };
    if (kind === "doughnut") {
      options.cutout = "58%";
      return { type: "doughnut", data: { labels: spec.labels, datasets }, options, plugins: [] };
    }
    const many = datasets.some((d) => (d.data || []).length > 2500);
    if (kind === "line" && many) {
      options.parsing = false; options.normalized = true;
      options.plugins.decimation = { enabled: true, algorithm: "min-max" };
    }
    const range = category ? null : (spec.xRange || extent(datasets));
    options.scales = {
      x: Object.assign({
        type: category ? "category" : (spec.xType || "linear"),
        title: { display: !!spec.xTitle, text: spec.xTitle || "", color: C.dim },
        grid: { color: C.grid }, ticks: { color: C.dim, maxRotation: spec.rotate ? 60 : 0, autoSkip: true },
      }, range ? { min: range[0], max: range[1] } : {}),
      y: Object.assign({
        type: spec.yType || "linear", beginAtZero: spec.beginAtZero !== false && kind === "bar",
        title: { display: !!spec.yTitle, text: spec.yTitle || "", color: C.dim },
        grid: { color: C.grid }, ticks: { color: C.dim },
      }, spec.yMin !== undefined ? { min: spec.yMin } : {}, spec.yMax !== undefined ? { max: spec.yMax } : {}),
    };
    if (spec.horizontal) options.indexAxis = "y";
    return { type: kind, data: { labels: spec.labels, datasets }, options, plugins: [bandsPlugin] };
  }

  function attachZoomPan(chart, canvas, tools, spec, category) {
    const state = { wheel: spec.wheelZoom !== false, full: null, drag: null };
    const s0 = () => chart.scales && chart.scales.x;
    function fullRange() {
      if (state.full) return state.full;
      const s = s0(); if (!s || !Number.isFinite(s.min) || !Number.isFinite(s.max)) return null;
      state.full = [s.min, s.max]; return state.full;
    }
    function readout() {
      const s = s0(), f = fullRange(); if (!s || !f || !tools.readout) return;
      const zoomed = s.min > f[0] + 1e-9 || s.max < f[1] - 1e-9;
      tools.readout.textContent = zoomed ? ("Showing " + fmt.n(category ? Math.round(s.min) + 1 : s.min, category ? 0 : 1) + " \u2013 " + fmt.n(category ? Math.round(s.max) + 1 : s.max, category ? 0 : 1)) : "Full range";
      canvas.style.cursor = zoomed ? "grab" : "default";
    }
    function setRange(lo, hi) {
      const f = fullRange(); if (!f) return;
      const minSpan = category ? 4 : (f[1] - f[0]) / 500;
      let span = Math.max(hi - lo, minSpan);
      span = Math.min(span, f[1] - f[0]);
      const mid = (lo + hi) / 2;
      lo = mid - span / 2; hi = mid + span / 2;
      if (lo < f[0]) { hi += f[0] - lo; lo = f[0]; }
      if (hi > f[1]) { lo -= hi - f[1]; hi = f[1]; }
      chart.options.scales.x.min = Math.max(lo, f[0]); chart.options.scales.x.max = Math.min(hi, f[1]);
      chart.update("none"); readout();
    }
    function zoomBy(factor, centre) {
      const s = s0(); if (!s) return;
      const c = centre === undefined ? (s.min + s.max) / 2 : centre;
      setRange(c - (c - s.min) * factor, c + (s.max - c) * factor);
    }
    function reset() {
      const f = fullRange(); if (!f) return;
      chart.options.scales.x.min = f[0]; chart.options.scales.x.max = f[1];
      chart.update("none"); readout();
    }
    canvas.style.touchAction = "pan-y";
    canvas.addEventListener("wheel", (e) => {
      if (!state.wheel) return;                     // wheel lock: the page scrolls normally
      const s = s0(); if (!s) return;
      e.preventDefault();
      const dy = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY;
      const rect = canvas.getBoundingClientRect();
      const centre = s.getValueForPixel(e.clientX - rect.left);
      zoomBy(Math.exp(dy * 0.0016), Number.isFinite(centre) ? centre : undefined);
    }, { passive: false });
    canvas.addEventListener("pointerdown", (e) => {
      const s = s0(); if (!s || e.button > 0) return;
      state.drag = { x: e.clientX, min: s.min, max: s.max };
      try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* not all browsers */ }
    });
    canvas.addEventListener("pointermove", (e) => {
      const d = state.drag, s = s0(); if (!d || !s) return;
      const f = fullRange(); if (!f || (d.min <= f[0] + 1e-9 && d.max >= f[1] - 1e-9)) return;
      const width = (chart.chartArea ? chart.chartArea.right - chart.chartArea.left : 0) || 1;
      const shift = -((e.clientX - d.x) * (d.max - d.min)) / width;
      setRange(d.min + shift, d.max + shift);
      canvas.style.cursor = "grabbing";
    });
    const end = (e) => { if (state.drag) { state.drag = null; try { canvas.releasePointerCapture(e.pointerId); } catch (err) { /* ignore */ } readout(); } };
    canvas.addEventListener("pointerup", end); canvas.addEventListener("pointercancel", end);
    canvas.addEventListener("dblclick", reset);
    if (tools.zoomIn) tools.zoomIn.addEventListener("click", () => zoomBy(0.6));
    if (tools.zoomOut) tools.zoomOut.addEventListener("click", () => zoomBy(1 / 0.6));
    if (tools.reset) tools.reset.addEventListener("click", reset);
    if (tools.wheel) tools.wheel.addEventListener("click", () => {
      state.wheel = !state.wheel;
      tools.wheel.textContent = state.wheel ? "Wheel zoom: on" : "Wheel zoom: off";
      tools.wheel.classList.toggle("active", state.wheel);
    });
    return { zoomBy, reset, state, setRange, readout };
  }

  /* kind: 'line' | 'bar' | 'scatter' | 'doughnut'.  Returns a handle; the chart is created right after
     the caller has put `container` in the page (next tick) so Chart.js can measure it. */
  function mount(container, kind, spec) {
    spec = spec || {};
    if (!available()) {
      container.appendChild(Atlas.states.stateBlock("\ud83d\udcca", "Chart library unavailable",
        "The Chart.js library could not be downloaded, so this graph cannot be drawn. The numbers are still listed in the table. Check your internet connection or ad-blocker and refresh with Ctrl+F5."));
      return { chart: null, destroy() {}, unavailable: true };
    }
    setDefaults();
    const interactive = spec.zoom !== false && kind !== "doughnut";
    const canvas = el("canvas", { role: "img", "aria-label": spec.ariaLabel || "Interactive chart" });
    const tools = {};
    let toolbar = null;
    if (interactive) {
      tools.zoomIn = el("button", { class: "btn btn-ghost chart-tool", type: "button", title: "Zoom in" }, "+");
      tools.zoomOut = el("button", { class: "btn btn-ghost chart-tool", type: "button", title: "Zoom out" }, "\u2212");
      tools.reset = el("button", { class: "btn btn-ghost chart-tool", type: "button", title: "Reset zoom (or double-click the chart)" }, "Reset");
      tools.wheel = el("button", { class: "btn btn-ghost chart-tool active", type: "button", title: "When on, the mouse wheel zooms while the pointer is over this chart" }, "Wheel zoom: on");
      tools.readout = el("span", { class: "muted chart-readout" }, "Full range");
      toolbar = el("div", { class: "chart-tools" }, [tools.zoomIn, tools.zoomOut, tools.reset, tools.wheel, tools.readout]);
    }
    const box = el("div", { class: "chart-box " + (spec.tall ? "tall" : "") }, canvas);
    const wrapper = el("div", { class: "chart-wrap" }, [toolbar, box,
      interactive ? el("div", { class: "chart-legend-note" }, "Hover for exact values \u00b7 scroll over the chart to zoom \u00b7 drag to pan \u00b7 double-click to reset") : null,
      spec.note ? el("div", { class: "chart-legend-note" }, spec.note) : null]);
    container.appendChild(wrapper);
    const handle = { canvas, chart: null, zoom: null, destroy() { if (handle.chart) { try { handle.chart.destroy(); } catch (e) { /* already gone */ } handle.chart = null; } }, spec, kind };
    const timer = setTimeout(() => {
      try {
        handle.chart = new global.Chart(canvas, buildConfig(kind, spec));
        if (interactive) handle.zoom = attachZoomPan(handle.chart, canvas, tools, spec, kind === "bar");
      } catch (error) {
        Atlas.diagnostics.errors.push("chart: " + (error && error.message));
        clearChildren(box); box.appendChild(Atlas.states.stateBlock("\ud83d\udcca", "This chart could not be drawn", "The data is listed in the table."));
      }
    }, 0);
    Atlas.charts.registry.push(handle);
    Atlas.onCleanup(() => { clearTimeout(timer); handle.destroy(); const i = Atlas.charts.registry.indexOf(handle); if (i >= 0) Atlas.charts.registry.splice(i, 1); });
    return handle;
  }
  function clearChildren(node) { while (node.firstChild) node.removeChild(node.firstChild); }

  Atlas.charts = {
    registry: [], available, mount, buildConfig,
    line: (c, s) => mount(c, "line", s), bar: (c, s) => mount(c, "bar", s),
    scatter: (c, s) => mount(c, "scatter", s), doughnut: (c, s) => mount(c, "doughnut", s),
    /* engine histogram {edges, labels, counts} -> bar chart */
    histogram(container, hist, spec) {
      const labels = (hist && hist.labels) || [];
      const counts = (hist && hist.counts) || [];
      if (!labels.length || !counts.length) { container.appendChild(Atlas.states.noDataBlock("No distribution is available for this quantity.")); return null; }
      return mount(container, "bar", Object.assign({
        labels, rotate: true, legend: false,
        datasets: [{ label: (spec && spec.label) || "Proteins", data: counts, backgroundColor: (spec && spec.color) || C.teal, borderWidth: 0, barPercentage: 1, categoryPercentage: 0.95 }],
        yTitle: "Number of proteins",
      }, spec || {}));
    },
  };
})(typeof window !== "undefined" ? window : globalThis);
