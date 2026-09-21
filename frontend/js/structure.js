/* Protein Atlas - interactive molecular viewer (3Dmol.js).
   * MolViewer      : thin wrapper (styles, colours, picking, deep-zoom "detail lens", screenshots)
   * Atlas.structure: full-screen workspace for proteins (experimental -> AlphaFold -> SWISS-MODEL -> ESMFold)
                      and for the 20 amino acids.  Open/closed state lives in the URL, so browser Back closes it. */
(function (global) {
  "use strict";
  const Atlas = global.Atlas;
  const { el, $, clear, fmt, api } = Atlas;
  const mol3d = () => (typeof global.$3Dmol !== "undefined" ? global.$3Dmol : null);
  const WATER = ["HOH", "WAT", "H2O", "DOD"];
  const ATOM_LIMIT = 25000;          // above this, "All atoms" and the detail lens are switched off to keep the viewer responsive

  function plddtColor(atom) {
    const b = atom.b;
    if (!(b >= 0)) return 0x9aa5a1;
    return b > 90 ? 0x0053d6 : b > 70 ? 0x65cbf3 : b > 50 ? 0xffdb13 : 0xff7d45;
  }

  class MolViewer {
    constructor(host, options) {
      const lib = mol3d();
      if (!lib) throw new Error("3Dmol library missing");
      this.host = host; this.small = !!(options && options.small);
      this.viewer = lib.createViewer(host, { backgroundColor: "#060b09", antialias: true });
      this.style = this.small ? "ballstick" : "cartoon";
      this.color = this.small ? "element" : "spectrum";
      this.showAtoms = false; this.showWater = false; this.showLabels = false; this.spinning = false;
      this.lens = true; this.lensOn = false; this.lensPinned = false; this.picked = null; this.lensLabels = []; this.pickLabel = null; this.labels = [];
      this.onPick = null; this.onLens = null; this.z0 = null; this.model = null; this.atoms = [];
      this.waterCount = 0; this.destroyed = false; this._lensBound = [];
      this.backbone = new Set();
      this._bindLens();
    }
    load(text, format) {
      const v = this.viewer;
      v.clear();
      this._clearLensLabels();
      this.model = v.addModel(text, format);
      this.atoms = this.model.selectedAtoms({});
      this.format = format;
      this.hasBackbone = this.atoms.some((a) => a.atom === "CA");
      this.waterCount = this.atoms.reduce((n, a) => n + (WATER.indexOf(a.resn) >= 0 ? 1 : 0), 0);
      if (!this.hasBackbone && (this.style === "cartoon")) this.style = "stick";
      this.picked = null; this.lensOn = false; this.lensPinned = false;
      this._indexBackbone();
      this.applyStyle();
      v.zoomTo(); v.render();
      this._rememberZoom();
      if (v.setClickable) {
        v.setClickable({}, true, (atom) => { if (atom) this.pick(atom); });
      }
      if (this.small && v.setHoverable) {
        v.setHoverable({}, true, (atom) => { if (this.onHover) this.onHover(atom); }, () => { if (this.onHover) this.onHover(null); });
      }
      return this.atoms.length;
    }
    _indexBackbone() {
      this.backbone = new Set();
      if (!this.small) return;
      const carboxyl = new Set(); // amino acid model: backbone = N, CA-like C alpha, carboxyl C and both O (+ H on them)
      const byIndex = this.atoms;
      const bonds = {};
      byIndex.forEach((a, i) => { (a.bonds || []).forEach((j) => { (bonds[i] = bonds[i] || []).push(j); }); });
      byIndex.forEach((a, i) => {
        if (a.elem !== "C") return;
        const oxy = (bonds[i] || []).filter((j) => byIndex[j] && byIndex[j].elem === "O");
        if (oxy.length === 2) { carboxyl.add(i); oxy.forEach((j) => carboxyl.add(j)); }
      });
      carboxyl.forEach((i) => {
        (bonds[i] || []).forEach((j) => {
          const n = byIndex[j];
          if (n && n.elem === "C" && !carboxyl.has(j)) {          // alpha carbon
            addAlpha(j);
          }
        });
      });
      function addAlpha(alpha) {
        carboxyl.add(alpha);
        (bonds[alpha] || []).forEach((k) => { const n = byIndex[k]; if (n && n.elem === "N") carboxyl.add(k); });
      }
      carboxyl.forEach((i) => { this.backbone.add(i); (bonds[i] || []).forEach((j) => { if (byIndex[j] && byIndex[j].elem === "H") this.backbone.add(j); }); });
    }
    _colorSpec(scheme) {
      switch (scheme) {
        case "chain": return { colorscheme: "chain" };
        case "ss": return { colorscheme: "ssPyMol" };
        case "confidence": return { colorfunc: plddtColor };
        case "element": return { colorscheme: "Jmol" };
        case "part": return { colorfunc: (atom) => (this.backbone.has(atom.index) ? 0x57d9c4 : 0xe8a33d) };
        default: return { color: "spectrum" };
      }
    }
    _atomColor() {           // colouring for atom-level representations (never "spectrum")
      if (this.color === "spectrum") return { colorscheme: "Jmol" };
      if (this.color === "ss") return { colorscheme: "Jmol" };
      return this._colorSpec(this.color);
    }
    applyStyle() {
      const v = this.viewer; if (!v || !this.model) return;
      v.removeAllSurfaces();
      const atomColor = this._atomColor();
      const big = this.atoms.length > ATOM_LIMIT;
      let spec;
      switch (this.style) {
        case "cartoon": spec = { cartoon: Object.assign({}, this._colorSpec(this.color)) }; break;
        case "stick": spec = { stick: Object.assign({ radius: 0.16 }, atomColor) }; break;
        case "sphere": spec = { sphere: Object.assign({}, atomColor) }; break;
        case "line": spec = { line: Object.assign({}, atomColor) }; break;
        case "surface": spec = { cartoon: Object.assign({ opacity: 0.6 }, this._colorSpec(this.color)) }; break;
        default: spec = { stick: Object.assign({ radius: 0.13 }, atomColor), sphere: Object.assign({ scale: 0.26 }, atomColor) };
      }
      v.setStyle({}, spec);
      if (!this.small) {
        v.setStyle({ resn: WATER }, this.showWater ? { sphere: { radius: 0.28, color: 0x5b9ee2 } } : {});
        v.addStyle({ hetflag: true, not: { resn: WATER } }, { stick: { radius: 0.2, colorscheme: "Jmol" } });
      }
      if ((this.showAtoms || this.lensOn) && this.style !== "stick" && this.style !== "ballstick" && !big) {
        v.addStyle({ not: { resn: WATER } }, { stick: Object.assign({ radius: 0.12 }, atomColor) });
      }
      if (this.style === "surface" && !this.small && this.atoms.length < 60000 && mol3d().SurfaceType) {
        const surfaceColor = this.color === "confidence" ? { colorfunc: plddtColor } : this.color === "chain" ? { colorscheme: "chain" } : { color: "white" };
        v.addSurface(mol3d().SurfaceType.VDW, Object.assign({ opacity: 0.8 }, surfaceColor), { hetflag: false });
      }
      if (this.picked && !this.small) {
        const sel = { chain: this.picked.chain, resi: this.picked.resi };
        v.addStyle(sel, { stick: { radius: 0.24, colorscheme: "Jmol" } });
      }
      this._paintLabels();
      v.render();
    }
    _paintLabels() {
      const v = this.viewer;
      this.labels.forEach((l) => { try { v.removeLabel(l); } catch (e) { /* ignore */ } });
      this.labels = [];
      if (this.showLabels && this.small) {
        const counters = {};
        this.atoms.forEach((a) => {
          if (a.elem === "H") return;
          counters[a.elem] = (counters[a.elem] || 0) + 1;
          this.labels.push(v.addLabel(a.elem + counters[a.elem], { position: { x: a.x, y: a.y, z: a.z }, fontSize: 12, fontColor: "#eaf1ec", backgroundColor: "#0c1210", backgroundOpacity: 0.7, showBackground: true, inFront: true }));
        });
      }
      if (this.pickLabel) { try { v.removeLabel(this.pickLabel); } catch (e) { /* ignore */ } this.pickLabel = null; }
      if (this.picked && !this.small) {
        this.pickLabel = v.addLabel(this.picked.resn + " " + this.picked.resi + (this.picked.chain ? " (" + this.picked.chain + ")" : ""),
          { position: { x: this.picked.x, y: this.picked.y, z: this.picked.z }, fontSize: 13, fontColor: "#eaf1ec", backgroundColor: "#0c1210", backgroundOpacity: 0.85, showBackground: true, inFront: true });
      }
    }
    setStyle(name) { this.style = name; this.applyStyle(); }
    setColor(name) { this.color = name; this.applyStyle(); }
    toggle(name, value) { this[name] = value === undefined ? !this[name] : !!value; if (name === "spinning") this.viewer.spin(this.spinning); else this.applyStyle(); return this[name]; }
    // ---- All atoms / Waters / Detail lens / zoom: these act on the structure that is already loaded (no network)
    toggleAtoms() {
      if (!this.showAtoms && this.atoms.length > ATOM_LIMIT) {
        Atlas.toast("This structure has " + fmt.int(this.atoms.length) + " atoms; \u201cAll atoms\u201d is limited to " + fmt.int(ATOM_LIMIT) + " to keep the viewer responsive.");
        return false;
      }
      this.showAtoms = !this.showAtoms; this.applyStyle(); return this.showAtoms;
    }
    toggleWater() {
      if (!this.showWater && !this.waterCount) { Atlas.toast("This structure contains no water molecules to show."); return false; }
      this.showWater = !this.showWater; this.applyStyle(); return this.showWater;
    }
    toggleLens() {
      if (this.small) return false;
      if (this.lensOn) {                              // details are showing -> switch the lens off
        this.lens = false; this.lensPinned = false; this.lensOn = false;
        this._clearLensLabels(); this.applyStyle(); if (this.onLens) this.onLens(false);
        return false;
      }
      if (this.atoms.length > ATOM_LIMIT) { Atlas.toast("The detail lens is limited to structures of up to " + fmt.int(ATOM_LIMIT) + " atoms."); return false; }
      this.lens = true; this.lensPinned = true;       // not showing -> show atoms, bonds and labels now, close to the region of interest
      this._focusDetail();
      return true;
    }
    _focusTarget() {
      if (this.picked) return this.picked;
      const ca = this.atoms.filter((a) => a.atom === "CA");
      if (!ca.length) return null;
      let view = null; try { view = this.viewer.getView(); } catch (e) { view = null; }
      if (!Array.isArray(view) || !view.slice(0, 3).every(Number.isFinite)) return ca[Math.floor(ca.length / 2)];
      const c = { x: -view[0], y: -view[1], z: -view[2] };
      let best = ca[0], bestD = Infinity;
      ca.forEach((a) => { const d = (a.x - c.x) * (a.x - c.x) + (a.y - c.y) * (a.y - c.y) + (a.z - c.z) * (a.z - c.z); if (d < bestD) { bestD = d; best = a; } });
      return best;
    }
    _focusDetail() {
      this.lensOn = true; this.applyStyle(); if (this.onLens) this.onLens(true);
      let far = true; try { const view = this.viewer.getView(); far = !this.z0 || Math.abs(view[3]) / this.z0 > 0.6; } catch (e) { far = true; }
      const t = far ? this._focusTarget() : null;
      if (t) {
        const sel = { resi: (t.resi - 4) + "-" + (t.resi + 4) };
        if (t.chain && t.chain !== " ") sel.chain = t.chain;
        try { this.viewer.zoomTo(sel, 500); } catch (e) { /* keep the current view */ }
      }
      this._afterCamera(far ? 700 : 0);
    }
    zoomBy(factor) {
      try { this.viewer.zoom(factor, 200); } catch (e) { return; }
      this._afterCamera(280);
    }
    _afterCamera(ms) {
      clearTimeout(this._camTimer);
      this._camTimer = setTimeout(() => { if (!this.destroyed) { this.updateLens(); try { this.viewer.render(); } catch (e) { /* hidden */ } } }, ms);
    }
    _clearLensLabels() {
      (this.lensLabels || []).forEach((l) => { try { this.viewer.removeLabel(l); } catch (e) { /* ignore */ } });
      this.lensLabels = [];
    }
    pick(atom) {
      this.picked = atom; this.applyStyle();
      if (this.onPick) this.onPick(atom);
    }
    gotoResidue(resi, chain) {
      const sel = chain ? { resi, chain } : { resi };
      const found = this.model.selectedAtoms(sel);
      if (!found.length) return false;
      const ca = found.find((a) => a.atom === "CA") || found[0];
      this.picked = ca; this.applyStyle();
      this.viewer.zoomTo(sel, 600);
      if (this.onPick) this.onPick(ca);
      this._afterCamera(700);                          // let the detail lens react once the camera has arrived
      return true;
    }
    resetCamera() {
      this.picked = null; this.lensOn = false; this.lensPinned = false;
      this._clearLensLabels();
      this.viewer.zoomTo(); this._rememberZoom(); this.applyStyle();
      if (this.onLens) this.onLens(false);
    }
    _rememberZoom() {
      try { const view = this.viewer.getView(); this.z0 = Array.isArray(view) ? Math.abs(view[3]) : null; } catch (e) { this.z0 = null; }
    }
    png() { try { return this.viewer.pngURI(); } catch (e) { return null; } }
    resize() { try { this.viewer.resize(); this.viewer.render(); } catch (e) { /* hidden */ } }
    destroy() {
      this.destroyed = true; clearTimeout(this._camTimer);
      this._lensBound.forEach(([type, fn, opts]) => this.host.removeEventListener(type, fn, opts));
      this._lensBound = [];
      try { this.viewer.spin(false); this.viewer.clear(); } catch (e) { /* ignore */ }
      clear(this.host);
    }

    // ---- deep-zoom "detail lens": when the camera is close, show atoms/bonds and nearby residue labels
    _bindLens() {
      const schedule = Atlas.debounce(() => { if (!this.destroyed) this.updateLens(); }, 160);
      [["wheel", { passive: true }], ["pointerup", undefined], ["touchend", { passive: true }]].forEach(([type, opts]) => {
        this.host.addEventListener(type, schedule, opts);
        this._lensBound.push([type, schedule, opts]);
      });
    }
    updateLens() {
      if (this.destroyed || !this.lens || !this.model || this.small || !this.z0) return;
      let view; try { view = this.viewer.getView(); } catch (e) { return; }
      if (!Array.isArray(view) || view.length < 4 || !view.slice(0, 4).every(Number.isFinite)) return;
      const ratio = Math.abs(view[3]) / this.z0;
      const on = this.lensPinned || (ratio < 0.5 && this.atoms.length <= ATOM_LIMIT);
      if (on !== this.lensOn) { this.lensOn = on; this.applyStyle(); if (this.onLens) this.onLens(on); }
      this._clearLensLabels();
      if (on) {
        const centre = { x: -view[0], y: -view[1], z: -view[2] };
        const near = this.atoms.filter((a) => a.atom === "CA" && Math.hypot(a.x - centre.x, a.y - centre.y, a.z - centre.z) < 11).slice(0, 30);
        near.forEach((a) => this.lensLabels.push(this.viewer.addLabel(a.resn + a.resi, { position: { x: a.x, y: a.y, z: a.z }, fontSize: 11, fontColor: "#eaf1ec", backgroundColor: "#0c1210", backgroundOpacity: 0.6, showBackground: true })));
        this.viewer.render();
      }
    }
  }
  Atlas.MolViewer = MolViewer;

  // ------------------------------------------------------------------ workspace (full screen)
  // W.token identifies the newest open/load request; any older request that finishes later is ignored.
  // W.current is the "source|id" whose coordinates are on screen; W.timeouts are how long a request may take before it is reported as failed.
  const W = { viewer: null, mode: null, listing: null, protein: null, token: 0, info: {}, current: null,
              timeouts: { listing: 60000, file: 150000, predict: 180000 } };
  // Coordinates already downloaded this session (reused instead of downloading again) and requests currently in flight (shared, never duplicated).
  const store = { files: new Map(), listings: new Map(), pending: new Map(), limit: 4 };
  // The overlays (#structure-loading, #structure-empty) have no ".hidden" rule in style.css, so toggling the class alone never hid them.
  // The inline display makes visibility follow the real state whatever the stylesheet says; "" hands control back to the stylesheet when shown.
  const show = (id, on) => { const n = $(id); if (!n) return; n.classList.toggle("hidden", !on); n.style.display = on ? "" : "none"; };
  function setTitle(name, sub) { $("#structure-modal-name").textContent = name; $("#structure-modal-sub").textContent = sub || ""; }
  function setLoading(text) { show("#structure-loading", !!text); if (text) $("#structure-loading-text").textContent = text; }
  function remember(map, key, value) { map.delete(key); map.set(key, value); while (map.size > store.limit) map.delete(map.keys().next().value); }
  function withTimeout(promise, ms, message) {
    let timer;
    const limit = new Promise((_, reject) => { timer = setTimeout(() => reject(new Atlas.AtlasError("timeout", message)), ms); });
    return Promise.race([promise, limit]).finally(() => clearTimeout(timer));
  }
  function once(key, make) {            // one request per key at a time; everybody waiting on it gets the same answer
    if (store.pending.has(key)) return store.pending.get(key);
    const request = make().finally(() => store.pending.delete(key));
    store.pending.set(key, request);
    return request;
  }
  function getFile(source, ident) {
    const key = source + "|" + ident;
    if (store.files.has(key)) return Promise.resolve(store.files.get(key));
    const limit = source === "esmfold" ? W.timeouts.predict : W.timeouts.file;
    return once("file:" + key, () => {
      const request = api("/structures/" + encodeURIComponent(source) + "/" + encodeURIComponent(ident) + "/file");
      request.then((file) => remember(store.files, key, file), () => {});      // keep it even if the caller already gave up: "Try again" is then instant
      return withTimeout(request, limit, "The structure server did not answer in time. Try again in a moment.");
    });
  }
  function getListing(protein) {
    if (store.listings.has(protein.id)) return Promise.resolve(store.listings.get(protein.id));
    return once("listing:" + protein.id, async () => {
      const listing = await withTimeout(api("/proteins/" + encodeURIComponent(protein.id) + "/structures"), W.timeouts.listing,
        "The structure databases did not answer in time. Try again in a moment.");
      if (!(listing.sources || []).slice(0, 3).some((s) => s.status === "unavailable")) remember(store.listings, protein.id, listing);
      return listing;
    });
  }

  function seg(items, current, onChange) {
    const group = el("div", { class: "seg-group" });
    items.forEach(([value, label]) => {
      const b = el("button", { class: "seg" + (value === current ? " active" : ""), type: "button", onclick: () => {
        Array.from(group.children).forEach((c) => c.classList.toggle("active", c === b)); onChange(value); } }, label);
      group.appendChild(b);
    });
    return group;
  }
  function toggleBtn(label, initial, onChange, title) {
    const b = el("button", { class: "seg toggle" + (initial ? " active" : ""), type: "button", title: title || "" }, label);
    b.addEventListener("click", () => { const on = onChange(); b.classList.toggle("active", !!on); });
    return b;
  }
  function buildControls(small) {
    const bar = $("#structure-controls"); clear(bar);
    const V = W.viewer;
    const styles = small ? [["ballstick", "Ball & stick"], ["stick", "Stick"], ["sphere", "Spacefill"], ["line", "Line"]]
      : [["cartoon", "Cartoon"], ["stick", "Stick"], ["ballstick", "Ball & stick"], ["sphere", "Spacefill"], ["surface", "Surface"]];
    const colors = small ? [["element", "Element"], ["part", "Backbone / side chain"]]
      : [["spectrum", "N\u2192C"], ["chain", "Chain"], ["ss", "Structure"], ["confidence", "Confidence"]];
    bar.appendChild(seg(styles, V.style, (s) => V.setStyle(s)));
    bar.appendChild(seg(colors, V.color, (c) => V.setColor(c)));
    if (small) {
      bar.appendChild(toggleBtn("Atom labels", false, () => V.toggle("showLabels"), "Label every heavy atom"));
    } else {
      bar.appendChild(toggleBtn("All atoms", false, () => V.toggleAtoms(), "Show every atom and bond as sticks"));
      bar.appendChild(toggleBtn("Waters", false, () => V.toggleWater(), "Show crystallographic waters"));
      bar.appendChild(toggleBtn("Detail lens", true, () => V.toggleLens(), "Show atoms, bonds and residue labels close to the region of interest (also switches on automatically when you zoom in close)"));
    }
    bar.appendChild(toggleBtn("\u27f2 Spin", false, () => V.toggle("spinning"), "Auto-rotate"));
    bar.appendChild(el("button", { class: "icon-btn", type: "button", title: "Zoom in", onclick: () => V.zoomBy(1.5) }, "+"));
    bar.appendChild(el("button", { class: "icon-btn", type: "button", title: "Zoom out", onclick: () => V.zoomBy(1 / 1.5) }, "\u2212"));
    bar.appendChild(el("button", { class: "icon-btn", type: "button", title: "Reset camera", onclick: () => V.resetCamera() }, "\u2922"));
    bar.appendChild(el("button", { class: "icon-btn", type: "button", title: "Save image (PNG)", onclick: () => {
      const uri = V.png(); if (!uri) { Atlas.toast("Image export is not available in this browser.", "error"); return; }
      const a = el("a", { href: uri, download: (W.protein ? W.protein.id : "molecule") + ".png" }); document.body.appendChild(a); a.click(); a.remove(); } }, "\ud83d\udcf7"));
    if (!small) {
      const input = el("input", { type: "text", placeholder: "Residue #", style: "width:96px;padding:6px 10px;border-radius:999px", "aria-label": "Jump to residue number" });
      input.addEventListener("keydown", (e) => {
        if (e.key !== "Enter") return;
        const n = parseInt(input.value, 10);
        if (!Number.isFinite(n) || !V.gotoResidue(n)) Atlas.toast("Residue " + fmt.text(input.value) + " is not in this structure.", "error");
      });
      bar.appendChild(input);
    }
  }
  function inspectorEmpty(small) {
    const box = $("#structure-inspector"); clear(box);
    box.appendChild(el("h4", {}, small ? "Atom inspector" : "Residue inspector"));
    box.appendChild(el("p", { class: "muted" }, small ? "Hover an atom to see its element and bonds. Drag to rotate, scroll to zoom, right-drag to pan." : "Click any atom or residue in the model. Zoom in close and atoms, bonds and residue labels appear automatically."));
  }
  function inspectorAtom(atom) {
    const box = $("#structure-inspector"); clear(box);
    if (!atom) { inspectorEmpty(W.mode === "molecule"); return; }
    const rows = W.mode === "molecule" ? [
      ["Element", atom.elem], ["Atom #", atom.serial], ["Bonded to", (atom.bonds || []).length + " atoms"],
      ["x, y, z (\u00c5)", [atom.x, atom.y, atom.z].map((v) => fmt.n(v, 2)).join(", ")],
    ] : [
      ["Residue", atom.resn + " " + atom.resi], ["Chain", atom.chain], ["Atom", atom.atom + " (" + atom.elem + ")"],
      [W.info.kind === "predicted" && W.info.provider === "AlphaFold DB" ? "pLDDT confidence" : "B-factor", fmt.n(atom.b, 1)],
      ["x, y, z (\u00c5)", [atom.x, atom.y, atom.z].map((v) => fmt.n(v, 2)).join(", ")],
    ];
    box.appendChild(el("h4", {}, W.mode === "molecule" ? "Atom" : "Residue"));
    rows.forEach(([k, v]) => box.appendChild(el("div", { class: "insp-row" }, [el("span", {}, k), el("b", {}, fmt.text(v))])));
  }

  function ensureViewer(small) {
    const host = $("#structure-canvas");
    if (W.viewer) { W.viewer.destroy(); W.viewer = null; }
    clear(host);
    if (!mol3d()) throw new Atlas.AtlasError("database_unavailable", "The 3D viewer library (3Dmol.js) could not be downloaded. Check your internet connection or ad-blocker, then refresh with Ctrl+F5.");
    W.viewer = new MolViewer(host, { small });
    W.viewer.onPick = inspectorAtom;
    W.viewer.onHover = (atom) => { if (W.mode === "molecule") inspectorAtom(atom); };
    W.viewer.onLens = (on) => { const n = $("#structure-lens-note"); if (n) n.textContent = on ? "Detail lens active: atoms, bonds and nearby residue labels are shown." : ""; };
    return W.viewer;
  }

  function open() {
    show("#structure-modal", true); document.body.style.overflow = "hidden";
    show("#structure-empty", false);
  }
  function close() {
    W.token++;
    setLoading(null); show("#structure-empty", false);
    show("#structure-modal", false); document.body.style.overflow = "";
    if (W.viewer) { W.viewer.destroy(); W.viewer = null; }
    W.mode = null; W.protein = null; W.current = null; W.listing = null;
  }

  function sourceOptions(listing) {
    const out = [];
    (listing.sources || []).forEach((s) => {
      if (s.key === "esmfold") {
        if (s.status === "available") out.push({ value: "esmfold|" + s.entries[0].id, label: "\u2728 " + s.entries[0].label, group: "On-demand prediction" });
        return;
      }
      (s.entries || []).forEach((e) => out.push({ value: s.key + "|" + e.id, label: (e.label || e.id), group: s.label }));
    });
    return out;
  }
  function statusList(listing) {
    return el("ul", { class: "db-status" }, (listing.sources || []).map((s) => el("li", {}, [
      el("b", {}, s.label), " \u2014 ",
      el("span", { class: s.status === "found" ? "badge badge-ok" : s.status === "unavailable" ? "badge badge-err" : "badge" },
        s.status === "found" ? "found (" + fmt.int(s.count) + ")" : s.status === "unavailable" ? "Structure database unavailable" : s.status === "available" ? "can be predicted" : s.status === "too_long" ? "too long to predict" : "no entry"),
    ])));
  }
  function showNoStructure(listing) {
    setLoading(null);
    const box = $("#structure-empty"); clear(box);
    box.appendChild(el("h3", {}, "No 3D structure was located"));
    box.appendChild(el("p", {}, listing.summary || "No structure was found."));
    box.appendChild(statusList(listing));
    box.appendChild(el("p", { class: "muted" }, "\u201cNo AlphaFold model\u201d never means \u201cno structure exists\u201d: every source above was searched separately."));
    show("#structure-empty", true);
    $("#structure-source").classList.add("hidden");
  }

  function failBox(title, message, retry) {
    const box = $("#structure-empty"); clear(box);
    box.appendChild(el("h3", {}, title));
    box.appendChild(el("p", {}, message));
    if (retry) box.appendChild(el("button", { class: "btn", onclick: retry }, "Try again"));
    show("#structure-empty", true);
  }

  async function loadProteinSource(protein, value, token) {
    const [source, ident] = value.split("|");
    const predicting = source === "esmfold";
    show("#structure-empty", false);
    if (W.viewer && W.current === value && W.mode === "protein") { setLoading(null); return; }     // already on screen: nothing to fetch, nothing to redraw
    const cached = store.files.get(value);
    if (!cached) setLoading(predicting ? "Predicting the structure with ESMFold \u2014 this can take up to a minute\u2026" : "Downloading coordinates\u2026");
    else setLoading(null);
    let file = cached;
    if (!file) {
      try { file = await getFile(source, ident); }
      catch (error) {
        if (token !== W.token) return;
        setLoading(null);
        failBox(error.kind === "not_found" ? "That structure is no longer available"
          : error.kind === "timeout" ? "The structure database is not responding" : "Structure database unavailable",
          error.message, () => loadProteinSource(protein, value, ++W.token));
        return;
      }
      if (token !== W.token) return;
    }
    // From here on this request owns the workspace; whatever happens below, the loading overlay is cleared when it ends.
    try {
      const viewer = ensureViewer(false);
      const count = viewer.load(file.data, file.format);
      W.current = value;
      const meta = file.meta || {};
      W.info = { provider: meta.source, kind: source === "alphafold" || source === "esmfold" || source === "swissmodel" ? "predicted" : "experimental" };
      buildControls(false);
      // AlphaFold/ESM models store confidence in the B-factor column: default to confidence colouring there
      if (source === "alphafold" || source === "esmfold") { viewer.color = "confidence"; viewer.applyStyle(); buildControls(false); }
      inspectorEmpty(false);
      const sub = [meta.source, meta.id, count ? fmt.int(count) + " atoms" : null, meta.mean_plddt !== undefined && meta.mean_plddt !== null ? "mean pLDDT " + fmt.n(meta.mean_plddt, 1) : null].filter(Boolean).join(" \u00b7 ");
      setTitle(protein.id, sub);
      const note = $("#structure-lens-note");
      if (note) note.textContent = source === "esmfold" ? "Predicted on demand from the sequence: treat low-confidence regions with care." : "";
    } catch (error) {
      W.current = null;
      Atlas.diagnostics.errors.push("viewer: " + (error && error.message));
      failBox("3D viewer unavailable", (error && error.message) || "The structure could not be drawn.", () => loadProteinSource(protein, value, ++W.token));
    } finally {
      if (token === W.token) setLoading(null);
    }
  }

  async function openProtein(protein, srcValue, onSourceChange) {
    const token = ++W.token;
    // Same protein already open (e.g. the source was changed): keep the workspace and only swap the structure, no new listing request.
    const switching = W.mode === "protein" && W.protein && W.protein.id === protein.id && !!W.listing;
    W.mode = "protein"; W.protein = protein;
    open();
    if (!switching) { setTitle(protein.id, protein.name || ""); clear($("#structure-controls")); inspectorEmpty(false); }
    const select = $("#structure-source"); select.classList.remove("hidden");
    let listing = switching ? W.listing : store.listings.get(protein.id);
    if (!listing) {
      setLoading("Searching PDB, AlphaFold DB and SWISS-MODEL for a structure\u2026");
      try { listing = await getListing(protein); }
      catch (error) {
        if (token !== W.token) return;
        setLoading(null);
        failBox(error.kind === "network" ? "Server unreachable" : error.kind === "timeout" ? "The structure databases are not responding" : "Structure database unavailable",
          error.message, () => openProtein(protein, srcValue, onSourceChange));
        return;
      }
      if (token !== W.token) return;
    }
    W.listing = listing;
    const options = sourceOptions(listing);
    if (!options.length) { showNoStructure(listing); return; }
    clear(select);
    const groups = {};
    options.forEach((o) => {
      if (!groups[o.group]) { groups[o.group] = el("optgroup", { label: o.group }); select.appendChild(groups[o.group]); }
      groups[o.group].appendChild(el("option", { value: o.value }, o.label));
    });
    let choice = options.find((o) => o.value === srcValue);
    if (!choice && listing.recommended) choice = options.find((o) => o.value === listing.recommended.source + "|" + listing.recommended.id);
    choice = choice || options[0];
    select.value = choice.value;
    select.onchange = () => { if (onSourceChange) onSourceChange(select.value); else loadProteinSource(protein, select.value, ++W.token); };
    await loadProteinSource(protein, choice.value, token);
  }

  function openMolecule(spec) {
    const token = ++W.token;
    W.mode = "molecule"; W.protein = null; W.info = {}; W.current = null; W.listing = null;
    open(); setTitle(spec.title, spec.subtitle); $("#structure-source").classList.add("hidden"); setLoading(null);
    try {
      const viewer = ensureViewer(true);
      viewer.load(spec.data, spec.format);
      buildControls(true); inspectorEmpty(true);
    } catch (error) {
      const box = $("#structure-empty"); clear(box); box.appendChild(el("h3", {}, "3D viewer unavailable")); box.appendChild(el("p", {}, error.message)); show("#structure-empty", true);
    }
    return token;
  }

  function init(onClose) {
    const closeBtn = $("#structure-close");
    if (closeBtn) closeBtn.addEventListener("click", () => { if (onClose) onClose(); else close(); });
    global.addEventListener("keydown", (e) => { if (e.key === "Escape" && W.mode) { if (onClose) onClose(); else close(); } });
    global.addEventListener("resize", Atlas.debounce(() => { if (W.viewer) W.viewer.resize(); }, 150));
  }

  Atlas.structure = { init, openProtein, openMolecule, close, isOpen: () => !!W.mode, state: W, store };
})(typeof window !== "undefined" ? window : globalThis);
