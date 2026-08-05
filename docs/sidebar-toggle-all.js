(function () {
  // Adds an "Expand all / Collapse all" control to the top of the Quarto
  // sidebar and remembers each section's open/closed state across page loads.
  //
  // Quarto renders a fresh, fully-expanded sidebar on every page, so without
  // persistence a "Collapse all" would be undone the moment you click a link.
  // We store per-section state (keyed by the section's label) in localStorage
  // and re-apply it on load.
  //
  // Wired in via each section's _metadata.yml (header-includes) rather than
  // _quarto.yml: the include kept getting stripped out of _quarto.yml by an
  // external formatter. The src is root-relative ("/sidebar-toggle-all.js");
  // Quarto rewrites it to the correct depth per page.
  var STORE_KEY = "nexus.sidebar.sections";

  function init() {
    var nav = document.getElementById("quarto-sidebar");
    if (!nav) return; // pages without a sidebar (Home, FAQ)

    var sections = Array.prototype.slice.call(
      nav.querySelectorAll("ul.sidebar-section.collapse")
    );
    if (sections.length < 2) return; // nothing worth a bulk toggle

    function sectionLabel(ul) {
      var li = ul.closest("li.sidebar-item-section");
      if (!li) return null;
      var el = li.querySelector(".sidebar-item-container .menu-text");
      return el ? el.textContent.trim() : null;
    }
    function toggleFor(ul) {
      return nav.querySelector(
        'a.sidebar-item-toggle[data-bs-target="#' + ul.id + '"]'
      );
    }
    function setExpanded(ul, open) {
      ul.classList.toggle("show", open);
      // The title link and the chevron both control the section; keep both in
      // sync so the chevron icon (driven by aria-expanded) points correctly.
      nav
        .querySelectorAll('[data-bs-target="#' + ul.id + '"]')
        .forEach(function (a) {
          a.setAttribute("aria-expanded", open ? "true" : "false");
        });
    }
    function readStore() {
      try {
        return JSON.parse(localStorage.getItem(STORE_KEY)) || {};
      } catch (e) {
        return {};
      }
    }
    function writeStore() {
      var state = {};
      sections.forEach(function (ul) {
        var label = sectionLabel(ul);
        if (label) state[label] = ul.classList.contains("show") ? "open" : "closed";
      });
      try {
        localStorage.setItem(STORE_KEY, JSON.stringify(state));
      } catch (e) {}
    }

    // 1. Re-apply any previously saved state, instantly (no animation flash).
    var saved = readStore();
    sections.forEach(function (ul) {
      var label = sectionLabel(ul);
      if (label && label in saved) setExpanded(ul, saved[label] === "open");
    });

    // 2. Remember the state whenever a single section is toggled. Bootstrap's
    //    collapse events bubble up to the nav.
    function onToggle(e) {
      var ul = e.target;
      if (ul && ul.classList && ul.classList.contains("sidebar-section")) writeStore();
    }
    nav.addEventListener("shown.bs.collapse", onToggle, true);
    nav.addEventListener("hidden.bs.collapse", onToggle, true);

    // 3. Bulk control. Click each section's toggle so the change routes through
    //    Bootstrap (animation + aria), then persist the intended uniform state
    //    directly (don't depend on transition timing).
    function setAll(open) {
      sections.forEach(function (ul) {
        if (ul.classList.contains("show") === open) return;
        var t = toggleFor(ul);
        if (t) t.click();
      });
      var state = {};
      sections.forEach(function (ul) {
        var label = sectionLabel(ul);
        if (label) state[label] = open ? "open" : "closed";
      });
      try {
        localStorage.setItem(STORE_KEY, JSON.stringify(state));
      } catch (e) {}
    }

    var bar = document.createElement("div");
    bar.className = "sidebar-toggle-all";

    var expand = document.createElement("button");
    expand.type = "button";
    expand.textContent = "Expand all";
    expand.setAttribute("aria-label", "Expand all sidebar sections");
    expand.addEventListener("click", function () {
      setAll(true);
    });

    var collapse = document.createElement("button");
    collapse.type = "button";
    collapse.textContent = "Collapse all";
    collapse.setAttribute("aria-label", "Collapse all sidebar sections");
    collapse.addEventListener("click", function () {
      setAll(false);
    });

    var sep = document.createElement("span");
    sep.className = "sep";
    sep.setAttribute("aria-hidden", "true");
    sep.textContent = "·";

    bar.appendChild(expand);
    bar.appendChild(sep);
    bar.appendChild(collapse);

    var container = nav.querySelector(".sidebar-menu-container") || nav;
    var list = container.querySelector("ul.list-unstyled") || container.querySelector("ul");
    if (list) container.insertBefore(bar, list);
    else container.insertBefore(bar, container.firstChild);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
