// epykit report — interactivity (vanilla JS, no dependencies).
(function () {
  "use strict";

  // ---- theme toggle (persisted) ----
  var root = document.documentElement;
  var tbtn = document.getElementById("themeBtn");
  function setTheme(t) {
    root.setAttribute("data-theme", t);
    if (tbtn) tbtn.textContent = t === "dark" ? "☀ Light" : "\u{1F319} Dark";
    try { localStorage.setItem("epykit-theme", t); } catch (e) {}
  }
  var saved = "light";
  try { saved = localStorage.getItem("epykit-theme") || "light"; } catch (e) {}
  setTheme(saved);
  if (tbtn) {
    tbtn.addEventListener("click", function () {
      setTheme(root.getAttribute("data-theme") === "dark" ? "light" : "dark");
    });
  }

  // ---- scroll-spy: highlight the active sidebar entry ----
  var links = Array.prototype.slice.call(document.querySelectorAll("#toc a"));
  var map = {};
  links.forEach(function (a) {
    var id = a.getAttribute("href");
    if (id && id.charAt(0) === "#") map[id.slice(1)] = a;
  });
  if (window.IntersectionObserver) {
    var obs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) {
          links.forEach(function (l) { l.classList.remove("active"); });
          if (map[e.target.id]) map[e.target.id].classList.add("active");
        }
      });
    }, { rootMargin: "-20% 0px -70% 0px" });
    document.querySelectorAll("main section").forEach(function (s) { obs.observe(s); });
  }

  // ---- per-table search filter ----
  document.querySelectorAll(".search").forEach(function (inp) {
    inp.addEventListener("input", function () {
      var q = inp.value.toLowerCase();
      var t = document.getElementById(inp.getAttribute("data-for"));
      if (!t) return;
      t.querySelectorAll("tbody tr").forEach(function (r) {
        r.style.display = r.textContent.toLowerCase().indexOf(q) !== -1 ? "" : "none";
      });
    });
  });

  // ---- click-to-sort table headers (numeric-aware) ----
  document.querySelectorAll("table.df th").forEach(function (th) {
    th.addEventListener("click", function () {
      var table = th.closest("table");
      var tb = table.querySelector("tbody");
      if (!tb) return;
      var rows = Array.prototype.slice.call(tb.rows);
      var idx = Array.prototype.indexOf.call(th.parentNode.children, th);
      var asc = th.getAttribute("data-asc") !== "1";
      th.setAttribute("data-asc", asc ? "1" : "0");
      rows.sort(function (a, b) {
        var x = a.cells[idx] ? a.cells[idx].textContent.replace(/[,%×\s]/g, "") : "";
        var y = b.cells[idx] ? b.cells[idx].textContent.replace(/[,%×\s]/g, "") : "";
        var nx = parseFloat(x), ny = parseFloat(y);
        if (!isNaN(nx) && !isNaN(ny)) return asc ? nx - ny : ny - nx;
        return asc ? x.localeCompare(y) : y.localeCompare(x);
      });
      rows.forEach(function (r) { tb.appendChild(r); });
    });
  });

  // ---- client-side CSV download of a rendered table ----
  document.querySelectorAll("[data-csv]").forEach(function (b) {
    b.addEventListener("click", function () {
      var t = document.getElementById(b.getAttribute("data-csv"));
      if (!t) return;
      var lines = Array.prototype.slice.call(t.querySelectorAll("tr")).map(function (r) {
        return Array.prototype.slice.call(r.cells).map(function (c) {
          return '"' + c.textContent.trim().replace(/"/g, '""') + '"';
        }).join(",");
      });
      var blob = new Blob([lines.join("\n")], { type: "text/csv" });
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url; a.download = b.getAttribute("data-csv") + ".csv";
      a.click(); URL.revokeObjectURL(url);
    });
  });

  // ---- copy methods text ----
  var cm = document.getElementById("copyMethods");
  if (cm) {
    cm.addEventListener("click", function () {
      var el = document.getElementById("methodsText");
      if (!el || !navigator.clipboard) return;
      navigator.clipboard.writeText(el.innerText).then(function () {
        var prev = cm.textContent; cm.textContent = "✓ Copied";
        setTimeout(function () { cm.textContent = prev; }, 1500);
      });
    });
  }
})();
