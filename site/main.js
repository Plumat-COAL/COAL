// COAL website: copy buttons, the suggested download, old links, and the offline cache.
// All three downloads show, in the same order; the card for the visitor's system only gets a "For this computer"
// mark. Links to the old one-page site (index.html#windows and so on) go to the page that now holds that section. The Copy button keeps its own label.
// No external requests. Without JavaScript everything still shows; only the Copy buttons do nothing.
(function () {
  "use strict";

  // Sections that moved off the home page, and where they are now.
  var MOVED = {
    install: "install.html", windows: "install.html", "smart-app-control": "install.html", macos: "install.html",
    linux: "install.html", cli: "install.html", downloads: "install.html", files: "install.html",
    verify: "install.html", example: "guide.html", misfit: "guide.html", limits: "guide.html", docs: "guide.html",
    pal2nal: "differences.html", speed: "differences.html", cite: "cite.html"
  };
  var id = location.hash.slice(1);
  if (id && MOVED[id] && !document.getElementById(id)) {
    location.replace(MOVED[id] + "#" + id);
    return;
  }

  function copyText(text, done) {
    if (navigator.clipboard && window.isSecureContext) {
      navigator.clipboard.writeText(text).then(function () { done(true); }, function () { done(fallback(text)); });
    } else {
      done(fallback(text));
    }
  }
  function fallback(text) {
    var ta = document.createElement("textarea");
    ta.value = text;
    ta.setAttribute("readonly", "");
    ta.className = "visually-hidden";
    document.body.appendChild(ta);
    ta.select();
    var ok = false;
    try { ok = document.execCommand("copy"); } catch (e) { ok = false; }
    document.body.removeChild(ta);
    return ok;
  }

  function suggest() {
    // Phones and tablets get no mark: COAL runs on a desktop or laptop.
    var ua = navigator.userAgent || "";
    if (/Android|iPhone|iPad|iPod|Mobile/i.test(ua)) return;
    var p = ((navigator.userAgentData && navigator.userAgentData.platform) || navigator.platform || ua).toLowerCase();
    var card = null;
    if (p.indexOf("win") === 0) card = document.getElementById("card-windows");
    else if (p.indexOf("mac") === 0 || p.indexOf("linux") !== -1) card = document.getElementById("card-pip");
    if (!card) return;
    card.classList.add("is-suggested");
    var mark = card.querySelector(".suggest");
    if (mark) mark.hidden = false;
  }

  function init() {
    suggest();
    var status = document.getElementById("copy-status");
    var buttons = document.querySelectorAll("button[data-copy]");
    for (var k = 0; k < buttons.length; k++) {
      buttons[k].addEventListener("click", function (e) {
        var btn = e.currentTarget;
        var label = btn.getAttribute("data-label") || btn.textContent;
        btn.setAttribute("data-label", label);
        var src = document.getElementById(btn.getAttribute("data-copy"));
        copyText(src.textContent, function (ok) {
          if (ok) {
            btn.textContent = "Copied";
            btn.setAttribute("data-copied", "");
            if (status) status.textContent = "Command copied. Paste it into a terminal.";
            setTimeout(function () { btn.textContent = label; btn.removeAttribute("data-copied"); }, 2500);
          } else {
            // Select the text so Ctrl+C / Cmd+C works.
            var range = document.createRange();
            range.selectNodeContents(src);
            var sel = window.getSelection();
            sel.removeAllRanges();
            sel.addRange(range);
            if (status) status.textContent = "Could not copy automatically. The command is selected: press Ctrl+C (Cmd+C on a Mac).";
          }
        });
      });
    }

    // Offline: keep a copy of this site's own files (same origin only) for the next visit.
    var nav = window.navigator;
    if ("serviceWorker" in nav && window.isSecureContext) {
      nav.serviceWorker.register("sw.js").catch(function () { /* offline cache is optional */ });
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
