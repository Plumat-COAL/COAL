// COAL download site: copy buttons and the offline cache.
// WEB-3: new file. See overnight/integration/CHANGES.md
// CHANGED(WEB-11): the visitor's system is no longer detected; the page shows all three downloads at once.
// No external requests. Without JavaScript everything still shows; only the Copy buttons do nothing.
(function () {
  "use strict";

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

  function init() {
    var status = document.getElementById("copy-status");
    var buttons = document.querySelectorAll("button[data-copy]");
    for (var k = 0; k < buttons.length; k++) {
      buttons[k].addEventListener("click", function (e) {
        var btn = e.currentTarget;
        var src = document.getElementById(btn.getAttribute("data-copy"));
        copyText(src.textContent, function (ok) {
          if (ok) {
            btn.textContent = "Copied";
            btn.setAttribute("data-copied", "");
            if (status) status.textContent = "Command copied. Paste it into a terminal.";
            setTimeout(function () { btn.textContent = "Copy"; btn.removeAttribute("data-copied"); }, 2500);
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
