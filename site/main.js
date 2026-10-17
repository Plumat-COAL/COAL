// COAL download site: pick the main download for the visitor's system, copy buttons, offline cache.
// WEB-3: new file. See overnight/integration/CHANGES.md
// CHANGED(WEB-8): was "COAL runs on Windows, macOS and Linux computers." for other systems; now says a computer is needed.
// No external requests. Without JavaScript the page shows both options (Windows installer and pip).
(function () {
  "use strict";

  var NAMES = { windows: "Windows", macos: "macOS", linux: "Linux" };

  // Exposed for tests: detectOS(userAgent, platform) -> "windows" | "macos" | "linux" | "other"
  function detectOS(ua, platform, touchPoints) {
    ua = ua || "";
    platform = platform || "";
    if (/Windows NT|Win64|WOW64/i.test(ua)) return "windows";
    if (/Android|iPhone|iPad|iPod/i.test(ua)) return "other";
    if (/Macintosh|Mac OS X/i.test(ua)) return touchPoints > 1 ? "other" : "macos";  // iPadOS asks for desktop sites
    if (/CrOS/i.test(ua)) return "other";
    if (/Linux|X11|FreeBSD/i.test(ua)) return "linux";
    if (/^Win/i.test(platform)) return "windows";
    if (/^Mac/i.test(platform)) return "macos";
    if (/Linux/i.test(platform)) return "linux";
    return "other";
  }
  window.coalDetectOS = detectOS;

  function setOS(os, detected) {
    var get = document.getElementById("get");
    if (!get) return;
    get.setAttribute("data-os", os);
    var line = document.getElementById("os-line");
    var label = document.getElementById("os-detected");
    line.hidden = false;
    if (NAMES[os]) {
      label.textContent = (detected ? "For " + NAMES[os] + " (detected)." : "For " + NAMES[os] + ".");
    } else {
      label.textContent = "COAL needs a Windows, macOS or Linux computer.";
    }
    var picks = document.querySelectorAll("[data-os-pick]");
    for (var i = 0; i < picks.length; i++) {
      picks[i].setAttribute("aria-pressed", picks[i].getAttribute("data-os-pick") === os ? "true" : "false");
    }
    var yours = document.querySelectorAll("[data-yours]");
    for (var j = 0; j < yours.length; j++) {
      yours[j].hidden = !(detected && yours[j].getAttribute("data-yours") === os);
    }
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

  function init() {
    var nav = window.navigator;
    var platform = (nav.userAgentData && nav.userAgentData.platform) || nav.platform || "";
    var os = detectOS(nav.userAgent, platform, nav.maxTouchPoints || 0);
    setOS(os, os !== "other");

    var picks = document.querySelectorAll("[data-os-pick]");
    for (var i = 0; i < picks.length; i++) {
      picks[i].addEventListener("click", function (e) { setOS(e.currentTarget.getAttribute("data-os-pick"), false); });
    }

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
            if (status) status.textContent = "Couldn't copy automatically: the command is selected, press Ctrl+C (Cmd+C on a Mac).";
          }
        });
      });
    }

    // Offline: keep a copy of this site's own files (same origin only) for the next visit.
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
