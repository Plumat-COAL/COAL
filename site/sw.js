// Offline copy of this site's own files. WEB-3: new file. See overnight/integration/CHANGES.md
// Network first (so a new version shows at once), the cached copy when offline. Same-origin GET only: downloads
// from github.com are never touched.
"use strict";
// CHANGED(WEB-8): was coal-site-v1 with img/broken-gene.png; now v2 with the logo assets and the WebP screenshot.
var CACHE = "coal-site-v2";
var FILES = ["./", "index.html", "style.css", "main.js", "favicon.svg", "explainer.html", "img/broken-gene.webp",
  "assets/coal-mark.svg", "assets/coal-logo.png", "assets/coal-logo@2x.png", "assets/coal-logo-dark.png",
  "assets/coal-logo-dark@2x.png", "assets/favicon-32.png", "assets/apple-touch-icon.png"];

self.addEventListener("install", function (event) {
  event.waitUntil(caches.open(CACHE).then(function (c) { return c.addAll(FILES); }).then(function () {
    return self.skipWaiting();
  }));
});

self.addEventListener("activate", function (event) {
  event.waitUntil(caches.keys().then(function (keys) {
    return Promise.all(keys.filter(function (k) { return k !== CACHE; }).map(function (k) { return caches.delete(k); }));
  }).then(function () { return self.clients.claim(); }));
});

self.addEventListener("fetch", function (event) {
  var req = event.request;
  if (req.method !== "GET" || new URL(req.url).origin !== self.location.origin) return;
  event.respondWith(fetch(req).then(function (res) {
    if (res.ok) {
      var copy = res.clone();
      caches.open(CACHE).then(function (c) { c.put(req, copy); });
    }
    return res;
  }).catch(function () {
    return caches.match(req, { ignoreSearch: true }).then(function (hit) {
      return hit || (req.mode === "navigate" ? caches.match("./") : Response.error());
    });
  }));
});
