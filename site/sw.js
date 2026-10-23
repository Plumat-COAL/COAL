// Offline copy of this site's own files.
// Network first (so a new version shows at once), the cached copy when offline. Same-origin GET only: downloads
// from github.com are never touched.
"use strict";
// Bump the cache name whenever FILES or their content changes.
var CACHE = "coal-site-v8";
var FILES = ["./", "index.html", "install.html", "guide.html", "differences.html", "cite.html", "style.css", "main.js",
  "favicon.svg", "explainer.html", "img/broken-gene.webp", "assets/fonts/Literata-Regular.woff2", "assets/fonts/Literata-Semibold.woff2",
  "assets/fonts/Literata-Italic.woff2",
  "assets/coal-mark-dark.svg", "assets/favicon-32.png", "assets/apple-touch-icon.png"];

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
