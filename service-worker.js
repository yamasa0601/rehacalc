const CACHE_NAME = "rehacalc-v71";

const ASSETS = [
  "./",
  "./index.html",
  "./evaluation.html",
  "./stroke-record.html",
  "./stroke-evaluation.html",
  "./record-destination.html",
  "./manifest.json",
  "./service-worker.js",
  "./icon-192.png",
  "./icon-512.png",
  "./assets/evaluation-form.pdf",
  "./assets/hdsr-mmse-sheet.pdf",
  "./assets/pdf-lib.min.js",
  "./docs/rehacalc_usage_guide.html",
  "./docs/rehacalc_usage_guide.pdf",
  "./calc/index.html",
  "./bbs/index.html",
  "./fac.html",
  "./sixmwt.html",
  "./dgi.html",
  "./ges.html",
  "./mmse.html",
  "./fab.html",
  "./fma.html",
  "./tis.html",
  "./mal.html",
  "./sara.html",
  "./minibestest.html",
  "./karvonen.html",
  "./keyform.html",
  "./gait/index.html",
  "./gait/gaitanalyze.js",
  "./gait/style.css",
  "./emg/index.html",
  "./emg/app.js",
  "./emg/style.css"
];

self.addEventListener("install", (event) => {
  event.waitUntil(
    caches
      .open(CACHE_NAME)
      .then((cache) => cache.addAll(ASSETS))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener("activate", (event) => {
  event.waitUntil(
    caches
      .keys()
      .then((keys) =>
        Promise.all(keys.filter((key) => key !== CACHE_NAME).map((key) => caches.delete(key)))
      )
      .then(() => self.clients.claim())
  );
});

self.addEventListener("fetch", (event) => {
  const request = event.request;
  if (request.method !== "GET") return;

  if (request.mode === "navigate") {
    event.respondWith(
      fetch(request).catch(() => caches.match("./index.html"))
    );
    return;
  }

  event.respondWith(
    caches.match(request).then((cached) => {
      if (cached) return cached;
      return fetch(request).then((response) => {
        const copy = response.clone();
        caches.open(CACHE_NAME).then((cache) => cache.put(request, copy));
        return response;
      });
    })
  );
});
