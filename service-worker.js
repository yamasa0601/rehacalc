const CACHE_NAME = "rehacalc-v82";

const ASSETS = [
  "./",
  "./index.html",
  "./evaluation.html",
  "./orthopedic-record.html",
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
  "./assessment-context.js",
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

self.addEventListener("message", (event) => {
  if (event.data?.type === "SKIP_WAITING") self.skipWaiting();
});

self.addEventListener("fetch", (event) => {
  const request = event.request;
  if (request.method !== "GET") return;

  event.respondWith(
    fetch(request, { cache: "no-store" }).then((response) => {
      const url = new URL(request.url);
      if (url.origin === self.location.origin && response.ok) {
        const copy = response.clone();
        caches.open(CACHE_NAME).then((cache) => cache.put(request, copy));
      }
      return response;
    }).catch(() => {
      return caches.match(request).then((cached) => {
        if (cached) return cached;
        if (request.mode === "navigate") return caches.match("./index.html");
        throw new Error("offline");
      });
    })
  );
});
