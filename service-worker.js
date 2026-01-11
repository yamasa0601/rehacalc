const CACHE_NAME = "rehacalc-v2";

// 初期キャッシュ（必要なものだけ）
// ※ emg をオフラインでも開きたいので追加
const ASSETS = [
  "./",
  "./index.html",
  "./manifest.json",
  "./service-worker.js",
  "./icon-192.png",
  "./icon-512.png",

  "./emg/index.html",
  "./emg/app.js",
  "./emg/style.css",
];

self.addEventListener("install", (event) => {
  event.waitUntil(
    caches.open(CACHE_NAME)
      .then((cache) => cache.addAll(ASSETS))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener("activate", (event) => {
  event.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(keys.filter((k) => k !== CACHE_NAME).map((k) => caches.delete(k))))
      .then(() => self.clients.claim())
  );
});

self.addEventListener("fetch", (event) => {
  const req = event.request;
  const url = new URL(req.url);

  // 画面遷移は network-first（古い404/古いindexに引っ張られない）
  if (req.mode === "navigate") {
    event.respondWith(
      fetch(req).catch(async () => {
        // オフライン時のフォールバック
        if (url.pathname.includes("/rehacalc/emg/")) {
          return caches.match("./emg/index.html");
        }
        return caches.match("./index.html");
      })
    );
    return;
  }

  // それ以外は cache-first（静的ファイル向け）
  event.respondWith(
    caches.match(req).then((cached) => {
      if (cached) return cached;

      return fetch(req).then((resp) => {
        // ★重要：成功レスポンスだけキャッシュ（404等はキャッシュしない）
        if (
          url.origin === self.location.origin &&
          req.method === "GET" &&
          resp.ok
        ) {
          const copy = resp.clone();
          caches.open(CACHE_NAME).then((cache) => cache.put(req, copy));
        }
        return resp;
      });
    })
  );
});

