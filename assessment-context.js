(() => {
  "use strict";

  const HIP_RECORDS_KEY = "rehacalc_evaluation_records_v1";
  const STROKE_RECORDS_KEY = "rehacalc_stroke_records_v1";
  const APP_BUILD = "20260701-tug-v80";
  const CONTEXT_KEY = "rehacalc_assessment_context_v1";
  const SNAPSHOT_PREFIX = "rehacalc_assessment_snapshot_v1";
  const TARGET_PREF_KEY = "rehacalc_assessment_target_v1";
  const TRANSFER_KEYS = new Set([
    "rehacalc_record_transfer_v1",
    "rehacalc_evaluation_draft_v1",
    "rehacalc_stroke_draft_v1",
    "rehacalc_keyform_transfer_v1"
  ]);

  const STAGES = [
    { id: "admission", label: "入棟時" },
    { id: "fac3", label: "FAC3" },
    { id: "month1", label: "1M" },
    { id: "month2", label: "2M" },
    { id: "month3", label: "3M" },
    { id: "discharge", label: "退院時" }
  ];

  const PAGE_CONFIGS = [
    { match: "fab.html", label: "FAB", target: "choice", reader: readFab, detailKey: "fab" },
    { match: "fac.html", label: "FAC", target: "choice", reader: readFac, detailKey: "fac" },
    { match: "dgi.html", label: "DGI", target: "choice", reader: readDgi, detailKey: "dgi" },
    { match: "ges.html", label: "GES", target: "choice", reader: readGes, detailKey: "ges" },
    { match: "mmse.html", label: "MMSE", target: "choice", reader: readMmse, detailKey: "mmse" },
    { match: "sixmwt.html", label: "6分間歩行テスト", target: "choice", reader: readSixMwt },
    { match: "fma.html", label: "FMA", target: "stroke", reader: readFma, detailKey: "fma" },
    { match: "tis.html", label: "TIS", target: "stroke", reader: readTis, detailKey: "tis" },
    { match: "mal.html", label: "MAL", target: "stroke", reader: readMal, detailKey: "mal" },
    { match: "sara.html", label: "SARA", target: "stroke", reader: readSara, detailKey: "sara" },
    { match: "minibestest.html", label: "Mini-BESTest", target: "stroke", reader: readMiniBestest, detailKey: "miniBestest" },
    { match: "bbs/index.html", label: "BBS", target: "bbs", reader: readBbs, customSnapshot: snapshotBbs, restoreCustomSnapshot: restoreBbs },
    { match: "calc/index.html", label: "HHD", target: "choice", reader: readCalc, detailKey: "hhd" },
    { match: "gait/index.html", label: "歩行加速度解析", target: "choice", reader: readGait }
  ];

  const scriptTarget = document.currentScript?.dataset?.recordTarget || "";
  const config = PAGE_CONFIGS.find((item) => normalizedPath().endsWith(item.match));
  if (!config) return;
  if (scriptTarget) config.target = scriptTarget;

  const rawSetItem = Storage.prototype.setItem;
  let panel;
  let statusEl;
  let lastSnapshotKey = "";
  let restoring = false;
  let saveTimer = null;
  let targetBeforeDatasetChange = "";
  let hasAssessmentInteraction = false;

  patchDraftTransfers();

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }

  function init() {
    refreshServiceWorkers();
    injectStyles();
    panel = buildPanel();
    statusEl = panel.querySelector("[data-rac-status]");
    insertPanel(panel);
    applyStoredContext();
    renderPatientOptions();
    lastSnapshotKey = snapshotKey();
    restoreSnapshot(lastSnapshotKey);
    bindPanelEvents();
    bindAssessmentEvents();
    setStatus("IDまたは氏名と時期を選ぶと、この画面の入力を自動保存します。", "waiting");
    renderPreviousComparison();
  }

  function refreshServiceWorkers() {
    if (!("serviceWorker" in navigator)) return;
    let reloaded = false;
    navigator.serviceWorker.addEventListener("controllerchange", () => {
      if (reloaded || sessionStorage.getItem(APP_BUILD)) return;
      reloaded = true;
      sessionStorage.setItem(APP_BUILD, "1");
      location.reload();
    });
    navigator.serviceWorker.getRegistrations().then((registrations) => {
      registrations.forEach((registration) => {
        registration.update();
        if (registration.waiting) {
          registration.waiting.postMessage({ type: "SKIP_WAITING" });
        }
        registration.addEventListener("updatefound", () => {
          const worker = registration.installing;
          if (!worker) return;
          worker.addEventListener("statechange", () => {
            if (worker.state === "installed" && navigator.serviceWorker.controller) {
              worker.postMessage({ type: "SKIP_WAITING" });
            }
          });
        });
      });
    }).catch(() => {});
  }

  function normalizedPath() {
    let path = location.pathname;
    if (path.endsWith("/")) path += "index.html";
    return path;
  }

  function today() {
    const date = new Date();
    date.setMinutes(date.getMinutes() - date.getTimezoneOffset());
    return date.toISOString().slice(0, 10);
  }

  function $(id) {
    return document.getElementById(id);
  }

  function clean(value) {
    return String(value ?? "").trim();
  }

  function numberText(value) {
    const text = clean(value).replaceAll(",", "");
    if (!text || text === "-") return null;
    const number = Number(text);
    return Number.isFinite(number) ? number : null;
  }

  function escapeHtml(value) {
    return clean(value)
      .replaceAll("&", "&amp;")
      .replaceAll("<", "&lt;")
      .replaceAll(">", "&gt;")
      .replaceAll('"', "&quot;");
  }

  function targetLabel(target) {
    return target === "stroke" ? "脳卒中・臨床用" : "大腿骨近位部骨折・研究用";
  }

  function stageLabel(stage) {
    return STAGES.find((item) => item.id === stage)?.label || clean(stage) || "-";
  }

  function buildPanel() {
    const root = document.createElement("section");
    root.id = "rehacalcAssessmentContext";
    root.className = "rac-panel";
    root.innerHTML = `
      <div class="rac-head">
        <div>
          <div class="rac-title">対象者・時期・保存</div>
          <div class="rac-sub">${escapeHtml(config.label)}の入力を、選択したIDと時期へ自動保存します。</div>
        </div>
        <button type="button" class="rac-save" data-rac-save>保存</button>
      </div>
      <div class="rac-grid">
        <label data-rac-target-wrap>
          <span>保存先</span>
          <select data-rac-target>
            <option value="hip">大腿骨近位部骨折・研究用</option>
            <option value="stroke">脳卒中・臨床用</option>
          </select>
        </label>
        <label>
          <span>保存済み対象者</span>
          <select data-rac-patient-select></select>
        </label>
        <label>
          <span>氏名</span>
          <input data-rac-patient-name autocomplete="off" placeholder="例：山田 太郎">
        </label>
        <label>
          <span>ID</span>
          <input data-rac-patient-id autocomplete="off" placeholder="例：12345">
        </label>
        <label>
          <span>評価日</span>
          <input data-rac-record-date type="date">
        </label>
        <label>
          <span>時期</span>
          <select data-rac-stage>
            ${STAGES.map((stage) => `<option value="${stage.id}">${stage.label}</option>`).join("")}
          </select>
        </label>
      </div>
      <div class="rac-status waiting" data-rac-status></div>
      <div class="rac-compare" data-rac-compare hidden>
        <div class="rac-compare-head">
          <div>
            <div class="rac-compare-title">前回下位項目</div>
            <div class="rac-compare-sub" data-rac-compare-sub></div>
          </div>
          <div class="rac-compare-actions">
            <button type="button" data-rac-apply-prev-all>全て反映</button>
            <button type="button" data-rac-copy-feedback>コピー</button>
          </div>
        </div>
        <div class="rac-compare-body" data-rac-compare-body></div>
      </div>
    `;
    configureTargetControl(root);
    return root;
  }

  function configureTargetControl(root) {
    const select = root.querySelector("[data-rac-target]");
    const wrap = root.querySelector("[data-rac-target-wrap]");
    if (config.target === "stroke" || config.target === "hip") {
      select.value = config.target;
      select.disabled = true;
      wrap.classList.add("rac-locked");
      return;
    }
    if (config.target === "bbs") {
      select.disabled = true;
      wrap.classList.add("rac-locked");
      return;
    }
    const stored = localStorage.getItem(TARGET_PREF_KEY);
    select.value = stored === "stroke" ? "stroke" : "hip";
  }

  function insertPanel(root) {
    const header = document.querySelector("body > header");
    if (header?.parentNode) {
      header.insertAdjacentElement("afterend", root);
      return;
    }
    document.body.insertBefore(root, document.body.firstElementChild);
  }

  function injectStyles() {
    if ($("rehacalcAssessmentContextStyles")) return;
    const style = document.createElement("style");
    style.id = "rehacalcAssessmentContextStyles";
    style.textContent = `
      .rac-panel {
        width: min(980px, calc(100% - 24px));
        margin: 12px auto;
        padding: 14px;
        border: 1px solid #d9e3e5;
        border-radius: 8px;
        background: #ffffff;
        color: #111827;
        box-shadow: 0 2px 10px rgba(15, 23, 42, .06);
        box-sizing: border-box;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Hiragino Sans", "Noto Sans JP", sans-serif;
      }
      .rac-head {
        display: flex;
        align-items: flex-start;
        justify-content: space-between;
        gap: 12px;
        margin-bottom: 12px;
      }
      .rac-title { font-size: 16px; font-weight: 900; line-height: 1.4; }
      .rac-sub { margin-top: 2px; color: #64748b; font-size: 12px; line-height: 1.55; }
      .rac-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 10px;
      }
      .rac-panel label { display: block; min-width: 0; }
      .rac-panel label span {
        display: block;
        margin: 0 0 4px;
        color: #475569;
        font-size: 12px;
        font-weight: 700;
      }
      .rac-panel input,
      .rac-panel select {
        width: 100%;
        min-height: 42px;
        padding: 9px 10px;
        border: 1px solid #cbd5e1;
        border-radius: 7px;
        background: #fff;
        color: #0f172a;
        font: inherit;
        font-size: 15px;
        box-sizing: border-box;
      }
      .rac-panel select:disabled { color: #475569; background: #f8fafc; }
      .rac-save {
        flex: 0 0 auto;
        min-width: 92px;
        min-height: 42px;
        padding: 10px 14px;
        border: 0;
        border-radius: 7px;
        background: #0f766e;
        color: #fff;
        font: inherit;
        font-size: 14px;
        font-weight: 900;
        cursor: pointer;
      }
      .rac-status {
        margin-top: 10px;
        padding: 9px 10px;
        border-radius: 7px;
        font-size: 12px;
        line-height: 1.55;
      }
      .rac-status.waiting { color: #475569; background: #f8fafc; border: 1px solid #e2e8f0; }
      .rac-status.saved { color: #0f5132; background: #edf8f1; border: 1px solid #b7e4ca; }
      .rac-status.warn { color: #7c2d12; background: #fff7ed; border: 1px solid #fed7aa; }
      .rac-compare {
        margin-top: 10px;
        padding: 10px;
        border: 1px solid #dbeafe;
        border-radius: 8px;
        background: #f8fbff;
      }
      .rac-compare-head {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 10px;
      }
      .rac-compare-title { font-size: 13px; font-weight: 900; }
      .rac-compare-sub { margin-top: 2px; color: #64748b; font-size: 12px; line-height: 1.45; }
      .rac-compare-actions { display: flex; gap: 6px; flex-wrap: wrap; justify-content: flex-end; }
      .rac-compare-actions button,
      .rac-prev-btn {
        border: 1px solid #bfdbfe;
        border-radius: 999px;
        background: #eff6ff;
        color: #1d4ed8;
        padding: 6px 9px;
        font: inherit;
        font-size: 12px;
        font-weight: 800;
        cursor: pointer;
      }
      .rac-compare-body {
        display: grid;
        gap: 6px;
        margin-top: 9px;
      }
      .rac-compare-item {
        display: grid;
        grid-template-columns: minmax(0, 1fr) auto auto;
        gap: 8px;
        align-items: center;
        padding: 7px;
        border: 1px solid #e2e8f0;
        border-radius: 7px;
        background: #fff;
        font-size: 12px;
      }
      .rac-compare-name {
        min-width: 0;
        color: #0f172a;
        font-weight: 800;
        line-height: 1.35;
      }
      .rac-delta {
        display: inline-block;
        min-width: 48px;
        padding: 3px 7px;
        border-radius: 999px;
        text-align: center;
        font-size: 11px;
        font-weight: 900;
        color: #334155;
        background: #f1f5f9;
      }
      .rac-delta.up { color: #166534; background: #dcfce7; }
      .rac-delta.down { color: #991b1b; background: #fee2e2; }
      @media (max-width: 760px) {
        .rac-head { flex-direction: column; }
        .rac-save { width: 100%; }
        .rac-grid { grid-template-columns: 1fr; }
        .rac-compare-head { flex-direction: column; }
        .rac-compare-actions { width: 100%; justify-content: flex-start; }
        .rac-compare-item { grid-template-columns: 1fr; }
      }
    `;
    document.head.appendChild(style);
  }

  function bindPanelEvents() {
    panel.querySelector("[data-rac-save]").addEventListener("click", () => {
      saveSnapshot(snapshotKey());
      saveAssessmentRecord({ manual: true });
      renderPreviousComparison();
    });

    panel.querySelector("[data-rac-apply-prev-all]")?.addEventListener("click", () => {
      applyPreviousItems();
    });

    panel.querySelector("[data-rac-copy-feedback]")?.addEventListener("click", () => {
      copyComparisonFeedback();
    });

    panel.querySelector("[data-rac-target]").addEventListener("change", () => {
      saveSnapshot(lastSnapshotKey);
      localStorage.setItem(TARGET_PREF_KEY, currentTarget());
      persistContext();
      renderPatientOptions();
      lastSnapshotKey = snapshotKey();
      restoreSnapshot(lastSnapshotKey);
      saveAssessmentRecord({ reason: "context" });
      renderPreviousComparison();
    });

    panel.querySelector("[data-rac-patient-select]").addEventListener("change", (event) => {
      saveSnapshot(lastSnapshotKey);
      applyPatientOption(event.target.selectedOptions[0]);
      persistContext();
      lastSnapshotKey = snapshotKey();
      restoreSnapshot(lastSnapshotKey);
      saveAssessmentRecord({ reason: "context" });
      renderPreviousComparison();
    });

    panel.querySelector("[data-rac-stage]").addEventListener("change", () => {
      saveSnapshot(lastSnapshotKey);
      persistContext();
      lastSnapshotKey = snapshotKey();
      restoreSnapshot(lastSnapshotKey);
      saveAssessmentRecord({ reason: "context" });
      renderPreviousComparison();
    });

    panel.querySelector("[data-rac-record-date]").addEventListener("change", () => {
      persistContext();
      saveSnapshot(snapshotKey());
      saveAssessmentRecord({ reason: "context" });
      renderPreviousComparison();
    });

    ["[data-rac-patient-name]", "[data-rac-patient-id]"].forEach((selector) => {
      panel.querySelector(selector).addEventListener("input", () => {
        persistContext();
        renderPatientOptions(false);
        lastSnapshotKey = snapshotKey();
        saveSnapshot(lastSnapshotKey);
        saveAssessmentRecord({ reason: "context" });
        renderPreviousComparison();
      });
    });
  }

  function bindAssessmentEvents() {
    targetBeforeDatasetChange = currentTarget();
    ["input", "change", "click"].forEach((eventName) => {
      document.addEventListener(eventName, (event) => {
        if (restoring) return;
        const target = event.target;
        if (!(target instanceof Element)) return;
        if (target.closest("#rehacalcAssessmentContext")) return;
        if (target.closest("a")) return;
        if (target.matches('input[type="file"]')) return;

        if (config.target === "bbs" && eventName === "change") {
          const newTarget = currentTarget();
          if (newTarget !== targetBeforeDatasetChange) {
            targetBeforeDatasetChange = newTarget;
            renderPatientOptions();
          }
        }

        if (!isMeasurementEventTarget(target)) return;
        hasAssessmentInteraction = true;

        clearTimeout(saveTimer);
        saveTimer = setTimeout(() => {
          saveSnapshot(snapshotKey());
          saveAssessmentRecord({ reason: "auto" });
          renderPreviousComparison();
        }, 250);
      }, true);
    });
  }

  function isMeasurementEventTarget(target) {
    if (target.closest("#items")) return true;
    if (target.matches("#dataset, #showZero, #prevTotal, #prevRecordSearch, #prevRecordSelect")) return false;
    if (target.matches("input, select, textarea")) return true;
    if (target.closest("button")) {
      const button = target.closest("button");
      const passiveButtons = new Set([
        "copy", "sendEval", "openKeyform", "csvExport", "reset", "expandAll",
        "send", "keyform", "g-copy", "h-bulk-open", "strengthModalClose",
        "bulk-clear", "bulk-send", "sendRecord", "createPdf", "clearScores"
      ]);
      return !passiveButtons.has(button.id);
    }
    return false;
  }

  function applyStoredContext() {
    let saved = {};
    try {
      saved = JSON.parse(localStorage.getItem(CONTEXT_KEY) || "{}");
    } catch {
      saved = {};
    }
    setPanelValue("[data-rac-patient-name]", saved.patientName || "");
    setPanelValue("[data-rac-patient-id]", saved.patientId || "");
    setPanelValue("[data-rac-record-date]", saved.recordDate || today());
    setPanelValue("[data-rac-stage]", saved.stage || "admission");
    if (config.target === "choice") {
      const target = saved.target === "stroke" || saved.target === "hip"
        ? saved.target
        : (localStorage.getItem(TARGET_PREF_KEY) || "hip");
      setPanelValue("[data-rac-target]", target);
    }
    if (config.target === "bbs") setPanelValue("[data-rac-target]", currentTarget());
  }

  function setPanelValue(selector, value) {
    const element = panel.querySelector(selector);
    if (element) element.value = value;
  }

  function persistContext() {
    const context = currentContext();
    rawSetItem.call(localStorage, CONTEXT_KEY, JSON.stringify(context));
    if (context.target === "hip" || context.target === "stroke") {
      rawSetItem.call(localStorage, TARGET_PREF_KEY, context.target);
    }
  }

  function currentTarget() {
    if (config.target === "stroke" || config.target === "hip") return config.target;
    if (config.target === "bbs") {
      const dataset = $("dataset")?.value;
      return dataset === "stroke" ? "stroke" : "hip";
    }
    return panel?.querySelector("[data-rac-target]")?.value === "stroke" ? "stroke" : "hip";
  }

  function currentContext() {
    const target = currentTarget();
    if (panel && config.target === "bbs") {
      const targetSelect = panel.querySelector("[data-rac-target]");
      if (targetSelect) targetSelect.value = target;
    }
    return {
      target,
      patientName: clean(panel?.querySelector("[data-rac-patient-name]")?.value),
      patientId: clean(panel?.querySelector("[data-rac-patient-id]")?.value),
      recordDate: panel?.querySelector("[data-rac-record-date]")?.value || today(),
      stage: panel?.querySelector("[data-rac-stage]")?.value || "admission"
    };
  }

  function contextPatientKey(context = currentContext()) {
    return clean(context.patientId) || clean(context.patientName);
  }

  function loadRecords(target = currentTarget()) {
    try {
      const records = JSON.parse(localStorage.getItem(target === "stroke" ? STROKE_RECORDS_KEY : HIP_RECORDS_KEY) || "[]");
      return Array.isArray(records) ? records : [];
    } catch {
      return [];
    }
  }

  function writeRecords(target, records) {
    rawSetItem.call(localStorage, target === "stroke" ? STROKE_RECORDS_KEY : HIP_RECORDS_KEY, JSON.stringify(records));
  }

  function patientKey(record) {
    return clean(record?.patientId) || clean(record?.patientName);
  }

  function patientLabel(record) {
    const name = clean(record?.patientName) || "氏名未入力";
    const id = clean(record?.patientId);
    return id ? `${name} (${id})` : name;
  }

  function renderPatientOptions(allowValueKeep = true) {
    const select = panel?.querySelector("[data-rac-patient-select]");
    if (!select) return;
    const selectedKey = allowValueKeep ? contextPatientKey() : "";
    const seen = new Map();
    loadRecords(currentTarget()).forEach((record) => {
      const key = patientKey(record);
      if (!key) return;
      const prev = seen.get(key);
      if (!prev || compareRecords(record, prev) > 0) seen.set(key, record);
    });
    select.innerHTML = "";
    const blank = document.createElement("option");
    blank.value = "";
    blank.textContent = seen.size ? "新規/手入力" : "保存済み対象者なし";
    select.appendChild(blank);
    [...seen.values()]
      .sort((a, b) => patientLabel(a).localeCompare(patientLabel(b), "ja"))
      .forEach((record) => {
        const option = document.createElement("option");
        option.value = patientKey(record);
        option.textContent = patientLabel(record);
        option.dataset.patientName = clean(record.patientName);
        option.dataset.patientId = clean(record.patientId);
        if (option.value === selectedKey) option.selected = true;
        select.appendChild(option);
      });
  }

  function applyPatientOption(option) {
    if (!option || !option.value) return;
    setPanelValue("[data-rac-patient-name]", option.dataset.patientName || "");
    setPanelValue("[data-rac-patient-id]", option.dataset.patientId || "");
  }

  function compareRecords(a, b) {
    const aKey = `${clean(a.recordDate)}|${clean(a.updatedAt || a.createdAt)}|${clean(a.id)}`;
    const bKey = `${clean(b.recordDate)}|${clean(b.updatedAt || b.createdAt)}|${clean(b.id)}`;
    return aKey.localeCompare(bKey);
  }

  function snapshotKey(context = currentContext()) {
    const patient = contextPatientKey(context) || "__no_patient__";
    return `${SNAPSHOT_PREFIX}:${normalizedPath()}:${context.target}:${patient}:${context.stage}`;
  }

  function getControls() {
    return [...document.querySelectorAll("input, select, textarea")].filter((element) => {
      if (element.closest("#rehacalcAssessmentContext")) return false;
      if (element.matches('input[type="button"], input[type="submit"], input[type="reset"], input[type="file"], input[type="hidden"]')) return false;
      if (element.matches("button")) return false;
      return true;
    });
  }

  function controlKey(element, index) {
    if (element.name) return `name:${element.name}`;
    if (element.id) return `id:${element.id}`;
    if (element.dataset?.field) return `data-field:${element.dataset.field}`;
    if (element.dataset?.key) return `data-key:${element.dataset.key}`;
    return `index:${index}`;
  }

  function saveSnapshot(key) {
    if (!key) return;
    const controls = {};
    getControls().forEach((element, index) => {
      const keyName = controlKey(element, index);
      if (element.type === "radio") {
        if (!controls[keyName]) controls[keyName] = { type: "radio", value: null };
        if (element.checked) controls[keyName].value = element.value;
        return;
      }
      if (element.type === "checkbox") {
        controls[keyName] = { type: "checkbox", checked: element.checked };
        return;
      }
      controls[keyName] = { type: element.tagName.toLowerCase(), value: element.value };
    });
    const snapshot = {
      savedAt: new Date().toISOString(),
      context: currentContext(),
      controls,
      custom: config.customSnapshot ? config.customSnapshot() : null
    };
    rawSetItem.call(localStorage, key, JSON.stringify(snapshot));
  }

  function restoreSnapshot(key) {
    let snapshot = null;
    try {
      snapshot = JSON.parse(localStorage.getItem(key) || "null");
    } catch {
      snapshot = null;
    }
    if (!snapshot) {
      setStatus(`${targetLabel(currentTarget())}：${stageLabel(currentContext().stage)}の下書きはまだありません。入力すると保存します。`, "waiting");
      return;
    }
    restoring = true;
    try {
      const controls = getControls();
      controls.forEach((element, index) => {
        const data = snapshot.controls?.[controlKey(element, index)];
        if (!data) return;
        if (element.type === "radio") {
          element.checked = data.value !== null && String(element.value) === String(data.value);
        } else if (element.type === "checkbox") {
          element.checked = !!data.checked;
        } else {
          element.value = data.value ?? "";
        }
        element.dispatchEvent(new Event("input", { bubbles: true }));
        element.dispatchEvent(new Event("change", { bubbles: true }));
      });
      if (snapshot.custom && config.restoreCustomSnapshot) {
        config.restoreCustomSnapshot(snapshot.custom);
      }
    } finally {
      restoring = false;
    }
    setStatus(`前回の入力下書きを読み込みました（${timeLabel(snapshot.savedAt)}）。`, "saved");
  }

  function timeLabel(iso) {
    const date = iso ? new Date(iso) : new Date();
    if (Number.isNaN(date.getTime())) return "";
    return date.toLocaleTimeString("ja-JP", { hour: "2-digit", minute: "2-digit" });
  }

  function saveAssessmentRecord(options = {}) {
    const payload = config.reader();
    const context = currentContext();
    if (options.reason === "context" && !hasAssessmentInteraction) return false;
    const patient = contextPatientKey(context);
    if (!patient) {
      setStatus("IDまたは氏名を入力すると、評価記録まで自動保存できます。", "waiting");
      return false;
    }
    if (!payload || !hasMeasurements(payload.measurements)) {
      if (options.manual) {
        setStatus("入力下書きを保存しました。評価記録に入れる点数はまだありません。", "warn");
      }
      return false;
    }
    savePayloadToRecords(payload, context);
    const mode = options.manual ? "保存しました" : "自動保存しました";
    setStatus(`${targetLabel(context.target)}：${patient} / ${stageLabel(context.stage)}へ${mode}。`, "saved");
    renderPatientOptions(true);
    return true;
  }

  function hasMeasurements(measurements) {
    return Object.values(measurements || {}).some((value) => clean(value) !== "");
  }

  function savePayloadToRecords(payload, context) {
    const records = loadRecords(context.target);
    const patient = contextPatientKey(context);
    const matching = records
      .filter((record) => patientKey(record) === patient && clean(record.stage) === clean(context.stage))
      .sort(compareRecords);
    const previous = matching.at(-1) || {};
    const record = {
      ...previous,
      id: previous.id || `as_${Date.now()}_${Math.random().toString(36).slice(2, 8)}`,
      patientName: context.patientName || previous.patientName || "",
      patientId: context.patientId || previous.patientId || "",
      recordDate: context.recordDate || previous.recordDate || today(),
      stage: context.stage,
      measurements: {
        ...(previous.measurements || {}),
        ...compactObject(payload.measurements)
      },
      assessmentDetails: {
        ...(previous.assessmentDetails || {}),
        ...(payload.assessmentDetails || {})
      },
      updatedAt: new Date().toISOString()
    };
    if (!record.createdAt) record.createdAt = previous.createdAt || new Date().toISOString();
    if (clean(payload.weight)) record.weight = clean(payload.weight);
    if (clean(payload.hhdArmCm)) {
      record.measurements.hhdArmCm = clean(payload.hhdArmCm);
    }
    const next = records.filter((item) => !(patientKey(item) === patient && clean(item.stage) === clean(context.stage)));
    next.push(record);
    writeRecords(context.target, next);
  }

  function compactObject(object) {
    const result = {};
    Object.entries(object || {}).forEach(([key, value]) => {
      if (clean(value) !== "") result[key] = clean(value);
    });
    return result;
  }

  function setStatus(message, state = "waiting") {
    if (!statusEl) return;
    statusEl.className = `rac-status ${state}`;
    statusEl.textContent = message;
  }

  function detailKey() {
    return config.detailKey || "";
  }

  function stageOrder(stage) {
    const index = STAGES.findIndex((item) => item.id === stage);
    return index >= 0 ? index : Number.MAX_SAFE_INTEGER;
  }

  function previousDetailRecord() {
    const key = detailKey();
    if (!key) return null;
    const context = currentContext();
    const patient = contextPatientKey(context);
    if (!patient) return null;
    const candidates = loadRecords(context.target)
      .filter((record) => patientKey(record) === patient)
      .filter((record) => clean(record.stage) !== clean(context.stage))
      .filter((record) => Array.isArray(record.assessmentDetails?.[key]?.items))
      .sort((a, b) => {
        const aEarlier = stageOrder(a.stage) < stageOrder(context.stage);
        const bEarlier = stageOrder(b.stage) < stageOrder(context.stage);
        if (aEarlier !== bEarlier) return aEarlier ? 1 : -1;
        const stageDiff = stageOrder(a.stage) - stageOrder(b.stage);
        if (stageDiff !== 0) return stageDiff;
        return compareRecords(a, b);
      });
    return candidates.at(-1) || null;
  }

  function currentAssessmentDetail() {
    const key = detailKey();
    if (!key) return null;
    try {
      const payload = config.reader();
      return payload?.assessmentDetails?.[key] || null;
    } catch {
      return null;
    }
  }

  function comparableItems() {
    const record = previousDetailRecord();
    if (!record) return null;
    const key = detailKey();
    const previousDetail = record.assessmentDetails?.[key];
    const currentDetail = currentAssessmentDetail();
    const currentMap = new Map();
    (currentDetail?.items || []).forEach((item, index) => {
      currentMap.set(itemIdentity(item, index), item);
    });
    return {
      record,
      previousDetail,
      currentDetail,
      items: (previousDetail.items || []).map((previous, index) => {
        const current = currentMap.get(itemIdentity(previous, index));
        const previousScore = numericScore(previous);
        const currentScore = current ? numericScore(current) : null;
        return {
          index,
          previous,
          current,
          previousScore,
          currentScore,
          delta: previousScore !== null && currentScore !== null ? currentScore - previousScore : null
        };
      }).filter((item) => item.previousScore !== null)
    };
  }

  function itemIdentity(item, index) {
    return clean(item.id) || clean(item.name) || `index:${index}`;
  }

  function numericScore(item) {
    const value = numberText(item?.score);
    return value === null ? null : value;
  }

  function itemLabel(item, index) {
    return clean(item?.name) || clean(item?.label) || `項目${index + 1}`;
  }

  function formatScore(value) {
    if (value === null || value === undefined) return "-";
    return Number.isInteger(Number(value)) ? String(Number(value)) : String(Number(value).toFixed(2)).replace(/\.?0+$/, "");
  }

  function renderPreviousComparison() {
    const box = panel?.querySelector("[data-rac-compare]");
    if (!box || !detailKey()) return;
    const body = box.querySelector("[data-rac-compare-body]");
    const sub = box.querySelector("[data-rac-compare-sub]");
    const comparison = comparableItems();
    if (!comparison?.items?.length) {
      box.hidden = true;
      return;
    }
    box.hidden = false;
    const label = `${patientLabel(comparison.record)} / ${clean(comparison.record.recordDate) || "-"} ${stageLabel(comparison.record.stage)}`;
    const improved = comparison.items.filter(isImprovedItem).length;
    const declined = comparison.items.filter(isDeclinedItem).length;
    sub.textContent = `前回：${label}。向上 ${improved}項目 / 低下 ${declined}項目`;
    body.innerHTML = comparison.items.map((item, arrayIndex) => {
      const deltaClass = isImprovedItem(item) ? "up" : isDeclinedItem(item) ? "down" : "";
      const unit = scoreUnit(item.previous);
      const deltaText = item.delta === null ? "未入力" : item.delta > 0 ? `+${formatScore(item.delta)}${unit}` : item.delta < 0 ? `${formatScore(item.delta)}${unit}` : "変化なし";
      const currentText = item.currentScore === null ? "今回 -" : `今回 ${formatScore(item.currentScore)}${unit}`;
      return `
        <div class="rac-compare-item">
          <div class="rac-compare-name">${escapeHtml(item.index + 1)}. ${escapeHtml(itemLabel(item.previous, item.index))}<br><span style="color:#64748b;font-weight:700;">前回 ${formatScore(item.previousScore)}${unit} / ${currentText}</span></div>
          <span class="rac-delta ${deltaClass}">${escapeHtml(deltaText)}</span>
          <button type="button" class="rac-prev-btn" data-rac-prev-index="${arrayIndex}">前回値を反映</button>
        </div>
      `;
    }).join("");
    body.querySelectorAll("[data-rac-prev-index]").forEach((button) => {
      button.addEventListener("click", () => {
        const index = Number(button.dataset.racPrevIndex);
        const item = comparison.items[index];
        if (!item) return;
        applyPreviousItem(item.previous, item.index);
      });
    });
  }

  function comparisonFeedbackText() {
    const comparison = comparableItems();
    if (!comparison?.items?.length) return "前回下位項目がありません。";
    const improved = comparison.items.filter(isImprovedItem);
    const declined = comparison.items.filter(isDeclinedItem);
    const unchangedMax = comparison.items.filter((item) => itemIsHigherBetter(item) && item.delta === 0 && numericScore(item.previous) === Number(item.previous?.max));
    const currentTotal = numberText(comparison.currentDetail?.total);
    const previousTotal = numberText(comparison.previousDetail?.total);
    const totalDelta = currentTotal !== null && previousTotal !== null ? currentTotal - previousTotal : null;
    const lines = [`${config.label} 下位項目フィードバック${Number.isFinite(totalDelta) ? `（前回比 ${totalDelta >= 0 ? "+" : ""}${formatScore(totalDelta)}点）` : ""}`];
    lines.push(`向上：${improved.length}項目`);
    improved.forEach((item) => lines.push(`・${item.index + 1}. ${itemLabel(item.previous, item.index)}：${formatScore(item.previousScore)}→${formatScore(item.currentScore)}${scoreUnit(item.previous)}（${formatDelta(item.delta)}）`));
    if (declined.length) {
      lines.push(`低下：${declined.length}項目`);
      declined.forEach((item) => lines.push(`・${item.index + 1}. ${itemLabel(item.previous, item.index)}：${formatScore(item.previousScore)}→${formatScore(item.currentScore)}${scoreUnit(item.previous)}（${formatDelta(item.delta)}）`));
    }
    if (unchangedMax.length) {
      lines.push(`満点維持：${unchangedMax.length}項目`);
      unchangedMax.forEach((item) => lines.push(`・${item.index + 1}. ${itemLabel(item.previous, item.index)}`));
    }
    return lines.join("\n");
  }

  function scoreUnit(item) {
    return clean(item?.unit) || "点";
  }

  function formatDelta(delta) {
    return delta > 0 ? `+${formatScore(delta)}` : formatScore(delta);
  }

  function itemIsHigherBetter(item) {
    return item.previous?.compare !== "neutral" && item.previous?.higherIsBetter !== false;
  }

  function isImprovedItem(item) {
    if (item.delta === null || item.previous?.compare === "neutral") return false;
    return item.previous?.higherIsBetter === false ? item.delta < 0 : item.delta > 0;
  }

  function isDeclinedItem(item) {
    if (item.delta === null || item.previous?.compare === "neutral") return false;
    return item.previous?.higherIsBetter === false ? item.delta > 0 : item.delta < 0;
  }

  function copyComparisonFeedback() {
    const text = comparisonFeedbackText();
    navigator.clipboard?.writeText(text)
      .then(() => setStatus("下位項目フィードバックをコピーしました。", "saved"))
      .catch(() => {
        const textarea = document.createElement("textarea");
        textarea.value = text;
        document.body.appendChild(textarea);
        textarea.select();
        document.execCommand("copy");
        textarea.remove();
        setStatus("下位項目フィードバックをコピーしました。", "saved");
      });
  }

  function applyPreviousItems() {
    const comparison = comparableItems();
    if (!comparison?.items?.length) return;
    comparison.items.forEach((item) => applyPreviousItem(item.previous, item.index, { silent: true }));
    afterApplyingPreviousItems("前回下位項目を全て反映しました。");
  }

  function applyPreviousItem(item, index, options = {}) {
    const applied = applyControlValue(item, index);
    if (!applied) return false;
    if (!options.silent) afterApplyingPreviousItems("前回値を反映しました。");
    return true;
  }

  function afterApplyingPreviousItems(message) {
    triggerCalculations();
    saveSnapshot(snapshotKey());
    saveAssessmentRecord({ reason: "auto" });
    renderPreviousComparison();
    setStatus(message, "saved");
  }

  function applyControlValue(item, index) {
    const control = item?.control || fallbackControl(item, index);
    const score = clean(item?.score);
    if (!control || score === "") return false;
    if (control.type === "radio") {
      const input = [...document.querySelectorAll('input[type="radio"]')]
        .find((candidate) => candidate.name === control.name && clean(candidate.value) === score);
      if (!input) return false;
      input.click();
      input.dispatchEvent(new Event("input", { bubbles: true }));
      input.dispatchEvent(new Event("change", { bubbles: true }));
      return true;
    }
    if (control.type === "select") {
      const select = control.selector
        ? document.querySelector(control.selector)
        : control.dataKey
          ? document.querySelector(`select[data-key="${cssAttr(control.dataKey)}"]`)
          : control.dataItem
            ? selectByDataItem(control)
            : null;
      if (!select) return false;
      select.value = score;
      select.dispatchEvent(new Event("input", { bubbles: true }));
      select.dispatchEvent(new Event("change", { bubbles: true }));
      return true;
    }
    if (control.type === "input") {
      const input = control.selector
        ? document.querySelector(control.selector)
        : control.dataField
          ? document.querySelector(`input[data-field="${cssAttr(control.dataField)}"]`)
          : control.dataArmField
            ? document.querySelector(`input[data-arm-field="${cssAttr(control.dataArmField)}"]`)
            : control.id
              ? $(control.id)
              : null;
      if (!input) return false;
      input.value = score;
      input.dispatchEvent(new Event("input", { bubbles: true }));
      input.dispatchEvent(new Event("change", { bubbles: true }));
      return true;
    }
    if (control.type === "checkbox") {
      const checkbox = document.querySelector(control.selector);
      if (!checkbox) return false;
      checkbox.checked = score === "1" || score === "true";
      checkbox.dispatchEvent(new Event("input", { bubbles: true }));
      checkbox.dispatchEvent(new Event("change", { bubbles: true }));
      return true;
    }
    return false;
  }

  function selectByDataItem(control) {
    const selector = `select[data-item="${cssAttr(control.dataItem)}"]`;
    const matches = [...document.querySelectorAll(selector)];
    if (control.side) return matches.find((select) => select.dataset.side === control.side) || null;
    return matches.find((select) => !select.dataset.side) || matches[0] || null;
  }

  function cssAttr(value) {
    return String(value ?? "").replaceAll("\\", "\\\\").replaceAll('"', '\\"');
  }

  function fallbackControl(item, index) {
    const key = detailKey();
    if (key === "fab") return { type: "radio", name: `q${index + 1}` };
    if (key === "fac") return { type: "radio", name: "facScore" };
    if (key === "dgi" || key === "ges") return { type: "radio", name: `q${index}` };
    if (key === "miniBestest") return { type: "radio", name: `item${index + 1}` };
    if (key === "mmse") return { type: "select", dataKey: item.id };
    if (key === "tis" || key === "fma") return { type: "radio", name: item.id };
    if (key === "sara") return { type: "select", dataItem: item.id, side: item.side || "" };
    if (key === "hhd") {
      if (clean(item.id).endsWith("ArmCm")) return { type: "input", dataArmField: item.id };
      return { type: "input", dataField: item.id };
    }
    return null;
  }

  function triggerCalculations() {
    ["calculate", "calc", "update"].forEach((name) => {
      if (typeof window[name] === "function") {
        try { window[name](); } catch {}
      }
    });
  }

  function patchDraftTransfers() {
    if (window.__rehacalcAssessmentContextPatched) return;
    window.__rehacalcAssessmentContextPatched = true;
    Storage.prototype.setItem = function patchedSetItem(key, value) {
      if (this === localStorage && TRANSFER_KEYS.has(String(key)) && typeof value === "string") {
        try {
          const draft = JSON.parse(value);
          const context = panel ? currentContext() : {};
          const decorated = decorateDraft(draft, context, String(key));
          if (panel && decorated?.measurements) {
            const target = key === "rehacalc_stroke_draft_v1"
              ? "stroke"
              : key === "rehacalc_evaluation_draft_v1"
                ? "hip"
                : context.target;
            if ((target === "hip" || target === "stroke") && contextPatientKey({ ...context, target })) {
              savePayloadToRecords(decorated, { ...context, target });
            }
          }
          return rawSetItem.call(this, key, JSON.stringify(decorated));
        } catch {
          return rawSetItem.call(this, key, value);
        }
      }
      return rawSetItem.call(this, key, value);
    };
  }

  function decorateDraft(draft, context, key) {
    const target = key === "rehacalc_stroke_draft_v1"
      ? "stroke"
      : key === "rehacalc_evaluation_draft_v1"
        ? "hip"
        : context.target;
    const next = {
      ...draft,
      preferredRecordTarget: target,
      recordDate: context.recordDate || draft.recordDate || today(),
      stage: context.stage || draft.stage || "admission"
    };
    if (context.patientName) next.patientName = context.patientName;
    if (context.patientId) next.patientId = context.patientId;
    try {
      const livePayload = config?.reader?.();
      if (livePayload) {
        next.measurements = {
          ...(livePayload.measurements || {}),
          ...(next.measurements || {})
        };
        next.assessmentDetails = {
          ...(livePayload.assessmentDetails || {}),
          ...(next.assessmentDetails || {})
        };
        if (!clean(next.weight) && clean(livePayload.weight)) next.weight = clean(livePayload.weight);
        if (!clean(next.hhdArmCm) && clean(livePayload.hhdArmCm)) next.hhdArmCm = clean(livePayload.hhdArmCm);
      }
    } catch {}
    return next;
  }

  function checkedValue(name) {
    return document.querySelector(`input[name="${name}"]:checked`)?.value ?? "";
  }

  function countChecked(names) {
    return names.filter((name) => clean(checkedValue(name)) !== "").length;
  }

  function radioInputs(name) {
    return [...document.querySelectorAll('input[type="radio"]')].filter((input) => input.name === name);
  }

  function maxRadioValue(name) {
    const values = radioInputs(name).map((input) => numberText(input.value)).filter((value) => value !== null);
    return values.length ? Math.max(...values) : null;
  }

  function radioGroupNames(scope = document) {
    const names = [];
    scope.querySelectorAll('input[type="radio"]').forEach((input) => {
      if (input.name && !names.includes(input.name)) names.push(input.name);
    });
    return names;
  }

  function radioItemLabel(name, index) {
    const input = radioInputs(name)[0];
    const item = input?.closest(".item, article, tr, section, .question");
    const label = item?.querySelector(".item-name, .itemname, .item-title, .title, h3, h2, td:nth-child(2)")?.textContent;
    return clean(label).replace(/\s+/g, " ") || `項目${index + 1}`;
  }

  function radioDetail(names, options = {}) {
    return names.map((name, index) => ({
      id: options.idMode === "name" ? name : index + 1,
      name: options.labels?.[index] || radioItemLabel(name, index),
      score: clean(checkedValue(name)),
      max: maxRadioValue(name),
      control: { type: "radio", name }
    }));
  }

  function selectMax(select) {
    const values = [...select.options].map((option) => numberText(option.value)).filter((value) => value !== null);
    return values.length ? Math.max(...values) : null;
  }

  function selectItemLabel(select, index) {
    const item = select.closest(".item, article, tr, section");
    const label = item?.querySelector(".item-title, .title, .item-name, td:nth-child(2), h3, h2")?.textContent;
    return clean(label).replace(/\s+/g, " ") || clean(select.getAttribute("aria-label")) || `項目${index + 1}`;
  }

  function readFab() {
    const names = ["q1", "q2", "q3", "q4", "q5", "q6"];
    if (!countChecked(names)) return null;
    return {
      sourceLabel: "FAB",
      measurements: { fab: $("total-score")?.textContent || "" },
      assessmentDetails: { fab: { total: numberText($("total-score")?.textContent), max: 18, items: radioDetail(names), updatedAt: new Date().toISOString() } }
    };
  }

  function readFac() {
    const score = checkedValue("facScore");
    if (!score) return null;
    return {
      sourceLabel: "FAC",
      measurements: { fac: score },
      assessmentDetails: {
        fac: {
          total: numberText(score),
          max: 5,
          items: [{
            id: "fac",
            name: "歩行自立度（FAC）",
            score: clean(score),
            max: 5,
            control: { type: "radio", name: "facScore" }
          }],
          updatedAt: new Date().toISOString()
        }
      }
    };
  }

  function readDgi() {
    const names = Array.from({ length: 8 }, (_, index) => `q${index}`);
    if (countChecked(names) < names.length) return null;
    return {
      sourceLabel: "DGI",
      measurements: { dgi: $("totalScore")?.textContent || "" },
      assessmentDetails: { dgi: { total: numberText($("totalScore")?.textContent), max: 24, items: radioDetail(names), updatedAt: new Date().toISOString() } }
    };
  }

  function readGes() {
    const names = Array.from({ length: 10 }, (_, index) => `q${index}`);
    if (!countChecked(names)) return null;
    return {
      sourceLabel: "GES",
      measurements: { mfes: $("total-score")?.textContent || "" },
      assessmentDetails: { ges: { total: numberText($("total-score")?.textContent), max: 100, items: radioDetail(names), updatedAt: new Date().toISOString() } }
    };
  }

  function readMmse() {
    const mmseSelects = [...document.querySelectorAll('.score-select[data-scale="mmse"]')];
    if (!mmseSelects.some((select) => clean(select.value) !== "")) return null;
    return {
      sourceLabel: "MMSE",
      measurements: { mmse: clean($("mmseTotal")?.textContent) },
      assessmentDetails: {
        mmse: {
          total: numberText($("mmseTotal")?.textContent),
          max: 30,
          items: mmseSelects.map((select, index) => ({
            id: select.dataset.key,
            name: selectItemLabel(select, index),
            score: clean(select.value),
            max: selectMax(select),
            control: { type: "select", dataKey: select.dataset.key }
          })),
          updatedAt: new Date().toISOString()
        }
      }
    };
  }

  function readSixMwt() {
    const distance = numberText($("totalDistance")?.textContent);
    if (!distance || distance <= 0) return null;
    return {
      sourceLabel: "6分間歩行テスト",
      measurements: { sixMinuteWalkDistance: String(distance).replace(/\.0$/, "") }
    };
  }

  function readFma() {
    const measurements = {};
    if (document.querySelector('#content-ue input[type="radio"]:checked')) {
      measurements.fmaUpper = clean($("ue-motor")?.textContent).split("/")[0];
    }
    if (document.querySelector('#content-le input[type="radio"]:checked')) {
      measurements.fmaLower = clean($("le-motor")?.textContent).split("/")[0];
    }
    if (!hasMeasurements(measurements)) return null;
    const items = radioDetail(radioGroupNames(document), { idMode: "name" });
    const upper = numberText(measurements.fmaUpper) ?? 0;
    const lower = numberText(measurements.fmaLower) ?? 0;
    return {
      sourceLabel: "FMA",
      measurements,
      assessmentDetails: {
        fma: {
          total: upper + lower,
          max: 100,
          items,
          updatedAt: new Date().toISOString()
        }
      }
    };
  }

  function readTis() {
    if (!document.querySelector('input[type="radio"]:checked')) return null;
    return {
      sourceLabel: "TIS",
      measurements: { tis: $("total-score")?.textContent || "" },
      assessmentDetails: {
        tis: {
          total: numberText($("total-score")?.textContent),
          max: 23,
          items: radioDetail(radioGroupNames(document), { idMode: "name" }),
          updatedAt: new Date().toISOString()
        }
      }
    };
  }

  function readMal() {
    const count = Number($("includedCount")?.textContent || 0);
    if (!count) return null;
    const items = [];
    document.querySelectorAll("#malBody tr").forEach((row, index) => {
      if (row.querySelector(".exclude")?.checked) return;
      const name = clean(row.children[1]?.textContent) || `項目${index + 1}`;
      const aou = row.querySelector(".aou");
      const qom = row.querySelector(".qom");
      items.push({
        id: `aou_${index + 1}`,
        name: `${name} AOU`,
        score: clean(aou?.value),
        max: 5,
        control: { type: "select", selector: `#malBody tr[data-row="${index}"] .aou` }
      });
      items.push({
        id: `qom_${index + 1}`,
        name: `${name} QOM`,
        score: clean(qom?.value),
        max: 5,
        control: { type: "select", selector: `#malBody tr[data-row="${index}"] .qom` }
      });
    });
    return {
      sourceLabel: "MAL",
      measurements: {
        malAou: $("aouAverage")?.textContent || "",
        malQom: $("qomAverage")?.textContent || ""
      },
      assessmentDetails: {
        mal: {
          total: null,
          max: null,
          items,
          updatedAt: new Date().toISOString()
        }
      }
    };
  }

  function readSara() {
    if (!$("ataxiaConfirmed")?.checked) return null;
    const selects = [...document.querySelectorAll("select[data-item]")];
    const sideLabel = { right: "右", left: "左" };
    const items = selects.map((select, index) => {
      const itemId = clean(select.dataset.item);
      const side = clean(select.dataset.side);
      return {
        id: side ? `${itemId}_${side}` : itemId,
        name: `${selectItemLabel(select, index)}${side ? ` ${sideLabel[side] || side}` : ""}`,
        score: clean(select.value),
        max: selectMax(select),
        side,
        higherIsBetter: false,
        control: { type: "select", dataItem: itemId, side }
      };
    });
    return {
      sourceLabel: "SARA",
      measurements: { sara: $("totalScore")?.textContent || "" },
      assessmentDetails: {
        sara: {
          total: numberText($("totalScore")?.textContent),
          max: 40,
          items,
          updatedAt: new Date().toISOString()
        }
      }
    };
  }

  function readMiniBestest() {
    if (typeof window.buildMiniDetail !== "function") return null;
    const detail = window.buildMiniDetail();
    if (detail.items.some((item) => item.score === null)) return null;
    detail.items = detail.items.map((item, index) => ({
      ...item,
      control: { type: "radio", name: `item${index + 1}` }
    }));
    return {
      sourceLabel: "Mini-BESTest",
      measurements: { miniBestest: String(detail.total) },
      assessmentDetails: { miniBestest: detail }
    };
  }

  function readBbs() {
    const total = clean($("totalScore")?.textContent);
    const staticScore = numberText($("staticScore")?.textContent) ?? 0;
    const dynamicScore = numberText($("dynamicScore")?.textContent) ?? 0;
    const totals = { total: numberText(total) ?? 0, staticScore, dynamicScore };
    const detail = typeof window.buildBbsDetail === "function"
      ? window.buildBbsDetail(totals)
      : buildBbsDetailFromDom(totals);
    return {
      sourceLabel: "BBS",
      measurements: { bbs: String(totals.total) },
      assessmentDetails: { bbs: detail }
    };
  }

  function buildBbsDetailFromDom(totals) {
    return {
      assessmentKey: "bbs",
      total: totals.total,
      max: 56,
      items: [...document.querySelectorAll("#items .item")].map((item, index) => ({
        id: index + 1,
        name: clean(item.querySelector(".itemname")?.textContent).replace(/^\d+\.\s*/, ""),
        score: numberText(item.querySelector(".btn.active")?.textContent),
        max: 4
      })),
      updatedAt: new Date().toISOString()
    };
  }

  function readCalc() {
    const measurements = {};
    const speed = clean($("res-mps")?.textContent);
    if (speed && speed !== "-") measurements.comfortable10mSpeed = speed;
    const steps = clean($("g-steps")?.value);
    if (steps) measurements.comfortable10mSteps = steps;
    const tug = clean($("res-tug")?.textContent);
    if (tug && tug !== "-") measurements.tugT = tug;

    const singleTarget = $("h-target")?.value;
    if (singleTarget) {
      const right = clean($("h-rn")?.value);
      const left = clean($("h-ln")?.value);
      if (right) measurements[`${singleTarget}Right`] = right;
      if (left) measurements[`${singleTarget}Left`] = left;
    }
    document.querySelectorAll(".strength-input").forEach((input) => {
      const value = clean(input.value);
      if (value && input.dataset.field) measurements[input.dataset.field] = value;
    });
    document.querySelectorAll(".strength-arm-input").forEach((input) => {
      const value = clean(input.value);
      if (value && input.dataset.armField) measurements[input.dataset.armField] = value;
    });
    const arm = clean($("h-arm")?.value);
    if (arm) measurements.hhdArmCm = arm;
    if (!hasMeasurements(measurements)) return null;
    const hhdKeys = Object.keys(measurements).filter((key) =>
      !["comfortable10mSpeed", "comfortable10mSteps", "tugT"].includes(key)
    );
    const sourceParts = [];
    if (speed && speed !== "-") sourceParts.push("10m歩行");
    if (tug && tug !== "-") sourceParts.push("TUG");
    if (hhdKeys.length) sourceParts.push("HHD");
    const sourceLabel = sourceParts.join("・") || "HHD";
    const hhdItems = hhdDetailItems();
    return {
      sourceLabel,
      weight: $("h-bw")?.value || "",
      hhdArmCm: arm,
      measurements,
      assessmentDetails: hhdItems.length ? {
        hhd: {
          total: null,
          max: null,
          weight: $("h-bw")?.value || "",
          items: hhdItems,
          updatedAt: new Date().toISOString()
        }
      } : {}
    };
  }

  function hhdDetailItems() {
    const labels = {
      hipFlexor: "股屈筋",
      kneeExtensor: "膝伸展",
      ankleDorsiflexor: "足背屈",
      hipExtensor: "股伸展",
      hipAbductor: "股外転",
      anklePlantarflexor: "足底屈"
    };
    const sideLabels = { Right: "右HHD", Left: "左HHD" };
    const items = [];
    document.querySelectorAll(".strength-arm-input").forEach((input) => {
      const value = clean(input.value);
      if (!value) return;
      const muscle = clean(input.dataset.muscle);
      const field = clean(input.dataset.armField);
      items.push({
        id: field,
        name: `${labels[muscle] || muscle} アーム`,
        score: value,
        unit: "cm",
        compare: "neutral",
        control: { type: "input", dataArmField: field }
      });
    });
    document.querySelectorAll(".strength-input").forEach((input) => {
      const value = clean(input.value);
      if (!value) return;
      const field = clean(input.dataset.field);
      const side = field.endsWith("Right") ? "Right" : field.endsWith("Left") ? "Left" : "";
      const muscle = side ? field.slice(0, -side.length) : field;
      items.push({
        id: field,
        name: `${labels[muscle] || muscle} ${sideLabels[side] || ""}`.trim(),
        score: value,
        unit: "N",
        control: { type: "input", dataField: field }
      });
    });
    return items;
  }

  function readGait() {
    const metrics = metricMapFromDom();
    const strideCv = metrics.get("ストライド変動係数");
    if (!strideCv) return null;
    return {
      sourceLabel: "歩行加速度解析",
      measurements: { strideCv }
    };
  }

  function metricMapFromDom() {
    const map = new Map();
    document.querySelectorAll("#metrics .metric").forEach((item) => {
      const label = clean(item.querySelector(".k")?.textContent);
      const value = clean(item.querySelector(".v")?.textContent).replace(/[^\d.+-]/g, "");
      if (label && value) map.set(label, value);
    });
    return map;
  }

  function snapshotBbs() {
    return {
      scores: [...document.querySelectorAll("#items .item")].map((item) => clean(item.querySelector(".btn.active")?.textContent))
    };
  }

  function restoreBbs(custom) {
    if (!Array.isArray(custom.scores)) return;
    custom.scores.forEach((score, index) => {
      const item = document.querySelectorAll("#items .item")[index];
      const button = [...(item?.querySelectorAll(".btn") || [])].find((candidate) => clean(candidate.textContent) === clean(score));
      if (button && !button.classList.contains("active")) button.click();
    });
  }
})();
