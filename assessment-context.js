(() => {
  "use strict";

  const HIP_RECORDS_KEY = "rehacalc_evaluation_records_v1";
  const STROKE_RECORDS_KEY = "rehacalc_stroke_records_v1";
  const APP_BUILD = "20260627-hhd-v76";
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
    { match: "fab.html", label: "FAB", target: "choice", reader: readFab },
    { match: "fac.html", label: "FAC", target: "choice", reader: readFac },
    { match: "dgi.html", label: "DGI", target: "choice", reader: readDgi },
    { match: "ges.html", label: "GES", target: "choice", reader: readGes },
    { match: "mmse.html", label: "MMSE", target: "choice", reader: readMmse },
    { match: "sixmwt.html", label: "6分間歩行テスト", target: "choice", reader: readSixMwt },
    { match: "fma.html", label: "FMA", target: "stroke", reader: readFma },
    { match: "tis.html", label: "TIS", target: "stroke", reader: readTis },
    { match: "mal.html", label: "MAL", target: "stroke", reader: readMal },
    { match: "sara.html", label: "SARA", target: "stroke", reader: readSara },
    { match: "minibestest.html", label: "Mini-BESTest", target: "stroke", reader: readMiniBestest },
    { match: "bbs/index.html", label: "BBS", target: "bbs", reader: readBbs, customSnapshot: snapshotBbs, restoreCustomSnapshot: restoreBbs },
    { match: "calc/index.html", label: "HHD", target: "choice", reader: readCalc },
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
      @media (max-width: 760px) {
        .rac-head { flex-direction: column; }
        .rac-save { width: 100%; }
        .rac-grid { grid-template-columns: 1fr; }
      }
    `;
    document.head.appendChild(style);
  }

  function bindPanelEvents() {
    panel.querySelector("[data-rac-save]").addEventListener("click", () => {
      saveSnapshot(snapshotKey());
      saveAssessmentRecord({ manual: true });
    });

    panel.querySelector("[data-rac-target]").addEventListener("change", () => {
      saveSnapshot(lastSnapshotKey);
      localStorage.setItem(TARGET_PREF_KEY, currentTarget());
      persistContext();
      renderPatientOptions();
      lastSnapshotKey = snapshotKey();
      restoreSnapshot(lastSnapshotKey);
      saveAssessmentRecord({ reason: "context" });
    });

    panel.querySelector("[data-rac-patient-select]").addEventListener("change", (event) => {
      saveSnapshot(lastSnapshotKey);
      applyPatientOption(event.target.selectedOptions[0]);
      persistContext();
      lastSnapshotKey = snapshotKey();
      restoreSnapshot(lastSnapshotKey);
      saveAssessmentRecord({ reason: "context" });
    });

    panel.querySelector("[data-rac-stage]").addEventListener("change", () => {
      saveSnapshot(lastSnapshotKey);
      persistContext();
      lastSnapshotKey = snapshotKey();
      restoreSnapshot(lastSnapshotKey);
      saveAssessmentRecord({ reason: "context" });
    });

    panel.querySelector("[data-rac-record-date]").addEventListener("change", () => {
      persistContext();
      saveSnapshot(snapshotKey());
      saveAssessmentRecord({ reason: "context" });
    });

    ["[data-rac-patient-name]", "[data-rac-patient-id]"].forEach((selector) => {
      panel.querySelector(selector).addEventListener("input", () => {
        persistContext();
        renderPatientOptions(false);
        lastSnapshotKey = snapshotKey();
        saveSnapshot(lastSnapshotKey);
        saveAssessmentRecord({ reason: "context" });
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
    return next;
  }

  function checkedValue(name) {
    return document.querySelector(`input[name="${name}"]:checked`)?.value ?? "";
  }

  function countChecked(names) {
    return names.filter((name) => clean(checkedValue(name)) !== "").length;
  }

  function radioDetail(names) {
    return names.map((name, index) => ({ id: index + 1, score: clean(checkedValue(name)) }));
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
      measurements: { fac: score }
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
          items: mmseSelects.map((select) => ({ id: select.dataset.key, score: clean(select.value) })),
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
    return { sourceLabel: "FMA", measurements };
  }

  function readTis() {
    if (!document.querySelector('input[type="radio"]:checked')) return null;
    return {
      sourceLabel: "TIS",
      measurements: { tis: $("total-score")?.textContent || "" }
    };
  }

  function readMal() {
    const count = Number($("includedCount")?.textContent || 0);
    if (!count) return null;
    return {
      sourceLabel: "MAL",
      measurements: {
        malAou: $("aouAverage")?.textContent || "",
        malQom: $("qomAverage")?.textContent || ""
      }
    };
  }

  function readSara() {
    if (!$("ataxiaConfirmed")?.checked) return null;
    return {
      sourceLabel: "SARA",
      measurements: { sara: $("totalScore")?.textContent || "" }
    };
  }

  function readMiniBestest() {
    if (typeof window.buildMiniDetail !== "function") return null;
    const detail = window.buildMiniDetail();
    if (detail.items.some((item) => item.score === null)) return null;
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
    const hhdKeys = Object.keys(measurements).filter((key) => key !== "comfortable10mSpeed");
    const sourceLabel = speed && speed !== "-"
      ? (hhdKeys.length ? "10m歩行・HHD" : "10m歩行")
      : "HHD";
    return {
      sourceLabel,
      weight: $("h-bw")?.value || "",
      hhdArmCm: arm,
      measurements
    };
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
