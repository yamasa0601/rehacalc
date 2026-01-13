/**
 * EMG walking app (iPhone robust)
 * - CSV parse hardened (BOM, CR, semicolon, tabs)
 * - NaN removal / interpolate tiny gaps
 * - Safe filtering (skip if too short / too many NaN)
 * - HS detection: MAD onset + peak fallback
 * - TA HS -> segment both channels
 * - Debug: shows "valid ratio" and env stats, plus HS reason
 *
 * Note: iPhone Safari + PWA caches are aggressive. Bump CACHE_NAME in service-worker.js.
 */

const $ = (id) => document.getElementById(id);

const state = {
  file1: null,
  file2: null,
  parsed1: null,
  parsed2: null,
  result1: null,
  result2: null,

  // video (optional)
  videoFile: null,
  roi: null,
  videoEvents: null,
  videoEnergy: null,
};

function setMsg(s, isErr=false){
  const el = $("msg");
  el.textContent = s;
  el.style.color = isErr ? "#ffb4b4" : "rgba(233,238,255,.72)";
}
function setVideoStatus(s, isErr=false){
  const el = $("videoStatus");
  if(!el) return;
  el.textContent = s;
  el.style.color = isErr ? "#ffb4b4" : "rgba(233,238,255,.72)";
}

/* =======================
   small utils
======================= */
function median(arr){
  const a = Array.from(arr).filter(Number.isFinite);
  if(a.length===0) return NaN;
  a.sort((x,y)=>x-y);
  const m = Math.floor(a.length/2);
  return (a.length%2===0) ? (a[m-1]+a[m])/2 : a[m];
}
function mad(arr){
  const med = median(arr);
  const dev = arr.map(v=>Math.abs(v-med));
  return median(dev);
}
function percentile(arr,p){
  const a = Array.from(arr).filter(Number.isFinite).sort((x,y)=>x-y);
  if(a.length===0) return NaN;
  const pos = (p/100)*(a.length-1);
  const i0 = Math.floor(pos);
  const i1 = Math.min(a.length-1,i0+1);
  const w = pos-i0;
  return a[i0]*(1-w)+a[i1]*w;
}
function linspace(a,b,n){
  const out = new Float64Array(n);
  const step = (b-a)/(n-1);
  for(let i=0;i<n;i++) out[i]=a+step*i;
  return out;
}
function reverseArray(x){
  const y=new Float64Array(x.length);
  for(let i=0;i<x.length;i++) y[i]=x[x.length-1-i];
  return y;
}
function absArray(x){
  const y=new Float64Array(x.length);
  for(let i=0;i<x.length;i++) y[i]=Math.abs(x[i]);
  return y;
}
function mean(arr){
  let s=0,c=0;
  for(const v of arr){ if(Number.isFinite(v)){ s+=v; c++; } }
  return c>0 ? s/c : NaN;
}
function varRatioFinite(x){
  let c=0;
  for(const v of x) if(Number.isFinite(v)) c++;
  return x.length ? c/x.length : 0;
}
function sanitizeArray(x){
  // remove non-finite by linear interpolation if short gaps; otherwise drop to 0
  const y = new Float64Array(x.length);
  for(let i=0;i<x.length;i++) y[i]=Number.isFinite(x[i])?x[i]:NaN;

  // forward fill start
  let first = -1;
  for(let i=0;i<y.length;i++){ if(Number.isFinite(y[i])){ first=i; break; } }
  if(first<0){
    // all NaN
    return new Float64Array(x.length).fill(0);
  }
  for(let i=0;i<first;i++) y[i]=y[first];

  // interpolate gaps
  let i=first;
  while(i<y.length){
    if(Number.isFinite(y[i])){ i++; continue; }
    const s=i-1;
    let e=i;
    while(e<y.length && !Number.isFinite(y[e])) e++;
    const left=y[s];
    const right = (e<y.length)?y[e]:left;
    const gap = e - s;
    // if huge gap, just fill with left
    for(let k=1;k<gap;k++){
      const w = k/gap;
      y[s+k] = left*(1-w)+right*w;
    }
    i=e;
  }
  return y;
}

/* =======================
   RMS envelope
======================= */
function movingRMS(x, win){
  const n=x.length;
  const half=Math.floor(win/2);
  const sq=new Float64Array(n);
  for(let i=0;i<n;i++) sq[i]=x[i]*x[i];

  const pref=new Float64Array(n+1);
  for(let i=0;i<n;i++) pref[i+1]=pref[i]+sq[i];

  const out=new Float64Array(n);
  for(let i=0;i<n;i++){
    const s=Math.max(0,i-half);
    const e=Math.min(n,i+half+1);
    const m=(pref[e]-pref[s])/(e-s);
    out[i]=Math.sqrt(m);
  }
  return out;
}

/* =======================
   Biquad filter (RBJ) + filtfilt
======================= */
function biquadCoeffs(type, fc, fs, Q=1/Math.SQRT2){
  const w0=(2*Math.PI*fc)/fs;
  const cosw0=Math.cos(w0);
  const sinw0=Math.sin(w0);
  const alpha=sinw0/(2*Q);
  let b0,b1,b2,a0,a1,a2;
  if(type==="lowpass"){
    b0=(1-cosw0)/2; b1=1-cosw0; b2=(1-cosw0)/2;
    a0=1+alpha; a1=-2*cosw0; a2=1-alpha;
  }else if(type==="highpass"){
    b0=(1+cosw0)/2; b1=-(1+cosw0); b2=(1+cosw0)/2;
    a0=1+alpha; a1=-2*cosw0; a2=1-alpha;
  }else throw new Error("unknown biquad type");
  return { b0:b0/a0, b1:b1/a0, b2:b2/a0, a1:a1/a0, a2:a2/a0 };
}
function biquadFilter(x,c){
  const y=new Float64Array(x.length);
  let x1=0,x2=0,y1=0,y2=0;
  for(let i=0;i<x.length;i++){
    const x0=x[i];
    const y0=c.b0*x0 + c.b1*x1 + c.b2*x2 - c.a1*y1 - c.a2*y2;
    y[i]=y0;
    x2=x1; x1=x0;
    y2=y1; y1=y0;
  }
  return y;
}
function filtfilt(x, coeffs){
  let y=biquadFilter(x,coeffs);
  y=reverseArray(y);
  y=biquadFilter(y,coeffs);
  y=reverseArray(y);
  return y;
}
function butter2x(x,type,fc,fs){
  const c=biquadCoeffs(type,fc,fs,1/Math.SQRT2);
  let y=filtfilt(x,c);
  y=filtfilt(y,c);
  return y;
}

/* =======================
   HS detection
======================= */
function risingEdges(mask){
  const starts=[], ends=[];
  for(let i=1;i<mask.length;i++){
    if(!mask[i-1] && mask[i]) starts.push(i);
    if(mask[i-1] && !mask[i]) ends.push(i);
  }
  if(mask[0]) starts.unshift(0);
  if(mask[mask.length-1]) ends.push(mask.length);
  const segs=[];
  for(let i=0;i<Math.min(starts.length,ends.length);i++) segs.push([starts[i],ends[i]]);
  return segs;
}
function localMaxima(y){
  const peaks=[];
  for(let i=1;i<y.length-1;i++) if(y[i]>y[i-1] && y[i]>=y[i+1]) peaks.push(i);
  return peaks;
}
function pickPeaksGreedy(y,minDist,thr,maxCount=null){
  const peaks=localMaxima(y).filter(i=>y[i]>thr);
  peaks.sort((i,j)=>y[j]-y[i]);
  const chosen=[];
  for(const idx of peaks){
    let ok=true;
    for(const c of chosen){ if(Math.abs(idx-c)<minDist){ ok=false; break; } }
    if(ok){
      chosen.push(idx);
      if(maxCount && chosen.length>=maxCount) break;
    }
  }
  chosen.sort((a,b)=>a-b);
  return chosen;
}
function refineToOnset(y,peakIdx,thrOff){
  let i=peakIdx;
  while(i>0 && y[i]>thrOff) i--;
  return i;
}
function applyOffsetAndMerge(hsRaw, offset, minGap, nMax){
  if(!hsRaw || hsRaw.length===0) return [];
  let hs=hsRaw.map(v=>v+offset).filter(v=>v>=0 && v<nMax).sort((a,b)=>a-b);
  const merged=[hs[0]];
  for(let i=1;i<hs.length;i++){
    if(hs[i]-merged[merged.length-1]>=minGap) merged.push(hs[i]);
  }
  return merged;
}
function scoreHS(hs, fs, plausible=[0.35,1.6], expected=null){
  if(!hs || hs.length<2) return Infinity;
  const ints=[];
  for(let i=1;i<hs.length;i++) ints.push((hs[i]-hs[i-1])/fs);
  const meanInt = mean(ints);
  const medInt = median(ints);
  const sd = Math.sqrt(mean(ints.map(v=>(v-meanInt)*(v-meanInt))) || 0);
  const cv = sd/(meanInt+1e-12);
  const [lo,hi]=plausible;
  const bad = ints.filter(v=>v<lo || v>hi).length/(ints.length||1);
  let s = 5*bad + 2*cv + 0.5*Math.abs(medInt-0.7);
  if(expected!==null) s += 0.8*Math.abs(hs.length-expected);
  return s;
}

function detectHS(env, fs, opts){
  const {
    expectedCount=null,
    minBurstMs=30,
    minGapMs=450,
    offsetMs=80,
    plausibleStep=[0.35,1.6],
  } = opts;

  const minBurst=Math.round(fs*minBurstMs/1000);
  const minGap=Math.round(fs*minGapMs/1000);
  const offset=Math.round(fs*offsetMs/1000);

  const med=median(env);
  let m=mad(env);

  if(!Number.isFinite(m) || m<1e-12){
    const p75=percentile(env,75), p25=percentile(env,25);
    m = (p75-p25)/1.349;
  }

  // MAD onset
  if(Number.isFinite(m) && m>1e-12){
    const ks=[];
    for(let k=1.2;k<=4.0+1e-9;k+=0.1) ks.push(Math.round(k*100)/100);

    let best=null;
    for(const k of ks){
      const thr = med + k*m;
      const mask=new Array(env.length);
      for(let i=0;i<env.length;i++) mask[i]=env[i]>thr;

      let segs=risingEdges(mask).filter(([s,e]) => (e-s)>=minBurst);
      if(segs.length===0) continue;

      let hs=segs.map(([s,_])=>s);
      hs=applyOffsetAndMerge(hs,offset,minGap,env.length);
      if(hs.length<2) continue;

      const sc=scoreHS(hs,fs,plausibleStep,expectedCount);
      if(best===null || sc<best.score) best={k,thr,hs,score:sc,mode:"MAD"};
    }
    if(best) return best;
  }

  // peak fallback
  const thrPeak = percentile(env, 80);
  const thrOff  = percentile(env, 60);
  let peaks = pickPeaksGreedy(env, Math.round(fs*minGapMs/1000), thrPeak, expectedCount);

  if(peaks.length>=2){
    const onsets = peaks.map(p=>refineToOnset(env,p,thrOff));
    const hs = applyOffsetAndMerge(onsets, Math.round(fs*offsetMs/1000), Math.round(fs*minGapMs/1000), env.length);
    if(hs.length>=2) return {k:"peak-fallback",thr:thrPeak,hs,score:scoreHS(hs,fs,plausibleStep,expectedCount),mode:"PEAK"};
  }

  throw new Error("HS推定に失敗しました（有効データ不足/NaN混入/信号弱い/閾値不適）。minGap/HP/offset/歩数を調整してください。");
}

function percentPeak(env, hs){
  let ref = NaN;
  if(hs && hs.length>=2){
    const s=hs[0], e=hs[hs.length-1];
    let mx=-Infinity;
    for(let i=s;i<e;i++) if(env[i]>mx) mx=env[i];
    ref = mx;
  }else{
    ref = Math.max(...env);
  }
  if(!Number.isFinite(ref) || ref<=0) return new Float64Array(env.length).fill(NaN);
  const out=new Float64Array(env.length);
  for(let i=0;i<env.length;i++) out[i]=(env[i]/ref)*100;
  return out;
}
function normalizeCycles(y, hs, nPoints){
  const grid=linspace(0,100,nPoints);
  const cycles=[];
  for(let i=0;i<hs.length-1;i++){
    const s=hs[i], e=hs[i+1];
    if(e<=s+2) continue;
    const segLen=e-s;
    const yseg=y.subarray(s,e);
    const out=new Float64Array(nPoints);
    const last=segLen-1;
    for(let j=0;j<nPoints;j++){
      const pos=(j/(nPoints-1))*last;
      const i0=Math.floor(pos);
      const i1=Math.min(last,i0+1);
      const w=pos-i0;
      out[j]=yseg[i0]*(1-w)+yseg[i1]*w;
    }
    cycles.push(out);
  }
  return {grid, cycles};
}
function meanAndSD(cycles){
  const n=cycles.length;
  if(n===0) return {mean:[], sd:[]};
  const m=cycles[0].length;
  const meanArr=new Float64Array(m);
  const sdArr=new Float64Array(m);

  for(let j=0;j<m;j++){
    let s=0,c=0;
    for(let i=0;i<n;i++){
      const v=cycles[i][j];
      if(Number.isFinite(v)){ s+=v; c++; }
    }
    meanArr[j]=c>0?s/c:NaN;
  }
  for(let j=0;j<m;j++){
    let s=0,c=0;
    for(let i=0;i<n;i++){
      const v=cycles[i][j];
      if(Number.isFinite(v)){ s+=(v-meanArr[j])*(v-meanArr[j]); c++; }
    }
    sdArr[j]=c>1?Math.sqrt(s/c):NaN;
  }
  return {mean:meanArr, sd:sdArr};
}
function applyHsToExistingResult(res, hs, nPoints){
  const hs2=hs.filter(i=>i>=0 && i<res.envRaw.length);
  const envPct=percentPeak(res.envRaw, hs2);
  const cyc=normalizeCycles(envPct, hs2, nPoints);
  const ms=meanAndSD(cyc.cycles);
  return {...res, hs:hs2, envPct, grid:cyc.grid, cycles:cyc.cycles, mean:ms.mean, sd:ms.sd};
}
function chooseEveryOtherIfDouble(hs, fs, expected=null){
  if(!hs || hs.length<6) return hs;
  const ints=[];
  for(let i=1;i<hs.length;i++) ints.push((hs[i]-hs[i-1])/fs);
  const medInt=median(ints);
  const looksDouble = (medInt<0.50) || (expected!==null && hs.length>=expected*1.6);
  if(!looksDouble) return hs;
  const a=hs.filter((_,i)=>i%2===0);
  const b=hs.filter((_,i)=>i%2===1);
  return (scoreHS(a,fs,[0.35,1.6],expected) <= scoreHS(b,fs,[0.35,1.6],expected)) ? a : b;
}

/* =======================
   CSV reading hardened
======================= */
async function readFileData(file){
  const buf = await file.arrayBuffer();
  let text="";
  try{
    text = new TextDecoder("utf-8", {fatal:false}).decode(buf);
  }catch{ text=""; }
  return {buf, text};
}
function bytesToLatin1String(buf){
  try{ return new TextDecoder("iso-8859-1",{fatal:false}).decode(buf); }
  catch{
    const bytes=new Uint8Array(buf);
    let out="";
    const CHUNK=4096;
    for(let i=0;i<bytes.length;i+=CHUNK){
      out += String.fromCharCode.apply(null, bytes.subarray(i,i+CHUNK));
    }
    return out;
  }
}
function normalizeCSVText(s){
  if(!s) return "";
  // remove BOM
  s = s.replace(/^\uFEFF/, "");
  // normalize line breaks
  s = s.replace(/\r\n/g,"\n").replace(/\r/g,"\n");
  // sometimes semicolon-separated
  // keep as-is; we will detect delimiter
  return s;
}
function detectDelimiter(line){
  const cComma = (line.match(/,/g)||[]).length;
  const cSemi  = (line.match(/;/g)||[]).length;
  const cTab   = (line.match(/\t/g)||[]).length;
  if(cTab>cComma && cTab>cSemi) return "\t";
  if(cSemi>cComma) return ";";
  return ",";
}
function parseCSV(buf, text, filename, fsDefault){
  const name = filename.replace(/\.csv$/i,"");
  const raw = bytesToLatin1String(buf);
  const rawNorm = normalizeCSVText(raw);

  // sensor format: contains "No,"
  const rawLines = rawNorm.split("\n");
  let start=-1, headerLine="";
  for(let i=0;i<Math.min(rawLines.length,5000);i++){
    const line = rawLines[i]||"";
    if(line.includes("No,")){
      start=i; headerLine=line.slice(line.indexOf("No,")); break;
    }
  }
  if(start>=0){
    const header = headerLine.split(",");
    const idxNo = header.findIndex(h=>(h||"").trim()==="No");
    let idxV = header.findIndex(h => (h||"").includes("mV") || (h||"").toLowerCase().includes("voltage") || (h||"").includes("電圧"));
    if(idxV<0) idxV=1;

    const no=[], x=[];
    for(let i=start+1;i<rawLines.length;i++){
      const line=rawLines[i];
      if(!line) continue;
      const parts=line.split(",");
      if(parts.length<Math.max(idxNo,idxV)+1) continue;
      const n=parseInt(parts[idxNo],10);
      const v=parseFloat(parts[idxV]);
      if(Number.isFinite(n) && Number.isFinite(v)){ no.push(n); x.push(v); }
    }
    const t=new Float64Array(no.length);
    for(let i=0;i<no.length;i++) t[i]=(no[i]-1)/fsDefault;

    return {name, fs:fsDefault, t, x:sanitizeArray(Float64Array.from(x)),
      columns:["No","value"], guessedTimeCol:"No", guessedSigCol:"value",
      table:null, format:"sensor"};
  }

  // table format
  let s = normalizeCSVText(text || rawNorm);
  let lines = s.split("\n").filter(l => (l||"").trim().length>0);
  if(lines.length<2) lines = rawLines.filter(l => (l||"").includes(","));

  const delim = detectDelimiter(lines[0]||"");
  const header = (lines[0]||"").split(delim).map(x=>(x||"").trim()).filter(Boolean);

  const rows=[];
  for(let i=1;i<lines.length;i++){
    const parts=(lines[i]||"").split(delim);
    if(parts.length < header.length) continue;
    rows.push(parts.slice(0, header.length));
  }

  const cols={};
  for(let j=0;j<header.length;j++){
    cols[header[j]] = new Float64Array(rows.length);
    for(let i=0;i<rows.length;i++){
      const v = parseFloat((rows[i][j]||"").replace(/"/g,""));
      cols[header[j]][i] = Number.isFinite(v) ? v : NaN;
    }
    cols[header[j]] = sanitizeArray(cols[header[j]]);
  }

  let timeCol = header.find(h=>["time","time_s","t","sec","seconds"].includes((h||"").toLowerCase())) || header[0];
  let sigCol  = header.find(h=>h!==timeCol) || header[0];

  return {name, fs:fsDefault, t:cols[timeCol], x:cols[sigCol],
    columns:header, guessedTimeCol:timeCol, guessedSigCol:sigCol,
    table:cols, format:"table"};
}
function populateSelect(sel, options, selected){
  if(!sel) return;
  sel.innerHTML="";
  for(const opt of options){
    const o=document.createElement("option");
    o.value=opt; o.textContent=opt;
    if(opt===selected) o.selected=true;
    sel.appendChild(o);
  }
}
function getSelectedSignal(parsed, sigSel, timeSel, fs){
  if(parsed.format==="sensor") return {t:parsed.t, x:parsed.x};
  const sig = sigSel?.value || parsed.guessedSigCol;
  const tim = timeSel?.value || parsed.guessedTimeCol;

  const x = sanitizeArray(parsed.table[sig]);
  const t0 = sanitizeArray(parsed.table[tim]);

  // detect if t0 is seconds
  const maxT = Math.max(...t0);
  const estDur = t0.length/fs;
  const isSeconds = Number.isFinite(maxT) && maxT > estDur*0.8;
  if(isSeconds) return {t:t0, x};

  // else: use sample-based seconds
  const tSec=new Float64Array(t0.length);
  for(let i=0;i<t0.length;i++) tSec[i]=(t0[i]-t0[0])/fs;
  return {t:tSec, x};
}

/* =======================
   compute pipeline (safe)
======================= */
function computePipeline(t, x, fs, params, expectedCount){
  const {hp, lp, rmsMs, minBurstMs, minGapMs, offsetMs, nPoints} = params;

  let y = sanitizeArray(x);
  const finiteRatio = varRatioFinite(y);

  if(finiteRatio < 0.95){
    // if too many non-finite, analysis will be unstable
    throw new Error(`CSVに欠損が多いです（有効率 ${(finiteRatio*100).toFixed(1)}%）。列選択やCSV形式を確認してください。`);
  }

  // demean
  const m = mean(y);
  for(let i=0;i<y.length;i++) y[i] = y[i]-m;

  // safe filter
  if(y.length < fs*2){
    // too short
    throw new Error("データが短すぎます（2秒未満）。");
  }

  // iPhoneでHP50が厳しい場合があるので、lp/hpの妥当性チェック
  if(hp >= fs/2) throw new Error("HighPassが高すぎます。");
  if(lp >= fs/2) throw new Error("LowPassが高すぎます。");

  if(hp > 0) y = butter2x(y, "highpass", hp, fs);
  if(lp > 0 && lp < fs/2) y = butter2x(y, "lowpass", lp, fs);

  y = absArray(y);
  const win = Math.max(1, Math.round(fs*rmsMs/1000));
  const env = movingRMS(y, win);

  const hsRes = detectHS(env, fs, {expectedCount, minBurstMs, minGapMs, offsetMs});
  const envPct = percentPeak(env, hsRes.hs);

  const cyc = normalizeCycles(envPct, hsRes.hs, nPoints);
  const ms = meanAndSD(cyc.cycles);

  return {
    t, envRaw: env, envPct,
    hs: hsRes.hs, k: hsRes.k, thr: hsRes.thr, mode: hsRes.mode || "MAD",
    grid: cyc.grid, cycles: cyc.cycles, mean: ms.mean, sd: ms.sd,
    debug: {
      finiteRatio,
      envMed: median(env),
      envMad: mad(env),
      envP95: percentile(env,95),
      envMax: Math.max(...env),
    }
  };
}

/* =======================
   Plot (Plotly)
======================= */
function plotCycle(targetId, res, showIndividual){
  const traces=[];
  if(showIndividual){
    for(const c of res.cycles){
      traces.push({
        x:Array.from(res.grid),
        y:Array.from(c),
        mode:"lines",
        line:{width:1,color:"rgba(200,200,200,0.35)"},
        hoverinfo:"skip",
        showlegend:false
      });
    }
  }
  if(res.cycles.length>0){
    const x=Array.from(res.grid);
    const meanArr=Array.from(res.mean);
    const sdArr=Array.from(res.sd);
    const upper=meanArr.map((v,i)=>v+sdArr[i]);
    const lower=meanArr.map((v,i)=>v-sdArr[i]);

    traces.push({x,y:upper,mode:"lines",line:{width:0},hoverinfo:"skip",showlegend:false});
    traces.push({x,y:lower,mode:"lines",fill:"tonexty",fillcolor:"rgba(59,130,246,0.18)",line:{width:0},hoverinfo:"skip",showlegend:false});
    traces.push({x,y:meanArr,mode:"lines",line:{width:3,color:"#3b82f6"},name:"mean"});
  }
  const layout={
    margin:{l:55,r:20,t:30,b:45},
    xaxis:{title:"gait cycle (%)",range:[0,100]},
    yaxis:{title:"%RMS (%peak)"},
    paper_bgcolor:"rgba(0,0,0,0)",
    plot_bgcolor:"rgba(0,0,0,0)",
    font:{color:"#e9eeff"}
  };
  return Plotly.newPlot(targetId,traces,layout,{displayModeBar:false,responsive:true});
}
function plotDebug(targetId, res){
  const shapes=(res.hs||[]).map(idx=>({
    type:"line",
    x0:res.t[idx], x1:res.t[idx],
    y0:0, y1:1,
    xref:"x", yref:"paper",
    line:{width:1,color:"rgba(255,255,255,0.25)"}
  }));
  const traces=[{
    x:Array.from(res.t),
    y:Array.from(res.envPct),
    mode:"lines",
    line:{width:2,color:"#22c55e"},
    name:"%peak env"
  }];
  const layout={
    margin:{l:55,r:20,t:30,b:45},
    xaxis:{title:"time (s)"},
    yaxis:{title:"%RMS (%peak)"},
    shapes,
    paper_bgcolor:"rgba(0,0,0,0)",
    plot_bgcolor:"rgba(0,0,0,0)",
    font:{color:"#e9eeff"},
    title:`valid ${(res.debug.finiteRatio*100).toFixed(1)}% | env max ${res.debug.envMax.toFixed(3)} | mode ${res.mode}`
  };
  return Plotly.newPlot(targetId,traces,layout,{displayModeBar:false,responsive:true});
}

/* =======================
   Downloads
======================= */
function makeHSCSV(t, hs){
  let s="hs_index,No,time_s\n";
  for(let i=0;i<hs.length;i++){
    s += `${i+1},${hs[i]+1},${t[hs[i]]}\n`;
  }
  return s;
}
function downloadText(filename,text){
  const blob=new Blob([text],{type:"text/plain;charset=utf-8"});
  const url=URL.createObjectURL(blob);
  const a=document.createElement("a");
  a.href=url; a.download=filename;
  document.body.appendChild(a); a.click(); a.remove();
  URL.revokeObjectURL(url);
}
async function downloadPlotPNG(plotId, filename){
  const div=document.getElementById(plotId);
  const url=await Plotly.toImage(div,{format:"png",height:800,width:1200,scale:2});
  const a=document.createElement("a");
  a.href=url; a.download=filename;
  document.body.appendChild(a); a.click(); a.remove();
}

/* =======================
   Main run
======================= */
function clearPlots(){
  ["plot1","plot2","debug1","debug2"].forEach(id=>{
    const el=document.getElementById(id);
    if(el) el.innerHTML="";
  });
  ["hsn1","hsn2","kmad1","kmad2","cycn1","cycn2","hssrc","hssrc2"].forEach(id=>{
    const el=$(id); if(el) el.textContent="-";
  });
}

async function handleFile(file, sigSel, timeSel, fs){
  if(!file) return null;
  const {buf,text}=await readFileData(file);
  const parsed=parseCSV(buf,text,file.name,fs);

  if(parsed.format==="table"){
    populateSelect(sigSel, parsed.columns, parsed.guessedSigCol);
    populateSelect(timeSel, parsed.columns, parsed.guessedTimeCol);
  }else{
    populateSelect(sigSel, ["value"], "value");
    populateSelect(timeSel, ["No"], "No");
  }
  return parsed;
}

async function run(){
  try{
    setMsg("読み込み中…");
    clearPlots();

    const fs=parseFloat($("fs").value);
    const hp=parseFloat($("hp").value);
    const lp=parseFloat($("lp").value);
    const rmsMs=parseFloat($("rmsms").value);
    const minBurstMs=parseFloat($("minburst").value);
    const minGapMs=parseFloat($("mingap").value);
    const offsetMs=parseFloat($("offset").value);
    const steps=parseInt($("steps").value,10);
    const nPoints=parseInt($("npoints").value,10);
    const showIndividual = $("showind").value==="1";
    const baselineSec = parseFloat($("baselineSec")?.value || "3");

    if(!state.file1 || !state.file2){
      setMsg("CSVを2つ選んでください。", true);
      return;
    }

    state.parsed1 = await handleFile(state.file1, $("sigcol1"), $("timecol1"), fs);
    state.parsed2 = await handleFile(state.file2, $("sigcol2"), $("timecol2"), fs);

    const expTA = (steps>0) ? Math.floor(steps/2) : null;

    const s1 = getSelectedSignal(state.parsed1, $("sigcol1"), $("timecol1"), fs);
    const s2 = getSelectedSignal(state.parsed2, $("sigcol2"), $("timecol2"), fs);

    // iPhoneで落ちやすいので、まず短いエラーチェック
    setMsg("解析中…");

    // まずTAを解析（ここで失敗するならCSV/パラメータ問題）
    const params = {hp, lp, rmsMs, minBurstMs, minGapMs, offsetMs, nPoints};
    state.result1 = computePipeline(s1.t, s1.x, fs, params, expTA);

    // TA HS refine: baselineSec以前削除
    let hsRef = state.result1.hs.filter(i => state.result1.t[i] >= baselineSec);

    // 倍検出なら1つおき
    hsRef = chooseEveryOtherIfDouble(hsRef, fs, expTA);

    // CSV2解析（HSは後で上書き）
    state.result2 = computePipeline(s2.t, s2.x, fs, params, null);

    // apply TA HS to both
    state.result1 = applyHsToExistingResult(state.result1, hsRef, nPoints);
    state.result2 = applyHsToExistingResult(state.result2, hsRef, nPoints);

    // UI
    $("name1").textContent = state.parsed1.name || "CSV1";
    $("name2").textContent = state.parsed2.name || "CSV2";
    $("hsn1").textContent = String(state.result1.hs.length);
    $("hsn2").textContent = String(state.result2.hs.length);
    $("kmad1").textContent = String(state.result1.k);
    $("kmad2").textContent = "(TA参照)";
    $("cycn1").textContent = String(Math.max(0, state.result1.hs.length-1));
    $("cycn2").textContent = String(Math.max(0, state.result2.hs.length-1));
    $("hssrc").textContent = `TA(${state.result1.mode}) + baseline(${baselineSec}s) + decimate`;
    $("hssrc2").textContent = `TA(${state.result1.mode})`;

    await plotCycle("plot1", state.result1, showIndividual);
    await plotDebug("debug1", state.result1);
    await plotCycle("plot2", state.result2, showIndividual);
    await plotDebug("debug2", state.result2);

    setMsg("完了。iPhoneでもHSが出るはずです（debugにvalid%/envが出ます）。");
  }catch(e){
    console.error(e);
    setMsg(String(e?.message || e), true);
  }
}

/* =======================
   UI wiring
======================= */
$("file1")?.addEventListener("change",(ev)=>{
  state.file1 = ev.target.files?.[0] || null;
  if(state.file1) setMsg(`CSV1 選択: ${state.file1.name}`);
});
$("file2")?.addEventListener("change",(ev)=>{
  state.file2 = ev.target.files?.[0] || null;
  if(state.file2) setMsg(`CSV2 選択: ${state.file2.name}`);
});
$("run")?.addEventListener("click", run);

$("clear")?.addEventListener("click", ()=>{
  state.file1=null; state.file2=null;
  if($("file1")) $("file1").value="";
  if($("file2")) $("file2").value="";
  clearPlots();
  setMsg("クリアしました。");
});

$("dlhs1")?.addEventListener("click",()=>{
  if(!state.result1) return;
  downloadText(`${state.parsed1?.name||"CSV1"}_HS.csv`, makeHSCSV(state.result1.t, state.result1.hs));
});
$("dlhs2")?.addEventListener("click",()=>{
  if(!state.result2) return;
  downloadText(`${state.parsed2?.name||"CSV2"}_HS_byTA.csv`, makeHSCSV(state.result2.t, state.result2.hs));
});
$("dlpng1")?.addEventListener("click", async ()=>{
  if(!state.result1) return;
  await downloadPlotPNG("plot1", `${state.parsed1?.name||"CSV1"}_cycle.png`);
});
$("dlpng2")?.addEventListener("click", async ()=>{
  if(!state.result2) return;
  await downloadPlotPNG("plot2", `${state.parsed2?.name||"CSV2"}_cycle_byTA.png`);
});

setVideoStatus("動画は任意（HSはCSVだけでも推定します）");
setMsg("CSVを2つ選択して「解析してプロット」を押してください。");
