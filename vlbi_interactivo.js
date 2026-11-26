(function () {
  "use strict";
  var C = 299792458, UAS = Math.PI / 180 / 3600 / 1e6, FS = 1000, BW = 200, A1 = 1, SEG = 1024;

  // ---------- utilidades numéricas ----------
  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(rng) {
    var u = 0; while (u === 0) u = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
  }
  function fft(re, im, inv) {
    var n = re.length, i, j, k, len, bit, t;
    for (i = 1, j = 0; i < n; i++) {
      bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (len = 2; len <= n; len <<= 1) {
      var ang = 2 * Math.PI / len * (inv ? 1 : -1), wr = Math.cos(ang), wi = Math.sin(ang), half = len >> 1;
      for (i = 0; i < n; i += len) {
        var cr = 1, ci = 0;
        for (k = 0; k < half; k++) {
          var a = i + k, b = a + half;
          var xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
          re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
          t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
    if (inv) for (i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }
  function sinc(x) { return x === 0 ? 1 : Math.sin(Math.PI * x) / (Math.PI * x); }

  function bandNoise(n, rng) {
    var re = new Float64Array(n), im = new Float64Array(n), k;
    for (k = 0; k < n; k++) re[k] = gauss(rng);
    fft(re, im, false);
    for (k = 0; k < n; k++) {
      var kk = k <= n / 2 ? k : k - n;
      if (Math.abs(kk) * FS / n > BW) { re[k] = 0; im[k] = 0; }
    }
    fft(re, im, true);
    var m = 0, v = 0;
    for (k = 0; k < n; k++) m += re[k];
    m /= n;
    for (k = 0; k < n; k++) v += (re[k] - m) * (re[k] - m);
    var sd = Math.sqrt(v / n);
    for (k = 0; k < n; k++) re[k] = (re[k] - m) / sd;
    return re;
  }
  function delayFFT(x, tau) { // x(t - tau), retardo posiblemente fraccionario
    var n = x.length, re = Float64Array.from(x), im = new Float64Array(n), k;
    fft(re, im, false);
    for (k = 0; k < n; k++) {
      var kk = k <= n / 2 ? k : k - n, ph = -2 * Math.PI * (kk * FS / n) * tau;
      var c = Math.cos(ph), s = Math.sin(ph), r = re[k], i2 = im[k];
      re[k] = r * c - i2 * s; im[k] = r * s + i2 * c;
    }
    fft(re, im, true);
    return re;
  }
  function xcorr(x1, x2, L) { // R12[k] = (1/N) sum x1[m] x2[m+k]
    var n = x1.length, nf = 2 * n, k;
    var ar = new Float64Array(nf), ai = new Float64Array(nf), br = new Float64Array(nf), bi = new Float64Array(nf);
    ar.set(x1); br.set(x2);
    fft(ar, ai, false); fft(br, bi, false);
    for (k = 0; k < nf; k++) {
      var re = ar[k] * br[k] + ai[k] * bi[k], im = ar[k] * bi[k] - ai[k] * br[k];
      ar[k] = re; ai[k] = im;
    }
    fft(ar, ai, true);
    var lags = [], r = [];
    for (k = -L; k <= L; k++) { lags.push(k); r.push(ar[(k + nf) % nf] / n); }
    return { lags: lags, r: r };
  }

  // ---------- parte A: señales, correlación y fase ----------
  function runA(p) {
    var n = 1 << p.nexp, rng = mulberry32(p.seed * 7919 + 13), i, k;
    var s = bandNoise(n, rng);
    var tt = p.tau + p.eps;               // ms; con FS = 1 kHz equivale a muestras
    var sd = delayFFT(s, tt / 1000);
    var x1 = new Float64Array(n), x2 = new Float64Array(n);
    for (i = 0; i < n; i++) {
      x1[i] = A1 * s[i] + p.sigma * gauss(rng);
      x2[i] = p.a2 * sd[i] + p.sigma * gauss(rng);
    }
    var kt = tt / 1000 * FS;
    var cc = xcorr(x1, x2, 600), lags = cc.lags, r = cc.r;
    var kmax = 0;
    for (k = 1; k < r.length; k++) if (r[k] > r[kmax]) kmax = k;
    var est = lags[kmax];
    if (kmax > 0 && kmax < r.length - 1) {
      var y0 = r[kmax - 1], y1 = r[kmax], y2 = r[kmax + 1], den = y0 - 2 * y1 + y2;
      if (den !== 0) est += 0.5 * (y0 - y2) / den;
    }
    var peak = -Infinity, sum = 0, sum2 = 0, cnt = 0;
    for (k = 0; k < r.length; k++) {
      if (Math.abs(lags[k] - kt) <= 3 + 0.5 && r[k] > peak) peak = r[k];
      if (Math.abs(lags[k] - kt) > 60) { sum += r[k]; sum2 += r[k] * r[k]; cnt++; }
    }
    var mean = sum / cnt, sd0 = Math.sqrt(Math.max(sum2 / cnt - mean * mean, 1e-30));
    var snrPeak = peak / sd0;
    var nearest = Math.round(kt);
    var expected = A1 * p.a2 * sinc(2 * BW * (nearest - kt) / FS);

    // N = 256 (comparación)
    var N0 = Math.min(256, n), short = [], sl = [];
    for (var lag = -150; lag <= 150; lag++) {
      var acc = 0;
      for (var m = 0; m < N0; m++) { var q = m + lag; if (q >= 0 && q < N0) acc += x1[m] * x2[q]; }
      sl.push(lag); short.push(acc / N0);
    }

    // Espectro cruzado por segmentos (Welch, ventana de Hann)
    var M = Math.floor(n / SEG), w = new Float64Array(SEG), sr = new Float64Array(SEG / 2 + 1), si = new Float64Array(SEG / 2 + 1);
    for (i = 0; i < SEG; i++) w[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / SEG);
    var a_r = new Float64Array(SEG), a_i = new Float64Array(SEG), b_r = new Float64Array(SEG), b_i = new Float64Array(SEG);
    for (var sg = 0; sg < M; sg++) {
      for (i = 0; i < SEG; i++) { a_r[i] = x1[sg * SEG + i] * w[i]; a_i[i] = 0; b_r[i] = x2[sg * SEG + i] * w[i]; b_i[i] = 0; }
      fft(a_r, a_i, false); fft(b_r, b_i, false);
      for (k = 0; k <= SEG / 2; k++) {
        sr[k] += a_r[k] * b_r[k] + a_i[k] * b_i[k];
        si[k] += a_r[k] * b_i[k] - a_i[k] * b_r[k];
      }
    }
    var k0 = lags[kmax], fx = [], ph = [], wt = [];
    for (k = 0; k <= SEG / 2; k++) {
      var f = k * FS / SEG;
      if (f >= 5 && f <= 0.95 * BW) {
        var ang = 2 * Math.PI * f * k0 / FS, c = Math.cos(ang), s2 = Math.sin(ang);
        var rr = sr[k] * c - si[k] * s2, ii = sr[k] * s2 + si[k] * c;
        fx.push(f); ph.push(Math.atan2(ii, rr)); wt.push(rr * rr + ii * ii);
      }
    }
    var sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
    for (i = 0; i < fx.length; i++) {
      sw += wt[i]; swx += wt[i] * fx[i]; swy += wt[i] * ph[i];
      swxx += wt[i] * fx[i] * fx[i]; swxy += wt[i] * fx[i] * ph[i];
    }
    var slope = 0, icpt = 0;
    if (fx.length > 1) {
      slope = (sw * swxy - swx * swy) / (sw * swxx - swx * swx);
      icpt = (swy - slope * swx) / sw;
    }
    var tauFit = (k0 / FS - slope / (2 * Math.PI)) * 1000; // ms

    var i0 = Math.min(2000, n - 250), t = [], c1 = [], c2 = [], n1 = [], n2 = [];
    for (i = 0; i < 250; i++) {
      t.push((i0 + i) / FS); c1.push(A1 * s[i0 + i]); c2.push(p.a2 * sd[i0 + i]);
      n1.push(x1[i0 + i]); n2.push(x2[i0 + i]);
    }
    var P1 = A1 * A1 + p.sigma * p.sigma, P2 = p.a2 * p.a2 + p.sigma * p.sigma;
    return {
      n: n, kt: kt, tt: tt, lags: lags, r: r, shortLags: sl, short: short, peak: peak, expected: expected, snrPeak: snrPeak,
      estPeak: est / FS * 1000, estPhase: tauFit, fx: fx, ph: ph, slope: slope, icpt: icpt,
      t: t, c1: c1, c2: c2, n1: n1, n2: n2,
      snr1: 10 * Math.log10(A1 * A1 / (p.sigma * p.sigma)), snr2: 10 * Math.log10(p.a2 * p.a2 / (p.sigma * p.sigma)),
      theory: A1 * p.a2 / Math.sqrt(P1 * P2) * Math.sqrt(n)
    };
  }

  // ---------- parte B: cobertura uv, visibilidades, imagen sucia ----------
  var STATIONS = [
    { name: "ALMA", lat: -23.02, lon: -67.75 },
    { name: "LMT", lat: 18.99, lon: -97.31 },
    { name: "SMT", lat: 32.70, lon: -109.89 },
    { name: "SMA", lat: 19.82, lon: -155.48 },
    { name: "IRAM-30m", lat: 37.07, lon: -3.39 }
  ];
  function uvCoverage(on, decDeg, lam) {
    var R = 6371e3, d = decDeg * Math.PI / 180, pos = [], u = [], v = [], i, j, s;
    for (s = 0; s < STATIONS.length; s++) {
      var la = STATIONS[s].lat * Math.PI / 180, lo = STATIONS[s].lon * Math.PI / 180;
      pos.push([R * Math.cos(la) * Math.cos(lo), R * Math.cos(la) * Math.sin(lo), R * Math.sin(la)]);
    }
    for (var step = 0; step < 144; step++) {
      var H = (15 * step / 6) * Math.PI / 180, up = [];
      for (s = 0; s < STATIONS.length; s++) {
        if (!on[s]) continue;
        var la2 = STATIONS[s].lat * Math.PI / 180, lo2 = STATIONS[s].lon * Math.PI / 180;
        var sinEl = Math.sin(la2) * Math.sin(d) + Math.cos(la2) * Math.cos(d) * Math.cos(H + lo2);
        if (Math.asin(sinEl) * 180 / Math.PI > 10) up.push(s);
      }
      for (i = 0; i < up.length; i++) for (j = i + 1; j < up.length; j++) {
        var b = [pos[up[j]][0] - pos[up[i]][0], pos[up[j]][1] - pos[up[i]][1], pos[up[j]][2] - pos[up[i]][2]];
        u.push((Math.sin(H) * b[0] + Math.cos(H) * b[1]) / lam);
        v.push((-Math.sin(d) * Math.cos(H) * b[0] + Math.sin(d) * Math.sin(H) * b[1] + Math.cos(d) * b[2]) / lam);
      }
    }
    return { u: u, v: v };
  }
  var GN = 73, GSTEP = 2, GEXT = 72; // rejilla de la imagen en microsegundos de arco
  function runB(p) {
    var lam = C / 230e9, uv = uvCoverage(p.stations, p.dec, lam), nu = uv.u.length, i, x, y, q;
    var coord = []; for (x = 0; x < GN; x++) coord.push((-GEXT + GSTEP * x));
    var r0 = p.ring / 2, wid = 3.5, img = new Float64Array(GN * GN), sym = new Float64Array(GN * GN), sI = 0, sS = 0;
    for (y = 0; y < GN; y++) for (x = 0; x < GN; x++) {
      var L = coord[x], Mm = coord[y], rr = Math.hypot(L, Mm), phi = Math.atan2(Mm, L);
      var g = Math.exp(-(rr - r0) * (rr - r0) / (2 * wid * wid));
      sym[y * GN + x] = g; img[y * GN + x] = g * (1 + p.asym * Math.cos(phi + Math.PI / 2));
      sI += img[y * GN + x]; sS += g;
    }
    for (q = 0; q < img.length; q++) { img[q] /= sI; sym[q] /= sS; }
    var thr = 1e-4 * Math.max.apply(null, img), pix = [];
    for (q = 0; q < img.length; q++) if (img[q] > thr) pix.push(q);
    // trigonometría por eje: cos/sen(2 pi u l) y (2 pi v m)
    var cu = [], su = [], cv = [], sv = [];
    for (i = 0; i < nu; i++) {
      var a = new Float64Array(GN), b = new Float64Array(GN), c = new Float64Array(GN), d2 = new Float64Array(GN);
      for (x = 0; x < GN; x++) {
        var ang1 = 2 * Math.PI * uv.u[i] * coord[x] * UAS, ang2 = 2 * Math.PI * uv.v[i] * coord[x] * UAS;
        a[x] = Math.cos(ang1); b[x] = Math.sin(ang1); c[x] = Math.cos(ang2); d2[x] = Math.sin(ang2);
      }
      cu.push(a); su.push(b); cv.push(c); sv.push(d2);
    }
    var rng = mulberry32(p.seed * 104729 + 7), sig = p.vnoise / 100, Vr = new Float64Array(nu), Vi = new Float64Array(nu), amp = [], rho = [];
    for (i = 0; i < nu; i++) {
      var re = 0, im = 0;
      for (var z = 0; z < pix.length; z++) {
        var px = pix[z] % GN, py = (pix[z] - px) / GN;
        var co = cu[i][px] * cv[i][py] - su[i][px] * sv[i][py], si = su[i][px] * cv[i][py] + cu[i][px] * sv[i][py];
        re += img[pix[z]] * co; im -= img[pix[z]] * si;
      }
      Vr[i] = re + sig * gauss(rng); Vi[i] = im + sig * gauss(rng);
      amp.push(Math.hypot(Vr[i], Vi[i])); rho.push(Math.hypot(uv.u[i], uv.v[i]) / 1e9);
    }
    var dirty = new Float64Array(GN * GN), beam = new Float64Array(GN * GN);
    if (nu > 0) {
      for (y = 0; y < GN; y++) for (x = 0; x < GN; x++) {
        var sd = 0, sb = 0;
        for (i = 0; i < nu; i++) {
          var co2 = cu[i][x] * cv[i][y] - su[i][x] * sv[i][y], si2 = su[i][x] * cv[i][y] + cu[i][x] * sv[i][y];
          sd += Vr[i] * co2 - Vi[i] * si2; sb += co2;
        }
        dirty[y * GN + x] = sd / nu; beam[y * GN + x] = sb / nu;
      }
    }
    var bmax = 0; for (i = 0; i < nu; i++) bmax = Math.max(bmax, Math.hypot(uv.u[i], uv.v[i]));
    var rmax = bmax > 0 ? 1.15 * bmax : 10e9, curveX = [], curveY = [], sPix = [];
    for (q = 0; q < sym.length; q++) if (sym[q] > 1e-4) sPix.push(q);
    for (i = 0; i <= 240; i++) {
      var rh = rmax * i / 240, acc = 0;
      for (var zz = 0; zz < sPix.length; zz++) acc += sym[sPix[zz]] * Math.cos(2 * Math.PI * rh * coord[sPix[zz] % GN] * UAS);
      curveX.push(rh / 1e9); curveY.push(Math.abs(acc));
    }
    var minPos = null;
    for (i = 1; i < curveY.length - 1; i++) if (curveY[i] < curveY[i - 1] && curveY[i] < curveY[i + 1] && curveY[i] < 0.5) { minPos = curveX[i]; break; }
    var side = 0, cen = beam[Math.floor(GN / 2) * GN + Math.floor(GN / 2)] || 1;
    for (y = 0; y < GN; y++) for (x = 0; x < GN; x++) if (Math.hypot(coord[x], coord[y]) > 15) side = Math.max(side, beam[y * GN + x] / cen);
    return {
      lam: lam, u: uv.u, v: uv.v, nu: nu, bmax: bmax, img: img, dirty: dirty, beam: beam, amp: amp, rho: rho,
      curveX: curveX, curveY: curveY, minPos: minPos, side: side, ring: p.ring,
      resol: bmax > 0 ? 1 / bmax / UAS : null, thin: 2.405 / (Math.PI * p.ring * UAS) / 1e9
    };
  }

  var core = { runA: runA, runB: runB, STATIONS: STATIONS };
  if (typeof module !== "undefined") module.exports = core;
  if (typeof document === "undefined") return;

  // ---------- interfaz ----------
  var DEF = { tau: 37.5, eps: 0, sigma: 3, a2: 0.6, nexp: 16, seed: 7, ring: 42, asym: 0.4, vnoise: 2, dec: 12.39 };
  var P = Object.assign({}, DEF), ON = STATIONS.map(function () { return true; });
  var SPEC = [
    {
      g: "Señales (supuestos ilustrativos)",
      badge: "correlación y fase",
      sub: "Parametros que definen el tipo de señal que llega del espacio"
    },
    {
      id: "tau",
      t: "Retardo geométrico τ_g",
      desc: "Diferencia de tiempo con la que llega cada señal a cada antena",
      min: 5, max: 60, step: 0.5, u: " ms", grp: "A"
    },
    {
      id: "eps",
      t: "Error de reloj ε",
      desc: "Desajuste entre los relojes internos de las antenas",
      min: -10, max: 10, step: 0.5, u: " ms", grp: "A"
    },
    {
      id: "sigma",
      t: "Ruido por antena σ_n",
      desc: "Interferencia térmica y atmosférica que oculta la señal débil del agujero negro.",
      min: 0.5, max: 8, step: 0.5, u: "", grp: "A"
    },
    {
      id: "a2",
      t: "Ganancia estación 2 A₂",
      desc: "Potencia de la señal recibida por la segunda antena y comparada por la primera/anterior",
      min: 0.2, max: 1.5, step: 0.1, u: "", grp: "A"
    },
    {
      id: "nexp",
      t: "Muestras integradas N",
      desc: "Canrtidad de mediciiones combinadas, Mientras más haya mayor calidad habrá",
      min: 10, max: 18, step: 1, grp: "A", fmt: function (v) { return (1 << v).toLocaleString("es-MX"); }
    },
    {
      g: "Cielo y estaciones",
      badge: "imagen y plano uv",
      sub: "Se define la forma del anillo de luz en base a los tipos de telescopios/antenas"
    },
    {
      id: "ring",
      t: "Diámetro del anillo",
      desc: "Tamaño angular de la sombra del agujero negro en microsegundos de arco (μas).",
      min: 20, max: 80, step: 1, u: " μas", grp: "B"
    },
    {
      id: "asym",
      t: "Asimetría de brillo",
      desc: "Efecto Doppler relativista que hace brillar más un lado del disco que el otro.",
      min: 0, max: 0.8, step: 0.05, u: "", grp: "B"
    },
    {
      id: "vnoise",
      t: "Ruido en visibilidades",
      desc: "Imperfecciones en los datos que agregan manchas o distorsión a la foto final.",
      min: 0, max: 10, step: 0.5, u: " %", grp: "B"
    },
    {
      id: "dec",
      t: "Declinación de la fuente",
      desc: "Qué tan alto está el agujero negro en el cielo, lo que cambia cómo lo captan las antenas.",
      min: -20, max: 60, step: 0.5, u: "°", grp: "B"
    }
  ];
  var dirtyA = true, dirtyB = true, timer = null, state = { A: null, B: null };
  var $ = function (id) { return document.getElementById(id); };

  function buildControls() {
    var host = $("controls"), html = "";
    SPEC.forEach(function (s) {
      if (s.g) {
        html += '<div class="group"><div class="group-head"><h2>' + s.g + '</h2><span class="badge">' + s.badge + '</span></div>' +
          '<p class="group-sub">' + s.sub + '</p></div>';
        return;
      }
      html += '<div class="ctl"><label for="' + s.id + '"><span>' + s.t + '</span><output id="o-' + s.id + '"></output></label>' +
        '<span class="ctl-sub">' + s.desc + '</span>' +
        '<input type="range" id="' + s.id + '" min="' + s.min + '" max="' + s.max + '" step="' + s.step + '" value="' + P[s.id] + '"></div>';
      if (s.id === "nexp") {
        html += '<div class="stations" id="stations"><div class="stations-head"><span class="badge">Estaciones activas</span>' +
          '<span class="ctl-sub">Enciende o apaga telescopios para ver cómo cambia la nitidez de la foto.</span></div></div>';
      }
    });
    html += '<div class="actions"><button class="primary" id="btn-work">Valores del trabajo (N = 262 144)</button><button id="btn-seed">Otra realización de ruido</button><button id="btn-reset">Restablecer</button><span class="status" id="status"></span></div>';
    host.innerHTML = html;
    // estaciones: se mueven al final del bloque de cielo
    var st = $("stations"), ringCtl = $("ring").closest(".ctl");
    STATIONS.forEach(function (s, i) {
      var l = document.createElement("label");
      l.innerHTML = '<input type="checkbox" id="st' + i + '" checked> ' + s.name;
      st.appendChild(l);
      l.firstChild.addEventListener("change", function (e) { ON[i] = e.target.checked; dirtyB = true; schedule(); });
    });
    host.insertBefore(st, ringCtl);
    SPEC.forEach(function (s) {
      if (s.g) return;
      var el = $(s.id);
      el.addEventListener("input", function () {
        P[s.id] = parseFloat(el.value); showOut(s); if (s.grp === "A") dirtyA = true; else dirtyB = true; schedule();
      });
      showOut(s);
    });
    $("btn-work").addEventListener("click", function () { P = Object.assign({}, DEF, { nexp: 18 }); syncControls(); dirtyA = dirtyB = true; schedule(); });
    $("btn-reset").addEventListener("click", function () { P = Object.assign({}, DEF); syncControls(); dirtyA = dirtyB = true; schedule(); });
    $("btn-seed").addEventListener("click", function () { P.seed += 1; dirtyA = dirtyB = true; schedule(); });
    $("show256").addEventListener("change", render);
  }
  function showOut(s) {
    var v = P[s.id], txt = s.fmt ? s.fmt(v) : (Math.abs(v) < 10 && s.step < 1 ? v.toFixed(s.step < 0.1 ? 2 : 1) : String(v));
    $("o-" + s.id).textContent = txt + (s.u || "");
  }
  function syncControls() {
    SPEC.forEach(function (s) { if (s.g) return; $(s.id).value = P[s.id]; showOut(s); });
    ON = STATIONS.map(function () { return true; });
    STATIONS.forEach(function (s, i) { $("st" + i).checked = true; });
  }
  function schedule() { $("status").textContent = "Calculando…"; clearTimeout(timer); timer = setTimeout(compute, 40); }
  function compute() {
    if (dirtyA) { state.A = runA(P); dirtyA = false; }
    if (dirtyB) { state.B = runB({ stations: ON, dec: P.dec, ring: P.ring, asym: P.asym, vnoise: P.vnoise, seed: P.seed }); dirtyB = false; }
    render(); $("status").textContent = "";
  }

  // ---------- dibujo ----------
  function colors() {
    var cs = getComputedStyle(document.documentElement), g = function (n) { return cs.getPropertyValue(n).trim(); };
    return {
      fg: g("--fg"), muted: g("--muted"), grid: g("--grid"),
      l1: g("--l1"), l2: g("--l2"), ref: g("--ref"), acc: g("--acc"), panel: g("--panel"),
      mono: g("--font-mono") || "monospace"
    };
  }
  function prep(cv) {
    var dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
    if (!(w > 0 && h > 0)) return null;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
    var g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, h);
    return { g: g, w: w, h: h };
  }
  function niceTicks(a, b, n) {
    var span = b - a; if (!(span > 0)) return [a];
    var raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
    var st = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p, out = [];
    for (var v = Math.ceil(a / st - 1e-9) * st; v <= b + 1e-9 * st; v += st) out.push(+v.toFixed(10));
    return out;
  }
  function fmt(v) { return String(Number(v.toPrecision(4))); }
  function plot(cv, o) {
    var pr = prep(cv); if (!pr) return;
    var g = pr.g, w = pr.w, h = pr.h, col = colors();
    var m = { l: o.ml || 44, r: 10, t: 8, b: o.mb || 30 };
    var x0 = o.x[0], x1 = o.x[1], y0 = o.y[0], y1 = o.y[1];
    var sx = function (v) { return m.l + (v - x0) / (x1 - x0) * (w - m.l - m.r); };
    var sy = function (v) { return h - m.b - (v - y0) / (y1 - y0) * (h - m.t - m.b); };
    g.font = "10.5px " + col.mono;
    g.lineWidth = 1; g.strokeStyle = col.grid; g.fillStyle = col.muted;
    g.textAlign = "center"; g.textBaseline = "top";
    niceTicks(x0, x1, w < 420 ? 4 : 7).forEach(function (v) {
      var X = Math.round(sx(v)) + 0.5; g.beginPath(); g.moveTo(X, m.t); g.lineTo(X, h - m.b); g.stroke(); g.fillText(fmt(v), X, h - m.b + 4);
    });
    g.textAlign = "right"; g.textBaseline = "middle";
    niceTicks(y0, y1, 5).forEach(function (v) {
      var Y = Math.round(sy(v)) + 0.5; g.beginPath(); g.moveTo(m.l, Y); g.lineTo(w - m.r, Y); g.stroke(); g.fillText(fmt(v), m.l - 5, Y);
    });
    g.save(); g.beginPath(); g.rect(m.l, m.t, w - m.l - m.r, h - m.t - m.b); g.clip();
    (o.vlines || []).forEach(function (vl) {
      g.strokeStyle = vl.color; g.setLineDash([4, 3]); g.lineWidth = 1.2; g.beginPath(); g.moveTo(sx(vl.x), m.t); g.lineTo(sx(vl.x), h - m.b); g.stroke(); g.setLineDash([]);
    });
    o.series.forEach(function (s) {
      g.globalAlpha = s.alpha || 1;
      if (s.type === "dots") {
        g.fillStyle = s.color;
        for (var i = 0; i < s.x.length; i++) { g.beginPath(); g.arc(sx(s.x[i]), sy(s.y[i]), s.r || 2, 0, 6.2832); g.fill(); }
      } else {
        g.strokeStyle = s.color; g.lineWidth = s.lw || 1.4; g.beginPath();
        for (var j = 0; j < s.x.length; j++) { if (j === 0) g.moveTo(sx(s.x[j]), sy(s.y[j])); else g.lineTo(sx(s.x[j]), sy(s.y[j])); }
        g.stroke();
      }
      g.globalAlpha = 1;
    });
    g.restore();
    g.fillStyle = col.muted; g.textAlign = "center"; g.textBaseline = "bottom";
    if (o.xl) g.fillText(o.xl, m.l + (w - m.l - m.r) / 2, h - 1);
    if (o.yl) { g.save(); g.translate(10, m.t + (h - m.t - m.b) / 2); g.rotate(-Math.PI / 2); g.textBaseline = "top"; g.fillText(o.yl, 0, -4); g.restore(); }
  }
  function range(arrs, pad) {
    var lo = Infinity, hi = -Infinity;
    arrs.forEach(function (a) { for (var i = 0; i < a.length; i++) { if (a[i] < lo) lo = a[i]; if (a[i] > hi) hi = a[i]; } });
    var d = (hi - lo) * (pad == null ? 0.08 : pad) || 1; return [lo - d, hi + d];
  }
  var STOPS = [[0, 0, 4], [59, 15, 112], [140, 41, 129], [222, 73, 104], [254, 159, 109], [252, 253, 191]], STOPX = [0, 0.25, 0.5, 0.75, 0.9, 1];
  function cmap(t) {
    t = Math.max(0, Math.min(1, t));
    for (var i = 1; i < STOPX.length; i++) if (t <= STOPX[i]) {
      var f = (t - STOPX[i - 1]) / (STOPX[i] - STOPX[i - 1]), a = STOPS[i - 1], b = STOPS[i];
      return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
    }
    return STOPS[STOPS.length - 1];
  }
  function heat(cv, data, circleDiam) {
    var pr = prep(cv); if (!pr) return;
    var g = pr.g, w = pr.w, h = pr.h, col = colors(), ml = 34, mb = 24, mt = 6, mr = 6;
    var side = Math.min(w - ml - mr, h - mt - mb), x0 = ml + (w - ml - mr - side) / 2, y0 = mt;
    var lo = Infinity, hi = -Infinity, i;
    for (i = 0; i < data.length; i++) { if (data[i] < lo) lo = data[i]; if (data[i] > hi) hi = data[i]; }
    var off = document.createElement("canvas"); off.width = GN; off.height = GN;
    var og = off.getContext("2d"), im = og.createImageData(GN, GN);
    for (var y = 0; y < GN; y++) for (var x = 0; x < GN; x++) {
      var c = cmap((data[(GN - 1 - y) * GN + x] - lo) / ((hi - lo) || 1)), k = (y * GN + x) * 4;
      im.data[k] = c[0]; im.data[k + 1] = c[1]; im.data[k + 2] = c[2]; im.data[k + 3] = 255;
    }
    og.putImageData(im, 0, 0);
    g.imageSmoothingEnabled = true; g.drawImage(off, x0, y0, side, side);
    g.strokeStyle = col.grid; g.strokeRect(x0 + 0.5, y0 + 0.5, side - 1, side - 1);
    g.font = "10.5px " + col.mono; g.fillStyle = col.muted;
    g.textAlign = "center"; g.textBaseline = "top";
    [-60, 0, 60].forEach(function (v) { var X = x0 + (v + GEXT) / (2 * GEXT) * side; g.fillText(String(v), X, y0 + side + 4); });
    g.textAlign = "right"; g.textBaseline = "middle";
    [-60, 0, 60].forEach(function (v) { var Y = y0 + side - (v + GEXT) / (2 * GEXT) * side; g.fillText(String(v), x0 - 5, Y); });
    g.textAlign = "center"; g.textBaseline = "bottom"; g.fillText("μas", x0 + side / 2, h);
    if (circleDiam) {
      g.strokeStyle = "rgba(255,255,255,0.85)"; g.lineWidth = 1.2; g.setLineDash([4, 3]); g.beginPath();
      g.arc(x0 + side / 2, y0 + side / 2, circleDiam / 2 / (2 * GEXT) * side, 0, 6.2832); g.stroke(); g.setLineDash([]);
    }
  }

  function render() {
    var A = state.A, B = state.B; if (!A || !B) return;
    var col = colors();
    // 1. señales
    var tx = [A.t[0], A.t[A.t.length - 1]];
    plot($("cv-clean"), { x: tx, y: range([A.c1, A.c2]), series: [
      { type: "line", x: A.t, y: A.c1, color: col.l1 }, { type: "line", x: A.t, y: A.c2, color: col.l2 }], yl: "amplitud", mb: 4 });
    plot($("cv-noisy"), { x: tx, y: range([A.n1, A.n2]), series: [
      { type: "dots", x: A.t, y: A.n1, color: col.l1, r: 1.8, alpha: 0.85 }, { type: "dots", x: A.t, y: A.n2, color: col.l2, r: 1.8, alpha: 0.85 }], xl: "tiempo [s]", yl: "amplitud", mb: 28 });
    // 2. correlación
    var lx = [], ly = [], i;
    for (i = 0; i < A.lags.length; i++) if (Math.abs(A.lags[i]) <= 150) { lx.push(A.lags[i]); ly.push(A.r[i]); }
    var show = $("show256").checked, arr = show ? [ly, A.short] : [ly];
    plot($("cv-corr"), { x: [-150, 150], y: range(arr), series: (show ? [{ type: "line", x: A.shortLags, y: A.short, color: col.muted, lw: 1, alpha: 0.7 }] : []).concat([{ type: "line", x: lx, y: ly, color: col.l1, lw: 1.5 }]),
      vlines: [{ x: A.kt, color: col.ref }], xl: "retardo de prueba k [ms]", yl: "R̂₁₂[k]" });
    // 3. fase
    var fy = A.fx.map(function (f) { return A.slope * f + A.icpt; });
    plot($("cv-phase"), { x: [0, 200], y: range([A.ph, fy], 0.1), series: [
      { type: "dots", x: A.fx, y: A.ph, color: col.l1, r: 2.2, alpha: 0.8 }, { type: "line", x: A.fx, y: fy, color: col.ref, lw: 2 }], xl: "frecuencia en banda base [Hz]", yl: "fase residual [rad]" });
    // 4. uv
    var uvc = $("cv-uv"); uvc.style.height = Math.min(uvc.clientWidth, 420) + "px";
    var ux = B.u.map(function (v) { return v / 1e9; }), vy = B.v.map(function (v) { return v / 1e9; });
    var ax = ux.concat(ux.map(function (v) { return -v; })), ay = vy.concat(vy.map(function (v) { return -v; }));
    plot(uvc, { x: [-10, 10], y: [-10, 10], series: [{ type: "dots", x: ax, y: ay, color: col.l1, r: 1.5 }], xl: "u [Gλ]", yl: "v [Gλ]" });
    // 5. visibilidades
    var xmax = B.bmax > 0 ? Math.max(B.curveX[B.curveX.length - 1], 1) : 10;
    plot($("cv-vis"), { x: [0, xmax], y: [0, 1.05], series: [
      { type: "line", x: B.curveX, y: B.curveY, color: col.ref, lw: 2 }, { type: "dots", x: B.rho, y: B.amp, color: col.l1, r: 2, alpha: 0.55 }],
      xl: "longitud de línea de base [Gλ]", yl: "|V| / V₀" });
    // 6. imágenes
    heat($("cv-model"), B.img, 0); heat($("cv-beam"), B.beam, 0); heat($("cv-dirty"), B.dirty, B.ring);
    readout();
  }
  function kv(l, v, ok) { return '<div class="kv"><span>' + l + '</span><b' + (ok ? ' class="ok"' : '') + '>' + v + '</b></div>'; }
  function readout() {
    var A = state.A, B = state.B, f2 = function (v) { return v.toFixed(2); }, f1 = function (v) { return v.toFixed(1); };
    var h = '<div class="head"><h2>Lectura de la simulación</h2><span class="badge">señales</span></div>';
    h += kv("Razón señal/ruido por muestra", f1(A.snr1) + " dB · " + f1(A.snr2) + " dB");
    h += kv("Retardo real τ_g + ε", f2(A.tt) + " ms");
    h += kv("Pico de correlación (esperado)", f2(A.peak) + " (" + f2(A.expected) + ")");
    h += kv("SNR del pico (fórmula ≈ " + f1(A.theory) + ")", f1(A.snrPeak), true);
    h += kv("Retardo estimado por el pico", f2(A.estPeak) + " ms", true);
    h += kv("Retardo estimado por la fase", f2(A.estPhase) + " ms", true);
    h += '<div class="head"><h2>Cobertura e imagen</h2><span class="badge">plano uv</span></div>';
    if (B.nu === 0) { h += '<div class="kv" style="grid-column:1/-1"><span>Sin líneas de base: activa al menos dos estaciones que vean la fuente a la vez.</span></div>'; }
    else {
      h += kv("Puntos uv (sin conjugados)", String(B.nu));
      h += kv("Línea de base máxima", f2(B.bmax / 1e9) + " Gλ · " + Math.round(B.bmax * B.lam / 1e3).toLocaleString("es-MX") + " km");
      h += kv("Resolución λ / B_max", f1(B.resol) + " μas", true);
      h += kv("Lóbulo lateral máx. del haz", f2(B.side) + " del pico");
      h += kv("Primer mínimo de |V| (anillo delgado)", (B.minPos == null ? "fuera de rango" : f2(B.minPos) + " Gλ") + " (" + f2(B.thin) + ")");
    }
    $("readout").innerHTML = h;
  }

  buildControls();
  compute();
  window.addEventListener("resize", function () { render(); });
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  if (mq.addEventListener) mq.addEventListener("change", render);
  new MutationObserver(render).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
})();
