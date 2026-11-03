#set page(paper: "us-letter", margin: 2cm, numbering: "1")
#set text(size: 10pt)
#set par(justify: true)
#set heading(numbering: "1.")
#show link: underline
#show figure.caption: set text(size: 9pt)
#let s = json("out/summary.json")
#let pct(x, d: 0) = str(calc.round(100 * x, digits: d)) + "%"
#let num(x, d: 2) = str(calc.round(x, digits: d))
#let sw(amin, nmin) = s.sweep.find(r => r.amin == amin and r.nmin == nmin)
#let P = s.params
#let n(m) = m.values().sum()
#let frac(m, g) = [#m.at(g + "/det") / #(m.at(g + "/det") + m.at(g + "/nodet"))]
#let prec(m) = m.at("TDS/det") / (m.at("TDS/det") + m.at("neg/det"))
#let rec(m) = m.at("TDS/det") / (m.at("TDS/det") + m.at("TDS/nodet"))
#let lprec(m) = m.at("TDS/lang") / (m.at("TDS/lang") + m.at("neg/lang"))
#let lrec(m) = m.at("TDS/lang") / (m.at("TDS/lang") + m.at("TDS/nolang"))
#let ntds = s.class_trains.tds_bipolar + s.class_trains.tds_unipolar + s.class_trains.tds_mixed

#align(center)[
  #text(16pt, weight: "bold")[Automatic detection of Time Domain Structures in THEMIS-A EFW bursts]
  #v(2pt)
  #datetime.today().display()
]

= Aim

Build a detector that selects Time Domain Structures (TDS) in THEMIS-A `efw` burst waveforms by their physical signature, separates them from waves and instrumental artifacts, and labels their type. Mozer et al. (2015, GRL 42, 3627) describe TDS as trains of ms-scale electric-field structures, mostly parallel to $bold(B)$: electron holes (bipolar $E_parallel$), double layers (unipolar $E_parallel$, net potential drop) and non-linear spiky waves. Lang Qin's hand-picked event list is a reference, not ground truth; disagreements with it are adjudicated by eye (@sec:val).

= Data

Search space: every `efw` burst (≥ 4096 samples, split at 0.5 s gaps; 8192 or 16384 S/s) on the #s.n_days days that carry at least one Lang event: #s.n_bursts bursts, #num(s.burst_hours, d: 0) h of waveform; #s.n_lang_bursts bursts overlap a Lang window (his times are 1 s floors, so windows are extended by 1 s). Calibration to DSL mV/m uses `efi_dsl` from `scripts/efw.jl`. Spin phase $phi$ comes from FGM SSL↔DSL, and $hat(bold(b))$ from FGM `fgh` (128 S/s; `fgl` fallback) interpolated to each sample. #s.nob_bursts bursts have no $hat(bold(b))$. Non-event days were not downloaded. Each burst carries the region of the majority of its day's Lang events.

= Method <sec:method>

Code: `detect/detect.jl` (detector; thresholds in `PARAMS`), `detect/run.jl` (pipeline, `summary.json`, figures), `detect/plots.jl`.

*Pre-processing.* Each component is high-passed by subtracting a twice-applied #num(1000 * P.at("τ"), d: 0) ms boxcar (cut-off ≈ 150 Hz), per contiguous chunk, because bursts contain ≈30 ms packet gaps whose edges ring; edge samples are not searched. This removes DC, spin tone and slow wake/photoelectron excursions, and keeps ms pulses. Detection runs on $E_parallel = bold(E) dot hat(bold(b))$ and on $E_z$ (spin axis, 5.6 m boom pair).

*Noise floor from quiet gaps.* $sigma_q$ is the #num(P.qq) quantile of the rms in #num(1000 * P.qblk, d: 0) ms blocks. In a TDS burst it is set by the quiet intervals between pulses. In a continuous wave it is the wave amplitude itself, so wave crests rarely reach the threshold.

*Pulses.* Local maxima of $|E|$ (#num(1000 * P.nms, d: 1) ms non-maximum suppression) with $|E| >= max(#P.snr sigma_q, #P.amin "mV/m")$, half-max lobe width between #P.wmin samples and #num(1000 * P.wmax, d: 1) ms, and peak / rms within ±#num(1000 * P.isow, d: 0) ms (pulse excluded) ≥ #P.iso. Each pulse gets:

- *shape fit*: least-squares $c + a G + b G'$ with $G$ a unit-peak Gaussian (unipolar, double-layer-like) and $G'$ its unit-peak derivative (bipolar, electron-hole-like), on a grid of widths and centres. The window is ±max(10 widths, #num(1000 * P.fitw, d: 0) ms), so the flanks count: an isolated pulse leaves a residual fraction $R approx 0$, while a wave packet leaves $R gt.tilde 0.5$. The shape is _bipolar_ if $|b| >= 2|a|$, _unipolar_ if $|a| >= 2|b|$, otherwise _mixed_;
- *intermittency*: local kurtosis over ±#num(1000 * P.kw, d: 0) ms (Gaussian noise 3, sinusoid 1.5, isolated pulses ≫ 3), and contrast = $|E|$ / (1.4826 × median $|E|$) within ±#num(1000 * P.isow, d: 0) ms excluding the pulse (Gaussian extremes reach ≈ 3–4);
- *in-band polarization*: the fitted template is regressed on all three filtered components, which gives $|E_parallel|\/|E_perp|$ in the pulse's own band;
- *periodicity*: maximum normalized autocorrelation `ac` within ±#num(1000 * P.isow, d: 0) ms at lags beyond 2 widths (≈ 1 for a quasi-monochromatic packet).

A pulse is *TDS-like* if it is not an artifact, has no clipped sample within ±#num(1000 * P.satw, d: 0) ms, has $|E_parallel|\/|E_perp| >= #P.pol$, and passes at least #P.votes of 3 tests: kurtosis ≥ #P.kurt, contrast ≥ #P.ctr, $R <= #P.rfit$.

*Trains and classes.* Pulses ≤ #num(1000 * P.gap, d: 0) ms apart form a train. Its class is chosen in this order:

+ `artifact` if ≥ #pct(P.art) of its pulses are artifacts;
+ `tds_bipolar` / `tds_unipolar` / `tds_mixed` if it has ≥ #P.nmin TDS-like pulses, by the shape of ≥ 2/3 of them (mixed otherwise);
+ `no_b` if $hat(bold(b))$ is unknown;
+ `wave_coherent` if the median `ac` ≥ #P.acmax (whistler / ion-acoustic packets);
+ `wave_broadband` otherwise.

A burst is detector-positive if it holds at least one `tds_*` train. Trains also record kurtosis, crest factor, duty cycle above $3 sigma_q$, width spread, $f_"ce"$ and clipped fraction (`detect/out/detections.csv`); the class does not use these.

*Physical caveat: boom response.* A structure of scale $L$ comparable to the 50 m spin-plane booms is attenuated there but not on the 5.6 m axial booms, which biases $|E_parallel|\/|E_perp|$ by the orientation of $hat(bold(b))$. $L = v times "width"$ needs the unknown structure speed. For $v gt.tilde 300$ km/s and widths ≳ 0.3 ms, $L gt.tilde 100$ m and the spin-plane response is ≳ 0.8, so no correction is applied.

#figure(image("figs/features.png", width: 88%), caption: [Pulse features (non-artifact pulses) by the class of their train; red: thresholds. The classes are built from these features, so the separation is partly by construction. What the figure shows is _which_ features separate which class: intermittency (kurtosis, contrast) separates TDS from broadband waves, which share TDS-like residuals and $E_parallel$ dominance; shape residual and $E_perp$ dominance separate coherent packets.]) <fig:feat>

= Artifacts <sec:art>

*Sun shadow.* Before any cut, pulse phases pile up at $phi = 45° + k 90°$ (@fig:phase), where a spin-plane boom is shadowed by the spacecraft. In $E_z$ and $E_parallel$ these are clean few-ms bipolar pulses of up to tens of mV/m, morphologically identical to electron holes (@fig:exa row 3, @fig:exb row 4). Three cuts remove them:

+ *spin recurrence*: a pulse recurs at the same phase (±1°) in at least half of the other spins of the burst, well above the chance count of a dense train (#s.cuts.rec pulses);
+ *single-pair shadow*: the pulse lies within #P.shadoww° of 45° mod 90° and ≥ #pct(P.pair) of its spin-plane template energy sits on one boom pair (#s.cuts.shadow1 pulses);
+ *fixed phase mask* (#pct(s.frac_masked_phase, d: 1) of the spin, red in @fig:phase) for the residual excess.

The shadow geometry changes around 2016 Q3: the 45°/225° peaks split to ±4.5°, so the mask is the union of both epochs. After the cuts, TDS-like pulses are flat in phase. In total #s.cuts.art of #s.n_pulses pulses are flagged as artifacts.

*ADC clipping.* Raw count −32768 is both the CDF fill value and negative full scale, so it reads as NaN; +32767 stays finite. Both are flagged as rail samples. NaN samples are interpolated for filtering, and pulses within ±#num(1000 * P.satw, d: 0) ms of a rail sample (#s.cuts.sat) are never TDS-like. Clipping is common: #pct(s.burst_clip.all) of all bursts, #pct(s.burst_clip.lang) of Lang bursts and #pct(s.burst_clip.det) of detector-positive bursts contain rail samples. Within ±#num(1000 * P.isow, d: 0) ms of a TDS train, #pct(s.tds_clip_frac) of trains have at least one rail sample, and #pct(s.tds_clip_frac_1pct) have more than 1%. Their peak amplitudes (`peak` in `detections.csv`) are lower bounds, and `clip` gives the fraction.

#figure(image("figs/phase.png", width: 88%), caption: [Spin phase of all pulses (grey), after the recurrence and single-pair shadow cuts (blue), and TDS-like pulses (black), before/after 2016-07. Red: fixed mask.]) <fig:phase>

= Validation <sec:val>

*Audit sample.* 90 bursts were labelled by eye from overviews of $E_parallel$ and $E_z$ plus ±20 ms zooms of all components (`detect/audit.csv`). Labels are tds, weak (TDS-like but ≲ 2–3 mV/m), wave, ambiguous (bursty activity where pulse vs wave is undecidable), artifact (sun-shadow / spin-plane only) and none. They are grouped as TDS = tds + weak, amb = ambiguous + none, neg = wave + artifact. The sample is stratified, not random:

- 25 bursts flagged by an earlier spike-only version of the detector but absent from Lang's list (_detector-only_);
- 25 Lang bursts that version missed (_Lang-only_);
- 10 bursts in neither (_neither_);
- 30 random sheath/SW bursts, 15 in Lang's list and 15 not (_sheath_), the hardest region for wave/TDS separation.

*Blind* labels were made without the current detector's verdict. *Adjudicated* labels re-examined every disagreement with ±12 ms pulse zooms; this moved #(s.audit_adj_shift.toward + s.audit_adj_shift.away) labels, #s.audit_adj_shift.toward towards the detector's verdict and #s.audit_adj_shift.away away. Adjudication therefore had sight of the detector and is optimistic; the blind numbers are the conservative estimate.

#let row(name, b, a) = (name, [#n(b)], frac(b, "TDS"), frac(b, "amb"), frac(b, "neg"), frac(a, "TDS"), frac(a, "amb"), frac(a, "neg"))
#figure(table(columns: 8, align: (left,) + (center,) * 7, stroke: none,
  table.hline(),
  table.header([], [], table.cell(colspan: 3)[blind], table.cell(colspan: 3)[adjudicated]),
  table.hline(start: 2, end: 5, stroke: 0.5pt), table.hline(start: 5, end: 8, stroke: 0.5pt),
  table.header([sample], [n], [TDS], [amb], [neg], [TDS], [amb], [neg]),
  table.hline(stroke: 0.5pt),
  ..row([detector-only], s.audit_blind.det, s.audit_adj.det),
  ..row([Lang-only], s.audit_blind.lang, s.audit_adj.lang),
  ..row([neither], s.audit_blind.tn, s.audit_adj.tn),
  ..row([sheath], s.audit_blind.sheath, s.audit_adj.sheath),
  table.hline(stroke: 0.5pt),
  ..row([*all*], s.audit_blind_all, s.audit_adj_all),
  table.hline(),
), caption: [Detector vs eye: bursts the detector calls TDS / bursts with that eye label.]) <tab:conf>

Ignoring amb, the detector reaches precision #pct(prec(s.audit_blind_all)) and recall #pct(rec(s.audit_blind_all)) against the blind labels, and #pct(prec(s.audit_adj_all)) and #pct(rec(s.audit_adj_all)) against the adjudicated ones (@tab:conf). Lang's list scores #pct(lprec(s.audit_lang_blind)) / #pct(lrec(s.audit_lang_blind)) (blind) and #pct(lprec(s.audit_lang_adj)) / #pct(lrec(s.audit_lang_adj)) (adjudicated) on the same bursts. These are rates on a sample enriched for disagreement, so they say who is right where the two differ, not population rates. On ambiguous bursts the detector splits evenly (#s.audit_adj_all.at("amb/det") positive of #(s.audit_adj_all.at("amb/det") + s.audit_adj_all.at("amb/nodet")) adjudicated).

#let J(g, k) = [#s.audit_joint_adj.at(g).at(k) (#s.audit_joint_blind.at(g).at(k))]
#figure(table(columns: 5, align: (left,) + (center,) * 4, stroke: none,
  table.hline(),
  table.header([eye label], [both], [detector only], [Lang only], [neither]),
  table.hline(stroke: 0.5pt),
  [TDS], J("TDS", "both"), J("TDS", "det"), J("TDS", "lang"), J("TDS", "none"),
  [amb], J("amb", "both"), J("amb", "det"), J("amb", "lang"), J("amb", "none"),
  [neg], J("neg", "both"), J("neg", "det"), J("neg", "lang"), J("neg", "none"),
  table.hline(),
), caption: [Per-burst agreement of detector and Lang with the adjudicated (blind) eye label, all 90 audited bursts.]) <tab:joint>

*Where Lang fails* (@tab:joint). Of the #(s.audit_joint_adj.TDS.det + s.audit_joint_adj.TDS.none) adjudicated TDS bursts missing from his list, the detector finds #s.audit_joint_adj.TDS.det. They include clear tens-of-mV/m $E_parallel$ bipolar trains (@fig:exa row 2). Of the #(s.audit_joint_adj.neg.lang + s.audit_joint_adj.neg.both) adjudicated-negative bursts he lists, #s.audit_lang_fp_eye.artifact contain only sun-shadow pulses (@fig:exb row 4) and #(s.audit_lang_fp_eye.at("wave", default: 0) + s.audit_lang_fp_eye.at("wave_coherent", default: 0)) are wave packets (@fig:exb row 2). His selection is also inconsistent: on 2020-09-08 the consecutive bursts \#3 and \#4 carry identical shadow pulses and nothing else, and only \#4 is listed (@fig:exa row 3).

*Where the detector fails.* The adjudicated errors are listed in @tab:err.

- *False positives* are long broadband trains (≈ 150–230 pulses) in which a handful of pulses pass the TDS-like votes. The `n_tds` ≥ #P.nmin rule is absolute, not a fraction of the train.
- *Misses* are weak (≈ 1–2 mV/m), sparse TDS with fewer than #P.nmin TDS-like pulses per train. Some of these trains are also dominated by shadow pulses and labelled `artifact`.
- *Undecided*: bursts without $hat(bold(b))$ (`no_b`) are never positive. Both such audited bursts are adjudicated ambiguous.

#figure(table(columns: 6, align: (left, left, left, center, center, left), stroke: none,
  table.hline(),
  table.header([burst], [sample], [eye (blind → adj.)], [Lang], [detector], [detector trains]),
  table.hline(stroke: 0.5pt),
  ..s.audit_det_errors.map(r => ([#r.date \##r.burst], r.set, [#r.class → #r.adjudicated], if r.lang [yes] else [no], if r.det [TDS] else [no], text(8pt, r.classes))).flatten(),
  table.hline(),
), caption: [Audited bursts where the detector contradicts the adjudicated label.]) <tab:err>

*Agreement with Lang over all bursts.* The detector and Lang agree on #s.all.tp bursts. The detector adds #s.all.fp bursts and misses #s.all.fn of his, so recall is #pct(s.all.recall) and precision #pct(s.all.precision) _against Lang_. #s.n_events_hit of #s.n_events Lang events (#pct(s.n_events_hit / s.n_events)) overlap a detected TDS train. By region, recall ranges from #pct(s.region.at("sheath/SW").recall) in the sheath to #pct(s.region.lobe.recall) in the lobe (@fig:agree). Raising the amplitude floor or the train size trades recall for precision but never lifts precision above ≈ #pct(sw(2.0, 10).precision). Given @tab:joint, the residual disagreement is dominated by TDS Lang skipped and by his artifact and wave bursts, not by detector noise.

#figure(image("figs/agreement.png", width: 95%), caption: [Left: burst-level precision/recall vs Lang when re-thresholding the TDS-like pulses by amplitude (labels) and number per train. Right: default detector per region.]) <fig:agree>

#figure(image("figs/examples_a.png", width: 100%), caption: [Burst overviews (high-passed $E_parallel$, $E_z$; Lang windows shaded; red ticks TDS-like pulses, grey other pulses) and ±20 ms zoom at the largest TDS-like pulse (or the largest unmasked $|E|$). Title: region, Lang, detector verdict and train classes.]) <fig:exa>

#figure(image("figs/examples_b.png", width: 100%), caption: [Harder cases. Row 1: sheath broadband bursts with large isolated pulses, detected as TDS (an ambiguous case). Row 2: coherent kHz packet, rejected as wave. Row 3: sheath TDS on a broadband background. Row 4: spin-plane shadow pulses, rejected.]) <fig:exb>

= Class statistics

#let cls = ("tds_bipolar", "tds_unipolar", "tds_mixed", "wave_broadband", "wave_coherent", "artifact", "no_b")
#let regs = ("PS", "PSBL", "lobe", "sheath/SW", "dayside/inner")
#figure(table(columns: 9, align: (left,) + (right,) * 8, stroke: none,
  table.hline(),
  table.header([class], [trains], [in Lang], [clipped], ..regs.map(r => [#r])),
  table.hline(stroke: 0.5pt),
  ..cls.map(c => (raw(c), [#s.class_trains.at(c)], pct(s.class_trains_lang.at(c)), pct(s.class_clip.at(c)),
    ..regs.map(r => [#s.class_region.at(r).at(c, default: 0)]))).flatten(),
  table.hline(),
), caption: [Trains by class: total, fraction lying in Lang bursts, fraction with rail samples, and counts per region.]) <tab:class>

Of #s.n_trains trains, #ntds (#pct(ntds / s.n_trains)) are TDS, which puts #s.n_det_bursts bursts in the detector-positive set (@tab:class). Most trains are `artifact`: shadow pulses form one-pulse trains every quarter spin. Bipolar trains outnumber unipolar ones about #num(s.class_trains.tds_bipolar / s.class_trains.tds_unipolar, d: 0) : 1, consistent with electron holes being the dominant TDS. `tds_mixed` trains, which hold both shapes or pulses fitting neither template cleanly, are nearly as common as bipolar ones. The sheath carries the most TDS trains and also the most broadband-wave trains. Lang bursts make up #pct(s.n_lang_bursts / s.n_bursts) of all bursts. Relative to that baseline, TDS trains are over-represented in them and coherent waves and artifacts under-represented, while #pct(s.class_trains_lang.wave_broadband) of broadband-wave trains also fall in his bursts.

TDS-like pulses (#s.n_tdslike) have a median half-max lobe width of #num(s.pulse_w_med) ms (4 samples; the sampling quantum is 0.12 ms). #pct(s.pulse_bip) fit bipolar and #pct(s.pulse_uni) unipolar. Filtered peak amplitudes have median #num(s.pulse_amp_q.at(0), d: 1), 90th percentile #num(s.pulse_amp_q.at(1), d: 1) and 99th percentile #num(s.pulse_amp_q.at(2), d: 0) mV/m. TDS trains last a median #num(s.train_dur_med) s with #num(s.train_n_med, d: 0) TDS-like pulses. Sheath pulses have the heaviest amplitude tail and the strongest $E_parallel$ dominance (@fig:pulses). Amplitudes come after a 150 Hz high-pass, without boom-response correction, and are clipped in some trains, so they are lower bounds for broad or large pulses.

#figure(image("figs/pulses.png", width: 100%), caption: [TDS-like pulses by region: filtered peak amplitude, half-max lobe width, in-band $|E_parallel|\/|E_perp|$.]) <fig:pulses>

= Limitations

- *Validation.* There is one labeller, the samples are small and stratified, and adjudication was done with sight of the detector. No population precision has been measured for the current detector.
- *Broadband waves.* Long sheath/dayside broadband trains can still be classified TDS through a few TDS-like pulses. The intermittency criteria act per pulse, not per train.
- *Weak, sparse TDS* (≲ 2 mV/m, fewer than #P.nmin per 0.25 s train) are missed by design.
- *Blind spots*: the #pct(s.frac_masked_phase, d: 1) phase mask, ±#num(1000 * P.satw, d: 0) ms around rail samples, and #s.nob_bursts bursts without $hat(bold(b))$.
- *Shape types.* Bipolar vs unipolar is a template ratio. Separating double layers from electron holes properly needs the net potential, hence the structure speed, which is not estimated here. The same unknown speed prevents a boom-response correction of amplitudes and of $|E_parallel|\/|E_perp|$.
- *Occurrence rates.* Only Lang event days are searched.

= Next steps

+ Audit a fresh random sample of the current disagreements (detector-only and Lang-only, stratified by region and class) blind to the detector, ideally with Lang as a second labeller.
+ Add train-level intermittency (fraction of TDS-like pulses, duty cycle, width spread; already in `detections.csv`) to suppress broadband false positives. Then fit thresholds or a small classifier on the audited labels.
+ Use SCM $delta bold(B)$ to separate electromagnetic (whistler) packets from electrostatic structures in the wave classes.
+ Estimate structure speed (e.g. timing between probes) to obtain scale, potential and boom-response-corrected amplitudes.
+ Download non-event days for occurrence rates per region and dwell time.

= Files

`julia --project=. -t 8 detect/run.jl` runs the detection when `detect/out/detections.csv` is absent, then rewrites `summary.json` (every number in this report) and `detect/figs/`.

- `detect/out/bursts.csv`: per burst (Lang overlap, region, verdict, TDS train/pulse counts, clipped fraction).
- `detect/out/detections.csv`: per train (times, pulse counts, peaks, features, class).
- `detect/out/pulses.csv` (≈ 0.85 GB): per pulse (features, artifact flags, shape, TDS-like flag).
- `detect/audit.csv`: eye labels (`class` blind, `adjudicated`).
