#set page(paper: "us-letter", margin: 2cm, numbering: "1")
#set text(size: 10pt)
#set par(justify: true)
#set heading(numbering: "1.")
#show link: underline
#let s = json("summary.json")
#let lc = json("lang_check.json")
#let pct(x) = str(calc.round(100 * x, digits: 1)) + "%"
#let num(x, d: 2) = if x == none { "–" } else { str(calc.round(x, digits: d)) }

#align(center)[
  #text(16pt, weight: "bold")[THEMIS-A Time Domain Structures: event dataset and statistics]
  #v(2pt)
  #datetime.today().display() · pipeline #raw(s.pipeline)
]

This draft covers the full list (#s.n_events unique intervals on #s.n_days days, 2010–2022). Dataset: `data/tds_tha.csv`, and `data/tds_tha.parquet` (same typed table with column dictionary; one row per list interval.

Column definitions: @tab:cols. Coverage: FGM #pct(s.frac_fgm), GMOM moments #pct(s.frac_moments) (burst-mode #pct(s.frac_burst_moments)), EFW waveform #pct(s.frac_efw) of events.

= Data and methods <sec:methods>
Three time windows are used, and the column table names which one each column refers to:
- *list interval*: Lang's [`start_time`, `end_time`], a wave-burst interval of typically 6 s that contains the TDS (@sec:lang); $t_c$ is its center. Field, moments and OMNI values are taken over it (moments: spin samples within ±1.5 s, half a spin, of it), and the background field over $t_c$ ± 5 min or ± 30 s.
- *TDS window*: [`tds_start`, `tds_end`], the span of the TDS trains the automatic detector finds inside the list interval (@sec:tdswin).
- *efw window* (`efw_window`): the TDS window where one was found, else the list interval extended by 1 s (as list interval are floored to whole seconds). All `efw_*` statistics use it.

- *Magnetic field*: FGM L2 `fgl` (4 S/s until mid-2015, 16 S/s after), falling back to spin-fit `fgs`. Component medians rather than means (@sec:quality). Background $hat(b)$ is the `fgs` median over ±5 min (window test: @sec:window).
- *Plasma moments*: THEMIS L2 *GMOM*, the ground-computed ESA+SST moments with spacecraft-potential correction. Burst-mode moments (`ptibb`/`ptebb`, full angular resolution, spin cadence) are used when present and consistent with reduced mode (density within ×3), which holds for #pct(s.frac_burst_moments) of events; other events use reduced mode (`ptirf`/`pterf`).
  - Reduced mode can't measure electron anisotropy because $T_(e perp)\/T_(e parallel) approx 1$ always, so `Te_aniso` is valid only for burst events.
  - Sparse full-mode (`ptiff`) moments look interpolated at event times, so they were not used.
- *Waveforms*: L1 `efw` (8192 or 16384 S/s). Spin phase comes from FGM SSL↔DSL. The calibration departs from SPEDAS in three ways:
  - each AC-coupled burst has its mean removed, because the hardware high-pass already removes DC, whereas SPEDAS adds a 2.5 mV/m DSL offset;
  - NaN-safe smoothing;
  - samples near the spin phases where a spin-plane boom points at the Sun (shadow spikes up to ~30 mV/m) are masked: ±7° around 45° and 225°, where the spike splits to ±4.5–6° from 2016 Q3 on, and ±4° around 135° and 315°; also ±4° around an empirically found second pair (165°, 346°).
  All components are then band-passed to 100–3000 Hz (4th-order Butterworth, zero-phase). The band holds ms-scale structures, removes the spin tone and the slow wake/photoelectron excursions that otherwise set the maxima, matches Lang's `max_Ez`, and its fixed upper edge makes 8192 and 16384 S/s bursts comparable. Outputs:
  - signed max $|E_z|$ (DSL spin axis) and signed max $|E_parallel|$, with $hat(b)$ from the window-median DSL field;
  - in-band $|E_perp|$ at the same sample, and the ratio;
  - in-band rms values;
  - the spin phase of each maximum, for artifact screening;
  - $e|E_parallel|lambda_D\/T_e$, which approximates $e phi\/T_e$ for a Debye-scale structure.
- *Derived*: $beta = 2 mu_0 (n_i T_i + n_e T_e)\/B^2$; $lambda_D = sqrt(epsilon_0 T_e\/n_e e^2)$; $f_(p e)\/f_(c e)$; $V_A$; $V_parallel = bold(V) dot hat(b)$ and $bold(V)_perp$; elevation $theta_B = arctan(B_z\/sqrt(B_x^2+B_y^2))$; $Delta B_z$ (mean over 1 min after the event center minus 1 min before); peak ion $V_x$ within ±5 min; $r\/r_"mp"$ against the Shue et al. (1998) magnetopause.
- *Solar wind / indices*: OMNI HRO 1-min: IMF, $V_"sw"$, $P_"dyn"$, AE, AL, SYM-H. Also the prior-1-h IMF $B_z$ mean, AL minimum and AE maximum. Kp from OMNI2 hourly.
- *Region label* (heuristic): sheath/SW if $r > r_"mp"$ or ($n > 5$ cm#super[-3] and $T_i < 1$ keV); dayside/inner if $X > -5 R_E$; otherwise PS ($beta > 0.5$), PSBL ($0.05 < beta < 0.5$) or lobe. The lobe class probably includes outer-PSBL intervals.

= Data quality screening <sec:quality>
Three defects handled in the input products produce extreme values if left in.
- *FGM spikes*: L2 `fgl`/`fgs` contain isolated spikes up to several thousand nT. On some 2022 days near $R approx 4 R_E$ they make up more than 10% of the spin fits in a 10-min window, which turned window means into $Delta B_z approx 5000$ nT and corrupted $hat(b)$. All field averages (event B, background B, $Delta B_z$, and the waveform $hat(b)$) are now component medians.
- *Broken plasma moments* (`mom_mask`):
  - Burst moments without an ESA contribution (SST-only sums with $n tilde 10^(-4)$ cm#super[-3], $T tilde 25$ keV) are rejected when their density disagrees with reduced mode by more than ×3; remaining $n < 10^(-3)$ cm#super[-3] moments are set to NaN (ions #pct(s.mom_mask.ion_lowN), electrons #pct(s.mom_mask.ele_lowN) of events).
  - Some reduced-mode ion moments, mostly in the sheath, carry no angular information: $|bold(V)| approx 10^(-3)$ km/s and $T_perp\/T_parallel equiv 1$. Their velocity and anisotropy are set to NaN (#pct(s.mom_mask.ion_isotropic)). Density and temperature are kept.
  - Single-spin velocity spikes reach several 1000 km/s; `Vx_peak_5min` uses a 3-spin running median.
  - `flag_ion`/`flag_ele` keep the GMOM quality bits (bit 128: $n_i$/$n_e$ mismatch) for further screening.
- *ADC saturation* (`efw_clip`): spin-plane AC channels hit the 16-bit rail near 70 mV/m. −32768 equals the CDF fill value and reads as missing, so saturated samples silently drop out of the maximum. #pct(s.frac_clip) of waveform events contain saturated samples (lobe #pct(s.region.lobe.clip), PS #pct(s.region.PS.clip)). *Their peak amplitudes are lower bounds*.
- *Broken waveform peaks* (`efw_phase_flag`): #pct(s.frac_phase_flag) of waveform events have their peak $|E_parallel|$ or $|E_z|$ at an over-populated spin phase, i.e. a spin-synchronous artifact rather than a TDS (@fig:efwq). Before band-passing and widening the 45°/225° mask this was 32%: broad shadow features set the maxima and masking only moved them to the mask edge. The residue sits just outside the 45° mask. Exclude flagged events from amplitude statistics.
#figure(
  image("figs/efw_quality.png", width: 100%),
  caption: [Spin phase of the efw maxima, all events and those not flagged by `efw_phase_flag`; clusters are spin-synchronous artifacts.],
) <fig:efwq>

= TDS windows (automatic detection) <sec:tdswin>

The automatic detector (`detect/report.typ`) searches every burst for trains of ms-scale, $E_parallel$-dominated, isolated pulses after a 150 Hz high-pass, rejects sun-shadow and spin-recurrent artifacts and wave packets, and labels trains bipolar, unipolar or mixed. `tds_start`/`tds_end` span the detector's TDS trains inside each list interval (extended by 1 s for Lang's floored times). Status (`det_status`): TDS found in #pct(s.det_status.tds) of events, none in #pct(s.det_status.no_tds), no $hat(b)$ #pct(s.det_status.no_b), no searchable burst #pct(s.det_status.no_burst). The median span is #num(s.tds_span_med, d: 1) s, against about 6 s for the list window. The found fraction ranges from #pct(s.tds_frac_by_region.at("sheath/SW")) in the sheath to #pct(s.tds_frac_by_region.lobe) in the lobe.

The detector misses weak (≲ 2 mV/m), sparse TDS by design and occasionally calls long broadband sheath trains TDS; on an audit sample enriched for disagreements it was right more often than the list (`detect/report.pdf`, Validation). Treat `no_tds` events as "not confirmed", not as false. `tds_peak` (150 Hz high-pass, single TDS-like pulses only, lower bound if `tds_clip` > 0) is the detector's own amplitude and agrees with `efw_Epar` (Spearman $rho$ = #num(s.tds_efw_corr)).

= Background-window test <sec:window>
@fig:bgwin varies the half-window $W$ of the background field from 15 s to 20 min (`scripts/bgwindow.jl`). The three diagnostics are:
- the angle between $hat(b)_W$ and the list-interval field;
- the angle between the medians before and after the event center, which is large when the window straddles a boundary or flapping;
- the change in $V_parallel$ relative to ±5 min.

In the magnetosphere (PS, PSBL, lobe, dayside) the field is steady on these scales:
- At ±5 min, $hat(b)$ is within a median 3–8° of the event field.
- Shortening the window to ±1 min changes $V_parallel$ by only a few km/s.
- ±5 min is therefore adequate there. ±1–2 min would follow the local field slightly better without changing any result.

The sheath is different:
- The field rotates on sub-minute scales, so $hat(b)$ at ±5 min is a median #num(s.bgwindow.at("300").ang_event, d: 0)° (all regions) and about 28° (sheath) from the event field.
- For sheath events the list-interval field (or $W lt.eq 30$ s) should be used.

Across all events, #pct(s.frac_bg_unsteady) have more than 30° between the before and after halves at ±5 min (`bg_rot`): PS #pct(s.region.PS.bg_unsteady), PSBL #pct(s.region.PSBL.bg_unsteady), dayside #pct(s.region.at("dayside/inner").bg_unsteady), sheath #pct(s.region.at("sheath/SW").bg_unsteady). For those events $V_parallel$/$V_perp$ depend on the window.

The table therefore carries both windows, and the choice is left to the analysis: `V_par_5min`/`V_perp_5min` (with `Vperp_5min_*`) for the steady magnetosphere, checked with `bg_rot`, and `V_par_30s`/`V_perp_30s` for the sheath and for events with large `bg_rot`. ±30 s is also the window of Lang's B (@sec:lang).
#figure(
  image("figs/bgwindow.png", width: 100%),
  caption: [Median sensitivity to the background half-window $W$, by region.],
) <fig:bgwin>

= Lang's list: definitions and consistency <sec:lang>
The list was produced by `idl/findthmbrst_v2.pro` (interactive selection in SPEDAS). `scripts/lang_check.jl` reproduces its quantities from our data on #lc.n_days random days (#lc.n_events events). Lang notes that some rows may be corrupted during the interactive session, so isolated disagreements with his code are expected; the systematic ones are listed here.
- *Event window*: start/end are a wave-burst interval from `thm_jbt_get_btrange`. Every list window lies inside one continuous burst segment of the L1 file (#pct(lc.window_inside)); #pct(lc.window_exact) coincide with its edges, the rest are consecutive pieces of a longer segment. The sub-window zoomed into for `max_Ez` is not recorded. Our list-interval values are therefore burst averages (typically 6 s, two spins).
- *B* (`lang_B*`): mean `fgs` GSM over list start ±30 s; reproduced to within 1 nT for #pct(lc.B_within1) of events (max #num(lc.B_maxdiff, d: 1) nT). It agrees with our list-interval B (Spearman $rho$ = #num(s.lang_corr_B_z) for $B_z$) up to the window difference.
- *V* (`lang_V*`, not in the table): ESA reduced ions (`peir` through `thm_part_products`, no spacecraft-potential correction), the sample with the largest $|V_x|$ within list start ±30 s. GMOM reduced moments with the same selection reproduce it ($rho$ = #num(lc.rho_V), median difference #num(lc.V_med_absdiff, d: 0) km/s, #lc.V_signflips sign flips in #lc.n_V). Lang's magnitudes are typically smaller (median ratio #num(lc.V_ratio_q.at(1)), 5th percentile #num(lc.V_ratio_q.at(0))), as expected for ESA-only moments. It is a 1-min peak, comparable to `Vx_peak_5min`, not to the event mean `Vi_x`: for 2011-01-02 15:35 UT, #num(lc.example.lang_Vx, d: 0) km/s (Lang) vs #num(lc.example.Vx30, d: 0) (GMOM, same selection) vs #num(lc.example.Vi_x, d: 0) (event mean). In #pct(lc.frac_VxVyVz_plain) of rows of the plain files $V_y$ and $V_z$ are copies of $V_x$; $V_x$ is still valid there. Missing velocity is written as 0.
- *`max_Ez`* is $max|E_parallel|$, not $E_z$: the parallel component from `thm_efi_clean_efw` (response deconvolution, spike removal, high-pass, field-aligned frame from `fgh`), band-passed to 100–3000 Hz, within the zoomed sub-window. Band-passing our $E_parallel$ the same way raises the rank correlation from #num(lc.rho_Ez_raw) to #num(lc.rho_Ez_bp), with median ratio #num(lc.Ez_ratio_bp, d: 1) (ours spans the whole burst, so ≥ 1 is expected; #lc.n_Ez events without `efw_phase_flag`). The rest of the scatter comes from the unrecorded sub-window, the deconvolution and the spike removal. The `*-QLP.txt` files (5 days) re-process existing days and replace the plain files; their `max_Ez` differs, consistent with re-picked sub-windows.
- *`rms_Ew`*: FBK L2 `fb_edc12` (`fb_edc34` if absent), bands 36–572 Hz, $sqrt(3 sum_k pi^2 A_k^2 \/ 8)$ ($pi^2\/8$ converts mean $|$amplitude$|$ to rms for a sinusoid; the factor 3 assumes isotropy), maximum over list start ±30 s. Reproduced to the printed precision for #pct(lc.rms_edc12_match) of events on days with `fb_edc12` (#lc.n_rms_edc12). On days without it (#lc.n_rms_edc34 events), #pct(lc.rms_edc34_daymax) carry the *whole-day* maximum instead (e.g. 2646 mV/m for every event of 2018-09-15). In the full list one value covers all events of the day for #pct(lc.rms_const_by_year.at("2018")) of events in 2018, #pct(lc.rms_const_by_year.at("2022")) in 2022 and #pct(lc.rms_const_by_year.at("2016")) in 2016 (#pct(lc.rms_const_all) overall). Do not use `rms_Ew` on those days.
#figure(
  image("figs/lang_compare.png", width: 80%),
  caption: [Lang's B and `max_Ez` vs our list-interval $B_z$ and band-passed $|E_parallel|$ (`efw_Epar`).],
) <fig:lang>

= Caveats
- *Boom length vs structure size*: in the sheath $lambda_D approx$ #num(s.debye_sheath, d: 0) m, smaller than the spin-plane booms (~35 m effective) and much smaller than their ~50 m tip-to-tip span. Sheath amplitudes there are therefore lower limits. The axial $E_z$ (2.8 m effective) is uncertain to about a factor of 2.
- *AC coupling*: EFW hardware high-pass suppresses unipolar (double-layer) signatures relative to bipolar electron holes.
- *Day boundaries*: only the event's start-day files are read, so windows crossing midnight are truncated.


= Column definitions
#let cols = csv("../data/columns.csv")
#show figure: set block(breakable: true)
#set par(justify: false)
#figure(
  table(
    columns: (28%, auto, auto, 1fr),
    stroke: 0.3pt,
    inset: 3pt,
    align: left,
    fill: (_, y) => if y > 0 and cols.at(y).at(2) == "Lang" { luma(225) },
    table.header(..cols.first().map(x => text(8pt, weight: "bold", x))),
    ..cols.slice(1).flatten().map(x => text(8pt, x)),
  ),
  caption: [Columns of `data/tds_tha.csv`, grouped. Shaded (source "Lang"): copied from Lang's list; all others are computed here. "Interval" is the list interval and $t_c$ its center (@sec:methods).],
) <tab:cols>
