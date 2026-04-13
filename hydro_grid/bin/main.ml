(** Electrical grid cascading-failure simulator.

    A 2D grid of generator nodes connected to a single shared bus (one
    system frequency).  The goal is to reproduce the macro dynamics of
    real grid collapses: rising demand pushes generators toward their
    limits, a weak unit trips, its load spills onto neighbors, and the
    failure propagates spatially.

    === Physics modeled ===

    Frequency.  A single system frequency f(t) evolves via the swing
    equation, normalized by capacity-weighted aggregate inertia H_sys.
    Inertia resists rapid frequency change -- grids with more spinning
    mass (high H) tolerate larger shocks before frequency deviates
    dangerously.

    Load damping (beta).  Real loads are frequency-sensitive: motors slow
    when frequency drops, naturally reducing demand.  Modeled as a linear
    term (1.5% of load per Hz) subtracted from the generation-load
    imbalance in the swing equation.  At 0.4 Hz below nominal this
    reduces demand by ~0.6%, partially offsetting any generation loss
    and slowing cascade propagation.

    Primary control (droop).  Each generator adjusts output proportionally
    to frequency error: if frequency drops, output increases.  5% droop
    means a 5%-of-nominal frequency deviation calls for 100% output change.
    This is fast (seconds) but leaves a steady-state frequency offset -- it
    arrests deviations, it does not eliminate them.

    Secondary control (AGC).  A slow integral controller (~120s time
    constant) that nudges each unit's reference setpoint p_ref to drive
    frequency back to nominal.  Real grids rely on AGC to track slow
    load changes (the daily curve).  Without it, droop alone cannot
    absorb a 20% trough-to-peak demand swing within the trip band --
    the required steady-state frequency offset would exceed the 0.5 Hz
    protection threshold.

    Governor dynamics.  Generators cannot jump to a new output instantly.
    A first-order lag (tau_gov ~11s) models the physical inertia of
    turbine governors (valve/gate travel).

    Trip protection.  Three independent triggers per unit:
    - Frequency band:  asymmetric under/over-frequency relay with
      inverse-time characteristic.  Under-frequency uses an
      accumulating dwell timer: each tick the frequency is below a
      node's df_fail threshold, the timer accumulates at a rate that
      increases with severity (1 + k * excess).  The unit trips when
      the accumulated dwell reaches the base limit (300s).  This models
      real ride-through requirements (IEEE C37.117): mild dips (just
      below threshold) take minutes to trip, severe excursions (2+ Hz
      below nominal) trip in seconds.  Per-node df_fail spread means
      each unit's timer starts at a different frequency, naturally
      staggering trips across many ticks.  The timer resets when
      frequency recovers above the node's threshold.
      Over-frequency uses the same inverse-time approach with its own
      dwell accumulator, threshold (1.5 Hz above nominal), and milder
      severity scaling (k=30 vs 60).  Generators tolerate brief
      over-frequency excursions (e.g. UFLS overshoot) without tripping;
      only sustained over-speed causes a trip.  Without this, UFLS
      load shedding that successfully arrests an under-frequency
      cascade can overshoot and kill the grid from above.
    - RoCoF:  rate of change of frequency exceeds threshold (inertia-
      based relay, protects against islanding / loss of synchronism).
    - Overload dwell:  unit is asked to exceed p_max for longer than
      tover_s (thermal / mechanical protection).
    Each node's trip thresholds (df_fail, rho, tover_s) are randomized
    +/- trip_spread around the nominal values.  In real grids, relay
    settings vary across units due to different equipment, calibration,
    and protection philosophies.  This spread means units trip
    individually rather than all at once, producing realistic gradual
    cascades where H_sys erodes as units drop off one by one.

    Cascading.  When a unit trips, its last output is redistributed as
    "proximal load" ℓ to its Von Neumann (N/S/E/W) neighbors.  This is
    local, not system-wide -- in a real grid, transmission constraints
    mean nearby units absorb most of the lost generation.  ℓ decays
    exponentially (tau_ell ~2h) and feeds into the droop setpoint, so
    neighbors ramp up to compensate.  If they hit their own limits, they
    trip too, and the failure spreads spatially.

    UFLS (under-frequency load shedding).  The last automated defense
    before total blackout.  Six stages shed load at progressively lower
    thresholds, defined as offsets from nominal frequency (f0-0.3,
    f0-0.5, f0-0.7, f0-0.9, f0-1.1, f0-1.3 Hz -- 36% total shed).
    For Portugal (f0=50 Hz): 49.7, 49.5, 49.3, 49.1, 48.9, 48.7 Hz.
    The first stage fires early to give the system maximum runway
    before exhausting UFLS.  Some cascades stabilize after one or two
    stages, others burn through all six and collapse.  After frequency
    recovers above f0-0.2 Hz for 60 seconds, one stage reconnects at
    a time -- modeling the slow, deliberate load restoration that
    operators perform after a disturbance.

    One-tick delay.  Redistributed load from a trip at tick t enters
    neighbors' ℓ buckets at t+1.  This prevents instantaneous positive
    feedback within a single timestep and models the finite propagation
    speed of electrical disturbances through the network.

    Init scaling.  At t=0 (midnight), load is at its daily minimum
    (l_mid - l_amp).  We scale each node's p and p_ref to match this
    load so the system starts balanced.  AGC then gradually raises p_ref
    as load rises toward the noon peak.  Without this scaling, the
    initial surplus would trip every unit on the first tick.

    === Plant types ===

    Mix modeled on an Iberian-style grid (Spain/Portugal) with high
    solar penetration.  Four archetypes assigned randomly:
    - Large baseload (10%): high inertia (8s), high capacity (3 MW),
      high utilization (85%).  Anchors the grid.
    - Medium (35%): moderate inertia (5s), moderate capacity (1.5 MW),
      75% utilization.  Bulk of the fleet.
    - Small peaker (25%): low inertia (2.5s), low capacity (0.6 MW),
      65% utilization.  First synchronous unit to trip under stress.
    - Solar (30%): zero inertia, no frequency response (r_pu ~ inf),
      0.5 MW capacity, 20% utilization (average capacity factor).
      Contributes generation but no grid stability services.  Inverter-
      based plants have tighter protection: lower df_fail (0.3 Hz vs
      0.5), more sensitive RoCoF (0.5 vs 1.0 Hz/s), and shorter ride-
      through (60s base dwell vs 300s).  This makes solar the first to
      disconnect during frequency disturbances, modeling real anti-
      islanding behavior.  As solar trips, generation drops but H_sys
      is unaffected (solar contributed no inertia), so the frequency
      decline accelerates less than the generation loss suggests --
      until synchronous units start tripping too.

    === Load model ===

    Daily sinusoid: load(t) = l_mid + l_amp * cos(2pi(t - t_peak)/86400).
    Peak at noon, trough at midnight.  l_mid = load_factor * total_p_ref,
    l_amp = amp_fraction * total_p_ref.  Additional step and ramp terms
    are available for stress-testing but default to zero.

    Gaussian noise (load_noise_pu * l_mid per tick) models the stochastic
    component of real demand -- people switching loads, industrial cycles,
    weather.  Without it, AGC tracks the smooth sinusoid almost perfectly
    and gen-load mismatch is unrealistically small.

    load_factor is the primary knob: at 0.90 the system is comfortable,
    at 1.0 peak demand equals reference generation (no headroom), above
    1.0 cascading failures become likely.

    === Geographic profiles ===

    Grid layout is controlled by a profile (default: Random_uniform).
    Random_uniform fills every cell with a randomly-typed plant.

    Portugal uses a 48x22 grid following the country's actual contour.
    Row 0 = north (Minho), row 47 = south (Algarve).  Atlantic coast
    on the left, Spain border on the right.  Geographic features
    modeled: narrow far north, wider around Porto, coastal bulge at
    Peniche/Nazaré, indentation at the Lisbon/Tagus estuary, widest
    in Alentejo, Algarve narrowing and shifting east, Cape São Vicente
    notch at the southwest corner.

    Plant mix has a strong latitude gradient mirroring the real grid:
    north is hydro-heavy (30% large baseload, only 5% solar), central
    is balanced (30% solar), south is solar-dominant (65%), Algarve is
    almost entirely solar (80%).  Solar plants have tighter protection
    (anti-islanding) with wider threshold spread (±40%) modeling the
    heterogeneous inverter fleet -- older installations trip quickly,
    modern fault-ride-through plants last longer.

    Cascades propagate differently than Random_uniform: the narrow
    north can be cut in two by a few trips, the solar-heavy south has
    low inertia making it fragile once synchronous anchors trip, and
    the Lisbon indentation creates a bottleneck between north and
    south.

    === Display ===

    Grid rendered every 1800 ticks (30 simulated minutes), switching to
    every tick once a cascade starts.  Logs every tick to <name>_sim.csv
    alongside a <name>_map.csv capturing the initial grid layout.

    === CLI ===

    ./grid_sim run <name> [tick_sleep]
        Run simulation, save to <name>_sim.csv and <name>_map.csv.
        Default tick_sleep = 0.5s.

    ./grid_sim replay <name> [start_tick] [tick_sleep]
        Replay from saved files.  No physics computed -- reads CSV and
        renders the display.  start_tick = 0 means show from the
        beginning (30-min snapshots then tick-by-tick on cascade).
        Non-zero start_tick skips display until that tick, then shows
        every tick.  Both <name>_sim.csv and <name>_map.csv must exist.
*)

[@@@warning "-32-37-69"]

(* ======================================================================
   Time constants
   ====================================================================== *)

let secs_per_min = 60
let secs_per_hour = 3600
let secs_per_day = 86400
let secs_per_hour_f = 3600.0
let secs_per_day_f = 86400.0

(* ======================================================================
   Emoji / Unicode glyphs
   ====================================================================== *)

let glyph_hydro   = "\xf0\x9f\x92\xa7"  (* 💧 water drop *)
let glyph_medium   = "\xf0\x9f\x8f\xad"  (* 🏭 factory *)
let glyph_small    = "\xf0\x9f\x94\xa7"  (* 🔧 wrench/turbine *)
let glyph_solar    = "\xf0\x9f\x8c\xbb"  (* 🌻 sunflower *)
let glyph_tripped  = "\xf0\x9f\x92\x80"  (* 💀 skull *)
let glyph_diamond  = "\xe2\x97\x86"      (* ◆ filled diamond *)
let glyph_middot   = "\xc2\xb7"          (* · middle dot *)
let glyph_vline    = "\xe2\x94\x82"      (* │ vertical line *)
let glyph_rtee     = "\xe2\x94\xa4"      (* ┤ right tee *)
let glyph_hline    = "\xe2\x94\x80"      (* ─ horizontal line *)
let glyph_tl       = "\xe2\x94\x8c"      (* ┌ top-left corner *)
let glyph_tr       = "\xe2\x94\x90"      (* ┐ top-right corner *)
let glyph_bl       = "\xe2\x94\x94"      (* └ bottom-left corner *)
let glyph_br       = "\xe2\x94\x98"      (* ┘ bottom-right corner *)
let glyph_ltee     = "\xe2\x94\x9c"      (* ├ left tee *)
let glyph_ttee     = "\xe2\x94\xac"      (* ┬ top tee *)
let glyph_btee     = "\xe2\x94\xb4"      (* ┴ bottom tee *)
let glyph_cross    = "\xe2\x94\xbc"      (* ┼ cross *)
let glyph_block_full  = "\xe2\x96\x88"  (* █ full block *)
let glyph_block_light = "\xe2\x96\x91"  (* ░ light shade *)

(* ======================================================================
   Color theme -- Catppuccin Mocha
   256-color ANSI approximations of the Catppuccin Mocha palette.
   Keep this block self-contained so it can move to config later.
   ====================================================================== *)

let ansi_reset      = "\027[0m"
let ansi_dim        = "\027[38;5;243m"  (* Overlay0 #6c7086 *)
let color_hydro     = "\027[38;5;111m"  (* Blue #89b4fa *)
let color_med       = "\027[38;5;147m"  (* Lavender #b4befe *)
let color_small     = "\027[38;5;216m"  (* Peach #fab387 *)
let color_solar     = "\027[38;5;223m"  (* Yellow #f9e2af *)
let color_freq_ok   = "\027[38;5;151m"  (* Green #a6e3a1 *)
let color_freq_warn = "\027[38;5;223m"  (* Yellow #f9e2af *)
let color_freq_crit = "\027[38;5;211m"  (* Red #f38ba8 *)
let color_alert     = "\027[38;5;211m"  (* Red #f38ba8 *)
let color_pause     = "\027[38;5;223m"  (* Yellow #f9e2af *)
let stress_bg_50    = "\027[48;5;23m"   (* Surface tinted green *)
let stress_bg_70    = "\027[48;5;58m"   (* Surface tinted yellow *)
let stress_bg_85    = "\027[48;5;52m"   (* Surface tinted red *)
let stress_bg_100   = "\027[48;5;88m"   (* Maroon dark *)

(* ======================================================================
   Key bindings
   Keep this block self-contained so it can move to config later.
   ====================================================================== *)

let key_quit      = 'q'
let key_pause     = ' '
let key_fast      = 'f'
let key_slow      = 's'
let key_rewind    = 'r'  (* replay only, not yet implemented *)
let key_term_size = 't'

let speed_min = 0.125
let speed_max = 16.0
let rewind_base = 10  (* frames per rewind step, not yet implemented *)

(* ======================================================================
   Configuration constants
   ====================================================================== *)

type plant_config = {
  name : string;
  h_s : float;        (* inertia constant (seconds) *)
  p_max : float;      (* nameplate capacity (MW) *)
  p_ref_ratio : float; (* utilization: p_ref = ratio * p_max *)
  prob : float;        (* probability of placement on grid *)
  r_pu : float;        (* droop in per-unit *)
  tau_gov_s : float;   (* governor time constant (seconds) *)
}

let large_baseload =
  { name = "large"; h_s = 8.0; p_max = 3.0;
    p_ref_ratio = 0.85; prob = 0.10;
    r_pu = 0.05; tau_gov_s = 11.0 }

let medium_plant =
  { name = "medium"; h_s = 5.0; p_max = 1.5;
    p_ref_ratio = 0.75; prob = 0.35;
    r_pu = 0.05; tau_gov_s = 11.0 }

let small_peaker =
  { name = "small"; h_s = 2.5; p_max = 0.6;
    p_ref_ratio = 0.65; prob = 0.25;
    r_pu = 0.05; tau_gov_s = 11.0 }

let solar =
  { name = "solar"; h_s = 0.0; p_max = 0.5;
    p_ref_ratio = 0.20; prob = 0.30;
    r_pu = 1e6; tau_gov_s = 1.0 }

let plant_types = [ large_baseload; medium_plant; small_peaker; solar ]

let default_rows = 36
let default_cols = 40
let default_seed = 42
let default_f0 = 50.0              (* nominal frequency (Hz), overridden by profile_f0 *)
let default_dt_s = 1.0             (* seconds per tick *)
let default_load_factor = 1.15     (* l_mid / total_p_ref *)
let default_amp_fraction = 0.25    (* l_amp / total_p_ref *)
let default_t_peak_s = 12 * secs_per_hour (* peak demand at noon *)
let default_load_step_mw = 0.0     (* sudden load step (MW), unused *)
let default_load_step_t = 0        (* tick when step kicks in *)
let default_load_ramp_mw_per_s = 0.0 (* linear ramp (MW/s), unused *)
let default_beta_pu = 0.015        (* load damping: 1.5% of load per Hz *)
let default_tau_ell_s = 2.0 *. secs_per_hour_f (* proximal load decay (seconds) *)
let default_df_fail = 0.5          (* under-frequency trip (Hz below nominal) *)
let default_df_fail_solar = 0.3    (* solar: tighter anti-islanding threshold *)
let default_df_fail_over = 1.5     (* over-frequency trip (Hz above nominal) *)
let default_freq_dwell_s = 300.0   (* inverse-time base dwell for under-freq trip (seconds) *)
let default_freq_dwell_solar_s = 60.0 (* solar: shorter ride-through *)
let default_freq_dwell_k = 60.0    (* under-freq severity scaling *)
let default_freq_dwell_over_s = 300.0 (* inverse-time base dwell for over-freq trip (seconds) *)
let default_freq_dwell_over_k = 30.0  (* over-freq severity scaling (less aggressive) *)
let default_rho = 1.0              (* RoCoF trip threshold (Hz/s) *)
let default_rho_solar = 0.5        (* solar: more sensitive to RoCoF *)
let default_tover_s = 120.0        (* overload dwell trip (seconds) *)
let default_tau_agc_s = 120.0      (* AGC integration time constant (seconds) *)
let default_load_noise_pu = 0.003  (* load noise std dev as fraction of l_mid *)
let default_trip_spread = 0.20     (* +/- fractional spread on per-node trip thresholds *)
let default_trip_spread_solar = 0.40 (* solar: wider spread, heterogeneous inverter fleet *)

type ufls_offset = { df : float; shed_frac : float }

let ufls_offsets = [|
  { df = -0.3; shed_frac = 0.05 };
  { df = -0.5; shed_frac = 0.05 };
  { df = -0.7; shed_frac = 0.06 };
  { df = -0.9; shed_frac = 0.06 };
  { df = -1.1; shed_frac = 0.07 };
  { df = -1.3; shed_frac = 0.07 };
|]

let ufls_recovery_df = -0.2       (* reconnect when f > f0 + this offset *)
let ufls_recovery_delay_s = 60.0  (* seconds above threshold before reconnecting one stage *)

let default_sim_hours = 24

(* ======================================================================
   Node and grid types
   ====================================================================== *)

type grid_profile = Random_uniform | Portugal

let default_profile = Portugal

let profile_dims = function
  | Random_uniform -> (default_rows, default_cols)
  | Portugal -> (44, 22)

let profile_f0 = function
  | Random_uniform -> 60.0
  | Portugal -> 50.0

let active_f0 = ref default_f0

type status = On | Off | Empty

type node = {
  status : status;
  h_s : float;
  p : float;
  p_ref : float;
  p_max : float;
  r_pu : float;
  tau_gov_s : float;
  ell : float;
  theta_s : float;
  freq_dwell : float;
  freq_dwell_over : float;
  freq_dwell_base : float;
  df_fail : float;
  rho : float;
  tover_s : float;
}

type grid = node array array

type sim = {
  g : grid;

  t : int;
  dt_s : float;

  f : float;
  f_prev : float;
  f0 : float;

  beta : float;

  l_mid : float;
  l_amp : float;
  t_peak_s : int;

  load_step_mw : float;
  load_step_t : int;
  load_ramp_mw_per_s : float;

  s_base : float;

  tau_ell_s : float;

  tau_agc_s : float;

  load_noise_pu : float;
  load_rng : Random.State.t;

  ufls_stages_fired : int;
  ufls_shed : float;
  ufls_recovery_s : float;
}

(* ======================================================================
   Utility
   ====================================================================== *)

let clip (x : float) (lo : float) (hi : float) : float =
  if x < lo then lo else if x > hi then hi else x

let clear () =
  print_string "\027[2J\027[H"  (* clear screen + cursor home -- no constant, terminal control *)

let cursor_home () =
  print_string "\027[H"

let cursor_hide () =
  print_string "\027[?25l"

let cursor_show () =
  print_string "\027[?25h"

let dim (g : grid) : int * int =
  let rows = Array.length g in
  if rows = 0 then (0, 0)
  else (rows, Array.length g.(0))

let in_bounds (g : grid) (r : int) (c : int) : bool =
  let rows, cols = dim g in
  r >= 0 && r < rows && c >= 0 && c < cols

let neighbors4_coords (g : grid) (r : int) (c : int) : (int * int) list =
  [ (r-1,c); (r+1,c); (r,c-1); (r,c+1) ]
  |> List.filter (fun (rr, cc) -> in_bounds g rr cc)

let exp_decay (dt_s : float) (tau_s : float) : float =
  if tau_s <= 0.0 then 0.0 else exp (-. dt_s /. tau_s)

let gaussian (rng : Random.State.t) : float =
  let u1 = Random.State.float rng 1.0 in
  let u2 = Random.State.float rng 1.0 in
  let u1 = Float.max u1 1e-30 in
  sqrt (-2.0 *. log u1) *. cos (2.0 *. Float.pi *. u2)

let saved_term_attr : Unix.terminal_io option ref = ref None

let restore_terminal () =
  print_string "\027[?25h";
  (match !saved_term_attr with
  | Some attr -> Unix.tcsetattr Unix.stdin Unix.TCSANOW attr; saved_term_attr := None
  | None -> ())

let with_raw_mode fn =
  let fd = Unix.stdin in
  let old_attr = Unix.tcgetattr fd in
  saved_term_attr := Some old_attr;
  at_exit restore_terminal;
  let raw = { old_attr with Unix.c_icanon = false; c_echo = false; c_vmin = 1; c_vtime = 0 } in
  Unix.tcsetattr fd Unix.TCSANOW raw;
  Fun.protect ~finally:(fun () -> Unix.tcsetattr fd Unix.TCSANOW old_attr; saved_term_attr := None) fn

let speed_mult = ref 1.0
let show_term_size = ref false
let term_size_cache = ref ""

let read_char_if_ready () =
  let ready, _, _ = Unix.select [Unix.stdin] [] [] 0.0 in
  if ready <> [] then begin
    let b = Bytes.create 1 in
    ignore (Unix.read Unix.stdin b 0 1);
    Some (Bytes.get b 0)
  end else None

let read_char_blocking () =
  let b = Bytes.create 1 in
  ignore (Unix.read Unix.stdin b 0 1);
  Bytes.get b 0

let query_term_size () =
  try
    let ic = Unix.open_process_in "stty size 2>/dev/null" in
    let line = input_line ic in
    ignore (Unix.close_process_in ic);
    match String.split_on_char ' ' line with
    | [h; w] -> term_size_cache := Printf.sprintf "%sx%s" w h
    | _ -> ()
  with _ -> ()

let handle_common_key (c : char) : unit =
  if c = key_quit then (Printf.printf "\n"; exit 0)
  else if c = key_fast then
    speed_mult := Float.min speed_max (!speed_mult *. 2.0)
  else if c = key_slow then
    speed_mult := Float.max speed_min (!speed_mult /. 2.0)
  else if c = key_term_size then begin
    show_term_size := not !show_term_size;
    if !show_term_size then query_term_size ()
  end

let handle_pause () =
  Printf.printf "%sPAUSED (space=resume q=quit f/s=speed t=size)%s%!" color_pause ansi_reset;
  let rec wait () =
    let c = read_char_blocking () in
    if c = key_pause then ()
    else begin handle_common_key c; wait () end
  in
  wait ();
  Printf.printf "\r\027[K%!" (* clear the pause message line *)

(* Sleep for secs, return key pressed during sleep (if any) *)
let interruptible_sleep (secs : float) : char option =
  if secs <= 0.0 then read_char_if_ready ()
  else begin
    let ready, _, _ = Unix.select [Unix.stdin] [] [] secs in
    if ready <> [] then begin
      let c = read_char_blocking () in
      Some c
    end else None
  end

(* ======================================================================
   Grid initialization
   ====================================================================== *)

let pick_plant_type (rng : Random.State.t) : plant_config =
  let r = Random.State.float rng 1.0 in
  let rec walk cum = function
    | [] -> small_peaker
    | [pt] -> ignore cum; pt
    | pt :: rest ->
      let cum' = cum +. pt.prob in
      if r < cum' then pt else walk cum' rest
  in
  walk 0.0 plant_types

let spread (rng : Random.State.t) (nominal : float) (frac : float) : float =
  let lo = nominal *. (1.0 -. frac) in
  let hi = nominal *. (1.0 +. frac) in
  lo +. Random.State.float rng (hi -. lo)

let empty_node = {
  status = Empty; h_s = 0.0; p = 0.0; p_ref = 0.0; p_max = 0.0;
  r_pu = 1e6; tau_gov_s = 1.0; ell = 0.0; theta_s = 0.0;
  freq_dwell = 0.0; freq_dwell_over = 0.0; freq_dwell_base = 1e6;
  df_fail = 1e6; rho = 1e6; tover_s = 1e6;
}

let node_of_plant (rng : Random.State.t) (pc : plant_config) : node =
  let p_ref = pc.p_ref_ratio *. pc.p_max in
  let is_solar = pc.h_s <= 0.0 in
  let df = if is_solar then default_df_fail_solar else default_df_fail in
  let rho_nom = if is_solar then default_rho_solar else default_rho in
  let dwell_base = if is_solar then default_freq_dwell_solar_s else default_freq_dwell_s in
  let sp = if is_solar then default_trip_spread_solar else default_trip_spread in
  { status = On;
    h_s = pc.h_s;
    p = p_ref;
    p_ref;
    p_max = pc.p_max;
    r_pu = pc.r_pu;
    tau_gov_s = pc.tau_gov_s;
    ell = 0.0;
    theta_s = 0.0;
    freq_dwell = 0.0;
    freq_dwell_over = 0.0;
    freq_dwell_base = dwell_base;
    df_fail = spread rng df sp;
    rho = spread rng rho_nom sp;
    tover_s = spread rng default_tover_s sp }

let init_grid_heterogeneous (rows : int) (cols : int) (seed : int) : grid =
  let rng = Random.State.make [| seed |] in
  Array.init rows (fun _ ->
    Array.init cols (fun _ ->
      let pc = pick_plant_type rng in
      node_of_plant rng pc
    )
  )

(* Portugal land mask following the country's contour on a 48x22 grid.
   Row 0 = north (Minho), row 47 = south (Algarve).
   Atlantic coast on the left, Spain border on the right.
   Key geographic features modeled:
   - Narrow north (Minho/Trás-os-Montes)
   - Wider around Porto/Douro
   - Slight coastal bulge at Peniche/Nazaré
   - Indentation at Lisbon/Tagus
   - Widest in Alentejo
   - Algarve narrows and shifts east
   - Cape São Vicente notch at southwest corner *)
let portugal_land_bounds (row : int) (rows : int) (cols : int) : (int * int) option =
  let frac = float_of_int row /. float_of_int (rows - 1) in
  let fcols = float_of_int cols in
  (* West coast offset from col 0: higher = more ocean on left *)
  let coast =
    if frac < 0.08 then 8.0                            (* far north Minho: narrow *)
    else if frac < 0.18 then 8.0 -. (frac -. 0.08) /. 0.10 *. 3.0  (* widens toward Porto *)
    else if frac < 0.28 then 5.0 -. (frac -. 0.18) /. 0.10 *. 1.0  (* Porto/Aveiro *)
    else if frac < 0.35 then 4.0 -. (frac -. 0.28) /. 0.07 *. 2.0  (* Peniche bulge *)
    else if frac < 0.42 then 2.0 +. (frac -. 0.35) /. 0.07 *. 2.0  (* Lisbon indentation *)
    else if frac < 0.52 then 4.0 -. (frac -. 0.42) /. 0.10 *. 1.0  (* Setúbal *)
    else if frac < 0.75 then 3.0 -. (frac -. 0.52) /. 0.23 *. 3.0  (* Alentejo coast: widest *)
    else if frac < 0.88 then 0.0 +. (frac -. 0.75) /. 0.13 *. 2.0  (* SW corner narrows *)
    else 2.0 +. (frac -. 0.88) /. 0.12 *. 6.0                      (* Algarve shifts east *)
  in
  (* East border: mostly at right edge, indented in Algarve *)
  let east =
    if frac < 0.08 then fcols -. 1.0 -. (0.08 -. frac) /. 0.08 *. 4.0 (* NE corner taper *)
    else if frac < 0.88 then fcols -. 1.0                              (* Spain border: full *)
    else fcols -. 1.0 -. (frac -. 0.88) /. 0.12 *. 4.0                (* Algarve east taper *)
  in
  let c0 = max 0 (int_of_float (Float.round coast)) in
  let c1 = min (cols - 1) (int_of_float (Float.round east)) in
  if c1 < c0 then None else Some (c0, c1)

let pick_plant_weighted (rng : Random.State.t)
    (weights : (plant_config * float) list) : plant_config =
  let r = Random.State.float rng 1.0 in
  let rec walk cum = function
    | [] -> small_peaker
    | [(pc, _)] -> ignore cum; pc
    | (pc, w) :: rest ->
      let cum' = cum +. w in
      if r < cum' then pc else walk cum' rest
  in
  walk 0.0 weights

let portugal_plant_weights (row : int) (rows : int)
    : (plant_config * float) list =
  let frac = float_of_int row /. float_of_int (rows - 1) in
  if frac < 0.30 then        (* North: hydro-heavy, minimal solar *)
    [(large_baseload, 0.30); (medium_plant, 0.40);
     (small_peaker, 0.25); (solar, 0.05)]
  else if frac < 0.55 then   (* Central: balanced *)
    [(large_baseload, 0.10); (medium_plant, 0.35);
     (small_peaker, 0.25); (solar, 0.30)]
  else if frac < 0.85 then   (* South: solar-dominant *)
    [(large_baseload, 0.05); (medium_plant, 0.15);
     (small_peaker, 0.15); (solar, 0.65)]
  else                        (* Algarve: almost all solar *)
    [(large_baseload, 0.02); (medium_plant, 0.08);
     (small_peaker, 0.10); (solar, 0.80)]

let init_grid_portugal (rows : int) (cols : int) (seed : int) : grid =
  let rng = Random.State.make [| seed |] in
  Array.init rows (fun r ->
    Array.init cols (fun c ->
      match portugal_land_bounds r rows cols with
      | None -> empty_node
      | Some (c0, c1) ->
        if c < c0 || c > c1 then empty_node
        else
          let weights = portugal_plant_weights r rows in
          let pc = pick_plant_weighted rng weights in
          node_of_plant rng pc
    )
  )

let init_grid (profile : grid_profile) (rows : int) (cols : int) (seed : int) : grid =
  match profile with
  | Random_uniform -> init_grid_heterogeneous rows cols seed
  | Portugal -> init_grid_portugal rows cols seed

(* ======================================================================
   Grid aggregates
   ====================================================================== *)

let gen_total (g : grid) : float =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n ->
      match n.status with On -> acc +. n.p | Off | Empty -> acc
    ) acc row
  ) 0.0 g

let gen_capacity (g : grid) : float =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n ->
      match n.status with On -> acc +. n.p_max | Off | Empty -> acc
    ) acc row
  ) 0.0 g

let s_base_of_grid (g : grid) : float =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n -> acc +. n.p_max) acc row
  ) 0.0 g

let total_p_ref (g : grid) : float =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n ->
      match n.status with On -> acc +. n.p_ref | Off | Empty -> acc
    ) acc row
  ) 0.0 g

let h_sys (g : grid) ~(s_base : float) : float =
  if s_base <= 0.0 then 0.0
  else
    let num =
      Array.fold_left (fun acc row ->
        Array.fold_left (fun acc n ->
          match n.status with
          | On -> acc +. (n.h_s *. n.p_max)
          | Off | Empty -> acc
        ) acc row
      ) 0.0 g
    in
    num /. s_base

let count_on (g : grid) : int =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n ->
      match n.status with On -> acc + 1 | Off | Empty -> acc
    ) acc row
  ) 0 g

let count_plants (g : grid) : int =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n ->
      match n.status with Empty -> acc | On | Off -> acc + 1
    ) acc row
  ) 0 g

let count_on_by_type (g : grid) : int * int * int * int =
  Array.fold_left (fun (hy, md, sm, so) row ->
    Array.fold_left (fun (hy, md, sm, so) n ->
      match n.status with
      | On ->
        if n.h_s <= 0.0 then (hy, md, sm, so + 1)
        else if n.p_max >= 2.5 then (hy + 1, md, sm, so)
        else if n.p_max >= 1.0 then (hy, md + 1, sm, so)
        else (hy, md, sm + 1, so)
      | Off | Empty -> (hy, md, sm, so)
    ) (hy, md, sm, so) row
  ) (0, 0, 0, 0) g

(* ======================================================================
   Load model
   ====================================================================== *)

let load_at (s : sim) : float =
  let period_s = secs_per_day_f in
  let x = 2.0 *. Float.pi *. (float_of_int (s.t - s.t_peak_s)) /. period_s in
  let base = s.l_mid +. s.l_amp *. cos x in

  let step =
    if s.t >= s.load_step_t then s.load_step_mw else 0.0
  in

  let ramp = s.load_ramp_mw_per_s *. float_of_int s.t in

  let noise =
    if s.load_noise_pu <= 0.0 then 0.0
    else s.load_noise_pu *. s.l_mid *. gaussian s.load_rng
  in

  (base +. step +. ramp +. noise) *. (1.0 -. s.ufls_shed)

(* ======================================================================
   Droop setpoint
   ====================================================================== *)

let droop_setpoint_raw (f : float) (f0 : float) (n : node) : float =
  let k_mw_per_hz =
    if n.r_pu <= 0.0 then 0.0 else n.p_ref /. (n.r_pu *. f0)
  in
  n.p_ref +. k_mw_per_hz *. (f0 -. f) +. n.ell

(* ======================================================================
   Scope (oscillograph) -- rolling chart of fleet composition
   ====================================================================== *)

let scope_len = 20
let scope_height = 21

type scope_entry = { s_hy : int; s_md : int; s_sm : int; s_so : int }

let scope_buf = Array.make scope_len { s_hy = 0; s_md = 0; s_sm = 0; s_so = 0 }
let scope_idx = ref 0
let scope_count = ref 0

let orig_total = ref 1



let scope_push (e : scope_entry) =
  scope_buf.(!scope_idx) <- e;
  scope_idx := (!scope_idx + 1) mod scope_len;
  if !scope_count < scope_len then incr scope_count

let scope_get (i : int) : scope_entry =
  let pos = (!scope_idx - !scope_count + i + scope_len * 2) mod scope_len in
  scope_buf.(pos)

let scope_row_for_pct (p : int) : int =
  let clamped = min p 100 in
  scope_height - 1 - (clamped * (scope_height - 1) / 100)

(* Frequency trace *)
let freq_height = 11
let freq_lo () = !active_f0 -. 3.0
let freq_hi () = !active_f0 +. 1.0
let freq_buf = Array.make scope_len 50.0
let freq_idx = ref 0
let freq_count = ref 0

let freq_push (f : float) =
  freq_buf.(!freq_idx) <- f;
  freq_idx := (!freq_idx + 1) mod scope_len;
  if !freq_count < scope_len then incr freq_count

let freq_get (i : int) : float =
  let pos = (!freq_idx - !freq_count + i + scope_len * 2) mod scope_len in
  freq_buf.(pos)

let freq_row_for_val (f : float) : int =
  let lo = freq_lo () in
  let hi = freq_hi () in
  let clamped = max lo (min hi f) in
  let frac = (hi -. clamped) /. (hi -. lo) in
  int_of_float (frac *. float_of_int (freq_height - 1) +. 0.5)

(* Event log *)
let event_log_len = 50
let event_log = Array.make event_log_len ""
let event_log_idx = ref 0
let event_log_count = ref 0

let event_push (msg : string) =
  event_log.(!event_log_idx) <- msg;
  event_log_idx := (!event_log_idx + 1) mod event_log_len;
  if !event_log_count < event_log_len then incr event_log_count

let event_get (i : int) : string =
  let pos = (!event_log_idx - !event_log_count + i + event_log_len * 2) mod event_log_len in
  event_log.(pos)

(* Stress string for CSV *)
let grid_stress_string (g : grid) : string =
  let buf = Buffer.create 1024 in
  Array.iter (fun row ->
    Array.iter (fun n ->
      match n.status with
      | Off | Empty -> Buffer.add_char buf '0'
      | On ->
        let ratio = if n.p_max <= 0.0 then 0.0 else n.p /. n.p_max in
        let d = min 9 (max 1 (int_of_float (ratio *. 9.0 +. 0.5))) in
        Buffer.add_char buf (Char.chr (d + Char.code '0'))
    ) row
  ) g;
  Buffer.contents buf

(* Heatmap: wrap glyph in ANSI 256-color background based on stress digit 0-9 *)
let stress_bg (stress : int) (glyph : string) : string =
  if stress <= 4 then glyph
  else if stress <= 6 then stress_bg_50 ^ glyph ^ ansi_reset
  else if stress <= 7 then stress_bg_70 ^ glyph ^ ansi_reset
  else if stress <= 8 then stress_bg_85 ^ glyph ^ ansi_reset
  else stress_bg_100 ^ glyph ^ ansi_reset

(* ======================================================================
   Display -- framed layout
   ====================================================================== *)

(* Repeat a UTF-8 glyph n times *)
let repeat_glyph (g : string) (n : int) : string =
  let buf = Buffer.create (n * String.length g) in
  for _ = 1 to n do Buffer.add_string buf g done;
  Buffer.contents buf

(* Sidebar panel row offsets (relative to grid row 0) *)
let scope_start = 0
let freq_start_off = scope_height + 1
let event_start_off = freq_start_off + freq_height + 1

(* Build sidebar content string for a given grid row *)
let sidebar_content (r : int) (buf : Buffer.t) : unit =
  let b = buf in
  (* Fleet scope *)
  let cr = r - scope_start in
  if cr >= 0 && cr < scope_height then begin
    let is_quarter =
      cr > 0 && cr < scope_height - 1 &&
      let axis_val = 100 * (scope_height - 1 - cr) / (scope_height - 1) in
      axis_val = 75 || axis_val = 50 || axis_val = 25
    in
    let label = if cr = 0 then "100"
      else if is_quarter then
        Printf.sprintf "%3d" (100 * (scope_height - 1 - cr) / (scope_height - 1))
      else "   " in
    let tick =
      if cr = 0 || cr = scope_height - 1 || is_quarter then glyph_rtee
      else glyph_vline
    in
    Buffer.add_string b (Printf.sprintf "%s%s" label tick);
    let dot = glyph_diamond in
    for i = 0 to !scope_count - 1 do
      let e = scope_get i in
      let r_so = scope_row_for_pct e.s_so in
      let r_hy = scope_row_for_pct e.s_hy in
      let r_md = scope_row_for_pct e.s_md in
      let r_sm = scope_row_for_pct e.s_sm in
      if cr = r_so then
        Buffer.add_string b (color_solar ^ dot ^ ansi_reset)
      else if cr = r_hy then
        Buffer.add_string b (color_hydro ^ dot ^ ansi_reset)
      else if cr = r_md then
        Buffer.add_string b (color_med ^ dot ^ ansi_reset)
      else if cr = r_sm then
        Buffer.add_string b (color_small ^ dot ^ ansi_reset)
      else
        Buffer.add_string b glyph_middot
    done;
    for _ = !scope_count to scope_len - 1 do
      Buffer.add_string b glyph_middot
    done
  end;
  (* Freq trace *)
  let fr = r - freq_start_off in
  if fr >= 0 && fr < freq_height then begin
    let hi = freq_hi () in
    let lo = freq_lo () in
    let f0 = !active_f0 in
    let nominal_row = int_of_float ((hi -. f0) /. (hi -. lo) *. float_of_int (freq_height - 1) +. 0.5) in
    let is_nominal = fr = nominal_row in
    let label =
      if fr = 0 then Printf.sprintf "%4.0f" hi
      else if fr = freq_height - 1 then Printf.sprintf "%4.0f" lo
      else if is_nominal then Printf.sprintf " %2.0f " f0
      else "    "
    in
    let tick =
      if fr = 0 || fr = freq_height - 1 || is_nominal then glyph_rtee
      else glyph_vline
    in
    Buffer.add_string b (Printf.sprintf "%s%s" label tick);
    let dot = glyph_diamond in
    for i = 0 to !freq_count - 1 do
      let f = freq_get i in
      let target_row = freq_row_for_val f in
      if fr = target_row then begin
        if f > f0 -. 0.5 then Buffer.add_string b (color_freq_ok ^ dot ^ ansi_reset)
        else if f > f0 -. 1.0 then Buffer.add_string b (color_freq_warn ^ dot ^ ansi_reset)
        else Buffer.add_string b (color_freq_crit ^ dot ^ ansi_reset)
      end else if is_nominal then
        Buffer.add_string b (ansi_dim ^ "-" ^ ansi_reset)
      else
        Buffer.add_string b glyph_middot
    done;
    for _ = !freq_count to scope_len - 1 do
      if is_nominal then Buffer.add_string b (ansi_dim ^ "-" ^ ansi_reset)
      else Buffer.add_string b glyph_middot
    done
  end;
  (* Event log *)
  let er = r - event_start_off in
  if er >= 0 && er < event_log_len then begin
    if er < !event_log_count then
      Buffer.add_string b (event_get er)
  end

(* Count visible (non-ANSI) character width -- assumes single-width chars + some 2-wide emoji *)
let visible_len (s : string) : int =
  let len = String.length s in
  let vis = ref 0 in
  let i = ref 0 in
  let in_esc = ref false in
  while !i < len do
    let c = Char.code s.[!i] in
    if !in_esc then begin
      if (c >= Char.code 'a' && c <= Char.code 'z') ||
         (c >= Char.code 'A' && c <= Char.code 'Z') then
        in_esc := false;
      incr i
    end else if c = 0x1b then begin
      in_esc := true; incr i
    end else if c < 0x80 then begin
      incr vis; incr i
    end else if c >= 0xf0 then begin
      (* 4-byte UTF-8 = emoji, typically 2 display columns *)
      vis := !vis + 2;
      i := !i + 4
    end else if c >= 0xe0 then begin
      (* 3-byte UTF-8 = misc symbols, assume 1 display column *)
      incr vis;
      i := !i + 3
    end else if c >= 0xc0 then begin
      incr vis;
      i := !i + 2
    end else begin
      incr vis; incr i
    end
  done;
  !vis

let pad_to (s : string) (width : int) : string =
  let vl = visible_len s in
  if vl >= width then s
  else s ^ String.make (width - vl) ' '

let hbar (w : int) : string = repeat_glyph glyph_hline w

(* Sidebar logical width in display columns: label(3-4) + tick(1) + data(scope_len) *)
let sidebar_w = 25

let draw ?(unicode = true) (s : sim)
    ~(rocof_prev : float) ~(rocof_next : float) : unit =
  let ob = Buffer.create 8192 in
  let pr fmt = Printf.bprintf ob fmt in
  Buffer.add_string ob "\027[H\027[?25l";
  let hs = h_sys s.g ~s_base:s.s_base in
  let l = load_at s in
  let gmw = gen_total s.g in
  let dp = gmw -. l in
  let on = count_on s.g in
  let total = count_plants s.g in
  let (n_hy, n_md, n_sm, n_so) = count_on_by_type s.g in
  let pct x = if total = 0 then 0 else x * 100 / total in
  let opct x = if !orig_total = 0 then 0 else x * 100 / !orig_total in
  let s_hy = opct n_hy in
  let s_md = opct n_md in
  let s_sm = opct n_sm in
  let s_so = opct n_so in
  scope_push { s_hy; s_md; s_sm; s_so };
  freq_push s.f;

  let hh = s.t / secs_per_hour in
  let mm = (s.t mod secs_per_hour) / secs_per_min in
  let ss = s.t mod 60 in

  let _rows, cols = dim s.g in
  let map_w = cols * 2 in

  let sym (n : node) : string =
    match n.status with
    | On ->
      let stress =
        if n.p_max <= 0.0 then 0
        else min 9 (max 1 (int_of_float (n.p /. n.p_max *. 9.0 +. 0.5)))
      in
      let glyph =
        if unicode then
          (if n.h_s <= 0.0 then glyph_solar
           else if n.p_max >= 2.5 then glyph_hydro
           else if n.p_max >= 1.0 then glyph_medium
           else glyph_small)
        else
          (if n.h_s <= 0.0 then "*"
           else if n.p_max >= 2.5 then "B"
           else if n.p_max >= 1.0 then "M"
           else "S")
      in
      stress_bg stress glyph
    | Off -> if unicode then glyph_tripped else "."
    | Empty -> if unicode then "  " else " "
  in

  (* Dashboard *)
  let cap = gen_capacity s.g in
  let reserve_pct = if l <= 0.0 then 0.0 else (cap -. l) /. l *. 100.0 in
  let bar_w = 12 in
  let filled = max 0 (min bar_w (int_of_float (reserve_pct /. 100.0 *. float_of_int bar_w +. 0.5))) in
  let reserve_color = if reserve_pct > 15.0 then color_freq_ok
    else if reserve_pct > 5.0 then color_freq_warn else color_freq_crit in
  let d1 = Printf.sprintf
    "t=%-6d %02d:%02d:%02d  f=%6.2f Hz  RoCoF %+6.2f / %+6.2f Hz/s  H_sys=%5.2f s"
    s.t hh mm ss s.f rocof_prev rocof_next hs in
  let d2 = Printf.sprintf
    "Gen=%7.1f MW  Load=%7.1f MW  DP=%+7.1f MW  On=%4d/%-4d  Rsv %s%s%s%s%s%+.0f%%%s"
    gmw l dp on total
    reserve_color (repeat_glyph glyph_block_full filled)
    (repeat_glyph glyph_block_light (bar_w - filled)) ansi_reset
    reserve_color reserve_pct ansi_reset in
  let d3 =
    if unicode then
      Printf.sprintf "%s%sHy=%3d%%%s  %s%sMd=%3d%%%s  %s%sSm=%3d%%%s  %s%sSo=%3d%%%s  %sTrip"
        glyph_hydro color_hydro (pct n_hy) ansi_reset
        glyph_medium color_med (pct n_md) ansi_reset
        glyph_small color_small (pct n_sm) ansi_reset
        glyph_solar color_solar (pct n_so) ansi_reset
        glyph_tripped
    else
      Printf.sprintf "Hy=%3d%%  Md=%3d%%  Sm=%3d%%  So=%3d%%"
        (pct n_hy) (pct n_md) (pct n_sm) (pct n_so)
  in
  let n_offsets = Array.length ufls_offsets in
  let d4 =
    if s.ufls_stages_fired > 0 then
      Printf.sprintf "%sUFLS: %d/%d stages fired  shed=%4.1f%%  recovery=%3.0f/%3.0fs%s"
        color_alert s.ufls_stages_fired n_offsets (s.ufls_shed *. 100.0)
        s.ufls_recovery_s ufls_recovery_delay_s ansi_reset
    else ""
  in

  let dash_lines = List.filter (fun s -> s <> "") [d1; d2; d3; d4] in
  let frame_w = map_w + 1 + sidebar_w in
  let max_dash = List.fold_left (fun acc l -> max acc (visible_len l)) 0 dash_lines in
  let total_w = max frame_w max_dash in
  let extra = total_w - frame_w in
  pr "%s%s%s\n" glyph_tl (hbar total_w) glyph_tr;
  List.iter (fun line ->
    pr "%s%s%s\n" glyph_vline (pad_to line total_w) glyph_vline
  ) dash_lines;

  (* Section headers for map + sidebar panels *)
  let fleet_hdr = " Fleet " in
  let fleet_hdr_bar = hbar (sidebar_w + extra - String.length fleet_hdr) in
  pr "%s%s%s%s%s%s\n"
    glyph_ltee (hbar map_w) glyph_ttee fleet_hdr fleet_hdr_bar glyph_rtee;

  (* Grid rows with sidebar *)
  Array.iteri (fun r row ->
    let map_buf = Buffer.create 64 in
    Array.iter (fun n -> Buffer.add_string map_buf (sym n)) row;
    let map_s = Buffer.contents map_buf in

    let sb_buf = Buffer.create 64 in
    (* Panel separator rows *)
    let sb_w = sidebar_w + extra in
    if r = scope_height then begin
      let freq_hdr = " Freq " in
      let freq_bar = hbar (sb_w - String.length freq_hdr) in
      pr "%s%s%s%s%s%s\n"
        glyph_vline (pad_to map_s map_w) glyph_ltee freq_hdr freq_bar glyph_rtee
    end else if r = freq_start_off + freq_height then begin
      let ev_hdr = " Events " in
      let ev_bar = hbar (sb_w - String.length ev_hdr) in
      pr "%s%s%s%s%s%s\n"
        glyph_vline (pad_to map_s map_w) glyph_ltee ev_hdr ev_bar glyph_rtee
    end else begin
      sidebar_content r sb_buf;
      let sb_s = Buffer.contents sb_buf in
      pr "%s%s%s%s%s\n"
        glyph_vline (pad_to map_s map_w) glyph_vline (pad_to sb_s sb_w) glyph_vline
    end
  ) s.g;

  (* Bottom border with status overlay *)
  let status_parts = ref [] in
  if !show_term_size && !term_size_cache <> "" then
    status_parts := !term_size_cache :: !status_parts;
  if !speed_mult <> 1.0 then
    status_parts := Printf.sprintf "%.3gx" !speed_mult :: !status_parts;
  let status = String.concat " " !status_parts in
  let border_inner = map_w + 1 + sidebar_w + extra in
  if status = "" then
    pr "%s%s%s\n" glyph_bl (hbar border_inner) glyph_br
  else begin
    let bar_left = border_inner - String.length status - 1 in
    pr "%s%s %s%s\n" glyph_bl (hbar bar_left) status glyph_br
  end;
  Buffer.add_string ob "\027[?25h";
  output_string stdout (Buffer.contents ob);
  flush stdout

let draw_replay ?(unicode = true)
    ~(type_grid : char array array) ~(grid_status : string) ~(stress : string)
    ~(total : int)
    ~(t : int) ~(f : float) ~(rocof_prev : float) ~(rocof_next : float)
    ~(h_sys : float) ~(gen : float) ~(load : float) ~(dp : float)
    ~(on : int)
    ~(ufls_stages_fired : int) ~(ufls_shed : float)
    ~(n_hy : int) ~(n_md : int) ~(n_sm : int) ~(n_so : int) () : unit =
  let ob = Buffer.create 8192 in
  let pr fmt = Printf.bprintf ob fmt in
  Buffer.add_string ob "\027[H\027[?25l";
  let pct x = if total = 0 then 0 else x * 100 / total in
  let opct x = if !orig_total = 0 then 0 else x * 100 / !orig_total in
  let s_hy = opct n_hy in
  let s_md = opct n_md in
  let s_sm = opct n_sm in
  let s_so = opct n_so in
  scope_push { s_hy; s_md; s_sm; s_so };
  freq_push f;

  let hh = t / secs_per_hour in
  let mm = (t mod secs_per_hour) / secs_per_min in
  let ss = t mod 60 in

  let rows = Array.length type_grid in
  let cols = if rows > 0 then Array.length type_grid.(0) else 0 in
  let map_w = cols * 2 in

  let sym_replay r c =
    let idx = r * cols + c in
    let st = if idx < String.length grid_status then grid_status.[idx] else 'E' in
    match st with
    | 'E' -> if unicode then "  " else " "
    | 'X' -> if unicode then glyph_tripped else "."
    | _ ->
      let tc = type_grid.(r).(c) in
      let sv = if idx < String.length stress then Char.code stress.[idx] - Char.code '0' else 0 in
      let glyph =
        if unicode then
          (match tc with
           | 'O' -> glyph_solar | 'H' -> glyph_hydro
           | 'M' -> glyph_medium | 'S' -> glyph_small | _ -> "  ")
        else
          (match tc with
           | 'O' -> "*" | 'H' -> "B" | 'M' -> "M" | 'S' -> "S" | _ -> " ")
      in
      stress_bg sv glyph
  in

  (* Replay reserve margin: sum p_max of online cells from type_grid *)
  let cap = ref 0.0 in
  for ri = 0 to rows - 1 do
    for ci = 0 to cols - 1 do
      let idx = ri * cols + ci in
      let st = if idx < String.length grid_status then grid_status.[idx] else 'E' in
      if st = 'O' then begin
        let pm = match type_grid.(ri).(ci) with
          | 'H' -> large_baseload.p_max | 'M' -> medium_plant.p_max
          | 'S' -> small_peaker.p_max | 'O' -> solar.p_max | _ -> 0.0 in
        cap := !cap +. pm
      end
    done
  done;
  let cap = !cap in

  (* Dashboard *)
  let reserve_pct = if load <= 0.0 then 0.0 else (cap -. load) /. load *. 100.0 in
  let bar_w = 12 in
  let filled = max 0 (min bar_w (int_of_float (reserve_pct /. 100.0 *. float_of_int bar_w +. 0.5))) in
  let reserve_color = if reserve_pct > 15.0 then color_freq_ok
    else if reserve_pct > 5.0 then color_freq_warn else color_freq_crit in
  let d1 = Printf.sprintf
    "t=%-6d %02d:%02d:%02d  f=%6.2f Hz  RoCoF %+6.2f / %+6.2f Hz/s  H_sys=%5.2f s"
    t hh mm ss f rocof_prev rocof_next h_sys in
  let d2 = Printf.sprintf
    "Gen=%7.1f MW  Load=%7.1f MW  DP=%+7.1f MW  On=%4d/%-4d  Rsv %s%s%s%s%s%+.0f%%%s"
    gen load dp on total
    reserve_color (repeat_glyph glyph_block_full filled)
    (repeat_glyph glyph_block_light (bar_w - filled)) ansi_reset
    reserve_color reserve_pct ansi_reset in
  let d3 =
    if unicode then
      Printf.sprintf "%s%sHy=%3d%%%s  %s%sMd=%3d%%%s  %s%sSm=%3d%%%s  %s%sSo=%3d%%%s  %sTrip"
        glyph_hydro color_hydro (pct n_hy) ansi_reset
        glyph_medium color_med (pct n_md) ansi_reset
        glyph_small color_small (pct n_sm) ansi_reset
        glyph_solar color_solar (pct n_so) ansi_reset
        glyph_tripped
    else
      Printf.sprintf "Hy=%3d%%  Md=%3d%%  Sm=%3d%%  So=%3d%%"
        (pct n_hy) (pct n_md) (pct n_sm) (pct n_so)
  in
  let n_offsets = Array.length ufls_offsets in
  let d4 =
    if ufls_stages_fired > 0 then
      Printf.sprintf "%sUFLS: %d/%d stages fired  shed=%4.1f%%%s"
        color_alert ufls_stages_fired n_offsets (ufls_shed *. 100.0) ansi_reset
    else ""
  in

  let dash_lines = List.filter (fun s -> s <> "") [d1; d2; d3; d4] in
  let frame_w = map_w + 1 + sidebar_w in
  let max_dash = List.fold_left (fun acc l -> max acc (visible_len l)) 0 dash_lines in
  let total_w = max frame_w max_dash in
  let extra = total_w - frame_w in
  pr "%s%s%s\n" glyph_tl (hbar total_w) glyph_tr;
  List.iter (fun line ->
    pr "%s%s%s\n" glyph_vline (pad_to line total_w) glyph_vline
  ) dash_lines;

  let fleet_hdr = " Fleet " in
  let fleet_hdr_bar = hbar (sidebar_w + extra - String.length fleet_hdr) in
  pr "%s%s%s%s%s%s\n"
    glyph_ltee (hbar map_w) glyph_ttee fleet_hdr fleet_hdr_bar glyph_rtee;

  for r = 0 to rows - 1 do
    let map_buf = Buffer.create 64 in
    for c = 0 to cols - 1 do
      Buffer.add_string map_buf (sym_replay r c)
    done;
    let map_s = Buffer.contents map_buf in

    let sb_w = sidebar_w + extra in
    if r = scope_height then begin
      let freq_hdr = " Freq " in
      let freq_bar = hbar (sb_w - String.length freq_hdr) in
      pr "%s%s%s%s%s%s\n"
        glyph_vline (pad_to map_s map_w) glyph_ltee freq_hdr freq_bar glyph_rtee
    end else if r = freq_start_off + freq_height then begin
      let ev_hdr = " Events " in
      let ev_bar = hbar (sb_w - String.length ev_hdr) in
      pr "%s%s%s%s%s%s\n"
        glyph_vline (pad_to map_s map_w) glyph_ltee ev_hdr ev_bar glyph_rtee
    end else begin
      let sb_buf = Buffer.create 64 in
      sidebar_content r sb_buf;
      let sb_s = Buffer.contents sb_buf in
      pr "%s%s%s%s%s\n"
        glyph_vline (pad_to map_s map_w) glyph_vline (pad_to sb_s sb_w) glyph_vline
    end
  done;

  (* Bottom border with status overlay *)
  let status_parts = ref [] in
  if !show_term_size && !term_size_cache <> "" then
    status_parts := !term_size_cache :: !status_parts;
  if !speed_mult <> 1.0 then
    status_parts := Printf.sprintf "%.3gx" !speed_mult :: !status_parts;
  let status = String.concat " " !status_parts in
  let border_inner = map_w + 1 + sidebar_w + extra in
  if status = "" then
    pr "%s%s%s\n" glyph_bl (hbar border_inner) glyph_br
  else begin
    let bar_left = border_inner - String.length status - 1 in
    pr "%s%s %s%s\n" glyph_bl (hbar bar_left) status glyph_br
  end;
  Buffer.add_string ob "\027[?25h";
  output_string stdout (Buffer.contents ob);
  flush stdout

(* ======================================================================
   Simulation step
   ====================================================================== *)

let step (s : sim) : sim * float * float * int * (int * int) list =
  let g = s.g in
  let rows, cols = dim g in

  let n_offsets = Array.length ufls_offsets in
  let ufls_stages_fired = ref s.ufls_stages_fired in
  let ufls_shed = ref s.ufls_shed in

  (* 0. UFLS recovery: reconnect one stage if frequency stable above threshold *)
  let ufls_recovery_s =
    if !ufls_stages_fired > 0 && s.f > s.f0 +. ufls_recovery_df then
      s.ufls_recovery_s +. s.dt_s
    else
      0.0
  in
  let ufls_recovery_s =
    if ufls_recovery_s >= ufls_recovery_delay_s && !ufls_stages_fired > 0 then begin
      ufls_stages_fired := !ufls_stages_fired - 1;
      ufls_shed := !ufls_shed -. ufls_offsets.(!ufls_stages_fired).shed_frac;
      0.0
    end else
      ufls_recovery_s
  in

  let s = { s with
    ufls_stages_fired = !ufls_stages_fired;
    ufls_shed = !ufls_shed;
    ufls_recovery_s } in

  (* 1. Frequency update *)
  let hs = h_sys g ~s_base:s.s_base in
  let l = load_at s in
  let gen = gen_total g in
  let dp = gen -. l in

  let f_next =
    if hs <= 0.0 || s.s_base <= 0.0 then s.f
    else
      let num_mw = dp -. s.beta *. (s.f -. s.f0) in
      let df_dt = (s.f0 /. (2.0 *. hs)) *. (num_mw /. s.s_base) in
      s.f +. s.dt_s *. df_dt
  in

  (* 1b. UFLS: if f_next crossed a threshold, fire stages and recompute *)
  let ufls_stages_fired = ref s.ufls_stages_fired in
  let ufls_shed = ref s.ufls_shed in

  while !ufls_stages_fired < n_offsets
        && f_next < s.f0 +. ufls_offsets.(!ufls_stages_fired).df do
    ufls_shed := !ufls_shed +. ufls_offsets.(!ufls_stages_fired).shed_frac;
    ufls_stages_fired := !ufls_stages_fired + 1
  done;

  let s, f_next =
    if !ufls_stages_fired > s.ufls_stages_fired then
      let s' = { s with
        ufls_stages_fired = !ufls_stages_fired;
        ufls_shed = !ufls_shed } in
      let l' = load_at s' in
      let dp' = gen -. l' in
      let f' =
        if hs <= 0.0 || s.s_base <= 0.0 then s.f
        else
          let num_mw = dp' -. s.beta *. (s.f -. s.f0) in
          let df_dt = (s.f0 /. (2.0 *. hs)) *. (num_mw /. s.s_base) in
          s.f +. s.dt_s *. df_dt
      in
      (s', f')
    else
      (s, f_next)
  in

  (* 2. RoCoF *)
  let rocof_prev = abs_float (s.f -. s.f_prev) /. s.dt_s in
  let rocof_next = abs_float (f_next -. s.f) /. s.dt_s in

  (* 2b. AGC: adjust p_ref toward nominal frequency *)
  let g =
    if s.tau_agc_s <= 0.0 then g
    else
      Array.map (fun row ->
        Array.map (fun n ->
          match n.status with
          | Off | Empty -> n
          | On ->
            let gain = n.p_max /. (n.r_pu *. s.f0 *. s.tau_agc_s) in
            let p_ref' = n.p_ref +. gain *. (s.f0 -. f_next) *. s.dt_s in
            let p_ref' = clip p_ref' 0.0 n.p_max in
            { n with p_ref = p_ref' }
        ) row
      ) g
  in

  (* 3. Trip logic -- per-node thresholds *)
  let freq_under = s.f0 -. f_next in
  let freq_over = f_next -. s.f0 in

  let trip = Array.make_matrix rows cols false in
  let ell_add = Array.make_matrix rows cols 0.0 in

  let g_theta =
    Array.mapi (fun r row ->
      Array.mapi (fun c n ->
        match n.status with
        | Off | Empty -> { n with theta_s = 0.0; freq_dwell = 0.0; freq_dwell_over = 0.0 }
        | On ->
          let p_star_raw = droop_setpoint_raw f_next s.f0 n in
          let theta_s =
            if p_star_raw > n.p_max then n.theta_s +. s.dt_s else 0.0
          in
          let excess_under = freq_under -. n.df_fail in
          let freq_dwell =
            if excess_under <= 0.0 then 0.0
            else n.freq_dwell +. s.dt_s *. (1.0 +. default_freq_dwell_k *. excess_under)
          in
          let excess_over = freq_over -. default_df_fail_over in
          let freq_dwell_over =
            if excess_over <= 0.0 then 0.0
            else n.freq_dwell_over +. s.dt_s *. (1.0 +. default_freq_dwell_over_k *. excess_over)
          in
          let freq_trip = freq_dwell >= n.freq_dwell_base
                          || freq_dwell_over >= default_freq_dwell_over_s in
          let rocof_trip = rocof_next > n.rho in
          let overload_trip = theta_s > n.tover_s in
          let will_trip = freq_trip || rocof_trip || overload_trip in
          trip.(r).(c) <- will_trip;
          { n with theta_s; freq_dwell; freq_dwell_over }
      ) row
    ) g
  in

  (* 4. Redistribute tripped output to neighbors *)
  let tripped_count = ref 0 in
  let tripped_cells = ref [] in
  for r = 0 to rows - 1 do
    for c = 0 to cols - 1 do
      if trip.(r).(c) then begin
        let n = g.(r).(c) in
        match n.status with
        | Off | Empty -> ()
        | On ->
          incr tripped_count;
          tripped_cells := (r, c) :: !tripped_cells;
          let nbrs = neighbors4_coords g r c
            |> List.filter (fun (rr, cc) -> g.(rr).(cc).status <> Empty) in
          let deg = List.length nbrs in
          if deg > 0 then
            let share = n.p /. float_of_int deg in
            List.iter (fun (rr, cc) ->
              ell_add.(rr).(cc) <- ell_add.(rr).(cc) +. share
            ) nbrs
      end
    done
  done;

  (* 5. Update outputs via droop + governor *)
  let g_p =
    Array.mapi (fun r row ->
      Array.mapi (fun c n ->
        if n.status = Empty then n
        else if trip.(r).(c) then
          { n with status = Off; p = 0.0; theta_s = 0.0; freq_dwell = 0.0; freq_dwell_over = 0.0 }
        else
          match n.status with
          | Off | Empty -> { n with p = 0.0; theta_s = 0.0; freq_dwell = 0.0; freq_dwell_over = 0.0 }
          | On ->
            let p_star_raw = droop_setpoint_raw f_next s.f0 n in
            let p_star = clip p_star_raw 0.0 n.p_max in
            let alpha = s.dt_s /. n.tau_gov_s in
            let p' = n.p +. alpha *. (p_star -. n.p) in
            let p' = clip p' 0.0 n.p_max in
            { n with p = p' }
      ) row
    ) g_theta
  in

  (* 6. Update ell buckets for next tick *)
  let decay = exp_decay s.dt_s s.tau_ell_s in
  let g_next =
    Array.mapi (fun r row ->
      Array.mapi (fun c n ->
        match n.status with
        | Off | Empty -> { n with ell = 0.0; theta_s = 0.0; freq_dwell = 0.0; freq_dwell_over = 0.0 }
        | On ->
          let ell' = decay *. n.ell +. ell_add.(r).(c) in
          { n with ell = ell' }
      ) row
    ) g_p
  in

  let s_next =
    { s with
      g = g_next;
      t = s.t + 1;
      f_prev = s.f;
      f = f_next }
  in
  (s_next, rocof_prev, rocof_next, !tripped_count, List.rev !tripped_cells)

(* ======================================================================
   Logging and run loop
   ====================================================================== *)

let grid_status_string (g : grid) : string =
  let buf = Buffer.create 1024 in
  Array.iter (fun row ->
    Array.iter (fun n ->
      Buffer.add_char buf (match n.status with On -> 'O' | Off -> 'X' | Empty -> 'E')
    ) row
  ) g;
  Buffer.contents buf

let cell_type_char (n : node) : char =
  match n.status with
  | Empty -> 'E'
  | On | Off ->
    if n.h_s <= 0.0 then 'O'
    else if n.p_max >= 2.5 then 'H'
    else if n.p_max >= 1.0 then 'M'
    else 'S'

let write_map (filename : string) (g : grid) ~(f0 : float) : unit =
  let oc = open_out filename in
  let rows, cols = dim g in
  let total = count_plants g in
  Printf.fprintf oc "%d,%d,%d,%.1f\n" rows cols total f0;
  Array.iter (fun row ->
    let codes = Array.to_list (Array.map (fun n ->
      String.make 1 (cell_type_char n)
    ) row) in
    Printf.fprintf oc "%s\n" (String.concat "," codes)
  ) g;
  close_out oc

let read_map (filename : string) : int * int * int * float * char array array =
  let ic = open_in filename in
  let header = input_line ic in
  let parts = String.split_on_char ',' header in
  let rows, cols, total, f0 = match parts with
    | [r; c; t; f] -> (int_of_string r, int_of_string c, int_of_string t, float_of_string f)
    | [r; c; t] -> (int_of_string r, int_of_string c, int_of_string t, 60.0)
    | _ -> failwith "bad map header"
  in
  let type_grid = Array.init rows (fun _ ->
    let line = input_line ic in
    let cells = String.split_on_char ',' line in
    Array.of_list (List.map (fun s -> s.[0]) cells)
  ) in
  close_in ic;
  (rows, cols, total, f0, type_grid)

let check_terminal_size ~(req_w : int) ~(req_h : int) : unit =
  try
    let ic = Unix.open_process_in "stty size 2>/dev/null" in
    let line = input_line ic in
    ignore (Unix.close_process_in ic);
    let parts = String.split_on_char ' ' line in
    match parts with
    | [h_s; w_s] ->
      let term_h = int_of_string h_s in
      let term_w = int_of_string w_s in
      if term_w < req_w || term_h < req_h then begin
        Printf.eprintf "Terminal too small: %dx%d, need at least %dx%d\n"
          term_w term_h req_w req_h;
        exit 1
      end
    | _ -> ()
  with _ -> ()

(* Dashboard line 1 is the widest fixed-format line: 76 visible chars + 2 borders *)
let dash_min_w = 78

let required_dims ~(rows : int) ~(cols : int) : int * int =
  let map_w = cols * 2 in
  let frame_w = map_w + 1 + sidebar_w + 2 in
  let w = max frame_w dash_min_w in
  let h = rows + 6 in
  (w, h)

let log_line (oc : out_channel) (s : sim)
    ~(rocof_prev : float) ~(rocof_next : float) ~(trips : int) : unit =
  let hs = h_sys s.g ~s_base:s.s_base in
  let l = load_at s in
  let gen = gen_total s.g in
  let dp = gen -. l in
  let on = count_on s.g in
  let (n_hy, n_md, n_sm, n_so) = count_on_by_type s.g in
  let gs = grid_status_string s.g in
  let ss = grid_stress_string s.g in
  Printf.fprintf oc "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%.6f,%d,%d,%d,%d,%s,%s\n"
    s.t s.f rocof_prev rocof_next hs gen l dp trips on
    s.ufls_stages_fired s.ufls_shed n_hy n_md n_sm n_so gs ss;
  flush oc

let run (s0 : sim) (steps : int) ~(name : string) ~(tick_sleep : float) : unit =
  with_raw_mode @@ fun () ->
  let sim_file = name ^ "_sim.csv" in
  let map_file = name ^ "_map.csv" in
  let oc = open_out sim_file in
  Printf.fprintf oc "t,f,rocof_prev,rocof_next,h_sys,gen,load,dp,trips,online,ufls_stages,ufls_shed,n_hy,n_md,n_sm,n_so,grid_status,stress\n";
  flush oc;
  write_map map_file s0.g ~f0:s0.f0;

  active_f0 := s0.f0;
  orig_total := count_plants s0.g;
  scope_count := 0; scope_idx := 0;
  freq_count := 0; freq_idx := 0;
  Array.fill freq_buf 0 scope_len s0.f0;
  event_log_count := 0; event_log_idx := 0;

  speed_mult := 1.0;
  draw s0 ~rocof_prev:0.0 ~rocof_next:0.0;

  let num_plants = count_plants s0.g in

  let cascade_started = ref false in

  let dispatch_key c =
    if c = key_pause then handle_pause ()
    else handle_common_key c
  in

  let rec aux s =
    if s.t > steps then ()
    else
      let (s_next, rocof_prev, rocof_next, trips, tripped_cells) = step s in

      log_line oc s_next ~rocof_prev ~rocof_next ~trips;

      if trips > 0 then cascade_started := true;

      if tripped_cells <> [] then begin
        let th = ref 0 and tm = ref 0 and ts = ref 0 and to_ = ref 0 in
        List.iter (fun (r, c) ->
          let n = s.g.(r).(c) in
          if n.h_s <= 0.0 then incr to_
          else if n.p_max >= 2.5 then incr th
          else if n.p_max >= 1.0 then incr tm
          else incr ts
        ) tripped_cells;
        let parts = ref [] in
        if !to_ > 0 then parts := Printf.sprintf "%dSo" !to_ :: !parts;
        if !ts > 0 then parts := Printf.sprintf "%dSm" !ts :: !parts;
        if !tm > 0 then parts := Printf.sprintf "%dMd" !tm :: !parts;
        if !th > 0 then parts := Printf.sprintf "%dHy" !th :: !parts;
        event_push (Printf.sprintf "t=%d -%s" s_next.t (String.concat " " !parts))
      end;

      let tick_by_tick = !cascade_started in
      let should_draw = tick_by_tick || s_next.t mod 1800 = 0 in
      let sleep = if tick_by_tick then tick_sleep /. !speed_mult else 0.0 in
      let end_now = (count_on s_next.g = 0) || (trips = num_plants) in
      if should_draw || end_now then begin
        draw s_next ~rocof_prev ~rocof_next;
        if tick_by_tick then begin
          match interruptible_sleep sleep with
          | Some c -> dispatch_key c
          | None -> ()
        end
      end;

      if end_now then () else aux s_next
  in

  (try aux s0 with e -> close_out_noerr oc; raise e);
  close_out oc

let run_replay (name : string) ~(start_tick : int) ~(tick_sleep : float) : unit =
  with_raw_mode @@ fun () ->
  let sim_file = name ^ "_sim.csv" in
  let map_file = name ^ "_map.csv" in
  if not (Sys.file_exists sim_file) then begin
    Printf.eprintf "Error: %s not found\n" sim_file; exit 1 end;
  if not (Sys.file_exists map_file) then begin
    Printf.eprintf "Error: %s not found\n" map_file; exit 1 end;
  let (rows, cols, total, f0, type_grid) = read_map map_file in
  let (req_w, req_h) = required_dims ~rows ~cols in
  check_terminal_size ~req_w ~req_h;
  active_f0 := f0;
  orig_total := total;
  scope_count := 0; scope_idx := 0;
  freq_count := 0; freq_idx := 0;
  Array.fill freq_buf 0 scope_len f0;
  event_log_count := 0; event_log_idx := 0;
  speed_mult := 1.0;
  let ic = open_in sim_file in
  ignore (input_line ic);
  clear ();
  let cascade_started = ref false in
  let prev_gs = ref "" in

  let dispatch_key c =
    if c = key_pause then handle_pause ()
    else handle_common_key c
  in

  (try while true do
    let line = input_line ic in
    let parts = String.split_on_char ',' line in
    match parts with
    | t_s :: f_s :: rp_s :: rn_s :: hs_s :: gen_s :: load_s :: dp_s
      :: trips_s :: on_s :: ufls_st_s :: ufls_sh_s
      :: nhy_s :: nmd_s :: nsm_s :: nso_s :: gs_s :: stress_rest ->
      let t = int_of_string t_s in
      let f = float_of_string f_s in
      let rocof_prev = float_of_string rp_s in
      let rocof_next = float_of_string rn_s in
      let h_sys = float_of_string hs_s in
      let gen = float_of_string gen_s in
      let load = float_of_string load_s in
      let dp = float_of_string dp_s in
      let trips = int_of_string trips_s in
      let on = int_of_string on_s in
      let ufls_stages_fired = int_of_string ufls_st_s in
      let ufls_shed = float_of_string ufls_sh_s in
      let n_hy = int_of_string nhy_s in
      let n_md = int_of_string nmd_s in
      let n_sm = int_of_string nsm_s in
      let n_so = int_of_string nso_s in
      let grid_status = gs_s in
      let stress = String.concat "," stress_rest in

      (* Generate trip events by diffing grid_status with previous *)
      if !prev_gs <> "" && String.length !prev_gs = String.length grid_status then begin
        let th = ref 0 and tm = ref 0 and ts = ref 0 and to_ = ref 0 in
        for idx = 0 to String.length grid_status - 1 do
          if !prev_gs.[idx] = 'O' && grid_status.[idx] = 'X' then begin
            let r = idx / cols in
            let c = idx mod cols in
            match type_grid.(r).(c) with
            | 'O' -> incr to_ | 'H' -> incr th
            | 'M' -> incr tm | 'S' -> incr ts | _ -> ()
          end
        done;
        if !th + !tm + !ts + !to_ > 0 then begin
          let parts = ref [] in
          if !to_ > 0 then parts := Printf.sprintf "%dSo" !to_ :: !parts;
          if !ts > 0 then parts := Printf.sprintf "%dSm" !ts :: !parts;
          if !tm > 0 then parts := Printf.sprintf "%dMd" !tm :: !parts;
          if !th > 0 then parts := Printf.sprintf "%dHy" !th :: !parts;
          event_push (Printf.sprintf "t=%d -%s" t (String.concat " " !parts))
        end
      end;
      prev_gs := grid_status;

      if trips > 0 then cascade_started := true;
      let past_start = t >= start_tick in
      let tick_by_tick =
        if start_tick > 0 then past_start else !cascade_started
      in
      let should_draw =
        tick_by_tick || (start_tick = 0 && t mod 1800 = 0)
      in
      let sleep = if tick_by_tick then tick_sleep /. !speed_mult else 0.0 in
      if should_draw then begin
        draw_replay
          ~type_grid ~grid_status ~stress ~total
          ~t ~f ~rocof_prev ~rocof_next ~h_sys ~gen ~load ~dp
          ~on ~ufls_stages_fired ~ufls_shed
          ~n_hy ~n_md ~n_sm ~n_so ();
        if tick_by_tick then begin
          match interruptible_sleep sleep with
          | Some c -> dispatch_key c
          | None -> ()
        end
      end
    | _ -> ()
  done with End_of_file -> ());
  close_in ic;
  ignore rows

(* ======================================================================
   Main
   ====================================================================== *)

let () =
  let argc = Array.length Sys.argv in
  if argc < 3 then begin
    Printf.eprintf "Usage: %s run <name> [tick_sleep]\n" Sys.argv.(0);
    Printf.eprintf "       %s replay <name> [start_tick] [tick_sleep]\n" Sys.argv.(0);
    exit 1
  end;
  match Sys.argv.(1) with
  | "replay" ->
    let name = Sys.argv.(2) in
    let start_tick = if argc > 3 then int_of_string Sys.argv.(3) else 0 in
    let tick_sleep = if argc > 4 then float_of_string Sys.argv.(4) else 0.5 in
    run_replay name ~start_tick ~tick_sleep
  | "run" ->
    let name = Sys.argv.(2) in
    let tick_sleep = if argc > 3 then float_of_string Sys.argv.(3) else 0.5 in
    let rows, cols = profile_dims default_profile in
    let (req_w, req_h) = required_dims ~rows ~cols in
    check_terminal_size ~req_w ~req_h;
    let f0 = profile_f0 default_profile in

    let g = init_grid default_profile rows cols default_seed in

    let s_base = s_base_of_grid g in
    let pref_total = total_p_ref g in

    let l_mid = default_load_factor *. pref_total in
    let l_amp = default_amp_fraction *. pref_total in

    let period_s = secs_per_day_f in
    let load_t0 =
      l_mid +. l_amp *. cos (2.0 *. Float.pi *. (float_of_int (0 - default_t_peak_s)) /. period_s)
    in
    let scale = load_t0 /. pref_total in
    let g = Array.map (fun row ->
      Array.map (fun n ->
        match n.status with
        | On -> { n with p = n.p_ref *. scale; p_ref = n.p_ref *. scale }
        | Off | Empty -> n
      ) row
    ) g in

    clear ();

    let s0 = {
      g;
      t = 0;
      dt_s = default_dt_s;

      f = f0;
      f_prev = f0;
      f0;

      beta = default_beta_pu *. l_mid;

      l_mid;
      l_amp;
      t_peak_s = default_t_peak_s;

      load_step_mw = default_load_step_mw;
      load_step_t = default_load_step_t;
      load_ramp_mw_per_s = default_load_ramp_mw_per_s;

      s_base;

      tau_ell_s = default_tau_ell_s;
      tau_agc_s = default_tau_agc_s;

      load_noise_pu = default_load_noise_pu;
      load_rng = Random.State.make [| default_seed + 1 |];

      ufls_stages_fired = 0;
      ufls_shed = 0.0;
      ufls_recovery_s = 0.0;
    } in

    run s0 (default_sim_hours * secs_per_hour) ~name ~tick_sleep
  | mode ->
    Printf.eprintf "Unknown mode: %s\n" mode;
    exit 1
