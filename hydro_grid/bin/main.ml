(** 2D grid of hydro-turbine nodes with:
    - global frequency f(t) (single-bus approximation),
    - per-node inertia H_i (seconds) and capacity P_max,i (MW),
    - aggregate inertia computed capacity-weighted,
    - primary control: droop + 1st-order governor + saturation,
    - trip rules (freq / RoCoF / overload dwell),
    - local redistribution of tripped output as a decaying "proximal load" ℓ
      to Von Neumann (N/S/E/W) neighbors.

    Time is discrete: one tick = 1 second (dt_s = 1.0).
    Load varies slowly as a daily sinusoid with period 24 hours (86400 ticks).

    One-tick delay (local stress propagation):
    - When a unit trips at tick t, we redistribute its last output into its
      neighbors' proximal-load buckets as ℓ(t+Δt).
    - That redistributed ℓ does NOT influence droop targets until tick (t+1),
      because we compute droop setpoints P⋆(t) using ℓ(t) (the value entering
      the tick), update governor outputs, and only then update ℓ(t+Δt).

    Printing: we render the grid every 60 ticks (once per simulation minute).
    If draw sleep_time = 1.0, you'll see one simulated minute per real second.
*)

[@@@warning "-32-37"]

type status = On | Off

type node = {
  status : status;

  (* Inertia constant H_i (seconds). *)
  h_s : float;

  (* Electrical power state (MW). *)
  p : float;      (* current output *)
  p_ref : float;  (* reference at nominal frequency (MW) *)
  p_max : float;  (* max output / saturation (MW) *)

  (* Primary control parameters. *)
  r_pu : float;       (* droop in per-unit, e.g. 0.05 (5% droop) *)
  tau_gov_s : float;  (* governor time constant (seconds) *)

  (* Cascading/coupling states. *)
  ell : float;      (* proximal load bucket ℓ_i (MW) *)
  theta_s : float;  (* overload dwell accumulator θ_i (seconds) *)
}

type grid = node array array

type sim = {
  g : grid;

  (* Time *)
  t : int;       (* tick index; 1 tick = 1 second *)
  dt_s : float;  (* seconds per tick (1.0) *)

  (* Frequency *)
  f : float;      (* Hz *)
  f_prev : float; (* Hz, previous tick *)
  f0 : float;     (* Hz, nominal *)

  (* Load damping coefficient (MW/Hz). *)
  beta : float;

  (* Daily load model (MW):
     load(t) = l_mid + l_amp*cos(2π(t - t_peak_s)/(24h))
     where t is in seconds.
  *)
  l_mid : float;
  l_amp : float;
  t_peak_s : int;  (* seconds into the day of peak load, e.g. 12*3600 *)

  (* System base (MW). Used to normalize ΔP into per-unit. *)
  s_base : float;

  (* Proximal load decay constant (seconds). *)
  tau_ell_s : float;

  (* Trip thresholds *)
  df_fail : float;  (* Hz band trip *)
  rho : float;      (* RoCoF threshold (Hz/s) *)
  tover_s : float;  (* overload dwell threshold (seconds) *)
}

let clip (x : float) (lo : float) (hi : float) : float =
  if x < lo then lo else if x > hi then hi else x

let clear () =
  print_string "\027[2J\027[H"

let in_bounds (g : grid) (r : int) (c : int) : bool =
  let rows = Array.length g in
  let cols = if rows = 0 then 0 else Array.length g.(0) in
  r >= 0 && r < rows && c >= 0 && c < cols

(* Von Neumann (4-neighborhood) coordinates. *)
let neighbors4_coords (g : grid) (r : int) (c : int) : (int * int) list =
  let candidates = [ (r-1,c); (r+1,c); (r,c-1); (r,c+1) ] in
  candidates |> List.filter (fun (rr,cc) -> in_bounds g rr cc)

let exp_decay (dt_s : float) (tau_s : float) : float =
  if tau_s <= 0.0 then 0.0 else exp (-. dt_s /. tau_s)

(* Daily sinusoidal load: max at t_peak_s, min 12 hours later. *)
let load_at (s : sim) : float =
  let period_s = 24.0 *. 3600.0 in
  let x = 2.0 *. Float.pi *. (float_of_int (s.t - s.t_peak_s)) /. period_s in
  s.l_mid +. s.l_amp *. cos x

let init_node
  ?(status=On)
  ?(h_s=5.0)
  ?(p_max=1.2)
  ?(p_ref=1.0)
  ?(r_pu=0.05)
  ?(tau_gov_s=5.0)
  () : node =
  let p_ref = Float.min p_ref p_max in
  let p = match status with On -> p_ref | Off -> 0.0 in
  { status; h_s; p; p_ref; p_max; r_pu; tau_gov_s; ell = 0.0; theta_s = 0.0 }

let init_grid rows cols ~(status:status) ~(h_s:float) ~(p_max:float) ~(p_ref:float) : grid =
  Array.init rows (fun _ ->
    Array.init cols (fun _ -> init_node ~status ~h_s ~p_max ~p_ref ())
  )

let with_off_positions (g : grid) (positions : (int*int) list) : grid =
  let off_tbl = Hashtbl.create (List.length positions) in
  List.iter (fun pos -> Hashtbl.replace off_tbl pos ()) positions;
  Array.mapi (fun i row ->
    Array.mapi (fun j n ->
      if Hashtbl.mem off_tbl (i,j) then { n with status = Off; p = 0.0 } else n
    ) row
  ) g

(* Total generation (MW). *)
let gen_total (g : grid) : float =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n ->
      match n.status with
      | On  -> acc +. n.p
      | Off -> acc
    ) acc row
  ) 0.0 g

(* System base (MW): sum of all nameplate capacities. *)
let s_base_of_grid (g : grid) : float =
  Array.fold_left (fun acc row ->
    Array.fold_left (fun acc n -> acc +. n.p_max) acc row
  ) 0.0 g

(* Capacity-weighted aggregate inertia:
     H_sys = (Σ_on H_i * S_i) / S_base
   Here S_i is approximated by p_max (MW), and S_base is a fixed system base.
*)
let h_sys (g : grid) ~(s_base:float) : float =
  if s_base <= 0.0 then 0.0
  else
    let num =
      Array.fold_left (fun acc row ->
        Array.fold_left (fun acc n ->
          match n.status with
          | On  -> acc +. (n.h_s *. n.p_max)
          | Off -> acc
        ) acc row
      ) 0.0 g
    in
    num /. s_base

(* Droop setpoint (MW):
   We interpret r_pu as a fractional frequency deviation of f0 that would call
   for ~1× change in power relative to p_ref. This yields a gain:
     k ≈ p_ref / (r_pu * f0)  [MW/Hz]
   and
     P⋆ = p_ref + k*(f0 - f) + ℓ

   Note: ℓ here is the bucket value entering the tick (one-tick delay).
*)
let droop_setpoint_raw (f : float) (f0 : float) (n : node) : float =
  let k_mw_per_hz =
    if n.r_pu <= 0.0 then 0.0 else n.p_ref /. (n.r_pu *. f0)
  in
  n.p_ref +. k_mw_per_hz *. (f0 -. f) +. n.ell

let draw ?(unicode=true) ?(sleep_time=1.0) (s : sim)
  ~(rocof_prev:float) ~(rocof_next:float) ~(trips:int) : unit =
  clear ();
  let hs = h_sys s.g ~s_base:s.s_base in
  let l  = load_at s in
  let gmw = gen_total s.g in
  let dp = gmw -. l in
  let mm = s.t / 60 in
  let ss = s.t mod 60 in
  Printf.printf
    "t=%d (%02d:%02d)  f=%.4f Hz  RoCoF(prev)=%.4f Hz/s  RoCoF(next)=%.4f Hz/s  H_sys=%.3f s  Gen=%.3f MW  Load=%.3f MW  ΔP=%.3f MW  trips=%d\n\n"
    s.t mm ss s.f rocof_prev rocof_next hs gmw l dp trips;

  let sym = function
    | {status=On; _}  -> if unicode then "🌀" else "T"
    | {status=Off; _} -> if unicode then "❌" else "X"
  in
  Array.iter (fun row ->
    Array.iter (fun n -> print_string (sym n)) row;
    print_newline ()
  ) s.g;
  print_newline ();
  flush stdout;
  Unix.sleepf sleep_time

(* One simulation tick (dt = 1 second):
   1) Update frequency via a normalized swing step.
   2) Compute RoCoF.
   3) Trip logic (freq / RoCoF / overload dwell).
   4) Redistribute tripped output into neighbors' ell_add (N/S/E/W).
   5) Update P via droop + 1st-order governor using ℓ(t) only.
   6) Update ℓ(t+Δt) = decay*ℓ(t) + ell_add (one-tick delay).
*)
let step (s : sim) : sim * float * float * int =
  let g = s.g in
  let rows = Array.length g in
  let cols = if rows = 0 then 0 else Array.length g.(0) in

  (* ---- (1) Frequency update (normalized by system base) ---- *)
  let hs = h_sys g ~s_base:s.s_base in
  let l  = load_at s in
  let gen = gen_total g in
  let dp = gen -. l in

  (* df/dt ≈ (f0 / (2H_sys)) * ((ΔP - β(f-f0)) / S_base) *)
  let f_next =
    if hs <= 0.0 || s.s_base <= 0.0 then s.f
    else
      let num_mw = dp -. s.beta *. (s.f -. s.f0) in
      let df_dt = (s.f0 /. (2.0 *. hs)) *. (num_mw /. s.s_base) in
      s.f +. s.dt_s *. df_dt
  in

  (* ---- (2) RoCoF (prev uses f_prev; next uses the step we just computed) ---- *)
  let rocof_prev = abs_float (s.f -. s.f_prev) /. s.dt_s in
  let rocof_next = abs_float (f_next -. s.f) /. s.dt_s in

  (* ---- (3) Trip logic and overload dwell update ---- *)
  let freq_trip_all = abs_float (f_next -. s.f0) > s.df_fail in
  let rocof_trip_all = rocof_next > s.rho in

  let trip = Array.make_matrix rows cols false in
  let ell_add = Array.make_matrix rows cols 0.0 in

  (* We accumulate θ when the unit is being "asked" to exceed its max.
     This uses the *raw* droop setpoint before clipping. *)
  let g_theta =
    Array.mapi (fun r row ->
      Array.mapi (fun c n ->
        match n.status with
        | Off -> { n with theta_s = 0.0 }
        | On ->
          let p_star_raw = droop_setpoint_raw f_next s.f0 n in
          let theta_s =
            if p_star_raw > n.p_max then n.theta_s +. s.dt_s else 0.0
          in
          let overload_trip = theta_s > s.tover_s in
          let will_trip = freq_trip_all || rocof_trip_all || overload_trip in
          trip.(r).(c) <- will_trip;
          { n with theta_s }
      ) row
    ) g
  in

  (* ---- (4) Redistribute tripped output into ell_add (one-tick delay) ---- *)
  let tripped_count = ref 0 in
  for r = 0 to rows - 1 do
    for c = 0 to cols - 1 do
      if trip.(r).(c) then begin
        let n = g.(r).(c) in
        match n.status with
        | Off -> ()
        | On ->
          incr tripped_count;
          let nbrs = neighbors4_coords g r c in
          let deg = List.length nbrs in
          if deg > 0 then
            let share = n.p /. float_of_int deg in
            List.iter (fun (rr,cc) ->
              ell_add.(rr).(cc) <- ell_add.(rr).(cc) +. share
            ) nbrs
      end
    done
  done;

  (* ---- (5) Update outputs with droop+governor using ℓ(t) only ---- *)
  let g_p =
    Array.mapi (fun r row ->
      Array.mapi (fun c n ->
        if trip.(r).(c) then
          { n with status = Off; p = 0.0; theta_s = 0.0 }
        else
          match n.status with
          | Off -> { n with p = 0.0; theta_s = 0.0 }
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

  (* ---- (6) Update ℓ buckets for NEXT tick: ℓ(t+Δt) ---- *)
  let decay = exp_decay s.dt_s s.tau_ell_s in
  let g_next =
    Array.mapi (fun r row ->
      Array.mapi (fun c n ->
        match n.status with
        | Off ->
          { n with ell = 0.0; theta_s = 0.0 }
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
      f = f_next;
    }
  in
  (s_next, rocof_prev, rocof_next, !tripped_count)

let log_line (oc : out_channel) (s : sim) ~(rocof:float) ~(trips:int) : unit =
  let hs = h_sys s.g ~s_base:s.s_base in
  let l  = load_at s in
  let gen = gen_total s.g in
  let dp = gen -. l in
  Printf.fprintf oc "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n"
    s.t s.f rocof hs gen l dp trips

let run (s0 : sim) (steps : int) : unit =
  let oc = open_out "sim_log.csv" in
  Printf.fprintf oc "t,f,rocof,h_sys,gen,load,dp,trips\n";
  flush oc;

  let rec aux s =
    if s.t > steps then ()
    else
      let (s_next, rocof_prev, rocof_next, trips) = step s in

      (* Log every tick (1 second). Use rocof_next because it corresponds to
         the transition that produced s_next (|f_next - f| / dt). *)
      log_line oc s_next ~rocof:rocof_next ~trips;

      (* Render only once per simulation minute (every 60 ticks). *)
      if s_next.t mod 60 = 0 then
        draw s_next ~rocof_prev ~rocof_next ~trips;

      aux s_next
  in

  (try aux s0 with e -> close_out_noerr oc; raise e);
  close_out oc

let () =
  (* Demo: 5x5 grid. *)
  let g =
    init_grid 5 5 
      ~status:On
      ~h_s:5.0
      ~p_max:1.5
      ~p_ref:1.0
  in

  (* Switch off diagonal as initial outages. *)
  let diag = List.init 5 (fun i -> (i, i)) in
  let g = with_off_positions g diag in

  let s_base = s_base_of_grid g in

  let s0 = {
    g;

    t = 5 * 3600 + 180;
    dt_s = 1.0;          (* 1 tick = 1 second *)

    f = 60.0;
    f_prev = 60.0;
    f0 = 60.0;

    beta = 0.0;          (* MW/Hz; tune later *)

    (* Load chosen to be near initial generation to avoid huge initial shocks.
       With 20 online nodes at ~1 MW each, gen ~20 MW at start. *)
    l_mid = 22.0;
    l_amp = 4.0;
    t_peak_s = 12 * 3600;

    s_base;

    tau_ell_s = 2.0 *. 3600.0;  (* proximal load decays over ~2 hours *)

    df_fail = 0.5;       (* Hz band trip *)
    rho = 1.0;           (* Hz/s RoCoF trip *)
    tover_s = 1.0;      (* seconds at "asked-above-max" before trip *)
  } in

  (* Run for 6 simulation hours (21600 seconds). *)
  run s0 (12 * 3600)
