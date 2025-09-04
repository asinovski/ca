(** 2D grid of hydro-turbine nodes.

    This file defines:
    - a turbine on/off [status],
    - a cell type [node],
    - a rectangular [grid] as [node array array],
    - constructors for nodes and grids,
    - helpers to switch selected coordinates to [Off],
    - bound checking and 8-neighborhood lookup,
    - a simple renderer to stdout.
*)

[@@@warning "-32-37"]

(** Turbine operating state. *)
type status = On | Off

(** A single grid cell. Extend later with kind/inertia/etc. *)
type node = { status : status }

(** Rectangular grid, indexed as [g.(row).(col)]. *)
type grid = node array array

(** [init_node ?status ()] creates a node with the given [status]
    (default = [On]). *)
let init_node ?(status=On) () = { status }

(** [init_grid rows cols ~status] builds a [rows × cols] grid where every cell
    is initialized to a node with [status].

    @param rows number of rows (≥ 0)
    @param cols number of columns (≥ 0)
    @param status initial status for all cells
    @return a fresh rectangular grid
    @raise Invalid_argument if [rows] or [cols] is negative
    @note Cells share the same immutable node value; sharing is harmless
          because [node] has no mutable fields. *)
let init_grid rows cols ~status =
  Array.make_matrix rows cols (init_node ~status:status ())

(** [with_off_positions g positions] returns a new grid where any cell at a
    coordinate in [positions] is set to [{status=Off}]. All other cells are
    copied from [g] unchanged.

    Out-of-bounds coordinates are ignored (no exception is raised).

    @param g source grid
    @param positions list of [(row, col)] pairs to switch off
    @return a new grid with requested cells off
    @since 0.1 *)
let with_off_positions (g : grid) (positions : (int*int) list) : grid =
  let off_tbl = Hashtbl.create (List.length positions) in
  List.iter (fun pos -> Hashtbl.replace off_tbl pos ()) positions;
  Array.mapi (fun i row ->
    Array.mapi (fun j n ->
      if Hashtbl.mem off_tbl (i,j) then {status=Off} else n
    ) row
  ) g

(** [in_bounds g r k] tells whether [(r, k)] is a valid index into [g]
    (i.e. [0 ≤ r < rows] and [0 ≤ k < cols]). *)
let in_bounds g r k =
  let rows = Array.length g in
  let cols = if rows = 0 then 0 else Array.length g.(0) in
  r >= 0 && r < rows && k >= 0 && k < cols

(** [neighbors g r k] returns the *nodes* in the 8-neighborhood (Moore
    neighborhood) around [(r,k)], filtered to in-bounds positions.

    The order is:
    [(r-1,k-1); (r-1,k); (r-1,k+1); (r,k-1); (r,k+1); (r+1,k-1); (r+1,k); (r+1,k+1)].

    Edge/corner cells naturally yield fewer neighbors.

    @param g grid
    @param r row index
    @param k column index
    @return a list of neighboring nodes (not their coordinates) *)
let neighbors g r k =
  let candidates = [
    (r - 1, k - 1); (r - 1, k); (r - 1, k + 1);
    (r, k - 1);                 (r, k + 1);
    (r + 1, k - 1); (r + 1, k); (r + 1, k + 1)
  ] in
  candidates
  |> List.filter (fun (r, k) -> in_bounds g r k)
  |> List.map (fun (r, k) -> g.(r).(k))

(** [draw ?unicode g] prints [g] to stdout as a compact grid.

    When [unicode] is [true] (default), online cells are rendered as ["🌀"]
    and offline cells as ["❌"]. When [unicode=false], they render as ["T"]
    and ["X"] respectively. Each row ends with a newline.

    @param unicode toggle Unicode vs ASCII glyphs
    @param g grid to render
    @return unit (prints to stdout) *)
let draw ?(unicode=true) (g : grid) =
  let sym = function
    | {status=On; _}  -> if unicode then "🌀" else "T"
    | {status=Off; _} -> if unicode then "❌" else "X"
  in
  Array.iter (fun row ->
    Array.iter (fun n -> print_string (sym n ^ "")) row;
    print_newline ()
  ) g

let () =
  let g = init_grid 5 5 ~status:On in
  draw g;
  print_newline ();
  let diag = List.init 5 (fun i -> (i, i)) in
  let g = with_off_positions g diag in
  draw g;
  List.iter
    (fun (r,c) ->
       let len = List.length (neighbors g r c) in
       Printf.printf "Len of neighbors of (%d,%d): %d\n" r c len)
    diag
