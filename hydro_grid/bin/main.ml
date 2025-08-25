[@@@warning "-32-37"]

type status = On | Off
type node = { status : status }
type grid = node array array

let init_node ?(status=On) () = { status }
let init_grid rows cols ~status =
  Array.make_matrix rows cols (init_node ~status:status ())

let with_off_positions (g : grid) (positions : (int*int) list) : grid =
  let off_tbl = Hashtbl.create (List.length positions) in
  List.iter (fun pos -> Hashtbl.replace off_tbl pos ()) positions;
  Array.mapi (fun i row ->
    Array.mapi (fun j n ->
      if Hashtbl.mem off_tbl (i,j) then {status=Off} else n
    ) row
  ) g

let in_bounds g r k =
  let rows = Array.length g in
  let cols = if rows = 0 then 0 else Array.length g.(0) in
  r >= 0 && r < rows && k >= 0 && k < cols

let neighbors g r k =
  let canditates = [
    (r - 1, k - 1); (r - 1, k); (r - 1, k + 1);
    (r, k - 1);                 (r, k + 1);
    (r + 1, k - 1); (r + 1, k); (r + 1, k + 1)
  ] in

  List.filter (fun (r, k) -> in_bounds g r k) canditates

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

  let g = with_off_positions g [(0, 0); (1, 1); (2, 2); (3, 3); (4, 4);] in
  draw g;

