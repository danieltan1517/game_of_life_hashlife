// hashlife_console.jai — ANSI console front-end.
// The simulation engine lives in hashlife_core.jai (shared with the UI and the
// test harness); this file only adds the terminal renderer and the CLI main().

#load "hashlife_core.jai";

#import "Basic";    // print, sleep_milliseconds, get_command_line_arguments, string_to_int
#import "String";   // ends_with (used to detect *.lif arguments)

// Render a VIEW_W × VIEW_H window automatically centred on the live cells.
// HashLife shifts the pattern's coordinates within the node on every step
// (centre/crop change the origin), so a fixed window would lose the pattern
// even though the population is intact. Re-deriving the centre each frame
// from the live-cell bounding box keeps the pattern in view.

VIEW_W  :: 64;
VIEW_H  :: 32;

render :: (t: *Node_Table, id: node_id, gen: u64) {
    // Clear screen with ANSI escape (works on most terminals)
    print("\x1b[2J\x1b[H");

    print("HashLife — Conway's Game of Life   gen: %   pop: %   level: %\n\n",
          gen, lookup(t, id).pop, level_of(id));

    b := find_bounds(t, id);
    // Centre the view on the pattern's bounding box (fall back to 0,0 if empty).
    cx : s64 = 0;
    cy : s64 = 0;
    if b.found {
        cx = (b.min_x + b.max_x) / 2;
        cy = (b.min_y + b.max_y) / 2;
    }

    half_w := VIEW_W / 2;
    half_h := VIEW_H / 2;
    size   := cast(u64) (1 << level_of(id));

    for row : 0 .. VIEW_H - 1 {
        for col : 0 .. VIEW_W - 1 {
            world_x := cx - half_w + col;
            world_y := cy - half_h + row;

            alive := false;
            if world_x >= 0 && world_y >= 0 {
                ux := cast(u64) world_x;
                uy := cast(u64) world_y;
                if ux < size && uy < size {
                    alive = get_cell(t, id, ux, uy, 0) > 0.0;
                }
            }
            print(ifx alive then "█" else "·");
        }
        print("\n");
    }
    print("\n");
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

print_usage :: () {
    print("Usage: hashlife [options] [file.lif]\n");
    print("\n");
    print("  file.lif         Life 1.05 pattern file to load\n");
    print("\n");
    print("  --glider         Built-in glider (default if no file given)\n");
    print("  --gun            Built-in Gosper glider gun\n");
    print("  --rpentomino     Built-in R-pentomino\n");
    print("\n");
    print("  --frames  <n>    Total frames to simulate (default 120)\n");
    print("  --sleep   <ms>   Milliseconds between frames (default 80)\n");
    print("\n");
    print("  --help           Show this message\n");
}

main :: () {
    t := create_table();
    defer free_table(t);

    // --- Parse command-line arguments ---
    // Jai stores argv in __command_line_arguments ([] string, index 0 is the exe).
    args := get_command_line_arguments();

    lif_file    : string;
    builtin     : string = "glider";   // default
    total_frames : int    = 120;
    sleep_ms     : s32    = 80;

    i := 1;   // skip argv[0] (executable name)
    while i < args.count {
        arg := args[i];
        if arg == {
        case "--help";
            print_usage();
            return;
        case "--glider";
            builtin  = "glider";
            lif_file = "";
        case "--gun";
            builtin  = "gun";
            lif_file = "";
        case "--rpentomino";
            builtin  = "rpentomino";
            lif_file = "";
        case "--frames";
            i += 1;
            if i < args.count {
                v, ok := string_to_int(args[i]);
                if ok && v > 0  total_frames = v;
            }
        case "--sleep";
            i += 1;
            if i < args.count {
                v, ok := string_to_int(args[i]);
                if ok && v >= 0  sleep_ms = cast(s32) v;
            }
        case;
            // Treat unrecognised arguments that end in .lif/.LIF as file paths.
            if ends_with(arg, ".lif") || ends_with(arg, ".LIF") {
                lif_file = arg;
            } else {
                print("Unknown argument: '%'. Run with --help for usage.\n", arg);
            }
        }
        i += 1;
    }

    // --- Load pattern ---
    id : node_id;

    if lif_file.count > 0 {
        ok : bool;
        id, ok = load_lif_file(t, lif_file);
        if !ok  return;
        print("Loaded: %   (initial pop: %)\n", lif_file, lookup(t, id).pop);
    } else {
        if builtin == {
        case "glider";
            print("Pattern: Glider\n");
            id = load_glider(t);
        case "gun";
            print("Pattern: Gosper Glider Gun\n");
            id = load_glider_gun(t);
        case "rpentomino";
            print("Pattern: R-pentomino\n");
            id = load_r_pentomino(t);
        }
    }

    // --- Simulation loop ---
    // advance(t, id, 1) with j starting at 0: bit 0 → successor(id, 0) = 1 generation.
    GENS_PER_ADVANCE :: 1;   // value passed to advance()
    STEP_SIZE        :: 1;   // actual generations per call (now correct)

    gen : u64 = 0;

    print("Starting simulation. Press Ctrl+C to stop.\n");
    sleep_milliseconds(800);

    for frame : 0 .. total_frames - 1 {
        render(t, id, gen);
        sleep_milliseconds(sleep_ms);

        id   = advance(t, id, GENS_PER_ADVANCE);
        gen += STEP_SIZE;

        if lookup(t, id).pop == 0 {
            render(t, id, gen);
            print("Pattern died at generation %.\n", gen);
            break;
        }
    }

    print("Simulation complete.  Final population: %  generations: %\n",
          lookup(t, id).pop, gen);
}
