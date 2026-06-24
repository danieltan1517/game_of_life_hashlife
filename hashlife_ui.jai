/*
    Conway's Game of Life — HashLife engine + interactive UI
    =========================================================

    This file combines two originals:

      * gol.jai             — a Simp/GetRect windowed UI that previously drove a
                              fixed 128x128 toroidal [M][M] s8 grid with a naive
                              O(M^2) generation step.
      * hashlife_console.jai — Gosper's HashLife algorithm (quadtree of hash-consed
                              macro-cells with memoised successors), previously
                              wrapped in a console renderer.

    The UI now drives the HashLife engine instead of the flat array, so the world
    is effectively unbounded and repetitive patterns advance in near-O(log t).

    Integration notes (read these — HashLife is not a drop-in for an array):

      1. State.  `gol_grid`/`gol_grid_next` are gone. The world is a single
         `node_id` (`gol_id`) interned in a `Node_Table` (`t`). "Clear" is just
         `gol_id = t.off`; one generation is `advance(t, gol_id, 1)`.

      2. Shifting coordinate frame.  `advance()` calls `centre`/`crop`, which
         re-root the quadtree, so a live cell's *node coordinates* change between
         generations. A fixed camera would lose the pattern. We therefore derive
         the camera from the live-cell bounding box (`find_bounds`) every frame
         while running — the same trick the console renderer used. This needs no
         changes to the algorithm itself.

      3. Camera.  `cam_x`/`cam_y` are the node-frame coordinates of the top-left
         visible cell. While playing we auto-follow the pattern's bbox centre.
         While paused the camera is stable and the arrow keys pan it.

      4. Editing.  The quadtree here only grows toward +x/+y (set_cell appends
         zero quadrants in b/c/d), so negative node coordinates aren't directly
         representable. `centre` shifts all live cells by +2^(level-1); the click
         handler uses that to pad the tree until the clicked cell is representable,
         adjusting the camera by the same amount so nothing appears to move.
         Editing is enabled only while paused.

      5. Step size.  HashLife's advantage is jumping many generations cheaply, so
         a "Step ... gens" slider feeds `advance(t, gol_id, n)` directly.

    Caveat: while running, the camera follows the bounding-box *centre*. For an
    emitter like the Gosper gun the bbox grows as gliders stream off, so the view
    drifts toward the centroid — pan with the arrow keys (after pausing) or pick a
    bounded pattern to watch it sit still.

    Command-line options (ported from the original console application):

        hashlife_ui [options] [file.lif]

        file.lif         Load a Life 1.05 pattern file on startup.
        --glider         Start with the built-in glider (default).
        --gun            Start with the Gosper glider gun.
        --rpentomino     Start with the R-pentomino.
        --frames  <n>    Auto-run for n generations on startup, then pause.
        --sleep   <ms>   Milliseconds between auto-steps (sets the "Update" speed).
        --help           Print usage and exit without opening a window.

    Unlike the console build, --frames/--sleep don't terminate the program: the
    window stays open and fully interactive once the startup run finishes.
*/

#import "Basic";
#import "Simp";
#import "Window_Creation";
#import "GetRect";
Input :: #import "Input";
#import "Math";
#import "String";   // split, begins_with, ends_with — used by the Life 1.05 loader
#import "File";     // read_entire_file — used by the Life 1.05 loader

// The HashLife engine, pattern loaders, and Life 1.05 parser now live in the
// shared core (also used by the console build and the test harness). The UI
// state, input glue, renderer, and main() follow below.
#load "hashlife_core.jai";

// ===========================================================================
// UI state + glue
// ===========================================================================

#scope_file

// The simulation world: a hash-consed quadtree node, advanced by the engine.
t:      *Node_Table;
gol_id:  node_id;

// View. cam_x/cam_y are the node-frame coordinates of the top-left visible cell.
// While running we auto-follow the pattern; while paused the arrow keys pan.
cam_x: s64 = 0;
cam_y: s64 = 0;
sq_width := 35.0;

// Where to place the camera over an empty world so freehand-drawn cells land at
// comfortably positive node coordinates (the tree only grows toward +x/+y).
BASE :: 1 << 16;

// Cap on cells drawn per axis so a fully zoomed-out window can't make the
// per-frame get_cell sweep pathologically expensive.
MAX_VIEW :: 512;

play := false;
max_time:  float = 1.0;   // seconds between auto-steps
step_gens: float = 1.0;   // generations advanced per step (HashLife jump size)

gen: u64 = 0;             // total generations elapsed (for the --frames limit)
frames_limit: int = 0;    // 0 = run indefinitely; set by --frames, cleared once reached

color_options :: string.[ "yellow", "green", "red", "blue" ];
color_option_choice: s32;
color_theme := Vector4.{1, 1, 0, 1};
prev_open := false;

pattern_options :: string.[ "Glider", "Glider Gun", "R-Pentomino", "Empty" ];
pattern_choice: s32;          // 0 = Glider (matches the pattern loaded in main)
pattern_prev_open := false;

getrect_theme: Overall_Theme;

// Number of whole cells that fit across `extent` pixels (clamped to MAX_VIEW).
visible_count :: (extent: float) -> int {
    n := cast(int)(extent / sq_width) + 1;
    if n < 1         n = 1;
    if n > MAX_VIEW  n = MAX_VIEW;
    return n;
}

// Position the camera so the pattern's bounding box is centred in the viewport.
// On an empty world, centre on BASE so click-to-draw lands at positive coords.
recenter_on_pattern :: (cols: int, rows: int) {
    b := find_bounds(t, gol_id);
    if b.found {
        cx := (b.min_x + b.max_x) / 2;
        cy := (b.min_y + b.max_y) / 2;
        cam_x = cx - cast(s64)(cols / 2);
        cam_y = cy - cast(s64)(rows / 2);
    } else {
        cam_x = cast(s64)BASE - cast(s64)(cols / 2);
        cam_y = cast(s64)BASE - cast(s64)(rows / 2);
    }
}

// Is the world cell (wx, wy) alive? Cells outside the node footprint are dead.
cell_alive :: (wx: s64, wy: s64) -> bool {
    if wx < 0 || wy < 0  return false;
    size := cast(s64)(cast(u64)1 << level_of(gol_id));
    if wx >= size || wy >= size  return false;
    return get_cell(t, gol_id, cast(u64) wx, cast(u64) wy, 0) > 0.0;
}

// Advance the simulation by the chosen number of generations.
step :: () {
    n := cast(u64) step_gens;
    if n < 1  n = 1;
    if lookup(t, gol_id).pop == 0 { play = false; return; }
    gol_id = advance(t, gol_id, n);
    gen   += n;
    if lookup(t, gol_id).pop == 0  play = false;   // stop spinning on a dead board
}

// Toggle the world cell currently displayed at screen (col, row).
// `centre` shifts all live cells by +2^(level-1); we use it to pad the tree
// until the target is representable, moving the camera by the same amount so
// nothing appears to shift.
toggle_cell :: (col: int, row: int) {
    while (cam_x + cast(s64) col < 0 || cam_y + cast(s64) row < 0) && level_of(gol_id) >= 1 {
        shift := cast(s64)(cast(u64)1 << (level_of(gol_id) - 1));
        gol_id = centre(t, gol_id);
        cam_x += shift;
        cam_y += shift;
    }
    nx := cam_x + cast(s64) col;
    ny := cam_y + cast(s64) row;
    if nx < 0 || ny < 0  return;
    alive := cell_alive(nx, ny);
    gol_id = set_cell(t, gol_id, cast(u64) nx, cast(u64) ny, !alive);
}

render_gol :: (win: Window_Type, current_time: float64, mouse_pressed: bool) {
    clear_render_target(.15, .15, .15, 1);
    defer swap_buffers(win);
    width, height := get_render_dimensions(win);

    set_shader_for_color();

    fwidth:  float = xx width;
    fheight: float = xx height;
    offset     := sq_width * 0.025;
    line_width := clamp(sq_width / 30.0, 0.55, 1.0);
    line_color := Vector4.{.35, .35, .35, 1};

    // Draw vertical grid lines.
    sq_x: float = 0.0;
    while sq_x < fwidth {
        immediate_quad(sq_x, fheight, sq_x + line_width, 0, line_color);
        sq_x += sq_width;
    }

    // Draw horizontal grid lines.
    sq_y: float = fheight;
    while sq_y > 0 {
        immediate_quad(0, sq_y, fwidth, sq_y + line_width, line_color);
        sq_y -= sq_width;
    }

    cols := visible_count(fwidth);
    rows := visible_count(fheight);

    // While running, follow the pattern (HashLife re-roots the tree each step,
    // so a fixed camera would lose it).
    if play  recenter_on_pattern(cols, rows);

    // Draw live cells by querying the HashLife quadtree over the visible window.
    for row: 0..rows-1 {
        row_top := fheight - (cast(float) row) * sq_width;   // top edge of this row
        wy := cam_y + cast(s64) row;
        for col: 0..cols-1 {
            col_left := (cast(float) col) * sq_width;
            wx := cam_x + cast(s64) col;
            if cell_alive(wx, wy) {
                immediate_quad(col_left + offset,            row_top - offset,
                               col_left + sq_width - offset, row_top - sq_width + offset,
                               color_theme);
            }
        }
    }

    ui_per_frame_update(win, width, height, current_time);

    ui_interaction := false;
    k := clamp(height * .05, 17.0, 30.0);
    r := get_rect(5, (xx height) - 5.0 - k, 8.5*k, k);

    play_text := ifx play "Stop" else "Simulate";
    if button(r, play_text) {
        play = !play;
        ui_interaction = true;
    }

    r.y -= k + 5.0;
    if button(r, "Next") {
        step();
        ui_interaction = true;
    }

    r.y -= k + 5.0;
    if button(r, "Zoom In") {
        sq_width = min(sq_width + 1.0, xx min(width, height));
        ui_interaction = true;
    }

    r.y -= k + 5.0;
    if button(r, "Zoom Out") {
        sq_width = max(sq_width - 1.0, 1.0);
        ui_interaction = true;
    }

    r.y -= k + 5.0;
    if button(r, "Clear") {
        gol_id = t.off;
        gen = 0;
        recenter_on_pattern(cols, rows);
        ui_interaction = true;
    }

    r.y -= k + 5.0;
    if slider(r, *max_time, 0.01, 3.0, 0.01, *getrect_theme.slider_theme, "Update ", " secs") {
        ui_interaction = true;
    }

    r.y -= k + 5.0;
    if slider(r, *step_gens, 1.0, 512.0, 1.0, *getrect_theme.slider_theme, "Step ", " gens") {
        ui_interaction = true;
    }

    r.y -= k + 5.0;
    choice := color_option_choice;
    changed, dropdown_state := dropdown(r, color_options, *choice);
    if (prev_open == true && dropdown_state.open == false) || (dropdown_state.open == true) {
        ui_interaction = true;
    }
    prev_open = dropdown_state.open;

    r.y -= k + 5.0;
    pchoice := pattern_choice;
    pchanged, pattern_state := dropdown(r, pattern_options, *pchoice);
    if (pattern_prev_open == true && pattern_state.open == false) || (pattern_state.open == true) {
        ui_interaction = true;
    }
    pattern_prev_open = pattern_state.open;

    defer {
        draw_popups();

        if choice != color_option_choice {
            color_option_choice = choice;
            if color_option_choice == 0 {
                color_theme = Vector4.{1, 1, 0, 1};   // yellow
            } else if color_option_choice == 1 {
                color_theme = Vector4.{0, 1, 0, 1};   // green
            } else if color_option_choice == 2 {
                color_theme = Vector4.{1, 0, 0, 1};   // red
            } else if color_option_choice == 3 {
                color_theme = Vector4.{.25, .25, 1, 1}; // blue
            }
        }

        if pchoice != pattern_choice {
            pattern_choice = pchoice;
            if pattern_choice == 0 {
                gol_id = load_glider(t);
            } else if pattern_choice == 1 {
                gol_id = load_glider_gun(t);
            } else if pattern_choice == 2 {
                gol_id = load_r_pentomino(t);
            } else {
                gol_id = t.off;
            }
            play = false;
            gen = 0;
            recenter_on_pattern(cols, rows);
        }
    }

    // Toggle a cell on left-click — only while paused, and only when the click
    // did not land on a UI widget.
    if !ui_interaction && mouse_pressed && !play {
        mouse_x, mouse_y := get_mouse_pointer_position(win, true);
        col := cast(int)(mouse_x / sq_width);
        row := cast(int)((fheight - mouse_y) / sq_width);
        if col >= 0 && col < cols && row >= 0 && row < rows {
            toggle_cell(col, row);
        }
    }
}

setup_getrect_theme :: (theme: Default_Themes) #expand {
    proc := default_theme_procs[theme];
    getrect_theme = proc();
    button_theme := *getrect_theme.button_theme;
    button_theme.label_theme.alignment = .Left;

    slider_theme := *getrect_theme.slider_theme;
    slider_theme.foreground.alignment = .Left;
    set_default_theme(getrect_theme);
}

#scope_export

print_usage :: () {
    print("John Conway's Game of Life — HashLife UI\n");
    print("\n");
    print("Usage: hashlife_ui [options] [file.lif]\n");
    print("\n");
    print("  file.lif         Life 1.05 pattern file to load on startup\n");
    print("\n");
    print("  --glider         Start with the built-in glider (default)\n");
    print("  --gun            Start with the Gosper glider gun\n");
    print("  --rpentomino     Start with the R-pentomino\n");
    print("\n");
    print("  --frames  <n>    Auto-run for n generations on startup, then pause\n");
    print("  --sleep   <ms>   Milliseconds between auto-steps (sets the Update speed)\n");
    print("\n");
    print("  --help           Show this message and exit\n");
    print("\n");
    print("In-window controls: Simulate/Stop, Next, Zoom In/Out, Clear, Update (speed),\n");
    print("Step (generations per step), and colour/pattern dropdowns. Arrow keys pan and\n");
    print("left-click toggles cells while paused.\n");
}

main :: () {
    t = create_table();
    defer free_table(t);

    // --- Parse command-line arguments (ported from the console build) ---
    args := get_command_line_arguments();

    lif_file : string;
    builtin  : string = "glider";   // default when no file is given

    i := 1;   // skip argv[0] (executable name)
    while i < args.count {
        arg := args[i];
        if arg == {
        case "--help";
            print_usage();
            return;                       // exit before opening a window
        case "--glider";
            builtin  = "glider";    lif_file = "";
        case "--gun";
            builtin  = "gun";       lif_file = "";
        case "--rpentomino";
            builtin  = "rpentomino"; lif_file = "";
        case "--frames";
            i += 1;
            if i < args.count {
                v, ok := string_to_int(args[i]);
                if ok && v > 0  frames_limit = v;
            }
        case "--sleep";
            i += 1;
            if i < args.count {
                v, ok := string_to_int(args[i]);
                if ok && v >= 0  max_time = cast(float) v / 1000.0;  // ms -> seconds
            }
        case;
            if ends_with(arg, ".lif") || ends_with(arg, ".LIF") {
                lif_file = arg;
            } else {
                print("Unknown argument: '%'. Run with --help for usage.\n", arg);
            }
        }
        i += 1;
    }

    // --- Window + UI init ---
    win := create_window(800, 600, "John Conway's Game of Life");
    window_width, window_height := get_render_dimensions(win);
    set_render_target(win);
    ui_init();
    setup_getrect_theme(.Grayscale);

    // --- Load the initial pattern ---
    if lif_file.count > 0 {
        id, ok := load_lif_file(t, lif_file);
        if ok {
            gol_id = id;
            print("Loaded: %   (initial pop: %)\n", lif_file, lookup(t, gol_id).pop);
            // Leave pattern_choice at its default; the dropdown is only a runtime
            // switcher for the built-ins and there is no "file" entry.
        } else {
            gol_id = load_glider(t);     // fall back to a glider on load failure
        }
    } else {
        if builtin == {
        case "glider";     gol_id = load_glider(t);      pattern_choice = 0;
        case "gun";        gol_id = load_glider_gun(t);  pattern_choice = 1;
        case "rpentomino"; gol_id = load_r_pentomino(t); pattern_choice = 2;
        }
    }

    recenter_on_pattern(visible_count(xx window_width), visible_count(xx window_height));

    // --frames means "run on startup, then pause": begin playing immediately.
    if frames_limit > 0  play = true;

    dt := 0.01;
    timer := 0.0;
    while eventloop := true {
        Input.update_window_events();
        for Input.get_window_resizes() {
            update_window(it.window);
            if it.window == win {
                window_width  = it.width;
                window_height = it.height;
            }
        }

        mouse_pressed := false;
        for event: Input.events_this_frame {
            if event.type == .QUIT then
                break eventloop;
            getrect_handle_event(event);
            if event.type == .KEYBOARD then {
                if event.key_pressed {
                    if event.key_code == {
                    case .MOUSE_BUTTON_LEFT;
                        mouse_pressed = true;
                    // Arrow keys pan the camera (effective while paused; while
                    // running, auto-follow overwrites the camera each frame).
                    case .ARROW_UP;
                        cam_y -= 1;
                    case .ARROW_DOWN;
                        cam_y += 1;
                    case .ARROW_LEFT;
                        cam_x -= 1;
                    case .ARROW_RIGHT;
                        cam_x += 1;
                    }
                }
            }
        }

        if timer >= max_time then {
            // Respect the --frames startup limit (no-op once it's been cleared).
            if frames_limit == 0 || gen < cast(u64) frames_limit  step();
            timer = 0.0;
        }

        // Once the startup run reaches the frame limit, pause and hand control
        // back to the user. Clearing frames_limit means manual play afterwards
        // isn't immediately re-paused.
        if play && frames_limit > 0 && gen >= cast(u64) frames_limit {
            play = false;
            frames_limit = 0;
            print("Reached frame limit at generation %.\n", gen);
        }

        current_time := seconds_since_init();
        render_gol(win, current_time, mouse_pressed);
        sleep_milliseconds(10);

        if play == true then {
            timer += dt;
        }

        reset_temporary_storage();
    }
}
