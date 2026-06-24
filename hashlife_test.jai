/*
    hashlife_test.jai — robustness, differential, and soak testing for HashLife.

    HashLife is unusually easy to get subtly wrong: hash-consing (join/merge),
    the successor memo cache with one-shot eviction, the binary decomposition in
    advance(), and the pad/crop boundary handling are all classic bug sites. A
    bug there typically does NOT crash — it silently drops or duplicates cells,
    or clips a fast pattern at the boundary. So the strategy here is layered:

      1. DIFFERENTIAL TESTING vs a naive O(n^2) Life simulator.
         The naive sim is obviously correct, so any disagreement is a HashLife
         bug. To make the naive (bounded, dead-border) grid a sound oracle for
         HashLife's infinite plane, we pad it by > number-of-generations on
         every side: information in Life propagates at most 1 cell/generation,
         so the border stays dead and the interior is identical to the plane.

      2. ORACLE-FREE INVARIANTS.
           - lookup(id).pop must equal the actual counted live-cell population.
           - advance(id, N) must equal N applications of advance(id, 1).
         (2) exercises completely different code paths — the big-power-of-two
         successors + memo cache vs the life_4x4 base case — that MUST agree.
         This catches binary-decomposition / clamping / padding bugs with no
         reference implementation at all.

      3. RESIZE / PROBE STRESS.
         Running the same simulation with a tiny initial table (16 slots, forces
         many resize_table + rehash cycles) and a large one must give identical
         results — proving rehashing and linear probing never corrupt content.

      4. RANDOM-SOUP FUZZING. Differential testing on random initial states.

      5. SOAK. A wall-clock-bounded loop that runs millions of generations across
         many fresh patterns, asserting invariants every iteration and doing a
         short-window oracle re-validation periodically, while tracking peak node
         count and table size to surface runaway memory.

    Comparisons are TRANSLATION-INVARIANT: HashLife re-centres the pattern inside
    the quadtree on every step (centre/crop shift the origin), and the naive grid
    has its own fixed origin, so absolute coordinates are meaningless across the
    two. We extract the live-cell set, subtract the min corner, encode each cell
    as a single sorted integer key, and compare the canonical key lists.

    Build:   jai build.jai - test
    Run:     ./gol-test                         (fast suite)
             ./gol-test --fuzz 40                (more soup cases)
             ./gol-test --soak 30                (30-second soak after the suite)
             ./gol-test --seed 12345 --soak 60
*/

#load "hashlife_core.jai";

#import "Basic";
#import "Random";
#import "File";    // write_entire_file (used by the .lif parser tests)

// ---------------------------------------------------------------------------
// Test bookkeeping
// ---------------------------------------------------------------------------

g_pass := 0;
g_fail := 0;

report :: (name: string, ok: bool) {
    if ok { g_pass += 1; print("  [ ok ] %\n", name); }
    else  { g_fail += 1; print("  [FAIL] %\n", name); }
}

iabs :: inline (x: s64) -> s64 { return ifx x < 0 then -x else x; }

Cell :: struct { x: s64; y: s64; }

// ---------------------------------------------------------------------------
// Canonical (translation-invariant) key extraction + comparison
//
// We support a normalized span up to KEY_W in each axis; assert if exceeded so
// a too-large pattern fails loudly rather than aliasing two cells to one key.
// ---------------------------------------------------------------------------

KEY_W :: 1 << 24;

encode_key :: inline (nx: s64, ny: s64) -> s64 { return ny * KEY_W + nx; }

// In-place heapsort: O(n log n) worst case, no recursion, no module dependency.
sift_down :: (a: [] s64, start: s64, n: s64) {
    root := start;
    while true {
        child := 2 * root + 1;
        if child >= n  break;
        if child + 1 < n && a[child] < a[child + 1]  child += 1;
        if a[root] < a[child] {
            a[root], a[child] = a[child], a[root];
            root = child;
        } else {
            break;
        }
    }
}

heapsort_s64 :: (a: [] s64) {
    n := a.count;
    if n < 2  return;
    i := n / 2 - 1;
    while i >= 0 { sift_down(a, i, n); i -= 1; }
    end := n - 1;
    while end > 0 {
        a[0], a[end] = a[end], a[0];
        sift_down(a, 0, end);
        end -= 1;
    }
}

keys_equal :: (a: [] s64, b: [] s64) -> bool {
    if a.count != b.count  return false;
    for i : 0 .. a.count - 1  if a[i] != b[i]  return false;
    return true;
}

// Recursively gather every live cell's (x, y) in the node's own coordinate frame.
collect_cells_rec :: (t: *Node_Table, id: node_id, ox: s64, oy: s64, out: *[..] Cell) {
    n := lookup(t, id);
    if n.pop == 0  return;              // empty subtree — skip (this is the speed win)
    if level_of(id) == 0 {
        array_add(out, .{ox, oy});
        return;
    }
    nn   := n.*;                        // copy before recursing (defensive; no mutation here)
    half := cast(s64)(cast(u64)1 << (level_of(id) - 1));
    collect_cells_rec(t, nn.a, ox,        oy,        out);
    collect_cells_rec(t, nn.b, ox + half, oy,        out);
    collect_cells_rec(t, nn.c, ox,        oy + half, out);
    collect_cells_rec(t, nn.d, ox + half, oy + half, out);
}

node_keys :: (t: *Node_Table, id: node_id) -> [..] s64 {
    cells: [..] Cell;
    defer array_free(cells);
    collect_cells_rec(t, id, 0, 0, *cells);

    keys: [..] s64;
    if cells.count == 0  return keys;

    mnx := cells[0].x;  mny := cells[0].y;
    for c : cells { if c.x < mnx mnx = c.x; if c.y < mny mny = c.y; }
    for c : cells {
        nx := c.x - mnx;  ny := c.y - mny;
        assert(nx >= 0 && ny >= 0 && nx < KEY_W && ny < KEY_W,
               "normalized span exceeds KEY_W; raise KEY_W");
        array_add(*keys, encode_key(nx, ny));
    }
    heapsort_s64(keys);
    return keys;
}

// ---------------------------------------------------------------------------
// Naive reference Life — a dense, double-buffered, dead-border grid.
// Obviously correct; used as the differential oracle.
// ---------------------------------------------------------------------------

Grid :: struct {
    w:     s64;
    h:     s64;
    cells: [] u8;
}

make_grid :: (w: s64, h: s64) -> Grid {
    g: Grid;
    g.w = w;  g.h = h;
    g.cells = NewArray(w * h, u8);     // zero-initialised
    return g;
}

free_grid :: (g: *Grid) {
    array_free(g.cells);
    g.cells.data  = null;
    g.cells.count = 0;
}

grid_get :: inline (g: *Grid, x: s64, y: s64) -> int {
    if x < 0 || y < 0 || x >= g.w || y >= g.h  return 0;   // dead border
    return cast(int) g.cells[y * g.w + x];
}

grid_set :: inline (g: *Grid, x: s64, y: s64, v: u8) {
    if x < 0 || y < 0 || x >= g.w || y >= g.h  return;
    g.cells[y * g.w + x] = v;
}

// One generation, writing into `scratch`, then swapping buffers.
grid_step :: (g: *Grid, scratch: *Grid) {
    for y : 0 .. g.h - 1 {
        for x : 0 .. g.w - 1 {
            n := grid_get(g, x-1, y-1) + grid_get(g, x, y-1) + grid_get(g, x+1, y-1)
               + grid_get(g, x-1, y)                         + grid_get(g, x+1, y)
               + grid_get(g, x-1, y+1) + grid_get(g, x, y+1) + grid_get(g, x+1, y+1);
            alive := g.cells[y * g.w + x] != 0;
            ns : u8 = 0;
            if n == 3 || (n == 2 && alive)  ns = 1;     // B3/S23
            scratch.cells[y * g.w + x] = ns;
        }
    }
    g.cells, scratch.cells = scratch.cells, g.cells;
}

grid_keys :: (g: *Grid) -> [..] s64 {
    keys: [..] s64;
    found := false;
    mnx, mny : s64 = 0, 0;
    for y : 0 .. g.h - 1 {
        for x : 0 .. g.w - 1 {
            if g.cells[y * g.w + x] != 0 {
                if !found { found = true; mnx = x; mny = y; }
                else      { if x < mnx mnx = x; if y < mny mny = y; }
            }
        }
    }
    if !found  return keys;
    for y : 0 .. g.h - 1 {
        for x : 0 .. g.w - 1 {
            if g.cells[y * g.w + x] != 0 {
                nx := x - mnx;  ny := y - mny;
                assert(nx < KEY_W && ny < KEY_W, "grid span exceeds KEY_W");
                array_add(*keys, encode_key(nx, ny));
            }
        }
    }
    heapsort_s64(keys);
    return keys;
}

grid_bbox_min :: (g: *Grid) -> s64, s64, bool {
    found := false;
    mnx, mny : s64 = 0, 0;
    for y : 0 .. g.h - 1 {
        for x : 0 .. g.w - 1 {
            if g.cells[y * g.w + x] != 0 {
                if !found { found = true; mnx = x; mny = y; }
                else      { if x < mnx mnx = x; if y < mny mny = y; }
            }
        }
    }
    return mnx, mny, found;
}

// ---------------------------------------------------------------------------
// Building both representations from the same cell list
// ---------------------------------------------------------------------------

PLACE :: 64;   // keep all HashLife coordinates comfortably positive

build_node_from_cells :: (t: *Node_Table, cells: [] Cell) -> node_id {
    id : node_id = t.off;
    for c : cells  id = set_cell(t, id, cast(u64)(c.x + PLACE), cast(u64)(c.y + PLACE), true);
    return id;
}

// Naive grid sized so the dead border is never reached within `gens` steps.
build_grid_from_cells :: (cells: [] Cell, gens: s64) -> Grid {
    if cells.count == 0  return make_grid(8, 8);
    mnx := cells[0].x;  mny := cells[0].y;  mxx := mnx;  mxy := mny;
    for c : cells {
        if c.x < mnx mnx = c.x;  if c.y < mny mny = c.y;
        if c.x > mxx mxx = c.x;  if c.y > mxy mxy = c.y;
    }
    pad := gens + 2;                                   // light cone grows <= 1 cell/gen
    g := make_grid((mxx - mnx + 1) + 2 * pad, (mxy - mny + 1) + 2 * pad);
    for c : cells  grid_set(*g, (c.x - mnx) + pad, (c.y - mny) + pad, 1);
    return g;
}

// Build a normalized cell list from a pattern-loader (load_glider, etc).
cells_from_builder :: (build: (*Node_Table) -> node_id) -> [..] Cell {
    tt := create_table();
    defer free_table(tt);
    id := build(tt);

    out: [..] Cell;
    collect_cells_rec(tt, id, 0, 0, *out);
    if out.count > 0 {
        mnx := out[0].x;  mny := out[0].y;
        for c : out { if c.x < mnx mnx = c.x; if c.y < mny mny = c.y; }
        for *c : out { c.x -= mnx;  c.y -= mny; }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Layer 1: differential test against the naive oracle
// ---------------------------------------------------------------------------

differential :: (name: string, cells: [] Cell, gens: s64, check_every: s64) -> bool {
    t := create_table();
    defer free_table(t);

    id      := build_node_from_cells(t, cells);
    g       := build_grid_from_cells(cells, gens);
    scratch := make_grid(g.w, g.h);
    defer free_grid(*g);
    defer free_grid(*scratch);

    ok  := true;
    gen : s64 = 0;
    while gen <= gens {
        if gen % check_every == 0 || gen == gens {
            ka := node_keys(t, id);
            kb := grid_keys(*g);
            defer array_free(ka);
            defer array_free(kb);

            hp := lookup(t, id).pop;             // invariant: stored pop == counted pop
            if cast(s64) hp != ka.count {
                print("    [%] pop invariant broken at gen %: node.pop=% counted=%\n",
                      name, gen, hp, ka.count);
                ok = false;
                break;
            }
            if !keys_equal(ka, kb) {
                print("    [%] DIVERGENCE at gen %: hashlife pop=% naive pop=%\n",
                      name, gen, ka.count, kb.count);
                ok = false;
                break;
            }
        }
        if gen == gens  break;
        id = advance(t, id, 1);
        grid_step(*g, *scratch);
        gen += 1;
    }
    return ok;
}

// ---------------------------------------------------------------------------
// Layer 2a: still lifes — invariant under any number of steps, including huge
//           power-of-two jumps (exercises deep successor memoisation).
// ---------------------------------------------------------------------------

still_life_test :: (name: string, cells: [] Cell, expect_pop: s64) {
    t := create_table();
    defer free_table(t);

    id := build_node_from_cells(t, cells);
    k0 := node_keys(t, id);
    defer array_free(k0);

    ok := cast(s64) lookup(t, id).pop == expect_pop && k0.count == expect_pop;

    steps := u64.[1, 2, 3, 7, 50, 1 << 10, 1 << 20, 1 << 30];
    for s : steps {
        r  := advance(t, id, s);
        kr := node_keys(t, r);
        defer array_free(kr);
        if !keys_equal(kr, k0)                       ok = false;
        if cast(s64) lookup(t, r).pop != expect_pop  ok = false;
    }
    report(tprint("still life: % (pop %)", name, expect_pop), ok);
}

// ---------------------------------------------------------------------------
// Layer 2b: oscillators — return to the same phase after `period` steps, but
//           must actually change after a single step (proves it isn't frozen).
// ---------------------------------------------------------------------------

oscillator_test :: (name: string, cells: [] Cell, period: s64) {
    t := create_table();
    defer free_table(t);

    id := build_node_from_cells(t, cells);
    k0 := node_keys(t, id);
    defer array_free(k0);

    p1 := advance(t, id, cast(u64) period);
    k1 := node_keys(t, p1);
    defer array_free(k1);

    p2 := advance(t, id, cast(u64)(period * 2));
    k2 := node_keys(t, p2);
    defer array_free(k2);

    s1 := advance(t, id, 1);
    ks := node_keys(t, s1);
    defer array_free(ks);

    same_after_period   := keys_equal(k0, k1) && keys_equal(k0, k2);
    evolves_after_1_step := period == 1 || !keys_equal(k0, ks);

    report(tprint("oscillator: % (period %)", name, period),
           same_after_period && evolves_after_1_step);
}

// ---------------------------------------------------------------------------
// Layer 2c: spaceships — shape returns after `period`, and the displacement
//           (measured in the naive grid, the trusted fixed lattice) matches the
//           ship's known signature. We first confirm the pattern really is a
//           clean spaceship in the *naive* sim, so a typo in our pattern can't
//           masquerade as a HashLife bug. kind 0 = diagonal c/4 (|dx|=|dy|=1),
//           kind 1 = orthogonal c/2 (one axis 2, other 0).
// ---------------------------------------------------------------------------

spaceship_test :: (name: string, cells: [] Cell, period: s64, kind: int) {
    // --- naive: validate the pattern AND measure ground-truth displacement ---
    g  := build_grid_from_cells(cells, period + 2);
    sc := make_grid(g.w, g.h);
    defer free_grid(*g);
    defer free_grid(*sc);

    x0, y0, f0 := grid_bbox_min(*g);
    kn0 := grid_keys(*g);
    defer array_free(kn0);

    for 1 .. period  grid_step(*g, *sc);

    x1, y1, f1 := grid_bbox_min(*g);
    kn1 := grid_keys(*g);
    defer array_free(kn1);

    dx := iabs(x1 - x0);  dy := iabs(y1 - y0);
    disp_sig := ifx kind == 0
        then (dx == 1 && dy == 1)
        else ((dx == 2 && dy == 0) || (dx == 0 && dy == 2));
    naive_clean := f0 && f1 && keys_equal(kn0, kn1) && disp_sig;

    // --- hashlife: shape returns after one period ---
    t := create_table();
    defer free_table(t);
    id := build_node_from_cells(t, cells);
    kh0 := node_keys(t, id);
    defer array_free(kh0);
    php := advance(t, id, cast(u64) period);
    khp := node_keys(t, php);
    defer array_free(khp);
    hashlife_shape_ok := keys_equal(kh0, khp);

    // --- the two systems agree generation-by-generation ---
    diff_ok := differential(name, cells, period * 4, 1);

    if !naive_clean
        print("    [%] note: pattern failed naive spaceship self-check (dx=%, dy=%)\n", name, dx, dy);

    report(tprint("spaceship: % (period %)", name, period),
           naive_clean && hashlife_shape_ok && diff_ok);
}

// ---------------------------------------------------------------------------
// Layer 2d: advance(id, N) must equal N applications of advance(id, 1).
// Oracle-free; binds the big-power-of-two successor path to the base case.
// ---------------------------------------------------------------------------

advance_consistency_test :: (name: string, cells: [] Cell) {
    t := create_table();
    defer free_table(t);
    id := build_node_from_cells(t, cells);

    Ns := s64.[1, 2, 3, 4, 5, 7, 8, 15, 16, 31, 64, 100, 257];
    all_ok := true;
    for N : Ns {
        big := advance(t, id, cast(u64) N);
        kb  := node_keys(t, big);

        step := id;
        for 1 .. N  step = advance(t, step, 1);
        ks := node_keys(t, step);

        eq := keys_equal(kb, ks);
        array_free(kb);
        array_free(ks);
        if !eq {
            print("    [%] advance(%) != % x advance(1)\n", name, N, N);
            all_ok = false;
        }
    }
    report(tprint("advance == repeated-step: %", name), all_ok);
}

// ---------------------------------------------------------------------------
// Layer 3: determinism + resize/probe stress.
// Identical results from huge and tiny initial tables prove rehashing and
// linear probing never corrupt or lose content.
// ---------------------------------------------------------------------------

run_to_keys :: (cells: [] Cell, gens: s64, init_size: u64) -> [..] s64 {
    t := create_table(init_size);
    defer free_table(t);
    id := build_node_from_cells(t, cells);
    for 1 .. gens  id = advance(t, id, 1);
    return node_keys(t, id);            // computed before the deferred free
}

determinism_test :: (name: string, cells: [] Cell, gens: s64) {
    a := run_to_keys(cells, gens, 1 << 16);
    b := run_to_keys(cells, gens, 1 << 16);   // same config, must repeat exactly
    c := run_to_keys(cells, gens, 16);        // forces many resize_table cycles
    defer array_free(a);
    defer array_free(b);
    defer array_free(c);
    report(tprint("determinism + resize stress: %", name),
           keys_equal(a, b) && keys_equal(a, c));
}

// ---------------------------------------------------------------------------
// set_cell / get_cell round-trip (quadtree growth + addressing)
// ---------------------------------------------------------------------------

set_get_roundtrip_test :: () {
    t := create_table();
    defer free_table(t);

    seen: [..] Cell;   defer array_free(seen);
    used: [..] s64;    defer array_free(used);   // encoded keys, for distinctness
    id : node_id = t.off;

    COUNT :: 300;
    n := 0;
    while n < COUNT {
        x   := cast(s64)(random_get() % 256);
        y   := cast(s64)(random_get() % 256);
        key := y * 256 + x;
        dup := false;
        for k : used  if k == key { dup = true; break; }
        if dup  continue;
        array_add(*used, key);
        array_add(*seen, .{x, y});
        id = set_cell(t, id, cast(u64) x, cast(u64) y, true);
        n += 1;
    }

    ok := cast(s64) lookup(t, id).pop == COUNT;
    for c : seen  if get_cell(t, id, cast(u64) c.x, cast(u64) c.y, 0) <= 0.0 { ok = false; break; }

    // Confirm coordinates we never set read back dead.
    misses := 0;
    while misses < 200 {
        x   := cast(s64)(random_get() % 256);
        y   := cast(s64)(random_get() % 256);
        key := y * 256 + x;
        present := false;
        for k : used  if k == key { present = true; break; }
        if present  continue;
        if get_cell(t, id, cast(u64) x, cast(u64) y, 0) > 0.0  ok = false;
        misses += 1;
    }

    report("set/get round-trip (300 distinct cells)", ok);
}

// ---------------------------------------------------------------------------
// Layer 4: random-soup fuzzing
// ---------------------------------------------------------------------------

soup_cells :: (w: s64, h: s64, density_pct: int) -> [..] Cell {
    cells: [..] Cell;
    for y : 0 .. h - 1 {
        for x : 0 .. w - 1 {
            if cast(int)(random_get() % 100) < density_pct  array_add(*cells, .{x, y});
        }
    }
    return cells;
}

fuzz_differential :: (count: int, gens: s64) {
    fails := 0;
    for iter : 0 .. count - 1 {
        w    := cast(s64)(16 + random_get() % 32);
        h    := cast(s64)(16 + random_get() % 32);
        dens := cast(int)(15 + random_get() % 35);
        cells := soup_cells(w, h, dens);
        defer array_free(cells);

        every := ifx gens / 4 < 1 then cast(s64) 1 else gens / 4;
        ok := differential(tprint("soup#%", iter), cells, gens, every);
        if !ok {
            fails += 1;
            print("    soup#% FAILED  (%x%, density %%%, %% gens) — rerun with --seed to reproduce\n",
                  iter, w, h, dens, gens);
        }
        reset_temporary_storage();
    }
    report(tprint("random-soup fuzz: % cases x % gens", count, gens), fails == 0);
}

// ---------------------------------------------------------------------------
// LIF parser robustness (malformed / edge-case input must not crash)
// ---------------------------------------------------------------------------

lif_parser_test :: () {
    t := create_table();
    defer free_table(t);
    tmp :: "____hashlife_test_tmp.lif";

    write_entire_file(tmp, "#Life 1.05\n#D test\n#P 0 0\n.*.\n.*.\n.*.\n");
    _,  valid_ok := load_lif_file(t, tmp);

    write_entire_file(tmp, "not a life file\n***\n");
    _,  ok_bad := load_lif_file(t, tmp);

    write_entire_file(tmp, "#Life 1.05\n");
    e3, ok_empty := load_lif_file(t, tmp);

    write_entire_file(tmp, "#Life 1.05\n#P -5 -7\n*\n");
    e4, ok_neg := load_lif_file(t, tmp);

    write_entire_file(tmp, "#Life 1.05\n#P 3\n*\n");          // truncated #P
    _,  _ = load_lif_file(t, tmp);                             // must merely survive

    write_entire_file(tmp, "#Life 1.05\n#P 100000 -250000\n*\n");
    e6, ok_big := load_lif_file(t, tmp);

    report("lif: valid pattern loads",        valid_ok);
    report("lif: rejects missing header",     !ok_bad);
    report("lif: empty pattern -> pop 0",     ok_empty && lookup(t, e3).pop == 0);
    report("lif: negative coords -> pop 1",   ok_neg   && lookup(t, e4).pop == 1);
    report("lif: truncated #P survived",      true);   // reaching here == no crash
    report("lif: large coords -> pop 1",      ok_big   && lookup(t, e6).pop == 1);
}

// ---------------------------------------------------------------------------
// Layer 5: SOAK
//
// Long-running loop bounded by wall clock. Each iteration: a fresh table, a
// random soup, advanced for many generations using random (1..8) step sizes
// (which mixes base-case and power-of-two successor paths). Every iteration we
// assert the pop invariant; periodically we export the current live state into
// a small dense grid and re-validate a short window against the naive oracle.
// Tracks peak node count and table size to surface runaway memory growth.
// ---------------------------------------------------------------------------

soak :: (seconds: float, base_seed: u32) {
    print("\n=== SOAK (target %.0f s) ===\n", seconds);
    start := seconds_since_init();

    iters         := 0;
    revalidations := 0;
    deaths        := 0;
    peak_nodes : u64 = 0;
    peak_slots : s64 = 0;
    seed := base_seed;

    while (seconds_since_init() - start) < seconds {
        random_seed(seed);
        seed += 1;

        t := create_table();

        w    := cast(s64)(16 + random_get() % 48);
        h    := cast(s64)(16 + random_get() % 48);
        dens := cast(int)(10 + random_get() % 40);
        cells := soup_cells(w, h, dens);
        id := build_node_from_cells(t, cells);
        array_free(cells);

        total := cast(s64)(200 + random_get() % 1200);
        done  : s64 = 0;
        while done < total {
            step := cast(u64)(1 + random_get() % 8);
            if cast(s64) step > total - done  step = cast(u64)(total - done);
            id = advance(t, id, step);
            done += cast(s64) step;
            if lookup(t, id).pop == 0 { deaths += 1; break; }
        }

        // Invariant: stored population equals the actual counted population.
        snap: [..] Cell;
        collect_cells_rec(t, id, 0, 0, *snap);
        if cast(s64) lookup(t, id).pop != snap.count {
            print("SOAK FAIL: pop invariant broken (iter %, seed %, node.pop=%, counted=%)\n",
                  iters, seed - 1, lookup(t, id).pop, snap.count);
            array_free(snap);
            free_table(t);
            g_fail += 1;
            return;
        }

        // Periodic short-window oracle re-validation (only when span is small
        // enough to fit a dense grid cheaply).
        if iters % 16 == 0 && snap.count > 0 {
            mnx := snap[0].x;  mny := snap[0].y;  mxx := mnx;  mxy := mny;
            for c : snap {
                if c.x < mnx mnx = c.x;  if c.y < mny mny = c.y;
                if c.x > mxx mxx = c.x;  if c.y > mxy mxy = c.y;
            }
            span := mxx - mnx;
            if (mxy - mny) > span  span = mxy - mny;
            if span < 400 {
                for *c : snap { c.x -= mnx;  c.y -= mny; }
                ok := differential(tprint("soak-revalidate#%", iters), snap, 16, 4);
                revalidations += 1;
                if !ok {
                    print("SOAK FAIL: oracle mismatch (iter %, seed %)\n", iters, seed - 1);
                    array_free(snap);
                    free_table(t);
                    g_fail += 1;
                    return;
                }
            }
        }
        array_free(snap);

        if t.count       > peak_nodes  peak_nodes = t.count;
        if t.index.count > peak_slots  peak_slots = t.index.count;

        free_table(t);
        iters += 1;
        reset_temporary_storage();
    }

    elapsed := seconds_since_init() - start;
    print("soak: % iterations, % oracle re-validations, % die-outs, %.1f s\n",
          iters, revalidations, deaths, elapsed);
    print("soak: peak node count = %, peak table slots = % (~%.1f MB index)\n",
          peak_nodes, peak_slots, cast(float)(peak_slots * size_of(Node)) / (1024.0 * 1024.0));
    print("soak: completed with no crashes and no invariant violations.\n");
    g_pass += 1;
}

// ---------------------------------------------------------------------------
// Test patterns (small ones hand-encoded; gun/r-pentomino reuse the loaders)
// ---------------------------------------------------------------------------

BLOCK   :: Cell.[ .{0,0}, .{1,0}, .{0,1}, .{1,1} ];
BEEHIVE :: Cell.[ .{1,0}, .{2,0}, .{0,1}, .{3,1}, .{1,2}, .{2,2} ];
LOAF    :: Cell.[ .{1,0}, .{2,0}, .{0,1}, .{3,1}, .{1,2}, .{3,2}, .{2,3} ];
BOAT    :: Cell.[ .{0,0}, .{1,0}, .{0,1}, .{2,1}, .{1,2} ];
TUB     :: Cell.[ .{1,0}, .{0,1}, .{2,1}, .{1,2} ];

BLINKER :: Cell.[ .{0,0}, .{1,0}, .{2,0} ];
TOAD    :: Cell.[ .{1,0}, .{2,0}, .{3,0}, .{0,1}, .{1,1}, .{2,1} ];
BEACON  :: Cell.[ .{0,0}, .{1,0}, .{0,1}, .{1,1}, .{2,2}, .{3,2}, .{2,3}, .{3,3} ];

GLIDER  :: Cell.[ .{1,0}, .{2,1}, .{0,2}, .{1,2}, .{2,2} ];
LWSS    :: Cell.[ .{1,0}, .{4,0}, .{0,1}, .{0,2}, .{4,2}, .{0,3}, .{1,3}, .{2,3}, .{3,3} ];

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

main :: () {
    args := get_command_line_arguments();

    soak_seconds : float = 0.0;
    fuzz_count   : int   = 12;
    seed         : u32   = 1;

    i := 1;
    while i < args.count {
        a := args[i];
        if a == {
        case "--soak"; i += 1; if i < args.count { v, ok := string_to_int(args[i]); if ok soak_seconds = cast(float) v; }
        case "--fuzz"; i += 1; if i < args.count { v, ok := string_to_int(args[i]); if ok fuzz_count   = v; }
        case "--seed"; i += 1; if i < args.count { v, ok := string_to_int(args[i]); if ok seed         = cast(u32) v; }
        case "--help";
            print("Usage: gol-test [--fuzz N] [--soak SECONDS] [--seed N]\n");
            return;
        case;
            print("Unknown argument: '%'\n", a);
        }
        i += 1;
    }

    random_seed(seed);

    print("HashLife test & robustness suite  (seed %)\n", seed);
    print("==========================================\n\n");

    print("-- Quadtree primitives --\n");
    set_get_roundtrip_test();

    print("\n-- Still lifes --\n");
    still_life_test("block",   BLOCK,   4);
    still_life_test("beehive", BEEHIVE, 6);
    still_life_test("loaf",    LOAF,    7);
    still_life_test("boat",    BOAT,    5);
    still_life_test("tub",     TUB,     4);

    print("\n-- Oscillators --\n");
    oscillator_test("blinker", BLINKER, 2);
    oscillator_test("toad",    TOAD,    2);
    oscillator_test("beacon",  BEACON,  2);

    print("\n-- Spaceships --\n");
    spaceship_test("glider", GLIDER, 4, 0);
    spaceship_test("LWSS",   LWSS,   4, 1);

    print("\n-- Differential vs naive oracle (known patterns) --\n");
    {
        gl := cells_from_builder(load_glider);       defer array_free(gl);
        gn := cells_from_builder(load_glider_gun);   defer array_free(gn);
        rp := cells_from_builder(load_r_pentomino);  defer array_free(rp);
        report("glider, 40 gens",            differential("glider",     gl, 40,  1));
        report("Gosper glider gun, 90 gens", differential("gun",        gn, 90,  1));
        report("R-pentomino, 200 gens",      differential("r-pentomino", rp, 200, 1));
    }

    print("\n-- Internal consistency (oracle-free) --\n");
    advance_consistency_test("blinker", BLINKER);
    {
        rp := cells_from_builder(load_r_pentomino);  defer array_free(rp);
        gn := cells_from_builder(load_glider_gun);   defer array_free(gn);
        advance_consistency_test("r-pentomino", rp);
        determinism_test("r-pentomino, 120 gens", rp, 120);
        determinism_test("glider gun, 90 gens",    gn, 90);
    }

    print("\n-- Random-soup fuzzing --\n");
    fuzz_differential(fuzz_count, 60);

    print("\n-- LIF parser robustness --\n");
    lif_parser_test();

    if soak_seconds > 0.0  soak(soak_seconds, seed + 1000);

    print("\n==========================================\n");
    print("RESULT: % passed, % failed\n", g_pass, g_fail);
    if g_fail > 0  exit(1);
}
