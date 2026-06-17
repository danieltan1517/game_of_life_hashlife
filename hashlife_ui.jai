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

// ===========================================================================
// HashLife engine (from hashlife_console.jai; the console renderer and the
// command-line main have been removed — the Life 1.05 loader is kept below
// so the UI can honour a file.lif argument).
// ===========================================================================

// ---------------------------------------------------------------------------
// Bit-packing helpers for node_id
// ---------------------------------------------------------------------------

HASH_BITS  :: 46;
HASH_MASK  :: (0xFFFF_FFFF_FFFF >> (64 - HASH_BITS));   // lower 46 bits

node_id :: u64;
UNUSED  : node_id : 0;  // sentinel: slot is empty in the open-address table

level_of :: inline (id: node_id) -> u64 {
    return (id >> 46) & 0xFFFF;
}

is_zero :: inline (id: node_id) -> bool {
    return (id >> 62) & 1 != 0;
}

is_marked :: inline (id: node_id) -> bool {
    return (id >> 63) & 1 != 0;
}

mark :: inline (id: node_id) -> node_id {
    return id | (cast(u64)1 << 63);
}

unmark :: inline (id: node_id) -> node_id {
    return id & ~(cast(u64)1 << 63);
}

// ---------------------------------------------------------------------------
// Node struct: one entry in the open-address hash table
// ---------------------------------------------------------------------------

Node :: struct {
    id:  node_id;   // the canonical key; UNUSED (0) means empty slot
    a:   node_id;   // top-left  quadrant (or 0 for leaf)
    b:   node_id;   // top-right quadrant
    c:   node_id;   // bot-left  quadrant
    d:   node_id;   // bot-right quadrant
    pop: u64;       // total living cells in this subtree
}

Next_Entry :: struct {
    from: node_id;
    to:   node_id;
    j:    u64;
}

// ---------------------------------------------------------------------------
// Node table
// ---------------------------------------------------------------------------

Node_Table :: struct {
    index: [] Node;       // open-address table, power-of-two size
    next:  [] Next_Entry; // successor cache, same size
    count: u64;
    off:   node_id;       // canonical "dead leaf" id
    on:    node_id;       // canonical "live leaf" id
}

// ---------------------------------------------------------------------------
// SplitMix64 — fast, good mixing function (same as C reference)
// ---------------------------------------------------------------------------

mix64 :: inline (x: u64) -> u64 {
    x += cast(u64)0x9e3779b97f4a7c15;
    x  = (x ^ (x >> 30)) * cast(u64)0xbf58476d1ce4e5b9;
    x  = (x ^ (x >> 27)) * cast(u64)0x94d049bb133111eb;
    return x ^ (x >> 31);
}

hash_quad :: inline (a: u64, b: u64, c: u64, d: u64) -> u64 {
    h : u64 = cast(u64)0x243f6a8885a308d3;
    h ^= mix64(a + cast(u64)0x9e3779b97f4a7c15 * 1);
    h  = (h * cast(u64)0x9e3779b97f4a7c15) ^ (h >> 32);
    h ^= mix64(b + cast(u64)0x9e3779b97f4a7c15 * 2);
    h  = (h * cast(u64)0x9e3779b97f4a7c15) ^ (h >> 32);
    h ^= mix64(c + cast(u64)0x9e3779b97f4a7c15 * 3);
    h  = (h * cast(u64)0x9e3779b97f4a7c15) ^ (h >> 32);
    h ^= mix64(d + cast(u64)0x9e3779b97f4a7c15 * 4);
    h  = (h * cast(u64)0x9e3779b97f4a7c15) ^ (h >> 32);
    return mix64(h);
}

// ---------------------------------------------------------------------------
// Table creation / teardown
// ---------------------------------------------------------------------------

create_table :: (initial_size: u64 = 1 << 16) -> *Node_Table {
    t := New(Node_Table);
    sz := initial_size;
    if sz < 16 then sz = 16;
    t.index = NewArray(cast(s64) sz, Node);
    t.next  = NewArray(cast(s64) sz, Next_Entry);

    t.off = (cast(u64)1 << 62) | (cast(u64)0 << 46) | (mix64(0) & HASH_MASK);
    t.on  = (cast(u64)0 << 62) | (cast(u64)0 << 46) | (mix64(1) & HASH_MASK);

    mask := cast(u64) t.index.count - 1;
    t.index[cast(s64)(t.off & mask)] = .{id=t.off, a=0, b=0, c=0, d=0, pop=0};
    t.index[cast(s64)(t.on  & mask)] = .{id=t.on,  a=0, b=0, c=0, d=0, pop=1};
    t.count = 2;
    return t;
}

free_table :: (t: *Node_Table) {
    array_free(t.index);
    array_free(t.next);
    free(t);
}

// ---------------------------------------------------------------------------
// Open-address lookup (linear probing)
// ---------------------------------------------------------------------------

lookup :: (t: *Node_Table, id: node_id) -> *Node {
    assert(t != null);
    mask  := cast(u64) t.index.count - 1;
    index := id;
    n     := *t.index[cast(s64)(index & mask)];
    while n.id != UNUSED && n.id != id {
        index += 1;
        n = *t.index[cast(s64)(index & mask)];
    }
    return n;
}

// ---------------------------------------------------------------------------
// Successor cache (one-shot eviction)
// ---------------------------------------------------------------------------

lookup_next :: (t: *Node_Table, from: node_id, j: u64) -> node_id {
    assert(t != null);
    mask := cast(u64) t.next.count - 1;
    h    := hash_quad(from, j, from, j) & mask;
    e    := *t.next[cast(s64)h];
    if e.from == from && e.j == j  return e.to;
    return UNUSED;
}

cache_next :: (t: *Node_Table, from: node_id, to: node_id, j: u64) {
    assert(t != null);
    mask := cast(u64) t.next.count - 1;
    h    := hash_quad(from, j, from, j) & mask;
    e    := *t.next[cast(s64)h];
    e.from = from;
    e.to   = to;
    e.j    = j;
}

// ---------------------------------------------------------------------------
// Compute the canonical id for the parent of four quadrants
// ---------------------------------------------------------------------------

merge :: (a: node_id, b: node_id, c: node_id, d: node_id) -> node_id {
    h       := hash_quad(a, b, c, d);
    level   := ((a >> 46) + 1) & 0xFFFF;
    all_zero := is_zero(a) && is_zero(b) && is_zero(c) && is_zero(d);
    if all_zero {
        return (cast(u64)1 << 62) | (cast(u64)(level << 46)) | (mix64(level) & HASH_MASK);
    } else {
        return (cast(u64)0 << 62) | (cast(u64)(level << 46)) | (h & HASH_MASK);
    }
}

// ---------------------------------------------------------------------------
// Table resize (double, rehash)
// ---------------------------------------------------------------------------

resize_table :: (t: *Node_Table) {
    assert(t != null);
    new_sz   := cast(u64) t.index.count * 2;
    new_idx  := NewArray(cast(s64) new_sz, Node);
    new_nxt  := NewArray(cast(s64) new_sz, Next_Entry);
    new_mask := new_sz - 1;

    for i : 0 .. t.index.count - 1 {
        old := t.index[i];
        if old.id != UNUSED {
            pos := old.id & new_mask;
            while new_idx[cast(s64)pos].id != 0  pos = (pos + 1) & new_mask;
            new_idx[cast(s64)pos] = old;
        }
    }

    array_free(t.index);
    array_free(t.next);
    t.index = new_idx;
    t.next  = new_nxt;
}

// ---------------------------------------------------------------------------
// join: intern four child nodes, returning the parent node_id
// ---------------------------------------------------------------------------

join :: (t: *Node_Table, a: node_id, b: node_id, c: node_id, d: node_id) -> node_id {
    assert(t != null);
    hash := merge(a, b, c, d);
    n    := lookup(t, hash);

    // Linear probe past "dopplegangers" (different content, same hash bits)
    while n.id != UNUSED {
        if n.a == a && n.b == b && n.c == c && n.d == d  return n.id;
        new_hash := mix64(hash);
        hash    ^= new_hash & HASH_MASK;
        n        = lookup(t, hash);
    }

    // Insert new node
    na := lookup(t, a);
    nb := lookup(t, b);
    nc := lookup(t, c);
    nd := lookup(t, d);

    n.id  = hash;
    n.a   = a;
    n.b   = b;
    n.c   = c;
    n.d   = d;
    n.pop = na.pop + nb.pop + nc.pop + nd.pop;
    t.count += 1;

    if t.count * 4 >= cast(u64) t.index.count  resize_table(t);
    return hash;
}

// ---------------------------------------------------------------------------
// Zero-filled node of given level
// ---------------------------------------------------------------------------

get_zero :: (t: *Node_Table, k: u64) -> node_id {
    if k == 0  return t.off;

    z := (cast(u64)1 << 62) | (k << 46) | (mix64(k) & HASH_MASK);
    if lookup(t, z).id == z  return z;
    zk := get_zero(t, k - 1);
    return join(t, zk, zk, zk, zk);
}

// ---------------------------------------------------------------------------
// Expand a node to the next level by padding with zeros around it
// ---------------------------------------------------------------------------

centre :: (t: *Node_Table, m_h: node_id) -> node_id {
    assert(level_of(m_h) >= 1, "centre called on a leaf node (level 0) — tree was not grown first");
    m := lookup(t, m_h).*;
    z := get_zero(t, level_of(m_h) - 1);
    na := join(t, z,   z,   z,   m.a);
    nb := join(t, z,   z,   m.b, z  );
    nc := join(t, z,   m.c, z,   z  );
    nd := join(t, m.d, z,   z,   z  );
    return join(t, na, nb, nc, nd);
}

// ---------------------------------------------------------------------------
// Base-case: apply Life rule to the 4x4 block, returning the central 2x2
// ---------------------------------------------------------------------------

base_life :: inline (
    a: node_id, b: node_id, c: node_id,
    d: node_id, e: node_id, f: node_id,
    g: node_id, h: node_id, i: node_id,
    ON: node_id, OFF: node_id) -> node_id
{
    pop := cast(int)(a == ON) + cast(int)(b == ON) + cast(int)(c == ON)
         + cast(int)(d == ON)                      + cast(int)(f == ON)
         + cast(int)(g == ON) + cast(int)(h == ON) + cast(int)(i == ON);
    if pop == 3 || (pop == 2 && e == ON)  return ON;
    return OFF;
}

life_4x4 :: (t: *Node_Table, id: node_id) -> node_id {
    n := lookup(t, id).*;
    a := lookup(t, n.a).*;
    b := lookup(t, n.b).*;
    c := lookup(t, n.c).*;
    d := lookup(t, n.d).*;

    on  := t.on;
    off := t.off;

    na := base_life(a.a, a.b, b.a, a.c, a.d, b.c, c.a, c.b, d.a, on, off);
    nb := base_life(a.b, b.a, b.b, a.d, b.c, b.d, c.b, d.a, d.b, on, off);
    nc := base_life(a.c, a.d, b.c, c.a, c.b, d.a, c.c, c.d, d.c, on, off);
    nd := base_life(a.d, b.c, b.d, c.b, d.a, d.b, c.d, d.c, d.d, on, off);

    return join(t, na, nb, nc, nd);
}

// ---------------------------------------------------------------------------
// successor: advance node by 2^j steps (clamped to max 2^(level-2))
// ---------------------------------------------------------------------------

sucjoin :: (t: *Node_Table, a: node_id, b: node_id, c: node_id, d: node_id, j: u64) -> node_id {
    return successor(t, join(t, a, b, c, d), j);
}

successor :: (t: *Node_Table, id: node_id, j: u64) -> node_id {
    level  := level_of(id);
    slevel := cast(s64) level;
    jj     := j;
    if cast(s64)jj > slevel - 2  jj = cast(u64)(slevel - 2);

    if is_zero(id)  return lookup(t, id).a;

    next := lookup_next(t, id, jj);
    if next != UNUSED  return next;

    if level == 2 {
        next = life_4x4(t, id);
        cache_next(t, id, next, jj);
        return next;
    }

    n := lookup(t, id).*;
    a := lookup(t, n.a).*;
    b := lookup(t, n.b).*;
    c := lookup(t, n.c).*;
    d := lookup(t, n.d).*;

    c1 := successor(t, a.id, jj);
    c2 := sucjoin  (t, a.b, b.a, a.d, b.c, jj);
    c3 := successor(t, b.id, jj);
    c4 := sucjoin  (t, a.c, a.d, c.a, c.b, jj);
    c5 := sucjoin  (t, a.d, b.c, c.b, d.a, jj);
    c6 := sucjoin  (t, b.c, b.d, d.a, d.b, jj);
    c7 := successor(t, c.id, jj);
    c8 := sucjoin  (t, c.b, d.a, c.d, d.c, jj);
    c9 := successor(t, d.id, jj);

    if jj < cast(u64)(slevel - 2) {
        c1n := lookup(t, c1).*;
        c2n := lookup(t, c2).*;
        c3n := lookup(t, c3).*;
        c4n := lookup(t, c4).*;
        c5n := lookup(t, c5).*;
        c6n := lookup(t, c6).*;
        c7n := lookup(t, c7).*;
        c8n := lookup(t, c8).*;
        c9n := lookup(t, c9).*;

        next = join(t,
                    join(t, c1n.d, c2n.c, c4n.b, c5n.a),
                    join(t, c2n.d, c3n.c, c5n.b, c6n.a),
                    join(t, c4n.d, c5n.c, c7n.b, c8n.a),
                    join(t, c5n.d, c6n.c, c8n.b, c9n.a));
    } else {
        next = join(t,
                    sucjoin(t, c1, c2, c4, c5, jj),
                    sucjoin(t, c2, c3, c5, c6, jj),
                    sucjoin(t, c4, c5, c7, c8, jj),
                    sucjoin(t, c5, c6, c8, c9, jj));
    }

    cache_next(t, id, next, jj);
    return next;
}

// ---------------------------------------------------------------------------
// Trimming helpers
// ---------------------------------------------------------------------------

inner :: (t: *Node_Table, id: node_id) -> node_id {
    n := lookup(t, id).*;
    a := lookup(t, n.a).*;
    b := lookup(t, n.b).*;
    c := lookup(t, n.c).*;
    d := lookup(t, n.d).*;
    return join(t, a.d, b.c, c.b, d.a);
}

is_padded :: (t: *Node_Table, id: node_id) -> bool {
    n  := lookup(t, id);
    na := lookup(t, n.a);
    nb := lookup(t, n.b);
    nc := lookup(t, n.c);
    nd := lookup(t, n.d);
    return na.pop == lookup(t, na.d).pop
        && nb.pop == lookup(t, nb.c).pop
        && nc.pop == lookup(t, nc.b).pop
        && nd.pop == lookup(t, nd.a).pop;
}

crop :: (t: *Node_Table, id: node_id) -> node_id {
    while level_of(id) > 3 && is_padded(t, id)  id = inner(t, id);
    return id;
}

pad :: (t: *Node_Table, id: node_id) -> node_id {
    while level_of(id) < 3 || is_padded(t, id)  id = centre(t, id);
    return id;
}

// ---------------------------------------------------------------------------
// advance: simulate forward by exactly `steps` generations.
// ---------------------------------------------------------------------------

advance :: (t: *Node_Table, id: node_id, steps: u64) -> node_id {
    j   : u64 = 0;
    rem := steps;
    while rem > 0 {
        if rem & 1 {
            while level_of(id) < j + 2                                      id = centre(t, id);
            while (cast(u64)1 << (level_of(id) - 2)) < (cast(u64)1 << j)  id = centre(t, id);
            id = centre(t, id);   // zero-border ring
            id = centre(t, id);   // extra room for translation
            id = successor(t, id, j);
        }
        j   += 1;
        rem >>= 1;
    }
    return crop(t, id);
}

// ---------------------------------------------------------------------------
// Cell set / get
// ---------------------------------------------------------------------------

set_cell :: (t: *Node_Table, id: node_id, x: u64, y: u64, state: bool) -> node_id {
    if level_of(id) == 0 && x == 0 && y == 0 {
        return ifx state then t.on else t.off;
    }
    while x >= (cast(u64)1 << level_of(id)) || y >= (cast(u64)1 << level_of(id)) {
        z  := get_zero(t, level_of(id));
        id  = join(t, id, z, z, z);
    }
    n      := lookup(t, id).*;
    offset := cast(u64)1 << (level_of(id) - 1);
    na, nb, nc, nd := n.a, n.b, n.c, n.d;

    if      x < offset && y < offset   na = set_cell(t, na, x,          y,          state);
    else if x >= offset && y < offset  nb = set_cell(t, nb, x - offset, y,          state);
    else if x < offset && y >= offset  nc = set_cell(t, nc, x,          y - offset, state);
    else                               nd = set_cell(t, nd, x - offset, y - offset, state);

    return join(t, na, nb, nc, nd);
}

get_cell :: (t: *Node_Table, id: node_id, x: u64, y: u64, level: u64) -> float {
    n := lookup(t, id);
    if level_of(id) == 0 || level_of(id) == level {
        return cast(float) n.pop / cast(float)(1 << (2 * level_of(id)));
    }
    size := cast(u64)1 << level_of(id);
    if x >= size || y >= size  return 0.0;
    offset := cast(u64)1 << (level_of(id) - 1);
    if      x < offset  && y < offset   return get_cell(t, n.a, x,          y,          level);
    else if x >= offset && y < offset   return get_cell(t, n.b, x - offset, y,          level);
    else if x < offset  && y >= offset  return get_cell(t, n.c, x,          y - offset, level);
    else                                return get_cell(t, n.d, x - offset, y - offset, level);
}

// ---------------------------------------------------------------------------
// Bounding box of live cells (descends only into non-empty quadrants)
// ---------------------------------------------------------------------------

Bounds :: struct {
    found: bool;
    min_x, min_y, max_x, max_y: s64;
}

find_bounds_rec :: (t: *Node_Table, id: node_id, ox: s64, oy: s64, b: *Bounds) {
    n := lookup(t, id);
    if n.pop == 0  return;

    if level_of(id) == 0 {
        if !b.found {
            b.found = true;
            b.min_x = ox; b.max_x = ox;
            b.min_y = oy; b.max_y = oy;
        } else {
            if ox < b.min_x  b.min_x = ox;
            if ox > b.max_x  b.max_x = ox;
            if oy < b.min_y  b.min_y = oy;
            if oy > b.max_y  b.max_y = oy;
        }
        return;
    }

    half := cast(s64)(cast(u64)1 << (level_of(id) - 1));
    find_bounds_rec(t, n.a, ox,        oy,        b);
    find_bounds_rec(t, n.b, ox + half, oy,        b);
    find_bounds_rec(t, n.c, ox,        oy + half, b);
    find_bounds_rec(t, n.d, ox + half, oy + half, b);
}

find_bounds :: (t: *Node_Table, id: node_id) -> Bounds {
    b: Bounds;
    find_bounds_rec(t, id, 0, 0, *b);
    return b;
}

// ---------------------------------------------------------------------------
// Built-in patterns (each returns a node_id with cells placed near the origin)
// ---------------------------------------------------------------------------

load_glider :: (t: *Node_Table) -> node_id {
    cells : [5][2] s64 : .[
        .[1, 0],
        .[2, 1],
        .[0, 2],
        .[1, 2],
        .[2, 2],
    ];
    id : node_id = t.off;
    for c : cells  id = set_cell(t, id, cast(u64)c[0] + 4, cast(u64)c[1] + 4, true);
    return id;
}

load_glider_gun :: (t: *Node_Table) -> node_id {
    cells : [36][2] s64 : .[
        .[24, 0],
        .[22, 1], .[24, 1],
        .[12, 2], .[13, 2], .[20, 2], .[21, 2], .[34, 2], .[35, 2],
        .[11, 3], .[15, 3], .[20, 3], .[21, 3], .[34, 3], .[35, 3],
        .[0,  4], .[1,  4], .[10, 4], .[16, 4], .[20, 4], .[21, 4],
        .[0,  5], .[1,  5], .[10, 5], .[14, 5], .[16, 5], .[17, 5], .[22, 5], .[24, 5],
        .[10, 6], .[16, 6], .[24, 6],
        .[11, 7], .[15, 7],
        .[12, 8], .[13, 8],
    ];
    id : node_id = t.off;
    for c : cells  id = set_cell(t, id, cast(u64)c[0] + 2, cast(u64)c[1] + 2, true);
    return id;
}

load_r_pentomino :: (t: *Node_Table) -> node_id {
    cells : [5][2] s64 : .[
        .[1, 0], .[2, 0],
        .[0, 1], .[1, 1],
                 .[1, 2],
    ];
    id : node_id = t.off;
    for c : cells  id = set_cell(t, id, cast(u64)c[0] + 30, cast(u64)c[1] + 30, true);
    return id;
}

// ---------------------------------------------------------------------------
// Life 1.05 file parser (verbatim from the console build)
//
// Format summary:
//   #Life 1.05          — required first line
//   #D <text>           — description (ignored)
//   #N                  — use normal (B3/S23) rule (we always use this)
//   #P <x> <y>          — move the pen to world coordinate (x, y)
//   .* rows             — pattern rows: '.' = dead, '*' = alive
//
// A file may contain any number of #P blocks; coordinates can be negative.
// We offset everything so the minimum (x, y) is placed at (16, 16).
// ---------------------------------------------------------------------------

// Parse a signed decimal integer from the start of `s`, advancing `s` past any
// leading whitespace and the digits consumed. Returns (value, ok).
parse_lif_int :: (s: *string) -> int, bool {
    // skip whitespace
    while s.count > 0 && (s.data[0] == #char " " || s.data[0] == #char "\t") {
        s.data  += 1;
        s.count -= 1;
    }
    if s.count == 0  return 0, false;

    negative := false;
    if s.data[0] == #char "-" {
        negative = true;
        s.data  += 1;
        s.count -= 1;
    } else if s.data[0] == #char "+" {
        s.data  += 1;
        s.count -= 1;
    }

    if s.count == 0 || s.data[0] < #char "0" || s.data[0] > #char "9"  return 0, false;

    value := 0;
    while s.count > 0 && s.data[0] >= #char "0" && s.data[0] <= #char "9" {
        value = value * 10 + cast(int)(s.data[0] - #char "0");
        s.data  += 1;
        s.count -= 1;
    }
    return ifx negative then -value else value, true;
}

// Load a Life 1.05 file and return the resulting node_id.
// Returns (id, true) on success, (t.off, false) on error.
load_lif_file :: (t: *Node_Table, filename: string) -> node_id, bool {
    text, ok := read_entire_file(filename);
    if !ok {
        print("Error: cannot read file '%'\n", filename);
        return t.off, false;
    }
    defer free(text.data);

    lines := split(text, "\n");
    defer array_free(lines);

    // --- Validate header ---
    if lines.count == 0 || !begins_with(lines[0], "#Life 1.05") {
        print("Error: '%' is not a Life 1.05 file\n", filename);
        return t.off, false;
    }

    // --- First pass: find the bounding box so we can shift to positive coords ---
    min_x : int = 0x7FFF_FFFF;
    min_y : int = 0x7FFF_FFFF;

    pen_x, pen_y : int = 0, 0;
    row_offset   : int = 0;

    for line : lines {
        // Strip carriage return if present (Windows line endings)
        trimmed := line;
        if trimmed.count > 0 && trimmed.data[trimmed.count - 1] == #char "\r"
            trimmed.count -= 1;

        if begins_with(trimmed, "#P") {
            rest := trimmed;
            rest.data  += 2;
            rest.count -= 2;
            x, xok := parse_lif_int(*rest);
            y, yok := parse_lif_int(*rest);
            if xok && yok {
                pen_x      = x;
                pen_y      = y;
                row_offset = 0;
            }
        } else if trimmed.count > 0 && trimmed[0] != #char "#" {
            for col : 0 .. trimmed.count - 1 {
                if trimmed[col] == #char "*" {
                    wx := pen_x + col;
                    wy := pen_y + row_offset;
                    if wx < min_x  min_x = wx;
                    if wy < min_y  min_y = wy;
                }
            }
            row_offset += 1;
        }
    }

    if min_x == 0x7FFF_FFFF {
        print("Warning: '%' contains no live cells\n", filename);
        return t.off, true;
    }

    // Offset so the topmost-leftmost live cell sits at (16, 16)
    MARGIN :: 16;
    offset_x := MARGIN - min_x;
    offset_y := MARGIN - min_y;

    // --- Second pass: place cells ---
    id           : node_id = t.off;
    pen_x        = 0;
    pen_y        = 0;
    row_offset   = 0;

    for line : lines {
        trimmed := line;
        if trimmed.count > 0 && trimmed.data[trimmed.count - 1] == #char "\r"
            trimmed.count -= 1;

        if begins_with(trimmed, "#P") {
            rest := trimmed;
            rest.data  += 2;
            rest.count -= 2;
            x, xok := parse_lif_int(*rest);
            y, yok := parse_lif_int(*rest);
            if xok && yok {
                pen_x      = x;
                pen_y      = y;
                row_offset = 0;
            }
        } else if trimmed.count > 0 && trimmed[0] != #char "#" {
            for col : 0 .. trimmed.count - 1 {
                if trimmed[col] == #char "*" {
                    wx := cast(u64)(pen_x + col + offset_x);
                    wy := cast(u64)(pen_y + row_offset + offset_y);
                    id = set_cell(t, id, wx, wy, true);
                }
            }
            row_offset += 1;
        }
    }

    return id, true;
}

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
