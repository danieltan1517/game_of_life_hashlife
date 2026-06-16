/*
    HashLife — Conway's Game of Life using Gosper's HashLife algorithm
    Jai implementation based on Gosper (1984) "Exploiting Regularities in Large Cellular Spaces"
    and the C reference implementation by John H. Williamson.

    Key idea:
      - The infinite grid is represented as a quadtree of "macro-cells".
      - A macro-cell at level n covers a 2^n × 2^n region.
      - Identical sub-trees are shared via hash-consing (one entry per unique pattern).
      - The RESULT of a macro-cell (its center after 2^(n-2) steps) is memoised,
        so repeated patterns are advanced for free.
      - This gives super-linear (often near-O(log t)) simulation of repetitive Life patterns.

    Layout of a node_id (u64 bit-packed identifier):
        [bit 63]   : MARKED flag (used during GC / vacuum)
        [bit 62]   : IS_ZERO flag (node represents all-dead cells)
        [bits 61-46]: level (16 bits, 0 = leaf cell)
        [bits 45-0] : hash / disambiguation bits

    Node quadrant convention (matches Gosper's paper):
        a | b
        -----
        c | d
    where a is top-left, b top-right, c bottom-left, d bottom-right.
*/

#import "Basic";
#import "String";
#import "File";

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

// Separate "next" cache: maps (from, j) -> to
// One-shot: collisions evict the existing entry (no chaining needed;
// a miss just means recomputation, which is always correct).
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

    // Build the two leaf ids using the same bit-packing as the C reference.
    // OFF leaf: IS_ZERO flag set, level=0, hash derived from 0
    t.off = (cast(u64)1 << 62) | (cast(u64)0 << 46) | (mix64(0) & HASH_MASK);
    // ON  leaf: IS_ZERO clear,  level=0, hash derived from 1
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
    
    // Check if the canonical zero-node for level k already exists
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
    z := get_zero(t, level_of(m_h) - 1);  // safe: level >= 1 guaranteed above
    na := join(t, z,   z,   z,   m.a);
    nb := join(t, z,   z,   m.b, z  );
    nc := join(t, z,   m.c, z,   z  );
    nd := join(t, m.d, z,   z,   z  );
    return join(t, na, nb, nc, nd);
}

// ---------------------------------------------------------------------------
// Base-case: apply Life rule to the 4×4 block, returning the central 2×2
// (This is the only place where the actual cellular automaton rule lives.)
//
// Quadrant layout of a 4×4 node (level=2):
//   Each quadrant is a 2×2 level-1 node, itself with four 1×1 leaves.
//   Leaf positions:
//     a.a a.b | b.a b.b
//     a.c a.d | b.c b.d
//     --------+--------
//     c.a c.b | d.a d.b
//     c.c c.d | d.c d.d
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
    // Conway's Life rule: survive on 2 or 3; born on 3
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
    slevel := cast(s64) level;   // signed so subtraction never underflows
    jj     := j;
    // j is the log2 of the step size: j=0 means 1 generation (2^0),
    // j=1 means 2 generations (2^1), etc. Clamp ONLY when j exceeds the
    // node's capacity (level-2). Do NOT special-case j==0 — that would make
    // a single-generation step impossible, which freezes period-2 oscillators
    // (the step would always be even and land on the same phase).
    if cast(s64)jj > slevel - 2  jj = cast(u64)(slevel - 2);

    // All-dead nodes: result is just a quadrant (also all-dead at level-1)
    if is_zero(id)  return lookup(t, id).a;

    next := lookup_next(t, id, jj);
    if next != UNUSED  return next;

    if level == 2 {
        next = life_4x4(t, id);
        cache_next(t, id, next, jj);
        return next;
    }

    // Copy node fields before any table mutations (rehash could move things)
    n := lookup(t, id).*;
    a := lookup(t, n.a).*;
    b := lookup(t, n.b).*;
    c := lookup(t, n.c).*;
    d := lookup(t, n.d).*;

    // Nine intermediate results covering the 3×3 overlapping sub-grid
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
        // Sub-natural: combine the inner corners of the nine results
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
        // Natural successor: recurse one more level on the nine results
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
    // Peel outer zero rings while we remain above level 3
    while level_of(id) > 3 && is_padded(t, id)  id = inner(t, id);
    return id;
}

pad :: (t: *Node_Table, id: node_id) -> node_id {
    // Expand until the outer ring is non-zero (or level >= 3)
    while level_of(id) < 3 || is_padded(t, id)  id = centre(t, id);
    return id;
}

// ---------------------------------------------------------------------------
// advance: simulate forward by exactly `steps` generations.
//
// HashLife's `successor(id, j)` computes the state after 2^j steps (clamped
// to 2^(level-2)). To advance an arbitrary number of generations we use binary
// decomposition: each set bit i in `steps` contributes 2^i generations, which
// we handle with successor(id, j=i). The loop therefore starts j at 0, not 1:
//   steps=1  → j=0 → successor(id, 0) → 1 generation  (base case: life_4x4)
//   steps=2  → j=1 → successor(id, 1) → 2 generations
//   steps=3  → j=0 then j=1 → 1 + 2 = 3 generations, etc.
//
// Padding: successor's result is the inner half of the node, so a translating
// pattern walks off the boundary without vacuum around it. Before every step
// we grow to the required level and add two centre rings for breathing room.
// crop() at the end keeps the tree compact between calls.
// ---------------------------------------------------------------------------

advance :: (t: *Node_Table, id: node_id, steps: u64) -> node_id {
    j   : u64 = 0;
    rem := steps;
    while rem > 0 {
        if rem & 1 {
            // Grow until the node can represent a 2^j-generation step,
            // then add padding rings so the result is never clipped.
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
    // A level-0 node is a single 1×1 cell at coordinate (0,0).
    // We can only set it directly when (x,y) == (0,0).
    // For any other coordinate we must grow into a larger tree first.
    if level_of(id) == 0 && x == 0 && y == 0 {
        return ifx state then t.on else t.off;
    }
    // Grow the tree until (x, y) fits inside its 2^level × 2^level footprint.
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
// Console rendering
// ---------------------------------------------------------------------------

// Find the bounding box of live cells by recursively descending only into
// non-empty quadrants. Coordinates are in the node's own frame.
// Returns whether any live cell was found, plus min/max x and y.
Bounds :: struct {
    found: bool;
    min_x, min_y, max_x, max_y: s64;
}

find_bounds_rec :: (t: *Node_Table, id: node_id, ox: s64, oy: s64, b: *Bounds) {
    n := lookup(t, id);
    if n.pop == 0  return;   // empty subtree — skip entirely (this is the speed win)

    if level_of(id) == 0 {
        // single live cell at (ox, oy)
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
// Patterns: each returns a node_id with the pattern placed near the origin.
// The renderer auto-centres on the live cells, so no coordinates are needed.
// ---------------------------------------------------------------------------

// Glider (period 4, translates diagonally)
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

// Gosper Glider Gun (period 30, emits gliders)
load_glider_gun :: (t: *Node_Table) -> node_id {
    // RLE-decoded live cells of the Gosper Glider Gun
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

// R-pentomino (chaotic 1100-generation methuselah)
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
// Life 1.05 file parser
//
// Format summary:
//   #Life 1.05          — required first line
//   #D <text>           — description (ignored)
//   #N                  — use normal (B3/S23) rule (we always use this)
//   #P <x> <y>          — move the pen to world coordinate (x, y)
//   .* rows             — pattern rows: '.' = dead, '*' = alive
//                         blank lines count as a dead row within a block
//
// A file may contain any number of #P blocks; coordinates can be negative.
// We offset everything so that the minimum (x, y) is placed at (16, 16),
// keeping the pattern away from (0, 0) and well inside the quadtree.
// ---------------------------------------------------------------------------

// Parse a signed decimal integer from the start of `s`, advancing `s` past
// any leading whitespace and the digits consumed. Returns (value, ok).
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

    // --- First pass: collect all live cell coordinates to find the bounding box ---
    // We need the min (x, y) so we can shift the pattern to positive coordinates.
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
            // Skip past "#P" and parse the two integers
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
            // Pattern row
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

    // If no live cells were found, return the empty node
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
