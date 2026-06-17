# Conway's Game of Life in Jai — HashLife

This project is an implementation of John Conway's Game of Life in [Jai](https://en.wikipedia.org/wiki/Jai_(programming_language)).

It began as a small didactic implementation: a fixed-size, power-of-two **toroidal** grid (the board wraps around at the edges) with a straightforward generation step that visits every cell. That version is the one shown in the demo video below, and it is still the easiest place to start if you just want to read a simple Game of Life.

The current code replaces that engine with **HashLife**, Bill Gosper's algorithm for simulating Life on an *unbounded* plane. The grid is no longer a fixed torus — the universe is effectively infinite, and large, regular patterns can be advanced enormous numbers of generations almost for free. The same engine powers two front-ends: an interactive windowed UI and a text-mode console application.

[Game of Life in Jai YouTube Demo Video](https://www.youtube.com/watch?v=wJJ6XU2Fv1w) (original simple toroidal version)

---

## On John Conway's Game of Life

John Conway's Game of Life is a cellular automaton created by the British mathematician John Conway. The world is a grid of cells, each either alive or dead. Every generation, each cell's next state is decided from its eight neighbours:

* Any live cell with fewer than two live neighbours dies, as if by underpopulation.
* Any live cell with two or three live neighbours lives on to the next generation.
* Any live cell with more than three live neighbours dies, as if by overpopulation.
* Any dead cell with exactly three live neighbours becomes a live cell, as if by reproduction.

This is the standard **B3/S23** rule. From these four simple rules, cells produce surprisingly intricate, chaotic, and aesthetically pleasing patterns.

[John Conway's Game of Life (Wikipedia)](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life)

---

## The HashLife algorithm

HashLife (Gosper, 1984, *"Exploiting Regularities in Large Cellular Spaces"*) is a way to simulate the Game of Life that exploits the heavy repetition found in most patterns — in both space and time. Instead of storing a flat array of cells and stepping every cell each generation, it represents the world as a **quadtree** and reuses work aggressively.

### How it works

* **Quadtree of macro-cells.** The infinite grid is stored as a tree of square nodes. A node at *level n* covers a `2^n × 2^n` region. Level 0 is a single cell. Every non-leaf node has four children — `a` (top-left), `b` (top-right), `c` (bottom-left), and `d` (bottom-right):

  ```
  a | b
  -----
  c | d
  ```

* **Hash-consing (canonical sharing).** Every distinct sub-pattern is stored exactly once in a hash table. Two regions of the board that look identical share the *same* node, no matter where they appear or how often. A large empty area, for example, is a single shared "all-dead" node rather than millions of dead cells. This is what keeps memory small even for huge boards.

* **Memoised successors.** For a node of level *n* (with *n ≥ 2*), HashLife can compute the future of its central region and **cache** that result, keyed on the node. Because identical nodes are shared, a result computed once is instantly available everywhere that pattern recurs. The only place the actual Life rule is evaluated is the base case: a level-2 (4×4) node, whose central 2×2 is produced by applying B3/S23 directly.

* **Stepping in powers of two.** The core operation, `successor(node, j)`, advances a node by `2^j` generations (capped at the largest jump the node is big enough to represent). Advancing by an arbitrary number of generations is done by binary decomposition — each set bit of the step count contributes one `successor` call. This is why HashLife can leap millions of generations in a single, cheap operation when a pattern is regular enough.

* **Node identifiers.** Internally each node is referred to by a 64-bit id that packs a flag for all-dead nodes, the node's level, and a hash. Ids are stored in an open-addressed hash table with linear probing, and a separate one-shot cache memoises successor results.

### Why it's fast

Because both space (shared subtrees) and time (cached successors) are deduplicated, regular and repetitive patterns are simulated in time that grows roughly with the *logarithm* of the number of generations, rather than linearly with the number of cells per step. A glider gun or a highway can run for billions of generations in moments. (Highly chaotic patterns share less and so benefit less — HashLife is fastest exactly when the pattern is regular.)

### How it differs from the simple version

* The original engine used a **fixed, power-of-two toroidal grid** that wrapped around at the edges. HashLife uses an **unbounded plane** — there are no edges and no wrap-around. Patterns may grow without limit, and the view follows them.
* The simple engine stepped one generation at a time over every cell. HashLife steps a whole quadtree at once and can advance by many generations per step.

---

## The console application

The console application is a text-mode front-end to the same HashLife engine. It loads a pattern, then repeatedly renders an ASCII view of the board to the terminal and advances the simulation, sleeping briefly between frames. The view is automatically re-centred on the living cells each frame (since the pattern drifts as it evolves), so a translating pattern such as a glider stays visible. If the pattern dies out, the run stops early.

It is driven entirely by command-line arguments:

```
hashlife [options] [file.lif]
```

| Option            | Description                                                          |
| ----------------- | -------------------------------------------------------------------- |
| `file.lif`        | Load a Life 1.05 pattern file (see below).                           |
| `--glider`        | Start with the built-in glider. *(Default if no file is given.)*     |
| `--gun`           | Start with the Gosper glider gun.                                    |
| `--rpentomino`    | Start with the R-pentomino.                                          |
| `--frames <n>`    | Total number of frames to simulate. *(Default: 120.)*                |
| `--sleep <ms>`    | Milliseconds to pause between frames. *(Default: 80.)*               |
| `--help`          | Print usage and exit.                                                |

The windowed UI accepts the same arguments. In the UI, the two timing-related options are adapted to the interactive setting: `--sleep <ms>` sets the initial simulation speed (the **Update** control), and `--frames <n>` auto-runs the simulation for *n* generations on startup and then pauses, leaving the window fully interactive afterwards. Passing a `.lif` file or one of the built-in flags chooses the pattern shown when the window opens.

---

## Loading patterns: Life 1.05 (`.lif`) files

Both front-ends can load patterns from **Life 1.05** files, a simple plain-text format for describing Game of Life boards. Give the path to a `.lif` (or `.LIF`) file on the command line.

The format is line-oriented:

| Line              | Meaning                                                                            |
| ----------------- | ---------------------------------------------------------------------------------- |
| `#Life 1.05`      | Required first line; identifies the file.                                          |
| `#D <text>`       | A human-readable description. Ignored by the loader.                               |
| `#N`              | Selects the normal **B3/S23** rule (the only rule supported here).                 |
| `#P <x> <y>`      | Moves the "pen" to board coordinate `(x, y)`. Subsequent rows are placed from here. |
| pattern rows      | A run of cells: `.` is a dead cell, `*` is a live cell.                             |

A file may contain any number of `#P` blocks, and coordinates may be negative. When a file is loaded, the parser scans every live cell to find the pattern's bounding box, then shifts the whole pattern so that its top-left-most live cell sits at `(16, 16)`. This keeps all coordinates positive and away from the origin, which suits the quadtree layout.

A glider, for example, can be written as:

```
#Life 1.05
#D Glider
#N
#P 0 0
.*.
..*
***
```

---

## Controls (windowed UI)

* **Click a grid square** to toggle that cell on or off. *(Editing is enabled while the simulation is paused.)*
* **Simulate / Stop** — start or stop the running simulation.
* **Next** — advance the board by the current step size (see *Step* below).
* **Zoom In / Zoom Out** — change how large each cell is drawn.
* **Clear** — empty the board.
* **Update** slider — the time between automatic steps; smaller is faster.
* **Step** slider — the number of generations advanced per step. Increasing this is where HashLife shows its strength: large jumps cost little more than small ones for regular patterns.
* **Colour** dropdown — choose the colour of live cells.
* **Pattern** dropdown — load a built-in pattern: *Glider*, *Glider Gun*, *R-Pentomino*, or *Empty*.
* **Arrow keys** — pan the camera up, down, left, or right.

Because the universe is unbounded, the view automatically **follows the pattern** while the simulation is running, so it never simply scrolls off the edge. Panning with the arrow keys is therefore most useful while paused. Note that for an open-ended emitter such as the Gosper gun, the stream of escaping gliders pulls the auto-centred view along with it — pause and pan, or watch a bounded pattern, to keep a fixed frame.

---

## Building

The project is written in Jai and built with the Jai compiler, for example:

```
jai hashlife_ui.jai        # the windowed UI
jai hashlife_console.jai   # the text-mode console application
```

(Adjust the file names to match your checkout.) The UI depends on Jai's `Simp`, `Window_Creation`, `GetRect`, `Input`, and `Math` modules; both front-ends use `Basic`, and the `.lif` loader additionally uses the `String` and `File` modules.
