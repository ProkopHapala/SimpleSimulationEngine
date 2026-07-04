# dataStructures

Header-only data structure library. No external dependencies beyond standard C++ and `macroUtils.h`.

- **datatypes.h** — Bare POD structs (float4, int4, double8, ...) for GPU buffers and C-style interchange. Guaranteed layout, no methods — maps directly to OpenCL types and raw binary.
- **datatypes2cl.h** — Bridge between double-precision CPU types (Vec3d, Quat4d) and float-precision GPU types (Quat4f). Provides pack/unpack/copy with index remapping (add_except) for mesh compaction on GPU upload.
- **Slots.h** — Fixed-capacity stack array for bounded-size associations. Born from edge-loop sorting where each vertex has at most 2 incident loop edges — Slots<int,2> avoids heap allocation entirely.
- **Buckets.h** — CSR (counting-sort) layout for static bucketing. Two-pass: count then scatter. Used as edgesOfVerts in Builder2 for O(1) "which edges touch vertex i?" queries. Predecessor to NeighChunks.
- **NeighChunks.h** — Cache-line-aligned neighbor lists with overflow chaining. Each vertex gets one 64-byte chunk (14 neighbors inline); overflow allocates extension chunks at the end of the flat array. Designed for dynamic mesh editing where vertex degree changes.
- **HashMap.h** — Open-addressing hash map for spatial bucketing with integer box-index keys. Knuth multiplicative hashing, power-of-2 capacity, hits[] array for early termination in getAllInBox(). Auto-resizes at 50% load.
- **HashMap64.h** — 64-bit key variant of HashMap for large-scale spatial hashing. Box IDs are 64-bit (21 bits/axis = 2M km per axis at 1km resolution). Notes discuss key width vs. spatial range tradeoffs.
- **HashMapT64.h** — Templated 64-bit key hash map. Same open-addressing design as HashMap64 but stores arbitrary type T instead of raw pointers. Caller sets EMPTY_O sentinel.
- **HashMap_temp.h** — Experimental hash map variants with bucket-based storage (HashMapField). Tests multiple hash functions; Knuth's method confirmed best. Historical/experimental — prefer HashMap.h or HashMapT64.h.
- **BatchBuff.h** — Chunked sparse array for O(1) random access to huge index ranges without allocating all slots. Batches allocated lazily on first write. BatchBuffPow2 uses power-of-2 batch size for bit-shift/bitmask instead of division.
- **Table.h** — Column-oriented table over a flat byte buffer. Records field offsets at bind time, provides runtime type-dispatched toStr/fromStr. Bridges cache-efficient POD struct arrays with named-field I/O (CSV, debugging).
- **Tree.h** — CRTP-based recursive tree template. TreeGen<T,TreeT> solves the incomplete-type problem; Tree<T> is value-based, PTree<T> is pointer-based with parent back-links. For scene graphs, BVH, nested mesh groups.
- **ListKeeper.h** — Simple growable int array with golden-ratio growth factor. Tracks free slots, supports shrink when emptiness exceeds threshold. Minimal slot allocator.
- **parsing.h** — Shared parser primitives: ParserItem struct (token span, branch count), charTypes[] LUT for ASCII classification. Foundation for all parsers in this folder.
- **LispParser.h** — S-expression parser inspired by jsmn. Tokenizes parenthesized lists into ParserItem array. Configurable open/close/separator characters.
- **StructParser.h** — Brace-delimited struct parser. Parses `{ field; field; }` syntax into ParserItem tree. Used for reading structured config files.
- **TreeParser.h** — Expression tree parser with 4 coupling levels: name-binding > void-binding > operator > tuple separator. Handles operator precedence and nested expressions. Most complex parser here.
- **AST.h** — Abstract syntax tree built from parser tokens. Flat array of Nodes with parent/child indices (no pointers). Nodes reference ParserItem tokens.
- **AST2Lang.h** — Symbol table and type system on top of AST. Defines Symbol, Type, Variable, Function structs for semantic analysis. Bridges parse tree to typed program representation.
- **ProgramGraph.h** — Scoped program graph with hashed names. Supports variable/function/type definitions with scope-aware lookup. Designed for a scripting language with C++ interop (typed callbacks + typeless macros).
- **CommandParser.h** — Simple command-line parser. Commands separated by `;`, first word is command name, arguments are primitive constants. Application provides its own call table via function pointers.
- **CommandParser2.h** — Extended scripting language parser with variable assignment (`$out = Func(in)`), multiple return values, typed callbacks and typeless macros. Three notation styles supported.
