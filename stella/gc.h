#ifndef STELLA_GC_H
#define STELLA_GC_H

#include <stdlib.h>
#include <stdio.h>

/// Heap size of the youngest generation of objects in bytes (only from-space)
/// If objects stay in 0th gen for more than one iteration,
/// then it is the size of the **first segment**
#define MAX_ALLOC_SIZE (1LU << 11)

/// Number of generations, including the youngest one
#define NGENS 3

/// Segment size (bytes) for each generation, given segment size of the generation before this one
#define GEN_SIZE(genSizeBefore) ((genSizeBefore) * 2U)
/// Number of iterations the object can live in the generation `0 <= gen < NGENS`
/// i.e. number of segments in each generation (segments have the same size)
#define NSEGMENTS(gen) ((gen) == 0 ? 1 : 3)

/// Granularity of cells, by which mutated data will be tracked (byte size of one cell)
#define CELL_SIZE 16LU

/** This macro is used whenever the runtime wants to READ a heap object's field.
 */
#define GC_READ_BARRIER(object, field_index, read_code) (void *)(gc_read_barrier(object, field_index), read_code) // NO BARRIER
/** This macro is used whenever the runtime wants to OVERWRITE a heap object's field.
 * This is NOT used when initializing object fields.
 */
#define GC_WRITE_BARRIER(object, field_index, contents, write_code) (gc_write_barrier(object, field_index, contents), write_code) // NO BARRIER

/** Allocate an object on the heap of AT LEAST size_in_bytes bytes.
 * If necessary, this should start/continue garbage collection.
 * Returns a pointer to the newly allocated object.
 */
void* gc_alloc(size_t size_in_bytes);

/** GC-specific code which must be executed on each READ operation.
 */
void gc_read_barrier(void *object, int field_index);
/** GC-specific code which must be executed on each WRITE operation
 * (except object field initialization).
 */
void gc_write_barrier(void *object, int field_index, void *contents);

/** Push a reference to a root (variable) on the GC's stack of roots.
 */
void gc_push_root(void **object);
/** Pop a reference to a root (variable) on the GC's stack of roots.
 * The argument must be at the top of the stack.
 */
void gc_pop_root(void **object);

/** Print GC statistics. Output must include at least:
 *
 * 1. Total allocated memory (bytes and objects).
 * 2. Maximum residency (bytes and object).
 * 3. Memory usage (number of reads and writes).
 * 4. Number of read/write barrier triggers.
 * 5. Total number of GC cycles (for each generation, if applicable).
 */
void print_gc_alloc_stats(void);

/** Print GC state. Output must include at least:
 *
 * 1. Heap state.
 * 2. Set of roots (or heap pointers immediately available from roots).
 * 3. Current allocated memory (bytes and objects).
 * 4. GC variable values (e.g. scan/next/limit variables in a copying collector).
 */
void print_gc_state(void);

/** Print current GC roots (addresses).
 * May be useful for debugging.
 */
void print_gc_roots(void);

#endif
