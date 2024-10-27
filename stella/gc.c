#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "runtime.h"
#include "gc.h"

/** Total allocated number of bytes (over the entire duration of the program). */
int total_allocated_bytes = 0;

/** Total allocated number of objects (over the entire duration of the program). */
int total_allocated_objects = 0;

int max_allocated_bytes = 0;
int max_allocated_objects = 0;

int total_reads = 0;
int total_writes = 0;

#define MAX_GC_ROOTS 1024

int gc_roots_max_size = 0;
int gc_roots_top = 0;
void **gc_roots[MAX_GC_ROOTS];

#define HEAP (MAX_ALLOC_SIZE / 2)

struct {
    char* toSpace;
    char* fromSpace;
    size_t scan;
    size_t next;
    size_t limit;
    bool active;
} gc;

static void collect(void);
static void initialize(void);
static stella_object* forward(stella_object* p);

size_t allocSize(void) {
    return MAX_ALLOC_SIZE;
}

void* gc_alloc(size_t size_in_bytes) {
    if (!gc.toSpace) initialize();

    assert(size_in_bytes % 8 == 0);

    if (gc.next + size_in_bytes >= HEAP)
        collect();

    void* dst = gc.fromSpace + gc.next;
    assert(gc.next + size_in_bytes < HEAP);
    gc.next += size_in_bytes;

    total_allocated_bytes += size_in_bytes;
    total_allocated_objects += 1;
    max_allocated_bytes = total_allocated_bytes;
    max_allocated_objects = total_allocated_objects;
//    return malloc(size_in_bytes);
    return dst;
}

void print_gc_roots(void) {
    printf("ROOTS: ");
    for (int i = 0; i < gc_roots_top; i++) {
        printf("%p ", gc_roots[i]);
    }
    printf("\n");
}

void print_gc_alloc_stats(void) {
    printf("Total memory allocation: %'d bytes (%'d objects)\n", total_allocated_bytes, total_allocated_objects);
    printf("Maximum residency:       %'d bytes (%'d objects)\n", max_allocated_bytes, max_allocated_objects);
    printf("Total memory use:        %'d reads and %'d writes\n", total_reads, total_writes);
    printf("Max GC roots stack size: %'d roots\n", gc_roots_max_size);
}

void print_gc_state(void) {
    // TODO: not implemented
}

#define FC(ptr) STELLA_OBJECT_HEADER_FIELD_COUNT((ptr)->object_header)
#define IN(space, ptr) (gc.space <= (char*)(ptr) && (char*)(ptr) < gc.space + HEAP)

void gc_read_barrier(void *object, int field_index) {
    if (IN(toSpace, object)) {
        stella_object* p = object;
        p->object_fields[field_index] = forward(p->object_fields[field_index]);
    }
    total_reads += 1;
}

void gc_write_barrier(void *object, int field_index, void *contents) {
    total_writes += 1;
}

void gc_push_root(void **ptr){
    gc_roots[gc_roots_top++] = ptr;
    if (gc_roots_top > gc_roots_max_size) { gc_roots_max_size = gc_roots_top; }
}

void gc_pop_root(void **ptr){
    assert(gc_roots[gc_roots_top - 1] == ptr);
    gc_roots_top--;
}

int debugFC(stella_object* p) {
    return FC(p);
}
bool debugDone(stella_object* p) {
    return IN(toSpace, p);
}
bool fromSpace(stella_object* p) {
    return gc.fromSpace <= (char*)(p) && (char*)(p) < gc.fromSpace + HEAP;
}

static void initialize(void) {
    gc.toSpace = malloc(HEAP);
    gc.fromSpace = malloc(HEAP);
}

static int getSize(stella_object* p) {
    return (1 + FC(p)) * sizeof(void*);
}

static void copyObject(void* dst, stella_object* src) {
    memcpy(dst, src, getSize(src));
}

static stella_object* forward(stella_object* p) {
    if (IN(toSpace, p)) return p; // already in to-space
    if (!IN(fromSpace, p)) return p; // unknown object

    assert(FC(p) > 0);

    stella_object* field0 = p->object_fields[0];
    if (IN(toSpace, field0)) return field0;

    int size = getSize(p);
    void* dst = gc.toSpace + gc.next;
    memcpy(dst, p, size);
    gc.next += size;
    assert(gc.next < HEAP);

    p->object_fields[0] = dst;

    return dst;
}

static void collect(void) {
    gc.active = true;
    gc.scan = gc.next = 0;
    for (int i = 0; i < gc_roots_top; ++i) {
        *gc_roots[i] = forward(*gc_roots[i]);
    }
    while (gc.scan < gc.next) {
        stella_object* p = (stella_object*)(gc.toSpace + gc.scan);
        for (int f = 0; f < FC(p); ++f)
            p->object_fields[f] = forward(p->object_fields[f]);
        gc.scan += getSize(p);
    }
    assert(gc.scan == gc.next);
    char* tmp = gc.toSpace;
    gc.toSpace = gc.fromSpace;
    gc.fromSpace = tmp;
    gc.active = false;
}
