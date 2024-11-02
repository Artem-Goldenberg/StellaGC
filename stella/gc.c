#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <stddef.h>

#include "runtime.h"
#include "gc.h"

int total_reads = 0;
int total_writes = 0;

#define MAX_GC_ROOTS 1024

int gc_roots_max_size = 0;
int gc_roots_top = 0;
void **gc_roots[MAX_GC_ROOTS];

#define HEAP (MAX_ALLOC_SIZE / 2)

typedef unsigned char byte;

typedef struct {
    char* fromSpace;
    char* toSpace;
    size_t segmentSize;
    /// Next pointer for each segment of the generation
    size_t* next;
    /// Scan pointers as well for convinience
    size_t* scan;

    /// Amount of allocated objects in each segment
    size_t* objects;
    /// Maximum amount of allocated objects ever in this generation
    size_t maxObjects;
    /// Maximum memory ever allocated in this gen
    size_t maxAlloc;
} Generation;

typedef struct {
    int gen;
    size_t segment;
    enum { FromSpace, ToSpace, NoSpace } space;
} Location;

static struct {
    Generation generations[NGENS];

    /// Generation, up to which we are currently collecting, -1 if not collecting
    int collectGen;

    /// `next` pointer of the generation which goes after the last one we're collecting
    /// It is set before each garbage collection to properly know,
    /// whether the object is forwarded to the last generations's `fromSpace`
    size_t barrier;

    /// Cumulative sizes for each generation backwards (see `initialize`)
    size_t offsets[NGENS];

    /// Cells to mark when user writes to memory, for each generation
    byte* marks[NGENS - 1];

    /// Amount of times each generation has been collected
    int timesCollected[NGENS];
    /// Maximum total bytes ever allocated
    size_t maxAlloc;
    /// Maximum total amount of objects
    size_t maxObjects;
} gc;

static Location getLocationFor(const void* p);

static void collect(void);
static void initialize(void);
static bool markCellAt(const void* ptr, int toGen);
static stella_object* forward(stella_object* p);
static size_t getSize(const stella_object* p);

#define FC(ptr) STELLA_OBJECT_HEADER_FIELD_COUNT((ptr)->object_header)

#define MarksToCover(size) ((size) / (CELL_SIZE) / 8u)

#define GetTop(gen, space, seg) ((gen).space + (gen).segmentSize * (seg) + (gen).next[seg])

void* gc_alloc(size_t size_in_bytes) {
    Generation* new = gc.generations; // first generation
    if (!new->fromSpace) initialize();

    assert(size_in_bytes % 8 == 0);

    if (new->next[0] + size_in_bytes > MAX_ALLOC_SIZE)
        collect();

    assert(new->next[0] + size_in_bytes <= MAX_ALLOC_SIZE);

    void* dst = new->fromSpace + new->next[0]; // first segment
    new->next[0] += size_in_bytes;
    new->objects[0]++;

    new->maxAlloc += size_in_bytes;
    new->maxObjects++;

    gc.maxAlloc += size_in_bytes;
    gc.maxObjects++;

    return dst;
}

void print_gc_roots(void) {
    printf("ROOTS: ");
    for (int i = 0; i < gc_roots_top; i++) {
        printf("%p ", gc_roots[i]);
    }
    printf("\n");
}

static const char genTableHeader[] =
"| Generation | Bytes    | MaxBytes | Objects | MaxObjects | Times Collected |";
static const char segTableHeader[] =
"| Segment | Bytes     | Objects |";

static void printSegRowSep(void) {
    for (unsigned i = 1; i < sizeof segTableHeader; ++i) putchar('-');
    printf("\n");
}
static void printGenRowSep(void) {
    for (unsigned i = 1; i < sizeof genTableHeader; ++i) putchar('-');
    printf("\n");
}

void print_gc_alloc_stats(void) {
    size_t totalBytes = 0;
    size_t totalObject = 0;
    int totalCollections = 0;

    // amount of allocated bytes and objects for each generation
    static size_t bytes[NGENS];
    static size_t objects[NGENS];

    printf("\nDetails about each generations:\n");

    for (int gen = 0; gen < NGENS; ++gen) {
        size_t genBytes = 0, genObjects = 0;

        Generation* g = gc.generations + gen;

        printf("\nGeneration %d (segments size = %zu bytes):\n", gen, g->segmentSize);

        printSegRowSep();
        puts(segTableHeader);
        printSegRowSep();

        for (int seg = 0; seg < NSEGMENTS(gen); ++seg) {
            genBytes += g->next[seg];
            genObjects += g->objects[seg];
            printf("| %-7d | %-9zu | %-7zu |\n", seg, g->next[seg], g->objects[seg]);
        }
        printSegRowSep();

        totalBytes += genBytes;
        totalObject += genObjects;
        totalCollections += gc.timesCollected[gen];

        bytes[gen] = genBytes;
        objects[gen] = genObjects;
    }

    printf("\nOverview Table:\n");

    printGenRowSep();
    puts(genTableHeader);
    printGenRowSep();

    for (int gen = 0; gen < NGENS; ++gen) {
        Generation* g = gc.generations + gen;
        printf("| %-10d | %-8zu | %-8zu | %-7zu | %-10zu | %-15d | \n",
               gen, bytes[gen], g->maxAlloc,
               objects[gen], g->maxObjects,
               gc.timesCollected[gen]);
    }

    printGenRowSep();

    printf("| Total      | %-8zu | %-8zu | %-7zu | %-10zu | %-15d | \n",
           totalBytes, gc.maxAlloc,
           totalObject, gc.maxObjects,
           totalCollections);

    printGenRowSep();

    printf("\n");
    printf("Total memory use:        %'d reads and %'d writes\n", total_reads, total_writes);
    printf("Max GC roots stack size: %d roots\n", gc_roots_max_size);
}

void print_gc_state(void) {

    /* 
     <0xFF32AC42>: {<0xAABC4342134>, <0xEF12C3221AC>}
    */

    printf("\nGarbage collector state:\n\n");

    if (gc.collectGen == -1)
        printf("Currently idle\n");
    else
        printf("Currently collecting up to Generation %d, printing info about to-spaces...\n",
               gc.collectGen);

    for (int gen = 0; gen < NGENS; ++gen) {
        Generation* g = gc.generations + gen;
        char* space = (gc.collectGen == -1 ? g->fromSpace : g->toSpace);

        printf("Generation %d:\n", gen);

        size_t genSize = NSEGMENTS(gen) * g->segmentSize;
        printf("\tFromSpace: %zu bytes: <%p>...<%p>\n",
               genSize, g->fromSpace, g->fromSpace + genSize);
        printf("\tToSpace: %zu bytes: <%p>...<%p>\n",
               genSize, g->toSpace, g->toSpace + genSize);

        printf("\tWith %d segments of %zu bytes:\n", NSEGMENTS(gen), g->segmentSize);

        for (int seg = 0; seg < NSEGMENTS(gen); ++seg, space += g->segmentSize) {

            printf("\tSegment %d: <%p>...<%p>:\n", seg, space, space + g->segmentSize);
            printf("\t\tAllocated: %zu bytes\n", g->next[seg]);
            printf("\t\tAvailable: %zu bytes\n", g->segmentSize - g->next[seg]);

            printf("\t\tTotal of %zu objects:\n", g->objects[seg]);

            for (size_t offset = 0; offset < g->next[seg];) {
                const stella_object* p = (stella_object*)(space + offset);
                if (FC(p) == 0) {
                    printf("Incomplete or corrupted object at %p, aborting...\n", p);
                    printf("Try printing later...");
                    return;
                }
                assert(FC(p) > 0);
                printf("\t\t\t<%p>: {<%p>", p, p->object_fields[0]);
                for (int f = 1; f < FC(p); ++f) {
                    printf(", <%p>", p->object_fields[f]);
                }
                printf("}\n");
                offset += getSize(p);
            }
        }
    }

    print_gc_roots();
}

void gc_read_barrier(void *object, int field_index) {
    total_reads += 1;
}

void gc_write_barrier(void *object, int field_index, void *contents) {
    total_writes += 1;
    Location loc = getLocationFor(contents);
    if (loc.space == FromSpace) {
        void* field = ((stella_object*)object)->object_fields + field_index;
        markCellAt(field, loc.gen);
    }
}

void gc_push_root(void **ptr){
    gc_roots[gc_roots_top++] = ptr;
    if (gc_roots_top > gc_roots_max_size) { gc_roots_max_size = gc_roots_top; }
}

void gc_pop_root(void **ptr){
    assert(gc_roots[gc_roots_top - 1] == ptr);
    gc_roots_top--;
}

#define In(gen, space, p) ( \
    InRange(\
        p, \
        space == FromSpace ? gc.generations[gen].fromSpace : gc.generations[gen].toSpace, \
        gc.generations[gen].segmentSize * NSEGMENTS(gen) \
    ) \
)
#define InRange(p, ptr, size) \
    ((char*)(ptr) <= (char*)(p) && (char*)(p) < (char*)(ptr) + (size))
#define Segment(p, ptr, segSize) \
    (((char*)(p) - (char*)(ptr)) / (segSize))

// TODO: debug
int fc(stella_object* p) {
    return FC(p);
}
int nsegs(int gen) {
    return NSEGMENTS(gen);
}
void * getToSpace(Generation* gen, int seg) {
    return GetTop(*gen, toSpace, seg);
}
void * getFromSpace(Generation* gen, int seg) {
    return GetTop(*gen, fromSpace, seg);
}
size_t allocSize(void) {
    return MAX_ALLOC_SIZE;
}

static int init = 0;
static void initialize(void) {
    init += 1;
    assert(init == 1);
    gc.collectGen = -1;
    for (size_t gen = 0, size = MAX_ALLOC_SIZE; gen < NGENS; ++gen, size = GEN_SIZE(size)) {
        gc.generations[gen] = (Generation) {
            .toSpace = malloc(size * NSEGMENTS(gen)),
            .fromSpace = malloc(size * NSEGMENTS(gen)),
            .segmentSize = size,
            .next = calloc(NSEGMENTS(gen), sizeof(size_t)),
            .scan = calloc(NSEGMENTS(gen), sizeof(size_t)),
            .objects = calloc(NSEGMENTS(gen), sizeof(size_t)),
        };
        assert(gc.generations[gen].toSpace);
        assert(gc.generations[gen].fromSpace);
        assert(gc.generations[gen].next);
        assert(gc.generations[gen].scan);
        assert(gc.generations[gen].objects);
    }
    gc.offsets[NGENS - 1] = 0;
    for (int gen = NGENS - 2; gen >= 0; --gen)
        gc.offsets[gen] = 
            gc.offsets[gen + 1] + gc.generations[gen + 1].segmentSize * NSEGMENTS(gen + 1);
    for (int gen = 0; gen < NGENS - 1; ++gen) {
        gc.marks[gen] = calloc(MarksToCover(gc.offsets[gen]), sizeof(byte));
        assert(gc.marks[gen]);
    }
}

static bool markCellAt(const void* ptr, int toGen) {
    Location loc = getLocationFor(ptr);
    assert(loc.space == FromSpace);
    if (loc.gen <= toGen) return false;

    assert(gc.offsets[toGen] - gc.offsets[loc.gen - 1] >= 0LL);

    size_t skipSize = gc.offsets[toGen] - gc.offsets[loc.gen - 1];
    size_t nskipBytes = MarksToCover(skipSize);

    ptrdiff_t offset = (char*)ptr - gc.generations[loc.gen].fromSpace;
    assert(offset >= 0);

    size_t cellIndex = offset / CELL_SIZE;
    size_t markIndex = cellIndex / 8;
    int bit = cellIndex % 8;

    gc.marks[toGen][nskipBytes + markIndex] |= 0x1 << bit;
    return true;
}

static size_t getSize(const stella_object* p) {
    return (1 + FC(p)) * sizeof(void*);
}

static Location getLocationFor(const void* p) {
    for (int gen = 0; gen < NGENS; ++gen) {
        Generation* g = gc.generations + gen;
        if (In(gen, FromSpace, p)) return (Location) {
            .gen = gen,
            .segment = Segment(p, g->fromSpace, g->segmentSize),
            .space = FromSpace
        };
        if (In(gen, ToSpace, p)) return (Location) {
            .gen = gen,
            .segment = Segment(p, g->toSpace, g->segmentSize),
            .space = ToSpace
        };
    }
    return (Location) { .space = NoSpace };
}

static bool atLocation(const void* p, const Location* loc) {
    Generation* gen = gc.generations + loc->gen;
    switch (loc->space) {
        case FromSpace: {
            size_t seg = Segment(p, gen->fromSpace, gen->segmentSize);
            return In(loc->gen, FromSpace, p) && seg == loc->segment;
        }
        case ToSpace: {
            size_t seg = Segment(p, gen->toSpace, gen->segmentSize);
            return In(loc->gen, ToSpace, p) && seg == loc->segment;
        }
        case NoSpace: return false;
    }
}

/// Accepts location of an object `loc` and mutates it to be the location of destination where to copy the object
/// Returns pointer to copy destination, or `NULL` if cannot determine (doesn't change `loc` in this case)
static void* transformToDestination(Location* loc) {
    Generation* gen = gc.generations + loc->gen;
    loc->space = ToSpace;

    if (loc->segment < NSEGMENTS(loc->gen) - 1) {
        // not the last segment
        loc->segment += 1;
        return GetTop(*gen, toSpace, loc->segment);
    }

    // last segment

    if (loc->gen == NGENS - 1) // in last generation
        return GetTop(*gen, toSpace, loc->segment);

    // exists next generation
    loc->gen += 1;
    loc->segment = 0;
    gen++;
    if (loc->gen > gc.collectGen) {
        // we just copy to next gen in it's from space (to youngest segment)
        loc->space = FromSpace; // only here we return FromSpace
        return GetTop(*gen, fromSpace, loc->segment);
    }

    // we are collecting next generation as well, copy to to space
    return GetTop(*gen, toSpace, loc->segment);
}

/// Checks if the object is already forwarded to the destination location `toLoc`
/// Returns pointer to the forwarded object if found, else `NULL`
static stella_object* isForwarded(const stella_object* p, const Location* toLoc) {
    assert(FC(p) > 0);
    void* field0 = p->object_fields[0];
    if (!atLocation(field0, toLoc)) return NULL;

    // check for this case where field0 points to the dest location
    // but it does not point to copied object, just to some old data
    if (gc.collectGen != NGENS - 1 && toLoc->space == FromSpace) {
        assert(toLoc->gen == gc.collectGen + 1);
        Generation* gen = gc.generations + toLoc->gen;
        assert(Segment(field0, gen->fromSpace, gen->segmentSize) == 0);

        size_t offset = ((char*)field0 - gc.generations[toLoc->gen].fromSpace);

        if (offset < gc.barrier) return NULL;
    }

    return field0;
}

static stella_object* forward(stella_object* p) {
    stella_object* current = p;
    stella_object* next = NULL;

    Location loc = getLocationFor(p), destLoc = loc;

    if (loc.gen > gc.collectGen) return p; // not collecting this gen
    if (loc.space != FromSpace) return p; // unknown object

    stella_object* dst = transformToDestination(&destLoc);
    assert(dst);
    assert(FC(current) > 0);

    if ((next = isForwarded(current, &destLoc))) return next; // already forwarded

    while (current) {
        assert(FC(current) > 0);
        next = NULL;
        Generation* destGen = gc.generations + destLoc.gen;

        size_t objectSize = getSize(current);
        destGen->next[destLoc.segment] += objectSize;
        destGen->objects[destLoc.segment]++;

        // out of memory
        if (destGen->next[destLoc.segment] > destGen->segmentSize) {
            printf("Fatal Error: not enough memory to continue\n");
            exit(EXIT_FAILURE);
        }

        dst->object_header = current->object_header;
        for (int f = 0; f < FC(current); ++f) {
            void* field = current->object_fields[f];
            dst->object_fields[f] = field;

            if (next) continue;

            if (field == current) continue; // skip self reference case

            // if field is poining to the same gen and segment as `loc`
            // and is not yet forwarded, then copy this field next
            if (atLocation(field, &loc) && !isForwarded(field, &destLoc))
                next = field;
        }
        current->object_fields[0] = dst;
        dst = (stella_object*)((char*)dst + objectSize);
        current = next;
    }
    return p->object_fields[0];
}

static void forwardInCell(void** cellStart) {
    for (size_t scan = 0; scan < CELL_SIZE / sizeof(void*); ++scan)
        cellStart[scan] = forward(cellStart[scan]);
}

static void forwardInCells(int toGen) {
    // no further points to the oldest ones
    if (toGen == NGENS - 1) return;

    byte* cells = gc.marks[toGen];
    int fromGen = toGen + 1;

    size_t covered = 0;
    for (unsigned i = 0; i < MarksToCover(gc.offsets[toGen]); ++i, covered += CELL_SIZE * 8) {
        if (covered >= gc.generations[fromGen].segmentSize * NSEGMENTS(fromGen)) {
            assert(covered == gc.generations[fromGen].segmentSize * NSEGMENTS(fromGen));
            fromGen++;
            covered = 0;
        }
        while (cells[i] != 0) {
            int cellIndex = __builtin_clz(cells[i]) - 24;
            void* cellStart =
                gc.generations[fromGen].fromSpace + covered + (7 - cellIndex) * CELL_SIZE;
            forwardInCell(cellStart);
            cells[i] &= 0xFF >> (cellIndex + 1);
        }
    }
}

static void forwardInLocation(Location destLoc) {
    Generation* g = gc.generations + destLoc.gen;

    size_t* scan = g->scan + destLoc.segment;

    char* space = destLoc.space == FromSpace ? g->fromSpace : g->toSpace;
    space += destLoc.segment * g->segmentSize;

    while (*scan < g->next[destLoc.segment]) {
        stella_object* p = (stella_object*)(space + *scan);
        assert(FC(p) > 0);
        for (int f = 0; f < FC(p); ++f)
            p->object_fields[f] = forward(p->object_fields[f]);
        *scan += getSize(p);
    }

    assert(*scan == g->next[destLoc.segment]);
}

static int generationToCollect(void) {
    int gen = 0;
    for (; gen < NGENS - 1; ++gen) {
        Generation* currentGen = gc.generations + gen;
        Generation* nextGen = gc.generations + gen + 1;
        size_t spaceLeft = nextGen->segmentSize - nextGen->next[0];
        // last segment can be stored entirely in the first segment of next gen
        if (currentGen->next[NSEGMENTS(gen) - 1] <= spaceLeft) break;
    }
    return gen;
}

#define Swap(type, a, b) do {type tmp = a; a = b; b = tmp;} while (0)

static void updateStats(void) {
    size_t totalBytes = 0, totalObjects = 0;
    for (int gen = 0; gen < NGENS; ++gen) {
        size_t genBytes = 0, genObjects = 0;
        for (int seg = 0; seg < NSEGMENTS(gen); ++seg) {
            genBytes += gc.generations[gen].next[seg];
            genObjects += gc.generations[gen].objects[seg];
        }
        if (genBytes > gc.generations[gen].maxAlloc)
            gc.generations[gen].maxAlloc = genBytes;
        if (genObjects > gc.generations[gen].maxObjects)
            gc.generations[gen].maxObjects = genObjects;
        totalBytes += genBytes;
        totalObjects += genObjects;
    }
    if (totalBytes > gc.maxAlloc) gc.maxAlloc = totalBytes;
    if (totalObjects > gc.maxObjects) gc.maxObjects = totalObjects;
}

static void collect(void) {
    // determine up to what generation we need to collect
    gc.collectGen = generationToCollect();

    // do not need barrier if collecting last gen
    if (gc.collectGen != NGENS - 1) {
        gc.barrier = gc.generations[gc.collectGen + 1].next[0];
    }

    // set `next` and `scan` to 0 in generations we are about to collect
    for (int gen = 0; gen <= gc.collectGen; ++gen) {
        memset(gc.generations[gen].next, 0, NSEGMENTS(gen) * sizeof(size_t));
        memset(gc.generations[gen].scan, 0, NSEGMENTS(gen) * sizeof(size_t));
        // stats
        memset(gc.generations[gen].objects, 0, NSEGMENTS(gen) * sizeof(size_t));
        gc.timesCollected[gen]++;
    }

    // forward all roots
    for (int i = 0; i < gc_roots_top; ++i) {
        *gc_roots[i] = forward(*gc_roots[i]);
    }

    // forward and correct pointers at marked cells
    for (int gen = 0; gen <= gc.collectGen; ++gen)
        forwardInCells(gen);

    // bfs scan
    Location destLoc;
    int cmp = 0;
    do { // until scans in all generations won't be equal to nexts
        for (int gen = 0; gen <= gc.collectGen; ++gen) {
            for (int seg = 0; seg < NSEGMENTS(gen); ++seg) {
                // find out where we copied this segment
                destLoc.gen = gen;
                destLoc.segment = seg;
                transformToDestination(&destLoc);
                forwardInLocation(destLoc);
            }
        }
        cmp = 0;
        for (int gen = 0; gen <= gc.collectGen && cmp == 0; ++gen) {
            Generation* g = gc.generations + gen;
            cmp = memcmp(g->next, g->scan, NSEGMENTS(gen) * sizeof(size_t));
        }
    } while (cmp != 0);

    updateStats();

    // to space and from space swap
    for (int gen = 1; gen <= gc.collectGen; ++gen)
        Swap(char*, gc.generations[gen].fromSpace, gc.generations[gen].toSpace);

    gc.collectGen = -1;
}
