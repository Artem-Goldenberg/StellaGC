src = stella

compiler = docker run -i fizruk/stella compile
sources = $(src)/runtime.c $(src)/gc.c

targets = factorial

CFLAGS += -std=c11
CFLAGS += -Wall -Wextra
CFLAGS += -Wno-unused-parameter -Wno-zero-length-array
CFLAGS += -Wno-self-assign

ifneq ($(DEBUG),)
CFLAGS += -DSTELLA_DEBUG
CFLAGS += -DSTELLA_GC_STATS
CFLAGS += -DSTELLA_RUNTIME_STATS
endif

$(targets): %: %.c $(sources)
	$(LINK.c) $^ -o $@

%.c: Examples/%.stella
	$(compiler) < $< > $@

.PHONY: clean
clean:
	$(RM) Examples/*.c $(targets) $(targets:%=%.c)
