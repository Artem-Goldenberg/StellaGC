src = stella

compiler = docker run -i fizruk/stella compile
sources = $(src)/runtime.c $(src)/gc.c

targets = factorial squarenums logicops power2

CFLAGS += -std=c11
CFLAGS += -Wall -Wextra
CFLAGS += -Wno-unused-parameter -Wno-zero-length-array
CFLAGS += -Wno-self-assign

ifeq ($(DEBUG),1)
CFLAGS += -DSTELLA_GC_STATS
CFLAGS += -DSTELLA_RUNTIME_STATS
else ifeq ($(DEBUG),2)
CFLAGS += -DSTELLA_DEBUG
CFLAGS += -DSTELLA_GC_STATS
CFLAGS += -DSTELLA_RUNTIME_STATS
endif

help:
	@echo "Available targets: $(targets) clean"

$(targets): %: %.c $(sources) $(sources:%.c=%.h)
	$(LINK.c) $< $(sources) -o $@

%.c: Examples/%.stella
	$(compiler) < $< > $@

.PHONY: clean
clean:
	$(RM) Examples/*.c $(targets) $(targets:%=%.c)
