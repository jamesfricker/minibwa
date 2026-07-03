CC ?= gcc
AR ?= ar

BUILD_DIR ?= build
OBJ_DIR := $(BUILD_DIR)/obj

PROG := minibwa
LIB_TARGET := libminibwa.a

CFLAGS ?= -std=c99 -g -Wall -O3
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lpthread -lz -lm
AUTO_CFLAGS :=
AUTO_CPPFLAGS :=
AUTO_LDFLAGS :=
AUTO_LDLIBS :=
OPT_CFLAGS :=
OPT_LDFLAGS :=
MIMALLOC_CFLAGS ?= -std=gnu11 -O3 -Wall -Wextra -DNDEBUG -DMI_MALLOC_OVERRIDE -DMI_OSX_INTERPOSE=1 -DMI_OSX_ZONE=1

PGO_PROFILE_DIR ?= $(BUILD_DIR)/pgo-data
PGO_PROFILE_DATA ?= $(PGO_PROFILE_DIR)/default.profdata
PGO_TRAIN_CMD ?=
LLVM_PROFDATA ?= $(shell command -v llvm-profdata 2>/dev/null || { command -v xcrun >/dev/null 2>&1 && xcrun -f llvm-profdata 2>/dev/null; })

CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -n 1)
PGO_GENERATE_FLAGS ?= -fprofile-generate=$(PGO_PROFILE_DIR)
ifneq ($(findstring clang,$(CC_VERSION)),)
PGO_USE_FLAGS ?= -fprofile-use=$(PGO_PROFILE_DATA)
PGO_MERGE_PROFILE := 1
else
PGO_USE_FLAGS ?= -fprofile-use=$(PGO_PROFILE_DIR) -fprofile-correction
PGO_MERGE_PROFILE := 0
endif

INCLUDES := -Iinclude -Isrc -Ithird_party/libsais
ARCH := $(shell uname -m)
omp ?= $(shell printf '\043include <omp.h>\nint main(){return 0;}' | $(CC) -x c -fopenmp -o /dev/null - 2>/dev/null && echo "1" || echo "0")

LIB_SRCS := \
	src/kommon.c \
	src/kalloc.c \
	src/bwt.c \
	src/l2bit.c \
	src/options.c \
	src/unmap.c \
	src/seed.c \
	src/par.c \
	src/map-algo.c \
	src/lchain.c \
	src/align.c \
	src/pe.c \
	src/sv_blacklist.c \
	src/mappability.c \
	src/cs.c \
	src/format.c \
	src/human_resources.c \
	src/ksw2_extz2_sse.c \
	src/ksw2_extd2_sse.c \
	src/ksw2_ll_sse.c

APP_SRCS := \
	src/kthread.c \
	third_party/libsais/libsais.c \
	third_party/libsais/libsais64.c \
	src/index.c \
	src/bseq.c \
	src/map-main.c \
	src/fastmap.c

GPL_SRCS := \
	third_party/bwtgen/QSufSort.c \
	third_party/bwtgen/bwtgen.c

MAIN_SRC := src/main.c
MIMALLOC_SRC := third_party/mimalloc/static.c

LIB_OBJS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(LIB_SRCS))
APP_OBJS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(APP_SRCS))
MAIN_OBJ := $(patsubst %.c,$(OBJ_DIR)/%.o,$(MAIN_SRC))
GPL_OBJS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(GPL_SRCS))
MIMALLOC_OBJ := $(patsubst %.c,$(OBJ_DIR)/%.o,$(MIMALLOC_SRC))

ifneq ($(asan),)
	override AUTO_CFLAGS += -fsanitize=address
	override AUTO_LDFLAGS += -fsanitize=address
	override AUTO_LDLIBS += -ldl
endif

ifeq ($(omp),1)
	override AUTO_CPPFLAGS += -DLIBSAIS_OPENMP
	override AUTO_CFLAGS += -fopenmp
	override AUTO_LDLIBS += -fopenmp
endif

ifneq ($(gpl),0)
	APP_OBJS += $(GPL_OBJS)
	override AUTO_CPPFLAGS += -DUSE_GPL
endif

ifeq ($(mimalloc),0)
	MIMALLOC_OBJ :=
	override AUTO_CPPFLAGS += -DHAVE_KALLOC
endif

ifeq ($(ARCH),x86_64)
	override AUTO_CFLAGS += -msse4.2 -mpopcnt
endif

ifeq ($(lto),1)
	override OPT_CFLAGS += -flto
	override OPT_LDFLAGS += -flto
endif

ifeq ($(pgo),generate)
	override OPT_CFLAGS += $(PGO_GENERATE_FLAGS)
	override OPT_LDFLAGS += $(PGO_GENERATE_FLAGS)
else ifeq ($(pgo),use)
	override OPT_CFLAGS += $(PGO_USE_FLAGS)
	override OPT_LDFLAGS += $(PGO_USE_FLAGS)
else ifneq ($(pgo),)
$(error pgo must be empty, generate, or use)
endif

ALL_CPPFLAGS := $(CPPFLAGS) $(AUTO_CPPFLAGS)
ALL_CFLAGS := $(CFLAGS) $(AUTO_CFLAGS) $(OPT_CFLAGS)
ALL_LDFLAGS := $(LDFLAGS) $(AUTO_LDFLAGS) $(OPT_LDFLAGS)
ALL_LINKFLAGS := $(CFLAGS) $(AUTO_CFLAGS) $(LDFLAGS) $(AUTO_LDFLAGS) $(OPT_LDFLAGS)
ALL_LDLIBS := $(LDLIBS) $(AUTO_LDLIBS)

DEPS := $(LIB_OBJS:.o=.d) $(APP_OBJS:.o=.d) $(MAIN_OBJ:.o=.d) $(MIMALLOC_OBJ:.o=.d)

OPT_STAMP := $(OBJ_DIR)/.opt-flags
OPT_FINGERPRINT := $(OPT_CFLAGS) $(OPT_LDFLAGS)

.PHONY: all clean examples test upstream-compare release-lto pgo pgo-generate pgo-train pgo-use FORCE

all: $(PROG)

$(LIB_TARGET): $(LIB_OBJS)
	$(AR) rcs $@ $^

$(PROG): $(LIB_TARGET) $(MIMALLOC_OBJ) $(APP_OBJS) $(MAIN_OBJ)
	$(CC) $(ALL_LINKFLAGS) $(MIMALLOC_OBJ) $(APP_OBJS) $(MAIN_OBJ) $(LIB_TARGET) -o $@ $(ALL_LDLIBS)

examples: $(LIB_TARGET)
	$(MAKE) -C examples

test: $(PROG) $(BUILD_DIR)/test-unmap-regions
	@mkdir -p $(BUILD_DIR)
	$(CC) $(ALL_CFLAGS) $(ALL_CPPFLAGS) $(INCLUDES) tests/test_sam_alt_hla.c $(LIB_TARGET) -o $(BUILD_DIR)/test_sam_alt_hla $(ALL_LDLIBS)
	$(BUILD_DIR)/test_sam_alt_hla
	$(CC) $(ALL_CFLAGS) $(ALL_CPPFLAGS) $(INCLUDES) tests/test-human-tiebreak.c $(LIB_TARGET) -o $(BUILD_DIR)/test-human-tiebreak $(ALL_LDLIBS)
	$(BUILD_DIR)/test-human-tiebreak
	$(BUILD_DIR)/test-unmap-regions tests/data/unmap_regions.38.tsv
	./bench/run-human-benchmark.py --out-dir .context/human-benchmark
	sh tests/test-problematic-mask.sh
	sh tests/test-single-end-fast.sh

upstream-compare:
	./bench/compare-upstream.sh

$(BUILD_DIR)/test-unmap-regions: tests/test-unmap-regions.c src/unmap.c src/kommon.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDES) -MMD -MP $^ -o $@ -lz -lm

release-lto:
	$(MAKE) clean
	$(MAKE) lto=1 all

pgo:
	@if [ -z "$(PGO_TRAIN_CMD)" ]; then \
		echo "PGO_TRAIN_CMD is required, for example:"; \
		echo "  PGO_TRAIN_CMD='./minibwa map -t1 -P -o /dev/null ref reads.fastq' make pgo"; \
		exit 2; \
	fi
	$(MAKE) pgo-generate
	$(MAKE) pgo-train
	$(MAKE) pgo-use

pgo-generate:
	$(MAKE) clean
	$(MAKE) pgo=generate all

pgo-train:
	@if [ -z "$(PGO_TRAIN_CMD)" ]; then \
		echo "PGO_TRAIN_CMD is required, for example:"; \
		echo "  PGO_TRAIN_CMD='./minibwa map -t1 -P -o /dev/null ref reads.fastq' make pgo-train"; \
		exit 2; \
	fi
	$(PGO_TRAIN_CMD)

pgo-use:
	@if [ "$(PGO_MERGE_PROFILE)" = "1" ]; then \
		if ! command -v "$(LLVM_PROFDATA)" >/dev/null 2>&1; then \
			echo "llvm-profdata is required to merge clang PGO profiles"; \
			exit 2; \
		fi; \
		if ! ls "$(PGO_PROFILE_DIR)"/*.profraw >/dev/null 2>&1; then \
			echo "no clang PGO profiles found in $(PGO_PROFILE_DIR)"; \
			exit 2; \
		fi; \
		"$(LLVM_PROFDATA)" merge -output="$(PGO_PROFILE_DATA)" "$(PGO_PROFILE_DIR)"/*.profraw; \
	fi
	@find "$(BUILD_DIR)" \( -name '*.o' -o -name '*.d' \) -delete 2>/dev/null || true
	rm -f $(PROG) $(LIB_TARGET)
	$(MAKE) pgo=use all

$(OPT_STAMP): FORCE
	@mkdir -p $(dir $@)
	@if [ ! -f "$@" ] || [ "$$(cat "$@" 2>/dev/null)" != "$(OPT_FINGERPRINT)" ]; then \
		printf '%s' "$(OPT_FINGERPRINT)" > "$@"; \
	fi

$(OBJ_DIR)/third_party/mimalloc/static.o: third_party/mimalloc/static.c $(OPT_STAMP)
	@mkdir -p $(dir $@)
	$(CC) -c $(MIMALLOC_CFLAGS) $(OPT_CFLAGS) -MMD -MP -Ithird_party/mimalloc $< -o $@

$(OBJ_DIR)/%.o: %.c $(OPT_STAMP)
	@mkdir -p $(dir $@)
	$(CC) -c $(ALL_CFLAGS) $(ALL_CPPFLAGS) $(INCLUDES) -MMD -MP $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(PROG) $(LIB_TARGET) a.out *.o *.a *.dSYM *~
	$(MAKE) -C examples clean

-include $(DEPS)
