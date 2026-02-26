PROGRAM_NAME := stroma_world

# Default simulation config.
CONFIG_FILE ?= config/PhysiCell_settings.xml

# Build mode: release or debug.
BUILD ?= release

# Path to a PhysiCell install tree.  Override on the command line or in the
# environment if your checkout lives somewhere else.
PHYSICELL_DIR ?= /work/$(USER)/PhysiCell

# Tell Make to search $(PHYSICELL_DIR) for any .cpp/.h it cannot find locally.
vpath %.cpp $(PHYSICELL_DIR)
vpath %.h   $(PHYSICELL_DIR)

# Compiler/toolchain.
CXX ?= g++
ARCH ?= native
OPENMP_FLAGS ?= -fopenmp

# JSON parser include: BioFVM bundles nlohmann/json.hpp as BioFVM/json.hpp.
JSON_INCLUDE ?= -I$(PHYSICELL_DIR)/BioFVM

CPPFLAGS := -I. -I$(PHYSICELL_DIR) $(JSON_INCLUDE)
STD_FLAGS := -std=c++17
COMMON_FLAGS := -march=$(ARCH) -m64 -fomit-frame-pointer -mfpmath=both $(OPENMP_FLAGS)

ifeq ($(BUILD),debug)
OPT_FLAGS := -O0 -g
else
OPT_FLAGS := -O2
endif

CXXFLAGS := $(STD_FLAGS) $(COMMON_FLAGS) $(OPT_FLAGS)
LDFLAGS :=
LDLIBS := $(OPENMP_FLAGS)

OBJ_DIR = build/obj/$(BUILD)
BUILD_BIN_DIR = build/bin/$(BUILD)
BUILD_PROGRAM = $(BUILD_BIN_DIR)/$(PROGRAM_NAME)
GENERATED_DIR := build/generated
GENERATED_MAIN := $(GENERATED_DIR)/main.cpp
TEST_CONFIG := build/config/PhysiCell_settings.test.xml

# PhysiCell/BioFVM source groups (standard project pattern).
BIOFVM_SOURCES := \
	BioFVM/BioFVM_vector.cpp \
	BioFVM/BioFVM_mesh.cpp \
	BioFVM/BioFVM_microenvironment.cpp \
	BioFVM/BioFVM_solvers.cpp \
	BioFVM/BioFVM_matlab.cpp \
	BioFVM/BioFVM_utilities.cpp \
	BioFVM/BioFVM_basic_agent.cpp \
	BioFVM/BioFVM_MultiCellDS.cpp \
	BioFVM/BioFVM_agent_container.cpp \
	BioFVM/pugixml.cpp

PHYSICELL_CORE_SOURCES := \
	core/PhysiCell_phenotype.cpp \
	core/PhysiCell_cell_container.cpp \
	core/PhysiCell_standard_models.cpp \
	core/PhysiCell_cell.cpp \
	core/PhysiCell_custom.cpp \
	core/PhysiCell_utilities.cpp \
	core/PhysiCell_constants.cpp \
	core/PhysiCell_basic_signaling.cpp \
	core/PhysiCell_signal_behavior.cpp \
	core/PhysiCell_rules.cpp

PHYSICELL_MODULE_SOURCES := \
	modules/PhysiCell_SVG.cpp \
	modules/PhysiCell_pathology.cpp \
	modules/PhysiCell_MultiCellDS.cpp \
	modules/PhysiCell_various_outputs.cpp \
	modules/PhysiCell_pugixml.cpp \
	modules/PhysiCell_settings.cpp \
	modules/PhysiCell_geometry.cpp

# Compile all custom module sources.
CUSTOM_SOURCES := $(sort $(wildcard custom_modules/*.cpp))

TEMPLATE_MAIN := $(PHYSICELL_DIR)/sample_projects/template/main.cpp
MAIN_SOURCE := $(GENERATED_MAIN)
ALL_SOURCES := $(BIOFVM_SOURCES) $(PHYSICELL_CORE_SOURCES) $(PHYSICELL_MODULE_SOURCES) $(CUSTOM_SOURCES) $(MAIN_SOURCE)
ALL_OBJECTS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(ALL_SOURCES))

.PHONY: all clean debug release test dirs build-mode

all: release

release:
	$(MAKE) BUILD=release build-mode

debug:
	$(MAKE) BUILD=debug build-mode

build-mode: $(BUILD_PROGRAM)
	cp $(BUILD_PROGRAM) $(PROGRAM_NAME)

$(PROGRAM_NAME): $(BUILD_PROGRAM)
	cp $< $@

$(BUILD_PROGRAM): $(ALL_OBJECTS) | dirs
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(ALL_OBJECTS) $(LDLIBS)

$(OBJ_DIR)/%.o: %.cpp | dirs
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(GENERATED_MAIN): $(TEMPLATE_MAIN) | dirs
	@mkdir -p $(dir $@)
	cp $< $@

dirs:
	@mkdir -p $(OBJ_DIR) $(GENERATED_DIR) build/config

clean:
	rm -rf build
	rm -f $(PROGRAM_NAME)

test: release
	@dt_diffusion=$$(awk -F'[<>]' '/<dt_diffusion[[:space:]]+units="min">/{print $$3; exit}' $(CONFIG_FILE)); \
	if [ -z "$$dt_diffusion" ]; then dt_diffusion=0.01; fi; \
	test_max_time=$$(awk -v dt="$$dt_diffusion" 'BEGIN{printf "%.10g", dt*100.0}'); \
	output_dir=build/test_output; \
	rm -rf "$$output_dir"; \
	mkdir -p "$$output_dir"; \
	awk -v max_time="$$test_max_time" -v out_dir="$$output_dir" '\
	BEGIN { replaced_folder = 0; } \
	/<max_time[[:space:]]+units="min">/ { sub(/>[^<]+</, ">" max_time "<"); } \
	(!replaced_folder && /<folder>[^<]*<\/folder>/) { \
		sub(/<folder>[^<]*<\/folder>/, "<folder>" out_dir "</folder>"); \
		replaced_folder = 1; \
	} \
	{ print }' $(CONFIG_FILE) > $(TEST_CONFIG); \
	echo "Running 100-step smoke test (dt_diffusion=$$dt_diffusion, max_time=$$test_max_time, output_dir=$$output_dir)"; \
	./$(PROGRAM_NAME) $(TEST_CONFIG); \
	if [ -z "$$(find "$$output_dir" -mindepth 1 -type f -print -quit)" ]; then \
		echo "ERROR: test failed; output directory '$$output_dir' is empty."; \
		exit 1; \
	fi; \
	echo "Test passed: output directory '$$output_dir' contains simulation outputs."
