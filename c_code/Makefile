###############################################################################
#    Copyright (c) 2020 Ben Cerjan
#
#    This file is part of simpleFDTD.
#
#    simpleFDTD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    simpleFDTD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with simpleFDTD.  If not, see <https://www.gnu.org/licenses/>.
###############################################################################

CXX = emcc

TARGET_DIR = ./../website/
EXTENSION = .html
_EXE = simulation
EXE = $(TARGET_DIR)$(_EXE)$(EXTENSION)
WASM = $(TARGET_DIR)$(_EXE).wasm
JS = $(TARGET_DIR)$(_EXE).js


SRC_DIR = src
OBJ_DIR = obj

SRC = $(widlcard $(SRC_DIR)/*.c)

# We don't need all of the objects as some of them are only needed
# for the binary version of the file
_OBJ = emscripten_out.o field_updates.o grid_free.o grid_initialization.o \
       imageShow.o memory_allocation.o ricker.o structure_functions.o \
       array_functions.o
OBJ = $(patsubst %, $(OBJ_DIR)/%,$(_OBJ))

HTML_TEMPLATE = ../html_templates/simulation_template.html

EXPORTS = -s "EXPORTED_FUNCTIONS=['_fdtdSim']"

CXXFLAGS += -s USE_SDL=2 -s WASM=1 -s NO_EXIT_RUNTIME=1\
           -s ERROR_ON_UNDEFINED_SYMBOLS=0 -s SAFE_HEAP=1 \
           -I ./include/fdtd

EMSDKFLAGS += $(CXXFLAGS) \
           -s "EXTRA_EXPORTED_RUNTIME_METHODS=['cwrap']"\
           --shell-file $(HTML_TEMPLATE) \
           -s ALLOW_MEMORY_GROWTH=1 \
           --output_eol linux \
           $(EXPORTS) \
           --profiling-funcs \
					 --js-library ./src/custom_js.js


#OPT = -O2

CFLAGS += -Wall $(OPT)

#LDFLAGS += -Llibs -DARMA_DONT_USE_WRAPPER

LDLIBS += -lm


.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(EMSDKFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@


#$(CXX) $(EMSDKFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CXX) $(CXXFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ)
	$(RM) $(EXE)
	$(RM) $(WASM)
	$(RM) $(JS)

clean_html:
	$(RM) $(EXE)
	$(RM) $(WASM)
	$(RM) $(JS)
