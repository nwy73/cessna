# Minimal Makefile for a single-file disting NT plugin

PROJECT = cessna_awg
DISTINGNT_API ?= distingNT_API
SRC = Cessna_AWG.cpp

CXX = arm-none-eabi-g++
CXXFLAGS = -std=c++11 -Wall -Wextra -fno-rtti -fno-exceptions \
  -mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard \
  -O3 -ffast-math -funroll-loops -fdata-sections -ffunction-sections \
  -I$(DISTINGNT_API)/include

BUILD_DIR = build
OUT_DIR = plugins

all: $(OUT_DIR)/$(PROJECT).o

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

$(BUILD_DIR)/$(PROJECT).o: $(SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $@

# Pack as a relocatable .o (what the NT expects)
$(OUT_DIR)/$(PROJECT).o: $(BUILD_DIR)/$(PROJECT).o | $(OUT_DIR)
	$(CXX) -r $(BUILD_DIR)/$(PROJECT).o -o $@

clean:
	rm -rf $(BUILD_DIR) $(OUT_DIR)
