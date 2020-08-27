.PHONY: csim
csim:
	rm -f out.ppm
	@echo "set csim 1" > gen.tcl
	@echo "set csynth 0" >> gen.tcl
	@echo "source hls.tcl" >> gen.tcl
	@AP_GCC_PATH=/usr/bin vitis_hls -f gen.tcl

.PHONY: csynth
csynth:
	@echo "set csim 0" > gen.tcl
	@echo "set csynth 1" >> gen.tcl
	@echo "source hls.tcl" >> gen.tcl
	@vitis_hls -f gen.tcl

.PHONY: compile
compile:
	g++ -o tb \
		tb.cpp \
		../kernel/rt.cpp \
		-O0 \
		-I../kernel -I../external/hlslib/include -I$(XILINX_VITIS)/include \
		-std=c++17 \
		-L$(XILINX_VITIS)/lnx64/lib/csim \
		-Wl,-rpath=$(XILINX_VITIS)/lnx64/lib/csim \
		-lhlsmc++-CLANG39 \
		-L$(XILINX_VITIS)/lnx64/tools/fpo_v7_0 \
		-Wl,-rpath=$(XILINX_VITIS)/lnx64/tools/fpo_v7_0 \
		-lIp_floating_point_v7_0_bitacc_cmodel \
		-lgmp -lmpfr \
		-pthread
	#LD_LIBRARY_PATH=$(LD_LIBARY_PATH):$(XILINX_VITIS)/lnx64/tools/fpo_v7_0 ./tb

.PHONY: clean
clean:
	rm -rf prj gen.tcl out.ppm *.log
