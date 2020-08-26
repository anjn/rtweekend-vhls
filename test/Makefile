.PHONY: csim
csim:
	@echo "set csim 1" > gen.tcl
	@echo "set csynth 0" >> gen.tcl
	@echo "source hls.tcl" >> gen.tcl
	@vitis_hls -f gen.tcl

.PHONY: csynth
csynth:
	@echo "set csim 0" > gen.tcl
	@echo "set csynth 1" >> gen.tcl
	@echo "source hls.tcl" >> gen.tcl
	@vitis_hls -f gen.tcl

.PHONY: compile
compile:
	g++ -o /dev/null \
		tb.cpp \
		../kernel/rt.cpp \
		-I../kernel -I../external/hlslib/include -I$(XILINX_VITIS)/include \
		-std=c++17 \
		-L$(XILINX_VITIS)/lnx64/lib/csim \
		-lhlsmc++-CLANG39 \
		-pthread

.PHONY: clean
clean:
	rm -rf prj gen.tcl out.ppm *.log
