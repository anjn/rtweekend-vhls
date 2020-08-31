.PHONY: csim
csim:
	rm -f out.ppm
	@echo "set csim 1" > gen.tcl
	@echo "set csyn 0" >> gen.tcl
	@echo "source hls.tcl" >> gen.tcl
	@AP_GCC_PATH=/usr/bin vitis_hls -f gen.tcl
	@gnome-open out.ppm

.PHONY: csyn
csyn:
	@echo "set csim 0" > gen.tcl
	@echo "set csyn 1" >> gen.tcl
	@echo "source hls.tcl" >> gen.tcl
	@vitis_hls -f gen.tcl \
		| sed -ru "s/^ERROR.*\$$/\o033[1;31m\0\o033[0m/g" \
		| sed -ru "s/^.*Unable to pipeline.*\$$/\o033[1;31m\0\o033[0m/g" \
		| sed -ru "s/^.*II Violation.*\$$/\o033[1;31m\0\o033[0m/g" \
		| sed -ru "s/^WARNING.*\$$/\o033[1;33m\0\o033[0m/g" \
		| sed -ru "s/^.*Final II.*\$$/\o033[1;34m\0\o033[0m/g" \
		| sed -ru "s/^.*Burst read of variable length and bit width.*\$$/\o033[1;34m\0\o033[0m/g" \
		| sed -ru "s/^.*Burst write of variable length and bit width.*\$$/\o033[1;34m\0\o033[0m/g" \
		| sed -ru "s/^.*Inferring multiple bus burst.*\$$/\o033[1;34m\0\o033[0m/g"

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
