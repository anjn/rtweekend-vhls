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
