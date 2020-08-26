# Host
CXX := $(XILINX_VITIS)/bin/xcpp
CXXFLAGS += -I$(XILINX_XRT)/include -Wall -fmessage-length=0 -std=c++14
CXXFLAGS += -Wno-unknown-pragmas -Wno-unused-label
LDFLAGS += -L$(XILINX_XRT)/lib -lOpenCL -lrt -lstdc++

# Kernel
xclbin_dir ?= xclbin
xclobj_dir ?= xclobj

dsa := $(basename $(notdir $(platform)))

VPP := $(XILINX_VITIS)/bin/v++
vpp_flags += -t $(build) --platform $(platform)
vpp_compile_flags +=
vpp_link_flags +=

XCLBIN += $(xclbin_dir)/$(xclbin_name).$(build).$(dsa).xclbin
XOBJS += $(patsubst %, $(xclobj_dir)/%.$(build).$(dsa).xo, $(kernels))

# for clean build
XCLBIN_ALL += $(xclbin_dir)/$(xclbin_name).*.$(dsa).xclbin*
XOBJS_ALL += $(patsubst %, $(xclobj_dir)/%.*.$(dsa).xo, $(kernels))

$(XOBJS): $(xclobj_dir)/%.$(build).$(dsa).xo: %.cpp
	@mkdir -p $(xclobj_dir)
	$(VPP) -c -k $* -o'$@' '$<' $(vpp_flags) $(vpp_compile_flags) --xp prop:kernel.$*.kernel_flags="$(kernel_flags)"

$(xclbin_dir)/$(xclbin_name).$(build).$(dsa).xclbin: $(XOBJS) $(xclbin_depends) $(rtl_kernels)
	@mkdir -p $(xclbin_dir)
	LANG= $(VPP) -l -o'$@' $(XOBJS) $(rtl_kernels) $(vpp_flags) $(vpp_link_flags)

emconfig.json:
	emconfigutil --platform $(platform)

.PHONY: device
device: $(XCLBIN)

.PHONY: host
host: all

clean_targets += clean-vitis

.PHONY: clean-vitis
clean-vitis:
	rm -rf _x .Xil .ipcache .run
	rm -rf emconfig.json sample_compile.ini sample_link.ini
	rm -rf *.log *.jou *.wcfg *.wdb
	rm -f $(XCLBIN_ALL) $(XOBJS_ALL)
	test $(xclbin_dir) != "." && rm -rf $(xclbin_dir) || true
	test $(xclobj_dir) != "." && rm -rf $(xclobj_dir) || true

.PHONY: help
help:
	@echo "make host"
	@echo "make device"
	@echo "make clean"

