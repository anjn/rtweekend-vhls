# Host
targets := main
objects := ocl_common xcl2
include_dirs := host kernel
bin_dir := .

#CXX := g++
#CXXFLAGS += -std=c++17
CXXFLAGS += -DOPENCL
LDFLAGS += -lxrt_core -pthread -luuid -lxilinxopencl

vpath %.cpp host kernel

# Kernel
xclbin_name := rt

platform := xilinx_u200_xdma_201830_2
kernels := rt
rtl_kernels := 

build ?= hw

kernel_flags += -std=c++17
vpp_flags += --kernel_frequency 400 --save-temps
vpp_flags += -Ikernel -Iexternal/hlslib/include
vpp_flags += -DHLSLIB_SYNTHESIS
vpp_flags += -g
vpp_flags += --profile_kernel data:all:all:all
vpp_compile_flags +=
vpp_link_flags += --report_level 1
vpp_link_flags += --config link.ini
#vpp_link_flags += --optimize 3
#vpp_link_flags += --advanced.prop prop:solution.kernel_compiler_margin=0.625

#vpp_link_flags += --xp vivado_prop:run.impl_1.strategy=Flow_RunPostRoutePhysOpt
#vpp_link_flags += --xp vivado_prop:run.impl_1.{STEPS.OPT_DESIGN.ARGS.DIRECTIVE}={Explore}
#vpp_link_flags += --xp vivado_prop:run.impl_1.{STEPS.PLACE_DESIGN.ARGS.DIRECTIVE}={Explore}
#vpp_link_flags += --xp vivado_prop:run.impl_1.{STEPS.PHYS_OPT_DESIGN.ARGS.DIRECTIVE}={Explore}
#vpp_link_flags += --xp vivado_prop:run.impl_1.{STEPS.ROUTE_DESIGN.ARGS.DIRECTIVE}={Explore}

include vitis.mk
include host.mk

