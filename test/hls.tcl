open_project -reset prj
set cflags "-std=c++17 -I../kernel -I../external/hlslib/include"
if {$csynth == 1} {
  set cflags "$cflags -DHLSLIB_SYNTHESIS"
}
add_files -cflags $cflags ../kernel/rt.cpp
add_files -cflags $cflags -tb tb.cpp
set_top rt
open_solution -flow_target vitis solution
set_part xcu200-fsgd2104-2-e
create_clock -period 400MHz -name default
if {$csim == 1} {
  csim_design
}
if {$csynth == 1} {
  csynth_design
}
exit
