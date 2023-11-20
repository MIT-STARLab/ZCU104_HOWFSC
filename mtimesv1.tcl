
# https://docs.xilinx.com/r/en-US/ug1399-vitis-hls/set_directive_bind_op
# Create a project
open_project -reset mtimesv1

# Add design files
add_files mtimesv1.cpp
# Add test bench & files
add_files -tb mtimesv_test.cpp

# Set the top-level function
set_top matrixTimesVector1

# ########################################################
# Create a solution, Define technology and clock rate
open_solution -reset sln1 -flow_target vitis
set_part  {xczu7ev-ffvc1156-2-e}
create_clock -period 25

# csim_design

# Set any optimization directives
# cyclic is 2/3 the time of block 
# 32 is 1/2 the time of 16
# 128 x 32 takes about 2100 clock cycles
# pipelining LOOP_M also pipelines LOOP_N
set_directive_array_partition -factor 16 -type cyclic matrixTimesVector1 A
set_directive_array_partition -factor 16 -type cyclic matrixTimesVector1 B
set_directive_array_partition -factor 16 -type cyclic matrixTimesVector1 C
set_directive_pipeline matrixTimesVector1/LOOP_N
# ap_memory or s_axilite
set_directive_interface -mode ap_memory matrixTimesVector1 A 
set_directive_interface -mode ap_memory matrixTimesVector1 B
set_directive_interface -mode ap_memory matrixTimesVector1 C
# End of directives
# seems to get stuck at "Starting scheduling". Maybe ap_memory doesn't get stuck?
csynth_design
cosim_design
export_design

exit
