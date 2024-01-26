call vitis_hls run_hls.tcl
find "#define MATRIX_" matrixvector.h
jsdb results.js mv_kernel solution1 matrixvector
rem get jsdb from jsdb.org -- sorry, Windows only, doesn't compile on 64-bit linux
rem type mv_kernel\solution1\sim\report\matrixvector_cosim.rpt
