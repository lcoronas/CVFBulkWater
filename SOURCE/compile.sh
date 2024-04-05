nvcc -arch=sm_86 --ptxas-options=--verbose --use_fast_math -O2 -o gpu_water3D gpu_water3D.cu mersenne_inline.cu -lm
