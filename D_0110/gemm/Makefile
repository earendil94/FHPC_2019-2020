###
MKL= -L${MKLROOT}/lib/intel64  -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

all: gemm.x gemm_gpu.x

cpu: gemm.x

gpu: gemm_gpu.x

#This -m64 is a suggestion for the mkl version to use
gemm.x: gemm.c
	gcc -DUSE_FLOAT $^ -m64 -I${MKLROOT}/include $(MKL)  -o $@

gemm_gpu.x: gemm_gpu.c
	nvcc -DUSE_FLOAT $^ -m64 -I${MKLROOT}/include $(MKL)  -lcublas  -o $@

clean:
	rm -rf *.x
