LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/cuda/lib64

CUDA_SDK_PATH = /usr/local/cuda
NVCC = $(CUDA_SDK_PATH)/bin/nvcc -ccbin g++-6

all: initialise runInit simulate runSim clean move

energies: energies.cpp
	g++ -o $@ $+ -I../

runEnergies: energies
	./energies

initialise.o:initialise.cu
	$(NVCC) -o $@ -c $< -I../

initialise: initialise.o
	$(NVCC) -o $@ $+ -lcufft -I../

runInit: initialise
	LD_LIBRARY_PATH=$(CUDA_SDK_PATH)/lib64 ./initialise

simulate.o:simulate.cu
	$(NVCC) -o $@ -c $< -I../

simulate: simulate.o
	$(NVCC) -o $@ $+ -I../

runSim: simulate
	./simulate

clean:
	rm -f simulate initialise simulate.o initialise.o energies

move:
	mv ./*.bin ./DataFiles
