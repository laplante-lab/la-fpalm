GCC = g++
CC = gcc
NVCC = nvcc
OPT = -O3 -s -DNDEBUG
WARN = -Wall 
IGNORE_WARN = -Wno-attributes -Wno-unused-result 
ERR = #-Werror

USER = jjravi

MATIO_DIR = -I/home/$(USER)/matio/include/

INCLUDES = -I/home/$(USER)/dip/include/ -I/usr/local/include/ -I/usr/local/cuda/include/ $(MATIO_DIR) -I./libs/nd2

SER_VAR	= -D SERIAL
OMP_VAR	= -D OMP 
OMP_OPT = -fopenmp

SER_DIR = SERIAL
OMP_DIR = OMP

OMPFLAGS = $(OPT) $(WARN) $(ERR) $(OMP_OPT) $(OMP_VAR)
SERFLAGS = $(OPT) $(WARN) $(ERR) $(SER_VAR)
#NV1FLGS = $(OPT) $(WARN) $(ERR) $(NVC_VAR)
#NVFLAGS	= -Xcompiler "-g" --gpu-architecture=compute_50 --gpu-code=compute_50,sm_50 -O3
NVFLAGS	= -Xcompiler "-s -DNDEBUG" --gpu-architecture=compute_50 --gpu-code=compute_50,sm_50 -O3
#NVFLAGS	= -gencode arch=compute_50,code=sm_50
#NVFLAGS	= -gencode arch=compute_50,code=sm_50 -O3
#NVFLAGS	= -arch=sm_62 -g -G -O0
NVCFLGS	= -I/usr/local/cuda/include -lcuda $(NVFLAGS)

CPPVER = -std=c++17

CPPFLAGS = $(CPPVER) ${INCLUDES}
CFLAGS =  ${INCLUDES}

SER_SRC = main.cpp ComputeManager.cpp DataIO.cpp noise_map/NoiseMap.cpp tools/plot.cpp tools/srhist_color.c tools/convolution.cpp tools/cHistRecon.c tools/cMakeSubregions.c gaussmle/kernelCall.cpp gaussmle/wrapper.cu gaussmle/kernel.o SRscmos.cpp tools/savedata.c
SER_OBJ	= $(addprefix $(SER_DIR)/, main.o ComputeManager.o DataIO.o NoiseMap.o plot.o srhist_color.o convolution.o cHistRecon.o cMakeSubregions.o kernelCall.o wrapper.o kernel.o SRscmos.o savedata.o)

OMP_SRC = main.cpp ComputeManager.cpp DataIO.cpp noise_map/NoiseMap.cpp tools/plot.cpp tools/srhist_color.c tools/convolution.cpp tools/cHistRecon.c tools/cMakeSubregions.c gaussmle/kernelCall.cpp gaussmle/wrapper.cu gaussmle/kernel.o SRscmos.cpp savedata.c
OMP_OBJ	= $(addprefix $(OMP_DIR)/, main.o ComputeManager.o DataIO.o NoiseMap.o plot.o srhist_color.o convolution.o cHistRecon.o cMakeSubregions.o kernelCall.o wrapper.o kernel.o SRscmos.o savedata.o)

all: serial omp
	@echo "Compilation Done ---> nothing else to make :) "

serial: CPPFLAGS += $(SERFLAGS)
serial: CFLAGS += $(SERFLAGS)
serial: serial_temp

#omp: CFLAGS = $(OMPFLAGS)
omp: CPPFLAGS += $(OMPFLAGS)
omp: CFLAGS += $(OMPFLAGS)
omp: omp_temp

serial_temp: $(SER_OBJ)
	#$(GCC) $(SERFLAGS) -o srscmos_serial $(SER_OBJ) -lm
	g++ -fPIC -o serial -fopenmp $(SER_OBJ) "-Wl,-rpath=/home/$(USER)/dip/lib" -L/home/$(USER)/dip/lib/ -lDIP -lstdc++fs -L/usr/local/cuda/lib64 -lcudart -ltiff "-Wl,-rpath,/home/$(USER)/matio/lib" -L/home/$(USER)/matio/lib -lmatio "-Wl,-rpath,./libs/nd2" -L./libs/nd2 -lNd2ReadSdk -lNd2File -pthread
	@echo "----------------------------------------------------------"
	@echo "---------C/C++ SRSCMOS SERIAL VERSION - John Ravi---------"
	@echo "----------------------------------------------------------"
	
omp_temp: $(OMP_OBJ)
	g++ -fPIC -o la-fpalm -fopenmp $(OMP_OBJ) "-Wl,-rpath=/home/$(USER)/dip/lib" -L/home/$(USER)/dip/lib/ -lDIP -lstdc++fs -L/usr/local/cuda/lib64 -lcudart -ltiff "-Wl,-rpath,/home/$(USER)/matio/lib" -L/home/$(USER)/matio/lib -lmatio "-Wl,-rpath,./libs/nd2" -L./libs/nd2 -lNd2ReadSdk -lNd2File -pthread
	@echo "----------------------------------------------------------"
	@echo "---------C/C++ SRSCMOS OPENMP VERSION - John Ravi---------"
	@echo "----------------------------------------------------------"
#$(GCC) $(OMPFLAGS) -o ocean_sim_omp $(OMP_OBJ) -lm

$(SER_DIR)/%.o: %.cpp
	@mkdir -p $(SER_DIR)
	$(GCC) $(IGNORE_WARN) $(CPPFLAGS)  -c $< -o $@

#SERIAL/%.o: noise_map/NoiseMap.cpp tools/mode.cpp
SERIAL/SRscmos.o: SRscmos.cpp 
	@mkdir -p $(SER_DIR)
	$(GCC) $(IGNORE_WARN) -Wno-unused-but-set-variable $(CPPFLAGS)  -c $< -o $@
	
SERIAL/NoiseMap.o: noise_map/NoiseMap.cpp
	@mkdir -p $(SER_DIR)
	$(GCC) $(CPPFLAGS)  -c $< -o $@

SERIAL/plot.o: tools/plot.cpp
	@mkdir -p $(SER_DIR)
	$(GCC) $(CPPFLAGS)  -c $< -o $@
	
SERIAL/savedata.o: tools/savadata.c
	@mkdir -p $(SER_DIR)
	$(CC) $(IGNORE_WARN) $(CFLAGS) -c $< -o $@

SERIAL/srhist_color.o: tools/srhist_color.c
	@mkdir -p $(SER_DIR)
	$(CC) $(IGNORE_WARN) $(CFLAGS) -c $< -o $@
	
SERIAL/convolution.o: tools/convolution.cpp
	@mkdir -p $(SER_DIR)
	$(GCC) $(CPPFLAGS)  -c $< -o $@
	
SERIAL/cHistRecon.o: tools/cHistRecon.c
	@mkdir -p $(SER_DIR)
	$(CC) $(CFLAGS)  -c $< -o $@

SERIAL/cMakeSubregions.o: tools/cMakeSubregions.c
	@mkdir -p $(SER_DIR)
	$(CC) $(CFLAGS)  -c $< -o $@
	
SERIAL/kernelCall.o: gaussmle/kernelCall.cpp 
	@mkdir -p $(SER_DIR)
	$(NVCC) $(NVCFLGS)  -c $< -o $@

SERIAL/wrapper.o: gaussmle/wrapper.cu
	@mkdir -p $(SER_DIR)
	$(NVCC) $(NVCFLGS)  -c $< -o $@

SERIAL/kernel.o: gaussmle/kernel.cu
	@mkdir -p $(SER_DIR)
	$(NVCC) $(NVCFLGS)  -c $< -o $@

##############

$(OMP_DIR)/main.o: main.cpp
	@mkdir -p $(OMP_DIR)
	$(GCC) $(CPPFLAGS)  -c $< -o $@

$(OMP_DIR)/ComputeManager.o: ComputeManager.cpp
	$(GCC) $(CPPFLAGS)  -c $< -o $@

$(OMP_DIR)/DataIO.o: DataIO.cpp
	$(GCC) $(CPPFLAGS)  -c $< -o $@

$(OMP_DIR)/SRscmos.o: SRscmos.cpp
	@mkdir -p $(OMP_DIR)
	$(GCC) $(IGNORE_WARN) -Wno-unused-but-set-variable $(CPPFLAGS)  -c $< -o $@

OMP/NoiseMap.o: noise_map/NoiseMap.cpp
	@mkdir -p $(SER_DIR)
	$(GCC) $(CPPVER) $(OPT) $(WARN) $(ERR) $(MATIO_DIR) -c $< -o $@

OMP/plot.o: tools/plot.cpp
	@mkdir -p $(SER_DIR)
	$(GCC) $(WARN) $(CPPVER) $(OPT) $(MATIO_DIR) $(INCLUDES)  -c $< -o $@
	
OMP/savedata.o: tools/savedata.c
	@mkdir -p $(SER_DIR)
	$(CC) $(IGNORE_WARN) $(CFLAGS) -c $< -o $@

OMP/srhist_color.o: tools/srhist_color.c
	@mkdir -p $(SER_DIR)
	$(CC) $(IGNORE_WARN) $(CFLAGS) -c $< -o $@

OMP/convolution.o: tools/convolution.cpp
	@mkdir -p $(SER_DIR)
	$(GCC) $(CPPVER) $(WARN) $(OPT)  -c $< -o $@
	
OMP/cHistRecon.o: tools/cHistRecon.c
	@mkdir -p $(SER_DIR)
	$(CC) $(WARN) $(OPT)  -c $< -o $@

OMP/cMakeSubregions.o: tools/cMakeSubregions.c
	@mkdir -p $(SER_DIR)
	$(CC) $(WARN) $(OPT)  -c $< -o $@
	
OMP/kernelCall.o: gaussmle/kernelCall.cpp 
	@mkdir -p $(SER_DIR)
	$(NVCC) $(NVCFLGS)  -c $< -o $@

OMP/wrapper.o: gaussmle/wrapper.cu
	@mkdir -p $(SER_DIR)
	$(NVCC) $(NVCFLGS)  -c $< -o $@

OMP/kernel.o: gaussmle/kernel.cu
	@mkdir -p $(SER_DIR)
	$(NVCC) $(NVCFLGS)  -c $< -o $@

clean:
	rm -rf SERIAL OMP
	rm -rf serial omp la-fpalm *.o
	rm -rf noise_map/*.o
	rm -rf tools/*.o
	rm -rf gaussmle/*.o
	rm -rf output/*
