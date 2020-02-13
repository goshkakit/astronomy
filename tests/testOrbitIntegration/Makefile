 #------------------------------------
# compile
CC = g++
NCC = nvcc

#------------------------------------
# flags
CFLAGS = -c -Wall -lm
NFLAGS = -c -arch=compute_20 -code=sm_20

#------------------------------------
#lib
#INCLUDE = -I/usr/local/cuda-5.0/include
#INCLUDE = -I/usr/local/cuda/include
INCLUDE = -I/common/cuda-5.0/include

#------------------------------------

SOURCESC = AdditionalFunction.cpp \
			Influence_Test.cpp \
			InfluenceAtmosphere.cpp \
			InfluenceEarthRotation.cpp \
			InfluenceEGM96.cpp \
			InfluenceForce.cpp \
			InfluenceNutationEarth.cpp \
			InfluencePlanet.cpp \
			InfluencePoleEarth.cpp \
			InfluenceSun.cpp \
			InfluenceTime.cpp \
			InfluiencePrecessionEarth.cpp \
			PredictListSetupAndCheck.cpp \
			PredictOrbitConst.cpp \
			PredictOrbitIntegrator.cpp \
			PredictOrbitSat.cpp \
			SatCloseApproach.cpp \
			PredictOrbitMod.cpp \
			PredictOrbitModTest.cpp \
			cutilNP.cpp \
			main_linux.cpp

SOURCESCU = InfluenceEGM96GPU.cu \
			InfluenceNutationEarthGPU.cu \
			InfluencePlanetGPU.cu \
			InfluencePoleEarthGPU.cu \
			PredictOrbitIntegratorGPU.cu

OBJECTSC = AdditionalFunction.o \
			Influence_Test.o \
			InfluenceAtmosphere.o \
			InfluenceEarthRotation.o \
			InfluenceEGM96.o \
			InfluenceForce.o \
			InfluenceNutationEarth.o \
			InfluencePlanet.o \
			InfluencePoleEarth.o \
			InfluenceSun.o \
			InfluenceTime.o \
			InfluiencePrecessionEarth.o \
			PredictListSetupAndCheck.o \
			PredictOrbitConst.o \
			PredictOrbitIntegrator.o \
			PredictOrbitSat.o \
			SatCloseApproach.o \
			PredictOrbitMod.o \
			PredictOrbitModTest.o \
			cutilNP.o \
			main_linux.o
			
OBJECTSCU = InfluenceEGM96GPU.o \
			InfluenceNutationEarthGPU.o \
			InfluencePlanetGPU.o \
			InfluencePoleEarthGPU.o \
			PredictOrbitIntegratorGPU.o

all: hello
    
hello: $(OBJECTSC) $(OBJECTSCU)
	$(NCC) $(OBJECTSC) $(OBJECTSCU) -o hello

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $<
	
%.o: %.cu
	$(NCC) $(NFLAGS) $<
