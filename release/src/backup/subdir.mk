################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/backup/exafmm.cpp \
../src/backup/load_balancer.cpp \
../src/backup/main.cpp \
../src/backup/new.cpp 

OBJS += \
./src/backup/exafmm.o \
./src/backup/load_balancer.o \
./src/backup/main.o \
./src/backup/new.o 

CPP_DEPS += \
./src/backup/exafmm.d \
./src/backup/load_balancer.d \
./src/backup/main.d \
./src/backup/new.d 


# Each subdirectory must supply rules for building sources it contributes
src/backup/%.o: ../src/backup/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -O3 -tbb -DNDEBUG -vec-report6 -xHost `pkg-config --cflags hpx_application` -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


