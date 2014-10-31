################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/exafmm.cpp \
../src/load_balancer.cpp \
../src/main.cpp \
../src/new.cpp 

OBJS += \
./src/exafmm.o \
./src/load_balancer.o \
./src/main.o \
./src/new.o 

CPP_DEPS += \
./src/exafmm.d \
./src/load_balancer.d \
./src/main.d \
./src/new.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -O3 -tbb -DNDEBUG -vec-report6 -xHost `pkg-config --cflags hpx_application` -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


