################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/load_balancer.cpp \
../src/main.cpp 

OBJS += \
./src/load_balancer.o \
./src/main.o 

CPP_DEPS += \
./src/load_balancer.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	icpc -DHPX_APPLICATION_EXPORTS -DNDEBUG -DHPX_ENABLE_ASSERT_HANDLER -I/home/dmarce1/include/hpx/external -O3 -c  -std=c++11 `pkg-config --cflags hpx_application` -vec-report=1 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


