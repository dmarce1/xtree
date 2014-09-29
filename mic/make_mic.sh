export PKG_CONFIG_PATH=/home/dmarce1/mic/lib/pkgconfig
#echo mpiicpc -O3 -mmic -mmic -w -DNO_OUTPUT -c -fmessage-length=0 `pkg-config --cflags hpx_application` -MMD -MP -MF"src/main.d" -MT"src/main.d" -o "src/main.o" "../src/main.cpp"
#     mpiicpc -O3 -mmic -mmic -w -DNO_OUTPUT -c -fmessage-length=0 `pkg-config --cflags hpx_application` -MMD -MP -MF"src/main.d" -MT"src/main.d" -o "src/main.o" "../src/main.cpp"
echo mpiicpc -O3 -mmic -mmic -w -c -fmessage-length=0 `pkg-config --cflags hpx_application` -MMD -MP -MF"src/load_balancer.d" -MT"src/load_balancer.d" -o "src/load_balancer.o" "../src/load_balancer.cpp"
     mpiicpc -O3 -mmic -mmic -w -c -fmessage-length=0 `pkg-config --cflags hpx_application` -MMD -MP -MF"src/load_balancer.d" -MT"src/load_balancer.d" -o "src/load_balancer.o" "../src/load_balancer.cpp"
echo mpiicpc  ./src/load_balancer.o ./src/main.o -mmic -L/home/dmarce1/lib `pkg-config --libs hpx_application`  -o  "xtree"
     mpiicpc  ./src/load_balancer.o ./src/main.o -mmic -L/home/dmarce1/lib `pkg-config --libs hpx_application`  -o  "xtree"

