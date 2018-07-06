g++ -Wall -I /usr/local/include/ -L/usr/local/lib -lgsl -lgslcblas -lm ./src/main.cpp ./src/Interp_Test.cpp
./a.out > ./src/data.out
gnuplot ./src/Solution.sh
open *.eps
