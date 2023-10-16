CFLAGS = -O3

all: evm

evm: main.o classes.o functions.o
		g++ -g $^ $(CFLAGS) -o $@

main.o: main.cpp
		g++ -c $(CFLAGS) $^

classes.o: classes.cpp
		g++ -c $(CFLAGS) $^

functions.o: functions.cpp
		g++ -c $(CFLAGS) $^

clear:
		rm -rf *.o evm
clear_evm:
		rm -rf *.o