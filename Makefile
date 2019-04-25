bsh-denovo: main.c htslib/libhts.a
	gcc -O2 -g -Wall $^ -lpthread -lz -lm -o $@

submodules:
	cd htslib && make libhts.a ; cd ..
