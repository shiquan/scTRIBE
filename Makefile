all:
	gcc -o TRIBE misc/TRIBE.c ../PISA/src/liba.a ../PISA/third_party/htslib-1.10.2/libhts.a -I../PISA/src -I../PISA/third_party/htslib-1.10.2 -lz -lcurl -lbz2 -llzma -pthread

