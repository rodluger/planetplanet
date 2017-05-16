# -*- makefile -*-

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
GCC_FLAGS1 = -fPIC -Wl,-Bsymbolic-functions -c -O3
GCC_FLAGS2 = -shared -O3 -Wl,-Bsymbolic-functions,-soname,ppo.so
endif
ifeq ($(UNAME_S),Darwin)
GCC_FLAGS1 = -fPIC -c
GCC_FLAGS2 = -shared -Wl,-install_name,ppo.so
endif

GCC = gcc

.PHONY: all
.SILENT: all

all:
	echo "[planetplanet] Compiling C source code..."
	${GCC} ${GCC_FLAGS1} orbit.c complex.c roots.c eyeball.c ppo.c
	echo "[planetplanet] Generating shared library..."
	gcc ${GCC_FLAGS2} -o ppo.so orbit.o complex.o roots.o eyeball.o ppo.o -lc
	rm orbit.o
	rm complex.o
	rm roots.o
	rm eyeball.o
	rm ppo.o
	echo "[planetplanet] Install successful."