# -*- makefile -*-
export OPENGL=0
include ../../rebound/src/Makefile.defs

all: librebound
	@echo ""
	@echo "[planetplanet] Compiling C source..."
	$(CC) -fPIC -I../../src/ $(OPT) $(PREDEF) -c orbit.c complex.c roots.c eyeball.c ppo.c
	@echo "[planetplanet] Generating shared library..."
	$(CC) -shared -Wl,-install_name,libppo.so -o libppo.so orbit.o complex.o roots.o eyeball.o ppo.o -L. -lrebound $(LIB)
	@-rm orbit.o
	@-rm complex.o
	@-rm roots.o
	@-rm eyeball.o
	@-rm ppo.o
	@echo ""
	@echo "[planetplanet] Installation successful."

librebound: 
	@echo "[planetplanet] Compiling REBOUND..."
	$(MAKE) -C ../../rebound/src/
	@-rm -f librebound.so
	@ln -s ../../rebound/src/librebound.so .

clean:
	@echo "[planetplanet] Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../rebound/src/ clean
	@echo "[planetplanet] Cleaning up shared library libppo.so ..."
	@-rm -vf libppo.so