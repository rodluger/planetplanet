# -*- makefile -*-
export OPENGL=0
include rebound/Makefile.defs

all: librebound
	@echo ""
	@echo "[planetplanet] Compiling C source..."
	$(CC) -fPIC -Irebound/ $(OPT) $(PREDEF) -c orbit.c complex.c roots.c eyeball.c ppo.c progress.c
	@echo "[planetplanet] Generating shared library..."
	$(CC) -shared -Wl,-install_name,libppo.so -o libppo.so orbit.o complex.o roots.o eyeball.o ppo.o progress.o -L. -lrebound librebound.so
	@-rm orbit.o
	@-rm complex.o
	@-rm roots.o
	@-rm eyeball.o
	@-rm ppo.o
	@-rm progress.o
	@echo ""
	@echo "[planetplanet] Installation successful."

librebound: 
	@echo "[planetplanet] Compiling REBOUND..."
	$(MAKE) -C rebound/
	@-rm -f librebound.so
	@ln -s rebound/librebound.so .

clean:
	@echo "[planetplanet] Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C rebound/ clean
	@echo "[planetplanet] Cleaning up shared library libppo.so ..."
	@-rm -vf libppo.so