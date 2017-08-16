# -*- makefile -*-
export OPENGL=0
include rebound/Makefile.defs

all: librebound
	@echo ""
	@echo "[planetplanet] Compiling C source..."
	$(CC) $(FLAGS1) -Irebound/ $(OPT) $(PREDEF) orbit.c numerical/complex.c numerical/roots.c eyeball.c ppo.c progress/progress.c
	@echo "[planetplanet] Generating shared library..."
	$(CC) $(FLAGS2)libppo.so -o libppo.so orbit.o complex.o roots.o eyeball.o ppo.o progress.o rebound/rebound.o rebound/tree.o rebound/particle.o rebound/gravity.o rebound/integrator.o rebound/integrator_whfast.o rebound/integrator_whfasthelio.o rebound/integrator_ias15.o rebound/integrator_sei.o rebound/integrator_leapfrog.o rebound/integrator_hermes.o rebound/boundary.o rebound/input.o rebound/output.o rebound/collision.o rebound/communication_mpi.o rebound/display.o rebound/tools.o rebound/derivatives.o rebound/simulationarchive.o rebound/glad.o rebound/integrator_janus.o rebound/transformations.o $(REALTIME)
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

clean:
	@echo "[planetplanet] Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C rebound/ clean
	@echo "[planetplanet] Cleaning up shared library libppo.so ..."
	@-rm -vf libppo.so