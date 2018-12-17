OBJ   = nrg.o nrgmedid.o nrgupdat.o nrginici.o nrgio.o
#CFLAGS= -g
CFLAGS= -O3 -funroll-loops -finline-functions 

iso: $(OBJ)
	gcc $(CFLAGS) $(OBJ) -lm -o u1
.c.o:
	gcc $(DEFINES)  -c $(CFLAGS) $<

clean: 
	/bin/rm -f $(OBJ) is

#		*Individual File Dependencies*
nrg.o: nrg.c nrg.h Makefile

nrgmedid.o: nrgmedid.c nrg.h Makefile

nrgupdat.o: nrgupdat.c nrg.h Makefile

nrginici.o: nrginici.c nrg.h Makefile

nrgio.o: nrgio.c nrg.h Makefile












