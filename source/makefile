objects = moult_prop.o moult_main.o moult_grid.o moult_gold.o moult_funs.o moult_sim.o rand.o
mode = -g  -Wall -ffloat-store
#macro = -DUNDRAW
compiler = gcc


moult	 : $(objects)
	$(compiler) $(mode) -o moult $(objects) -lm
	if [ -d ../build ]; then \
	   mv moult ../build ; \
	else \
	     mv moult .. ;\
	fi
	rm -f $(objects)

moult_prop.o : moult_prop.c moult.h 
	$(compiler) $(mode) -c moult_prop.c 
moult_main.o : moult_main.c moult.h 
	$(compiler) $(mode) -c moult_main.c 
moult_grid.o : moult_grid.c moult.h 
	$(compiler) $(mode) -o moult_grid.o -c moult_grid.c 
moult_gold.o : moult_gold.c moult.h 
	$(compiler) $(mode) -c moult_gold.c 
moult_funs.o : moult_funs.c moult.h 
	$(compiler) $(mode) -c moult_funs.c 
moult_sim.o : moult_sim.c moult.h 
	$(compiler) $(mode) -c moult_sim.c 
rand.o : rand.c
	$(compiler) $(mode) -c rand.c 


.PHONY : clean
clean :
	-rm $(objects)
