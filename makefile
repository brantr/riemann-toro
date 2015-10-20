EXEC   = riemann

OPTIMIZE =  -O2

OPT   = $(OPTIMIZE)

OBJS   = main.o routines.o riemann_toro.o

CC     = clang++

INCL   = routines.h riemann_toro.h

LIBS   = -lm -lgsl -lgslcblas

CFLAGS = $(OPT)

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

