IDIR =../include
CFLAGS=-I$(IDIR)
LDFLAGS= -shared
CPPFLAGS=-fPIC
ODIR=obj
LDIR =../lib
# BIN = ../bin
LTARGET=$(LDIR)/jetmodel.so
# TARGET=$(BIN)/jetmodel
LIBS=-lm -lgsl -lgslcblas

_DEPS = jet_fns.h mtwister.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = func_model.o jet_fns.o mtwister.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(LTARGET)
debug: CFLAGS += -DDEBUG -g
debug: all

parallel: CFLAGS += -fopenmp
parallel: all

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(CPPFLAGS) 

$(LTARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDFLAGS) 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
