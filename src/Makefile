IDIR =../include
CFLAGS=-I$(IDIR)
LDFLAGS= -shared
CPPFLAGS=-fPIC
ODIR=obj
# ODIR_FULL = obj_full
LDIR =../lib
BIN = ../bin
# LTARGET_FULL=$(LDIR)/fullmodel.so
# TARGET_FULL=$(BIN)/fullmodel
LTARGET=$(LDIR)/jetmodel.so
TARGET=$(BIN)/jetmodel
LIBS=-lm -lgsl -lgslcblas

_DEPS = jet_fns.h mtwister.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = func_model.o jet_fns.o mtwister.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# _OBJ_FULL = full_model.o jet_fns.o mtwister.o 
# OBJ_FULL = $(patsubst %,$(ODIR_FULL)/%,$(_OBJ_FULL))

all: $(LTARGET) $(TARGET) #$(LTARGET_FULL) $(TARGET_FULL)
debug: CFLAGS += -DDEBUG -g
debug: all

parallel: CFLAGS += -fopenmp
parallel: all

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(CPPFLAGS) 

# $(ODIR_FULL)/%.o: %.c $(DEPS)
# 	$(CC) -c -o $@ $< $(CFLAGS) $(CPPFLAGS) 

$(LTARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDFLAGS) 

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) 

# $(LTARGET_FULL): $(OBJ_FULL)
# 	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(LDFLAGS) 

# $(TARGET_FULL): $(OBJ_FULL)
# 	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 