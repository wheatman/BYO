include ../makefile.variables

INCLUDE_DIRS = -I../ -I../external/parlaylib/include/ -I../external/ -I../external/abseil-cpp/ -I../external/SSTGraph/include/ -I../external/ParallelTools/ -I../external/StructOfArrays/include/ -I../external/dhb/include/ -I../external/reducer/ -I../external/terrace/include/

OBJDIR = ../bin/gbbs/

ALL_PRE = graph_io io
ALL= $(addprefix $(OBJDIR), $(addsuffix .o, $(ALL_PRE))) $(addprefix $(OBJDIR), $(addsuffix .a, $(ALL_PRE)))

all: $(ALL)

$(OBJDIR)graph_io.o : graph_io.cc $(addprefix $(OBJDIR), io.o)
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.o : %.cc
	$(CC) $(INCLUDE_DIRS) $(CFLAGS) $(PFLAGS) -c $< -o $@

$(OBJDIR)%.a : $(OBJDIR)%.o
	ar -rcs $@ $<

.PHONY : clean

clean :
	rm -f *.o $(ALL)
