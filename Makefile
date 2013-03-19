#CC=/Users/xiaocai/software/mpich-3.0.2/bin/mpic++
CC=mpic++ -g
#CC=CC
CFLAGS=-DDEBUG
#CFLAGS=-std=c++11
#INCDIR=-I/Users/xiaocai/software/mpich-3.0.2/include/ \
       -I/Users/xiaocai/software/parallel-netcdf-1.3.1/include      
#LIBDIR=-L/Users/xiaocai/software/mpich-3.0.2/lib/ \
       -L/Users/xiaocai/software/parallel-netcdf-1.3.1/lib

INCDIR = -I/global/homes/w/wkliao/PnetCDF/include
LIBDIR = -L/global/homes/w/wkliao/PnetCDF/lib -lpnetcdf

LIB=-lpnetcdf -lmpich

all:graph_parallel 

graph_parallel:GraphBase.cpp BBox.cpp Grid.cpp Tree.cpp admm.cpp
	$(CC) $(INCDIR) $(CFLAGS) $^ -o $@ $(LIBDIR) $(LIB)

#Grid.o:Graphbase.cpp BBox.cpp Grid.cpp
#	$(CC) $(INCDIR) $(CFLAGS) $^ $(LIBDIR) $(LIB)

clean:
	rm -f graph_parallel *.o

