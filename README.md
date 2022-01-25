# input format
a connected undirect graph, with each edge appears only once in the file

# run
* a simple usage is given in main.c, change the parameter in the head of main.c before running
* make && cat n10000.inp | ./main_treem
* the output contains information for preprocessing and max-flow calculation

# how to generate input file
you can refer to igraph lib in python, for example, uniformly random graph (Erdos_Renyi), scale-free graph (Barabasi), unit disk graph (GRG).
