#include "hi_treem.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    //no mode 2, so it is just randomly picking up children
    pd->P = 100; // P% percent of total passes in mode 1, the remaining in mode 2
    pd->total = 1; //number of total passes
    pd->LEVEL_NUM = 5; 

    loadGraphData(pd); //load graph data from standard input

    initPreprocData(pd); //init data structure

    preProc(pd); // preproc by traversing the graph for pd->total times

    calcuRandomPairs(200,pd); // randomly choose 100 node pairs and calcu their min-cut and output
  
    exit(0);

}
