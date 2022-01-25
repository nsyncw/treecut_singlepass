/* defs.h */

//#ifdef EXCESS_TYPE_LONG
///typedef unsigned long excessType;
//#else
typedef unsigned long long int excessType; /* change to double if not supported */
//#endif

typedef unsigned long cType;
typedef unsigned int sType;
typedef unsigned char aType;

typedef 
   struct edgeProp
{

   cType endNode;
   cType cap;
   cType w; /*weight*/
   cType avgCV; /*average CV of cut set where this edge belongs to */
   long tmp;

   struct edgeProp* rev; //reverse

}edgeP;


typedef  /* node */
   struct nodeProp
{
   edgeP* edges;
   cType maxEdges;
   cType nIdx;
   cType totalCap;
   
   sType* orderedEdges;

} nodeP;

typedef
struct NodePropExtra_
{
   cType fa;
   cType dep;
   long cv;
   short s;


} NodePropExtra;


//这个是记录所有节点的一个数据结构，这一个记录所有节点
typedef
struct NodePropArr_{

   cType * pfa;
   cType * pdep;
   long * pcv;
   long * poh; //when traversing: record the edges' capacity connected to f by tree of a child n
               //when traversed: record the new updated cv
   long * pcof; //record the remove of n, bring to the change of father' cv
   short * ps;
   cType * pacc_upid; // up node id
   long* pacc_upmincv; // min of [currnet node, upnode)
   cType * pacc_pos_upmincv; // the diff of depth of nearest min value
   long ** pacc_cut_adjust; // for multiple level acc
   aType ** pacc_cut_adjust_sign; // for multiple level acc
   // cType * pacc_jointid; //the nearest joint id;
   // long* pacc_jointmcv; //the min of [cur, joint node)
   
} NodePropArr;

typedef
struct GraphData_{

  long N, M;
  nodeP *nodes;

} GraphData;

/////////////////////////////////////////
///////////////////////////////////////// The definition of data structure
/////////////////////////////////////////

//data structure that holds graph data

//data structure that holds data for randomrization
typedef
struct RandomData_{
  cType *randNums;
  cType randNumIdx;
} RandomData;


//data structure that holds preprocessing data
typedef
struct PreprocData_{

  //graph data
  GraphData *gd;

  //holds data in all passes
  NodePropArr* allResults;

  //pre-generated data
  RandomData* rd;

  //<BEING> hold the hot data in current pass
  cType *gpfa;
  cType *gpdep;
  long *gpcv;
  long *gpoh;
  long *gpcof;
  short *gps;
  cType *gpaccup;
  long *gpaccmcv;
  cType *gpaccposmcv;
  long **gcutadj;
  aType **gcutadjsign;

  //<END>

  cType *roots; //records root in each pass

  int mode; //the traversing mode, as explained in the paper

  int P; // P% percent of total passes in mode 1, the remaining in mode 2
  int total; //number of total passes
  int SPAN_LEN; //the length of a ancestor span, used for acceleration in traversal trees
  int LEVEL_NUM;

} PreprocData;

