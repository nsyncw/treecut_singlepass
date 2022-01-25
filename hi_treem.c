/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <math.h>

#include "types_treem.h"  /* type definitions */
#include "parser_treem.c" /* parser */
#include "timer.c"        /* timing routine */

#define MAX_LONG LONG_MAX

#define min(s, t) (s) < (t) ? (s) : (t)
#define max(s, t) (s) > (t) ? (s) : (t)







///////////////////////////////////////////
///////////////////////////////////////////The definition of functions
///////////////////////////////////////////


//The function for allocation
void *walloc(unsigned int num, unsigned int size)
{
  void *ptr = calloc(num, size);
  assert(ptr != NULL);
  return ptr;
}

//The function for randomization
cType mrand(RandomData* rd)
{
  return rd->randNums[rd->randNumIdx++];
}

RandomData* initrand(cType len)
{
  RandomData *rd = walloc(1, sizeof(RandomData));

  srand((int)(timer() * 1000));

  rd->randNums = (cType *)walloc(len, sizeof(cType));
  rd->randNumIdx = 0;

  for (int i = 0; i < len; i++)
  {
    rd->randNums[i] = (cType)rand();
  }

  return rd;
}

/////////////////////////The two function for heap sorting edges 
void HeapAdjustDown(sType *idx, edgeP * edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    // assert(idx[0] != idx[3]);
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]].tmp > edges[idx[i]].tmp )    
            i++;  

        if(edges[idx[i]].tmp <= edges[tempIdx].tmp )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void HeapSort(sType *idx, edgeP * edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        HeapAdjustDown(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        // printf("swap 0 with %d \n",i);
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        HeapAdjustDown(idx,edges,0,i-1);  
    }  

}  
  
/////////////////////The function to sort edges using capacity
void deOrderEdgeByRandomCap(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    // pedges[i].tmp = -1*mrand(pd->rd) % (pedges[i].cap+1); //+1防止0的情况
    pedges[i].tmp = 1000-pedges[i].cap+1; //+1防止0的情况
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    HeapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}

///////////////////The function to sort edges using the value of currently minimal cut one edge belongs to 
void aOrderEdgeByAvgCV(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    long acv = pp->avgCV;
    if(acv == 0 ){
      pp->tmp = MAX_LONG;
    }
    else{
      pedges[i].tmp = mrand(pd->rd) % acv; 
    }
  }

    HeapSort(idxs,pedges,cnt);
}

///////////////////按度数升序访问 
void aOrderEdgeByDegree(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;
  nodeP* nodes = pd->gd->nodes;
  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    cType zn = pp->endNode;
    pedges[i].tmp =  (nodes+zn)->nIdx;//mrand(pd->rd) % acv; 
  }

    HeapSort(idxs,pedges,cnt);
}

// //根据度数，但是同时
// void aOrderEdgeByDegreeAver(nodeP *np,PreprocData *pd)
// {
//   if (np->nIdx == 0)
//   {
//     return;
//   }
//   int cnt = np->nIdx;
//   nodeP* nodes = pd->gd->nodes;
//   sType *idxs = np->orderedEdges;
//   edgeP *pedges = np->edges;

//   for (int i = 0; i < cnt; i++)
//   {
//     edgeP *pp = pedges +i;
//     cType zn = pp->endNode;
//     pedges[i].tmp =  (nodes+zn)->nIdx;//mrand(pd->rd) % acv; 
//   }

//     HeapSort(idxs,pedges,cnt);
// }



/*
  the function to add more data to a traversal tree to accelerate the searching in the tree
  The idea is to precalcuate minimal cv value of a span of nodes in a traversal tree, e.g., when SPAN_LEN = 100, and a node has depth of 200, then the algorithm will pre-calculate the minimal cv value of the nodes between the node (dep=200) and an acestor(dep=101)
  upid is the id of the last SPAN node, mcv is the min cv among all previous nodes in the recent SPAN
  lastDepMCV is the depth of the node depth that has the minimal cv in the span
  lastJointNodeId is the last ancestor node id that has more than one child nodes
  lastJointMCV is the cv of lastJoineNodeId

  改版后预处理算法要有变化：
      对每个节点的值，都要改进下考虑当前出发节点z的考虑下游分支的更小的cv'=cv2+-cof

  (1)处理时：
      在非段头节点中还要考虑父亲cv'的值，如果更小，则更新指向父亲
      跨段的cv2+-cof直接放到下段mcv中，因为solve中确认跨段了才使用当前段的mcv
  (2)求解时：
      看solve的注释


*/

long gRoot = 0;
void buildAcc(PreprocData *pd, cType curN, cType upid, cType upid_LC, long mcv, long* applyAdj, cType* depMap)
{
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);
  int cnt = (nodes + curN)->nIdx;
  edgeP *pedges = (nodes + curN)->edges;

  long curCV = pd->gpcv[curN];
  cType curDep = pd->gpdep[curN];

  assert(curDep == 0 || curCV >0);

  depMap[curDep] = curN;
  pd->gpaccup[curN] = upid; //the previous segment tail

  int SPAN_LEN = pd->SPAN_LEN;
  int LN = pd->LEVEL_NUM;
  /******************************************/
  //作为祖先节点更新的点
  pd->gpaccup[curN] = upid; //the previous segment tail
  pd->gpaccposmcv[curN] = upid_LC; //the LN Child of previous segment tail
  
  if(curDep % SPAN_LEN == 0){
    upid = curN;
    upid_LC = 0; //节点id最小是1，所以0可以表示无效
  }
  else{
    //不需要更新什么，或者在后面逻辑中包含了

  }


  /****************************************/
  cType ances; 
  for (int i = 1; i <= pd->LEVEL_NUM; i++)
  {
    if (curDep  < (unsigned int)i)
    {
      break;
    }

    ances = depMap[curDep - i];
    assert(ances != 0);
    if (applyAdj[ances] == 0) //ances的applyAdj是0，说明未被祖先设置过
    {
      if (pd->gcutadjsign[curN][i] == 1 ) // 
      {
        // printf("update ApplyAdj level i %ld, curN %ld dep %ld, to set applyAdj[%ld] dep %ld to %ld\n",i,curN,pd->gpdep[curN], ances,pd->gpdep[ances],pd->gcutadj[curN][i]);
        applyAdj[ances] = pd->gcutadj[curN][i];
        pd->gcutadj[curN][i] = -1 * pd->gcutadj[curN][i];
      }
    }

  }

  //更新完，此时curN的LN祖先(最上面祖先)值，在这个支线上肯定不会再变了，可以尝试更新minCV了
  if(curDep >= LN ){
      if(curDep % SPAN_LEN == LN){
        upid_LC = curN; //更新作为祖先的upid_LC值
        if(LN > 0){
          //如果LN是0，这就是上个段尾，不适用此更新规则
          pd->gpaccposmcv[curN] = upid_LC; //把当前节点LC也设置为最新，后面好处理 
        }
      }
      else if(curDep % SPAN_LEN == (LN+1)){
        mcv = MAX_LONG;
      }

      assert(applyAdj[ances] <= 0);
      assert(pd->gcutadj[ances][0] <= 0);
      assert(pd->gpcv[ances]+pd->gcutadj[ances][0] > 0);
      //例行更新祖先值，这个和段无关,而且是循环
      ances = depMap[curDep-LN];
      // printf("applyAdj[ances] %ld\n",applyAdj[ances]);
      //curN的poh装的是LN祖先的在包含curN及其小祖先的基础上的树的最小割，因为applyAdj是负值或0
      //需要反向补充回去的目标值，应该不是祖先的pcv,而是祖先pcv还要小的，树被裁剪后的最小值，应此时cv+[ances][0]
      pd->gpoh[curN] = pd->gpcv[ances] + pd->gcutadj[ances][0] - applyAdj[ances]; 
      mcv = min(mcv, pd->gpoh[curN]);
      pd->gpaccmcv[curN] = mcv; //这个值记录的是LN祖先的该段的从段首到该祖先的所有割值的最小

  }
  else{
    //nothing to update

  }


  while (cnt > 0)
  {
    cType zn = pedges->endNode;
    if (pd->gpfa[zn] == curN)
    {
      buildAcc(pd,zn, upid, upid_LC, mcv, applyAdj, depMap);
    }
    pedges++;
    cnt--;
  }

  //恢复添加到祖先身上的值，恢复成0
  for (int i = 1; i <= pd->LEVEL_NUM; i++)
  {
    if (curDep  < (unsigned int)i)
    {
      break;
    }
    cType ances = depMap[curDep - i];
    // printf("curN %ld, dep %ld, ances %ld dep %ld\n",curN, pd->gpdep[curN], ances, pd->gpdep[ances]);
    // printf("level i %ld sign %ld, adj %ld  --  apladj %ld\n",i,pd->gcutadjsign[curN][i],pd->gcutadj[curN][i],applyAdj[ances] );
    assert(ances != 0);
    //这里是还结果的地方，之前赋值给applyAdj[ances]，这里还回去
    if (pd->gcutadjsign[curN][i] == 1 && pd->gcutadj[curN][i] > 0)
    {
      pd->gcutadj[curN][i] = -1 * pd->gcutadj[curN][i];
      assert(applyAdj[ances] == pd->gcutadj[curN][i]);
      applyAdj[ances] = 0;
    }
  }

}

void calcuTotalCap(PreprocData *pd)
{
  nodeP *nodes = pd->gd->nodes;

  for (cType curN = 1; curN <= pd->gd->N; curN++)
  {
    nodeP *np = nodes + curN;
    edgeP *pedges = np->edges;
    int cnt = np->nIdx;
    np->totalCap = 0;
    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + ni;
      np->totalCap += eh->cap;
    }
  }
}

/*

算法整体步骤：相对之前treem算法的改进
(1)预处理
	//按之前计算，加上一个值放到n上，代表对于n的father,如果去掉n这一支，f的mc变化
	sum_f = 0
	对f的每个为访问的child n:
		在遍历n前，在f上设置oh_f并置0
		在遍历计算中，每次访问祖先，除了计算mc，还更新oh_f (加上就可以)
		n返回后
			此时知道mc_n, oh_f(oh_f就是n这一支连到f的边的值，要包含f_n父子的边)
			这时计算n这一支去掉对f的mc的影响cof_n = oh_f - (mc_n-oh_f) = 2* oh_f - mc_n
				//oh_f是增加的割，mc_n-oh_f是之间到非f的割
			如果cof_n是负值 //说明割值降低了，值得庆贺
				就加到sum_f上,即sum_f 加上 负值，sum_f指的是f的所有减少割值的子n去掉，总共减少的割值
				如果不减少，这个n就不去掉
	//所有child遍历完之后计算mc2_f，就是f的最优的偏特树，此时这个树包含某些子树，而且mc是最小的
	得到mc2_f = mc_f + sum_f 
				
(2)计算时: solve 和 build的时候
	//因为mc2_f是已经去掉cof_n为负值的n的，正值就不管，还在mc_f中
	每次向上回溯，n回溯到f时，如果cof_n是负值，说明如果要保留n这一支，目前f最优割就包含n，即此时f用于计算的mc应该取(mc2_f + -cof_n)，即加上n这一支减掉的值
  现在问题来了：
    f的mc2和访问哪一支有关系，buildAcc预处理咋做？
    这样就意味着，不同的底层上来，每个节点的mc还不一样，导致 节点段 中最小值还不一样
    buildAcc记录的是向上的，所以可以记录的呀

*/
// long level_cumsum[110]; //保存子的每个level的负值的和
//////////////////////////////////The function to traverse the graph for one pass, i.e. the checkNode function of Algorithm 2 in the paper
void markCut(cType curN, PreprocData *pd)
{
  // printf("call markcut in %ld\n",curN);
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);

  // printf("******curN is %ld\n",curN);
  assert((nodes + curN)->nIdx > 0);

  short *curS = pd->gps + curN;

  assert(*curS == 0);

  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedEdges == NULL)
  {
    np->orderedEdges = (sType *)walloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedEdges[i] = i;
    }
  }

  long cap;
  sType *idxs = np->orderedEdges;

  if (pd->mode == 1)
  {
    deOrderEdgeByRandomCap(np,pd); //
  }
  else if (pd->mode == 2)
  {
    aOrderEdgeByAvgCV(np,pd);
  }
  else if (pd->mode == 3){
    //安装度数最小最先访问
    aOrderEdgeByDegree(np,pd);
  }

  // long sum_f = 0;

  //多层遍历算法步骤：每次遍历前，记录上面几位祖先的值到cut_adjust中
  //问题：还没有遍历成功
  cType fa = curN;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    if(pd->gpdep[fa] == 0){
      break;
    }
    // printf("fa is %ld, fasfa is %ld,  dep is %ld\n",fa,pd->gpfa[fa],pd->gpdep[fa]);
    fa = pd->gpfa[fa];
    // printf("   -> fa is %ld, dep is %ld\n",fa,pd->gpdep[fa]);
    
    pd->gcutadj[curN][i] = pd->gpoh[fa]; 
  }


  for (int ni = 0; ni < cnt; ni++)
  {
    // nodeP* znp = nodes+eh->endNode;

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;

    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn];

    assert(!(pd->gpfa[zn] == curN && zs == 2));

    // printf("zn is %ld (curN %ld) \n",zn,curN);
    if (zs == 1) // zn is an ancestor node or father node
    {
      // printf("\t zs is 1, zn %ld , *cruCV %ld += %ld\n",zn,*curCV,cap);
      cap = eh->cap;
      *curCV += cap;
      pd->gpcv[zn] -= cap;
      pd->gpoh[zn] += cap;
    }
    else if (zs == 0) // zn is not accessed, i.e., a child node
    {
      //poh临时使用，后面会作为cv2使用：ph记录curN的当前儿子zn子树遍历中，连到curN的边的容量的和，就是直接和curN连接的割值。所以需要置0
      // pd->gpoh[curN] = 0;

      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;

      // printf("\t zs is 0, zn %ld , *cruCV %ld before markcut\n",zn, *curCV);
      markCut(zn,pd);

      //assert(pd->gpoh[curN] > 0);
      // printf("----marCut return \n");
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      // printf("\t zs is %ld , zn %ld , *cruCV %ld += %ld after markcut\n",pd->gps[zn],zn, *curCV,pd->gpcv[zn]);
      *curCV += pd->gpcv[zn];

//       //此时：gpoh[curN]是zn子树连到curN的边总容量，即zn树所有节点脱离curN树后，curN树多出来的割
//         //zn树不连其他curN的子的树，只连curN或上面祖先
//       // pd->gpcv[zn] - pd->gpoh[curN]是zn子树在curN之上的割值, 即zn子树去掉后，curN减少的割值
		 // pd->gpoh[curN] 就是zn子树去掉后，增加的割值
         //cof_zn 就是zn树去掉后，整体增加的割值，如果小于0，就可以去掉
      // long cof_zn = pd->gpoh[curN] - ( pd->gpcv[zn] - pd->gpoh[curN]);
      // if(cof_zn < 0){ //如果去掉能进一步优化割值(减少)，则记录
      //   sum_f += cof_zn;
      // }
      // assert(pd->gpcv[zn] >= pd->gpoh[curN]);
      // pd->gpcof[zn] = cof_zn;
      // pd->gpcof[zn] = pd->gpoh[curN];

    }
    else
    {
      //这种情况就是,zn是curN的子树的叶子，正好连到curN有条边
      // printf("zs is %ld\n",zs);
      assert(pd->gpdep[curN] < pd->gpdep[zn]);
      
    }

  }

  // if (pd->gpoh[curN] == 0)
  // {
  //   printf("poh error: (cnt is %d ) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n", cnt, curN, gRoot, *curCV, pd->gpdep[curN]);
  //   assert(1 == 2);
  // }
  //自己有没算进去？遍历所有邻居就是把自己算进去了
  //------------多层遍历算法步骤：每次遍历后，更新祖先增加的值到cut_adjust中
  fa = curN;
  cType faBelow = 0; //连fa到curN之前的容量值
  int actual_ln = -1;
  for(int i=1; i<=pd->LEVEL_NUM; i++ ){
    // printf("fa is %ld, dep is %ld\n",fa,pd->gpdep[fa]);
  
    if(pd->gpdep[fa] == 0){
      actual_ln = i-1;
      break;
    }

    fa = pd->gpfa[fa];
    //此时gcutadj保留的是curN中到fa这个祖先节点的边容量总和
    pd->gcutadj[curN][i] = pd->gpoh[fa] - pd->gcutadj[curN][i];
    //(1) 下面计算后的faBelow是curN子树所有cut到fa内部(包含fa)到curN前的(不包含curN)的边容量总和
    faBelow += pd->gcutadj[curN][i];
    //(2) 而curN到fa之外(不包含fa)的cut值为  *curV - faBelow;
    //则：(1)-(2)的值是curN的子树全部去除后，对相应祖先fa的cut的增加值，负值最好
    pd->gcutadj[curN][i] = faBelow - (*curCV - faBelow);    
  
  }

  //actual_ln说明curN最大可以作为哪一层？ 0就是没有(curN就是root)
  if(actual_ln < 0){
    actual_ln = pd->LEVEL_NUM;
  }

  //如果level_num=5但是actual_ln=2,说明当前节点最多只能作为第2层，而且就是底层了
  //作为底层

  //作为最底层
  // if(pd->gcutadj[curN][actual_ln] < 0){
  //   pd->gcutadjsign[curN][actual_ln] = 1;
  // }

  //----------多层遍历算法步骤：


  // assert(actual_ln < 20);//给个限制，最多20层
  // memset(level_cumsum,0,(pd->LEVEL_NUM+1)*sizeof(long));
  //现在比较curN的i层和子的i+1层之间的值，两层的这两个值对应的祖先是同一个
  long tempCumSum = 0;
  pd->gcutadj[curN][0] = 0; //0层记录的是自己对自己的cv值的改变，或者说整个树中通过裁剪能达到的最大的可减少负值(保留自己为根的前提下)，该负值针对curN为根的树
  for (int i = (actual_ln < pd->LEVEL_NUM ? actual_ln : actual_ln - 1); i >= 0; i--) //只能从curN的level-1开始，因为子是要level开始；i可以为0，因为对应子时i+1层
  {
    //子对同个fa的原始cv的负值的总和，存在这里
    tempCumSum = 0;

    for (int ni = 0; ni < cnt; ni++)
    {
      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      assert(zn != 0);
      assert(zn != curN);

      if (pd->gpfa[zn] != curN)
      {
        continue;
      }

      if (pd->gcutadj[zn][i + 1] < 0)
      {
        tempCumSum += pd->gcutadj[zn][i + 1];
      }
    }

      // printf("i %ld, level_cumsum[i + 1] %ld\n",i,tempCumSum);

    if (pd->gcutadj[curN][i] < tempCumSum) //tempCumSum是全部子的负值的总和(对祖先可减少的总负值)
    {
      // if(curN == 9904 && i == 4){
      //   printf("curN 9904 i=4 set sign self 1");
      // }
      pd->gcutadjsign[curN][i] = 1; //说明curN子树中，全树全部移除反而减少的更多
      // pd->gcutadjsign[curN][i+1] = 0; //设置的不是curN，应该是zn
    }
    else
    {


      // if(curN == 9904 && i == 4){
      //   printf("curN 9904 i=4 set sign children 1");
      // }
      pd->gcutadj[curN][i] = tempCumSum;
      pd->gcutadjsign[curN][i] = 0; //虽然用了子孙中所有负值，但是设0表示没有用自己全树参与
      // pd->gcutadjsign[curN][i+1] = 1;
    }

  }


  // if(pd->gpoh[curN] == 0 && curN != gRoot){ //当邻居全是祖先是有可能的，已经在上面加了判断
  //   printf("poh error: (cnt is %d ) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
  //   assert(1==2);
  // }
  if(*curCV == 0 && curN != gRoot){
    printf("cv=0 error: (cnt is %d) curN is %ld (gRoot is %ld), cv is %ld, dep is %ld \n",cnt,curN, gRoot,*curCV,pd->gpdep[curN]);
    assert(1==2);
  }

  // assert(pd->gpoh[curN] >0);

  if(pd->mode == 1){
//update ver 2 w,according to
    for (int ni = 0; ni < cnt; ni++)
    {
      // nodeP* znp = nodes+eh->endNode;

      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      // nodeP *znp = nodes+zn;

      assert(zn != 0);
      assert(zn != curN);
      short zs = pd->gps[zn];
    

      //progate weight to curN's edges
      if (zs == 1 && pd->gpdep[zn] != *curDep - 1)
      {
          cType weight = eh-> w;
          if(eh->avgCV == 0){
            eh->avgCV = MAX_LONG;
          }
          eh->avgCV = min(eh->avgCV, *curCV);//((eh->avgCV) * weight + *curCV)/(weight+1);
          eh->w = weight+1;
          
          edgeP *reh = eh->rev;

          if(reh->avgCV == 0){
            reh->avgCV = MAX_LONG;
          }

          weight = reh-> w;
          reh->avgCV = min(reh->avgCV, *curCV);//((reh->avgCV) * weight + *curCV)/(weight+1);
          reh->w = weight+1;

      }
      
    }
  }

    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  // printf("set curN %ld to stat 2 \n",curN);
  *curS = 2;

}

////////////////////////////////The function to obtain min-cut value of given node pair, i.e., Algorithm 3 in the paper
/*
 (2)求解时：
      深度大的出发节点直接cv2，
      另外一个不是这个的祖先时，另外一个也可以用cv2
      [这个先不优化有点复杂]如果t是s祖先，
        如果cv2正好割掉下面，直接可以用cv2
        如果不是，也可以计算割掉对应树后这个t的割

      PS: t是s祖先，可以用cv2,t不是s的祖先也可以cv2，两个cv2都可以用
      问题是，如果t是s祖先，t的cv2对应的割并不能切断s这条支线，这样就不对了
        如果t是s祖先，其实我们算的是去掉这个支后的t的最小割
      那就简化：
        只要t不是s祖先，就可以用t的cv2

      后面的逐个逼近，需要用cv'
      
*/
//对给定节点，更新祖先值，并记录到按深度差为key的数组中，数组也就LN长度
void traceUp2LN(NodePropArr* pnp, cType startDep, cType curN, cType curDep, long *adj, cType LN)
{

  for (int i = 1; i <= LN; i++)
  {
    if ( curDep < (unsigned int)i || (curDep+LN) <= (startDep+i) )
    {
      break; //祖先已经不存在了，或已经到startDep的LN祖先了，没必要计算了
    }

    if (pnp->pacc_cut_adjust_sign[curN][i] == 1)
    {
      adj[startDep+i-(curDep)] = pnp->pacc_cut_adjust[curN][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
    }
  }
}

//函数只要检查，LN_v的最优值，是否包含v; 如果不包含，区分是v_top(不含)以下点最优不包含，还是v_top(含)以上
//
//v到LN_v的情况，每个点
long isOptimalExclude(NodePropArr* pnp, cType v, cType v_top, cType LN_v, cType *pDep, cType *pFa, cType LN){
  
  cType minDep = 0;
  // printf("isOptimalExclude: v %ld dep %ld, LN_v %ld dep %ld\n",v,pDep[v],LN_v,pDep[LN_v]);
  while(v != LN_v){
    //如果有节点对应LN_v标志是1，说明肯定排除在最优之外，v的祖先有没有，v都被排除了。
    //如果v没有但祖先有，也可以通过再Fa回溯检测到
    if(pDep[v] - pDep[LN_v] <= LN && pnp->pacc_cut_adjust_sign[v][pDep[v] - pDep[LN_v]]==1){
      //如果不包含，还要体现到哪不包含
      minDep = pDep[v];    
    }
    v = pFa[v];
  }

  // printf("isOptimalExclude DONE \n");

  //除了v到v_top(不含)被包含或不包含的情况，还有一种是根本不知道，因为不在这个深度裁剪，这个时候什么情况？
  //这种咋处理?
  //从共同LN_v开始每个祖先，有几种可能：
  //(1)LN_v的最优包含一方，不包含另外一方
  //具体情况: 
  //(1) LN_s(不含)之下的祖先，单独看就行，其实就只有LN_s，其他也没意义
  /*(2) LN_s(含)及之上的祖先，则需要看不含一方的最高点，如果不含的最高点在LN_s(不含)以下，可以认为不含，因为去掉t不影响s。
        如果在之上就麻烦了，去掉t也去掉s了，这个就不好办了

  */
  //默认LN_v最优不需要去掉v
  //如果包含，那就是都包含，这个不用继续看了
  if(minDep > 0){
    if(minDep <= pDep[v_top]){
      return 2; //在v_top及之上断开了
    }
    else{
      return 1;
    }
  }

  return 0;
}

//其实还可以考察共同点v上面的祖先点W，只要v到w都没被丢弃，则s和t单独被丢弃的情况，都可以作为最优解来处理
//因为w不确定，其实是找最低的w, w的最小值
long solveMaxFlowAccVER4(PreprocData *pd, cType root, long minCandi, NodePropArr* pnp, cType s, cType t, int SPAN_LEN, aType LN)
{
  cType *pDep = pnp->pdep;
  long *pCV = pnp->pcv;
  cType *pFa = pnp->pfa;
  // long *poh = pnp->poh;
  // cType *paccup = pnp->pacc_upid;
  // cType *paccup_LC = pnp->pacc_pos_upmincv;
  // long *paccmcv = pnp->pacc_upmincv;

  // cType os = s, ot = t;
  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  printf("\nstart: s %ld dep %ld, t %ld dep %ld\n",s,pDep[s],t,pDep[t]);

  assert(pDep[s] >= pDep[t]);

  // long *adj = (long *)walloc(pd->gd->N+1, sizeof(long));
  long *adjs = (long *)walloc(pd->gd->N+1, sizeof(long));
  long *adjt = (long *)walloc(pd->gd->N+1, sizeof(long));
  // memset(adj, 0, (pd->gd->N+1) *  sizeof(long));
  memset(adjs, 0, (pd->gd->N+1) *  sizeof(long));
  memset(adjt, 0, (pd->gd->N+1) *  sizeof(long));

  long mcv = MAX_LONG;

  long startDeps = pDep[s];
  long startDept = pDep[t];

  while(pDep[s] > pDep[t]){

    //先更新对祖先的补回值
    for (int i = 1; i <= LN; i++)
    {
      if (pnp->pacc_cut_adjust_sign[s][i] == 1)
      {
        adjs[ startDeps - (pDep[s]-i) ] = pnp->pacc_cut_adjust[s][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
      }
    }

    //调用祖先优化
    mcv = min(mcv, pCV[s] + pnp->pacc_cut_adjust[s][0] - adjs[startDeps - pDep[s]]); 
    s = pFa[s];

  }

  assert(pDep[s] == pDep[t]);

  if(s == t){
    goto end;
  }

  //此时s和t不相等，一起上溯
  while(s != t){
    for (int i = 1; i <= LN; i++)
    {
      if (pnp->pacc_cut_adjust_sign[s][i] == 1)
      {
        adjs[ startDeps - (pDep[s]-i) ] = pnp->pacc_cut_adjust[s][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
      }

      if (pnp->pacc_cut_adjust_sign[t][i] == 1)
      {
        adjs[ startDept - (pDep[t]-i) ] = pnp->pacc_cut_adjust[t][i]; //这个是curN整个树去掉相对i祖先的cv的影响负值
      }


    }

    //调用祖先优化
    mcv = min(mcv, pCV[s] + pnp->pacc_cut_adjust[s][0] - adjs[startDeps - pDep[s]]); 
    mcv = min(mcv, pCV[t] + pnp->pacc_cut_adjust[t][0] - adjs[startDept - pDep[t]]); 

    s = pFa[s];
    t = pFa[t];
  }

end:
  free(adjs);
  free(adjt);
  return mcv;


}

// long solveMaxFlowAccVER5(cType root, long minCandi, NodePropArr* pnp, cType s, cType t, int SPAN_LEN, aType LN)
// {
//   cType *pDep = pnp->pdep;
//   long *pCV = pnp->pcv;
//   cType *pFa = pnp->pfa;
//   long *poh = pnp->poh;
//   cType *paccup = pnp->pacc_upid;
//   // cType *paccup_LC = pnp->pacc_pos_upmincv;
//   // long *paccmcv = pnp->pacc_upmincv;

//   assert(s != t);
//   if (pDep[s] < pDep[t])
//   {
//     cType tmp = s;
//     s = t;
//     t = tmp;
//   }

//   assert(pDep[s] >= pDep[t]);

//   long *adj = (long *)walloc(LN+1, sizeof(long));
//   long *adj2 = (long *)walloc(LN+1, sizeof(long));
//   memset(adj, 0, (LN+1) *  sizeof(long));
//   memset(adj2, 0, (LN+1) *  sizeof(long));

//   long mcv = MAX_LONG;
//   cType LN_s = s , LN_t = t ;
//   int alsoT = 0;

//   if(pDep[LN_s] == pDep[LN_t]){
//     alsoT = 1;
//   }

//   assert(pDep[LN_t] <= pDep[LN_s]);

//   //先找到s和t的LN祖先，如果s到LN包含t，直接就可以退出了
//   while (pDep[s] - pDep[LN_s] < LN) //差整个LN的不用纳入循环计算
//   {

//     // printf("check: (LN_s is %ld, dep is %ld ) LN_t is %ld, dep is %ld  \n", LN_s, pDep[LN_s], LN_t, pDep[LN_t]);
//     assert(pDep[LN_t] <= pDep[LN_s]);

//     if (alsoT > 0)
//     {
//       assert(pDep[LN_t] == pDep[LN_s]);
//       //这里mcv到底用什么值？按理说应该用对应节点的cv，去掉LN子减掉的值；还是最优树，加上子去掉的值
//       mcv = min(mcv, pCV[LN_t] + pnp->pacc_cut_adjust[LN_t][0] - adj2[pDep[t] - pDep[LN_t]]); 
//       traceUp2LN(pnp, pDep[t], LN_t, pDep[LN_t], adj2, LN);
//       LN_t = pFa[LN_t];
//     }

//     mcv = min(mcv, pCV[LN_s] + pnp->pacc_cut_adjust[LN_s][0] - adj[pDep[s] - pDep[LN_s]]);
//     traceUp2LN(pnp, pDep[s], LN_s, pDep[LN_s], adj, LN);
//     LN_s = pFa[LN_s];

//     if (pDep[LN_s] == pDep[LN_t])
//     {
//       if (LN_s == LN_t)
//       {
//         //遇到t了，说明t就是s的祖先，而在LN之内，则：直接返回
//         goto end;
//       }

//       //到这里说明不是t,但是两者深度一样了，此时需要将t也纳入一起上溯的过程
//       alsoT = 1;
//     }

//     // if(pDep[LN_s] == 0){
//     //   //这种情况不可能存在！最多t和s汇总到root，在上面就退出了
//     //   LN_s = 0; //0说明还没到祖先就退出了
//     //   break;
//     // }
//   }
//     //这里，如果LN_s是深度为0的情况，说明s还没算到LN就到root了，那就不可能存在这种情况
//     //所以LN_s在这里必须深度不是0，而且>0


//     assert(pDep[LN_s] > 0);
//     assert(pDep[s] - pDep[LN_s] == LN);
//     assert(pDep[LN_t] > 0 || LN_t == root);
//     assert(pDep[LN_t] <= pDep[LN_s]);
//     assert(pDep[LN_t] == pDep[LN_s] || LN_t == t); //要么LN_t没动和t一样，要么深度和LN_s一样

//   //虽然LN_s是正常的LN外的祖先，但不排除，LN_s还没到正常的LN_t，就已经和t分支融合了
//   //此时，需要将LN_s赶到和LN_t平齐，再一直上升到LN_t(可能没到LN_t就结束了)
//   //用LN_t是因为t分支可能也被更新过

//   // printf("1mcv: %ld\n",mcv);
//   while (pDep[LN_s] > pDep[LN_t])
//   {

//     // printf("check3: (LN_s is %ld, dep is %ld, pup is %ld ) s is %ld, dep is %ld.  LN_t is %ld, dep is %ld (root is %ld)  \n",LN_s,pDep[LN_s], paccup[LN_s], s, pDep[s], LN_t, pDep[LN_t],root);
//     assert(pDep[LN_s] > 0 || LN_s == root);
//     assert(LN_s > 0);

//     mcv = min(mcv, poh[s]);
//     s = pFa[s];
//     LN_s = pFa[LN_s];

//   // printf("check4: (LN_s is %ld, dep is %ld ) s is %ld, dep is %ld ,pacc[LN_s] is %ld, dep is %ld  \n",LN_s,pDep[LN_s], s, pDep[s],paccup_LC[LN_s],pDep[paccup_LC[LN_s]]);
//     assert(pDep[s] - pDep[LN_s] == LN);

//   }

//   if (LN_s == LN_t)
//   {
//     goto end;
//   }

//   // printf("2mcv: %ld\n",mcv);
//   assert(pDep[LN_t] == pDep[LN_s]);
//   assert(pDep[paccup[LN_s]] < pDep[LN_t]);

//   //此时LN_s所在段的上个段尾已经在t之上了，即t肯定和LN_s在同段了，而且在t下面
//   // printf("check2: (LN_s is %ld, dep is %ld ) LN_t is %ld, dep is %ld  \n",LN_s,pDep[LN_s], LN_t, pDep[LN_t]);

//   //此时和t一起

//   // printf("3mcv: %ld\n",mcv);
//   //此时如果LN_t还不是最终LN_t，即没有达到和t距离是LN(因为前面可能s上溯到LN_s时同时带动了t分支)
//   //需要同时LN_s一起上溯
//   while(pDep[t] - pDep[LN_t] < LN)
//   {
//       assert(pDep[LN_t] == pDep[LN_s]);
//       if (LN_s == LN_t)
//       {
//         goto end;
//       }
//       mcv = min(mcv, pCV[LN_t] + pnp->pacc_cut_adjust[LN_t][0] - adj2[pDep[t] - pDep[LN_t]]);
//       traceUp2LN(pnp, pDep[t], LN_t, pDep[LN_t], adj2, LN);
//       LN_t = pFa[LN_t];      

//       mcv = min(mcv, poh[s]);
//       s = pFa[s];
//       LN_s = pFa[LN_s];
//   }

//   if (LN_s == LN_t)
//   {
//     goto end;
//   }

//   assert(pDep[LN_t] == pDep[LN_s]);
//   assert(pDep[t] - pDep[LN_t] == LN);
//   // printf("check5: (LN_s is %ld, dep is %ld ) s is %ld, dep is %ld  \n",LN_s,pDep[LN_s], s, pDep[s]);
//   assert(pDep[s] - pDep[LN_s] == LN);

//   // printf("5mcv: %ld\n",mcv);
//   //此时LN_s和LN_t在一个段里了
//   while(LN_s != LN_t){
//     mcv = min(mcv, poh[s]);
//     s = pFa[s];
//     LN_s = pFa[LN_s];

//     mcv = min(mcv, poh[t]);
//     t = pFa[t];
//     LN_t = pFa[LN_t];        
//   }

// end:
//   //VER55555555555555555555555555555555
//   // printf("check: (LN_s is %ld, dep is %ld ) LN_t is %ld, dep is %ld  \n",LN_s,pDep[LN_s], LN_t, pDep[LN_t]);
//   assert(LN_s == LN_t);
//   //最优是最优了，关键是比如LN_s取最优一定不包含s吗？怎么判断？
//   //只要从s到LN_s中有一个对应LN_s的cutsign是1就可以
//   // printf("mcv : %ld, final root: (LN_s/LN_t is %ld, dep is %ld ) \n",mcv,LN_s,pDep[LN_s]);
//   free(adj);
//   free(adj2);
//   return mcv;


// }
//////////////////////function to load graph data
void loadGraphData(PreprocData *pd){
  pd->gd = walloc(1,sizeof(GraphData));  
  printf("c\nc hi_treem version 1.1\n");
  printf("c Copyright C by nsyncw, nsyncw@gmail.com\nc\n");

  parse(&(pd->gd->N), &(pd->gd->M), &(pd->gd->nodes));

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", pd->gd->N, pd->gd->M);
}

///////////////////function to initialize data structure
void initPreprocData(PreprocData *pd){
  pd->rd = NULL;
  pd->SPAN_LEN = (int)(sqrt(pd->gd->N));

  pd->roots = (cType *)walloc(pd->total + 2, sizeof(cType));
  pd->allResults = walloc(pd->total+2, sizeof(NodePropArr));

  NodePropArr * allResults = pd->allResults;
  cType len = pd->gd->N + 2;

  calcuTotalCap(pd);

  int LN = pd->LEVEL_NUM+1;
  for (int i = 0; i < pd->total; i++)
  {
    allResults[i].pfa = (cType *)walloc(len, sizeof(cType));
    allResults[i].pdep = (cType *)walloc(len, sizeof(cType));
    allResults[i].pcv = (long *)walloc(len, sizeof(long));
    allResults[i].poh = (long *)walloc(len, sizeof(long));
    allResults[i].pcof = (long *)walloc(len, sizeof(long));
    allResults[i].ps = (short *)walloc(len, sizeof(short));
    allResults[i].pacc_upid = (cType *)walloc(len, sizeof(cType));
    allResults[i].pacc_upmincv = (long *)walloc(len, sizeof(long));
    allResults[i].pacc_pos_upmincv = (cType *)walloc(len, sizeof(cType));

    long* ptr = (long *)walloc(len*LN, sizeof(long));
    memset(ptr, 0, len * LN * sizeof(long));
    allResults[i].pacc_cut_adjust = (long **)walloc(len, sizeof(long*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust[j] = ptr+j*LN;
    }


    aType* ptr2 = (aType *)walloc(len*LN, sizeof(aType));
    memset(ptr2, 0, len * LN * sizeof(aType));
    allResults[i].pacc_cut_adjust_sign = (aType **)walloc(len, sizeof(aType*));
    for(int j = 0; j<len; j++){
      allResults[i].pacc_cut_adjust_sign[j] = ptr2+j*LN;
    }    

    memset(allResults[i].pfa, 0, len * sizeof(cType));
    memset(allResults[i].pdep, 0, len * sizeof(cType));
    memset(allResults[i].pcv, 0, len * sizeof(long));
    memset(allResults[i].poh, 0, len * sizeof(long));
    memset(allResults[i].pcof, 0, len * sizeof(long));
    memset(allResults[i].ps, 0, len * sizeof(short));
    memset(allResults[i].pacc_upid, 0, len * sizeof(cType));
    memset(allResults[i].pacc_upmincv, 0, len * sizeof(long));
    // memset(allResults[i].pacc_cut_adjust, 0, len * sizeof(cType*)); //存的指针不能初始化，指针再指向区域已经初始化了
    // memset(allResults[i].pacc_cut_adjust_sign, 0, len * sizeof(aType*));

  }  
}


/////////////////////function to traverse the graph data for multiple times, i.e., Algorithm 2 in the paper
void preProc(PreprocData *pd){
  double tm;
  double totalProcTime = 0;
  NodePropArr *allResults = pd->allResults;
  //calculate total cap of one node
  cType root;

  cType len = pd->gd->N + 2;
  long *apply_adj = (long *)walloc(len, sizeof(long));
  cType *depth_map = (cType *)walloc(len, sizeof(cType));

  for (int ipass = 0; ipass < pd->total; ipass++)
  {
    if(pd->rd != NULL){
      free(pd->rd);
      pd->rd = NULL;
    }
    pd->rd = initrand(pd->gd->M*2);    
    // printf("the %d times\n",i);
    pd->gpfa = allResults[ipass].pfa;
    pd->gpdep = allResults[ipass].pdep;
    pd->gpcv = allResults[ipass].pcv;
    pd->gpoh = allResults[ipass].poh;
    pd->gpcof = allResults[ipass].pcof;
    pd->gps = allResults[ipass].ps;
    pd->gpaccup = allResults[ipass].pacc_upid;
    pd->gpaccmcv = allResults[ipass].pacc_upmincv;
    pd->gpaccposmcv = allResults[ipass].pacc_pos_upmincv;
    pd->gcutadj = allResults[ipass].pacc_cut_adjust;
    pd->gcutadjsign = allResults[ipass].pacc_cut_adjust_sign;

    if (pd->P == 300)
    {
      pd->mode = 3;
    }
    else
    {
      pd->mode = ipass < pd->P * pd->total / 100 ? 1 : 2;
    }
    //----------------准备深度遍历   
    // if(pd->mode == 3){
    //   //此时预期是sf图，而根据我们生成sf图的方式，1是度数最大的节点
    //   root = 90000;
    // }
    // else{
    root = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % pd->gd->N);
    // }
    pd->roots[ipass] = root;
    pd->gpdep[root] = 0;
    printf("pass %d before markCut: root is %ld\n",ipass,root);
    fflush(stdout);
    //printf("pass %d, randidx %ld, root is %ld\n",i, randNumIdx,root);
    // printf("root fa %ld\n",allResults[i].pfa[root]);
    tm = timer();
    gRoot = root;
    
    markCut(root,pd);
    pd->gpcv[root] = MAX_LONG;
    pd->gpoh[root] = MAX_LONG;


    //----------------准备构建加速结构
    /*
    printf("c before buildAcc\n");
	  fflush(stdout);

    memset(apply_adj, 0, len *  sizeof(long));
    memset(depth_map, 0, len *  sizeof(cType));
    buildAcc(pd, root, root, 0,MAX_LONG, apply_adj,depth_map);
    printf("c after buildAcc\n");
    fflush(stdout);
    */

    totalProcTime += timer() - tm;
    printf("c proctime for onepass: %10.06f\n", timer() - tm);
    if (ipass % 10 == 0)
    {
      printf("c the %d passes\n", ipass);
    }
  }

  free(apply_adj);
  free(depth_map);

  printf("c preprocess times %10.6f\n", totalProcTime);

}


//////////////////////////function to calculate multiple random node pairs, i.e., the calling of Algoirthm 3 for multiple times
void calcuRandomPairs(int numOfPairs, PreprocData *pd){
  double totalTime = 0;
  long mv = MAX_LONG;

  double curTime = 0;
  cType ns, nt;

  if(pd->rd != NULL){
    free(pd->rd);
    pd->rd = NULL;
  }
  pd->rd = initrand(pd->gd->M*2);  

  for (int ipair = 0; ipair < numOfPairs;)
  {

    // printf("%d\n",i);
    ns = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
    nt = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
    if (ns != nt)
    {
      // mv = MAX_LONG;
      mv = min((pd->gd->nodes+ns)->totalCap, (pd->gd->nodes+nt)->totalCap);
      
      curTime = timer();
      for (int j = 0; j < pd->total; j++)
      {
        cType root = pd->roots[j];
        // memset(apply_adj, 0, len *  sizeof(long));
        long tmp = solveMaxFlowAccVER4(pd, root,mv, pd->allResults+j, ns, nt,pd->SPAN_LEN,pd->LEVEL_NUM);
        // printf("--solve return %ld\n",tmp);
        if (mv > tmp)
        {
          mv = tmp;
        }
        // printf("--mv is %ld\n",mv);
      }
      curTime = timer() - curTime;
      totalTime += curTime;
      ipair++;
      printf("c hi_treem_res(n,s,mflow,tm) %lu %lu %12.01f %12.06f1\n", ns, nt, 1.0 * mv, curTime);
    }
  }

  printf("c run ok! average time %10.6f\n", totalTime / numOfPairs);


}


