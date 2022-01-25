#!/bin/bash
CDIR=`pwd`
cd /mnt/disk2/progs/graphLab/maxflow/hipr35
CNT=0
CNT2=0
INP=$1
while read nop comm ns nt mv tm
do
	echo "ns nt is $ns $nt"
    bash $CDIR/solveGeneral.sh $INP $ns $nt | ./hi_pr | grep "c flow:" > tmpcmp.txt
	cat tmpcmp.txt
    while read comm1 comm2 mv2 tm2
    do
        VAL=0
        if [ "$mv" == "$mv2" ] 
        then
            VAL=1
            let CNT2+=1
        fi
        let CNT+=1   

        echo $CNT2/$CNT $VAL `echo "scale=2; $tm2/$tm" | bc` mflow $mv $mv2, tm $tm $tm2 

    done <tmpcmp.txt
done <$CDIR/$2

