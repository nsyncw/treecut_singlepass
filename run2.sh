make
if [ "$?" != 0 ]
then
	echo "make error"
    exit 1 #参数正确，退出状态为0
fi
echo "start to run treem"
ulimit -s unlimited
ulimit -v unlimited
NNUM=1000000
TP="sf"
UTP="u$TP"
cat /mnt/E/mydir/myproj/research/graphLab/PCSP/lab_for_graph/data/n$NNUM$TP.inp | ./main_treem | tee tmp$NNUM$TP.txt
cat tmp$NNUM$TP.txt | grep "hi_treem_res" > out$NNUM$TP.txt
bash comphi.sh /mnt/E/mydir/myproj/research/graphLab/PCSP/lab_for_graph/data/n$NNUM$UTP.inp out$NNUM$TP.txt
