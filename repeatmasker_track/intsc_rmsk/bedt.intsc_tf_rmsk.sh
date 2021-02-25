for f in ../../peak/*.bed
do
	bn=`basename $f`
	rn=`echo $bn | cut -d "_" -f 1,2`
	bedtools intersect -a ../rmsk.bed -b ${f} -f 0.40 -wa  > ${rn}_rmsk.bed


done
