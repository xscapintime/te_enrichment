for f in ../../resize_peak/resized_peaks/*.bed
do
	bn=`basename $f`
	rn=`echo $bn | sed 's/.bed//g'`

	pct=0.70

	echo $bn "===================>\n""bedtools intersecting...\nfraction of overlap: " ${pct}
	bedtools intersect \
		-a ../del_qmarkline_rmsk.bed -b ${f} \
		-f ${pct} -wa  > ${rn}_rmsk.bed

	#echo "bedtools intersecting...\n fraction of overlap: " ${pct}
done
