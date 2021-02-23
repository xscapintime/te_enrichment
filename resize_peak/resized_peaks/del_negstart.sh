
#### Delete tracks with negative number for start site after resizing
for f in *.bed
do
	srt=`cat $f | awk '$2 < 0' | wc -l`
	if [ ${srt} -gt 0 ]
	then
		echo "ACTION: del" ${srt} "line(s) from" $f
		cat $f | sed '/-[0-9]/d' > tmp
		mv tmp ${f}
	fi
done

		
