for f in ../../peak/*.bed
do
	bn=`basename ${f}`
	rn=`echo ${bn} | cut -d "_" -f 1`
	echo ${rn} `cat ${f} | awk 'OFS="\t" {SUM += $3-$2} END {print SUM, SUM/NR}'` >> len_of_peak.txt
done

sleep 0.2

for f in ../../peak/ctcf.rad21/*.bed
do
	bn=`basename ${f}`
	rn=`echo ${bn} | cut -d "_" -f 1,2`
	echo ${rn} `cat ${f} | awk 'OFS="\t" {SUM += $3-$2} END {print SUM, SUM/NR}'` >> len_of_peak.txt
done

sleep 0.2

for f in ../../peak/rand_backgd/*.bed.gz
do
	bn=`basename ${f}`
	rn=`echo ${bn} | cut -d "." -f 1`
	echo ${rn} `zcat ${f} | awk 'OFS="\t" {SUM += $3-$2} END {print SUM, SUM/NR}'` >> rand_len_of_peak.txt
done
