# biostars, https://www.biostars.org/p/278914/#278915
cat rmsk.txt | awk -v OFS="\t" '{ print $6, $7, $8, $12, $13, $10 }' \
	> rmsk.bed
