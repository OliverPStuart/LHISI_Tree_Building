
FASTA=$1
OUT=$(echo ${FASTA} | cut -d "." -f1 | sed 's/$/\.nex/g')


N_TAXA=$(grep ">" ${FASTA} | wc -l)
N_CHAR=$(awk 'NR==2 {print length}' ${FASTA})

# Print taxa block

echo "#NEXUS" > $OUT
echo "BEGIN TAXA;" >> $OUT
echo " DIMENSIONS NTAX=${N_TAXA};" >> $OUT
echo " TAXLABELS" >> $OUT
bioawk -c fastx '{print $name}' ${FASTA} | sed 's/^/  /g' >> $OUT
echo "  ;" >>  $OUT
echo "END;" >> $OUT

# Print data block

echo "BEGIN CHARACTERS;" >> $OUT
echo " DIMENSIONS NCHAR=${N_CHAR};" >> $OUT
echo " FORMAT" >> $OUT
echo "  DATATYPE=protein" >> $OUT
echo "  MISSING=?"  >> $OUT
echo "  GAP=-;"  >> $OUT
echo " MATRIX" >> $OUT
bioawk -c fastx '{print "  ",$name,"\t",$seq}' ${FASTA} >> $OUT
echo "  ;" >> $OUT
echo "END;" >> $OUT
