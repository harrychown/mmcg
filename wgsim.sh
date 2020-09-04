# Run wgsim and seqkt for MMCG

# Input:
# $1 length of first read 300 (-1)
# $2 length of second read 300 (-2)
# $3 error rate 0.002 (-e)
# $4 number of sim reads (-N)
# $5 number of seqkt reads
# $6 input file
# $7 output file prefix
# $8 read label
# $9 output directory

# -S set seed 662491

# Create output file names
out_1="${7}_1.fastq"
out_2="${7}_2.fastq"

wgsim -S 662491 -r 0 -R 0 -X 0 -1 $1 -2 $2 -e $3 -N $4 $6 $out_1 $out_2

# Count simulator reads
echo $(cat $out_1|wc -l)/4|bc
echo $(cat $out_2|wc -l)/4|bc

# Pass reads into seqkt for resampling
seqtk sample -s662491 $out_1 $5 > sub1.fastq 
seqtk sample -s662491 $out_2 $5 > sub2.fastq 

rm $out_1 $out_2

# Count subsampled reads
echo $(cat sub1.fastq|wc -l)/4|bc
echo $(cat sub2.fastq|wc -l)/4|bc

# Add accession ID to each read
sed -i "s/@/@$8_/" sub1.fastq 
cat sub1.fastq >> "${9}${out_1}"

sed -i "s/@/@$8_/" sub2.fastq 
cat sub2.fastq >> "${9}${out_2}"

rm sub1.fastq sub2.fastq


