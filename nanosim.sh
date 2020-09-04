# Run nanosim and seqkt for MMCG

# Input:
# $1 number of sim reads (-N)
# $2 number of seqkt reads
# $3 input file
# $4 output file prefix
# $5 read label
# $6 output directory


# Create output file names
out="${4}_nanopore.fasta"

cd ./
simulator.py genome -dna_type linear --seed 662491 -rg $3 -n $1 > std_out.txt 2>&1
rm std_out.txt

# Count simulator reads
echo $(grep -c ">" simulated_aligned_reads.fasta)

# Pass reads into seqkt for resampling
seqtk sample -s662491 simulated_aligned_reads.fasta $2 > sub.fasta 


rm simulated_aligned_reads.fasta

# Count subsampled reads
echo $(grep -c ">" sub.fasta)



sed -i "s/>/>$5_/" sub.fasta 
cat sub.fasta >> "${6}${out}"


rm sub.fasta







