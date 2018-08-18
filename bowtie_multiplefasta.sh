for sample in `ls /media/sample/fastqfiles/*R1.fastq`
do
dir="/media/sample/fastqfiles"
base=$(basename $sample "_R1.fastq")
bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq -2 ${dir}/${base}_R2.fastq -S ${dir}/${base}.sam
done
echo "bowtie2 -x path_to_my_index -1 ${dir}/${base}_R1.fastq -2 ${dir}/${base}_R2.fastq -S ${dir}/${base}.sam"

