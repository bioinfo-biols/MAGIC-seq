#!/bin/bash
sample=$1

sampath=`dirname ${sample}`
ifile1=`basename ${sample}`
array2=(${ifile1//_R1/ })
sample=${array2[0]}

fastq1=${sampath}/${sample}_R1.fastq.gz
fastq2=${sampath}/${sample}_R2.fastq.gz

project_path=/Project/MAGIC_seq/Mouse_Adult_Brain_10X
#Barcode-info
whitelist=${project_path}/data/Barcode-T3-70/whitelist
ID=${project_path}/data/Barcode-T3-70/T3-ids-barcode.txt
#result_path
st_path=${project_path}/data/result_STARsolo/${sample}
#ref-info
MAP=/database/ref_mm10_M23_release98/mm10_filtered_star_150bp_index
ANN=/database/ref_mm10_M23_release98/mm10_filtered.gtf

#Other-info
log_file=${st_path}/${sample}_st_log.txt
t_num=10
node=batch


echo "

source  /.bashrc
source activate st_pipeline

ulimit -n 4096


echo 'Step1-----seqkit------'
echo \$(date +%T)

seqkit split2 -1 ${fastq1} -2 ${fastq2} -p 10 -O ${sampath}/${sample}_split/out -f -j ${t_num}

mkdir -p ${st_path}

cd ${st_path}

mkdir ./STARsolo ./tmp

files=(001 002 003 004 005 006 007 008 009 010)


for file in \"\${files[@]}\"
do
    seqkit subseq -j 10 -r 1:28  ${sampath}/${sample}_split/out/${sample}_R1.part_\${file}.fastq.gz -o ${sample}_\${file}_reformat_R1.fastq.gz

    rm test*gz
done


cat *_reformat_R1.fastq.gz > ${sample}_reformat_R1.fastq.gz

rm ${sample}_0*_reformat_R1.fastq.gz

cat ${sampath}/${sample}_split/out/${sample}_R2.part_0*.fastq.gz > ${sample}_reformat_R2.fastq.gz

rm -r ${sampath}/${sample}_split

echo 'Step2-----st_pipeline------'
echo \$(date +%T)

/software/st_pipeline/scripts/st_pipeline_run.py \\
  --output-folder ./ \\
  --temp-folder ./tmp \\
  --umi-start-position 16 \\
  --umi-end-position 28 \\
  --ids $ID \\
  --expName ${sample} \\
  --ref-map ${MAP} \\
  --ref-annotation ${ANN} \\
  --verbose \\
  --threads ${t_num} \\
  --log-file ${log_file} \\
  --min-length-qual-trimming 30 \\
  --disable-mapping \\
  --disable-annotation \\
  --disable-barcode \\
  ${sample}_reformat_R1.fastq.gz ${sample}_reformat_R2.fastq.gz

echo 'Step3-----STARsolo------'
echo \$(date +%T)

/software/STAR/source/STAR  --genomeDir ${MAP} \\
  --outFileNamePrefix ${st_path}/STARsolo/${sample}_ \\
  --readFilesCommand cat \\
  --readFilesIn st_R2.fastq st_R1.fastq \\
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM sF \\
  --outSAMtype BAM SortedByCoordinate \\
  --limitBAMsortRAM 121539607552 \\
  --soloType CB_UMI_Simple \\
  --soloCBwhitelist ${whitelist} \\
  --soloCBstart 1 \\
  --soloCBlen 16 \\
  --soloUMIstart 17 \\
  --soloUMIlen 12 \\
  --soloFeatures Gene GeneFull SJ Velocyto \\
  --soloMultiMappers EM \\
  --soloUMIdedup 1MM_All \\
  --soloCellFilter EmptyDrops_CR \\
  --soloCellReadStats Standard \\
  --clipAdapterType CellRanger4 \\
  --outReadsUnmapped Fastx \\
  --runThreadN ${t_num} 



samtools index -@ ${t_num} ${st_path}/STARsolo/${sample}_Aligned.sortedByCoord.out.bam
time pigz -p ${t_num} ${st_path}/STARsolo/*_Solo.out/*/*/*mtx



echo 'Step4----rename-gzip------'
echo \$(date +%T)

rm ${sample}_reformat_*.fastq.gz 

mv st_R2.fastq ${sample}_st_R2.fastq 
mv st_R1.fastq ${sample}_st_R1.fastq 

time pigz -p ${t_num} *fastq 


mv STARsolo/${sample}_Unmapped.out.mate1 STARsolo/${sample}_Unmapped_R2.fastq 
mv STARsolo/${sample}_Unmapped.out.mate2 STARsolo/${sample}_Unmapped_R1.fastq 
time pigz -p ${t_num} STARsolo/*_Unmapped_*.fastq


echo 'End-----------'8
echo \$(date +%T)

chmod -R 775 ${st_path}

"  > ${project_path}/qsub/qsub_solo/${sample}.St.sh

echo "#PBS -N ${sample}
#PBS -l nodes=1:ppn=${t_num}
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -q ${node}

cd ${project_path}/qsub/qsub_solo/

bash ${sample}.St.sh

"  > ${project_path}/qsub/qsub_solo/${sample}.St.qsub
