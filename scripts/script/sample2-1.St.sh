
#-------Starting-Log-----#
set -e
echo "[$(date)] @$(hostname) ${PBS_JOBNAME} started"
start=$(date +%s)
#-------Starting-Log-----#

source  /histor/zhao/pangkun/.bashrc
source activate st_pipeline

ulimit -n 4096

mkdir -p /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1


echo 'Step0-----QC-reports and Extract barcode and UMI sequence---'
echo $(date +%T)

/histor/public/software/fastp/fastp \
    -i /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/rawdata/sample2-1/sample2-1_R1.fastq.gz -I /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/rawdata/sample2-1/sample2-1_R2.fastq.gz \
    -o /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_R1.cleaned.fq.gz \
    -O /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_R2.cleaned.fq.gz \
    -Q -L -G -A \
    -m --merged_out /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_se.cleaned.fq.gz \
    -w 20 --detect_adapter_for_pe \
    -j /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1.raw.fastp.json \
    -h /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1.raw.fastp.html
rm /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/*fq.gz

seqkit split2 -1 /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/rawdata/sample2-1/sample2-1_R1.fastq.gz -2 /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/rawdata/sample2-1/sample2-1_R2.fastq.gz -p 10 -O /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_split/out -f -j 20

cd /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1

# 定义文件列表
files=(001 002 003 004 005 006 007 008 009 010)

# 遍历文件列表
for file in "${files[@]}"
do
    seqkit subseq -j 20 -r 1:8  /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_split/out/sample2-1_R1.part_${file}.fastq.gz -o test1-8.fastq.gz
    seqkit subseq -j 20 -r 27:46  /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_split/out/sample2-1_R1.part_${file}.fastq.gz -o test27-46.fastq.gz
    seqkit concat test1-8.fastq.gz test27-46.fastq.gz -j 20 -o sample2-1_${file}_reformat_R1.fastq.gz
    rm test*gz
done

cat *_reformat_R1.fastq.gz > sample2-1_reformat_R1.fastq.gz
rm sample2-1_0*_reformat_R1.fastq.gz
cat /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_split/out/sample2-1_R2.part_0*.fastq.gz > sample2-1_reformat_R2.fastq.gz
rm -r /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_split

echo 'Step1-----fastp-QC------'
echo $(date +%T)

# TA 建库需要运行下面一行
cutadapt -g AAGCAGTGGTATCAACGCAGAGTGAATGGG -e 0.01 -j 20 -o sample2-1_reformat_cutadapt_R2.fastq.gz sample2-1_reformat_R2.fastq.gz

/histor/public/software/fastp/fastp \
    -i sample2-1_reformat_R1.fastq.gz -I sample2-1_reformat_cutadapt_R2.fastq.gz \
    -o /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_R1.cleaned.fq.gz \
    -O /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_R2.cleaned.fq.gz \
    -l 28 -x -g \
    -w 20 --detect_adapter_for_pe \
    -j /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1.barcode.fastp.json \
    -h /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1.barcode.fastp.html


echo 'Step2-----STARsolo------'
echo $(date +%T)

/histor/zhao/pangkun/software/STAR/source/STAR  --genomeDir /histor/zhao/pangkun/database/ref_mm10_M23_release98/ref_mm10_st_pipeline_150_index \
  --outFileNamePrefix /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/STARsolo_se/sample2-1_ \
  --readFilesCommand zcat \
  --readFilesIn /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_R2.cleaned.fq.gz /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/sample2-1_R1.cleaned.fq.gz \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM sF \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 121539607552 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/Barcode-T3-70/whitelist.txt \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloMultiMappers EM \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR \
  --soloCellReadStats Standard \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --outReadsUnmapped Fastx \
  --runThreadN 20 


samtools index -@ 20 /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/STARsolo_se/sample2-1_Aligned.sortedByCoord.out.bam

time pigz -p 20 /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/STARsol*/*_Solo.out/*/*/*mtx


echo 'Step3----rename-gzip------'
echo $(date +%T)

rm STARsolo_se/*_Unmapped.out.mate*

echo 'Step4----get-info------'
echo $(date +%T)

conda deactivate

cd /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1/STARsolo_se

samtools view -h sample2-1_Aligned.sortedByCoord.out.bam | grep -v 'CB:Z:-' > sample2-1_Aligned_AB.sortedByCoord.out.sam
samtools view -@ 20 -b sample2-1_Aligned_AB.sortedByCoord.out.sam > sample2-1_Aligned_AB.sortedByCoord.out.bam
rm sample2-1_Aligned_AB.sortedByCoord.out.sam
samtools index -@ 20 sample2-1_Aligned_AB.sortedByCoord.out.bam
samtools flagstat -@ 20 sample2-1_Aligned_AB.sortedByCoord.out.bam > sample2-1_flagstat.txt

echo 'End-----------'
echo $(date +%T)

chmod -R 775 /histor/zhao/pangkun/Zhaolab/Hubeiyu/20230712-T3-70-32um-ele/data/result_STARsolo/sample2-1

#-------Ending-Log-----#
echo "[$(date)] ${PBS_JOBID}@$(hostname) $PBS_JOBNAME finished"
end=$(date +%s)
secs=$(($end-$start))
printf 'Elapsed Time: %02dh:%02dm:%02ds\n' $((secs/3600)) $((secs%3600/60)) $((secs%60))
#-------Ending-Log-----#


