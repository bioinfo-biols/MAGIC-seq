# MAGIC-seq Pipeline

This is a public repository for all code connected to MAGIC-seq.


The following files/parameters are commonly required :
- FASTQ files (Read 1 containing the spatial information and the UMI and read 2 containing the genomic sequence) 
- A genome index generated with STAR 
- An annotation file in GTF or GFF3 format (optional when using a transcriptome)
- The file containing the barcodes and array coordinates (look at the folder "ids" to use it as a reference). 


For the QC step, the stpipeline/common/filterInputReads.py in st_pipeline must be replaced with the supplied filterInputReads.py to generate the QC post-fastq for the subsequent steps.


## Preprocessing of MAGIC-seq raw data:

1. Triple/Nine grid chip data
```
#raw data
fastq1=${sampath}/${sample}_R1.fastq.gz
fastq2=${sampath}/${sample}_R2.fastq.gz
#Barcode-info
whitelist=${project_path}/data/Barcode-T9-70/whitelist
ID=${project_path}/data/Barcode-T9-70/T9-ids-barcode.txt
#result_path
st_path=${project_path}/data/result_STARsolo/${sample}
#ref-info
MAP=/database/ref_mm10_M23_release98/mm10_filtered_star_index
ANN=/database/ref_mm10_M23_release98/mm10_filtered.gtf
#Other-info
log_file=${st_path}/${sample}_st_log.txt
t_num=32

#1.Data preprocessing: extracting barcode and UMI sequences from read1
seqkit split2 -1 ${fastq1} -2 ${fastq2} -p 10 -O ${sampath}/${sample}_split/out -f -j ${t_num}
mkdir -p ${st_path}
cd ${st_path}
mkdir ./STARsolo ./tmp
files=(001 002 003 004 005 006 007 008 009 010)

for file in \"\${files[@]}\"
do
    seqkit subseq -j ${t_num} -r 1:8  ${sampath}/${sample}_split/out/${sample}_R1.part_\${file}.fastq.gz -o test1-8.fastq.gz
    seqkit subseq -j ${t_num} -r 27:34  ${sampath}/${sample}_split/out/${sample}_R1.part_\${file}.fastq.gz -o test27-34.fastq.gz
    seqkit subseq -j ${t_num} -r 53:72  ${sampath}/${sample}_split/out/${sample}_R1.part_\${file}.fastq.gz -o test53-72.fastq.gz
    seqkit concat test1-8.fastq.gz test27-34.fastq.gz -j ${t_num} -o test.fastq.gz
    seqkit concat test.fastq.gz test53-72.fastq.gz -j ${t_num} -o ${sample}_\${file}_reformat_R1.fastq.gz
    rm test*gz
done

cat *_reformat_R1.fastq.gz > ${sample}_reformat_R1.fastq.gz
rm ${sample}_0*_reformat_R1.fastq.gz
cat ${sampath}/${sample}_split/out/${sample}_R2.part_0*.fastq.gz > ${sample}_reformat_R2.fastq.gz
rm -r ${sampath}/${sample}_split

#2.Quality control of sequencing data
/software/st_pipeline/scripts/st_pipeline_run.py \
  --output-folder ./ \
  --temp-folder ./tmp \
  --umi-start-position 16 \
  --umi-end-position 28 \
  --ids $ID \
  --expName ${sample} \
  --ref-map ${MAP} \
  --ref-annotation ${ANN} \
  --verbose \
  --threads ${t_num} \
  --log-file ${log_file} \
  --min-length-qual-trimming 30 \
  --disable-mapping \
  --disable-annotation \
  --disable-barcode \
  ${sample}_reformat_R1.fastq.gz ${sample}_reformat_R2.fastq.gz
#3.Generate gene expression matrix
/software/STAR/source/STAR  --genomeDir ${MAP} \
  --outFileNamePrefix ${st_path}/STARsolo/${sample}_ \
  --readFilesCommand cat \\
  --readFilesIn st_R2.fastq st_R1.fastq \\
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM sF \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 121539607552 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist ${whitelist} \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloMultiMappers EM \
  --soloUMIdedup 1MM_All \
  --soloCellFilter EmptyDrops_CR \
  --soloCellReadStats Standard \
  --clipAdapterType CellRanger4 \
  --outReadsUnmapped Fastx \
  --runThreadN ${t_num} 


```
2. Splicing grid chip data
```
#raw data
fastq1=${sampath}/${sample}_R1.fastq.gz
fastq2=${sampath}/${sample}_R2.fastq.gz
#Barcode-info
whitelist=${project_path}/data/Barcode-M9-70/whitelist
ID=${project_path}/data/Barcode-M9-70/M9_ST_ids_barcode_chip1_C18.txt
#result_path
st_path=${project_path}/data/result_STARsolo/${sample}
#ref-info
MAP=/database/ref_mm10_M23_release98/mm10_filtered_star_index
ANN=/database/ref_mm10_M23_release98/mm10_filtered.gtf
#Other-info
log_file=${st_path}/${sample}_st_log.txt
t_num=32

#1.Data preprocessing: extracting barcode and UMI sequences from read1
seqkit split2 -1 ${fastq1} -2 ${fastq2} -p 10 -O ${sampath}/${sample}_split/out -f -j ${t_num}
mkdir -p ${st_path}
cd ${st_path}
mkdir ./STARsolo ./tmp

files=(001 002 003 004 005 006 007 008 009 010)
for file in \"\${files[@]}\"
do
    seqkit subseq -j ${t_num} -r 1:8  ${sampath}/${sample}_split/out/${sample}_R1.part_\${file}.fastq.gz -o test1-8.fastq.gz
    seqkit subseq -j ${t_num} -r 27:34  ${sampath}/${sample}_split/out/${sample}_R1.part_\${file}.fastq.gz -o test27-34.fastq.gz
    seqkit subseq -j ${t_num} -r 53:72  ${sampath}/${sample}_split/out/${sample}_R1.part_\${file}.fastq.gz -o test53-72.fastq.gz
    seqkit concat test1-8.fastq.gz test27-34.fastq.gz -j ${t_num} -o test.fastq.gz
    seqkit concat test.fastq.gz test53-72.fastq.gz -j ${t_num} -o ${sample}_\${file}_reformat_R1.fastq.gz
    rm test*gz
done

cat *_reformat_R1.fastq.gz > ${sample}_reformat_R1.fastq.gz
rm ${sample}_0*_reformat_R1.fastq.gz
cat ${sampath}/${sample}_split/out/${sample}_R2.part_0*.fastq.gz > ${sample}_reformat_R2.fastq.gz
rm -r ${sampath}/${sample}_split

#2.Quality control of sequencing data
/software/st_pipeline/scripts/st_pipeline_run.py \
  --output-folder ./ \
  --temp-folder ./tmp \
  --umi-start-position 24 \
  --umi-end-position 36 \
  --ids $ID \
  --expName ${sample} \
  --ref-map ${MAP} \
  --ref-annotation ${ANN} \
  --verbose \
  --threads ${t_num} \
  --log-file ${log_file} \
  --min-length-qual-trimming 30 \
  --disable-mapping \
  --disable-annotation \
  --disable-barcode \
  ${sample}_reformat_R1.fastq.gz ${sample}_reformat_R2.fastq.gz
#3.Generate gene expression matrix
/software/STAR/source/STAR  --genomeDir ${MAP} \
--outFileNamePrefix ${st_path}/STARsolo/${sample}_ \
--readFilesCommand cat \
--readFilesIn st_R2.fastq st_R1.fastq \\
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM sF \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 121539607552 \
--soloType CB_UMI_Simple \
--soloCBwhitelist ${whitelist} \
--soloCBstart 1 \
--soloCBlen 24 \
--soloUMIstart 25 \
--soloUMIlen 12 \
--soloFeatures Gene GeneFull SJ Velocyto \
--soloMultiMappers EM \
--soloUMIdedup 1MM_All \
--soloCellFilter EmptyDrops_CR \
--soloCellReadStats Standard \
--clipAdapterType CellRanger4 \
--outReadsUnmapped Fastx \
--runThreadN ${t_num}

```

##  Convert  gene expression matrix to anndata

1. Triple/Nine grid chip data (T9-50um-70barcode)
```
from st_processing import get_adata_STARsolo
sample='sample1-1'
file_path='/MAGIC_seq/Mouse_Adult_Organ_T9_70_50um/data/result_STARsolo/'+sample+'/STARsolo/'
image_file_path='/MAGIC_seq/Mouse_Adult_Organ_T9_70_50um/data/Image/'
reg='reg7'
HE_point1=(384,614)
Spot_point1=(369,556)
HE_point2=(4623,4742)
Spot_point2=(4634, 4634)
line_point_r1c1=(746,676) 
line_point_r1c70=(773,4194)
line_point_r70c1=(4266,632)
adata=get_adata_STARsolo(sample,chip_type,reg,channels_num,barcode_num,res_um=res_um,only_HE=False,
                        img_file=True,EM=True,Velocyto=False,soloFeatures='GeneFull',raw=True,species='mouse',
                        file_path=file_path,image_file_path=image_file_path,Barcode_file_path=Barcode_file_path,geneinfo_path=geneinfo_path,
                        HE_point1=HE_point1,Spot_point1=Spot_point1,HE_point2=HE_point2,Spot_point2=Spot_point2,
                        line_point_r1c1=line_point_r1c1,line_point_r1c70=line_point_r1c70,line_point_r70c1=line_point_r70c1)
```
2. Splicing grid chip data(M9-15um-70barcode)
```
from st_processing import get_adata_STARsolo
sample='ZOO9-3'
file_path='/MAGIC_seq/Mouse_Adult_Brain_M9_70_15um/data/result_STARsolo/'+sample+'/STARsolo/'
image_file_path='/MAGIC_seq/Mouse_Adult_Brain_M9_70_15um/data/Image/'
Barcode_file_path='/MAGIC_seq/Mouse_Adult_Brain_M9_70_15um/data/Barcode-M9-70/'

res_um=15 
channels_num=210 
barcode_num=70 
reg='reg1' 
chip_type='M9'

HE_point1=(138,169)
Spot_point1=(112,227)
HE_point2=(4840,4691)
Spot_point2=(4827, 4720)
line_point_r1c1=(735,694) 
line_point_r1c70=(740,4426)
line_point_r70c1=(4461,669)

adata=get_adata_STARsolo(sample,chip_type,reg,channels_num,barcode_num,res_um=res_um,add_M9_15=True,
                        img_file=True,EM=True,Velocyto=True,soloFeatures='GeneFull',raw=True,species='mouse',
                        file_path=file_path,image_file_path=image_file_path,Barcode_file_path=Barcode_file_path,geneinfo_path=geneinfo_path,
                        HE_point1=HE_point1,Spot_point1=Spot_point1,HE_point2=HE_point2,Spot_point2=Spot_point2,
                        line_point_r1c1=line_point_r1c1,line_point_r1c70=line_point_r1c70,line_point_r70c1=line_point_r70c1)
```


