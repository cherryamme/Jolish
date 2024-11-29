### 
file=/home/jiangchen/project/Jolish/CYP21A2_Cyc/MCGD074-11-cyclone_correct.fq.gz

reference=/home/jiangchen/project/lrs_thal_pipeline/source/hg38_thal.mmi

bam=${file}.bam
# minimap2 -x map-ont -t 8 $reference $file > ${file}.paf
minimap2  \
            --MD -ax map-ont -t 4 -L --secondary=no -Y \
            $reference \
            $file | samtools sort -@2 -m 256Mb -T $(dirname $file)/temp -o ${bam} - && samtools index ${bam}


### 
### minimap2 hg38
file=/home/jiangchen/project/Jolish/CYP21A2_Cyc/MCGD074-11-cyclone_correct.fq.gz
file=/home/jiangchen/project/Jolish/CYP21A2_Cyc/MCGD092-06-ONT.fa

reference=/home/jiangchen/00_software/reference/GRCh38/hg38.mmi

bam=${file}.bam
# minimap2 -x map-ont -t 8 $reference $file > ${file}.paf
minimap2  \
            --MD -ax map-ont -t 4 -L --secondary=no -Y \
            $reference \
            $file | samtools sort -@2 -m 256Mb -T $(dirname $file)/temp -o ${bam} - && samtools index ${bam}


### 


### 
file=/home/jiangchen/project/Jolish/test_aa.fq
reference=/home/jiangchen/project/lrs_thal_pipeline/source/hg38_thal.mmi
reference_fa=/home/jiangchen/project/lrs_thal_pipeline/source/hg38_thal.fa
# 1. 生成比对文件
minimap2 -ax map-ont $reference $file > $file.sam

# 2. 运行 Racon 进行校正
racon-v1.4.21/build/bin/racon -m 8 -x -6 -g -8 -w 500 $file $file.sam $reference_fa -u > $file.racon.fq


### 


### 
file=/home/jiangchen/project/Jolish/test_aa.fq
THREADS=4
minimap2 -x ava-ont -w 20 -K 2g -f 0.005 -t $THREADS $file $file | awk '{ if($4 - $3 >= 0.6 * $2 || $9 - $8 >= 0.6 * $7) print $0}' > ovlp.paf
Fec -x 1 -t $THREADS -r 0.6 -a 400 -c 0 -l 1000 -m 0.005 -f 0.2 ovlp.paf $file corrected1.fasta
minimap2 -x ava-ont -w 20 -K 2g -f 0.005 -t $THREADS corrected1.fasta corrected1.fasta | awk '{ if($4 - $3 >= 0.2 * $2 || $9 - $8 >= 0.2 * $7) print $0}' > ovlp.paf
Fec -x 1 -R -t $THREADS -r 0.6 -a 1000 -c 4 -l 2000 -m 0.005 -f 0.2 ovlp.paf corrected1.fasta corrected2.fasta

### 



### 
cd /home/jiangchen/project/Jolish/benchmark
for i in 2000
do
name=benchmark_$i
rg -z '%-$' -A3 '/home/jiangchen/project/Jolish/example/3.7_aaD/alpha_THAL226-34_hg38_hap1.fq.gz' |head -4000 > ${name}_raw.fq
seqtk trimfq -b 2000 -L $i ${name}_raw.fq > $name.fq
done
### 




### hyperfine
name=jolish_10w
cd /home/jiangchen/project/Jolish/benchmark
data=/home/jiangchen/project/Jolish/benchmark/benchmark_
outfile=/home/jiangchen/project/Jolish/benchmark/result_
/home/jiangchen/00_software/hyperfine/hyperfine --show-output --export-csv ${name}.tsv --export-json ${name}.json --export-markdown ${name}.md --runs 3 -L length 2000,4000,8000 -L threads 4,10,20 "/home/jiangchen/project/Jolish/target/release/Jolish -i $data{length}.fq -t {threads} -o ${outfile}_t{threads}_{length}.fq.gz"

### 



### 

/home/jiangchen/project/Jolish/target/release/Jolish -i /home/jiangchen/project/Jolish/CYP21A2_Cyc/MCGD074-11-ont_5000.fq -c 100 -t 20 -m 500 -o /home/jiangchen/project/Jolish/CYP21A2_Cyc/MCGD074-11-cyclone_correct.fq.gz

### 


### 
file=/home/jiangchen/project/Jolish/benchmark/benchmark_1000.fq.bam
file=/home/jiangchen/project/Jolish/benchmark/result__t20_1000_c100.fq.gz.bam
file=/home/jiangchen/project/Jolish/benchmark/benchmark_1000_abpoa.fq.gz.bam

reference=/home/jiangchen/project/lrs_thal_pipeline/source/hg38_thal.fa

outdir=/home/jiangchen/project/Jolish/benchmark/counterr_abpoa

counterr -bam $file -genome $reference -output_dir $outdir -len_min_read 10 -mapq_thres 7 -len_min_contig 10 -len_min_aln 10

### 


### abpoachunk
abpoa_handle_py=/home/jiangchen/project/lrs_thal_pipeline/modules/abpoa_handle.py
abpoa=/home/jiangchen/project/lrs_thal_pipeline/modules/abpoa

python $abpoa_handle_py $abpoa /home/jiangchen/project/Jolish/benchmark/benchmark_1000.fq.gz 5 50 2> /dev/null | gzip > /home/jiangchen/project/Jolish/benchmark/benchmark_1000_abpoa.fq.gz

### 