#### map all longreads to reference

reference="Pant_hap.fasta"

for file in all #NZ_102_02 NZ_120_02 NZ_123_05 NZ_124_02 NZ_126_03 NZ_134_01 NZ_143_01 NZ_146_04 NZ_148_01 NZ_155_01 NZ_159_03
do
ngmlr -t 64 -r $reference -q $file.fastq -o $file.n.sam
sambamba view -t 32 -S -f bam -o $file.n.bam $file.n.sam
sambamba sort -t 32 -o $file.sorted.n.bam $file.n.bam
sambamba index -t 32 $file.sorted.n.bam
done


#### run sniffles to detect structural variants

sniffles -t 32 -m all.sorted.n.bam -v all_sniffles.vcf



### new basecalling with guppy

./ont-guppy-cpu/bin/guppy_basecaller --compress_fastq -i rawdata/NZ_102_02/fast5_pass -s newbasecalls_2 \
--cpu_threads_per_caller 32 --num_callers 1 -c ont-guppy-cpu/data/dna_r10.3_450bps_sup.cfg


### using nanopolish


for chneck in NZ_102_02 NZ_120_02 NZ_123_05 NZ_124_02 NZ_126_03 NZ_134_01 NZ_143_01 NZ_146_04 NZ_148_01 NZ_155_01 NZ_159_03
do
nanopolish index -d rawdata/$chneck/fast5_pass/ reads_fastq/$chneck.fastq
nanopolish call-methylation -t 16 -r reads_fastq/$chneck.fastq -b mappings/$chneck.sorted.bam -g Pant_hap.fasta > $chneck.meth.tsv
done

### using nanopolish and MCaller

for READS in NZ_102_02 NZ_120_02 NZ_124_02 NZ_126_03 NZ_134_01 NZ_143_01 NZ_146_04 NZ_148_01 NZ_155_01 NZ_159_03
do
nanopolish eventalign -t 24 --scale-events -n -r reads_fastq/$READS.fastq -b mappings/$READS.sorted.bam -g Pant_hap.fasta > $READS.eventalign.tsv
#python ./mCaller/mCaller.py -m GATC -r Pant_hap.fasta -d ./mCaller/r95_twobase_model_NN_6_m6A.pkl -e $READS.eventalign.tsv -f reads_fastq/$READS.fastq -b A
done



############################################################################

### basecalling using megalodon:

megalodon \
    NZ_102_02/fast5_pass \
    --guppy-server-path /usr/bin/guppy_basecall_server \
    --guppy-config dna_r10.4_e8.1_modbases_5hmc_5mc_cg_fast.cfg \
    --remora-modified-bases dna_r10.4_e8.1 fast v3.3 5mc CG 1 \
    --outputs basecalls mappings mod_mappings mods \
    --reference Pant_hap.fasta \
    --processes 12


### results from modified base calling_   _5mC.bed
### intersect methylation with transcript mapping

for BED in all12 NZ_102_02 NZ_120_02 NZ_123_05 NZ_124_02 NZ_126_03 NZ_132_01 NZ_134_01 NZ_143_01 NZ_146_04 NZ_148_01 NZ_155_01 NZ_159_03
do
bedtools intersect -b $BED\_5mC.bed -a trans.bed > intersect_trans_$BED.bed
done

### contigs_count

for BED in NZ_102_02 NZ_120_02 NZ_123_05 NZ_124_02 NZ_126_03 NZ_132_01 NZ_134_01 NZ_143_01 NZ_146_04 NZ_148_01 NZ_155_01 NZ_159_03
do
cut -f 4 intersect_trans_$BED.bed | sort | uniq -c > trans_contigs_count_$BED.txt
done

###contigs_names

for BED in all12 NZ_102_02 NZ_120_02 NZ_123_05 NZ_124_02 NZ_126_03 NZ_132_01 NZ_134_01 NZ_143_01 NZ_146_04 NZ_148_01 NZ_155_01 NZ_159_03
do
cut -f 4 intersect_trans_$BED.bed | sort | uniq > trans_contigs_meth_$BED.txt
done


