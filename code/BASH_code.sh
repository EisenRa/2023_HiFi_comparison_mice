###

################################################################################
#Downsampled short read analysis:
module load seqtk/1.3 pigz/2.3.4
for i in *_1.fastq.gz; do seqtk sample -s 1337 $i 5700000 > ${i/_1.fastq.gz/_ss_1.fastq}; done
for i in *_2.fastq.gz; do seqtk sample -s 1337 $i 5700000 > ${i/_2.fastq.gz/_ss_2.fastq}; done
for i in *_1.fastq.gz; do seqtk sample -s 1337 $i 6300000 > ${i/_1.fastq.gz/_ss_1.fastq}; done
for i in *_2.fastq.gz; do seqtk sample -s 1337 $i 6300000 > ${i/_2.fastq.gz/_ss_2.fastq}; done
for i in C04M3/*_1.fastq.gz; do seqtk sample -s 1337 $i 11400000 > $(basename ${i/_1.fastq.gz/_ss40gb_1.fastq}); done
for i in C04M3/*_2.fastq.gz; do seqtk sample -s 1337 $i 11400000 > $(basename ${i/_2.fastq.gz/_ss40gb_2.fastq}); done
for i in C13F5/*_1.fastq.gz; do seqtk sample -s 1337 $i 12600000 > $(basename ${i/_1.fastq.gz/_ss40gb_1.fastq}); done
for i in C13F5/*_2.fastq.gz; do seqtk sample -s 1337 $i 12600000 > $(basename ${i/_2.fastq.gz/_ss40gb_2.fastq}); done

pigz -p 40 */*.fastq

snakemake -s 0_Code/2_Coassembly_Binning.snakefile --use-conda --conda-frontend conda --configfile 0_Code/2_Assembly_Binning_config.yaml --cores 40

#Map short reads to long-read assemblies
for i in ../../Short_reads/Coassembly_binning/2_Reads/3_Host_removed/C04M3_subsampled/*_1.fastq.gz; do bowtie2 \
--threads 40 \
-x C04M3_LongReadAssembly.fa \
-1 $i \
-2 ${i/_1.fastq.gz/_2.fastq.gz} \
| samtools sort -@ 40 -o $(basename ${i/_1.fastq.gz/.bam});
done

for i in ../../Short_reads/Coassembly_binning/2_Reads/3_Host_removed/C13F5_subsampled/*_1.fastq.gz; do bowtie2 \
--threads 40 \
-x C13F5_LongReadAssembly.fa \
-1 $i \
-2 ${i/_1.fastq.gz/_2.fastq.gz} \
| samtools sort -@ 40 -o $(basename ${i/_1.fastq.gz/.bam});
done



#HybridSPAdes

metaspades.py -m 185 -t 40 --pacbio /home/projects/ku-cbd/people/antalb/mice_pacbio/m64241e_211217_170900.hifi_reads.fastq.gz \
-1 C04M3_20Gbp_1.fastq.gz -2 C04M3_20Gbp_2.fastq.gz -k 21,33,55,77,99 -o C04M3_20Gbp_HybridSPAdes

metaspades.py -m 185 -t 40 --pacbio /home/projects/ku-cbd/people/antalb/mice_pacbio/m64241e_211219_040338.hifi_reads.fastq.gz \
-1 C13F5_20Gbp_1.fastq.gz -2 C13F5_20Gbp_2.fastq.gz -k 21,33,55,77,99 -o C13F5_20Gbp_HybridSPAdes



#dRep
for i in MAGs/*/; do dRep dereplicate $(basename $i) -p 40 -comp 70 -sa 0.98 --genomeInfo MAGs/genome_info.csv -g $i/*.fa; done

#dRep all
dRep dereplicate ALL_longread_META_redo -comp 70 -sa 0.98 --genomeInfo MAGs_redo_longread_META/genome_info_new.csv -p 40 -g MAGs_redo_longread_META/*/*.fa.gz

#Flye has non-speicifc contig headers, rename them
for i in *.fa.gz; do rename.sh in=$i out=${i/.fa.gz/_renamed.fa.gz} prefix=${i/.fa.gz/-}; done

#Map to drepped genomes
cat dereplicated_genomes/*_renamed.fa.gz > catted_ref.fna.gz

bowtie2-build --threads 40 catted_ref.fna.gz catted_ref.fna.gz

for i in ../../../Short_reads/Coassembly_binning/2_Reads/3_Host_removed/*subsampled/*_1.fastq.gz;
do bowtie2 --threads 40 -x catted_ref.fna.gz -1 $i -2 ${i/_1/_2} | samtools sort -@ 40 -o $(basename ${i/_1.fastq.gz/.bam});
done

coverm genome -b *.bam -s - -m count covered_fraction length -t 40 --min-covered-fraction 0 > CoverM_MAG_table_META.tsv


#Map long reads to dereplicated MAG catalogue
module load minimap2/2.24r1122
module load samtools/1.14

minimap2 -ax map-hifi \
catted_ref.fna.gz \
-t 40 \
m64241e_211219_040338.hifi_reads.fastq.gz \
| samtools sort -o long_read_mapping/C13M5.bam

minimap2 -ax map-hifi \
catted_ref.fna.gz \
-t 40 \
m64241e_211217_170900.hifi_reads.fastq.gz \
| samtools sort -o long_read_mapping/C04M3_host.bam



#Gtdbtk
export GTDBTK_DATA_PATH=/home/projects/ku-cbd/people/rapeis/0_DBs/release202

for i in *Gb*; do gtdbtk classify_wf --genome_dir $i/dereplicated_genomes --extension .fa --prefix $i --cpus 40 --scratch_dir . --out_dir "$i"_gtdbtk; done

for i in */classify/*tsv; do sed '1d;' $i | cat; done >> gtdbtk_merged_temp.tsv
head -1 20Gb_long_gtdbtk/classify/20Gb_long.bac120.summary.tsv > header.tsv
cat header.tsv gtdbtk_merged_temp.tsv > gtdbtk_merged.tsv



#DRAM
module load tools miniconda3/4.10.1
conda activate DRAM

for i in `seq 1 10`; do mkdir -p group$i; find dereplicated_genomes/ -type f | head -n 23 | xargs -i mv "{}" group$i; done
for i in `seq 2 10`; do sed -i'' "s/group1/group$i/" DRAM-batch$i.sh; done
for i in `seq 2 10`; do sed -i'' "s/batch1/batch$i/" DRAM-batch$i.sh; done
for i in `seq 2 10`; do qsub180 DRAM-batch$i.sh; done

DRAM.py merge_annotations -i batch* -o merged_annotations
DRAM.py distill -i merged_annotations/annotations.tsv --rrna_path merged_annotations/rrnas.tsv --trna_path merged_annotations/trnas.tsv -o DRAM_distillate




# hifiasm analyses

## assembly
#!/bin/sh
#SBATCH -c 32 --mem 372G --time 24:00:00 # number of cores
/projects/mjolnir1/people/ncl550/0_software/hifiasm-meta/hifiasm_meta -t32 -o m64241e_211217_170900 m64241e_211217_170900.hifi_reads.fastq.gz 2>asm_m6
4241e_211217_170900.log

/projects/mjolnir1/people/ncl550/0_software/hifiasm-meta/hifiasm_meta -t32 -o m64241e_211219_040338 m64241e_211219_040338.hifi_reads.fastq.gz 2>asm_m6
4241e_211219_040338.log


for i in *p_ctg.gfa; do cat $i | awk '$1=="S" && ($2 ~ /.c$/) {printf ">%s\n%s\n", $2, $3} ' | gzip -1 > "$i"_CIRC_contigs.fa.gz; done
for i in *p_ctg.gfa; do cat $i | awk '$1=="S" {printf ">%s\n%s\n", $2, $3} ' | gzip -1 > "$i"_ALL_contigs.fa.gz; done



##renaming contigs
rename.sh in=m64241e_211217_170900.rescue.fa out=C04M3_longread_hifiasm_CIRCULAR.fa prefix=C04M3_longread_hifiasm_CIRCULAR addprefix=t
rename.sh in=m64241e_211219_040338.rescue.fa out=C13F5_longread_hifiasm_CIRCULAR.fa prefix=C13F5_longread_hifiasm_CIRCULAR addprefix=t

rename.sh in=m64241e_211217_170900.rescue.all.fa out=C04M3_longread_hifiasm_ALL.fa prefix=C04M3_longread_hifiasm_ALL addprefix=t
rename.sh in=m64241e_211219_040338.rescue.all.fa out=C13F5_longread_hifiasm_ALL.fa prefix=C13F5_longread_hifiasm_ALL addprefix=t


##Splitting circular bins
seqkit split C04M3_longread_hifiasm_CIRCULAR.fa -i
seqkit split C13M5_longread_hifiasm_CIRCULAR.fa -i
seqkit split C04M3_longread_hifiasm_ALL.fa -i
seqkit split C13M5_longread_hifiasm_ALL.fa  -i




##CheckM 
checkm lineage_wf -x fa C04M3_longread_hifiasm_CIRCULAR.fa.split/ C04M3_longread_hifiasm_CIRCULAR_checkm --tab_table -f C04M3_longread_hifiasm_CIRCULA
R_checkm.tsv -t 32 --pplacer_threads 8

for i in *checkm.tsv; do cut -f1,12,13 $i --output-delimiter=',' > ${i/.tsv/.csv}; done
for i in *.csv; do sed '1d;' $i >> genomes.csv; done
sed -i'' 's/,/.fa,/' genomes.csv
echo -e "genome,completeness,contamination" > header.csv
cat header.csv genomes.csv > checkm_for_drep.csv

#get stats for contigs:
for i in C04M3_longread_hifiasm_ALL.fa.split/*.fa; do echo $(basename $i) >> names.tsv && grep '>' $i | wc -l >> ncontigs.tsv && grep -v '>' $i | wc -c >> length.tsv; done




# PACBIO HIFI PIPELINE:
cut -f1,7,8 C04M3.HiFi-MAG.summary.txt > C04M3_metawrap_70_10.stats
cut -f1,7,8 C13F5.HiFi-MAG.summary.txt > C13F5_metawrap_70_10.stats



#MAG mapping rate 
for i in MAGs/*.fa.gz; do bowtie2-build --threads 32 $i $i; done


for i in MAGs/*.fa.gz; do for r in Reads/*_1.fastq.gz; do bowtie2 --threads 32 -x $i -1 $r -2 ${r/_1.fa/_2.fa} | samtools view -b -@ 32 - | samtools sort -@ 32 -o BAMs/$(basename "${r/_1.fastq.gz/}")_"$(basename ${i/.fa.gz/.bam})"; done; done


# long read renaming
rename.sh in=C04M3_hifi_mags.fa.gz out=C04M3_hifi_mags_renamed.fa.gz zl=9 prefix=contig^


619 mags pre dRep

Barnap 16S
for i in *.fa; do barrnap --threads 4 --incseq < $i > ${i/.fa/_RNAs.fa} 2> ${i/.fa/_STATS.tsv}; done

for i in *STATS.tsv; do grep 'Found: 16S_rRNA' $i >> ${i/_STATS/_16SrRNAs}; done

for file in *_16SrRNAs.tsv; do name=${file/_16SrRNAs.tsv/}; temp_file=$(mktemp); while IFS=$'\t' read -r -a line; do echo -e "${line[*]}\t$name" >> "$temp_file"; done < "$file"; mv "$temp_file" "$file"; done
cat *_16SrRNAs.tsv > barrnap_results.tsv

ANTISMASH
