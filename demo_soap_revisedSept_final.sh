#!/bin/sh
# this script runs soap denovo for 
#$ -N soapdn
#$ -S /bin/sh
#$ -M username@niaid.nih.gov
#$ -j y
#$ -m e
#$ -pe threaded 2
#$ -l h_vmem=10G

# for larger assembly try:
## -l h_vmem=10G,mem_free=200G,h_rt=48:00:00 
## -pe threaded 10

# full parameters for SOAPdenovo can be found in /usr/local/bio_apps/SOAPdenovo-1.05/SOAPdenovo-63mer all

cd ~/denovo/
export PATH=$PATH:/usr/local/bio_apps/SOAPdenovo-1.05/
export PATH=$PATH:/usr/local/bio_apps/bwa-0.6.1/
export PATH=$PATH:/usr/local/bio_apps/BEDTools-2.16.2/bin/


#### RUN SOAP DENOVO #######
SOAPdenovo-63mer all -s soap_config_sept2012.txt -d 10 -D 10 -L 200 -p 2 -K 27 -o soap_out/vcho27


#### GET N50 STATS FOR CONTIGS ####
/usr/local/bio_apps/abyss-1.3.2/bin/abyss-fac -t 100 soap_out/vcho27.contig > contig_stat.txt


#### DETERMINE HOW MUCH OF REFERENCE GENOME WAS COVERED ####
# get Contigs that are more than 1000 bp
perl get_contigs_ends.pl -e 0 -c 1k soap_out/vcho27.contig > soap_out/vcho27.contig_1k.fa

# get coverage
bwa bwasw -f soap_vibrio_bwa.sam vibrio/637000333.fna soap_out/vcho27.contig_1k.fa




# samtools to create index and generate bam file for viewing in browser
#/usr/local/bio_apps/samtools-0.1.16/samtools faidx vibrio/637000333.fna
samtools view -bt vibrio/637000333.fna.fai soap_vibrio_bwa.sam > soap_vibrio_bwa.bam
samtools sort soap_vibrio_bwa.bam soap_vibrio_bwa_sorted
samtools index soap_vibrio_bwa_sorted.bam


############

#Part C: compute coverage of genome by contigs
bedtools coverage -abam soap_vibrio_bwa_sorted.bam -b vibrio/genome_size.bed > soap_to_genome_coverage.txt
