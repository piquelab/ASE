# Pileup makefile
bamFolder=../bwa-mem

##aiScript=/nfs/hpnfs/piquelab/gmb/AI/results/ai_v3.5_noImpute.R
aiScript=/nfs/hpnfs/groups/piquelab/gmb/AI/results/ai_v3.6_noImpute.R
1KGSnps=/wsu/home/groups/piquelab/data/SNPs/1KGSnps.wg.cnvFilt1.uni.bed.gz
1KGSnpsAs=/wsu/home/groups/piquelab/data/SNPs/1KGSnps.wg.cnvFilt1.uni.as1k.bed.gz
genome=/wsu/home/groups/piquelab/data/RefGenome/hg19.fa
mappFile=/wsu/home/groups/piquelab/DNaseQTLs/mappability.x8b


filtStr='$$4>0 && $$5 !~ /[^\^][<>]/ && $$5 !~ /\+[0-9]+[ACGTNacgtn]+/ && $$5 !~ /-[0-9]+[ACGTNacgtn]+/ && $$5 !~ /[^\^]\*/'
bedStr='{ print $$1,$$2-1,$$2,$$3,$$4,$$5,$$6}'

numThreads=2
Qsub.ppn=2
Qsub.q=mmtxq
Qsub.N=pileup

# Rule to submit it as a job to the cluster 
%.Qsub:
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q $(Qsub.q) -l nodes=1:ppn=$(Qsub.ppn) -N $(Qsub.N) -o $@ -e $@.e

# Pileup reads from the bam file using 1KG SNPs
%.pileup.gz: $(bamFolder)/%.bam $(genome) $(1KGSnps)
	samtools rmdup $< - \
		| samtools view -h - \
		| tr ' ' _ | samFilterMappXbFile ${mappFile} stdin stdout 	\
		| samtools view -Sbu -q10 - 					\
		| samtools mpileup -f $(genome) -l $(1KGSnps) - | bgzip > $@

# Filter and intersect the pileup with the SNP alleles for reference
%.pileup.bed.gz: %.pileup.gz $(1KGSnpsAs)
	less $< | awk $(filtStr) | awk -v OFS='\t' $(bedStr) | sortBed -i stdin | intersectBed -a stdin -b $(1KGSnpsAs) -wo | gzip > $@

# Run the AI script and
%.pileup.clean.bed.gz: %.pileup.bed.gz
	-R --vanilla --args $^ $@ $*.pileup.loci.gz < /nfs/hpnfs/groups/piquelab/charvey/ASE/aiPreProcessing.v0.R  > $@ 2> $@.e

# combine all reads 
%




firstrule: $(patsubst %.pileup.bed.gz,%.pileup.clean.bed.gz.Qsub,$(wildcard *.pileup.bed.gz))

pileup_inCluster: $(patsubst $(bamFolder)/%.bam,%.pileup.bed.gz.Qsub,$(wildcard $(bamFolder)/*.bam))

pileup_all: $(patsubst $(bamFolder)/%.bam,%.pileup.bed.gz,$(wildcard $(bamFolder)/*.bam))



clean:
	rm -r $(patsubst $(bamFolder)/%.bam,%*,$(wildcard $(bamFolder)/*.bam))
