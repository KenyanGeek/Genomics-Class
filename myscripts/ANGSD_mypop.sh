myrepo="/users/h/c/hchege/Genomics-Class"

mkdir ${myrepo}/myresults/ANGSD

output="${myrepo}/myresults/ANGSD"

mypop="XCS"

ls /data/project_data/RS_ExomeSeq/mapping/BWA/${mypop}_*sorted.rm*.bam >${output}/${mypop}_bam.list

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# ANGSD -b ${output}/${mypop}_bam.list \
# -ref ${REF} -anc ${REF} \
# -out ${output}/${mypop}_allsites \
# -nThreads 1 \
# -remove_bads 1 \
# -C 50 \
# -baq 1 \
# -minMapQ 20 \
# -minQ 20 \
# -setMinDepth 3 \
# -minInd 2 \
# -setMinDepthInd 1 \
# -setMaxDepthInd 17 \
# -skipTriallelic 1 \
# -GL 1 \
# -doCounts 1 \
# -doMajorMinor 1 \
# -doMaf 1 \
# -doSaf 1 \
# -doHWE 1 \
# # -SNP_pval 1e-6

# Now let's calculate 
ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} -anc ${REF} \
-out ${output}/${mypop}_folded_allsites \
-nThreads 1 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-setMinDepth 3 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-fold 1
# -SNP_pval 1e-6

#Get a rough first estimate of the sfs

realSFS ${output}/${mypop}_folded_allsites.saf.idx \
-maxIter 1000 -tole 1e-6 -P 1 \
> ${output}/${mypop}_outFold.sfs

#Get refined estimate of sfs and do

ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} -anc ${REF} \
-out ${output}/${mypop}_folded_allsites \
-nThreads 1 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-setMinDepth 3 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-fold 1 \
-pest ${output}/${mypop}_outFold.sfs \
-doThetas 1

#Use the doTheta output from above to estimate nucleotide diversity

thetaStat do_stat ${output}/${mypop}_folded_allsites.thetas.idx