We ran several of our bioinformatics analyses  in the Terra cloud computing environment from the Broaad Institute https://app.terra.bio/, using the following wdl workflows.


For sequencing run quality assessment, demultiplexing, and alignment we ran the following publically available workflows on Terra
FastQC  https://portal.firecloud.org/?return=terra#methods/Constantine_Workflows/FastQC_Constantine/18
bcl2fastq - https://portal.firecloud.org/?return=terra#methods/cumulus/bcl2fastq/5
Starsolo - https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/STARsolo:2.0.0?tab=info
Citeseq-count - https://portal.firecloud.org/?return=terra#methods/bemead/CITEseq-count/7

To run cNMF we used this publically available workflow on Terra
https://portal.firecloud.org/?return=terra#methods/mparikh/cNMF/20
