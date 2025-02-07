
3.1. `LRAA - Long Read Alignment Assembler`
============================================

LRAA is a isoform discovery and quantification tool developed in-house at the Methods Development Lab 
to perform accurate and efficient ID and isoform quantification based on long read isoform sequencing alignments. 
It supports both PacBio and ONT reads and can be run 3 modes as below:

   - De-novo (reference annotation-free) isoform identification and quantification
   - Reference-guided isoform detection and quantification
   - Isoform expression quantification only

For analyzing the single cell dataset the recommended mode is Ref-Guided ID and quantification.


      | Dockstore : `LRAA.wdl <https://dockstore.org/workflows/github.com/MethodsDev/LongReadAlignmentAssembler/LRAA>`_
      | GitHub: `LRAA:main <https://github.com/MethodsDev/LongReadAlignmentAssembler>`_
      | Test Data can be found here (public, requester-pays) : `add file path`



**Example of input arguments for LRAA workflow for Ref Guided ID and quantifaction**

.. code:: bash

  :number-lines: 
  
  {
    "LRAA_wf.sample_id": "PBMC_BioIVT_10x3p_complete_LRAA",
    "LRAA_wf.main_chromosomes": "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",
    "LRAA_wf.inputBAM": "${this.minimap2_bam}",
    "LRAA_wf.numThreads": "${}",
    "LRAA_wf.referenceGenome": "gs://mdl-refs/GRCh38/GRCh38_no_alt.fa",
    "LRAA_wf.annot_gtf": "gs://mdl-refs/GRCh38/GRCh38.gencode.v39.annotation.gtf"
  }
