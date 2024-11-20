Full Length Bulk Downstream Processing
=======================================

After performing QC on the aligned reads using LongRNAseqQC, 
we can proceed to perform transcript analysis using an array of tools put togther by folks using Long read sequencing.

There are multiple applications one can use for many of these sub-steps, below is the currenty recommended workflow.

High-level Workflow for Bulk:
-----------------------------
The current recommended high-level workflow is as below:

 - Performing Sample QC and PCA with `Ptr` utility from `Trinity toolkit <https://github.com/trinityrnaseq/trinityrnaseq/wiki>`_.
 - `minimap2 <https://lh3.github.io/minimap2/minimap2.html>`_ for aligning reads
 - `Isoquant <https://github.com/ablab/IsoQuant>`_ for Reference guided Isoform Discovery.
 - `Stringtie <https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual>`_ merge to merge reconstructed assemblies to generate new Ref
 - Quantification with `Isoquant <https://github.com/ablab/IsoQuant>`_ against reconstructed ref
 - `Gffcompare <https://github.com/gpertea/gffcompare>`_ to compare against reference and fetch Ref Ids
 - Differential Expression analysis using DEseq2 implementation from the Trinity
 - Exploring isoform switches and isoform visualization with isoformSwitchAnalysisR
 - MDLgenomeview for inspecting read support


`Test Data`
~~~~~~~~~~~
The test dataset used for demonstarting downstream processing worklfows for Full Length correpsond to HG002 from the `Genome-in-a-Bottle project <https://www.nist.gov/programs-projects/genome-bottle>`_ 
and `UHRR Universal Human Reference RNA <https://www.thermofisher.com/order/catalog/product/QS0639>`_ cell lines. 

(Add information on the kits used here)

`minimap2`
~~~~~~~~~~
The minimap2 workflow, as detailed below, maps long reads or their assemblies to a reference genome.
The details and parameters can be found in the `minimap2 manual <https://lh3.github.io/minimap2/minimap2.html>`_ 

Workflow configuration for runnning the minimap2 workflow over cloud platforms supporting Cromwell like Terra can be found here:-
The workflow is designed to be run on a sample. Technical replicates must be merged before.


      | Dockstore : `Minimap2_LR.wdl <https://dockstore.org/workflows/github.com/broadinstitute/MDL-workflows/Minimap2_LR>`_
      | Github: `minimap2_LR <https://github.com/broadinstitute/MDL-workflows/blob/main/LR-tools/minimap2_LR/minimap2_LR.wdl>`_
      | Test Data can be found here (public, requester-pays) : `add file path`


.. csv-table:: minimap2
   :file: ../_subpages/tables/minimap2.csv
   :header-rows: 1


**Example of input arguments for minimap2 workflow for alignment with Human Ref genome**

.. code:: bash
  :number-lines: 


   {
    "Minimap2_LR.Minimap2Task.cpu" : "${8}",
    "Minimap2_LR.Minimap2Task.diskSizeGB" : "${500}",
    "Minimap2_LR.inputReads" : "${this.ubam}",
    "Minimap2_LR.referenceGenome" : "gs://mdl-refs/GRCh38/GRCh38_no_alt.fa",
    "Minimap2_LR.juncBED" : "gs://mdl-refs/GRCh38/GRCh38.gencode.v39.annotation.sorted.bed",
    "Minimap2_LR.sampleName" : "${this.sample_id}",
    "Minimap2_LR.readType" : "PacBioIsoSeq",
    "Minimap2_LR.customArguments" : "-G 1250k",
    "Minimap2_LR.keepUnmapped" : "true",
    "Minimap2_LR.allowSecondary" : "false",
    "Minimap2_LR.preemptible_tries" : "${3}"
  }

`Isoquant - Isoform Discovery`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. csv-table:: Isoquant
  :file: ../_subpages/tables/isoquant.csv
  :widths: 20,25,55
  :header-rows: 1

`Stringtie`
~~~~~~~~~~~
`StringTie tool <https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual>`_  is run with the --merge option, with a list of GTF/GFF files as inpout. It merges/assembles the input transcripts into a non-redundant set of transcripts. 
In the workflow stringtie is used to merge the reconstructed GTFs from Isoquant ID to create a new referrence.

Workflow configuration for runnning GffCompare on Terra can be found here:-

      | Dockstore : `stringtie_merge:main <https://dockstore.org/workflows/github.com/broadinstitute/MDL-workflows/StringTieMerge>`_
      | Github: `stringtie_merge_and_reestimate.wdl <https://github.com/broadinstitute/MDL-workflows/blob/main/LR-tools/stringtie_merge/stringtie_merge_and_reestimate.wdl>`_
      | Test Data can be found here (public, requester-pays) : `add file path` 


`Differential Expression with DEseq2`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For performing differential expression analysis, `DEseq2 <http://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_  implementation from Trinity RNA-seq project is detailed here.
Trinity provides support for other differential expression analysis tools namely, edgeR, limma/voom and ROTS.

More details and in-depth instructions on using the Trinity package can be found on the 
`Trinity wiki <https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Differential-Expression>`_

Trinity Toolkit also supports analysis of short read data. The overall workflow is modified to leverage short-read data if paired short read data for the samples is present.


.. code:: bash
  :number-lines:

  docker run -it -v /Users/usename/local_data_dir:/mnt/data trinityrnaseq/trinityrnaseq:latest
  /usr/local/bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
  --matrix /mnt/data/combined_transcript_counts_matrix.tsv \
  --method DESeq2 \
  --samples_file samples_desc.txt


`Gffcompare`
~~~~~~~~~~~~
`GffCompare is a utility <https://ccb.jhu.edu/software/stringtie/gffcompare.shtml>`_ used to compare two GTF/GFF files, which in referrnce based ID, is a reconstructed GTF resulted from merging individual GTFs from Isoquant ID with stringtie to the referrence annotation GTF.
The quick command is as below:


.. code:: bash
  :number-lines:
  
  conda create -n gffcompare bioconda::gffutils
  conda activate gffcompare
  gffcompare stringtie_merged.gtf -r gencode.vM32.annotation.gtf  

The tracking file generated in results conatins the matching trasncripts between samples. 
As GffCompare here is run with `-r` option, the 3rd column contains information about the reference annotation transcript.

`isoformSwitchAnalysisR`
~~~~~~~~~~~~~~~~~~~~~~~~~

The `isoformSwitchAnalyzeR <https://www.bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html>`_ is an R package developed to enable statistical identification of isoform switching and alternative isoform usage.
The vignete here is based on the `isoformAnlayzeR tutorial provided <https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html>`_

.. toctree::
   :maxdepth: 1 
   :caption: isoformSwitchAnalysisR:

   ./isoformSwitch_bulk


`Generating functional annotations`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To generate switch plots with reference annotation with isoformSwitchAnlysisR we can supply the annotations generated 
using various tools listed below to the switchObject generated by isoformSwitchAnalysisPart1(). 

Pfam annotations are required, in addition we can provide annotations generated with the tools below.
CPC2 Coding Potential Calculator : https://cpc2.gao-lab.org/

Pfam - domain annotation `pfam_scan.pl -as -dir isoformSwitchAnalysisPart1_results -fasta isoformSwitchAnalyzeR_isoform_AA_complete.fasta -cpu 4 -e_seq 10.0 -e_dom 10.0 > Pfam_result.txt`

IUPred Intrinsically disordered proteins (IDPs) : https://iupred2a.elte.hu/

SignalP Signal peptide and cleavage sites in gram+, gram- and eukaryotic amino acid sequences (signal pipetide at N terminus) : https://services.healthtech.dtu.dk/services/SignalP-5.0/

TMHMM - transmembrane domain https://dtu.biolib.com/DeepTMHMM


      | Dockstore : `func_anno_main.wdl <https://dockstore.org/workflows/github.com/MethodsDev/IsoFuncAnnot/func_anno_main>`_
      | Github: `func_anno_main <https://github.com/MethodsDev/IsoFuncAnnot/blob/main/anno_main.wdl>`_


.. code:: bash
  :number-lines: 

  {
    "anno_main.inputAAfasta":"${this.inputAAfasta}",
    "anno_main.inputNTfasta":"${this.inputNTfasta}",
    "anno_main.pfam.pfamArgs":" -e_seq 10.0 -e_dom 10.0"
  }
    
Output :

.. code:: bash
  :number-lines: 

  {
  "anno_main.cpc2Out":"${this.cpc2Out}",
  "anno_main.iupredErrlogfile":"${this.iupredErrlogfile}",
  "anno_main.iupredOut":"${this.iupredOut}",
  "anno_main.pfamOut":"${this.pfamOut}",
  "anno_main.signalPOut":"${this.signalPOut}",
  "anno_main.tmhmmOut":"${this.tmhmmOut}"
  }

