2. Full Length Bulk Secondary Processing
=========================================

After performing QC on the aligned reads using LongRNAseqQC, 
we can proceed to perform transcript analysis using an array of tools put togther by folks using Long read sequencing.

There are multiple applications one can use for many of these sub-steps, below is the currently recommended workflow.

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
The test dataset used for demonstrating downstream processing workflows for Full Length correspond to HG002 from the `Genome-in-a-Bottle project <https://www.nist.gov/programs-projects/genome-bottle>`_ 
and K562 cell lines. 


2.1 `minimap2`
~~~~~~~~~~~~~~~~~
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