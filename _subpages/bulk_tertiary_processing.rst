Tertiary Processing
===================

After performing QC on the aligned reads using LongRNAseqQC, 
we can proceed to perform transcript analysis using an array of tools put togther by folks using Long read sequencing.

There are multiple applications one can use for many of these sub-steps, below is the currenty recommended workflow.

High-level Workflow:
--------------------
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


`minimap2`
~~~~~~~~~~


`Isoquant - Isoform Discovery`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Stringtie`
~~~~~~~~~~~

`Isoquant - Quantification`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Differential Expression with DEseq2`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


`isoformSwitchAnalysisR`
~~~~~~~~~~~~~~~~~~~~~~~~~

