Tertiary Processing
===================

Sub-heading1
---------------

Option1 for sub1:
~~~~~~~~~~~~~~~~

Option2 for sub1:
~~~~~~~~~~~~~~~~~

High-level Workflow:
--------------------
   
   - Performing Sample QC and PCA with `Ptr` utility from `Trinity toolkit <https://github.com/trinityrnaseq/trinityrnaseq/wiki>`_.
   - `minimap2 <https://lh3.github.io/minimap2/minimap2.html>`_ for aligning reads
   - `Isoquant <https://github.com/ablab/IsoQuant>`_ for Reference guided Isoform Discovery.
   - `Stringtie <https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual>`_ merge to merge reconstructed assemblies to generate new Ref
   - Quantification with `Isoquant <https://github.com/ablab/IsoQuant>`_ against reconstructed ref
   - `Gffcompare <https://github.com/gpertea/gffcompare>`_ to compare against reference and fetch Ref Ids
   - Differential Expression analysis using DEseq2 implementation from the Trinity
   - Exploring isoform switches and isoform visualization with isoformSwitchAnalysisR
   - MDLgenomeview for inspecting read support