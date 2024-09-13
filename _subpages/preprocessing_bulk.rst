Analysis workflows
++++++++++++++++++
The analysis workflows for Kinnex products are similar for pre-processing with a few tweaks in parameters, and diverge there after. 
The preliminary processing mainly leverages various applications from PacBio Analysis Toolkit. 
Most workflows have WDL wrappers publicly accessible via Dockstore with test data and reference files to follow along placed in public requester-pays Google Storage bucket.

a. Kinnex Full Length 
=====================

Overall Workflow in a nutshell
------------------------------

.. figure:: ../_images/bulk_workflow.png
   :alt: Novel Methods and R&D
   :align: left

Preliminary analysis
~~~~~~~~~~~~~~~~~~~~
The pre-processing workflows extract clean s-reads using 3 tools as below which can then be provided to the alignment applications and other downstream workflows similar to those used to analyze Isoseq data.

   - `skera <https://skera.how/>`_ for de-concatenating the MAS arrays into individual cDNA molecules and generate segmented reads (s-reads),
   - `lima <https://lima.how/>`_ to de-multiplex a bulk pool and remove unwanted combinations while orienting sequences in 5’ → 3’ orientation and 
   - `isoseq refine <https://isoseq.how/getting-started.html>`_ for trimming poly(A) tails and extracting Full length non-concatemer reads (FLNC) from s-reads.

The pbskera workflow processes raw HiFi reads generated with Revio LR seqeuncers.(for previous versions, kinetic info is presumed to be removed). The HiFi reads are a current default, and can be plugged in directly into the workflow. The direct command executed here is:

Downstream analysis
~~~~~~~~~~~~~~~~~~~
   
   - Performing Sample QC and PCA with Ptr utility from Trinity toolkit.
   - Minimap2 for aligning reads
   - Isoquant for Reference guided Isoform Discovery.
   - Stringtie merge to merge reconstructed assemblies to generate new Ref
   - Quantification with Isoquant against reconstructed ref
   - Gffcompare to compare against reference and fetch Ref Ids
   - Differential Expression analysis using DEseq2 implementation from the Trinity
   - Exploring isoform switches and isoform visualization with isoformSwitchAnalysisR
   - MDLgenomeview for inspecting read support


b. Kinnex Single Cell
=====================