.. kinnex-docs-external documentation master file, created by
   sphinx-quickstart on Fri Sep  6 01:29:45 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html
  :file: _images/MDL_RTD.svg

Resource for long-read RNA isoform sequencing analysis
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This website serves as an updated resource for bulk and single-cell RNA isoform sequencing analysis from Kinnex data. 
Using latest methods and best practices, we have compiled a series of workflows and notebook vignettes to facilitate data processing and analysis. 
We hope this resource will serve as a useful guide for this new and feature rich data type.

The Kinnex protocol is based on the MAS-ISO-seq method, developed at the Broad. For information on the approach, please reference the original publication in Nature Biotechnology:
   `High-throughput RNA isoform sequencing using programmed cDNA concatenation <https://doi.org/10.1038/s41587-023-01815-7>`_

From receiving sequencing data from `PacBio's Revio <https://www.pacb.com/revio/>`_ sequencing platform, 
the document steps through various pre-processing workflows for obtaining cleaned s-reads suitable for downstream processing and tertiary analysis workflows and vignettes developed by the `MDL team <https://methodsdevlab.org/>`_ at `Broad Clinical Labs <https://broadclinicallabs.org/>`_ 
to explore Kinnex datatypes.

`Repository of Public Datasets for Kinnex, Sequel2e and previous MAS-ISO-seq versions <https://downloads.pacbcloud.com/public/dataset/>`_

.. grid:: 2 

    .. grid-item-card::  Kinnex Full Length 
      :link: _subpages/bulk_landing_page
      :link-type: doc

      ``Overview and workflow`` 

    .. grid-item-card::  Kinnex Single Cell
      :link: _subpages/sc_landing_page
      :link-type: doc

      ``Overview and workflow``



.. toctree::
   :maxdepth: 2
   :caption: Index : 

   _subpages/pb_guidelines


.. toctree::
   :maxdepth: 4
   :caption: Kinnex Full Length :

   _subpages/bulk_landing_page
   _subpages/preprocessing_bulk
   _subpages/bulk_secondary_processing
   _subpages/readQC
   jupyter-notebooks/demo-ITV
   _subpages/bulk_tertiary_processing
   jupyter-notebooks/id_mapping_from_gffcompare
   _subpages/isoformSwitch_bulk


.. toctree::
   :maxdepth: 4
   :caption: Kinnex Single Cell :

   _subpages/sc_landing_page
   _subpages/preprocessing_sc
   _subpages/secondary_processing
   _subpages/readQC
   jupyter-notebooks/demo-ITV
   _subpages/single_cell_tertiary_processing_500cells
   _subpages/sc_subsampling.rst


.. toctree::
   :maxdepth: 1
   :caption: Citations :

   _subpages/citations.rst
   _subpages/todo.rst


Doc set-up Notes
=================
To be removed during cleanup!!

1. Link to best practices for Kinnex

   a. Guidelines - PacBio SMRT link v13 guide for Kinnex product documentation - Troubleshooting guide - Wetlab perspective

2. Processing the reads - Divide by Product type

   a. PacBio documentation 
   b. Our workflows with Dockstore links
   c. Cleaned S-reads → downstream bioinformatic analyses 

3. Read-level Quality Controls : 
   
   a. Read level SQANTI with RNAQC+

4. Analysis Vignettes for Kinnex Full Length product type: 

   a. Aligning to the genome - minimap2
   b. Isoform ID with IsoQuant 
   c. isoformQC
   d. Isoform Quantification
   e. Isoform Annotation - gffcompare, mdl- functional annotations attr vignette
   f. Visualization - IGV, Genomeview
   g. Differential Expression DEseq2
   h. isoformSwitchAnlysisR

5. Analysis Vignettes for Kinnex Single Cell Isoform analysis



Indices and tables
==================

* :ref:`genindex`
* :ref:`search`


Foot Note
==========

This is a dev-test version of rtd-website to house MDL documentation. First we are trying to create the scanpy theme.
The theme for scanpy is sphinx-book-theme, with patches for readthedocs-sphinx-search

Next, we'll create the layout for the documentation and all the subsequent .rst pages. 
Once done with the the landing index page and we'll add/edit the relevant information to corresponding pages.

