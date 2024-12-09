
To-do Tasks/ Feedback
=====================

`All available emoticons <https://sphinxemojicodes.readthedocs.io/en/stable/>`_

|:bulb:| - discussion

|:green_circle:| - completed

|:salad:| - done but review before publishing

|:o:| - open/currently working

|:stop_sign:| - blocked


Outstanding tasks as of Dec 9th:
---------------------------------
1. single cell vignette update static images with Full dataset
2. add in pseudo-bulk vignette (RMD plus html for static viewing)
3. isoformSwitch for single cell
4. add static plots in ipynb for ITV below respective code blocks
5. revise toctree structure 

Nov 20th:
----------
1. Integrate ITV - git lfs
2. isoformSwitchAnalysisR - for bulk debug pfam
3. isoformSwitchAnalyzeR - for single-cell - pseudo-bulk approach


Sep4th:
-------
1. Can we add/Should add toctree in Index left side panel? Yes we can we'd like. |:salad:|
2. Use seperate test datasets for each of the 3 sub-workflows pre, sec and tert. |:salad:|
3. Add info on bulk multiplexing - preferrably schametic |:salad:|

Sep 25th:
---------
4. Use PacBio dataset for single cell analysis, from cleaned s-reads but push it through the R-vignette we have |:salad:|
5. for bulk use cell lines from tech transfer  - start from Isoquant counts - ref guided |:salad:|
6. Ok to use google collab - find out approx. mem usage and $ required to run our vigenettes - using personal collab account |:bulb:| 
7. updated version of MDL genomeview ITV - add static page |:stop_sign:|
8. Add RNA seq QC+ page |:salad:|

Oct 24th:
---------
9. remove secondary processing or add minimap2 there and a toctree for ipynbs |:bulb:| 
10. Wet-lab best practices - revise and complete the page |:stop_sign:| 
11. Clearly mark vigenttes vs wdls |:salad:|
12. swap out MDL genomeview to ITV in workflows and other |:salad:|
13. can we have overarching workflow for all secondary analysis steps? merge- minimap2 - LongRNAseqQC? |:bulb:|
14. Reference to other wiki’s that our not a vignette
15. Add all links referrences glossary at the end

Oct25th:
--------

.. figure:: ../_images/test_datasets.png
   :alt: Novel Methods and R&D
   :align: left



Nov15th:
--------
16. Add the script to combine isoquant tanscript and gene counts |:o
17. long reads vigentte - debug collab env
18. literature review for isoform switching tools publichsed for single cell 
19. isoformSwitchAnlysisR for single cell



Other possibilities for Full Length:
FLNC - Full Length Non-concatemer reads - post-refine available for below datasets


.. figure:: ../_images/pb_bulk_datasets.png
   :alt: Novel Methods and R&D
   :align: left

HG002
~~~~~~

HG002 cells were purchased from Coriell and grown in RPMI1640 using Glutamax 
media with 16% FBS and 0.5% Penicillin-Streptomycin. RNA was isolated from 
10x10^6 HG002 cells using Trizol reagent and Phasemaker tubes. RNA quality was 
assessed using a Bioanalyzer.

UHRR (Universal Human Reference RNA)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Vendor – Agilent
Part No - 740000
UHRR total RNA was purchased from Agilent and directly used for cDNA generation.

Heart and Cerebellum samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Collaborator – Seattle Children’s Research Institute (SCRI)

Two heart samples from prenatal specimens with trisomy 21 were obtained from the 
Birth Defects Research Laboratory tissue repository. Total RNA was isolated from 
50 mg of fresh frozen tissue using the Promega Maxwell Kit.

Two heart samples from prenatal specimens were obtained from the Birth Defects 
Research Laboratory tissue repository. Total RNA was isolated from 50 mg of fresh 
frozen tissue using the Promega Maxwell RNA Extraction Kit.

Three cerebellum samples were obtained from the Birth Defects Research Laboratory 
tissue repository. Total RNA was isolated from fresh frozen tissue sections or 
following laser capture microdissection using the Qiagen RNeasy Micro Kit as 
described in PMID:34140698.

The collaborator has granted permission to release this dataset.
