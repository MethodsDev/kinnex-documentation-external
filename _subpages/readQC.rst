Read-level Quality Controls
+++++++++++++++++++++++++++

LongRNAqc+ is a workflow intended to look at base reads quality rather than collapsed reads, and to easily compare multiple samples on Terra. 
Taking an aligned BAM as input for a sample, it runs a long-read adapted version of RNAqc<>, `SQANTI3 <https://www.nature.com/articles/s41592-024-02229-2>`_, and optionally IsoQuant. 

Once multiple samples have been processed, a downstream workflow can be run on a selection of samples to generate a comparison report based on the Sqanti3 outputs.

Main workflow: `LongRNAqcPlusFromBAM`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Workflow configuration for runnning the main workflow over cloud platforms supporting Cromwell like Terra can be found here:-

      | Dockstore : `LongRNAqcPlusFromBAM.wdl <https://dockstore.org/workflows/github.com/broadinstitute/MDL-workflows/LongRNAqcPlusFromBAM:main>`_
      | Github: `LongRNAqcPlusFromBAM <https://github.com/broadinstitute/MDL-workflows/blob/main/LR-tools/LongRNAqc/LongRNAqcPlusFromBAM.wdl>`_
      | Test Data can be found here (public, requester-pays) : (add download link from github or gs here)

**Input arguments for LongRNAqcPlusFromBAM**

.. csv-table:: `LongRNAqcPlusFromBAM`
   :file: ../_subpages/tables/longRNAqcmain.csv
   :widths: 20,25,55
   :header-rows: 1


`LongRNAqc+ Plotting`
~~~~~~~~~~~~~~~~~~~~~
Workflow configuration for runnning LongRNAqc+ on Terra can be found here:-

      | Dockstore : `LongRNAqcPlotting:main <https://dockstore.org/workflows/github.com/broadinstitute/MDL-workflows/LongRNAqcPlotting:main>`_
      | Github: `LongRNAqcPlotting.wdl <https://github.com/broadinstitute/MDL-workflows/blob/main/LR-tools/LongRNAqc/LongRNAqcPlotting.wdl>`_
      | Test Data can be found here (public, requester-pays) : `add file path` 


**Example of input arguments for LongRNAqc+ workflow for alignment with human ref genome**

.. code:: bash
  :number-lines: 
  
  {
    "LongRNAqcPlotting.classificationFile":"${this.sqantiClassificationTSV}",
    "LongRNAqcPlotting.junctionFile":"${this.sqantiJunctionsTSV}",
    "LongRNAqcPlotting.includeSaturation":"${false}",
    "LongRNAqcPlotting.outputPrefix":"LongRNAqcPlots",
    "LongRNAqcPlotting.sampleName":"${this.samples.sample_id}"
    }



**Example of plots generated as a part of the output_report.pdf populated in `QC_plots` field**

.. list-table:: 
    :widths: 50 50

    * - .. figure:: ../_images/longRNAqc.1.png
           :alt: Kinnex Full Length - FSM

           Normalized transcript length distriubtions - FS

      - .. figure:: ../_images/longRNAqc.2.png
           :alt: Kinnex Full Length - NIC

           Normalized transcript length distriubtions - NIC 


.. list-table:: 
    :widths: 50 50

    * - .. figure:: ../_images/longRNAqc.5.png
           :alt: Kinnex Full Length - Category wise

           Isoform Distribution across SQANTI Structural Categories

      - .. figure:: ../_images/longRNAqc.4.png
           :alt: Kinnex Full Length - FSM

           Isoform Distribution Across FSM's and Samples


.. figure:: ../_images/longRNAqc.3.png
   :scale: 45%
   :align: left