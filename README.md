# *_Comparing metagenomics and total RNA-Seq pipelines for microbial diversity assessments_*
This project represents the first chapter of my PhD, in which I compare metagenomics and total RNA-Seq data-processing pipelines to identify which sequencing technique and which pipeline give the most accurate estimate of the diversity of a microbial mock community.
The repo contains scripts to run pipelines on compute canada clusters or locally, as well as scripts for statistics.

The workflow from sample preparation to accuracy evaluation is shown in the following:
<img src="https://github.com/hempelc/metagenomics-vs-totalRNASeq/blob/master/workflow.png" alt="workflow" width="400"/>

The processed pipelines are combinations of typically utilized data-processing tools, and the following figure gives an overview of the tested tools per processing step. ALl combinations of tools were tested, resulting in 3,064 tested pipelines.
<img src="https://github.com/hempelc/metagenomics-vs-totalRNASeq/blob/master/pipeline_steps.png" alt="pipelines" width="400"/>

The results are to be published in a paper.
