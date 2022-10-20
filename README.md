# *_Comparing metagenomics and total RNA-Seq pipelines for microbial diversity assessments_*
This project represents the first chapter of my PhD, in which I compare metagenomics and total RNA-Seq data-processing pipelines to identify which sequencing technique and which pipeline give the most accurate estimate of the diversity of a microbial mock community.
The repo contains scripts to run all pipelines and scripts for statistics.

The project is part of the following publication:<br>
<b>Hempel C. A.</b>, Wright N., Harvie J., Hleap J. S., Adamowicz S. J., and Steinke D. (2022): Metagenomics versus total RNA sequencing: most accurate data-processing tools, microbial identification accuracy and perspectives for ecological assessments. Nucleic Acids Research 50(16):9279–9293. [https://doi.org/10.1093/nar/gkac689](https://doi.org/10.1093/nar/gkac689)

The workflow from sample preparation to accuracy evaluation is shown in the following:

<img src="https://github.com/hempelc/metagenomics-vs-totalRNASeq/blob/master/workflow.png" alt="workflow" width="250"/>

The processed pipelines are combinations of typically utilized data-processing tools, and the following figure gives an overview of the tested tools per processing step. We tested all combinations of tools (1,536 combinations in total).

<img src="https://github.com/hempelc/metagenomics-vs-totalRNASeq/blob/master/pipeline_steps.png" alt="pipelines" width="800"/>
