![logo](./images/XSVLEN.PNG)




## Introduction

Benchmarking SVs is not always easy so here we provide a new method for verifying the presence of an SV using haplotype-resolved assemblies.


# What is XSVLen?

XSVLen takes cuteSV or Sniffles VCF and using reference coordinates will produce reference sequences having included inserted sequences or deleted sequences. By creating this variant sequences we could check for the presence of the predicted variants in haplotype-resolved assemblies.

The variant sequences will be mapped into the assemblies using minimap2 and length of alignment will be compared to variant queries length. 

A report summary will be produced using an R script. 

# Workflow diagram

![logo](./images/Workflow.png)

# How to use XSVLen
  
A Nextflow workflow is produced

