# GNRS (Graph nonreference sequences)
Pangenome of human nonreference sequences from population-scale long-read sequencing   

## Description 
In order to make sure the results are reproduceable, the pipeline is performed using framework [**Snakemake**](https://snakemake.readthedocs.io/en/stable/) coupled with the environment conducted by [**Anoconda**](https://www.anaconda.com/). And the pipeline can be used in other cohort with long-read sequencing.

The workflow of GNRS on the population-scale long-read sequencing are below:
![image](https://user-images.githubusercontent.com/42490165/187616551-c578eb18-95d2-4a84-82c1-b0e7290c5fd5.png)

### Schematic representation of GraphNRS
* a, Long-read sequencing data from different platforms are de novo assembled and polished. 
* b, The NRSs are anchored to GRCh38. Placed NRSs are clustered to select the representative NRSs, and unplaced NRSs are clustered after filtering out contaminants and centromeric repeats. Then, we merge the placed and the unplaced NRSs to obtain the nonredundant NRSs of the whole population. 
* c, vg is used to construct the graph pangenome, and NRS genotyping is performed for each NRS of the individual.


## Requirements


## Configure the environment


## Quick start for the pipeline


## Introduction of pipeline



