# Large-scale genomic surveillance reveals immunosuppression drives mutation dynamics in persistent SARS-CoV-2 infections

Code for "Large-scale genomic surveillance reveals immunosuppression drives mutation dynamics in persistent SARS-CoV-2 infections"

Link to paper: To be added

## Overview

The presented code aims to reproduce the analyses from "Large-scale genomic surveillance reveals immunosuppression drives mutation dynamics in persistent SARS-CoV-2 infections". In brief, the goal of the code is to identify persistent infections, identify risk factors for persistent infections, analyze associated mutational patterns, and link health conditions to accelerated selective genetic changes in the virus.

## Folders

* ```case_control_code```: Code used for preparing data and running the case control study.
* ```defining_snp_thresholds```: Contains code to identify the false-positive rates at different rare SNP thresholds.
* ```dnds```: Code to conduct dnds/kaks analyses.
* ```figures```: Contains code that leads to the four main figures in the article.
* ```identifying_persistent_infections```: Code used to identify persistent infections from the full set of sequences.
* ```identifying_substitutions_and_mutations```: Code used to identify the placement of mutations and non-synonymous changes, as well as analyzing the association between different risk factors and diagnoses with the rate of (non)synonymous change.
* ```identical_sequences```: Code to analyze the distance between the final persistent infection sequences and all other sequences in the database within 11 days.
* ```rebounding```: Code to analyze the occurrence of rebounding among persistent infections.
* ```recurrent_mutations```: Code to identify independently-occurring and recurrent mutations among persistent infections.
* ```tip_length_analysis```: Code to analyze tip lengths for persistent versus control infections.

All code is numbered in the intended order of use.


## General Data

Due to the sensitive nature of the data, we cannot provide individual-level metadata or further metadata for the sequences.


## Software overview

R (4.2.3), Python (3.10.10), brms (2.20.3), Guppy, Medaka, iVar (1.4.3), BCFtools (1.18), VADR, Trim Galore (0.6.10).

## Authors

* Mark Khurana (<mark.khurana@sund.ku.dk>)
* Alexandros Katsiferis (<alexandros.katsiferis@sund.ku.dk>)
* Neil Scheidwasser (<neil.clow@sund.ku.dk>)


## License

Apache 2.0 License

## Citation

To be added.
