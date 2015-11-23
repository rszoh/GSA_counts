#!/bin/bash

module load R-latest
/share_home/rszoh/test_R_script/GSA_counts/Test_based_on_clustering_Power_analysis_wrt_Sig.R $SGE_TASK_ID
module unload R-latest