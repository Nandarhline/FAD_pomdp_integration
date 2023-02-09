FILE              OUTPUT VARIABLES              OUTPUT FILE NAME            REMARK
SIF_FAD           MCS crack samples, pF         Deterioration.mat, pf_FAD   Change FM model / failure criteria here 
dbn_dr_trans      dr_env                        dr_ENV.mat                  Play with discretization here 
dbn_dr_uncond                                                               To compare pF between MCS and DBN 
MCS_updat         pF_MCS, pf_MCS, beta_MCS      mcs_upd_out.mat             Change POD model here, need to specify inspection years 
dbn_dr_cond       dr_env (with POD model)       dr_ENV.mat                  Input the same POD model and inspection years, to compare updated pF between MCS and DBN 
parsing_infinite  dr_env (with rewards)         IFH_input.mat               Change the cost model here 
parsing_finite    dr_env (FH models)            FN_input.mat                Transforming state, transition, observation, reward models  
headingSarsop                                   Name.pomdp                  Writing the pomdp input file 

## Citation
If you use this DBN-POMDP implementation for complex deterioration models and failure criteira in your own work, please cite our paper:

```
@article{Hlaing2022,
author = {Nandar Hlaing and Pablo G. Morato and Jannie S. Nielsen and Peyman Amirafshari and Athanasios Kolios and Philippe Rigo},
title = {{Inspection and maintenance planning for offshore wind structural components: integrating fatigue failure criteria with Bayesian networks and Markov decision processes}},
journal = {Structure and Infrastructure Engineering},
volume = {18},
number = {7},
pages = {983-1001},
year  = {2022},
publisher = {Taylor & Francis},
doi = {10.1080/15732479.2022.2037667},
}

```
