Description of the files: 
----------------------------------------------------------------------
Linear inversion to get initial samples: 
  1. state0.priors.m : Main script to run the linear inversions
  2. makemesh_full_inv.m : making triangular mesh 
  3. grn_func.m = calculates greens function with triangular mesh as input
  4. prior_samples.m = function that does linear inversion 
  5. lin_inv.m : inside prior_samples script to compute inversion
----------------------------------------------------------------------

----------------------------------------------------------------------
The ATMIP algorithm is divided into 3 parts of job: 
1. Initializes with the prior samples at stage 0
  a) firsthalf_run.m : Main script to run 
  b) ATMIP_remaining1.m : runs under the main script 
2. Continues at the current stage (without changing stage)
  a) nexthalf_run.m : Main script to run
  b) ATMIP_remaining2.m : runs under the main script 
3. Changes the stage: 
  a) change_stagerun.m : Main script to run 
  b) ATMIP_change.m : runs under the main script
----------------------------------------------------------------------

----------------------------------------------------------------------
Other scripts : 
1. AMH.m : Adaptive Metropolis Hastings sampling 
2. CalcTriDisps.m : used in grn_fnc to calculate displacement due to triangular dislocations
3. deterministicR.m : used in ATMIP_remaining during bootstrapping
4. example_code.m : contains scripts to run in the cluster
5. findpoint.m : used in mesh generation 
6. grn_func.m : calculates greens function with triangular mesh as input
7. laplacian.m : used with grn_func for slip smoothness 
8. likelinew.m : calculates likelihood function 
9. posteriorTohoku.m : calculates the posterior 
10. priorTohoku.m : calculates the geometrical prior
11. reversefindpoint.m : used in mesh generation
12. slipprior.m : calculates slip smoothness prior 
----------------------------------------------------------------------

