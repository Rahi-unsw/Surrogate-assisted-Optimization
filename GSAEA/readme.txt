%%%%%%%%% Algorithm Name: GSAEA (Generalized Surrogate-Assisted Evolutionary Algorithm) %%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Algorithm Name: SASSEA (Surrogate-Assisted Steady-State Evolutionary Algorithm) %%%%%%%%%%%%%%%%%%%%

Please cite the following papers if you use this algorithm.

%------------------------------- Copyright ------------------------------------
% Copyright (c) MDO Group, UNSW, Australia. You are free to use the GSAEA and SASSEA for
% research purposes. All publications which use this code should acknowledge the 
% use of "SASSEA", "GSAEA" and references: 
% "K.H.Rahi, H.K.Singh, T.Ray, A steady-state algorithm for solving expensive multi-objective 
% optimization problems with non-parallelizable evaluations, IEEE Transactions on Evolutionary Computation, 2022" and
% "K.H.Rahi, H.K.Singh, T.Ray, A generalized surrogate-assisted evolutionary
% algorithm for expensive multi-objective Optimization, IEEE Congress on Evolutionary Computation, 2023 (Accepted for publication)".
%------------------------------------------------------------------------------ 

%%%%% Steps to run the code %%%%%%%%%%%%

1. Go to Params.m script. You can specify all the parameters. Please, write the problem scripts (follow TR1.m, TR2.m files etc.). All the parameters are kept
   in the 'param' structure.
2. Go to the ND_Sort folder and in the matlab command prompt, write mex('E_NDSort.c') to built the mex file successfully.
3. Run the Multirun.m script now. If you don't want to fire parallel runs, then comment line 62, 63 and 71 and write 'parfor' instead of 'for'.
4. All the data will be saved in a newly generated 'Data' folder. The format to store data is arranged as follows:
   Data->Variant-(Variant number)>Problem name->Trial number->Stored data (ProblemName_Archive.m, Parameters.m).

 