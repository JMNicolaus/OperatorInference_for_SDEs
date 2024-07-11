# OperatorInference_for_SDEs
This repository contains the MATLAB implementation of operator inference (OpInf) approach for stochastic differential equations (SDEs) as described in:

[1] M. A. Freitag, J. M. Nicolaus, M. Redmann
 	[Nonintrusive model order reduction for stochastic differential equations](https://doi.org/10.48550/arXiv.2407.05724)<details><summary>BibTex</summary><pre>
@misc{freitag2024nonintrusivemodelorderreduction,
      title={Nonintrusive model order reduction for stochastic differential equations}, 
      author={M. A. Freitag and J. M. Nicolaus and M. Redmann},
      year={2024},
      eprint={2407.05724},
      archivePrefix={arXiv},
      primaryClass={math.NA},
      url={https://arxiv.org/abs/2407.05724}, 
}</pre></details>

The scripts 'genPlots.m' and 'genPlotSubsspace.m' produce the figures and the table in [1].
To obtain these figures, one needs to specify which example and which snapshot method to use. 
The predefined examples are available by setting the variable FOM.eqtype to 'Heat', '2dHeat' or 'ConvectionReaction'.
Specifying the variable 'snapshotType' to be 'moment' or 'state' results in the moment-snapshot matrix or the state-snapshot matrix to be used, respectively. 
The parameters 'L','s' and 'h' correspond to the number of sampled trajectories, the number of time-steps and their size. 
