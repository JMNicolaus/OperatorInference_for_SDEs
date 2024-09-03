# OperatorInference_for_SDEs
This repository contains the MATLAB(2022a) implementation of operator inference (OpInf) approach for stochastic differential equations (SDEs) as described in:

[1] M. A. Freitag, J. M. Nicolaus, M. Redmann
 	[Learning Stochastic Reduced Models from Data: A Nonintrusive Approach](https://doi.org/10.48550/arXiv.2407.05724)<details><summary>BibTex</summary><pre>
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
The Subfigures of the Figures 5.2 and 5.3 are produced by 'genPlots.m'.
The Figure 5.1 and Table 5.1 are produced by 'genPlotSubspace.m'.


The 'ConvectionReaction' example uses the matrices provided by the file '[pde.mat](https://www.slicot.org/objects/software/shared/bench-data/pde.zip)', which is part of the benchmark examples of
[SLICOT - Subroutine Library in Systems and Control Theory](https://www.slicot.org/20-site/126-benchmark-examples-for-model-reduction).<details><summary>BibTex</summary><pre>
@MANUAL{slicot_pde,
 title =        {{SLICOT} - Subroutine Library in Systems and Control Theory},
 organization = {Niconet e.V.},
 address =      {\url{http://www.slicot.org}},
 key =          {SLICOT}
}</pre></details>
