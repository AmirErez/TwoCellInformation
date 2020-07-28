# TwoCellInformation
Code and data for the journal article "Cell-to-cell information at a feedback-induced bifurcation point"

Phys. Rev. Lett. 125, 048103 – Published 22 July 2020

https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.048103


* Data - Data in tabular form, summarized using Collect.m from Gillespie simulation results.
  Read Gillespie simulations .mat files are too large to store here, but were generated using
  exec_scan_thetax_thetay_h0.cmd and similar SLURM scripts.
 
* Gillespie - Code for Gillespie simulations. Simulations were written in MATLAB and compiled to mex.
  The compilation scripts are SimulateSchlogl2cell_codegenscript.m and SimulateHill2cell_codegenscript.m

* Plot - Code for plotting the figures
