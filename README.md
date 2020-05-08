# Algorithms in Computational Biology (CSE 494/598) - Spring 2020 - Project
## Proposal

To generate a pdf from the tex file:

`make CSE494_598_project_proposal.pdf`

To remove the intermediate files created by generating the pdf:

`make clean`

## Simulations

To run a simulation (may take very long the first time you run this):

`bash ./scripts/run_simulation.sh n` (where n is the number of desired sample contributors)

To tally all simulation results:

`bash ./scripts/hla_mtdna_tally.sh`

Results may be found in `simulations/`.
