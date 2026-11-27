# Smith reproduction code

All Smith main scripts, helpers and method notes are consolidated here.
Open main scripts in RStudio and execute numbered sections interactively.
Paths and settings are at the top; no command-line arguments are required.

| Main script | Role | Output location under Results/ |
| --- | --- | --- |
| ERA5_Smith_PNW_GEV_MLE.R | Annual GEV with/without 2021; Table 0.1 | ERA5_Smith_PNW_GEV_MLE/ |
| ERA5_Smith_PNW_GEV_MCMC.R | GEV posterior risks/endpoints; Section 0.6.1 | ERA5_Smith_PNW_GEV_MCMC/ |
| ERA5_Smith_PNW_POT_MLE.R | Seasonal PP MLE; Section 0.6.2 | ERA5_Smith_PNW_POT_MLE/ |
| ERA5_Smith_PNW_POT_MCMC.R | Seasonal PP posterior; Table 0.2 | ERA5_Smith_PNW_POT_MCMC/ |

Shared helpers are in helpers/. NOTES.md files explain method choices and limits.
All source paths use Code/smith_reproduce. Inputs and saved results/checkpoints
were not relocated, preserving existing MLE hashes and PP resume compatibility.
GEV MCMC needs the new GEV MLE snapshot; PP MCMC uses its original PP snapshot.
Do not source an MLE main script just to load functions; source helpers only.

For texmex 2.4.9, attach library(texmex) before PP declustering if evmReal is not
found: that package dynamically resolves the function. The move did not alter
PP model or sampling algorithms. Other non-Smith analysis scripts are unchanged.

GEV MLE and its four-chain pilot ran on 2026-11-18; see WORKFLOW.md for review.
Existing PP MLE/pilot outputs
retain their review limitations. Long-chain settings are not launched automatically.

