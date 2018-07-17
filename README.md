# Pandora &nu;<sub>e</sub> selection
This module is a [LArSoft](http://www.larsoft.org) analyzer that builds a ROOT `TTree` with information from &nu;<sub>e</sub> candidates, reconstructed with the [Pandora framework](https://github.com/PandoraPFA).
It looks for neutrino PFParticles with at least one daughter shower and at least one daughter track, or at least two daughter showers in the full Pandora PFParticle hierarchy.

There are three FCL files, two for data (`run_PandoraOnly_data_bnb.fcl`, `run_PandoraOnly_data_extbnb.fcl`) and one for Monte Carlo (`run_PandoraOnly.fcl`).

## Requirements
- `uboonecode v06_26_01_17` (MCC8.7)
- Checkout `uboonecode` feature `feature/alister1_EventWeightTreeUtility`
- Checkout `uboonecode` feature `wvdp_lightcharge_v6_26`
