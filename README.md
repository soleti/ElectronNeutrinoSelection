[![Build Status](https://travis-ci.org/soleti/ElectronNeutrinoSelection.svg?branch=master)](https://travis-ci.org/soleti/ElectronNeutrinoSelection)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/16286/badge.svg)](https://scan.coverity.com/projects/soleti-electronneutrinoselection)

# WORK IN PROGRESS - DOES NOT HAVE ALL THE FEATURES OF THE MCC8 VERSION

# Pandora &nu;<sub>e</sub> selection

### [Doxygen Documentation](https://soleti.github.io/ElectronNeutrinoSelection/html/annotated.html)

This module is a [LArSoft](http://www.larsoft.org) analyzer that builds a ROOT `TTree` with information from &nu;<sub>e</sub> candidates, reconstructed with the [Pandora framework](https://github.com/PandoraPFA).
It looks for neutrino PFParticles with at least one daughter shower and at least one daughter track, or at least two daughter showers in the full Pandora PFParticle hierarchy.

There are three FCL files, two for data (`run_PandoraOnly_data_bnb.fcl`, `run_PandoraOnly_data_extbnb.fcl`) and one for Monte Carlo (`run_PandoraOnly.fcl`).

## Requirements
- `uboonecode v07_05_00_02` (MCC9)

# License

See the [LICENSE](LICENSE) file for license rights and limitations (MIT).
