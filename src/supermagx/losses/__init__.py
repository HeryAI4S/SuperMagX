"""Independent-current transport and induced-current loss models."""

from .joule import JouleResult, joule_energy
from .eddy import EddyResult, eddy_energy, eddy_energy_slab_modes
from .network import (CoupledLossOperator, CoupledLossResult,
                      LowFrequencyLossCoefficients, compile_coupled_loss,
                      coupled_network_loss)

__all__ = ["JouleResult", "joule_energy", "EddyResult", "eddy_energy",
           "eddy_energy_slab_modes", "CoupledLossOperator", "CoupledLossResult",
           "LowFrequencyLossCoefficients", "compile_coupled_loss",
           "coupled_network_loss"]
