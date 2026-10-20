"""Independent-current transport and induced-current loss models."""

from .joule import JouleResult, joule_energy
from .eddy import EddyResult, eddy_energy, eddy_energy_slab_modes
from .network import CoupledLossResult, coupled_network_loss

__all__ = ["JouleResult", "joule_energy", "EddyResult", "eddy_energy",
           "eddy_energy_slab_modes", "CoupledLossResult", "coupled_network_loss"]
