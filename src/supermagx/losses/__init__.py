"""Independent-current transport and induced-current loss models."""

from .joule import JouleResult, joule_energy
from .eddy import EddyResult, eddy_energy

__all__ = ["JouleResult", "joule_energy", "EddyResult", "eddy_energy"]
