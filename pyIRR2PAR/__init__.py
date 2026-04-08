from .IRR2PAR import IRR2PAR
from .JPitarch_Ed_380_443_490_555 import JPitarch_PAR_from_Ed
from .JTan2025_Ed_380_443_490_555 import JTan_PAR_from_Ed

from importlib.metadata import version as _version

try:
    __version__ = _version("pyIRR2PAR")
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "9999"

__all__ = ("IRR2PAR", "JPitarch_PAR_from_Ed", "JTan_PAR_from_Ed")
