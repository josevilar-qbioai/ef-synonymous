"""EF-Synonymous — glass-box thermodynamic scoring of synonymous variants (σ + G>A).

Motor self-contained (numpy-only). Comparte exactamente la lógica del prototipo
web y de core/synonymous_sigma_predictor.py, verificada a 1e-6.
"""
__version__ = "0.2.0"

from .engine import sigma_signed, is_g_to_a, stacking_profile, synonymy, STACKING_TURNER
from .predictor import SynonymousSigmaPredictor

__all__ = ["SynonymousSigmaPredictor", "sigma_signed", "is_g_to_a",
           "stacking_profile", "synonymy", "STACKING_TURNER"]
