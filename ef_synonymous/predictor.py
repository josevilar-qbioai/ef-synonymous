"""
Predictor glass-box de patogenicidad de variantes SINÓNIMAS (σ + G>A).

Regresión logística de 2 features físicas que reproduce a la CNN dedicada del
Paper 1 (AUC 0.709, Turner RNA; a la par con la CNN en el set gene-matched, 0.659). Inferencia solo-numpy. El método `score` devuelve el dossier
completo, auditable y trazable variante a variante (argumento SaMD/ACMG).
"""
import json
import os

import numpy as np

from .engine import sigma_signed, is_g_to_a, synonymy

_MODEL_PATH = os.path.join(os.path.dirname(__file__), "model.json")


def _pkg_version():
    try:
        from . import __version__
        return __version__
    except Exception:
        return "unknown"

# Umbrales ACMG ILUSTRATIVOS (sin calibrar clínicamente).
ACMG_BENIGN_LE = 0.35
ACMG_PATHOGENIC_GE = 0.65


def _acmg(prob):
    if prob >= ACMG_PATHOGENIC_GE:
        return {"code": "PP3", "interpretation": "supporting pathogenic"}
    if prob <= ACMG_BENIGN_LE:
        return {"code": "BP4", "interpretation": "supporting benign"}
    return {"code": "—", "interpretation": "indeterminate"}


class SynonymousSigmaPredictor:
    FEATURES = ["sigma_signed", "is_GtoA"]

    def __init__(self, coef, intercept, mean, std, meta=None):
        self.coef = np.asarray(coef, dtype=float)
        self.intercept = float(intercept)
        self.mean = np.asarray(mean, dtype=float)
        self.std = np.asarray(std, dtype=float)
        self.meta = meta or {}

    @classmethod
    def load(cls, path=_MODEL_PATH):
        d = json.load(open(path))
        meta = {k: d[k] for k in ("params", "version", "training", "auc",
                                  "reference_cnn_auc", "window", "local_w") if k in d}
        return cls(d["coef"], d["intercept"], d["mean"], d["std"], meta)

    def featurize(self, cds, cds_pos, ref, alt):
        s = sigma_signed(cds, cds_pos, ref, alt)
        if s is None:
            return None
        return np.array([s, is_g_to_a(ref, alt)], dtype=float)

    def predict_proba(self, cds, cds_pos, ref, alt):
        x = self.featurize(cds, cds_pos, ref, alt)
        if x is None:
            return None
        z = (x - self.mean) / self.std
        logit = self.intercept + float(self.coef @ z)
        return float(1.0 / (1.0 + np.exp(-logit)))

    def score(self, cds, cds_pos, ref, alt, gene=None, transcript=None):
        """Dossier completo (glass-box). Devuelve dict, o None si mismatch de CDS."""
        x = self.featurize(cds, cds_pos, ref, alt)
        if x is None:
            return None
        z = (x - self.mean) / self.std
        contrib = {f: float(self.coef[i] * z[i]) for i, f in enumerate(self.FEATURES)}
        logit = self.intercept + sum(contrib.values())
        prob = float(1.0 / (1.0 + np.exp(-logit)))
        aa = len(cds) // 3
        return {
            "tool": "ef-synonymous", "tool_version": _pkg_version(),
            "use": "RESEARCH USE ONLY — not a medical device",
            "model": {
                "name": "synonymous sigma+G>A (glass-box logistic)",
                "params": self.meta.get("params", "Turner (RNA)"),
                "version": self.meta.get("version"),
                "coef": self.coef.tolist(), "intercept": self.intercept,
                "mean": self.mean.tolist(), "std": self.std.tolist(),
                "training": self.meta.get("training"), "auc": self.meta.get("auc"),
                "reference_cnn_auc": self.meta.get("reference_cnn_auc"),
                "patent": "P202630522 (OEPM)",
                "paper": "Zenodo 10.5281/zenodo.20275792",
            },
            "input": {
                "gene": gene, "transcript": transcript,
                "cds_pos": cds_pos, "ref": ref.upper(), "alt": alt.upper(),
                "cds_length": len(cds),
            },
            "features": {"sigma_signed_kcal_mol": float(x[0]), "is_GtoA": int(x[1])},
            "zscores": {"sigma_signed": float(z[0]), "is_GtoA": float(z[1])},
            "contributions_to_logit": dict({"intercept": self.intercept}, **contrib),
            "logit": logit, "score_sigma": prob,
            "acmg": dict(_acmg(prob), thresholds={"benign_le": ACMG_BENIGN_LE,
                        "pathogenic_ge": ACMG_PATHOGENIC_GE}, calibrated=False,
                        note="illustrative thresholds, not clinically calibrated"),
            "synonymy": synonymy(cds, cds_pos, ref, alt),
            "applicability": {"protein_aa": aa, "reduced_reliability": aa > 2500,
                              "basis": "P-L1b length"},
            "disclaimer": "Supporting computational evidence (AUC ~0.68). "
                          "Do not use for clinical decisions.",
        }
