"""
Ground truth (CACNA1C c.213G>A, GCG->GCA Ala, sinónima) — el mismo valor que
verifica la implementación JS del prototipo web a 1e-6.
σ = 1.67 ; P(patogénica) = 0.688068.

El CDS completo (marco de lectura desde el ATG) está en cacna1c_cds.txt: la
comprobación de sinonimia requiere el frame real (posición 1 = inicio).
"""
import os
from ef_synonymous import SynonymousSigmaPredictor, sigma_signed, synonymy

CDS = open(os.path.join(os.path.dirname(__file__), "cacna1c_cds.txt")).read().strip()
POS, REF, ALT = 213, "G", "A"


def test_sigma_signed():
    assert abs(sigma_signed(CDS, POS, REF, ALT) - 1.67) < 1e-6


def test_probability():
    clf = SynonymousSigmaPredictor.load()
    assert abs(clf.predict_proba(CDS, POS, REF, ALT) - 0.688068) < 1e-5


def test_synonymy():
    syn = synonymy(CDS, POS, REF, ALT)
    assert syn["synonymous"] is True
    assert syn["wt_codon"] == "GCG" and syn["mut_codon"] == "GCA"
    assert syn["wt_aa"] == syn["mut_aa"] == "A"


def test_mismatch_returns_none():
    assert sigma_signed(CDS, POS, "A", "G") is None


def test_dossier_shape():
    clf = SynonymousSigmaPredictor.load()
    d = clf.score(CDS, POS, REF, ALT, gene="CACNA1C")
    assert d["acmg"]["code"] == "PP3"        # 0.688 >= 0.65
    assert d["model"]["params"] == "Turner (RNA)"
    assert set(d["contributions_to_logit"]) == {"intercept", "sigma_signed", "is_GtoA"}


if __name__ == "__main__":
    for fn in [test_sigma_signed, test_probability, test_synonymy,
               test_mismatch_returns_none, test_dossier_shape]:
        fn(); print("ok:", fn.__name__)
    print("TODOS OK")
