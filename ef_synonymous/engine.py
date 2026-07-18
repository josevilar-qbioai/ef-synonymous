"""
Motor termodinámico self-contained (numpy-only) para EF-Synonymous.

Portado 1:1 de core/energy.py + core/synonymous_sigma_predictor.py, pero SIN
dependencias del proyecto grande: aquí van embebidos los parámetros Turner y la
lógica de featurización. El paquete corre solo con numpy — sin PyTorch, sin GPU.

σ(x) = suma con signo del ΔΔG de apilamiento nearest-neighbor (Turner, RNA) en
±10 nt alrededor de la variante. Verificado idéntico a la referencia Python y a
la implementación JS del prototipo web (tolerancia 1e-6).

Autor: Jose Antonio Vilar Sanchez — QMetrika Labs.
"""
import numpy as np

# Parámetros de apilamiento nearest-neighbor, Turner (RNA), kcal/mol.
# Claves en alfabeto DNA (T); la secuencia RNA se normaliza U->T.
STACKING_TURNER = {
    "AA": -0.93, "AT": -1.10, "AG": -2.08, "AC": -2.24,
    "TA": -1.33, "TT": -0.93, "TG": -2.35, "TC": -2.11,
    "GA": -2.35, "GT": -2.24, "GG": -3.26, "GC": -3.42,
    "CA": -2.11, "CT": -2.08, "CG": -2.36, "CC": -3.26,
}

WINDOW = 128
LOCAL_W = 10


def stacking_profile(sequence, params=STACKING_TURNER):
    """Perfil de energía de apilamiento nucleótido a nucleótido (idéntico a core)."""
    seq = sequence.upper().replace("U", "T")
    n = len(seq)
    if n < 2:
        return np.zeros(n)
    bonds = np.array([params.get(seq[i:i + 2], 0.0) for i in range(n - 1)])
    profile = np.zeros(n)
    profile[0] = bonds[0]
    profile[-1] = bonds[-1]
    for i in range(1, n - 1):
        profile[i] = (bonds[i - 1] + bonds[i]) / 2
    return profile


def sigma_signed(cds, cds_pos, ref, alt, params=STACKING_TURNER, ws=WINDOW, w=LOCAL_W):
    """σ signed local. Devuelve float, o None si el ref no casa con el CDS.

    Args:
        cds: secuencia CDS del transcrito (str, A/C/G/T)
        cds_pos: posición 1-based en el CDS
        ref, alt: alelos a NIVEL DE CDS (hebra codificante). Para hebra − usar el
                  complemento, o parsear el HGVS c.{pos}{REF}>{ALT}.
    """
    p = cds_pos - 1
    if p < 0 or p >= len(cds) or cds[p].upper() != ref.upper():
        return None
    h = ws // 2
    start = max(0, p - h)
    end = min(len(cds), p + h)
    wt = cds[start:end]
    mut = (cds[:p] + alt + cds[p + 1:])[start:end]
    c = p - start
    fw = stacking_profile(wt, params)
    fm = stacking_profile(mut, params)
    n = min(len(fw), len(fm))
    d = fm[:n] - fw[:n]
    lo, hi = max(0, c - w), min(n, c + w + 1)
    return float(d[lo:hi].sum())


def is_g_to_a(ref, alt):
    """Transición G>A a nivel de CDS (o C>U tras normalizar T->U)."""
    r = str(ref).upper().replace("T", "U")
    a = str(alt).upper().replace("T", "U")
    return int(r == "G" and a == "A")


# Código genético estándar (chequeo de sinonimia).
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def synonymy(cds, cds_pos, ref, alt):
    """Devuelve dict con codón WT/mut, aminoácidos y si es sinónima, o None."""
    p = cds_pos - 1
    cs = (p // 3) * 3
    if cs + 3 > len(cds):
        return None
    wt_c = cds[cs:cs + 3].upper()
    off = p - cs
    mut_c = wt_c[:off] + alt.upper() + wt_c[off + 1:]
    a_wt, a_mut = CODON_TABLE.get(wt_c), CODON_TABLE.get(mut_c)
    return {
        "wt_codon": wt_c, "mut_codon": mut_c, "wt_aa": a_wt, "mut_aa": a_mut,
        "synonymous": bool(a_wt and a_mut and a_wt == a_mut),
    }
