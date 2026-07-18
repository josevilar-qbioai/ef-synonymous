"""
CLI de EF-Synonymous.

Ejemplos:
    ef-syn score --cds-file NM_000719_cds.fasta --hgvs c.213G>A --gene CACNA1C
    ef-syn score --cds ATGGTC... --pos 213 --ref G --alt A --json
    ef-syn score --transcript ENST00000399655 --hgvs c.213G>A   # fetch Ensembl

Inferencia solo-numpy. El mismo motor (Turner + σ) del prototipo web, verificado
a 1e-6. Con --json imprime el dossier auditable completo.
"""
import argparse
import json
import re
import sys

from .predictor import SynonymousSigmaPredictor

_HGVS = re.compile(r"c\.(\d+)([ACGTUacgtu])>([ACGTUacgtu])")


def parse_hgvs(s):
    m = _HGVS.search(s or "")
    if not m:
        return None
    return int(m.group(1)), m.group(2).upper().replace("U", "T"), m.group(3).upper().replace("U", "T")


def read_cds(literal=None, path=None):
    """Devuelve la secuencia CDS desde --cds (literal) o --cds-file (fasta/raw)."""
    if literal:
        seq = literal
    elif path:
        seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    else:
        return None
    seq = re.sub(r"\s", "", seq).upper().replace("U", "T")
    return seq if re.fullmatch(r"[ACGT]+", seq or "") else None


def fetch_cds(transcript_id, timeout=20):
    """Descarga el CDS de un transcrito desde Ensembl REST (runtime, red del usuario)."""
    import urllib.request
    url = "https://rest.ensembl.org/sequence/id/%s?type=cds;content-type=application/json" % transcript_id
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)["seq"].upper()


def _human(d):
    g = d["input"]["gene"] or "?"
    inp = d["input"]
    print("EF-Synonymous  ·  c.%d%s>%s  (%s)" % (inp["cds_pos"], inp["ref"], inp["alt"], g))
    print("  score_sigma (P patogénica): %.3f" % d["score_sigma"])
    a = d["acmg"]
    print("  ACMG: %s (%s)  [umbrales ilustrativos, sin calibrar]" % (a["code"], a["interpretation"]))
    syn = d["synonymy"]
    if syn:
        tag = "sinónima ✓" if syn["synonymous"] else "NO sinónima ⚠ (fuera de dominio)"
        print("  codón: %s(%s) -> %s(%s)  — %s" % (syn["wt_codon"], syn["wt_aa"],
              syn["mut_codon"], syn["mut_aa"], tag))
    ap = d["applicability"]
    if ap["reduced_reliability"]:
        print("  applicability: proteína grande (~%d aa) → fiabilidad reducida (P-L1b)" % ap["protein_aa"])
    print("  glass-box (contribuciones al logit):")
    c = d["contributions_to_logit"]
    print("    intercepto            %+.3f" % c["intercept"])
    print("    sigma_signed  %6.2f  %+.3f" % (d["features"]["sigma_signed_kcal_mol"], c["sigma_signed"]))
    print("    is_GtoA       %6d  %+.3f" % (d["features"]["is_GtoA"], c["is_GtoA"]))
    print("    logit = %.3f  ->  P = %.3f" % (d["logit"], d["score_sigma"]))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ef-syn", description="Glass-box scoring de variantes sinónimas (σ + G>A).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score", help="puntuar una variante sinónima")
    src = s.add_mutually_exclusive_group(required=True)
    src.add_argument("--cds", help="secuencia CDS literal (A/C/G/T)")
    src.add_argument("--cds-file", help="fichero FASTA o texto con el CDS")
    src.add_argument("--transcript", help="ID de transcrito Ensembl (descarga el CDS)")
    s.add_argument("--hgvs", help="variante c.{pos}{REF}>{ALT}")
    s.add_argument("--pos", type=int, help="posición 1-based en el CDS")
    s.add_argument("--ref", help="alelo de referencia (nivel CDS)")
    s.add_argument("--alt", help="alelo alternativo (nivel CDS)")
    s.add_argument("--gene", help="nombre del gen (metadato)")
    s.add_argument("--json", action="store_true", help="imprime el dossier auditable JSON")
    args = ap.parse_args(argv)

    # CDS
    if args.transcript:
        try:
            cds = fetch_cds(args.transcript)
        except Exception as e:
            sys.exit("Error al descargar el CDS de %s: %s" % (args.transcript, e))
    else:
        cds = read_cds(args.cds, args.cds_file)
    if not cds:
        sys.exit("CDS inválido o vacío (solo A/C/G/T).")

    # variante
    if args.hgvs:
        p = parse_hgvs(args.hgvs)
        if not p:
            sys.exit("HGVS no reconocida. Usa c.213G>A.")
        pos, ref, alt = p
    elif args.pos and args.ref and args.alt:
        pos, ref, alt = args.pos, args.ref.upper(), args.alt.upper()
    else:
        sys.exit("Especifica la variante con --hgvs o con --pos/--ref/--alt.")

    if pos < 1 or pos > len(cds):
        sys.exit("La posición c.%d está fuera del CDS (longitud %d)." % (pos, len(cds)))
    if cds[pos - 1] != ref:
        sys.exit("Mismatch: en c.%d el CDS tiene %s, no %s. ¿Transcrito/hebra correctos?"
                 % (pos, cds[pos - 1], ref))

    clf = SynonymousSigmaPredictor.load()
    d = clf.score(cds, pos, ref, alt, gene=args.gene, transcript=args.transcript)
    if args.json:
        print(json.dumps(d, indent=2))
    else:
        _human(d)


if __name__ == "__main__":
    main()
