# ef-synonymous

Scoring **glass-box** de patogenicidad de variantes **sinónimas** a partir de un
único observable físico interpretable: el ΔΔG de apilamiento *nearest-neighbor*
con signo (parámetros Turner de RNA) más el sesgo mutacional G>A.

Un modelo de **2 variables** que reproduce a la CNN dedicada del Paper 1
(AUC 0.671 pooled / 0.680 leave-genes-out ≈ CNN 0.683). Inferencia **solo con
numpy** — sin PyTorch, sin GPU, sin llamadas a servicios externos. Es la misma
lógica del prototipo web y de `core/synonymous_sigma_predictor.py`, verificada a
1e-6.

> ⚠ **Research use only.** No es un dispositivo médico. Los umbrales de mapeo
> ACMG son ilustrativos y no están calibrados clínicamente.

## Instalación

```bash
pip install -e ef-synonymous-tool     # desde el repo
# o, publicado:  pip install ef-synonymous
```

## CLI

```bash
# CDS literal + HGVS
ef-syn score --cds ATGGTC...TAG --hgvs c.213G>A --gene CACNA1C

# desde un FASTA local
ef-syn score --cds-file NM_000719.cds.fasta --pos 213 --ref G --alt A

# descargando el CDS de un transcrito (Ensembl REST, red del usuario)
ef-syn score --transcript ENST00000399655 --hgvs c.213G>A

# dossier auditable JSON (base del argumento SaMD)
ef-syn score --cds-file cds.fasta --hgvs c.213G>A --json
```

Salida humana:

```
EF-Synonymous  ·  c.213G>A  (CACNA1C)
  score_sigma (P patogénica): 0.688
  ACMG: PP3 (supporting pathogenic)  [umbrales ilustrativos, sin calibrar]
  codón: GCG(A) -> GCA(A)  — sinónima ✓
  glass-box (contribuciones al logit):
    intercepto            -0.062
    sigma_signed    1.67  +0.071
    is_GtoA            1  +0.782
    logit = 0.791  ->  P = 0.688
```

## Librería

```python
from ef_synonymous import SynonymousSigmaPredictor
clf = SynonymousSigmaPredictor.load()
clf.predict_proba(cds, cds_pos=213, ref="G", alt="A")   # 0.688...
clf.score(cds, 213, "G", "A", gene="CACNA1C")           # dossier completo (dict)
```

## Plugin Ensembl VEP

```bash
cp vep_plugin/EFSynonymous.pm ~/.vep/Plugins/
vep -i input.vcf --plugin EFSynonymous
```

Añade `EF_sigma`, `EF_score` y `EF_acmg` a las variantes `synonymous_variant`.
El plugin delega en la CLI (mismo motor). Para escala, la vía de producción es
una cache precomputada por transcrito o un endpoint REST — este plugin es la
implementación de referencia de la integración.

## Verificación

```bash
python -m pytest tests/          # o:  python tests/test_ground_truth.py
```

Reproduce el valor validado (CACNA1C c.213G>A → σ=1.67, P=0.688068), el mismo que
verifica la implementación JavaScript del prototipo web.

## Método y licencia

Código bajo licencia MIT. El **método** de huella termodinámica del mRNA está
cubierto por la patente **P202630522** (OEPM); su uso comercial requiere licencia.
Paper 1: [Zenodo 10.5281/zenodo.20275792](https://doi.org/10.5281/zenodo.20275792).

QMetrika Labs · Jose Antonio Vilar Sánchez.
