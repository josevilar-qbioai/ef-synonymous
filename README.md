# ef-synonymous

[![CI](https://github.com/josevilar-qbioai/ef-synonymous/actions/workflows/ci.yml/badge.svg)](https://github.com/josevilar-qbioai/ef-synonymous/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/ef-synonymous.svg)](https://pypi.org/project/ef-synonymous/)
[![Python](https://img.shields.io/pypi/pyversions/ef-synonymous.svg)](https://pypi.org/project/ef-synonymous/)
[![License: MIT](https://img.shields.io/badge/License-MIT-teal.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/Paper-10.5281%2Fzenodo.20275792-blue)](https://doi.org/10.5281/zenodo.20275792)

**English** · [Español](#español)

---

## English

**Glass-box** pathogenicity scoring for **synonymous** variants from a single
interpretable physical observable: the signed nearest-neighbor stacking ΔΔG
(Turner RNA parameters) plus the G>A mutational bias.

A **2-variable** model that reproduces the dedicated CNN of Paper 1
(AUC 0.671 pooled / 0.680 leave-genes-out; dedicated CNN 0.709, on par on the gene-matched set). Inference is
**numpy-only** — no PyTorch, no GPU, no external services. It is the same logic
as the web prototype and `core/synonymous_sigma_predictor.py`, verified to 1e-6.

> ⚠ **Research use only.** Not a medical device. The ACMG mapping thresholds are
> illustrative and not clinically calibrated.

### Install

```bash
pip install ef-synonymous
```

### CLI

```bash
# literal CDS + HGVS
ef-syn score --cds ATGGTC...TAG --hgvs c.213G>A --gene CACNA1C

# from a local FASTA
ef-syn score --cds-file NM_000719.cds.fasta --pos 213 --ref G --alt A

# fetching the CDS of a transcript (Ensembl REST, user's network)
ef-syn score --transcript ENST00000399655 --hgvs c.213G>A

# auditable JSON dossier (the SaMD argument)
ef-syn score --cds-file cds.fasta --hgvs c.213G>A --json
```

Human-readable output:

```
EF-Synonymous  ·  c.213G>A  (CACNA1C)
  score_sigma (P pathogenic): 0.688
  ACMG: PP3 (supporting pathogenic)  [illustrative thresholds, uncalibrated]
  codon: GCG(A) -> GCA(A)  — synonymous ✓
  glass-box (contributions to the logit):
    intercept             -0.062
    sigma_signed    1.67  +0.071
    is_GtoA            1  +0.782
    logit = 0.791  ->  P = 0.688
```

### Library

```python
from ef_synonymous import SynonymousSigmaPredictor
clf = SynonymousSigmaPredictor.load()
clf.predict_proba(cds, cds_pos=213, ref="G", alt="A")   # 0.688...
clf.score(cds, 213, "G", "A", gene="CACNA1C")           # full dossier (dict)
```

### Ensembl VEP plugin

```bash
cp vep_plugin/EFSynonymous.pm ~/.vep/Plugins/
vep -i input.vcf --plugin EFSynonymous
```

Adds `EF_sigma`, `EF_score` and `EF_acmg` to `synonymous_variant` records. The
plugin delegates to the CLI (same engine). At scale, the production path is a
precomputed per-transcript cache or a REST endpoint — this plugin is the
reference implementation of the integration.

### How it works (and what it is not)

The model is a logistic regression over **two physical observables** that matches
the dedicated CNN of Paper 1 (AUC 0.709, Turner RNA):

- **σ (sigma_signed):** signed sum of nearest-neighbor stacking ΔΔG (Turner RNA
  parameters) over ±10 nt around the variant — local thermodynamic
  destabilization of the mRNA.
- **is_GtoA:** G>A transition (CpG mutational bias, 2.60× enriched in pathogenic).

Everything runs locally with numpy. Each prediction is reproducible, auditable
and traceable variant by variant — the **glass-box argument for SaMD/ACMG**.
AUC ~0.68 = supporting evidence, not a verdict; wide gray zone; illustrative,
uncalibrated ACMG thresholds. Do not use for clinical decisions.

### Method and license

Code under the MIT license. The **method** (mRNA nearest-neighbor thermodynamic
fingerprinting) is covered by patent **P202630522** (OEPM, Spain); commercial use
of the method may require a separate license.
Paper 1: [Zenodo 10.5281/zenodo.20275792](https://doi.org/10.5281/zenodo.20275792).

QMetrika Labs · Jose Antonio Vilar Sánchez.

---

## Español

Scoring **glass-box** de patogenicidad de variantes **sinónimas** a partir de un
único observable físico interpretable: el ΔΔG de apilamiento *nearest-neighbor*
con signo (parámetros Turner de RNA) más el sesgo mutacional G>A.

Un modelo de **2 variables** que reproduce a la CNN dedicada del Paper 1
(AUC 0.671 pooled / 0.680 leave-genes-out; CNN dedicada 0.709, a la par en el set gene-matched). Inferencia **solo con
numpy** — sin PyTorch, sin GPU, sin servicios externos. Es la misma lógica del
prototipo web y de `core/synonymous_sigma_predictor.py`, verificada a 1e-6.

> ⚠ **Research use only.** No es un dispositivo médico. Los umbrales de mapeo
> ACMG son ilustrativos y no están calibrados clínicamente.

### Instalación

```bash
pip install ef-synonymous
```

### CLI

```bash
# CDS literal + HGVS
ef-syn score --cds ATGGTC...TAG --hgvs c.213G>A --gene CACNA1C

# desde un FASTA local
ef-syn score --cds-file NM_000719.cds.fasta --pos 213 --ref G --alt A

# descargando el CDS de un transcrito (Ensembl REST, red del usuario)
ef-syn score --transcript ENST00000399655 --hgvs c.213G>A

# dossier auditable JSON (el argumento SaMD)
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

### Librería

```python
from ef_synonymous import SynonymousSigmaPredictor
clf = SynonymousSigmaPredictor.load()
clf.predict_proba(cds, cds_pos=213, ref="G", alt="A")   # 0.688...
clf.score(cds, 213, "G", "A", gene="CACNA1C")           # dossier completo (dict)
```

### Plugin Ensembl VEP

```bash
cp vep_plugin/EFSynonymous.pm ~/.vep/Plugins/
vep -i input.vcf --plugin EFSynonymous
```

Añade `EF_sigma`, `EF_score` y `EF_acmg` a las variantes `synonymous_variant`. El
plugin delega en la CLI (mismo motor). Para escala, la vía de producción es una
cache precomputada por transcrito o un endpoint REST — este plugin es la
implementación de referencia de la integración.

### Cómo funciona (y qué NO es)

El modelo es una regresión logística de **dos observables físicos** que iguala a
la CNN dedicada del Paper 1 (AUC 0.709, Turner RNA):

- **σ (sigma_signed):** suma con signo del ΔΔG de apilamiento nearest-neighbor
  (parámetros Turner de RNA) en ±10 nt alrededor de la variante — desestabilización
  termodinámica local del mRNA.
- **is_GtoA:** transición G>A (sesgo mutacional CpG, enriquecido 2.60× en patogénicas).

Todo corre en local con numpy. Cada predicción es reproducible, auditable y
trazable variante a variante — el **argumento glass-box para SaMD/ACMG**.
AUC ~0.68 = evidencia de apoyo, no un veredicto; zona gris amplia; umbrales ACMG
ilustrativos y sin calibrar. No usar para decisiones clínicas.

### Método y licencia

Código bajo licencia MIT. El **método** (huella termodinámica nearest-neighbor del
mRNA) está cubierto por la patente **P202630522** (OEPM); su uso comercial puede
requerir licencia.
Paper 1: [Zenodo 10.5281/zenodo.20275792](https://doi.org/10.5281/zenodo.20275792).

QMetrika Labs · Jose Antonio Vilar Sánchez.
