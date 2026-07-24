# Changelog

Todas las versiones notables de `ef-synonymous`. Formato basado en
[Keep a Changelog](https://keepachangelog.com/); versionado
[SemVer](https://semver.org/).

## [0.2.2] — 2026-07-24

### Changed
- Actualizada la referencia de la CNN dedicada del Paper 1 tras reentrenar con
  parámetros Turner (RNA): **AUC 0.683 → 0.709** (out-of-fold 0.693). El modelo
  σ+G>A de 2 variables (0.671 pooled / 0.680 leave-genes-out) recupera la mayor
  parte y queda a la par con la CNN en el set gene-matched (CNN 0.659).
- `model.json`: `reference_cnn_auc` 0.683 → 0.709 (aparece en el dossier auditable).
- README y docstrings actualizados. El motor σ+G>A no cambia — solo la cifra de referencia.

## [0.2.1] — 2026-07-18

### Changed
- README bilingüe (English + Español) para la ficha de PyPI.
- `pip install ef-synonymous` como instalación principal (ya publicado en PyPI).

## [0.2.0] — 2026-07-18

### Added
- Motor self-contained `engine.py` (numpy-only): parámetros Turner embebidos,
  `stacking_profile`, `sigma_signed`, `is_g_to_a`, `synonymy`. Sin dependencia de
  `core/` ni PyTorch.
- `SynonymousSigmaPredictor` con `predict_proba` y `score` (dossier auditable
  completo: features, z-scores, contribuciones al logit, ACMG, sinonimia,
  applicability, versión y coeficientes del modelo).
- CLI `ef-syn score`: entrada por CDS literal, FASTA o descarga por transcrito
  (Ensembl REST); salida humana o `--json`.
- Plugin Ensembl VEP `EFSynonymous.pm` (`EF_sigma`, `EF_score`, `EF_acmg` para
  variantes `synonymous_variant`).
- Tests de ground truth (CACNA1C c.213G>A → σ=1.67, P=0.688068), verificados
  también contra la implementación JS del prototipo web (tolerancia 1e-6).
- Empaquetado `pyproject.toml` (entry point `ef-syn`, `model.json` como
  package-data), CI GitHub Actions (Python 3.8/3.10/3.12).

### Model
- σ + G>A (regresión logística de 2 variables, parámetros Turner de RNA).
  Entrenado sobre ClinVar sinónimas v2 (n=1957; 657P/1300B; 456 genes).
  AUC 0.671 pooled / 0.680 leave-genes-out ≈ CNN dedicada del Paper 1 (0.683).

### Notes
- Research use only. Umbrales de mapeo ACMG ilustrativos, sin calibrar.
- Método cubierto por patente P202630522 (OEPM); código bajo MIT.
