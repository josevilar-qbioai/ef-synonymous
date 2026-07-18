#!/usr/bin/env bash
#
# release.sh — publica una versión de ef-synonymous en un solo comando.
#
#   ./release.sh 0.2.1            # bump + test + build + check + subir a PyPI
#   ./release.sh 0.2.1 --test     # igual, pero sube a TestPyPI (ensayo)
#   ./release.sh 0.2.1 --dry      # bump + test + build + check, SIN subir
#
# Requisitos: python3, y `pip install build twine`. Tokens en ~/.pypirc
# (usa pypirc.template). Requiere 2FA en la cuenta PyPI/TestPyPI.
#
set -euo pipefail

NEW="${1:-}"
MODE="${2:-}"
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

if [ -z "$NEW" ]; then
  echo "uso: ./release.sh X.Y.Z [--test | --dry]"; exit 1
fi
if ! echo "$NEW" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$'; then
  echo "✗ versión no válida (esperado semver X.Y.Z): $NEW"; exit 1
fi

# 0) Preflight: herramientas necesarias en este intérprete
missing=""
for mod in build twine pytest; do
  python3 -c "import $mod" 2>/dev/null || missing="$missing $mod"
done
if [ -n "$missing" ]; then
  echo "✗ Faltan herramientas en $(python3 -c 'import sys;print(sys.executable)'):$missing"
  echo "  Instálalas con:  pip install$missing"
  exit 1
fi

echo "▶ Preparando release $NEW …"

# 1) Bump de versión (single-source: pyproject + __init__)
python3 - "$NEW" <<'PY'
import re, sys
v = sys.argv[1]
for path, pat in [("pyproject.toml", r'(?m)^version = ".*"'),
                  ("ef_synonymous/__init__.py", r'(?m)^__version__ = ".*"')]:
    s = open(path).read()
    repl = 'version = "%s"' % v if path.endswith("toml") else '__version__ = "%s"' % v
    s2, n = re.subn(pat, repl, s, count=1)
    assert n == 1, "no se encontró la línea de versión en " + path
    open(path, "w").write(s2)
print("  · versión fijada a", v, "en pyproject.toml y __init__.py")
PY

# 2) Recordatorio de CHANGELOG (no bloquea)
grep -q "\[$NEW\]" CHANGELOG.md 2>/dev/null || \
  echo "  ⚠ recuerda añadir la entrada [$NEW] en CHANGELOG.md"

# 3) Tests
echo "▶ Tests…"
python3 -m pytest tests/ -q

# 4) Build limpio
echo "▶ Build…"
rm -rf dist build ./*.egg-info
python3 -m build >/dev/null
python3 -m twine check dist/*

# 5) Subida (según modo)
case "$MODE" in
  --dry)
    echo "✓ Artefactos en dist/ (modo --dry, no se sube). Contenido:"; ls -1 dist/ ;;
  --test)
    echo "▶ Subiendo a TestPyPI…"
    python3 -m twine upload -r testpypi dist/*
    echo "✓ En TestPyPI. Prueba: pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ef-synonymous==$NEW" ;;
  "")
    echo "▶ Subiendo a PyPI (real)…"
    python3 -m twine upload dist/*
    echo "✓ Publicado. pip install ef-synonymous==$NEW"
    echo "  Recuerda: git commit -am 'release $NEW' && git tag v$NEW && git push --tags" ;;
  *)
    echo "✗ modo desconocido: $MODE (usa --test o --dry)"; exit 1 ;;
esac
