#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

PYTHON_VERSION="${PYTHON_VERSION:-3.14}"
PRODUCT_NAME="${PRODUCT_NAME:-MieShield}"
ARCH="${ARCH:-x64}"
BUILD_ENV="${BUILD_ENV:-$ROOT_DIR/.build-venv/linux}"
OUTPUT_DIR="${OUTPUT_DIR:-$ROOT_DIR/dist/linux}"
LTO="${LTO:-yes}"

export UV_CACHE_DIR="${UV_CACHE_DIR:-$ROOT_DIR/.uv-cache}"
export NUITKA_CACHE_DIR="${NUITKA_CACHE_DIR:-$ROOT_DIR/.nuitka-cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$ROOT_DIR/.mplconfig}"

echo "==> Creating build environment: $BUILD_ENV"
uv venv --clear --python "$PYTHON_VERSION" "$BUILD_ENV"

VERSION="${VERSION:-$("$BUILD_ENV/bin/python" -c 'import tomllib; print(tomllib.load(open("pyproject.toml", "rb"))["project"]["version"])')}"
DIST_PATH="$OUTPUT_DIR/mie_shield.dist"
TAR_PATH="$ROOT_DIR/dist/$PRODUCT_NAME-$VERSION-linux-$ARCH.tar.gz"

echo "==> Installing runtime dependencies into build environment"
env UV_PROJECT_ENVIRONMENT="$BUILD_ENV" uv sync --locked

echo "==> Installing build-only dependencies"
uv pip install --python "$BUILD_ENV/bin/python" "nuitka>=4,<5" ordered-set zstandard

echo "==> Removing previous Linux build output"
rm -rf "$DIST_PATH" "$OUTPUT_DIR/mie_shield.build" "$TAR_PATH"
mkdir -p "$OUTPUT_DIR"

NUITKA_ARGS=(
  --mode=standalone
  --enable-plugin=pyside6
  --include-package=PyMieScatt
  --nofollow-import-to=scipy.integrate
  --nofollow-import-to=scipy.integrate.*
  --nofollow-import-to=scipy.stats
  --nofollow-import-to=scipy.stats.*
  --assume-yes-for-downloads
  --lto="$LTO"
  --python-flag=-O
  --python-flag=no_docstrings
  --output-filename="$PRODUCT_NAME"
  --output-dir="$OUTPUT_DIR"
)

echo "==> Building Linux application with Nuitka"
"$BUILD_ENV/bin/python" -m nuitka "${NUITKA_ARGS[@]}" mie_shield.py

if [[ ! -x "$DIST_PATH/$PRODUCT_NAME" ]]; then
  echo "Expected executable was not created: $DIST_PATH/$PRODUCT_NAME" >&2
  exit 1
fi

echo "==> Creating TAR artifact"
tar -C "$OUTPUT_DIR" -czf "$TAR_PATH" mie_shield.dist

echo "==> Done"
echo "Dist: $DIST_PATH"
echo "TAR:  $TAR_PATH"
