#!/usr/bin/env bash
set -euo pipefail
portable_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export PATH="$portable_root/bin/x86_64-linux:$PATH"
export TEXMFCNF="$portable_root/texmf-dist/web2c"
export TEXMFROOT="$portable_root"
export TEXMFDIST="$portable_root/texmf-dist"
export TEXMFSYSVAR="$portable_root/texmf-var"
export TEXMFSYSCONFIG="$portable_root/texmf-config"
export TEXMFLOCAL="$portable_root/texmf-local"
export TEXMFHOME="$portable_root/texmf-home"
export TEXMFVAR="$portable_root/user-var"
export TEXMFCONFIG="$portable_root/user-config"
export TEXMF="{$TEXMFHOME,!!$TEXMFLOCAL,!!$TEXMFSYSCONFIG,!!$TEXMFSYSVAR,!!$TEXMFDIST}"
exec "$@"
