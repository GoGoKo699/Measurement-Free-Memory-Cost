#!/usr/bin/env bash
set -euo pipefail

# Export a portable public TeX toolchain, never manuscript content.
output_dir=${1:-export}
mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd)
portable="$output_dir/texlive2025"
test ! -e "$portable"
mkdir -p "$portable/bin/x86_64-linux" "$portable/texmf-dist" "$portable/texmf-var"

pdflatex --version | tee "$output_dir/compiler-version.txt"
grep -F 'TeX Live 2025' "$output_dir/compiler-version.txt"
kpsewhich revtex4-2.cls | tee "$output_dir/revtex-location.txt"
test -s "$output_dir/revtex-location.txt"
texroot=$(kpsewhich -var-value=TEXMFROOT)
texdist=$(kpsewhich -var-value=TEXMFDIST)
texvar=$(kpsewhich -var-value=TEXMFSYSVAR)
texroot=$(realpath "$texroot")
texdist=$(realpath "$texdist")
texvar=$(realpath "$texvar")
printf '%s\n' "$texroot" > "$output_dir/original-texmf-root.txt"
for notice in LICENSE.TL LICENSE.CTAN; do
  if test -f "$texroot/$notice"; then cp -a "$texroot/$notice" "$portable/$notice"; fi
done

for tool in pdftex kpsewhich; do
  cp -L "$(command -v "$tool")" "$portable/bin/x86_64-linux/$tool"
done
ln -s pdftex "$portable/bin/x86_64-linux/pdflatex"
ln -s pdftex "$portable/bin/x86_64-linux/latex"
ldd "$portable/bin/x86_64-linux/pdftex" > "$output_dir/pdftex-ldd.txt" || true

# Preserve complete package code, not a fragile list inferred from one paper.
# pdfLaTeX uses Type1/TFM/VF/encoding/map resources; omit OTF/TTF/doc/source bulk.
for relative in tex web2c dvips fonts/afm fonts/enc fonts/map fonts/tfm fonts/type1 fonts/vf; do
  if test -d "$texdist/$relative"; then
    mkdir -p "$portable/texmf-dist/$(dirname "$relative")"
    cp -aL "$texdist/$relative" "$portable/texmf-dist/$relative"
  fi
done
for relative in web2c/pdftex fonts; do
  if test -d "$texvar/$relative"; then
    mkdir -p "$portable/texmf-var/$(dirname "$relative")"
    cp -aL "$texvar/$relative" "$portable/texmf-var/$relative"
  fi
done
format=$(kpsewhich -engine=pdftex -format=fmt pdflatex.fmt)
test -s "$format"
mkdir -p "$portable/texmf-var/web2c/pdftex"
cp -L "$format" "$portable/texmf-var/web2c/pdftex/pdflatex.fmt"
mkdir -p "$portable/texmf-home" "$portable/texmf-local" "$portable/texmf-config" "$portable/user-var" "$portable/user-config"

# Use one relocatable wrapper with explicit tree ownership. No local2023 inputs.
cp manuscript-build/run-texlive2025.sh "$portable/run-texlive2025.sh"
chmod 755 "$portable/run-texlive2025.sh" "$portable/bin/x86_64-linux/pdftex" "$portable/bin/x86_64-linux/kpsewhich"
mktexlsr "$portable/texmf-dist" "$portable/texmf-var"
"$portable/run-texlive2025.sh" kpsewhich -all revtex4-2.cls > "$output_dir/exported-revtex-location.txt"
"$portable/run-texlive2025.sh" kpsewhich -engine=pdftex -format=fmt pdflatex.fmt > "$output_dir/exported-format-location.txt"
grep -F "$portable/" "$output_dir/exported-revtex-location.txt"
grep -F "$portable/" "$output_dir/exported-format-location.txt"

mkdir -p "$output_dir/smoke"
cp manuscript-build/toolchain-smoke.tex "$output_dir/smoke/toolchain-smoke.tex"
(
  cd "$output_dir/smoke"
  "$portable/run-texlive2025.sh" pdflatex -recorder -interaction=nonstopmode -halt-on-error toolchain-smoke.tex
  "$portable/run-texlive2025.sh" pdflatex -recorder -interaction=nonstopmode -halt-on-error toolchain-smoke.tex
) > "$output_dir/smoke-console.txt" 2>&1
test -s "$output_dir/smoke/toolchain-smoke.pdf"
grep -F 'TeX Live 2025' "$output_dir/smoke/toolchain-smoke.log"
cp "$output_dir/smoke/toolchain-smoke.log" "$output_dir/smoke.log"
cp "$output_dir/smoke/toolchain-smoke.fls" "$output_dir/smoke.fls"

# Every absolute input in the relocated build must lie inside its exported tree.
# The only other allowed inputs are the intentionally public smoke source and cwd.
while IFS= read -r line; do
  case "$line" in
    'INPUT /'*)
      input=${line#INPUT }
      case "$input" in "$portable/"*|"$output_dir/smoke/"*) ;; *)
        printf 'Nonportable input: %s\n' "$input" >&2
        exit 1
      esac
      ;;
  esac
done < "$output_dir/smoke.fls"

du -sk "$portable" > "$output_dir/uncompressed-size-kib.txt"
find "$portable" -type f | wc -l > "$output_dir/file-count.txt"
# The historic container omits xz. Archive here; the Ubuntu runner compresses it.
tar -C "$output_dir" -cf "$output_dir/texlive2025.tar" texlive2025
ls -lh "$output_dir/texlive2025.tar"
