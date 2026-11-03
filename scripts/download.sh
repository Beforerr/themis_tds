#!/usr/bin/env bash
# Mirror THEMIS-A daily CDFs (Berkeley server; SPDF via MIRROR=https://spdf.gsfc.nasa.gov/pub/data/themis is several times slower) for event dates. Usage: download.sh <level/product> (e.g. l2/fgm)
set -euo pipefail
cd "$(dirname "$0")/.."
prod=$1; name=$(basename "$prod"); lvl=$(dirname "$prod")
out=data/cdf/$name; mkdir -p "$out"
base=${MIRROR:-https://themis.ssl.berkeley.edu/data/themis}/tha/$prod
for y in $(cut -c1-4 data/dates.txt | sort -u); do
  curl -s "$base/$y/" | grep -o "tha_${lvl}_${name}_${y}[0-9]*_v[0-9]*\.cdf" | sort -u
done > "$out/.listing"
# latest version per event date -> "url path"
awk -F_ '{f[$4]=$0} END{for(d in f) print d, f[d]}' "$out/.listing" | sort | join - data/dates.txt \
  | awk -v b="$base" -v o="$out" '{print b "/" substr($1,1,4) "/" $2, o "/" $2}' \
  | while read -r url path; do [ -s "$path" ] || echo "$url $path"; done \
  | xargs -P 6 -n 2 sh -c 'curl -sf --retry 3 -o "$1.part" "$0" && mv "$1.part" "$1"' || true
echo "$prod: $(ls "$out"/*.cdf 2>/dev/null | wc -l) files"
