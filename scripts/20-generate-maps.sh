#!/bin/bash
set -euo pipefail

mkdir -p output

python3 src/build_hpa_hgnc_map.py -o output/hpa-hgnc-map.csv
python3 src/build_deepcell_types_hgnc_map.py -o output/deepcell-types-hgnc-map.csv
python3 src/build_asctb_hgnc_map.py -o output/asctb-hgnc-map.csv
python3 src/build_omap_hgnc_map.py -o output/omap-hgnc-map.csv
python3 src/build_avr_hgnc_map.py -o output/avr-hgnc-map.csv --omap output/omap-hgnc-map.csv
python3 src/build_ctann_hgnc_map.py -o output/ctann-hgnc-map.csv
python3 src/build_review_items.py -o output/review-items.csv
