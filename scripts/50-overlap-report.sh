#!/bin/bash
set -euo pipefail

# Build the HPA x HRA overlap report from the maps made by 20-generate-maps.sh.
python3 src/build_overlap_report.py --maps-dir output --report-dir overlap-report
