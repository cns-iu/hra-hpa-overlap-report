#!/bin/bash
set -euo pipefail

mkdir -p downloads

curl -fsSL https://raw.githubusercontent.com/vanvalenlab/deepcell-types/refs/heads/master/deepcell_types/channel_mapping.yaml -o downloads/deepcell-types-channels.yaml

curl -fsSL -H "Accept: text/turtle" https://purl.humanatlas.io/collection/hra -o downloads/hra.ttl
curl -fsSL -H "Accept: text/turtle" https://purl.humanatlas.io/graph/avr -o downloads/avr.ttl

curl -fsSL https://www.proteinatlas.org/download/proteinatlas.tsv.zip -o downloads/proteinatlas.tsv.zip
curl -fsSL https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt -o downloads/hgnc_complete_set.txt
curl -fsSL https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/withdrawn.txt -o downloads/hgnc_withdrawn.txt
