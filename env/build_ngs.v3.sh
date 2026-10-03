#!/bin/bash
# Build the ngs.v3 environment as an exact copy of ngs.v2 plus MACS3 and gffread.
# Run on the login node (needs internet). ngs.v2 itself is not modified.
#
# ngs.v2 mixes conda packages with pip upgrades (pypiper 0.15.1, pipestat 0.13.1, looper 2.1.1 on top of
# older conda metadata), so a plain clone is not exact:
#   1. explicit conda export of ngs.v2 -> new prefix (same builds, from the package cache)
#   2. pip packages that were pip-installed in ngs.v2, at the same versions, --no-deps
#   3. MACS3 and gffread (bioconda) with --freeze-installed (gffread: GFF3->GTF for genome resource builds)
#   4. R packages installed from within R in ngs.v2 (copied)
#   5. write env/ngs.v3.explicit.txt and env/ngs.v3.pip.txt (the pinned definition)
set -euo pipefail
MM=/home/gschotta/bin/micromamba
ROOT=/store24/project24/becgsc_001/micromamba
V2=$ROOT/envs/ngs.v2
V3=$ROOT/envs/ngs.v3
HERE=$(cd "$(dirname "$0")" && pwd)
export MAMBA_ROOT_PREFIX=$ROOT

[ -e "$V3" ] && { echo "ERROR: $V3 exists, remove it first"; exit 1; }

"$MM" env export -p "$V2" --explicit > "$HERE/ngs.v2.explicit.txt"
"$MM" create -y -p "$V3" -f "$HERE/ngs.v2.explicit.txt"

# pip-installed packages in ngs.v2 (INSTALLER == pip), same versions
PIPPKG="divvy==0.6.0 logmuse==0.3.0 looper==2.1.1 piper==0.15.1 pipestat==0.13.1 pydantic-argparse==0.10.0 \
        pydantic-settings==2.14.0 python-dotenv==1.2.2 ubiquerg==0.9.3 yacman==1.0.0"
"$V3/bin/python3" -m pip install --no-deps $PIPPKG

"$MM" install -y -p "$V3" -c conda-forge -c bioconda --freeze-installed macs3 gffread

# R packages that were installed from within R in ngs.v2 (install.packages / BiocManager, not conda):
# (or upgraded there over the conda version); same R 4.5.1 and system libraries -> copy the package folders
# (list: ngs.v3.R_packages_from_v2.txt; afterwards the R libraries of ngs.v2 and ngs.v3 are identical)
for p in $(cat "$HERE/ngs.v3.R_packages_from_v2.txt"); do rm -rf "$V3/lib/R/library/$p"; cp -a "$V2/lib/R/library/$p" "$V3/lib/R/library/"; done

"$MM" env export -p "$V3" --explicit > "$HERE/ngs.v3.explicit.txt"
"$V3/bin/python3" -m pip list --format=freeze > "$HERE/ngs.v3.pip.txt"
"$V3/bin/python3" -m pip check || true
for t in trimmomatic bowtie2 STAR samtools picard featureCounts bamCoverage rsem-calculate-expression macs3 gffread; do
    printf "%-28s %s\n" "$t" "$(ls "$V3/bin/$t" >/dev/null 2>&1 && echo ok || echo MISSING)"
done
"$V3/bin/python3" -c "import pypiper, looper; from importlib.metadata import version as v; print('pypiper', v('piper'), 'pipestat', v('pipestat'), 'looper', v('looper'))"
