#!/bin/bash
# Project ArterialAtlas.nii (MNI volume) to fsLR-32k surfaces.
# ArterialAtlas_level2.nii is intentionally not used.

set -euo pipefail

ATLAS_ROOT="/host/verges/tank/data/mo/software/Arterial_Atlas-main"
ATLAS_NII="${ATLAS_ROOT}/data/Atlas/ArterialAtlas.nii"
OUT_DIR="/host/verges/tank/data/mo/atlas/ArterialAtlas"

# These must be fsLR-32k midthickness surfaces whose coordinates are in MNI space.
# If your files have different names, change only these two variables.
SURF_L="/host/verges/tank/data/mo/software/CanlabCore-master/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii"

SURF_R="/host/verges/tank/data/mo/software/CanlabCore-master/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii"

mkdir -p "${OUT_DIR}"

for f in "${ATLAS_NII}" "${SURF_L}" "${SURF_R}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required file not found: ${f}" >&2
        exit 1
    fi
done

command -v wb_command >/dev/null 2>&1 || {
    echo "ERROR: wb_command is not available in PATH" >&2
    exit 1
}

# -enclosing assigns the value of the voxel containing each surface vertex.
# This preserves the integer arterial-territory labels.
wb_command -volume-to-surface-mapping \
    "${ATLAS_NII}" "${SURF_L}" \
    "${OUT_DIR}/lh.ArterialAtlas.32k_fsLR.func.gii" \
    -enclosing

wb_command -volume-to-surface-mapping \
    "${ATLAS_NII}" "${SURF_R}" \
    "${OUT_DIR}/rh.ArterialAtlas.32k_fsLR.func.gii" \
    -enclosing

cp -f "${ATLAS_ROOT}/data/ArterialAtlasLables.txt" \
      "${OUT_DIR}/ArterialAtlasLabels.txt"

echo "COMPLETED"
echo "Left : ${OUT_DIR}/lh.ArterialAtlas.32k_fsLR.func.gii"
echo "Right: ${OUT_DIR}/rh.ArterialAtlas.32k_fsLR.func.gii"

