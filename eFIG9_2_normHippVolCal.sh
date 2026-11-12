#!/bin/bash

BASE_DIR="/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/hippunfold_v1.4.1/hippunfold"
FS_DIR="/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives/fastsurfer"
SES="ses-pre"
OUTPUT_CSV="/host/verges/tank/data/mo/demo/52all_hippocampal_volumes_normalized.csv"

echo "subject,hemi,Sub,CA1,CA2,CA3,CA4,DG,SRLM,Cyst,total_hipp_vol,TIV,normalized_hipp_vol" > ${OUTPUT_CSV}

for SUBJ_PATH in "${BASE_DIR}"/sub-PX*; do
    SUBJ=$(basename ${SUBJ_PATH})
    TSV="${SUBJ_PATH}/${SES}/anat/${SUBJ}_${SES}_space-cropT1w_desc-subfields_atlas-multihist7_volumes.tsv"

    # ✅ Confirmed path format: sub-PX1801001_ses-pre
    ASEG="${FS_DIR}/${SUBJ}_${SES}/stats/aseg.stats"

    if [ ! -f "${ASEG}" ]; then
        echo "WARNING: No aseg.stats for ${SUBJ}, TIV set to NA" >&2
        TIV="NA"
    else
        # ✅ Confirmed format: field $4 contains the number
        TIV=$(grep "Estimated Total Intracranial Volume" ${ASEG} | awk -F',' '{print $4}' | tr -d ' ')
        echo "DEBUG: ${SUBJ} TIV=${TIV}" >&2
    fi

    if [ -f "${TSV}" ]; then
        python3 << EOF >> ${OUTPUT_CSV}
import csv

tsv_file = "${TSV}"
tiv_str  = "${TIV}"
subj     = "${SUBJ}"

with open(tsv_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        row = {k.strip(): v.strip() for k, v in row.items()}
        hemi = row['hemi']
        Sub  = float(row['Sub'])
        CA1  = float(row['CA1'])
        CA2  = float(row['CA2'])
        CA3  = float(row['CA3'])
        CA4  = float(row['CA4'])
        DG   = float(row['DG'])
        SRLM = float(row['SRLM'])
        Cyst = float(row['Cyst'])
        total = Sub + CA1 + CA2 + CA3 + CA4 + DG + SRLM

        if tiv_str != "NA":
            tiv = float(tiv_str)
            norm_vol = total / tiv
        else:
            tiv = "NA"
            norm_vol = "NA"

        print(f"{subj},{hemi},{Sub},{CA1},{CA2},{CA3},{CA4},{DG},{SRLM},{Cyst},{total},{tiv},{norm_vol}")
EOF
    else
        echo "WARNING: Missing TSV for ${SUBJ}" >&2
    fi
done

echo "Done! Results saved to ${OUTPUT_CSV}"
