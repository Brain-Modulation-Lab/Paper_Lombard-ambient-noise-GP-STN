#!/bin/bash

SUBJECT_FS=DM10XX

CT_FIXED="/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/ecogloc/rpostop_ct.nii"
CT_MOVING="/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/ecogloc/preop_ct.nii"
CT_MOVING_COREG="/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/ecogloc/rpreop_ct_SLICER.nii"

PY_SCRIPT="/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/preproc/sub_${SUBJECT_FS}_ses_intraop_B05a_coreg_ct2ct.py"

/Applications/Slicer.app/Contents/MacOS/Slicer \
  --no-main-window \
  --python-script "$PY_SCRIPT" \
  "$CT_FIXED" "$CT_MOVING" "$CT_MOVING_COREG"




# pop up Slicer GUI for manual verification
/Applications/Slicer.app/Contents/MacOS/Slicer


# --------------------------

# SUBJECT=DM10XX

# # RPREOP_CT=/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/leaddbs/rpreop_ct.nii
# CT_FIXED=/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/ecogloc/rpostop_ct.nii
# CT_MOVING=/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/leaddbs/preop_ct.nii
# CT_MOVING_COREG=/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/leaddbs/rpreop_ct.nii

# # /Applications/Slicer.app/Contents/MacOS/Slicer --no-main-window --python-script ./coregister_ct.py
# /Applications/Slicer.app/Contents/MacOS/Slicer --python-script ./coregister_ct.py --CT_FIXED --CT_MOVING --CT_MOVING_COREG


# sys.exit()