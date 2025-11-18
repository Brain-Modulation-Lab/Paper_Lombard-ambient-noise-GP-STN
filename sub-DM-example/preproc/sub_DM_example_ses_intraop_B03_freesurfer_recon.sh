# run by YOUR NAME 2023 MM DD



# get subject name from user
read -p 'SUBJECT (eg, DM1025): ' SUBJECT_FS
# export SUBJECT_FS=DM1029

# Display the prompt to the user
read -p "Is this the correct subject?: ${SUBJECT_FS} (y/n) " response

# Check the user's response
if [ "$response" = "y" ]; then
    echo "Continuing with the script..."
    # Add your script logic here
else
    echo "Aborting the script. Incorrect subject."
    exit 1
fi



export FREESURFER_HOME=/Applications/freesurfer/7.4.1

export SUBJECTS_DIR=/Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer

source $FREESURFER_HOME/SetUpFreeSurfer.sh


# ensure we are running thes from the right folder
cd /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/


recon-all -subject ${SUBJECTS_DIR} -i /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/anat_t1.nii -autorecon-all -notal-check


chmod -R 777 ./freesurfer/.



mri_convert -i /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/T1.mgz -o /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/T1.nii -it mgz -ot nii

mri_convert -i /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/lh.ribbon.mgz -o /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/gray_left.nii -it mgz -ot nii

mri_convert -i /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/rh.ribbon.mgz -o /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/gray_right.nii -it mgz -ot nii

mri_convert -i /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/ribbon.mgz -o /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/t1_class.nii -it mgz -ot nii

mri_convert -i /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/ribbon.mgz -o /Volumes/Nexus4/DBS/derivatives/sub-${SUBJECT_FS}/freesurfer/freesurfer/mri/t1_class.nii -it mgz -ot nii





# # view volumes and percellations surface 
# freeview -v \
# freesurfer/mri/T1.mgz \
# freesurfer/mri/wm.mgz \
# freesurfer/mri/brainmask.mgz \
# freesurfer/mri/aseg.mgz:colormap=lut:opacity=0.2 \
# -f freesurfer/surf/lh.white:edgecolor=blue \
# freesurfer/surf/lh.pial:edgecolor=red \
# freesurfer/surf/rh.white:edgecolor=blue \
# freesurfer/surf/rh.pial:edgecolor=red



# view cortical surface 
freeview \
-v \
freesurfer/mri/T1.mgz \
-f \
freesurfer/surf/lh.pial:annot=aparc.annot:name=pial_aparc:visible=0 \
freesurfer/surf/lh.pial:annot=aparc.a2009s.annot:name=pial_aparc_des:visible=0 \
freesurfer/surf/lh.inflated:overlay=lh.thickness:overlay_threshold=0.1,3::name=inflated_thickness:visible=0 \
freesurfer/surf/lh.inflated:visible=0 \
freesurfer/surf/lh.white:visible=0 \
freesurfer/surf/lh.pial \
--viewport 3d
