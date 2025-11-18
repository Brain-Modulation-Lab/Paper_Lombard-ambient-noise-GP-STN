import slicer
import sys


# Parse command-line arguments dynamically
CT_FIXED = sys.argv[1]
CT_MOVING = sys.argv[2]
CT_MOVING_COREG = sys.argv[3]

# ct1_path = '/Volumes/Nexus4/DBS/derivatives/sub-DM1046/ecogloc/ct1.nii'  # fixed volume (rpostop_ct_2)
# ct2_path = '/Volumes/Nexus4/DBS/derivatives/sub-DM1046/ecogloc/ct2.nii'  # moving volume (preop_ct_1)
# ct2_coregistered_path = '/Volumes/Nexus4/DBS/derivatives/sub-DM1046/ecogloc/ct2_coregistered.nii'

# Load volumes
fixedVolume = slicer.util.loadVolume(CT_FIXED)  # fixed volume (rpostop_ct_2)
movingVolume = slicer.util.loadVolume(CT_MOVING)  # moving volume (preop_ct_1)

# Create an empty output volume node in the Slicer scene
outputVolumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode', 'ct2_coregistered')

# Set parameters including Advanced Optimization settings
parameters = {
    "fixedVolume": fixedVolume.GetID(),
    "movingVolume": movingVolume.GetID(),
    "outputVolume": outputVolumeNode.GetID(),

    # Registration phases (from your GUI screenshot)
    "useRigid": False,
    "useScaleVersor3D": True,            # Rigid+Scale(7 DOF)
    "useScaleSkewVersor3D": False,       # Rigid+Scale+Skew(10 DOF)
    "useAffine": True,                   # Affine(12 DOF)

    # General parameters
    "samplingPercentage": 0.01,
    "splineGridSize": '14,10,12',
    "initializeTransformMode": 'useGeometryAlign',

    # Advanced Optimization parameters (from your second screenshot)
    "numberOfIterations": 1500,
    "maximumStepLength": 0.05,
    "minimumStepLength": 0.001,
    "relaxationFactor": 0.5,
    "translationScale": 1000.0,
    "reproportionScale": 1.0,
    "skewScale": 1.0,
    "maxBSplineDisplacement": 0.0,

    # Explicitly set transforms to None if not needed
    "linearTransform": None,
    "bsplineTransform": None,
}

# Run registration synchronously
cliNode = slicer.cli.runSync(slicer.modules.brainsfit, None, parameters)

# Check for errors
if cliNode.GetStatus() & cliNode.ErrorsMask:
    print(f"Registration failed: {cliNode.GetErrorText()}")
    sys.exit(1)
else:
    print("Registration successful!")

# Save the registered volume to file explicitly
slicer.util.saveNode(outputVolumeNode, CT_MOVING_COREG)

sys.exit()