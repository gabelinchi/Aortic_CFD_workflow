import sys #used to edit Python Runtime Environment (edit python run environment)
import os #operating system functionality
import os.path as osp #makes the software able to manipulate file paths
import numpy as np #we know this one
from glob import glob #returns a possibly empty list of path names that match the input pathname
                      #which has to be a string containing a path specification in absolute or relative style
import pyvista as pv  #3D plotting and mesh analysis for Visualisation Toolkit (which does 3D computer graphics)
                      #Install in Anaconda prompt: pip install pyvista  
import utils as ut    #This is a seperate python file from research

#NOTE: .vtp files are Paraview file formats, which is a post processing visualisation
# tool, visualisation is for now not necessary, but important to know. check bro

#-----------------------------------------------------------------------------------------------------------------------
## Options

#output .vtp files (do we need to automate this? or call this code as a custom library/function)
outputDir = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Mapping_output' # path for saving resampled .vtp files
saveName = 'Mapped_velocity_profile_test'   # filename of resamples .vtp files

#Input files
source_profile_dir = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Velocity_profiles' # directory containing .vtp files associated to a 4D profile of the flow profile
target_profile_fn = r''  # can be a .stl, .vtk or .vtp file (mesh file containing the inlet)


flip_normals = True # usually set to True, but might have to change depending on target plane orientation (how to check this?)
leftmost_idx_on_target = 000      # index of the leftmost point in the target plane w.r.t the subject

#Inputs for the interpolation
intp_options = {
    'zero_boundary_dist': 0.2, # percentage of border with zero velocity (smooth damping at the border)
    'zero_backflow': True, # set all backflow components to zero
    'kernel': 'linear', # RBF interpolation kernel (linear is recommended)
    'smoothing': 0.01, # interpolation smoothing, range recommended [0, 2]
    'degree': 0,
    'hard_noslip': False} # degree of polynomial added to the RBF interpolation matrix


#-----------------------------------------------------------------------------------------------------------------------
## Read data
# This reads all the flow profile files and reads them with pyvista so that they can be used further. Notice that source_profiles
# is a list and target_plane is a single mesh
source_profiles = [pv.read(i) for i in sorted(glob(osp.join(source_profile_dir, '*.vtp')))]
target_plane = pv.read(target_profile_fn)

num_frames = len(source_profiles) #just takes the amount of velocity profiles
source_pts = [source_profiles[k].points for k in range(num_frames)] # Collects for every flow profile mesh the points
source_coms = [source_pts[k].mean(0) for k in range(num_frames)] #Collects for every flow profile mesh the COM based on the points
target_pts = target_plane.points #same as above for target mesh
target_com = target_pts.mean(0) #COM of target mesh
target_normal = target_plane.compute_normals()['Normals'].mean(0) #Calculates the normal in the COM of the target_plane
normals = [source_profiles[k].compute_normals()['Normals'].mean(0) for k in range(num_frames)] #Calculates the normal in the COM of the flow profiles
if flip_normals: normals = [normals[k] * -1 for k in range(num_frames)] #Flips the normals if flip_normals is true


#-----------------------------------------------------------------------------------------------------------------------
## Align source to target

# center at origin for simplicity (centers everything at their respective origins)
target_pts -= target_com
source_pts = [source_pts[k] - source_coms[k] for k in range(num_frames)]

# normalize w.r.t. max coordinate norm
targetmax = np.max(np.sqrt(np.sum(target_pts ** 2, axis=1))) #Euclidian norm to calculate the maximum distance to the origin
pts = [source_pts[k] * targetmax for k in range(num_frames)] #Normalises the points with respect to the maximum distance to the target origin

# rotate to align normals
Rots = [ut.rotation_matrix_from_vectors(normals[k], target_normal) for k in range(num_frames)] #Uses utils to calculate the rotation matrices between the flow profile
                                                                                               #normals and the target normal 
pts = [Rots[k].dot(pts[k].T).T for k in range(num_frames)] #rotates the normalised velocity profiles points to the target coordinate systems
vel = [Rots[k].dot(source_profiles[k]['Velocity'].T).T for k in range(num_frames)] #rotates the velocity vectors to the target coordinate system

# second rotation to ensure consistent in-plane alignment
lm_ids = [np.argmax(source_pts[k][:, 0]) for k in range(num_frames)] #Grabs the indeces of the points that are the furthest away (as stuff is now in plane this represents the angle on which source and target are rotated)
Rots_final = [ut.rotation_matrix_from_vectors(pts[k][lm_ids[k], :], target_pts[leftmost_idx_on_target, :]) for k in range(num_frames)] #Calculates the final rotation matrix for in plane rotation
pts = [Rots_final[k].dot(pts[k].T).T for k in range(num_frames)] #Does the final rotation of the velocity profile points
vel = [Rots_final[k].dot(vel[k].T).T for k in range(num_frames)] #Does the final rotation of the velocity profile matrices


# create new polydatas
aligned_planes = [source_profiles[k].copy() for k in range(num_frames)] #Creates new same shaped profiles from the velocity profiles
for k in range(num_frames): #Replaces the points and velocity vectors with the new rotated points and velocity vectors
    aligned_planes[k].points = pts[k]
    aligned_planes[k]['Velocity'] = vel[k]

# spatial interpolation (google this but basically makes a new intepolated velocity profile based on the target points)
interp_planes = ut.interpolate_profiles(aligned_planes, target_pts, intp_options)
interp_planes[2].plot()

# recenters the velocity profiles at the target profile origin (this is nice for further modifications)
for k in range(num_frames):
    interp_planes[k].points += target_com


#-----------------------------------------------------------------------------------------------------------------------
## Save profiles to .vtp
os.makedirs(outputDir, exist_ok=True) #Makes the output directory according to the path. If this one already exists there is no error raised.
for k in range(num_frames):
    interp_planes[k].save(osp.join(outputDir, saveName + '_{:02d}.vtp'.format(k))) #Saves the mapped velocity profile in the output directory with the right nummering


