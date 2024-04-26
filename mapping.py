import sys 
import os 
import os.path as osp 
import numpy as np 
from glob import glob 
import pyvista as pv   
import utils as ut   

#NOTE: .vtp files are Paraview file formats
#NOTE: This code based on the code from 'Data-driven generation of 4D velocity profiles in the aneurysmal ascending aorta' Saitta et al.

#-----------------------------------------------------------------------------------------------------------------------
## Options

saveName = 'Mapped_velocity_profile'     # filename of mapped .vtp files
flip_normals = True                      # usually set to True, but might have to change depending on target plane orientation (how to check this?)
leftmost_idx_on_target = 000             # index of the leftmost point in the target plane w.r.t the subject

#Inputs for the interpolation
intp_options = {
    'zero_boundary_dist': 0.2,  # percentage of border with zero velocity (smooth damping at the border)
    'zero_backflow': True,      # set all backflow components to zero
    'kernel': 'linear',         # RBF interpolation kernel (linear is recommended)
    'smoothing': 0.01,          # interpolation smoothing, range recommended [0, 2]
    'degree': 0,
    'hard_noslip': False}       # degree of polynomial added to the RBF interpolation matrix


#Mapping function that stores the profiles in the output folder and returns the mappend Polydata as well as the velocities
def mapping(source_profile_dir, target_profile_fn, outputDir):

    ## Read data
    # Reads all PolyData, extracts the points, calculates the COM based on these points and calculates the normals in the COG
    source_profiles = [pv.read(i) for i in sorted(glob(osp.join(source_profile_dir, '*.vtp')))]
    target_plane = pv.read(target_profile_fn)
    num_frames = len(source_profiles)
    source_pts = [source_profiles[k].points for k in range(num_frames)] 
    source_coms = [source_pts[k].mean(0) for k in range(num_frames)] 
    target_pts = target_plane.points
    target_com = target_pts.mean(0)
    target_normal = target_plane.compute_normals()['Normals'].mean(0) 
    normals = [source_profiles[k].compute_normals()['Normals'].mean(0) for k in range(num_frames)]
    if flip_normals: normals = [normals[k] * -1 for k in range(num_frames)] #Flips the normals if flip_normals is true

    ## Align source to target

    # center at origin for simplicity (centers everything at their respective origins)
    target_pts -= target_com
    source_pts = [source_pts[k] - source_coms[k] for k in range(num_frames)]

    # normalize w.r.t. max coordinate norm
    targetmax = np.max(np.sqrt(np.sum(target_pts ** 2, axis=1))) #Euclidian norm to calculate the maximum distance to the origin
    pts = [source_pts[k] * targetmax for k in range(num_frames)] #Normalises the points with respect to the maximum distance to the target origin


    # rotate to align normals
    Rots = [ut.rotation_matrix_from_vectors(normals[k], target_normal) for k in range(num_frames)] 
    pts = [Rots[k].dot(pts[k].T).T for k in range(num_frames)]                         #rotates source points to target
    vel = [Rots[k].dot(source_profiles[k]['Velocity'].T).T for k in range(num_frames)] #rotates velocity vectors to target

    # second rotation to ensure consistent in-plane alignment
    lm_ids = [np.argmax(source_pts[k][:, 0]) for k in range(num_frames)] #Grabs the indeces of the points that are the furthest away (as stuff is now in plane this represents the angle on which source and target are rotated)
    Rots_final = [ut.rotation_matrix_from_vectors(pts[k][lm_ids[k], :], target_pts[leftmost_idx_on_target, :]) for k in range(num_frames)] #Calculates the second (final) rotation matrix for in plane rotation
    pts = [Rots_final[k].dot(pts[k].T).T for k in range(num_frames)] #Does the final rotation of points
    vel = [Rots_final[k].dot(vel[k].T).T for k in range(num_frames)] #Does the final rotation of velocity

    # create new polydatas and inserts the rotated points/velocity
    aligned_planes = [source_profiles[k].copy() for k in range(num_frames)] #Creates new same shaped profiles from the source profiles
    for k in range(num_frames): 
        aligned_planes[k].points = pts[k]
        aligned_planes[k]['Velocity'] = vel[k]

    # spatial interpolation 
    interp_planes = ut.interpolate_profiles(aligned_planes, target_pts, intp_options)

    # recenters the velocity profiles at the target profile origin for further modifications
    for k in range(num_frames):
        interp_planes[k].points += target_com

    vel_final = [interp_planes[k]['Velocity'] for k in range(num_frames)]

    ## Save profiles to .vtp (may be removed from final)
    os.makedirs(outputDir, exist_ok=True) #Makes the output directory according to the path. If this one already exists there is no error raised.
    for k in range(num_frames):
        interp_planes[k].save(osp.join(outputDir, saveName + '_{:02d}.vtp'.format(k)))

    return interp_planes, vel_final