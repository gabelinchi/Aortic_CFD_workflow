import sys 
import os 
import os.path as osp 
import numpy as np 
from glob import glob 
import pyvista as pv   
import utils as ut   

intp_options = {
    'zero_boundary_dist': 0.2,  # percentage of border with zero velocity (smooth damping at the border)
    'zero_backflow': False,      # set all backflow components to zero
    'kernel': 'linear',         # RBF interpolation kernel (linear is recommended)
    'smoothing': 0.5,          # interpolation smoothing, range recommended [0, 2]
    'degree': 0,                # degree of polynomial added to the RBF interpolation matrix
    'hard_noslip': False}       


#NOTE: .vtp files are Paraview file formats
#NOTE: This code based on the code from 'Data-driven generation of 4D velocity profiles in the aneurysmal ascending aorta' Saitta et al.

#Mapping function that stores the profiles in the output folder and returns the mappend Polydata as well as the velocities
def vel_mapping(source_profile_dir, target_plane, outputDir, intp_options):

    ## Options
    saveName = 'Mapped_velocity_profile'     # filename of mapped .vtp files
    flip_normals = False                      # usually set to True, but might have to change depending on target plane orientation (how to check this?)
    vel_outputDir = osp.join(outputDir, r'Mapped_Velocity_Profiles')

    print('Start velocity profile mapping')
    ## Read data
    # Reads all PolyData, extracts the points, calculates the COM based on these points and calculates the normals in the COG
    source_profiles = [pv.read(i) for i in sorted(glob(osp.join(source_profile_dir, '*.vtp')))]
    num_frames = len(source_profiles)
    source_pts = [source_profiles[k].points for k in range(num_frames)] 
    source_coms = [source_pts[k].mean(0) for k in range(num_frames)]
    target_plane = target_plane.extract_surface()
    target_pts = target_plane.points
    leftmost_idx_on_target = min(range(len(target_pts[: ,0])), key = target_pts[: ,1].__getitem__)             # index of the leftmost point (most negative in y direction) in the target plane w.r.t the subject
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

    print('com after interpolation', interp_planes[0].points.mean(0) + target_com)
    # recenters the velocity profiles at the target profile origin for further modifications
    for k in range(num_frames):
        interp_planes[k].points = interp_planes[k].points + target_com #+= target_com

    print('com final', interp_planes[0].points.mean(0))
    vel_final = [interp_planes[k]['Velocity'] for k in range(num_frames)]

    ## Save profiles to .vtp (may be removed from final)
    os.makedirs(vel_outputDir, exist_ok=True) #Makes the output directory according to the path. If this one already exists there is no error raised.
    for k in range(num_frames):
        interp_planes[k].save(osp.join(vel_outputDir, saveName + '_{:02d}.vtp'.format(k)))

    print('Velocity profile mapping done')

    return interp_planes, num_frames, source_profiles

velocity_map, n_maps, source_profiles = vel_mapping(r'C:\Users\lmorr\Documents\TU\23-24\BEP\Velocity_profiles', pv.read('test_inlet.vtk'), r'C:\Users\lmorr\Documents\TU\23-24\BEP\Git_repository\Aortic_CFD_workflow-3', intp_options)

print(velocity_map[0]['Velocity'])
n = 0
for i in velocity_map:
    i = i.extract_surface()
    source_profiles[n].plot()
    if True:
        plt = pv.Plotter()
        #plt.add_mesh(pv.read('test_inlet.vtk'), show_edges = True, color = 'black')
        plt.add_arrows(i.points, 20 * i['Velocity'], color = 'black')
        plt.add_mesh(i)
        plt.show()
    n = n + 1

