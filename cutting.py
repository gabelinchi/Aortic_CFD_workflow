import numpy as np
import pyvista as pv
import utils as ut
import matplotlib.pyplot as pyplt


#main function that runs the cutter script
def main_cutter(inlet, wall, plot=False):
    '''
    Function that is the main for the cutting of the aorta geometry. It calls centerline, areaselection and cutting functions to perform the cut
    :arg1 inlet: pyvista Polydata
    :arg2 wall: pyvista PolyData

    returns pyvista PolyData of the wall after the cut
    '''
    print('Start cutting')
    #Calculates areas along the centerline of the aorta (only first 40 mm). Also outputs nodes/normals and edge profiles for the final cut and visualisation
    centernodes, centernormals, edgeprofiles, slice_areas = centerline(inlet, wall)

    #Plots the centerline
    if plot:
        plt = pv.Plotter()
        plt.add_mesh(wall, style ='wireframe')
        plt.add_points(centernodes, color = 'red')
        plt.show()
        edgeprofiles.plot()

    #Calculates the smallest area across the first 40mm of the inlet
    smallest_area, smallest_area_index = areaselection(slice_areas)

    print('Smallest cross-section calculated')
    print('Smallest cross-section: ',smallest_area)

    #Calculates the normal in the final cutting plane centerpoint
    normal_final = ut.average_normal(centernormals, 4, smallest_area_index)
    #Final centernode
    center_final = centernodes[smallest_area_index].flatten()

    new_geometry = cut(center_final, normal_final, wall, plot=plot)

    print('Cutting done')    
    return new_geometry

def centerline(inlet, wall, dist=40, flip_norm=False):
    '''
    Function that calculates an approximation of the centerline of the wall geometry
    :arg1 inlet: pyvista Polydata
    :arg2 wall: pyvista PolyData
    :opt arg3: the distance from the inlet at which the function stops calculating

    returns:
    centernodes: the points of which the centerline is build upon
    centernormals: the normals in every centernode, perpendicular to the edgeprofiles
    edgeprofiles: profile from the wall when you cut in the centernode along the centernormal
    slice_ares: internal area of each edgeprofile (used in areaselection function)
    '''
    #From the inlet surface mesh extract the boundary edges
    inlet_boundary = inlet.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)

    #Calculate the centernode of the inlet mesh and the normal of the centernode
    inlet_points = inlet_boundary.points
    inlet_centerpoint = inlet_points.mean(0)
    inlet_normal = inlet.compute_normals()['Normals'].mean(0)
    if flip_norm: inlet_normal=-inlet_normal

    #Calculate the area of the inlet edgeprofile
    inlet_area = inlet.area

    #Create empty storage variables, for centernodes/centernormals/areas
    centernodes = np.empty((0,3))
    centernormals = np.empty((0,3))
    slice_areas = np.array([])


    #Use the inlet information calculated above as initial values in the storage variables (append) and and make the centernode and normal equal to a calculation variable
    centernodes = np.vstack([centernodes, inlet_centerpoint])
    centernormals = np.vstack([centernormals, inlet_normal])
    edgeprofiles = inlet_boundary
    slice_areas = np.append(slice_areas, inlet_area)

    #Calculation variables
    center = inlet_centerpoint
    normal = inlet_normal

    #Creat count variable starting at 0
    count = 0

    #Create a while loop with a max count around 200 or something
    while count < 495:
        #Create a new point by adding the directional vector to the centernode of the previous edgeprofile (previous iteration)
        inter_center = np.add(center, normal)
        #Use this new point and the previous directional vector to make a cut of the wall mesh
        inter_profile = get_clip_perimeter(inter_center, normal, wall) #Extract the edge profile (temporary) and isolate the profile (with connectivity) which we want to continue working with
        #From the isolated profile (temporary), calculate the new center point
        inter_profile_points = inter_profile.points
        new_center = inter_profile_points.mean(0)    
        #Calculate the normalised relative vector from the centernode of the previous edgeprofile (previous iteration) and the new center point
        A = np.subtract(inter_center, center)
        B =  0.1 * np.subtract(new_center, inter_center)
        new_normal = ut.normalise(np.add(A, B))

        #Make a new cut of the wall mesh with the new center point and the new normalised directional vector (also isolate this edgeprofile again)
        new_profile = get_clip_perimeter(new_center, new_normal, wall)

        #Calculate the area of the new profile
        new_profile_closed = new_profile.delaunay_2d()
        new_area = new_profile_closed.area

        #Store the new center point, the normalised directional vector, the last (second) edge profile and the area of the last edge profile in the storage variables
        centernodes = np.vstack([centernodes, new_center])
        centernormals = np.vstack([centernormals, new_normal])
        edgeprofiles += new_profile
        slice_areas = np.append(slice_areas, new_area)

        #Calculate the distance of the new center point and the center point of the inlet, if the distance is larger than a certain threshold. We break the while loop
        if np.linalg.norm(new_center - inlet_centerpoint) >= dist:
            break

        #Make the calculation variables the new center point and the new normalised directional vector
        center = new_center
        normal = new_normal

        count += 1

    print('Centerline generated')
    
    return centernodes, centernormals, edgeprofiles, slice_areas

def areaselection(areas):
    '''
    Function that selects the minimal area after the aortic root along the centerline
    :arg1 areas: all areas of the slices along the centerline

    returns
    smallest area: the minimal area after the aortic root
    smallest area index: index of smallest area in areas
    '''
    # Find the indices where the tube transitions from wide to narrow and narrow to wide
    diff_areas = np.diff(areas)
    transition_indices = np.where(np.logical_or(diff_areas < 0, diff_areas == 0))[0]

    #Find the areas between these transition points
    area_segments = np.split(areas, transition_indices + 1)

    #Select the segments that have a narrow part and are larger than 1. This result in multiple arrays that go from a smaller to a larger area
    narrow_segments = [segment for segment in area_segments if len(segment) > 1 and segment[0] < segment[-1]]

    #Select the smallest area from the second small segment. This should be the constriction after the aortic root
    if len(narrow_segments) > 1:
        smallest_area = min(narrow_segments[1]) #Gives out of bound error when there is no constriction
    else:
        counter = 0 
        window_size = 3
        moving_averages = []

        while counter < len(diff_areas) - window_size + 1:
   
            # Store elements from i to i+window_size
            # in list to get the current window
            window = diff_areas[counter : counter + window_size]
        
            # Calculate the average of current window
            window_average = round(sum(window) / window_size, 2)
            
            # Store the average of current
            # window in moving average list
            moving_averages.append(window_average)
            
            # Shift window to right by one position
            counter += 1
        transition_index = np.argmin(moving_averages)
        smallest_area = areas[transition_index]

    #Grab the index of the smallest_area
    smallest_index_calc = np.where(smallest_area == areas)
    smallest_area_index = smallest_index_calc[0]

    return smallest_area, smallest_area_index

def post_cutter(inlet, wall, plot=False, flip_norm=False):
    '''
    Cuts a vessel at the point where the centerline is horizontal, returns downstream geomtetry.
    '''
    print('Start cutting')
    #Calculates areas along the centerline of the aorta (only first 40 mm). Also outputs nodes/normals and edge profiles for the final cut and visualisation
    centernodes, centernormals, edgeprofiles, slice_areas = centerline(inlet, wall, dist = 150, flip_norm=flip_norm)

    #Plots the centerline
    if plot:
        plt = pv.Plotter()
        plt.add_mesh(wall, style ='wireframe')
        plt.add_points(centernodes, color = 'red')
        plt.show()
        edgeprofiles.plot()

    #Finds most horizontal part of the centerline
    horiz_index = np.argmin(np.abs(centernormals @ np.array([0,0,1])))
    print('Horizontal normal found')

    #Calculates the normal in the final cutting plane centerpoint
    normal_final = ut.average_normal(centernormals, 2, horiz_index)

    #Final centernode
    center_final = centernodes[horiz_index].flatten()
    new_geometry = cut(center_final, -normal_final, wall, plot=plot)

    print('Cutting done')    
    return new_geometry

def cut(point, normal, wall, plot=False):
    """
    Function that cuts a vessel geometry along a plane defined by a point & vector, and keeps all regions upstream
    of the point.
    :arg1 point: numpy array, xyz
    :arg2 normal: vector
    :arg3 wall: pyvista PolyData

    returns clipped geometry as PolyData
    """
    # Clip geometry
    reg1, reg2 = wall.clip(normal=normal, origin=point, return_clipped=True) #reg2 is in direction of vector

    if plot==True:
        plt = pv.Plotter()
        plt.add_mesh(reg1, style= 'wireframe', color='green')
        plt.add_mesh(reg2, style= 'wireframe', color='red')
        plt.add_points(point)
        plt.show()

    # Extract geometry to keep from reg2
    reg2 = reg2.connectivity('all')
    keep_id = reg2.point_data['RegionId'][reg2.find_closest_point(point)]
    reg2 = reg2.extract_cells(np.where(reg2.cell_data['RegionId'] == keep_id)[0])
    
    # Extract geometry to keep from reg1
    reg1 = reg1.connectivity('all')
    del_id = reg1.point_data['RegionId'][reg1.find_closest_point(point)]
    reg1 = reg1.extract_cells(np.where(reg1.cell_data['RegionId'] != del_id)[0])

    # Recombine regions
    clipped = reg1.merge(reg2).clean()
    #clipped.clear_data() #Commented to fix postproc, breaks main

    if plot==True:
        clipped.plot()
    return(clipped)

def get_clip_perimeter(point, normal, wall, plot=False):
    '''
    Function that cuts a vessel geometry along a plane defined by a point & vector, and returns the resulting perimeter
    :arg1 point: numpy array, xyz
    :arg2 normal: vector
    :arg3 wall: pyvista PolyData

    returns pyvista PolyData of the perimeter resulting from the cut
    '''
    # Clip geometry
    normal=np.array(normal)
    clipped = wall.clip(normal=(normal * -1), origin=point) 

    # Extract edges
    edges = clipped.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)

    # Extract edge closest to point
    edge = edges.connectivity('closest', point)

    if plot==True:
        edges.plot()
        edge.plot()
    return(edge)

#Development and testing code
#---------------------------#
def ref_geom():
    def make_points():
        """Helper to make XYZ points"""
        theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
        z = np.linspace(-2, 2, 100)
        r = z**2 + 1
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        return (np.column_stack((x, y, z))*4, r*4)

    points, r = make_points()
    spline = pv.Spline(points, 1000)
    spline["theta"] = 0.4 * np.arange(len(spline.points))
    spline["radius"] = np.abs(np.sin(spline["theta"]))*0.3 + 1
    tube = spline.tube(scalars="radius", absolute=True)
    tube.plot(scalars="theta", smooth_shading=True)

    scale = 5
    tube=tube.scale(scale)
    points=spline.points*scale

    """ tubestart = cut(points[0]+0.02*(points[1]-points[0]), -(points[1]-points[0]), tube)
    tubestart = tubestart.connectivity('closest', points[0]).extract_surface()
    tubecenternodes = centerline(tubestart, tube, 1000, flip_norm=True)[0]
    pv.lines_from_points(tubecenternodes).plot()
    pv.lines_from_points(tubecenternodes).save('C:/Users/Yarran/OneDrive/Uni_Documenten/2023-24/BEP/VMTK/ourcenterline.vtp')
    tube.save('C:/Users/Yarran/OneDrive/Uni_Documenten/2023-24/BEP/VMTK/tube.vtp')

    plt = pv.Plotter()
    plt.add_mesh(tube, style='wireframe')
    plt.add_points(tubecenternodes, color='red')
    plt.show() """

    vmtkcenterline = pv.read('C:/Users/Yarran/OneDrive/Uni_Documenten/2023-24/BEP/VMTK/tubecenterline.vtp')
    ourcenterline = pv.read('C:/Users/Yarran/OneDrive/Uni_Documenten/2023-24/BEP/VMTK/ourcenterline.vtp')

    vmtkcenterline.plot()

    plt = pv.Plotter()
    plt.add_mesh(ourcenterline, color='red')
    plt.add_mesh(vmtkcenterline, color = 'green')
    plt.add_mesh(tube, style='wireframe')
    plt.show()

    dist_vmtk = np.linalg.norm(vmtkcenterline.find_closest_cell(points, return_closest_point=True)[1]-points, axis=1)
    dist_slicing = np.linalg.norm(ourcenterline.find_closest_cell(points, return_closest_point=True)[1]-points, axis=1)

    rel_vmtk = dist_vmtk/(spline['radius']*scale)
    rel_slicing = dist_slicing/(spline['radius']*scale)

    # Original indices (0 to 99)
    original_indices = np.arange(100)
    # New indices (0 to 999)
    new_indices = np.linspace(0, 99, 1000)
    # Interpolated array of length 1000
    r_interp = np.interp(new_indices, original_indices, r)*scale
    print(tube.bounds)

    pyplt.plot(rel_vmtk[10:-10]* 100, label = 'vmtk error [% of radius]')
    pyplt.plot(rel_slicing[10:-10]* 100, label = 'custom alg. error [% of radius]')
    pyplt.plot(r_interp/10, label = 'turn radius [cm]')
    pyplt.xlabel('Datapoints along centerline')
    pyplt.legend()
    pyplt.show()
    return

def aorta_geom():
    wall = pv.read("Aortic_CFD_workflow/geometries/case_0002/meshes/wall.stl")
    inlet = pv.read("Aortic_CFD_workflow/geometries/case_0002/meshes/inlet.stl")
    vmtk_centerline= pv.read('C:/Users/Yarran/OneDrive/Uni_Documenten/2023-24/BEP/VMTK/centerline.vtp')
    centerlinedata = centerline(inlet, wall, dist=1000)
    centernodes = centerlinedata[0]
    areas = centerlinedata[3]
    plt = pv.Plotter()
    plt.add_mesh(wall, style='wireframe')
    plt.add_mesh(vmtk_centerline, color='green')
    plt.add_points(centernodes, color='red')
    plt.show()

    dist = np.linalg.norm(vmtk_centerline.find_closest_cell(centernodes, return_closest_point=True)[1]-centernodes, axis=1)
    dist_rel = dist/((areas/np.pi)**0.5)
    pyplt.plot(dist_rel[15:-10]*100, label = 'custom alg.')
    pyplt.xlabel('Datapoints along centerline')
    pyplt.ylabel('Error [% of radius]')
    pyplt.legend()
    pyplt.show()
    return

ref_geom()
aorta_geom()