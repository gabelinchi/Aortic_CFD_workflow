import numpy as np
import pyvista as pv
import utils as ut
import cutting_with_vector as cwv

#import the inlet and the wall with pyvista
inlet = pv.read("geometries\input\inlet.stl")
wall = pv.read("geometries\input\wall.stl")

#main function that runs the cutter script
def main_cutter(inlet, wall, plot=False):

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

    print('Smallest cross-section: ',smallest_area)

    #Calculates the normal in the final cutting plane centerpoint
    normal_final = ut.average_normal(centernormals, 4, smallest_area_index)
    #Final centernode
    center_final = centernodes[smallest_area_index].flatten()

    print(normal_final)
    print(center_final)
    new_geometry = cwv.cut(center_final, normal_final, wall, plot=True)    
    return new_geometry

def centerline(inlet, wall):
    #From the inlet surface mesh extract the boundary edges
    inlet_boundary = inlet.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)

    #Calculate the centernode of the inlet mesh and the normal of the centernode
    inlet_points = inlet_boundary.points
    inlet_centerpoint = inlet_points.mean(0)
    inlet_normal = inlet.compute_normals()['Normals'].mean(0)

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
    while count < 200:

        #Create a new point by adding the directional vector to the centernode of the previous edgeprofile (previous iteration)
        inter_center = np.add(center, normal)
        #Use this new point and the previous directional vector to make a cut of the wall mesh

        inter_profile = cwv.get_clip_perimeter(inter_center, normal, wall) #Extract the edge profile (temporary) and isolate the profile (with connectivity) which we want to continue working with
        #From the isolated profile (temporary), calculate the new center point
        inter_profile_points = inter_profile.points
        new_center = inter_profile_points.mean(0)    
        #Calculate the normalised relative vector from the centernode of the previous edgeprofile (previous iteration) and the new center point
        A = np.subtract(inter_center, center)
        B =  0.1 * np.subtract(new_center, inter_center)
        new_normal = ut.normalise(np.add(A, B))

        #Make a new cut of the wall mesh with the new center point and the new normalised directional vector (also isolate this edgeprofile again)
        new_profile = cwv.get_clip_perimeter(new_center, new_normal, wall)

        #Calculate the area of the new profile
        new_profile_closed = new_profile.delaunay_2d()
        new_area = new_profile_closed.area

        #Store the new center point, the normalised directional vector, the last (second) edge profile and the area of the last edge profile in the storage variables
        centernodes = np.vstack([centernodes, new_center])
        centernormals = np.vstack([centernormals, new_normal])
        edgeprofiles += new_profile
        slice_areas = np.append(slice_areas, new_area)

        #Calculate the distance of the new center point and the center point of the inlet, if the distance is larger than a certain threshold. We break the while loop
        if np.linalg.norm(new_center - inlet_centerpoint) >= 40:
            break

        #Make the calculation variables the new center point and the new normalised directional vector
        center = new_center
        normal = new_normal

        count += 1
    
    return centernodes, centernormals, edgeprofiles, slice_areas

def areaselection(areas):
    # Find the indices where the tube transitions from wide to narrow and narrow to wide
    diff_areas = np.diff(areas)
    transition_indices = np.where(np.logical_or(diff_areas < 0, diff_areas == 0))[0]

    #Find the areas between these transition points
    area_segments = np.split(areas, transition_indices + 1)

    #Select the segments that have a narrow part and are larger than 1. This result in multiple arrays that go from a smaller to a larger area
    narrow_segments = [segment for segment in area_segments if len(segment) > 1 and segment[0] < segment[-1]]

    #Select the smallest area from the last small segment (might want to find a smarter solution for this, but for now this works)
    smallest_area = min(narrow_segments[1]) #Gives out of bound error when there is no constriction

    smallest_index_calc = np.where(smallest_area == areas)
    smallest_area_index = smallest_index_calc[0]
    return smallest_area, smallest_area_index


#Select the centernode and the directional vector for the final cut, based on area criteria (this is more refined, but do this after prove of concept for the previous stage)

geometry = main_cutter(inlet, wall, plot=True)

geometry.plot()


#Do the final cut of the wall geometry



