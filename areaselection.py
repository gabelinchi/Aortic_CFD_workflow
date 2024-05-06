import numpy as np

#Function that return the area and the index of the smallest area in the  slice_area array. Written with help of chatGPT
#Possible new feature: Selection of specific small part is now always the last one that occurs, might want to automatically define this
areas = np.array([10, 20, 30, 40, 50, 30, 20, 10, 15, 20, 25, 30, 40])

def areaselection(areas):
    # Find the indices where the tube transitions from wide to narrow and narrow to wide
    diff_areas = np.diff(areas)
    transition_indices = np.where(np.logical_or(diff_areas < 0, diff_areas == 0))[0]

    #Find the areas between these transition points
    area_segments = np.split(areas, transition_indices + 1)

    #Select the segments that have a narrow part and are larger than 1. This result in multiple arrays that go from a smaller to a larger area
    narrow_segments = [segment for segment in area_segments if len(segment) > 1 and segment[0] < segment[-1]]

    #Select the smallest area from the last small segment (might want to find a smarter solution for this, but for now this works)
    smallest_area = min(narrow_segments[-1])

    smallest_index_calc = np.where(smallest_area == areas)
    smallest_area_index = smallest_index_calc[0]
    return smallest_area, smallest_area_index
