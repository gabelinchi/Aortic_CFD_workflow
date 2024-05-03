import numpy as np
import pyvista as pv

#import the inlet and the wall with pyvista

#From the inlet surface mesh extract the boundary edges

#Calculate the centernode of the inlet mesh and the normal of the centernode

#Calculate the area of the inlet edgeprofile

#Create empty storage variables, for centernodes/centernormals/edgeprofiles/areas

#Use the inlet information calculated above as initial values in the storage variables (append) and and make the centernode and normal equal to a calculation variable

#Creat count variable starting at 0

#Create a while loop with a max count around 100 or something
    #In while loop
    #Create a new point by adding the directional vector to the centernode of the previous edgeprofile (previous iteration)
    #Use this new point and the previous directional vector to make a cut of the wall mesh
    #Extract the edge profile (temporary) and isolate the profile (with connectivity) which we want to continue working with
    #From the isolated profile (temporary), calculate the new center point
    #Calculate the relative vector from the centernode of the previous edgeprofile (previous iteration) and the new center point
    #Normalise this vector and make a new cut of the wall mesh with the new center point and the new normalised directional vector (also isolate this edgeprofile again)
    #Store the new center point, the normalised directional vector, the last (second) edge profile and the area of the last edge profile in the storage variables
    #Make the calculation variables the new center point and the new normalised directional vector
    #Calculate the distance of the new center point and the center point of the inlet, if the distance is larger than a certain threshold. We break the while loop

#Select the centernode and the directional vector for the final cut, based on area criteria (this is more refined, but do this after prove of concept for the previous stage)
#Do the final cut of the wall geometry



