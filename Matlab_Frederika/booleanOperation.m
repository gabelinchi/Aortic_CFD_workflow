function [F_out,V_out] = booleanOperation(FA_in,VA_in,FB_in,VB_in,operation)
%% 
% Eric Trudel (2023). Surface Booleans (https://www.mathworks.com/matlabcentral/fileexchange/122502-surface-booleans), MATLAB Central File Exchange. Retrieved May 16, 2023.
% This function performs boolean operation for two close surfaces.
% input:
% FA_in: faces of 1st surface. 
% VB_in: vertices of 1st surface. 
% FB_in: faces of 2nd surface.
% VB_in: vertices of 2nd surface.
% operation: option of boolean operation.
% subtraction (union - 2nd Shape) --> operation = 1
% subtraction (union - 1nd Shape) --> operation = 2
% intersection --> operation = 3
% union --> operation = 4
%
% output:
% F_out,V_out: faces and vertices data of output structure.
%
% notes:
% install all function in Surface Booleans before using this function.
% https://www.mathworks.com/matlabcentral/fileexchange/122502-surface-booleans
%%
%turn off warnings
warning('off','all')
%tolerance for nodes
tolerance = 10^-8; %floating point tolerance
visual=0; %plot the floodfill algorithm? (your choice 1=yes,0=no)

%contruct the structure to be used for surface booleans
FV_A.faces = FA_in;
FV_A.vertices = VA_in;

FV_B.faces = FB_in;
FV_B.vertices = VB_in;

% plot the shapes together, ensure normals are pointing in the same direction
% figure(3)
% clf
% hold on
% view(3)
% axis equal
% Triang_A = triangulation(FV_A.faces,FV_A.vertices);
% trisurf(Triang_A,'edgecolor','k','Facecolor','y')
% F = featureEdges(Triang_A,pi/10);
% F2 = faceNormal(Triang_A);
% centroids = (Triang_A.Points(Triang_A.ConnectivityList(:,1),:)+Triang_A.Points(Triang_A.ConnectivityList(:,2),:)+Triang_A.Points(Triang_A.ConnectivityList(:,3),:))/3;
% quiver3(centroids(:,1),centroids(:,2),centroids(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'Color','b'); %check the face normals, must be consistent
% 
% % axis equal
% Triang_B = triangulation(FV_B.faces,FV_B.vertices);
% trisurf(Triang_B,'edgecolor','k','Facecolor','g')
% F = featureEdges(Triang_B,pi/10);
% F2 = faceNormal(Triang_B);
% centroids = (Triang_B.Points(Triang_B.ConnectivityList(:,1),:)+Triang_B.Points(Triang_B.ConnectivityList(:,2),:)+Triang_B.Points(Triang_B.ConnectivityList(:,3),:))/3;
% quiver3(centroids(:,1),centroids(:,2),centroids(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'Color','b'); %check the face normals, must be consistent
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check Mesh Quality %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%test if the input meshes are closed or not

logC = test_mesh_watertight(FV_B,1);
logS = test_mesh_watertight(FV_A,1);

%must loop for each free volume in "Sphere"
if ~logC || ~logS
    fprintf('Warning! Initial Meshes are not Watertight!\n');
end

[intSurface_tmp,FV_A_tmp, FV_B_tmp,sur_A_int,sur_B_int,edges_co_planar_tmp ,co_planar_switch] = sectioned_surface_intersection(FV_A,FV_B,tolerance);

%cleaned intersection line segments
edge_list_init = intSurface_tmp.edges; %edges seem to be directed properly if all the triangles have the same normals

%%
%%%%%%%%%%%%%%%%%%%
%%%Clean Meshes %%%
%%%%%%%%%%%%%%%%%%%

%clean meshes before connecting intersection edges together
[surfaceA, surfaceB] = pre_process_mesh(intSurface_tmp,FV_A, FV_A_tmp ,FV_B, FV_B_tmp,sur_A_int,sur_B_int, edge_list_init,edges_co_planar_tmp,tolerance);

%%
%FROMING INTERSECTION LOOPS%
fprintf('Forming Intersection Loops\n');
%for soft contour loops (intersecting contours)
%the loops must be split into two parts
%project lines to have and 'x' shape
all_edges = unique([surfaceA.edges; surfaceB.edges],'rows','stable');
[edge_loops,edge_loops_directed]  = form_intersection_loops_boolean(all_edges);


%%
%CREATING SUB SURFACES%
[m_aBlocks, orientation_matrix_tot ] = create_sub_blocks(surfaceA,surfaceB,edge_loops,co_planar_switch,visual);

%%
%Update the boolean when co-planar is turned on
if co_planar_switch
     
     m_aBlocks{1}.iso_bool(sum(orientation_matrix_tot,2)==0,:) = 0; 
     m_aBlocks{2}.iso_bool(sum(orientation_matrix_tot,2)==0,:) = 0;
end
%%

%Determine if the subblocks are union, subtraction or intersection%

%subtraction (union - 2nd Shape) --> operation = 1
%subtraction (union - 1nd Shape) --> operation = 2
%intersection --> operation = 3
%union --> operation = 4


new_nodes = uniquetol(surfaceA.vertices,tolerance,'Byrows',true);

    
tris_block = m_aBlocks{operation}.trias;
nodes_block = surfaceA.vertices(tris_block(:),:);
[~,new_ids] = ismembertol(nodes_block,new_nodes,tolerance,'Byrows',true );
tris_block(:) = new_ids;
F_out = tris_block;
V_out = surfaceA.vertices;
[F_out,V_out]=patchCleanUnused(F_out,V_out);

% figure
% upper = triangulation(F_out,V_out);
% trisurf(upper,'facecolor','r','edgecolor','k')
% view(45,25)
% axis equal


end