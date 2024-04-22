%% From arterial lumen .stl to fluid and solid mesh to FEBio input file

% Starting point: binary .stl files of the fluid domain, wall domain and calcifications 
%The surfaces (.stl) of the components must NOT intersect, be sure this does not
%happen, and eventually correct this in the software used for segmentation. 

%%
clear; close all; clc;
%path with all functions
addpath('C:\Users\fffontana\OneDrive - Delft University of Technology\Comp modelling and simulations\GIBBON implementation\ICA\functions');

%%
% Plot settings
markerSize=50;
edgeThickness=3;
fontSize=25; 

%% Import a binary type STL file as patch data and create the struct for lumen and calcifications

%Set main folder
path = ['C:\Users\fffontana\OneDrive - Delft University of Technology\Comp modelling and simulations\GIBBON implementation\ICA\'];
mkdir(path)
fileName='EMC029_LICA_lumen.stl'; %could be changed in order to perform batch analyses
[stlStruct] = import_STL_bin(fileName); %careful, stl must be in binary format (not ASCII) - another function 


mesh_struct.lumen.F=stlStruct.solidFaces{1};
mesh_struct.lumen.V=stlStruct.solidVertices{1};
%mesh_struct.lumen.N=stlStruct.solidNormals{1}; %normals to the faces

fileName='EMC029_LICA_wall.stl'; %could be changed in order to perform batch analyses
[stlStruct] = import_STL_bin(fileName);
mesh_struct.wall.F=stlStruct.solidFaces{1};
mesh_struct.wall.V=stlStruct.solidVertices{1};


% !!! add a way to find in advance how many calcium there are !!!
num_calcium = {'1'};

for i=1:size(num_calcium,2)
 
    fileName=['EMC029_LICA_calcium',char(num_calcium(i)),'.stl'];
    [stlStruct] = import_STL_bin(fileName); %careful, stl must be in binary format (not ASCII)
    mesh_struct.calcium{i}.F=stlStruct.solidFaces{1};
    mesh_struct.calcium{i}.V=stlStruct.solidVertices{1};
    %mesh_struct.calcium{i}.N=stlStruct.solidNormals{1};

end


%% Merging nodes and translating geometry

[mesh_struct.lumen.F,mesh_struct.lumen.V]=mergeVertices(mesh_struct.lumen.F,mesh_struct.lumen.V);
[mesh_struct.wall.F,mesh_struct.wall.V]=mergeVertices(mesh_struct.wall.F,mesh_struct.wall.V);

for i=1:size(num_calcium,2)

    [mesh_struct.calcium{i}.F,mesh_struct.calcium{i}.V]=mergeVertices(mesh_struct.calcium{i}.F,mesh_struct.calcium{i}.V);

end

Vc = triSurfCentroid(mesh_struct.lumen.F,mesh_struct.lumen.V); %centroid of the geometry
mesh_struct.lumen.V = mesh_struct.lumen.V - Vc; %origin is now coincient with the centroid of the structure
mesh_struct.wall.V = mesh_struct.wall.V - Vc;
for i=1:size(num_calcium,2)

   mesh_struct.calcium{i}.V = mesh_struct.calcium{i}.V - Vc;

end

Vc = triSurfCentroid(mesh_struct.lumen.F,mesh_struct.lumen.V); %getting the new coordinates

%% Plotting the .stl 
cFigure;
title('Imported patch data from STL','fontSize',fontSize);
gpatch(mesh_struct.lumen.F,mesh_struct.lumen.V,'r');
gpatch(mesh_struct.wall.F,mesh_struct.wall.V,'g');
for i=1:size(num_calcium,2)
  gpatch(mesh_struct.calcium{i}.F, mesh_struct.calcium{i}.V,'b');
end
axisGeom(gca,fontSize);
alpha(0.3); %transparency to check eventual intersection between components
camlight('headlight');
lighting phong; %axis off; 
drawnow;

%% Identify the surface
%Angular threshold in radians - for surfaces detection
a=(30/180)*pi;
%Detect surface features
G=patchFeatureDetect(mesh_struct.lumen.F,mesh_struct.lumen.V,a); %through this function, the number of surfaces is detected (and each surface is labeled as 1,2,3...). Also, each face is assigned to a certain surface, being associated with the label 1, 2, 3, ...     
num = zeros(max(G),1); %coloumn vector with number entries = number of surfaces (= max value you can find in G)
for i = 1:max(G)
    num(i) = sum(ismember(G,i)); %storing in num(i) how many faces belong to the i-th surface
end
num_sort = sort(num,'descend');

cFigure;
gpatch(mesh_struct.lumen.F,mesh_struct.lumen.V,G);hold on
icolorbar;
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
view(45,45)
axis off;
%through this visualization we learn that label 1 is assigned to the wall, while 2 and 3 to the inlets

F_fs = mesh_struct.lumen.F(G==2,:);

%don't need it, may cancel these lines
% F_outlet = mesh_struct.lumen.F(G==2,:);
% F_inlet = mesh_struct.lumen.F(G==3,:);

%I think it's not needed:
[C]=patchConnectivity(mesh_struct.lumen.F,mesh_struct.lumen.V); %info on connectivity between edges, vertices and faces - C is a struct


%% Remesh the fluid-structure surface and compute new surfaces for inlet and outlet

%get the vertices associated with F_inner only
[F_fs_old,V_fs_old]=patchCleanUnused(F_fs,mesh_struct.lumen.V); %function that removes unused vertices from patch data - We have info on the faces in the internal surface (F_inner) but not the corresponding vertixes

F_fs = F_fs_old;
V_fs = V_fs_old;

% % !!! this remeshing could be avoided !!! - 
% optionStruct.pointSpacing=0.5; %Set desired point spacing
% optionStruct.disp_on=0; % Turn off command window text display
% % optionStruct.pre.max_hole_area=100; %Max hole area for pre-processing step
% % optionStruct.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
% % optionStruct.anisotropy=1; %Use anisotropy (~=0) to capture geometry or favour isotropic triangles (=0)
% [F_fs,V_fs]=ggremesh(F_fs_old,V_fs_old,optionStruct); 
% 
F_set={F_fs}; %F_fs in cell format
V_set={V_fs}; %V_fs in cell format
% 
% cFigure;
% subplot(1,2,1); hold on;
% title('Input mesh');
% gpatch(F_fs_old,V_fs_old,'w','k');
% axisGeom; axis off; 
% camlight headlight;
% 
% subplot(1,2,2); hold on;
% title('Geogram remeshed');
% gpatch(F_fs,V_fs,'gw','k',1,1);
% axisGeom; axis off; 
% camlight headlight;
% 
% gdrawnow;

%compute inlet and outlet surfaces:
Eb = detect_mesh_holes_and_boundary(F_fs); %detects the vertices that belong to the boundary in coincidence to holes
cFigure;
F_set_section={};
V_set_section={};
for i = 1:3 %creating the inlets and outlets for the fluid (no solid yet), storing them in F_set and V_set (for cycle with i=1:3 because the operation is repeated for the 3 sections) 
pointSpacing=0.3; %Desired point spacing
resampleCurveOpt=0;
interpMethod='linear'; %or 'natural'
boundaries=V_fs(Eb{i},:);
[F_support,V_support]=regionTriMesh3D({boundaries},pointSpacing,resampleCurveOpt,interpMethod);
F_set{end+1}=F_support;
V_set{end+1}=V_support;
F_set_section{end+1}=F_support;
V_set_section{end+1}=V_support;
plotV(boundaries,'.b','MarkerSize',10);hold on
gpatch(F_fs,V_fs,'r');hold on
gpatch(F_set_section{i},V_set_section{i},'g');hold on
end
axis equal
view(45,45)
%Fset and Vset are structs in which - until now - are stored the vertices
%and faces of the f-s surface the 3 sections separately, while in
%Fset_section and Vset_section the sections only (to facilitate the
%numbering (i goes from 1 to 3, not from 2 to 4 as it should considering
%Fset and Vset)

%%
% Identify which one is the inlet and which is the outlet -> criterion:
% outlet has higher z coordinate
min_z = inf; %inf = positive infinite
for i=1:3 %in 1 there's the f-s surface
    [~,V_support] = patchCleanUnused(F_set_section{i},V_set_section{i}); %process repeated for each cell (= each section)
    z_coords = V_support(:, 3); %isolate the z coordinates  
    
    min_z_part = min(z_coords);    
    % Update minimum z-coordinate if necessary 
    if min_z_part < min_z
        min_z = min_z_part;
        min_z_vertex_index = i;
    end
end
F_inlet = F_set_section{min_z_vertex_index};
V_inlet = V_set_section{min_z_vertex_index};

ind_outlet=1:3;
ind_outlet=ind_outlet(ind_outlet~=min_z_vertex_index);

%identify which outlet is ACA and which outlet is MCA: ACA has the max x coordinate
max_x = 0;
for i=1:2 %because we have 2 outlets
[~,V_support] = patchCleanUnused(F_set_section{ind_outlet(i)},V_set_section{ind_outlet(i)});
x_coords = V_support(:, 1); %isolate the x coordinates  
max_x_part = max(x_coords);    
    % Update minimum z-coordinate if necessary 
    if max_x_part > max_x
        max_x = max_x_part;
        max_x_vertex_index = ind_outlet(i);
    end

end

F_outlet_ACA = F_set_section{max_x_vertex_index};
V_outlet_ACA = V_set_section{max_x_vertex_index};

ind_MCA=1:3;
ind_MCA=ind_MCA(ind_MCA~=min_z_vertex_index);
ind_MCA=ind_MCA(ind_MCA~=max_x_vertex_index);


F_outlet_MCA = F_set_section{ind_MCA};
V_outlet_MCA = V_set_section{ind_MCA};

cFigure;
gpatch(F_fs,V_fs,'kw','none',0.2);hold on
gpatch(F_inlet,V_inlet,'r');hold on
gpatch(F_outlet_ACA,V_outlet_ACA,'b');hold on
gpatch(F_outlet_MCA,V_outlet_MCA,'g');hold on
title('Surfaces before cutting','fontSize',fontSize);
axis equal
view(45,45)
legend('Solid-Fluid Interface','Inlet','ACA','MCA','FontSize',15)

%[~,V_inlet] = patchCleanUnused(F_inlet,V_inlet); - %not useful, doesn't clean any vertix


%% Cut the sections - inlet and outlet

V_wall = mesh_struct.wall.V;
F_wall = mesh_struct.wall.F;

[F_new,V_new,F_new_wall,V_new_wall] = cutVessel_bifurcation_wall(F_fs,V_fs,F_wall,V_wall,F_inlet,V_inlet,F_outlet_MCA,V_outlet_MCA,F_outlet_ACA,V_outlet_ACA);

%just for backup of what comes out of the cutVessel function
F_backup_lumen = F_new; 
V_backup_lumen = V_new;
F_backup_wall = F_new_wall;
V_backup_wall = V_new_wall;

%%
%already remeshed, in the cutVessel function
F_fs = F_new;
V_fs = V_new;

F_wall = F_new_wall;
V_wall = V_new_wall;

cFigure;
gpatch(F_fs,V_fs);
gpatch(F_wall,V_wall,'kw',0.3)
title('lumen and wall cut and remeshed','fontSize',25);
view(45,25)
% axis off;
axis equal;
%close(gcf)

%% Add extensions

optionStruct2.extendDistance=5; %Extend distance
optionStruct2.extendMethod=4; %Method to use for extruding - see gdoc for details

Eb=patchBoundary(F_fs,V_fs);
VE=patchCentre(Eb,V_fs); %Edge centre coordinates

Eb_wall=patchBoundary(F_wall,V_wall);
VE_wall=patchCentre(Eb_wall,V_wall);

%attempt
% boundaries = [];
% for i=1:3
% boundaries = [boundaries; V_fs(Eb{i},:)];
% end

%Questa strategia non va bene, per il tipo di dato che si ottiene da detect_mesh_holes_and_boundary, ovvero gli indici dei nodi che
%appartengono agli edges - mentre a noi per fare l'extrusion serve una matrice n x 2 con gli indici delle coppie di nodi che formano gli edges
% Eb = detect_mesh_holes_and_boundary(F_fs); %All boundary edges (inlet and outlets)
% VE = [V_fs(Eb{1,1},:); V_fs(Eb{2,1},:); V_fs(Eb{3,1},:)];
% Eb = [Eb{1,1}'; Eb{2,1}'; Eb{3,1}'];

cFigure
plotV(VE,'.b','MarkerSize',10);hold on
plotV(VE_wall,'.b','MarkerSize',10);hold on
axisGeom; axis off; 
camlight headlight;


%logicTop_out=VE(:,3)>0; %Logic for outlet edges (Z-coordinate above 0)
logicTop_in=VE(:,3)<0; %Logic for inlet edges (Z-coordinate below 0)
%Eb_out=Eb(logicTop_out,:); %Set of edges at the outlet
Eb_in=Eb(logicTop_in,:); %Set of edges at the inlet
logicTop_ACA= VE(:,1)>0 & VE(:,3)>0; %Logic for outlets is z>0, for ACA x>0
Eb_ACA=Eb(logicTop_ACA,:);
logicTop_MCA=VE(:,1)<0 & VE(:,3)>0; %Logic for outlets is z>0, for MCA x<0
Eb_MCA=Eb(logicTop_MCA,:);


[F_long_in,V_long_in]=patchExtend(F_fs,V_fs,Eb_in,optionStruct2);
[F_long_ACA,V_long_ACA]=patchExtend(F_fs,V_fs,Eb_ACA,optionStruct2);
[F_long_MCA,V_long_MCA]=patchExtend(F_fs,V_fs,Eb_MCA,optionStruct2);

cFigure
title(['Extensions added']);
gpatch(F_fs,V_fs,'rw','k',1); hold on;
gpatch(F_long_ACA,V_long_ACA,'bw','k',1);
gpatch(F_long_MCA,V_long_MCA,'bw','k',1);
gpatch(F_long_in,V_long_in,'bw','k',1);
axisGeom; axis off; 
camlight headlight;

%the 3 separate components are put in the same cell array so to be joined through the joinElementSets function
F_fs = {F_fs, F_long_ACA, F_long_MCA, F_long_in};
V_fs = {V_fs, V_long_ACA, V_long_MCA V_long_in};

[F_fs, V_fs, C_fs] = joinElementSets(F_fs,V_fs);
[F_fs,V_fs]=patchCleanUnused(F_fs,V_fs); %just to be sure no issue from the previous joining function


%Extensions of the solid wall:
logicTop_in=VE_wall(:,3)<0; %Logic for inlet edges (Z-coordinate below 0)
%Eb_out=Eb(logicTop_out,:); %Set of edges at the outlet
Eb_in=Eb_wall(logicTop_in,:); %Set of edges at the inlet
logicTop_ACA= VE_wall(:,1)>0 & VE_wall(:,3)>0; %Logic for outlets is z>0, for ACA x>0
Eb_ACA=Eb_wall(logicTop_ACA,:);
logicTop_MCA=VE_wall(:,1)<0 & VE_wall(:,3)>0; %Logic for outlets is z>0, for MCA x<0
Eb_MCA=Eb_wall(logicTop_MCA,:);


[F_long_in,V_long_in]=patchExtend(F_wall,V_wall,Eb_in,optionStruct2);
[F_long_ACA,V_long_ACA]=patchExtend(F_wall,V_wall,Eb_ACA,optionStruct2);
[F_long_MCA,V_long_MCA]=patchExtend(F_wall,V_wall,Eb_MCA,optionStruct2);

cFigure
title(['Extensions added']);
gpatch(F_wall,V_wall,'rw','k',1); hold on;
gpatch(F_long_ACA,V_long_ACA,'bw','k',1);
gpatch(F_long_MCA,V_long_MCA,'bw','k',1);
gpatch(F_long_in,V_long_in,'bw','k',1);
axisGeom; axis off; 
camlight headlight;

%the 3 separate components are put in the same cell array so to be joined through the joinElementSets function
F_wall = {F_wall, F_long_ACA, F_long_MCA, F_long_in};
V_wall = {V_wall, V_long_ACA, V_long_MCA V_long_in};

[F_wall, V_wall, C_wall] = joinElementSets(F_wall,V_wall);
[F_wall,V_wall]=patchCleanUnused(F_wall,V_wall); %just to be sure no issue from the previous joining function


% once the extensions are added, remesh before adding the BL -> otherwise it's like there's a discontinuity where 
% the extensions were added and after the BL are added, in that region holes in the mesh form.
optionStruct.pointSpacing=0.3; %Set desired point spacing
optionStruct.disp_on=0; % Turn off command window text display
[F_fs,V_fs]=ggremesh(F_fs,V_fs,optionStruct); %remesh after having added the extensions, to avoid issue in the mesh when adding boundary layers (the elements detatch otherwise)
[F_wall,V_wall]=ggremesh(F_wall,V_wall,optionStruct);


cFigure
gpatch(F_fs,V_fs,'w','k',1); hold on; %check the result of the joining operation
gpatch(F_wall,V_wall,'w','k',1); hold on;
axisGeom;
camlight headlight;


%% Obtain the arterial wall

%since it's built with patchThick, the solid domain will have a penta6 mesh, as the boundary layers
% [E_solid,V_solid,Fq1_solid,Fq2_solid] = patchThick(F_fs,V_fs,1,1.2,4); %thickness of each element = 1.5 - elements along the thickness = 8
% F_solid = element2patch(E_solid,[],'penta6'); %get the faces for patches visualization
% 
% cFigure;
% gpatch(F_solid,V_solid,'bw','k',0.5); hold on 
% gpatch(F_fs,V_fs,'r','k'); hold on
% gpatch(Fq2_solid, V_solid, 'g'); hold on
% legend('Vessel Wall','','Inner surface', 'External surface') %the second component of the legend is omitted because it correspondes to the vessel wall again, since FE_w is a cell array of 2 cells (one for tri and the other for quad faces), and the legend sees it as twice
% axis off; 
% axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
% %view(45,45)


cFigure;
gpatch(F_fs, V_fs, 'r'); hold on
gpatch(F_wall, V_wall, 'g'); hold on
for i=1:size(num_calcium,2)
  gpatch(mesh_struct.calcium{i}.F, mesh_struct.calcium{i}.V,'b');
end
legend('f-s surface','External surface', 'calcifications');
alpha(0.3);
axis off; 
axisGeom(gca, fontSize);


%% Mesh solid domain and calcifications - tetgen

F_solid_sup = {F_wall, F_fs};
V_solid_sup = {V_wall, V_fs};
[F_solid_sup, V_solid_sup] = joinElementSets(F_solid_sup,V_solid_sup);  %surface of the wall, both internal and external side
[F_solid_sup,V_solid_sup]=patchCleanUnused(F_solid_sup,V_solid_sup);

% optionStruct.pointSpacing=0.3; %Set desired point spacing
% optionStruct.disp_on=0; % Turn off command window text display
% [F_solid_sup,V_solid_sup]=ggremesh(F_solid_sup,V_solid_sup,optionStruct);

cFigure
gpatch(F_solid_sup, V_solid_sup, 'r'); hold on
axis off; 
axisGeom(gca, fontSize);

%Compute inlet and outlet surfaces:
Eb = detect_mesh_holes_and_boundary(F_solid_sup); %detects the vertices that belong to the boundary in coincidence to holes

boundaries={};
for i=1:6
boundaries{end+1} = V_solid_sup(Eb{i},:);
end

boundary_inlet={};
boundary_outlet1={};
boundary_outlet2={};


for i=1:6
if boundaries{i}(:,3) < 0
    boundary_inlet{end+1} = boundaries{i};
elseif (boundaries{i}(:,3) > 0) & (boundaries{i}(:,1) > 0)
        boundary_outlet1{end+1} =  boundaries{i};
else 
      boundary_outlet2{end+1} =boundaries{i};
end
end

%check identification is correct
cFigure
plotV(boundary_inlet{1},'.b','MarkerSize',10);hold on
plotV(boundary_inlet{2},'.g','MarkerSize',10);hold on
plotV(boundary_outlet1{1},'.b','MarkerSize',10);hold on
plotV(boundary_outlet1{2},'.g','MarkerSize',10);hold on
plotV(boundary_outlet2{1},'.b','MarkerSize',10);hold on
plotV(boundary_outlet2{2},'.g','MarkerSize',10);hold on
axis equal
view(45,45)

boundaries = {boundary_inlet, boundary_outlet1, boundary_outlet2};

cFigure;
F_set={F_solid_sup};
V_set={V_solid_sup};
for i = 1:3 %creating the inlets and outlets for the solid, storing them in F_set and V_set (for cycle with i=1:3 because the operation is repeated for the 3 sections) 
pointSpacing=0.3; %Desired point spacing
resampleCurveOpt=0;
interpMethod='linear'; %linear or 'natural'
%boundaries=[V_solid_sup(Eb{i},:);V_solid_sup(Eb{i+1},:)];
[F_support,V_support]=regionTriMesh3D(boundaries{i},pointSpacing,resampleCurveOpt,interpMethod);
F_set{end+1}=F_support;
V_set{end+1}=V_support;
plotV(boundaries{i},'.b','MarkerSize',10);hold on
gpatch(F_solid_sup,V_solid_sup,'r');hold on
gpatch(F_set{i+1},V_set{i+1},'g');hold on
end
axis equal
view(45,45)

%check solid surface mesh is complete and capped:
cFigure
gpatch(F_set, V_set, 'r'); hold on
axis off; 
axisGeom(gca, fontSize);

[F_solid_sup,V_solid_sup,C] = joinElementSets(F_set,V_set);
[F_solid_sup,V_solid_sup]=patchCleanUnused(F_solid_sup,V_solid_sup);

cFigure
gpatch(F_solid_sup, V_solid_sup, 'r'); hold on
axis off; 
axisGeom(gca, fontSize);

F_set2 = {F_solid_sup};
V_set2 = {V_solid_sup};
V_regions(1,:)= getInnerPoint(F_solid_sup,V_solid_sup);
regionTetVolumes(1) = tetVolMeanEst(F_solid_sup,V_solid_sup)*3;


%remeshing calcifications to obtain an homogeneous mesh of the calcif
optionStruct.pointSpacing=0.2; %Set desired point spacing
optionStruct.disp_on=0; % Turn off command window text display
cFigure;
for i=1:size(num_calcium,2)
  [mesh_struct.calcium{i}.F,mesh_struct.calcium{i}.V]=ggremesh(mesh_struct.calcium{i}.F,mesh_struct.calcium{i}.V,optionStruct);
  gpatch(mesh_struct.calcium{i}.F, mesh_struct.calcium{i}.V,'b');
end
legend('calcifications');
axis off; 
axisGeom(gca, fontSize);

for i = 1:size(num_calcium,2)
F = mesh_struct.calcium{i}.F;
V = mesh_struct.calcium{i}.V;
V_regions(1+i,:) = getInnerPoint(F,V);
regionTetVolumes(1+i) = tetVolMeanEst(F,V)*2;
F_set2{end+1} = F;
V_set2{end+1} = V;
end

%%
stringOpt='-pqAaY'; %Options for tetgen
[F,V,C] = joinElementSets(F_set2,V_set2);
[F,V]=patchCleanUnused(F,V);

%Create tetgen input structure
inputStruct.stringOpt=stringOpt; %Tetgen options
inputStruct.Faces=F; %Boundary faces
inputStruct.Nodes=V; %Nodes of boundary
inputStruct.faceBoundaryMarker=C;
inputStruct.regionPoints=V_regions; %Interior points for regions
inputStruct.holePoints=[]; %Interior points for holes
inputStruct.regionA=regionTetVolumes; %Desired tetrahedral volume for each region
inputStruct.minRegionMarker=1;
[meshOutput]=runTetGen(inputStruct); %Run tetGen
% temp = meshOutput.elementMaterialID;
% temp(meshOutput.elementMaterialID == 1)= -3*ones(sum(meshOutput.elementMaterialID==1),1);
% meshOutput.elementMaterialID = temp;
meshOutput.elementMaterialID = -meshOutput.elementMaterialID;
meshView(meshOutput,[]);
axis equal
view(0,0)

Fb = meshOutput.facesBoundary;
ind_Fb = meshOutput.boundaryMarker;
V_solid = meshOutput.nodes; %IMPORTANT
%faces:
F_solid=Fb(ind_Fb==1,:); 
F_calcium1 = Fb(ind_Fb==2,:);
% F_calcium2 = Fb(ind_Fb==3,:);
% F_calcium3 = Fb(ind_Fb==4,:);
% F_calcium4 = Fb(ind_Fb==5,:);

cFigure
gpatch(F_calcium1, V_solid, 'r'); hold on %for visualization, changed F_calcium2 with F_calcium1 and F_solid
axis off; 
axisGeom(gca, fontSize);

ind_El = abs(meshOutput.elementMaterialID);
E_solid = meshOutput.elements(ind_El==1,:); %IMPORTANT
E_calcium = meshOutput.elements(ind_El~=1,:); %IMPORTANT

%identify f-s surface:
a=(30/180)*pi;
%Detect surface features
G=patchFeatureDetect(F_solid,V_solid,a);     
num = zeros(max(G),1); %coloumn vector with number entries = number of surfaces (= max value you can find in G)
for i = 1:max(G)
    num(i) = sum(ismember(G,i)); %storing in num(i) how many faces belong to the i-th surface
end
num_sort = sort(num,'descend');

cFigure;
gpatch(F_solid,V_solid,G);hold on
icolorbar;
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
view(45,45)
axis off;
%through this visualization we learn how labels were assigned

F_fs = F_solid(G==2,:);


%% Fluid domain: Add boundary layers

% bias = 2 and with 5 segments (segment = layer)
n = 31; % 1+2+4+8+16
layerThickness = 0.2; %Yanjing used 0.5 - I guess it depends on the dimension of the sections, mine are very small and using 0.5 was too much
unit=layerThickness/n;

E_BL1 = {};
V_BL1 = {};
F_inner_temp = F_fs; %created temp matrixes for the internal surface, because as we create each BL, the internal surface changes (at the beginning, internal surface = fluid-structure surface)
V_inner_temp = V_solid;

% loop for each segment
for i = 1:5
    
    segmentThickness = (2^(i-1))*unit; %thickness of each segment is different -> smaller as closer the outer surface
    [E,V,Fq1,Fq2]=patchThick(F_inner_temp,V_inner_temp,-1,segmentThickness,1); %elements along the thickness = 1
    
    E_BL1{6-i} = E; %cell arrays to then join the sets
    V_BL1{6-i} = V; %cell arrays to then join the sets
    F_inner_temp = Fq2;
    V_inner_temp = V;
    
    %visualize for each cycle what is going on (so the layers/segments that are being added) 
    cFigure;
    F=element2patch(E,[],'penta6');
    gpatch(F_fs,V_fs,'r','none',0.3);hold on %plotting the starting point = f-s surface
    gpatch(F,V);hold on %plotting the layer added
    axis equal
    view(45,45)
    close(gcf)
end

[E_BL1,V_BL1,C] = joinElementSets(E_BL1,V_BL1);
[E_BL1,V_BL1]=patchCleanUnused(E_BL1,V_BL1);
E_temp = E_BL1;
% invert the element to avoid negative volume
E_BL1(:,1:3) = E_temp(:,4:6);
E_BL1(:,4:6) = E_temp(:,1:3);
F_BL1=element2patch(E_BL1,[],'penta6'); %FE_BL1 is a cell array because the faces of penta6 elements are some triangles, other quadrilaterals, so the matrixes have a different size and can't be combined in one only

cFigure;
%gpatch(F_fs,V_fs,'kw','none',0.3);hold on %f-s surface = outer surface of the fluid domain / inner surface of the solid domain
gpatch(F_BL1,V_BL1);hold on %vertices/faces of the whole structure we built with the boundary layers
gpatch(F_inner_temp,V_inner_temp,'b');hold on %inner surface
legend('Boundary Layers','','Inner surface')
axis off; 
axis equal
view(45,45)

cFigure;
gpatch(F_inner_temp,V_inner_temp,'b');
axis off; 
axis equal
view(45,45)


%% Cap inlet / outlet fluid domain

%Once again, compute inlet and outlet surfaces:
Eb = detect_mesh_holes_and_boundary(F_inner_temp); %detects the vertices that belong to the boundary in coincidence to holes
cFigure;
F_set={F_inner_temp};
V_set={V_inner_temp};
for i = 1:3 %creating the inlets and outlets for the fluid (no solid yet), storing them in F_set and V_set (for cycle with i=1:2 because the operation is repeated for the 2 sections) 
pointSpacing=0.3; %Desired point spacing
resampleCurveOpt=0;
interpMethod='linear'; %or 'natural'
boundaries=V_inner_temp(Eb{i},:);
[F_support,V_support]=regionTriMesh3D({boundaries},pointSpacing,resampleCurveOpt,interpMethod);
F_set{end+1}=F_support;
V_set{end+1}=V_support;
plotV(boundaries,'.b','MarkerSize',10);hold on
gpatch(F_fs,V_fs,'r');hold on
gpatch(F_set{i+1},V_set{i+1},'g');hold on
end
axis equal
view(45,45)

F_fluid_sup = F_set;
V_fluid_sup = V_set;

cFigure;
gpatch(F_fluid_sup,V_fluid_sup,'r'); hold on
axis off; 
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
%view(45,45)

%% Mesh fluid domain - tetgen

[F_fluid,V_fluid,C]=joinElementSets(F_set,V_set);
[F_fluid,V_fluid]=patchCleanUnused(F_fluid,V_fluid);
V_regions= getInnerPoint(F_fluid,V_fluid);
regionTetVolumes = tetVolMeanEst(F_fluid,V_fluid)*3;

stringOpt='-pqAaY'; %Options for tetgen

%Create tetgen input structure
inputStruct.stringOpt=stringOpt; %Tetgen options
inputStruct.Faces=F_fluid; %Boundary faces
inputStruct.Nodes=V_fluid; %Nodes of boundary
inputStruct.faceBoundaryMarker=C;
inputStruct.regionPoints=V_regions; %Interior points for regions
inputStruct.holePoints=[]; %Interior points for holes
inputStruct.regionA=regionTetVolumes; %Desired tetrahedral volume for each region
inputStruct.minRegionMarker=1;
[meshOutput1]=runTetGen(inputStruct); %Run tetGen
% temp = meshOutput.elementMaterialID;
% temp(meshOutput.elementMaterialID == 1)= -3*ones(sum(meshOutput.elementMaterialID==1),1);
% meshOutput.elementMaterialID = temp;
meshOutput1.elementMaterialID = -meshOutput1.elementMaterialID;
meshView(meshOutput1,[]);
axis equal
view(0,0)

E_fluid = meshOutput1.elements;
V_fluid = meshOutput1.nodes;

%generated for future visualization
F_fluid = element2patch(E_fluid,[],'tet4');

[E_fluid_set,V_fluid_set,C_fluid_set]=joinElementSets({E_fluid,E_BL1},{V_fluid,V_BL1}); %generated cell array with both the tetrahedral fluid mesh and the pentahedral fluid mesh (the one in the BL) - Can't use a single matrix because the number of vertices is different for the 2 mesh types

%% Show the whole domain - fluid+solid

cFigure;
gpatch(F_solid,V_solid,'rw'); hold on
gpatch(F_fluid,V_fluid,'r'); hold on
gpatch(F_BL1,V_BL1,'r'); hold on
gpatch(F_calcium1, V_solid, 'g'); hold on
% gpatch(F_calcium2, V_solid, 'g'); hold on
% gpatch(F_calcium3, V_solid, 'g'); hold on
% gpatch(F_calcium4, V_solid, 'g'); hold on
axis off; 
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
view(45,45)
legend('Vessel wall','Blood','','','Calcifications','FontSize',fontSize)


% Combine all parts
[E_set,V_set,C_set]=joinElementSets({E_solid,E_calcium,E_fluid,E_BL1},{V_solid,V_solid,V_fluid,V_BL1});
%this way we have a common numbering for the vertices !! which is necessary for the input file definition and FE modelling. 
%Note that V_set is one big matrix with all the vertices of fluid/solid, while E_set is a cell array, keeping the division 
% of the 4 parts (solid, calcium, fluid and BL) - fluid and BL are separated because of missing consistency between elements type
[E_set,V_set]=mergeVertices(E_set,V_set); %to avoid having double nodes, for example in the f-s interface

%% Indentify f-s suface (FS interface / Non-slip Condition)

F = element2patch(E_set{4},[],'penta6'); %faces of the boundary layers domain
indBoundary=tesBoundary(F);
Fb=F{1}(indBoundary{1},:);

%Angular threshold in radians
a=(30/180)*pi;
%Detect surface features
G=patchFeatureDetect(Fb,V_set,a);
num = zeros(max(G),1);
for i = 1:max(G)
    num(i) = sum(ismember(G,i));
end
num_sort = sort(num,'descend');

% Used to check the numbering of the surfaces for inlet and outlet
cFigure;
gpatch(Fb,V_set,G);hold on
icolorbar;
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
view(45,45)

F_fs = Fb(G==2,:);

%invert the nodes order (facet winding) - issue encountered when importing in FEBio:
%nodes order for each face as performed by GIBBON/MatLab is not compliant with what FEBio requires
F_support = F_fs;
F_fs(:,1) = F_support(:,3);
F_fs(:,2) = F_support(:,2);
F_fs(:,3) = F_support(:,1);

cFigure;
gpatch(F,V_set,'kw');hold on
gpatch(F_fs,V_set,'r');hold on
axis equal
view(45,25)

% F = element2patch(E_set{3},[],'penta6'); %faces of the boundary layers domain
%[G]=tesgroup(F{1});
%F_fs=F{1}(G(:,end),:);
%alternative method to get F_fs, but I'm not sure of how it works -> in
%F{1} are stored all the tri elements of the BL domain, which are 2
%surfaces. How does G detect only the fs surface ones?

%F_fs is the matrix storing the indexes of the vertices of the f-s surface,
%using the numbering of the V_set (so the set that keeps into account the
%vertices of the whole domain) - of course F-fs has indexes for each face
%different than the global one, but that's ok since ids for element sets in
%the FEBio input file have their own numbering that start from 1 for each
%set - this logic holds true also for the following element sets

%on the solid side:
% F = element2patch(E_set{1},[],'penta6');
% 
% indBoundary=tesBoundary(F);
% Fb=F{1}(indBoundary{1},:);
% 
% %Angular threshold in radians
% a=(30/180)*pi;
% %Detect surface features
% G=patchFeatureDetect(Fb,V_set,a);
% num = zeros(max(G),1);
% for i = 1:max(G)
%     num(i) = sum(ismember(G,i));
% end
% num_sort = sort(num,'descend');
% 
% % Used to check the numbering of the surfaces for inlet and outlet
% cFigure;
% gpatch(Fb,V_set,G);hold on
% icolorbar;
% axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
% view(45,45)
% 
% F_fs_solid = Fb(G==1,:);


%% Indentify Inlet and Outlet solid domain (Fix Displacement)

FE = element2patch(E_set{1},[],'tet4');
indBoundary=tesBoundary(FE);
Fb=FE(indBoundary,:); %using the indexes info to isolate the boundary faces from FE

%Angular threshold in radians
a=(30/180)*pi;
%Detect surface feature
G=patchFeatureDetect(Fb,V_set,a);
num = zeros(max(G),1);
for i = 1:max(G)
    num(i) = sum(ismember(G,i));
end

cFigure;
gpatch(Fb,V_set,G);hold on
icolorbar;
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
view(45,45)

F_fix = [Fb(G==3,:);Fb(G==5,:);Fb(G==6,:)];

[~,V_fix] = patchCleanUnused(F_fix,V_set);
ind_fix = find(ismember(V_set,V_fix,'rows')); %in ins_fix are stored the indexes of the vertices that belong to F_fix, with the numbering of V_set
%Notice that V_fix of course has a local numbering, doesn't take into
%consideration the global numbering of vertices (of course, numbering goes
%from 1 to dim of the matrix) and if we compute F_fix through the
%patchCleanUnused, it would consider the local numbering inside V_fix,
%which is useless for the FEBio input file - may be useful only for
%visualization purposes in MatLab. 

cFigure;
gpatch(FE,V_set,'kw');hold on
gpatch(F_fix,V_set,'r');hold on
plotV(V_set(ind_fix,:),'.b','MarkerSize',30);hold on
axis equal
view(45,25)
legend('Vessel Wall','','Fix Displacement')


%% Identify Inlet and Outlet and BL fluid domain (Fix Displacement and fluid BC)

%Boundary layers:
FE = element2patch(E_set{4},[],'penta6');
indBoundary=tesBoundary(FE);
Fb=FE{2}(indBoundary{2},:);

%Angular threshold in radians
a=(30/180)*pi;
%Detect surface features
G=patchFeatureDetect(Fb,V_set,a);
num = zeros(max(G),1);
for i = 1:max(G)
    num(i) = sum(ismember(G,i));
end
num_sort = sort(num,'descend');

cFigure;
gpatch(Fb,V_set,G);hold on
icolorbar;
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
view(45,45)
%Notice from this figure whether each BL is recognized as a different surface (could be that elements are not in plane)

F_BL = Fb;
[~,V_BL] = patchCleanUnused(F_BL,V_set);
ind_BL = find(ismember(V_set,V_BL,'rows'));

%F_out_BL = [Fb(G==1,:);Fb(G==3,:);Fb(G==5,:);Fb(G==7,:);Fb(G==9,:)];
%F_in_BL = [Fb(G==2,:);Fb(G==4,:);Fb(G==6,:);Fb(G==8,:);Fb(G==10,:)];
F_ACA_BL = [Fb(G==3,:)];
F_MCA_BL = [Fb(G==2,:)];
F_out_BL = [F_ACA_BL;F_MCA_BL];
F_in_BL = [Fb(G==1,:)];
ind_fix = [ind_fix; ind_BL];
F_fix = {F_fix, F_BL};

%To verify in and out BL are identified correctly
% cFigure;
% gpatch(F_out_BL,V_set,'r');hold on
% gpatch(F_in_BL,V_set,'b');hold on
% axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
% view(45,45)

cFigure;
gpatch(FE,V_set,'kw');hold on
gpatch(F_fix,V_set,'r');hold on
plotV(V_set(ind_fix,:),'.b','MarkerSize',30);hold on
axis equal
view(45,25)
legend('Blood domain','','Fix Displacement')


%Fluid domain (excluded BL):
FE = element2patch(E_set{3},[],'tet4');
indBoundary=tesBoundary(FE);
Fb=FE(indBoundary,:);

%Angular threshold in radians
a=(30/180)*pi;
%Detect surface features
G=patchFeatureDetect(Fb,V_set,a);
num = zeros(max(G),1);
for i = 1:max(G)
    num(i) = sum(ismember(G,i));
end
num_sort = sort(num,'descend');

% Used to check the numbering of the surfaces for inlet and outlet
cFigure;
gpatch(Fb,V_set,G);hold on
icolorbar;
axisGeom(gca, fontSize); %scales the axis and axis don't scale as you rotate
view(45,45)

F_fluid_BC = [Fb(G==2,:);Fb(G==3,:);Fb(G==4,:)];
F_out_fluid = [Fb(G==3,:);Fb(G==4,:)];
F_in_fluid = Fb(G==2,:);
[~,V_fluid_BC] = patchCleanUnused(F_fluid_BC,V_set);
ind_fluid_BC = find(ismember(V_set,V_fluid_BC,'rows'));

ind_fix = [ind_fix; ind_fluid_BC];
ind_fix = ind_fix';


cFigure;
gpatch(FE,V_set,'kw');hold on
gpatch(F_fix,V_set,'r');hold on
gpatch(F_fluid_BC,V_set,'r');hold on
plotV(V_set(ind_fix,:),'.b','MarkerSize',30);hold on
axis equal
view(45,25)
legend('Vessel Wall','','Fix Displacement')


%% Create FSI model - input file for FEBio
%ref, type on Command Window "open HELP_febioStruct2xml"

% FEA control settings
numTimeSteps=1550; %Number of time steps desired
max_refs=5; %Max reforms
max_ups=50; %Set to zero to use full-Newton iterations
opt_iter=50; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=0; %Minimum time step size
dtmax=0.005; %Maximum time step size

clear febio_spec
[febio_spec]=febioStructTemplate;
%febio_spec version
febio_spec.ATTR.version='4.0';

%Module section
febio_spec.Module.ATTR.type='fluid-FSI';

% febio_spec.Loads.surface_load{2}.velocity.ATTR.lc='1';
% febio_spec.Loads.surface_load{2}.velocity.VAL=-250;

%Control section
febio_spec.Control.analysis='DYNAMIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=0.001;
febio_spec.Control.plot_zero_state=0;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.time_stepper.aggressiveness=1;
febio_spec.Control.solver.ATTR.type = 'fluid-FSI';
febio_spec.Control.solver.symmetric_stiffness='non-symmetric';
febio_spec.Control.solver.reform_each_time_step=0;
febio_spec.Control.solver.diverge_reform=0;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.max_residual=1e20;
febio_spec.Control.solver.min_residual=1e-10;
febio_spec.Control.solver.rtol=0.001;
febio_spec.Control.solver.vtol=0.001;
febio_spec.Control.solver.ftol=0.001;
febio_spec.Control.solver.predictor=0;
febio_spec.Control.solver.min_volume_ratio=0;
febio_spec.Control.solver.order=2;
febio_spec.Control.solver.qn_method.ATTR.type = 'Broyden';
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.plot_level='PLOT_MUST_POINTS';
febio_spec.Control.output_level='OUTPUT_MUST_POINTS';
%remove unused field, belonging to GIBBBON .feb file template, which is for the SOLID module (not FSI)
febio_spec.Control.solver = rmfield(febio_spec.Control.solver,'arc_length_scale');
febio_spec.Control.solver = rmfield(febio_spec.Control.solver,'arc_length');
febio_spec.Control.solver = rmfield(febio_spec.Control.solver,'logSolve');
febio_spec.Control.solver = rmfield(febio_spec.Control.solver,'gamma');
febio_spec.Control.solver = rmfield(febio_spec.Control.solver,'beta');
febio_spec.Control.solver = rmfield(febio_spec.Control.solver,'alpha');
%Steps to use a load curve for the dtmax -> must points
febio_spec.Control.time_stepper = rmfield(febio_spec.Control.time_stepper,'dtmax');
febio_spec.Control.time_stepper.dtmax.ATTR.lc='2';
febio_spec.Control.time_stepper.dtmax.VAL=dtmax;

%Constants - need to be in consistant units
% febio_spec.Globals.Constants.T=0;
% febio_spec.Globals.Constants.P=0;
% febio_spec.Globals.Constants.R=8.31446;
% febio_spec.Globals.Constants.Fc=96485.3;

%Material section
materialName1 = 'VesselWall';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Mooney-Rivlin';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.density=1*10^(-9);
febio_spec.Material.material{1}.k=4.16;
febio_spec.Material.material{1}.pressure_model='default';
febio_spec.Material.material{1}.c1=0.083;
febio_spec.Material.material{1}.c2=0;

materialName2 = 'Calcification';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Mooney-Rivlin';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1*10^(-9);
febio_spec.Material.material{2}.k=13.3;
febio_spec.Material.material{2}.pressure_model='default';
febio_spec.Material.material{2}.c1=1.6;
febio_spec.Material.material{2}.c2=0;

materialName3 = 'Fluid';
febio_spec.Material.material{3}.ATTR.name=materialName3;
febio_spec.Material.material{3}.ATTR.type='fluid-FSI';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.solid.ATTR.type='neo-Hookean';
febio_spec.Material.material{3}.solid.density=0;
febio_spec.Material.material{3}.solid.E=1*10^(-15);
febio_spec.Material.material{3}.solid.v=0;
febio_spec.Material.material{3}.fluid.density=1.06*10^(-9);
febio_spec.Material.material{3}.fluid.k=2200;
febio_spec.Material.material{3}.fluid.viscous.ATTR.type='Newtonian fluid';
febio_spec.Material.material{3}.fluid.viscous.kappa=0;
febio_spec.Material.material{3}.fluid.viscous.mu=3.5*10^(-9);

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_set,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V_set; %The nodel coordinates

% -> Elements
partName='VesselWall';
febio_spec.Mesh.Elements{1}.ATTR.name=partName; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{1}.ATTR.id=1;
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E_set{1},1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E_set{1}; %The element matrix
eleID_start=size(E_set{1},1);

partName='Calcification';
febio_spec.Mesh.Elements{2}.ATTR.name=partName; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{2}.ATTR.id=2;
febio_spec.Mesh.Elements{2}.elem.ATTR.id=eleID_start+(1:1:size(E_set{2},1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E_set{2}; %The element matrix
eleID_start = eleID_start+size(E_set{2},1);

partName='Fluid';
febio_spec.Mesh.Elements{3}.ATTR.name=partName; %Name of this part
febio_spec.Mesh.Elements{3}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{3}.ATTR.id=3;
febio_spec.Mesh.Elements{3}.elem.ATTR.id=eleID_start+(1:1:size(E_set{3},1))'; %Element id's
febio_spec.Mesh.Elements{3}.elem.VAL=E_set{3}; %The element matrix
eleID_start = eleID_start+size(E_set{3},1);

partName='Fluid_BL';
febio_spec.Mesh.Elements{4}.ATTR.name=partName; %Name of this part
febio_spec.Mesh.Elements{4}.ATTR.type='penta6'; %Element type
febio_spec.Mesh.Elements{4}.ATTR.id=3; %same ID between Fluid and Fluid BL
febio_spec.Mesh.Elements{4}.elem.ATTR.id=eleID_start+(1:1:size(E_set{4},1))'; %Element id's
febio_spec.Mesh.Elements{4}.elem.VAL=E_set{4}; %The element matrix
eleID_start = eleID_start+size(E_set{4},1);

% -> Node Sets
nodeSetName1='FixDisplacement';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=ind_fix; %must be a row vector !!!

% -> Element sets (surfaces)
surfaceName1 = 'fsSurface_fluid';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.tri3.ATTR.id=(1:1:size(F_fs,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL=F_fs;

surfaceName2 = 'fluidOutSurface';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.tri3.ATTR.id=(1:1:size(F_out_fluid,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL=F_out_fluid;
fluid_out_ID_start=size(F_out_fluid,1);
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=fluid_out_ID_start+(1:1:size(F_out_BL,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=F_out_BL;

surfaceName3 = 'fluidInSurface';
febio_spec.Mesh.Surface{3}.ATTR.name=surfaceName3;
febio_spec.Mesh.Surface{3}.tri3.ATTR.id=(1:1:size(F_in_fluid,1))';
febio_spec.Mesh.Surface{3}.tri3.VAL=F_in_fluid;
fluid_in_ID_start=size(F_in_fluid,1);
febio_spec.Mesh.Surface{3}.quad4.ATTR.id=fluid_in_ID_start+(1:1:size(F_in_BL,1))';
febio_spec.Mesh.Surface{3}.quad4.VAL=F_in_BL;

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name='VesselWall';

febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.name='Calcification';

febio_spec.MeshDomains.SolidDomain{3}.ATTR.mat=materialName3;
febio_spec.MeshDomains.SolidDomain{3}.ATTR.name='Fluid';

febio_spec.MeshDomains.SolidDomain{4}.ATTR.mat=materialName3;
febio_spec.MeshDomains.SolidDomain{4}.ATTR.name='Fluid_BL';

%Boundary conditions section
% -> Fix displacement
febio_spec.Boundary.bc{1}.ATTR.name='zero_displ_xyz';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;

% -> No-slip condition
febio_spec.Boundary.bc{2}.ATTR.name='No_slip';
febio_spec.Boundary.bc{2}.ATTR.type='zero fluid velocity';
febio_spec.Boundary.bc{2}.ATTR.node_set='@surface:fsSurface_fluid';
febio_spec.Boundary.bc{2}.wx_dof=1;
febio_spec.Boundary.bc{2}.wy_dof=1;
febio_spec.Boundary.bc{2}.wz_dof=1;

% -> Outlet fluid resistance
febio_spec.Boundary.bc{3}.ATTR.name='Outlet_R';
febio_spec.Boundary.bc{3}.ATTR.type='fluid resistance';
febio_spec.Boundary.bc{3}.ATTR.surface='fluidOutSurface';
febio_spec.Boundary.bc{3}.R.ATTR.lc='1';
febio_spec.Boundary.bc{3}.R.VAL=1.334e-06; %double of the value I used without the bifurcation
febio_spec.Boundary.bc{3}.pressure_offset.ATTR.lc='1';
febio_spec.Boundary.bc{3}.pressure_offset.VAL=0.01;

%Loads section
% -> Fluid fsi tractions
febio_spec.Loads.surface_load{1}.ATTR.name='FluidFSITraction';
febio_spec.Loads.surface_load{1}.ATTR.type='fluid-FSI traction';
febio_spec.Loads.surface_load{1}.ATTR.surface='fsSurface_fluid';
febio_spec.Loads.surface_load{1}.shell_bottom=0;

% -> Inlet fluid normal velocity
%surface = fluidInSurface
febio_spec.Loads.surface_load{2}.ATTR.name='Inlet_Vel';
febio_spec.Loads.surface_load{2}.ATTR.type='fluid normal velocity';
febio_spec.Loads.surface_load{2}.ATTR.surface='fluidInSurface';
febio_spec.Loads.surface_load{2}.velocity.ATTR.lc='3';
febio_spec.Loads.surface_load{2}.velocity.VAL=-0.5; %-0.5 with a loading curve which is not normalized, otherwise -250 if max 500 
febio_spec.Loads.surface_load{2}.prescribe_nodal_velocities=1;
febio_spec.Loads.surface_load{2}.parabolic=1;
febio_spec.Loads.surface_load{2}.prescribe_rim_pressure=1;

%Load controllers section
% -> Sigmoid LC
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_sigmoid';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='CONTROL POINTS';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0,0 ; 0.124,0 ; 0.126,1 ; 0.25,1];
%febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0,0 ; 0.49,0 ; 0.51,1 ; 1,1];

% -> Must point LC
febio_spec.LoadData.load_controller{2}.ATTR.name='LC_must_points';
febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='STEP';
febio_spec.LoadData.load_controller{2}.points.pt.VAL=[0,0 ; 0.05,0.005 ; 0.1,0.005 ; 0.15,0.005 ; 0.2,0.005 ; 0.25,0.005 ; ...
                                                      0.3,0.005 ; 0.35,0.005 ; 0.4,0.005 ; 0.45,0.005 ; 0.5,0.005 ; 0.55,0.005 ; ...
                                                      0.6,0.005 ; 0.65,0.005 ; 0.7,0.005 ; 0.75,0.005 ; 0.8,0.005 ; 0.85,0.005 ; ...
                                                      0.9,0.005 ; 0.95,0.005 ; 1,0.005 ; 1.05,0.005 ; 1.10,0.005 ; 1.15,0.005 ; ...
                                                      1.20,0.005 ; 1.25,0.005 ; 1.30,0.005 ; 1.35,0.005 ; 1.40,0.005 ; 1.45,0.005 ; ...
                                                      1.50,0.005 ; 1.55,0.005];


% -> Pulsatile LC
Pulsatile_curve = importdata('3_pulsatile_max_500.txt','\t');
febio_spec.LoadData.load_controller{3}.ATTR.name='LC_pulsatile';
febio_spec.LoadData.load_controller{3}.ATTR.id=3;
febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{3}.extend='CONSTANT';
febio_spec.LoadData.load_controller{3}.points.pt.VAL=Pulsatile_curve;

%Output section
febioFebFileName=[path,'EMC028_RICA.feb']; %FEB file name
febioLogFileName=[path,'EMC028_RICA.txt']; %FEBio log file name

%Plotfile
febio_spec.Output.plotfile.ATTR.type='febio';
febio_spec.Output.plotfile.var{1}.ATTR.type='displacement';
febio_spec.Output.plotfile.var{2}.ATTR.type='stress';
febio_spec.Output.plotfile.var{3}.ATTR.type='reaction forces';
febio_spec.Output.plotfile.var{4}.ATTR.type='fluid pressure';
febio_spec.Output.plotfile.var{5}.ATTR.type='fluid velocity';
% ... add all wanted variables

%Logfile
% febio_spec.Output.logfile.node_data{1}.ATTR.data='x;y;z';
% febio_spec.Output.logfile.node_data{1}.ATTR.file='nodes_position';
% febio_spec.Output.logfile.node_data{1}.ATTR.node_set='@surface:fsSurface'; %node set chosen just as example to try this option
% % 
% febio_spec.Output.logfile.node_data{2}.ATTR.data='ux;uy;uz';
% febio_spec.Output.logfile.node_data{2}.ATTR.file='nodes_velocity';
% febio_spec.Output.logfile.node_data{2}.ATTR.node_set='@surface:fsSurface';
% ... add all the data to be saved in the logfile - or in a separate file (to be specified, as done here)

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
save('EMC028_RICA.mat', '-struct', 'febio_spec');

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; % = 1, Display information on the command window
febioAnalysis.runMode='external';
[runFlag]=runMonitorFEBio(febioAnalysis); %START FEBio NOW!!!!!!!!
