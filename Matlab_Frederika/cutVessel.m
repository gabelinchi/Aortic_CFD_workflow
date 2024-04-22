function [F_new,V_new] = cutVessel(F_fs,V_fs,F_inlet,V_inlet,F_outlet,V_outlet) 

[F_new,V_new] = Subtraction(F_fs,V_fs,F_inlet,V_inlet,1);
[F_new,V_new] = Subtraction(F_new,V_new,F_outlet,V_outlet,2);

optionStruct.pre.min_comp_area = 0;
optionStruct.post.min_comp_area = 0;
optionStruct.pointSpacing=0.4;
[F_new,V_new]=ggremesh(F_new,V_new,optionStruct); %remesh after the cutting = the subtraction

end

%% function Subtraction

function [F_out,V_out] = Subtraction(F_in,V_in,F_ref,V_ref,type) %type = 1 or 2, depending on that the bias for the cylinder is set tu -bias (inverted) - this quantity is used for translation of the cylinder used for boolean operation (in function subtraction)

[~,V_ref] = patchCleanUnused(F_ref,V_ref);

rcoor = sqrt( (V_ref(:,1) - mean(V_ref(:,1)) ).^2 + ...
          (V_ref(:,2) - mean(V_ref(:,2)) ).^2 + ...
          (V_ref(:,3) - mean(V_ref(:,3)) ).^2 ); %raggi per ogni nodo - con mean calcoliamo le coord x,y,z del centroide della sezione
%il centroide coincide con il centro del cerchio solo se i nodi sono
%equamente distribuiti (e in questo caso ne siamo sicuri?) pero' piu o meno
%sono distribuiti equamente, quindi eventuale errore sara' piccolo

r_ref = max(rcoor); %identifichiamo il raggio piu grande
bias = 2;
r_cylinder = 1.5*r_ref; %il cilindro deve avere raggio 1.5 volte piu grande del massimo raggio della sezione

if type == 2
    bias = -bias;
end
% if min(V_in(:,3))<0
%     bias = -bias;
% end

[F_box,V_box] = CreateCylinder(r_cylinder,V_ref,type); %created cylinder for subtraction -> another function

T  = [1 0 0 0;...
      0 1 0 0;...
      0 0 1 bias;...
      0 0 0 1];
V_box=tform(T,V_box); %cilindro traslato secondo la matrice T, cosi da intersecare il vaso

figure('WindowState', 'maximized');
gpatch(F_box,V_box,'r',0.3); hold on
gpatch(F_in,V_in,'kw','none',0.3)
view(45,25)
axis equal
%close(gcf);

[F_new,V_new] = booleanOperation(F_in,V_in,F_box,V_box,1); %function from another toolbox, that requires also other functions from the same toolbox

[G] = tesgroup(F_new); %assign index 1 or 2 to the faces F_new depending on which side they belong to (after the cutting, we have 2 parts/sides)

%%Use this if the cutting was not done properly ?? - ex. cylinder used was too small - r_cylinder is 1.5 r_ref right now 
% while size(G,2) == 1
%     [F_box,V_box] = CreateCylinder(1.1*r_cylinder,V_ref,bias);
%     [F_new,V_new] = booleanOperation(F_in,V_in,F_box,V_box,1);
%     [G] = tesgroup(F_new);
% end

%Process to get rid of the extremes generated after the boolean operation with the cylinder
 num = sum(G);
 [~,ind_save] = max(num);
 F_new = F_new(G(:,ind_save),:);

%Process to get rid of the inlet/outlet surfaces 
a=(30/180)*pi;
G=patchFeatureDetect(F_new,V_new,a); %Detect surface features
clear num
for i = 1:max(G)
    num(i) = sum(ismember(G,i));
end
num_sort = sort(num,'descend');
F_test = F_new(G==find(num==num_sort(2)),:);
figure;
gpatch(F_new,V_new,'kw','none',0.3);
gpatch(F_test,V_new);
view(45,25)
axis equal
F_out = F_new(G~=find(num==num_sort(2)),:);
V_out = V_new;
[~,~] = patchCleanUnused(F_out,V_out);
close(gcf)
 
 figure('WindowState', 'maximized');
 gpatch(F_out,V_out);
 gpatch(F_in,V_in,'kw','none',0.3)
 title('Result of cutting operation', num2str(type), 'fontSize',25);
 view(45,25)
 axis equal
 %close(gcf)

end

%% function CreateCylinder

function [F_box,V_box] = CreateCylinder(r_cylinder,V_ref,bias)

pointSpacing=0.3; %Desired point spacing
resampleCurveOpt=0;
interpMethod='linear'; %or 'natural'
[F_ref,V_ref]=regionTriMesh3D({V_ref},pointSpacing,resampleCurveOpt,interpMethod);

[~,~,N]=patchNormal(F_ref,V_ref);
N = mean(N);
% Normalize the normal vector
N = N / norm(N);
V_mean = mean(V_ref);
pointSpacing=0.5;
cylRadius=r_cylinder;
cylLength=0.3;
[meshStruct]=hexMeshCylinder(cylRadius,cylLength,pointSpacing);
F_box = meshStruct.facesBoundary;
V_box = meshStruct.nodes;
[F_box,V_box] = patchCleanUnused(F_box,V_box);
[F_box,V_box]=quad2tri(F_box,V_box,'f');
OR = V_mean / norm(V_mean)*(norm(V_mean)+bias);
T  = [1 0 0 OR(1);...
      0 1 0 OR(2);...
      0 0 1 OR(3);...
      0 0 0 1];

rotation_axis = cross([0 0 1], N);
rotation_angle = acos(dot([0 0 1], N));

% Create the rotation matrix using the axis-angle representation
u = rotation_axis / norm(rotation_axis); % Normalize the rotation axis
cos_theta = cos(rotation_angle);
sin_theta = sin(rotation_angle);
one_minus_cos_theta = 1 - cos_theta;

R = [cos_theta + u(1)^2 * one_minus_cos_theta,      u(1) * u(2) * one_minus_cos_theta - u(3) * sin_theta, u(1) * u(3) * one_minus_cos_theta + u(2) * sin_theta;
     u(2) * u(1) * one_minus_cos_theta + u(3) * sin_theta, cos_theta + u(2)^2 * one_minus_cos_theta,      u(2) * u(3) * one_minus_cos_theta - u(1) * sin_theta;
     u(3) * u(1) * one_minus_cos_theta - u(2) * sin_theta, u(3) * u(2) * one_minus_cos_theta + u(1) * sin_theta, cos_theta + u(3)^2 * one_minus_cos_theta];

[V_box]=rotate_vertices(V_box,R);
V_box=tform(T,V_box);

end