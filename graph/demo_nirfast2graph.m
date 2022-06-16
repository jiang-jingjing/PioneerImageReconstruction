%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert nirfast mesh to graph
% Author: Jingjing Jiang jing.jing.jiang@outlook.com
%  
% created on 2022.03.25 
% modified  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add paths for required packages
% define nirfast paths
% Nirfaster: GPU-fascilitated model
% Nirfast8: old CPU version
pathNirfaster = '/media/jiang/WD10T/Data//SoftwarePackages/NIRFASTer';
pathNirfast8 = '/media/jiang/WD10T/Data/SoftwarePackages/nirfast8';
addpath(genpath(pathNirfast8))
pathPioneerIR = '/media/jiang/WD10T/Data/Projects/PioneerImageReconstruction';
addpath(genpath(pathPioneerIR))
%% load mesh  / create meshfldr_data = '/media/jiang/WD10T/Data/Projects/PioneerImageReconstruction/demo_nirfast_2d';

fldr_data = '/media/jiang/WD10T/Data/Projects/PioneerImageReconstruction/demo_nirfast_2d';
mesh = load_mesh([fldr_data '/mesh/Rectangle-stnd-mesh']);
%% conver to Graph
NN = size(mesh.nodes,1);
MatA = zeros(NN,NN);
[NE1 NE2] = size(mesh.elements);
for ie = 1:NE1
    if NE2 == 3
    v1 = mesh.elements(ie,1);
    v2 = mesh.elements(ie,2);
    v3 = mesh.elements(ie,3);
    MatA(v1,v2) = 1;
    MatA(v2,v1) = 1;
    MatA(v2,v3) = 1;
    MatA(v3,v2) = 1;
    MatA(v1,v3) = 1;
    MatA(v3,v1) = 1;
    end
end
figure, imagesc(MatA)
% Mat_Adj = mesh.elements;
G_mesh = graph(MatA);

%% visualization
MuaLabel = num2str(mesh.mua);
MuaLabel = string(MuaLabel(:,4));
% Plot graph.
figure
% plot(G_mesh,'NodeLabel',MuaLabel)
plot(G_mesh)
 