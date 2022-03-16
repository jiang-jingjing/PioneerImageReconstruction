function img2d=mesh2image2D(mesh,dens ,np)
%>@ function img2d=mesh2image2D(mesh,dens ,np)
%>
%> convert a 2D nirfast mesh to an image
%>  
%> author: jingjing jiang jing.jing.jiang@outlook.com
%> created: 2022.03.15
% adapted from % function [img2d_coord, conc]=CrossSecMesh2(mesh,img2d_type,coord,dens,F,np)
% 
%> mesh is a NIRFAST mesh
%> dens is a size of pixels in mm, i.e. the mesh size of 20 mm will result 
%>      in 20 pixels of desn=1
%> np is a number of points in each dimention, it is 2 values vector,
%>       each value reperesents its own axis, e.g. for 'XY' type the
%>       np=[100, 50] means 100 pints along X axis and 50 points 
%>      along Y axis

if nargin==3
    flg=1;
else 
    flg=0;
end

%% Create cross-section
% define limits of the field of view
% xlim=[min(mesh.nodes(:,1))-3, max(mesh.nodes(:,1))+3];
% ylim=[min(mesh.nodes(:,2))-3, max(mesh.nodes(:,2))+3];
xlim=[min(mesh.nodes(:,1)), max(mesh.nodes(:,1))];
ylim=[min(mesh.nodes(:,2)), max(mesh.nodes(:,2))];
%  get mua and mus
F{1}=scatteredInterpolant(mesh.nodes(:,1:2),mesh.mua(:),...
        'linear','none');
F{2}=scatteredInterpolant(mesh.nodes(:,1:2),mesh.mus(:),...
        'linear','none');
% define step of discretization
if flg==1
    xstep=(xlim(2)-xlim(1))/np(1);
    ystep=(ylim(2)-ylim(1))/np(2);
elseif flg==0
    xstep=dens;
    ystep=dens;
    np(1)=ceil((xlim(2)-xlim(1))/xstep);
    np(2)=ceil((ylim(2)-ylim(1))/ystep);
end
% create coordinates of the cross-section points
tx=xlim(1)-0.5*xstep:xstep:xlim(2)-0.5*xstep;
ty=ylim(1)-0.5*ystep:ystep:ylim(2)-0.5*ystep;
[xq,yq]=meshgrid(tx,ty);

%% Calculate values

if strcmp(mesh.type,'spec')
    % get concentration values
    for j=1:size(mesh.conc,2)
        conc(:,:,j)=F{j}(xq,yq);
    end
elseif strcmp(mesh.type,'stnd')
    % get absorption and scattering values
    conc(:,:,1)=F{1}(xq,yq);
    conc(:,:,2)=F{2}(xq,yq);
else
    error('Unknown mesh type');
end

% make coordinates output
img2d_coord(:,:,1)=xq;
img2d_coord(:,:,2)=yq;
img2d.coord = img2d_coord;
img2d.mua = squeeze(conc(:,:,1));
img2d.mus = squeeze(conc(:,:,2));
end