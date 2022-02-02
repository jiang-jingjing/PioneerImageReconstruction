function [cs_coord, conc]=CrossSecMesh2(mesh,cs_type,coord,dens,F,np)
%CrossSecMesh   Calculates values of NIRFAST mesh along cross section plane
% mesh is a NIRFAST mesh
% cs_type is a type of cross-section: 'XY', 'XZ', 'YZ'
% coord is a coordinate of the cross-section, e.g. for 'XY' type the coord
%       is the coordinate along Z axis
% dens is a size of pixels in mm, i.e. the mesh size of 20 mm will result 
%       in 20 pixels of desn=1
% np is a number of points in each dimention, it is 2 values vector,
%       each value reperesents its own axis, e.g. for 'XY' type the
%       np=[100, 50] means 100 pints along X axis and 50 points 
%       along Y axis

if nargin==6
    flg=1;
else 
    flg=0;
end

%% Create cross-section
% define limits of the field of view
xlim=[min(mesh.nodes(:,1))-3, max(mesh.nodes(:,1))+3];
ylim=[min(mesh.nodes(:,2))-3, max(mesh.nodes(:,2))+3];
zlim=[min(mesh.nodes(:,3))-1, max(mesh.nodes(:,3))+1];

% create the section points
if isequal(cs_type,'XY')
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
    zq=coord.*ones(size(xq));
elseif isequal(cs_type,'XZ')
    % define step of discretization
    if flg==1
        xstep=(xlim(2)-xlim(1))/np(1);
        zstep=(zlim(2)-zlim(1))/np(2);
    elseif flg==0
        xstep=dens;
        zstep=dens;
        np(1)=ceil((xlim(2)-xlim(1))/xstep);
        np(2)=ceil((zlim(2)-zlim(1))/zstep);
    end
    % create coordinates of the cross-section points
    tx=xlim(1)-0.5*xstep:xstep:xlim(2)-0.5*xstep;
    tz=zlim(1)-0.5*zstep:zstep:zlim(2)-0.5*zstep;
    [xq,zq]=meshgrid(tx,tz);
    yq=coord.*ones(size(xq));
elseif isequal(cs_type,'YZ')
    % define step of discretization
    if flg==1
        ystep=(ylim(2)-ylim(1))/np(1);
        zstep=(zlim(2)-zlim(1))/np(2);
    elseif flg==0
        ystep=dens;
        zstep=dens;
        np(1)=ceil((ylim(2)-ylim(1))/zstep);
        np(2)=ceil((zlim(2)-zlim(1))/zstep);
    end
    % create coordinates of the cross-section points
    ty=ylim(1)-0.5*ystep:ystep:ylim(2)-0.5*ystep;
    tz=zlim(1)-0.5*zstep:zstep:zlim(2)-0.5*zstep;
    [yq,zq]=meshgrid(ty,tz);
    xq=coord.*ones(size(yq));
else
    error('Undefined type of cross-section');
end

%% Calculate values

if strcmp(mesh.type,'spec')
    % get concentration values
    for j=1:size(mesh.conc,2)
        conc(:,:,j)=F{j}(xq,yq,zq);
    end
elseif strcmp(mesh.type,'stnd')
    % get absorption and scattering values
    conc(:,:,1)=F{1}(xq,yq,zq);
    conc(:,:,2)=F{2}(xq,yq,zq);
else
    error('Unknown mesh type');
end

% make coordinates output
cs_coord(:,:,1)=xq;
cs_coord(:,:,2)=yq;
cs_coord(:,:,3)=zq;

end