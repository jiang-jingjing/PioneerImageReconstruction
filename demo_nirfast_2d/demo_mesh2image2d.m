% demo of using mesh2image2D.m
% conver mesh to 2D image
% It requires the preloaded mesh_anom / or running demo_2d_nirfast.m first
%% plot original mesh
plotimage(mesh_anom,mesh_anom.mua)
title('original mesh')
%% convert 2D mesh to 2D image
img=mesh2image2D(mesh_anom,0.1)
%% plot converted 2D image
figure
imagesc(img.mua)
title('converted image')
axis xy
axis equal; 
axis off;
colormap hot;    
colorbar('horiz');
axis equal; 
axis off;
colormap hot;