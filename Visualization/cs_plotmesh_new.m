function CS=cs_plotmesh_new(mesh,coord,dens,flg,JS)
%CS_PLOTMESH    Plots cross sections of 3D NIRFAST mesh
%   The function is designed to supplement the standard NIRFAST plot tools.
%   It plots cross sections of the mesh with an access to MATLAB figure
%   features, such as colormap and zoom.
%   CS_PLOTMESH is not designed to replace the plotmesh() function from
%   NIRFAST since it has neither sliding features no 3D view, and it is 
%   significantly slower.
%
%   CS = CS_PLOTMESH(MESH,COORD,DENS) plots pixelated cross sections of 
%   MESH at position COORD with pixel size DENS in mm. By default tissue 
%   oxygen saturation StO2 is plotted for spectral mesh and absorption mua 
%   for standard mesh, the plots are grouped in one window. CS is an output
%   with all the values of the 3 cross sections in a structure format 
%   CS.XX.aaaaa, where XX is a cross section type (XY, XZ, YZ), and aaaaa 
%   is the data type ('coord', 'conc', 'StO2', etc.)
%
%   CS = CS_PLOTMESH(MESH,COORD,DENS,FLG,JS) plots cross sections in a mode 
%   set by FLG in one window (JS='jnt') or in separate windows (JS='sep').
%   By default FLG is equal to 1 and JS is equal to 'jnt'. 
%       FLG = 0 -- no plot, just CS output is generated
%       FLG = 1 -- default view, StO2 or mua plots
%       FLG = 2,3,.. -- generate plots of chromophores and optical 
%           properties. Value FLG-1 defines the number of a chromophore in
%           the order it is stored in the MESH. In case of standard MESH
%           the FLG-1 value refers to either mua (FLG=2) or mus (FLG=3)
%       NOTE: FLG=1 and FLG=2 for 'stnd' mesh results in the same output
%
%       Negative FLG value results in the same output as abs(FLG) but the
%       XZ cross section is flipped around Z axes. In this case the cross
%       section appears on a screen in the same way as the standard
%       plotmesh() output. Since the right-handed coordinate system is used
%       in NIRFAST and this representation is kept through all the views,
%       the directions of the axis are as following:
%           view 1 (top left corner) -- 3D view, X: lr (left to right), 
%               Y: bt (bottom to top), Z: out (out of the screen)
%           view 2 (top right) -- YZ section, X: out, Y: lr, Z: bt
%           view 3 (bottom left) -- XZ section, X: rl, Y: out, Z: bt
%           view 4 (bottom right) -- XY section, X: lr, Y: bt, Z: out
%       With negative FLG values the cross section views are set to mimic
%       the standard representation, even though there is no meaning of 3D
%       rotation in cross section view.
%
% subfunctions  CrossSecMesh.m
%
% Author: Alexander Kalyanov
% e-mail: Alexander.Kalyanov@usz.ch
% Release: 0.9.1
% Release date: 28/04/2017
%
% modified: 22/11/2019 by Jingjing Jiang jjiang@ethz.ch 
% 0.9.1 release notes:
% Interpolation algorithm is employed instead of search of closest nodes.
% By this the speed of the function was significantly increased, it is now
% up to 6 times faster.

%% Check the input and set the mode
if ~strcmp(mesh.type,'spec') && ~strcmp(mesh.type,'stnd')
    % only specteral and standard 3D meshes are supported 
    error('Wrong mesh type, only standard and spectral 3D meshes are supported');
end

flg_NFV=0;              % standard non NIRFAST view
if nargin==3            % check if the flg variable exist
    flg=1;              % if not, set flg to 1: standard plotting
    JS='jnt';           % joint/separate set to joint
elseif flg<0            % negative flg, NIRFAST view mode
    flg=abs(flg);       % set flg to positive
    flg_NFV=1;          % set flg_NFV to 1, NIRFAST veiw
end
if flg>1                % if flg exist and referes to nonstandard mode
    if strcmp(mesh.type,'spec')
        % if mesh is spectral, check the number of chromophores,
        % the requiered chromophor (flg-1) should exist in the list
        if (flg-1)>(size(mesh.chromscattlist,1)-2)
            error('The chromophore does not exist');
        end
    elseif strcmp(mesh.type,'stnd')
        % if mesh is standard, check that the requiered optical parameter
        % (flg-1) is less than 2 (there are only 2 optical parameters to be
        % plotted: mua and mus
        if (flg-1)>2
            error('Reference to the wrong optical parameter');
        end
    end
end

%% Calculate cross sections
if strcmp(mesh.type,'spec')     % calculate chromophore concentrations
    for i=1:size(mesh.conc,2)
        F{i}=scatteredInterpolant(mesh.nodes,mesh.conc(:,i),...
            'natural','none');
    end
    % XY
    [CS.XY.coord, CS.XY.conc]=CrossSecMesh2(mesh,'XY',coord(3),dens,F);
    % XZ
    [CS.XZ.coord, CS.XZ.conc]=CrossSecMesh2(mesh,'XZ',coord(2),dens,F);
    % YZ
    [CS.YZ.coord, CS.YZ.conc]=CrossSecMesh2(mesh,'YZ',coord(1),dens,F);
    % Calculate StO2 values
    CS.XY.StO2=CS.XY.conc(:,:,strcmp(mesh.chromscattlist,'HbO'))./ ...
        (CS.XY.conc(:,:,strcmp(mesh.chromscattlist,'HbO'))+...
        CS.XY.conc(:,:,strcmp(mesh.chromscattlist,'deoxyHb')));
    CS.XZ.StO2=CS.XZ.conc(:,:,strcmp(mesh.chromscattlist,'HbO'))./ ...
        (CS.XZ.conc(:,:,strcmp(mesh.chromscattlist,'HbO'))+...
        CS.XZ.conc(:,:,strcmp(mesh.chromscattlist,'deoxyHb')));
    CS.YZ.StO2=CS.YZ.conc(:,:,strcmp(mesh.chromscattlist,'HbO'))./ ...
        (CS.YZ.conc(:,:,strcmp(mesh.chromscattlist,'HbO'))+...
        CS.YZ.conc(:,:,strcmp(mesh.chromscattlist,'deoxyHb')));
elseif strcmp(mesh.type,'stnd') % calcualte optical properties
%     F{1}=scatteredInterpolant(mesh.nodes,mesh.mua(:),...
%         'natural','none');
%     F{2}=scatteredInterpolant(mesh.nodes,mesh.mus(:),...
%         'natural','none'); % COMMENTED BY JINGJING 2019.06.13
% ----------------------- ADDED BY JINGJING 2019.06.13
    F{1}=scatteredInterpolant(mesh.nodes,mesh.mua(:),...
        'linear','none');
    F{2}=scatteredInterpolant(mesh.nodes,mesh.mus(:),...
        'linear','none');
% ------------------------ END
    % XY
    [CS.XY.coord, CS.XY.op]=CrossSecMesh2(mesh,'XY',coord(3),dens,F);
    % XZ
    [CS.XZ.coord, CS.XZ.op]=CrossSecMesh2(mesh,'XZ',coord(2),dens,F);
    % YZ
    [CS.YZ.coord, CS.YZ.op]=CrossSecMesh2(mesh,'YZ',coord(1),dens,F);
end

%% Plot data
if flg==1       % standard plot
    if strcmp(mesh.type,'spec')     % plot StO2
        
        dscr={'StO_2'; ...
            ['Position: (' num2str(coord(1)) ', ' num2str(coord(2)) ...
            ', ' num2str(coord(3)) ')']};
        PlotCrsSec(CS,CS.XY.StO2,CS.XZ.StO2,CS.YZ.StO2, ...
            dscr,flg_NFV,[0 1],JS,1)
    elseif strcmp(mesh.type,'stnd') % plot mua (absorption)
        dscr={'Absorption mua'; ...
            ['Position: (' num2str(coord(1)) ', ' num2str(coord(2)) ...
            ', ' num2str(coord(3)) ')']};
        PlotCrsSec(CS,CS.XY.op(:,:,1),CS.XZ.op(:,:,1),CS.YZ.op(:,:,1),...
            dscr,flg_NFV,[],JS)
    end
elseif flg>1
    u=flg-1; % plot mode referes to cromophore index, or to mua/mus values
    if strcmp(mesh.type,'spec')
        dscr{1}=mesh.chromscattlist{u};
        dscr{2}=['Position: (' num2str(coord(1)) ', ' num2str(coord(2)) ...
            ', ' num2str(coord(3)) ')'];
        PlotCrsSec(CS,CS.XY.conc(:,:,u),CS.XZ.conc(:,:,u), ...
            CS.YZ.conc(:,:,u),dscr,flg_NFV,[],JS)
    elseif strcmp(mesh.type,'stnd')
        foo={'Mua', 'Mus'};
        dscr{1}=foo{u};
        dscr{2}=['Position: (' num2str(coord(1)) ', ' num2str(coord(2)) ...
            ', ' num2str(coord(3)) ')'];
        PlotCrsSec(CS,CS.XY.op(:,:,u),CS.XZ.op(:,:,u), ...
            CS.YZ.op(:,:,u),dscr,flg_NFV,[],JS)
    else
        error('Unknown mesh type');
    end
end

end

function PlotCrsSec(CS,XY,XZ,YZ,dscr,flg_NFV,cAxLm,JS,flg)

if nargin==8
    flg=0;
end

if strcmp(JS,'jnt')                         % plot the sections together
   hf= figure;                                 % create a new image
end

% X-Z
V=GetFoV(CS.XZ.coord);                      % get FoV
if strcmp(JS,'jnt')
%     subplot(2,2,3);
dx = 0.45;
dz = 0.25;
% pos1 = [0.55 0.45 dz dx];
pos1 = [0.1 0.1 dx dz];
% h_xz=subplot('Position',pos1);
h_xz=axes('Position',pos1);

%     imagesc(XZ','XData',V(3,:),'YData',V(1,:));% plot XZ section in subplot
    imagesc(XZ,'XData',V(1,:),'YData',V(3,:));% plot XZ section in subplot

elseif strcmp(JS,'sep')
    rect=get(0,'defaultfigureposition');
    w=round(0.5*rect(3))+5;
    h=rect(4)+80;
    figure('Name','X-Z cross section',...
        'Position',[rect(1)-w rect(2)-h rect(3) rect(4)]);
    imagesc(XZ,'XData',V(1,:),'YData',V(3,:));% plot XZ section in new fig
end
% imagesc(XZ,'XData',V(1,:),'YData',V(3,:));  % plot XZ section
title('X-Z cross section','FontSize',10);
hcb = colorbar; hcb.Visible = 'off';
axis xy; axis image
ax = gca;

% ax.YTick = 0:10:60;
[ZM XM] = size(XZ)  ;
ax.XTick = 0:10:(XM-6);
% ax.YTick = 0:10:(ZM-2);
% ax.YTick = 0:10:(ZM-6);
% ax.XTick = 0:10:60;
set(hcb, 'visible','off')
% xlabel('X, mm'); ylabel('Z, mm');
% xlabel('Z, mm','FontSize',9); 
xlabel('X, mm','FontSize',9); 
ylabel('Z, mm','FontSize',9); 
if ~isempty(cAxLm)                          % set the caxis if provided
    caxis(cAxLm);
else
    ma_xz=max(XZ(:));
    mi_xz=min(XZ(:));
    % adjust the caxis if all the values are the same
    if ma_xz==mi_xz
        ma_xz=ma_xz+0.01*abs(ma_xz);
        mi_xz=mi_xz-0.01*abs(mi_xz);
%         caxis([mi_xz ma_xz]);
    end
end
if flg_NFV==1                               % flip X-axis if needed
    set(gca,'xdir','reverse');
end
if flg==1
    colormap(cm_sto2);
end

% Y-Z
V=GetFoV(CS.YZ.coord);                      % get FoV
if strcmp(JS,'jnt')
%     subplot(2,2,2);
% pos1 = [0.1 0.1 dx dz];
pos1 = [0.55 0.45 dz dx];
% h_yz = subplot('Position',pos1);
h_yz = axes('Position',pos1);
% imagesc(YZ,'XData',V(2,:),'YData',V(3,:));% plot YZ section in subplot
imagesc(YZ','XData',V(3,:),'YData',V(2,:));% plot YZ section in subplot

elseif strcmp(JS,'sep')
    figure('Name','Y-Z cross section', ...
        'Position',[rect(1)+w rect(2) rect(3) rect(4)]);
    imagesc(YZ,'XData',V(2,:),'YData',V(3,:));% plot YZ section in new fig
end
title('Y-Z cross section','FontSize',10);
% colorbar;
axis xy; axis image
[YM ZM] = size(YZ) ;
ax.YTick = 0:10:(YM-6);
ax.XTick = 0:10:(ZM-6);
% xlabel('Y, mm','FontSize',9); ylabel('Z, mm','FontSize',9);
xlabel('Z, mm','FontSize',9); 

if ~isempty(cAxLm)                          % set the caxis if provided
    caxis(cAxLm);
else
    ma_yz=max(YZ(:));
    mi_yz=min(YZ(:));
    % adjust the caxis if all the values are the same
    if ma_yz==mi_yz
        ma_yz=ma_yz+0.01*abs(ma_yz);
        mi_yz=mi_yz-0.01*abs(mi_yz);
%         caxis([mi_yz ma_yz]);
    end
end
if flg==1
    colormap(cm_sto2);
end

% X-Y
V=GetFoV(CS.XY.coord);                      % get FoV
if strcmp(JS,'jnt')
%     subplot(2,2,4);
%     imagesc(XY,'XData',V(1,:),'YData',V(2,:));% plot XY section
    pos1 = [0.1 0.45 dx  dx ];
%     h_xy = subplot('Position',pos1);
    h_xy = axes('Position',pos1);

    imagesc(XY,'XData',V(1,:),'YData',V(2,:));% plot XY section
%     imagesc(XY','XData',V(2,:),'YData',V(1,:));% plot XY section
elseif strcmp(JS,'sep')
    figure('Name','X-Y cross section',...
        'Position',[rect(1)+w rect(2)-h rect(3) rect(4)]);
    imagesc(XY,'XData',V(1,:),'YData',V(2,:));% plot XY section in new fig
end    
title('X-Y cross section','FontSize',10);
% colorbar; 
axis xy; axis image
ax = gca;
[XM YM] = size(XY) ;
ax.YTick = 0:10:(YM-6);
ax.XTick = 0:10:(XM-6);
ax.FontSize = 8;
% xlabel('Y, mm'); 
ylabel('Y, mm', 'FontSize',9);
if ~isempty(cAxLm)                          % set the caxis if provided
    caxis(cAxLm);
else
    ma_xy=max(XY(:));
    mi_xy=min(XY(:));
    % adjust the caxis if all the values are the same
    if ma_xy==mi_xy
        ma_xy=ma_xy+0.01*abs(ma_xy);
        mi_xy=mi_xy-0.01*abs(mi_xy);
%         caxis([mi_xy ma_xy]);
    end
end

% uniq color range
mi = min([mi_xz mi_yz mi_xy]);
ma = max([ma_xz ma_yz ma_xy]);
axes(h_yz)
caxis([mi  ma ]);
axes(h_xy)
caxis([mi  ma ]);
axes(h_xz)
caxis([mi  ma ]);
if flg==1
    colormap(cm_sto2);
end

hcb.Position = [0.8  0.1 0.03 0.8];
hcb.Visible = 'on';

% add discription
if strcmp(JS,'jnt')
    ax = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax);
%     text(0.15, 0.85, dscr, 'FontSize', 12)
    text(0.55, 0.25, dscr, 'FontSize', 8)

elseif strcmp(JS,'sep')
    f=round(0.4*rect(3));
    rect1=[rect(1)-w rect(2)+rect(4)-f f f];
    figure('Name','Details', 'Position',rect1);
    ax = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax);
    text(0.15, 0.85, dscr, 'FontSize', 12);
end

end

function V=GetFoV(coord)
% minimum and maximum along X axis
V(1,1)=min(min(coord(:,:,1)));
V(1,2)=max(max(coord(:,:,1)));
% along Y axis
V(2,1)=min(min(coord(:,:,2)));
V(2,2)=max(max(coord(:,:,2)));
% along Z axis
V(3,1)=min(min(coord(:,:,3)));
V(3,2)=max(max(coord(:,:,3)));

end

function cmap=cm_sto2

cmap=[
         0         0         0
    0.0392    0.0392    0.5882
    0.0426    0.0426    0.6061
    0.0460    0.0460    0.6240
    0.0494    0.0494    0.6419
    0.0529    0.0529    0.6598
    0.0563    0.0563    0.6777
    0.0597    0.0597    0.6957
    0.0631    0.0631    0.7136
    0.0665    0.0665    0.7315
    0.0699    0.0699    0.7494
    0.0733    0.0733    0.7673
    0.0767    0.0767    0.7852
    0.0801    0.0801    0.8031
    0.0835    0.0835    0.8210
    0.0870    0.0870    0.8389
    0.0904    0.0904    0.8568
    0.0938    0.0938    0.8747
    0.0972    0.0972    0.8926
    0.1006    0.1006    0.9105
    0.1040    0.1040    0.9284
    0.1074    0.1074    0.9463
    0.1108    0.1108    0.9642
    0.1142    0.1142    0.9821
    0.1176    0.1176    1.0000
    0.1147    0.1336    0.9912
    0.1117    0.1496    0.9824
    0.1087    0.1656    0.9736
    0.1057    0.1816    0.9649
    0.1027    0.1976    0.9561
    0.0998    0.2136    0.9473
    0.0968    0.2296    0.9385
    0.0938    0.2456    0.9297
    0.0908    0.2616    0.9209
    0.0878    0.2776    0.9122
    0.0849    0.2936    0.9034
    0.0819    0.3096    0.8946
    0.0789    0.3256    0.8858
    0.0759    0.3416    0.8770
    0.0729    0.3576    0.8682
    0.0700    0.3736    0.8595
    0.0670    0.3896    0.8507
    0.0640    0.4056    0.8419
    0.0610    0.4216    0.8331
    0.0580    0.4376    0.8243
    0.0551    0.4536    0.8155
    0.0521    0.4696    0.8067
    0.0491    0.4856    0.7980
    0.0461    0.5016    0.7892
    0.0431    0.5176    0.7804
    0.0790    0.5312    0.7806
    0.1148    0.5447    0.7807
    0.1506    0.5583    0.7809
    0.1865    0.5718    0.7811
    0.2223    0.5854    0.7813
    0.2581    0.5989    0.7815
    0.2939    0.6125    0.7816
    0.3298    0.6260    0.7818
    0.3656    0.6396    0.7820
    0.4014    0.6531    0.7822
    0.4373    0.6667    0.7823
    0.4731    0.6802    0.7825
    0.5089    0.6938    0.7827
    0.5448    0.7073    0.7829
    0.5806    0.7209    0.7831
    0.6164    0.7344    0.7832
    0.6522    0.7480    0.7834
    0.6881    0.7615    0.7836
    0.7239    0.7751    0.7838
    0.7597    0.7886    0.7839
    0.7956    0.8022    0.7841
    0.8314    0.8157    0.7843
    0.8314    0.8157    0.7843
    0.8314    0.8157    0.7843
    0.8314    0.8157    0.7843
    0.8314    0.8157    0.7843
    0.8314    0.8157    0.7843
    0.8314    0.8157    0.7843
    0.8428    0.7873    0.7578
    0.8543    0.7588    0.7314
    0.8657    0.7304    0.7049
    0.8771    0.7020    0.6784
    0.8886    0.6735    0.6520
    0.9000    0.6451    0.6255
    0.9114    0.6167    0.5990
    0.9229    0.5882    0.5725
    0.9343    0.5598    0.5461
    0.9458    0.5314    0.5196
    0.9572    0.5029    0.4931
    0.9686    0.4745    0.4667
    0.9718    0.4769    0.5200
    0.9749    0.4792    0.5733
    0.9780    0.4816    0.6267
    0.9812    0.4839    0.6800
    0.9843    0.4863    0.7333
    0.9875    0.4886    0.7867
    0.9906    0.4910    0.8400
    0.9937    0.4933    0.8933
    0.9969    0.4957    0.9467
    1.0000    0.4980    1.0000];
    
end

%%
% function cmap=cm_sto2
% 
% cmap=[
%          0         0         0
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%          0         0    1.0000
%     0.0177    0.0174    0.9954
%     0.0354    0.0347    0.9908
%     0.0531    0.0521    0.9862
%     0.0708    0.0694    0.9816
%     0.0884    0.0868    0.9771
%     0.1061    0.1041    0.9725
%     0.1238    0.1215    0.9679
%     0.1415    0.1388    0.9633
%     0.1592    0.1562    0.9587
%     0.1769    0.1736    0.9541
%     0.1946    0.1909    0.9495
%     0.2123    0.2083    0.9449
%     0.2300    0.2256    0.9403
%     0.2477    0.2430    0.9357
%     0.2653    0.2603    0.9312
%     0.2830    0.2777    0.9266
%     0.3007    0.2950    0.9220
%     0.3184    0.3124    0.9174
%     0.3361    0.3298    0.9128
%     0.3538    0.3471    0.9082
%     0.3715    0.3645    0.9036
%     0.3892    0.3818    0.8990
%     0.4069    0.3992    0.8944
%     0.4245    0.4165    0.8899
%     0.4422    0.4339    0.8853
%     0.4599    0.4512    0.8807
%     0.4776    0.4686    0.8761
%     0.4953    0.4859    0.8715
%     0.5130    0.5033    0.8669
%     0.5307    0.5207    0.8623
%     0.5484    0.5380    0.8577
%     0.5661    0.5554    0.8531
%     0.5837    0.5727    0.8486
%     0.6014    0.5901    0.8440
%     0.6191    0.6074    0.8394
%     0.6368    0.6248    0.8348
%     0.6545    0.6421    0.8302
%     0.6722    0.6595    0.8256
%     0.6899    0.6769    0.8210
%     0.7076    0.6942    0.8164
%     0.7253    0.7116    0.8118
%     0.7430    0.7289    0.8072
%     0.7606    0.7463    0.8027
%     0.7783    0.7636    0.7981
%     0.7960    0.7810    0.7935
%     0.8137    0.7983    0.7889
%     0.8314    0.8157    0.7843
%     0.8314    0.8157    0.7843
%     0.8314    0.8157    0.7843
%     0.8314    0.8157    0.7843
%     0.8314    0.8157    0.7843
%     0.8314    0.8157    0.7843
%     0.8314    0.8157    0.7843
%     0.8390    0.8012    0.7713
%     0.8467    0.7868    0.7583
%     0.8544    0.7724    0.7453
%     0.8620    0.7579    0.7323
%     0.8697    0.7435    0.7192
%     0.8774    0.7290    0.7062
%     0.8850    0.7146    0.6932
%     0.8927    0.7002    0.6802
%     0.9004    0.6857    0.6672
%     0.9080    0.6713    0.6542
%     0.9157    0.6568    0.6412
%     0.9234    0.6424    0.6281
%     0.9310    0.6280    0.6151
%     0.9387    0.6135    0.6021
%     0.9463    0.5991    0.5891
%     0.9540    0.5846    0.5761
%     0.9617    0.5702    0.5631
%     0.9693    0.5558    0.5501
%     0.9770    0.5413    0.5370
%     0.9847    0.5269    0.5240
%     0.9923    0.5124    0.5110
%     1.0000    0.4980    0.4980];
% 
% end

%%
function cmap=cm_hot

cmap=[
             0         0         0
    0.2980         0         0
    0.3181         0         0
    0.3382         0         0
    0.3582         0         0
    0.3783         0         0
    0.3983         0         0
    0.4184         0         0
    0.4384         0         0
    0.4585         0         0
    0.4785         0         0
    0.4986         0         0
    0.5187         0         0
    0.5387         0         0
    0.5588         0         0
    0.5788         0         0
    0.5989         0         0
    0.6189         0         0
    0.6390         0         0
    0.6590         0         0
    0.6791         0         0
    0.6992         0         0
    0.7192         0         0
    0.7393         0         0
    0.7593         0         0
    0.7794         0         0
    0.7994         0         0
    0.8195         0         0
    0.8396         0         0
    0.8596         0         0
    0.8797         0         0
    0.8997         0         0
    0.9198         0         0
    0.9398         0         0
    0.9599         0         0
    0.9799         0         0
    1.0000         0         0
    1.0000    0.0270         0
    1.0000    0.0541         0
    1.0000    0.0811         0
    1.0000    0.1081         0
    1.0000    0.1351         0
    1.0000    0.1622         0
    1.0000    0.1892         0
    1.0000    0.2162         0
    1.0000    0.2432         0
    1.0000    0.2703         0
    1.0000    0.2973         0
    1.0000    0.3243         0
    1.0000    0.3514         0
    1.0000    0.3784         0
    1.0000    0.4054         0
    1.0000    0.4324         0
    1.0000    0.4595         0
    1.0000    0.4865         0
    1.0000    0.5135         0
    1.0000    0.5405         0
    1.0000    0.5676         0
    1.0000    0.5946         0
    1.0000    0.6216         0
    1.0000    0.6486         0
    1.0000    0.6757         0
    1.0000    0.7027         0
    1.0000    0.7297         0
    1.0000    0.7568         0
    1.0000    0.7838         0
    1.0000    0.8108         0
    1.0000    0.8378         0
    1.0000    0.8649         0
    1.0000    0.8919         0
    1.0000    0.9189         0
    1.0000    0.9459         0
    1.0000    0.9730         0
    1.0000    1.0000         0
    1.0000    1.0000    0.0385
    1.0000    1.0000    0.0769
    1.0000    1.0000    0.1154
    1.0000    1.0000    0.1538
    1.0000    1.0000    0.1923
    1.0000    1.0000    0.2308
    1.0000    1.0000    0.2692
    1.0000    1.0000    0.3077
    1.0000    1.0000    0.3462
    1.0000    1.0000    0.3846
    1.0000    1.0000    0.4231
    1.0000    1.0000    0.4615
    1.0000    1.0000    0.5000
    1.0000    1.0000    0.5385
    1.0000    1.0000    0.5769
    1.0000    1.0000    0.6154
    1.0000    1.0000    0.6538
    1.0000    1.0000    0.6923
    1.0000    1.0000    0.7308
    1.0000    1.0000    0.7692
    1.0000    1.0000    0.8077
    1.0000    1.0000    0.8462
    1.0000    1.0000    0.8846
    1.0000    1.0000    0.9231
    1.0000    1.0000    0.9615
    1.0000    1.0000    1.0000];

end