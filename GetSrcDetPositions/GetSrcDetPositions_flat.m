 function pos = GetSrcDetPosisions_flat(fldr,  flnm_h, fldr_pos, pixel_size,...
    ListSrc, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: GetSrcDetPositions_flat
% # determine model positions of sources and detectors from a measurment
% # for Pioneer system with the flat probe
% for NIRFAST
% # for Nirfast, the voxel size is 1 mm, 
% # the default model center is (45, 45)mm
% # input: fldr: path for the measured data
%          flnm_h:  path for timing_data
%          fldr_pos: path for the folder to save source detector positions
%          pixel_size: in mm
%          ListSrc: sources [1:11]
%          
% created: 2019.11.06 jingjing jiang jing.jing.jiang@outlook.ch
% modified: 2020.11.27 Jingjing Jiang 
%           2021.05.10 Jingjing Jiang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varagin
isRep = 0;
ratioR = 0.87;
fovC = [16 16];
modelCenterRound = [45 45];
if ~isempty(varargin)
    if length(varargin) >= 1
        % log file name
        if isnumeric(varargin{1})  
           if varargin{1}>=0
            isRep = 1;
            rep_id = varargin{1};
            end
        else
            error('Bad 6th argument value. number expected for repetition number.')
        end
    end
    if length(varargin) >= 2
        % log file name
        if isnumeric(varargin{2}) 
            ratioR = varargin{2};
        else
            error('Bad 7th argument value. number expected for radius ratio.')
        end
    end
    if length(varargin) >= 3
        % log file name
        if isnumeric(varargin{3}) 
            fovC = varargin{3};
        else
            error('Bad 8th argument value. number expected for the center of fov.')
        end
    end
    if length(varargin) >= 4
        % log file name
        if isnumeric(varargin{4}) 
            modelCenterRound = varargin{4};
        else
            error('Bad 9th argument value. number expected for the center of fov.')
        end
    end    
end
% get detector source positions  
if exist(fldr_pos,'dir')~=7
    mkdir(fldr_pos);
end
cw_allsrc = zeros(32);
nei = [-1 -1
             -1  0 
             -1  1
              0 -1
              0  1
              1  -1
              1  0
              1  1
             ];
unitmm = pixel_size;
% modelCenterRound = [50, 50] ;
% modelCenterRound = [45, 45] ;

for isrc = ListSrc
    % load timing responses 
    if isRep
        load([fldr flnm_h{1} num2str(isrc) '_' num2str(rep_id) flnm_h{2}]);
    else
        load([fldr flnm_h{1} num2str(isrc) flnm_h{2}]);
    end
    
    % get CW data
     if iscell(timing_response) % old data type before 2021.04

%         [nc nr] = size(timing_response);
%         [~,Y] = meshgrid(0.5:nc-0.5, 0.5:nr-0.5);
%         kk = 0.25;
%         [Xq,Yq] = meshgrid(0.5:kk:nc-0.5, 0.5:kk:nr-0.5);
             for ic = 1:nc
                for ir = 1:nr
                    CW_raw(ir, ic) = sum(timing_response{ir,ic});
                end
             end
             CW = CW_raw;
        for ic = 1:nc
            for ir = 1:nr 
                if((ir-16)^2+(ic-16 )^2 < 11.5*11.5)  & (CW(ir,ic)==0 ||isnan(CW(ir,ic))) 
                    c_n = 0;
                    v_n = 0;
                    for nn = 1:8
                        if  CW_raw(ir+nei(nn,1), ic+nei(nn,2)) &  ...
                                ~ isnan( CW_raw(ir+nei(nn,1), ic+nei(nn,2)))
                            v_n = v_n +  CW_raw(ir+nei(nn,1), ic+nei(nn,2)) ;
                            c_n = c_n+1;
                        end
                    end
                    CW(ir, ic) = v_n / c_n;
                end
            end
        end
     else % new data type: matrix,  after 2021.04
        [ng nc nr] = size(timing_response);
%         [~,Y] = meshgrid(0.5:nc-0.5, 0.5:nr-0.5);
%         kk = 0.25;
%         [Xq,Yq] = meshgrid(0.5:kk:nc-0.5, 0.5:kk:nr-0.5);
        CW_raw  = squeeze(sum(timing_response,1));      
        CW = CW_raw;
        for ic = 1:nc
            for ir = 1:nr 
%                 if((ir-16)^2+(ic-16 )^2 < 11.5*11.5)  & ...
%                         ( ~ CW(ir,ic)==0)
                if((ir-fovC(1))^2+(ic-fovC(2))^2 < 12*12)  & ...
                        ( CW(ir,ic)==0)
             
                    c_n = 0;
                    v_n = 0;
                    for nn = 1:8
                        if  CW_raw(ir+nei(nn,1), ic+nei(nn,2)) &  ...
                                 ( CW_raw(ir+nei(nn,1), ic+nei(nn,2)))
                            v_n = v_n +  CW_raw(ir+nei(nn,1), ic+nei(nn,2)) ;
                            c_n = c_n+1;
                        end
                    end
                    CW(ir, ic) = v_n / c_n;
                end
            end
        end                       
     
    end
  %  interpolation
%     Vq = interp2(X,Y,CW,Xq,Yq,'spline');      
%     [rr, cc] = find(ismember(Vq, max(CW(:))));
%     row_max(isrc) = Xq(rr,cc);
%     col_max(isrc) = Yq(rr,cc);
    threshFoV = 98;
    thres =   prctile(CW(:),threshFoV); % select region of FOV, 95%
    mask = CW > thres;
    CW_masked = CW .* mask;
    stats = regionprops(mask);
    centroid = stats.Centroid;
    row_max(isrc) = centroid(1);
    col_max(isrc) = centroid(2);
    figure(101)
    subplot(331)
    imagesc(CW_raw)
    ax = gca;
    ax.YDir = 'normal'
    title('original')
    subplot(332)
    imagesc( CW)
    ax = gca;
    ax.YDir = 'normal'
    title('fill in gaps')
    
   
%     subplot(323)
%     imagesc(  Vq)
%     title('interpolation')
    subplot(333)
    imagesc(CW_masked)
    ax = gca;
    ax.YDir = 'normal'
    title('selected region')
 
end

% linear regression to fit a circle
[R,XC,YC,ERR] = circfit(row_max,col_max);
  

figure(101),subplot(334)
plot(row_max, col_max,'ro');
hold on
plot(row_max(1), col_max(1),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6])
plot(row_max(2), col_max(2),'ro','MarkerSize',8,'MarkerFaceColor',[1 .6 .6])
axis equal

plot(XC, YC,'mo');
viscircles([XC, YC], R)
title('fitted relative source positions')
%% process source data
% load relative source positions
% load('METRICS.mat','PosList')
PosList =[39    54
    26    50
    17    41
    15    28
    20    16
    31     8
    44     8
    56    15
    61    26
    60    39
    53    49
    38    31]; %% the last is for the central position
%%
row2 = PosList(:,1);
col2 = PosList(:,2);

row2N = row2(1:end-1) -row2(end);
col2N = col2(1:end-1) -col2(end);

figure(101),subplot(335)
hold on
plot(row2N, col2N, 'bo')
plot(row2N(1), col2N(1),'bo','MarkerSize',10,'MarkerFaceColor',[.6 1  .6])
plot(row2N(2), col2N(2),'bo','MarkerSize',8,'MarkerFaceColor',[.6 1  .6])
 title('measured relative source positions')

%% apply rotation to source positions 
[R2,XC2,YC2,ERR2] = circfit(row2N,col2N); %[unit: mm]
for isrc = ListSrc
    vector_det = [row_max(isrc)-XC col_max(isrc)-YC];
    vector_2 = [row2N(isrc)-XC2 col2N(isrc)-YC2];
    ang(isrc) = acos(dot(vector_2,vector_det)/norm(vector_2)/norm(vector_det));
end
ang
theta = median(ang) /pi * 180  % changed by jingjing 2019.11.06
% theta = 180;
% Create rotation matrix
A = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% Rotate 
PosScrRot= (A*[row2N col2N]')';
figure(101),
subplot(336)
hold on
plot(PosScrRot(:,1), PosScrRot(:,2),'yo','MarkerSize',4,'MarkerFaceColor',[.6 .6 1  ])
plot(PosScrRot(1,1), PosScrRot(1,2),'yo','MarkerSize',10,'MarkerFaceColor',[.6 .6 1  ])
plot(PosScrRot(2,1), PosScrRot(2,2),'yo','MarkerSize',8,'MarkerFaceColor',[.6 .6 1  ])
% center 
[R2,XC2,YC2,ERR2] = circfit(PosScrRot(:,1),PosScrRot(:,2));
figure(101),subplot(336)
viscircles([XC2, YC2], R2,'EdgeColor','b');


figure(101),subplot(336)
hold on
plot(PosScrRot(:,1), PosScrRot(:,2),'ko')
plot(PosScrRot(1,1), PosScrRot(1,2),'ko','MarkerSize',10,'MarkerFaceColor',[.6 .6 1  ])
plot(PosScrRot(2,1), PosScrRot(2,2),'ko','MarkerSize',8,'MarkerFaceColor',[.6 .6 1  ])
for isrc = ListSrc
    plot([PosScrRot(isrc,1) XC2], [PosScrRot(isrc,2) YC2],'-g')
end
%% get model coordinates of detectors 
% select detectors
% origin of Piccolo coordinate = [0 0], the first pixel coordinate is [0.5 0.5]
x = [1:32]  ;
[X,Y] = meshgrid(x,x)  ;
R0 = R * ratioR;
% IndexInFOVPiccolo =  (XC - X).^2 + (YC - Y).^2 <= R0^2; 
IndexInFOVPiccolo =  (YC - X).^2 + (XC - Y).^2 <= R0^2; % changed by
% jingjing 2021.06.03
[xFOV, yFOV] = find(IndexInFOVPiccolo);
PosDetCam = [xFOV yFOV];
xFOVP = xFOV - 0.5;
yFOVP = yFOV - 0.5;
% transform from Piccolo coordinate to model coordinate
% introduce a parameter: pixel size 

% modelCenter = modelCenterRound + [XC YC] - floor([XC YC]);
modelCenter = modelCenterRound + [YC XC] - floor([YC XC]); % changed by
% % jingjing 2021.06.03
[xFOV, yFOV] = find(IndexInFOVPiccolo);

% select detectors
 
%% move to [0 0 ]
xDet0 = xFOVP - XC;
yDet0 = yFOVP - YC;% changed by
% jingjing 2021.06.03
% xDet0 = xFOVP - YC; % changed on 2019.11.06 by jingjing 
% yDet0 = yFOVP - XC;
% figure(31),hold on, plot(xDet0, yDet0,'o')

% %% pos x 2 for MC model
% xDet0_2x = xDet0*2;
% yDet0_2x = yDet0.*2;
% %% get model coordinates of detectors
% PosDetModel = [xDet0_2x  yDet0_2x ] + ...
%     repmat(2.* modelCenter, size(xFOV ,1),1 );
% %% get model coordinates of sources
% PosScrRot_2x = PosScrRot ./unitmm;
% PosSrcModel = PosScrRot_2x + ...
%     repmat(2.* modelCenter, size(PosScrRot,1),1 ); % [unit: pixel/2]
% 
% % dis_src  =sqrt(sum( (PosSrcModel - repmat(2.* modelCenter, size(PosScrRot,1),1 )).^2,2));
% [XC YC] % ceter of ring in piccolo coordinate [32x32 pixels]
% % source positions are easily determined
% % detector potions are very difficult because the FOV is not fixed and 
% % pixel size is unknown 
%  
% %% plot source and detector positions for MC
% figure(101),subplot(337)
% plot(PosSrcModel(:,1), PosSrcModel(:,2),'mo')
% hold on
% plot(PosDetModel(:,1), PosDetModel(:,2),'o')
%  
% plot(PosSrcModel(1,1), PosSrcModel(1,2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6   ])
% plot(PosSrcModel(2,1), PosSrcModel(2,2),'ro','MarkerSize',7,'MarkerFaceColor',[1 .6 .6   ])
% plot(PosDetModel(1,1), PosDetModel(1,2),'ro','MarkerSize',10,'MarkerFaceColor',[  .6 .6  1 ])
% plot(PosDetModel(2,1), PosDetModel(2,2),'ro','MarkerSize',7,'MarkerFaceColor',[ .6 .6   1])
% legend('sources','detectors', 'source 1', 'source 2', 'detector 2','detector 2')
% 
% xlabel(['x [' num2str(pixel_size) '/2 mm]'])
% ylabel(['y [' num2str(pixel_size) '/2 mm]'])
% grid on
% title('MC model source / detector positions')
% %% save data
% fln_out = 'detSrcPosModel.mat';
% save([fldr_pos '/' fln_out], 'PosSrcModel', 'PosDetModel', ...
%     'pixel_size','modelCenter','unitmm','IndexInFOVPiccolo', 'PosDetCam');
% disp(['MC model positions are saved to ' [fldr_pos '/' fln_out]])

%% source detector positions for NIRFAST 
% the model center is (45, 45) by default
% for modeling in NIFAST, the corresponding phantom size is 90 x 90 x z
% shift them when it is needed to create models of different sizes
PosSrcFEM  = PosScrRot + ...
    repmat(modelCenterRound, size(PosScrRot,1),1 ); % [mm]
PosDetFEM = [xDet0  yDet0 ] .*pixel_size+ ...
    repmat( modelCenterRound, size(xFOV ,1),1 );

% %% x,y -> y, x % commented out by jingjing 2020.10.24
% temp = PosSrcFEM(:,1);
% PosSrcFEM(:,1) = PosSrcFEM(:,2);
% PosSrcFEM(:,2) = temp;

%% plot source and detector positions
figure(101),subplot(338)
plot(PosSrcFEM(:,1), PosSrcFEM(:,2),'mo')
hold on
plot(PosDetFEM(:,1), PosDetFEM(:,2),'o')
 
plot(PosSrcFEM(1,1), PosSrcFEM(1,2),'ro','MarkerSize',10,'MarkerFaceColor',[1 .6 .6   ])
plot(PosSrcFEM(2,1), PosSrcFEM(2,2),'ro','MarkerSize',7,'MarkerFaceColor',[1 .6 .6   ])
plot(PosDetFEM(1,1), PosDetFEM(1,2),'ro','MarkerSize',10,'MarkerFaceColor',[  .6 .6  1 ])
plot(PosDetFEM(2,1), PosDetFEM(2,2),'ro','MarkerSize',7,'MarkerFaceColor',[ .6 .6   1])
legend('sources','detectors', 'source 1', 'source 2', 'detector 2','detector 2')

xlabel(['x  [mm]'])
ylabel(['y [mm]'])
grid on
title('model source / detector positions')

% save to fldr_pos

fln_out = 'detSrcPosFEM.mat';
save([fldr_pos '/' fln_out],  'PosSrcFEM', 'PosDetFEM',...
    'pixel_size','modelCenter','unitmm','IndexInFOVPiccolo', ...
    'PosDetCam');
disp(['model positions are saved to ' [fldr_pos '/' fln_out]])
 
pos.det.coordinates = PosDetFEM;
pos.src.coordinates = PosSrcFEM;
pos.unitmm = unitmm;
pos.modelCenter = modelCenter;
pos.camera = PosDetCam;
pos.det.pixel = pixel_size;
fln_out = 'posSrcDet_Flat.mat';
save([fldr_pos '/' fln_out], 'pos')

end