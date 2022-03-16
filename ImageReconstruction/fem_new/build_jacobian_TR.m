function [J] = build_jacobian_TR(mesh,data,varargin)
% created by jingjing 2019.09.04
%
%% check in/out

narginchk(2,3);
nargoutchk(0,1);

%% get optional inputs, handle variable input

% default
isAll = false;

if ~isempty(varargin)
    if length(varargin) == 1
        if ischar(varargin{1}) || isstring(varargin{1})
            if strcmp(varargin{1},'all')
                isAll = true;
            end
        else
            error('Bad 3rd argument value. Text expected. Please see help for details on how to use this function.')
        end
    else
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% BODY

% indexes of active pairs
index_pairs = find(logical(mesh.link(:,3)));
% [labels of active sources, index of source per active source-detector pair, ~]
[index_source_num, index_sources, ~] = sources_active_index(mesh);
% [labels of active detectors, index of detector per active source-detector pair, ~]
[index_detector_num, index_detectors, ~] = detectors_active_index(mesh);

% required as phi and aphi columns exist only for active and unique sources
% and detectors! Forward and adjoint fields are calculated for sources and
% detectors as calculated by 'sources_active_index' and
% 'detectors_active_index'
% indexes of columns in phi (sources permutatin matrix)
source_index = zeros(length(index_pairs),1);
% indexes of columns in aphi (detectors permutatin matrix)
detector_index = zeros(length(index_pairs),1);

% loop through active source-detectors pairs
for ind_active_pair = 1:length(index_pairs)
    % required as phi columns exist only for active and unique sources! 
    mask_source = index_source_num == index_sources(index_pairs(ind_active_pair));
    mask_detector = index_detector_num == index_detectors(index_pairs(ind_active_pair));
    source_index(ind_active_pair) = find(mask_source);
    detector_index(ind_active_pair) = find(mask_detector);
end

num_dt = size(data.time,2);
for i = 1 : 2*num_dt-1
%     % spatial convolution of fields for all enabled source-dtector pairs
%     J.mua(:,:,i) = -IntFG_tet4_CPU(mesh,...
%                                    squeeze(data.phifft(:,:,i)),...
%                                    squeeze(data.aphifft(:,:,i)),...
%                                    source_index,detector_index)';
%     % spatial convolution of gradients of fields for all enabled source-dtector pairs
%     J.mus(:,:,i) = -IntgradFgradG_tet4_CPU(mesh,...
%                                    squeeze(data.phifft(:,:,i)),...
%                                    squeeze(data.aphifft(:,:,i)),...
%                                    source_index,detector_index)';
    % spatial convolution of fields for all enabled source-dtector pairs
    J.mua(:,:,i) = IntFG_tet4_CPU(mesh,...
                                   squeeze(data.phifft(:,:,i)),...
                                   squeeze(data.aphifft(:,:,i)),...
                                   source_index,detector_index)';
    % spatial convolution of gradients of fields for all enabled source-dtector pairs
    J.mus(:,:,i) = IntgradFgradG_tet4_CPU(mesh,...
                                   squeeze(data.phifft(:,:,i)),...
                                   squeeze(data.aphifft(:,:,i)),...
                                   source_index,detector_index)';
end