function [J,data] = jacobian_stnd_TR_allBins(mesh, t_range, t_step, varargin)
% created by Hamid on 2019.05
% modified by jingjing on 2019.09.04
% 
narginchk(1,5);
nargoutchk(0,2);

%% get optional inputs, handle variable input

% default
 
second_mesh_basis = [];
solver = get_solver;
OPTIONS = solver_options;
isAll = false;
 

if ~isempty(varargin)

    if length(varargin) >= 1
        % second mesh basis
        if isstruct(varargin{1}) || ~isempty(varargin{1})
            second_mesh_basis = varargin{1};
        else
            warning('Bad 4nd argument value. A mesh structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 2
        % solver
        if ischar(varargin{2}) || isstring(varargin{2})
            % user specified, sanity check
            solver = get_solver(varargin{2});
        elseif isstruct(varargin{2})
            OPTIONS = varargin{2};
        elseif isempty(varargin{2})
            solver = get_solver;
        else
            warning('Bad 5th argument value. Solver name or solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 3
        % solver options
        if isstruct(varargin{3})
            OPTIONS = varargin{3};
        elseif isempty(varargin{3})
            OPTIONS = solver_options;
        else
            warning('Bad 6th argument value. Solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 4
        % if scattering for CW data as well
        if ischar(varargin{4}) || isstring(varargin{4})
            if strcmp(varargin{4},'all')
                isAll = true;
            else
                warning('Bad 7th argument value. Text expected. Please see help for details on how to use this function.')
            end
        end
    end
     
    if length(varargin) > 5
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end
% JACOBIAN_STND_TR Calculates spatial distributions of sensitivity for TR
%% ------- added by jingjing 2021.03.17 --------------------------%
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;

[data] = femdata_stnd_TR(mesh,t_range,t_step,'field', 'BiCGStab_GPU',OPTIONS);
% [data] = femdata_stnd_TR(mesh,t_range,t_step,'field', 'BiCGStab_CPU',OPTIONS);

%% comment out the original 
%[data] = femdata_stnd_TR(mesh,t_range,t_step,'field');  
%% !MELLIN! Now calculate adjoint sources and data (detectors are new sources)

% make sources copy and swap with detectors
sources_copy = mesh.source;
detectors_copy = mesh.meas;
link_copy = mesh.link;

% swap sources
mesh.source.fixed = 0;
mesh.source.num = detectors_copy.num;
mesh.source.coord = detectors_copy.coord;
mesh.source.fwhm = zeros(size(detectors_copy.num));
mesh.source.int_func = detectors_copy.int_func;
%swap detector
mesh.meas.fixed = 0;
mesh.meas.num = sources_copy.num;
mesh.meas.coord = sources_copy.coord;
% a missing detectors int functions 'mesh.meas.int_func' will calculate
% itself in 'get_boundary_data' if needed (however, we don't need/use
% adjoint boundary data and 'mesh.meas.int_func' is not needed) 

%swap link
mesh.link(:,1) = link_copy(:,2);
mesh.link(:,2) = link_copy(:,1);

% calculate adjoint field for all new sources (old detectors)
% data_detector = femdata_TR(mesh,t_range,t_step,'field',solver_name_matlab_iterative);
% changed by jingjing 2021.08.03
data_detector = femdata_stnd_TR(mesh,t_range,t_step,...
    'field', 'BiCGStab_GPU',OPTIONS);

data.aphi = data_detector.phi;
clear data_detector
% swap back sources, detectors and link
mesh.source = sources_copy;
mesh.meas = detectors_copy;
mesh.link = link_copy;

%do bins
% for i = 1 : size(bins,1)
%     ind(i,1) = find(data.time==bins(i,1));
%     ind(i,2) = find(data.time==bins(i,2));
% end
%% spatially convolve fileds for sources and detectors
% copy boundary data with possible NaN values for disabled pairs
boundary_data_copy = data.tpsf;

% jacobian build requires data for enabled pairs only
data.tpsf = data.tpsf(logical(data.link(:,3)),:);

num_dt = size(0:t_step:t_range,2)-1;
for i = 1 : size(data.phi,2)
    data.phifft(:,i,:) = fft(squeeze(data.phi(:,i,:)),2*num_dt-1,2);
end
for i = 1 : size(data.aphi,2)
    data.aphifft(:,i,:) = fft(squeeze(data.aphi(:,i,:)),2*num_dt-1,2);
end

% J = build_jacobian_TR(mesh,data,'all');
if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
    % build jacobian on data interpoladed onto the coarse mesh
    
        % absorption and scattering
         J  = build_jacobian_TR(second_mesh_basis,interpolatef2r(mesh,...
             second_mesh_basis,data),'all');
 
else
    
     J = build_jacobian_TR(mesh,data,'all');
end
% interpolatef2r(mesh,second_mesh_basis,data);

for i = 1 : size(J.mua,1)
    J.mua(i,:,:) = ifft(squeeze(J.mua(i,:,:))',2*num_dt-1)';
    J.mus(i,:,:) = ifft(squeeze(J.mus(i,:,:))',2*num_dt-1)';
end

Jmua(:,:,:) = J.mua(:,:,1:num_dt);
Jmus(:,:,:) = J.mus(:,:,1:num_dt);

J.mua = Jmua;
J.mus = Jmus;
clear Jmu*

% mus is infact kappa

%% jacobian ln(tpsf) 

[pp nn tt] = size(J.mua);
J.mualnI = zeros(pp,nn,tt);
J.muslnI = J.mualnI;
for p = 1:pp
    for t_range = 1:tt
        if data.tpsf(p,t_range) == 0
        else
            J.mualnI(p,:,t_range) = J.mua(p,:,t_range)./data.tpsf(p,t_range);
            J.muslnI(p,:,t_range) = J.mus(p,:,t_range)./data.tpsf(p,t_range);
        end
    end
end
end


function [data_recon] = interpolatef2r(fwd_mesh,recon_mesh,data_fwd)
% This function interpolates fwd_mesh data onto recon_mesh
% Used to calculate the Jacobian on second mesh
% added by jingjing 2019.09.04
    % copy boundary data
    data_recon.tpsf = data_fwd.tpsf;
    data_recon.time = data_fwd.time;
    % calculate interpolation functions if needed
    if ~isfield(fwd_mesh,'fine2coarse')
        fwd_mesh.fine2coarse = second_mesh_basis(fwd_mesh,recon_mesh);
    end

    % loop through nodes of recon_mesh
    for i = 1:size(recon_mesh.nodes,1)
        % loop throurg moments
%         for ind_tpsf = 1:size(data_fwd.phi,3)
        for ind_tpsf = 1:size(data_fwd.phifft,3)

            % if the recon_mesh node belongs to the fwd_mesh, interpolate
            if fwd_mesh.fine2coarse(i,1) > 0
%                 data_recon.phi(i,:,ind_tpsf) = fwd_mesh.fine2coarse(i,2:end) * ...
%                     squeeze(data_fwd.phi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:,ind_tpsf));
%                 data_recon.aphi(i,:,ind_tpsf) = fwd_mesh.fine2coarse(i,2:end) * ...
%                     squeeze(data_fwd.aphi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:,ind_tpsf));
                data_recon.phifft(i,:,ind_tpsf) = fwd_mesh.fine2coarse(i,2:end) * ...
                    squeeze(data_fwd.phifft(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:,ind_tpsf));
                data_recon.aphifft(i,:,ind_tpsf) = fwd_mesh.fine2coarse(i,2:end) * ...
                    squeeze(data_fwd.aphifft(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:,ind_tpsf));
            end
        end
    end

end
