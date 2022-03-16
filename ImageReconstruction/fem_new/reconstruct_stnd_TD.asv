function [mesh, pj_error] = reconstruct_stnd_TD(mesh, t_cfg, ...
    data_measured, varargin)
% reconstruct_stnd_TD Recovers spatial distributions of optical properties
% using temporal data
% author: Jingjing Jiang jing.jing.jiang@outlook.com
% created  2022.03.16 
%
% adapted from RECONSTRUCT_STND_FD
% 
% SYNTAX:
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED)
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER)
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL)
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH)
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX)
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX,
%                              REG_THIKONOV)
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX,
%                              REG_THIKONOV, LOG_FILENAME) 
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX,
%                              REG_THIKONOV, LOG_FILENAME, FILTER_FLAG)
%  [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, [], [], [], [], [], [], [])
%  [MESH, PJ_ERROR] = RECONSTRUCT_STND_FD(___)
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED) Recovers for
%   optical properties within the MESH (FEM mesh in NIRFAST format or mesh
%   name within the MATLAB search path) considering sources modulated at
%   FREQUENCY (>=0Hz) in Hz and boundary data DATA_MEASURED. FREQUENCY=0Hz
%   enables Continuous Wave (CW) analysis.
% WARNING:
%   * DATA_MEASURED should be calibrated in respect to diffusion theory (or
%     MC) and expressed in photons per mm^-2s^-1. See help for calibration
%     procedures to learn how it is done.
%   **By default, this syntax will trigger absolute recovery of absorption
%     parameter. Regardles of the FREQUENCY value.
% INPUT:
%   DATA_MEASURED should be formated as data structure returned by
%   'femdata_stnd_FD' with following required entries:
%    - DATA_MEASURED.LINK - shows how sources and detectors are linked
%       together into pairs. Formated as 'mesh.link' field. See the
%       'load_mesh', 'save_mesh' functions help for information on the
%       'link' field format. 'DATA_MEASURED.LINK' should match the
%       'MESH.link'.
%    - DATA_MEASURED.AMPLITUDE - m by p matrix with absolute value of the
%       measured photon fluence rate. Calibrated using theory and expressed
%       in mm^-2s^-1. m is number of measurements (source-detector pairs)
%       as enabled in the 'DATA_MEASURED.LINK' and p is number of
%       experiments. The experiment is a dataset used for reconstruction.
%       The reconstruction procedure will bw carried out p times.
%    - DATA_MEASURED.PHASE - m by p matrix with angle (polar coordinate) of
%       the measured photon fluence rate and expressed in degrees. Ignored
%       for CW sources (if FREQUENCY=0Hz). m and p size the same as for
%       'DATA_MEASURED.AMPLITUDE'
%   DATA_MEASURED can also be a name of data file in NIRFAST data file
%   format. The file shoud exist in the MATLAB search path.
% OUTPUT:
%   The returned MESH structure is a copy of the input MESH except the
%   returned structure have the 'MESH.mua' and/or 'MESH.mus' and/or
%   'MESH.kappa' filled with the recovered parameters. the recovered fields
%   have n by p size, where n is number of nodes and p is the number of
%   experiments as specified in the DATA_MEASURED.
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER)
%   The PARAMETER is used to force a recovery for particular optical
%   parameter(s). PARAMETER have follwoing minigfull text values:
%    - 'mua' - (default) recover for absorption only
%    - 'mus' -  recover for scattering only
%    - 'all' - (default) absorption
%   Other values are ignored. The measurement type (CW or FD) does not
%   determine the PARAMETER option and 'mus' or 'all' can be used for CW
%   sources as well. However, the 'all' option for the CW should be rather
%   used to illustrate how the CW measurements can not distinguish between
%   absorption and scattering at the same time.
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL)
%   ISDIFFERENTIAL is a text parameter that have only one meaningfull
%   value:
%    - 'differential' - provided data (DATA_MEASURED) are for differen
%       imagingand as such the differential imaging reconstruction will be
%       carried out
%    - all other entries are ignored and absolute recovery is performed
%   By default, the absolute recovery is carried out.
% WARNING:
%   In the differential mode, the DATA_MEASURED should be specified as
%   differential data in the following manner:
%    - DATA_MEASURED.LINK - as in the absolute imaging.
%    - DATA_MEASURED.AMPLITUDE - m by p matrix with attenuations A between
%       some rest and stimulation responses. Assuming that light intensities
%       I_rest and I_stim were measured during the rest and stimulation
%       periods repsectively, the attenuation will be defined as follows:
%       A=ln(I_rest/I_stim). m is number of measurements (source-detector
%       pairs) as enabled in the 'DATA_MEASURED.LINK' and p is number of
%       experiments. The experiment is a dataset used for reconstruction.
%       The reconstruction procedure will bw carried out p times.
%    - DATA_MEASURED.PHASE - m by p matrix with changes in angle d_phi
%       (polar coordinate) between some rest and stimulation periods,
%       expressed in degrees. Assuming that the phase shift phi_rest and
%       phi_stim were measured during the rest and stimulation periods
%       repsectively, the change in the phase shift will be defined as
%       follows: d_phi=phi_stim-phi_rest. Ignored for CW sources (if
%       FREQUENCY=0Hz). m and p sizes are the same as for
%       'DATA_MEASURED.AMPLITUDE' 
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH)
%   RECON_MESH is the second mesh basis used for reconstruction. The idea
%   is to carry out the reconstruction on the RECON_MESH and interpolate
%   the result back onto the MESH. I ngeneral the RECON_MESH is designed to
%   be a coarse version of the MESH. The RECON_MESH can be specified in 4
%   ways:
%    - as a mesh structure in NIRFAST format
%    - as a name of a '*.mat' file that holds the mesh structure in NIRFAST
%      format where the mesh should be named 'recon_mesh' as only this
%      variable will be loaded from this file.
%    - as a root file name of mesh saved in NIRFAST format and loaded using
%      the 'load_mesh' function.
%    - as a pixel/voxel basis specified as 3 or 2 elements vector: [Nx Ny]
%      for 2D meshes or [Nx Ny Nz] for 3D meshes. Where Nx, Ny and Nz
%      specify number of evenly spaced nodes within the MESH bounding-box.
%      the RECON_MESH will be created by interpolating the MESH into
%      Nx*Ny*(Nz) evely spaced nodes.
%   Empty by default and the reconstruction is carried out on the MESH.
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX)
%   ITER_MAX (default 30) is the maximum iteration steps.
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX,
%                              REG_THIKONOV)
%   REG_THIKONOV (default 10) is the initial regularization value (usually
%   between 1 and 1000). 10 seems to be a good choice for adult human brain
%   imaging. Unitless.
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX,
%                              REG_THIKONOV, LOG_FILENAME)
%   LOG_FILENAME (empty by default) is a root file name to save the
%   recovery log information as displayed on the MATLAB Command Window. The
%   Command Window output will be saved in the 'LOG_FILENAME.log' file.
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX,
%                              REG_THIKONOV, LOG_FILENAME, FILTER_FLAG)
%   FILTER_FLAG (default 0) is the flag enabling spatial filtering of
%   recovered optical properties. The filtering is carried out at each
%   reconstruction step. FILTER_FLAG is a scalar and is interpreted as
%   follows:
%    - 'FILTER_FLAG=0' - (default) do not filter
%    - 'FILTER_FLAG<0' - apply median spatial filtering as follows: for
%      each node set median value of all connected neighbours. Repeat
%      abs(FILTER_FLAG) times.
%    - 'FILTER_FLAG>0' - apply mean spatial filtering as follows: for each
%      node set mean value of all connected neighbours. Repeat
%      abs(FILTER_FLAG) times.
% 
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, PARAMETER,
%                              ISDIFFERENTIAL, RECON_MESH, ITER_MAX,
%                              REG_THIKONOV, LOG_FILENAME, FILTER_FLAG, ...
%                              IsMellin)
%      - IsMellin = true , mellin is required
% [MESH] = RECONSTRUCT_STND_FD(MESH, FREQUENCY, DATA_MEASURED, [],
%                              [], [], [], [], [], FILTER_FLAG)
%   Optinal parameters can be empty. Empty entries are neglected, default
%   values are used.
% 
% [MESH, PJ_ERROR] = RECONSTRUCT_STND_FD(___) Returns the projection error
%   (PJ_ERROR) expressed as total square error of recovery. The PJ_ERROR is
%   a matrix of Ni by Ne size, where Ni is number of finished iterations up 
%   to the ITER_MAX and Ne is number of experiments (data setups as
%   specified in DATA_MEASURED).
% 
% ALGORITHM
%   This functions solves iteratively the following pertuirbation problem
%   with Rytov approach. In this example, for simplicity, we assume
%   FREQUENCY=0:
%    1: start with guess: mua = quess
%    2: while error change < -2%
%    3:   calculate sensitivity for current mua: J = J_mua^A
%    4:   calculate boundary data for current mua (I_update)
%    5:   calculate mismatch in boundary data (-A) and the error change
%    6:   if error change >= -2%
%    7:     break
%    8:   end
%    9:   solve for update d_mua: -A = J*d_mua' => d_mua' = J'*(J*J'+lambda*I)\(-A)
%   10:   apply update: mua = mua + d_mua;
%   11: end
%   Where: 'guess' is the initial abosrption from 'MESH.mua', the Jacobian
%     J and current boundary data (I_update) are returned by 'jacobian_FD',
%     the data mishamtch is -A=ln(I_update/I_measure) and the I_measure
%     is the function input 'DATA_MEASURED', the error is sum((-A).^2)
%     and the error change is the change between consecutive steps, I is
%     the identity matrix and lambda is the Thikonov regularization
%     parameter that reduces at each iteration by a factor of 10^0.25.
%     '-A' notation is used to emphasize resemblance with attenuation
%     defined as ln(reference/update) and the minus sign swaps numerator
%     with denominator under the logarithm to ln(update/reference) as
%     required in the Rytov perturbation approach.
% 
% READ: Arridge S. and M. Schweiger Photon-measurement density functions.
%       Part 2: Finite-element-method calculations. APPLIED OPTICS, 34,
%       1995, p. 8026-8037
%       Arridge S.Optical tomography in medical imaging. Inverse Problems,
%       15, 1999, p. R41�R93
%       Dehghani H. et al. Near infrared optical tomography using NIRFAST: algorithm
%       for numerical model and image reconstruction. Communications in Numerical
%       Methods in Engineering, 25, 2008, p. 711�732.
%  
% See also RECONSTRUCT_FD, RECONSTRUCT_SPEC_FD, JACOBIAN_FD, FEMDATA_FD.
% 
%   Part of NIRFAST package.
%   H. Dehghani and S. Wojtkiewicz 2018

%% check in/out

narginchk(3,11);
nargoutchk(0,2);

%% get optional inputs, handle variable input

% default
parameter = 'mua';
isDifferential = false;
isReconMesh = false;
recon_mesh = [];
iter_max = 30;
reg_thikonov = 10;
log_filename = [];
filter_flag = 0;
isMellin = false;

if ~isempty(varargin)
    if length(varargin) >= 1
        % log file name
        if ischar(varargin{1}) || isstring(varargin{1})
            if strcmp(varargin{1},'mua') || strcmp(varargin{1},'mus') || strcmp(varargin{1},'all')
                parameter = varargin{1};
            else
                warning(['Unknown parameter to reconstruct name ''' varargin{1} '''. Default ''mua'' used.']);
            end
        else
            error('Bad 4th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 2
        % if differential imaging
        if ischar(varargin{2}) || isstring(varargin{2}) || isempty(varargin{2}) 
            if strcmp(varargin{2},'differential')
                isDifferential = true;
            end
        else
            error('Bad 5th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 3
        % second mesh basis
        if isstruct(varargin{3}) || isempty(varargin{3})
            recon_mesh = varargin{3};
        elseif (numel(varargin{3}) >= 2) && (numel(varargin{3}) <= 3)
            recon_mesh = varargin{3};
        else
            error('Bad 6th argument value. A mesh structure or pixel/voxel size expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 4
        % max number of iterations
        if ~ischar(varargin{4}) && ~isstring(varargin{4}) && (numel(varargin{4})<=1)
            % if not set to empty on purpose
            if ~isempty(varargin{4})
                % sanity check
                if varargin{4} <= 0
                    error('Nonpositive numer of iterations. Please see help for details on how to use this function.')
                else 
                    iter_max = varargin{4};
                end
            end    
        else
            error('Bad 7th argument value. A scalar expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 5
        % Thikonov regularizatin parameter
        if ~ischar(varargin{5}) && ~isstring(varargin{5}) && (numel(varargin{5})<=1)
            % if not set to empty on purpose
            if ~isempty(varargin{5})
                % sanity check
                if varargin{5} <= 0
                    error('Nonpositive Thikonov regularization parameter. Please see help for details on how to use this function.')
                else 
                    reg_thikonov = varargin{5};
                end
            end    
        else
            error('Bad 8th argument value. A scalar expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 6
        % log file name
        if ischar(varargin{6}) || isstring(varargin{6})
            if ~canwrite(varargin{6})
                warning(['No write permission to save the log ''' varargin{6} '.log'' file. Reconstruction will not be logged.']);
            else
                log_filename = varargin{6};
            end
        else
            error('Bad 9th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) == 7
        % filter flag
        if ~ischar(varargin{7}) && ~isstring(varargin{7}) && (numel(varargin{7})<=1)
            % if not set to empty on purpose
            if ~isempty(varargin{7})
                filter_flag = varargin{7};
            end    
        else
            error('Bad 10th argument value. A scalar expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) > 8
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
        if length(varargin) >= 8
        % log file name
        if ischar(varargin{8}) || isstring(varargin{8})
            if ~strcmp(varargin{8},'mellin')
                warning(['please give a string for Mellin']);
                isMellin = 0;
            else
                isMellin = 1;
            end
        else
            error('Bad 10th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
end

if ~isempty(recon_mesh)
    isReconMesh = true;
end


%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end
% and check the type
if ~strcmp(mesh.type,'stnd')
    error(['Incorrect mesh type ''' mesh.type '''. The ''stnd'' type expected.']);
end
% now check link match
if any(any(data_measured.link ~= mesh.link))
    warning('The ''data.link'' and ''mesh.link'' do not match. Please check your inputs. The result might be faulty or errors might follow.')
end


% %% check frequency
% 
% % check if CW measurement
% isCW = true;
% if frequency > 0
%     isCW = false;
% elseif frequency < 0
%     error('Negative frequency value. Please see help for details on how to use this function.')
% end
% %% max_order
% if max_order < 0
%     error('Negative moment order value. Please see help for details on how to use this function.')
% end
%% load data if needed, check data size
% 'data' is the calibrated experimental data, differential experimental data or simulated data

% load from file if needed
if ischar(data_measured) || isstring(data_measured)
    data_measured = load_data(data_measured);
end

% enabled source-detector pairs
% mask_data = logical(data_measured.link(:,3));

% % number of active measurements: CW mode (amplitude), FD mode (amplitude + phase)
% meas_no = sum(isfinite(data_measured.amplitude(:,1)));
% if ~isCW
%     % double number of measurements for FD (amnplitude and phase)
%     meas_no = meas_no * 2;
% end
% ------------- changed by jingjing jiang 2019.05.22 START ---------------%
% % number of active measurements: 
meas_no = size(data_measured.tpsf,1);
meas_no = meas_no * t_cfg.num_gates;

% -------------- changed by jingjing jiang 2019.05.22 END ----------------%
% numer of experiments
% experiments_no = size(data_measured.moments,2);
experiments_no = 1;
%% PREPARE RECONSTRUCTION BASIS
% Set, load or calculate mesh for reconstruction basis

if isReconMesh
    %first check if load from file
    if ischar(recon_mesh) || isstring(recon_mesh)
        [~,~,fileEext] = fileparts(recon_mesh);
        % check if MATLAB '.mat' or NIRFAST mesh format
        if strcmp(fileEext,'.mat')
            load(recon_mesh,'recon_mesh');
        else
            recon_mesh = load_mesh(recon_mesh);
        end
    elseif size(recon_mesh,2) == mesh.dimension
        %pixel basis
        [mesh.fine2coarse,recon_mesh] = pixel_basis(recon_mesh,mesh);
    end
    % show error if something went wrong....
    if ~isstruct(recon_mesh)
        error('Bad second basis reconstruction mesh. Please see help on how to use the ''RECON_MESH'' parameter.')
    end
    % check if we have all needed interpolation functions
    if ~isfield(recon_mesh,'coarse2fine')
        [mesh.fine2coarse, recon_mesh.coarse2fine] = second_mesh_basis(mesh,recon_mesh);
    end
    if ~isfield(mesh,'fine2coarse')
        mesh.fine2coarse = second_mesh_basis(mesh,recon_mesh);
    end

    % set other parameters for the second mesh.
    recon_mesh.type = mesh.type;
    recon_mesh.link = mesh.link;
    recon_mesh.source = mesh.source;
    recon_mesh.meas = mesh.meas;
    recon_mesh.dimension = mesh.dimension;
    recon_mesh.element_area = ele_area_c(recon_mesh.nodes(:,1:recon_mesh.dimension),recon_mesh.elements);
else
    % no second mesh basis specified
    recon_mesh = mesh;
end

% get number of nodes used for reconstruction
nodes_size = size(recon_mesh.nodes,1);

% allocate reconstructed parameters based on user requests
if strcmp(parameter,'all')
    % absorption and scattering
    mua_recon = zeros(size(mesh.nodes,1),experiments_no);
    mus_recon = zeros(size(mua_recon));
elseif strcmp(parameter,'mus')
    % scattering only
    mus_recon = zeros(size(mesh.nodes,1),experiments_no);
else
    % absorption only
    mua_recon = zeros(size(mesh.nodes,1),experiments_no);
end

%% check for input regularization

lambda.value = reg_thikonov;
% determine regularization type based on the problem to solve size
% JJt -> J*J'; JtJ -> J'*J
if meas_no < nodes_size
    lambda.type = 'JJt';
else
    lambda.type = 'JtJ';
end


%% Initiate projection error and log file

pj_error = Inf(iter_max,experiments_no);

if ~isempty(log_filename)
    % Initiate log file
    fid_log = fopen([log_filename '.log'],'w');
    fprintf(fid_log,'---------------------------------\n');
    fprintf(fid_log,'Forward Mesh                    = %s\n',mesh.name);
    if isReconMesh
        if isfield(recon_mesh,'name')
            fprintf(fid_log,'Second Recon Basis              = yes, ''%s''\n', recon_mesh.name);
        else
            fprintf(fid_log,'Second Recon Basis              = yes\n');
        end
    else
        fprintf(fid_log,'Second Recon Basis              = no\n');
    end
    fprintf(fid_log,'maximum order of moments        = %f Hz\n',max_order);
    fprintf(fid_log,'Initial Regularization          = %d\n',lambda.value);
    fprintf(fid_log,'Filter                          = %d\n',filter_flag);
end

%% RECONSTRUCTION

% loop through experiments
for ind_experiment = 1:experiments_no
    % copy initial optical properties, used at the end to reset the mesh for the next experiment
    mua_copy = mesh.mua;
    mus_copy = mesh.mus;

    % converge for this experiment
    for ind_iter = 1:iter_max

        % Calculate jacobian and data for this update
%         [J,data_update] = calculate_jacobian(mesh,frequency,recon_mesh,isReconMesh,parameter);
% ------------ changed by jingjing jiang 2019.05.22 START --------------- %
       [J,data_update] = jacobian_stnd_TR_allBins(mesh, ...
           t_cfg.range, t_cfg.step,recon_mesh);

%        [J,data_update] = convert2gates_TR(J,data_update,...
%            t_cfg.range, t_cfg.step, t_cfg.num_gates); % get gated data and jacobian
        % reshape J into n_moments x nodes
% ------------ changed by jingjing jiang 2019.05.22 END ----------------- %
        % Interpolate onto recon mesh if needed
        if isReconMesh
            % use already calculated interpolation functions (mesh.fine2coarse)
            recon_mesh = interpolate_mesh2mesh(mesh,recon_mesh,mesh.fine2coarse);
        end

        % save the guess data for differential imaging
        if isDifferential
            if ind_iter == 1
%                 data_guess.amplitude = data_update.amplitude;
%                 data_guess.phase = data_update.phase;
                data_guess.tpsf = data_update.tpsf;
            end
        else
            % set as empty as not needed for absolute imaging
            data_guess = [];
        end

        % difference in measured and simulated data
        data_diff = data_difference_rytov(data_update, data_measured, ind_experiment, isDifferential, data_guess);

        % square error
%         pj_error(ind_iter,ind_experiment) = sum(abs(data_diff.^2));
        pj_error(ind_iter,ind_experiment) = nansum(abs(data_diff.^2));

        % display and log if needed
        disp('---------------------------------');
        disp(['Experiment/Iteration Number     = ' num2str(ind_experiment) '(of ' num2str(experiments_no) ')/' num2str(ind_iter)]);
        disp(['Projection error                = ' num2str(pj_error(ind_iter,ind_experiment))]);

        if ~isempty(log_filename)
            fprintf(fid_log,'---------------------------------\n');
            fprintf(fid_log,'Experiment/Iteration Number     = %d(of %d)/%d\n',ind_experiment,experiments_no,ind_iter);
            fprintf(fid_log,'Projection error                = %f\n',pj_error(ind_iter,ind_experiment));
        end

        % stopping criteria
        if ind_iter > 1
            % percentage of error change as compared to last step
            error_change = diff(pj_error(ind_iter-1:ind_iter,ind_experiment))/pj_error(ind_iter-1,ind_experiment) * 100;
            disp(['Projection error change         = ' num2str(error_change,'%4.2f') '%']);
            if ~isempty(log_filename)
                fprintf(fid_log,'Projection error change         = %f %%\n',error_change);
            end
            % check for stopping criteria
            if error_change >= -2
                disp('---------------------------------');
                disp('STOPPING CRITERIA REACHED       : error change >= -2%');
                if ~isempty(log_filename)
                    fprintf(fid_log,'---------------------------------\n');
                    fprintf(fid_log,'STOPPING CRITERIA REACHED       : error change >= -2%%\n');
                end
                break;
            end
        end

        % Normalize Jacobian by optical properties (get absorption and scatterin parts into similar amplitude level)
        J = normalize_jacobian(J,recon_mesh,parameter);

        % Set/update Thikonov regularization
        if ind_iter > 1
            reg_thikonov = reg_thikonov./10^0.25;
%             reg_thikonov = reg_thikonov./10^0.5;
%             reg_thikonov = reg_thikonov./10;

        else
            reg_thikonov = lambda.value;
        end

        % build Hessian
        if strcmp(lambda.type, 'JJt')
            Hess = (J*J');
        else
            Hess = (J'*J);
        end

        % calculate and apply the Thikonov regularization for amplitude and phase
%         Hess = regularize_thikonov(Hess, reg_thikonov, isCW);
        num_dt = t_cfg.range/ t_cfg.step;
        Hess = regularize_thikonov(Hess, reg_thikonov,num_dt);
        % Solve! Calculate the update.
        if strcmp(lambda.type, 'JJt')
            recon_update = J'*(Hess\data_diff);
        else
            recon_update = Hess\J'*data_diff;
        end
        
        %% 
        clear J data_update
        
        % apply the recovered update (normalize back the absorption and scatterin parts and add to the previous solution)
        recon_mesh = apply_update(recon_mesh, recon_update, parameter);

        % Interpolate back the recovered optical properties to the fine mesh, if needed
        if isReconMesh
            mesh = interpolate_mesh2mesh(recon_mesh,mesh,recon_mesh.coarse2fine,'mua','mus','kappa');
        else
            % simply apply new values if no second mesh basis is used
            mesh.mua = recon_mesh.mua;
            mesh.mus = recon_mesh.mus;
            mesh.kappa = recon_mesh.kappa;
        end

        % Filter if needed, see the help for details
        if filter_flag > 1
            mesh = mean_filter(mesh,abs(filter_flag));
        elseif filter_flag < 0
            mesh = median_filter(mesh,abs(filter_flag));
        end

    end

    % show message if maximum number of iterations reached
    if ind_iter >= iter_max
        disp('---------------------------------');
        disp('STOPPING CRITERIA REACHED       : max iterations');
        if ~isempty(log_filename)
            fprintf(fid_log,'---------------------------------\n');
            fprintf(fid_log,'STOPPING CRITERIA REACHED       : max iterations\n');
        end
    end

    % store this experiment recovery
    if strcmp(parameter,'all')
        mua_recon(:,ind_experiment) = mesh.mua;
        mus_recon(:,ind_experiment) = mesh.mus;
    elseif strcmp(parameter,'mus')
        mus_recon(:,ind_experiment) = mesh.mus;
    else
        mua_recon(:,ind_experiment) = mesh.mua;
    end

    % reset the mesh for the next experiment recovery
    mesh.mua = mua_copy;
    mesh.mus = mus_copy;
    mesh.kappa = 1./(3*(mesh.mua + mesh.mus));

end

% close the log file if needed
if ~isempty(log_filename)
    % close log file!
    fclose(fid_log);
end

% set the recovery for all experiments in the output mesh
if strcmp(parameter,'all')
    if isDifferential
        % substract the quess if differential
        mesh.mua = mua_recon - mua_copy;
        mesh.mus = mus_recon - mus_copy;
    else
        mesh.mua = mua_recon;
        mesh.mus = mus_recon;
    end
    mesh.kappa = 1./(3*(mesh.mua + mesh.mus));
elseif strcmp(parameter,'mus')
    if isDifferential
        % substract the quess if differential
        mesh.mus = mus_recon - mus_copy;
    else
        mesh.mus = mus_recon;
    end
    mesh.kappa = 1./(3*(mesh.mua + mesh.mus));
else
    if isDifferential
        % substract the quess if differential
        mesh.mua = mua_recon - mua_copy;
    else
        mesh.mua = mua_recon;
    end
    mesh.kappa = 1./(3*(mesh.mua + mesh.mus));
end

% if error returned
if nargout == 2
    % return only finished itrations, up to the max iterations
    pj_error = pj_error(isfinite(pj_error(:,1)),:);
end


end



% %% Calculate jacobians
% function [J,data] = calculate_jacobian(mesh,frequency,recon_mesh,isReconMesh,parameter)
% % See help for 'jacobian_FD'
%     if ~strcmp(parameter,'mua')
%         % always absorption and scattering sensitivity for CW and FD
%         if isReconMesh
%             [J,data] = jacobian_FD(mesh,frequency,recon_mesh,[],[],'all');
%         else
%             [J,data] = jacobian_FD(mesh,frequency,[],[],[],'all');
%         end
%     else
%         % absorption and scattering sensitivity for FD and absorption only for CW
%         if isReconMesh
%             [J,data] = jacobian_FD(mesh,frequency,recon_mesh);
%         else
%             [J,data] = jacobian_FD(mesh,frequency);
%         end
%     end
% 
%     % Set jacobian as Rytov approximation values (attenuation and phase parts)
%     J = J.complete;
% end

%% difference data - Rytov approach
function [data_diff] = data_difference_rytov(data_update, data_measured, ind_experiment, isDifferential, data_guess)
% data_update and data_measured should be vectors of the same size
% 'reference' is the simulated data for current solution, Rytov approach, thus log of amplitude
% 'anomaly' is the measured data reformated to match Jacobians formats, Rytov approach, thus log of amplitude

    % enabled source-detector pairs
%     mask_data = logical(data_measured.link(:,3))-1;
    mask_data = logical(data_measured.link(:,3));
% %     if isCW
%         % if Continuous Wave (zero frequency), only log of amplitude (attenuation denominator)
%         if isDifferential
%             % attenuation at this step
%             reference = log(data_guess.amplitude(mask_data,1)./data_update.amplitude(mask_data,1));
%             % measured attenuation
%             anomaly = data_measured.amplitude(mask_data,ind_experiment);
%             % attenuation change (measured - simulated)
%             data_diff = anomaly - reference;
%         else
%             % ln of update
%             reference = log(data_update.amplitude(mask_data,1));
%             % ln of reference
%             anomaly = log(data_measured.amplitude(mask_data,ind_experiment));
%             % -(attenuation), -ln(measured/current_solution)
%             data_diff = reference - anomaly;
%         end
%     else
%         % convert phase to radians
%         data_update.phase = data_update.phase/180*pi;
%         % make sure all phase shift is between 0 and 2*pi
%         data_update.phase(mask_data,data_update.phase(mask_data,:)<0) = data_update.phase(mask_data,data_update.phase(mask_data,:)<0) + 2*pi;
%         % if Frequency Domain, log of amplitude (attenuation denominator)
        % and phase in radians
        num_gates = size(data_measured.tpsf,2);
        if isDifferential
            error('Differential version has not been implemented for moments- based reconstruction, set isDifferential = 1')
%             data_guess.phase = data_guess.phase/180*pi;
%             data_guess.phase(mask_data,data_guess.phase(mask_data,:)<0) = data_guess.phase(mask_data,data_guess.phase(mask_data,:)<0) + 2*pi;
%             % attenuation + phase change at this step
%             reference = [log(data_guess.amplitude(mask_data,1)./data_update.amplitude(mask_data,1)); ...
%                 data_update.phase(mask_data,1)-data_guess.phase(mask_data,1)];
%             
%             %measured attenuation and phase change
%             anomaly = [data_measured.amplitude(mask_data,ind_experiment); data_measured.phase(mask_data,ind_experiment)];
%             
%             % split to match the Jacobian interlaced format
%             data_diff = zeros(size(reference));
%             % attenuation change
%             data_diff(1:2:end) = anomaly(1:end/2) - reference(1:end/2);
% %             data_diff(1:2:end) = reference(1:end/2) - anomaly(1:end/2);
%             % phase change difference
%             data_diff(2:2:end) = anomaly(end/2+1:end) - reference(end/2+1:end);
% %             data_diff(2:2:end) = reference(end/2+1:end) - anomaly(end/2+1:end);
        else
%             data_measured.phase = data_measured.phase/180*pi;
%             data_measured.phase(mask_data,data_measured.phase(mask_data,:)<0) = data_measured.phase(mask_data,data_measured.phase(mask_data,:)<0) + 2*pi;
%             reference = [log(data_update.amplitude(mask_data,1)); data_update.phase(mask_data,1)];
%             anomaly = [log(data_measured.amplitude(mask_data,ind_experiment)); data_measured.phase(mask_data,ind_experiment)];
%             
%             % split to match the Jacobian interlaced format
%             data_diff = zeros(size(reference));
%             % -(attenuation)
%             data_diff(1:2:end) = reference(1:end/2) - anomaly(1:end/2);
%             % -(phase difference)
%             data_diff(2:2:end) = reference(end/2+1:end) - anomaly(end/2+1:end);
%             reference = log( data_update.tpsf(mask_data, 1));
%             anomaly = log(data_measured.tpsf(mask_data, 1));
            reference =   log(data_update.tpsf(mask_data, 1)+...
                (data_update.tpsf(mask_data, 1)<=0));
            anomaly = log(data_measured.tpsf(mask_data, 1)+...
                (data_measured.tpsf(mask_data, 1)<=0));
            if num_gates > 1
                for ior = 2:num_gates
%                     reference = [reference;- data_update.tpsf(mask_data,ior)];
%                     anomaly = [anomaly;- data_measured.tpsf(mask_data, ior)];
                    reference = [reference; log(data_update.tpsf(mask_data, ior)+...
                (data_update.tpsf(mask_data, ior)<=0));];
                     anomaly = [anomaly; log(data_measured.tpsf(mask_data, ior)+...
                (data_measured.tpsf(mask_data, ior)<=0))];
                end
            end
            
% split to match the Jacobian interlaced format
            len = size(reference,1);
            data_diff = zeros(len,1);
            mm = len/num_gates;
            for ior = 1 : num_gates
                data_diff(ior:num_gates:len) = ...
                    reference(1+mm*(ior-1):mm*ior) ...
                    - anomaly(1+mm*(ior-1):mm*ior);
            end
        end
%     end
% remove nan values
% data_diff(isnan(data_diff))=0;

end


%% Normalize Jacobian by optical properties
function [J] = normalize_jacobian(J,recon_mesh,parameter)
    if strcmp(parameter,'mua')
        J_temp = J.mualnI;
%         [nn, pp, mm] = size(J_temp);
        [pp, nn, mm] = size(J_temp);
        J_norm = zeros(pp*mm,nn);
        for ii = 1:mm
%             J_m = squeeze(J_temp(:,:,ii)');
            J_m = squeeze(J_temp(:,:,ii));
            J_norm(ii:mm:pp*mm,:) = J_m;
        end
        J = J_norm .* recon_mesh.mua';
    elseif strcmp(parameter,'mus')
        error('has not implement normalize_jacobian for mus')
    elseif strcmp(parameter,'all')
        error('has not implement normalize_jacobian for all')
    end
%     if size(J,2) == nodes_size
%         % if only mua part of Jacobian present
%         J = J .* recon_mesh.mua';
%     elseif size(J,2) == 2*nodes_size
%         % if kappa and mua parts of Jacobian present
%         if strcmp(parameter,'all')
%             J = J .* [recon_mesh.kappa' recon_mesh.mua'];
%         elseif strcmp(parameter,'mus')
%             J = J(:,1:end/2) .* recon_mesh.kappa';
%         else
%             J = J(:,end/2+1:end) .* recon_mesh.mua';
%         end
%     else
%         error('Jacobian and mesh sizes do not match.')
%     end
end

%% regularize Thikonov
            
function [Hess] = regularize_thikonov(Hess, reg_thikonov,num_gates)
% calculate regularization for amplitude and phase
    reg = ones(size(Hess,1),1);
    nn =  num_gates;
    for ior = 1: nn
        value_reg = max(diag(Hess(ior:nn:end,ior:nn:end)));
        reg_temp = reg_thikonov*value_reg;
        reg(ior:nn:end)=reg(ior:nn:end).*reg_temp;
    end
%     if isCW
%         % scale by max of the Hessian diagonal
%         reg_amp = reg_thikonov*max(diag(Hess));
%         reg = reg.*reg_amp;
%     else
%         % scale by max of the amplitude-Hessian diagonal
%         reg_amp = reg_thikonov*max(diag(Hess(1:2:end,1:2:end)));
%         % scale by max of the phase-Hessian diagonal
%         reg_phs = reg_thikonov*max(diag(Hess(2:2:end,2:2:end)));
%         reg(1:2:end) = reg(1:2:end).*reg_amp;
%         reg(2:2:end) = reg(2:2:end).*reg_phs;
%     end
% 
%     % Add regularisation to Hessian diagonal
    for ind = 1:size(Hess,1)
        Hess(ind,ind) = Hess(ind,ind) + reg(ind);
    end
    
end

%% apply the recovered update
function recon_mesh = apply_update(recon_mesh, recon_update, parameter)
% simply add the update to existing values
    if strcmp(parameter,'all')
        % first, normalize back the update as the Jacobian was first
        % normalized by optical properties
        recon_update = recon_update.*[recon_mesh.kappa; recon_mesh.mua];
        recon_mesh.kappa = recon_mesh.kappa + (recon_update(1:end/2));
        recon_mesh.mua = recon_mesh.mua + (recon_update(end/2+1:end));
        recon_mesh.mus = (1./(3.*recon_mesh.kappa))-recon_mesh.mua;
    elseif strcmp(parameter,'mus')
        % first, normalize back the update
        recon_update = recon_update.*recon_mesh.kappa;
        recon_mesh.kappa = recon_mesh.kappa + recon_update;
        recon_mesh.mus = (1./(3.*recon_mesh.kappa))-recon_mesh.mua;
    else
        % first, normalize back the update
        recon_update = recon_update.*recon_mesh.mua;
        recon_mesh.mua = recon_mesh.mua + recon_update;
        recon_mesh.kappa = 1./(3*(recon_mesh.mua + recon_mesh.mus));
    end
end