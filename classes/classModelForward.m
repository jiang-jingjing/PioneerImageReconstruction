%> @file classForwardModel.m
%> @brief class description
% ======================================================================
%> @brief Here we have a brief description of the class.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
classdef ModelForward
    properties 
        Type
    end
    properties 
        Path
    end
    methods
        function model = getType(str)
            if strcmp(str, 'MC')
                model.Type = str;
                Path = []; %MCXLAB
                addpath(genpath(Path))
            elseif strcmp(str, 'FEM')
                model.Type = str;
                Path = []; % Nirfaster
                addpath(genpath(Path))
            else
                error('forward model is MC or FEM')
            end
        end
    end
end
   