%>@brief Brief description of the function
%>
%> image reconstruction for fourier domain modality 
%>
%>@param volIn tissue volume
%>@param paras parameters  
%>
%> @retval volOut reconstructed volume
  
function [volOut, varargout] = reconstructionFD(volIn, paras, varargin)
en