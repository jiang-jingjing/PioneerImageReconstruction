%>@brief Brief description of the function
%>
%> Shannon Entropy (SE) 
%>
%>@param sig input signal (propability distribution)
%>
%> @retval se shannon entropy
%>
%> author: jingjing jiang jing.jing.jiang@outlook.com
%> created: 2021.08.11
function se = ShannonEntropy(sig)
% normalization of the input signal
% p = (sig - mean(sig))/std(sig);
% 
% p = sig ./ sum(sig(:));
% se = (-1) * ((p(p>0)).*(log(p(p>0))));