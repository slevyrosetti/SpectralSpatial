function f = freq_vector(minf,maxf,samp_f)
% function f = freq_vector(minf,maxf,samp_f)
% Builds frequency vector at particular sampling rate centered at 0 Hz
% Sydney Williams, University of Michigan 2016
%
% INPUTS:
% minf      [1]         - lowest frequency in vector in Hz
% maxf      [1]         - highest frequency in vector in Hz
% samp_f    [1]         - frequency sampling rate in Hz

% OUTPUTS:
% f         [Nf]        - frequency samples vector

nf=floor((maxf-minf)/samp_f)+1;             % number of frequency samples in simulation
f=[-nf/2:(nf/2-1)]*samp_f;                  % vector of frequencies
shift=round(minf-f(1));                     % shifts and centers simulation frequencies
f=f+shift;
if isempty(find(f==0))
    posneg=find(f(1:end-1).*f(2:end)<=0);	% finds where sign change is for zero insertion
    f=[f(1:posneg) 0 f(posneg+1:end)];      % inserts zero
end

end