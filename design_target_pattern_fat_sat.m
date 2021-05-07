function [d,f,minf,maxf,nobj_f] = design_target_pattern_fat_sat(nom_fa,b0map,roi,b0fieldStrength)
%design_target_pattern_fat_sat Designs the target pattern for optimization
%of a fat saturation pulse taking into account the frequency inhomogeneity
%offsets.
%   
% INPUTS:
% (required)
% nom_fa           [1]          - nominal flip angle in degrees
% b0map            [nx ny nz]   - B0 fieldmap in Hz
% roi              [nx ny nz]   - logical mask of the object
% b0fieldStrength  [1]          - B0 field strength
%
% OUTPUTS:
% d                [nf nx ny]   - target pattern
% f                [nf]         - vector of the frequencies included in the
%                               design patten (in Hz)
% minf             [1]          - minimum frequency included in the design
%                               (in Hz)
% maxf             [1]          - maximum frequency included in the design
%                               (in Hz)
% nobj_f           [1]          - number of sampling frequencies included
%                               in the design (length(f))

excBW=50;                                           % excitation (for fat) and non-excitation (for water) BW of 50Hz
f_fat=-3.3*42.58*b0fieldStrength;                   % fat frequency
minf=(min(b0map(roi))+f_fat-excBW/2)*1.05;          % include shifted fat frequency and take margins in addition to the excitation BW
maxf=(max([max(b0map(roi)) 0])+excBW/2)*1.05;       % include shifted water frequency and take margins in addition to the excitation BW
df=10;                                              % sampling rate of 10Hz
f=minf:df:maxf;                                     % frequency axis
nobj_f=length(f);

% start implementing the target pattern
d=zeros([nobj_f, size(b0map)]);
for x=1:size(b0map,1)
    for y=1:size(b0map,2)
        if roi(x,y)==1
            fatVoxFreq=b0map(x,y)+f_fat;                  % fat frequency of the voxel given the B0 offset
            [~, fatVoxFreqIndex]=min(abs(f-fatVoxFreq));  % closest position of this frequency on the target design frequency grid (first dimension)
            d(fatVoxFreqIndex,x,y)=sind(nom_fa);          % apply pulse at this frequency with nominal flip angle
            half_excBW_index_span=round(excBW/2/df);      % half the span of the excitation BW on the frequency grid
            d((fatVoxFreqIndex-half_excBW_index_span):(fatVoxFreqIndex+half_excBW_index_span),x,y)=sind(nom_fa);  % apply flip angle also around the fat frequency with excitation BW
        end
    end
end

end

