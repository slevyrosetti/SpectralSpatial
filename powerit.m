function [vec,value]=powerit(start,A,toler)
% function [vec,value]=powerit(start,A,toler)
% Power method for computing eigenvalues
% Susan Holmes, Stanford, 2002
%
% INPUTS:
% start            [npts]       - initialization vector of Eigenvalues
% A                [npts npts]  - either A'A or A'WA, where A is system matrix
% toler            [1]          - stopping tolerance, normally 1e-5
%
% OUTPUTS:
% vec              [npts]       - Eigenvector associated with largest Eigenvalue
% value            [1]          - largest Eigenvalue

    dd=1;
    x=start;
    n=10;
    iter=1;
    while dd> toler
        y=A*x;
        dd=abs(norm(x)-n);
        n=norm(x);
        x=y/n;
        %pause
        if mod(iter,5)==0
            fprintf('no. iter.: %3.0f\n', iter); 
        else
        end
        iter=iter+1;
    end
    vec=x;
    value=n;
    fprintf('total iter.: %3.0f\n', iter);
end
