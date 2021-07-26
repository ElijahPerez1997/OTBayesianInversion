function dS=Sinkhorn2D(f,g,Qcell,lambda)
% Compute debiased-Sinkhorn distance
%Written by M. Motamed at UNM, 2020

%---------------------------
% Required Inputs:
%---------------------------
% f,g = two nx1 column vectors in the probability simplex 
%       (nonnegative & summing to one). 
%
% Q= nxn matrix, equal to exp(-lambda*C), where C is nxn distance matrix 
%    (C has zeros on the diagonal and such that c_ij < c_ik + c_kj)
%    It is given by Qcell and lambda
%---------------------------
tolerance=1e-2; %tol in stopping criterion 
maxIter=500; % max number of Sinkhorn fixed point iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compt=0;
n=length(f);
u=ones(n,1)/n;

s=MVM_Q(Qcell,u); %s=Q*u;
v=g./s;
Qv=MVM_Q(Qcell,v); %s=Q*v;
while compt<maxIter
    
    u=f./Qv;
    Qu=MVM_Q(Qcell,u); %s=Q*u;
    v=g./Qu;
    Qv=MVM_Q(Qcell,v); %s=Q*v;
    compt=compt+1;
    
    % check the stopping criterion every 10 iterations
    if mod(compt,10)==1 || compt==maxIter
        Criterion=norm(u.*Qv-f,1);
        
        if Criterion<tolerance || isnan(Criterion) 
           break;
        end
    end
end

%new way
dS=(f'* log(u) + g'* log(v))/lambda - 1/lambda;

%Old way
%dS=sum(u.*(QC*v));
%dS= -sum(u.*((Q.*log(Q))*v))/lambda; %this is the same as dS above

