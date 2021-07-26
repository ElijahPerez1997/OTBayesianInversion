function s=MVM_Q(Qcell,v) 
%Inputs:
%Qcell = a cell array with d cells where the 1st cell is Q1 and the last
%cell is Qd
%We have Q=kron(Qd,...,Q1)
%v = an n x 1 column vector, where n=n1.n2...nd 

%Output:
%s = Q*v
              
%After each recurssion loop, we get rid of the data (Qd) we already used 

d=length(Qcell);
Qd=Qcell{d};
nd=length(Qd);    
n=length(v); 
nn=n/nd;

if d==1
    s=Qd*v;
else
   W=reshape(v,[nn,nd]);
   V=zeros(nn,nd);
   for i=1:nn
     ww=W(i,:)';
     V(i,:)=(Qd*ww)';     
   end
   
   S=zeros(nn,nd);
   for j=1:nd
       S(:,j)=MVM_Q(Qcell(1:end-1),V(:,j));
   end
   s=reshape(S,[n,1]);  
end
