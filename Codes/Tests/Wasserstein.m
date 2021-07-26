function W2=Wasserstein(f,g,t)

% f,g = two nx1 column vectors in the probability simplex 
%       (nonnegative & summing to one).
% f, g are normalized signals, not the original signals, to ensure
% positivity and unit mass

%t (or x) is the independent variable for f and g

%Written by Mohammad Motamed at UNM in 2017
%=======================================================
dt=t(2)-t(1);
n=length(t);

%Step 1. Compute CDF inverse of f and g

Fint=zeros(size(f));
Gint=Fint;
for k=1:n
    Fint(k)=sum(f(1:k));
    Gint(k)=sum(g(1:k));
end

%Step 2. Compute the map tG = G^{-1}oF by a search algorithm (here: built-in find) and interpolation

tG=zeros(size(t)); %pre-allocation of memory
tG(1)=t(1);
for k=2:n-1
    val=Fint(k);
    
    kL=find(Gint<=val); 
    kR=find(Gint>val);  
    
    if isempty(kL)
        tG(k)=0;
    elseif isempty(kR)
        tG(k)=t(end);
    else
        kL=kL(end); 
        tL=t(kL); 
        GL=Gint(kL);
        kR=kR(1); tR=t(kR); GR=Gint(kR);
        tG(k)=tL+(tR-tL)*(val-GL)/(GR-GL);
        %plot([tG(k) tG(k)],[0 1],'m')
        %plot([t(k) tG(k)],[val val],'k')
    end
    
end
tG(n)=t(n);

%Step 3. Compute the W2 distance 

dt=1;
W2=dt*sum(((t-tG).^2).*f); %This is in fact the square of W2 
                           %(should take its square root if needed)


