function c0 = fitc0(oxd,caloxd)

c = ones(1,size(caloxd,1)); 
c = c./sum(max(0,c));

tol   = 1e-2;
maxit = 1e+4;

c0   = c;
res  = 1;
res0 = 1;
it   = 1;
while res>tol && it<maxit

    c = max(0.001,c0 + max(0.05,c0)/2.*randn(size(c)).*res0^0.25);
    c = c./sum(c);

    res = norm(abs(c*caloxd - oxd)./oxd.^0.5)/sqrt(length(oxd));

    if it==1; res0 = res; end

    if res<1.01*res0
        c0   = c;
        res0 = res;
    end

    it = it+1;
end