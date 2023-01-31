function c0 = fitc0(oxd,cal)

c = ones(1,cal.nc); 
c = c./sum(max(0,c));

tol   = 1e-3;
maxit = 1e3;

c0   = c;
res  = 1;
res0 = 1;
it   = 1;
while res>tol && it<maxit

    c = max(0,c0 + max(0.01,c0)/10.*randn(size(c)).*res0^0.5);
    c = c./sum(c);

    res = norm((abs(c*cal.oxd - oxd)/oxd).^0.5)/sqrt(cal.nc);

    if res<1.01*res0
        c0 = c;
        res0 = res;
    end

    it = it+1;
end