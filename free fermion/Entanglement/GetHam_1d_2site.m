function [H] = GetHam_1d_2site(nsite, t1, ta, tb, tab, ua, ub)
H = zeros(2*nsite, 2*nsite);
for i1 = 1:nsite
    j1 = i1 + 1;
    if j1>nsite
        j1 = 1;
    end
    ia = 2*i1-1;
    ib = 2*i1-0;
    ja = 2*j1-1;
    jb = 2*j1-0;
    
    H(ia,ib) = t1;
    
    H(ia,ja) = ta;
    H(ib,jb) = tb;
    
    H(ib,ja) = tab;
end
H = H + H';

for i1 = 1:nsite
    H(2*i1-1,2*i1-1) = ua;
    H(2*i1-0,2*i1-0) = ub;
end
end