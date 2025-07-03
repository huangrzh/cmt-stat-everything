function [H] = GetHam_1d_2site4(nsite,miu)
H = zeros(2*nsite, 2*nsite);
for i1 = 1:nsite
    j1 = i1 + 1;
    if j1>nsite
        j1 = 1;
    end
    
    mm = i1-1;
    if mm<1
        mm = mm + nsite;
    end
    
    ia = 2*i1-1;
    ib = 2*i1-0;
    jam = 2*j1-1;
    jbp = 2*j1-0;
    jbm = 2*mm-0;
    
    H(ib,jbp) = 0.5; % (2,2)
    
    H(ia,ib) = 1; % (1,2)
    H(ia,jbp) = -0.5;
    H(ia,jbm) = -0.5;
end
H = H + H';

for i1 = 1:nsite
    ia = 2*i1-1;
    H(ia,ia) = 1; % (1,1) 
end
H = H - miu*eye(2*nsite,2*nsite);
end