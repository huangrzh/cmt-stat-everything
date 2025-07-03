function [H] = GetHam_1d_2site3(nsite)
H = zeros(2*nsite, 2*nsite);
for i1 = 1:nsite
    j1 = i1 + 1;
    if j1>nsite
        j1 = 1;
    end
    pp = j1 + 1;
    if pp>nsite
        pp = 1;
    end
    mm = i1-2;
    if mm<1
        mm = mm + nsite;
    end
    
    ia = 2*i1-1;
    ib = 2*i1-0;
    ja = 2*j1-1;
    jb = 2*j1-0;
    ka = 2*pp-1;
    kp = 2*pp-0;    
    km = 2*mm-0;
    
    %H(ia,ja) = 0; % (1,1)
    
    H(ib,jb) = 0.5; % (2,2)
    
    H(ia,ib) = 1; % (1,2)
    H(ia,kp) = -0.5;
    H(ia,km) = -0.5;
end
H = H + H';

for i1 = 1:nsite
    ia = 2*i1-1;
    H(ia,ia) = 1; % (1,1) 
end
H = H - eye(2*nsite,2*nsite);
end