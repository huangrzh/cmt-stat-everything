function [H] = GetHam_1d(t1, t2, miu, nsite)
H = zeros(nsite, nsite);
for i1 = 1:nsite
    j1 = i1 + 1;
    if j1>nsite
        j1 = 1;
    end
    j2 = j1 + 1;
    if j2>nsite
        j2 = 1;
    end
    H(i1,j1) = -t1;
    H(j1,i1) = -t1;
    H(i1,j2) = -t2;
    H(j2,i1) = -t2;
    H(i1,i1) = -miu;
end
end