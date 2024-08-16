% Particle-hole pair excitations in a single-band free fermion
% 25/02/2022, HRZ


close all;
N = 100;
k = (1:N)/N;
Ek = -cos(k*2*pi);
figure; hold on;

for k1 = 1:N
    for k2 = 1:N
        if and(Ek(k1)<=0, Ek(k2)>=0)
            q = k(k2) - k(k1);
            if q<0
                q = q + 1;
            end
            plot(q, Ek(k2)-Ek(k1), '*');
        end
    end
end

xlabel('q/2\pi')
ylabel('E(q)')
title('Particle-hole pair excitation')