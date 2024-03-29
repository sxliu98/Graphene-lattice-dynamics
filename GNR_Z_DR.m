%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        ZGNR_Dispersion                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
newcolors = [234, 32, 39; 0, 98, 102; 27, 20, 100; 87, 88, 187; 111, 30, 81;
             238, 90, 36; 0, 148, 50; 6, 82, 221; 153, 128, 250; 131, 52, 113;
             247, 159, 31; 163, 203, 56; 18, 137, 167; 217, 128, 250; 181, 52, 113;
             255, 195, 18; 196, 229, 56; 18, 203, 196; 253, 167, 223; 237, 76, 103]./255; 

Nx = [4 8 16];
a = 1.42e-10; % distance between neighboring atoms
phi = diag([36.50, 24.50, 9.82, 8.80, -3.23, -0.40, 3.00, -5.25, 0.15, -1.92, 2.29, -0.58])*10;
m = 1.99e-26; % mass of carbon atom

%% Interatomic force constants
[KAB_I, KBA_I] = rotate(2/3*pi, phi(1:3, 1:3));
[KAA_II, KBB_II] = rotate(1/3*pi, K(1/6*pi, phi(4:6, 4:6)));
[KAB_III, KBA_III] = rotate(2/3*pi, K(1/3*pi, phi(7:9, 7:9)));
[KAB_IV1, KBA_IV1] = rotate(2/3*pi, K(acos(2.5/sqrt(7)), phi(10:12, 10:12)));
[KAB_IV2, KBA_IV2] = rotate(2/3*pi, K(2*pi-acos(2.5/sqrt(7)), phi(10:12, 10:12)));
KAB_IV = cat(3,KAB_IV1, KAB_IV2);
KBA_IV = cat(3,KBA_IV1, KBA_IV2);

%% Interatomic force constants
Nk = 1000;
k_full = linspace(0,2*pi/3/a,Nk);
figure('OuterPosition',[0 0 1200 600])
for n = 1:3
    N = Nx(n);
    len = 12*N;
    omega = zeros(len,Nk);

for i = 1:length(k_full)
    kx = k_full(i);
    D = zeros(len,len); 
    D00 = zeros(3,3);
    Dxx = sum(KBA_I,3) + sum(KBB_II,3) + sum(KBA_III,3) + sum(KBA_IV,3);

%%% D5
    D52 = KBA_IV(:,:,3)*exp(1i * kx * (sqrt(3)/2*a)) + KBA_IV(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D53 = KBB_II(:,:,5)*exp(1i * kx * (-sqrt(3)/2*a)) + KBB_II(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));
    D54 = KBA_I(:,:,3)*exp(1i * kx * 0) + KBA_III(:,:,2)*exp(1i * kx * (-sqrt(3)*a)) + KBA_III(:,:,3)*exp(1i * kx * (sqrt(3)*a));
    D55 = KBB_II(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBB_II(:,:,4)*exp(1i * kx * (-sqrt(3)*a));
    D56 = KBA_I(:,:,1)*exp(1i * kx * (sqrt(3)/2*a)) + KBA_I(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KBA_IV(:,:,2)*exp(1i * kx * (-3*sqrt(3)/2*a)) + KBA_IV(:,:,4)*exp(1i * kx * (3*sqrt(3)/2*a));
    D57 = KBB_II(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KBB_II(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a));
    D58 = KBA_III(:,:,1)*exp(1i * kx * 0) + KBA_IV(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBA_IV(:,:,5)*exp(1i * kx * (-sqrt(3)*a));

%%% D6
    D63 = KAB_III(:,:,1)*exp(1i * kx * 0) + KAB_IV(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAB_IV(:,:,5)*exp(1i * kx * (sqrt(3)*a));
    D64 = KAA_II(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KAA_II(:,:,3)*exp(1i * kx * (sqrt(3)/2*a));
    D65 = KAB_I(:,:,1)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB_I(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KAB_IV(:,:,2)*exp(1i * kx * (3*sqrt(3)/2*a)) + KAB_IV(:,:,4)*exp(1i * kx * (-3*sqrt(3)/2*a));
    D66 = KAA_II(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAA_II(:,:,4)*exp(1i * kx * sqrt(3)*a);
    D67 = KAB_I(:,:,3)*exp(1i * kx * 0) + KAB_III(:,:,2)*exp(1i * kx * (sqrt(3)*a)) + KAB_III(:,:,3)*exp(1i * kx * (-sqrt(3)*a));
    D68 = KAA_II(:,:,5)*exp(1i * kx * (sqrt(3)/2*a)) + KAA_II(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D69 = KAB_IV(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB_IV(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));

%%% D7
    D74 = KBA_IV(:,:,3)*exp(1i * kx * (sqrt(3)/2*a)) + KBA_IV(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D75 = KBB_II(:,:,5)*exp(1i * kx * (-sqrt(3)/2*a)) + KBB_II(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));
    D76 = KBA_I(:,:,3)*exp(1i * kx * 0) + KBA_III(:,:,2)*exp(1i * kx * (-sqrt(3)*a)) + KBA_III(:,:,3)*exp(1i * kx * (sqrt(3)*a));
    D77 = KBB_II(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBB_II(:,:,4)*exp(1i * kx * (-sqrt(3)*a));
    D78 = KBA_I(:,:,1)*exp(1i * kx * (sqrt(3)/2*a)) + KBA_I(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KBA_IV(:,:,2)*exp(1i * kx * (-3*sqrt(3)/2*a)) + KBA_IV(:,:,4)*exp(1i * kx * (3*sqrt(3)/2*a));
    D79 = KBB_II(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KBB_II(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a));
    D710 = KBA_III(:,:,1)*exp(1i * kx * 0) + KBA_IV(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBA_IV(:,:,5)*exp(1i * kx * (-sqrt(3)*a));

%%% D8
    D85 = KAB_III(:,:,1)*exp(1i * kx * 0) + KAB_IV(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAB_IV(:,:,5)*exp(1i * kx * (sqrt(3)*a));
    D86 = KAA_II(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KAA_II(:,:,3)*exp(1i * kx * (sqrt(3)/2*a));
    D87 = KAB_I(:,:,1)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB_I(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KAB_IV(:,:,2)*exp(1i * kx * (3*sqrt(3)/2*a)) + KAB_IV(:,:,4)*exp(1i * kx * (-3*sqrt(3)/2*a));
    D88 = KAA_II(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAA_II(:,:,4)*exp(1i * kx * (sqrt(3)*a));
    D89 = KAB_I(:,:,3)*exp(1i * kx * 0) + KAB_III(:,:,3)*exp(1i * kx * (-sqrt(3)*a)) + KAB_III(:,:,2)*exp(1i * kx * (sqrt(3)*a));
    D810 = KAA_II(:,:,5)*exp(1i * kx * (sqrt(3)/2*a)) + KAA_II(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D811 = KAB_IV(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB_IV(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));

   
     D(1:3,[1:12 12*N-8:12*N]) = [Dxx-D55  -D56 -D57 -D58 -D52 -D53 -D54];
     D(4:6,[1:15 12*N-5:12*N]) = [-D65 Dxx-D66 -D67  -D68 -D69 -D63 -D64];
     D(7:9,[1:18 12*N-2:12*N]) = [-D75 -D76 Dxx-D77 -D78 -D79  -D710 -D74];
     D(10:12,1:21) = [-D85 -D86 -D87 Dxx-D88 -D89 -D810 -D811];
    
    for j = 0:N-3
     D(12*j+13:12*j+15,12*j+4:12*j+24) = [-D52 -D53 -D54 Dxx-D55  -D56 -D57 -D58];
     D(12*j+16:12*j+18,12*j+7:12*j+27) = [-D63 -D64 -D65 Dxx-D66 -D67  -D68 -D69];
     D(12*j+19:12*j+21,12*j+10:12*j+30) = [-D74 -D75 -D76 Dxx-D77 -D78 -D79  -D710];
     D(12*j+22:12*j+24,12*j+13:12*j+33) = [-D85 -D86 -D87 Dxx-D88 -D89 -D810 -D811];
    end

     D(12*N-11:12*N-9,12*N-20:12*N) = [-D52 -D53 -D54 Dxx-D55 -D56 -D57 -D58];
     D(12*N-8:12*N-6,[1:3 12*N-17:12*N]) = [-D69 -D63 -D64 -D65 Dxx-D66 -D67  -D68];
     D(12*N-5:12*N-3,[1:6 12*N-14:12*N]) = [-D79  -D710 -D74 -D75 -D76 Dxx-D77 -D78];     
     D(12*N-2:12*N,[1:9 12*N-11:12*N]) = [-D89 -D810 -D811 -D85 -D86 -D87 Dxx-D88];

    e = eig(D);
    w = sort(real(sqrt(e)));
    omega(:,i) = w/sqrt(m);
end

subplot(1,3,n)
plot(k_full/(pi/sqrt(3)/a),omega/1e14,'b-','LineWidth',1.5)
set(gca,'linewidth',1.5,'FontSize',14);
xlabel('$k_x/k_{max}$','Interpreter','latex','FontSize',20,'FontWeight','bold')
ylabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold')
xlim([0 1])
ylim([0 3.2])
title([num2str(N*4),'-ZGNR'])

end
%% Rotation function
function [R, RB] = rotate(theta, r)
    Um = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    inv_Um = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    U = [cos(pi), sin(pi), 0; -sin(pi), cos(pi), 0; 0, 0, 1];
    inv_U = [cos(pi), -sin(pi), 0; sin(pi), cos(pi), 0; 0, 0, 1];
    n = 2*pi/theta;
    if numel(r) == 3
        R = zeros(3,n);      
        RB = zeros(3,n);
    else
        R = zeros(3,3,n);      
        RB = zeros(3,3,n);
    end
    for i = 1:n
        if numel(r) == 3
            r = inv_Um * r;
            rb = r*(-1);
            R(:, i) = r;
            RB(:, i) = rb;
        else
            r = inv_Um * r * Um;
            rb = inv_U * r * U;
            R(:, :, i) = r;
            RB(:, :, i) = rb;
        end

    end
end

function k = K(theta, k)
    U = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    inv_U = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    k = inv_U * k * U;
end