%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        AGNR_Dispersion                         %%%
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

%% Plotting dispersion relationship
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
    
    %%% Atom 9
        D95 = KBB_II(:,:,4)*exp(1i * kx * 0);
        D96 = KBA_IV(:,:,2)*exp(1i * kx * (-1/2*a));
        D98 = KBA_III(:,:,2)*exp(1i * kx * a) + KBA_IV(:,:,5)*exp(1i * kx * (-2*a));
        D910 = KBA_I(:,:,2)*exp(1i * kx * (-1/2*a)) + KBA_IV(:,:,6)*exp(1i * kx * (5/2*a));
        D911 = KBB_II(:,:,3)*exp(1i * kx * (-3/2*a)) + KBB_II(:,:,5)*exp(1i * kx * (3/2*a));
        D912 = KBA_I(:,:,3)*exp(1i * kx * a) + KBA_III(:,:,1)*exp(1i * kx * (-2*a));  
        D913 = KBB_II(:,:,1)*exp(1i * kx * 0);
        D914 = KBA_I(:,:,1)*exp(1i * kx * (-1/2*a)) + KBA_IV(:,:,3)*exp(1i * kx * (5/2*a));  
        D915 = KBB_II(:,:,2)*exp(1i * kx * (-3/2*a)) + KBB_II(:,:,6)*exp(1i * kx * (3/2*a));
        D916 = KBA_III(:,:,3)*exp(1i * kx * (a)) + KBA_IV(:,:,1)*exp(1i * kx * (-2*a));
        D918 = KBA_IV(:,:,4)*exp(1i * kx * (-1/2*a));
    DA = [-D00  -D00 -D00 -D00 -D95 -D96 -D00 -D98  Dxx -D910 -D911 -D912 -D913 -D914 -D915 -D916 -D00 -D918 -D00 -D00];
    
    
    %%% Atom 10
        D101 = KAB_IV(:,:,4)*exp(1i * kx * (1/2*a));
        D105 = KAB_I(:,:,1)*exp(1i * kx * (1/2*a)) + KAB_IV(:,:,3)*exp(1i * kx * (-5/2*a));
        D106 = KAA_II(:,:,1)*exp(1i * kx * 0);
        D107 = KAB_III(:,:,3)*exp(1i * kx * (-a)) + KAB_IV(:,:,1)*exp(1i * kx * (2*a));
        D108 = KAA_II(:,:,2)*exp(1i * kx * (3/2*a)) + KAA_II(:,:,6)*exp(1i * kx * (-3/2*a));
        D109 = KAB_I(:,:,2)*exp(1i * kx * (1/2*a)) + KAB_IV(:,:,6)*exp(1i * kx * (-5/2*a));
        D1011 = KAB_I(:,:,3)*exp(1i * kx * (-a)) + KAB_III(:,:,1)*exp(1i * kx * (2*a));
        D1012 = KAA_II(:,:,3)*exp(1i * kx * (3/2*a)) + KAA_II(:,:,5)*exp(1i * kx * (-3/2*a));
        D1013 = KAB_IV(:,:,2)*exp(1i * kx * (1/2*a));
        D1014 = KAA_II(:,:,4)*exp(1i * kx * (0));
        D1015 = KAB_III(:,:,2)*exp(1i * kx * (-a)) + KAB_IV(:,:,5)*exp(1i * kx * (2*a));
    DB = [-D101 -D00 -D00 -D00 -D105 -D106 -D107 -D108 -D109 Dxx -D1011 -D1012 -D1013 -D1014 -D1015 -D00 -D00 -D00 -D00 -D00];
    
    %%% Atom 11
        D114 = KBA_IV(:,:,2)*exp(1i * kx * (-1/2*a));
        D115 = KBB_II(:,:,5)*exp(1i * kx * (3/2*a)) + KBB_II(:,:,3)*exp(1i * kx * (-3/2*a));
        D116 = KBA_III(:,:,2)*exp(1i * kx * a) + KBA_IV(:,:,5)*exp(1i * kx * (-2*a));
        D117 = KBB_II(:,:,4)*exp(1i * kx * 0);
        D118 = KBA_I(:,:,2)*exp(1i * kx * (-1/2*a)) + KBA_IV(:,:,6)*exp(1i * kx * (5/2*a));
        D119 = KBB_II(:,:,6)*exp(1i * kx * (3/2*a)) + KBB_II(:,:,2)*exp(1i * kx * (-3/2*a));
        D1110 = KBA_I(:,:,3)*exp(1i * kx * (a)) + KBA_III(:,:,1)*exp(1i * kx * (-2*a));
        D1112 = KBA_I(:,:,1)*exp(1i * kx * (-1/2*a)) + KBA_IV(:,:,3)*exp(1i * kx * (5/2*a));
        D1114 = KBA_III(:,:,3)*exp(1i * kx * (a)) + KBA_IV(:,:,1)*exp(1i * kx * (-2*a));
        D1115 = KBB_II(:,:,1)*exp(1i * kx * 0);
        D1116 = KBA_IV(:,:,4)*exp(1i * kx * (-1/2*a));
    DC = [-D00 -D00 -D00 -D114 -D115 -D116 -D117 -D118 -D119 -D1110 Dxx -D1112 -D00 -D1114 -D1115 -D1116 -D00 -D00 -D00 -D00];
    
    %%% Atom 12
        D125 = KAB_III(:,:,3)*exp(1i * kx * (-a)) + KAB_IV(:,:,1)*exp(1i * kx * (2*a));
        D127 = KAB_IV(:,:,4)*exp(1i * kx * (1/2*a));
        D128 = KAA_II(:,:,1)*exp(1i * kx * 0);
        D129 = KAB_I(:,:,3)*exp(1i * kx * (-a)) + KAB_III(:,:,1)*exp(1i * kx * (2*a));
        D1210 = KAA_II(:,:,2)*exp(1i * kx * (3/2*a)) + KAA_II(:,:,6)*exp(1i * kx * (-3/2*a));
        D1211 = KAB_I(:,:,1)*exp(1i * kx * (1/2*a)) + KAB_IV(:,:,3)*exp(1i * kx * (-5/2*a));
        D1213 = KAB_III(:,:,2)*exp(1i * kx * (-a)) + KAB_IV(:,:,5)*exp(1i * kx * (2*a));
        D1214 = KAA_II(:,:,3)*exp(1i * kx * (3/2*a)) + KAA_II(:,:,5)*exp(1i * kx * (-3/2*a));
        D1215 = KAB_I(:,:,2)*exp(1i * kx * (1/2*a)) + KAB_IV(:,:,6)*exp(1i * kx * (-5/2*a));
        D1216 = KAA_II(:,:,4)*exp(1i * kx * 0);
        D1219 = KAB_IV(:,:,2)*exp(1i * kx * (1/2*a));
    DD = [-D00 -D00 -D00 -D00 -D125 -D00 -D127 -D128 -D129 -D1210 -D1211 Dxx -D1213 -D1214 -D1215 -D1216 -D00 -D00 -D1219 -D00];
    
         D(1:3,[1:36 12*N-23:12*N]) = [DA(:,25:60) DA(:,1:24)];
         D(4:6,[1:36 12*N-23:12*N]) = [DB(:,25:60) DB(:,1:24)];
         D(7:9,[1:36 12*N-23:12*N]) = [DC(:,25:60) DC(:,1:24)];
         D(10:12,[1:36 12*N-23:12*N]) = [DD(:,25:60) DD(:,1:24)];
    
         D(13:15,[1:48 12*N-11:12*N]) = [DA(:,13:60) DA(:,1:12)];
         D(16:18,[1:48 12*N-11:12*N]) = [DB(:,13:60) DB(:,1:12)];
         D(19:21,[1:48 12*N-11:12*N]) = [DC(:,13:60) DC(:,1:12)];
         D(22:24,[1:48 12*N-11:12*N]) = [DD(:,13:60) DD(:,1:12)];
    
        if N >= 5
        for j = 0:N-5
         D(12*j+25:12*j+27,12*j+1:12*j+60) = DA;
         D(12*j+28:12*j+30,12*j+1:12*j+60) = DB;
         D(12*j+31:12*j+33,12*j+1:12*j+60) = DC;
         D(12*j+34:12*j+36,12*j+1:12*j+60) = DD;
        end
        end
    
         D(12*N-23:12*N-21,[1:12 12*N-47:12*N]) = [DA(:,49:60) DA(:,1:48)];
         D(12*N-20:12*N-18,[1:12 12*N-47:12*N]) = [DB(:,49:60) DB(:,1:48)];
         D(12*N-17:12*N-15,[1:12 12*N-47:12*N]) = [DC(:,49:60) DC(:,1:48)];
         D(12*N-14:12*N-12,[1:12 12*N-47:12*N]) = [DD(:,49:60) DD(:,1:48)];
    
         D(12*N-11:12*N-9,[1:24 12*N-35:12*N]) = [DA(:,37:60) DA(:,1:36)];
         D(12*N-8:12*N-6,[1:24 12*N-35:12*N]) = [DB(:,37:60) DB(:,1:36)];
         D(12*N-5:12*N-3,[1:24 12*N-35:12*N]) = [DC(:,37:60) DC(:,1:36)];
         D(12*N-2:12*N,[1:24 12*N-35:12*N]) = [DD(:,37:60) DD(:,1:36)];
    
         if N == 4
         D(1:3,:) = [DA(:,25:48) DA(:,1:12)+DA(:,49:60) DA(:,13:24)];
         D(4:6,:) = [DB(:,25:48) DB(:,1:12)+DB(:,49:60) DB(:,13:24)];
         D(7:9,:) = [DC(:,25:48) DC(:,1:12)+DC(:,49:60) DC(:,13:24)];
         D(10:12,:) = [DD(:,25:48) DD(:,1:12)+DD(:,49:60) DD(:,13:24)];
    
         D(13:15,:) = [DA(:,13:48) DA(:,1:12)+DA(:,49:60)];
         D(16:18,:) = [DB(:,13:48) DB(:,1:12)+DB(:,49:60)];
         D(19:21,:) = [DC(:,13:48) DC(:,1:12)+DC(:,49:60)];
         D(22:24,:) = [DD(:,13:48) DD(:,1:12)+DD(:,49:60)];
    
         D(25:27,:) = [DA(:,1:12)+DA(:,49:60) DA(:,13:48)];
         D(28:30,:) = [DB(:,1:12)+DB(:,49:60) DB(:,13:48)];
         D(31:33,:) = [DC(:,1:12)+DC(:,49:60) DC(:,13:48)];
         D(34:36,:) = [DD(:,1:12)+DD(:,49:60) DD(:,13:48)];
    
         D(37:39,:) = [DA(:,37:48) DA(:,1:12)+DA(:,49:60) DA(:,13:36)];
         D(40:42,:) = [DB(:,37:48) DB(:,1:12)+DB(:,49:60) DB(:,13:36)];
         D(43:45,:) = [DC(:,37:48) DC(:,1:12)+DC(:,49:60) DC(:,13:36)];
         D(46:48,:) = [DD(:,37:48) DD(:,1:12)+DD(:,49:60) DD(:,13:36)];
         end
    
        e = eig(D);
        w = sort(real(sqrt(e)));
        omega(:,i) = w/sqrt(m);
    
    end

    subplot(1,3,n)
    plot(k_full/(pi/3/a),omega/1e14,'b-','LineWidth',1.5)
    set(gca,'linewidth',1.5,'FontSize',14);
    xlabel('$k_x/k_{max}$','Interpreter','latex','FontSize',20,'FontWeight','bold')
    ylabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold')
    xlim([0 1])
    ylim([0 3.2])
    title([num2str(N*2),'-AGNR'])

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