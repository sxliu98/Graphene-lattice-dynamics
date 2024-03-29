%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Graphene_5NNFC_2D                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
a = 1.42e-10; % distance between neighboring atoms
m = 1.99e-26; % mass of carbon atom
h_bar = 1.0546e-34; % reduced Planck constant
kB = 1.38065e-23; % Boltzmann constant

pi = 3.141592653589793;
phi = diag([25.88, 8.42, 6.183, 4.037, -3.044, -0.492, -3.016, 3.948, 0.516, 0.564, 0.129, -0.521, 1.035, 0.166, 0.110])*16.0217662;

%% Coordinates of atoms
[A_I, B_I] = rotate(2/3*pi, [1; 0; 0]*a);
[A_II, B_II] = rotate(1/3*pi, [3/2; sqrt(3)/2; 0]*a);
[A_III, B_III] = rotate(2/3*pi, [1; sqrt(3); 0]*a);
[A_IV1, B_IV1] = rotate(2/3*pi, [2.5; sqrt(3)/2; 0]*a);
[A_IV2, B_IV2] = rotate(2/3*pi, [2.5; -sqrt(3)/2; 0]*a);
[A_V, B_V] = rotate(1/3*pi,[3; 0; 0]*a);
A_IV = [A_IV1, A_IV2];                  
B_IV = [B_IV1, B_IV2];

%% Interatomic force constants
[KAB_I, KBA_I] = rotate(2/3*pi, phi(1:3, 1:3));
[KAA_II, KBB_II] = rotate(1/3*pi, K(1/6*pi, phi(4:6, 4:6)));
[KAB_III, KBA_III] = rotate(2/3*pi, K(1/3*pi, phi(7:9, 7:9)));
[KAB_IV1, KBA_IV1] = rotate(2/3*pi, K(acos(2.5/sqrt(7)), phi(10:12, 10:12)));
[KAB_IV2, KBA_IV2] = rotate(2/3*pi, K(2*pi-acos(2.5/sqrt(7)), phi(10:12, 10:12)));
[KAA_V, KBB_V] = rotate(1/3*pi, phi(13:15, 13:15));
KAB_IV = cat(3,KAB_IV1, KAB_IV2);
KBA_IV = cat(3,KBA_IV1, KBA_IV2);

%% Three directions of k-space
n = 100;
k_full = zeros((30 + round(sqrt(3) * 10)) * n, 3); 
for i = 1:((30 + round(sqrt(3) * 10)) * n)
    if i < n * round(10 * sqrt(3)) 
        kx = (i-1) * 2 * pi / 3 / a / (n * round(10 * sqrt(3)));
        ky = 0;
    elseif i < (10 + round(10 * sqrt(3))) * n
        kx = 2 * pi / 3 / a;
        ky = (i - n * round(10 * sqrt(3))) / (10 * n - 1) * 2 * pi / 3 / a / sqrt(3);
    else
        kx = 2 * pi / 3 / a - (i - (10 + round(10 * sqrt(3))) * n) / (n * 20 - 1) * 2 * pi / 3 / a;
        ky = kx / sqrt(3);
    end
    k = [kx, ky, 0]; 
    k_full(i, :) = k;
end

%% Solving the characteristic equation
N_k = size(k_full,1);
omega = zeros(6,N_k);
for l = 1:N_k
    k = k_full(l,:);

    D = zeros(6, 6);
    DAA = zeros(3,3);
    DBB = zeros(3,3);
    DBA = zeros(3,3);
    DAB = zeros(3,3);

    for i = 1:3
        % dynamics matrix
        DAA = DAA + KAA_II(:,:,i) * exp(1i * k * (-A_II(:, i))) + KAA_II(:,:,i+3) * exp(1i * k * (-A_II(:, i+3))) + ...
               KAA_V(:,:,i) * exp(1i * k * (-A_V(:, i))) + KAA_V(:,:,i+3) * exp(1i * k * (-A_V(:, i+3)));
        DBB = DBB + KBB_II(:,:,i) * exp(1i * k * (-B_II(:, i))) + KBB_II(:,:,i+3) * exp(1i * k * (-B_II(:, i+3))) + ...
               KBB_V(:,:,i) * exp(1i * k * (-B_V(:, i))) + KBB_V(:,:,i+3) * exp(1i * k * (-B_V(:, i+3)));
        DAB = DAB + KAB_I(:,:,i) * exp(1i * k * (-A_I(:, i))) + KAB_III(:,:,i) * exp(1i * k * (-A_III(:, i))) +...
               KAB_IV(:,:,i) * exp(1i * k * (-A_IV(:, i))) + KAB_IV(:,:,i + 3) * exp(1i * k * (-A_IV(:, i + 3)));
        DBA = DBA + KBA_I(:,:,i) * exp(1i * k * (-B_I(:, i))) + KBA_III(:,:,i) * exp(1i * k * (-B_III(:, i))) + ...
               KBA_IV(:,:,i) * exp(1i * k * (-B_IV(:, i))) + KBA_IV(:,:,i + 3) * exp(1i * k * (-B_IV(:, i + 3)));
    end
    D(1:3, 4:6) = -DAB;
    D(4:6, 1:3) = -DBA;
    D(1:3, 1:3) = sum(KAB_I,3) + sum(KAA_II,3) + sum(KAB_III,3) + sum(KAB_IV1,3) + sum(KAB_IV2,3) + sum(KAA_V,3) - DAA;
    D(4:6, 4:6) = sum(KAB_I,3) + sum(KAA_II,3) + sum(KAB_III,3) + sum(KAB_IV1,3) + sum(KAB_IV2,3) + sum(KAA_V,3)- DBB;
    e = eig(D);
    w = sort(e);
    omega(:,l) = real(sqrt(w))/sqrt(m);
end

%% Correction of dispersion relationship
idx0 = (1:6)';
idx = (1:6)';
w_full_b = omega;
slope = zeros(6,2);
for i = 1:N_k
    if i>1 && i<N_k-1
        for j = 1:5
            if abs(w_full_b(idx0(j),i-1)-(w_full_b(idx0(j+1),i-1))) < 5e11
                if abs(slope(j,1)-slope(j,2))-abs(slope(j,1)-slope(j+1,2))>0
                    if(j==4 || j==1)
                        continue;
                    end
                    pos_1 = find(idx==j);
                    pos_2 = find(idx==(j+1));
                    tmp = idx(pos_1);
                    idx(pos_1) = idx(pos_2);
                    idx(pos_2) = tmp;
                end
            end
        end
    end
    if i>1
        slope(:,1) = w_full_b(idx0,i)-w_full_b(idx0,i-1);
    end
    if(i<N_k-1)
        slope(:,2) = w_full_b(idx0,i+1)-w_full_b(idx0,i);
    end
    omega(idx0,i) = w_full_b(idx,i);
end

%% Plotting dispersion relationship
newcolors = [234, 32, 39; 0, 98, 102; 27, 20, 100; 87, 88, 187; 111, 30, 81;
             238, 90, 36; 0, 148, 50; 6, 82, 221; 153, 128, 250; 131, 52, 113;
             247, 159, 31; 163, 203, 56; 18, 137, 167; 217, 128, 250; 181, 52, 113;
             255, 195, 18; 196, 229, 56; 18, 203, 196; 253, 167, 223; 237, 76, 103]./255; 
figure('OuterPosition',[0 0 600 450])
s = sqrt(3);

xk = [0, s, s + 1, s + 3];
kk = linspace(0, 4.7, N_k); 
colororder(newcolors)
hold on
for i = 1:6
    plot(kk,omega(i,:),'LineWidth',2);
end
set(gca,'linewidth',1.5,'FontSize',14);
xticks(xk);
xticklabels({'Γ', 'M', 'K', 'Γ'});
box on
xlim([0, 4.75]);
ylim([0, 3.1e14]);
xlabel('k', 'FontSize', 16)
ylabel('ω', 'FontSize', 16);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';

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