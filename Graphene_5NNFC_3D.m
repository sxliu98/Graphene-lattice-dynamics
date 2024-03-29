%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Graphene_5NNFC_3D                         %%%
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

%% k-space grid settings
kpoints = 200;
kx =  linspace(-2*pi/(a*3),2*pi/(a*3),kpoints);
ky = linspace(-pi/(a*sqrt(3)),pi/(a*sqrt(3)),kpoints);
[X,Y] = meshgrid(kx,ky);
for i = 1:kpoints
    for l = 1:kpoints
        if Y(i,l) < sqrt(3)*X(i,l)-2*sqrt(3)/3*pi/a || Y(i,l) < -sqrt(3)*X(i,l)-2*sqrt(3)/3*pi/a ...
                ||Y(i,l) > -sqrt(3)*X(i,l)+2*sqrt(3)/3*pi/a || Y(i,l) > sqrt(3)*X(i,l)+2*sqrt(3)/3*pi/a
            X(i,l) = NaN;
            Y(i,l) = NaN;
        end
    end
end
k_full = zeros(3,kpoints,kpoints);
for i = 1:kpoints
    for l = 1:kpoints
        k_full(:,i,l)  = [X(i,l), Y(i,l), 0];
    end
end

%% %% Solving the characteristic equation
k_num = kpoints*kpoints;
omega = zeros(kpoints,kpoints,6);
for i = 1:kpoints
    for j = 1:kpoints
        k = k_full(:,i,j)';
        D = zeros(6, 6);
        DAA = zeros(3,3);
        DBB = zeros(3,3);
        DBA = zeros(3,3);
        DAB = zeros(3,3);
    
        for l = 1:3
            % dynamics matrix
            DAA = DAA + KAA_II(:,:,l) * exp(1i * k * (-A_II(:, l))) + KAA_II(:,:,l+3) * exp(1i * k * (-A_II(:, l+3))) + ...
                   KAA_V(:,:,l) * exp(1i * k * (-A_V(:, l))) + KAA_V(:,:,l+3) * exp(1i * k * (-A_V(:, l+3)));
            DBB = DBB + KBB_II(:,:,l) * exp(1i * k * (-B_II(:, l))) + KBB_II(:,:,l+3) * exp(1i * k * (-B_II(:, l+3))) + ...
                   KBB_V(:,:,l) * exp(1i * k * (-B_V(:, l))) + KBB_V(:,:,l+3) * exp(1i * k * (-B_V(:, l+3)));
            DAB = DAB + KAB_I(:,:,l) * exp(1i * k * (-A_I(:, l))) + KAB_III(:,:,l) * exp(1i * k * (-A_III(:, l))) +...
                   KAB_IV(:,:,l) * exp(1i * k * (-A_IV(:, l))) + KAB_IV(:,:,l + 3) * exp(1i * k * (-A_IV(:, l + 3)));
            DBA = DBA + KBA_I(:,:,l) * exp(1i * k * (-B_I(:, l))) + KBA_III(:,:,l) * exp(1i * k * (-B_III(:, l))) + ...
                   KBA_IV(:,:,l) * exp(1i * k * (-B_IV(:, l))) + KBA_IV(:,:,l + 3) * exp(1i * k * (-B_IV(:, l + 3)));
        end
        D(1:3, 4:6) = -DAB;
        D(4:6, 1:3) = -DBA;
        D(1:3, 1:3) = sum(KAB_I,3) + sum(KAA_II,3) + sum(KAB_III,3) + sum(KAB_IV1,3) + sum(KAB_IV2,3) + sum(KAA_V,3) - DAA;
        D(4:6, 4:6) = sum(KAB_I,3) + sum(KAA_II,3) + sum(KAB_III,3) + sum(KAB_IV1,3) + sum(KAB_IV2,3) + sum(KAA_V,3)- DBB;
        e = eig(D);
        w = sort(e);
        omega(i,j,:) = real(sqrt(w))/sqrt(m);
    end
end

%% Plotting dispersion relationship
figure('OuterPosition',[100 100 1200 800])   
for i=1:6
    subplot(2,3,i); 
    surfc(X,Y,omega(:,:,i),'EdgeAlpha',0.1),shading flat;
    set(gca,'linewidth',1.5,'FontSize',14);
    zlabel('ω', 'FontSize', 14);
    xlabel("k_x",'FontSize', 14);
    ylabel("k_y",'FontSize', 14);
    hold on
    t=0:2*pi/6:2*pi;
    y0=sin(t)*2*pi/3;x0=cos(t)*2*pi/3;
    plot3(x0.*(max(max(X))/max(x0)),y0.*(max(max(Y))/max(y0)),ones(1,7),'b','linewidth',1)
%
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