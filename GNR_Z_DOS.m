%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          ZGNR_DOS                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;

load DOS.mat frequency pdf
w_grap = frequency; g_grap = pdf;

newcolors = [234, 32, 39; 0, 98, 102; 27, 20, 100; 87, 88, 187; 111, 30, 81;
             238, 90, 36; 0, 148, 50; 6, 82, 221; 153, 128, 250; 131, 52, 113;
             247, 159, 31; 163, 203, 56; 18, 137, 167; 217, 128, 250; 181, 52, 113;
             255, 195, 18; 196, 229, 56; 18, 203, 196; 253, 167, 223; 237, 76, 103]./255; 

Nx = [2 4 8];

a = 1.42e-10;
f = diag([36.50, 24.50, 9.82, 8.80, -3.23, -0.40, 3.00, -5.25, 0.15, -1.92, 2.29, -0.58])*10;
m = 1.99e-26;

[KAB1, KBA1] = rotate(2/3*pi, f(1:3, 1:3));
[KAA, KBB] = rotate(1/3*pi, K(1/6*pi, f(4:6, 4:6)));
[KAB3, KBA3] = rotate(2/3*pi, K(1/3*pi, f(7:9, 7:9)));
[KAB4f, KBA4f] = rotate(2/3*pi, K(acos(2.5/sqrt(7)), f(10:12, 10:12)));
[KAB4s, KBA4s] = rotate(2/3*pi, K(2*pi-acos(2.5/sqrt(7)), f(10:12, 10:12)));
KAB4 = cat(3,KAB4f, KAB4s);
KBA4 = cat(3,KBA4f, KBA4s);

dlist = 1000;

klist = linspace(0,2*pi/3/a,dlist);

figure('OuterPosition',[0 0 1200 600])
for n = 1:3
    N = Nx(n);
    len = 12*N;
    omega = zeros(len,dlist);

for i = 1:length(klist)
    kx = klist(i);
    D = zeros(len,len); 
    D00 = zeros(3,3);
    Dxx = sum(KBA1,3) + sum(KBB,3) + sum(KBA3,3) + sum(KBA4,3);

%%% D5
    D52 = KBA4(:,:,3)*exp(1i * kx * (sqrt(3)/2*a)) + KBA4(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D53 = KBB(:,:,5)*exp(1i * kx * (-sqrt(3)/2*a)) + KBB(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));
    D54 = KBA1(:,:,3)*exp(1i * kx * 0) + KBA3(:,:,2)*exp(1i * kx * (-sqrt(3)*a)) + KBA3(:,:,3)*exp(1i * kx * (sqrt(3)*a));
    D55 = KBB(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBB(:,:,4)*exp(1i * kx * (-sqrt(3)*a));
    D56 = KBA1(:,:,1)*exp(1i * kx * (sqrt(3)/2*a)) + KBA1(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KBA4(:,:,2)*exp(1i * kx * (-3*sqrt(3)/2*a)) + KBA4(:,:,4)*exp(1i * kx * (3*sqrt(3)/2*a));
    D57 = KBB(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KBB(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a));
    D58 = KBA3(:,:,1)*exp(1i * kx * 0) + KBA4(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBA4(:,:,5)*exp(1i * kx * (-sqrt(3)*a));

%%% D6
    D63 = KAB3(:,:,1)*exp(1i * kx * 0) + KAB4(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAB4(:,:,5)*exp(1i * kx * (sqrt(3)*a));
    D64 = KAA(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KAA(:,:,3)*exp(1i * kx * (sqrt(3)/2*a));
    D65 = KAB1(:,:,1)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB1(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KAB4(:,:,2)*exp(1i * kx * (3*sqrt(3)/2*a)) + KAB4(:,:,4)*exp(1i * kx * (-3*sqrt(3)/2*a));
    D66 = KAA(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAA(:,:,4)*exp(1i * kx * sqrt(3)*a);
    D67 = KAB1(:,:,3)*exp(1i * kx * 0) + KAB3(:,:,2)*exp(1i * kx * (sqrt(3)*a)) + KAB3(:,:,3)*exp(1i * kx * (-sqrt(3)*a));
    D68 = KAA(:,:,5)*exp(1i * kx * (sqrt(3)/2*a)) + KAA(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D69 = KAB4(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB4(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));

%%% D7
    D74 = KBA4(:,:,3)*exp(1i * kx * (sqrt(3)/2*a)) + KBA4(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D75 = KBB(:,:,5)*exp(1i * kx * (-sqrt(3)/2*a)) + KBB(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));
    D76 = KBA1(:,:,3)*exp(1i * kx * 0) + KBA3(:,:,2)*exp(1i * kx * (-sqrt(3)*a)) + KBA3(:,:,3)*exp(1i * kx * (sqrt(3)*a));
    D77 = KBB(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBB(:,:,4)*exp(1i * kx * (-sqrt(3)*a));
    D78 = KBA1(:,:,1)*exp(1i * kx * (sqrt(3)/2*a)) + KBA1(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KBA4(:,:,2)*exp(1i * kx * (-3*sqrt(3)/2*a)) + KBA4(:,:,4)*exp(1i * kx * (3*sqrt(3)/2*a));
    D79 = KBB(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KBB(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a));
    D710 = KBA3(:,:,1)*exp(1i * kx * 0) + KBA4(:,:,1)*exp(1i * kx * (sqrt(3)*a)) + KBA4(:,:,5)*exp(1i * kx * (-sqrt(3)*a));

%%% D8
    D85 = KAB3(:,:,1)*exp(1i * kx * 0) + KAB4(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAB4(:,:,5)*exp(1i * kx * (sqrt(3)*a));
    D86 = KAA(:,:,2)*exp(1i * kx * (-sqrt(3)/2*a)) + KAA(:,:,3)*exp(1i * kx * (sqrt(3)/2*a));
    D87 = KAB1(:,:,1)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB1(:,:,2)*exp(1i * kx * (sqrt(3)/2*a)) + KAB4(:,:,2)*exp(1i * kx * (3*sqrt(3)/2*a)) + KAB4(:,:,4)*exp(1i * kx * (-3*sqrt(3)/2*a));
    D88 = KAA(:,:,1)*exp(1i * kx * (-sqrt(3)*a)) + KAA(:,:,4)*exp(1i * kx * (sqrt(3)*a));
    D89 = KAB1(:,:,3)*exp(1i * kx * 0) + KAB3(:,:,3)*exp(1i * kx * (-sqrt(3)*a)) + KAB3(:,:,2)*exp(1i * kx * (sqrt(3)*a));
    D810 = KAA(:,:,5)*exp(1i * kx * (sqrt(3)/2*a)) + KAA(:,:,6)*exp(1i * kx * (-sqrt(3)/2*a));
    D811 = KAB4(:,:,3)*exp(1i * kx * (-sqrt(3)/2*a)) + KAB4(:,:,6)*exp(1i * kx * (sqrt(3)/2*a));

   
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

    flatten_w = omega(:);
    bins = 0.01e14:1e12:3.05e14;
    figure
    h = histogram(flatten_w,bins,'Normalization','pdf');
    set(gca,'linewidth',1.5,'FontSize',14);
    set(h,'edgecolor','none');
    xlabel('Frequency');
    ylabel('DOS');
    xlim([0 3.05e14])
    ylim([0 2.5e-14])
    
    pdf(n,:) = h.Values;
    frequency = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
    
end

%%
save DOS_ZGNR-8-16-32.mat  frequency pdf
%%
figure('OuterPosition',[100 100 600 450])
plot(frequency,pdf(1,:),'r-','LineWidth',2);
hold on
plot(frequency,pdf(2,:)+2.2e-14,'b-','LineWidth',2);
plot(frequency,pdf(3,:)+4.4e-14,'g-','LineWidth',2);
plot(w_grap,sum(g_grap)/5+6.6e-14,'k-','LineWidth',2);
set(gca,'linewidth',1.5,'FontSize',14);
xlabel('Frequency');
ylabel('DOS');
xlim([0 3.1e14])
ylim([0 10e-14])
xlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
set(gca,'yticklabel',[])
ylabel('$g(\omega)$ [a.u.]','Interpreter','latex','FontSize',20,'FontWeight','bold');
legend('8-ZGNR','16-ZGNR','32-ZGNR','graphene')
legend boxoff
