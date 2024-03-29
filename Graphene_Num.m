%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        Graphene_Num                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load DOS.mat frequency pdf
newcolors = [234, 32, 39; 0, 98, 102; 27, 20, 100; 87, 88, 187; 111, 30, 81;
             238, 90, 36; 0, 148, 50; 6, 82, 221; 153, 128, 250; 131, 52, 113;
             247, 159, 31; 163, 203, 56; 18, 137, 167; 217, 128, 250; 181, 52, 113;
             255, 195, 18; 196, 229, 56; 18, 203, 196; 253, 167, 223; 237, 76, 103]./255; 

TX = [300 400 500 600 700 800];
omega = frequency;
h_bar = 1.0546e-34;
kB = 1.38065e-23;
%%

n0 = zeros(length(TX),length(omega));
for i = 1:length(TX)
    for j = 1:length(omega)
        n0(i,j) = 1/(exp(h_bar*omega(j)/kB/TX(i))-1);
    end
end

figure('OuterPosition',[100 100 600 450])
semilogy(omega/1e14,n0,'linewidth',2)
hold on
yline(1,'linewidth',1.5)
set(gca,'linewidth',1.5,'FontSize',14);
colororder(newcolors)
legend('300K','400K','500K','600K','700K','800K')
legend boxoff
xlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
ylabel('$ \overline n (\omega,T) $','Interpreter','latex','FontSize',20,'FontWeight','bold');
ylim([1e-2 100])
xlim([0 3])
%%
temp = ["300K";"400K";"500K";"600K";"700K";"800K"];
figure('OuterPosition',[100 100 1300 800])
tiledlayout(2,3)
for i = 1:6
    nexttile
    N = n0(i,:).*pdf*3.6484e25;
    semilogy(omega/1e14,N,'linewidth',2)
    title(temp(i));
    set(gca,'linewidth',1.5,'FontSize',14);
    colororder(newcolors)
    legend('ZA','TA','LA','ZO','TO','LO')
    legend boxoff
    ylim([1e8 1e14])
    xlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
    ylabel('$ N (\omega,T) $','Interpreter','latex','FontSize',20,'FontWeight','bold');
end
%%