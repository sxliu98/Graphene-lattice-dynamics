%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Graphene_DOS                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
load  data.mat
a0 = 1.42e-10;
newcolors = [234, 32, 39; 0, 98, 102; 27, 20, 100; 87, 88, 187; 111, 30, 81;
             238, 90, 36; 0, 148, 50; 6, 82, 221; 153, 128, 250; 131, 52, 113;
             247, 159, 31; 163, 203, 56; 18, 137, 167; 217, 128, 250; 181, 52, 113;
             255, 195, 18; 196, 229, 56; 18, 203, 196; 253, 167, 223; 237, 76, 103]./255; 
%%
[n1,n2,n3]=size(W);
center=[(n1+1)/2,(n2+1)/2];
A=[];
W_new = W;
for i=1:6
    w=W_new(:,:,i);
    z1=mean(mean(w(fix(center(1)):ceil(center(1)),fix(center(2)):ceil(center(2)))));%中心点z值
    z2=mean(mean(w));
    if z1>=z2
        A(i)=1;
    else
        A(i)=2;
    end
end

for t=1:5
    a1=find(A==1);
    a2=find(A==2);
    if isempty(a1)==0
        for i=1:n1
            for j=1:n2
                [~,b]=sort(W_new(i,j,a1));
                b=b(:);
                W_new(i,j,a1)=W_new(i,j,a1(b));
            end
        end
    end
    if isempty(a2)==0
        for i=1:n1
            for j=1:n2
                [~,b]=sort(W_new(i,j,a2));
                b=b(:);
                W_new(i,j,a2)=W_new(i,j,a2(b));
            end
        end
    end
    for i=1:length(a1)
        for j=1:length(a2)
            w1=W_new(:,:,a1(i));
            w2=W_new(:,:,a2(j));
            a=find(isnan(X)==0);
            b=find(abs(w1(a)-w2(a))<1e11);
            if isempty(b)==0
                xyz=[X(a(b)),Y(a(b)),w1(a(b))];
                theta=[];
                for k=1:size(xyz,1)
                    if xyz(k,1)>0 && xyz(k,2)>0
                        theta(k,1)=atan(xyz(k,2)/xyz(k,1));
                    elseif xyz(k,1)>0 && xyz(k,2)<0
                        theta(k,1)=atan(xyz(k,2)/xyz(k,1))+2*pi;
                    elseif xyz(k,1)<0 && xyz(k,2)<0
                        theta(k,1)=atan(xyz(k,2)/xyz(k,1))+pi;
                    elseif xyz(k,1)<0 && xyz(k,2)>0
                        theta(k,1)=atan(xyz(k,2)/xyz(k,1))+pi;
                    elseif xyz(k,1)>0 && xyz(k,2)==0
                        theta(k,1)=0;
                    elseif xyz(k,1)==0 && xyz(k,2)<0
                        theta(k,1)=1.5*pi;
                    elseif xyz(k,1)<0 && xyz(k,2)==0
                        theta(k,1)=pi;
                    elseif xyz(k,1)==0 && xyz(k,2)>0
                        theta(k,1)=0.5*pi;
                    end 
                end
                [~,index]=sort(theta);
                xyz=xyz(index,:);
                xyz=[xyz;xyz(1,:)];
                [in,on]=inpolygon(X(a),Y(a),xyz(:,1),xyz(:,2));
                c1=[find(in==1);find(on==1)];
                c1=unique(c1);
                aa=a(c1);
                c=setdiff(a,aa);
                ww=[w1(c),w2(c)]';
                w1(c)=min(ww);
                w2(c)=max(ww);
                W_new(:,:,a1(i))=w1;
                W_new(:,:,a2(j))=w2;
            end
        end
    end
end
%%
mode = ["ZA";"TA";"LA";"ZO";"TO";"LO"];
figure('OuterPosition',[100 100 1200 800])   
tiledlayout(2,3)
for i=1:6
    nexttile
    surfc(X*a0,Y*a0,W_new(:,:,i)/1e14,'EdgeAlpha',0.1),shading flat;
    title(mode(i));
    set(gca,'linewidth',1.5,'FontSize',14);
    zlim([0 max(max(max(W_new)))/1e14]);
    xlim([-2*pi/3-0.2,2*pi/3+0.2]);
    ylim([-pi/sqrt(3)-0.2,pi/sqrt(3)+0.2]);
    xlabel('$k_x a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    ylabel('$k_y a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    zlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
    hold on
    t=0:2*pi/6:2*pi;
    y0=sin(t)*2*pi/3;x0=cos(t)*2*pi/3;
    plot3(x0,y0,zeros(1,7),'b','linewidth',1)
end
cb = colorbar;
cb.Layout.Tile = 'east';

save frequency_mode.mat X Y W_new
%%
figure('OuterPosition',[100 100 1300 800])   
tiledlayout(2,3)
for i=1:6
    nexttile
    [M,c]=contour(X*a0,Y*a0,W_new(:,:,i)/1e14,18);
    c.LineWidth = 2;
    hold on
    nn = 9;
    [DX,DY] = gradient(W_new(1:nn:200,1:nn:200,i),2*pi/300/a0*nn,pi/sqrt(3)/100/a0*nn);
    if i == 2 || i ==3
        DX(abs(DX)>3e4) = NaN; DY(abs(DY)>3e4) = NaN;
    else
        DX(abs(DX)>1e4) = NaN; DY(abs(DY)>1e4) = NaN;
    end
    DX([1,end],:) = NaN; DY([1,end],:) = NaN;
    quiver(X(1:nn:200,1:nn:200)*a0,Y(1:nn:200,1:nn:200)*a0,DX,DY,'k');
    title(mode(i));
    set(gca,'linewidth',1.5,'FontSize',14);
    zlim([0 max(max(max(W_new)))/1e14]);
    xlim([-2*pi/3-0.2,2*pi/3+0.2]);
    ylim([-pi/sqrt(3)-0.2,pi/sqrt(3)+0.2]);
    xlabel('$k_x a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    ylabel('$k_y a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    zlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
    hold on
    t=0:2*pi/6:2*pi;
    y0=sin(t)*2*pi/3;x0=cos(t)*2*pi/3;
    plot3(x0,y0,zeros(1,7),'b','linewidth',1)
end

%%
figure('OuterPosition',[100 100 1300 800])   
tiledlayout(2,3)
for i=1:6
    nexttile
    [DX,DY] = gradient(W_new(:,:,i),2*pi/300/a0,pi/sqrt(3)/100/a0);
    if i == 2 || i ==3
        DX(abs(DX)>3e4) = NaN; DY(abs(DY)>3e4) = NaN;
    else
        DX(abs(DX)>1e4) = NaN; DY(abs(DY)>1e4) = NaN;
    end
    DX([1,end],:) = NaN; DY([1,end],:) = NaN;
    quiver(X*a0,Y*a0,DX,DY,'k');        
    title(mode(i));
    set(gca,'linewidth',1.5,'FontSize',14);
    xlim([-2*pi/3-0.2,2*pi/3+0.2]);
    ylim([-pi/sqrt(3)-0.2,pi/sqrt(3)+0.2]);
    xlabel('$k_x a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    ylabel('$k_y a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    zlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
    hold on
    t=0:2*pi/6:2*pi;
    y0=sin(t)*2*pi/3;x0=cos(t)*2*pi/3;
    plot3(x0,y0,zeros(1,7),'b','linewidth',1)
end

%%
figure('OuterPosition',[100 100 600 450])   
for i=1:6
    [DX,DY] = gradient(W_new(:,:,i),2*pi/300/a0,pi/sqrt(3)/100/a0);
    if i == 2 || i ==3
        DX(abs(DX)>3e4) = NaN; DY(abs(DY)>3e4) = NaN;
    else
        DX(abs(DX)>1e4) = NaN; DY(abs(DY)>1e4) = NaN;
    end
    DX([1:10,end-10:end],:) = NaN; DY([1:10,end-10:end],:) = NaN;

    vg = sqrt(DX.^2+DY.^2);
    omega_i = W_new(1:16:200,1:16:200,i);
    omega = omega_i(:);
    vel = vg(1:16:200,1:16:200);
    velocity = vel(:);
    plot(omega/1e14,velocity/1e3,'o','LineWidth',1.5)
    colororder(newcolors)
    hold on
end
    set(gca,'linewidth',1.5,'FontSize',14);    
    xlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
    ylabel('$ v_g(\omega),10^3 $ m/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
    legend('ZA','TA','LA','ZO','TO','LO')
    legend boxoff
%%

%%
W_new([1:10,end-10:end],:) = NaN;
figure
for i =1:6
    Wi = W_new(:,:,i);
    flatten_w = Wi(:);
    bins = 0.01e14:1e12:3.05e14;
    h = histogram(flatten_w,bins,'Normalization','pdf');
    hold on
    pdf(i,:) = h.Values;
    frequency = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
end
set(gca,'linewidth',1.5,'FontSize',14);
set(h,'edgecolor','none');
xlabel('Frequency');
ylabel('DOS');
xlim([0 3e14])
ylim([0 8e-14])
%%
newcolors = [234, 32, 39; 0, 98, 102; 27, 20, 100; 87, 88, 187; 111, 30, 81;
             238, 90, 36; 0, 148, 50; 6, 82, 221; 153, 128, 250; 131, 52, 113;
             247, 159, 31; 163, 203, 56; 18, 137, 167; 217, 128, 250; 181, 52, 113;
             255, 195, 18; 196, 229, 56; 18, 203, 196; 253, 167, 223; 237, 76, 103]./255; 
figure('OuterPosition',[100 100 600 450])
plot(frequency*1e-14,pdf,'-','LineWidth',1.5);
hold on
plot(frequency*1e-14,sum(pdf,1),'k-','LineWidth',2);
colororder(newcolors)
set(gca,'linewidth',1.5,'FontSize',14);
xlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
set(gca,'yticklabel',[])
ylabel('$g(\omega)$ [a.u.]','Interpreter','latex','FontSize',20,'FontWeight','bold');
legend('ZA','TA','LA','ZO','TO','LO','Total')
legend boxoff
xlim([0 3.5])
ylim([0 5e-14])

%%
save DOS.mat frequency pdf
%%
h_bar = 1.0546e-34;
kB = 1.38065e-23;
w_width = h.BinWidth;
capacity_per_x = @(x) kB*(x.^2).*exp(x)./((exp(x)-1).^2);
T_limit = 1e4;
x_limit = (h_bar/kB)*frequency' * (1./T_limit);
capacity_limit = capacity_per_x(x_limit).*sum(pdf)';
capacity_limit = sum(capacity_limit)*w_width;
c_lim_true = 24.943*1000/12;
coff = c_lim_true./capacity_limit;
T_list1 = (1:1:3000);
x_mat1 = (h_bar/kB)*frequency' * (1./T_list1);
capacity_1 = capacity_per_x(x_mat1).*(pdf(1,:))';
capacity_1 = sum(capacity_1)*w_width*coff;
capacity_2 = capacity_per_x(x_mat1).*(pdf(2,:))';
capacity_2 = sum(capacity_2)*w_width*coff;
capacity_3 = capacity_per_x(x_mat1).*(pdf(3,:))';
capacity_3 = sum(capacity_3)*w_width*coff;
capacity_4 = capacity_per_x(x_mat1).*(pdf(4,:))';
capacity_4 = sum(capacity_4)*w_width*coff;
capacity_5 = capacity_per_x(x_mat1).*(pdf(5,:))';
capacity_5 = sum(capacity_5)*w_width*coff;
capacity_6 = capacity_per_x(x_mat1).*(pdf(6,:))';
capacity_6 = sum(capacity_6)*w_width*coff;

capacity_tol = capacity_per_x(x_mat1).*sum(pdf)';
capacity_tol = sum(capacity_tol)*w_width*coff;

figure('OuterPosition',[100 100 600 450])
semilogy(T_list1,capacity_1,"LineStyle","-",'LineWidth',1.5)
hold on
plot(T_list1,capacity_2,"LineStyle","-",'LineWidth',1.5)
plot(T_list1,capacity_3,"LineStyle","-",'LineWidth',1.5)
plot(T_list1,capacity_4,"LineStyle","-",'LineWidth',1.5)
plot(T_list1,capacity_5,"LineStyle","-",'LineWidth',1.5)
plot(T_list1,capacity_6,"LineStyle","-",'LineWidth',2)
plot(T_list1,capacity_tol,'k',"LineStyle","-",'LineWidth',2)
colororder(newcolors)
T_exp = [10.79705 59.59274 96.57802 142.57651 214.06065 294.55801 397.43302 497.82162 598.21022 698.28802 798.67662 893.16 997.89981 1098.91001 1191.5286 1404.42702];
Cv_exp = [14.35108 82.31417 173.42559 287.24368 482.75724 725.13698 999.63935 1236.76266 1403.80027 1554.48455 1661.9493 1732.03501 1790.43977 1853.51691 1889.14381 1956.30928];
plot(T_exp,Cv_exp,'rs','MarkerFacecolor','r')
set(gca,'linewidth',1.5,'FontSize',14);
xlim([0 800]);
ylim([1 3000])
xlabel('$T$, K','Interpreter','latex','FontSize',20,'FontWeight','bold')
ylabel('$C_v(T)$, J/(kg.K)','Interpreter','latex','FontSize',20,'FontWeight','bold')
legend('ZA','TA','LA','ZO','TO','LO','Total','Experiment')
legend boxoff
hold off