clc,clear,close all;        
x1 = [8 16 32 48 64 80 96];
y1 = [0.267783 0.277327 0.499198 0.639440 0.899527 1.011029 1.228623];
x2 = [ 0 32 32.01            64 64.01          96 96.01];
y2 = [ 6.79 6.79 9.38         9.38 11.97         11.97 14.56];
AX=plotyy(x1,y1,x2,y2);
        pos3=axis(AX(1));
        pos4=axis(AX(2));
        xlabel('Population size','FontName','Times New Roman','FontSize',8);
        set(get(AX(1),'Ylabel'),'string','time/s','FontName','Times New Roman','FontSize',8);

        set(get(AX(2),'Ylabel'),'string','time/us','FontName','Times New Roman','FontSize',8);
        set(AX(1),'xlim',[16 96]);
        set(AX(2),'xlim',[16 96]);
        title('(c) Temperature and Pressure','FontSize',21,'position',[0.194 26.252  0]);
        

        set(AX(1),'FontSize',12) ;
        set(AX(2),'FontSize',12) ;
 
        hold on;
        c = polyfit(x1, y1, 2); %进行拟合，c为2次拟合后的系数
        d = polyval(c, x1, 1); %拟合后，每一个横坐标对应的值即为d
        line1 = plot(x1, d, 'b'); %拟合后的曲线
        hold on;
        line2= plot(x2, y2, 'g'); %拟合后的曲线
        hold on;
        h3=legend([line1 line2],'Software','Hardware'); 
        set(h3,'FontName','Times New Roman','FontSize',10,'FontWeight','normal');

      

% clc,clear,close all;
% a = [16 32 48 64 80 96]; %横坐标
% b = [0.277327 0.499198 0.639440 0.899527 1.011029 1.228623]; %纵坐标
% plot(a, b, 'b'); %自然状态的画图效果
% hold on;
% %第一种，画平滑曲线的方法
% c = polyfit(a, b, 2); %进行拟合，c为2次拟合后的系数
% d = polyval(c, a, 1); %拟合后，每一个横坐标对应的值即为d
% plot(a, d, 'r'); %拟合后的曲线
% 
% plot(a, b, '*'); %将每个点 用画出来
% hold on;