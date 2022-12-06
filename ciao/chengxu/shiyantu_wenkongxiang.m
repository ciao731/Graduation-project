clc;
close all;
clear;

w1=importdata('C:\Users\1\Desktop\项目\温度补偿\温度实验数据\对比测量2\20211207-t.csv');
w2=importdata('C:\Users\1\Desktop\项目\温度补偿\温度实验数据\对比测量2\20211207-m.csv');
w3=importdata('C:\Users\1\Desktop\项目\温度补偿\温度实验数据\对比测量2\20211207-p.csv');%  毛刺很大0.05，补偿效果正常

displacement=w2.data(:,1);
temperature(:,1)=w1.data(:,10);
pressure=w3.data(:,1);
[average_tem,average_dis_ch,average_pre] =data_handle(displacement,temperature(:,1),pressure);
time(1,:) = 0:10/60/60:(length(average_dis_ch)-1)*1/6/60;

font_size = 12;
data_num  = 4500;

hfigure=figure('numbertitle','off','name','分段补偿效果对比');
subplot(2,1,1)
plot(time(1:data_num),average_dis_ch(1:data_num))
xlabel('time/h','FontName','Times New Roman','FontSize',font_size);
ylabel('displacement/nm','FontName','Times New Roman','FontSize',font_size);
xlim([time(1) time(data_num)])
set(gca,'FontSize',12) ;

subplot(2,1,2)
        AX=plotyy(time,average_tem,time,(average_pre/1000));
        pos3=axis(AX(1));
        pos4=axis(AX(2));
        xlabel('time/h','FontName','Times New Roman','FontSize',font_size);
        set(get(AX(1),'Ylabel'),'string','temperature/℃','FontName','Times New Roman','FontSize',font_size);
                xlim([time(1) time(data_num)]);
        set(get(AX(2),'Ylabel'),'string','pressure/kPa','FontName','Times New Roman','FontSize',font_size);
        set(AX(1),'xlim',[time(1) time(data_num)]);
        set(AX(2),'xlim',[time(1) time(data_num)]);

        

        set(AX(1),'FontSize',font_size) ;
        set(AX(2),'FontSize',font_size) ;
        hold on;
        line31=plot(time(1:data_num),average_tem(1:data_num),'Color',[0 0.447 0.741]);
        hold on;
       
        line32=plot(time(1:data_num),(average_pre(1:data_num)/1000),'Color',[0.85 0.325 0.098]);
        hold on;
        h3=legend([line31 line32],'Temperature','Pressure'); 
        set(h3,'FontName','Times New Roman','FontSize',8,'FontWeight','normal');
        
set(hfigure,'position',[100 100 800 500]);