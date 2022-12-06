clc;
close all;
clear;

w1=importdata('C:\Users\1\Desktop\项目\温度补偿\温度实验数据\7.21 21点-次日9点/20210721.csv');
w2=importdata('C:/Users/1/Desktop/项目/温度补偿/温度实验数据/7.21 21点-次日9点/measure_2021-07-21212230.csv');

displacement=w2.data(:,1);
temperature(:,1)=w1.data(:,10);
[average_tem,average_dis_ch,average_pre] =data_handle(displacement,temperature(:,1),displacement);
time(1,:) = 0:10/60/60:(length(average_tem)-1)*1/6/60;

font_size = 12;
data_num  = 2520;
subplot(2,1,1)
plot(time(1:data_num),average_dis_ch(1:data_num))
xlabel('time/h','FontName','Times New Roman','FontSize',font_size);
ylabel('displacement/nm','FontName','Times New Roman','FontSize',font_size);
set(gca,'FontSize',12) ;
subplot(2,1,2)
plot(time(1:data_num),average_tem(1:data_num))
xlabel('time/h','FontName','Times New Roman','FontSize',font_size);
ylabel('temperature/℃','FontName','Times New Roman','FontSize',font_size);
set(gca,'FontSize',12) ;