clc;
close all;
clear;
tic
% ��20 06��ǰ�룩 ����30 04
for i = 1:1 %%%�ԱȲ���2 С�¶ȱ仯
  
w1=importdata('C:\Users\1\Desktop\��Ŀ\�¶Ȳ���\�¶�ʵ������\�ԱȲ���2\20211117-d-t.csv');
w2=importdata('C:\Users\1\Desktop\��Ŀ\�¶Ȳ���\�¶�ʵ������\�ԱȲ���2\20211117-d-m.csv');
w3=importdata('C:\Users\1\Desktop\��Ŀ\�¶Ȳ���\�¶�ʵ������\�ԱȲ���2\20211117-d-p.csv');
end

temperature(:,1)=w1.data(:,10);





channel_num = 12;
landa0 =632.991354;

for i=1:1 %��������
    if exist('w3')
        pressure=w3.data(:,1);
    else
        for i=1:length(temperature(:,1))
            pressure(i)=1.013*10^3;
        end
    end
    if (channel_num == 1)
        o_length = 90;
        o_length2 = 90;
        displacement=w2.data(:,1);
        displacement2=w2.data(:,1);
        two_channel_com = 0;
    elseif (channel_num == 2)
        o_length = 45;
        o_length2 = 45;
        displacement=w2.data(:,2);
        displacement2=w2.data(:,2);
        two_channel_com = 0;
    elseif (channel_num == 12)
        o_length = 90;
        o_length2 = 45;
        displacement=w2.data(:,1);
        displacement2=w2.data(:,2);
        two_channel_com = 1;
    else
        error('���벻����Ҫ��');
    end
end
for i=1:1
xishu = 1.4;
d_end = ceil(length(displacement)/xishu);
t_end = ceil(length(temperature)/xishu);
p_end = ceil(length(pressure)/xishu);

[average_tem,average_dis_ch1,average_pre] =data_handle(displacement,temperature,pressure);
[average_tem,average_dis_ch2,average_pre] =data_handle(displacement2,temperature,pressure);
% [average_tem,average_dis_ch1,average_pre] =data_handle(displacement(1:d_end),temperature(1:t_end,1),pressure(1:p_end));
% [average_tem,average_dis_ch2,average_pre] =data_handle(displacement2(1:d_end),temperature(1:t_end,1),pressure(1:p_end));


time = 0:10/60/60:(length(average_tem)-1)*1/6/60;
plot_begin=1;
plot_len=length(average_dis_ch1);

average_data_plot = 0;                          %λ�ơ��¶ȡ���ѹ�������ݾ���data_handle������ͼ

zheshelv_compare_with_tem=0;                    %pso�Ż�ǰ�������ʶԱ��Լ��¶�
pso_result_improve=0;                           %pso�Ż�ǰ�󲹳�Ч���Աȣ����¶ȣ�
edlen_data_compare_improve=0;                   %pso�Ż�ǰ��edlen����ֵ�Ա�
edpso_burst_compare       =0;                   %����pso�ֶβ�����δ�ֶβ�����Ч���Ա�

all_method_result_fliter_compare_with_tem=0;    %�������Ĳ���������¶ȵĶԱȣ��˲���
all_method_data_and_result=0;                   %�������������ݺͽ���Ա�
all_order_pso_result=0;                         %����pso����Ч���Ա�
lunwentu_daiyouhua = 1;
lunwentu_wuyouhua  = 1;

end


train_num = 6;
train_index = ceil(length(average_tem)/train_num);
gbest_edpso = [2.68*10^-9, -9.36*10^-7];
% [fitnessgbest_e1,gbest_e1]=PSO(average_tem(train_index*2:train_index*3),average_dis_ch1(train_index*2:train_index*3),average_pre(train_index*2:train_index*3),o_length,10,1); %1204
[fitnessgbest_e1,gbest_edpso]=PSO(average_tem,average_dis_ch1,average_pre,o_length,10,1); %1204
gbest_0=[-4.182298472036175,-239.4007268006262,-3916.842551820620,47.616565483641246,-7.352948152555437];%01-0.2014
gbest_2=[-8.949034637700064,-102.0346592028632,-8.963643789028142];%01-0.2022







gbest_edlen=[2.68*10^-9, -9.36*10^-7];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i = 1: length(average_dis_ch1)
data_pso_2(i) = gbest_2(1)*(average_tem(i)-average_tem(1))^2+gbest_2(2)*(average_tem(i)-average_tem(1))+gbest_2(3);
data_pso_edlen(i) = gbest_0(1)*(average_tem(i)-average_tem(1))^2+gbest_0(2)*(average_tem(i)-average_tem(1))+gbest_0(3)*10^3/(1+gbest_0(4)*average_tem(i))-gbest_0(3)*10^3/(1+gbest_0(4)*average_tem(1))+gbest_0(5);


n_edlen(i)                = gbest_edlen(1)  * average_pre(i) + gbest_edlen(2)  * average_tem(i);
n_edpso(i)                = gbest_edpso(1) * average_pre(i) + gbest_edpso(2) * average_tem(i);
compensation_edlen_ch1(i) = (o_length  / 1000 * (n_edlen(i)-n_edlen(1)) * 10^9);
compensation_edlen_ch2(i) = (o_length2 / 1000 * (n_edlen(i)-n_edlen(1)) * 10^9);
compensation_edpso_ch1(i) = (o_length  / 1000 * (n_edpso(i)-n_edpso(1)) * 10^9);
compensation_edpso_ch2(i) = (o_length2 / 1000 * (n_edpso(i)-n_edpso(1)) * 10^9);
end

burst_width = 80;
%%��burst_widthΪ���ڴ�С���л����ִ�����ѵ��
for i = 1 : fix(length(average_tem)/burst_width)
    [fitnessgbest_e_burst_ch1{i},gbest_e_burst_ch1{i}] =PSO([average_tem((i-1)*burst_width+1 : i*burst_width)],[average_dis_ch1((i-1)*burst_width+1 : i*burst_width)],[average_pre((i-1)*burst_width+1 : i*burst_width)],o_length,10,0);
    [fitnessgbest_e_burst_ch2{i},gbest_e_burst_ch2{i}] =PSO([average_tem((i-1)*burst_width+1 : i*burst_width)],[average_dis_ch2((i-1)*burst_width+1 : i*burst_width)],[average_pre((i-1)*burst_width+1 : i*burst_width)],o_length2,10,0);
    for t = 1:burst_width
        index = (i-1)*burst_width+t;
        n_edpso_burst_ch1(index)                = gbest_e_burst_ch1{i}(1)  * average_pre(index) + gbest_e_burst_ch1{i}(2)  * average_tem(index);
        n_edpso_burst_ch2(index)                = gbest_e_burst_ch2{i}(1)  * average_pre(index) + gbest_e_burst_ch2{i}(2)  * average_tem(index);
        compensation_edpso_burst_ch1(index) = (o_length  / 1000 * (n_edpso_burst_ch1(index)-n_edpso_burst_ch1((i-1)*burst_width+1)) * 10^9);
        compensation_edpso_burst_ch2(index) = (o_length2 / 1000 * (n_edpso_burst_ch2(index)-n_edpso_burst_ch2((i-1)*burst_width+1)) * 10^9);
    end
    
end
figure();
subplot(2,1,1)
plot(compensation_edpso_burst_ch1)
subplot(2,1,2)
plot(compensation_edpso_burst_ch2)
%%��ÿһ��ѵ���Ľ����ȥ��һ�ε��ۼ�����߲���һ��1*10�Ĺ̶����ڶ���һ�ε�ѵ�������ƽ���Ӷ������ۼ����
for i = 1:fix(length(average_tem)/burst_width)
    if(i~=1)
        burst_error_ch1 = mean(compensation_edpso_burst_ch1((i-1)*burst_width+1 : (i-1)*burst_width+10)) - mean(compensation_edpso_burst_ch1((i-1)*burst_width-10+1 : (i-1)*burst_width));
        compensation_edpso_burst_ch1((i-1)*burst_width+1 : i*burst_width) = compensation_edpso_burst_ch1((i-1)*burst_width+1 : i*burst_width) - burst_error_ch1;
        burst_error_ch2 = mean(compensation_edpso_burst_ch2((i-1)*burst_width+1 : (i-1)*burst_width+10)) - mean(compensation_edpso_burst_ch2((i-1)*burst_width-10+1 : (i-1)*burst_width));
        compensation_edpso_burst_ch2((i-1)*burst_width+1 : i*burst_width) = compensation_edpso_burst_ch2((i-1)*burst_width+1 : i*burst_width) - burst_error_ch2;       
    end
end

%%��ĩ�˲��㻬��ѵ�����ڵ�ʣ�����ݣ�ʹ�����һ�εĻ����ִ�ѵ��������в��������Ҳ����ۼ���
for i = index : length(average_tem)
    n_edpso_burst_ch1(i)                = gbest_e_burst_ch1{fix(length(average_tem)/burst_width)}(1)  * average_pre(i) + gbest_e_burst_ch1{fix(length(average_tem)/burst_width)}(2)  * average_tem(i);
    n_edpso_burst_ch2(i)                = gbest_e_burst_ch2{fix(length(average_tem)/burst_width)}(1)  * average_pre(i) + gbest_e_burst_ch2{fix(length(average_tem)/burst_width)}(2)  * average_tem(i);
    compensation_edpso_burst_ch1(i) = (o_length  / 1000 * (n_edpso_burst_ch1(i)-n_edpso_burst_ch1(index - burst_width +1)) * 10^9) - burst_error_ch1;
    compensation_edpso_burst_ch2(i) = (o_length2 / 1000 * (n_edpso_burst_ch2(i)-n_edpso_burst_ch2(index - burst_width +1)) * 10^9) - burst_error_ch2;
end

figure();
plot(compensation_edpso_burst_ch1)

data_pso                = data_pso_2;
result_pso              = average_dis_ch1 + data_pso;
result_pso_edlen        = average_dis_ch1 + data_pso_edlen;
result_edlen_ch1        = average_dis_ch1 + compensation_edlen_ch1;
result_edlen_ch2        = average_dis_ch2 + compensation_edlen_ch2;
result_edpso_ch1        = average_dis_ch1 + compensation_edpso_ch1;
result_edpso_ch2        = average_dis_ch2 + compensation_edpso_ch2;
result_edpso_burst_ch1  = average_dis_ch1 + compensation_edpso_burst_ch1;
result_edpso_burst_ch2  = average_dis_ch2 + compensation_edpso_burst_ch2;

data_pso_filter                 = lowp(data_pso               , 0.5 , 80 , 0.1 , 30 , 200);
result_pso_filter               = lowp(result_pso             , 0.5 , 80 , 0.1 , 30 , 200);
compensation_edlen_ch1_filter   = lowp(compensation_edlen_ch1 , 0.5 , 80 , 0.1 , 30 , 200);
data_pso_edlen_filter           = lowp(data_pso_edlen         , 0.5 , 80 , 0.1 , 30 , 200);
result_pso_edlen_filter         = lowp(result_pso_edlen       , 0.5 , 80 , 0.1 , 30 , 200);

result_edlen_ch1_filter         = lowp(result_edlen_ch1       , 0.5 , 80 , 0.1 , 30 , 200);
result_edlen_ch2_filter         = lowp(result_edlen_ch2       , 0.5 , 80 , 0.1 , 30 , 200);
result_edpso_ch1_filter         = lowp(result_edpso_ch1       , 0.5 , 80 , 0.1 , 30 , 200);
result_edpso_ch2_filter         = lowp(result_edpso_ch2       , 0.5 , 80 , 0.1 , 30 , 200);

result_edpso_burst_ch1_filter   = lowp(result_edpso_burst_ch1 , 0.5 , 80 , 0.1 , 30 , 200);
result_edpso_burst_ch2_filter   = lowp(result_edpso_burst_ch2 , 0.5 , 80 , 0.1 , 30 , 200);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error(1:8)=0;
for i = 1:length(average_dis_ch1)
    error(1) = error(1) + (average_dis_ch2(1,i)-0)^2;
    error(2) = error(2) + (average_dis_ch1(1,i)-0)^2;
    error(3) = error(3) + (result_edlen_ch2_filter(1,i)-0)^2;
    error(4) = error(4) + (result_edlen_ch1_filter(1,i)-0)^2;
    error(5) = error(5) + (result_edpso_ch2_filter(1,i)-0)^2;
    error(6) = error(6) + (result_edpso_ch1_filter(1,i)-0)^2;
    error(7) = error(7) + (result_edpso_burst_ch2_filter(1,i)-0)^2;
    error(8) = error(8) + (result_edpso_burst_ch1_filter(1,i)-0)^2;
end
for i = 1:8
    
error(i) =sqrt(error(i)/length(average_dis_ch1))

end


%%%%%%%%%%%%%%%%%%%%%%%%%%��ͼ����
for i = 1:1     %%λ�ơ��¶ȡ���ѹ�������ݾ���data_handle������ͼ
    if (average_data_plot == 1)
        figure('numbertitle','off','name','ƽ��������');
        subplot(3,1,1);
        plot(average_tem);
        title('�¶�');
        subplot(3,1,2);
        plot(average_dis_ch1);
        title('λ��');
        subplot(3,1,3);
        title('��ѹ');
        plot(average_pre)
    end
end
for i = 1:1     %%�������Ĳ���������¶ȵĶԱȣ��˲���
    if(all_method_result_fliter_compare_with_tem == 1)
        figure('numbertitle','off','name','����������Ч��');
        subplot(2,1,1);
        % plot(average_dis_ch1);
        % hold on;
        plot(result_pso_filter);
        hold on;
        plot(result_edlen_ch1_filter);
        hold on;
        plot(result_edpso_ch1_filter);
        hold on;
        plot(result_pso_edlen_filter);
        hold on;
        plot(average_dis_ch1);
        xlim([plot_begin plot_begin+plot_len]);
        xlabel('����'); %%���ú�����
        ylabel('λ��/nm'); %%����������
        title('λ������');
        legend('PSO','edlen','�Ż���edlen','pso-edlen','ԭʼ����');
        subplot(2,1,2);
        plot(average_tem);
        xlim([plot_begin plot_begin+plot_len]);
        xlabel('����'); %%���ú�����
        ylabel('�¶�/��'); %%����������
        title('�¶�����');
        legend('ͨ��1','ͨ��3','ͨ��5','ͨ��7');
    end
end
for i = 1:1     %%PSO�Ż�ǰ�������ʶԱ��Լ��¶�
    if (zheshelv_compare_with_tem == 1)
        figure('NumberTitle', 'off','name','pso�Ż�ǰ��edlen������')
        subplot(3,1,1);
        plot(n_edpso);
        title('pso�Ż�ǰedlen������Ч��');
        subplot(3,1,2);
        plot(n_edlen);
        title('pso�Ż���edlen������Ч��');
        subplot(3,1,3);
        plot(average_tem);
        title('�¶�����');
    end
end
for i = 1:1     %%pso�Ż�ǰ��ellen����Ч���Ա�
if(pso_result_improve==1)%%pso�Ż�ǰ�󲹳�Ч���Ա�
    figure('NumberTitle', 'off','name','pso�Ż�ǰ�����Ա�')
    subplot(2,1,1);
    plot(result_edpso_ch1_filter);
    hold on ;
    plot(result_edlen_ch1_filter);
    hold on;
    plot(average_dis_ch1);
    title('pso�Ż�ǰ��edlen����Ч��');
    legend('pso�Ż���','pso�Ż�ǰ','ԭʼ����');
    subplot(2,1,2);
    plot(average_tem);
    title('�¶�����');
end
end
for i = 1:1     %%�������������ݺͽ���Ա�
    if(all_method_data_and_result==1)
        figure('NumberTitle', 'off', 'Name', '�������������ݺͽ���Ա�');
        plot(average_dis_ch1);
        hold on
        plot(data_pso);
        hold on;
        plot(result_pso);
        hold on;
        plot(compensation_edlen_ch1);
        hold on;
        plot(result_edlen_ch1)
        hold on;
        plot(compensation_edpso_ch1);
        hold on;
        plot(result_edpso_ch1);
        hold on;
        plot(data_pso_edlen);
        hold on;
        plot(result_pso_edlen);
        xlim([plot_begin plot_begin+plot_len]);
        legend('ԭ����','PSO��������','PSO�������','edlen��������','edlen�������','�Ż���edlen��������','�Ż���edlen�������','PSO_edlen��������','PSO_edlen�������');
        title('�˲�ǰ�����ݶԱ�');
    end
end
for i = 1:1     %%����pso����Ч��
    if(all_order_pso_result==1)
        result_pso1_filter(:,1)=lowp(result_pso1,0.5,80,0.1,30,200);
        result_pso3_filter(:,1)=lowp(result_pso3,0.5,80,0.1,30,200);
        result_pso4_filter(:,1)=lowp(result_pso4,0.5,80,0.1,30,200);
        result_pso5_filter(:,1)=lowp(result_pso5,0.5,80,0.1,30,200);
        figure;
        subplot(2,1,1);
        plot(result_pso1_filter);
        hold on;
        plot(result_pso_filter);
        hold on;
        plot(result_pso3_filter);
        hold on;
        plot(result_pso4_filter);
        hold on;
        plot(result_pso5_filter);
        legend('һ��','����','����','�Ľ�','���','ԭ����');
        title('���׵�PSO�������');
        hold on;
        subplot(2,1,2);
        plot(average_dis_ch1);
        title('ԭʼ����')
    end
end
for i = 1:1     %%�Ż�ǰ��edlen����ֵ�Ա�
    if(edlen_data_compare_improve==1)
        figure('numbertitle','off','name','�Ż�ǰedlen����ֵ/�Ż���edlen����ֵ')
        plot(compensation_edlen_ch1);
        hold on;
        plot(result_edpso_burst_ch1_filter);
        legend('ǰ','��');
    end
end

for i = 1:1     %%����ͨ���Ա�
    if(two_channel_com==1)
        figure('numbertitle','off','name','����ͨ������Ч���Ա�')
        subplot(2,1,1);
        plot(average_dis_ch1);
        hold on;
        plot(result_edpso_ch1);
        hold on;
        legend('ͨ��1����ǰ','ͨ��1������');
        subplot(2,1,2);
        plot(average_dis_ch2);
        hold on;
        plot(result_edpso_ch2);        
        legend('ͨ��2����ǰ','ͨ��2������');       
    end
end

for i = 1:1     %%����ͼ
    market_index = 1:length(average_dis_ch2)/10-1:length(average_dis_ch2);
        if(lunwentu_daiyouhua==1)
        hfigure=figure('numbertitle','off','name','����ͨ������Ч���Ա�');
        subplot(3,1,1,'Position',[0.1 0.74 0.8 0.2]);
        line1=plot(time(1),average_dis_ch2(1),'-o','Color',[0 0.447 0.741]);
        hold on;
        plot(time(ceil(market_index)),average_dis_ch2(ceil(market_index)),'o','Color',[0 0.447 0.741]);
        hold on;
        plot(time,average_dis_ch2,'LineWidth',1.5,'Color',[0 0.447 0.741]);
        hold on;
        
        line2=plot(time(1),result_edlen_ch2_filter(1),'-s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time(ceil(market_index)),result_edlen_ch2_filter(ceil(market_index)),'s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time,result_edlen_ch2_filter,'LineWidth',1.5,'Color',[0.929 0.694 0.125]);
        
        line3=plot(time(1),result_edpso_ch2_filter(1),'-d','Color',[0.85 0.325 0.098]);
        hold on;
        plot(time(ceil(market_index)),result_edpso_ch2_filter(ceil(market_index)),'d','Color',[0.85 0.325 0.098]);
        hold on;
        plot(time,result_edpso_ch2_filter,'LineWidth',1.5,'Color',[0.85 0.325 0.098]);

        pos1=axis;
        xlabel('time/h','FontName','Times New Roman','FontSize',8);
        ylabel('displacement/nm','FontName','Times New Roman','FontSize',8);
        h1=legend([line1 line2 line3],'Raw data','Unoptimized Edlen','Optimized Edlen'); 
        xlim([time(1) 7.522]);
        title('(a) 45mm','FontSize',21,'position',[0.06 20.69  0]);
        set(gca,'FontSize',12) ;
        
        
        
        
        subplot(3,1,2,'Position',[0.1 0.43 0.8 0.2]);
        line21=plot(time(1),average_dis_ch1(1),'-o','Color',[0 0.447 0.741]);
        hold on;
        plot(time(ceil(market_index)),average_dis_ch1(ceil(market_index)),'o','Color',[0 0.447 0.741]);
        hold on;
        plot(time,average_dis_ch1,'LineWidth',1.5,'Color',[0 0.447 0.741]);
        hold on;
        
        line22=plot(time(1),result_edlen_ch1_filter(1),'-s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time(ceil(market_index)),result_edlen_ch1_filter(ceil(market_index)),'s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time,result_edlen_ch1_filter,'LineWidth',1.5,'Color',[0.929 0.694 0.125]);
        
        line23=plot(time(1),result_edpso_ch1_filter(1),'-d','Color',[0.85 0.325 0.098]);
        hold on;
        plot(time(ceil(market_index)),result_edpso_ch1_filter(ceil(market_index)),'d','Color',[0.85 0.325 0.098]);
        hold on;
        plot(time,result_edpso_ch1_filter,'LineWidth',1.5,'Color',[0.85 0.325 0.098]);
        
        pos2=axis;
        xlabel('time/h','FontName','Times New Roman','FontSize',8');
        ylabel('displacement/nm','FontName','Times New Roman','FontSize',8)
        h2=legend([line21 line22 line23],'Raw data','Unoptimized Edlen','Optimized Edlen'); 
%         set(h2,'FontName','Times New Roman','FontSize',8,'FontWeight','normal')
        xlim([time(1) 7.522]);
        title('(b) 90mm','FontSize',21,'position',[0.061 51.626  0]);
        set(gca,'FontSize',12) ;
        
        
        
        
        subplot(3,1,3,'Position',[0.1 0.12 0.8 0.2]);
        AX=plotyy(time,average_tem,time,(average_pre/1000));
        pos3=axis(AX(1));
        pos4=axis(AX(2));
        xlabel('time/h','FontName','Times New Roman','FontSize',8);
        set(get(AX(1),'Ylabel'),'string','temperature/��','FontName','Times New Roman','FontSize',8);
                xlim([time(1) 7.522]);
        set(get(AX(2),'Ylabel'),'string','pressure/kPa','FontName','Times New Roman','FontSize',8);
        set(AX(1),'xlim',[time(1) 7.522]);
        set(AX(2),'xlim',[time(1) 7.522]);
        title('(c) Temperature and Pressure','FontSize',21,'position',[0.194 26.252  0]);
        

        set(AX(1),'FontSize',12) ;
        set(AX(2),'FontSize',12) ;
        set(hfigure,'position',[100 100 800 500]);
        hold on;
        line31=plot(time(1),average_tem(1),'-o','Color',[0 0.447 0.741]);
        hold on;
        plot(time(ceil(market_index)),average_tem(ceil(market_index)),'o','Color',[0 0.447 0.741]);
        hold on;
        line32=plot(time(1),(average_pre(1)/1000),'Color',[0.85 0.325 0.098]);
        hold on;
        line33=plot(time(ceil(market_index)),(average_pre(ceil(market_index))/1000),'o','Color',[0.85 0.325 0.098]);
        hold on;
        h3=legend([line31 line32],'Temperature','Pressure'); 
        set(h3,'FontName','Times New Roman','FontSize',8,'FontWeight','normal');

        
        
        
        
        set(h1,'FontName','Times New Roman','FontSize',8,'FontWeight','normal','Orientation','horizon','position',[0.565 0.91 0.27 0.025])
        set(h2,'FontName','Times New Roman','FontSize',8,'FontWeight','normal','Orientation','horizon','position',[0.565 0.599 0.27 0.025])
        set(h3,'FontName','Times New Roman','FontSize',8,'FontWeight','normal','Orientation','horizon','position',[0.718 0.289 0.147 0.025]);
    end
    
    
    
    
    
    
 if(lunwentu_wuyouhua==1)
        hfigure=figure('numbertitle','off','name','����ͨ������Ч���Ա�');
        subplot(3,1,1,'Position',[0.1 0.74 0.8 0.2]);
        line1=plot(time(1),average_dis_ch2(1),'-o','Color',[0 0.447 0.741]);
        hold on;
        plot(time(ceil(market_index)),average_dis_ch2(ceil(market_index)),'o','Color',[0 0.447 0.741]);
        hold on;
        plot(time,average_dis_ch2,'LineWidth',1.5,'Color',[0 0.447 0.741]);
        hold on;
        
        line2=plot(time(1),result_edlen_ch2_filter(1),'-s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time(ceil(market_index)),result_edlen_ch2_filter(ceil(market_index)),'s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time,result_edlen_ch2_filter,'LineWidth',1.5,'Color',[0.929 0.694 0.125]);
        

        pos1=axis;
        xlabel('time/h','FontName','Times New Roman','FontSize',8);
        ylabel('displacement/nm','FontName','Times New Roman','FontSize',8);
        h1=legend([line1 line2],'Raw data','Unoptimized Edlen'); 
        xlim([time(1) time(length(time))]);
        title('(a) 45mm','FontSize',21,'position',[0.06 20.69  0]);
        set(gca,'FontSize',12) ;
        
        
        
        
        subplot(3,1,2,'Position',[0.1 0.43 0.8 0.2]);
        line21=plot(time(1),average_dis_ch1(1),'-o','Color',[0 0.447 0.741]);
        hold on;
        plot(time(ceil(market_index)),average_dis_ch1(ceil(market_index)),'o','Color',[0 0.447 0.741]);
        hold on;
        plot(time,average_dis_ch1,'LineWidth',1.5,'Color',[0 0.447 0.741]);
        hold on;
        
        line22=plot(time(1),result_edlen_ch1_filter(1),'-s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time(ceil(market_index)),result_edlen_ch1_filter(ceil(market_index)),'s','Color',[0.929 0.694 0.125]);
        hold on;
        plot(time,result_edlen_ch1_filter,'LineWidth',1.5,'Color',[0.929 0.694 0.125]);
        

        
        pos2=axis;
        xlabel('time/h','FontName','Times New Roman','FontSize',8');
        ylabel('displacement/nm','FontName','Times New Roman','FontSize',8)
        h2=legend([line21 line22],'Raw data','Unoptimized Edlen'); 
        xlim([time(1) time(length(time))]);
        title('(b) 90mm','FontSize',21,'position',[0.061 40.626  0]);
        set(gca,'FontSize',12) ;
        
        
        
        
        subplot(3,1,3,'Position',[0.1 0.12 0.8 0.2]);
        AX=plotyy(time,average_tem,time,(average_pre/1000));
        pos3=axis(AX(1));
        pos4=axis(AX(2));
        xlabel('time/h','FontName','Times New Roman','FontSize',8);
        set(get(AX(1),'Ylabel'),'string','temperature/��','FontName','Times New Roman','FontSize',8);
                xlim([time(1) time(length(time))]);
        set(get(AX(2),'Ylabel'),'string','pressure/kPa','FontName','Times New Roman','FontSize',8);
        set(AX(1),'xlim',[time(1) time(length(time))]);
        set(AX(2),'xlim',[time(1) time(length(time))]);
        title('(c) Temperature and Pressure','FontSize',21,'position',[0.194 26.252  0]);
        

        set(AX(1),'FontSize',12) ;
        set(AX(2),'FontSize',12) ;
        set(hfigure,'position',[100 100 800 500]);
        hold on;
        line31=plot(time(1),average_tem(1),'-o','Color',[0 0.447 0.741]);
        hold on;
        plot(time(ceil(market_index)),average_tem(ceil(market_index)),'o','Color',[0 0.447 0.741]);
        hold on;
        line32=plot(time(1),(average_pre(1)/1000),'Color',[0.85 0.325 0.098]);
        hold on;
        line33=plot(time(ceil(market_index)),(average_pre(ceil(market_index))/1000),'o','Color',[0.85 0.325 0.098]);
        hold on;
        h3=legend([line31 line32],'Temperature','Pressure'); 

        
        set(h1,'FontName','Times New Roman','FontSize',8,'FontWeight','normal','Orientation','horizon','position',[0.64118 0.901 0.259 0.041])
        set(h2,'FontName','Times New Roman','FontSize',8,'FontWeight','normal','Orientation','horizon','position',[0.64118 0.5930 0.2587 0.0412])
        set(h3,'FontName','Times New Roman','FontSize',8,'FontWeight','normal','Orientation','horizon','position',[0.683 0.284 0.217 0.041]);
    end

    


%     if(lunwentu_daiyouhua==1)
%         hfigure=figure('numbertitle','off','name','����ͨ������Ч���Ա�');
%         subplot(3,1,1,'Position',[0.1 0.72 0.8 0.22]);
%         plot(time,average_dis_ch2,'LineWidth',2);
%         hold on;
%         plot(time,result_edlen_ch2_filter,'LineWidth',2);
%         hold on;
%         plot(time,result_edpso_ch1_filter,'LineWidth',2)
%         pos1=axis;
%         xlabel('time/h','FontName','Times New Roman','FontSize',18,'position',[1*pos1(2) 1.06*pos1(3)]);
%         ylabel('displacement/nm','FontName','Times New Roman','FontSize',18,'rotation',0,'position',[pos1(1)+1.2 1.1*pos1(4)]);
%         h1=legend('Raw data','Unoptimized Edlen','Optimized Edlen'); 
%         xlim([time(1) time(length(time))]);
%         title('45mm','FontSize',22);
%         set(gca,'FontSize',12) ;
%         
%         
%         
%         
%         subplot(3,1,2,'Position',[0.1 0.39 0.8 0.22]);
%         plot(time,average_dis_ch1,'LineWidth',2);
%         hold on;
%         plot(time,result_edlen_ch1_filter,'LineWidth',2);
%         hold on;
%         plot(time,result_edpso_ch1_filter,'LineWidth',2);
%         pos2=axis;
%         xlabel('time/h','FontName','Times New Roman','FontSize',18,'position',[1*pos2(2) 1.3*pos2(3)]);
%         ylabel('displacement/nm','FontName','Times New Roman','FontSize',18,'rotation',0,'position',[pos2(1)+1.2 2.75*pos2(4)])
%         h2=legend('Raw data','Unoptimized Edlen','Optimized Edlen'); 
% %         set(h2,'FontName','Times New Roman','FontSize',10,'FontWeight','normal')
%         xlim([time(1) time(length(time))]);
%         title('90mm','FontSize',22);
%         set(gca,'FontSize',12) ;
%         
%         
%         
%         
%         subplot(3,1,3,'Position',[0.1 0.06 0.8 0.22]);
%         AX=plotyy(time,average_tem,time,(average_pre/1000));
%         pos3=axis(AX(1));
%         pos4=axis(AX(2));
%         xlabel('time/h','FontName','Times New Roman','FontSize',18,'position',[1*pos3(2) 0.9994*pos3(3)]);
%         set(get(AX(1),'Ylabel'),'string','temperature/��','FontName','Times New Roman','FontSize',19,'rotation',0,'position',[pos3(1)+1.04 1.0001*pos3(4)]);
%                 xlim([time(1) time(length(time))]);
%         set(get(AX(2),'Ylabel'),'string','pressure/kPa','FontName','Times New Roman','FontSize',19,'rotation',0,'position',[pos4(2)*0.90 1.00115*pos4(4)]);
%         set(AX(1),'xlim',[time(1) time(length(time))]);
%         set(AX(2),'xlim',[time(1) time(length(time))]);
%         h3=legend('Temperature','Pressure'); 
%         title('temperature and pressure','FontSize',22);
%         set(h3,'FontName','Times New Roman','FontSize',11,'FontWeight','normal');
%         set(AX(1),'FontSize',12) ;
%         set(AX(2),'FontSize',12) ;
%         set(hfigure,'position',[100 100 800 500]);
%         set(h1,'FontName','Times New Roman','FontSize',10,'FontWeight','normal','position',[0.75 0.85 0.08 0.08])
%         set(h2,'FontName','Times New Roman','FontSize',10,'FontWeight','normal','position',[0.75 0.52 0.08 0.08])
%         set(h3,'FontName','Times New Roman','FontSize',11,'FontWeight','normal','position',[0.71 0.19 0.15 0.09]);
%     end
end


toc



