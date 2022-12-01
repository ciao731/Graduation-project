clc;
close all;
clear;
t=20:0.05:30;
p=98:0.1:102;
t0 = 20;
p0 = 101.325;

y_p = 1+p/(93214.60*(1+0.0036610*t0))*(1+10^-8*(0.5953-0.009876*t0));
for i = 1:600
    t(i) = 20+(i-1)*0.05;
y_t(i) = 1+p0/93214.60*(1+10^-8*(0.5953-0.009876*t(i)))/(1+0.0036610*t(i));
end
figure
plot(p,y_p)
xlabel('气压/kPa','FontSize',10);
ylabel('折射率','FontSize',10);
set(gcf,'Position',[100,100,600,360])
figure()
plot(t,y_t)
xlabel('温度/°C','FontSize',10);
ylabel('折射率','FontSize',10);
set(gcf,'Position',[100,100,600,360])
