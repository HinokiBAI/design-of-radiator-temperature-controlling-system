MATLAB 代码1：气气传热

clc;
clear;
% 定义变量
C_w = 4200;             % 定义水的比热容
lou_w = 1000;           % 定义水的密度
V_w = 0.5;              % 定义水流量
m_w = lou_w*V_w;        % 定义水的质量

K_s = 50;               % 暖气片传热系数
F_s = 2;                % 暖气片对流换热面积
T_h = 40;               % 暖气水温

T_out = 0;              % 室外温度
sigma1 = 0.25;          % 墙体厚度
lamda1 = 1.2;           % 墙体传热系数
lamda2 = 5.7;           % 玻璃传热系数
h1 = 0.9;               % 内表面传热系数
h2 = 1.6;               % 外表面传热系数
dx = 0.05;              % 给定差分变量取值范围
Lx = 5;                 % 定义x轴距离
x = 0:dx:Lx;            % 差分密度
Lxx = 5;                % 同上
xx = 0:dx:Lxx;

dz = 0.05;
Lz = 2.8;               % 同上
z = 0:dz:Lz;                
dy = 0.05;
Ly = 10;                % 同上
y = 0:dy:Ly;
time = 360;             % *10即为实际时间
dt = 0.001;
t = 0:dt:time;

Ax = (-2*eye(length(x))+diag(ones(1, length(x)-1), 1)+diag(ones(1, length(x)-1), -1));  % 定义X方向的A矩阵
Ay = (-2*eye(length(y))+diag(ones(1, length(y)-1), 1)+diag(ones(1, length(y)-1), -1));  % 定义Y方向上的A矩阵
Az = (-2*eye(length(z))+diag(ones(1, length(z)-1), 1)+diag(ones(1, length(z)-1), -1));  % 定义Z方向上的A矩阵

[Z, X] = meshgrid(z, x);  
[Y, XX] = meshgrid(y, xx);

U0xoz = ones(length(Ax), length(Az))*0;  % 定义南北墙初始温度为0℃
U0xoy = ones(length(Ax), length(Ay))*0;  % 定义xoy面初始温度0℃

fxozs = ones(length(Ax), length(Az));
fxozs(20:90, 1:25) = 40;                 % 定义热源(仅南墙暖气作用)

fxozn = ones(length(Ax), length(Az))*3;
                % 定义热源(仅北墙暖气作用)

fxoy1 = 5*exp(-1/2*(((XX-2.5).^2)/1+((Y-8.75).^2)/1));
fxoy2 = 5*exp(-1/2*(((XX-2.5).^2)/1+((Y-1.25).^2)/1));

border_n = 8*exp(-1/2*(((z-1).^2)/1));       % 定义边界条件
border_nx = 10*exp(-1/2*(((x-2.5).^2)/1));   % 定义边界条件
Uxozs = U0xoz;
Uxozn = U0xoz;
Uxoy = U0xoy;
a = 0.22;                                    % 定义散热系数

xm = zeros(1, time*10);                      % 控制变量的定义
xm_in_s = zeros(1, time*10);
xm_in_n = zeros(1, time*10);
xm_in_xoy = zeros(1, time*10);

temp = 1;                                    % 控制变量的定义 
temp_in_s = 1;
temp_in_n = 1;
temp_in_xoy = 1;

Q_w = 0;                                     %  墙体导热初值设置
Q_s = 0;                                     % 对流换热初值设置

for n = 1:length(t)-1                        % 迭代时间,计算热传导偏微分方程
    Uxozs = Uxozs+a^2*(1/dx^2*Ax*Uxozs+1/dz^2*Uxozs*Az+fxozs )*dt;   % 南墙的温度分布
    min_xoz_s = mean(Uxozs,2);      		 % 计算南墙平均温度分布,压缩为一维向量作为xoy面的边界条件
    Uxozs(:, 1) = border_nx;             	 % 边界设置
    Uxozs(:, end) = Uxozs(:,end-1);
    Uxozs(1, :) =  border_n;
    Uxozs(end, :) = border_n;
    
    Uxozn = Uxozn+a^2*(1/dx^2*Ax*Uxozn+1/dz^2*Uxozn*Az+fxozn)*dt;   % 北墙的温度分布
    min_xoz_n = mean(Uxozn,2);      		% 计算北墙平均温度分布,压缩为一维向量
    Uxozn(:, 1) = border_nx;                % 边界设置
    Uxozn(:, end) = Uxozn(:,end-1);
    Uxozn(1, :) =  border_n;
    Uxozn(end, :) = border_n;

     % xoy墙的温度分布
    Uxoy = Uxoy+a^2*(1/dx^2*Ax*Uxoy+1/dy^2*Uxoy*Ay+fxoy1+fxoy2)*dt; 
    Uxoy(:, 1) =min_xoz_s;                	% 绝缘设置
    Uxoy(:, end) = min_xoz_n;
    Uxoy(1, :) =  Uxoy(2, :);
    Uxoy(end, :) = Uxoy(end-1, :);
    
    min_in_s = mean(Uxozs(:));     % 南墙的均值
    min_in_n = mean(Uxozn(:));     % 北墙的均值
    min_in_xoy = mean(Uxoy(:));    % xoy的均值
    %---------------回水温度的计算--
    Q_w = ((min_in_s'-T_out)/((1/h1)+(sigma1/lamda1)+(1/h2)))*5*2.8*2*3600;  % 墙体换热
    Q_s = K_s*F_s*(T_h-min_in_s')*3600;   % 对流换热第一个暖气片
    Q_t = Q_w+Q_s;                        % 计算总能量
    T_x = T_h-(Q_t/(C_w*m_w));            % 计算回水温度
    fxozs(20:90, 1:25) = T_x;             % 将回水温度设置为北墙的热源
    xm(1, temp) = T_x;                    % 将回水温度存入向量,之后绘图
    xm_in_s(1, temp_in_s) = min_in_s;     % 保存数据画图
    xm_in_n(1, temp_in_n) = min_in_n;
    xm_in_xoy(1, temp_in_xoy) = min_in_xoy;
    
    if mod(n,1000) == 1                   % 加速画图
        temp = temp + 1;
        temp_in_s = temp_in_s + 1;
        temp_in_n = temp_in_n + 1;
        temp_in_xoy = temp_in_xoy + 1;
        subplot(311);
        surf(X,Z,Uxozs);                  % 绘制南墙温度分布
        shading interp;
        axis([x(1) x(end) z(1) z(end) 0 40]);
        title('南墙温度分布')
        xlabel('x/m');
        ylabel('z/m');
        zlabel('T/℃');
        subplot(312);
        surf(X,Z,Uxozn);                  % 绘制北墙温度分布
        shading interp;
        axis([x(1) x(end) z(1) z(end) 0 40]);
        title('北墙温度分布');
        xlabel('x/m');
        ylabel('z/m');
        zlabel('T/℃');
        subplot(313);
        surf(XX,Y,Uxoy);                  % 绘制xoy墙温度分布
        shading interp;
        axis([xx(1) xx(end) y(1) y(end) 0 40]);
        title('xoy面温度分布');
        xlabel('x/m');
        ylabel('y/m');
        zlabel('T/℃');
%     view(2);  % 改变视角方向
        getframe;
    end
end
% 绘图
x = 1:1:time*10;
figure(4);
plot(x, xm(1, 1:time*10));
legend('出水口温度')
title('暖气回水温度与时间的关系');
xlabel('时间/s');
ylabel('回水温度/℃');
axis([x(1) x(time*10) 0 40]);
figure(5);
plot(x, xm_in_s(1, 1:time*10));
hold on
plot(x, xm_in_n(1, 1:time*10));
hold on
plot(x, xm_in_xoy(1, 1:time*10));
legend('南墙平均温度', '北墙平均温度',  '底面xoy面平均温度')
axis([x(1) x(time*10) 0 30]);
title('各个面温度随时间变化曲线');
xlabel('时间/s');
ylabel('温度/℃');

**************************************************************************************************************************************
**************************************************************************************************************************************
**************************************************************************************************************************************

Python 代码2：NTU

from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

c = 4.2*10**3
S = 2
U = 2.7
Tci = 0
def NTU(G, T):
    Cmin = G*c
    NTU = S*U/Cmin
    E = (1 - np.exp(-2*NTU))/2
    qmax = Cmin*(T - Tci)
    return E*qmax

q_list = []
G_list = np.array(np.linspace(0.1, 4, 1000))
T_list = np.array(np.linspace(40, 80, 1000))


X,Y = np.meshgrid(G_list, T_list)
Z = NTU(X, Y)

# 绘制X与Z的曲线
fig = plt.figure()

ax = plt.axes()
ax.plot(G_list, Z[0])
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
ax.set_xlabel('流量(Kg/h)')
ax.set_ylabel('热传率(W)')
# 设置DPI为300
fig.set_dpi(300)
# 设置图片大小
fig.set_size_inches(8.5, 5.5)
# 设置标题
plt.title('流量与热传率的关系图')
plt.savefig('G_q.png', dpi=1000)

# 绘制Y与Z的曲线
fig = plt.figure()
ax = plt.axes()
ax.plot(T_list, Z[:,0])
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
ax.set_xlabel('温度(℃)')
ax.set_ylabel('热传率(W)')
# 设置DPI为300
fig.set_dpi(300)
# 设置图片大小
fig.set_size_inches(8.5, 5.5)
# 设置标题
plt.title('温度与热传率的关系图')
plt.savefig('T_q.png', dpi=1000)

plt.show()


