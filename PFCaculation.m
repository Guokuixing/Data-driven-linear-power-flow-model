% 本算法为毕业设计《数据驱动的线性化潮流方程建模方法及开源实现》的研究内容，作者郭魁星，指导教师康重庆教授、张宁副教授。
% 使用Caculation.m文件可以实现对IEEE任意节点系统进行潮流分布计算，输入case及节点电压幅值和相角数据即可获得系统潮流分布。
% 经本研究验证，该算法在IEEE-5、30、33bw、57、118节点系统中均实现了较高的计算精度和较快的计算速度。

%% function[data_out]=...
    PFCaculation(case_name, V_in, Va_in, bus_branch)
% bus_branch: 0,caculate pf and qf;1,caculate p and q

%% define parameters
generate_data = 0;%1,data generation is needed; 0,data is already generated
upper_bound = 1.20;%upper bound of generated load
lower_bound = 0.80;%lower bound of generated load

G_range = 0.1; %range of power generation variations
Q_range = 0.25; %range of Q variations
Q_per = 0.2; %Q percentage on P
V_range = 0.01; %range of voltage magnitude variations of PV buses
L_range = 0.05; %range of load in different nodes
L_corr = 0.9; %covariance
Va_range = 7;%degree
Va_num = [];
dc_ac = 1; %0-dc;1-ac;
random_load = 1; %1,random 0,not random with bounder 2,not random with covariance

data_size = 500;% training data size
address = '';% address to read and save the data filess

%%  training data generation
data_name = [address case_name '_training_data'];
if (generate_data)
    mpc = ext2int(loadcase(case_name));
    [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
    DataGeneration(case_name, Q_per, data_name, dc_ac, G_range, ...
        upper_bound, lower_bound, Q_range, V_range, data_size, L_range, ...
        random_load, Va_range, ref, L_corr);      
end
load([data_name,'.mat']);

%%  linear regression
%  get bus index lists of each type of bus
    mpc = ext2int(loadcase(case_name));
    [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);

    [Xp, Xq, Xpf, Xqf, Xpt, Xqt, yanzheng] =...
        RegressionForward(num_load, num_branch, data, address, case_name);

%% caculation
size_in = size(V_in,1);
for i =1:size_in
    if(bus_branch)
        data.p(i, :) = [Va_in(i,:) * pi / 180 V_in(i,:).^2 1] * Xp';
        data.q(i, :) = [Va_in(i,:) * pi / 180 V_in(i,:).^2 1]*Xq';
    else
        data.pf(i, :) = [Va_in(i,:) * pi / 180 V_in(i,:).^2 1] * Xpf';
        data.qf(i, :) = [Va_in(i,:) * pi / 180 V_in(i,:).^2 1] * Xqf';
    end
end
    