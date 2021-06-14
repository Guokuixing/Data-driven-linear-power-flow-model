% 本算法为毕业设计《数据驱动的线性化潮流方程建模方法及开源实现》的研究内容，作者郭魁星，指导教师康重庆教授、张宁副教授。
% StartUp.m文件为线性化潮流方程建模并应用于N-1准则检验的算法实现过程。
clc;
clear;
%% define parameters
generate_data = 0;%1,data generation is needed; 0,data is already generated
generate_test_data = 0;%1,data generation is needed; 0,data is already generated
generate_test_n = 1;%1,data generation is needed; 0,data is already generated
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
data_size_test = 300;% testing data size
case_name = 'case5';
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

[Xp_dlpf, Xq_dlpf,~, ~, ~] = DLPF(mpc);
Xp_dlpf = full(Xp_dlpf);
Xq_dlpf = full(Xq_dlpf);

    [Xp, Xq, Xpf, Xqf, Xpt, Xqt, yanzheng] =...
        RegressionForward(num_load, num_branch, data, address, case_name);

%% generate testing data
upper_bound = 1.2;
lower_bound = 0.8;
data_name = [address case_name '_testing_data'];
if (generate_test_data)
    DataGeneration(case_name, Q_per, data_name, dc_ac, G_range,...
        upper_bound, lower_bound, Q_range, V_range, data_size_test, L_range, ...
        random_load, Va_range, ref, L_corr); 
end
load([data_name,'.mat']);
num_train = size(data.P, 1);

%% verify the accuracy

    [delta, test] = ...
        TestAccuracyForward(num_train, data, Xp, Xq, Xp_dlpf, Xq_dlpf, B);

%% n-1principle test
data_name = [address case_name '_testing_n_data'];
if (generate_test_n)
    mpc = ext2int(loadcase(case_name));
    [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
    [ntest] = nDataGeneration(case_name, Q_per, data_name, dc_ac, G_range, ...
        upper_bound, lower_bound, Q_range, V_range, L_range, ...
        random_load, Va_range, ref, L_corr,test); 
    load([data_name,'.mat']);
end
%% verify the n-1principle accuracy
if (generate_test_n)
    [delta, ntest, XXp, XXq] = ...
        nTestAccuracyForward(num_branch, num_load, data, Xp, Xq, Xpf, Xqf, Xpt, Xqt, Xp_dlpf, Xq_dlpf, B, case_name, ntest);
end
