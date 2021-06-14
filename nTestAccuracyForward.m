function [delta, ntest, XXp, XXq] = ...
    nTestAccuracyForward(num_branch, num_load, data, Xp, Xq,Xpf, Xqf,Xpt, Xqt, Xp_dlpf, Xq_dlpf, B,case_name, ntest)
% this function test the accuracy of forward regression
%% calculate the results by data-driven linearized equations
    A1 = Xpf;
    A2 = Xqf;
    A3 = Xpt;
    A4 = Xqt;
for m = 1:num_branch
    XXp = zeros(num_load,2*num_load+1);
    XXq = zeros(num_load,2*num_load+1);
    mpc = ext2int(loadcase(case_name));
    Connect = mpc.branch(:,1:2);
    Xpf = A1;
    Xqf = A2;
    Xpt = A3;
    Xqt = A4;
    Xpf(m,:) = 0;
    Xqf(m,:) = 0;
    Xpt(m,:) = 0;
    Xqt(m,:) = 0;
    for i = 1:num_load
        id_f = find(Connect(:,1) == i);
        id_t = find(Connect(:,2) == i);
        num_idf = size(id_f,1);
        num_idt = size(id_t,1);
        for k =1:num_idf
            j = Connect(id_f(k,:),2);
            XXp(i,i) = XXp(i,i)+Xpf(id_f(k,:),i);
            XXp(i,j) = XXp(i,j)+Xpf(id_f(k,:),j);
            XXp(i,num_load+i) = XXp(i,num_load+i)+Xpf(id_f(k,:),num_load+i);
            XXp(i,num_load+j) = XXp(i,num_load+j)+Xpf(id_f(k,:),num_load+j);
            XXp(i,2*num_load+1) = XXp(i,2*num_load+1)+Xpf(id_f(k,:),2*num_load+1);
            XXq(i,i) = XXq(i,i)+Xqf(id_f(k,:),i);
            XXq(i,j) = XXq(i,j)+Xqf(id_f(k,:),j);
            XXq(i,num_load+i) = XXq(i,num_load+i)+Xqf(id_f(k,:),num_load+i);
            XXq(i,num_load+j) = XXq(i,num_load+j)+Xqf(id_f(k,:),num_load+j);
            XXq(i,2*num_load+1) = XXq(i,2*num_load+1)+Xqf(id_f(k,:),2*num_load+1);
        end
        for k =1:num_idt
            j = Connect(id_t(k,:),1);
            XXp(i,i) = XXp(i,i)+Xpt(id_t(k,:),i);
            XXp(i,j) = XXp(i,j)+Xpt(id_t(k,:),j);
            XXp(i,num_load+i) = XXp(i,num_load+i)+Xpt(id_t(k,:),num_load+i);
            XXp(i,num_load+j) = XXp(i,num_load+j)+Xpt(id_t(k,:),num_load+j);
            XXp(i,2*num_load+1) = XXp(i,2*num_load+1)+Xpt(id_t(k,:),2*num_load+1);
            XXq(i,i) = XXq(i,i)+Xqt(id_t(k,:),i);
            XXq(i,j) = XXq(i,j)+Xqt(id_t(k,:),j);
            XXq(i,num_load+i) = XXq(i,num_load+i)+Xqt(id_t(k,:),num_load+i);
            XXq(i,num_load+j) = XXq(i,num_load+j)+Xqt(id_t(k,:),num_load+j);
            XXq(i,2*num_load+1) = XXq(i,2*num_load+1)+Xqt(id_t(k,:),2*num_load+1);
        end
    end
    ntest.p.fitting(m, :) = [data.Va(m,:) * pi / 180 data.V(m,:).^2 1] * XXp';
    ntest.q.fitting(m, :) = [data.Va(m,:) * pi / 180 data.V(m,:).^2 1] * XXq';
%     ntest.p.dcpf(m, :) = B * data.Va(m, :)' * pi / 180;
%     ntest.p.dlpf(m, :) = [data.Va(m,:) * pi / 180 data.V(m,:)]*Xp_dlpf';
%     ntest.q.dlpf(m, :) = [data.Va(m,:) * pi / 180 data.V(m,:)]*Xq_dlpf';
end
%% calculate the errors, note that the value of nan or inf is removed
temp = abs((data.P - ntest.p.fitting)./data.P);
delta.p.fitting = temp * 100;
temp(find(isnan(temp)==1)) = [];
temp(find(isinf(temp)==1)) = [];
delta.p.find = mean(mean(temp)) * 100;

temp = abs((data.Q - ntest.q.fitting)./data.Q);
delta.q.fitting = temp* 100;
temp(find(isnan(temp)==1)) = [];
temp(find(isinf(temp)==1)) = [];
delta.q.find = mean(mean(temp)) * 100;

temp = abs((data.P - ntest.p.dcpf)./data.P);
ntest.dcpf.p = temp;
temp(find(isnan(temp)==1)) = [];
temp(find(isinf(temp)==1)) = [];
delta.p.dcpf = mean(mean(temp)) * 100;

temp = abs((data.P - ntest.p.dlpf)./data.P);
ntest.dlpf.p = temp;
temp(find(isnan(temp)==1)) = [];
temp(find(isinf(temp)==1)) = [];
delta.p.dlpf = mean(mean(temp)) * 100;

temp = abs((data.Q - ntest.q.dlpf)./data.Q);
ntest.dlpf.q = temp;
temp(find(isnan(temp)==1)) = [];
temp(find(isinf(temp)==1)) = [];
delta.q.dlpf = mean(mean(temp)) * 100;
