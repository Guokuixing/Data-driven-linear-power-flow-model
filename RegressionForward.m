function [Xp, Xq, Xpf, Xqf, Xpt, Xqt,yanzheng] =...
    RegressionForward(num_load, num_branch, data, address, case_name)
% this function conduct the forward regression
    % ordinary least squares
        Xp = zeros(num_load,2*num_load+1);
        Xq = zeros(num_load,2*num_load+1);
        Xpf = zeros(num_branch,2*num_load+1);
        Xqf = zeros(num_branch,2*num_load+1);
        Xpt = zeros(num_branch,2*num_load+1);
        Xqt = zeros(num_branch,2*num_load+1);
        mpc = ext2int(loadcase(case_name));
        Connect = mpc.branch(:,1:2);
        for i = 1:num_load
            id_f = find(Connect(:,1) == i);
            id_t = find(Connect(:,2) == i);
            num_idf = size(id_f,1);
            num_idt = size(id_t,1);
            for k =1:num_idf
                pf = [data.PF(:,id_f(k,:))/mpc.baseMVA];
                j = Connect(id_f(k,:),2);
                V_Va_pf = [data.Va(:,i) * pi / 180 data.Va(:,j) * pi / 180 data.V(:,i).^2 data.V(:,j).^2 ones(size(data.V, 1), 1)];
                b = regress(pf, V_Va_pf);
                Xpf(id_f(k,:),i) = Xpf(id_f(k,:),i)+b(1);
                Xpf(id_f(k,:),j) = Xpf(id_f(k,:),j)+b(2);
                Xpf(id_f(k,:),num_load+i) = Xpf(id_f(k,:),num_load+i)+b(3);
                Xpf(id_f(k,:),num_load+j) = Xpf(id_f(k,:),num_load+j)+b(4);
                Xpf(id_f(k,:),2*num_load+1) = Xpf(id_f(k,:),2*num_load+1)+b(5);
                Xp(i,i) = Xp(i,i)+b(1);
                Xp(i,j) = Xp(i,j)+b(2);
                Xp(i,num_load+i) = Xp(i,num_load+i)+b(3);
                Xp(i,num_load+j) = Xp(i,num_load+j)+b(4);
                Xp(i,2*num_load+1) = Xp(i,2*num_load+1)+b(5);
                qf = [data.QF(:,id_f(k,:))/mpc.baseMVA];
                b = regress(qf, V_Va_pf);
                Xqf(id_f(k,:),i) = Xqf(id_f(k,:),i)+b(1);
                Xqf(id_f(k,:),j) = Xqf(id_f(k,:),j)+b(2);
                Xqf(id_f(k,:),num_load+i) = Xqf(id_f(k,:),num_load+i)+b(3);
                Xqf(id_f(k,:),num_load+j) = Xqf(id_f(k,:),num_load+j)+b(4);
                Xqf(id_f(k,:),2*num_load+1) = Xqf(id_f(k,:),2*num_load+1)+b(5);
                Xq(i,i) = Xq(i,i)+b(1);
                Xq(i,j) = Xq(i,j)+b(2);
                Xq(i,num_load+i) = Xq(i,num_load+i)+b(3);
                Xq(i,num_load+j) = Xq(i,num_load+j)+b(4);
                Xq(i,2*num_load+1) = Xq(i,2*num_load+1)+b(5);
            end
            for k =1:num_idt
                pf = [data.PT(:,id_t(k,:))/mpc.baseMVA];
                j = Connect(id_t(k,:),1);
                V_Va_pf = [data.Va(:,i) * pi / 180 data.Va(:,j) * pi / 180 data.V(:,i).^2 data.V(:,j).^2 ones(size(data.V, 1), 1)];
                b = regress(pf, V_Va_pf);
                Xpt(id_t(k,:),i) = Xpt(id_t(k,:),i)+b(1);
                Xpt(id_t(k,:),j) = Xpt(id_t(k,:),j)+b(2);
                Xpt(id_t(k,:),num_load+i) = Xpt(id_t(k,:),num_load+i)+b(3);
                Xpt(id_t(k,:),num_load+j) = Xpt(id_t(k,:),num_load+j)+b(4);
                Xpt(id_t(k,:),2*num_load+1) = Xpt(id_t(k,:),2*num_load+1)+b(5);
                Xp(i,i) = Xp(i,i)+b(1);
                Xp(i,j) = Xp(i,j)+b(2);
                Xp(i,num_load+i) = Xp(i,num_load+i)+b(3);
                Xp(i,num_load+j) = Xp(i,num_load+j)+b(4);
                Xp(i,2*num_load+1) = Xp(i,2*num_load+1)+b(5);
                qf = [data.QT(:,id_t(k,:))/mpc.baseMVA];
                b = regress(qf, V_Va_pf);
                Xqt(id_t(k,:),i) = Xqt(id_t(k,:),i)+b(1);
                Xqt(id_t(k,:),j) = Xqt(id_t(k,:),j)+b(2);
                Xqt(id_t(k,:),num_load+i) = Xqt(id_t(k,:),num_load+i)+b(3);
                Xqt(id_t(k,:),num_load+j) = Xqt(id_t(k,:),num_load+j)+b(4);
                Xqt(id_t(k,:),2*num_load+1) = Xqt(id_t(k,:),2*num_load+1)+b(5);
                Xq(i,i) = Xq(i,i)+b(1);
                Xq(i,j) = Xq(i,j)+b(2);
                Xq(i,num_load+i) = Xq(i,num_load+i)+b(3);
                Xq(i,num_load+j) = Xq(i,num_load+j)+b(4);
                Xq(i,2*num_load+1) = Xq(i,2*num_load+1)+b(5);
            end            
        end
        for i = 1:500
            yanzheng.p(i,:) = [data.Va(i,:) * pi / 180 data.V(i,:).^2 1] * Xp';
            yanzheng.q(i,:) = [data.Va(i,:) * pi / 180 data.V(i,:).^2 1]*Xq';
        end
