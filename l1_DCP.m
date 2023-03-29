function [y_est_mtx, fail] = l1_DCP(Rx, Rg, dRg, Nmax, epsilon, wBR_mtx, wR_mtx)

%Inputs:
%Rx - receiver position matrix H \times L
%Rg - M \times L bistatic range matrix
%dRg - M \times L direct range matrix
%Nmax - user-specified maximum number of CCCP iterations
%epsilon - threshold to terminate the iterations
%wBR_mtx - weighting factors for indirect path delays
%wR_mtx - weighting factors for direct path delays

%Outputs:
%y_est_mtx - estimates matrix ? \times Ncccp
%fail - indicator for the feasibility of the run

fail = false;

[M, ~] = size(Rg);
[H, L] = size(Rx);

cnt = 0;

y_ini = 1000*(rand(H+H*M+2*M*L,1)-0.5);

y_est_mtx = [];

y_est_mtx = [y_est_mtx,y_ini];

y_kth = y_ini;

%main loop
while cnt <= Nmax
    
cvx_begin quiet

variables y_vec(H+H*M+2*M*L)

expression obj
obj = 0;
for m = 1:M
    for l = 1:L
        obj = obj + wBR_mtx(m,l)*y_vec(H+H*M+(l-1)*M+m) + wR_mtx(m,l)*y_vec(H+H*M+M*L+(l-1)*M+m);
    end
end

minimize obj

subject to

for m = 1:M
    for l = 1:L

        Rg(m,l) - y_vec(H+H*M+(l-1)*M+m) - (norm(y_kth(1:H) - y_kth(H+(m-1)*H+1:H+(m-1)*H+H)) + norm(y_kth(1:H) - Rx(:,l))) - ...
            ((y_kth(1:H) - y_kth(H+(m-1)*H+1:H+(m-1)*H+H))/norm(y_kth(1:H) - y_kth(H+(m-1)*H+1:H+(m-1)*H+H)) + (y_kth(1:H) - Rx(:,l))/norm(y_kth(1:H) - Rx(:,l)))'*(y_vec(1:H) - y_kth(1:H)) - ...
            ((y_kth(H+(m-1)*H+1:H+(m-1)*H+H) - y_kth(1:H))/norm(y_kth(H+(m-1)*H+1:H+(m-1)*H+H) - y_kth(1:H)))'*(y_vec(H+(m-1)*H+1:H+(m-1)*H+H) - y_kth(H+(m-1)*H+1:H+(m-1)*H+H)) <= 0;

        (norm(y_vec(1:H) - y_vec(H+(m-1)*H+1:H+(m-1)*H+H)) + norm(y_vec(1:H) - Rx(:,l))) - (Rg(m,l) + y_kth(H+H*M+(l-1)*M+m)) - (y_vec(H+H*M+(l-1)*M+m) - y_kth(H+H*M+(l-1)*M+m)) <= 0;

        dRg(m,l) - y_vec(H+H*M+M*L+(l-1)*M+m) - norm(y_kth(H+(m-1)*H+1:H+(m-1)*H+H) - Rx(:,l)) - ((y_kth(H+(m-1)*H+1:H+(m-1)*H+H) - Rx(:,l))/norm(y_kth(H+(m-1)*H+1:H+(m-1)*H+H) - Rx(:,l)))'* ...
            (y_vec(H+(m-1)*H+1:H+(m-1)*H+H) - y_kth(H+(m-1)*H+1:H+(m-1)*H+H)) <= 0;

        norm(y_vec(H+(m-1)*H+1:H+(m-1)*H+H) - Rx(:,l)) - (dRg(m,l) + y_kth(H+H*M+M*L+(l-1)*M+m)) - (y_vec(H+H*M+M*L+(l-1)*M+m) - y_kth(H+H*M+M*L+(l-1)*M+m)) <= 0;

    end
end

cvx_end

    y = y_vec;
    
    if (isnan(sum(y))) || (sum(y) == +Inf) || (sum(y) == -Inf)
        fail = true;
        return
    end
    
    y_est_mtx = [y_est_mtx,y];
    
    y_kth = y;
    
    cnt = cnt + 1;
    
    if (norm(y_est_mtx(1:H,end) - y_est_mtx(1:H,end-1))/(min(norm(y_est_mtx(1:H,end)), norm(y_est_mtx(1:H,end-1)))) <= epsilon)
        return
    end
    
end


end

