function [x_vec] = SDP_modified(Tx, Rx, Rg)
%If the convex hull issues arise,
%making the mixed semidefinite/second-order
%cone program infeasible, please proceed
%by omitting the second-order cone constraints.
[~, M] = size(Tx);
[~, L] = size(Rx);


w_mtx = zeros(M,L);

for mm = 1:M
    for ll = 1:L
        w_mtx(mm,ll) = 1 - (Rg(mm,ll)/sum(sum(Rg)));
    end
end


cvx_begin

variables dmt_vec(M) gmt_vec(M) hl_vec(L) c_mtx(M,L) x_vec(2) z_cons uml_mtx(M,L)

expression obj
    
obj = 0;

for m = 1:M
    for l = 1:L
        obj = obj + w_mtx(m,l)*(Rg(m,l)^2 - 2*Rg(m,l)*dmt_vec(m) + gmt_vec(m) - hl_vec(l) - c_mtx(m,l))^2 + 0.01*(c_mtx(m,l)^2 + dmt_vec(m)^2) + 1*uml_mtx(m,l)^2;
    end
end


minimize obj

subject to

[eye(2), x_vec; x_vec', z_cons] == semidefinite(3);

for m = 1:M
    for l = 1:L
        c_mtx(m,l) >= 0;
        Rg(m,l)^2 - 2*Rg(m,l)*dmt_vec(m) + gmt_vec(m) + uml_mtx(m,l) >= hl_vec(l);
        uml_mtx(m,l) >= 0;
    end
end

for m = 1:M
    gmt_vec(m) == [Tx(:,m);-1]'*[eye(2),x_vec;x_vec',z_cons]*[Tx(:,m);-1];
    norm(Tx(:,m) - x_vec) <= dmt_vec(m);
end

for l = 1:L
    hl_vec(l) == [Rx(:,l);-1]'*[eye(2),x_vec;x_vec',z_cons]*[Rx(:,l);-1];
end


cvx_end

end

