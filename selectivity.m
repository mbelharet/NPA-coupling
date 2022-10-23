function s =  selectivity( L, L_prey,a1, a2, rho1, rho2)
%     a1 = param(1); %5;
%     a2 = param(2); %0.05;
%     rho1 = param(3); %3;
%     rho2 = param(4); %100;

    for i = 1:length(L)
        for j = 1:length(L_prey)
            f = 1./(1 + exp(a1*(rho1 - L(i)/L_prey(j))));
            g = 1 - 1/(1 + exp(a2*(rho2 - L(i)/L_prey(j))));
            s(i,j) = f * g;
        end
    end
    
end