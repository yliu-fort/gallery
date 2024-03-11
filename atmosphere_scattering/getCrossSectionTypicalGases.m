function total_sigma_s = getCrossSectionTypicalGases(gas_fractions, lambdas)
% O2,N2,Ar,Ne,H2,CO2,CH4
% e.g. for earth, gas_fractions = [0.26,0.74,0,0,0,0,0];
% total_sigma_s (nl, 1)
Ns = 6.02e23; % Avogadro's number

gas_name = ["O2","N2","Ar","Ne","H2","CO2","CH4"];

measured_lambdas = [632.8 514.5 488.0 457.9 363.8];

unit_coeff = (1e-31*Ns)*((measured_lambdas*1e-3).^4);

% from NASA report
measured_sigma_s = unit_coeff .*[...
2.06 4.88 6.50 8.39 20.03 ;...
2.24 5.61 7.26 10.38 23.82 ;...
2.08 5.46 7.24 10.13 23.00 ;...
0.103 0.25 0.33 0.42 1.01 ;...
0.036 0.086 0.115 0.15 0.35 ;...
0.493 1.17 1.56 2.01 4.80 ;...
7.28 17.25 23.00 29.60 70.70 ;...
5.26 12.44 16.59 21.40 51.10 ];

gas_fractions = gas_fractions / sum(gas_fractions);

total_sigma_s = 0*lambdas;
for l = 1:numel(lambdas)
    total_sigma_s(l) = 0;
    for c = 1:numel(gas_fractions)
        total_sigma_s(l) = total_sigma_s(l) ...
            + gas_fractions(c) ...
            .* max(0, interp1(measured_lambdas, measured_sigma_s(c,:),lambdas(l),'nearest','extrap'));
    end
end


end