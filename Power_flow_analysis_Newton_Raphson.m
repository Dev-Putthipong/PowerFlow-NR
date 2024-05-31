% MATLAB code for power flow analysis by using Newton-Raphson method
% Created by Putthipong Niyomkitjakankul
% May 2024
% Disclaimer: This MATLAB Code was created out of the personal interest and 
%             is intended for educational purposes only. Please note that it may contain error.

clearvars -except BranchData BusData
BusData= sortrows(BusData);
N = height(BusData);
S_base = 100; % Adjustable
Y = zeros(N);

for k=1:height(BranchData)
    a = BranchData.Tap(k);
    shift = deg2rad(BranchData.Shift(k));
    [ax, ay] = pol2cart(shift, a);
    Y(BranchData.bus_i(k), BranchData.bus_j(k)) = Y(BranchData.bus_i(k), BranchData.bus_j(k)) - 1/(conj(ax+1i*ay)*(BranchData.R(k)+1i*BranchData.X(k)));
    Y(BranchData.bus_j(k), BranchData.bus_i(k)) = Y(BranchData.bus_j(k), BranchData.bus_i(k)) - 1/((ax+1i*ay)*(BranchData.R(k)+1i*BranchData.X(k)));
    Y(BranchData.bus_i(k), BranchData.bus_i(k)) = Y(BranchData.bus_i(k), BranchData.bus_i(k)) + 1/(a^2*(BranchData.R(k)+1i*BranchData.X(k))) + 1i*BranchData.half_B(k);
    Y(BranchData.bus_j(k), BranchData.bus_j(k)) = Y(BranchData.bus_j(k), BranchData.bus_j(k)) + 1/(BranchData.R(k)+1i*BranchData.X(k)) + 1i*BranchData.half_B(k);
end

for k=1:N
    Y(k, k) = Y(k, k) + BusData.G(k) + 1i*BusData.B(k);
end

upper_limit_violation = false(N, 1);
lower_limit_violation = false(N, 1);
upper_pv2pq = false(N, 1);
lower_pv2pq = false(N, 1);
BusData_copy = BusData;
while 1
    V = BusData_copy.V(:);
    Phase = BusData_copy.Phase(:);
    
    V(isnan(V)) = 1; % Adjustable
    Phase(isnan(Phase)) = 0; % Adjustable
    
    Log_load = BusData_copy.Type == "Load" | BusData_copy.Type == "PQ";
    Log_pv_load = BusData_copy.Type == "PV" | Log_load;
    
    P_variable = BusData_copy.Bus(Log_pv_load);
    Delta_variable = P_variable;
    Q_variable = BusData_copy.Bus(Log_load);
    V_variable = Q_variable;
    
    P_sch = (BusData_copy.Pgen(Log_pv_load) - BusData_copy.Pload(Log_pv_load)) / S_base;
    Q_sch = (BusData_copy.Qgen(Log_load) - BusData_copy.Qload(Log_load)) / S_base;
    
    tolerance = 1e-10; % Adjustable
    iteration = 0;
    
    while 1
        
        iteration = iteration + 1;
        J1_full = zeros(N);
        J2_full = zeros(N);
        J3_full = zeros(N);
        J4_full = zeros(N);
        for i = 1:N
            x=0;
            y=0;
            z=0;
            w=0;
            for j = 1:N
                if i ~= j
                    J1_full(i, j) = V(i)*V(j)*abs(Y(i, j))*sin(Phase(i) - Phase(j) - angle(Y(i, j)));
                    w = w + J1_full(i, j);
    
                    J2_full(i, j) = V(i)*abs(Y(i, j))*cos(Phase(i) - Phase(j) - angle(Y(i, j)));
                    x = x + (V(j)/V(i))*J2_full(i, j);
    
                    J3_full(i, j) = -V(i)*V(j)*abs(Y(i, j))*cos(Phase(i) - Phase(j) - angle(Y(i, j)));
                    y = y - J3_full(i, j);
    
                    J4_full(i, j) = V(i)*abs(Y(i, j))*sin(Phase(i) - Phase(j) - angle(Y(i, j)));
                    z = z + (V(j)/V(i))*J4_full(i, j);
                end
            end      
            J1_full(i, i) = -w;
            J2_full(i, i) = 2*V(i)*abs(Y(i, i))*cos(angle(Y(i, i))) + x;
            J3_full(i, i) = y;
            J4_full(i, i) = -2*V(i)*abs(Y(i, i))*sin(angle(Y(i, i))) + z;
        end
        
        P_size = length(P_variable);
        Delta_size = P_size;
        Q_size = length(Q_variable);
        V_size = Q_size;
    
        J1 = zeros(P_size, Delta_size);
        J2 = zeros(P_size, V_size);
        P_cal = [];
        for i = 1:P_size
            k = P_variable(i);
            for j = 1:Delta_size
                J1(i, j) = J1_full(k, Delta_variable(j));
            end        
            for j = 1:V_size
                J2(i, j) = J2_full(k, V_variable(j));
            end        
            pcal = 0;
            for j = 1:N
                pcal = pcal + V(k)*V(j)*abs(Y(k,j))*cos(Phase(k) - Phase(j) - angle(Y(k,j)));
            end
            P_cal = [P_cal; pcal];
        end
    
        J3 = zeros(Q_size, Delta_size);
        J4 = zeros(Q_size, V_size);
        Q_cal = [];
        for i = 1:Q_size
            k = Q_variable(i);
            for j = 1:Delta_size
                J3(i, j) = J3_full(Q_variable(i), Delta_variable(j));
            end
            for j = 1:V_size
                J4(i, j) = J4_full(Q_variable(i), V_variable(j));
            end
            qcal = 0;
            for j = 1:N
                qcal = qcal + V(k)*V(j)*abs(Y(k,j))*sin(Phase(k) - Phase(j) - angle(Y(k,j)));
            end
            Q_cal = [Q_cal; qcal];
        end
    
        Full_Jacobian = [J1_full J2_full; J3_full J4_full];
        Jacobian = [J1 J2; J3 J4];
    
        residual = [P_sch-P_cal; Q_sch-Q_cal];
        mismatch = Jacobian\residual;
        
        Phase(Delta_variable) = Phase(Delta_variable) + mismatch(1:Delta_size);
        V(V_variable) = V(V_variable) + mismatch(Delta_size+1:Delta_size+V_size);
        
        Newton(iteration) = struct('Jacobian', Jacobian, 'Pcal', P_cal, 'Qcal', Q_cal, 'Residual', residual, 'Mismatch', mismatch, 'V', V, 'Phase', Phase);
    
        if norm(mismatch) < tolerance
            break;
        end
        
    end
    
    I_line_flow = table();
    S_line_flow = table();
    S_network = zeros(N, 1);
    for k= 1:height(BranchData)
        i = BranchData.bus_i(k);
        j = BranchData.bus_j(k);
        a = BranchData.Tap(k);
        shift = deg2rad(BranchData.Shift(k));
        [ax, ay] = pol2cart(shift, a);
        [Vi_x, Vi_y] = pol2cart(Phase(i), V(i));
        [Vj_x, Vj_y] = pol2cart(Phase(j), V(j));
    
        R = BranchData.R(k);
        X = BranchData.X(k);
        y = 1/(R + 1i*X);
        
        if shift ~= 0
            Iij = -(y/conj(ax + 1i*ay))*(Vj_x+1i*Vj_y) + (y/a^2)*(Vi_x+1i*Vi_y);
            Iji = y*(Vj_x+1i*Vj_y) - (y/(ax + 1i*ay))*(Vi_x+1i*Vi_y);
            Sij = (Vi_x+1i*Vi_y) * conj(Iij)*S_base;
            Sji = (Vj_x+1i*Vj_y) * conj(Iji)*S_base;
            Ii0 = NaN;
            Ij0 = NaN;
            Iij_line = NaN;
            Iji_line = NaN;
            Si0 = NaN;
            Sj0 = NaN;
            Sij_line = NaN;
            Sji_line = NaN;
        else
            Mid = y/a;
            Half_i = 1i*BranchData.half_B(k) + y*(1-a)/a^2;
            Half_j = 1i*BranchData.half_B(k) + y*(a-1)/a;
        
            Iij_line = ((Vi_x+1i*Vi_y) - (Vj_x+1i*Vj_y)) * Mid;
            Ii0 = Half_i*(Vi_x + 1i*Vi_y);
            Iij = Iij_line + Ii0;
            
            Sij_line = (Vi_x+1i*Vi_y)*conj(Iij_line)*S_base;
            Si0 = (Vi_x+1i*Vi_y)*conj(Ii0)*S_base;
            Sij = Sij_line + Si0;
    
            Iji_line = ( (Vj_x+1i*Vj_y)-(Vi_x+1i*Vi_y) ) * Mid ;
            Ij0 = Half_j*(Vj_x + 1i*Vj_y);
            Iji = Iji_line + Ij0;
        
            Sji_line = (Vj_x + 1i*Vj_y)*conj(Iji_line)*S_base;
            Sj0 = (Vj_x + 1i*Vj_y)*conj(Ij0)*S_base;
            Sji = Sji_line + Sj0;
        end
        S_network(i) = S_network(i) + Sij;
        S_network(j) = S_network(j) + Sji;
    
        Total_power2ground = Si0 + Sj0;
        RX_loss = Sij_line + Sji_line;
        Total_loss = Sij + Sji;
    
        I_line_flow = [I_line_flow; {i j Iij_line Ii0 Iij Iji_line Ij0 Iji}];
        S_line_flow = [S_line_flow; {i j Sij_line Si0 Sij Sji_line Sj0 Sji Total_power2ground RX_loss Total_loss}];
    end
    
    I_line_flow.Properties.VariableNames = {'i', 'j' 'Iij_line', 'Ii0', 'Iij', 'Iji_line', 'Ij0', 'Iji'};
    S_line_flow.Properties.VariableNames = {'i', 'j' 'Sij_line', 'Si0', 'Sij', 'Sji_line', 'Sj0', 'Sji', 'Si0+Sj0', 'Sij_line+Sji_line', 'Total_loss'};
    
    Summary_line_flow = table(S_line_flow.i, S_line_flow.j, real(S_line_flow.Sij), imag(S_line_flow.Sij), real(S_line_flow.Sji), imag(S_line_flow.Sji), real(S_line_flow.Total_loss), imag(S_line_flow.Total_loss));
    Summary_line_flow.Properties.VariableNames = {'i', 'j', 'Pij', 'Qij', 'Pji', 'Qji', 'P_loss', 'Q_loss'};
    
    S_line_flow = [S_line_flow; {NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, sum(S_line_flow.Si0+Sj0(:)), sum(S_line_flow.Sij_line+Sji_line(:)), sum(S_line_flow.Total_loss(:))}];
    Summary_line_flow = [Summary_line_flow; {NaN,NaN,NaN,NaN,NaN,NaN, sum(Summary_line_flow.P_loss(:)), sum(Summary_line_flow.Q_loss(:))}];
    
    Q_injected = zeros(N, 1);
    P_bus_loss = zeros(N,1);
    for k=1:N
        Q_injected(k) = S_base*BusData_copy.B(k)*V(k)^2;
        P_bus_loss(k) = S_base*BusData_copy.G(k)*V(k)^2;
    end
    
    P_load = BusData_copy.Pload(:);
    Q_load = BusData_copy.Qload(:);
    P_network = real(S_network);
    Q_network = imag(S_network);
    P_gen = P_network + P_load + P_bus_loss;
    Q_gen = Q_network + Q_load - Q_injected;
    
    Phase = rad2deg(Phase);
    No_Bus = BusData_copy.Bus(:);
    Bus_information = table(No_Bus, V, Phase, P_gen, Q_gen, P_load, Q_load, Q_injected, P_bus_loss);
    
    Bus_information = [Bus_information; {NaN,NaN,NaN, sum(P_gen), sum(Q_gen), sum(P_load), sum(Q_load), sum(Q_injected), sum(P_bus_loss)}];
    
    disp(Summary_line_flow);
    disp(Bus_information);
    
    check_upper_voltage = (V > BusData.V) & (upper_pv2pq);
    check_lower_voltage = (V < BusData.V) & (lower_pv2pq);
    pass_upper_voltage = (V < BusData.V) & (upper_pv2pq);
    pass_lower_voltage = (V > BusData.V) & (lower_pv2pq);
    
    upper_limit_violation = round(Q_gen, 6) > BusData.Qmax;
    lower_limit_violation = round(Q_gen, 6) < BusData.Qmin;
    upper_limit_violation(BusData.Type ~= "PV") = false;
    lower_limit_violation(BusData.Type ~= "PV") = false;

    upper_pv2pq = pass_upper_voltage | upper_limit_violation;
    lower_pv2pq = pass_lower_voltage | lower_limit_violation;
    
    if any(upper_limit_violation) | any(lower_limit_violation) | any(check_lower_voltage) | any(check_upper_voltage)
        if any(upper_limit_violation) | any(lower_limit_violation) 
            disp('Q_violation');
            fprintf('upper Q violation at bus : %d\n', BusData.Bus(upper_limit_violation));
            fprintf('lower Q violation at bus : %d\n', BusData.Bus(lower_limit_violation));
            BusData_copy.Qgen(upper_limit_violation) = BusData.Qmax(upper_limit_violation);
            BusData_copy.Qgen(lower_limit_violation) = BusData.Qmin(lower_limit_violation);
            BusData_copy.Type(upper_limit_violation | lower_limit_violation) = "PQ";
            BusData_copy.V(upper_limit_violation | lower_limit_violation) = NaN;
        end
        if any(check_lower_voltage) | any(check_upper_voltage)
            disp('Out of voltage range');
            fprintf('out of voltage range (upper limit) at bus : %d\n', BusData.Bus(check_upper_voltage));
            fprintf('out of voltage range (lower limit) at bus : %d\n', BusData.Bus(check_lower_voltage));
            BusData_copy.V(check_upper_voltage) = BusData.V(check_upper_voltage);
            BusData_copy.V(check_lower_voltage) = BusData.V(check_lower_voltage);
            BusData_copy.Type(check_upper_voltage | check_lower_voltage)= "PV";
            BusData_copy.Qgen(check_upper_voltage | check_lower_voltage)= NaN;
        end
        if input("(1) Fix (2) Don't fix : ") == 2
            break;
        end
    else
        disp("The power flow analysis is completed");
        break;
    end
end

clearvars -except BranchData BusData Bus_information I_line_flow S_line_flow Summary_line_flow Newton Y P_sch Q_sch
