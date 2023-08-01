% To determine Q1 (upper)bound of maximum quantum violation of Noncontextual inequalities. 
%  Similar code for a different preparation equivalence is explained in https://quantum-journal.org/papers/q-2021-06-29-484/

%  Equivalence condition on preparation: 
%     1/4(X0+X1+X6+X7)=1/4(X2+X3+X4+X5)=1/4(X0+X2+X5+X7)=1/4(X0+X3+X4+X7)=4(X1+X3+X4+X6)=1/4(X1+X2+X5+X6)=1/4(X0+X3+X5+X6)=1/4(X1+X2+X4+X7)
% Equivalence condition on measurement: 1/3(M_{0/0}+M_{0/1}+M_{0/2}) = 1/3(M_{1/0}+M_{1/1}+M_{1/2}) M_{K/Y}

X = 8; % eight settings for the preparation device
Y = 3; % three settings for the measurement device
K = 2; % binary outcomes
O = 2*Y*(K-1)+1; % total number of operators in our list  % \mathcal{O} = {I, U_k^y, U{\dagger}_k^y} for all k,y
conP = []; % an empty list for constraints on the moment matrix level
conM = []; % an empty list for constraints on the substrate level
Prob = cell(X,Y,K-1); % a cell for probabilities
S = []; % an empty list for upper bounds
G = cell(X,1); % cell for X moment matrices

for x = 0:X-1
    G{x+1}=sdpvar(O,O,'hermitian','complex'); % declaration of our SDP variables
    conP = [conP;G{x+1} >= 0;]; % semi-definite constraints
end

conP = [conP; G{1} + G{2} + G{7} + G{8} == G{6} + G{3} + G{5} + G{4}; G{1} + G{2} + G{7} + G{8} == G{6} + G{3} + G{1} + G{8}; G{1} + G{2} + G{7} + G{8} == G{2} + G{7} + G{5} + G{4}; G{1} + G{2} + G{7} + G{8} == G{1} + G{8} + G{4} + G{5}; G{1} + G{2} + G{7} + G{8} == G{6} + G{3} + G{2} + G{7}; G{1} + G{2} + G{7} + G{8} == G{6} + G{1} + G{4} + G{7}; G{1} + G{2} + G{7} + G{8} == G{2} + G{3} + G{8} + G{5}]; % preparation equivalences

idx = @(y, k, u) 2*(K-1)*y + 2*k + u + 2; % function to return the position of the operators 

for x = 0:X-1
    for y = 0:Y-1
        for k = 0:K-2 
            conM = [conM; G{x+1}(1,idx(y,k,0)) == G{x+1}(idx(y,k,1),1)]; 
            conM = [conM; G{x+1}(idx(y,k,0),1) == G{x+1}(1,idx(y,k,1))];
        end
    end
    for j = 1:O
        conM = [conM; G{x+1}(j,j) == 1]; % unitarity constraints
    end
end

for x = 0:X-1
    for j = 1:O
        for k = 0:K-2 
            sum1 = 0; sum2 = 0; 
            for y = 0:Y-1
                sum1 = sum1 + G{x+1}(j,idx(y,k,0)) + G{x+1}(j,idx(y,k,1));
                sum2 = sum2 + G{x+1}(idx(y,k,0),j) + G{x+1}(idx(y,k,1),j);
            end
        conM = [conM; sum1==0; sum2==0]; % measurement equivalences     
        end
        
    end
end

for x = 0:X-1
    for y = 0:Y-1
        for k = 0:K-2
            Prob{x+1,y+1,k+1} = 0.5 + 0.25 * (G{x+1}(1, idx(y,k,0)) + G{x+1}(1, idx(y,k,1)));
        end
    end
end


S1 = real(4*Prob{1,1,1} - 4*Prob{2,1,1} - 2*Prob{5,1,1} - 5*Prob{1,2,1} + 4*Prob{3,2,1} + 3*Prob{5,2,1});

diagnostics = optimize([conP;conM], -S1, sdpsettings('solver', 'sdpt3'));
S1=value(S1)
