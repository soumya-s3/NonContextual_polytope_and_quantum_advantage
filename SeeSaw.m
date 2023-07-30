% Finding the lower bound of maximum quantum value in contextuality scenario 

% Scenario: X(Preparation)={0,1,2,3}; Y(Measurement)={0,1,2}; Z(Outcome)={0,1} 
% Equivalence condition on preparation: (1/3)(X0 + X1 + X2) = (1/2)(X0 + X3) 
% No nontrivial equivalence condition on measurement;    Optimization will only be over preparation variables
% First obtain the Noncontextual inequalities according to the method provided in ...

d=2; % dimension of preparation states   ## increasing the dimension may give better bounds in some scenarios

rho0_v = sdpvar(d,d,'hermitian','complex');  % density matrix corresponding to each preparation
rho1_v = sdpvar(d,d,'hermitian','complex'); 
rho2_v = sdpvar(d,d,'hermitian','complex'); 
rho3_v = sdpvar(d,d,'hermitian','complex'); 

ConStates = [rho0_v>=0;trace(rho0_v)==1;rho1_v>=0;trace(rho1_v)==1;rho2_v>=0;trace(rho2_v)==1;rho3_v>=0;trace(rho3_v)==1;(1/3)*(rho0_v+rho1_v+rho2_v)==(1/2)*(rho0_v+rho3_v)]; % Trivial Constraints
        
rho0 = RandomDensityMatrix(d);
rho1 = RandomDensityMatrix(d);
rho2 = RandomDensityMatrix(d);
rho3 = RandomDensityMatrix(d);
    
ops = sdpsettings('solver','sdpt3','verbose',0)   % Uses 'sdpt3' solver

for iter = 1:10  % increase iteration steps if it doesn't converge
    
%[V0,D0] = eig(-3*rho0-3*rho1+2*rho2+2*rho3)
[V0,D0] = eig(-rho0+2*rho1);[V1,D1]=eig(rho0-2*rho2); % choose according to objective state
%M0 = kron(V0(:,d),V0(:,d)');M1 = kron(V1(:,d),V1(:,d)');

M0 = zeros(d,d); M1 = zeros(d,d); 
for i = 1:d
   if D0(i,i) > 0
   M0 = M0 + kron(V0(:,i),V0(:,i)');   % taking the measurement operator as sum of projectors having +ve eigenvalue
   end                                    % eventually produces projection operator with largest eigenvalue
end 
for i = 1:d
   if D1(i,i) > 0
   M1 = M1 + kron(V1(:,i),V1(:,i)'); 
   end
end

    %objective_state = trace(M0*(-3*rho0_v-3*rho1_v+2*rho2_v+2*rho3_v));
    objective_state = trace(M0*(-rho0_v+2*rho1_v)+M1*(rho0_v-2*rho2_v));  % after each iteration update the measurement operators and optimize over states
    Q1_state=optimize(ConStates,-real(objective_state),ops);
    rho0=value(rho0_v);
    rho1=value(rho1_v);
    rho2=value(rho2_v);
    rho3=value(rho3_v);
    value(objective_state)    
end 

rho0
rho1
rho2
rho3
M0
M1