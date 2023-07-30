% Finding the lower bound of maximum quantum value in contextuality scenario 
% Scenario: X(Preparation)={0,1,2,3,4,5,6,7}; Y(Measurement)={0,1,2}; Z(Outcome)={0,1} 
% Equivalence condition on preparation:
% 1/4(X0+X1+X6+X7)=1/4(X2+X3+X4+X5)=1/4(X0+X2+X5+X7)=1/4(X0+X3+X4+X7)=4(X1+X3+X4+X6)=1/4(X1+X2+X5+X6)=1/4(X0+X3+X5+X6)=1/4(X1+X2+X4+X7)
% Equivalence condition on measurement: 1/3(M_{0/0}+M_{0/1}+M_{0/2}) = 1/3(M_{1/0}+M_{1/1}+M_{1/2}) M_{Z/Y} 
% Optimization over both preparation & measurement variable
% First obtain the Noncontextual inequalities according to the method provided in ...

d=2; % dimension of preparation states   ## increasing the dimension may give better bounds in some scenarios

    rho0_v = sdpvar(d,d,'hermitian','complex');  % density matrix corresponding to each preparation
    rho1_v = sdpvar(d,d,'hermitian','complex'); 
    rho2_v = sdpvar(d,d,'hermitian','complex'); 
    rho3_v = sdpvar(d,d,'hermitian','complex'); 
    rho4_v = sdpvar(d,d,'hermitian','complex'); 
    rho5_v = sdpvar(d,d,'hermitian','complex'); 
    rho6_v = sdpvar(d,d,'hermitian','complex'); 
    rho7_v = sdpvar(d,d,'hermitian','complex'); 

ConStates = [rho0_v>=0;trace(rho0_v)==1;rho1_v>=0;trace(rho1_v)==1;rho2_v>=0;trace(rho2_v)==1; ...  % Trivial constraints
    rho3_v>=0;trace(rho3_v)==1;rho4_v>=0;trace(rho4_v)==1;rho5_v>=0;trace(rho5_v)==1;...
    rho6_v>=0;trace(rho6_v)==1;rho7_v>=0;trace(rho7_v)==1];
ConStates = [ConStates; (1/4)*(rho0_v + rho1_v + rho6_v + rho7_v) == (1/4)*(rho2_v + rho3_v + rho4_v + rho5_v)]; % constraints from preparation equivalence

ConStates = [ConStates; (1/4)*(rho2_v + rho3_v + rho4_v + rho5_v) == (1/4)*(rho0_v + rho2_v + rho5_v + rho7_v)];

ConStates = [ConStates; (1/4)*(rho0_v + rho2_v + rho5_v + rho7_v) == (1/4)*(rho1_v + rho3_v + rho4_v + rho6_v)];

ConStates = [ConStates; (1/4)*(rho1_v + rho3_v + rho4_v + rho6_v) == (1/4)*(rho0_v + rho3_v + rho4_v + rho7_v)];

ConStates = [ConStates; (1/4)*(rho0_v + rho3_v + rho4_v + rho7_v) == (1/4)*(rho1_v + rho2_v + rho5_v + rho6_v)];

ConStates = [ConStates; (1/4)*(rho1_v + rho2_v + rho5_v + rho6_v) == (1/4)*(rho0_v + rho3_v + rho5_v + rho6_v)];

ConStates = [ConStates; (1/4)*(rho0_v + rho3_v + rho5_v + rho6_v) == (1/4)*(rho1_v + rho2_v + rho4_v + rho7_v)];

ConStates = [ConStates; (1/4)*(rho1_v + rho2_v + rho4_v + rho7_v) == (1/4)*(rho0_v + rho1_v + rho6_v + rho7_v)];

% Three measurements corresponding to three choices of Y. We converted all outcome 1 measurements to outcome 0 by M_{1/Y} = I - M_{0/Y}
 M0_v = sdpvar(d,d,'hermitian','complex'); 
 M1_v = sdpvar(d,d,'hermitian','complex'); 
 M2_v = sdpvar(d,d,'hermitian','complex'); 
 ConMea = [M0_v>=0;M1_v>=0;M2_v>=0;M0_v<=eye(d);M1_v<=eye(d);M2_v<=eye(d);2*M0_v+2*M1_v+2*M2_v==3*eye(d)]; % constraints from measurement equivalence
    
 
rho0 = RandomDensityMatrix(d);
rho1 = RandomDensityMatrix(d);
rho2 = RandomDensityMatrix(d);
rho3 = RandomDensityMatrix(d);
rho4 = RandomDensityMatrix(d);
rho5 = RandomDensityMatrix(d);
rho6 = RandomDensityMatrix(d);
rho7 = RandomDensityMatrix(d);


ops = sdpsettings('solver','sdpt3','verbose',0)  % Uses 'sdpt3' solver

for iter = 1:10  % increase iteration steps if it doesn't converge
    
    objective_mea = trace(M0_v*(-2*rho0+2*rho1+2*rho2) + M1_v*(-rho0+2*rho1));
   
    Q1_mea=optimize(ConMea,-real(objective_mea),ops);  % Optimize over measurement
    M0=value(M0_v);
    M1=value(M1_v);
    M2=value(M2_v);
    value(objective_mea);
    
    
    objective_state =  trace(M0*(-2*rho0_v+2*rho1_v+2*rho2_v) + M1*(-rho0_v+2*rho1_v)); 
    
    Q1_state=optimize(ConStates,-real(objective_state),ops); % Optimize over states
    rho0=value(rho0_v);
    rho1=value(rho1_v);
    rho2=value(rho2_v);
    rho3=value(rho3_v);
    rho4=value(rho4_v);
    rho5=value(rho5_v);
    rho6=value(rho6_v);
    rho7=value(rho7_v);
    value(objective_state)
    

    
end 
M0
M1