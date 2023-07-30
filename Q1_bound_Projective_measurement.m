% Written according to https://quantum-journal.org/papers/q-2021-06-29-484/

% Scenario : X ={0, 1, 2, 3}; y = {0, 1, 2}; k = {0, 1}  X: Preparation, y: Measurement, k: Outcome choice
% Equivalence condition on preparation: (1/3)(X0 + X1 + X2) = (1/2)(X0 + X3) 
% No nontrivial equivalence condition on measurement 
% The existence of a quantum model for given behaviour p(k/x,y) implies the existence of a moment matrix (Gamma) 

function SRAC = SRACQ1()         
% \mathcal{O} : sequence of monomials of the operators (M_k^y) 
% \mathcal{O} = {I, M_0^0, M_0^1, M_0^2}
% Define moment matrix corresponding to each preparation
Gamma0 = sdpvar(4,4,'hermitian','complex');           %  (4,4) complex, hermitian matrix 
Gamma1 = sdpvar(4,4,'hermitian','complex');
Gamma2 = sdpvar(4,4,'hermitian','complex');
Gamma3 = sdpvar(4,4,'hermitian','complex');

Cons = [Gamma0>=0;Gamma1>=0;Gamma2>=0;Gamma3>=0];     % positive semidefinite 
                                                      
Cons = [Cons; (1/3)*(Gamma0 + Gamma1 + Gamma2) == (1/2)*(Gamma0 + Gamma3)]; % equivalence relation corresponding to the equivalence on preparations 
 
Cons = [Cons; Gamma0(1,1)==1;Gamma1(1,1)==1;Gamma2(1,1)==1;Gamma3(1,1)==1]; % trivial constraints 

ops = sdpsettings('solver','sdpt3','verbose',0)       % sdpt3 used as solver to solve the sdp problem; 'verbose' argument is optional

% In general M_k^y are POVMs but for this scenario(no nontrivial equivalence condition on measurement)
% considering M_k^y as Projection operators is enough to get all quantum behaviours. 

Cons = [Cons; Gamma0(1,2)==Gamma0(2,2); Gamma0(1,3) == Gamma0(3,3); Gamma0(1,4)==Gamma0(4,4)];  % follows from the assumption of projective measurement
Cons = [Cons; Gamma1(1,2)==Gamma1(2,2); Gamma1(1,3) == Gamma1(3,3); Gamma1(1,4)==Gamma1(4,4)];
Cons = [Cons; Gamma2(1,2)==Gamma2(2,2); Gamma2(1,3) == Gamma2(3,3); Gamma2(1,4)==Gamma2(4,4)];
Cons = [Cons; Gamma3(1,2)==Gamma3(2,2); Gamma3(1,3) == Gamma3(3,3); Gamma3(1,4)==Gamma3(4,4)];

% - p(0/0,0) + 2p(0/1,0) + p(0/0,1) - 2p(0/2,1) <= 2    ## Noncontextual inequality obtained for given scenario 
Objective = real(-Gamma0(1,2) + 2*Gamma1(1,2) +Gamma0(1,3) -2*Gamma2(1,3))  % corresponding to the above Noncontextual inequality 
%Objective = real(-Gamma0(1,2) +Gamma1(1,2) +Gamma4(1,2))

optimize(Cons,-Objective, ops); % performing the optimization   

SRAC = value(Objective);
