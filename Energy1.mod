
var
C           (long_name='consumption')
L           % Labor
W           % Wage
R           % Nominal Interest rate
pai         % Inflation
lamba       % Lagrange multiplier
San         % 
Y           % Output
A_Y         % Company Technology
K_Y         % Company Capital
O           % Energy Aggregate
P_O         % Energy Price
R_Y         % Company Capital Return
Q_Y         % Company Asset Price
O_F         % Fossil Energy
O_G         % Green Energy
P_F         % Fossil Energy Price
P_G         % Green Energy Price
A_F         % Fossil Energy Technology
K_F         % Fossil Energy Capital
X           % Fossil resource
R_F         % Fossil Capital Return
Q_F         % Fossil Asset Price
A_G         % Green Energy Technology
K_G         % Green Energy Capital
R_G         % Green Energy Capital Return
Q_G         % Green Asset Price
I_Y         % Company Investment
I_F         % Fossil Energy Investment
I_G         % Green Energy Investment
S_Y         % Company Asset
S_F         % Fossil Energy Asset
S_G         % Green Energy Asset
N           % Bank Net Worth
B           % Deposits
N_E         % Existing Bank Net Worth
N_N         % New Bank Net Worth
Zeta        % Risk-weighted Leverage
V_N         %
I
V
spread_Y
spread_F
spread_G
PSI_F
PSI_G
PHI_SS
tao_X
Xs
F_S
K_S
lev
G
PXX
pai_1
X_1
X_2
MC
d
;
varexo
e_AY
e_AF
e_AG
e_X
e_Xs
e_FS
e_KS
e_PX
;
parameters
beta        % Discount Factor
a           % Internal Habit Parameter
varphi      % Frisch Elasticity Of Labor Supply
alpha       % Capital Share In Output
miu         % Labor Share In Output
theta       %
epsilon     %
delta_Y     % Depreciation Rate Of Company Capital
omega_O     % Share Of Fossil Energy
epsilon_O   % Elasticity
alpha_F     % Capital Share In Fossil Energy Output
delta_F     % Depreciation Rate Of Fossil Energy Capital
delta_G     % Depreciation Rate Of Green Energy Capital
tao_Y       % Investment adjustment costs of company
tao_F       % Investment adjustment costs of fossil energy
tao_G       % Investment adjustment costs of green energy
M           % 
omega_N     % 
lamba_Y     % Absconding rate
lamba_F     % Relative absconding rate sector fossil energy
lamba_G     % Relative absconding rate sector green energy
rho_R       % Smoothing of the Taylor rule
rho_Y       % Output coefficient of the Taylor rule
rho_pai     % Inflation coefficient of the Taylor rule
rho_AY      % Autoc. shock
rho_AF      % Autoc. shock
rho_AG      % Autoc. shock
AYss        % Steady State Of Company Technology
PX          % Steady State Of Fossil Resource Price
omega_G
theta_V
AYFss
AYGss
fai_L
taoXss
rho_X
rho_Xs
rho_FS
FSss
KSss
rho_KS
rho_PX
PXXs

Css           (long_name='consumption')
Lss           % Labor
Wss           % Wage
Rss           % Nominal Interest rate
paiss         % Inflation
lambass       % Lagrange multiplier
Sanss         % 
Yss           % Output
A_Yss         % Company Technology
K_Yss         % Company Capital
Oss           % Energy Aggregate
P_Oss         % Energy Price
R_Yss         % Company Capital Return
Q_Yss         % Company Asset Price
O_Fss         % Fossil Energy
O_Gss         % Green Energy
P_Fss         % Fossil Energy Price
P_Gss        % Green Energy Price
A_Fss         % Fossil Energy Technology
K_Fss         % Fossil Energy Capital
Xss           % Fossil resource
R_Fss         % Fossil Capital Return
Q_Fss         % Fossil Asset Price
A_Gss         % Green Energy Technology
K_Gss        % Green Energy Capital
R_Gss         % Green Energy Capital Return
Q_Gss         % Green Asset Price
I_Yss         % Company Investment
I_Fss         % Fossil Energy Investment
I_Gss         % Green Energy Investment
S_Yss         % Company Asset
S_Fss         % Fossil Energy Asset
S_Gss         % Green Energy Asset
Nss           % Bank Net Worth
Bss           % Deposits
N_Ess         % Existing Bank Net Worth
N_Nss        % New Bank Net Worth
Zetass        % Risk-weighted Leverage
V_Nss         %
Iss
Vss
spread_Yss
spread_Fss
spread_Gss
PSI_Fss
PSI_Gss
PHI_SSss
tao_Xss
Xsss
F_Sss
K_Sss
levss
Gss
PXXss
pai_1ss
X_1ss
X_2ss
MCss
dss
;
load para.mat;  % load mat file created in console
for jj=1:length(M_.param_names)
set_param_value(M_.param_names{jj},eval(M_.param_names{jj})); 
end;
model;
%Household
1/(C-a*C(-1))-a*beta/(C(+1)-a*C)=lamba;                                              % FOC Consumption 1
fai_L*(L)^(varphi)=lamba*W;                                                                % FOC Labor2
-lamba+beta*lamba(+1)*R/pai(+1)=0;                                                   % Euler Equation3
San=lamba(+1)/lamba;                                                                 %4
%Company
Y*d=A_Y*(K_Y(-1)*K_S)^(alpha)*L^(miu)*O^(1-alpha-miu);  
d=(1-theta)*(pai_1)^(-epsilon)*(pai)^(epsilon)+(pai)^(epsilon)*(theta)*d(-1);%7
pai^(1-epsilon)=(1-theta)*(pai_1)^(1-epsilon)+theta;%32
pai_1=epsilon*pai*X_1/(epsilon-1)/X_2;%8
X_1=lamba*MC*Y+theta*beta*X_1(+1)*pai(+1)^(epsilon);%9
X_2=lamba*Y+theta*beta*X_2(+1)*pai(+1)^(epsilon-1);%10                                     % Product Function5
W*L=MC*miu*Y;                                                                        % FOC Labor11
P_O*O=MC*(1-alpha-miu)*Y;                                                            % FOC Energy12                             
R_Y=(alpha*MC*Y/K_Y(-1)/K_S+(1-delta_Y)*Q_Y)*K_S/Q_Y(-1);                            % FOC Capital13
O=((omega_O)^(1/epsilon_O)*(O_F)^((epsilon_O -1)/epsilon_O)+(1-omega_O)^(1/epsilon_O)*(O_G)^((epsilon_O -1)/epsilon_O))^((epsilon_O)/(epsilon_O - 1)); % Energy CES14
O_F=omega_O*(P_F/P_O)^(-epsilon_O)*O;                                                % Fossil Energy Demands15   
O_G=(1-omega_O)*(P_G/P_O)^(-epsilon_O)*O;                                            % Green Energy Demands16
%Fossil Energy Product
O_F=A_F*(F_S*K_F(-1))^(alpha_F)*(X)^(1-alpha_F);                                             % Fossil Energy Product Function17
(1-alpha_F)*(P_F)*O_F/X=(PX*PXX+tao_X*Xs);                                                               % FOC Fossil Resource18
R_F=(alpha_F*(P_F)*O_F/K_F(-1)/F_S+(1-delta_F)*Q_F)*F_S/Q_F(-1); 
%Green Energy Product
O_G=A_G*K_G(-1);                                                                         % Green Energy Product Function20
R_G=(P_G*O_G/K_G(-1)+(1-delta_G)*Q_G)/Q_G(-1);                                     % FOC Green Capital21
%Capital Producers
K_Y=I_Y+(1-delta_Y)*K_Y(-1)*K_S;                                                         % Capital 22
K_F=I_F+(1-delta_F)*K_F(-1)*F_S;                                                         % Capital23
K_G=I_G+(1-delta_G)*K_G(-1);                                                         % Capital24
Q_Y=1+(tao_Y/2)*((I_Y/I_Y(-1))-1)^(2)+tao_Y*((I_Y/I_Y(-1))-1)*((I_Y/I_Y(-1)))-beta*San*tao_Y*((I_Y(+1)/I_Y)-1)*(I_Y(+1)/I_Y)^(2); % 25
Q_F=1+(tao_F/2)*((I_F/I_F(-1))-1)^(2)+tao_F*((I_F/I_F(-1))-1)*((I_F/I_F(-1)))-beta*San*tao_F*((I_F(+1)/I_F)-1)*(I_F(+1)/I_F)^(2); % 26
Q_G=1+(tao_G/2)*(I_G/I_G(-1)-1)^(2)+tao_G*(I_G/I_G(-1)-1)*(I_G/I_G(-1))-beta*San*tao_G*((I_G(+1)/I_G)-1)*(I_G(+1)/I_G)^(2); % 27
I=I_Y+I_F+I_G;%28
%Bank
Q_Y*K_Y=Q_Y*S_Y;%29
Q_F*K_F=Q_F*S_F;%30
Q_G*K_G=Q_G*S_G;%31
(Q_Y*S_Y)+Q_F*S_F+Q_G*S_G=N+B;
N=N_E+N_N;%33
N_N=omega_N*(Q_Y(-1)*S_Y(-1)+Q_F(-1)*S_F(-1)+Q_G(-1)*S_G(-1));%34
N_E=theta_V*((R_Y-R)*Q_Y(-1)*S_Y(-1)/N(-1)+(R_F-R)*Q_F(-1)*S_F(-1)/N(-1)+(R_G-R)*Q_G(-1)*S_G(-1)/N(-1)+R)*N(-1);%35
Zeta=(Q_Y*S_Y+lamba_F*Q_F*S_F+lamba_G*Q_G*S_G)/N;%36
V_N=lamba_Y*beta*San*((1-theta_V)+theta_V*V_N(+1))*R(+1)/(lamba_Y-beta*San*((1-theta_V)+theta_V*V_N(+1))*(R_Y(+1)-R(+1)));%37
V=V_N*N;%38
lamba_F*San*((1-theta_V)+theta_V*V_N(+1))*(R_Y(+1)-R(+1))=San*((1-theta_V)+theta_V*V_N(+1))*(R_F(+1)-R(+1));%39
lamba_G*San*((1-theta_V)+theta_V*V_N(+1))*(R_F(+1)-R(+1))=lamba_F*San*((1-theta_V)+theta_V*V_N(+1))*(R_G(+1)-R(+1));%40
Zeta=V_N/lamba_Y;%41
%Center Bank
(R/Rss)=(R(-1)/Rss)^(rho_R)*((Y/Y(-1))^(rho_Y)*(pai)^(rho_pai))^(1-rho_R); %44
G = omega_G *Y ;
Y=C+I_Y+I_F+I_G+(tao_Y/2)*(I_Y/I_Y(-1)-1)^(2)+(tao_F/2)*(I_F/I_F(-1)-1)^(2)+(tao_G/2)*(I_G/I_G(-1)-1)^(2) + G;%46
A_Y=rho_AY*A_Y(-1)+(1-rho_AY)*AYss+e_AY;
A_F=rho_AF*A_F(-1)+(1-rho_AF)*AYFss+e_AF;                                                                
A_G=rho_AG*A_G(-1)+(1-rho_AG)*AYGss+e_AG;                                                                                                                                         %
spread_Y=R_Y-R;
spread_F=R_F-R;
spread_G=R_G-R;
PSI_F=spread_F/spread_Y;
PSI_G=PSI_F*spread_G/spread_F;
PHI_SS=(S_Y+PSI_F*S_F+PSI_G*S_G)/N; 
tao_X=rho_X*tao_X(-1)+(1-rho_X)*tao_Xss+e_X;         
Xs = Xs(-1)^(rho_Xs)*exp(e_Xs); 
F_S=rho_FS*F_S(-1)+(1-rho_FS)*FSss-e_FS; 
K_S=rho_KS*K_S(-1)+(1-rho_KS)*KSss-e_KS; 
lev=O_G / O;
PXX = rho_PX*PXX(-1)+(1-rho_PX)*PXXs+e_PX;                                                       
end;
initval;
C=Css;
L=Lss;
W=Wss;
R=Rss;
pai=paiss;
lamba=lambass;
San=Sanss;
Y=Yss;
A_Y=A_Yss;
K_Y=K_Yss;
O=Oss;
P_O=P_Oss;
R_Y=R_Yss;
Q_Y=Q_Yss;
O_F=O_Fss;
O_G=O_Gss;
P_F=P_Fss;
P_G=P_Gss;
A_F=A_Fss;
K_F=K_Fss;
X=Xss;
R_F=R_Fss;
Q_F=Q_Fss;
A_G=A_Gss;
K_G=K_Gss;
R_G=R_Gss;
Q_G=Q_Gss;
I_Y=I_Yss;
I_F=I_Fss;
I_G=I_Gss;
S_Y=S_Yss;
S_F=S_Fss;
S_G=S_Gss;
N=Nss;
B=Bss;
N_E=N_Ess;
N_N=N_Nss;
Zeta=Zetass;
V_N=V_Nss;
I=Iss;
V=Vss;
spread_Y=spread_Yss;
spread_F=spread_Fss;
spread_G=spread_Gss;
PSI_F=PSI_Fss;
PSI_G=PSI_Gss;
PHI_SS=PHI_SSss;
tao_X=tao_Xss;
Xs=Xsss;
F_S=F_Sss;
K_S=K_Sss;
lev=levss;
G=Gss;
PXX=PXXss;
pai_1=pai_1ss;
X_1=X_1ss;
X_2=X_2ss;
MC=MCss;
d=dss;
end;
resid(1);
steady;
check;
model_diagnostics;
shocks;
var e_X = 0.0109^2;
end;
stoch_simul(irf=40,order=1) ;