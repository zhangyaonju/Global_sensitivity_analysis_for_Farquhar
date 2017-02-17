%%
%stomatal conductance model
%Farquhar model
%temperature adjustment
function Ac = Farquhar_model(para,I,Tk,D,Ca,LAI)
%fix parameter
Ox = 0.21;
R = 8.314;
%get parameter
aq = para(1);
Kc25 = para(2);
Ekc = para(3);
Ko25 = para(4);
Eko =para(5);
Vm25 = para(6);
Evm = para(7);
Gama25 = para(8);
EGama = para(9);
kn = para(10);
rjv = para(11);
m = para(12);
D0 = para(13);

Vm = Vm25*exp(Evm.*(Tk-298)/R./Tk/298);
Gama = Gama25*exp(EGama.*(Tk-298)/R./Tk/298);
KC = Kc25*exp(Ekc.*(Tk-298)/R./Tk/298);
KO = Ko25*exp(Eko.*(Tk-298)/R./Tk/298);
Jm = rjv*Vm;

%calculation of Ci using four models

%Ci = Ca*(1-1/m/hs);                  %%%% Ball-Berry model
Ci = Ca-(Ca-Gama).*(1+D/10/D0)/m;     %%%% Leuning model
%Ci = Ca-sqrt((Ca-Tau)*1.6*D)/lamda); %%%% Farquhar model
%Ci = Ca(1-1/(1+g1/sqrt(D)));         %%%% Medlyn model

%assimilation
Jc = Vm.*(Ci-Gama)./(Ci+KC.*(1+Ox./KO));
J = aq.*I.*Jm./sqrt(Jm.^2+(aq.*I).^2);
Je = J.*(Ci-Gama)/4./(Ci+2*Gama);

A = min(Jc,Je);
A(A<0)=0;

Ac = A.*(1-exp(-kn.*LAI))/kn;






