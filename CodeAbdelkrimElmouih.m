%Mini projet “Estimation des paramètres du modèles thermodynamique UNIFAC appliqué en transfert de matière”
% Question 1 : Détermination des paramètres a12,a21 dans le cas de diffusion de Méthanol (A) dans l’eau(B) à T=313,13K avec Xa=0,35.
R=[1.4311   0.92]; 
q=[1.432   1.4];
Dab0=[2.1*10^(-5)   2.67*10^(-5)]; %Dab0ab et Dab0ba.
DabExpab=1.33*10^-5; %la valeur expérimentale.
Lam=[R(1).^(1/3)   R(2).^(1/3)]; %Lambda 1 et Lambda 2.
X=[0.35   1-0.35]; %Fractions molaires methanol / eau.
T=313.13; 
TET= [(X(1)*q(1))/(X(1)*q(1)+X(2)*q(2))   (X(2)*q(2))/(X(1)*q(1)+X(2)*q(2))]; % TETA(i)
Phi=[(X(1)*Lam(1))/(X(1)*Lam(1)+X(2)*Lam(2))   (X(2)*Lam(2))/(X(1)*Lam(1)+X(2)*Lam(2))]; %PHI(i)
% TO(ij) et TETA(ij) sont cités directement dans l'equation ci dessous :
%les boucle sont utilises pour avoir la plus divergence possible
for i=0 : 100 : 1000
    for j=0 : 100 : 1000
fun=@(A) abs(exp(X(2)*log(Dab0(1))+X(1)*Dab0(2)+2*(X(1)*log(X(1)./Phi(1))+X(2)*log(X(2)./Phi(2)))+...
    2*X(1)*X(2)*((Phi(1)./X(1))*(1-(Lam(1)./Lam(2)))+(Phi(2)./X(2))*(1-(Lam(2)./Lam(1))))+...
    X(2)*q(1)*((1-((TET(2)*exp(-A(2)/T))/(TET(1)+TET(2)*exp(-A(2)/T))).^2)*log(exp(-A(2)/T))+...
    (1-((TET(2))/(TET(1)*exp(-A(1)/T)+TET(2))).^2)*exp(-A(1)/T)*log(exp(-A(1)/T)))+...
    X(1)*q(2)*((1-((TET(1)*exp(-A(1)/T))/(TET(1)*exp(-A(1)/T)+TET(2))).^2)*log(exp(-A(1)/T))+...
    (1-((TET(1))/(TET(1)+TET(2)*exp(-A(2)/T))).^2)*exp(-A(2)/T)*log(exp(-A(2)/T))))-DabExpab);

A0 = [j,i];
A = fminsearch(fun,A0);
if i==0 && j==0
    min=fun(A);
end
if fun(A)<min
    min=fun(A);
    Aopt=A;
end
    end
end
fprintf('les paramètres A12 et A21 estimés sont :')
A = Aopt
fprintf('L ecart entre la valeur expérimentale et théorique (en pourcentage) :')
Ecart = 100*min/DabExpab 
%ce programme est destiné de retrouver le petit ecart exist entre la
%valeur expérimentale et théorique, autrement dit la fonction "fminsearch"
%utilisé ici vas etre minimisé la valeur absolue entre ces dernieres.


%Traçage de la variation du coefficient de diffusion en fonction du fraction molaire jusqu'à : Xa=0,7.

TO=[exp(-A(1)/T) exp(-A(2)/T) 1 1]; 
TETA=[(TET(1)*TO(1))/(TET(1)*TO(1)+TET(2)*TO(4)) (TET(2)*TO(2))/(TET(1)*TO(3)+TET(2)*TO(2))...
    (TET(1)*TO(3))/(TET(1)*TO(3)+TET(2)*TO(2))  (TET(2)*TO(4))/(TET(1)*TO(1)+TET(2)*TO(4))]; 

X = 0:0.0001:0.7;
Y = exp((1-X).*log(Dab0(1))+X.*Dab0(2)+2.*(X.*log(X./Phi(1))+(1-X).*log((1-X)./Phi(2)))+...
    2.*X.*(1-X).*((Phi(1)./X).*(1-(Lam(1)./Lam(2)))+(Phi(2)./(1-X)).*(1-(Lam(2)./Lam(1))))+...
    (1-X).*q(1).*((1-TETA(2).^2).*log(TO(2))+(1-TETA(4).^2)*TO(1).*log(TO(1)))+X.*q(2).*((1-TETA(1).^2).*log(TO(1))+(1-TETA(3).^2).*TO(2).*log(TO(2))));
plot(X,Y,'r')
title('Dab=f(Xa)')
xlabel('=== Xa ===>')
ylabel('=== Dab ===>')
