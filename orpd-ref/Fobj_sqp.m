function [Fobj]=Fobj_sqp(L)
% =============================================
%PRORACUN MULTI-OBJEKTIVNE FUNKCIJE KOD OPTIMALNIH TOKOVA SNAGA. 
%Stott-ov raspregnuti postupak se koristi za proracun tokova snaga.
%Jordan Radosavljevic, jun, 2012., FTN KM
% =============================================
%Fobj - objektivna (multi-objektivna) funkcija
%L - vektor upravljackih promenljivih;
global testsistem Case

%Definicija test sistema:
switch testsistem
      case 'ts_wscc9'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_wscc9;
     case 'ts_ieee14'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee14;
     case 'ts_ieee30'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee30;
     case 'ts_ieee39'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee39;
     case 'ts_ieee57'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee57;
     case 'ts_ieee118'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee118;
     case 'ts_bus6'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_bus6;
end
%[mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax]=ts_ieee30;   %test sistem IEEE30 
%--------------------------------------------------------------------------
[m,n]=size(cvorovi);cvor=cvorovi(:,1);tip=cvorovi(:,2);
V0=cvorovi(:,3);teta0=cvorovi(:,4);Pg=cvorovi(:,5);Qg=cvorovi(:,6);
ppot=cvorovi(:,7);qpot=cvorovi(:,8);
%--------------------------------------------------------------------------
t=mreza(:,7);
%--------------------------------------------------------------------------
[ngen,nkolg]=size(generatori);gcvor=generatori(:,1);Vgmin=generatori(:,2);
Vgmax=generatori(:,3);Pgmin=generatori(:,4);Pgmax=generatori(:,5);
Qgmin=generatori(:,6);Qgmax=generatori(:,7);
a=generatori(:,8);b=generatori(:,9);c=generatori(:,10);d=generatori(:,11);e=generatori(:,12);
%--------------------------------------------------------------------------
[ntf,nkolt]=size(transformatori);redbrt=transformatori(:,1);
tmin=transformatori(:,2);tmax=transformatori(:,3);
%--------------------------------------------------------------------------
[nkom,nkolk]=size(kompenzatori);kcvor=kompenzatori(:,1);
Qgcmin=kompenzatori(:,2);Qgcmax=kompenzatori(:,3);
%--------------------------------------------------------------------------
% Sortiranje cvorova
npq=0; npv=0;
for k=1:m 
    if tip(k)==0
        sl=cvor(k);
    elseif tip(k)==2
        npq=npq+1;
        pq(npq)=cvor(k);
    elseif tip(k)==1
        npv=npv+1;
        pv(npv)=cvor(k);
    end
end
%--------------------------------------------------------------------------
%Definisanje upravljackih promenljivih (kontrolne promenljive cije vrednosti se optimiziraju) i njihovih granicnih vrednosti:
for k=1:ngen
V0(gcvor(k))=L(k);  %moduli napona generatorskih cvorova + BLR cvor
end

for k=(ngen+1):(ngen+ntf)
t(redbrt(k-ngen))=L(k); %prenosni odnosi regulacionih transformatora
end

for k=(ngen+ntf+1):(ngen+ntf+nkom)
Qg(kcvor(k-ngen-ntf))=L(k);            %snage kompenzatora su pri (proracunu tokova snaga predstavljene su kao reaktivne snage generisanja u cvorovima)
end
%--------------------------------------------------------------------------
%Snage injektiranja u cvorovima
pinj=Pg-ppot; 
qinj=Qg-qpot;
%--------------------------------------------------------------------------
%Izracunavanje matrice admitansi nezavisnih cvorova[Y] 
[Y,yg,yg0]=Ybus(mreza,cvorovi,t);
G=real(Y); B=imag(Y);
ggr=real(yg); bgr=imag(yg);
ggr0=real(yg0); bgr0=imag(yg0);
%Formiranje redukovanih matrica susceptansi
B1=B-bgr0;
for k=1:m
    if tip(k)==0
        B1(:,k)=[];
        B1(k,:)=[];
    end
end
B2rr=B-bgr0;
for k=1:m
    if tip(k)==0 || tip(k)==1
        B2rr(:,k)=0;
        B2rr(k,:)=0;
    end
    if tip(k)==2
        index(k)=k;
    end
end
index=nonzeros(index);
for k=1:length(index)
    B2r(k,:)=B2rr(index(k),:);
end
for k=1:length(index)
    B2(:,k)=B2r(:,index(k));
end
%--------------------------------------------------------------------------

%Iterativna petlja proracuna tokova snaga Stotovom metodom
%--------------------------------------------------------------------
%Potrebne deklaracije:
iter=0; itermax=100; epsilon=0.00001; maxraz=epsilon+1;
V=V0;teta=teta0;
while (iter < itermax) && (maxraz > epsilon) 
iter=iter+1;
tetanew=teta;Vnew=V;
    for k=1:m
        suma=0;
        for j=1:m
            suma=suma+Vnew(j)*(G(k,j)*cos(tetanew(k)-tetanew(j))...
                           +B(k,j)*sin(tetanew(k)-tetanew(j)));
        end
            deltaP(k)=pinj(k)-Vnew(k)*suma;
    end
    deltaP(sl)=[];
    Vnew(sl)=[];
    tetanew(sl)=[];
    teta(sl)=[];
    dteta=-inv(B1)*(deltaP'./Vnew);
    teta=tetanew+dteta;
    tetapomoc=zeros(m,1);
    tetapomoc(sl)=teta0(sl);
    tetapomoc(1:sl-1)=teta(1:sl-1);
    tetapomoc(sl+1:m)=teta(sl:m-1);
    teta=tetapomoc;
    Vnewpomoc=zeros(m,1);
    Vnewpomoc(sl)=V0(sl);
    Vnewpomoc(1:sl-1)=Vnew(1:sl-1);
    Vnewpomoc(sl+1:m)=Vnew(sl:m-1);
    Vnew=Vnewpomoc;
    for k=1:m
        suma=0;
        for j=1:m
            suma=suma+Vnew(j)*(G(k,j)*sin(teta(k)-teta(j))...
                           -B(k,j)*cos(teta(k)-teta(j)));
        end
        deltaQ(k)=qinj(k)-Vnew(k)*suma;    
    end
B2i=-inv(B2);
br1=0;
for k=1:m
    if tip(k)==2
        br1=br1+1;
        Vpq(br1)=Vnew(k);
        deltaQpq(br1)=deltaQ(k);
    end
end
for kk=1:npq
deltaV(kk)=B2i(kk,:)*(deltaQpq'./Vpq');
Vpqnew(kk)=Vpq(kk)+B2i(kk,:)*(deltaQpq'./Vpq');
end
br=0;
for k=1:m
    if tip(k)==2
    br=br+1;
    V(k)=Vpqnew(br);
    end
end
 
maxraz=max(deltaP);
end  %Kraj iterativne petlje
%--------------------------------------------------------------------------

% Proracun snaga injektiranja
for k=1:m
    sumap=0;
    sumaq=0;
    for l=1:m
        sumap=sumap+V(l)*(G(k,l)*cos(teta(k)-teta(l))...
                         +B(k,l)*sin(teta(k)-teta(l)));
        sumaq=sumaq+V(l)*(G(k,l)*sin(teta(k)-teta(l))...
                         -B(k,l)*cos(teta(k)-teta(l)));
    end
    pinj(k)=V(k)*sumap;
    qinj(k)=V(k)*sumaq;
    sinj(k)=pinj(k)+qinj(k)*sqrt(-1);
end

%Snage generisanja 
Pg=pinj+ppot; Qg=qinj+qpot;

%Proracun tokova aktivnih i reaktivnih snaga po granama EES
for i=1:m
for j=1:m
tokovip(i,j)=(V(i)^2)*ggr(i,j)-V(i)*V(j)*(ggr(i,j)*cos(teta(i)-teta(j))...
                                           +bgr(i,j)*sin(teta(i)-teta(j)));
tokoviq(i,j)=(-V(i)^2)*bgr(i,j)-(V(i)^2)*bgr0(i,j)-V(i)*V(j)...              %Reaktivna snaga po rednim granama sa otocnim admitansama
            *(ggr(i,j)*sin(teta(i)-teta(j))-bgr(i,j)*cos(teta(i)-teta(j)));
tokoviqX(i,j)=(-V(i)^2)*bgr(i,j)-V(i)*V(j)...                                %Reaktivna snaga po rednim granama bez otocnih admitansi !!!
            *(ggr(i,j)*sin(teta(i)-teta(j))-bgr(i,j)*cos(teta(i)-teta(j)));
end
end

%Gubici snaga u sistemu
pgub_sum=sum(pinj);
qgub_sum=sum(sum(tokoviqX));
%--------------------------------------------------------------------------
%Funkcija troskova generatora: a+b*P+c*P^2
Fcost=0;
for k=1:ngen
    Cg(k)=a(k)+b(k)*Pg(gcvor(k))+c(k)*Pg(gcvor(k))^2;
    Fcost=Fcost+Cg(k);
end
%--------------------------------------------------------------------------
%Funkcija troskova with valve point effect
Fcostvpe=0;
for k=1:ngen
    Cgvpe(k)=a(k)+b(k)*Pg(gcvor(k))+c(k)*Pg(gcvor(k))^2+abs(d(k)*sin(e(k)*(Pgmin(k)-Pg(gcvor(k)))));
    Fcostvpe=Fcostvpe+Cgvpe(k);
end
%--------------------------------------------------------------------------
%Odstupanje napona od nominalne vrednosti 1 p.u.
FdeltaVpq=0; %suma odstupanja napona u svim PQ cvorovima
for k=1:npq
    FdeltaVpq=FdeltaVpq+abs(V(pq(k))-1);
end
%--------------------------------------------------------------------------
%Penalne funkcije-pri narusavanju funcionalnih ogranicenja. 
%Ove penalne funkcije se ugradjuju u objektivnu funkciju
%==========================================================================
%lambda=mean(c); % za tezinski faktor je uzeta srednja vrednost koeficijenta c funkcije troskova goriva.
%Penalna funkcija za aktivnu snagu BLR sabirnicu (slack bus - oznaka sl) 
%lambdaPg=5*10^6;lambdaVpq=50000;lambdaQg=50000; %Tezinski faktori prema: M. Ghasemi, A novel hybrid algorithm.... (2014)
%lambdaPg=10;lambdaVpq=20;lambdaQg=0.02; %Tezinski faktori prema: O. Alsac, B. Stott, Optimal load flow with steady-state security (1974)
%lambdaPg=1000;lambdaVpq=5000;lambdaQg=1000; %Tezinski faktori prema: meni za 0.95<Vpq<1.10
lambdaPg=10^6;lambdaVpq=10^6;lambdaQg=10^4; %Tezinski faktori prema: meni za 0.95<Vpq<1.05

for k=1:ngen
    if gcvor(k)==sl
if Pg(sl)<Pgmin(k)
    Pgsllim=Pgmin(k);
elseif Pg(sl)>Pgmax(k)
    Pgsllim=Pgmax(k);
else
    Pgsllim=Pg(sl);
end
    end
end
Pf_Pgsl=lambdaPg*(Pg(sl)-Pgsllim)^2;

%Penalna funkcija za napona u PQ cvorovima: Vpqmin=0.95; Vpqmax=1.05;
for k=1:npq
    if V(pq(k))<Vpqmin
        Vpqlim(k)=Vpqmin;
    elseif V(pq(k))>Vpqmax
        Vpqlim(k)=Vpqmax;
    else
       Vpqlim(k)=V(pq(k)); 
    end
    deltaVpq(k)=(V(pq(k))-Vpqlim(k))^2;
end
Pf_Vpq=lambdaVpq*sum(deltaVpq);

%Penalne funkcije za reaktivne snage generatora:
for k=1:ngen
    if Qg(gcvor(k))<Qgmin(k)
        Qglim(k)=Qgmin(k);
    elseif Qg(gcvor(k))>Qgmax(k)
        Qglim(k)=Qgmax(k);
    else
        Qglim(k)=Qg(gcvor(k));
    end
    deltaQg(k)=(Qg(gcvor(k))-Qglim(k))^2;
end
Pf_Qg=lambdaQg*(sum(deltaQg));

%Objektivna funkcija
%==========================================================================
switch Case
      case 'Fcost'              %Fuel cost minimization
      Fobj=Fcost+Pf_Pgsl+Pf_Vpq+Pf_Qg;  
      case 'Ploss'              %Power loss minimization
      Fobj=pgub_sum+Pf_Pgsl+Pf_Vpq+Pf_Qg;
      case 'VD'               %Voltage profile improvement
      Fobj=FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
      case 'Fcost_Ploss'        %Fuel cost and power loss minimization
      Fobj=Fcost+1950*pgub_sum+Pf_Pgsl+Pf_Vpq+Pf_Qg;
      case 'Fcost_VD'         %Fuel cost minimization and voltage profile improvement
      Fobj=Fcost+200*FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
      case 'Fcost_Ploss_VD'   %Fuel cost, Ploss minimization and voltage profile improvement
      Fobj=Fcost+1950*pgub_sum+200*FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
      case 'Fcostvpe'           %Fuel cost minimiz. with valve point effect
      Fobj=Fcostvpe+Pf_Pgsl+Pf_Vpq+Pf_Qg;
end
%Case1  Fuel cost minimization
%Fobj=Fcost+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%--------------------------------------------------------------------------
%Case x1  Power loss minimization
%Fobj=pgub_sum+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%--------------------------------------------------------------------------
%Case x2  Voltage profile improvement
%Fobj=FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%--------------------------------------------------------------------------
%Case 2  Fuel cost minimization and voltage profile improvement
%Fobj=Fcost+200*FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%--------------------------------------------------------------------------
%Case 3  Fuel cost and power loss minimization
%Fobj=Fcost+1950*pgub_sum+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%--------------------------------------------------------------------------
%Case y1  Fuel cost, Ploss minimization and voltage profile improvement
%Fobj=Fcost+1950*pgub_sum+200*FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%--------------------------------------------------------------------------
%Case 4a  Fuel cost minimiz. with valve point effect
%Fobj=Fcostvpe+Pf_Pgsl+Pf_Vpq+Pf_Qg;
return
