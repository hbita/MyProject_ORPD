function [] = st(L,testsistem,Case,N,max_it,BRP,FBEST,TOC)
% =============================================
%PRORACUN MULTI-OBJEKTIVNE FUNKCIJE KOD OPTIMALNIH TOKOVA SNAGA. 
%Stott-ov raspregnuti postupak se koristi za proracun tokova snaga.
%Jordan Radosavljevic, jun, 2012., FTN KM
% =============================================
%Fobj - objektivna (multi-objektivna) funkcija
%L - vektor upravljackih promenljivih;
%% default arguments
                                
if nargin < 3
    Case = 'Ploss'; N=25;max_it=100;BRP=1;FBEST=1;TOC=1;
end
 

%Definicija test sistema:
switch testsistem
%      case 'ts_wscc9'
%        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_wscc9;
%     case 'ts_ieee14'
%        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee14;
     case 'ts_ieee30'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee30;
%     case 'ts_ieee39'
%        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee39;
%    case 'ts_ieee57'
%        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee57;
     case 'ts_ieee118'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_ieee118;
%     case 'ts_bus6'
%        [mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax,Sbase]=ts_bus6;
end
%[mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax]=ts_ieee30;   %test sistem IEEE30 
%--------------------------------------------------------------------------
[m,n]=size(cvorovi);cvor=cvorovi(:,1);tip=cvorovi(:,2);
V0=cvorovi(:,3);teta0=cvorovi(:,4);Pg=cvorovi(:,5);Qg=cvorovi(:,6);
ppot=cvorovi(:,7);qpot=cvorovi(:,8);Vbase=cvorovi(:,10);
%--------------------------------------------------------------------------
t=mreza(:,7);
%--------------------------------------------------------------------------
[ngen,nkolg]=size(generatori);gcvor=generatori(:,1);Vgmin=generatori(:,2);
Vgmax=generatori(:,3);Pgmin=generatori(:,4);Pgmax=generatori(:,5);
Qgmin=generatori(:,6);Qgmax=generatori(:,7);
a=generatori(:,8);b=generatori(:,9);c=generatori(:,10);d=generatori(:,11);e=generatori(:,12);
%--------------------------------------------------------------------------
[ntf,nkolt]=size(transformatori);redbrt=transformatori(:,1);
tmin=transformatori(:,2);tmax=transformatori(:,3);tfbusi=transformatori(:,4);tfbusj=transformatori(:,5);
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

%--------------------------------------------------------------------------
%Snage generisanja 
%Pg=pgen;
Pg(sl)=pinj(sl)+ppot(sl);
Qg=qinj+qpot;
%--------------------------------------------------------------------------
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
br=1;
for ii=1:m
for jj=1:m
if abs(tokovip(ii,jj)+sqrt(-1)*tokoviq(ii,jj))~=0
tokovis(br,1:2)=[ii jj];
tokovis(br,3)=tokovip(ii,jj)+sqrt(-1)*tokoviq(ii,jj);
tokovis(br,4)=tokovip(ii,jj)+sqrt(-1)*tokoviqX(ii,jj)+tokovip(jj,ii)+sqrt(-1)*tokoviqX(jj,ii);
br=br+1;
end;end;end
%--------------------------------------------------------------------------
%Tokovi snaga po granama mreze
Cvor_i=tokovis(:,1); 
Cvor_j=tokovis(:,2);
Pij=real(tokovis(:,3));
Qij=imag(tokovis(:,3));
%--------------------------------------------------------------------------
%Gubici snage po granama i ukupni gubici u sistemu
pgub_ij=real(tokovis(:,4)); 
qgub_ij=imag(tokovis(:,4));
pgub_sum=sum(pinj);          %Ukupni gubici aktivne snage
qgub_sum=sum(qinj);
qgub_sumX=sum(qgub_ij)/2; %Ukupni gubici reaktivne snage (suma gubitaka na rednim reaktansama grana mreze)
Qshunt=(qgub_sum-qgub_sumX); %Ukupna reaktivna snaga koju generisu otocne kapacitivnosti (vodova i/ili onih prikljucenih u cvorove mreze)
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
%Odstupanje napona u PQ cvorovima od nominalne vrednosti 1 p.u.
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
lambdaPg=1000;lambdaVpq=100000;lambdaQg=10000; %Tezinski faktori prema: meni za 0.95<Vpq<1.05

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
%      case 'Fcost'              %Fuel cost minimization
%      Fobj=Fcost+Pf_Pgsl+Pf_Vpq+Pf_Qg;  
      case 'Ploss'              %Power loss minimization
      Fobj=pgub_sum+Pf_Pgsl+Pf_Vpq+Pf_Qg;
      case 'VD'               %Voltage profile improvement
      Fobj=FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%      case 'Fcost_Ploss'        %Fuel cost and power loss minimization
%      Fobj=Fcost+1950*pgub_sum+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%      case 'Fcost_VD'         %Fuel cost minimization and voltage profile improvement
%      Fobj=Fcost+200*FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%      case 'Fcost_Ploss_VD'   %Fuel cost, Ploss minimization and voltage profile improvement
%      Fobj=Fcost+1950*pgub_sum+200*FdeltaVpq+Pf_Pgsl+Pf_Vpq+Pf_Qg;
%      case 'Fcostvpe'           %Fuel cost minimiz. with valve point effect
%      Fobj=Fcostvpe+Pf_Pgsl+Pf_Vpq+Pf_Qg;
end
%st(L,testsistem,Case,BRP,N,max_it,minFBEST,maxFBEST,meanFBEST,stdFBEST,meanTOC,method)
%--------------------------------------------------------------------------
%Priprema tabela za prikaz rezultata
%tabela_cvorovi=[cvor V teta*180/pi Pg*Sbase Qg*Sbase ppot*Sbase qpot*Sbase Vbase]; 
%tabela_grane=[Cvor_i Cvor_j Pij*Sbase Qij*Sbase pgub_ij*Sbase qgub_ij*Sbase];
%--------------------------------------------------------------------------
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
for k=1:ngen
Vgen(k)=L(k); %moduli napona generatorskih cvorova + BLR cvor
end
tabela_Vg=[gcvor Vgen'];
%--------------------------------------------------------------------------
for k=(ngen+1):(ngen+ntf)
ttf(k-ngen)=L(k); %prenosni odnosi regulacionih transformatora
end
tabela_T=[redbrt tfbusi tfbusj ttf'];
%--------------------------------------------------------------------------
for k=(ngen+ntf+1):(ngen+ntf+nkom)
Qkond(k-ngen-ntf)=L(k)*Sbase;            %snage kompenzatora su pri (proracunu tokova snaga predstavljene su kao reaktivne snage generisanja u cvorovima)
end
tabela_Qc=[kcvor Qkond'];
%--------------------------------------------------------------------------

disp('  R E S U L T S: ')
disp(' ');
fprintf('\n  ===========================')
fprintf('\n  OPTIMUM CONTROL VARIABLES: ')
fprintf('\n  ===========================')
%disp(' ')
%fprintf('\n  Generator active powers')
%fprintf('\n  -----------------------')
%fprintf('\n   Bus    Pg    ')
%fprintf('\n   No    (MW)   ')
%fprintf('\n  --------------')
%for i=1:npv
%fprintf('\n%5d%10.4f',tabela_Pg(i,:));
%end
disp(' ');
fprintf('\n  Generator voltages ')
fprintf('\n  -------------------')
fprintf('\n   Bus     Vg   ')
fprintf('\n   No    (p.u.) ')
fprintf('\n  --------------')
for i=1:ngen
fprintf('\n%5d%10.4f',tabela_Vg(i,:));
end
disp(' ');
fprintf('\n  Transformer tap settings   ')
fprintf('\n  ---------------------------')
fprintf('\n  Brnch   From   To      T   ')
fprintf('\n   No     Bus    Bus   (p.u.)')
fprintf('\n  ---------------------------')
for i=1:ntf
fprintf('\n%5d%7d%7d%10.4f',tabela_T(i,:));
end
disp(' ');
fprintf('\n  Shunt VAR compensations ')
fprintf('\n  ----------------------- ')
fprintf('\n   Bus     Qc   ')
fprintf('\n   No    (MVAr) ')
fprintf('\n  --------------')
for i=1:nkom
fprintf('\n%5d%10.4f',tabela_Qc(i,:));
end
disp(' ');
fprintf('\n  =================================================')
fprintf('\n  BEST RESULTS:                                    ')
fprintf('\n  =================================================')
fprintf('\n  Objective Function, OF: %7.5f (units)',Fobj);
%fprintf('\n  Fuel Cost, Fcost: %7.5f ($/h)',Fcost);
fprintf('\n  Power Loss, Ploss: %7.5f (MW)',pgub_sum*Sbase);
fprintf('\n  Voltage Deviation, VD: %7.5f (p.u.)',FdeltaVpq);
fprintf('\n  ------------------------------------------------')
disp(' ');
%--------------------------------------------------------------------------
%Statisticki pokazatelji rezultata proracuna:
minFBEST=Fobj;        %minimalna (najbolja) vrednsot obj. funkcije
maxFBEST=max(FBEST);   %maksimalna (najgora) vrednsot obj. funkcije
meanFBEST=mean(FBEST); %srednja vrednsot obj. funkcije
stdFBEST=std(FBEST);   %standardna devijacija vrednosti obj. funkcije u toku zadatok BRP
meanTOC=mean(TOC);     %srednje vreme trajanja proracuna
%--------------------------------------------------------------------------
if BRP > 1
fprintf('\n  ================================================')
fprintf('\n  STATISTICS:                                     ')
fprintf('\n  ================================================')
fprintf('\n  Number of Runs: %4d',BRP);
fprintf('\n  Minimum value of the OF: %7.5f (units)',minFBEST);
fprintf('\n  Maximum value of the OF: %7.5f (units)',maxFBEST);
fprintf('\n  Mean value of the OF: %7.5f (units)',meanFBEST);
fprintf('\n  Standard deviation of the OF: %7.5f (units)',stdFBEST);
fprintf('\n  Mean value of computational time: %7.5f (s)',meanTOC);
fprintf('\n  ------------------------------------------------')
end
%--------------------------------------------------------------------------
%Rezultati tokova snaga sa optimalnim kontrolnim promenljivama
Qc(1:m)=0;Vmin(1:m)=0;Vmax(1:m)=0;Qgmi(1:m)=0;Qgma(1:m)=0;
for k=1:m
    Pg(pq)=0;
    Qg(pq)=0;
    for kk=1:nkom
        if k==kcvor(kk)
            Qc(k)=Qkond(kk);
        end
    end
    if tip(k)==2
        Vmin(k)=Vpqmin;
        Vmax(k)=Vpqmax;
    else
    for kk=1:ngen
        if k==gcvor(kk)
            Vmin(k)=Vgmin(kk);
            Vmax(k)=Vgmax(kk);
            Qgmi(k)=Qgmin(kk);
            Qgma(k)=Qgmax(kk);
        end
    end
    end
end
Qgmin=Qgmi'*Sbase;Qgmax=Qgma'*Sbase;
disp(' ')
fprintf('\n  ========================================================')
fprintf('\n  BALANCE OF ACTIVE POWER:                                ')
fprintf('\n  ========================================================')
fprintf('\n  Load: %7.3f (MW)',sum(ppot)*Sbase);
fprintf('\n  Generation: %7.3f (MW)',sum(Pg)*Sbase);
fprintf('\n  Loss (R*I^2): %7.3f (MW)',pgub_sum*Sbase);
fprintf('\n  IMBALANCE: %7.3f (MW)',sum(Pg)*Sbase-sum(ppot)*Sbase-pgub_sum*Sbase);
fprintf('\n  --------------------------------------------------------')
disp(' ')
fprintf('\n  ========================================================')
fprintf('\n  BALANCE OF REACTIVE POWER:                              ')
fprintf('\n  ========================================================')
fprintf('\n  Load: %7.3f (MVAr)',sum(qpot)*Sbase);
fprintf('\n  Generation: %7.3f (MVAr)',sum(Qg)*Sbase);
fprintf('\n  Shunt VAR conmpensations: %7.3f (MVAr)',sum(Qkond));
fprintf('\n  Shunt admitances and Branch charging (inj): %7.3f (MVAr)',-Qshunt*Sbase);
fprintf('\n  Loss (X*I^2): %7.3f (MVAr)',qgub_sumX*Sbase);
fprintf('\n  IMBALANCE: %7.3f (MVAr)',sum(Qg)*Sbase+sum(Qkond)-Qshunt*Sbase-sum(qpot)*Sbase-qgub_sumX*Sbase);
fprintf('\n  --------------------------------------------------------')
disp(' ')
fprintf('\n  ========================================================')
fprintf('\n  VIOLATING CONSTRAINTS?                                  ')
fprintf('\n  ========================================================')
Vviolating(1:m)=0;Qgviolating(1:m)=0;
for k=1:m
    if V(k)<Vmin(k) || V(k)>Vmax(k)
        Vviolating(k)=k;
    end
    if Qg(k)<Qgmin(k) || Qg(k)>Qgmax(k)
        Qgviolating(k)=k;
    end
end
    Vviolat=find(Vviolating);
    Qgviolat=find(Qgviolating);
    if isempty(Vviolat) 
fprintf('\n  Voltages at all buses are within permissible limits                  ')
    else
fprintf('\n  Violating voltage limit at bus: %4d',Vviolat);        
    end
    if isempty(Qgviolat) 
fprintf('\n  Reactive power outputs of all generators are within permissible limits')
    else
fprintf('\n  Violating reactive power limit of generator at bus: %4d',Qgviolat);        
    end  
fprintf('\n  --------------------------------------------------------\n')

%----------------------------------------------------------------------------------------
%Priprema tabela za prikaz rezultata tokova snaga po cvorovima i granama
%tabela_cvorovi=[cvor V teta*180/pi Pg*Sbase Qg*Sbase Qc' ppot*Sbase qpot*Sbase Vbase Vmin' Vmax' Qgmin Qgmax];
tabela_cvorovi=[cvor V teta*180/pi Pg*Sbase Qg*Sbase Qc' ppot*Sbase qpot*Sbase]; 
tabela_grane=[Cvor_i Cvor_j Pij*Sbase Qij*Sbase pgub_ij*Sbase qgub_ij*Sbase];
%--------------------------------------------------------------------------
fprintf('\n  =========================================================================')
fprintf('\n  BUS VOLTAGES AND POWERS UNDER OPTIMAL CONTROL VARIABLES                  ')
fprintf('\n  =========================================================================')
fprintf('\n   Bus     V        teta       Pg       Qg        Qc       Pload     Qload ')
fprintf('\n   No    (p.u.)     (deg)     (MW)     (MVAr)    (MVAr)    (MW)     (MVAr) ')
fprintf('\n  -----  ----------------   ---------------------------    ----------------')
for i=1:m
fprintf('\n%5d%10.4f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f',tabela_cvorovi(i,:));
end
fprintf('\n  --------------------------------------------------------------------------')
disp(' ')
disp(' ')
fprintf('\n  ==========================================================')
fprintf('\n  BRANCH POWER FLOW AND LOSS UNDER OPTIMAL CONTROL VARIABLES')
fprintf('\n  ==========================================================')
fprintf('\n     From   To       P          Q        Ploss     Qloss    ')
fprintf('\n     Bus    Bus     (MW)      (MVAr)     (MW)      (MVAr)   ')
fprintf('\n     ----------    -----------------     ----------------   ')
for i=1:length(Cvor_i)
fprintf('\n%7d%7d%12.3f%10.3f%10.3f%10.3f',tabela_grane(i,:));
end
fprintf('\n  ----------------------------------------------------------\n')

L=L';       %vektor upravljackih promenljivih (optimalno resenje)
V;
Qg;
Pg;
return