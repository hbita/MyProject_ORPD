%Definisanje ogranicenja kontrolnih (upravljackih promenljvih) za OPF
function [down,up,dim]=ogranicenja(testsistem)

%If lower bounds of dimensions are the same, then 'down' is a value.
%Otherwise, 'down' is a vector that shows the lower bound of each dimension.
%This is also true for upper bounds of dimensions.

%Definicija test sistema:
switch testsistem
      case 'ts_wscc9'
        [~,cvorovi,generatori,transformatori,kompenzatori,~,~,~]=ts_wscc9;
     case 'ts_ieee14'
        [~,cvorovi,generatori,transformatori,kompenzatori,~,~,~]=ts_ieee14;
     case 'ts_ieee30'
        [~,cvorovi,generatori,transformatori,kompenzatori,~,~,~]=ts_ieee30;
     case 'ts_ieee39'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,~,Vpqmax,~]=ts_ieee39;
     case 'ts_ieee57'
        [~,cvorovi,generatori,transformatori,kompenzatori,~,~,Sbase]=ts_ieee57;
     case 'ts_ieee118'
        [~,cvorovi,generatori,transformatori,kompenzatori,~,~,Sbase]=ts_ieee118;
     case 'ts_bus6'
        [mreza,cvorovi,generatori,transformatori,kompenzatori,~,Vpqmax,Sbase]=ts_bus6;
end

%[mreza,cvorovi,generatori,transformatori,kompenzatori,Vpqmin,Vpqmax]=ts_ieee30;  %test sistem IEEE30 

[m,~]=size(cvorovi);cvor=cvorovi(:,1);tip=cvorovi(:,2);
[~,~]=size(generatori);gcvor=generatori(:,1);Vgmin=generatori(:,2);
Vgmax=generatori(:,3);Pgmin=generatori(:,4);Pgmax=generatori(:,5);
tmin=transformatori(:,2);tmax=transformatori(:,3);
Qgcmin=kompenzatori(:,2);Qgcmax=kompenzatori(:,3);

% Sortiranje cvorova (da bi se indentifikovao BLR cvor)
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
% Pgming=Pgmin;Pgmaxg=Pgmax;
% 
% for k=1:ngen
%     if gcvor(k)==sl
%         Pgming(k)=[]; Pgmaxg(k)=[]; %eliminasanje BLR (sl) cvora iz ogranicenja koja se odnose na upravljacke promenljive
%     end
% end
down=[Vgmin;tmin;Qgcmin]'; %donje ogranicenje (lower bounds) 
up=[Vgmax;tmax;Qgcmax]'; %gornje ogranicenje (upper bounds)
dim=length(down);


