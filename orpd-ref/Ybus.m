function [Y,yg,yg0]=Ybus(mreza,cvorovi,t) 
%===================================
%Funkcijski program Ybus
%Ulazni argumenti:
%mreza - Matrica u kojoj se zadaju konfiguracija i parametri elemenata sistema
%cvorovi - Matrica u kojoj se zadaju snage, naponi i admitanse u cvorovima sistema za proracun tokova snaga
%t - vektor ciji su elementi prenosni odnosi transformatora
%Izlazni argumenti:
%Y - matrica admitansi nezavisnih cvorova mreze,
%yg - matrice admitansi rednih grana mreze
%yg0 - matrica admitansi otocnih grana mreze 
%Jordan Radosavljevic , FTN KM
%====================================
[ngr,npod]=size(mreza);cvi=mreza(:,1);cvj=mreza(:,2);tipel=mreza(:,3);
rgr=mreza(:,4);xgr=mreza(:,5);bgrz=mreza(:,6);
ygrz=1i*bgrz;
zgr=rgr+1i*xgr;
[m,n]=size(cvorovi);Yot=cvorovi(1:m,9);
Yot=diag(1i*Yot); % otocne admitanse u cvorovima (kondezatori, prigusnice...) 
for k=1:ngr %izracunavanje elemenata matrica yg i yg0
   zg(cvi(k),cvj(k))=zgr(k);
   zg(cvj(k),cvi(k))=zgr(k);
   yg(cvi(k),cvj(k))=1/zg(cvi(k),cvj(k))/t(k);
   yg(cvj(k),cvi(k))=1/zg(cvj(k),cvi(k))/t(k); 
   if tipel(k)==1 
   yg0(cvi(k),cvj(k))=ygrz(k)/2;
   yg0(cvj(k),cvi(k))=ygrz(k)/2;
   elseif tipel(k)==2
   yg0(cvi(k),cvj(k))=(1-t(k))/t(k)^2*(1/zgr(k));
   yg0(cvj(k),cvi(k))=(t(k)-1)/t(k)*(1/zgr(k));
   end
end
yg0=yg0+Yot; 
nbus=max(max(cvi,cvj));
for ii=1:nbus %izracunavanje elemenata matrice Y
    for jj=1:nbus
        if ii==jj
       Y(ii,jj)=sum(yg(ii,:))+sum(yg0(ii,:));
        else
       Y(ii,jj)=-yg(ii,jj);
        end
end
end