%OPF by fmincon solver (sqp, interior-point,....) - klasicne metode su zakon!!!
global testsistem Case
testsistem='ts_ieee30';
Case='Ploss';

[low,up,dim]=ogranicenja(testsistem);
[X]=initialization(dim,1,up,low);
L0=X(1,:);
%%-------------------------------------------------------------------------
%RateLS=0.95;
%RateLS=1;
%if iteration>max_it*0.50
%if rand>RateLS % We decide to perform a Local Search.
    options=optimset('MaxFunEvals',1e4,'Display','off',...
    'algorithm','sqp','UseParallel','never');
    [LbestLS,FbestLS] = fmincon(@Fobj_sqp,L0,[],[],[],[],low,up,[],options);
%    iteration
%end
%end
%%------------------------------------------------------------------------- 
st(LbestLS,testsistem)