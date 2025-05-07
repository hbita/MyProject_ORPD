function [Fbest,Lbest,BestChart]=ABC(numAgents,maxIter,testSystem,Case)
    [lb,ub,D]=constraints(testSystem);
    FoodNumber=numAgents; 
    limit=numAgents*D;
    maxCycle=maxIter;
    Range = repmat((ub-lb),[FoodNumber 1]);
    Lower = repmat(lb, [FoodNumber 1]);
    Foods = rand(FoodNumber,D) .* Range + Lower;
    for i=1:FoodNumber
        ObjVal(i)=obj_fun(Foods(i,:),testSystem,Case);
        Fitness=calculateFitness(ObjVal);
    end
    %reset trial counters
    trial=zeros(1,FoodNumber);
    %/*The best food source is memorized*/
    BestInd=find(ObjVal==min(ObjVal));
    BestInd=BestInd(end);
    GlobalMin=ObjVal(BestInd);
    GlobalParams=Foods(BestInd,:);
    iter=1;
    while ((iter <= maxCycle))

     % EMPLOYED BEE PHASE 
        for i=1:(FoodNumber)
            
            %/*The parameter to be changed is determined randomly*/
            Param2Change=fix(rand*D)+1;
            
            %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
            neighbour=fix(rand*(FoodNumber))+1;
        
            %/*Randomly selected solution must be different from the solution i*/        
                while(neighbour==i)
                    neighbour=fix(rand*(FoodNumber))+1;
                end;
            
        sol=Foods(i,:);
        %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
        sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
            
        %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
            ind=find(sol<lb);
            sol(ind)=lb(ind);
            ind=find(sol>ub);
            sol(ind)=ub(ind);
            
            %evaluate new solution
            %ObjValSol=feval(objfun,sol);
            ObjValSol=obj_fun(sol,testSystem,Case);
            FitnessSol=calculateFitness(ObjValSol);
            
        % /*a greedy selection is applied between the current solution i and its mutant*/
        if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
                Foods(i,:)=sol;
                Fitness(i)=FitnessSol;
                ObjVal(i)=ObjValSol;
                trial(i)=0;
            else
                trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
        end
        end
        prob=(0.9.*Fitness./max(Fitness))+0.1;
        i=1;
        t=0;
        while(t<FoodNumber)
            if(rand<prob(i))
                t=t+1;
                %/*The parameter to be changed is determined randomly*/
                Param2Change=fix(rand*D)+1;
                
                %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
                neighbour=fix(rand*(FoodNumber))+1;
               
                %/*Randomly selected solution must be different from the solution i*/        
                    while(neighbour==i)
                        neighbour=fix(rand*(FoodNumber))+1;
                    end;
                
               sol=Foods(i,:);
               %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
               sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
                
               %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
                ind=find(sol<lb);
                sol(ind)=lb(ind);
                ind=find(sol>ub);
                sol(ind)=ub(ind);
                
                %evaluate new solution
                %ObjValSol=feval(objfun,sol);
                ObjValSol=obj_fun(sol,testSystem,Case);
                FitnessSol=calculateFitness(ObjValSol);
                
               % /*a greedy selection is applied between the current solution i and its mutant*/
               if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
                    Foods(i,:)=sol;
                    Fitness(i)=FitnessSol;
                    ObjVal(i)=ObjValSol;
                    trial(i)=0;
                else
                    trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
               end;
            end;
            
            i=i+1;
            if (i==(FoodNumber)+1) 
                i=1;
            end  
        end
        ind=find(ObjVal==min(ObjVal));
        ind=ind(end);
        if (ObjVal(ind)<GlobalMin)
        GlobalMin=ObjVal(ind);
        GlobalParams=Foods(ind,:);
        end
        ind=find(trial==max(trial));
        ind=ind(end);
        if (trial(ind)>limit)
            Bas(ind)=0;
            sol=(ub-lb).*rand(1,D)+lb;
            %ObjValSol=feval(objfun,sol);
            ObjValSol=Fobj(sol,testSystem,Case);
            FitnessSol=calculateFitness(ObjValSol);
            Foods(ind,:)=sol;
            Fitness(ind)=FitnessSol;
            ObjVal(ind)=ObjValSol;
        end
        GLM(iter)=GlobalMin; 
        iter=iter+1;
    end
    Fbest=GlobalMin;
    Lbest=GlobalParams;
    BestChart=GLM;
function fFitness=calculateFitness(fObjV)
    fFitness=zeros(size(fObjV));
    ind=find(fObjV>=0);
    fFitness(ind)=1./(fObjV(ind)+1);
    ind=find(fObjV<0);
    fFitness(ind)=1+abs(fObjV(ind));