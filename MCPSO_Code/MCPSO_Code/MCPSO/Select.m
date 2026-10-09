function [Population,FrontNo] = Select(Population,FrontNo,y)
% And one offspring to the population and delete the worst one

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Identify the solutions in the last front
Population = [Population,y];
PopObj     = Population.objs;
[N,M]      = size(PopObj);
FrontNo    = UpdateFront(PopObj,FrontNo);

refPoint=max(PopObj,[],1);
refPoint=refPoint+1;
%% Delete the worst solution
if max(FrontNo) > 1
    Distance  = pdist2(PopObj,PopObj);
    Distance(logical(eye(N))) = inf;
    LastFront = find(FrontNo==max(FrontNo));
    [~,worst] = min(min(Distance(LastFront,:),[],2));
    worst     = LastFront(worst);
else
    if M==2
        frontObjvs=PopObj;
        fitV=zeros(N,1);
        [frontObjvs,IX]=sortrows(frontObjvs,1);
        fitV(IX(1))=(frontObjvs(2,1)-frontObjvs(1,1)).* (refPoint(2)-frontObjvs(1,2));
        fitV(IX(2:N-1))=(frontObjvs(3:N,1)-frontObjvs(2:N-1,1)).* (frontObjvs(1:N-2,2)-frontObjvs(2:N-1,2));
        fitV(IX(N))=(refPoint(1)-frontObjvs(N,1)).*(frontObjvs(N-1,2)-frontObjvs(N,2));
        [~,worst]=min(fitV);
    else
        totalHV=hv(PopObj',refPoint);
        fitV=zeros(N,1);
        for i=1:N
            tmpObjvs=PopObj;
            tmpObjvs(i,:)=[];
            tmpHV=hv(tmpObjvs',refPoint);
            fitV(i,:)=totalHV-tmpHV;
        end
        [~,worst]=min(fitV);
    end
end
FrontNo = UpdateFront(PopObj,FrontNo,worst);
Population(worst) = [];
end