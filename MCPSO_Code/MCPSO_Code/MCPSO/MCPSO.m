function MCPSO(Global)
% <algorithm> <M>
% MOPSO with Dual-Space Cooperative Directions
% K --- 10 --- Number of Clustering

%------------------------------- Reference --------------------------------
% Dual-Space Cooperative Directions for Multiobjective Particle Swarm
% Optimization

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter Settings
     K = Global.ParameterSet(10);

     %% Generate random swarm
     Swarm = Global.Initialization();
     FrontNo    = NDSort(Swarm.objs,inf); % NSorting
     
    %% Optimization
    while Global.NotTermination(Swarm)
        %% Spectral Clustering
        PopDec = Swarm.decs;
        DecKNNG = SimGraph_NearestNeighbors(PopDec',6,1, 1);   
        [ClusterLab,~] = SpectralClustering(DecKNNG, K, 3);
        %% Regularity Modeling
        Model = ReModeling(Swarm,ClusterLab,K);
        %% Fitness Assigning
        Fitness    = CalFitness(Swarm.objs);
        %% Get Pbest and Gbest
        Pbest=[];
        Gbest=[];
        for i = 1:Global.N
            % Pbest
            k = ClusterLab(i);
            Current = find(ClusterLab==k);
            CurFitness = Fitness(Current);
            PbestIndex = TournamentSelection(2,1,CurFitness);
            PbestDec = Swarm(Current(PbestIndex)).decs;
            Pbest = [Pbest;PbestDec];
            
            % Gbest
            if ~isempty(Model(k).eVector)
                lower = Model(k).a - 0.25*(Model(k).b-Model(k).a);
                upper = Model(k).b + 0.25*(Model(k).b-Model(k).a);
                trial = rand(1,Model(k).L).*(upper-lower) + lower;
                GbestDec = Model(k).mean + trial*Model(k).eVector(:,1:Model(k).L)';
            else
                GbestDec = Model(k).mean + randn(1,Global.D);
            end
            Gbest = [Gbest;GbestDec];
            
        end
        %% Update Particle and Swarm
        NewParticle = UpdateParticle(Swarm,Pbest,Gbest);
        for i = 1 : Global.N
            y = NewParticle(i);
            [Swarm,FrontNo] = Select(Swarm,FrontNo,y); % S-Metric Slection
        end
        % NSGA-II Selection
        %[Population,~,~] = EnvironmentalSelection([Swarm,NewParticle],Global.N); 
    end
end