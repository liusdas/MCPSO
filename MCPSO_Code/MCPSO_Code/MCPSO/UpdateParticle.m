function Offspring = UpdateParticle(Swarm,Pbest,Gbest)
% Particle swarm optimization in MCPSO

    Global	= GLOBAL.GetObj();
    N       = Global.N;
    D       = Global.D;
    ParticleVel = Swarm.adds(zeros(N,D));
    ParticlePosition = Swarm.decs;

    %% Parameter Settings
    W = 0.4;
    B1     = repmat(rand(N,1),1,D);
    B2     = repmat(rand(N,1),1,D);

    %% Update of Particle
    OffVel    = W.*ParticleVel + B1.*(Pbest-ParticlePosition) + B2.*(Gbest-ParticlePosition);
    Offspring = ParticlePosition + OffVel;

    
    %% Polynomial mutation
    proM=1;
    disM=20;
    
    Lower  = repmat(Global.lower,N,1);
    Upper  = repmat(Global.upper,N,1);

    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    
    Offspring = INDIVIDUAL(Offspring,OffVel);
end