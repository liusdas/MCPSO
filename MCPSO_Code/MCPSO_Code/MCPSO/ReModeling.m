function Model = ReModeling(Population,ClusterLab,K)
% Regularity Modeling by PCA 

    PopDec = Population.decs;
    [N,D] = size(PopDec);
    M      = length(Population(1).obj);
    Model = struct('mean',   num2cell(PopDec(1:K,:),2),...  % The mean of the model
                   'PI',     eye(D),...                     % The matrix PI
                   'eVector',[],...                         % The eigenvectors
                   'eValue', [],...                         % The eigenvalues
                   'L',      [],...                         % The the dimensionality of the model
                   'a',      [],...                         % The lower bound of the projections
                   'b',      []);                           % The upper bound of the projections
               
    %% Modeling
    for k = 1 : K
        current = find(ClusterLab==k);
        if length(current)<2
            Model(k).mean    = PopDec(current,:);
            Model(k).PI      = eye(D);
            Model(k).eVector = [];
            Model(k).eValue  = [];
        else
            Model(k).mean    = mean(PopDec(current,:),1);
            [eVector,eValue] = eig(cov(PopDec(current,:)-repmat(Model(k).mean,length(current),1)));
            [eValue,rank]    = sort(diag(eValue),'descend');
            Model(k).eValue  = real(eValue);
            Model(k).eVector = real(eVector(:,rank));
            Model(k).PI      = Model(k).eVector(:,M:end)*Model(k).eVector(:,M:end)';
            
            % Estimation of dimensionality
            % S = Model.eValue;
            % QS=sum(S);
            % explained = QS./sum(QS);
            % sum_explained = 0;
            % L = 0;
            % while sum_explained <= 0.95  
            %     L = L + 1;
            %     sum_explained = sum_explained + explained(L);
            % end

            Model(k).L     = M-1; % (m-1)-D model
            
            hyperRectangle = (PopDec(current,:)-repmat(Model(k).mean,length(current),1))*Model(k).eVector(:,1:Model(k).L);
            Model(k).a     = min(hyperRectangle,[],1);
            Model(k).b     = max(hyperRectangle,[],1);
        end
    end
end