function [p, r, K, R] = compute_clinical_correlation(Group, Neuroscore, covariate, alpha, N_all_edges)
% COMPUTE_CLINICAL_CORRELATION
% Computes edgewise correlation of connectivity strength to clinical
% neuroscores with Spearman's rank correlation
% 27/11/2020 Pauliina Yrjölä, BABA Center, Finland
%
%   INPUT ARGUMENTS
%   Group: cell array {1 x N freq.} of connectivity matrices for all subjects [N parcels x N
%   parcels x N subj.]. Connectivity matrices per subject must be square matrices.
%   Neuroscore: vector array [1 x N subj.] of clinical neuroscores for each subject.
%   covariate: vector array [1 x N subj.] of variable (eg. age) to use as a
%   covariate. If no covariate, give [].
%   alpha: significance level
%   N_all_edges: number of edges for computing 
%       K = N_significant_edges/N_all_edges
%   
%   OUTPUT ARGUMENTS
%   p: cell array {1 x N freq.} of p-value matrices [N parcels x N parcels]
%   r: cell array {1 x N freq.} of r-value matrices [N parcels x N parcels]
%   K: vector array of fraction K as a function of frequency [2 x N freq.]
%       K(1,:) -> fraction of positive correlation (r >= 0)
%       K(2,:) -> fraction of negative correlation (r < 0)
%   R: vector array of the mean effect size of the significant network as a function of frequency [2 x N freq.]
%       R(1,:) -> effect size of positive correlation (r >= 0)
%       R(2,:) -> effect size of negative correlation (r < 0)


% Get parameters from Input arguments
N_parcels = size(Group{1,1},1);    % Number of parcels
N_Fc = size(Group,2);              % Number of frequency bands

% Initialize Output arguments
p{1,N_Fc} = [];
r{1,N_Fc} = [];
K = zeros(2,N_Fc);
R = zeros(2,N_Fc);

   
for f = 1:N_Fc
    for Parcel1 = 1:N_parcels
        for Parcel2 = 1:N_parcels   

            % Extract edge vector
            edge_vector_connectivity = abs(squeeze(Group{1,f}(Parcel1,Parcel2,:)));
            
            % Compute Spearman correlation
            if isempty(covariate)
                [r{1,f}(Parcel1,Parcel2), p{1,f}(Parcel1,Parcel2)] = corr(edge_vector_connectivity, Neuroscore, 'type', 'Spearman');
            else
                [r{1,f}(Parcel1,Parcel2), p{1,f}(Parcel1,Parcel2)] = partialcorr(edge_vector_connectivity, Neuroscore, covariate, 'type', 'Spearman');
            end
            
        end
    end 

    % Compute number of significant correlations and divide by
    % number of all connections            
    K(1,f) = nnz(p{1,f} < alpha & r{1,f} >= 0)/N_all_edges;
    K(2,f) = nnz(p{1,f} < alpha & r{1,f} < 0)/N_all_edges;
    
    % Compute mean effect size of the significant network
    R(1,f) = mean(r{1,f}(r{1,f} >= 0 & p{1,f} >= 0 & p{1,f} < alpha));
    R(2,f) = mean(r{1,f}(r{1,f} < 0 & p{1,f} >= 0 & p{1,f} < alpha));
end

end

