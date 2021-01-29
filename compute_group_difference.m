function [p, r, K, R] = compute_group_difference(Group1, Group2, alpha, direction, N_all_edges)
% COMPUTE_GROUP_DIFFERENCE 
% Computes statistical group difference by one-tailed Wilcoxon rank sum test
% 26/11/2020 Pauliina Yrjölä, BABA Center, Finland
%
%   INPUT ARGUMENTS
%   Group1: cell array {1 x N freq.} of connectivity matrices for all subjects [N parcels x N
%   parcels x N subj.]. Connectivity matrices per subject must be square matrices.
%   Group2: cell array {1 x N freq.} of connectivity matrices for all subjects [N parcels x N
%   parcels x N subj.]. Connectivity matrices per subject must be square matrices.
%   alpha: significance level
%   direction: direction of Wilcoxon rank sum test. Options: 
%       'right' (Group1 > Group2)
%       'left' (Group1 < Group2)
%   N_all_edges: number of edges for computing 
%       K = N_significant_edges/N_all_edges
%   
%   OUTPUT ARGUMENTS
%   p: cell array {1 x N freq.} of p-value matrices [N parcels x N parcels]
%   r: cell array {1 x N freq.} of r-value matrices [N parcels x N parcels]
%   K: vector array of fraction K as a function of frequency [1 x N freq.]
%   R: vector array of mean effect size [1 x N freq.]


% Get parameters from Input arguments
N_parcels = size(Group1{1,1},1);    % Number of parcels
N_Fc = size(Group1,2);              % Number of frequency bands

% Initialize Output arguments
p{1,N_Fc} = [];         
r{1,N_Fc} = [];         
K = zeros(1,N_Fc);
R = zeros(1,N_Fc);
S = zeros(1,N_Fc);

for f = 1:N_Fc
    f

    % Loop through connectivity matrix (edges) 
    for Parcel1 = 1:N_parcels
        for Parcel2 = 1:N_parcels
            edge_vector_Group1 = abs(squeeze(Group1{1,f}(Parcel1,Parcel2,:)));
            edge_vector_Group2 = abs(squeeze(Group2{1,f}(Parcel1,Parcel2,:)));

            % Wilcoxon rank sum test 
            [p{1,f}(Parcel1,Parcel2), ~, stats] = ranksum(edge_vector_Group1, edge_vector_Group2, 'Tail', direction);

            % Compute effect size
            N_1 = size(edge_vector_Group1,1);
            N_2 = size(edge_vector_Group2,1);
            
            % Convert W (Wilcoxon rank sum statistic) into Mann-Whitney U-test statistic
            W = stats.ranksum;
            U = W - (N_1*(N_1+1)/2);
            
            % Rank-biserial correlation
            r{1,f}(Parcel1,Parcel2) = 1-2*U/(N_1*N_2);
        end
    end    

    % Take number of significant edges in p matrix and divide by
    % number of all edges to get K
    K(f) = nnz(p{1,f} <= alpha)/N_all_edges;        
    
    % Compute mean effect size of the significant network
    R(f) = mean(r{1,f}(p{1,f} >= 0 & p{1,f} < alpha));
    
end

end

