function [p, K] = compute_group_difference(Group1, Group2, alpha, direction, N_all_edges)
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
%   K: vector array of fraction K as a function of frequency [1 x N freq.]


% Get parameters from Input arguments
N_parcels = size(Group1{1,1},1);    % Number of parcels
N_Fc = size(Group1,2);              % Number of frequency bands

% Initialize Output arguments
p{1,N_Fc} = [];         
K = zeros(1,N_Fc);


for f = 1:N_Fc
    f

    % Loop through connectivity matrix (edges) 
    for Parcel1 = 1:N_parcels
        for Parcel2 = 1:N_parcels
            edge_vector_Group1 = squeeze(Group1{1,f}(Parcel1,Parcel2,:));
            edge_vector_Group2 = squeeze(Group2{1,f}(Parcel1,Parcel2,:));

            % Wilcoxon rank sum test 
            p{1,f}(Parcel1,Parcel2) = ranksum(edge_vector_Group1, edge_vector_Group2, 'Tail', direction);

        end
    end    

    % Take number of significant edges in p matrix and divide by
    % number of all edges to get K
    K(f) = nnz(p{1,f} <= alpha)/N_all_edges;

end

end

