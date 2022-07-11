function [p, r, K, R] = compute_group_difference(Group1, Group2, alpha, direction, N_edges)
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
%   N_edges: number of upper triangle edges for computing 
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

for f = 1:N_Fc
    f

    % Loop through connectivity matrix (edges) 
    for Parcel1 = 1:N_parcels
        for Parcel2 = 1:N_parcels
            edge_vector_Group1 = squeeze(Group1{1,f}(Parcel1,Parcel2,:));
            edge_vector_Group2 = squeeze(Group2{1,f}(Parcel1,Parcel2,:));

            % Wilcoxon rank sum test 
            [p{1,f}(Parcel1,Parcel2), ~, stats] = ranksum(edge_vector_Group1, edge_vector_Group2, 'Tail', direction);

            % Compute effect size
            N_1 = size(edge_vector_Group1,1);
            N_2 = size(edge_vector_Group2,1);
            N = N_1+N_2;
            
            % Convert W into Mann-Whitney U-test statistic and take smaller U
            W_1 = stats.ranksum; % rank sum statistic
            W_2 = (N*(N+1)/2)- W_1; % rank sum statistic
            U_1 = W_1 - (N_1*(N_1+1)/2); % Convert W (Wilcoxon rank sum statistic) into Mann-Whitney U-test statistic   
            U_2 = W_2 - (N_2*(N_2+1)/2); % Convert W (Wilcoxon rank sum statistic) into Mann-Whitney U-test statistic 
            U = min([U_1, U_2]); % take smaller U
            
            % Rank-biserial correlation
            r{1,f}(Parcel1,Parcel2) = 1-2*U/(N_1*N_2);
        end
    end    

    % Take number of significant edges in p matrix and divide by
    % number of all edges to get K
    p_significant = p{1,f}.*(p{1,f} < alpha);
    K(f) = nnz(triu(p_significant))/N_edges; 
    
    % Compute mean effect size of the significant network
    r_significant = r{1,f}.*(p{1,f} < alpha);
    R(f) = sum(sum(triu(r_significant)))/nnz(triu(r_significant));
    
end

end
