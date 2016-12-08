%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Author: Xi PENG@milab.org, Sichuan University.
% pangsaai@gmail.com
% Date: 6, Apr. 2013
% Description:  This code is developed based on the code of SSC [3].

% IF your used any part of this code, PLEASE approximately cited the following works;

% Reference:
% [1] Xi Peng, Lei Zhang, Zhang Yi,
%     Scalable Sparse Subspace Clustering,
%     The 26th IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Portland, Oregon, USA, June, 2013.
% [2] Xi Peng, Lei Zhang, Zhang Yi, 
%     Constructing L2-Graph for Subspace Learning and Segmentation, 
%     arXiv:1209.0841.
% [3] E. Elhamifar and R. Vidal,
%     Sparse Subspace Clustering,
%     IEEE International Conference on Computer Vision and Pattern Recognition, 2009. 


% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL HOLDER AND CONTRIBUTORS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

%--------------------------------------------------------------------------
% This function takes a NxN coefficient matrix and returns a NxN adjacency
% matrix by choosing only the K strongest connections in the similarity
% graph
% CMat: NxN coefficient matrix
% K: number of strongest edges to keep; if K=0 use all the coefficients
% CKSym: NxN symmetric adjacency matrix
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function CKSym = BuildAdjacency(CMat,K)

N = size(CMat,1);
CAbs = abs(CMat);
for i = 1:N
    c = CAbs(:,i);
    [PSrt,PInd] = sort(c,'descend');
    if abs( c(PInd(1)) ) ~=0
        CAbs(:,i) = CAbs(:,i) ./ abs( c(PInd(1)) );
    else
        CAbs(:,i) = CAbs(:,i) ./ eps;
    end
end

CSym = CAbs + CAbs';

if (K ~= 0)
    [Srt,Ind] = sort( CSym,1,'descend' );
    CK = zeros(N,N);
    for i = 1:N
        for j = 1:K
            if CSym( Ind(1,i),i )~=0
                CK( Ind(j,i),i ) = CSym( Ind(j,i),i ) ./ eps;
            else
                CK( Ind(j,i),i ) = CSym( Ind(j,i),i ) ./ CSym( Ind(1,i),i );
            end
        end
    end
    CKSym = CK + CK';
else
    CKSym = CSym;
end