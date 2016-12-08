%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Author: Xi PENG@milab.org, Sichuan University.
% pangsaai@gmail.com
% Date: 6, Apr. 2013
% Description:  This code is developed for find the cluster assiginment of out-of-sample data using CRC [3]. 

% IF your used any part of this code, PLEASE approximately cited our works;

% Reference:
% [1] Xi Peng, Lei Zhang, Zhang Yi,
%     Scalable Sparse Subspace Clustering,
%     The 26th IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Portland, Oregon, USA, June, 2013.
% [2] Xi Peng, Lei Zhang, Zhang Yi, 
%     Constructing L2-Graph for Subspace Learning and Segmentation, 
%     arXiv:1209.0841.
% [3] Deng Cai, Xiaofei He, and Jiawei Han,
%     Document Clustering Using Locality Preserving Indexing, 
%     in IEEE TKDE, 2005.


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


function MIhat = MutualInfo(L1,L2)
%   mutual information
%
%   version 2.0 --May/2007
%   version 1.0 --November/2003
%
%   Written by ** (** AT gmail.com)
%===========    
L1 = L1(:);
L2 = L2(:);
if size(L1) ~= size(L2)
    error('size(L1) must == size(L2)');
end

Label = unique(L1);
nClass = length(Label);

Label2 = unique(L2);
nClass2 = length(Label2);
if nClass2 < nClass
     % smooth
     L1 = [L1; Label];
     L2 = [L2; Label];
elseif nClass2 > nClass
     % smooth
     L1 = [L1; Label2];
     L2 = [L2; Label2];
end


G = zeros(nClass);
for i=1:nClass
    for j=1:nClass
        G(i,j) = sum(L1 == Label(i) & L2 == Label(j));
    end
end
sumG = sum(G(:));

P1 = sum(G,2);  P1 = P1/sumG;
P2 = sum(G,1);  P2 = P2/sumG;
if sum(P1==0) > 0 || sum(P2==0) > 0
    % smooth
    error('Smooth fail!');
else
    H1 = sum(-P1.*log2(P1));
    H2 = sum(-P2.*log2(P2));
    P12 = G/sumG;
    PPP = P12./repmat(P2,nClass,1)./repmat(P1,1,nClass);
    PPP(abs(PPP) < 1e-12) = 1;
    MI = sum(P12(:) .* log2(PPP(:)));
    MIhat = MI / max(H1,H2);
    %%%%%%%%%%%%%   why complex ?       %%%%%%%%
    MIhat = real(MIhat);
end



