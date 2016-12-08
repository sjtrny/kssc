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


function [newL2] = bestMap(L1,L2)
%bestmap: permute labels of L2 to match L1 as good as possible
%   [newL2] = bestMap(L1,L2);
%
%   version 2.0 --May/2007
%   version 1.0 --November/2003
%
%   Written by Deng Cai (dengcai AT gmail.com)


%===========    

L1 = L1(:);
L2 = L2(:);
if size(L1) ~= size(L2)
    error('size(L1) must == size(L2)');
end

Label1 = unique(L1);
nClass1 = length(Label1);
Label2 = unique(L2);
nClass2 = length(Label2);

nClass = max(nClass1,nClass2);
G = zeros(nClass);
for i=1:nClass1
	for j=1:nClass2
		G(i,j) = length(find(L1 == Label1(i) & L2 == Label2(j)));
	end
end

[c,t] = hungarian(-G);
newL2 = zeros(size(L2));
for i=1:nClass2
    newL2(L2 == Label2(i)) = Label1(c(i));
end

