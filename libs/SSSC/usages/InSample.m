%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Author: Xi PENG@milab.org, Sichuan University.
% pangsaai@gmail.com
% Date: 6, Apr. 2013
% Description:  This code is developed for clustering in-sample data. 

% IF your used any part of this code, PLEASE approximately cited our works;

% Reference:
% [1] Xi Peng, Lei Zhang, Zhang Yi,
%     Scalable Sparse Subspace Clustering,
%     The 26th IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Portland, Oregon, USA, June, 2013.
% [2] Xi Peng, Lei Zhang, Zhang Yi, 
%     Constructing L2-Graph for Subspace Learning and Segmentation, 
%     arXiv:1209.0841.

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
function Predict_label=InSample(Tr_dat, lambda, tolerance, par, nClass)

for i = 1:size(Tr_dat,2)
    if i == 1
        [tmp_x, tmp_iter] = SolveHomotopy(Tr_dat(:,2:end), Tr_dat(:,i), ...
            'maxIteration', par.maxIteration,...
            'isNonnegative', par.isNonnegative, ...
            'lambda', lambda, ...
            'tolerance', tolerance);
    else
        [tmp_x, tmp_iter] = SolveHomotopy([Tr_dat(:,1:i-1) Tr_dat(:,i+1:end)], Tr_dat(:,i), ...
            'maxIteration', par.maxIteration,...
            'isNonnegative', par.isNonnegative, ...
            'lambda', lambda, ...
            'tolerance', tolerance);
    end;
    % -- get the coeffient matrix    
    if i == 1
        coef = [0; tmp_x];
    else
        coef = [coef [tmp_x(1:i-1); 0; tmp_x(i:end)] ];
    end
end;
% --- Building Adiacency graph
CKSym = BuildAdjacency(coef,0);
% --- perform Normalized Symmetric spectral clustering
Predict_label = SC(CKSym,nClass);
Predict_label = reshape(Predict_label,1,[]);
