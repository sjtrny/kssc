%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Author: Xi Peng@milab.org, Sichuan University.
% pangsaai@gmail.com
% Date: 6, Apr. 2013
% Description:  This code is developed for large scale data clustering [1]. 

% IF your used any part of this code, PLEASE approximately cited our works;
% In addition, we adopted some codes from other works, Please approximately
%              cited the works if you used the corresponding codes or databases.

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

close all;
clear all;
clc;

%% --------------------------------------------------------------------------

addpath ('data/');
addpath('usages/');
% --- loading data
CurData = 'YaleB_48_42';
load (CurData);  
load ([CurData '_landmarkID']); 
% --- parameters configuration
par.nClass             =   length(unique(labels));
par                    =   L1ParameterConfig(par);

% --- dimension reduction using PCA
dat = FeatureEx(DATA, par);
clear DATA;

% --- split the data into two parts for landmark and non-landmark
Tr_dat = dat(:,landmark_ID(1:par.landmarkNO)); % the landmark data;
Tt_dat = dat(:,landmark_ID(1+par.landmarkNO:end)); % the non-landmark data;
Tr_gnd = labels(landmark_ID(1:par.landmarkNO));
Tt_gnd = labels(landmark_ID(1+par.landmarkNO:end));
labels = labels(landmark_ID);
clear dat;


fprintf([' \n--------- Running the experiment when lambda = ' num2str(par.lambda) ' | tolerance = ' num2str(par.tolerance) ' ---------\n ']);
% clustering the in-sample data
tic;
Tr_plabel = InSample(Tr_dat, par.lambda, par.tolerance, par, length(unique(Tr_gnd)));
time=toc;
% clustering the out-of-sample data
tic;
Tt_plabel = OutSample(Tr_dat, Tt_dat, Tr_plabel);
time=time+toc;
% Calculate the clustering quality in Accuracy and NMI
P_label = [Tr_plabel Tt_plabel]';
clear Classifier_plabel;
P_label = bestMap(labels,P_label);
labels = reshape(labels,[],1);
accuracy = length(find(labels == P_label))/length(labels);
nmi = MutualInfo(labels,P_label);

fprintf([' *** the accuracy scores is: ' num2str(accuracy) '\n']);
fprintf([' + the normalized mutual information score is: ' num2str(nmi) '\n']);
fprintf([' - The time cost is : ' num2str(time) '\n\n']);

clear Tr_dat Tt_dat ans fid i j;
clear Tr_gnd Tr_plabel Tt_gnd Tt_plabel labels landmark_ID P_label;
