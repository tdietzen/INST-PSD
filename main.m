%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2020 Thomas Dietzen
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt).
%
% If you find it useful, please cite:
% [1] T. Dietzen, M. Moonen, and T. van Waterschoot, 'Instantaneous PSD
% estimation for speech enhancement based on generalized principal
% components,' ESAT-STADIUS Tech. Rep., KU Leuven, Belgium, submitted for
% publication, Mar 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PREAMBLE

%clear;
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
set(0,'DefaultFigureWindowStyle','docked');


%% CONFIGURATION

%%% STFT settings
%
% sample rate
fs = 16000;
% STFT parameters
N_STFT = 512;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));
N_STFT_half = floor(N_STFT/2)+1;
% frequency vector
f = linspace(0,fs/2,N_STFT_half);

%%% ACOUSTIC SETTINGS
%
% speed of sound
c = 340;
% microphone positions
micPos = [...
    0, 0;...
    0.08, 0;...
    0.16, 0;...
    0.24, 0;...
    0.32, 0;...
    ];
% number of microphones
M = size(micPos,1);
% source angle
sourceAng = 30;
% SNR
SNR = 5;
% audio signals
x_TD = audioread('audio/x.wav');
s_TD = audioread('audio/s.wav');
v_TD = audioread('audio/v_SNR_0dB.wav');
v_TD = adjustSNR(x_TD, v_TD, SNR);
y_TD = x_TD + v_TD;

%%% ALGORITHMIC SETTINGS
%
% time constant and forgetting factor
tau   = 1;
zeta  = tau2forget(tau, R_STFT, fs);
% initial RETFs
h_FT = squeeze(doa2steervec(micPos, sourceAng, N_STFT_half, fs, c));
% diffuse coherence matrix
Gamma_FT = calc_diffcoherence(micPos,N_STFT,fs,c,1e-3);

%%% FIGURE SETTINGS
%
% spectogram figure settings
xTickProp = [0, R_STFT/fs, fs/R_STFT];
yTickProp = [0, fs/(2000*R_STFT), R_STFT/2];
cRange    = [-45 15];


%% STFT PROCESSING

y_STFT   = calc_STFT(y_TD, fs, win, N_STFT, R_STFT, 'onesided');
s_STFT   = calc_STFT(s_TD, fs, win, N_STFT, R_STFT, 'onesided');

% plot
figure('Name','microphone signal');
plotSpec(y_STFT(:,:,1),  'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
%
figure('Name','source image');
plotSpec(s_STFT(:,:,1), 'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
drawnow;


%% ESTIMATE PSDs

fprintf(' * estimate PSDs...\n');

% correlation matrix of microphone signal
Psi_y_STFT = estim_corrmat(y_STFT, zeta);
%
% compute GEVD
[P_STFT, lambda_inst_STFT, lambda_sm_STFT] = desmooth_GEVD(Psi_y_STFT, Gamma_FT, y_STFT);
%
% compute instantaneous PSD estimates
[phi_s_inst_hat, phi_d_inst_hat] = estim_PSD(P_STFT, lambda_inst_STFT, Gamma_FT, h_FT);
%
% compute smooth PSD estimates
[phi_s_sm_hat, phi_d_sm_hat] = estim_PSD(P_STFT, lambda_sm_STFT, Gamma_FT, h_FT);

% plot
figure('Name','source image PSD estim., smooth');
plotSpec(phi_s_sm_hat(:,:,1), 'pow', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
%
figure('Name','diffuse PSD estim., smooth');
plotSpec(phi_d_sm_hat(:,:,1), 'pow', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
%
figure('Name','source image PSD estim., instantaneous');
plotSpec(phi_s_inst_hat(:,:,1), 'pow', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
%
figure('Name','diffuse PSD estim., instantaneous');
plotSpec(phi_d_inst_hat(:,:,1), 'pow', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
drawnow;


%% MWF

fprintf(' * run MWF...\n');

% number of frames
numFrames = size(y_STFT,2);
%
% init outputs
e_MVDR_STFT        = zeros(N_STFT_half, numFrames);
e_MWF_inst_STFT    = zeros(N_STFT_half, numFrames);
e_MWF_sm_STFT      = zeros(N_STFT_half, numFrames);

% run per frequency bin
for k = 2:N_STFT_half
    
    % reorganize data
    y_stack              = shiftdim(squeeze(y_STFT(k,:,:)),1);                % from (1 x numFrames x M)       to (M x numFrames)
    h                    = h_FT(k,:).';                                       % from (1 x M)                   to M
    Gamma                = squeeze(Gamma_FT(k,1,:,:));                        % from (1 x 1 x M x M)           to (M x M)
    phi_s_inst_stack     = shiftdim(squeeze(phi_s_inst_hat(k,:,1)),1);        % from (1 x numFrames x 1)       to (1 x numFrames)
    phi_d_inst_stack     = shiftdim(squeeze(phi_d_inst_hat(k,:,1)),1);        % from (1 x numFrames x 1)       to (1 x numFrames)
    phi_s_sm_stack       = shiftdim(squeeze(phi_s_sm_hat(k,:,1)),1);          % from (1 x numFrames x 1)       to (1 x numFrames)
    phi_d_sm_stack       = shiftdim(squeeze(phi_d_sm_hat(k,:,1)),1);          % from (1 x numFrames x 1)       to (1 x numFrames)
    %
    % run MWF
    [e_MVDR_stack,...
        e_MWF_inst_stack,...
        e_MWF_sm_stack ] = ...
        MWF(...
        y_stack,...
        h,...
        Gamma,...
        phi_s_inst_stack,...
        phi_d_inst_stack,...
        phi_s_sm_stack,...
        phi_d_sm_stack );
    %
    % save output
    e_MVDR_STFT(k,:)            = e_MVDR_stack;
    e_MWF_inst_STFT(k,:)        = e_MWF_inst_stack;
    e_MWF_sm_STFT(k,:)          = e_MWF_sm_stack;
    
end

% plot
figure('Name','MVDR');
plotSpec(e_MVDR_STFT(:,:,1),        'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
%
figure('Name','MWF, smooth PSDs');
plotSpec(e_MWF_sm_STFT(:,:,1),      'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
%
figure('Name','MWF, inst PSDs');
plotSpec(e_MWF_inst_STFT(:,:,1),    'mag', xTickProp, yTickProp, cRange, 0); ylabel('f/kHz');
drawnow;


%% ISTFT PROCESSING

e_MVDR_TD           = calc_ISTFT(e_MVDR_STFT, win, N_STFT, R_STFT,      'onesided');
e_MWF_inst_TD       = calc_ISTFT(e_MWF_inst_STFT, win, N_STFT, R_STFT,  'onesided');
e_MWF_sm_TD         = calc_ISTFT(e_MWF_sm_STFT, win, N_STFT, R_STFT,    'onesided');


%% WRITE AUDIO

audiowrite(['.' filesep 'audio' filesep 'v_SNR_'            num2str(SNR) 'dB.wav'], v_TD            ,fs);
audiowrite(['.' filesep 'audio' filesep 'y_SNR_'            num2str(SNR) 'dB.wav'], y_TD            ,fs);
audiowrite(['.' filesep 'audio' filesep 'e_MVDR_SNR_'       num2str(SNR) 'dB.wav'], e_MVDR_TD       ,fs);
audiowrite(['.' filesep 'audio' filesep 'e_MWF_inst_SNR_'   num2str(SNR) 'dB.wav'], e_MWF_inst_TD   ,fs);
audiowrite(['.' filesep 'audio' filesep 'e_MWF_sm_SNR_'     num2str(SNR) 'dB.wav'], e_MWF_sm_TD     ,fs);

fprintf('\nDONE.\n');


%% HELPER FUNCTIONS

function [x2_scaled, scaling] = adjustSNR(x1, x2, SNR_x1_x2)
%
% scale signals
SingalPower = norm(x1(:,1));
NoisePower  = norm(x2(:,1));
if isinf(SNR_x1_x2)
    scaling = 0;
else
    scaling     = ((SingalPower/NoisePower)/db2mag(SNR_x1_x2));
end
x2_scaled   = scaling*x2;
end 
