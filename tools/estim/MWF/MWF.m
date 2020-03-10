function [ e_MVDR_stack, e_MWF_inst_stack, e_MWF_sm_stack ] = MWF( y_stack, h, Gamma, phi_s_inst_stack, phi_d_inst_stack, phi_s_sm_stack, phi_d_sm_stack )
% runs the MWF in one frequency bin.
%
% IN:
% y_stack                 microphone signal - channels x frames
% h                       RETF - channels
% Gamma                   diffuse coherence matrix - channels x channels
% phi_s_inst_stack        instantaneous source image PSD
% phi_d_inst_stack        instantaneous diffuse PSD
% phi_s_sm_stack          smooth source image PSD
% phi_d_sm_stack          smooth diffuse PSD
%
% OUT:
% e_MVDR_stack            MVDR output - frames
% e_MWF_inst_stack        MWF output using instantaneous PSDs - frames
% e_MWF_sm_stack          MWF output using smooth PSDs - frames


% get dimensions
numFrames = size(y_stack,2);

% init output
e_MVDR_stack        = zeros(1, numFrames);
e_MWF_inst_stack    = zeros(1, numFrames);
e_MWF_sm_stack      = zeros(1, numFrames);

% MVDR
w_MVDR_num = Gamma\h;
w_MVDR_den = h'*w_MVDR_num;
w_MVDR = w_MVDR_num/w_MVDR_den;

for i_frame = 1:numFrames
    
    %%% Load Data
    %
    y = y_stack(:,i_frame);
    phi_s_inst = phi_s_inst_stack(i_frame);
    phi_d_inst = phi_d_inst_stack(i_frame);
    phi_s_sm = phi_s_sm_stack(i_frame);
    phi_d_sm = phi_d_sm_stack(i_frame);

    %%% Apply MVDR
    e_MVDR = w_MVDR'*y;
    
    %%% Spectral Post Processing
    %
    % compute gains
    gain_inst   = phi_s_inst/(phi_s_inst + phi_d_inst/w_MVDR_den);
    gain_sm     = phi_s_sm/(phi_s_sm + phi_d_sm/w_MVDR_den);

    %
    % apply gains
    e_MWF_inst    = gain_inst*e_MVDR;
    e_MWF_sm      = gain_sm*e_MVDR;

    %%% Save Data
    %
    % save
    e_MVDR_stack(:,i_frame)         = e_MVDR;
    e_MWF_inst_stack(:,i_frame)     = e_MWF_inst;
    e_MWF_sm_stack(:,i_frame)       = e_MWF_sm;

end

end
