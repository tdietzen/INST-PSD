function [P_STFT, lambda_inst_STFT, lambda_sm_STFT] = desmooth_GEVD(Psi_y_STFT, Gamma_FT, y_STFT, varargin) 
% [P_STFT, lambda_inst_STFT, lambda_sm_STFT] = desmooth_GEVD(Psi_x_STFT, Gamma_FT, varargin) 
% performs GEVD and computes instantaneous eigenvalues.
%
% IN:
% Psi_y_STFT                correlation matrix - freqbins x frames x channels x channels
% Gamma_FT                  diffuse coherence matrix - freqbins x 1 x channels x channels
% y_STFT                    microphone signal - freqbins x frames x channels
% 'lambdaMin', lambdaMin    eigenvalue threshold
%
% OUT:
% P_STFT                    eingevectors - freqbins x frames x channels x channels
% lambda_inst_STFT          instantaneous eigenvalues - freqbins x frames x channels
% lambda_sm_STFT            smooth eigenvalues - freqbins x frames x channels


% dimensions
[N_FT_half, L, M, ~]      = size(Psi_y_STFT);     % number of frequency bins, frames, microphones

% default options
lambdaMin  = 0;   % minimum power

% read options from input
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}                         
            case 'lambdaMin'
                lambdaMin     = varargin{i+1};              
        end
    end
end

% init
P_STFT                 = zeros(N_FT_half,L,M,M);
lambda_inst_STFT       = zeros(N_FT_half,L,M);
lambda_sm_STFT         = zeros(N_FT_half,L,M);

for k=2:N_FT_half
    
    % apply regularization
    Gamma = squeeze(Gamma_FT(k,1,:,:));
    for l = 1:L
        
        %%% GEVD %%%
        %
        Psi_y = squeeze(Psi_y_STFT(k,l,:,:));
        % generalized eigenvalue decomposition
        [P, Lambda_sm] = eig(Psi_y, Gamma);
        % ignore residual imaginary component and set negative values to zero
        lambda_sm = real(diag(Lambda_sm));
        % rescaling such that P'*Gamma*P = I and P'*Psi_y*P = Lambda
        P = P/diag(sqrt(real(diag(P'*Gamma*P))));
        %
        % get instantaneous eigenvalues
        y = squeeze(y_STFT(k,l,:));
        lambda_inst = abs(P'*y).^2;
        lambda_inst(lambda_inst < lambdaMin) = lambdaMin;
        %
        % save
        P_STFT(k,l,:,:)         = shiftdim(P, -2);
        lambda_sm_STFT(k,l,:)   = shiftdim(lambda_sm, -2);
        lambda_inst_STFT(k,l,:) = shiftdim(lambda_inst, -2);
        
    end               
end

end