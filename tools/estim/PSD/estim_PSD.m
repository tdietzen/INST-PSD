function [phi_s_STFT, phi_d_STFT] = estim_PSD(P_STFT, lambda_STFT, Gamma_FT, h_FT, varargin) 
% function [phi_s_STFT, phi_d_STFT] = estim_PSD(P_STFT, lambda_STFT, Gamma_FT, h_FT, varargin) 
% estimates PSDs.
%
% IN:
% P_STFT                    eigenvectors - freqbins x frames x channels x channels
% lambda_STFT               eigenvalues - freqbins x frames x channels
% Gamma_FT                  diffuse coherence matrix - freqbins x 1 x frames x channels
% h_FT                      RETF estimate - freqbins x channels
% 'phiMin', phiMin          power threshold
%
% OUT:
% phi_s_STFT                early PSD estimate - freqbins x frames
% phi_d_STFT                diffuse PSD estimate - freqbins x frames

% dimensions
[N_FT_half, L, ~]      = size(lambda_STFT);             % number of frequency bins, frames, microphones

% default values
phiMin                    = 0;                          % minimum power threshold

% parse name-value pairs
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        switch varargin{i} 
            case 'phiMin'
                phiMin      = varargin{i+1};          
        end
    end
end

% init
phi_s_STFT        = zeros(N_FT_half,L);
phi_d_STFT        = zeros(N_FT_half,L);


for k=2:N_FT_half
    
    % get h
    h  = h_FT(k,:).';                                   
    
    % get Gamma
    Gamma = squeeze(Gamma_FT(k,1,:,:));
    
    for l = 1:L
               
        % get P and lambda
        P = squeeze(P_STFT(k,l,:,:));
        lambda = squeeze(lambda_STFT(k,l,:));
        %
        % find index of maximum eigenvalues
        [~, maxIdx] = sort(lambda, 'descend');
        p1 = P(:,maxIdx(1)); 
        lambda_1 = lambda(maxIdx(1));
        lambda_m = lambda(maxIdx(2:end));
        %
        % compute late reverberant PSD
        phi_d = mean(lambda_m); 
        %
        % save
        phi_d_STFT(k,l) = phi_d;
        %
        % compute lambda_xe
        lambda_xe = lambda_1 - phi_d;

        % compute PSD
        phi_s_hat = lambda_xe*abs(h'*Gamma*p1/(h'*h))^2;
        phi_s_hat = max(phiMin, phi_s_hat);
        
        % save
        phi_s_STFT(k,l) = phi_s_hat;

    end   
end