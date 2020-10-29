classdef ModelCSCS < handle
    properties (Access = protected)
        isSonar
    end
    
    methods (Access = protected)

        function mps = waveSpeed(app)
            specificHeatAtSTP = 1.4;  % source NASA
            gasConstantForEarthsTroposphere = 286.0; % source NASA
            offsetFromCtoK = 273.15; % source NASA
            conversionFactorAtSTP=(specificHeatAtSTP*gasConstantForEarthsTroposphere)^.5;
            cMps=@(Tc) conversionFactorAtSTP*(Tc+offsetFromCtoK)^.5;  
            %by definition 1hr = 3600s and 1 nautical mile = 1852 m
            %cKts=@(Tc) cMps(Tc)*3600/1852;  
            if(app.isSonar) 
                mps = cMps(24.0);           % soundspeed at temperature 24-deg C
            else
                mps = 3.0e8;                % light speed
            end
        end

        function [scale, B, T, Fs, Ts, N, c, H] = getDesignParameters( ...
            app, ...
            modulationBandwidth, ...
            pulseDuration, ...
            samplingFrequency, ...
            arrayHeight ...
        )
            scale = (app.isSonar)*1e3 + (1-app.isSonar)*1e6;
            if modulationBandwidth == 0    % chirp frequency modulation bandwidth
                B = 1.0*scale;
            else
                B = modulationBandwidth*scale;
            end
            if pulseDuration == 0    % pulse duration
                T = 1/scale;
            else
                T = pulseDuration/scale;
            end
            if samplingFrequency>=2.0*B
                Fs = samplingFrequency;             % sampling frequency
            else
                Fs=2.0*B;
            end
            Ts = 1/Fs;              % sampling time per grid
            N = T/Ts;               % grid points
            c = waveSpeed(app);
            H = arrayHeight*1e3/((app.isSonar)*1e6 + (1-app.isSonar)*1e3); % array height [m]
        end

        %% Additive White Gaussian Noise
        function noise_awgn = awgn(app, St, SNR)
            % y = awgn(x, SNR) adds AWGN noise vector to signal 'x' to generate 
            % a resulting signal vector y of specified SNR in dB
            
            % Set the random generator seed to default (for comparison only)
            %rng('default');
            
            L = length(St);
                        
            % Calculate actual symbol energy
            Esym = sum(abs(St).^2)/L;
            
            % Find the noise spectral density
            N0 = Esym/SNR;
            if(isreal(St))
                % Standard deviation for AWGN Noise when x is real
                noiseSigma = sqrt(N0);
                % Computed noise
                noise_awgn = noiseSigma*randn(1, L);
            else
                % Standard deviation for AWGN Noise when x is complex
                noiseSigma = sqrt(N0/2);
                % Computed noise
                noise_awgn = noiseSigma*(randn(1, L) + 1i*randn(1, L));
            end
    
        end
        
        %% Kaiser Window
        function htw_ks = kaiser_window(app, ht, alpha)
            % Kaiser window function for matched filter
            % GW
            % 2016-7-22
            
            % N points in s(t)
            N = length(ht);
            Ns = linspace(-0.5*N, 0.5*N, N);
            
            % Kaiser window
            % wt_ks = I0u/I0d,
            
            % I0 is the zeroth-order modified Bessel function of the first kind,
            % I0u = upper term, I0d = lower term,
            
            % alpha is an arbitrary, non-negative real number that determines ...
            % the shape of the window. In the frequency domain, it determines ...
            % the trade-off between main-lobe width and side lobe level, ...
            % which is a central decision in window design.
            
            I0u = besseli(0, pi*alpha*sqrt(1-((2*Ns)/N).^2));
            I0d = besseli(0, pi*alpha);
            wt_ks = I0u./I0d;
            htw_ks = ht.*wt_ks;
        end
        
    end

    methods (Access = public)

        % Construct app
        function app = ModelCSCS(isSonar)
            if nargin == 1
                app.isSonar = isSonar;
            else
                app.isSonar = 0;
            end
            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
        end
    end
end