classdef ModelCSCS < HowToCSCS
    
    methods (Static)

        function mps = waveSpeed(isRadar)
            specificHeatAtSTP = 1.4;  % source NASA
            gasConstantForEarthsTroposphere = 286.0; % source NASA
            offsetFromCtoK = 273.15; % source NASA
            conversionFactorAtSTP=(specificHeatAtSTP*gasConstantForEarthsTroposphere)^.5;
            cMps=@(Tc) conversionFactorAtSTP*(Tc+offsetFromCtoK)^.5;  
            if(isRadar) 
                mps = 3.0e8;                % light speed
            else
                mps = cMps(24.0);           % soundspeed at temperature 24-deg C
            end
        end

        function [scale, B, T, Fs, Ts, N, c, H] = getDesignParameters( ...
            isRadar, ...
            modulationBandwidth, ...
            pulseDuration, ...
            frequencyMultiplier, ...
            arrayHeight ...
        )
            scale = 1e6;
            if(~isRadar)
                scale = 1e3;
            end
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
            if frequencyMultiplier>=2
                Fs = frequencyMultiplier*B;             % sampling frequency
            else
                Fs = 2.0*B;
            end
            Ts = 1/Fs;              % sampling time per grid
            N = T/Ts;               % grid points
            c = ModelCSCS.waveSpeed(isRadar);
            H = arrayHeight*1e-6*scale; % array height [m for radar, mm for sonar]
        end
        
        %% Additive White Gaussian Noise
        function noise_awgn = awgn(St, SNR)
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
        function htw_ks = kaiser_window(ht, alpha)
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
        
        %% Chirp set
        function [t, St, Sot, Sotdb,Ht] = chirp(K,T,N,doAWGN,SNR,doKaiser,alpha,chirpType)            
            A = 1;
            
            switch chirpType
                case -1
                    t = linspace(-T, 0, N);

                case 0
                    t = linspace(-0.5*T, 0.5*T, N);
                    
                case 1
                    t = linspace(0, T, N);                    
            end

            % -- General waveform -- %
            St = A*exp(1i*pi*K*t.^2);
            
            if doAWGN
                St = St + ModelCSCS.awgn(St, SNR);
            end
            
            % -- Matched Filter -- %
            t0 = linspace(-0.5*T, 0.5*T, N);
            Ht = exp(-1i*pi*K*t0.^2);
            
            if doKaiser
                Ht = ModelCSCS.kaiser_window(Ht, alpha);
            end
            
            % -- Convolution (s(t), h(t)) -- %
            switch chirpType
                case -1
                    Sot = conv(St, Ht);
                    Sot = Sot(end-length(St)+1:end);

                case 0
                    Sot = conv(St, Ht, 'Same');
                    
                case 1
                    Sot = conv(St, Ht);
                    Sot = Sot(1:length(St));
            end
            Sotn = Sot./max(Sot);
            Sotdb = 20*log10(abs(Sotn));
        end
        
        % Targets set
        function [tm, Srt, Stc, Sotc, Sotcdb] = targets(K,T,Ts,c,H,doAWGN,SNR,Ht,chirpType,angi,GR1)
            
            %% -- Multiple Targets -- %
            GR0 = H*tand(angi);            % ground range from Nadir to Near Range [m]
            
            if isempty(find(GR1, 1)) == 1
                GR1 = 0;
            else
                GR1 = GR1(GR1~=0);
            end
            
            % Ground Range between Nadir and Targets
            GR = GR0 + GR1;
            
            % Range between Radar and Targets.
            R = sqrt(GR.^2 + H^2);
            
            % range interval
            % time space before echo signals [s].
            switch chirpType
                case -1
                    Rmin = min(R) - c*T/1.2;
                    Rmax = max(R) + c*T/3;

                case 0
                    Rmin = min(R) - c*T/3;
                    Rmax = max(R) + c*T/3;
                    
                case 1
                    Rmin = min(R) - c*T/3;
                    Rmax = max(R) + c*T/1.2;                    
            end
            
            Rwid = Rmax - Rmin;    % receive window [meter]
            Twid = 2*Rwid/c;       % receive window [second]
            Nwid = ceil(Twid/Ts);  % receive window [number]
            
            % radar cross section
            
            RCS = ones(1, length(R));
            
            % -- Generate the echo -- %
            % receive window
            % open window when t = 2*Rmin/C
            % close window when t = 2*Rmax/C
            tm = linspace(2*Rmin/c, 2*Rmax/c, Nwid);
            
            % number of targets
            MN = length(R);
            
            % (total window time) - (received moment of each echo signals)
            % the received moment of each signals should be td = 0
            % for convenience, we set echo received time located at chirp center.
            t_echo_center = 2*R/c;
            td = repmat(tm, MN, 1) - repmat(t_echo_center', 1, Nwid);
            
            % radar echo from point targets
            switch chirpType
                case -1
                    Srt = exp(1i*pi*K*td.^2).*(td >= -T & td <= 0);

                case 0
                    Srt = exp(1i*pi*K*td.^2).*(abs(td)<=T/2);

                case 1
                    Srt = exp(1i*pi*K*td.^2).*(td >= 0 & td <= T);                    
            end

            if doAWGN
                size_srt = size(Srt);
                for ii = 1 : size_srt
                    Srt(ii, :) = Srt(ii,:) + ModelCSCS.awgn(Srt(ii, :),SNR);
                end
            end
            % combine all signals continuity
            Stc = RCS*Srt;
            
            % -- Convolution (s(t), h(t)) of Multiple Targets -- %
            switch chirpType
                case -1
                    Sotc = conv(Stc, Ht);
                    Sotc = Sotc(end-length(Stc)+1:end);

                case 0
                    Sotc = conv(Stc, Ht, 'Same');
                    
                case 1
                    Sotc = conv(Stc, Ht);
                    Sotc = Sotc(1:length(Stc));                    
            end
            Sotcn = Sotc/max(Sotc);
            Sotcdb = 20*log10(abs(Sotcn));
            
        end

        function [scale, Fsn, t, St, Sot, Sotdb, tm, Srt, Stc, Sotc, Sotcdb] = simulate(...
                isRadar, ...
                modulationBandwidth, ...
                pulseDuration, ...
                frequencyMultiplier, ...
                arrayHeight, ...
                doAWGN, ...
                SNR, ...
                doKaiser, ...
                alfa, ...
                chirpType, ...
                angi, ...
                GR1 ...
        )
            [scale, B, T, Fs, Ts, N, c, H] = ModelCSCS.getDesignParameters( ...
                isRadar,...
                modulationBandwidth, ...
                pulseDuration, ...
                frequencyMultiplier, ...
                arrayHeight ...
            );
            Fsn=Fs/B;
            K = B/T;
          
            [t, St, Sot, Sotdb, Ht] = ModelCSCS.chirp(K,T,N,doAWGN,SNR,doKaiser,alfa,chirpType);
            [tm, Srt, Stc, Sotc, Sotcdb] = ModelCSCS.targets(K,T,Ts,c,H,doAWGN,SNR,Ht,chirpType,angi,GR1);
        end            
    end

    methods (Access = public)

        % Construct model
        function model = ModelCSCS(isNotRadar)
            if nargin == 0
                isNotRadar = 0;
            end
            
            model.isRadar = ~isNotRadar;
            
            if nargout == 0
                clear model
            end
        end
    end
end