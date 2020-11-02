% Implementation of a HowToCSCS class
%
% Matlab's implementation of multiple inheritance requires that if any
% superclass be a subclass of handle, then all superclasses must be
% subclasses of handle.  This requirement is dictated by Mathworks,
% the developers of the Matlab language.  As of Nov 2, 2020, the full text
% documenting this requirement can be found at:
%
% https://www.mathworks.com/help/matlab/matlab_oop/supporting-both-handle-and-value-subclasses-handlecompatible.html
%
% The default behaviour of the HowToCSCS provides properties
% that are typically associated with radars.  This behaviour is asserted by
% the publicly accessible boolean value "isRadar".
% The boolean value "isRadar" is asserted unless altered in a subclass.
%
% HowToCSCS is constructed without any input arguments.
% It returns a sample instance that demonstrates a default implementation.
% The instance is deleted unless returned to the caller.
%
% An instance of HowToCSCS has no implied fitness for use except to serve as 
% an example of datatypes returned from a static function named "simulate".
% 
% These datatypes could represent a transmitted radar signal and idealized 
% signal processing results in order to produce a mockup of signals.
% 
% A sample implmentation is provided for the static method "simulate". 
% Assertions are provided that reflect the original author's concept of 
% parameters from which signals can be generated.
%

classdef (HandleCompatible) HowToCSCS
    properties (Access = public)
        isRadar
    end
    
    methods (Access = public)
        % Construct model
        function model = HowToCSCS()
            model.isRadar = true;
            if nargout == 0
                clear model
            end
        end
    end

    methods (Static)

        function [scale, Fsn, t, St, Sot, Sotdb, tm, Srt, Stc, Sotc, Sotcdb] = simulate(...
                isRadar, ...
                modulationBandwidth, ...
                pulseDuration, ...
                frequencyMultiplier, ...
                arrayHeight, ...
                doAWGN, ...
                SNR, ...
                doKaiser, ...
                alpha, ...
                chirpType, ...
                angi, ...
                GR1 ...
        )
            assert(isRadar, "The base implementation is for Radar; no other behaviour is implemented")
            scale = 1e6;
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
            Fsn = Fs/B;
            K = B/T;
            Ts = 1/Fs;              % sampling time per grid
            N = T/Ts;               % grid points
            assert(~doAWGN,"additive white noise is not yet implemented");
            if(doAWGN)
               assert(0.01<=SNR && SNR<=2.0,"The value for SNR is expected to be between 0.01 and 2.0");
            end
            assert(~doKaiser,"Kaiser window is not yet implemented");
            if(doKaiser)
               assert(0.0<=alpha && alpha<=3.0,"The value for alpha is expected to be between 0 and 3");
            end
            switch chirpType
                case -1 % DOWNCHIRP
                    assert(false,"The base implementation is for Symmetric chirp; Downchirp is not implemented");                    

                case 0 % SYMMETRIC CHIRP
                    t = linspace(-0.5*T, 0.5*T, N);
                    
                case 1 % UPCHIRP
                    assert(false,"The base implementation is for Symmetric chirp; Upchirp is not implemented");                    
            end
            St = exp(1i*pi*K*t.^2);
            % -- Matched Filter -- %
            t0 = linspace(-0.5*T, 0.5*T, N);
            Ht = exp(-1i*pi*K*t0.^2);
            % S o t where o is the convolution of St with its impulse response
            Sot = conv(St, Ht, 'Same'); 
            Sotn = Sot./max(Sot);
            Sotdb = 20*log10(abs(Sotn));
            H=arrayHeight*1e-6*scale; % array height [m for radar, mm for sonar]            
            GR0 = H*tand(angi);            % ground range from Nadir to Near Range [m]            
            %% -- Multiple Targets -- %
            if isempty(find(GR1, 1)) == 1
                GR1 = 0;
            else
                GR1 = GR1(GR1~=0);
            end            
            % Ground Range between Nadir and Targets
            GR = GR0 + GR1;            
            % Range between Radar and Targets.
            R = sqrt(GR.^2 + H^2);
            c = 300.0*scale; %speed of light when scale is 1e6
            Rmin = min(R) - c*T/3;
            Rmax = max(R) + c*T/3;
            Rwid = Rmax - Rmin;    % receive window [meter]
            Twid = 2*Rwid/c;       % receive window [second]
            Nwid = ceil(Twid/Ts);  % receive window [number]
            tm = linspace(2*Rmin/c, 2*Rmax/c, Nwid);
            t_echo_center = 2*R/c;
            MN = length(R);
            td = repmat(tm, MN, 1) - repmat(t_echo_center', 1, Nwid);
            Srt = exp(1i*pi*K*td.^2).*(abs(td)<=T/2);
            RCS = ones(1, length(R));
            Stc = RCS*Srt;
            Sotc = conv(Stc, Ht, 'Same');
            Sotcn = Sotc/max(Sotc);
            Sotcdb = 20*log10(abs(Sotcn));
        end            
    end

end