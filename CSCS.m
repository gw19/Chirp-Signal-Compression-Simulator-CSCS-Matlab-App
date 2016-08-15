classdef CSCS < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure                   % Single ...
        UISot                  matlab.ui.control.UIAxes           % Signal ...
        UISotdb                matlab.ui.control.UIAxes           % Signal ...
        UISt                   matlab.ui.control.UIAxes           % Radar E...
        UIStc                  matlab.ui.control.UIAxes           % Combina...
        TabGroup               matlab.ui.container.TabGroup       % Chirp S...
        Tab                    matlab.ui.container.Tab            % Chirp Set
        LabelSlider            matlab.ui.control.Label            % Bandwid...
        Bw                     matlab.ui.control.Slider           % [0 100]
        disp_B                 matlab.ui.control.NumericEditField % [1 100]
        Label                  matlab.ui.control.Label            % Pulse d...
        Tp                     matlab.ui.control.Slider           % [1 50]
        disp_T                 matlab.ui.control.NumericEditField % [1 50]
        Label2                 matlab.ui.control.Label            % Samplin...
        Fsp                    matlab.ui.control.Slider           % [2 12]
        disp_Fs                matlab.ui.control.NumericEditField % [2 12]
        Label3                 matlab.ui.control.Label            % x B
        Button_Run             matlab.ui.control.Button           % Run
        LabelKnob              matlab.ui.control.Label            % SNR = 
        KnobAWGN               matlab.ui.control.Knob             % [1.0000...
        dispAWGN               matlab.ui.control.NumericEditField % [1.0000...
        Label4                 matlab.ui.control.Label            % Kaiser ...
        swKaiser               matlab.ui.control.Switch           % Off, On
        Label5                 matlab.ui.control.Label            % Alpha = 
        KnobKaiser             matlab.ui.control.Knob             % [0 3]
        dispKaiser             matlab.ui.control.NumericEditField % [0 3]
        LampAutoUpdate         matlab.ui.control.Lamp            
        swAutoUpdate           matlab.ui.control.Switch           % Off, On
        Label6                 matlab.ui.control.Label            % Auto Run
        LabelSwitch            matlab.ui.control.Label            % White G...
        swAWGN                 matlab.ui.control.Switch           % Off, On
        LampAWGN               matlab.ui.control.Lamp            
        LampKaiser             matlab.ui.control.Lamp            
        LabelDropDown          matlab.ui.control.Label            % Waveform
        swChirp                matlab.ui.control.DropDown         % Up-Chir...
        Tab2                   matlab.ui.container.Tab            % Targets
        Label28                matlab.ui.control.Label            % Distanc...
        swTarget1              matlab.ui.control.ToggleSwitch     % Off, On
        LampTarget1            matlab.ui.control.Lamp            
        swTarget2              matlab.ui.control.ToggleSwitch     % Off, On
        LampTarget2            matlab.ui.control.Lamp            
        swTarget3              matlab.ui.control.ToggleSwitch     % Off, On
        LampTarget3            matlab.ui.control.Lamp            
        swTarget4              matlab.ui.control.ToggleSwitch     % Off, On
        LampTarget4            matlab.ui.control.Lamp            
        swTarget5              matlab.ui.control.ToggleSwitch     % Off, On
        LampTarget5            matlab.ui.control.Lamp            
        swTarget6              matlab.ui.control.ToggleSwitch     % Off, On
        LampTarget6            matlab.ui.control.Lamp            
        Panel                  matlab.ui.container.Panel         
        Label27                matlab.ui.control.Label            % m
        Label26                matlab.ui.control.Label            % m
        Label25                matlab.ui.control.Label            % m
        Label24                matlab.ui.control.Label            % m
        Label23                matlab.ui.control.Label            % m
        Label22                matlab.ui.control.Label            % m
        Label21                matlab.ui.control.Label            % Target 6
        InputTarget6           matlab.ui.control.NumericEditField % [1 Inf]
        Label20                matlab.ui.control.Label            % Target 5
        InputTarget5           matlab.ui.control.NumericEditField % [1 Inf]
        Label19                matlab.ui.control.Label            % Target 4
        InputTarget4           matlab.ui.control.NumericEditField % [1 Inf]
        Label18                matlab.ui.control.Label            % Target 3
        InputTarget3           matlab.ui.control.NumericEditField % [1 Inf]
        LabelNumericEditField2 matlab.ui.control.Label            % Target 2
        InputTarget2           matlab.ui.control.NumericEditField % [1 Inf]
        LabelNumericEditField  matlab.ui.control.Label            % Target 1
        InputTarget1           matlab.ui.control.NumericEditField % [1 Inf]
        LabelSlider2           matlab.ui.control.Label            % Nadir
        gr0                    matlab.ui.control.Slider           % [50 1000]
        Label29                matlab.ui.control.Label            % Near Range
        Label30                matlab.ui.control.Label            % Targets
        Lamp                   matlab.ui.control.Lamp            
        Lamp2                  matlab.ui.control.Lamp            
        Lamp3                  matlab.ui.control.Lamp            
        Lamp4                  matlab.ui.control.Lamp            
        Lamp5                  matlab.ui.control.Lamp            
        Lamp6                  matlab.ui.control.Lamp            
        LabelSlider3           matlab.ui.control.Label            % Height
        h                      matlab.ui.control.Slider           % [400 1000]
        LabelNumericEditField3 matlab.ui.control.Label            % Angle o...
        Theta                  matlab.ui.control.NumericEditField % [10 45]
        Label31                matlab.ui.control.Label            % [degree]
        disp_NearRange         matlab.ui.control.NumericEditField % [70 1000]
        Label32                matlab.ui.control.Label            % [km]
        Label33                matlab.ui.control.Label            % [km]
        disp_h                 matlab.ui.control.NumericEditField % [400 1000]
        KnobTheta              matlab.ui.control.Knob             % [10 45]
        Label34                matlab.ui.control.Label            % Set Rad...
        Tab3                   matlab.ui.container.Tab            % Hints
        Label13                matlab.ui.control.Label            % Author:...
        Label14                matlab.ui.control.Label            % Nationa...
        Label15                matlab.ui.control.Label            % Graduat...
        Label16                matlab.ui.control.Label            % gw.chen...
        Label17                matlab.ui.control.Label            % Version...
        Label35                matlab.ui.control.Label            % Step 1:...
        TextArea               matlab.ui.control.TextArea         % First, ...
        Label36                matlab.ui.control.Label            % Step 2:...
        TextArea2              matlab.ui.control.TextArea         % Second,...
        Label7                 matlab.ui.control.Label            % ¡X¡X¡X¡X¡X¡X ...
        UIStm                  matlab.ui.control.UIAxes           % Multipl...
    end

    % GW
    % 2016-8-5
    
    properties (Access = private)
        autoUpdate
    end

    methods (Access = private)
        
        function updatePlot(app)
            
            if app.Bw.Value == 0    % chirp frequency modulation bandwidth
                B = 1*1e6;
            else
                B = app.Bw.Value*1e6;
            end
            if app.Tp.Value == 0    % pulse duration
                T = 1*1e-6;
            else
                T = app.Tp.Value*1e-6;
            end
            Fsn = app.Fsp.Value;
            Fs = Fsn*B;             % sampling frequency
            Ts = 1/Fs;              % sampling time per grid
            N = T/Ts;               % grid points
            c = 3e8;                % light speed

            t = linspace(0, T, N);
            
            [St, Sot, Sotdb, tm, Srt, Stc, Sotc, Sotcdb] = chirp(app);
            
            % -- Plot figures -- %
            % s(t)
            plot(app.UISt, t*1e6, real(St));
            xlim(app.UISt, [t(1)*1e6, t(end)*1e6]);
            ylim(app.UISt, [-max(abs(St))-0.5, max(abs(St))+0.5]);
            
            % s_output(t)
            plot(app.UISot, t*1e6, real(Sot)./Fsn);
            axis(app.UISot, 'tight');

            % s_output(t) in [dB]
            tt1 = min(t(Sotdb>=-13));
            tt2 = max(t(Sotdb>=-13));
            ttd = tt2 - tt1;
            zoomin = 0.2;
            plot(app.UISotdb, t*1e6, Sotdb);
            ylabel(app.UISotdb, 'Sout(t)');
            axis(app.UISotdb, 'tight');
            xlim(app.UISot, [(tt1 - ttd/zoomin)*1e6, (tt2 + ttd/zoomin)*1e6]);
            xlim(app.UISotdb, [(tt1 - ttd/zoomin)*1e6, (tt2 + ttd/zoomin)*1e6]);
            %xlim(app.UISotdb, [(t(1)+zoomin*T)*1e6, (t(end)-zoomin*T)*1e6]);
            ylim(app.UISotdb, [-30, 0]);
            
            target_number = 6;
            swTarget_Value = zeros(1, target_number);
            for ii = 1 : target_number
                swTarget_Value_index = eval(['app.swTarget', num2str(ii), '.Value']);
                if strcmp(swTarget_Value_index, 'On')
                    swTarget_Value(ii) = 1;
                end
            end
            
            if sum(swTarget_Value) >= 2
                hold(app.UIStm, 'off');
                for ii = 1 : sum(swTarget_Value)
                    plot(app.UIStm, tm*1e6, real(Srt(ii,:)-2*(ii-1)));
                    hold(app.UIStm, 'on');
                end
                axis(app.UIStm, 'tight');
                %xlim(app.UIStm, [tm(1)*1e6, tm(end)*1e6]);
                
                plot(app.UIStc, tm*1e6, real(Stc));
                axis(app.UIStc, 'tight');
                %xlim(app.UIStc, [tm(1)*1e6, tm(end)*1e6]);
                ylim(app.UIStc, [-max(abs(Stc))-0.5, max(abs(Stc))+0.5]);
                xlabel(app.UIStc, 'Time [\mus]');
                
                % s_output(t)
                ttm1 = min(tm(Sotcdb>=-13));
                ttm2 = max(tm(Sotcdb>=-13));
                ttmd = ttm2 - ttm1;
                zoomin = 0.3;
                
                plot(app.UISot, tm*1e6, real(Sotc)./Fsn);
                axis(app.UISot, 'tight');
                xlim(app.UISotdb, [(ttm1 - ttmd/zoomin)*1e6, (ttm2 + ttmd/zoomin)*1e6]);
                xlabel(app.UISot, 'Time [\mus]');
                
                % s_output(t) in [dB]
                plot(app.UISotdb, tm*1e6, Sotcdb);
                ylabel(app.UISotdb, 'S_o_u_t(t)');
                axis(app.UISotdb, 'tight');
                xlim(app.UISot, [(ttm1 - ttmd/zoomin)*1e6, (ttm2 + ttmd/zoomin)*1e6]);
                xlim(app.UISotdb, [(ttm1 - ttmd/zoomin)*1e6, (ttm2 + ttmd/zoomin)*1e6]);
                ylim(app.UISotdb, [-30, 0]);
                xlabel(app.UISotdb, 'Time [\mus]');
                
            else
                hold(app.UIStm, 'off');
                % each s(t)
                plot(app.UIStm, t*1e6, real(St));
                xlim(app.UIStm, [t(1)*1e6, t(end)*1e6]);
                ylim(app.UIStm, [-max(abs(St))-0.5, max(abs(St))+0.5]);
                
                % combined s(t)
                plot(app.UIStc, t*1e6, real(St));
                xlim(app.UIStc, [t(1)*1e6, t(end)*1e6]);
                ylim(app.UIStc, [-max(abs(St))-0.5, max(abs(St))+0.5]);
            end
        end

        %% Chirp set
        function [St, Sot, Sotdb, tm, Srt, Stc, Sotc, Sotcdb] = chirp(app)
            
            if app.Bw.Value == 0    % chirp frequency modulation bandwidth
                B = 1*1e6;
            else
                B = app.Bw.Value*1e6;
            end
            if app.Tp.Value == 0    % pulse duration
                T = 0.1*1e-6;
            else
                T = app.Tp.Value*1e-6;
            end
            Fsn = app.Fsp.Value;
            Fs = Fsn*B;             % sampling frequency
            Ts = 1/Fs;              % sampling time per grid
            N = T/Ts;               % grid points
            c = 3e8;                % light speed
            K = B/T;
            A = 1;
            
            switch app.swChirp.Value
                case 'Symmetry-Chirp'
                    t = linspace(-0.5*T, 0.5*T, N);
                    
                case 'Up-Chirp'
                    t = linspace(0, T, N);
                    
                case 'Down-Chirp'
                    t = linspace(-T, 0, N);
            end

            % -- General waveform -- %
            St = A*exp(1i*pi*K*t.^2);
            
            if strcmp(app.swAWGN.Value, 'On')
                St = awgn(app, St);
            end
            
            % -- Matched Filter -- %
            t0 = linspace(-0.5*T, 0.5*T, N);
            Ht = exp(-1i*pi*K*t0.^2);
            
            if strcmp(app.swKaiser.Value, 'On')
                Ht = kaiser_window(app, Ht);
            end
            
            % -- Convolution (s(t), h(t)) -- %
            switch app.swChirp.Value
                case 'Symmetry-Chirp'
                    Sot = conv(St, Ht, 'Same');
                    
                case 'Up-Chirp'
                    Sot = conv(St, Ht);
                    Sot = Sot(1:length(St));
                    
                case 'Down-Chirp'
                    Sot = conv(St, Ht);
                    Sot = Sot(end-length(St)+1:end);
            end
            Sotn = Sot./max(Sot);
            Sotdb = 20*log10(abs(Sotn));
            
            
            %% -- Multiple Targets -- %
            H = app.h.Value*1e3;           % radar height [m]
            angi = app.Theta.Value;        % angle of incidence [degree]
            GR0 = H*tand(angi);            % ground range from Nadir to Near Range [m]
            target_number = 6;
            GR1 = zeros(1, target_number);
            for ii = 1 : target_number
                swTarget_Value = eval(['app.swTarget', num2str(ii), '.Value']);
                if strcmp(swTarget_Value, 'On')
                    GR1(ii) = eval(['app.InputTarget', num2str(ii), '.Value']);
                end
            end
            
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
            switch app.swChirp.Value
                case 'Symmetry-Chirp'
                    Rmin = min(R) - c*T/3;
                    Rmax = max(R) + c*T/3;
                    
                case 'Up-Chirp'
                    Rmin = min(R) - c*T/3;
                    Rmax = max(R) + c*T/1.2;
                    
                case 'Down-Chirp'
                    Rmin = min(R) - c*T/1.2;
                    Rmax = max(R) + c*T/3;
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
            switch app.swChirp.Value
                case 'Symmetry-Chirp'
                    Srt = exp(1i*pi*K*td.^2).*(abs(td)<=T/2);

                case 'Up-Chirp'
                    Srt = exp(1i*pi*K*td.^2).*(td >= 0 & td <= T);
                    
                case 'Down-Chirp'
                    Srt = exp(1i*pi*K*td.^2).*(td >= -T & td <= 0);
            end

            if strcmp(app.swAWGN.Value, 'On')
                size_srt = size(Srt);
                for ii = 1 : size_srt
                    Srt(ii, :) = awgn(app, Srt(ii, :));
                end
            end
            % combine all signals continuity
            Stc = RCS*Srt;
            
            % -- Convolution (s(t), h(t)) of Multiple Targets -- %
            switch app.swChirp.Value
                case 'Symmetry-Chirp'
                    Sotc = conv(Stc, Ht, 'Same');
                    
                case 'Up-Chirp'
                    Sotc = conv(Stc, Ht);
                    Sotc = Sotc(1:length(Stc));
                    
                case 'Down-Chirp'
                    Sotc = conv(Stc, Ht);
                    Sotc = Sotc(end-length(Stc)+1:end);
            end
            Sotcn = Sotc/max(Sotc);
            Sotcdb = 20*log10(abs(Sotcn));
            
        end
        
        %% Additive White Gaussian Noise
        function St_awgn = awgn(app, St)
            % y = awgn(x, SNR) adds AWGN noise vector to signal 'x' to generate 
            % a resulting signal vector y of specified SNR in dB
            
            % Set the random generator seed to default (for comparison only)
            %rng('default');
            
            L = length(St);
            
            % SNR to linear scale
            SNR = app.KnobAWGN.Value;
            
            % Calculate actual symbol energy
            Esym = sum(abs(St).^2)/L;
            
            % Find the noise spectral density
            N0 = Esym/SNR;
            if(isreal(St)), 
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
    
            % Received signal
            St_awgn = St + noise_awgn;
        end
        
        %% Kaiser Window
        function htw_ks = kaiser_window(app, ht)
            % Kaiser window function for matched filter
            % GW
            % 2016-7-22

            % Parameters
            % Alpha
            alfa = app.KnobKaiser.Value;
            
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
            
            I0u = besseli(0, pi*alfa*sqrt(1-((2*Ns)/N).^2));
            I0d = besseli(0, pi*alfa);
            wt_ks = I0u./I0d;
            htw_ks = ht.*wt_ks;
        end

    end


    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.swAutoUpdate.Value = 'On';
            app.swAWGN.Value = 'Off';
            app.LampAWGN.Color = [0.5, 0.5, 0.5];
            app.KnobAWGN.Enable = 'Off';
            app.dispAWGN.Enable = 'off';
            app.swKaiser.Value = 'Off';
            app.LampKaiser.Color = [0.5, 0.5, 0.5];
            app.KnobKaiser.Enable = 'Off';
            app.dispKaiser.Enable = 'off';
            app.Button_Run.Visible = 'Off';
            app.swChirp.Value = 'Symmetry-Chirp';
            app.InputTarget1.Enable = 'Off';
            app.InputTarget2.Enable = 'Off';
            app.InputTarget3.Enable = 'Off';
            app.InputTarget4.Enable = 'Off';
            app.InputTarget5.Enable = 'Off';
            app.InputTarget6.Enable = 'Off';
            app.LampTarget1.Color = [0.5, 0.5, 0.5];
            app.LampTarget2.Color = [0.5, 0.5, 0.5];
            app.LampTarget3.Color = [0.5, 0.5, 0.5];
            app.LampTarget4.Color = [0.5, 0.5, 0.5];
            app.LampTarget5.Color = [0.5, 0.5, 0.5];
            app.LampTarget6.Color = [0.5, 0.5, 0.5];
            app.autoUpdate = 1;
            updatePlot(app)
        end

        % Bw value changed function
        function BwValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % Bw value changing function
        function BwValueChanging(app, event)
            format('Bank');
            if event.Value == 0
                app.disp_B.Value = 1;
            else
                app.disp_B.Value = event.Value;
            end
        end

        % Button_Run button pushed function
        function Button_RunButtonPushed(app)
            updatePlot(app)
        end

        % Tp value changed function
        function TpValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % disp_B value changed function
        function disp_BValueChanged(app)
            app.Bw.Value = app.disp_B.Value;
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % disp_T value changed function
        function disp_TValueChanged(app)
            app.Tp.Value = app.disp_T.Value;
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % Tp value changing function
        function TpValueChanging(app, event)
            format('Bank')
            if event.Value == 0
                app.disp_T.Value = 1;
            else
                app.disp_T.Value = event.Value;
            end
        end

        % disp_Fs value changed function
        function disp_FsValueChanged(app)
            app.Fsp.Value = app.disp_Fs.Value;
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % Fsp value changing function
        function FspValueChanging(app, event)
            format('Bank')
            app.disp_Fs.Value = event.Value;
        end

        % swAWGN value changed function
        function swAWGNValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swAWGN.Value, 'On')
                app.KnobAWGN.Enable = 'on';
                app.dispAWGN.Enable = 'on';
                app.LampAWGN.Color = [0, 1, 0] ;
            else
                app.KnobAWGN.Enable = 'off';
                app.dispAWGN.Enable = 'off';
                app.LampAWGN.Color = [0.5, 0.5, 0.5] ;  
            end      
        end

        % swAutoUpdate value changed function
        function swAutoUpdateValueChanged(app)
            if strcmp(app.swAutoUpdate.Value, 'On')
                app.autoUpdate = 1;
                app.Button_Run.Visible = 'Off' ;
                app.LampAutoUpdate.Color = [0, 1, 0] ;
            else
                app.autoUpdate = 0;
                app.Button_Run.Visible = 'On' ;
                app.LampAutoUpdate.Color = [0.5, 0.5, 0.5] ;  
            end
            
        end

        % Fsp value changed function
        function FspValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % KnobAWGN value changed function
        function KnobAWGNValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % dispAWGN value changed function
        function dispAWGNValueChanged(app)
            app.KnobAWGN.Value = app.dispAWGN.Value;
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % swKaiser value changed function
        function swKaiserValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swKaiser.Value, 'On')
                app.KnobKaiser.Enable = 'on';
                app.dispKaiser.Enable = 'on';
                app.LampKaiser.Color = [0, 1, 0] ;
            else
                app.KnobKaiser.Enable = 'off';
                app.dispKaiser.Enable = 'off';
                app.LampKaiser.Color = [0.5, 0.5, 0.5] ;  
            end
        end

        % KnobKaiser value changed function
        function KnobKaiserValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % dispKaiser value changed function
        function dispKaiserValueChanged(app)
            app.KnobKaiser.Value = app.dispKaiser.Value;
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % KnobAWGN value changing function
        function KnobAWGNValueChanging(app, event)
            format('Bank')
            app.dispAWGN.Value = event.Value;
        end

        % KnobKaiser value changing function
        function KnobKaiserValueChanging(app, event)
            format('Bank')
            app.dispKaiser.Value = event.Value;
        end

        % swChirp value changed function
        function swChirpValueChanged(app, event)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % InputTarget1 value changed function
        function InputTarget1ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end    
        end

        % InputTarget2 value changed function
        function InputTarget2ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end           
        end

        % InputTarget3 value changed function
        function InputTarget3ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % InputTarget4 value changed function
        function InputTarget4ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end           
        end

        % InputTarget5 value changed function
        function InputTarget5ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end          
        end

        % InputTarget6 value changed function
        function InputTarget6ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end           
        end

        % swTarget1 value changed function
        function swTarget1ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swTarget1.Value, 'On')
                app.InputTarget1.Enable = 'on';
                app.LampTarget1.Color = [0, 1, 0];
            else
                app.InputTarget1.Enable = 'off';
                app.LampTarget1.Color = [0.5, 0.5, 0.5];
            end
        end

        % swTarget2 value changed function
        function swTarget2ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swTarget2.Value, 'On')
                app.InputTarget2.Enable = 'on';
                app.LampTarget2.Color = [0, 1, 0];
            else
                app.InputTarget2.Enable = 'off';
                app.LampTarget2.Color = [0.5, 0.5, 0.5];
            end
        end

        % swTarget3 value changed function
        function swTarget3ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swTarget3.Value, 'On')
                app.InputTarget3.Enable = 'on';
                app.LampTarget3.Color = [0, 1, 0];
            else
                app.InputTarget3.Enable = 'off';
                app.LampTarget3.Color = [0.5, 0.5, 0.5];
            end
        end

        % swTarget4 value changed function
        function swTarget4ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swTarget4.Value, 'On')
                app.InputTarget4.Enable = 'on';
                app.LampTarget4.Color = [0, 1, 0];
            else
                app.InputTarget4.Enable = 'off';
                app.LampTarget4.Color = [0.5, 0.5, 0.5];
            end
        end

        % swTarget5 value changed function
        function swTarget5ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swTarget5.Value, 'On')
                app.InputTarget5.Enable = 'on';
                app.LampTarget5.Color = [0, 1, 0];
            else
                app.InputTarget5.Enable = 'off';
                app.LampTarget5.Color = [0.5, 0.5, 0.5];
            end
        end

        % swTarget6 value changed function
        function swTarget6ValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
            if strcmp(app.swTarget6.Value, 'On')
                app.InputTarget6.Enable = 'on';
                app.LampTarget6.Color = [0, 1, 0];
            else
                app.InputTarget6.Enable = 'off';
                app.LampTarget6.Color = [0.5, 0.5, 0.5];
            end
        end

        % h value changed function
        function hValueChanged(app)
            format('Bank');
            app.disp_NearRange.Value = app.h.Value*tand(app.Theta.Value);
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % disp_h value changed function
        function disp_hValueChanged(app)
            format('Bank');
            app.h.Value = app.disp_h.Value;
            app.disp_NearRange.Value = app.h.Value*tand(app.Theta.Value);
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % h value changing function
        function hValueChanging(app, event)
            format('Bank');
            app.disp_h.Value = event.Value;
        end

        % Theta value changed function
        function ThetaValueChanged(app)
            app.disp_NearRange.Value = app.h.Value*tand(app.Theta.Value);
        end

        % KnobTheta value changed function
        function KnobThetaValueChanged(app)
            if app.autoUpdate
                updatePlot(app)
            end
        end

        % KnobTheta value changing function
        function KnobThetaValueChanging(app, event)
            format('Bank')
            app.Theta.Value = event.Value;
            app.disp_NearRange.Value = app.h.Value*tand(app.Theta.Value);
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 60 1151 617];
            app.UIFigure.Name = 'Single and Multi-target Chirp Signal Compression Simulator';
            setAutoResize(app, app.UIFigure, true)

            % Create UISot
            app.UISot = uiaxes(app.UIFigure);
            title(app.UISot, 'Signal Compressed after Matched Filter');
            xlabel(app.UISot, 'Time [us]');
            ylabel(app.UISot, 'Amplitude');
            app.UISot.Box = 'on';
            app.UISot.XGrid = 'on';
            app.UISot.YGrid = 'on';
            app.UISot.Position = [740 310 400 273];

            % Create UISotdb
            app.UISotdb = uiaxes(app.UIFigure);
            title(app.UISotdb, 'Signal Compressed after Matched Filter in [dB]');
            xlabel(app.UISotdb, 'Time [us]');
            ylabel(app.UISotdb, 'dB');
            app.UISotdb.Box = 'on';
            app.UISotdb.XGrid = 'on';
            app.UISotdb.YGrid = 'on';
            app.UISotdb.Position = [740 5 399 298];

            % Create UISt
            app.UISt = uiaxes(app.UIFigure);
            title(app.UISt, 'Radar Echo Signals');
            ylabel(app.UISt, 's(t)');
            app.UISt.Box = 'on';
            app.UISt.XGrid = 'on';
            app.UISt.YGrid = 'on';
            app.UISt.Position = [313 408 428 175];

            % Create UIStc
            app.UIStc = uiaxes(app.UIFigure);
            title(app.UIStc, 'Combination of all Signals');
            xlabel(app.UIStc, 'Time [us]');
            ylabel(app.UIStc, 's(t)');
            app.UIStc.Box = 'on';
            app.UIStc.XGrid = 'on';
            app.UIStc.YGrid = 'on';
            app.UIStc.Position = [313 5 428 200];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Units = 'pixels';
            app.TabGroup.Position = [8 8 298 567];

            % Create Tab
            app.Tab = uitab(app.TabGroup);
            app.Tab.Units = 'pixels';
            app.Tab.Title = 'Chirp Set';

            % Create LabelSlider
            app.LabelSlider = uilabel(app.Tab);
            app.LabelSlider.HorizontalAlignment = 'center';
            app.LabelSlider.VerticalAlignment = 'center';
            app.LabelSlider.FontSize = 16;
            app.LabelSlider.Position = [51 458 121 20];
            app.LabelSlider.Text = 'Bandwidth [MHz]';

            % Create Bw
            app.Bw = uislider(app.Tab);
            app.Bw.MajorTickLabels = {'1', '20', '40', '60', '80', '100'};
            app.Bw.ValueChangedFcn = createCallbackFcn(app, @BwValueChanged);
            app.Bw.ValueChangingFcn = createCallbackFcn(app, @BwValueChanging, true);
            app.Bw.Position = [28 449 232 3];
            app.Bw.Value = 30;

            % Create disp_B
            app.disp_B = uieditfield(app.Tab, 'numeric');
            app.disp_B.ValueChangedFcn = createCallbackFcn(app, @disp_BValueChanged);
            app.disp_B.Limits = [1 100];
            app.disp_B.RoundFractionalValues = 'on';
            app.disp_B.Position = [196 457 52 22];
            app.disp_B.Value = 30;

            % Create Label
            app.Label = uilabel(app.Tab);
            app.Label.HorizontalAlignment = 'center';
            app.Label.VerticalAlignment = 'center';
            app.Label.FontSize = 16;
            app.Label.Position = [51 382 133 20];
            app.Label.Text = 'Pulse duration [us]';

            % Create Tp
            app.Tp = uislider(app.Tab);
            app.Tp.Limits = [1 50];
            app.Tp.MajorTickLabels = {'1', '8', '15', '22', '29', '36', '43', '50'};
            app.Tp.ValueChangedFcn = createCallbackFcn(app, @TpValueChanged);
            app.Tp.ValueChangingFcn = createCallbackFcn(app, @TpValueChanging, true);
            app.Tp.Position = [28 372 232 3];
            app.Tp.Value = 5;

            % Create disp_T
            app.disp_T = uieditfield(app.Tab, 'numeric');
            app.disp_T.ValueChangedFcn = createCallbackFcn(app, @disp_TValueChanged);
            app.disp_T.Limits = [1 50];
            app.disp_T.RoundFractionalValues = 'on';
            app.disp_T.Position = [196 381 52 22];
            app.disp_T.Value = 5;

            % Create Label2
            app.Label2 = uilabel(app.Tab);
            app.Label2.HorizontalAlignment = 'center';
            app.Label2.VerticalAlignment = 'center';
            app.Label2.FontSize = 14;
            app.Label2.Position = [43 311 129 18];
            app.Label2.Text = 'Sampling Frequency';

            % Create Fsp
            app.Fsp = uislider(app.Tab);
            app.Fsp.Limits = [2 12];
            app.Fsp.ValueChangedFcn = createCallbackFcn(app, @FspValueChanged);
            app.Fsp.ValueChangingFcn = createCallbackFcn(app, @FspValueChanging, true);
            app.Fsp.Position = [28 301 232 3];
            app.Fsp.Value = 5;

            % Create disp_Fs
            app.disp_Fs = uieditfield(app.Tab, 'numeric');
            app.disp_Fs.ValueChangedFcn = createCallbackFcn(app, @disp_FsValueChanged);
            app.disp_Fs.Limits = [2 12];
            app.disp_Fs.RoundFractionalValues = 'on';
            app.disp_Fs.Position = [180 309 26 22];
            app.disp_Fs.Value = 5;

            % Create Label3
            app.Label3 = uilabel(app.Tab);
            app.Label3.FontSize = 14;
            app.Label3.Position = [211 311 21 18];
            app.Label3.Text = 'x B';

            % Create Button_Run
            app.Button_Run = uibutton(app.Tab, 'push');
            app.Button_Run.ButtonPushedFcn = createCallbackFcn(app, @Button_RunButtonPushed);
            app.Button_Run.FontSize = 28;
            app.Button_Run.FontColor = [0.0627 0 0.7686];
            app.Button_Run.Position = [127 6 151 42];
            app.Button_Run.Text = 'Run';

            % Create LabelKnob
            app.LabelKnob = uilabel(app.Tab);
            app.LabelKnob.HorizontalAlignment = 'center';
            app.LabelKnob.Position = [46.5 67 40 15];
            app.LabelKnob.Text = 'SNR = ';

            % Create KnobAWGN
            app.KnobAWGN = uiknob(app.Tab, 'continuous');
            app.KnobAWGN.Limits = [0.01 2];
            app.KnobAWGN.MajorTicks = [0.01 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2];
            app.KnobAWGN.ValueChangedFcn = createCallbackFcn(app, @KnobAWGNValueChanged);
            app.KnobAWGN.ValueChangingFcn = createCallbackFcn(app, @KnobAWGNValueChanging, true);
            app.KnobAWGN.Position = [52 108 60 60];
            app.KnobAWGN.Value = 1;

            % Create dispAWGN
            app.dispAWGN = uieditfield(app.Tab, 'numeric');
            app.dispAWGN.ValueChangedFcn = createCallbackFcn(app, @dispAWGNValueChanged);
            app.dispAWGN.Limits = [0.01 2];
            app.dispAWGN.Position = [86 63 26 22];
            app.dispAWGN.Value = 1;

            % Create Label4
            app.Label4 = uilabel(app.Tab);
            app.Label4.HorizontalAlignment = 'center';
            app.Label4.Position = [177 227 81 15];
            app.Label4.Text = 'Kaiser Window';

            % Create swKaiser
            app.swKaiser = uiswitch(app.Tab, 'slider');
            app.swKaiser.ValueChangedFcn = createCallbackFcn(app, @swKaiserValueChanged);
            app.swKaiser.Position = [187 204 45 20];

            % Create Label5
            app.Label5 = uilabel(app.Tab);
            app.Label5.HorizontalAlignment = 'center';
            app.Label5.Position = [180 67 45 15];
            app.Label5.Text = 'Alpha = ';

            % Create KnobKaiser
            app.KnobKaiser = uiknob(app.Tab, 'continuous');
            app.KnobKaiser.Limits = [0 3];
            app.KnobKaiser.ValueChangedFcn = createCallbackFcn(app, @KnobKaiserValueChanged);
            app.KnobKaiser.ValueChangingFcn = createCallbackFcn(app, @KnobKaiserValueChanging, true);
            app.KnobKaiser.Position = [188 108 60 60];
            app.KnobKaiser.Value = 0.8;

            % Create dispKaiser
            app.dispKaiser = uieditfield(app.Tab, 'numeric');
            app.dispKaiser.ValueChangedFcn = createCallbackFcn(app, @dispKaiserValueChanged);
            app.dispKaiser.Limits = [0 3];
            app.dispKaiser.Position = [227 63 26 22];
            app.dispKaiser.Value = 0.8;

            % Create LampAutoUpdate
            app.LampAutoUpdate = uilamp(app.Tab);
            app.LampAutoUpdate.Position = [92 10 15 15];

            % Create swAutoUpdate
            app.swAutoUpdate = uiswitch(app.Tab, 'slider');
            app.swAutoUpdate.ValueChangedFcn = createCallbackFcn(app, @swAutoUpdateValueChanged);
            app.swAutoUpdate.Position = [33 10 33.75 15];
            app.swAutoUpdate.Value = 'On';

            % Create Label6
            app.Label6 = uilabel(app.Tab);
            app.Label6.Position = [25 30 51 15];
            app.Label6.Text = 'Auto Run';

            % Create LabelSwitch
            app.LabelSwitch = uilabel(app.Tab);
            app.LabelSwitch.HorizontalAlignment = 'center';
            app.LabelSwitch.Position = [25 227 119 15];
            app.LabelSwitch.Text = 'White Gaussian Noise';

            % Create swAWGN
            app.swAWGN = uiswitch(app.Tab, 'slider');
            app.swAWGN.ValueChangedFcn = createCallbackFcn(app, @swAWGNValueChanged);
            app.swAWGN.Position = [48 204 45 20];

            % Create LampAWGN
            app.LampAWGN = uilamp(app.Tab);
            app.LampAWGN.Position = [119 207 15 15];

            % Create LampKaiser
            app.LampKaiser = uilamp(app.Tab);
            app.LampKaiser.Position = [257 207 15 15];

            % Create LabelDropDown
            app.LabelDropDown = uilabel(app.Tab);
            app.LabelDropDown.HorizontalAlignment = 'right';
            app.LabelDropDown.FontSize = 14;
            app.LabelDropDown.Position = [37 501 64 18];
            app.LabelDropDown.Text = 'Waveform';

            % Create swChirp
            app.swChirp = uidropdown(app.Tab);
            app.swChirp.Items = {'Up-Chirp', 'Down-Chirp', 'Symmetry-Chirp'};
            app.swChirp.ValueChangedFcn = createCallbackFcn(app, @swChirpValueChanged, true);
            app.swChirp.Position = [105 500 143 20];
            app.swChirp.Value = 'Symmetry-Chirp';

            % Create Tab2
            app.Tab2 = uitab(app.TabGroup);
            app.Tab2.Units = 'pixels';
            app.Tab2.Title = 'Targets';

            % Create Label28
            app.Label28 = uilabel(app.Tab2);
            app.Label28.FontSize = 14;
            app.Label28.FontWeight = 'bold';
            app.Label28.Position = [11 263 277 18];
            app.Label28.Text = 'Distance Between Target and Near Range';

            % Create swTarget1
            app.swTarget1 = uiswitch(app.Tab2, 'toggle');
            app.swTarget1.Orientation = 'horizontal';
            app.swTarget1.ValueChangedFcn = createCallbackFcn(app, @swTarget1ValueChanged);
            app.swTarget1.Position = [209 226 45 20];

            % Create LampTarget1
            app.LampTarget1 = uilamp(app.Tab2);
            app.LampTarget1.Position = [279 229 15 15];

            % Create swTarget2
            app.swTarget2 = uiswitch(app.Tab2, 'toggle');
            app.swTarget2.Orientation = 'horizontal';
            app.swTarget2.ValueChangedFcn = createCallbackFcn(app, @swTarget2ValueChanged);
            app.swTarget2.Position = [209 193 45 20];

            % Create LampTarget2
            app.LampTarget2 = uilamp(app.Tab2);
            app.LampTarget2.Position = [279 196 15 15];

            % Create swTarget3
            app.swTarget3 = uiswitch(app.Tab2, 'toggle');
            app.swTarget3.Orientation = 'horizontal';
            app.swTarget3.ValueChangedFcn = createCallbackFcn(app, @swTarget3ValueChanged);
            app.swTarget3.Position = [209 161 45 20];

            % Create LampTarget3
            app.LampTarget3 = uilamp(app.Tab2);
            app.LampTarget3.Position = [279 164 15 15];

            % Create swTarget4
            app.swTarget4 = uiswitch(app.Tab2, 'toggle');
            app.swTarget4.Orientation = 'horizontal';
            app.swTarget4.ValueChangedFcn = createCallbackFcn(app, @swTarget4ValueChanged);
            app.swTarget4.Position = [209 126 45 20];

            % Create LampTarget4
            app.LampTarget4 = uilamp(app.Tab2);
            app.LampTarget4.Position = [279 129 15 15];

            % Create swTarget5
            app.swTarget5 = uiswitch(app.Tab2, 'toggle');
            app.swTarget5.Orientation = 'horizontal';
            app.swTarget5.ValueChangedFcn = createCallbackFcn(app, @swTarget5ValueChanged);
            app.swTarget5.Position = [209 94 45 20];

            % Create LampTarget5
            app.LampTarget5 = uilamp(app.Tab2);
            app.LampTarget5.Position = [279 97 15 15];

            % Create swTarget6
            app.swTarget6 = uiswitch(app.Tab2, 'toggle');
            app.swTarget6.Orientation = 'horizontal';
            app.swTarget6.ValueChangedFcn = createCallbackFcn(app, @swTarget6ValueChanged);
            app.swTarget6.Position = [209 60 45 20];

            % Create LampTarget6
            app.LampTarget6 = uilamp(app.Tab2);
            app.LampTarget6.Position = [279 63 15 15];

            % Create Panel
            app.Panel = uipanel(app.Tab2);
            app.Panel.BorderType = 'line';
            app.Panel.BackgroundColor = [0.9373 0.9373 0.9373];
            app.Panel.FontName = 'Helvetica';
            app.Panel.FontUnits = 'pixels';
            app.Panel.FontSize = 12;
            app.Panel.Units = 'pixels';
            app.Panel.Position = [2 58 179 198];

            % Create Label27
            app.Label27 = uilabel(app.Panel);
            app.Label27.FontSize = 16;
            app.Label27.Position = [157 2 20 20];
            app.Label27.Text = 'm';

            % Create Label26
            app.Label26 = uilabel(app.Panel);
            app.Label26.FontSize = 16;
            app.Label26.Position = [157 35 20 20];
            app.Label26.Text = 'm';

            % Create Label25
            app.Label25 = uilabel(app.Panel);
            app.Label25.FontSize = 16;
            app.Label25.Position = [157 68 20 20];
            app.Label25.Text = 'm';

            % Create Label24
            app.Label24 = uilabel(app.Panel);
            app.Label24.FontSize = 16;
            app.Label24.Position = [157 102 20 20];
            app.Label24.Text = 'm';

            % Create Label23
            app.Label23 = uilabel(app.Panel);
            app.Label23.FontSize = 16;
            app.Label23.Position = [157 135 20 20];
            app.Label23.Text = 'm';

            % Create Label22
            app.Label22 = uilabel(app.Panel);
            app.Label22.FontSize = 16;
            app.Label22.Position = [157 168 20 20];
            app.Label22.Text = 'm';

            % Create Label21
            app.Label21 = uilabel(app.Panel);
            app.Label21.HorizontalAlignment = 'right';
            app.Label21.FontSize = 16;
            app.Label21.Position = [8 1 60 20];
            app.Label21.Text = 'Target 6';

            % Create InputTarget6
            app.InputTarget6 = uieditfield(app.Panel, 'numeric');
            app.InputTarget6.ValueChangedFcn = createCallbackFcn(app, @InputTarget6ValueChanged);
            app.InputTarget6.Limits = [1 Inf];
            app.InputTarget6.FontSize = 16;
            app.InputTarget6.Position = [75 2 75 20];
            app.InputTarget6.Value = 1;

            % Create Label20
            app.Label20 = uilabel(app.Panel);
            app.Label20.HorizontalAlignment = 'right';
            app.Label20.FontSize = 16;
            app.Label20.Position = [8 35 60 20];
            app.Label20.Text = 'Target 5';

            % Create InputTarget5
            app.InputTarget5 = uieditfield(app.Panel, 'numeric');
            app.InputTarget5.ValueChangedFcn = createCallbackFcn(app, @InputTarget5ValueChanged);
            app.InputTarget5.Limits = [1 Inf];
            app.InputTarget5.FontSize = 16;
            app.InputTarget5.Position = [75 35 75 20];
            app.InputTarget5.Value = 1;

            % Create Label19
            app.Label19 = uilabel(app.Panel);
            app.Label19.HorizontalAlignment = 'right';
            app.Label19.FontSize = 16;
            app.Label19.Position = [8 67 60 20];
            app.Label19.Text = 'Target 4';

            % Create InputTarget4
            app.InputTarget4 = uieditfield(app.Panel, 'numeric');
            app.InputTarget4.ValueChangedFcn = createCallbackFcn(app, @InputTarget4ValueChanged);
            app.InputTarget4.Limits = [1 Inf];
            app.InputTarget4.FontSize = 16;
            app.InputTarget4.Position = [75 68 75 20];
            app.InputTarget4.Value = 1;

            % Create Label18
            app.Label18 = uilabel(app.Panel);
            app.Label18.HorizontalAlignment = 'right';
            app.Label18.FontSize = 16;
            app.Label18.Position = [8 102 60 20];
            app.Label18.Text = 'Target 3';

            % Create InputTarget3
            app.InputTarget3 = uieditfield(app.Panel, 'numeric');
            app.InputTarget3.ValueChangedFcn = createCallbackFcn(app, @InputTarget3ValueChanged);
            app.InputTarget3.Limits = [1 Inf];
            app.InputTarget3.FontSize = 16;
            app.InputTarget3.Position = [75 102 75 20];
            app.InputTarget3.Value = 1;

            % Create LabelNumericEditField2
            app.LabelNumericEditField2 = uilabel(app.Panel);
            app.LabelNumericEditField2.HorizontalAlignment = 'right';
            app.LabelNumericEditField2.FontSize = 16;
            app.LabelNumericEditField2.Position = [8 134 60 20];
            app.LabelNumericEditField2.Text = 'Target 2';

            % Create InputTarget2
            app.InputTarget2 = uieditfield(app.Panel, 'numeric');
            app.InputTarget2.ValueChangedFcn = createCallbackFcn(app, @InputTarget2ValueChanged);
            app.InputTarget2.Limits = [1 Inf];
            app.InputTarget2.FontSize = 16;
            app.InputTarget2.Position = [75 135 75 20];
            app.InputTarget2.Value = 1;

            % Create LabelNumericEditField
            app.LabelNumericEditField = uilabel(app.Panel);
            app.LabelNumericEditField.HorizontalAlignment = 'right';
            app.LabelNumericEditField.FontSize = 16;
            app.LabelNumericEditField.Position = [8 168 60 20];
            app.LabelNumericEditField.Text = 'Target 1';

            % Create InputTarget1
            app.InputTarget1 = uieditfield(app.Panel, 'numeric');
            app.InputTarget1.ValueChangedFcn = createCallbackFcn(app, @InputTarget1ValueChanged);
            app.InputTarget1.Limits = [1 Inf];
            app.InputTarget1.FontSize = 16;
            app.InputTarget1.Position = [75 168 75 20];
            app.InputTarget1.Value = 1;

            % Create LabelSlider2
            app.LabelSlider2 = uilabel(app.Tab2);
            app.LabelSlider2.HorizontalAlignment = 'right';
            app.LabelSlider2.FontSize = 16;
            app.LabelSlider2.Position = [2 325 39 20];
            app.LabelSlider2.Text = 'Nadir';

            % Create gr0
            app.gr0 = uislider(app.Tab2);
            app.gr0.Limits = [50 1000];
            app.gr0.MajorTickLabels = {'50', '', '', '', '', '1000'};
            app.gr0.Enable = 'off';
            app.gr0.Position = [52 334 101 3];
            app.gr0.Value = 600;

            % Create Label29
            app.Label29 = uilabel(app.Tab2);
            app.Label29.HorizontalAlignment = 'right';
            app.Label29.FontSize = 14;
            app.Label29.Position = [110 344 76 18];
            app.Label29.Text = 'Near Range';

            % Create Label30
            app.Label30 = uilabel(app.Tab2);
            app.Label30.HorizontalAlignment = 'right';
            app.Label30.FontSize = 14;
            app.Label30.Position = [208 307 48 18];
            app.Label30.Text = 'Targets';

            % Create Lamp
            app.Lamp = uilamp(app.Tab2);
            app.Lamp.Position = [185 328 15 15];
            app.Lamp.Color = [0 0.4824 1];

            % Create Lamp2
            app.Lamp2 = uilamp(app.Tab2);
            app.Lamp2.Position = [199 328 15 15];
            app.Lamp2.Color = [1 0.498 0];

            % Create Lamp3
            app.Lamp3 = uilamp(app.Tab2);
            app.Lamp3.Position = [224 328 15 15];
            app.Lamp3.Color = [1 0.749 0];

            % Create Lamp4
            app.Lamp4 = uilamp(app.Tab2);
            app.Lamp4.Position = [250 328 15 15];
            app.Lamp4.Color = [0.6235 0 0.7725];

            % Create Lamp5
            app.Lamp5 = uilamp(app.Tab2);
            app.Lamp5.Position = [264 328 15 15];
            app.Lamp5.Color = [1 1 0];

            % Create Lamp6
            app.Lamp6 = uilamp(app.Tab2);
            app.Lamp6.Position = [278 328 15 15];
            app.Lamp6.Color = [0 0.8824 1];

            % Create LabelSlider3
            app.LabelSlider3 = uilabel(app.Tab2);
            app.LabelSlider3.HorizontalAlignment = 'right';
            app.LabelSlider3.FontSize = 14;
            app.LabelSlider3.Position = [2 480 41 18];
            app.LabelSlider3.Text = 'Height';

            % Create h
            app.h = uislider(app.Tab2);
            app.h.Limits = [400 1000];
            app.h.Orientation = 'vertical';
            app.h.ValueChangedFcn = createCallbackFcn(app, @hValueChanged);
            app.h.ValueChangingFcn = createCallbackFcn(app, @hValueChanging, true);
            app.h.Position = [39 352 3 111];
            app.h.Value = 700;

            % Create LabelNumericEditField3
            app.LabelNumericEditField3 = uilabel(app.Tab2);
            app.LabelNumericEditField3.HorizontalAlignment = 'right';
            app.LabelNumericEditField3.FontSize = 14;
            app.LabelNumericEditField3.Position = [152 480 116 18];
            app.LabelNumericEditField3.Text = 'Angle of Incidence';

            % Create Theta
            app.Theta = uieditfield(app.Tab2, 'numeric');
            app.Theta.ValueChangedFcn = createCallbackFcn(app, @ThetaValueChanged);
            app.Theta.Limits = [10 45];
            app.Theta.Position = [227 451 43 22];
            app.Theta.Value = 40;

            % Create Label31
            app.Label31 = uilabel(app.Tab2);
            app.Label31.FontSize = 14;
            app.Label31.Position = [227 431 52 18];
            app.Label31.Text = '[degree]';

            % Create disp_NearRange
            app.disp_NearRange = uieditfield(app.Tab2, 'numeric');
            app.disp_NearRange.Limits = [70 1000];
            app.disp_NearRange.Editable = 'off';
            app.disp_NearRange.BackgroundColor = [0.9373 0.9373 0.9373];
            app.disp_NearRange.Position = [118 361 43 22];
            app.disp_NearRange.Value = 600;

            % Create Label32
            app.Label32 = uilabel(app.Tab2);
            app.Label32.FontSize = 14;
            app.Label32.Position = [165 363 27 18];
            app.Label32.Text = '[km]';

            % Create Label33
            app.Label33 = uilabel(app.Tab2);
            app.Label33.FontSize = 14;
            app.Label33.Position = [92 480 27 18];
            app.Label33.Text = '[km]';

            % Create disp_h
            app.disp_h = uieditfield(app.Tab2, 'numeric');
            app.disp_h.ValueChangedFcn = createCallbackFcn(app, @disp_hValueChanged);
            app.disp_h.Limits = [400 1000];
            app.disp_h.Position = [46 478 43 22];
            app.disp_h.Value = 700;

            % Create KnobTheta
            app.KnobTheta = uiknob(app.Tab2, 'continuous');
            app.KnobTheta.Limits = [10 45];
            app.KnobTheta.ValueChangedFcn = createCallbackFcn(app, @KnobThetaValueChanged);
            app.KnobTheta.ValueChangingFcn = createCallbackFcn(app, @KnobThetaValueChanging, true);
            app.KnobTheta.FontSize = 11;
            app.KnobTheta.Position = [158 412 42 42];
            app.KnobTheta.Value = 40;

            % Create Label34
            app.Label34 = uilabel(app.Tab2);
            app.Label34.FontSize = 14;
            app.Label34.FontWeight = 'bold';
            app.Label34.Position = [15 512 270 18];
            app.Label34.Text = 'Set Radar Height and Angle of Incidence';

            % Create Tab3
            app.Tab3 = uitab(app.TabGroup);
            app.Tab3.Units = 'pixels';
            app.Tab3.Title = 'Hints';
            app.Tab3.ForegroundColor = [1 0.4824 0];

            % Create Label13
            app.Label13 = uilabel(app.Tab3);
            app.Label13.HorizontalAlignment = 'center';
            app.Label13.VerticalAlignment = 'center';
            app.Label13.FontSize = 16;
            app.Label13.FontWeight = 'bold';
            app.Label13.Position = [1 55 192 20];
            app.Label13.Text = 'Author: Guan-Wen CHEN';

            % Create Label14
            app.Label14 = uilabel(app.Tab3);
            app.Label14.HorizontalAlignment = 'center';
            app.Label14.VerticalAlignment = 'center';
            app.Label14.FontSize = 14;
            app.Label14.FontAngle = 'italic';
            app.Label14.Position = [1 20 219 18];
            app.Label14.Text = 'National Central University, Taiwan';

            % Create Label15
            app.Label15 = uilabel(app.Tab3);
            app.Label15.HorizontalAlignment = 'center';
            app.Label15.VerticalAlignment = 'center';
            app.Label15.FontSize = 14;
            app.Label15.FontAngle = 'italic';
            app.Label15.Position = [1 37 225 18];
            app.Label15.Text = 'Graduate Institute of Space Science';

            % Create Label16
            app.Label16 = uilabel(app.Tab3);
            app.Label16.HorizontalAlignment = 'center';
            app.Label16.VerticalAlignment = 'center';
            app.Label16.FontSize = 14;
            app.Label16.Position = [1 3 146 18];
            app.Label16.Text = 'gw.chen19@gmail.com';

            % Create Label17
            app.Label17 = uilabel(app.Tab3);
            app.Label17.VerticalAlignment = 'center';
            app.Label17.FontSize = 16;
            app.Label17.Position = [0.5 74 86 20];
            app.Label17.Text = 'Version: 1.0';

            % Create Label35
            app.Label35 = uilabel(app.Tab3);
            app.Label35.FontSize = 16;
            app.Label35.FontWeight = 'bold';
            app.Label35.Position = [3 515 290 20];
            app.Label35.Text = 'Step 1:  Chirp Set - Set the basic chirp';

            % Create TextArea
            app.TextArea = uitextarea(app.Tab3);
            app.TextArea.Editable = 'off';
            app.TextArea.FontName = 'Times New Roman';
            app.TextArea.FontSize = 16;
            app.TextArea.BackgroundColor = [0.9373 0.9373 0.9373];
            app.TextArea.Position = [-1 303 298 205];
            app.TextArea.Value = {'First, determine the bandwidth, pulse '; 'duration, sampling frequency and chirp '; 'form you want, it will display a simple '; 'basic chirp waveform s(t), and results '; 'of the output signal compressed after a '; 'default, corresponding matched filter h(t). '; 'It will auto update the results when you '; 'change the initial setting.'; ''; 'Then it will use this basic waveform for bothof the single and multiple targets simulation.'};

            % Create Label36
            app.Label36 = uilabel(app.Tab3);
            app.Label36.FontSize = 16;
            app.Label36.FontWeight = 'bold';
            app.Label36.Position = [2 274 294 20];
            app.Label36.Text = 'Step 2:  Targets - Set ranges of targets';

            % Create TextArea2
            app.TextArea2 = uitextarea(app.Tab3);
            app.TextArea2.Editable = 'off';
            app.TextArea2.FontName = 'Times New Roman';
            app.TextArea2.FontSize = 16;
            app.TextArea2.BackgroundColor = [0.9373 0.9373 0.9373];
            app.TextArea2.Position = [-1 111 298 156];
            app.TextArea2.Value = {'Second, assume that all of the targets are '; 'on the ground range, input the radar height, '; 'incidence angle, and ranges between the '; 'near range and targets, check if it can be '; 'identified by compressing the combination '; 'signal after the matched filter.'};

            % Create Label7
            app.Label7 = uilabel(app.UIFigure);
            app.Label7.HorizontalAlignment = 'center';
            app.Label7.VerticalAlignment = 'center';
            app.Label7.FontName = 'Times New Roman';
            app.Label7.FontSize = 20;
            app.Label7.FontWeight = 'bold';
            app.Label7.FontColor = [1 0.498 0];
            app.Label7.Position = [119.5 591 913 26];
            app.Label7.Text = '¡X¡X¡X¡X¡X¡X  Easily Simulate Chirp Signal Compression of Single/Multiple Range Targets  ¡X¡X¡X¡X¡X¡X';

            % Create UIStm
            app.UIStm = uiaxes(app.UIFigure);
            title(app.UIStm, 'Multiple-Targets  Radar Echo Signals');
            ylabel(app.UIStm, 's(t)');
            app.UIStm.Box = 'on';
            app.UIStm.XGrid = 'on';
            app.UIStm.YGrid = 'on';
            app.UIStm.Position = [313 206 428 204];
        end
    end

    methods (Access = public)

        % Construct app
        function app = CSCS()

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end