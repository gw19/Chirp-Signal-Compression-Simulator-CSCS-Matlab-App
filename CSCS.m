classdef CSCS < matlab.apps.AppBase & ViewCSCS

    methods (Access = private)
        
        function updateCSCS(app)
            
            updatePlot(app)

        end
        
    end


    methods (Access = private)

        % Bw value changed function
        function BwValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
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
            updateCSCS(app)
        end

        % Tp value changed function
        function TpValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end
        end

        % disp_B value changed function
        function disp_BValueChanged(app)
            app.Bw.Value = app.disp_B.Value;
            if app.autoUpdate
                updateCSCS(app)
            end
        end

        % disp_T value changed function
        function disp_TValueChanged(app)
            app.Tp.Value = app.disp_T.Value;
            if app.autoUpdate
                updateCSCS(app)
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
                updateCSCS(app)
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
                updateCSCS(app)
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
                updateCSCS(app)
            end
        end

        % KnobAWGN value changed function
        function KnobAWGNValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end
        end

        % dispAWGN value changed function
        function dispAWGNValueChanged(app)
            app.KnobAWGN.Value = app.dispAWGN.Value;
            if app.autoUpdate
                updateCSCS(app)
            end
        end

        % swKaiser value changed function
        function swKaiserValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
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
                updateCSCS(app)
            end
        end

        % dispKaiser value changed function
        function dispKaiserValueChanged(app)
            app.KnobKaiser.Value = app.dispKaiser.Value;
            if app.autoUpdate
                updateCSCS(app)
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
                updateCSCS(app)
            end
        end

        % InputTarget1 value changed function
        function InputTarget1ValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end    
        end

        % InputTarget2 value changed function
        function InputTarget2ValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end           
        end

        % InputTarget3 value changed function
        function InputTarget3ValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end
        end

        % InputTarget4 value changed function
        function InputTarget4ValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end           
        end

        % InputTarget5 value changed function
        function InputTarget5ValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end          
        end

        % InputTarget6 value changed function
        function InputTarget6ValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
            end           
        end

        % swTarget1 value changed function
        function swTarget1ValueChanged(app)
            if app.autoUpdate
                updateCSCS(app)
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
                updateCSCS(app)
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
                updateCSCS(app)
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
                updateCSCS(app)
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
                updateCSCS(app)
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
                updateCSCS(app)
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
                updateCSCS(app)
            end
        end

        % disp_h value changed function
        function disp_hValueChanged(app)
            format('Bank');
            app.h.Value = app.disp_h.Value;
            app.disp_NearRange.Value = app.h.Value*tand(app.Theta.Value);
            if app.autoUpdate
                updateCSCS(app)
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
                updateCSCS(app)
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
    methods (Access = protected)

        % Create UIFigure and components
        function createComponents(app)
            createComponents@ViewCSCS(app)
            setAutoResize(app, app.UIFigure, true)

            app.Bw.ValueChangedFcn = createCallbackFcn(app, @BwValueChanged);
            app.Bw.ValueChangingFcn = createCallbackFcn(app, @BwValueChanging, true);
            app.disp_B.ValueChangedFcn = createCallbackFcn(app, @disp_BValueChanged);
            app.Tp.ValueChangedFcn = createCallbackFcn(app, @TpValueChanged);
            app.Tp.ValueChangingFcn = createCallbackFcn(app, @TpValueChanging, true);
            app.disp_T.ValueChangedFcn = createCallbackFcn(app, @disp_TValueChanged);
            app.Fsp.ValueChangedFcn = createCallbackFcn(app, @FspValueChanged);
            app.Fsp.ValueChangingFcn = createCallbackFcn(app, @FspValueChanging, true);
            app.disp_Fs.ValueChangedFcn = createCallbackFcn(app, @disp_FsValueChanged);
            app.Button_Run.ButtonPushedFcn = createCallbackFcn(app, @Button_RunButtonPushed);
            app.KnobAWGN.ValueChangedFcn = createCallbackFcn(app, @KnobAWGNValueChanged);
            app.KnobAWGN.ValueChangingFcn = createCallbackFcn(app, @KnobAWGNValueChanging, true);
            app.dispAWGN.ValueChangedFcn = createCallbackFcn(app, @dispAWGNValueChanged);
            app.swKaiser.ValueChangedFcn = createCallbackFcn(app, @swKaiserValueChanged);
            app.KnobKaiser.ValueChangedFcn = createCallbackFcn(app, @KnobKaiserValueChanged);
            app.KnobKaiser.ValueChangingFcn = createCallbackFcn(app, @KnobKaiserValueChanging, true);
            app.dispKaiser.ValueChangedFcn = createCallbackFcn(app, @dispKaiserValueChanged);
            app.swAutoUpdate.ValueChangedFcn = createCallbackFcn(app, @swAutoUpdateValueChanged);
            app.swAWGN.ValueChangedFcn = createCallbackFcn(app, @swAWGNValueChanged);
            app.swChirp.ValueChangedFcn = createCallbackFcn(app, @swChirpValueChanged, true);
            app.swTarget1.ValueChangedFcn = createCallbackFcn(app, @swTarget1ValueChanged);
            app.swTarget2.ValueChangedFcn = createCallbackFcn(app, @swTarget2ValueChanged);
            app.swTarget3.ValueChangedFcn = createCallbackFcn(app, @swTarget3ValueChanged);
            app.swTarget4.ValueChangedFcn = createCallbackFcn(app, @swTarget4ValueChanged);
            app.swTarget5.ValueChangedFcn = createCallbackFcn(app, @swTarget5ValueChanged);
            app.swTarget6.ValueChangedFcn = createCallbackFcn(app, @swTarget6ValueChanged);
            app.InputTarget6.ValueChangedFcn = createCallbackFcn(app, @InputTarget6ValueChanged);
            app.InputTarget5.ValueChangedFcn = createCallbackFcn(app, @InputTarget5ValueChanged);
            app.InputTarget4.ValueChangedFcn = createCallbackFcn(app, @InputTarget4ValueChanged);
            app.InputTarget3.ValueChangedFcn = createCallbackFcn(app, @InputTarget3ValueChanged);
            app.InputTarget2.ValueChangedFcn = createCallbackFcn(app, @InputTarget2ValueChanged);
            app.InputTarget1.ValueChangedFcn = createCallbackFcn(app, @InputTarget1ValueChanged);
            app.h.ValueChangedFcn = createCallbackFcn(app, @hValueChanged);
            app.h.ValueChangingFcn = createCallbackFcn(app, @hValueChanging, true);
            app.Theta.ValueChangedFcn = createCallbackFcn(app, @ThetaValueChanged);
            app.disp_h.ValueChangedFcn = createCallbackFcn(app, @disp_hValueChanged);
            app.KnobTheta.ValueChangedFcn = createCallbackFcn(app, @KnobThetaValueChanged);
            app.KnobTheta.ValueChangingFcn = createCallbackFcn(app, @KnobThetaValueChanging, true);
        end
    end

    methods (Access = public)

        % Construct app
        function app = CSCS()

            % Create and configure components
            createComponents(app)

            if nargout == 0
                run(app)
                clear app
            end
        end

        function run(app)
            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            % delete(app.UIFigure) 
        end
    end
end