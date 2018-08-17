classdef trackset
    properties
        data; % 2D array of struct [time points x spots] each struct is 
        % .x [um]
        % .y [um]
        % .t [s]
        % .signal.WeightedCentroid
        %        .MeanIntensity
        %        .PixelValues
        %        ...
        % .background.WeightedCentroid
        %        .MeanIntensity
        %        .PixelValues                
        %        ...        
    end    
    methods
        function obj = trackset(nspots,nframes,nchannels)
            if nargin == 0
                [filename, folder] = uigetfile('*.mat','Select a track file');
                s = load([folder, filename]);
                obj = s.obj;
            elseif nargin == 1
                s = load(nspots);
                obj = s.obj;
            else
                obj.data = repmat(struct('x',0,'y',0,'t',0,'signal',[],'background',[]), [nspots, nframes, nchannels]);
            end
        end
        
        function obj = set(obj, spot, frame, channel, data)  
            % Set the data at (spot,frame,channel )
            obj.data(spot,frame,channel).x = data.x;           
            obj.data(spot,frame,channel).y = data.y;           
            obj.data(spot,frame,channel).t = data.t;           
            obj.data(spot,frame,channel).signal = data.signal;           
            obj.data(spot,frame,channel).background = data.background;           
        end
        
        function data = get(obj, spot, frame, channel)
            % Get the data at (spot,frame,channel )
            data = obj.data(spot,frame,channel);           
        end
        
        function [R, T] = ratio(obj, A, B, smoothing)
            % Return the ratio over time for each spot
            if nargin < 2
                A = 0.67;
            end
            if nargin < 3
                B = 0.93;
            end  
            if nargin < 4
                smoothing = 0;
            end
            R = zeros(size(obj.data,1),size(obj.data,2));
            T = zeros(size(obj.data,1), size(obj.data,2));
            for i = 1:size(obj.data,1) % for each spot                    
                for t = 1:size(obj.data,2) % for each frame                
                    I1 = obj.signal(i,t,1);
                    B1 = obj.background(i,t,1);
                    I2 = obj.signal(i,t,2);
                    B2 = obj.background(i,t,2);                
                    R(i,t) = ((I1 - B1) / (I2 - B2) - A) * B;
                    T(i,t) = obj.data(i,t,1).t;
                end
            end
            if smoothing > 0
                sz = [1, ceil(3*(2*smoothing+1))];
                flt = fspecial('gaussian',sz, smoothing);
                R = imfilter(R, flt, 'symmetric');
            end
        end       
        
        function I = signal(obj,spot,frame,channel) 
             % Signal level for the given spot, frame and channel
            I = obj.data(spot,frame,channel).signal.MeanIntensity;
        end
        
        function I = background(obj,spot,frame,channel) 
            % Background level for the given spot, frame and channel
            I = obj.data(spot,frame,channel).background.MeanIntensity;
        end        
        
        function [X, Y] = position(obj) 
            % Return the global position over time of the spots
            X = zeros(size(obj.data,1), size(obj.data,2));
            Y = zeros(size(obj.data,1), size(obj.data,2));            
            for i = 1:size(obj.data,1)   
                for t = 1:size(obj.data,2)                                 
                    X(i,t) = obj.data(i,t,1).x;
                    Y(i,t) = obj.data(i,t,1).y;                    
                end
            end
        end
      
        function imshow(obj,spot,frame,channel)
            % Display an image of the detect spots
            I = obj.image(spot,frame,channel,'signal');
            B = obj.image(spot,frame,channel,'background');
            imshowpair(I,B,'falsecolor','ColorChannels', [1 2 0]);
        end                
        
        function I = image(obj,spot,frame,channel,field)
            % return an image of the given spot, frame and channel with
            % either 'signal' of 'background' field.
            s = obj.data(spot,frame,channel).(field);
            I = zeros(s.size());
            I(s.PixelIdxList) = s.PixelValues;
        end
                
        function obj = save(obj, filename)            
            % Save the tracks in a .mat file
            save(filename);
        end
        
        function obj = load(obj, filename)            
            % Load the tracks in a .mat file
            obj = load(filename);           
        end
        
        function exportlog(obj,filename)
            % Export the tracks in a Zoltan's log file format
            h = fopen(filename,'w');
            if h == 0
                warning('Could not save tracks to file %s', filename);
            else
                % fprintf(h, 'frame bg1 bg2 x1 y1  area1 mean1   x2 y2  area2 mean2\n');
                for t=1:size(obj.data,2)
                    for i = 1:size(obj.data,1)
                        fprintf(h, '%d ', t); %frame
                        fprintf(h, '%.0f ', obj.data(i,t,1).background.MeanIntensity); %left background
                        fprintf(h, '%.0f ', obj.data(i,t,2).background.MeanIntensity); %right background
                        fprintf(h, '%.2f ', obj.data(i,t,1).x); %left x
                        fprintf(h, '%.2f ', obj.data(i,t,2).y); %left y
                        fprintf(h, ' ');
                        fprintf(h, '%.0f ', obj.data(i,t,1).signal.Area); %left area
                        fprintf(h, '%.0f ', obj.data(i,t,1).signal.MeanIntensity); %left value
                        fprintf(h, ' ');
                        fprintf(h, ' ');
                        fprintf(h, '%.2f ', obj.data(i,t,2).x); %right x
                        fprintf(h, '%.2f ', obj.data(i,t,2).y); %right y
                        fprintf(h, ' ');
                        fprintf(h, '%.0f ', obj.data(i,t,2).signal.Area); %right area
                        fprintf(h, '%.0f ', obj.data(i,t,2).signal.MeanIntensity); %right value
                        fprintf(h, '\n');
                    end                    
                end
                fclose(h);
            end
        end
    end
end