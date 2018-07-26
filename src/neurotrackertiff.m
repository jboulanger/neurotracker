classdef neurotrackertiff
    % neurotrackertiff
    %   Read images in a Tiff file acquired with neurotracker from Xavier
    %   Baels.        
    %
    %  nt = neurotrackertiff('test.tif');
    %  nt.playmovie();
    %  nt.closetiff();
    %
    %  Tiff tags:
    %  Tag 65000    TAG_TRACKING_POSITION_X,  TIFF_SHORT,  "Tracking Result X component [pixel]" ,
    %  Tag 65001    TAG_TRACKING_POSITION_Y, TIFF_SHORT,  "Tracking Result Y component [pixel]",
    %  Tag 65002    TAG_TARGET_POSITION_X,  TIFF_SHORT,  "Target Position X component [pixel]",
    %  Tag 65003    TAG_TARGET_POSITION_Y, TIFF_SHORT,  "Target Position Y component [pixel]",
    %  Tag 65004    TAG_TRAY_POSITION_X,  TIFF_LONG,  "Stage Position X component [micrometer]",
    %  Tag 65005    TAG_TRAY_POSITION_Y, TIFF_LONG,  "Stage Position Y component [micrometer]",
    %  Tag 65006    TAG_PIXEL_WIDTH,  TIFF_FLOAT,  "Pixel width [micrometer]",
    %  Tag 65007    TAG_PIXEL_HEIGHT, TIFF_FLOAT,  "Pixel height [micrometer]",
    %  Tag 65008    TAG_MAGNIFICATION,  TIFF_SHORT,  "Magnification",
    %  Tag 65009    TAG_TIMESTAMP, TIFF_LONG,  "Timestamp starting from 0 [ms]",
    %  Tag 65010    TAG_EXPOSURE_TIME,  TIFF_LONG,  "Exposure time [ms]"

    properties
        pathtofile; % a single string or a cell of string
        datasize = [1 1 1 2 1]; % size in all dimensions [xyzct]
        warningmode; % prints warning or just errors
        metadata; % the list of all metadata as a array of structs
        imagecount; % number of images for each file
        tif; % TIFF object or array of tiff objects
        stageindex; % Nimage x (x,y) array            
        magnification; % Microscope magnification
        pixelsize; % [dx,dy] Size of the pixel in um 
        timestamps; % array of time stamps for each image
        exposuretime; % Exposure time for each image
        flip; % Flip images X and Y directions        
        stageunitum; % stage unit in um
    end

    methods      
        function obj = neurotrackertiff(url, warningmode)
            % Setup the tiff reader and load metadata
            %
            if nargin < 2
                warningmode = 0;
            end
            obj.pathtofile = url;
            obj = obj.setwarningmode(warningmode);
            obj.flip = ~[true,false;false,true];
            obj.pixelsize =[6.5 6.5];     
            obj.stageunitum = 1;
            obj = obj.loadmetadata();
            obj = obj.opentiff();
            obj = obj.parsestageindex();
            obj = obj.parsetimestamps();
        end
        
        function printinfo(obj)
            fprintf('*** neurotrackertiff info ***\n');
            if ~iscell(obj.pathtofile)
                fprintf('file path: ''%s''\n', obj.pathtofile);
            else
                fprintf('file paths:');
                for i = numel(pathtofile)
                    fprintf( '''%d: %s''\n', i, obj.pathtofile{i});
                end
            end
            fprintf('pixel size: %.2fx%.2fum\n', obj.pixelsize);
            fprintf('magnification: %.2fX\n', obj.magnification);
            fprintf('stage unit: %fum\n', obj.stageunitum);
            fov = max(obj.range(),[],2) - min(obj.range(),[],2);
            fprintf('total field of view: %.2f x %.2f um\n', fov);
            D = max(obj.timestamps) - min(obj.timestamps);
            fprintf('duration %d frames / %.2fsec\n', obj.datasize(5), D);
        end

        function obj = setwarningmode(obj, mode)
            % Set the warning mode
            obj.warningmode = mode;
            if obj.warningmode < 2
                warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
            else
                warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
            end
        end
        
        function obj = setpixelsize(obj, camera_pixelsize_um)
            % set the physical pixel size of the camera in um
            obj.pixelsize = [camera_pixelsize_um camera_pixelsize_um];
        end
        
        function obj = setmagnification(obj, M)
            % set the total magnification of the microscope
            obj.magnification = M;
        end

        function obj = opentiff(obj)
            % Open the tif files
            if iscell(obj.pathtofile)
                for n = 1:numel(obj.pathtofile)
                    obj.debug(sprintf('opening tif file %s', obj.pathtofile{n}));
                    obj.tif{n} = Tiff(obj.pathtofile{n}, 'r');
                end
            else
                obj.debug(sprintf('opening tif file %s', obj.pathtofile));
                obj.tif = Tiff(obj.pathtofile, 'r');
            end
            im = obj.imread(1,1);
            obj.datasize(1:2) = [size(im,2) size(im,1)];
        end

        function obj = closetiff(obj)
            % close the tiff file
            if iscell(obj.tif)
                for n = 1:numel(obj.tif)
                    close(obj.tif{n});
                end
            else
                 close(obj.tif);
            end
        end

        function obj = loadmetadata(obj)
            % load all the metadata from the tiff file
            if iscell(obj.pathtofile)
                obj.metadata = [];
                obj.imagecount = zeros(1,length(obj.pathtofile));
                for n = 1:length(obj.pathtofile)
                    obj.debug(sprintf('Loading metadata from %s\n', obj.pathtofile{n}));
                    s = imfinfo(obj.pathtofile{n});
                    obj.metadata = [obj.metadata; s];
                    obj.imagecount(n) = length(s);
                end
            else
                obj.debug(sprintf('Loading metadata from %s\n', obj.pathtofile));
                obj.metadata = imfinfo(obj.pathtofile);
                obj.imagecount(1) = length(obj.metadata);
            end
            obj.datasize(5) = sum(obj.imagecount) / obj.datasize(4);
            obj.debug(sprintf(' - image count : %d\n', sum(obj.imagecount)));
            obj.debug(sprintf(' - number of time points : %d\n', obj.datasize(5)));            
            obj.pixelsize  = [obj.gettagvalue(1, 65006), obj.gettagvalue(1, 65007)];
            obj.magnification = obj.gettagvalue(1, 65008);  
            obj.exposuretime = obj.gettagvalue(1, 65010);  
        end
        
        function idx = gettagindex(obj, n, ID)
            idx = find([obj.metadata(n).UnknownTags(:).ID] == ID);
        end
        
        function value = gettagvalue(obj, n, ID)
            idx = obj.gettagindex(n, ID);
            %idx = ID - 65000 + 1;
            value = obj.metadata(n).UnknownTags(idx).Value;
        end
        
        function obj = parsestageindex(obj)           
           obj.debug('Loading stage XY positions\n');          
           obj.stageindex = zeros(length(obj.metadata),2);
           for n = 1:2:length(obj.metadata);               
               X = obj.gettagvalue(n, 65004);
               Y = obj.gettagvalue(n, 65005);
               obj.stageindex(n,:) =  [X,Y];
               obj.stageindex(n+1,:) =  [X,Y];
               %obj.stageindex(n-1,:) =  obj.stageindex(n,:); % often the value is present only for on of the two channels
           end           
           %obj.stageindex(:,1) = obj.stageindex(:,1) - min(obj.stageindex(:,1));
           %obj.stageindex(:,2) = obj.stageindex(:,2) - min(obj.stageindex(:,2));           
        end
        
        function value = stageposition(obj,frame,channel,axis)
            n = obj.offset(frame,channel);
            value = obj.stageindex(n,axis) * obj.stageunitum;
        end
        
        function R = stagerange(obj)
            R = [min(obj.stageindex);max(obj.stageindex)]' * obj.stageunitum;            
        end
        
        function R = range(obj)            
            R = obj.stagerange() + obj.pixelsize' / obj.magnification .* obj.datasize(1:2)' * [0 1];            
        end

        function obj = parsetimestamps(obj) 
           obj.debug('Loading time stamps\n');                     
           for n = 1:2:length(obj.metadata);               
               obj.timestamps(n) = obj.gettagvalue(n, 65009) / 1000;
               obj.timestamps(n + 1) = obj.timestamps(n);
           end
           obj.timestamps = obj.timestamps - min(obj.timestamps);                      
        end        

        function setstageaxis(obj)        
            % set the axis to match the all stage range
            R = obj.range()';
            axis(R(:));
        end

        function setimageaxis(obj, frame)            
            % set the axis to match the current frame
            n = obj.offset(frame,1);
            axis(obj.stageposition(n,1)' * [1 1] + obj.pixelsize' / obj.magnification .* obj.datasize(1:2)' * [0 1]);            
        end

        function n = length(obj)
            n = obj.datasize(5);
        end

        function n = offset(obj, frame, channel)
            % n = offset(obj, frame, channel)
            % compute the offset inside the tiff stack
            n = (frame - 1) * obj.datasize(4) + channel;
        end

        function im = imread(obj, frame, channel)
            n = obj.offset(frame, channel);
            if iscell(obj.pathtofile)
                k = 1;
                while n > obj.imagecount(k)
                    n = n - obj.imagecount(k);
                    k = k + 1;
                end
                obj.debug(sprintf('Loading image #%d from %s\n', n, obj.pathtofile{k}));
                setDirectory(obj.tif{k}, n);
                im = read(obj.tif{k});
            else
                obj.debug(sprintf('Loading image #%d from %s\n', n, obj.pathtofile));
                setDirectory(obj.tif, n);
                im = read(obj.tif);
            end            
            if obj.flip(channel,1) == true
                im = im(end:-1:1,:);
            end
            if obj.flip(channel,2) == true
                im = im(:,end:-1:1);
            end
        end

        function ref = imref2d(obj, frame, channel)            
            d = obj.stageposition(frame,channel,1:2);
            s =  obj.pixelsize ./ obj.magnification;      
            xl = [d(1) - s(1)*obj.datasize(1)/2, d(1) + s(1)*(obj.datasize(1)/2-1)];
            yl = [d(2) - s(2)*obj.datasize(2)/2, d(2) + s(2)*(obj.datasize(2)/2-1)];            
            ref = imref2d(obj.datasize(1:2), xl, yl);
        end
        
        function etime = elapsedtime(obj, frame, channel)
            etime = obj.timestamps(obj.offset(frame, channel));
        end        

        function imshowpair(obj, frame)
            % display the image at a given frame number in the stage
            % coordinates
            im1 = obj.imread(frame, 1);            
            ref1 = obj.imref2d(frame,1);
            im2 = obj.imread(frame, 2);            
            ref2 = obj.imref2d(frame,2);            
            imshowpair(im1,ref1,im2,ref2,'falsecolor','ColorChannels', [1 2 0]);
        end

        function playmovie(obj, usestageaxis, stepframe)
            if nargin < 2
                usestageaxis = true;
            end
            if nargin < 3
                stepframe = 1;
            end
            for frame = 1:stepframe:obj.length()
                obj.imshowpair(frame);
                %hold on;
                %plot(obj.stageposition(1:2:end,1), obj.stageposition(1:2:end,2));
                %plot(obj.stageposition(2*frame,1), obj.stageposition(2*frame,2),'ro');
                %hold off
                if usestageaxis == true
                    obj.setstageaxis();
                else
                    obj.setimageaxis(frame);
                end
                title(sprintf('Frame #%d / %d [%.2f / %.2f s]', ...
                    frame, obj.length(), obj.timestamps(frame)/1000, ...
                    max(obj.timestamps/1000)));
                drawnow
            end
        end
        
        function stitch(obj,stepframe)
            if nargin < 2
                stepframe = 1;
            end
            for frame = 1:stepframe:obj.length()
                obj.imshowpair(frame);                
                obj.setstageaxis();                
                hold on;
            end
            hold off
        end
        
        function P = pano(obj, channel, stepframe, blendingmode)            
            % P = nt.pano(channel, stepframe, blendingmode)
            %
            % Input:
            %  channel : color channel 
            %  stepframe : frame to skip
            %  blendingmode : if equal to 1, 'average projection' otherwise
            %                 'max projection'
            % Output 
            %  panoramic image for the 'channel' with 1 frame every
            %  'stepframe'
            %
            s = obj.pixelsize / obj.magnification; 
            % determine the size in pixel of covered by all images           
            R = obj.range() ./ (s' * [1 1]); % XY range                     
            sz = ceil(R(:,2) - R(:,1))';                      
            sz = sz(2:-1:1);
            % compute a zoom factor
            if max(sz) > 0.5*max(get(groot,'ScreenSize'))
                alpha = 0.5*max(get(groot,'ScreenSize')) / max(sz);
                sz = ceil(alpha * sz);
            end            
            fprintf(1,'pano size is %d x %d pixel (zoom:%f)\n', sz(1), sz(2), alpha);            
            img = obj.imread(1, channel);
            % two accumulators 
            P = mean(double(img(:))) .* ones(sz);
            N = ones(sz);            
            for frame = 1:stepframe:obj.length()
                % load and scale the image
                im = obj.imread(frame, channel);
                im = imresize(im, alpha);
                % compute the position in the destination image                
                d = alpha * (obj.stageposition(frame,channel,1:2) ./ s(:)' - R(:,1)');
                l = max(1, ceil(d(2:-1:1))); % lower bound
                u = l + size(im) - 1; % upper boud                         
                if blendingmode == 1
                    P(l(1):u(1),l(2):u(2)) = P(l(1):u(1),l(2):u(2)) + double(im);
                    N(l(1):u(1),l(2):u(2)) = N(l(1):u(1),l(2):u(2)) + 1;
                else
                    P(l(1):u(1),l(2):u(2)) = max(P(l(1):u(1),l(2):u(2)) , double(im));
                end            
            end
            if blendingmode == 1
                P(N>0) = P(N>0) ./ N(N>0);
            end
            %P = P - min(P(:));
        end

        function [x,y] = tostagecoords(obj, x, y, frame, channel)
            % [x,y] = tostagecoords(obj, x, y, frame, channel)
            % 
            % convert the coordinates X,Y in pixel in the image to the 
            % coordinates of the real world by taking into account the
            % stage position, the pixel size and the magnification
            %
            d = obj.stageposition(frame, channel,1:2);     
            s = obj.pixelsize / obj.magnification;
            c = obj.datasize(1:2) / 2;
            x = d(1) + s(1) * (x - c(1));
            y = d(2) + s(2) * (y - c(2));
        end
        
        
        function obj = calibratestageunit(obj,stepframe)
            % obj = calibratestageunit(obj,stepframe)
            % calibrate the stage unit 
            % note sequence must be a film of a moving stage only
            if nargin < 2
                stepframe = 10;
            end
            obj.stageunitum = 1;
            clf;
            k = 1;            
            for frame = 1:stepframe:obj.length() - 1;                
                im1 = obj.imread(frame, 1);
                im2 = obj.imread(frame + 1, 1);
                tform = imregcorr(im2single(im1),im2single(im2),'translation');
                I(k,:) = tform.T(3,1:2);
                S(k,:) = obj.stageposition(frame + 1,1,1:2) - obj.stageposition(frame,1,1:2);
                De(k) = norm(tform.T(3,1:2)) * obj.pixelsize(1) / obj.magnification;
                Ds(k) = norm(obj.stageposition(frame + 1,1,1:2) - obj.stageposition(frame,1,1:2));
                imshowpair(im1,circshift(im2,-round(tform.T(3,2:-1:1))),'falsecolor','ColorChannels', [1 2 0]);
                axis on
                title(sprintf('%d/%d',frame,obj.length()))
                drawnow                
                k = k + 1;
            end
            subplot(121)
            plot(De,Ds);
            xlabel('image [um]');
            ylabel('stage [um]');
            axis square
            subplot(122)
            c = mean(De/Ds);
            obj.stageunitum = c;
            plot(I(:,1),I(:,2));hold on; plot(c*S(:,1),c*S(:,2));hold off;
            legend('image','stage')
            xlabel('{\Delta}X [um]');
            ylabel('{\Delta}Y [um]');
            axis equal
            axis square
            fprintf('Stage unit is %.3fum\n', c);
        end

        function debug(obj, msg)
            if obj.warningmode > 0
                fprintf(1, msg);
            end
        end
    end

end
