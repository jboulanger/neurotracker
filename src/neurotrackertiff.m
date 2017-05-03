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
    %  Tag 65004    TAG_TRAY_POSITION_X,  TIFF_LONG,  "Stage Position X component [?]",
    %  Tag 6500596  TAG_TRAY_POSITION_Y, TIFF_LONG,  "Stage Position Y component [?]",
    %  Tag 6500695  TAG_PIXEL_WIDTH,  TIFF_FLOAT,  "Pixel width [micrometer]",
    %  Tag 65007    TAG_PIXEL_HEIGHT, TIFF_FLOAT,  "Pixel height [micrometer]",
    %  Tag 65008    TAG_MAGNIFICATION,  TIFF_SHORT,  "Magnification",
    %  Tag 65009    TAG_TIMESTAMP, TIFF_LONG,  "Timestamp starting from 0 [ms]",
    %  Tag 65010    TAG_EXPOSURE_TIME,  TIFF_LONG,  "Exposure time [ms]"

    properties
        pathtofile;
        datasize = [1 1 1 2 1];
        warningmode;
        metadata;
        imagecount;
        tif;
        stageposition; % Nx2
        limits;
        pixelsize;        
        magnification;
        timestamps;
        exposuretime;
        flip = [true, false];
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
            obj = obj.loadmetadata();
            obj = obj.opentiff();
            obj = obj.loadstagepositions();
            obj.flip = [true,false];
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
            if iscell(obj.tif)
                for n = 1:numel(obj.tif)
                    close(obj.tif{n});
                end
            else
                 close(obj.tif);
            end
        end

        function obj = loadmetadata(obj)
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
            value = obj.metadata(n).UnknownTags(idx).Value;
        end
        
        function obj = loadstagepositions(obj)           
           obj.debug('Loading stage positions\n');          
           obj.stageposition = zeros(length(obj.metadata),2);
           for n = 1:2:length(obj.metadata);
               obj.stageposition(n,:) =  [obj.gettagvalue(n, 65004), obj.gettagvalue(n, 65005)];
               obj.stageposition(n + 1,:) = obj.stageposition(n,:);
               obj.timestamps(n) = obj.gettagvalue(n, 65009);
               obj.timestamps(n + 1) = obj.timestamps(n);
           end
           obj.timestamps = obj.timestamps - min(obj.timestamps);           
           obj.stageposition(:,1) = obj.stageposition(:,1) - min(obj.stageposition(:,1));
           obj.stageposition(:,2) = obj.stageposition(:,2) - min(obj.stageposition(:,2));
           %center = repmat(mean(obj.stageposition), size(obj.stageposition,1), 1);
           %obj.stageposition = 2 * center - obj.stageposition;
        end

        function setstageaxis(obj)
            s = 1e6/obj.pixelsize(1); 
            axis([min(obj.stageposition(:,1)) - s*obj.datasize(1)/2, ...
                  max(obj.stageposition(:,1)) + s*obj.datasize(1)/2, ...
                  min(obj.stageposition(:,2)) - s*obj.datasize(2)/2, ...
                  max(obj.stageposition(:,2)) + s*obj.datasize(2)/2]);
        end

        function setimageaxis(obj, frame)
            n = obj.offset(frame,1);
            axis([obj.stageposition(n,1) - obj.datasize(1)/2, ...
                  obj.stageposition(n,1) + obj.datasize(1)/2, ...
                  obj.stageposition(n,2) - obj.datasize(2)/2, ...
                  obj.stageposition(n,2) + obj.datasize(2)/2]);
        end

        function n = length(obj)
            n = obj.datasize(5);
        end

        function n = offset(obj, frame, channel)
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
            if obj.flip(1) == true
                im = im(end:-1:1,:);
            end
            if obj.flip(2) == true
                im = im(:,end:-1:1);
            end
        end

        function ref = imref2d(obj, frame, channel)
            n = obj.offset(frame, channel);
            d = obj.stageposition(n,:);
            s = 1e6/obj.pixelsize(1); 
            xl = [d(1) - s*obj.datasize(1)/2, d(1) + s*(obj.datasize(1)/2 - 1)];
            yl = [d(2) - s*obj.datasize(2)/2, d(2) + s*(obj.datasize(2)/2 - 1)];            
            ref = imref2d(obj.datasize(1:2), xl, yl);
        end
        
        function etime = elapsedtime(obj, frame, channel)
            etime = obj.timestamps(obj.offset(frame, channel));
        end

        function imshowpair(obj, frame)
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
                hold on;
                plot(obj.stageposition(1:2:end,1), obj.stageposition(1:2:end,2));
                plot(obj.stageposition(2*frame,1), obj.stageposition(2*frame,2),'ro');
                hold off
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
            s = 1e6/obj.pixelsize(1); 
            mini = min(obj.stageposition)/s - obj.datasize(2:-1:1)/2;
            maxi = max(obj.stageposition)/s + obj.datasize(2:-1:1)/2;
            sz = ceil(maxi - mini);
            sz = sz(2:-1:1);            
            if max(sz) > 0.5*max(get(groot,'ScreenSize'))
                alpha = 0.5*max(get(groot,'ScreenSize')) / max(sz);
                sz = ceil(alpha * sz);
            end
            fprintf(1,'pano size is %d x %d pixel (zoom:%f)\n', sz(1), sz(2), alpha);
            img = obj.imread(1, channel);
            P = mean(double(img(:))) .* ones(sz); 
            N = ones(sz);          
            for frame = 1:stepframe:obj.length()
                im = obj.imread(frame, channel);
                im = imresize(im, alpha);
                n = obj.offset(frame, channel);
                d = alpha * obj.stageposition(n,:)/s;  
                l = max(1, ceil(d(2:-1:1)));
                u = l + size(im) - 1;
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
            d = obj.stageposition(obj.offset(frame, channel),:);
            s = 1e6/obj.pixelsize(1); 
            x = s*x + d(1) - s*obj.datasize(1) / 2;
            y = s*y + d(2) - s*obj.datasize(2) / 2;
        end

        function debug(obj, msg)
            if obj.warningmode > 0
                fprintf(1, msg);
            end
        end
    end

end
