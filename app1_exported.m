classdef app1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        ImageDataTab                    matlab.ui.container.Tab
        DicomDataLabel                  matlab.ui.control.Label
        ModalityEditField               matlab.ui.control.EditField
        ModalityEditFieldLabel          matlab.ui.control.Label
        PatientnameEditField            matlab.ui.control.EditField
        PatientnameEditFieldLabel       matlab.ui.control.Label
        PatientageEditField             matlab.ui.control.EditField
        PatientageEditFieldLabel        matlab.ui.control.Label
        BodypartEditField               matlab.ui.control.EditField
        BodypartEditFieldLabel          matlab.ui.control.Label
        ImageDataLabel                  matlab.ui.control.Label
        WidthEditField                  matlab.ui.control.EditField
        WidthEditFieldLabel             matlab.ui.control.Label
        HeightEditField                 matlab.ui.control.EditField
        HeightEditFieldLabel            matlab.ui.control.Label
        TSizeBitsEditField              matlab.ui.control.EditField
        TSizeBitsEditFieldLabel         matlab.ui.control.Label
        BitDepthEditField               matlab.ui.control.EditField
        BitDepthEditFieldLabel          matlab.ui.control.Label
        ImgColorEditField               matlab.ui.control.EditField
        ImgColorEditFieldLabel          matlab.ui.control.Label
        BrowseButton                    matlab.ui.control.Button
        ImagePathEditField              matlab.ui.control.EditField
        ImagePathEditFieldLabel         matlab.ui.control.Label
        UIAxes                          matlab.ui.control.UIAxes
        ImageCoordinatesTab             matlab.ui.container.Tab
        ColorPixelButton                matlab.ui.control.Button
        WhiteButton                     matlab.ui.control.Button
        UIAxes2                         matlab.ui.control.UIAxes
        InterpolationTab                matlab.ui.container.Tab
        Panel_2                         matlab.ui.container.Panel
        UIAxes4                         matlab.ui.control.UIAxes
        ZoomingFactorEditField          matlab.ui.control.NumericEditField
        OriginalImagePanel              matlab.ui.container.Panel
        UIAxes3                         matlab.ui.control.UIAxes
        ZoomingFactorEditFieldLabel     matlab.ui.control.Label
        ImageVersionDropDown            matlab.ui.control.DropDown
        ImageVersionDropDownLabel       matlab.ui.control.Label
        BrowseButton_2                  matlab.ui.control.Button
        ImagePathEditField_2            matlab.ui.control.EditField
        ImagePathEditField_2Label       matlab.ui.control.Label
        HistogramTab                    matlab.ui.container.Tab
        ChooseOperationPanel            matlab.ui.container.Panel
        BrowseButton_3                  matlab.ui.control.Button
        EqualizedImgsNormalizedHistogramButton  matlab.ui.control.Button
        EqualizationButton              matlab.ui.control.Button
        NormalizedHistogramButton       matlab.ui.control.Button
        HistogramEqualizationPanel      matlab.ui.container.Panel
        UIAxes8                         matlab.ui.control.UIAxes
        NormalizedHistogramPanel        matlab.ui.container.Panel
        UIAxes7                         matlab.ui.control.UIAxes
        EqualizedImagePanel             matlab.ui.container.Panel
        UIAxes6                         matlab.ui.control.UIAxes
        ImagePathEditField_3            matlab.ui.control.EditField
        ImagePathEditField_3Label       matlab.ui.control.Label
        ImagePanel                      matlab.ui.container.Panel
        UIAxes5                         matlab.ui.control.UIAxes
        SpatialFilteringTab             matlab.ui.container.Tab
        EnhancedImagePanel              matlab.ui.container.Panel
        UIAxes10                        matlab.ui.control.UIAxes
        OriginalImagePanel_2            matlab.ui.container.Panel
        UIAxes9                         matlab.ui.control.UIAxes
        BrowseButton_4                  matlab.ui.control.Button
        ImagePathEditField_4            matlab.ui.control.EditField
        ImagePathEditField_4Label       matlab.ui.control.Label
        EnterinputsforEnhancementPanel  matlab.ui.container.Panel
        ScaleButton                     matlab.ui.control.Button
        KFactorEditField                matlab.ui.control.NumericEditField
        KFactorEditFieldLabel           matlab.ui.control.Label
        KernelSizeEditField             matlab.ui.control.NumericEditField
        KernelSizeEditFieldLabel        matlab.ui.control.Label
        EnhanceButton                   matlab.ui.control.Button
        FourierTab                      matlab.ui.container.Tab
        ApplyFourierButton              matlab.ui.control.Button
        PhasePanel                      matlab.ui.container.Panel
        UIAxes13                        matlab.ui.control.UIAxes
        MagnitudePanel                  matlab.ui.container.Panel
        UIAxes12                        matlab.ui.control.UIAxes
        ImagePanel_2                    matlab.ui.container.Panel
        UIAxes11                        matlab.ui.control.UIAxes
        BrowseButton_5                  matlab.ui.control.Button
        ImagePathEditField_5            matlab.ui.control.EditField
        ImagePathEditField_5Label       matlab.ui.control.Label
        UIAxes12_2                      matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        whiteImage=255 * ones(20, 20, 3, 'uint8'); % Description
        Image_File;% Description
        GrayImage;
         % Description
        freq;
        normalized
        cdf % Description
        cummlative
        output % Description
        EqualizedImg % Description
        Eqfreq;
        Eqnormalized
        % Description
         % Description
        kernel % Description
        enhanced % Description
        subtraction % Description
        Kfactor % Description
         % Description
    end
    
    
    
    methods (Access = public)
        
        
        
    end
    
    methods (Access = private)
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)
            %Open dialog box for browsing
            [Filename, Pathname,indx] = uigetfile({'*.jpg';'*.bmp';'*.dcm'}, 'All Files');
            %if user canceled browsing 
            if Filename==0
                return;
            end
            %if extensions jpg or bmp
             if (indx==1)||(indx==4)
                fullname=[Pathname,Filename];
                app.ImagePathEditField.Value=fullname;
                app.DicomDataLabel.Visible='off';
                app.ModalityEditFieldLabel.Visible='off';
                app.ModalityEditField.Visible='off';
                app.PatientnameEditFieldLabel.Visible='off';
                app.PatientnameEditField.Visible='off';
                app.PatientageEditFieldLabel.Visible='off';
                app.PatientageEditField.Visible='off';
                app.BodypartEditFieldLabel.Visible='off';
                app.BodypartEditField.Visible='off';
                
                %try & catch error due to corrupted image
                %Display image with its info
                try
                    
                    ImageFile=imread(fullname);
                    imshow(ImageFile,[],'Parent' ,app.UIAxes);
                    [r,c,channels]=size(ImageFile);
                    ImageData=imfinfo(fullname);
                    app.WidthEditField.Value=num2str(ImageData.Width);
                    app.HeightEditField.Value=num2str(ImageData.Height);
                    app.BitDepthEditField.Value=num2str(ImageData.BitDepth);
                    SizeinBits=(ImageData.FileSize)*8;
                    
                    app.TSizeBitsEditField.Value=num2str(SizeinBits);
                    if (ImageData.BitDepth==1 && channels==1 )
                        app.ImgColorEditField.Value='Binary';
                    elseif (ImageData.BitDepth==8 && channels==1)
                        app.ImgColorEditField.Value='Greyscale';
                    else
                       app.ImgColorEditField.Value='RGB'; 
                    end
                    %app.ImgColorEditField.Value=(ImageData.ColorType);
                    

                
                catch ME
                    warndlg('Corrupted Image! Please rebrowse.','Warning');
                    
                    
                    
                
                end
             % if extension dcm (division of if conditions due to library used differently in dicom)    
             elseif (indx==3)
                DicomFullName=[Pathname,Filename];
                app.ImagePathEditField.Value=DicomFullName;
                %try & catch error due to corrupted dicom image
                %Display dicom image with its image and header info
                try
                    DicomImgFile=dicomread(DicomFullName); 
                    imshow(DicomImgFile,[],'Parent',app.UIAxes);
                    DicomImgData=dicominfo(DicomFullName);
                    app.DicomDataLabel.Visible='on';
                    app.ModalityEditFieldLabel.Visible='on';
                    app.ModalityEditField.Visible='on';
                    app.PatientnameEditFieldLabel.Visible='on';
                    app.PatientnameEditField.Visible='on';
                    app.PatientageEditFieldLabel.Visible='on';
                    app.PatientageEditField.Visible='on';
                    app.BodypartEditFieldLabel.Visible='on';
                    app.BodypartEditField.Visible='on';
                    app.WidthEditField.Value=num2str(DicomImgData.Width);
                    app.HeightEditField.Value=num2str(DicomImgData.Height);
                    app.BitDepthEditField.Value=num2str(DicomImgData.BitDepth);
                    SizeinBits=(DicomImgData.FileSize)*8;
                    app.TSizeBitsEditField.Value=num2str(SizeinBits);
                    app.ImgColorEditField.Value=(DicomImgData.ColorType);
                    app.ModalityEditField.Value=char(DicomImgData.Modality);
                    % try and catch some dicominfo that are missing in some
                    % dicomimages
                    try
                        app.PatientageEditField.Value=num2str(DicomImgData.PatientAge);
                        
                    catch ME
                        app.PatientageEditField.Value='Not Found';
                        
                    end
                    try
                        app.BodypartEditField.Value=char(DicomImgData.BodyPartExamined);
                        
                    catch ME
                        app.BodypartEditField.Value='Not Found';
                        
                    end
                    
 
                    try
                        app.PatientnameEditField.Value=struct2array(DicomImgData.PatientName);
                        
                    catch ME
                        app.PatientnameEditField.Value='Not Found';
                        
                    end
                    
                        
                   
                    
                    
                    
                    
                catch ME
                   warndlg('Corrupted Dicom Image! Please rebrowse.','Warning');
                  
                end
                
            
                
                
            end
            
            
                
           
            
            
            
     
       
        end

        % Button pushed function: ColorPixelButton
        function ColorPixelButtonPushed(app, event)

            %get pixel coordinates that need coloring
            app.whiteImage(2,20,:)=[255,0,0];
            app.whiteImage(2,19,:)=[255,0,0];
            app.whiteImage(2,18,:)=[255,0,0];
            app.whiteImage(2,17,:)=[255,0,0];
            app.whiteImage(20,2,:)=[0,0,255];
            app.whiteImage(19,2,:)=[0,0,255];
            app.whiteImage(18,2,:)=[0,0,255];
            app.whiteImage(17,2,:)=[0,0,255];
            imshow(app.whiteImage,[],'Parent',app.UIAxes2);
        end

        % Button pushed function: WhiteButton
        function WhiteButtonPushed(app, event)
          %display rgb white image
           app.whiteImage=255 * ones(20, 20, 3, 'uint8');
           imshow(app.whiteImage,[],'Parent',app.UIAxes2);
    
           
        end

        % Button pushed function: BrowseButton_2
        function BrowseButton_2Pushed(app, event)
            
            [Filename, Pathname] = uigetfile({'*.jpg';'*.bmp';'*.dcm'}, 'All Files');
            %app.UIFigure.Visible = 'off';
            %app.UIFigure.Visible = 'on';
            %if user canceled browsing 
            if Filename==0
                return;
            end
            
            fullname=[Pathname,Filename];
            app.ImagePathEditField_2.Value=fullname;
            
            try
                app.Image_File=imread(fullname);
            catch ME
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end
            
                [rows, columns, numberOfColorChannels] = size(app.Image_File);
          
            
                if numberOfColorChannels==3
                    
                 
                    app.GrayImage= rgb2gray(app.Image_File);
                    app.UIAxes3.Position=[0,0,rows,columns];
               
                  
                    imshow(app.GrayImage,[],'Parent',app.UIAxes3 );
                  
                
                    
 
                else
               
                    app.GrayImage=app.Image_File;
                    app.UIAxes3.Position=[0,0,rows,columns];
                    imshow(app.GrayImage,[],'Parent' ,app.UIAxes3);
                end
                
            
                
                
                
            
        end

        % Value changed function: ImageVersionDropDown
        function ImageVersionDropDownValueChanged(app, event)
            value = app.ImageVersionDropDown.Value;
    
                    
            [row, col] = size(app.GrayImage);
            zoomval=app.ZoomingFactorEditField.Value;
            
            
            
             if (strcmp(value,'Nearest Neighbor Interpolation')==1)
                disp('nearest');
                
                NewRow=round(zoomval*row);
                NewCol=round(zoomval*col);
                %Loop of row
                RowRange = 1 : NewRow;
                AlterRow = round(RowRange * (row - 1) / (NewRow - 1) + (NewRow - row) / (NewRow- 1) );
                %Loop of col
                ColRange = 1 : NewCol;
                AlterCol = round(ColRange * (col - 1) / (NewCol - 1) + (NewCol - col) / (NewCol - 1) );    
                nearestinterpol(RowRange,ColRange) = app.GrayImage(AlterRow,AlterCol);
                app.Panel_2.Title='Nearest Neighbor Interpolation';
                app.UIAxes4.Position=[0,0,NewRow,NewCol];
                imshow(nearestinterpol,[],'Parent' ,app.UIAxes4);
                [npolrow,npolcol]=size(nearestinterpol);
                disp('Old:');
                disp(row);
                disp(col);
                disp('New:');
                disp(npolrow); 
                disp(npolcol);
%             
             elseif (strcmp(value,'Linear Interpolation')==1)

                disp('bilinear');
                NewRow = (zoomval * row);
                NewCol = (zoomval * col);
                LinInterpol = zeros(NewRow,NewCol);
                rowscale = (row/NewRow);
                colscale = (col/NewCol);
                for i=1:NewRow
                    %row index of output 
                     y = (rowscale * i) + (0.5 * (1 - 1/zoomval));
                     for j=1:NewCol 
                        %col index of output
                        x = (colscale * j) + (0.5 * (1 - 1/zoomval));
                        %Any values out of acceptable range
                        x(x < 1) = 1;
                        x(x > row - 0.001) = row - 0.001;
                        x1 = floor(x);
                        x2 = x1 + 1;
                        y(y < 1) = 1;
                        y(y > col - 0.001) = col - 0.001;
                        y1 = floor(y);
                        y2 = y1 + 1;
                        %// 4 Neighboring Pixels
                        NP1 = app.GrayImage(y1,x1);
                        NP2 = app.GrayImage(y1,x2);
                        NP3 = app.GrayImage(y2,x1); 
                        NP4 = app.GrayImage(y2,x2);
                        %// 4 Pixels Weights
                        PW1 = (y2-y)*(x2-x);
                        PW2 = (y2-y)*(x-x1);
                        PW3 = (x2-x)*(y-y1);
                        PW4 = (y-y1)*(x-x1);
                        LinInterpol(i,j) = PW1 * NP1 + PW2 * NP2 + PW3 * NP3 + PW4 * NP4;
                      end
                end

           
%                 img = double(app.GrayImage);
%                 %New rows and columns,Pre allocated memory
%                 NewRow = round(zoomval* row);
%                 NewCol = round(zoomval * col);
%                 %Vectorization cycle,Prevent spillage
%                 RowRange = 1 : NewRow - 1;
%                 ColRange = 1 : NewCol - 1;
%                 %Projection back to position
%                 AlterRow = RowRange * (row - 1) / (NewRow - 1) + (NewRow - row) / (NewRow - 1);
%                 AlterCol = ColRange * (col - 1) / (NewCol - 1) + (NewCol - col) / (NewCol - 1);
%    
%                 %Rounding original coordinates
%                 %Integer part (i)
%                 i_row = floor(AlterRow);
%                 i_col = floor(AlterCol);
%                 %Fraction part(u) around distance of a point
%                 u_row = AlterRow - i_row;
%                 u_col = AlterCol- i_col;
%                 w_row = 1 - u_row;
%                 w_col = 1 - u_col;
% %Determine the weight according to the distance, and get the interpolated image
%                 LinInterpol(RowRange,ColRange) = w_row' * w_col .* img(i_row,i_col) + w_row' * u_col .* img(i_row,i_col + 1)...
%         + u_row' * w_col .* img(i_row + 1,i_col) + u_row' * u_col .* img(i_row + 1,i_col + 1);
%                 %LinInterpol = uint8(LinInterpol);
                
                   app.Panel_2.Title='Linear Interpolation';
                   app.UIAxes4.Position=[0,0,NewRow,NewCol];
                   imshow(LinInterpol,[],'Parent',app.UIAxes4);
                   [LInterRow,LInterCol]=size(LinInterpol);
                   disp('Old:');
                   disp(row);
                   disp(col);
                   disp('New:');
                   disp(LInterRow); 
                   disp(LInterCol);
   
             end
            
        
        
                 
            
                
                
            
            
        end

        % Button pushed function: BrowseButton_3
        function BrowseButton_3Pushed(app, event)
            %Open dialog box for browsing
            [Filename, Pathname] = uigetfile({'*.jpg';'*.bmp';'*.dcm';'*.tif'}, 'All Files');
            %if user canceled browsing 
            if Filename==0
                return;
            end
            fullname=[Pathname,Filename];
            app.ImagePathEditField_3.Value=fullname;
            try
                app.Image_File=imread(fullname);
            catch ME
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end
            
                [rows, columns, numberOfColorChannels] = size(app.Image_File);
          
            
                if numberOfColorChannels==3
                    
                 
                    app.GrayImage= rgb2gray(app.Image_File);
                    app.UIAxes5.Position=[0,0,rows,columns];
                    imshow(app.GrayImage,'Parent',app.UIAxes5 );
                          
 
                else
               
                    app.GrayImage=app.Image_File;
                    app.UIAxes5.Position=[0,0,rows,columns];
                    imshow(app.GrayImage,'Parent' ,app.UIAxes5);
                end
           
        end

        % Button pushed function: NormalizedHistogramButton
        function NormalizedHistogramButtonPushed(app, event)
            app.freq=zeros(256,1);
            app.normalized=zeros(256,1);
            [r,c]=size(app.GrayImage);
       
            for row=1:r
                for col=1:c
                   
        
	                pix=app.GrayImage(row,col);
	                app.freq(pix+1)=app.freq(pix+1)+1;                               %pix plus 1 because index in matlab arrays start at 1%array should have positive values
                    app.normalized(pix+1)=app.freq(pix+1)/(r*c);
                    
                end
            end
        
       
            
            bar(app.UIAxes7,app.normalized);

            
        end

        % Button pushed function: EqualizationButton
        function EqualizationButtonPushed(app, event)
            [r,c]=size(app.GrayImage);
            app.EqualizedImg = uint8(zeros(r,c));
            eqhist=zeros(256,1);
            app.cdf = zeros(256,1);
            app.cummlative = zeros(256,1);
            app.output = zeros(256,1);
            app.normalized=zeros(256,1);
            sum =0 ;

            intensityLevel = 255;
            %get cdf 
            for i = 1:size(app.GrayImage)
                sum =sum +app.freq(i);
                app.cummlative(i) = sum;
                app.cdf(i) = (app.cummlative(i))/(r*c);
                app.output(i) = round(app.cdf(i) * intensityLevel);
                eqhist(app.output(i)+1)= eqhist(app.output(i)+1)+app.normalized(i);
            end
            for i = 1:r
                for j = 1:c
                    app.EqualizedImg(i,j) = app.output(app.GrayImage(i,j) + 1);
                end
            end
           
            app.UIAxes6.Position=[0,0,r,c];
            imshow(app.EqualizedImg,'Parent' ,app.UIAxes6);
            app.HistogramEqualizationPanel.Title='Histogram Equalization';
            bar(app.UIAxes8,eqhist);
        end

        % Button pushed function: 
        % EqualizedImgsNormalizedHistogramButton
        function EqualizedImgsNormalizedHistogramButtonPushed(app, event)
           app.Eqfreq=zeros(1,256);
           app.Eqnormalized=zeros(1,256);
           [row,col]=size(app.EqualizedImg);
           for i=1:row
               for j=1:col
                   pixels=app.EqualizedImg(i,j);
                   app.Eqfreq(pixels+1)= app.Eqfreq(pixels+1)+1;
                   app.Eqnormalized(pixels+1)=(app.Eqfreq(pixels+1))/(row*col);
               end
           end
          
          app.HistogramEqualizationPanel.Title='Equalized Image Normalized Histogram';
          bar(app.UIAxes8,app.Eqnormalized);
         
        end

        % Button pushed function: BrowseButton_4
        function BrowseButton_4Pushed(app, event)
            %Open dialog box for browsing
            [Filename, Pathname] = uigetfile({'*.jpg';'*.bmp';'*.dcm';'*.tif'}, 'All Files');
            %if user canceled browsing 
            if Filename==0
                return;
            end
            fullname=[Pathname,Filename];
            app.ImagePathEditField_4.Value=fullname;
            try
                app.Image_File=imread(fullname);
            catch ME
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end
            
                [rows, columns, numberOfColorChannels] = size(app.Image_File);
          
            
                if numberOfColorChannels==3
                    
                 
                    app.GrayImage= rgb2gray(app.Image_File);
                    app.UIAxes9.Position=[0,0,rows,columns];
                    imshow(app.GrayImage,'Parent',app.UIAxes9 );
                          
 
                else
               
                    app.GrayImage=app.Image_File;
                    app.UIAxes9.Position=[0,0,rows,columns];
                    imshow(app.GrayImage,'Parent' ,app.UIAxes9);
                end
        end

        % Button pushed function: EnhanceButton
        function EnhanceButtonPushed(app, event)

            [row,col]=size(app.GrayImage);
            Output=uint8(zeros(row,col));
            kerValue=app.KernelSizeEditField.Value;
            if(mod(kerValue,2)==0)
                warndlg('Please enter an odd number for kernel size','Warning');
                return;
            end
            try
             app.kernel=(ones(kerValue))/(kerValue*kerValue);
            [krow,kcol]=size(app.kernel);
            if(size(app.kernel)>size(app.GrayImage))
                warndlg('kernel size greater than image size','Warning');
                return;
            end
            kerRowMid = ceil(krow/2); %middle of kernel is first pixel in new img
            kerColMid = ceil(kcol/2);
            zeropad=(kerValue-1)/2;
            Padding=zeros(row+krow-1,col+kcol-1);
            for i=1:row
                for j=1:col
                        Padding(i+zeropad,j+zeropad)=app.GrayImage(i,j);
                end
            end
            %Padding = padarray(app.GrayImage, [kerRowMid - 1, kerColMid - 1], 'replicate', 'both');
            for NewImgRow = kerRowMid : size(Padding,1) - (kerRowMid - 1)
                 for NewImgCol = kerColMid : size(Padding,2) - (kerColMid - 1)
                     sum = 0;
                     for kerRow = 1 : krow
                         for kerCol = 1 : kcol
                             multVal = app.kernel(kerRow, kerCol) * Padding((kerRow - kerRowMid + NewImgRow), (kerCol - kerColMid + NewImgCol)); 
                             sum = sum + multVal;
                         end
                     end
                     
                     Output((NewImgRow - kerRowMid + 1), (NewImgCol - kerColMid + 1)) = sum;
                 end
            end
%            
%             Padding=zeros(row+krow-1,col+kcol-1);
%             for i=1:row
%                 for j=1:col
%                     Padding(i+1,j+1)=app.GrayImage(i,j);
%                 end
%             end
%            
%             for i=1:size(Padding,1)-(krow-1) 
%                 for j=1:size(Padding,2)-(kcol-1)
%                     multVal=Padding(i:i+(krow-1),j:j+(kcol-1))*app.kernel(krow,kcol);
%                     Output(i,j)=sum(multVal(:));
%                 end
%             end
            
            
            app.Kfactor=app.KFactorEditField.Value;
%             en=(imfilter(app.GrayImage,app.kernel));
%             eh=(Kfactor*(app.GrayImage-en))+app.GrayImage;
            app.subtraction=app.GrayImage-Output;
            app.enhanced=(app.Kfactor*(app.GrayImage-Output))+app.GrayImage;
            app.EnhancedImagePanel.Title='Enhanced Image';
            app.UIAxes10.Position=[0,0,row,col];
            imshow(uint8(app.enhanced),'Parent' ,app.UIAxes10);

%             app.UIAxes9.Position=[0,0,row,col];
%             imshow(uint8(eh),'Parent' ,app.UIAxes9);
                
            catch ME
                warndlg('Please browse an image or re-enter odd integers in kernel','Warning');
                
            end
            
           
            
        end

        % Button pushed function: BrowseButton_5
        function BrowseButton_5Pushed(app, event)
            [Filename, Pathname] = uigetfile({'*.jpg';'*.bmp';'*.dcm';'*.tif'}, 'All Files');
            %if user canceled browsing 
            if Filename==0
                return;
            end
            fullname=[Pathname,Filename];
            app.ImagePathEditField_5.Value=fullname;
            try
                app.Image_File=imread(fullname);
            catch ME
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end
            
                [rows, columns, numberOfColorChannels] = size(app.Image_File);
          
            
                if numberOfColorChannels==3
                    
                 
                    app.GrayImage= rgb2gray(app.Image_File);
                    app.UIAxes11.Position=[0,0,rows,columns];
                    imshow(app.GrayImage,'Parent',app.UIAxes11 );
                          
 
                else
               
                    app.GrayImage=app.Image_File;
                    app.UIAxes11.Position=[0,0,rows,columns];
                    imshow(app.GrayImage,'Parent' ,app.UIAxes11);
                end
        end

        % Button pushed function: ApplyFourierButton
        function ApplyFourierButtonPushed(app, event)
            [rows, columns] = size(app.GrayImage);
            magnitude=log(1+abs(fftshift(fft2(app.GrayImage))));
            phase=angle(fftshift(fft2(app.GrayImage)));
            app.UIAxes12.Position=[0,0,rows,columns];
            imshow(magnitude,[],'Parent' ,app.UIAxes12);
            app.UIAxes13.Position=[0,0,rows,columns];
            imshow(phase,[],'Parent' ,app.UIAxes13);
        end

        % Button pushed function: ScaleButton
        function ScaleButtonPushed(app, event)
           [row,col]=size(app.GrayImage);
           %mini=min(app.enhanced);
           %maxi=max(app.enhanced);
           %scaled=floor(((app.enhanced-mini)./(maxi-mini))*255);
           scaled=app.subtraction-min(app.subtraction);
           scaled1=255*(scaled./max(scaled));
           scaled2=scaled1.*app.Kfactor+app.GrayImage;
           app.EnhancedImagePanel.Title='Scaled Enhanced Image';
           app.UIAxes10.Position=[0,0,row,col];
           imshow(scaled2,'Parent' ,app.UIAxes10);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.8 0.8 0.8];
            app.UIFigure.Position = [100 100 787 496];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes12_2
            app.UIAxes12_2 = uiaxes(app.UIFigure);
            title(app.UIAxes12_2, 'Title')
            xlabel(app.UIAxes12_2, 'X')
            ylabel(app.UIAxes12_2, 'Y')
            zlabel(app.UIAxes12_2, 'Z')
            app.UIAxes12_2.Position = [8 360 253 128];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 790 496];

            % Create ImageDataTab
            app.ImageDataTab = uitab(app.TabGroup);
            app.ImageDataTab.Title = 'Image Data';
            app.ImageDataTab.BackgroundColor = [0.9412 0.9412 0.9412];

            % Create UIAxes
            app.UIAxes = uiaxes(app.ImageDataTab);
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [72 158 496 243];

            % Create ImagePathEditFieldLabel
            app.ImagePathEditFieldLabel = uilabel(app.ImageDataTab);
            app.ImagePathEditFieldLabel.HorizontalAlignment = 'right';
            app.ImagePathEditFieldLabel.Position = [44 426 72 22];
            app.ImagePathEditFieldLabel.Text = 'Image Path';

            % Create ImagePathEditField
            app.ImagePathEditField = uieditfield(app.ImageDataTab, 'text');
            app.ImagePathEditField.Position = [119 426 414 23];

            % Create BrowseButton
            app.BrowseButton = uibutton(app.ImageDataTab, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.Position = [547 425 51 24];
            app.BrowseButton.Text = {'Browse'; ''};

            % Create ImgColorEditFieldLabel
            app.ImgColorEditFieldLabel = uilabel(app.ImageDataTab);
            app.ImgColorEditFieldLabel.HorizontalAlignment = 'right';
            app.ImgColorEditFieldLabel.Position = [14 14 58 22];
            app.ImgColorEditFieldLabel.Text = 'Img Color';

            % Create ImgColorEditField
            app.ImgColorEditField = uieditfield(app.ImageDataTab, 'text');
            app.ImgColorEditField.Position = [83 15 74 15];

            % Create BitDepthEditFieldLabel
            app.BitDepthEditFieldLabel = uilabel(app.ImageDataTab);
            app.BitDepthEditFieldLabel.HorizontalAlignment = 'right';
            app.BitDepthEditFieldLabel.Position = [13 35 55 22];
            app.BitDepthEditFieldLabel.Text = 'Bit Depth';

            % Create BitDepthEditField
            app.BitDepthEditField = uieditfield(app.ImageDataTab, 'text');
            app.BitDepthEditField.Position = [83 35 74 15];

            % Create TSizeBitsEditFieldLabel
            app.TSizeBitsEditFieldLabel = uilabel(app.ImageDataTab);
            app.TSizeBitsEditFieldLabel.HorizontalAlignment = 'right';
            app.TSizeBitsEditFieldLabel.Position = [15 56 59 22];
            app.TSizeBitsEditFieldLabel.Text = 'TSize Bits';

            % Create TSizeBitsEditField
            app.TSizeBitsEditField = uieditfield(app.ImageDataTab, 'text');
            app.TSizeBitsEditField.Position = [84 57 74 15];

            % Create HeightEditFieldLabel
            app.HeightEditFieldLabel = uilabel(app.ImageDataTab);
            app.HeightEditFieldLabel.HorizontalAlignment = 'right';
            app.HeightEditFieldLabel.Position = [28 77 40 22];
            app.HeightEditFieldLabel.Text = 'Height';

            % Create HeightEditField
            app.HeightEditField = uieditfield(app.ImageDataTab, 'text');
            app.HeightEditField.Position = [83 77 74 15];

            % Create WidthEditFieldLabel
            app.WidthEditFieldLabel = uilabel(app.ImageDataTab);
            app.WidthEditFieldLabel.HorizontalAlignment = 'right';
            app.WidthEditFieldLabel.Position = [33 98 36 22];
            app.WidthEditFieldLabel.Text = 'Width';

            % Create WidthEditField
            app.WidthEditField = uieditfield(app.ImageDataTab, 'text');
            app.WidthEditField.Position = [84 98 74 15];

            % Create ImageDataLabel
            app.ImageDataLabel = uilabel(app.ImageDataTab);
            app.ImageDataLabel.Position = [61 126 67 22];
            app.ImageDataLabel.Text = 'Image Data';

            % Create BodypartEditFieldLabel
            app.BodypartEditFieldLabel = uilabel(app.ImageDataTab);
            app.BodypartEditFieldLabel.HorizontalAlignment = 'right';
            app.BodypartEditFieldLabel.Visible = 'off';
            app.BodypartEditFieldLabel.Position = [410 28 57 22];
            app.BodypartEditFieldLabel.Text = 'Body part';

            % Create BodypartEditField
            app.BodypartEditField = uieditfield(app.ImageDataTab, 'text');
            app.BodypartEditField.Visible = 'off';
            app.BodypartEditField.Position = [482 28 71 14];

            % Create PatientageEditFieldLabel
            app.PatientageEditFieldLabel = uilabel(app.ImageDataTab);
            app.PatientageEditFieldLabel.HorizontalAlignment = 'right';
            app.PatientageEditFieldLabel.Visible = 'off';
            app.PatientageEditFieldLabel.Position = [400 53 66 22];
            app.PatientageEditFieldLabel.Text = 'Patient age';

            % Create PatientageEditField
            app.PatientageEditField = uieditfield(app.ImageDataTab, 'text');
            app.PatientageEditField.Visible = 'off';
            app.PatientageEditField.Position = [481 53 71 14];

            % Create PatientnameEditFieldLabel
            app.PatientnameEditFieldLabel = uilabel(app.ImageDataTab);
            app.PatientnameEditFieldLabel.HorizontalAlignment = 'right';
            app.PatientnameEditFieldLabel.Visible = 'off';
            app.PatientnameEditFieldLabel.Position = [390 77 76 22];
            app.PatientnameEditFieldLabel.Text = 'Patient name';

            % Create PatientnameEditField
            app.PatientnameEditField = uieditfield(app.ImageDataTab, 'text');
            app.PatientnameEditField.Visible = 'off';
            app.PatientnameEditField.Position = [481 77 71 14];

            % Create ModalityEditFieldLabel
            app.ModalityEditFieldLabel = uilabel(app.ImageDataTab);
            app.ModalityEditFieldLabel.HorizontalAlignment = 'right';
            app.ModalityEditFieldLabel.Visible = 'off';
            app.ModalityEditFieldLabel.Position = [419 98 50 22];
            app.ModalityEditFieldLabel.Text = 'Modality';

            % Create ModalityEditField
            app.ModalityEditField = uieditfield(app.ImageDataTab, 'text');
            app.ModalityEditField.Visible = 'off';
            app.ModalityEditField.Position = [482 102 71 14];

            % Create DicomDataLabel
            app.DicomDataLabel = uilabel(app.ImageDataTab);
            app.DicomDataLabel.Visible = 'off';
            app.DicomDataLabel.Position = [439 126 68 22];
            app.DicomDataLabel.Text = 'Dicom Data';

            % Create ImageCoordinatesTab
            app.ImageCoordinatesTab = uitab(app.TabGroup);
            app.ImageCoordinatesTab.Title = 'Image Coordinates';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.ImageCoordinatesTab);
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.DataAspectRatio = [1 1 1];
            app.UIAxes2.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes2.Position = [218 60 357 231];

            % Create WhiteButton
            app.WhiteButton = uibutton(app.ImageCoordinatesTab, 'push');
            app.WhiteButton.ButtonPushedFcn = createCallbackFcn(app, @WhiteButtonPushed, true);
            app.WhiteButton.Position = [125 336 182 38];
            app.WhiteButton.Text = 'White';

            % Create ColorPixelButton
            app.ColorPixelButton = uibutton(app.ImageCoordinatesTab, 'push');
            app.ColorPixelButton.ButtonPushedFcn = createCallbackFcn(app, @ColorPixelButtonPushed, true);
            app.ColorPixelButton.Position = [465 336 182 38];
            app.ColorPixelButton.Text = 'Color Pixel';

            % Create InterpolationTab
            app.InterpolationTab = uitab(app.TabGroup);
            app.InterpolationTab.Title = 'Interpolation';
            app.InterpolationTab.BackgroundColor = [0.9412 0.9412 0.9412];
            app.InterpolationTab.ForegroundColor = [0.149 0.149 0.149];
            app.InterpolationTab.Scrollable = 'on';

            % Create ImagePathEditField_2Label
            app.ImagePathEditField_2Label = uilabel(app.InterpolationTab);
            app.ImagePathEditField_2Label.HorizontalAlignment = 'right';
            app.ImagePathEditField_2Label.FontName = 'Buxton Sketch';
            app.ImagePathEditField_2Label.FontSize = 16;
            app.ImagePathEditField_2Label.FontWeight = 'bold';
            app.ImagePathEditField_2Label.FontColor = [0.149 0.149 0.149];
            app.ImagePathEditField_2Label.Position = [48 426 68 22];
            app.ImagePathEditField_2Label.Text = 'Image Path';

            % Create ImagePathEditField_2
            app.ImagePathEditField_2 = uieditfield(app.InterpolationTab, 'text');
            app.ImagePathEditField_2.FontName = 'Buxton Sketch';
            app.ImagePathEditField_2.FontSize = 16;
            app.ImagePathEditField_2.FontWeight = 'bold';
            app.ImagePathEditField_2.FontColor = [0.149 0.149 0.149];
            app.ImagePathEditField_2.Position = [119 424 414 25];

            % Create BrowseButton_2
            app.BrowseButton_2 = uibutton(app.InterpolationTab, 'push');
            app.BrowseButton_2.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_2Pushed, true);
            app.BrowseButton_2.BackgroundColor = [0.8 0.8 0.8];
            app.BrowseButton_2.FontName = 'Buxton Sketch';
            app.BrowseButton_2.FontSize = 16;
            app.BrowseButton_2.FontWeight = 'bold';
            app.BrowseButton_2.FontColor = [0.149 0.149 0.149];
            app.BrowseButton_2.Position = [544 425 71 25];
            app.BrowseButton_2.Text = {'Browse'; ''};

            % Create ImageVersionDropDownLabel
            app.ImageVersionDropDownLabel = uilabel(app.InterpolationTab);
            app.ImageVersionDropDownLabel.HorizontalAlignment = 'right';
            app.ImageVersionDropDownLabel.FontName = 'Buxton Sketch';
            app.ImageVersionDropDownLabel.FontSize = 16;
            app.ImageVersionDropDownLabel.FontWeight = 'bold';
            app.ImageVersionDropDownLabel.Position = [311 53 87 22];
            app.ImageVersionDropDownLabel.Text = 'Image Version';

            % Create ImageVersionDropDown
            app.ImageVersionDropDown = uidropdown(app.InterpolationTab);
            app.ImageVersionDropDown.Items = {'Nearest Neighbor Interpolation', 'Linear Interpolation', ''};
            app.ImageVersionDropDown.ValueChangedFcn = createCallbackFcn(app, @ImageVersionDropDownValueChanged, true);
            app.ImageVersionDropDown.BackgroundColor = [0.8 0.8 0.8];
            app.ImageVersionDropDown.Position = [413 53 139 22];
            app.ImageVersionDropDown.Value = '';

            % Create ZoomingFactorEditFieldLabel
            app.ZoomingFactorEditFieldLabel = uilabel(app.InterpolationTab);
            app.ZoomingFactorEditFieldLabel.HorizontalAlignment = 'right';
            app.ZoomingFactorEditFieldLabel.FontName = 'Buxton Sketch';
            app.ZoomingFactorEditFieldLabel.FontSize = 16;
            app.ZoomingFactorEditFieldLabel.FontWeight = 'bold';
            app.ZoomingFactorEditFieldLabel.Position = [45 53 99 22];
            app.ZoomingFactorEditFieldLabel.Text = 'Zooming Factor';

            % Create OriginalImagePanel
            app.OriginalImagePanel = uipanel(app.InterpolationTab);
            app.OriginalImagePanel.BorderType = 'none';
            app.OriginalImagePanel.TitlePosition = 'centertop';
            app.OriginalImagePanel.Title = 'Original Image';
            app.OriginalImagePanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.OriginalImagePanel.FontName = 'Buxton Sketch';
            app.OriginalImagePanel.FontWeight = 'bold';
            app.OriginalImagePanel.Scrollable = 'on';
            app.OriginalImagePanel.FontSize = 14;
            app.OriginalImagePanel.Position = [16 126 301 248];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.OriginalImagePanel);
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.AmbientLightColor = 'none';
            app.UIAxes3.Position = [28 16 230 175];

            % Create ZoomingFactorEditField
            app.ZoomingFactorEditField = uieditfield(app.InterpolationTab, 'numeric');
            app.ZoomingFactorEditField.FontName = 'Buxton Sketch';
            app.ZoomingFactorEditField.FontSize = 16;
            app.ZoomingFactorEditField.FontWeight = 'bold';
            app.ZoomingFactorEditField.BackgroundColor = [0.8 0.8 0.8];
            app.ZoomingFactorEditField.Position = [152 53 100 22];

            % Create Panel_2
            app.Panel_2 = uipanel(app.InterpolationTab);
            app.Panel_2.BorderType = 'none';
            app.Panel_2.TitlePosition = 'centertop';
            app.Panel_2.FontName = 'Buxton Sketch';
            app.Panel_2.FontWeight = 'bold';
            app.Panel_2.Scrollable = 'on';
            app.Panel_2.FontSize = 14;
            app.Panel_2.Position = [341 126 282 248];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.Panel_2);
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.Position = [13 20 232 167];

            % Create HistogramTab
            app.HistogramTab = uitab(app.TabGroup);
            app.HistogramTab.Title = 'Histogram';
            app.HistogramTab.BackgroundColor = [0.651 0.8588 0.7412];

            % Create ImagePanel
            app.ImagePanel = uipanel(app.HistogramTab);
            app.ImagePanel.BorderType = 'none';
            app.ImagePanel.TitlePosition = 'centertop';
            app.ImagePanel.Title = 'Image';
            app.ImagePanel.BackgroundColor = [0.8902 0.9804 0.949];
            app.ImagePanel.FontName = 'Lucida Console';
            app.ImagePanel.Scrollable = 'on';
            app.ImagePanel.Position = [154 226 294 188];

            % Create UIAxes5
            app.UIAxes5 = uiaxes(app.ImagePanel);
            zlabel(app.UIAxes5, 'Z')
            app.UIAxes5.Position = [53 39 166 110];

            % Create ImagePathEditField_3Label
            app.ImagePathEditField_3Label = uilabel(app.HistogramTab);
            app.ImagePathEditField_3Label.HorizontalAlignment = 'right';
            app.ImagePathEditField_3Label.FontName = 'Lucida Console';
            app.ImagePathEditField_3Label.FontSize = 14;
            app.ImagePathEditField_3Label.FontWeight = 'bold';
            app.ImagePathEditField_3Label.Position = [157 439 92 22];
            app.ImagePathEditField_3Label.Text = 'Image Path';

            % Create ImagePathEditField_3
            app.ImagePathEditField_3 = uieditfield(app.HistogramTab, 'text');
            app.ImagePathEditField_3.FontName = 'Lucida Console';
            app.ImagePathEditField_3.FontSize = 14;
            app.ImagePathEditField_3.FontWeight = 'bold';
            app.ImagePathEditField_3.Position = [252 439 512 23];

            % Create EqualizedImagePanel
            app.EqualizedImagePanel = uipanel(app.HistogramTab);
            app.EqualizedImagePanel.BorderType = 'none';
            app.EqualizedImagePanel.TitlePosition = 'centertop';
            app.EqualizedImagePanel.Title = 'Equalized Image';
            app.EqualizedImagePanel.BackgroundColor = [0.8902 0.9804 0.949];
            app.EqualizedImagePanel.FontName = 'Lucida Console';
            app.EqualizedImagePanel.Scrollable = 'on';
            app.EqualizedImagePanel.Position = [482 226 282 187];

            % Create UIAxes6
            app.UIAxes6 = uiaxes(app.EqualizedImagePanel);
            zlabel(app.UIAxes6, 'Z')
            app.UIAxes6.Position = [34 32 164 118];

            % Create NormalizedHistogramPanel
            app.NormalizedHistogramPanel = uipanel(app.HistogramTab);
            app.NormalizedHistogramPanel.BorderType = 'none';
            app.NormalizedHistogramPanel.TitlePosition = 'centertop';
            app.NormalizedHistogramPanel.Title = 'Normalized Histogram';
            app.NormalizedHistogramPanel.BackgroundColor = [0.8902 0.9804 0.949];
            app.NormalizedHistogramPanel.FontName = 'Lucida Console';
            app.NormalizedHistogramPanel.Scrollable = 'on';
            app.NormalizedHistogramPanel.Position = [155 15 293 189];

            % Create UIAxes7
            app.UIAxes7 = uiaxes(app.NormalizedHistogramPanel);
            zlabel(app.UIAxes7, 'Z')
            app.UIAxes7.Position = [13 16 245 155];

            % Create HistogramEqualizationPanel
            app.HistogramEqualizationPanel = uipanel(app.HistogramTab);
            app.HistogramEqualizationPanel.BorderType = 'none';
            app.HistogramEqualizationPanel.TitlePosition = 'centertop';
            app.HistogramEqualizationPanel.Title = 'Histogram Equalization';
            app.HistogramEqualizationPanel.BackgroundColor = [0.8902 0.9804 0.949];
            app.HistogramEqualizationPanel.FontName = 'Lucida Console';
            app.HistogramEqualizationPanel.Scrollable = 'on';
            app.HistogramEqualizationPanel.Position = [481 15 283 188];

            % Create UIAxes8
            app.UIAxes8 = uiaxes(app.HistogramEqualizationPanel);
            zlabel(app.UIAxes8, 'Z')
            app.UIAxes8.Position = [13 16 239 155];

            % Create ChooseOperationPanel
            app.ChooseOperationPanel = uipanel(app.HistogramTab);
            app.ChooseOperationPanel.Title = 'Choose Operation';
            app.ChooseOperationPanel.BackgroundColor = [0.8902 0.9804 0.949];
            app.ChooseOperationPanel.FontName = 'Comic Sans MS';
            app.ChooseOperationPanel.FontWeight = 'bold';
            app.ChooseOperationPanel.Position = [13 15 127 398];

            % Create NormalizedHistogramButton
            app.NormalizedHistogramButton = uibutton(app.ChooseOperationPanel, 'push');
            app.NormalizedHistogramButton.ButtonPushedFcn = createCallbackFcn(app, @NormalizedHistogramButtonPushed, true);
            app.NormalizedHistogramButton.Icon = 'images.png';
            app.NormalizedHistogramButton.BackgroundColor = [0.651 0.8588 0.7412];
            app.NormalizedHistogramButton.FontName = 'Lucida Console';
            app.NormalizedHistogramButton.FontWeight = 'bold';
            app.NormalizedHistogramButton.Position = [14 215 102 58];
            app.NormalizedHistogramButton.Text = {'Normalized'; 'Histogram'};

            % Create EqualizationButton
            app.EqualizationButton = uibutton(app.ChooseOperationPanel, 'push');
            app.EqualizationButton.ButtonPushedFcn = createCallbackFcn(app, @EqualizationButtonPushed, true);
            app.EqualizationButton.BackgroundColor = [0.651 0.8588 0.7412];
            app.EqualizationButton.FontName = 'Lucida Console';
            app.EqualizationButton.FontWeight = 'bold';
            app.EqualizationButton.Position = [12 130 104 58];
            app.EqualizationButton.Text = 'Equalization';

            % Create EqualizedImgsNormalizedHistogramButton
            app.EqualizedImgsNormalizedHistogramButton = uibutton(app.ChooseOperationPanel, 'push');
            app.EqualizedImgsNormalizedHistogramButton.ButtonPushedFcn = createCallbackFcn(app, @EqualizedImgsNormalizedHistogramButtonPushed, true);
            app.EqualizedImgsNormalizedHistogramButton.BackgroundColor = [0.651 0.8588 0.7412];
            app.EqualizedImgsNormalizedHistogramButton.FontName = 'Lucida Console';
            app.EqualizedImgsNormalizedHistogramButton.FontWeight = 'bold';
            app.EqualizedImgsNormalizedHistogramButton.Position = [12 41 104 63];
            app.EqualizedImgsNormalizedHistogramButton.Text = {'Equalized Img''s'; 'Normalized '; 'Histogram'; ''};

            % Create BrowseButton_3
            app.BrowseButton_3 = uibutton(app.ChooseOperationPanel, 'push');
            app.BrowseButton_3.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_3Pushed, true);
            app.BrowseButton_3.Icon = 'icin.jpg';
            app.BrowseButton_3.BackgroundColor = [0.651 0.8588 0.7412];
            app.BrowseButton_3.FontName = 'Lucida Console';
            app.BrowseButton_3.FontSize = 16;
            app.BrowseButton_3.FontWeight = 'bold';
            app.BrowseButton_3.Position = [16 298 100 58];
            app.BrowseButton_3.Text = {'Browse'; ''};

            % Create SpatialFilteringTab
            app.SpatialFilteringTab = uitab(app.TabGroup);
            app.SpatialFilteringTab.Title = 'Spatial Filtering';
            app.SpatialFilteringTab.BackgroundColor = [0.749 0.8784 0.9294];

            % Create EnterinputsforEnhancementPanel
            app.EnterinputsforEnhancementPanel = uipanel(app.SpatialFilteringTab);
            app.EnterinputsforEnhancementPanel.BorderType = 'none';
            app.EnterinputsforEnhancementPanel.TitlePosition = 'centertop';
            app.EnterinputsforEnhancementPanel.Title = 'Enter inputs for Enhancement';
            app.EnterinputsforEnhancementPanel.BackgroundColor = [0.749 0.8784 0.9294];
            app.EnterinputsforEnhancementPanel.FontName = 'SimSun';
            app.EnterinputsforEnhancementPanel.FontWeight = 'bold';
            app.EnterinputsforEnhancementPanel.FontSize = 14;
            app.EnterinputsforEnhancementPanel.Position = [127 12 560 104];

            % Create EnhanceButton
            app.EnhanceButton = uibutton(app.EnterinputsforEnhancementPanel, 'push');
            app.EnhanceButton.ButtonPushedFcn = createCallbackFcn(app, @EnhanceButtonPushed, true);
            app.EnhanceButton.Icon = '1480703.png';
            app.EnhanceButton.BackgroundColor = [0.8706 0.9412 0.9686];
            app.EnhanceButton.FontName = 'SimSun';
            app.EnhanceButton.FontSize = 14;
            app.EnhanceButton.FontWeight = 'bold';
            app.EnhanceButton.Position = [324 45 166 28];
            app.EnhanceButton.Text = 'Enhance';

            % Create KernelSizeEditFieldLabel
            app.KernelSizeEditFieldLabel = uilabel(app.EnterinputsforEnhancementPanel);
            app.KernelSizeEditFieldLabel.HorizontalAlignment = 'right';
            app.KernelSizeEditFieldLabel.FontName = 'SimSun';
            app.KernelSizeEditFieldLabel.FontSize = 14;
            app.KernelSizeEditFieldLabel.FontWeight = 'bold';
            app.KernelSizeEditFieldLabel.Position = [5 51 93 22];
            app.KernelSizeEditFieldLabel.Text = 'Kernel Size';

            % Create KernelSizeEditField
            app.KernelSizeEditField = uieditfield(app.EnterinputsforEnhancementPanel, 'numeric');
            app.KernelSizeEditField.HorizontalAlignment = 'center';
            app.KernelSizeEditField.FontName = 'SimSun';
            app.KernelSizeEditField.FontSize = 14;
            app.KernelSizeEditField.FontWeight = 'bold';
            app.KernelSizeEditField.Position = [136 51 160 22];

            % Create KFactorEditFieldLabel
            app.KFactorEditFieldLabel = uilabel(app.EnterinputsforEnhancementPanel);
            app.KFactorEditFieldLabel.HorizontalAlignment = 'right';
            app.KFactorEditFieldLabel.FontName = 'SimSun';
            app.KFactorEditFieldLabel.FontSize = 14;
            app.KFactorEditFieldLabel.FontWeight = 'bold';
            app.KFactorEditFieldLabel.Position = [7 14 69 22];
            app.KFactorEditFieldLabel.Text = 'K Factor';

            % Create KFactorEditField
            app.KFactorEditField = uieditfield(app.EnterinputsforEnhancementPanel, 'numeric');
            app.KFactorEditField.HorizontalAlignment = 'center';
            app.KFactorEditField.FontName = 'SimSun';
            app.KFactorEditField.FontSize = 14;
            app.KFactorEditField.FontWeight = 'bold';
            app.KFactorEditField.Position = [138 14 159 22];

            % Create ScaleButton
            app.ScaleButton = uibutton(app.EnterinputsforEnhancementPanel, 'push');
            app.ScaleButton.ButtonPushedFcn = createCallbackFcn(app, @ScaleButtonPushed, true);
            app.ScaleButton.Icon = '153-1530625_justice-scale-icon-weighing-scale.png';
            app.ScaleButton.BackgroundColor = [0.8706 0.9412 0.9686];
            app.ScaleButton.FontName = 'SimSun';
            app.ScaleButton.FontSize = 14;
            app.ScaleButton.FontWeight = 'bold';
            app.ScaleButton.Position = [325 14 166 29];
            app.ScaleButton.Text = 'Scale';

            % Create ImagePathEditField_4Label
            app.ImagePathEditField_4Label = uilabel(app.SpatialFilteringTab);
            app.ImagePathEditField_4Label.HorizontalAlignment = 'right';
            app.ImagePathEditField_4Label.FontName = 'SimSun';
            app.ImagePathEditField_4Label.FontSize = 14;
            app.ImagePathEditField_4Label.FontWeight = 'bold';
            app.ImagePathEditField_4Label.Position = [97 425 85 22];
            app.ImagePathEditField_4Label.Text = 'Image Path';

            % Create ImagePathEditField_4
            app.ImagePathEditField_4 = uieditfield(app.SpatialFilteringTab, 'text');
            app.ImagePathEditField_4.FontName = 'SimSun';
            app.ImagePathEditField_4.FontSize = 14;
            app.ImagePathEditField_4.FontWeight = 'bold';
            app.ImagePathEditField_4.Position = [185 425 414 23];

            % Create BrowseButton_4
            app.BrowseButton_4 = uibutton(app.SpatialFilteringTab, 'push');
            app.BrowseButton_4.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_4Pushed, true);
            app.BrowseButton_4.Icon = 'icin.jpg';
            app.BrowseButton_4.BackgroundColor = [0.8706 0.9412 0.9686];
            app.BrowseButton_4.FontName = 'SimSun';
            app.BrowseButton_4.FontSize = 14;
            app.BrowseButton_4.FontWeight = 'bold';
            app.BrowseButton_4.Position = [614 422 95 26];
            app.BrowseButton_4.Text = {'Browse'; ''};

            % Create OriginalImagePanel_2
            app.OriginalImagePanel_2 = uipanel(app.SpatialFilteringTab);
            app.OriginalImagePanel_2.BorderType = 'none';
            app.OriginalImagePanel_2.TitlePosition = 'centertop';
            app.OriginalImagePanel_2.Title = 'Original Image';
            app.OriginalImagePanel_2.BackgroundColor = [0.8706 0.9412 0.9686];
            app.OriginalImagePanel_2.FontName = 'SimSun';
            app.OriginalImagePanel_2.FontWeight = 'bold';
            app.OriginalImagePanel_2.Scrollable = 'on';
            app.OriginalImagePanel_2.FontSize = 14;
            app.OriginalImagePanel_2.Position = [61 143 312 247];

            % Create UIAxes9
            app.UIAxes9 = uiaxes(app.OriginalImagePanel_2);
            zlabel(app.UIAxes9, 'Z')
            app.UIAxes9.Position = [6 14 300 185];

            % Create EnhancedImagePanel
            app.EnhancedImagePanel = uipanel(app.SpatialFilteringTab);
            app.EnhancedImagePanel.BorderType = 'none';
            app.EnhancedImagePanel.TitlePosition = 'centertop';
            app.EnhancedImagePanel.Title = 'Enhanced Image';
            app.EnhancedImagePanel.BackgroundColor = [0.8706 0.9412 0.9686];
            app.EnhancedImagePanel.FontName = 'SimSun';
            app.EnhancedImagePanel.FontWeight = 'bold';
            app.EnhancedImagePanel.Scrollable = 'on';
            app.EnhancedImagePanel.FontSize = 14;
            app.EnhancedImagePanel.Position = [423 143 313 247];

            % Create UIAxes10
            app.UIAxes10 = uiaxes(app.EnhancedImagePanel);
            zlabel(app.UIAxes10, 'Z')
            app.UIAxes10.Position = [7 14 300 185];

            % Create FourierTab
            app.FourierTab = uitab(app.TabGroup);
            app.FourierTab.Title = 'Fourier';

            % Create ImagePathEditField_5Label
            app.ImagePathEditField_5Label = uilabel(app.FourierTab);
            app.ImagePathEditField_5Label.HorizontalAlignment = 'right';
            app.ImagePathEditField_5Label.Position = [44 426 72 22];
            app.ImagePathEditField_5Label.Text = 'Image Path';

            % Create ImagePathEditField_5
            app.ImagePathEditField_5 = uieditfield(app.FourierTab, 'text');
            app.ImagePathEditField_5.Position = [119 426 414 23];

            % Create BrowseButton_5
            app.BrowseButton_5 = uibutton(app.FourierTab, 'push');
            app.BrowseButton_5.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_5Pushed, true);
            app.BrowseButton_5.Position = [547 425 51 24];
            app.BrowseButton_5.Text = {'Browse'; ''};

            % Create ImagePanel_2
            app.ImagePanel_2 = uipanel(app.FourierTab);
            app.ImagePanel_2.TitlePosition = 'centertop';
            app.ImagePanel_2.Title = 'Image';
            app.ImagePanel_2.Scrollable = 'on';
            app.ImagePanel_2.Position = [38 30 356 360];

            % Create UIAxes11
            app.UIAxes11 = uiaxes(app.ImagePanel_2);
            zlabel(app.UIAxes11, 'Z')
            app.UIAxes11.Position = [51 68 255 202];

            % Create MagnitudePanel
            app.MagnitudePanel = uipanel(app.FourierTab);
            app.MagnitudePanel.TitlePosition = 'centertop';
            app.MagnitudePanel.Title = 'Magnitude';
            app.MagnitudePanel.Scrollable = 'on';
            app.MagnitudePanel.Position = [449 217 284 173];

            % Create UIAxes12
            app.UIAxes12 = uiaxes(app.MagnitudePanel);
            zlabel(app.UIAxes12, 'Z')
            app.UIAxes12.Position = [28 10 228 133];

            % Create PhasePanel
            app.PhasePanel = uipanel(app.FourierTab);
            app.PhasePanel.TitlePosition = 'centertop';
            app.PhasePanel.Title = 'Phase';
            app.PhasePanel.Scrollable = 'on';
            app.PhasePanel.Position = [448 29 284 173];

            % Create UIAxes13
            app.UIAxes13 = uiaxes(app.PhasePanel);
            zlabel(app.UIAxes13, 'Z')
            app.UIAxes13.Position = [33 18 228 133];

            % Create ApplyFourierButton
            app.ApplyFourierButton = uibutton(app.FourierTab, 'push');
            app.ApplyFourierButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyFourierButtonPushed, true);
            app.ApplyFourierButton.Position = [634 426 100 22];
            app.ApplyFourierButton.Text = 'Apply Fourier';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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