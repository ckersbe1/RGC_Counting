%% RGCcountingSizeSegmentation2

% Calvin J. Kersbergen, 20171117
% ckersbe1@jhmi.edu
% Calabresi Lab
% Johns Hopkins University School of Medicine, Department of Neurology
% Matlab R2017b

% This is a semi-automated cell counting algorithm optimized for 
% nuclear (Brn3a) staining of Retinal Ganglion Cells (RGCs) in flat mount 
% preparations.
% Selections can be from central, middle, or peripheral retina.
% It entails a low-pass filter to reduce noise, grayscale automated or 
% manual thresholding, and identification of cell boundaries. Based on the 
% average cell size, it divides up large clusters of overlapping cells into
% estimated numbers of cells within the cluster (Size Segmentation). 
% Output is the number of cells and cell density in the loaded .tif file. 
% Option for user to self-adjust contrast is available. 
% Can be modified for other applications, but currently optimized for the
% described analysis of RGCs only. 

% This script is published here:
% Jing, J., Smith, M. D., Kersbergen, C. J., Kam, T., Viswanathan, M., 
% Martin, K., Dawson, T. M., Dawson, V. L., Zack, D. J., Whartenby, K. A. 
% & Calabresi, P. A. Glial Pathology and Retinal Neurotoxicity in the 
% Anterior Visual Pathway in Experimental Autoimmune Encephalomyelitis.
% (2019). 
% 

% Additional details available upon request to the Calabresi Lab. 

clear all, clc
close all
%user selects file
[fname pname] = uigetfile({'*.tif';'*.TIF';'*.tiff';'*.TIFF'},'select tiff file for counting');
cd(pname)
rawImage = imread(fname); % opens image
figure
imshow(rawImage) % displays raw image to user in Figure 1
title ('raw image')
xlabel ('select if from center (c), middle (m), or periphery (p)')
fprintf ('select if from center (c), middle (m), or periphery (p)')
done = 0;
fignum = gcf;

while not(done) % selection of area chooses the size threshold the algorithm will use (only for automatic thesholding)
    waitforbuttonpress
    pressed = get(fignum, 'CurrentCharacter');
    if pressed == 'p'
        thresholdLow = 0.3;  % cells smaller than 30% of average cell eliminated because of high peripheral noise
        done = 1;
        fprintf ('\n selected periphery \n')
    elseif pressed == 'm'
        thresholdLow = 0.05; % cells smaller than 5% of average cell eliminated
        done = 1;
        fprintf ('\n selected middle \n')
    elseif pressed == 'c'
        thresholdLow = 0.05;
        done = 1;
        fprintf ('\n selected center \n')
    else
        beep
        disp ('\n not a valid key \n')
    end
end
% convert to grayscale
rawImage = rgb2gray(rawImage); % this may or may not need to be turned on or off depending on the image input (RGB or Uint8)

B = wiener2(rawImage, [7 7]); % low pass filter, this is adjustable. 7x7 provides a nice filter for brn3a RGCs
if pressed == 'p' %periphery selected for automatic threshold - pre-optimized with predicted cell size based on this thresholding
    B = imadjust(B,[0.125 0.9],[]); 
     thresHold2 = adaptthresh(B,0.001); 
     C = imbinarize(B, thresHold2); 
     avgCellSize2 = 230; 
else % center or middle selected for automatic threshold - pre-optimized with predicted cell size based on this thresholding
    thresHold = graythresh(B); 
    C = im2bw(B, thresHold); 
    avgCellSize2 = 180; 
end

hold on

figure %shows raw image next to automatic thresholding image
subplot(1,2,1)
imshow(rawImage)
title('raw image')
subplot(1,2,2)
imshow (C)
title('thresholded image')
xlabel('is this thresholding good? g = good, f = fail')

fprintf ('\n is this thresholding good? g = good, f = fail \n')
done = 0;
done1 = 0;
fignum = gcf;
while not(done)
    waitforbuttonpress
    pressed1 = get(fignum, 'CurrentCharacter');
    if pressed1 == 'g' %automatic thresholding is good, advances to counts
        fprintf( '\n automatic thresholding selected \n')
        done = 1;
    elseif pressed1 == 'f' %automatic thresholding failed, need manual thresholding and estimation of cell size
        
        thresholdLow = 0.3;
        while not(done1)
        fprintf(' \n manual thresholding selected \n')
        B1 = B;
        %manual thresholding
        figure;
        imshow(B)
        xlabel('adjust contrast. press any key once complete')
        hfig = imcontrast; %opens window to adjust contrast
        fprintf('\n press any button once complete \n')
        waitforbuttonpress
        
        window_min = str2double(get(findobj(hfig, 'tag', 'window min edit'), 'String')); %saves low threshold
        window_max = str2double(get(findobj(hfig, 'tag', 'window max edit'), 'String')); %saves high threshold
        Cnew = imadjust(B1,[window_min/255 window_max/255],[0 1]); % converts from 0-255 to 0-1
        thresHold2 = adaptthresh(Cnew,0.1); %thresholding algorithm
        Cnew = imbinarize(Cnew, thresHold2); %binarizes
        %convert modified image to workable form
        
        subplot(1,3,1) %plots raw, automatic, and manual side by side 
        imshow(rawImage)
        title('raw image')
        subplot(1,3,2)
        imshow(C)
        title('automatic threshold')
        subplot(1,3,3)
        imshow(Cnew)
        title('manual threshold')
        xlabel('Is this threshold better? Press "y" if good, "n" to repeat')
        
        fprintf('\n Is this threshold better? Press any button to advance \n')
        fignum1 = gcf;
        waitforbuttonpress
        pressed2 = get(fignum1, 'CurrentCharacter');
            if pressed2 == 'y' %thresholding good, advances to counts
                done1 = 1;
                lowThresh =  window_min %displays thresholds used to user
                highThresh = window_max %displays thresholds used to user
            elseif pressed2 == 'n' %thresholding not good, loops until user is satisfied
            
                done1 = 0;
            end
        
        end
        
        done = 1;
        C = Cnew;
    end
end

[B,L,N] = bwboundaries(C,'noholes'); % determines boundaries of cells B = boundary
imgStats = regionprops(C,'area','centroid'); % gives center and area of each enclosed boundary

figure; %plots cells with boundaries and labels, indicated below
imshow(C); 
title('counted cells')
xlabel('Red+ = 1 cell. Blue* = multiple cells. green* = removed for size. yellow* = removed for shape')
hold on;
       for k = 1:length(B)
          plot(imgStats(k).Centroid(1), imgStats(k).Centroid(2),'r+') %plots all cells with red + sign
           boundary = B{k};
           if(k > N)
               plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
           else
               plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
           end
       end
       
% size filter 

i = 1;
if pressed1 == 'f' % if manual thresholding, performs automatic determination of area because not optimized (optimized for automatic already)
    for k = 1:N 
    if imgStats(k).Area < 400 && imgStats(k).Area > 20 % uses range to determine average cell size - excluding very small or very large cells
       
    avgCellSize(i) = imgStats(k).Area;
     i = i + 1;
    else
    end
    end
    avgCellSize2 = mean(avgCellSize(:));
end

count = N; % preliminary cell count without size segmentation
sizeCellremoved = 0; % count of cells removed

 for k = 1:length (B)
      if pressed == 'p'
      numCell = 1;
   
     if imgStats(k).Area < thresholdLow*avgCellSize2 
         % cells smaller than threshold of 0.3 - optimized. Could also be
         % performed based on variance of cell size or distribution, but
         % not necessarily gaussian. 
         count = count - 1;
         sizeCellremoved = sizeCellremoved + 1;
         plot(imgStats(k).Centroid(1), imgStats(k).Centroid(2),'g*')
     elseif imgStats(k).Area > (1.75*avgCellSize2)
         % cells larger than 1.75 times an average cell. Could also be
         % performed based on variance of cell size or distribution, but
         % not necessarily gaussian.  
         numCell = round(imgStats(k).Area / avgCellSize2);
         if numCell > 1
             count = count + (numCell - 1);
             plot(imgStats(k).Centroid(1), imgStats(k).Centroid(2),'b*')
         else
         end

     end
      elseif pressed == 'c' || pressed == 'm'
      numCell = 1;
    
       if imgStats(k).Area < thresholdLow*avgCellSize2 
         % cells smaller than threshold of 0.05 - optimized. Could also be
         % performed based on variance of cell size or distribution, but
         % not necessarily gaussian. 
         count = count - 1;
         sizeCellremoved = sizeCellremoved + 1;
         plot(imgStats(k).Centroid(1), imgStats(k).Centroid(2),'g*') %plots cells eliminated
     elseif imgStats(k).Area > (1.75*avgCellSize2) 
         % cells larger than 1.75 times an average cell. Could also be
         % performed based on variance of cell size or distribution, but
         % not necessarily gaussian. 
         numCell = round(imgStats(k).Area / avgCellSize2);
         if numCell > 1
             count = count + (numCell - 1);
             plot(imgStats(k).Centroid(1), imgStats(k).Centroid(2),'b*') %plots cells split by size segmentation
         else
         end
     
        end  
    end

    
end
 
if (count/k > 1.5) || sizeCellremoved/N > 0.2 % this appears if there is a lot of noise removed or a lot of "large" cells in the image. 
    fprintf('\n WARNING: high noise in image  \n')
elseif avgCellSize2 < 160 || avgCellSize2 > 240 % if large differences from normal RGC cell size, thresholding may be inaccurate.
    fprintf('\n WARNING: check cell size segmentation - cells may be too small or too big\n')
else
end

% scaling of imported image
% 0.496 um/pixel in all images unless new microscope, new objective, etc. 
scaleF = 0.496;
scale1 = size(rawImage);
areaMM = ((scale1(1)*scaleF)/1000)*((scale1(2)*scaleF)/1000);

cellDensityMM = count/areaMM;
count
cellDensityMM