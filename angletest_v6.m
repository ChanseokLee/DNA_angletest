% Matlab code for included angle measurement
% This code can be used to the structures have one center hinge
% Before start, make an 'input' folder and put original and high-contrast AFM images

clear all
clc
if (~exist('output'))
    mkdir('output')
end

filename = '161013_Z Height_Forward_002';
mkdir(strcat('output\',filename));

I_original = imread(strcat('input\',filename,'.png')); % Original AFM image
I = imread(strcat('input\',filename,'_mod.png')); % High-contrast AFM image

% Make black and white(BW) image
level = graythresh(I);
bw = im2bw(I,level);
bw = bwareaopen(bw,20); % filtering too small particles
bw = padarray(bw,[5 5]);

set(figure, 'Position', [300 200 1500 700]);
sub1 = subplot(1,3,1);
title('Original image');
imshow(I_original);

sub2 = subplot(1,3,2);
title('Before filtering');
imageHandle = imshow(I, 'InitialMagnification', 'fit');
axesHandle = ancestor(imageHandle, 'axes');

% Get the extrema points for each labeled object.
s = regionprops(bw, 'Extrema');
a = regionprops(bw, 'Area');
for k=1:numel(a)
    Area(k)=a(k).Area;
end
% Get most frequent value of particle sizes. Different criteria can be used (e.g. median)
modeArea = mode(Area); 

% Plot the image after filtering
hold(axesHandle, 'on');
for k = 1:numel(s)
   e = s(k).Extrema;
   text(e(1,1), e(1,2), sprintf('%d\n', k), ...
     'Parent', axesHandle, ...
      'Clipping', 'on', ...
      'Color', 'r', ...
      'FontWeight', 'bold');  
end
hold(axesHandle, 'off');
text(10,10,strcat('\color{red}Objects Found:',num2str(numel(a))));

% Filtering large and small objects
LB = modeArea-40;
UB = modeArea+40;
L2 = xor(bwareaopen(bw,LB), bwareaopen(bw,UB));

sub3 = subplot(1,3,3);
title('After filtering');
imageHandle2 = imshow(L2, 'InitialMagnification', 'fit');
axesHandle2 = ancestor(imageHandle2, 'axes');
s2 = regionprops(L2, 'Extrema');
a2 = regionprops(L2, 'Area');
l2 = regionprops(L2, 'Image');
c2 = regionprops(L2, 'Centroid');
pixelList = regionprops(L2, 'PixelList');

for k=1:numel(a2)
    Area2(k)=a2(k).Area;    
end

hold(axesHandle2, 'on');
for k = 1:numel(s2)
   e2 = s2(k).Extrema;
   text(e2(1,1), e2(1,2), sprintf('%d\n', k), ...
     'Parent', axesHandle2, ...
      'Clipping', 'on', ...
      'Color', 'g', ...
      'FontWeight', 'bold');  
end
hold(axesHandle2, 'off');
text(10,10,strcat('\color{green}Objects Found:',num2str(numel(a2))));

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 16 7]);
saveas(gcf,strcat('output\',filename,'\result_overview.png'));
close

% Plot area information
figure
subplot(1,2,1)
histogram(Area, 100)
title('Before filtering')
subplot(1,2,2)
histogram(Area2,100)
title('After filtering')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 10 5]);
saveas(gcf,strcat('output\',filename,'\Area_histogram.png'));
close
fprintf('***** %d structures were found *****\n',numel(l2));

% Individual object analysis
for k=1:numel(l2)
    kk=0;
    idvParticleRaw = l2(k).Image;
    idvParticle = padarray(idvParticleRaw, [5,5]);
    idvPCentroid = regionprops(idvParticle, 'Centroid');
    idvPCentroidCoord = idvPCentroid.Centroid;
    idvPCorner = corner(idvParticle,3);
    
    % Cut individual structure
    idvPixelCoord=pixelList(k).PixelList;
    Max = max(idvPixelCoord);
    Min = min(idvPixelCoord);
    idvParticleOriginal=imcrop(I_original,[Min(1)-10 Min(2)-10 Max(1)-Min(1)+10 Max(2)-Min(2)+10]);
    imwrite(idvParticleOriginal, strcat('output\',filename,'\sample',num2str(k),'.png'));
           
    figure(k)    
    set(figure(k),'Position',[500 500 1000 500]);
    subplot(1,3,1)
    imshow(idvParticleOriginal);
    title(['Particle #',num2str(k),' - original']);
        
    subplot(1,3,2);
    imshow(idvParticle);
    hold on
    plot(idvPCentroidCoord(1),idvPCentroidCoord(2), 'b*')
    plot(idvPCorner(:,1), idvPCorner(:,2), 'r*')
    title(['Particle #',num2str(k)]);      
    hold off
        
    % 1st principal component analysis to find centerlines
    K=idvPixelCoord;
    minX=min(K(:,1));
    maxY=max(K(:,2));
    for i=1:size(K,1)
        K(i,1)=K(i,1)-minX+1;
        K(i,2)=-K(i,2)+maxY+1;
    end        
       
    meanK=mean(K);
    [coeff_K,score_K]=princomp(K);
    maxScoreK1=max(abs(score_K(:,1)));
    maxScoreK2=max(abs(score_K(:,2)));
    v1_pre=coeff_K(:,1);
    v2_pre=coeff_K(:,2);
    x_pca = linspace(-size(bw,2),size(bw,2));
        
    %Two principal axes
    y1_pre = v1_pre(2)/v1_pre(1)*(x_pca-meanK(1))+meanK(2);
    y2_pre = v2_pre(2)/v2_pre(1)*(x_pca-meanK(1))+meanK(2);
    
    %Classifying data
    for i=1:numel(K(:,1))
        if (K(i,2) > v1_pre(2)/v1_pre(1)*(K(i,1)-meanK(1))+meanK(2))
            K(i,3) = 1;            
        else
            K(i,3) = 2;
        end
        if (K(i,2) > v2_pre(2)/v2_pre(1)*(K(i,1)-meanK(1))+meanK(2))
            K(i,4) = 1;            
        else
            K(i,4) = 2;            
        end    
    end
    
    K_img1 = [K(K(:,3)==1,1), K(K(:,3)==1,2)];
    K_img2 = [K(K(:,3)==2,1), K(K(:,3)==2,2)];
    K_img3 = [K(K(:,4)==1,1), K(K(:,4)==1,2)];
    K_img4 = [K(K(:,4)==2,1), K(K(:,4)==2,2)];
        
    img1=[];img2=[];img3=[];img4=[];    
    for i=1:length(K_img1)
        img1(K_img1(i,1),K_img1(i,2))=1;
    end
    for i=1:length(K_img2)
        img2(K_img2(i,1),K_img2(i,2))=1;
    end
    for i=1:length(K_img3)
        img3(K_img3(i,1),K_img3(i,2))=1;
    end
    for i=1:length(K_img4)
        img4(K_img4(i,1),K_img4(i,2))=1;
    end
    
    img1=bwareaopen(img1,10);img2=bwareaopen(img2,10);
    img3=bwareaopen(img3,10);img4=bwareaopen(img4,10);    
    
    numel_img1=numel(regionprops(img1,'Image'));
    numel_img2=numel(regionprops(img2,'Image'));
    numel_img3=numel(regionprops(img3,'Image'));
    numel_img4=numel(regionprops(img4,'Image'));
    
    coord_img1=regionprops(img1,'pixelList');coord_img2=regionprops(img2,'pixelList');
    coord_img3=regionprops(img3,'pixelList');coord_img4=regionprops(img4,'pixelList');
    coord_img1=coord_img1.PixelList;coord_img2=coord_img2.PixelList;
    coord_img3=coord_img3.PixelList;coord_img4=coord_img4.PixelList;
    
    coord_img1=flip(coord_img1,2);coord_img2=flip(coord_img2,2);
    coord_img3=flip(coord_img3,2);coord_img4=flip(coord_img4,2);
    
    resultCheck = 0;
    
    % Choose the correct centerline
    if numel_img1+numel_img2 == 2
       if numel_img3+numel_img4 ~= 2
            K_1=coord_img1;
            K_2=coord_img2;
            ycenter=y1_pre;
       else
            fprintf('Check the result of particle #%d\n',k);            
            resultCheck = 1;
            %if maxScoreK1 < maxScoreK2  
            if maxScoreK1 > maxScoreK2  % Use this criteria when included angle is very small                
                K_1=coord_img1;
                K_2=coord_img2;
                ycenter=y1_pre;
            else
                K_1=coord_img3;
                K_2=coord_img4;
                ycenter=y2_pre;
            end
        end
    else
        if numel_img3+numel_img4 == 2
            K_1=coord_img3;
            K_2=coord_img4;
            ycenter=y2_pre;
        else
            fprintf('Incorrect particle #%d\n',k);
            resultCheck = 1;
            K_1=coord_img1;
            K_2=coord_img2;
            ycenter=y1_pre;
        end
    end              
                  
    % 2nd principal component analysis to calculate bending angle    
    meanK1 = mean(K_1);
    meanK2 = mean(K_2);
    [coeff_K1]=princomp(K_1);
    [coeff_K2]=princomp(K_2);
    v1=coeff_K1(:,1);
    v2=coeff_K2(:,1);
    y1_slope=v1(2)/v1(1);
    y2_slope=v2(2)/v2(1);
    if y1_slope == inf
        y1_slope=1000;
    end
    if y2_slope == inf
        y2_slope=1000;
    end
        
    y1 = y1_slope*(x_pca-meanK1(1))+meanK1(2);
    y2 = y2_slope*(x_pca-meanK2(1))+meanK2(2);
    [intersect_x,intersect_y]=polyxpoly(x_pca,y1,x_pca,y2);
    
    if ~isempty([intersect_x intersect_y])
        v1_mod = meanK1 - [intersect_x intersect_y];
        v1_mod = [v1_mod 0];
        v2_mod = meanK2 - [intersect_x intersect_y];
        v2_mod = [v2_mod 0];
        theta(k) = atan2(norm(cross(v1_mod,v2_mod)),dot(v1_mod,v2_mod))*180/pi();
    else
        theta(k) = 180;        
    end        
               
    subplot(1,3,3);
    plot(K_1(:,1),K_1(:,2),'r.','MarkerSize',12);
    hold on
    plot(K_2(:,1),K_2(:,2),'b.','MarkerSize',12);
    plot(x_pca, ycenter, 'black', 'Linewidth', 1);
    plot(x_pca, y1, 'r', 'Linewidth', 1.5);
    plot(x_pca, y2, 'b', 'Linewidth', 1.5); 
            
    axis equal
    axis ([-10 size(idvParticle,2)+5 -10 size(idvParticle,1)+5]);
    if resultCheck == 0
        title(['PCA analysis : ',num2str(round(theta(k)*100)/100),' deg']);
    else
        title(['PCA analysis : ',num2str(round(theta(k)*100)/100),' deg'],'Color','r');        
    end              
    hold off   
   
    saveas(gcf,strcat('output\',filename,'\particle',num2str(k),'_result.png'));
    fprintf('Particle #%d is finished\n', k)
    
    result(k,:)=[k theta(k) 0];
    if resultCheck == 1
        result(k,3) = 1;
    end              
    clear K
    close    
end

resultFN = strcat('output\',filename,'\result_',filename,'.xlsx');
xlswrite(resultFN, result, 1)
fprintf('***** Analysis completed *****\n');
   

