function afterProcessingImages=myImageProcessingFL(beforeProcessingImages)
%This function can segreate orignal images as you want and returen the mask images
%   beforeProcessingImages: Fluorescent images stack;
%   afterProcessingImages: Binary images stack;
imageType='uint16'; %here you can chenge your image type
gaussianFilter=fspecial('gaussian',[3, 3],10); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
grayThresh=150;
areaThreshold=10;
maxIntensity=600;
for iframe=1:size(beforeProcessingImages,3)
    % adjust the images contraast
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imbinarize(afterProcessingImages(:,:,iframe),grayThresh/(2^16-1));%1.2 for the image with few particles
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image=afterProcessingImages(:,:,iframe);
            image(cc(iCC).PixelIdxList)=0;
            afterProcessingImages(:,:,iframe)=image;
        end
    end
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes'); % fill holes process
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % find outLine process
end
afterProcessingImages=logical(afterProcessingImages);
end