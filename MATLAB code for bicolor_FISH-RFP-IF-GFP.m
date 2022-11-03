clc;
clear;
close all;
fclose('all');

dis_cutoff = 11;
RFP_cutoff = [5,Inf];
RFP_intensity_cutoff = 300;
RFP_z_series = 3;
GFP_cutoff = [20,Inf];           

scale_bar = 0.11;
height = 3;
se = strel('disk',5);

path = 'E:\HYP\20221101\';
imgs_RFP = dir([path,'*RFP*.tif']);
imgs_GFP = dir([path,'*GFP*.tif']);
imgs_RFP_rm = dir([path,'.*RFP*.tif']);
imgs_GFP_rm = dir([path,'.*GFP*.tif']);
imgs_RFP(ismember({imgs_RFP.name},{imgs_RFP_rm.name})) = [];
imgs_GFP(ismember({imgs_GFP.name},{imgs_GFP_rm.name})) = [];


if ~isequal(numel(imgs_RFP),numel(imgs_GFP))
    error('The number of pictures varies\n');
end

save_path = [path,'result\'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end

RFP_dis = [];

for id = 1:numel(imgs_RFP)
    fprintf('%d / %d\n',id,numel(imgs_RFP));
    tiffInfo = imfinfo([path,imgs_RFP(id).name]);
    z = numel(tiffInfo);
    x = tiffInfo.Width;
    y = tiffInfo.Height;
    step_length = height/(z-1);
    RFP_cutoff_real = [max(RFP_cutoff(1),1),min(RFP_cutoff(2),x*y)];
    GFP_cutoff_real = [max(GFP_cutoff(1),1),min(GFP_cutoff(2),x*y)];
    
    GFP_coor = [];
    GFP_map = zeros(x,y,z);
    for frame = 1:z
        data_GFP = double(imread([path,imgs_GFP(id).name],frame));
        data_GFP = mat2gray(data_GFP);
        data_GFP_bi = imbinarize(data_GFP);
        data_GFP_bi = double(bwareaopen(data_GFP_bi,GFP_cutoff_real(1),4));
        data_GFP_bi = imdilate(data_GFP_bi,se);
        data_GFP_bi = imfill(data_GFP_bi,'holes');
        data_GFP_bi = imerode(data_GFP_bi,se);
        data_GFP_bi = data_GFP_bi-double(bwareaopen(data_GFP_bi,GFP_cutoff_real(2),4));
%         CC = bwconncomp(data_GFP_bi,4);
%         numPixels = cellfun(@numel,CC.PixelIdxList);
%         for j = 1:numel(numPixels)
%             idx = cell2mat(CC.PixelIdxList(j));
%             xx = mod(idx-1,x)+1;
%             yy = ceil(idx/x);
%             if numPixels(j)>=GFP_cutoff(2)
%                 data_GFP_bi(xx,yy) = 0;
%                 continue;
%             end
%         end
        
        GFP_map(:,:,frame) = data_GFP_bi;
        [coor_x,coor_y] = find(data_GFP_bi);
        GFP_coor = [GFP_coor;[coor_x,coor_y,ones(numel(coor_x),1)*frame]];
    end
    GFP_coor(:,1:2) = GFP_coor(:,1:2)*scale_bar;
    GFP_coor(:,3) = GFP_coor(:,3)*step_length;
    imwritestack(GFP_map, [save_path,imgs_GFP(id).name(1:end-4),'-binary.tif'])
    
    RFP_map = zeros(x,y,z);
    for frame = 1:z
        data_RFP = double(imread([path,imgs_RFP(id).name],frame));
        s = sort(reshape(data_RFP,x*y,1),'descend');
        data_RFP(data_RFP>s(round(x*y*1e-5))) = s(round(x*y*1e-5));
        data_RFP_norm = mat2gray(data_RFP);
        data_RFP_bi = get_frontground(data_RFP_norm);
        data_RFP_bi = imfill(data_RFP_bi,'holes');
        data_RFP_bi(data_RFP<RFP_intensity_cutoff) = 0;
        data_RFP_bi = double(bwareaopen(data_RFP_bi,RFP_cutoff_real(1),4))-double(bwareaopen(data_RFP_bi,RFP_cutoff_real(2),4));
        RFP_map(:,:,frame) = data_RFP_bi;
%         [coor_x,coor_y] = find(data_Cy5_bi);
%         RFP_coor = [RFP_coor;[coor_x,coor_y,ones(numel(coor_x),1)*(frame-1)*step_length]];
    end
    imwritestack(RFP_map, [save_path,imgs_RFP(id).name(1:end-4),'-binary2.tif'])
    CC = bwconncomp(RFP_map,18);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    RFP_coor = [];
    for j = 1:numel(numPixels)
        idx = cell2mat(CC.PixelIdxList(j));
        zz = ceil(idx/(x*y));
        yy = ceil((mod(idx-1,x*y)+1)/x);
        xx = mod((mod(idx-1,x*y)+1)-1,x)+1;
%         RG = sqrt(sum(((xx-mean(xx))*scale_bar).^2+((yy-mean(yy))*scale_bar).^2+((zz-mean(zz))*step_length).^2)/numel(xx));
%         R = (3*numel(xx)*(step_length*scale_bar^2)/(4*pi))^(1/3);
        if max(zz)-min(zz)<RFP_z_series
            RFP_map(xx,yy,zz) = 0;
            continue;
        end
%         if R<RG
%             RFP_map(xx,yy,zz) = 0;
%             continue;
%         end
        RFP_coor = [RFP_coor;[mean(xx),mean(yy),mean(zz)]];
    end
    imwritestack(RFP_map, [save_path,imgs_RFP(id).name(1:end-4),'-binary.tif'])
    RFP_coor(:,1:2) = RFP_coor(:,1:2)*scale_bar;
    RFP_coor(:,3) = RFP_coor(:,3)*step_length;
%     scatter3(RFP_coor(:,2),RFP_coor(:,1),RFP_coor(:,3));
%     axis equal

    RFP_distance = zeros(1,size(RFP_coor,1));
    for j = 1:size(RFP_coor,1)
        RFP_distance(j) = sqrt(min(sum((GFP_coor-RFP_coor(j,:)).^2,2)));
    end
    RFP_distance(RFP_distance>dis_cutoff) = [];
    RFP_dis = [RFP_dis,RFP_distance];
end
save([save_path,'result.mat']);
figure
histogram(RFP_dis,'BinWidth',0.5)
xlabel('\mum');
ylabel('count');
legend({'RFP'})
saveas(gcf,[save_path,'result.png']);
if exist([save_path,'result.xlsx'],'file')
    delete([save_path,'result.xlsx']);
end
writecell({'RFP',RFP_dis},[save_path,'result.xlsx']);

function [data_thres] = get_frontground(data)
    kernal_meanfilter = ones(10)/10^2;
    data = mat2gray(data);
    data_dc = imgaussfilt(data,0.8);

    data_c = conv2(data_dc,kernal_meanfilter,'same');
    data_sub = max(data_dc-data_c,0);
    data_sub = mat2gray(data_sub);
    data_thres = data_dc>data_c;
    v = data_sub(data_thres>eps);
    thre = graythresh(v);
    data_thres = data_thres.*(data_sub>2*thre);
    thre = graythresh(data_dc);
    data_thres = data_thres.*(data_dc>2*thre);
end

function imwritestack(stack, filename)
	t = Tiff(filename, 'w');

	tagstruct.ImageLength = size(stack, 1);
	tagstruct.ImageWidth = size(stack, 2);
	tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
	tagstruct.BitsPerSample = 32;
	tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
	tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

	for k = 1:size(stack, 3)
		t.setTag(tagstruct)
		t.write(single(stack(:, :, k)));
		t.writeDirectory();
	end

	t.close();
end
