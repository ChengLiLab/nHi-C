clc;
clear;
close all;
fclose('all');

dis_cutoff = 11;
RFP_cutoff = [20,Inf];
Cy5_cutoff = [5,Inf];
Cy5_intensity_cutoff = 130;
Cy5_z_series = 3;
GFP_cutoff = [5,Inf];
GFP_intensity_cutoff = 150;
GFP_z_series = 3;           

se = strel('disk',5);

scale_bar = 0.11; %nm
height = 3;
path = 'E:\HYP\20221025\';
imgs_Cy5 = dir([path,'*Cy5*.tif']);
imgs_RFP = dir([path,'*RFP*.tif']);
imgs_GFP = dir([path,'*GFP*.tif']);
imgs_Cy5_rm = dir([path,'.*Cy5*.tif']);
imgs_RFP_rm = dir([path,'.*RFP*.tif']);
imgs_GFP_rm = dir([path,'.*GFP*.tif']);
imgs_Cy5(ismember({imgs_Cy5.name},{imgs_Cy5_rm.name})) = [];
imgs_RFP(ismember({imgs_RFP.name},{imgs_RFP_rm.name})) = [];
imgs_GFP(ismember({imgs_GFP.name},{imgs_GFP_rm.name})) = [];

if ~isequal(numel(imgs_Cy5),numel(imgs_RFP),numel(imgs_GFP))
    error('The number of pictures varies\n');
end

save_path = [path,'result\'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end

Cy5_dis = [];
GFP_dis = [];

for id = 1:numel(imgs_RFP)
    fprintf('%d / %d\n',id,numel(imgs_RFP));
    tiffInfo = imfinfo([path,imgs_RFP(id).name]);
    z = numel(tiffInfo);
    x = tiffInfo.Width;
    y = tiffInfo.Height;
    step_length = height/(z-1);
    RFP_cutoff_real = [max(RFP_cutoff(1),1),min(RFP_cutoff(2),x*y)];
    GFP_cutoff_real = [max(GFP_cutoff(1),1),min(GFP_cutoff(2),x*y)];
    Cy5_cutoff_real = [max(Cy5_cutoff(1),1),min(Cy5_cutoff(2),x*y)];
    
    RFP_coor = [];
    RFP_map = zeros(x,y,z);
    for frame = 1:z
        data_RFP = double(imread([path,imgs_RFP(id).name],frame));
        data_RFP = mat2gray(data_RFP);
        data_RFP_bi = imbinarize(data_RFP);
        data_RFP_bi = double(bwareaopen(data_RFP_bi,RFP_cutoff_real(1),4));
        data_RFP_bi = imdilate(data_RFP_bi,se);
        data_RFP_bi = imfill(data_RFP_bi,'holes');
        data_RFP_bi = imerode(data_RFP_bi,se);
        data_RFP_bi = data_RFP_bi-double(bwareaopen(data_RFP_bi,RFP_cutoff_real(2),4));
        
        RFP_map(:,:,frame) = data_RFP_bi;
        [coor_x,coor_y] = find(data_RFP_bi);
        RFP_coor = [RFP_coor;[coor_x,coor_y,ones(numel(coor_x),1)*frame]];
    end
    imwritestack(RFP_map, [save_path,imgs_RFP(id).name(1:end-4),'-binary.tif'])
    RFP_coor(:,1:2) = RFP_coor(:,1:2)*scale_bar;
    RFP_coor(:,3) = RFP_coor(:,3)*step_length;
    
    Cy5_map = zeros(x,y,z);
    for frame = 1:z
        data_Cy5 = double(imread([path,imgs_Cy5(id).name],frame));
        data_Cy5_norm = mat2gray(data_Cy5);
        data_Cy5_bi = get_frontground(data_Cy5_norm);
        data_Cy5_bi = imfill(data_Cy5_bi,'holes');
        data_Cy5_bi(data_Cy5<Cy5_intensity_cutoff) = 0;
        data_Cy5_bi = double(bwareaopen(data_Cy5_bi,Cy5_cutoff_real(1),4))-double(bwareaopen(data_Cy5_bi,Cy5_cutoff_real(2),4));
        
        Cy5_map(:,:,frame) = data_Cy5_bi;
%         [coor_x,coor_y] = find(data_Cy5_bi);
%         Cy5_coor = [Cy5_coor;[coor_x,coor_y,ones(numel(coor_x),1)*(frame-1)*step_length]];
    end
    CC = bwconncomp(Cy5_map,18);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    Cy5_coor = [];
    for j = 1:numel(numPixels)
        idx = cell2mat(CC.PixelIdxList(j));
        zz = ceil(idx/(x*y));
        yy = ceil((mod(idx-1,x*y)+1)/x);
        xx = mod((mod(idx-1,x*y)+1)-1,x)+1;
%         RG = sqrt(sum(((xx-mean(xx))*scale_bar).^2+((yy-mean(yy))*scale_bar).^2+((zz-mean(zz))*step_length).^2)/numel(xx));
%         R = (3*numel(xx)*(step_length*scale_bar^2)/(4*pi))^(1/3);
        if max(zz)-min(zz)<Cy5_z_series
            Cy5_map(xx,yy,zz) = 0;
            continue;
        end
%         if RG>R
%             Cy5_map(xx,yy,zz) = 0;
%             continue;
%         end
        Cy5_coor = [Cy5_coor;[mean(xx),mean(yy),mean(zz)]];
    end
    imwritestack(Cy5_map, [save_path,imgs_Cy5(id).name(1:end-4),'-binary.tif'])
    Cy5_coor(:,1:2) = Cy5_coor(:,1:2)*scale_bar;
    Cy5_coor(:,3) = Cy5_coor(:,3)*step_length;
%     scatter3(Cy5_coor(:,1),Cy5_coor(:,2),Cy5_coor(:,3));
%     axis equal
    
    GFP_map = zeros(x,y,z);
    for frame = 1:z
        data_GFP = double(imread([path,imgs_GFP(id).name],frame));
        data_GFP_norm = mat2gray(data_GFP);
        data_GFP_bi = get_frontground(data_GFP_norm);
        data_GFP_bi = imfill(data_GFP_bi,'holes');
        data_GFP_bi(data_GFP<GFP_intensity_cutoff) = 0;
        data_GFP_bi = double(bwareaopen(data_GFP_bi,GFP_cutoff_real(1),4))-double(bwareaopen(data_GFP_bi,GFP_cutoff_real(2),4));
        
        GFP_map(:,:,frame) = data_GFP_bi;
%         [coor_x,coor_y] = find(data_Cy5_bi);
%         GFP_coor = [GFP_coor;[coor_x,coor_y,ones(numel(coor_x),1)*(frame-1)*step_length]];
    end
    CC = bwconncomp(GFP_map,18);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    GFP_coor = [];
    for j = 1:numel(numPixels)
        idx = cell2mat(CC.PixelIdxList(j));
        zz = ceil(idx/(x*y));
        yy = ceil((mod(idx-1,x*y)+1)/x);
        xx = mod((mod(idx-1,x*y)+1)-1,x)+1;
%         RG = sqrt(sum(((xx-mean(xx))*scale_bar).^2+((yy-mean(yy))*scale_bar).^2+((zz-mean(zz))*step_length).^2)/numel(xx));
%         R = (3*numel(xx)*(step_length*scale_bar^2)/(4*pi))^(1/3);
        if max(zz)-min(zz)<GFP_z_series
            GFP_map(xx,yy,zz) = 0;
            continue;
        end
%         if RG>R
%             GFP_map(xx,yy,zz) = 0;
%             continue;
%         end
        GFP_coor = [GFP_coor;[mean(xx),mean(yy),mean(zz)]];
    end
    imwritestack(GFP_map, [save_path,imgs_GFP(id).name(1:end-4),'-binary.tif'])
    GFP_coor(:,1:2) = GFP_coor(:,1:2)*scale_bar;
    GFP_coor(:,3) = GFP_coor(:,3)*step_length;
%     scatter3(GFP_coor(:,2),GFP_coor(:,1),GFP_coor(:,3));
%     axis equal

    Cy5_distance = zeros(1,size(Cy5_coor,1));
    for j = 1:size(Cy5_coor,1)
        Cy5_distance(j) = sqrt(min(sum((RFP_coor-Cy5_coor(j,:)).^2,2)));
    end
    GFP_distance = zeros(1,size(GFP_coor,1));
    for j = 1:size(GFP_coor,1)
        GFP_distance(j) = sqrt(min(sum((RFP_coor-GFP_coor(j,:)).^2,2)));
    end
    Cy5_distance(Cy5_distance>dis_cutoff) = [];
    GFP_distance(GFP_distance>dis_cutoff) = [];
    Cy5_dis = [Cy5_dis,Cy5_distance];
    GFP_dis = [GFP_dis,GFP_distance];
end
save([save_path,'result.mat']);
figure
histogram(Cy5_dis,'BinWidth',0.5)
hold on
histogram(GFP_dis,'BinWidth',0.5)
xlabel('\mum');
ylabel('count');
legend({'Cy5','GFP'})
saveas(gcf,[save_path,'result.png']);
delete([save_path,'result.xlsx']);
writecell({'Cy5',Cy5_dis;'GFP',GFP_dis},[save_path,'result.xlsx']);

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
