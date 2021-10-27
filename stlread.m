function [points,Nvecs,nvecpoints] = stlread(file,angleomit)

    if ~exist(file,'file')
        error('file not found');
    end
    
    fid=fopen(file,'r');
    M = fread(fid,inf,'uint8=>uint8');
    fclose(fid);
    
    info = M(85:end);
    
    
    nFaces = typecast(M(81:84),'uint32');
    
    nvecs = zeros(nFaces,6);
    
    verts = zeros(nFaces*3,3);

    for i=0:nFaces-1
        facet = info(50*i+1:50*i+50);
        v1 = [typecast(facet(13:16),'single'),typecast(facet(17:20),'single'),typecast(facet(21:24),'single')];
        v2 = [typecast(facet(25:28),'single'),typecast(facet(29:32),'single'),typecast(facet(33:36),'single')];
        v3 = [typecast(facet(37:40),'single'),typecast(facet(41:44),'single'),typecast(facet(45:48),'single')];
        
        nvecs(i+1,1:3) = [typecast(facet(1:4),'single'),typecast(facet(5:8),'single'),typecast(facet(9:12),'single')];
        nvecs(i+1,4) = (v1(1)+v2(1)+v3(1))/3;
        nvecs(i+1,5) = (v1(2)+v2(2)+v3(2))/3;
        nvecs(i+1,6) = (v1(3)+v2(3)+v3(3))/3;
        verts(3*i+1:3*i+3,1) = [v1(1) v2(1) v3(1)];
        verts(3*i+1:3*i+3,2) = [v1(2) v2(2) v3(2)];
        verts(3*i+1:3*i+3,3) = [v1(3) v2(3) v3(3)];
    end
    nvecs = nvecs(nvecs(:,2)>0,:);
    nvecs = nvecs(atan(nvecs(:,2)./sqrt(nvecs(:,1).^2+nvecs(:,3).^2))>angleomit.*pi./180,:);
    Z = [nvecs(:,5);verts(:,2)];
    X = [nvecs(:,6);verts(:,3)];
    Y = [nvecs(:,4);verts(:,1)];
    nZ = nvecs(:,2);
    nX = nvecs(:,3);
    nY = nvecs(:,1);
    points = [X Y Z];
    Nvecs = [nX nY nZ];
    nvecpoints = [nvecs(:,6) nvecs(:,4) nvecs(:,5)];
end

