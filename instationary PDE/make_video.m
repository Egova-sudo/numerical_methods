function make_video(data)
v = VideoWriter('newfile','MPEG-4');
v.FrameRate = 30;
v.Quality = 100;
%v.CompressionRatio = 30;
open(v)
[s1,s2,s3] = size(data);
Nx = s1-2; hx = 1/(Nx+1);
Ny = s2-2; hy = 1/(Ny+1);
xx = 0:hx:1; yy = 0:hy:1;
for index = 1:s3
    surf(xx,yy,data(:,:,index));
    zlim([0 1])
    axis manual
    set(gca,'nextplot','replacechildren');
    frame = getframe(gcf);
    writeVideo(v,frame)
end
close(v)
end