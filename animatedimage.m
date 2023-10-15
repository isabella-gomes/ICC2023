function animatedimage(Xesttime,Yesttime,targetpos, userspos)

for t=1:1500
    xreta = linspace(0,targetpos(1,1),100);
    yreta = linspace(0,targetpos(1,2),100);
    plot(Xesttime(:,t),Yesttime(:,t),'ob',targetpos(1,1),targetpos(1,2),'or',xreta,yreta,'r-', userspos(:,1),userspos(:,2),'og');
    ylim([0 1000]);
    xlim([0 1000]);
    grid on;
 
    % gif utilities
    set(gcf,'color','w'); % set figure background to white
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    outfile = 'teste1.gif';
 
    % On the first loop, create the file. In subsequent loops, append.
    if t==1
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
    end
 
end