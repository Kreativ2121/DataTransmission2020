clc
close
clear

fileID = fopen('Time.txt','r');
x= fscanf(fileID,'%f')';
fileID = fopen('DecTTL.txt','r');
y= fscanf(fileID,'%f')';

w=plot(x, y,'m-','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.7,0,0.7],...
     'MarkerSize',10);
 hold on
title('Odkodowane TTL')
xlabel('Ox')
ylabel('Oy')