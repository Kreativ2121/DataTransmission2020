clc
close
clear

fileID = fopen('time.txt','r');
x= fscanf(fileID,'%f')';
fileID = fopen('PSK.txt','r');
y= fscanf(fileID,'%f');


w=plot(x,y,'m-','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.7,0,0.7],...
     'MarkerSize',10);
 hold on
title('Zakodowane PSK')
xlabel('Ox')
ylabel('Oy')