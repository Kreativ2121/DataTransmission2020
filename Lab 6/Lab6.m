clc
close
clear

fileID = fopen('zad3Mprim.txt','r');
y= fscanf(fileID,'%f')';

w=plot(y,'m-','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.7,0,0.7],...
     'MarkerSize',10);
 hold on
title('FSK - Zad4 - Widmo')
xlabel('Ox')
ylabel('Oy')