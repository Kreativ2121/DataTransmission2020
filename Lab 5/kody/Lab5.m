clc
close
clear

fileID = fopen('zad1OX.txt','r');
x= fscanf(fileID,'%f')';
fileID = fopen('zad1Zp.txt','r');
y= fscanf(fileID,'%f')';

w=plot(x,y,'m-','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.7,0,0.7],...
     'MarkerSize',10);
 hold on
title('Wykres Zad 3')
xlabel('Ox')
ylabel('Oy')