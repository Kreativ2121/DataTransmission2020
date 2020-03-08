clc
close
clear

x=-10:0.01:10;
fileID = fopen('data3.txt','r');
y= fscanf(fileID,'%f')';

w=plot(x,y,'m:','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.7,0,0.7],...
     'MarkerSize',10);
 hold on
 plot(0,1,'r*')
title('Wykres zad 1')
xlabel('Ox')
ylabel('Oy')