clc
close
clear

%fileID = fopen('Zad2ASK.txt','r');
fileID = fopen('DemASK.txt','r');
y= fscanf(fileID,'%f')';

w=plot(y,'m-','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.7,0,0.7],...
     'MarkerSize',10);
 hold on
title('Demodulacja PSK - Mprim')
xlabel('Ox')
ylabel('Oy')