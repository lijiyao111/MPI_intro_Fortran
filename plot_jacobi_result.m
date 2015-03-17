clear all

source=load('source.txt');
field=load('field.txt');

figure(1)
%surf(source)
imagesc(source)
%set(gca,'zdir','reverse')
axis equal
xlim([1 72])
title('Source (density anomaly)')

export_fig source.pdf -transparent

figure(2)
%surf(field)
imagesc(field)
axis equal
xlim([1 72])
title('Field (gravity anomaly)')

export_fig field.pdf -transparent

load fullfield

diff=fullfield-field;
figure(4)
surf(diff)
