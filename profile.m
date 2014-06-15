close all; clear; clc
Ts = 0.0001;

% 载波
c1 = load('f1.dat');
c2 = load('f2.dat');
subplot(221)
plot(c1(:,1), c1(:,2), 'r')
hold on
plot(c2(:,1), c2(:,2), 'b')
hold off
axis([0, Ts, -5, 5])
xlabel('Time')
ylabel('Magnitude')
title('Carrier')

% 已调信号
cx = load('tx.dat');
subplot(222)
plot(cx(:,1), cx(:,2))
xlabel('Time')
ylabel('Magnitude')
title('Signal Modulated')

% 信道输出信号 
cxpn = load('tx_pn.dat');
subplot(223)
plot(cxpn(:,1), cxpn(:,2))
xlabel('Time')
ylabel('Magnitude')
title('Channel Output')

% 已解调信号
ch_outp1 = load('ch_outp1.dat');
ch_outp2 = load('ch_outp2.dat');
subplot(224)
plot(ch_outp1(:,1), ch_outp1(:,2), 'r')
hold on
plot(ch_outp2(:,1), ch_outp2(:,2), 'b')
hold off
xlabel('Time')
ylabel('Magnitude')
title('Signal Demodulated')

print('profile.png', '-dpng')
close all;

