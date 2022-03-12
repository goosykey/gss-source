ebn0db1 = 0:20;

ebn0db = repmat(0:20,[4,1]);
esn0db = ebn0db+10.*log10(repmat([1 2 3 4]',[1,21]));

ebn0 = 10.^(0.1*ebn0db);
M = repmat([2 4 8 16]',[1,21]);
m = log2(M);

BER = 1./m .* erfc(sqrt(ebn0).*sin(pi./M));

[BER1,SER1] = berawgn(ebn0db(1,:),'psk',2,'nondiff');
[BER2,SER2] = berawgn(ebn0db(1,:),'psk',4,'nondiff');
[BER3,SER3] = berawgn(ebn0db(1,:),'psk',8,'nondiff');
[BER4,SER4] = berawgn(ebn0db(1,:),'psk',16,'nondiff');

figure;

% % BER MINE
% semilogy(ebn0db1',BER(1,:)','-c'); hold on; ylim([1e-08, 1]);
% semilogy(ebn0db1',BER(2,:)','-m')
% semilogy(ebn0db1',BER(3,:)','-y')
% semilogy(ebn0db1',BER(4,:)','-k')


% BERAWGN
semilogy(ebn0db1',BER1','-r'); hold on; ylim([1e-08, 1]);
semilogy(ebn0db1',BER2','-g')
semilogy(ebn0db1',BER3','-b')
semilogy(ebn0db1',BER4','-k')


% % SERAWGN
% semilogy(esn0db(1,:)',SER1',':c')
% semilogy(esn0db(2,:)',SER2',':m')
% semilogy(esn0db(3,:)',SER3',':y')
% semilogy(esn0db(4,:)',SER4',':k')

grid;


legend ('BPSK','QPSK','8PSK','16PSK');
