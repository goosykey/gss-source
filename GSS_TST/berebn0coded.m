ebn0db1 = 0:20;

ebn0db = repmat(0:20,[4,1]);
esn0db = ebn0db+10.*log10(repmat([1 2 3 4]',[1,21]));

ebn0 = 10.^(0.1*ebn0db);
M = repmat([2 4 8 16]',[1,21]);
m = log2(M);

%BER = 1./m .* erfc(sqrt(ebn0).*sin(pi./M));

trellis = poly2trellis([5 4],[23 35 0; 0 5 13]);
trellis1 = poly2trellis(7,{'1 + x^3 + x^4 + x^5 + x^6', ...
    '1 + x + x^3 + x^4 + x^6'});
spect = distspec(trellis,4);
spect1 = distspec(trellis1,4);
% berub = bercoding(1:10,'conv','hard',2/3,spect); % BER bound


[BER1,SER1] = berawgn(ebn0db(1,:),'psk',4,'nondiff');
BER2 = bercoding(ebn0db(1,:),'block','hard',23,18,3);
BER3 = bercoding(ebn0db(1,:),'block','hard',23,12,7);
BER4 = bercoding(ebn0db(1,:),'conv','hard',2/3,spect); % conv coded 2/3
BER5 = bercoding(ebn0db(1,:),'conv','hard',1/2,spect1); % conv coded 1/2



figure;

% % BER MINE
% semilogy(ebn0db1',BER(1,:)','-c'); hold on; ylim([1e-08, 1]);
% semilogy(ebn0db1',BER(2,:)','-m')
% semilogy(ebn0db1',BER(3,:)','-y')
% semilogy(ebn0db1',BER(4,:)','-k')


% BERAWGN
semilogy(ebn0db1',BER1','-r'); hold on; ylim([1e-08, 1]);
%semilogy(ebn0db1',BER2','-g')
%semilogy(ebn0db1',BER3','-b')
semilogy(ebn0db1',BER4','-g')
semilogy(ebn0db1',BER5','-b')


% % SERAWGN
% semilogy(esn0db(1,:)',SER1',':c')
% semilogy(esn0db(2,:)',SER2',':m')
% semilogy(esn0db(3,:)',SER3',':y')
% semilogy(esn0db(4,:)',SER4',':k')

grid;


legend ('QPSK_u_n_c_o_d_e_d','QPSK_c_o_d_e_d_2_/_3','QPSK_c_o_d_e_d_1_/_2');
