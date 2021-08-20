function statePrint(iter,iter_h,iter_d,objInfo,maxCons,dIndi,lamCri,sig)
disp('Iteration step')
disp([' It.: ' sprintf('%4i\t',iter) ' It_h.: ' sprintf('%4i\t',iter_h) 'It_d'...
 sprintf('%3i\t',iter_d)]);
disp('Current value')
disp([ ' Obj.: ' sprintf('%6.3f\t',objInfo.val)  ...
	' cons: ' sprintf('%4i\t',maxCons)]);
disp('Parameter state')
disp([ ' dIndi.: ' sprintf('%6.3f\t',dIndi)  ...
	' lamMax: ' sprintf('%4i\t',lamCri) ' Sig: ' sprintf('%4i\t',sig)]);