function E2=BPM_steg(E1,delta_z,N,a,lambda_noll,nmat,daempmat)
 k_noll = 2*pi/lambda_noll;
 n_PAS = sum(sum(nmat.*abs(E1).^2)*a)*a/(sum(sum(abs(E1).^2)*a)*a);
 E2_PAS = PAS(E1,delta_z,N,a,lambda_noll,n_PAS);
 faskorrektion = k_noll.*(nmat-n_PAS).*delta_z;
 if max(max(abs(faskorrektion)))>(2*pi*0.02) % max tillåten faskorrektion är 2% av 2pi
     disp('Steglängd för stor')
 return 
 else
 E2=E2_PAS.*exp(1i*faskorrektion).*daempmat;
 end
 % keyboard % avbryter körningen
end