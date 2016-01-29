function SharrockVal = Sharrock(KVkT,H0,R,f0,n)
    SharrockVal = H0 * (1 - ( (1/KVkT)*log(f0*H0/2 / KVkT ./ R)).^(1/n));  
end