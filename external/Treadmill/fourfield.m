function frfld = fourfield(caim)
if isfield(caim,'network') && isfield(caim,'fireprop')
   frfld = zeros(size(caim.fireprop.fire,1),6);
   for i = 1:size(caim.fireprop.fire,1)
       %%  network
       if caim.network.netprob(i) >0  
           % run 
           if caim.fireprop.fire(i,2) >0 && caim.fireprop.fire(i,3) == 0 
               frfld(i,1) = 1;
           end
           % run and not run
           if caim.fireprop.fire(i,2) >0 && caim.fireprop.fire(i,3) > 0 
               frfld(i,2) = 1;
           end
           % not run
           if caim.fireprop.fire(i,2) == 0 && caim.fireprop.fire(i,3) > 0 
               frfld(i,3) = 1;
           end
       else %no network
        % run 
           if caim.fireprop.fire(i,2) >0 && caim.fireprop.fire(i,3) == 0 
               frfld(i,4) = 1;
           end
           % run and not run
           if caim.fireprop.fire(i,2) >0 && caim.fireprop.fire(i,3) > 0 
               frfld(i,5) = 1;
           end
           % not run
           if caim.fireprop.fire(i,2) == 0 && caim.fireprop.fire(i,3) > 0 
               frfld(i,6) = 1;
           end
       end
   end
   frfld = sum(frfld,1)/size(frfld,1);
end