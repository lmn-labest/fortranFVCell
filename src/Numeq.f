c *********************************************************************      
c * NUMEQ: numero de equacoes                                         *
c *********************************************************************      
       subroutine numeq(id,num,elconn,numel,nshared,neq)
       implicit none
       integer id(nshared,*),elconn(nshared,*),num(*)
       integer i,j,aux,numel,nshared,nelviz,neq,n
c ...
       neq = numel
       do i = 1, numel
         n = num(i)
         do j = 1, nshared
           nelviz  = elconn(j,i)
           if( nelviz .ne.-1) then 
             id(j,n) = num(nelviz) 
           else
             id(j,n) = -1
           endif
         enddo        
       enddo
c .....................................................................
c
c ...
       return
       end
c *********************************************************************
