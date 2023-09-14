      program fileio
      implicit none
      integer i
      open(20,file="cubes.dat")

      do i =1,100
        write(20,*)i,i*i,i**3
      end do
      
      close(20)
      end program fileio
