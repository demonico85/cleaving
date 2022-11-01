program trapezoid2

  implicit none

  integer, parameter :: n=10

  integer :: i,nargs,system
  integer :: ios,nlines
  real*8 :: w1,w2,w3,w4,v1,v2,v3,v4,totvar,prefactor

  character(len=30) :: namefile
 character(len=50) :: arg


system = -10 ! must be specified

! system = 0 -> wAlls
! system = 1 -> wElls


    nargs = iargc()
    if(nargs .eq. 0) then
        stop '4 name file needed'
    end if

  write (*,*)

  call getarg( 1, arg )
  read(arg,*) namefile
  call trapezoid_integration(namefile,w1,v1)

  call getarg( 2, arg )
  read(arg,*) namefile
  call trapezoid_integration(namefile,w2,v2)

  call getarg( 3, arg )
  read(arg,*) namefile
  call trapezoid_integration(namefile,w3,v3)


  call getarg( 4, arg )
  read(arg,*) namefile
  call trapezoid_integration(namefile,w4,v4)



    i=1
    do 
        call getarg( i, arg )
        if(arg .ne. '-wa') then
            i=i+1
            if(i .gt. nargs)then
                exit
            end if
        else
            system = 0
            exit
        end if
    end do
    
    i=1
    do 
        call getarg( i, arg )
        if(arg .ne. '-we') then
            i=i+1
            if(i .gt. nargs)then
                exit
            end if
        else
            system = 1
            exit
        end if
    end do

if (system .lt. 0 ) stop "Input walls/wells missing"

if(system .eq. 0)then
    prefactor = -1.0
else
    prefactor = 1.0
end if

    w1=prefactor*w1
    w2=prefactor*w2
    w4=prefactor*w4

  write (*,'(A,2X,F7.4,1X,A,1X,F7.5)') 'Step1 = ',w1, " +- ",v1
  write (*,'(A,2X,F7.4,1X,A,1X,F7.5)') 'Step2 = ',w2, " +- ",v2
  write (*,'(A,2X,F9.6,1X,A,1X,F9.7)') 'Step3 = ',w3, " +- ",v3
  write (*,'(A,2X,F7.4,1X,A,1X,F7.5)') 'Step4 = ',w4, " +- ",v4

  totvar = v1**2 + v2**2 + v3**2 + v4**2
  totvar = sqrt(totvar)

  write (*,*)
  write (*,'(A,2X,F7.4,1X,A,1X,F7.5)')"The interfacial Free Energy is: ",w1+w2+w3+w4, " +- ",totvar
  write (*,*)

  contains

    subroutine trapezoid_integration(nf,integral,var)
      implicit none
      integer :: n
      real ::  end_val
      real*8,intent(out) :: integral, var
      real*8,pointer :: work(:),x(:),dx(:)
      integer :: i,jj
        character(len=30),intent(in) ::nf 

  open(unit=10,file=nf,status="old",iostat=ios)
  if(ios .ne. 0)stop 'File does not exist' 
    n = 0
    
    do 
        read(10,*,iostat=ios)
        if(ios .ne. 0)exit
        n=n+1
    end do
    rewind(10)

    allocate(work(n),x(n),dx(n))
    do i=1,n
        read(10,*,iostat=ios)x(i),work(i),dx(i)
    end do


      integral = 0.0
      do i=1,n-1
         integral = integral+(x(i+1)-x(i))*(work(i+1)+work(i))

      end do
    
    integral=0.5*integral


     var = 0.0
     do i=1,n-1
        var = var+((x(i+1)-x(i))*0.5)**2 * (dx(i+1)**2+dx(i)**2)
     end do

    var = sqrt(var)

    close(10)
    end subroutine trapezoid_integration

end program trapezoid2
