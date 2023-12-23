program calcw

implicit none

integer :: ios,ntypes,a,iexit,i,j,rows,matrixdim,jj,b,nargs,kk
integer, pointer :: indx(:)
real*8,pointer :: intermatrix(:,:)
real*8 :: AA,AB,totinteraction,N,lambda,delta,CC,DD,expn,lam,minlam
character(len=30) :: namefile,outfile,arg,serout,namelam,typepair


nargs = iargc()
if(nargs .ne. 2) then
      write(*,*) 'usage: s3work <delta> <n>'
      STOP
end if 

  call getarg( 1, arg )
  read(arg,*) typepair
  
  
  if (typepair .ne. 'C' .and. typepair .ne. 'L') then    
	write(*,*) 'Error: you must insert as argument C/L'
    STOP	
  end if



ntypes=4
matrixdim=ntypes+1 ! ntypes+zero col
rows=ntypes+2 ! ntypes+zero col + index
lambda=0.0



allocate(intermatrix(matrixdim,matrixdim))
allocate(indx(matrixdim))

outfile="work.dat"
    open(unit=11,file=outfile,status="replace",iostat=ios)


    lam=1.0
    minlam=1.0

do jj=1,500
    write(namefile,*)"inters3.",jj,".dat"
    call StripSpaces(namefile)
!write(*,*)jj,namefile
    open(unit=10,file=namefile,status="old",iostat=ios)
    write(namelam,*)"lambda.",jj,".dat"
    call StripSpaces(namelam)
!write(*,*)jj,namefile
    open(unit=30,file=namelam,status="old",iostat=ios)
    if(ios .eq. 0 )then 
    write(serout,*)"ave.F.",jj,".out"
    call StripSpaces(serout)
    open(unit=20,file=serout,status="replace")
    write(20,*)
    write(20,*)

    read(10,*)
    read(10,*)
    read(10,*)

    read(30,*)
    read(30,*)
    read(30,*)a,lambda

    iexit=1
    AA=0.0
    AB=0.0
    N=0.0
    kk=0
do 
    read(10,*,iostat=ios)b
    if(ios .ne. 0)then
        iexit = 0
        exit
    end if
    do i=1,matrixdim
        read(10,*,iostat=ios)a,(intermatrix(i,j),j=1,matrixdim)
        if(ios .ne. 0)then
            iexit = 0
            exit
        end if
    end do


    if(iexit .eq. 0)exit

if (typepair == 'C') then
    AA = AA + intermatrix(4,4) + intermatrix(4,5) + intermatrix(5,5) 
    
else 
 !AA = AA + intermatrix(4,4) + intermatrix(4,5) + intermatrix(5,5) 
end if




 close(10)
 close(20)
 close(30)
end do


contains

    subroutine StripSpaces(string)
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

    end subroutine



function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper



end program




