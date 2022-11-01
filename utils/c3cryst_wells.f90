program calcw

implicit none

integer :: Zzero,ios,ntypes,a,iexit,i,j,rows,matrixdim,jj,b,nargs,kk,filen
integer, pointer :: indx(:)
real*8,pointer :: intermatrix(:,:)
real*8 :: AA,AB,totinteraction,N,lambda,delta,CC,DD,expn,lam,minlam
character(len=30) :: namefile,outfile,arg,serout,namelam
character(len=8) :: namebit

!nargs = iargc()
!if(nargs .ne. 2) then
!      write(*,*) 'usage: s3work <delta> <n>'
!      STOP
!end if 

!  call getarg( 1, arg )
!  read(arg,*) delta

!call getarg( 2, arg )
!  read(arg,*) expn

!delta=0.01





ntypes=2
matrixdim=ntypes+1 ! ntypes+zero col
rows=ntypes+2 ! ntypes+zero col + index
lambda=0.0


namefile="inters3.1.dat"
open(unit=10,file=namefile,status="old",iostat=ios)
if(ios .ne. 0)then
   namefile="inters2.1.dat"
   open(unit=10,file=namefile,status="old",iostat=ios)
   if(ios .ne. 0)then
      namefile="inters1.1.dat"
      open(unit=10,file=namefile,status="old",iostat=ios)
      if(ios .ne. 0)then
        write(*,*)"ERROR: unkown .dat files"
        write(*,*)"Exiting...."
        stop
      else
         filen=1
      end if
   else
      filen=2
   end if
else
   filen=3
end if 


select case (filen)
case (1)
   namebit="inters1."
case (2)
   namebit="inters2."
case (3)
   namebit="inters3."
end select
!write(*,*)namebit

allocate(intermatrix(matrixdim,matrixdim))
allocate(indx(matrixdim))

outfile="work.dat"
    open(unit=11,file=outfile,status="replace",iostat=ios)

!write(*,*)int(expn),expn

!if (int(expn) .eq. 1)then
    lam=1.0
    minlam=1.0
!else if(int(expn) .eq. 2)then
!    lam=2.0*lambda
!    minlam=2.0*(1.0-lambda)
!else
!    write(*,*)"Error: exponent not implemented yet"
!    stop
!end if

do jj=1,500
    write(namefile,*)namebit,jj,".dat"
    call StripSpaces(namefile)
!write(*,*)jj,namefile
    open(unit=10,file=namefile,status="old",iostat=ios)
    if(ios .ne. 0)then
        write(namefile,*)namebit,jj-1,".dat"
        call StripSpaces(namefile)
        write(*,*)"Last read file: ",namefile
        exit
    end if
    write(namelam,*)"shiftz.",jj,".dat"
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

    DD = intermatrix(2,2) + intermatrix(2,3) + intermatrix(3,3) 


    kk=kk+1
   write(200,*)kk,kk,kk,DD,lambda

end do

rewind(200)

        do 
           read(200,*,iostat=ios)kk,kk,kk,DD,lambda
           if(ios .ne. 0)exit
           if(kk .eq. 1) Zzero=DD
           write(20,*) kk,kk,kk,DD-Zzero,lambda
        end do

        close(200,status="delete")




   end if
 close(20)
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






end program


