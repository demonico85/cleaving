program calcw

implicit none

integer :: ios,a,iexit,kk,jj
real*8 :: interaction,lambda,b
character(len=30) :: namefile,outfile,arg,serout,namelam
character(len=8) :: namebit


   namebit="inters3."
   namelam="zdir.dat"

    call StripSpaces(namelam)
open(unit=30,file=namelam,status="old",iostat=ios)
if(ios .ne. 0)then
    write(*,*)"ERROR: File zdir.dat not found"
    write(*,*)"Exiting..."
    stop 
end if

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
    write(serout,*)"ave.F.",jj,".out"
    call StripSpaces(serout)
    open(unit=20,file=serout,status="replace")
    write(20,*)
    write(20,*)

    read(10,*)
    read(10,*)
    read(10,*)

    read(30,*)lambda

    iexit=1
    kk=0
do 
    read(10,*,iostat=ios)b
    if(ios .ne. 0)then
        iexit = 0
        exit
    end if
    read(10,*,iostat=ios)
    read(10,*,iostat=ios)a, b, interaction


    if(iexit .eq. 0)exit

    kk=kk+1
   write(20,*)kk,kk,kk,interaction,lambda
end do


 close(10)
 close(20)
end do

close(30)

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


