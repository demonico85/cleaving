program calcw

implicit none

integer :: ios,ntypes,a,iexit,i,j,rows,matrixdim,jj,b,nargs,kk
integer, pointer :: indx(:)
real*8,pointer :: intermatrix(:)
real*8 :: totinteraction,N,lambda,delta,time
character(len=30) :: namefile,outfile,arg,serout
character(len=1) :: sfx


nargs = iargc()
if(nargs .ne. 1) then
      write(*,*) 'usage: s3int <delta>'
      STOP
end if 

  call getarg( 1, arg )
  read(arg,*) delta

!delta=0.01


        write(namefile,*)"wall_a.1.dat"
        call StripSpaces(namefile)

        open(unit=10,file=namefile,status="old",iostat=ios)

            read(10,*)
            read(10,*)
            read(10,*)


            read(10,*)time,a
            
            ntypes=a-1

        close(10)

write(*,*)ntypes

!ntypes=8
matrixdim=ntypes+1 ! ntypes+zero col
rows=ntypes+2 ! ntypes+zero col + index
lambda=0.0



allocate(intermatrix(matrixdim))
allocate(indx(matrixdim))






do kk=1,2
    if(kk == 1)then
        sfx="a"
        outfile="int_a.dat"
    else 
        sfx="b"
        outfile="int_b.dat"
    end if

    open(unit=11,file=outfile,status="replace",iostat=ios)

    do jj=1,500
        write(namefile,*)"wall_",sfx,".",jj,".dat"
        call StripSpaces(namefile)

        open(unit=10,file=namefile,status="old",iostat=ios)
        iexit=1
!write(*,*)jj,namefile,ios
        if(ios .eq. 0 )then 

            read(10,*)
            read(10,*)
            read(10,*)

            do 
                read(10,*,iostat=ios)time,a
                if(ios .ne. 0)then
                    iexit = 0
                    exit
                end if
                do i=1,matrixdim
                    read(10,*,iostat=ios)a,intermatrix(i)
                    if(ios .ne. 0)then
                        iexit = 0
                        exit
                    end if
                end do

            if(iexit .eq. 0)exit

!write(*,*)time
            write(11,*) time*0.001,(intermatrix(i),i=2,matrixdim)

        end do

    else 
        close(10)
        exit
    end if

    close(10)

    end do

    close(11)

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


