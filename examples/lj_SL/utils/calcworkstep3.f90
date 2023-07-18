program calcw

implicit none

integer :: ios,ntypes,a,iexit,i,j,rows,matrixdim,jj,b,nargs,kk
real*8,pointer :: intermatrix(:,:)
real*8 :: AA,AB,totinteraction,N,lambda,delta,CC,DD,expn,lam,minlam
character(len=30) :: namefile,outfile,arg,serout,namelam


ntypes=8
matrixdim=ntypes+1 ! ntypes+zero col
rows=ntypes+2 ! ntypes+zero col + index
lambda=0.0



allocate(intermatrix(matrixdim,matrixdim))


do jj=1,500
    write(namefile,*)"inters3.",jj,".dat"
    call StripSpaces(namefile)
    open(unit=10,file=namefile,status="old",iostat=ios)
    write(namelam,*)"lambda.",jj,".dat"
    call StripSpaces(namelam)
    open(unit=30,file=namelam,status="old",iostat=ios)
    if(ios .eq. 0 )then 
    write(serout,*)"ave.F.",jj,".out"
    call StripSpaces(serout)
    open(unit=20,file=serout,status="replace")
    write(20,*)'# Time-averaged data for fix'
    write(20,*)'# TimeStep c_thermo_temp c_thermo_pe f_totW v_lambda'

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

    CC = intermatrix(2,6) + intermatrix(3,8) + intermatrix(4,9)  + intermatrix(5,7) + &
        &  intermatrix(2,7) +intermatrix(3,9) + intermatrix(5,6) +  intermatrix(4,8)


    DD = intermatrix(2,3) + intermatrix(2,4) + intermatrix(3,5) + intermatrix(4,5) 


    kk=kk+1
    write(20,*)kk,kk,kk,2.0*DD+CC,lambda

end do


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






end program


