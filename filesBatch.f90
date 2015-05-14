module filesBatch
implicit none
contains

    subroutine list2file(filePath,outFilename)
    use dflib
    implicit none
    character(*), intent(in) :: filePath, outFilename
    character(len=255) :: cmd
    logical :: res
    
    ! type "dir /? " in command window to learn what's the mean of the following cmd
    cmd='dir /a-d/b/s '//trim(filePath)//'>'//trim(outFilename)
    !cmd='dir /a-d/b '//trim(filePath)//'>'//trim(outFilename)
    res = systemqq(cmd)
    end subroutine list2file
!BL
    subroutine getFileLines(filename,totalLines)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: totalLines

    ! input filename's open-status
    integer :: inStat  
    
    ! read-file status
    integer :: readStat

    totalLines = 0

    open(7,file=Filename,iostat=inStat,status='old',action='read')
    in_ok: if(inStat/=0)then
        write(*,*) 'Input file OPEN failed: iostat=', inStat
        write(*,*) 'Shutting Down ... ... '
    else    
        ! Just read nothing to count lines in the file
        do 
            read(7,*,iostat=readStat) 
            totalLines = totalLines + 1
            read_ok: if(readStat < 0) then
                !write(*,*)'End of file, iostat=',readStat
                !write(*,*)'<'//trim(filename)//'>',' file has total lines:',totalLines-1
                totalLines = totalLines-1
                close(7)
                return 
            elseif(readStat > 0)then  ! read bad formating data 
                !write(*,*)'Intel Fortran run-time library (forRTL) error number:',readStat
                !write(*,*)'>>>Check maybe bad record in file:', inFilename
            else ! read successfully
                !write(*,*)startDate, startTime, &
                !      stationID, pointID, diInster, fInster, pd, &
                !      dBaseline, hBaseline, zBaseline, fBaseline
            endif read_ok
        enddo
    endif in_ok
    close(7)                 
    end subroutine getFileLines
end module filesBatch
