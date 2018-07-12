      PROGRAM conv_F90
          character (len=512) :: line1,line2
          integer :: eof

          logical :: debug = .FALSE.


          CALL get_line(line1,eof,5)
          IF (debug) write(6,*) 'Old: ',trim(line1)
          DO WHILE(eof == 0)
            CALL get_line(line2,eof,5)
            IF (debug) write(6,*) 'Old: ',trim(line2)
            IF (line2(6:6) /= ' ' .AND. line2(1:1) == ' ') THEN
              line1(73:73) = '&'
              line2(6:6) = ' '
            END IF
            IF (line1(1:1) /= ' ') line1(1:1) = '!'
            IF (trim(adjustl(line1)) /= '&') write(6,'(a)') trim(line1)
            line1 = line2
          END DO
          IF (line1(1:1) /= ' ') line1(1:1) = '!'
          IF (trim(adjustl(line1)) /= '&') write(6,'(a)') trim(line1)
          
          END PROGRAM conv_F90
          SUBROUTINE get_line(line,eof,nio)
            character (len=512) :: line
            character :: c
            integer :: nio,eof
            integer :: i

            i=0
            eof = 0
            line(1:512) = ' '
            DO
              i = i + 1
              read(nio,'(a)',advance='no',eor=999,end=998) line(i:i)
            END DO
  998       eof = 1
  999       CONTINUE
            IF (i >= 511) 
     *           write(6,*) 'WARINNIG, the line can be incomplete'
          END SUBROUTINE get_line
