      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
C***BEGIN PROLOGUE  XERRWD
C***SUBSIDIARY
C***PURPOSE  Write error message with values.
C***LIBRARY   MATHLIB
C***CATEGORY  R3C
C***TYPE      REAL*8 (XERRWV-S, XERRWD-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
C  as given here, constitute a simplified version of the SLATEC error
C  handling package.
C
C  All arguments are input arguments.
C
C  MSG    = The message (character array).
C  NMES   = The length of MSG (number of characters).
C  NERR   = The error number (not used).
C  LEVEL  = The error level..
C           0 or 1 means recoverable (control returns to caller).
C           2 means fatal (run is aborted--see note below).
C  NI     = Number of integers (0, 1, or 2) to be printed with message.
C  I1,I2  = Integers to be printed, depending on NI.
C  NR     = Number of reals (0, 1, or 2) to be printed with message.
C  R1,R2  = Reals to be printed, depending on NR.
C
C  Note..  this routine is machine-dependent and specialized for use
C  in limited context, in the following ways..
C  1. The argument MSG is assumed to be of type CHARACTER, and
C     the message is printed with a format of (1X,A).
C  2. The message is assumed to take only one line.
C     Multi-line messages are generated by repeated calls.
C  3. If LEVEL = 2, control passes to the statement   STOP
C     to abort the run.  This statement may be machine-dependent.
C  4. R1 and R2 are assumed to be in real*8 and are printed
C     in D21.13 format.
C
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   920831  DATE WRITTEN
C   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
C   930329  Modified prologue to SLATEC format. (FNF)
C   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
C   930922  Minor cosmetic change. (FNF)
C***END PROLOGUE  XERRWD
C
C*Internal Notes:
C
C For a different default logical unit number, IXSAV (or a subsidiary
C routine that it calls) will need to be modified.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERRWD.. None
C Function routine called by XERRWD.. IXSAV
C-----------------------------------------------------------------------
C**End
C
C  Declare arguments.
C
      REAL*8 R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG
C
C  Declare local variables.
C
      INTEGER LUNIT, IXSAV, MESFLG
C
C  Get logical unit number and message print flag.
C
C***FIRST EXECUTABLE STATEMENT  XERRWD
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C
C  Write the message.
C
      WRITE (LUNIT,10)  MSG
 10   FORMAT(1X,A)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',D21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
C
C  Abort the run if LEVEL = 2.
C
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERRWD ----------------------
      END
