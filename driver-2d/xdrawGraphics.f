c initializeGraphics
c####################################################################
      subroutine  initializeGraphics 

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use parameters

      use timeStepping

      use variables

      use graphics

      use iosetup

      use mg_setup

      implicit none

c Call variables

c Local variables

      integer  :: i,nplots,prof_ivar(20)

c Begin program

c Set graphics plotting range

      imin = 1
      imax = nxd
      jmin = 0
      jmax = nyd+1

c Set graphics dumping interval

      dfreq = 8d0

      if (tmax.gt.0d0) then
        if (dstep.eq.0d0) dstep = dt*max(int((tmax-time)/dfreq/dt),1)
        rstep = min(dt*max(int((tmax-time)/dfreq/dt),1),dstep)
cc        rstep = dstep
        ndstep = -1
      else
        if (ndstep.eq.0) ndstep = max(numtime/int(dfreq),1)
        nrstep = min(max(numtime/int(dfreq),1),nrstep)
cc        nrstep = ndstep
        dstep = 1e30
      endif

c Initialize contour graph array

      call initializeContourPlots

c Initialize profile description array

      where (array_graph%descr /= '')
        prof_desc = 'Profile ' // array_graph%descr 
      elsewhere
        prof_desc = ''
      end where

c Initialize diagnostics description array (see evaluateDiagnostics)

      call initializeTimePlots

c Determine number of contour plots

      nqty = 0
      do i = 1,size(array_graph)
        if (len(trim(array_graph(i)%descr)) == 0 ) exit
        nqty = nqty + 1
      enddo

c Open graphics files

      call openGraphicsFiles(restart)

c Initialize contour plots

      if (.not.restart) then
        call contour_step(xl(imin,ngrd),yl(jmin,ngrd),ucontour,0)
      endif
    
c Create draw*.in files

c     Time histories

      nplots = count (sel_diag /= 0)

      call createDrawInGfile(nplots,lineplotfile,'Time histories'
     .           ,'time',sel_diag,diag_desc,diag_ivar,'drawgamma.in')

c     Contours

      nplots = count (sel_cont /= 0)

      call createDrawInMfile(nplots,contourfile,'Contour Plots'
     .        ,'time','x','y',sel_cont,array_graph%descr
     .        ,'-L50','drawcontours.in')

c     Profiles (selects same profiles as contour plots except vector plots)

      prof_ivar(:) = 0 !All profiles have same (default) independent variable: y

      call createDrawInGfile(nplots,profilefile,'Profiles','y'
     .          ,abs(sel_cont),prof_desc,prof_ivar,'drawprofiles.in')

c End program

      return
      end

c dumpTimeStepPlots
c####################################################################
      subroutine dumpTimeStepPlots(nplot,tmplot,itime)

c--------------------------------------------------------------------
c     Dumps time plots
c--------------------------------------------------------------------

      use timeStepping

      use graphics

      use iosetup

      use mg_setup

      implicit none

c Call variables

      integer*4   nplot,itime
      real*8      tmplot

c Local variables

      integer*4   i,j,ieq

c Begin program

      if (nplot.eq.ndstep.or.tmplot.ge.dstep.or.itime.eq.0) then

        nplot  = 0
        if (itime.gt.0) tmplot = tmplot - dstep

        call profplot(uprofile)

        call contour_step(xl(imin,ngrd),yl(jmin,ngrd),ucontour,1)

        call line_diagnostics(itime)

      endif

      if (itime > 0) write(ulineplot) real(time),real(diagnostics)

c End program

      return
      end

c profplot
c####################################################################
      subroutine profplot(unit)

c--------------------------------------------------------------------
c     Does time plots of averaged profiles
c--------------------------------------------------------------------

      use parameters

      use mg_setup

      use graphics

      implicit none

c Call variables

      integer*4   unit

c Local variables

      integer*4   j,igr

c Begin program

c Average in periodic direction and dump plots

      do j=0,nyd+1
        do igr =1,nqty
          profiles(igr) = sum(array_graph(igr)%array(1:nxd-1,j),1)
     .                   /(nxd-1)
        enddo
        write(unit) real(yl(j,ngrd)),real(profiles(1:nqty))
      enddo
      write(unit)

c End

      return
      end

c contour
c #########################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
      implicit none               !For safe fortran
c -------------------------------------------------------------------------
c     Contours arr in xdraw format
c     Notes:
c      put the next 2 lines in main
c      open(unit=nunit,file='contour.bin',form='unformatted') before
c      close(unit=nunit) after
c -------------------------------------------------------------------------

c Call variables

      integer*4      nx,ny,iopt,nunit
      real*8         arr(0:nx+1,0:ny+1),xmin,xmax,ymin,ymax

c Local variables

      integer*4      i,j

c Begin program

      if(iopt.eq.0) then
        write(nunit) nx-1,ny-1,0
        write(nunit) real(xmin,4),real(xmax,4),real(ymin,4),real(ymax,4) 
      endif
      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)

      return
      end

c contour_step
c #########################################################################
      subroutine contour_step(x,y,nunit,iopt)

c -------------------------------------------------------------------------
c     Contours arr in xdraw format
c     Notes:
c      put the next 2 lines in main
c      open(unit=nunit,file='contour.bin',form='unformatted') before
c      close(unit=nunit) after
c -------------------------------------------------------------------------

      use graphics

      implicit none

c Call variables

      integer*4 :: iopt,nunit
      real*8    :: x(imin:imax),y(jmin:jmax)

c Local variables

      integer*4      nsblk,nusblk
      integer*4      i,j,igr

c Begin program

      nsblk  = 1          ! Number of structured blocks
      nusblk = 0          ! Number of unestructured blocks

      if(iopt.eq.0) then

c     Write initialization headings

        write(nunit) nsblk,nusblk,nqty
        write(nunit) (imax-imin),(jmax-jmin)
        write(nunit) ((real(x(i),4),i=imin,imax),j=jmin,jmax),
     $               ((real(y(j),4),i=imin,imax),j=jmin,jmax)

      else

c     Write arrays

        do igr=1,nqty
          write(nunit) ((real(array_graph(igr)%array(i,j),4)
     .                  ,i=imin,imax),j=jmin,jmax)
        enddo

      endif

      return
      end

c finalizeGraphics
c####################################################################
      subroutine  finalizeGraphics 

c--------------------------------------------------------------------
c     Set graphics files and dumping intervals
c--------------------------------------------------------------------

      use iosetup

      implicit none

c Call variables

c Local variables

c Begin program

cc      call openGraphicsFiles(.true.)

c Finalize time traces output (separator)

      write(ulineplot)

c Close all files

      call closeGraphicsFiles

c End program

      return
      end

*deck createDrawInGfile
c#####################################################################
      subroutine createDrawInGfile(nvar,binf,graphtitle,ivar
     .                       ,dvariables,descr,ivariables,drawfile)

c---------------------------------------------------------------------
c     Creates the G-type (x-y curves) draw.in file for xdraw graphics.
c     On call sequence:
c       * nvar : number of variables in consideration for plotting
c       * binf : name of binary file to be read
c       * graphtitle : character array of arbitrary length
c       * ivar : independent variable axis name
c       * dvariables : integer array of size (nvar) identifying the 
c                      dependent variables (although they can also be 
c                      independent variables)
c       * descr      : character array containing dependent variables
c                      descriptions.
c       * ivariables : integer array specifying independent variables 
c                      for all the dependent variables.
c                      These include any of the dependent
c                      variable indeces plus 0 (indicating the true
c                      independent variable).
c       * drawfile : character variable specifying the name of the
c                    draw.in file.
c---------------------------------------------------------------------

      use iosetup

      implicit none

c Call variables

      integer     :: nvar,ivariables(20),dvariables(nvar)
      character(*):: binf,drawfile,graphtitle,descr(20),ivar

c Local variables

      integer      :: i,j
      character(50):: title
      logical      :: flag

c Begin program

      open (unit=mfile,file=drawfile,status='unknown')

      write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      G'
      write (mfile,'(a)') 'filename(s)'
      write (mfile,'(a,/)') trim(binf)

      write (mfile,'(2a,/)') 'graph title: ',trim(graphtitle)
      write (mfile,'(a)') 'variable names'
      write (mfile,'(2a)') ' 0       ',trim(ivar)
      do i=1,size(descr)
        if (len(trim(descr(i))) == 0) exit
        write (mfile,'(i2,2a)') i  ,'       ',trim(descr(i))
      enddo
      write (mfile,'(/,a)') 'ix       iy        window title'

      do i=1,size(descr)
        if (len(trim(descr(i))) == 0) exit
        flag = .false.
        do j = 1,nvar
          if (dvariables(j).eq.i) then
            flag = .true.
            exit
          endif
        enddo
        if (flag) then
          write (mfile,'(i2,a,i2,2a)')       ivariables(i),'      '
     .                                ,i,'         ',trim(descr(i))
        else
          write (mfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
     .                                ,i,'         ',trim(descr(i))
        endif
      enddo

      write (mfile,*)

      close (mfile)

c End program

      return
      end

*deck createDrawInMfile
c#####################################################################
      subroutine createDrawInMfile(nvar,binf,graphtitle,ivart,ivarx
     .                       ,ivary,variables,descr,options,drawfile)

c---------------------------------------------------------------------
c     Creates drawfile for xdraw graphics.
c---------------------------------------------------------------------

      use iosetup

      implicit none

c Call variables

      integer     :: nvar,variables(nvar)
      character(*):: binf,drawfile,graphtitle,options,descr(20)
      character(*):: ivart,ivarx,ivary

c Local variables

      integer      :: i,j
      character(50):: title
      logical      :: flag

c Begin program

      open (unit=mfile,file=drawfile,status='unknown')

      write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      M'
      write (mfile,'(a)') 'filename(s)'
      write (mfile,'(a,/)') trim(binf)


      write (mfile,'(2a,/)') 'comment:      ',graphtitle
      write (mfile,'(a)') 'independent variable names'
      write (mfile,'(2a)')   't       ',trim(ivart)
      write (mfile,'(2a)')   'x       ',trim(ivarx)
      write (mfile,'(2a,/)') 'y       ',trim(ivary)
      write (mfile,'(a)') 'dependent variable names'
      do i=1,size(descr)
        if (len(trim(descr(i))) == 0) exit
        write (mfile,'(i2,2a)') i  ,'       ',descr(i)
      enddo
      write (mfile,'(/,a)') 'iqty       options     window title'

      do i=1,size(variables)
        if (variables(i) > 0) then
          write (mfile,'(i2,4a)') variables(i)
     .         ,'         ',options
     .         ,'         ',trim(descr(variables(i)))
        elseif (variables(i) == 0) then
          cycle
        else
          variables(i) = - variables(i)
          title = 'Vector plot (' // trim(descr(variables(i)))
     .              // ',' // trim(descr(variables(i)+1))//')'
          write (mfile,'(i2,a,i2,4a)') variables(i),' ',variables(i)+1
     .         ,'       ',options,'         ',title
        endif
      enddo
      write (mfile,*)

      close (mfile)

c End program

      return
      end
