    program cosubs
! COrrelated Sampling Using Batch Statistics, using MCNP TFCs.
! Jeffrey A. Favorite, Los Alamos National Laboratory, XTD-RTA.
! An early Python version for testing was written by Michael Squire, U. Texas-San Antonio.
!
! Copyright 2023. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by Triad National
! Security, LLC for the U.S. Department of Energy/National Nuclear Security
! Administration. All rights in the program are reserved by Triad National Security,
! LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive,
! paid-up, irrevocable worldwide license in this material to reproduce, prepare
! derivative works, distribute copies to the public, perform publicly and display
! publicly, and to permit others to do so.
!
    implicit none
! compile with
!   ifort -o cosubs -r8 -i8 -check bounds -warn unused -traceback cosubs.F90
    integer iui,iuo
    parameter (iui=2,iuo=3)
! used to store information returned from the calc_* routines and possibly printed.
    integer maxinfo
    parameter (maxinfo=9)
! rtmp is for testing.
    real*8 diff(maxinfo),sens(maxinfo),c1,c2,relpert,rn0,rtimes(3) !,rtmp
    integer ntal,npert,nentry,nbatch,nfiles,nhash,ifmt,italratio,iratio,iunc, &
     ierr,ierr1,ierr2,npsf,iprod,nba,nthb,i,n0,it,itn,ip,ip1,npsb,nbatch1
    real*8,  allocatable, dimension(:) :: fbase1,fpos1,fneg1,fbase2,fpos2,fneg2
    real*8,  allocatable, dimension(:,:,:,:) :: tally_base,tally_pos,tally_neg
    integer, allocatable, dimension(:) :: nps,ital,ipert
    character bfile*120,pfile*120,nfile*120,ofile*120
    logical lsens,lpert,lpfile,lnfile,lofile
!
! iprod=0/1 not production version/production version
!  (a logical variable is not used to allow different levels later on if desired).
    iprod=1
!
! lpfile=.true. if there is a positive-pert. file
! lnfile=.true. if there is a negative-pert. file
! lofile=.true. if there is an output file
    call rd_input_line(iuo,bfile,pfile,nfile,ofile,lpfile,lnfile,lofile,nhash,relpert,italratio,iunc)
!
! nps=total nps up to batch as read from TFC
! npsf=total nps in problem
! ntal=number of tallies in problem
! npert=number of perturbations in problem
! nfiles=number of files to read (base files)
! nhash=number of # signs in base file name
! nentry=number of entries in TFC
! nbatch=number of batches to use in uncertainty
! nbatch1=nbatch+1
! ifmt=format of TFCs
! relpert=relative perturbation
! italratio=tally to use as denominator in ratios
    ierr=0 ; ierr1=0 ; ierr2=0
    call pass1TFC(iui,bfile,ntal,npert,nfiles,nhash,nentry,ifmt,ierr)
    if(ierr > 0)then
      if(lofile)close(iuo)
      stop
    end if
    if(lpfile)then
      call checkTFC(iui,pfile,bfile,ntal,npert,nfiles,nhash,nentry,ifmt,ierr1)
    end if
    if(lnfile)then
      call checkTFC(iui,nfile,bfile,ntal,npert,nfiles,nhash,nentry,ifmt,ierr2)
    end if
! maybe if method=-1, -2, or -3 there would be a reason to compute differences and
! sensitivities with perturbed tallies.
    if(npert > 0.and.(lpfile.or.lnfile))then
      write(*,'("warning. perturbations are unused because there is a -p and/or -n file.")')
      npert=0
    end if
    write(*,'("number of files: ",i0)')nfiles
    write(*,'("number of tallies: ",i0)')ntal
    write(*,'("number of perturbations: ",i0)')npert
    if(lofile)then
      write(iuo,'("number of files: ",i0)')nfiles
      write(iuo,'("number of tallies: ",i0)')ntal
      write(iuo,'("number of perturbations: ",i0)')npert
    end if
    if(ierr1 > 0.or.ierr2 > 0)then
      if(lofile)close(iuo)
      stop
    end if
!
! c1 is the constant for relative sensitivity using one-sided diff. or pert.
! c2 is the constant for relative sensitivity using central diff.
! lsens=.true. if relative sensitivity is calculated (otherwise difference or ratio only)
! lpert=.true. if there is one file with a pert card
! iratio=0/1/2 no ratio/ratio only/diff and sens
    if(italratio > 0)then
      if(.not.lpfile.and..not.lnfile)then
        iratio=1
      else
        iratio=2
      end if
    else
      iratio=0
    end if
    if(.not.lpfile.and..not.lnfile.and.iratio == 0)then
      lpert=.true.
      if(npert == 0)then
        write(*,'("error. base file only but no tally specified and no perturbations.")')
        ierr=1
      end if
      if(relpert == 0.d0)then
        write(*,'("error. base file only but no tally or relative perturbation specified.")')
        ierr=2
      else
        lsens=.true.
        c1=1.d0/relpert
      end if
    else
      lpert=.false.
      if(relpert /= 0.d0)then
        lsens=.true.
        c1=1.d0/relpert
        c2=0.5d0/relpert
      else
        lsens=.false.
        c2=0.d0
        c1=0.d0
      end if
    end if
!
    if(ierr > 0)then
      if(lofile)close(iuo)
      stop
    end if
!
    nbatch=nfiles*nentry
    nbatch1=nbatch+1
    write(*,'("number of batches: ",i0)')nbatch
    if(lofile)write(iuo,'("number of batches: ",i0)')nbatch
! allocate arrays.
! tally(1,...) tally value
! tally(2,...) relative uncertainty
    allocate(tally_base(1:2,1:ntal,0:npert,1:nbatch1))
    allocate(tally_pos(1:2,1:ntal,0:npert,1:nbatch1))
    allocate(tally_neg(1:2,1:ntal,0:npert,1:nbatch1))
    allocate(nps(0:nbatch))
    allocate(ital(1:ntal))
    allocate(ipert(1:npert))
    allocate(fbase1(1:nbatch))
    allocate(fpos1(1:nbatch))
    allocate(fneg1(1:nbatch))
    if(iratio == 2)then
      allocate(fbase2(1:nbatch))
      allocate(fpos2(1:nbatch))
      allocate(fneg2(1:nbatch))
    end if
!   call allocate_arrays1(nbatch,ntal,npert, &
!    nps,ital,ipert, &
!    fbase1,fpos1,fneg1, &
!    tally_base,tally_pos,tally_neg)
!
! get tallies.
    nps(0:nbatch)=0
    ital(1:ntal)=0
    ipert(1:npert)=0
    rtimes(1:3)=0.d0
    tally_base(1:2,1:ntal,0:npert,1:nbatch1)=0.d0
    tally_pos(1:2,1:ntal,0:npert,1:nbatch1)=0.d0
    tally_neg(1:2,1:ntal,0:npert,1:nbatch1)=0.d0
    call readTFC(iui,iuo,lofile,bfile,bfile,ifmt,ntal,npert,nentry,nbatch,nps,nfiles,nhash,ital,ipert,tally_base,rtimes(1))
    if(lpfile)then
      call readTFC(iui,iuo,lofile,pfile,bfile,ifmt,ntal,npert,nentry,nbatch,nps,nfiles,nhash,ital,ipert,tally_pos,rtimes(2))
    end if
    if(lnfile)then
      call readTFC(iui,iuo,lofile,nfile,bfile,ifmt,ntal,npert,nentry,nbatch,nps,nfiles,nhash,ital,ipert,tally_neg,rtimes(3))
    end if
!   do i=1,nbatch
!     write(*,'(i12,1p4e18.10)')nps(i),tally_base(1:2,1,0,i),tally_pos(1:2,1,0,i)
!   end do ! i
!
! check italratio. it changes from the problem tally name to the tally index here.
! itn is the tally name to pass to write_head2.
    itn=0
    if(iratio > 0)then
      do i=1,ntal
        if(ital(i) == italratio)then
          italratio=i
          itn=ital(italratio)
          exit
        end if
      end do ! i
      if(itn == 0)then
        write(*,'("error. tally ",i0," does not exist.")')italratio
        if(lofile)close(iuo)
        stop
      end if
    end if
! check nps intervals
    npsf=0
    npsb=nps(1)
    do i=1,nbatch
      if(nps(i)-nps(i-1) > 0 .and. nps(i)-nps(i-1) /= npsb)then
        write(*,'("error. nps intervals are not equal.",/, &
         "  batch ",i0," present interval ",i0, &
         "; previous interval ",i0,".")')i,nps(i)-nps(i-1),npsb
        ierr=1
      else if(nps(i)-nps(i-1) < 0 .and. nps(i) /= npsb)then
        write(*,'("error. nps intervals are not equal.",/, &
         "  batch ",i0," present interval ",i0, &
         "; previous interval ",i0,".")')i,nps(i),npsb
        ierr=1
      end if
      if(nps(i)-nps(i-1) > 0)then
        npsb=nps(i)-nps(i-1)
      else
        npsb=nps(i)
      end if
      if(mod(i,nentry) == 0)then
        npsf=npsf+nps(i)
      end if
    end do ! i
    write(*,'("number of histories: ",i0)')npsf
    if(lofile)write(iuo,'("number of histories: ",i0)')npsf
    if(iunc == 0)then
      write(*,'("uncertainties are absolute, not relative.")')
      if(lofile)write(iuo,'("uncertainties are absolute, not relative.")')
    else if(iunc == 1)then
      write(*,'("uncertainties are relative.")')
      if(lofile)write(iuo,'("uncertainties are relative.")')
    else if(iunc == 2)then
      write(*,'("uncertainties are relative, in percent.")')
      if(lofile)write(iuo,'("uncertainties are relative, in percent.")')
    end if
!
! ip is the perturbation number for most tallies. all perturbations are read,
! but only ip=1 is used, and only if the base file is the only one specified.
    ip=0
    do it=1,ntal
      if(iratio > 0.and.it == italratio)then
        if(iprod == 1)cycle ! don't do italratio/italratio
      end if
      do ip1=0,npert
        if(lpert.and.ip1 == 0)then
          cycle
        end if
! 6 is the unit number for standard out
        if(lofile)then
          call write_head2(iuo,iratio,lpert,lsens,ital(it),itn,ipert(ip1))
        else
          call write_head2(6,iratio,lpert,lsens,ital(it),itn,ipert(ip1))
        end if
! cycle through batches.
! 10 is the smallest number of batches allowed.
        do nba=10,nbatch
! this is for testing.
!       do nba=nbatch,nbatch
!         write(*,'("nbatch,nba=",2i6)')nbatch,nba
          if(mod(nbatch,nba) /= 0)cycle
! every nth batch.
          nthb=nbatch/nba
! use this until I figure out how to mix tallies from different runs.
!         write(*,'("nentry,nthb=",2i6)')nentry,nthb
          if(mod(nentry,nthb) /= 0)cycle
          n0=0
          do i=nthb,nbatch,nthb
            n0=n0+1
! why don't I use "if(i == nthb)then" or "if(n0 == 1)then" here?
            if(nps(i) == nps(nthb))then
!             write(*,'("A n0,i,nps=",2i6,3i12)')n0,i,nps(i)
              if(iratio == 0)then
                fbase1(n0)=tally_base(1,it,ip,i)
              else if(iratio == 1)then
                fbase1(n0)=tally_base(1,italratio,ip,i)
              else if(iratio == 2)then
                fbase1(n0)=tally_base(1,it,ip,i)
                fbase2(n0)=tally_base(1,italratio,ip,i)
                fpos2(n0)=tally_pos(1,italratio,ip,i)
                fneg2(n0)=tally_neg(1,italratio,ip,i)
              end if
              if(lpert)then
                fpos1(n0)=tally_base(1,it,ip1,i)
              else if(iratio == 1)then
                fpos1(n0)=tally_base(1,it,ip,i)
              else
                fpos1(n0)=tally_pos(1,it,ip,i)
              end if
              fneg1(n0)=tally_neg(1,it,ip,i)
            else
!             write(*,'("B n0,i,nps=",2i6,3i12)')n0,i,nps(i)-nps(i-nthb) !,nps(nthb)
              rn0=real(nps(i)-nps(i-nthb),8)
              if(iratio == 0)then
                fbase1(n0)=(tally_base(1,it,ip,i)*real(nps(i),8) &
                 -tally_base(1,it,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
              else if(iratio == 1)then
                fbase1(n0)=(tally_base(1,italratio,ip,i)*real(nps(i),8) &
                 -tally_base(1,italratio,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
              else if(iratio == 2)then
                fbase1(n0)=(tally_base(1,it,ip,i)*real(nps(i),8) &
                 -tally_base(1,it,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
                fbase2(n0)=(tally_base(1,italratio,ip,i)*real(nps(i),8) &
                 -tally_base(1,italratio,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
                fpos2(n0)=(tally_pos(1,italratio,ip,i)*real(nps(i),8) &
                 -tally_pos(1,italratio,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
                fneg2(n0)=(tally_neg(1,italratio,ip,i)*real(nps(i),8) &
                 -tally_neg(1,italratio,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
              end if
              if(lpert)then
                fpos1(n0)=(tally_base(1,it,ip1,i)*real(nps(i),8) &
                 -tally_base(1,it,ip1,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
              else if(iratio == 1)then
                fpos1(n0)=(tally_base(1,it,ip,i)*real(nps(i),8) &
                 -tally_base(1,it,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
              else
                fpos1(n0)=(tally_pos(1,it,ip,i)*real(nps(i),8) &
                 -tally_pos(1,it,ip,i-nthb)*real(nps(i-nthb),8)) &
                 /rn0
              end if
              fneg1(n0)=(tally_neg(1,it,ip,i)*real(nps(i),8) &
               -tally_neg(1,it,ip,i-nthb)*real(nps(i-nthb),8)) &
               /rn0
            end if
! this is for testing.
!           rtmp=fbase1(n0)*fpos1(n0)
!           if(rtmp < 0.5d0)rtmp=0.d0
!           write(*,'(2i12,1p3e14.6)')n0,i,fbase1(n0),fpos1(n0),rtmp
          end do ! i
          diff(1:maxinfo)=0.d0
          sens(1:maxinfo)=0.d0
! Index  Information
!   1      Difference or relative sensitivity
!   2      Uncertainty using batch statistics with covariances
!   3      Uncertainty using batch statistics without covariances
!   4      Uncertainty using the last entries in the TFC (no covariances)
!   5      Theoretical uncertainty including correlations
!   6      Theoretical uncertainty without correlations
! all uncertainties are absolute, not relative.
          if(lpert)then
! diff is first group by itself
! sens is first group divided by second group
            call calc_pert(n0,maxinfo, &
             tally_base(1:2,it,ip1,nbatch1),fpos1(1:n0), &
             tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
             c1,diff,sens)
            if(ip1 == 1)then
              sens(9)=sqrt(sens(4)**2 &
               -(tally_base(1,it,2,nbatch1)**2*tally_base(2,it,2,nbatch1)**2 &
               -tally_base(1,it,ip1,nbatch1)**2*tally_base(2,it,ip1,nbatch1)**2 &
               -tally_base(1,it,ip,nbatch1)**2*tally_base(2,it,ip,nbatch1)**2) &
               *sens(1)**2/tally_base(1,it,ip1,nbatch1)/tally_base(1,it,ip,nbatch1))
            end if
          else if(iratio == 1)then
! ratio is first group divided by second group
            call calc_ratio(n0,maxinfo, &
             tally_base(1:2,it,ip,nbatch1),fpos1(1:n0), &
             tally_base(1:2,italratio,ip,nbatch1),fbase1(1:n0), &
             diff,sens)
          else if(iratio == 0)then
            if(.not.lnfile)then
! diff is first group minus second group
              call calc_diff(n0,maxinfo, &
               tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
               tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
               diff)
! sens is (first group minus second group) divided by second group
              if(lsens)then
                call calc_sens1(n0,maxinfo,1, &
                 tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
                 tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
                 c1,sens)
              end if
            else if(.not.lpfile)then
! diff is first group minus second group
              call calc_diff(n0,maxinfo, &
               tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
               tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
               diff)
! sens is (second group minus first group) divided by second group
              if(lsens)then
                call calc_sens1(n0,maxinfo,2, &
                 tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
                 tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
                 c1,sens)
              end if
            else
! diff is first group minus second group
              call calc_diff(n0,maxinfo, &
               tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
               tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
               diff)
! sens is (first group minus second group) divided by third group
              if(lsens)then
                call calc_sens2(n0,maxinfo, &
                 tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
                 tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
                 tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
                 c2,sens)
              end if
            end if
          else if(iratio == 2)then
            if(.not.lnfile)then
! diff is (first group divided by second group) minus (third group divided by fourth group)
              call calc_diff_ratio(n0,maxinfo, &
               tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
               tally_pos(1:2,italratio,ip,nbatch1),fpos2(1:n0), &
               tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
               tally_base(1:2,italratio,ip,nbatch1),fbase2(1:n0), &
               diff)
! sens is [(first group divided by second group) minus (third group divided by fourth group)]
!  divided by (third group divided by fourthgroup)
              if(lsens)then
                call calc_sens1_ratio(n0,maxinfo,1, &
                 tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
                 tally_pos(1:2,italratio,ip,nbatch1),fpos2(1:n0), &
                 tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
                 tally_base(1:2,italratio,ip,nbatch1),fbase2(1:n0), &
                 c1,sens)
              end if
            else if(.not.lpfile)then
! diff is (first group divided by second group) minus (third group divided by fourth group)
              call calc_diff_ratio(n0,maxinfo, &
               tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
               tally_neg(1:2,italratio,ip,nbatch1),fneg2(1:n0), &
               tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
               tally_base(1:2,italratio,ip,nbatch1),fbase2(1:n0), &
               diff)
! sens is [(third group divided by fourth group) minus (first group divided by second group)]
!  divided by (third group divided by fourth group)
              if(lsens)then
                call calc_sens1_ratio(n0,maxinfo,2, &
                 tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
                 tally_neg(1:2,italratio,ip,nbatch1),fneg2(1:n0), &
                 tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
                 tally_base(1:2,italratio,ip,nbatch1),fbase2(1:n0), &
                 c1,sens)
              end if
            else
! diff is (first group divided by second group) minus (third group divided by fourth group)
              call calc_diff_ratio(n0,maxinfo, &
               tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
               tally_pos(1:2,italratio,ip,nbatch1),fpos2(1:n0), &
               tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
               tally_neg(1:2,italratio,ip,nbatch1),fneg2(1:n0), &
               diff)
! sens is [(first group divided by second group) minus (third group divided by fourth group)]
!  divided by (fifth group divided by sixth group)
              if(lsens)then
                call calc_sens2_ratio(n0,maxinfo, &
                 tally_pos(1:2,it,ip,nbatch1),fpos1(1:n0), &
                 tally_pos(1:2,italratio,ip,nbatch1),fpos2(1:n0), &
                 tally_neg(1:2,it,ip,nbatch1),fneg1(1:n0), &
                 tally_neg(1:2,italratio,ip,nbatch1),fneg2(1:n0), &
                 tally_base(1:2,it,ip,nbatch1),fbase1(1:n0), &
                 tally_base(1:2,italratio,ip,nbatch1),fbase2(1:n0), &
                 c2,sens)
              end if
            end if
          end if
! output
! iunc=0/1/2 absolute/relative/relative in percent
          if(iunc > 0)then
            diff(2:maxinfo)=diff(2:maxinfo)/abs(diff(1))
            sens(2:maxinfo)=sens(2:maxinfo)/abs(sens(1))
            if(iunc == 2)then
              diff(2:maxinfo)=diff(2:maxinfo)*100.d0
              sens(2:maxinfo)=sens(2:maxinfo)*100.d0
            end if
          end if
          if(lofile)then
            write(iuo,'(i8,1p4e14.6)',advance="no")n0,diff(1:4)
            if(lsens)then
              write(iuo,'(1p4e14.6)')sens(1:4)
            else
              write(iuo,'("")')
            end if
          else
            write(*,'(i8,1p4e14.6)',advance="no")n0,diff(1:4)
            if(lsens)then
              write(*,'(1p4e14.6)')sens(1:4)
            else
              write(*,'("")')
            end if
          end if
        end do ! nba
        if(lpert)then
          if(lofile)then
            write(iuo,'(2x,"warning. method=2 is assumed.")')
          else
            write(*,'(2x,"warning. method=2 is assumed.")')
          end if
        end if
! print theoretical uncertainties.
        if(iprod == 0.and..not.lpert)then
          diff(5:6)=diff(5:6)/sqrt(real((npsf-1),8))
          sens(5:6)=sens(5:6)/sqrt(real((npsf-1),8))
          if(lofile)then
            write(iuo,'(2x,"theoretical,corr.",1pe17.6)',advance="no")diff(5)
!           write(iuo,'(2x,"theoretical,corr.",1pe17.6)',advance="no")diff(7)
            if(sens(5) /= 0.d0)then
              write(iuo,'(42x,1pe14.6)')sens(5)
            else
              write(iuo,'("")')
            end if
            write(iuo,'(2x,"theoretical,uncorr.",1pe15.6)',advance="no")diff(6)
!           write(iuo,'(2x,"theoretical,uncorr.",1pe15.6)',advance="no")diff(8)
            if(sens(6) /= 0.d0)then
              write(iuo,'(42x,1pe14.6)')sens(6)
            else
              write(iuo,'("")')
            end if
          else
            write(*,'(2x,"theoretical,corr.",1pe17.6)',advance="no")diff(5)
!           write(*,'(2x,"theoretical,corr.",1pe17.6)',advance="no")diff(7)
            if(sens(5) /= 0.d0)then
              write(*,'(42x,1pe14.6)')sens(5)
            else
              write(*,'("")')
            end if
            write(*,'(2x,"theoretical,uncorr.",1pe15.6)',advance="no")diff(6)
!           write(*,'(2x,"theoretical,uncorr.",1pe15.6)',advance="no")diff(8)
            if(sens(6) /= 0.d0)then
              write(*,'(42x,1pe14.6)')sens(6)
            else
              write(*,'("")')
            end if
          end if
        else if(iprod == 0.and.lpert.and.ip1 == 1)then
          if(lofile)then
            write(iuo,'(2x,"exact",71x,1pe14.6)')sens(9)
          else
            write(*,'(2x,"exact",71x,1pe14.6)')sens(9)
          end if
        end if
      end do ! ip1
    end do ! it
! these integer fields will never be filled, not even in seconds (1 year = 3.1536E+07 s)
! if nfiles > 0 and the average time is desired, write with "1pe17.6".
    if(lofile)then
      write(iuo,'(/,2x,"time for base",i21," s")')nint(rtimes(1),8)
      if(lpfile)write(iuo,'(2x,"time for positive",i17," s")')nint(rtimes(2),8)
      if(lnfile)write(iuo,'(2x,"time for negative",i17," s")')nint(rtimes(3),8)
    else
      write(*,'(/,2x,"time for base",i21," s")')nint(rtimes(1),8)
      if(lpfile)write(*,'(2x,"time for positive",i17," s")')nint(rtimes(2),8)
      if(lnfile)write(*,'(2x,"time for negative",i17," s")')nint(rtimes(3),8)
    end if
!
    if(lofile)then
      close(iuo)
      write(*,'(/,"see output file ",a," for results.")')trim(ofile)
    end if
    end program

    subroutine rd_input_line(iuo,bfile,pfile,nfile,ofile,lpfile,lnfile,lofile,nhash,relpert,italratio,iunc)
! read command-line arguments. some of this is from
! https://cyber.dabamos.de/programming/modernfortran/command-line-arguments.html
    implicit none
    character version*8
    parameter (version="20240206")
    integer iuo,nhash,italratio,iunc
    character bfile*120,pfile*120,nfile*120,ofile*120
    logical lpfile,lnfile,lofile
    real*8 relpert
    integer i,icycle,ierr,ios,ii
    character arg*120,cunc*120,line*140
    character(len=:), allocatable :: command
    integer :: cmdlen
    logical lrelpert,lexist
!
    bfile=""
    pfile="<none>" ; nfile="<none>" ; ofile="<none>"
    lpfile=.false. ; lnfile=.false. ; lofile=.false.
    relpert=0.d0
    nhash=0
    italratio=0
    lrelpert=.false.
    iunc=0
    cunc=""
!
    icycle=0
    ierr=0
    do i=1,command_argument_count()
      if(icycle == 1)then
        icycle=0
        cycle
      end if
      call get_command_argument(i, arg)
      select case (arg)
        case ('-v')
          write(*,'("version ",a)')version
          stop
!
        case ('-h')
          call print_help()
!
        case ('-b')
          call get_command_argument(i+1, arg)
          icycle=1
          if(len_trim(bfile) > 0)then
            write(*,'("error. base file set again.")')
            ierr=ierr+1
          else
            if(len_trim(arg) == 0)then
              write(*,'("error reading -b file name.")')
              ierr=ierr+1
            end if
            bfile=arg
          end if
!
        case ('-n')
          call get_command_argument(i+1, arg)
          icycle=1
          if(lnfile)then
            write(*,'("error. negative-pert file set again.")')
            ierr=ierr+1
          else
            if(len_trim(arg) == 0)then
              write(*,'("error reading -n file name.")')
              ierr=ierr+1
            end if
            nfile=arg
            lnfile=.true.
          end if
!
        case ('-p')
          call get_command_argument(i+1, arg)
          icycle=1
          if(lpfile)then
            write(*,'("error. positive-pert file set again.")')
            ierr=ierr+1
          else
            if(len_trim(arg) == 0)then
              write(*,'("error reading -p file name.")')
              ierr=ierr+1
            end if
            pfile=arg
            lpfile=.true.
          end if
!
        case ('-r')
          call get_command_argument(i+1, arg)
          icycle=1
          if(lrelpert)then
            write(*,'("error. relative perturbation set again.")')
            ierr=ierr+1
          else
            read(arg,*,iostat=ios)relpert
            if(ios /= 0)then
              write(*,'("error reading relative perturbation.")')
              ierr=ierr+1
            end if
            lrelpert=.true.
          end if
!
        case ('-f')
          call get_command_argument(i+1, arg)
          icycle=1
          if(italratio /= 0)then
            write(*,'("error. tally set again.")')
            ierr=ierr+1
          else
            read(arg,*,iostat=ios)italratio
            if(ios /= 0)then
              write(*,'("error reading tally.")')
              ierr=ierr+1
            end if
          end if
!
        case ('-u')
          call get_command_argument(i+1, arg)
          icycle=1
          if(len_trim(cunc) > 0)then
            write(*,'("error. uncertainty output format (-u) set again.")')
            ierr=ierr+1
          else
            cunc=arg
            if(len_trim(arg) == 0)then
              write(*,'("error reading uncertainty output format (-u).")')
              ierr=ierr+1
            else if(trim(cunc) == "abs")then
              iunc=0
            else if(trim(cunc) == "rel")then
              iunc=1
            else if(trim(cunc) == "prc")then
              iunc=2
            else
              write(*,'("error. unknown uncertainty output format: ",a)')trim(cunc)
              write(*,'("  must be abs, rel, or prc.")')
              ierr=ierr+1
            end if
          end if
!
        case ('-o')
          call get_command_argument(i+1, arg)
          icycle=1
          if(lofile)then
            write(*,'("error. output file set again.")')
            ierr=ierr+1
          else
            if(len_trim(arg) == 0)then
              write(*,'("error reading -o file name.")')
              ierr=ierr+1
            end if
            ofile=arg
            lofile=.true.
          end if
!
        case default
          write(*,'(2a, /)')'unrecognized command-line option: ',trim(arg)
          call print_help()
      end select
    end do ! i
    if(ierr > 0)stop
!
! default output format for printing (iunc is already set)
    if(len_trim(cunc) == 0)then
      cunc="abs"
    end if
!
! https://fortranwiki.org/fortran/files/character_handling_in_Fortran.html
    call get_command(length=cmdlen)
    if(cmdlen > 0) then
      allocate(character(cmdlen) :: command)
      call get_command(command)
    end if
!
! 6 is the unit number for standard out
    call write_head1(6,relpert,italratio,cmdlen,version,command,bfile,nfile,pfile,cunc,ofile)
!
    ierr=0
    if(len_trim(bfile) == 0)then
      write(*,'("error. base file is required.")')
      ierr=1
    else
      ii=index(bfile,"#")
      if(ii /= 0)then
        nhash=1 ! will count them in pass1TFC
      end if
    end if
    if(relpert < 0.d0)then
      write(*,'("error. relative perturbation must be positive. ", &
       "sign is determined by -p or -n except for pert.")')
      ierr=1
    end if
    if(italratio < 0)then
      write(*,'("error. tally must be positive. ")')
      ierr=1
    end if
    if(ierr > 0)stop
!
! try to determine if ofile is an mcnp output file (by user accident)
    if(lofile)then
      inquire(file=ofile,exist=lexist)
      if(lexist)then
        open(unit=iuo,file=ofile,status="old") 
        read(iuo,'(a)')line
        if(index(line,"Code Name") /= 0)then
          ierr=1
        else
          read(iuo,'(a)')line
          if(index(line,"comment") /= 0)ierr=1
          if(index(line,"warning") /= 0)ierr=1
          if(index(line," 1- ") /= 0)ierr=1
          if(index(line," 2- ") /= 0)ierr=1
        end if
        close(iuo)
        if(ierr == 1)then
          write(*,'("error. output file ",a," appears to be an existing MCNP output file.", &
           " will not overwrite.")')trim(ofile)
          stop
        end if
      end if
      open(unit=iuo,file=ofile,status="unknown") 
      call write_head1(iuo,relpert,italratio,cmdlen,version,command,bfile,nfile,pfile,cunc,ofile)
    end if
    return
    end subroutine

    subroutine print_help()
    write(*,'("command-line options:")')
    write(*,'(4x,"-b       MCNP output file for base case")')
    write(*,'(4x,"-n       MCNP output file for negative perturbation")')
    write(*,'(4x,"-p       MCNP output file for positive perturbation")')
    write(*,'(4x,"-r       Relative perturbation used in MCNP calculations")')
    write(*,'(4x,"-f       Tally to divide other tallies by")')
    write(*,'(4x,"-u       abs/rel/prc Output absolute/relative/percent uncertainties")')
    write(*,'(4x,"-o       Correlated sampling output file to write")')
    write(*,'(4x,"-v       Print code version and exit")')
    write(*,'(4x,"-h       Print command-line options and exit")')
    stop
    end subroutine

    subroutine pass1TFC(iui,xfile,ntal,npert,nfiles,nhash,nentry,ifmt,ierr)
! get the number of tallies, perturbations, and batches;
! get the format of the TFC.
    implicit none
    integer iui,ntal,npert,nfiles,nhash,nentry,ifmt,ierr
    character xfile*120
    integer i1,i2,i3,it(3),ios,nl,i,ii,iptmp(2)
    character tfile*120,line*140,frmt*24
    logical lexist

!
    ntal=0
    npert=0
    nentry=0
    ifmt=0
!
    tfile=xfile
! for multiple files, this only checks the first one.
    if(nhash > 0)then
      nfiles=0 ! will count them at the end
! ii is where # starts, nf is the number of them
      ii=index(xfile,"#")
      nhash=ii-1
      do i=ii,len_trim(xfile)
        if(xfile(i:i) == "#")then
          nhash=nhash+1
          tfile(i:i)="0"
        else
          exit
        end if
      end do ! i
      tfile(nhash:nhash)="1"
      nhash=nhash-ii+1
    else
      nfiles=1
    end if
    open(unit=iui,file=tfile,status="old",iostat=ios) 
    if(ios /= 0)then
      write(*,'("error. output file ",a," not found.")')trim(tfile)
      ierr=1
      return
    end if
    nl=0 ! count lines for debugging if needed
 10 read(iui,'(a)',iostat=ios)line
    if(ios /= 0)then
      write(*,'("error. TFC in output file ",a," not found.")')trim(tfile)
      ierr=2
      return
    end if
    nl=nl+1
! reading this part of the file depends on output file format.
    if(line(1:25) /= "1tally fluctuation charts")go to 10
    read(iui,'(a)',iostat=ios)line ! blank
    nl=nl+1
 20 read(iui,'(a)',iostat=ios)line ! tally numbers
    nl=nl+1
    i1=0 ; i2=0 ; i3=0
    it(1:3)=0
    nentry=0
    i1=index(line,"tally")
! dbcn 49j 0
    if(i1 == 29)then
      ifmt=1
! dbcn 49j 2
    else if(i1 == 34)then
      ifmt=2
! dbcn 49j 1
    else if(i1 == 30)then
      ifmt=3
    else
      write(*,'("error. undefined TFC format.")')
      ifmt=-1
      ierr=3
      return
    end if
    read(line(i1+5:i1+13),*)it(1)
    ntal=ntal+1
    i2=index(line(i1+5:len_trim(line)),"tally")
    if(i2 > 0)then
      i2=i2+i1+4
      read(line(i2+5:i2+13),*)it(2)
      ntal=ntal+1
      i3=index(line(i2+5:len_trim(line)),"tally")
      if(i3 > 0)then
        i3=i3+i2+4
        read(line(i3+5:i3+13),*)it(3)
        ntal=ntal+1
      end if
    end if
    read(iui,'(a)',iostat=ios)line ! titles
    nl=nl+1
 30 read(iui,'(a)',iostat=ios)line
    nl=nl+1
!   write(*,'("DEBUG ",a)')trim(line)
    if(len_trim(line) == 0)then ! TFC block ends
      read(iui,'(a)',iostat=ios)line
      nl=nl+1
      if(line(2:11) == "**********")then ! entire TFC ends
        go to 90
      else ! another block of tallies
        go to 20
      end if
    else if(line(2:19) == "tally data written")then ! entire TFC ends
      go to 90
! sometimes this is repeated; why?
    else if(line(1:25) == "1tally fluctuation charts".and. &
     index(line,"perturbation") == 0)then
      read(iui,'(a)',iostat=ios)line ! blank
      nl=nl+1
      go to 20
    else if(line(1:25) == "1tally fluctuation charts".and. &
     index(line,"perturbation") > 0)then ! perturbations start
      read(line(45:len_trim(line)),*)iptmp(1)
      npert=npert+1
      go to 40
    else
      nentry=nentry+1
    end if
    go to 30

! looking for perturbations
 40 read(iui,'(a)',iostat=ios)line
    nl=nl+1
    if(line(1:25) == "1tally fluctuation charts")then
      read(line(45:len_trim(line)),*)iptmp(2)
      if(iptmp(2) /= iptmp(1))then
        npert=npert+1
      end if
      iptmp(1)=iptmp(2)
    else if(line(2:11) == "**********")then ! entire TFC ends
      go to 90
    end if
    go to 40
!
 90 close(iui)
!
! check how many files there are.
    if(nhash > 0)then
!     write(*,'("ii,nhash=",2i6,i12)')ii,nhash,10**nhash-1
      write(frmt,'("(i",i0,".",i0,")")')nhash,nhash
      do i=1,10**nhash-1
        write(tfile(ii:ii+nhash-1),frmt)i
        inquire(file=tfile,exist=lexist)
        if(lexist)then
          nfiles=nfiles+1
        else
          exit
        end if
      end do ! i
    end if
!
!   write(*,'("file=",a)')trim(xfile)
!   write(*,'("nfiles=",i0)')nfiles
!   write(*,'("nhash=",i0)')nhash
!   write(*,'("ntal=",i0)')ntal
!   write(*,'("npert=",i0)')npert
!   write(*,'("nentry=",i0)')nentry
    return
    end subroutine

    subroutine checkTFC(iui,xfile,bfile,ntal,npert,nfiles,nhash,nentry,ifmt,ierr)
! check the number of tallies, perturbations, entries in TFC, and # in file name;
! check the TFC format; check number of files and retain the smaller number.
    implicit none
    integer iui,ntal,npert,nfiles,nhash,nentry,ifmt,ierr
    character xfile*120,bfile*120
    integer ntmp(5),itmp
!
! these need to be set on entry to pass1TFC
    ntmp(3)=nfiles
    ntmp(4)=nhash
    call pass1TFC(iui,xfile,ntmp(1),ntmp(2),ntmp(3),ntmp(4),ntmp(5),itmp,ierr)
    if(ierr == 0)then ! don't check if file is wrong.
      if(ntmp(1) /= ntal)then
        write(*,'("error. number of tallies in ",a," not equal to number in ",a, &
         ": ",i0," /= ",i0)')trim(xfile),trim(bfile),ntmp(1),ntal
        ierr=1
      end if
      if(ntmp(2) /= npert)then
        write(*,'("error. number of perturbations in ",a," not equal to number in ",a, &
         ": ",i0," /= ",i0)')trim(xfile),trim(bfile),ntmp(2),npert
        ierr=2
      end if
      if(ntmp(3) /= nfiles)then
        write(*,'("warning. number of files in ",a," not equal to number in ",a, &
         ": ",i0," /= ",i0)')trim(xfile),trim(bfile),ntmp(3),nfiles
        nfiles=min(nfiles,ntmp(3))
        write(*,'(2x,"will use the smaller number, ",i0,".")')nfiles
      end if
      if(ntmp(4) /= nhash)then
        write(*,'("error. number of # in ",a," not equal to number of # in ",a, &
         ": ",i0," /= ",i0)')trim(xfile),trim(bfile),ntmp(4),nhash
        ierr=4
      end if
      if(ntmp(5) /= nentry)then
        write(*,'("error. number of TFC entries in ",a," not equal to number in ",a, &
         ": ",i0," /= ",i0)')trim(xfile),trim(bfile),ntmp(5),nentry
        ierr=5
      end if
      if(itmp /= ifmt)then
        write(*,'("error. TFC format in ",a," not equal to TFC format in ",a, &
         ".")')trim(xfile),trim(bfile)
        ierr=6
      end if
    end if
    return
    end subroutine

    subroutine readTFC(iui,iuo,lofile,xfile,bfile,ifmt,ntal,npert,nentry,nbatch,nps,nfiles,nhash,ital,ipert,tally,rtimes)
! set up to read all TFCs. the actual TFC reading is done in rd_tfc_lines.
! read timing information.
! compare tally names and nps numbers.
    implicit none
    integer iui,iuo,ifmt,ntal,npert,nentry,nbatch,nps(0:nbatch),nfiles,nhash,ital(ntal),ipert(npert)
    real*8 tally(2,ntal,0:npert,nbatch+1),rtimes
    character xfile*120,bfile*120
    logical lofile
    real*8 tally0(2,ntal,nentry)
    integer i,ii,ierr,ios,ip,ital0(ntal),ipert0(npert),nps0(nentry),itime(12),nf,ilast,itimes,nf1,nbatch1
    character tfile*120,line*140,frmt*24
!
    nbatch1=nbatch+1
    tfile=xfile
    if(nhash > 0)then
      ii=index(xfile,"#")
      write(frmt,'("(i",i0,".",i0,")")')nhash,nhash
    end if
    nf1=0
    ierr=0
    do nf=1,nfiles
      if(nhash > 0)write(tfile(ii:ii+nhash-1),frmt)nf
      tally0(1:2,1:ntal,1:nentry)=0.d0
      nps0(1:nentry)=0
      write(*,'("reading ",a,"...")')trim(tfile)
      open(unit=iui,file=tfile,status="old",iostat=ios) 
! this was already checked in pass1TFC but something could have happened.
      if(ios /= 0)then
        write(*,'("error in readTFC. output file ",a," not found.")')trim(tfile)
        return
      end if
      nf1=nf1+1
      ilast=(nf-1)*nentry
 10   read(iui,'(a)')line
      if(line(1:25) /= "1tally fluctuation charts")go to 10
      backspace(iui)
      ip=0
      ital0(1:ntal)=0
      ipert0(1:npert)=0
      call rd_tfc_lines(iui,ifmt,ntal,nentry,nps0,ital0,ip,tally0)
      tally(1:2,1:ntal,ip,ilast+1:ilast+nentry)=tally0(1:2,1:ntal,1:nentry)
! the ultimate tally
      tally(1,1:ntal,ip,nbatch1)=tally(1,1:ntal,ip,nbatch1)+tally0(1,1:ntal,nentry)
      tally(2,1:ntal,ip,nbatch1)=tally(2,1:ntal,ip,nbatch1)+(tally0(1,1:ntal,nentry)*tally0(2,1:ntal,nentry))**2
      if(ital(1) == 0)then
        ital(1:ntal)=ital0(1:ntal)
      else
! error checking
        do i=1,ntal
          if(ital0(i) /= ital(i))then
            write(*,'("error. tally numbers in ",a," not equal to ", &
             "tally numbers in ",a,".")')trim(tfile),trim(bfile)
            if(lofile)then
              write(iuo,'("error. tally numbers in ",a," not equal to ", &
               "tally numbers in ",a,".")')trim(tfile),trim(bfile)
            end if
            ierr=ierr+1
            exit
          end if
        end do ! i
      end if
      if(nps(ilast+nentry) == 0)then
        nps(ilast+1:ilast+nentry)=nps0(1:nentry)
      else
! error checking
        do i=1,nentry
          if(nps0(i) /= nps(i))then
            write(*,'("error. nps numbers in ",a," not equal to ", &
             "nps numbers in ",a,".")')trim(tfile),trim(bfile)
            if(lofile)then
              write(iuo,'("error. nps numbers in ",a," not equal to ", &
               "nps numbers in ",a,".")')trim(tfile),trim(bfile)
            end if
            ierr=ierr+1
            exit
          end if
        end do ! i
      end if
      do ip=1,npert
        call rd_tfc_lines(iui,ifmt,ntal,nentry,nps0,ital0,ipert0(ip),tally0)
        tally(1:2,1:ntal,ip,ilast+1:ilast+nentry)=tally0(1:2,1:ntal,1:nentry)
! the ultimate tally
        tally(1,1:ntal,ip,nbatch1)=tally(1,1:ntal,ip,nbatch1)+tally0(1,1:ntal,nentry)
        tally(2,1:ntal,ip,nbatch1)=tally(2,1:ntal,ip,nbatch1)+(tally0(1,1:ntal,nentry)*tally0(2,1:ntal,nentry))**2
        if(ipert(npert) == 0)then
          ipert(ip)=ipert0(ip)
        else
! error checking
          if(ipert0(ip) /= ipert(ip))then
            write(*,'("error. perturbation numbers in ",a," not equal to ", &
             "perturbation numbers in ",a,".")')trim(tfile),trim(bfile)
            if(lofile)then
              write(iuo,'("error. perturbation numbers in ",a," not equal to ", &
               "perturbation numbers in ",a,".")')trim(tfile),trim(bfile)
            end if
            ierr=ierr+1
            exit
          end if
        end if
      end do ! ip
!
! read timing information.
 20   read(iui,'(a)',iostat=ios)line
      if(ios /= 0)go to 90
      if(index(line,"version") == 0.and.index(line,"probid") == 0)go to 20
!     write(*,'(a)')trim(line)
      i=index(line,"/")
! this line depends on output file format. end date starts 27 characters
! after first "/" in loddat.
! indices in itime: end time             start time
!                   01/02/03 04:05:06    07/08/09 10:11:12
      read(line(i+27:len_trim(line)),'(i2,1x,i2,1x,i2,i3,1x,i2,1x,i2, &
       31x,i2,1x,i2,1x,i2,i3,1x,i2,1x,i2)')itime(1:12)
      call calc_time(itime(1:6),itime(7:12),itimes)
      rtimes=rtimes+real(itimes,8)
 90   close(iui)
    end do ! nf
    if(nf1 > 1)write(*,'(i0," files read.")')nf1
! the ultimate tally relative uncertainty and average
    do i=1,ntal
      do ip=0,npert
        if(tally(1,i,ip,nbatch1) /= 0.d0)then
          tally(2,i,ip,nbatch1)=sqrt(tally(2,i,ip,nbatch1))/tally(1,i,ip,nbatch1)
        end if
        tally(1,i,ip,nbatch1)=tally(1,i,ip,nbatch1)/real(nfiles,8)
      end do ! ip
    end do ! i
!
    if(ierr > 0)then
      if(lofile)close(iuo)
      stop
    end if
    return
    end subroutine

    subroutine rd_tfc_lines(iui,ifmt,ntal,nentry,nps0,ital0,ipert0,tfc)
! read the nps, tally mean, and standard deviation from the TFC.
! read tally names and perturbation names.
! reading this depends on output file format.
    implicit none
    integer iui,ifmt,ntal,nentry,nps0(nentry),ital0(ntal),ipert0,ierr1,ierr2
    real*8 tfc(2,ntal,nentry)
    integer i,ib,npb,n,n0
    character line*140
!
! npb=max number of tallies on a line in TFC
    if(ifmt == 1)then
      npb=3
    else if(ifmt == 2.or.ifmt == 3)then
      npb=2
    end if
    ierr1=0 ; ierr2=0
    n0=1
    do i=1,ntal,npb
!     write(*,'(i0)')i
      read(iui,'(a)')line ! 1tally fluctuation charts
      if(len_trim(line(45:54)) /= 0)then
        if(line(54:54) == "*")then
          if(ierr1 == 0)then
            write(*,'("error. pert number too large for field.")')
          end if
          ierr1=1
          ipert0=-1
        else
          read(line(45:54),*)ipert0
        end if
      end if
      read(iui,'(a)')line ! blank
      read(iui,'(a)')line ! tally...
      n=n0
      if(ifmt == 1)then
        read(line(34:42),*)ital0(n)
        if(n < ntal)then
          n=n+1
          read(line(74:82),*)ital0(n)
          if(n < ntal)then
            n=n+1
            read(line(114:122),*)ital0(n)
          end if
        end if
      else if(ifmt == 2)then
        read(line(39:47),*)ital0(n)
        if(n < ntal)then
          n=n+1
          read(line(89:97),*)ital0(n)
        end if
      else if(ifmt == 3)then
        read(line(35:43),*)ital0(n)
        if(n < ntal)then
          n=n+1
          read(line(81:89),*)ital0(n)
        end if
      end if
      read(iui,'(a)')line ! titles
      do ib=1,nentry
        read(iui,'(a)')line
!       write(*,'(a)')trim(line)
        if(line(13:13) == "*")then
          if(ierr2 == 0)then
            write(*,'("error. nps too large for field.")')
          end if
          ierr2=1
          nps0(ib)=-1
        else
          read(line(1:13),'(i13)')nps0(ib)
        end if
        n=n0
        if(ifmt == 1)then
          read(line(14:33),*)tfc(1,n,ib),tfc(2,n,ib)
          if(n < ntal)then
            n=n+1
            read(line(54:73),*)tfc(1,n,ib),tfc(2,n,ib)
            if(n < ntal)then
              n=n+1
              read(line(94:113),*)tfc(1,n,ib),tfc(2,n,ib)
            end if
          end if
        else if(ifmt == 2)then
          read(line(14:38),*)tfc(1,n,ib),tfc(2,n,ib)
          if(n < ntal)then
            n=n+1
            read(line(60:88),*)tfc(1,n,ib),tfc(2,n,ib)
          end if
        else if(ifmt == 3)then
          read(line(14:34),*)tfc(1,n,ib),tfc(2,n,ib)
          if(n < ntal)then
            n=n+1
            read(line(60:80),*)tfc(1,n,ib),tfc(2,n,ib)
          end if
        end if
      end do ! ib
      n0=n0+npb
    end do ! i
!   do ib=1,nentry
!     write(*,'(i12,1pe15.8,0pf8.5)')ib,tfc(1:2,ntal,ib)
!   end do ! ib
    return
    end subroutine

    subroutine calc_diff(n0,maxinfo,tally2,f2,tally0,f0,diff)
! calculate the difference: (2-0)
    implicit none
    integer n0,maxinfo
    real*8 tally2(2),f2(n0),tally0(2),f0(n0),diff(maxinfo)
    integer n
    real*8 v0,v2,t1,rn0,Q0,Q2
!
    diff(1)=tally2(1)-tally0(1)
    v0=(tally0(1)*tally0(2))**2
    v2=(tally2(1)*tally2(2))**2
! Eq. (2)
    diff(4)=sqrt(v0+v2)
    do n=1,n0
      t1=(f0(n)-tally0(1))**2+(f2(n)-tally2(1))**2
      diff(3)=diff(3)+t1
! Eq. (4)
      diff(2)=diff(2)+t1 &
       -2.d0*(f0(n)-tally0(1))*(f2(n)-tally2(1))
      diff(7)=diff(7)+abs(f0(n)-f2(n))-(f0(n)-f2(n))**2
      diff(8)=diff(8)+f2(n)+f0(n)-f2(n)**2-f0(n)**2
    end do ! n
    rn0=real(n0*(n0-1),8)
    diff(3)=sqrt(diff(3)/rn0)
    diff(2)=sqrt(diff(2)/rn0)
!   diff(7)=sqrt(diff(7)/real((n0-1),8))
!   diff(8)=sqrt(diff(8)/real((n0-1),8))
! theoretical uncertainty. need to divide by sqrt(nps-1).
    Q0=1.d0 ; Q2=1.d0
!   Q0=2.d0 ; Q2=3.d0
!   Q0=1.d0 ; Q2=0.954198473282d0
! correlated.
    diff(5)=(abs(diff(1))-diff(1)**2)
!   diff(5)=(Q2*tally2(1)+Q0*tally0(1) &
!    -2.d0*Q2*Q0*min(tally2(1)/Q2,tally0(1)/Q0)-diff(1)**2)
    if(diff(5) > 0.d0)diff(5)=sqrt(diff(5))
! uncorrelated.
    diff(6)=(tally2(1)+tally0(1)-tally2(1)**2-tally0(1)**2)
!   diff(6)=((Q2*tally2(1)+Q0*tally0(1)-tally2(1)**2-tally0(1)**2))
    if(diff(6) > 0.d0)diff(6)=sqrt(diff(6))
    return
    end subroutine

    subroutine calc_pert(n0,maxinfo,tally1,f1,tally0,f0,c1,diff,sens)
! calculate relative sensitivity using ratio from pert: 1/0
    implicit none
    integer n0,maxinfo
    real*8 tally1(2),f1(n0),tally0(2),f0(n0),c1,diff(maxinfo),sens(maxinfo)
    integer n
    real*8 v0,v1,t2,rn0
!
    diff(1)=tally1(1)
    diff(4)=abs(tally1(1)*tally1(2))
! Eq. (11)
    sens(1)=tally1(1)/tally0(1)*c1
    v0=tally0(2)**2
    v1=tally1(2)**2
! Eqs. (13) and (7)
    sens(4)=abs(sens(1))*sqrt(v0+v1)
    do n=1,n0
      t2=(f1(n)/tally1(1)-1.d0)**2+(f0(n)/tally0(1)-1.d0)**2
      sens(3)=sens(3)+t2
! Eqs. (13) and (9)
      sens(2)=sens(2)+t2 &
       -2.d0*(f1(n)/tally1(1)-1.d0)*(f0(n)/tally0(1)-1.d0)
    end do ! n
    rn0=real(n0*(n0-1),8)
    sens(3)=sqrt(sens(3)/rn0)*abs(sens(1))
    sens(2)=sqrt(sens(2)/rn0)*abs(sens(1))
    return
    end subroutine

    subroutine calc_ratio(n0,maxinfo,tally1,f1,tally0,f0,diff)
! calculate ratio: 1/0
    implicit none
    integer n0,maxinfo
    real*8 tally1(2),f1(n0),tally0(2),f0(n0),diff(maxinfo)
    integer n
    real*8 v0,v1,t2,rn0
!
! Eq. (6)
    diff(1)=tally1(1)/tally0(1)
    v0=tally0(2)**2
    v1=tally1(2)**2
! Eq. (7)
    diff(4)=abs(diff(1))*sqrt(v0+v1)
    do n=1,n0
      t2=(f1(n)/tally1(1)-1.d0)**2+(f0(n)/tally0(1)-1.d0)**2
      diff(3)=diff(3)+t2
! Eq. (9)
      diff(2)=diff(2)+t2 &
       -2.d0*(f1(n)/tally1(1)-1.d0)*(f0(n)/tally0(1)-1.d0)
    end do ! n
    rn0=real(n0*(n0-1),8)
    diff(3)=sqrt(diff(3)/rn0)*abs(diff(1))
    diff(2)=sqrt(diff(2)/rn0)*abs(diff(1))
    return
    end subroutine

    subroutine calc_sens1(n0,maxinfo,is,tally1,f1,tally0,f0,c1,sens)
! calculate a one-sided difference estimate of the relative sensitivity:
!  (1-0)/0 or (0-1)/0, depending on is
    implicit none
    integer n0,maxinfo,is
    real*8 tally1(2),f1(n0),tally0(2),f0(n0),c1,sens(maxinfo)
    integer n
    real*8 v0,v1,t2,rn0,sss
!
! Eqs. (14) and (15)
    if(is == 1)then
      sens(1)=(tally1(1)-tally0(1))/tally0(1)*c1
    else if(is == 2)then
      sens(1)=(tally0(1)-tally1(1))/tally0(1)*c1
    end if
    sss=abs(c1*tally1(1)/tally0(1))
    v0=tally0(2)**2
    v1=tally1(2)**2
! Eqs. (13) and (7)
    sens(4)=sss*sqrt(v0+v1)
    do n=1,n0
      t2=(f1(n)/tally1(1)-1.d0)**2+(f0(n)/tally0(1)-1.d0)**2
      sens(3)=sens(3)+t2
! Eqs. (13) and (9)
      sens(2)=sens(2)+t2 &
       -2.d0*(f1(n)/tally1(1)-1.d0)*(f0(n)/tally0(1)-1.d0)
    end do ! n
    rn0=real(n0*(n0-1),8)
    sens(3)=sqrt(sens(3)/rn0)*sss
    sens(2)=sqrt(sens(2)/rn0)*sss
    return
    end subroutine

    subroutine calc_sens2(n0,maxinfo,tally2,f2,tally1,f1,tally0,f0,c2,sens)
! calculate a central-difference estimate of the relative sensitivity: (2-1)/0
    implicit none
    integer n0,maxinfo
    real*8 tally2(2),f2(n0),tally1(2),f1(n0),tally0(2),f0(n0),c2,sens(maxinfo)
    integer n
    real*8 v0,v1,v2,t2,rn0,Q0,Q1,Q2
!
! Eq. (16)
    sens(1)=(tally2(1)-tally1(1))/tally0(1)*c2
    v0=(tally0(1)*tally0(2))**2
    v1=(tally1(1)*tally1(2))**2
    v2=(tally2(1)*tally2(2))**2
! Eq. (18)
    sens(4)=sqrt(c2**2*(v1+v2)+sens(1)**2*v0)/tally0(1)
    do n=1,n0
      t2=c2**2*((f2(n)-tally2(1))**2+(f1(n)-tally1(1))**2) &
       +sens(1)**2*(f0(n)-tally0(1))**2
      sens(3)=sens(3)+t2
! Eq. (20)
      sens(2)=sens(2)+t2 &
       -2.d0*c2**2*(f2(n)-tally2(1))*(f1(n)-tally1(1)) &
       -2.d0*c2*sens(1)*(f0(n)-tally0(1)) &
       *((f2(n)-tally2(1))-(f1(n)-tally1(1)))
    end do ! n
    rn0=real(n0*(n0-1),8)
    sens(3)=sqrt(sens(3)/rn0)/tally0(1)
    sens(2)=sqrt(sens(2)/rn0)/tally0(1)
! theoretical uncertainty. need to divide by sqrt(nps-1).
    Q0=1.d0 ; Q1=1.d0 ; Q2=1.d0
!   Q1=1.d0 ; Q0=2.d0 ; Q2=3.d0
! correlated.
    sens(5)=(c2/tally0(1))**2*abs(tally2(1)-tally1(1)) &
     +sens(1)**2/tally0(1) &
     -2.d0*c2*sens(1)/tally0(1)**2 &
     *(min(tally2(1),tally0(1))-min(tally1(1),tally0(1)))
!   sens(5)=(c2/tally0(1))**2*(Q2*tally2(1)+Q1*tally1(1) &
!    -2.d0*Q2*Q1*min(tally2(1)/Q2,tally1(1)/Q1)) &
!    +sens(1)**2*Q0/tally0(1) &
!    -2.d0*c2*sens(1)*Q0/tally0(1)**2 &
!    *(Q2*min(tally2(1)/Q2,tally0(1)/Q0) &
!    -Q1*min(tally1(1)/Q1,tally0(1)/Q0))
    if(sens(5) > 0.d0)sens(5)=sqrt(sens(5))
! uncorrelated.
    sens(6)=((c2/tally0(1))**2*(Q2*tally2(1)+Q1*tally1(1) &
     -tally2(1)**2-tally1(1)**2)+(sens(1)/tally0(1))**2 &
     *(Q0*tally0(1)-tally0(1)**2))
    if(sens(6) > 0.d0)sens(6)=sqrt(sens(6))
    return
    end subroutine

    subroutine calc_diff_ratio(n0,maxinfo,tally1,f1,tally2,f2,tally3,f3,tally4,f4,diff)
! calculate the difference: (1/2)-(3/4)
    implicit none
    integer n0,maxinfo
    real*8 tally1(2),f1(n0),tally2(2),f2(n0),tally3(2),f3(n0),tally4(2),f4(n0),diff(maxinfo)
    integer n
    real*8 rat1,rat2,v1,v2,v3,v4,t1,rn0
!
! Eq. (21)
    rat1=tally1(1)/tally2(1)
    rat2=tally3(1)/tally4(1)
    diff(1)=rat1-rat2
    v1=tally1(2)**2
    v2=tally2(2)**2
    v3=tally3(2)**2
    v4=tally4(2)**2
! Eq. (22)
    diff(4)=sqrt(rat1**2*(v1+v2)+rat2**2*(v3+v4))
    do n=1,n0
      t1=rat1**2*((f1(n)/tally1(1)-1.d0)**2+(f2(n)/tally2(1)-1.d0)**2) &
       +rat2**2*((f3(n)/tally3(1)-1.d0)**2+(f4(n)/tally4(1)-1.d0)**2)
      diff(3)=diff(3)+t1
! Eq. (24)
      diff(2)=diff(2)+t1 &
       -2.d0*rat1**2*(f1(n)/tally1(1)-1.d0)*(f2(n)/tally2(1)-1.d0) &
       -2.d0*rat2**2*(f3(n)/tally3(1)-1.d0)*(f4(n)/tally4(1)-1.d0) &
       +2.d0*rat1*rat2*((f1(n)/tally1(1)-1.d0)*(f4(n)/tally4(1)-1.d0) &
       +(f2(n)/tally2(1)-1.d0)*(f3(n)/tally3(1)-1.d0) &
       -(f1(n)/tally1(1)-1.d0)*(f3(n)/tally3(1)-1.d0) &
       -(f2(n)/tally2(1)-1.d0)*(f4(n)/tally4(1)-1.d0))
    end do ! n
    rn0=real(n0*(n0-1),8)
    diff(3)=sqrt(diff(3)/rn0)
    diff(2)=sqrt(diff(2)/rn0)
    return
    end subroutine

    subroutine calc_sens1_ratio(n0,maxinfo,is,tally1,f1,tally2,f2,tally3,f3,tally4,f4,c1,sens)
! calculate a one-sided difference estimate of the relative sensitivity:
!  [(1/2)-(3/4)]/(3/4) or [(3/4)-(1/2)]/(3/4), depending on is
    implicit none
    integer n0,maxinfo,is
    real*8 tally1(2),f1(n0),tally2(2),f2(n0),tally3(2),f3(n0),tally4(2),f4(n0),c1,sens(maxinfo)
    integer n
    real*8 rat0,rat1,v1,v2,v3,v4,t1,sss,rn0
!
! Eqs. (25) and (26)
    rat0=tally3(1)/tally4(1)
    rat1=tally1(1)/tally2(1)
    if(is == 1)then
      sens(1)=c1*(rat1-rat0)/rat0
    else if(is == 2)then
      sens(1)=c1*(rat0-rat1)/rat0
    end if
    v1=tally1(2)**2
    v2=tally2(2)**2
    v3=tally3(2)**2
    v4=tally4(2)**2
! Eq. (27)
    sss=c1*abs(rat1/rat0)
    sens(4)=sss*sqrt(v1+v2+v3+v4)
    do n=1,n0
      t1=((f1(n)/tally1(1)-1.d0)**2+(f2(n)/tally2(1)-1.d0)**2) &
       +((f3(n)/tally3(1)-1.d0)**2+(f4(n)/tally4(1)-1.d0)**2)
      sens(3)=sens(3)+t1
! Eq. (29)
      sens(2)=sens(2)+t1 &
       -2.d0*(f1(n)/tally1(1)-1.d0)*(f2(n)/tally2(1)-1.d0) &
       -2.d0*(f3(n)/tally3(1)-1.d0)*(f4(n)/tally4(1)-1.d0) &
       +2.d0*((f1(n)/tally1(1)-1.d0)*(f4(n)/tally4(1)-1.d0) &
       +(f2(n)/tally2(1)-1.d0)*(f3(n)/tally3(1)-1.d0) &
       -(f1(n)/tally1(1)-1.d0)*(f3(n)/tally3(1)-1.d0) &
       -(f2(n)/tally2(1)-1.d0)*(f4(n)/tally4(1)-1.d0))
    end do ! n
    rn0=real(n0*(n0-1),8)
    sens(3)=sqrt(sens(3)/rn0)*sss
    sens(2)=sqrt(sens(2)/rn0)*sss
    return
    end subroutine

    subroutine calc_sens2_ratio(n0,maxinfo,tally1,f1,tally2,f2,tally3,f3,tally4,f4, &
     tally5,f5,tally6,f6,c2,sens)
! calculate the relative sensitivity: [(1/2)-(3/4)]/(5/6)
    implicit none
    integer n0,maxinfo
    real*8 tally1(2),f1(n0),tally2(2),f2(n0),tally3(2),f3(n0),tally4(2),f4(n0), &
     tally5(2),f5(n0),tally6(2),f6(n0),c2,sens(maxinfo)
    integer n
    real*8 rat0,rat1,rat2,dnum,v1,v2,v3,v4,v5,v6,t2,sss,rn0
!
! Eq. (30)
    rat0=tally5(1)/tally6(1)
    rat1=tally1(1)/tally2(1)
    rat2=tally3(1)/tally4(1)
    dnum=rat1-rat2
    sens(1)=c2*dnum/rat0
    v1=tally1(2)**2
    v2=tally2(2)**2
    v3=tally3(2)**2
    v4=tally4(2)**2
    v5=tally5(2)**2
    v6=tally6(2)**2
! Eq. (31)
    sss=c2/abs(rat0)
    sens(4)=sss*sqrt(rat1**2*(v1+v2)+rat2**2*(v3+v4)+dnum**2*(v5+v6))
    do n=1,n0
      t2=rat1**2*((f1(n)/tally1(1)-1.d0)**2+(f2(n)/tally2(1)-1.d0)**2) &
       +rat2**2*((f3(n)/tally3(1)-1.d0)**2+(f4(n)/tally4(1)-1.d0)**2) &
       +dnum**2*((f5(n)/tally5(1)-1.d0)**2+(f6(n)/tally6(1)-1.d0)**2)
      sens(3)=sens(3)+t2
! Eq. (33)
      sens(2)=sens(2)+t2 &
       -2.d0*rat1**2*(f1(n)/tally1(1)-1.d0)*(f2(n)/tally2(1)-1.d0) &
       -2.d0*rat2**2*(f3(n)/tally3(1)-1.d0)*(f4(n)/tally4(1)-1.d0) &
       -2.d0*dnum**2*(f5(n)/tally5(1)-1.d0)*(f6(n)/tally6(1)-1.d0) &
       +2.d0*rat1*rat2*((f1(n)/tally1(1)-1.d0)*(f4(n)/tally4(1)-1.d0) &
       +(f2(n)/tally2(1)-1.d0)*(f3(n)/tally3(1)-1.d0) &
       -(f1(n)/tally1(1)-1.d0)*(f3(n)/tally3(1)-1.d0) &
       -(f2(n)/tally2(1)-1.d0)*(f4(n)/tally4(1)-1.d0)) &
       +2.d0*dnum*rat1*((f1(n)/tally1(1)-1.d0)*(f6(n)/tally6(1)-1.d0) &
       +(f2(n)/tally2(1)-1.d0)*(f5(n)/tally5(1)-1.d0) &
       -(f1(n)/tally1(1)-1.d0)*(f5(n)/tally5(1)-1.d0) &
       -(f2(n)/tally2(1)-1.d0)*(f6(n)/tally6(1)-1.d0)) &
       -2.d0*dnum*rat2*((f3(n)/tally3(1)-1.d0)*(f6(n)/tally6(1)-1.d0) &
       +(f4(n)/tally4(1)-1.d0)*(f5(n)/tally5(1)-1.d0) &
       -(f3(n)/tally3(1)-1.d0)*(f5(n)/tally5(1)-1.d0) &
       -(f4(n)/tally4(1)-1.d0)*(f6(n)/tally6(1)-1.d0))
    end do ! n
    rn0=real(n0*(n0-1),8)
    sens(3)=sqrt(sens(3)/rn0)*sss
    sens(2)=sqrt(sens(2)/rn0)*sss
    return
    end subroutine

    subroutine allocate_arrays1(nbatch,ntal,npert, &
     nps,ital,ipert, &
     fbase1,fpos1,fneg1, &
     tally_base,tally_pos,tally_neg)
! attempt to allocate arrays "properly" but this does not work.
    integer, allocatable, dimension(:) :: nps,ital,ipert
    real*8,  allocatable, dimension(:) :: fbase1,fpos1,fneg1
    real*8,  allocatable, dimension(:,:,:,:) :: tally_base,tally_pos,tally_neg
    logical :: force_alloc = .FALSE.

    ierr=0

    if (force_alloc) then
       deallocate(nps, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'nps(0:nbatch)'
          stop
       end if
    end if
    if (.NOT. allocated(nps)) &
         allocate(nps(0:nbatch), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'nps(0:nbatch)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(ital, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'ital(1:ntal)'
          stop
       end if
    end if
    if (.NOT. allocated(ital)) &
         allocate(ital(1:ntal), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'ital(1:ntal)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(ipert, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'ipert(1:npert)'
          stop
       end if
    end if
    if (.NOT. allocated(ipert)) &
         allocate(ipert(1:npert), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'ipert(1:npert)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(fbase1, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'fbase1(1:nbatch)'
          stop
       end if
    end if
    if (.NOT. allocated(fbase1)) &
         allocate(fbase1(1:nbatch), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'fbase1(1:nbatch)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(fpos1, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'fpos1(1:nbatch)'
          stop
       end if
    end if
    if (.NOT. allocated(fpos1)) &
         allocate(fpos1(1:nbatch), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'fpos1(1:nbatch)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(fneg1, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'fneg1(1:nbatch)'
          stop
       end if
    end if
    if (.NOT. allocated(fneg1)) &
         allocate(fneg1(1:nbatch), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'fneg1(1:nbatch)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(tally_base, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'tally_base(1:2,1:ntal,0:npert,1:nbatch+1)'
          stop
       end if
    end if
    if (.NOT. allocated(tally_base)) &
         allocate(tally_base(1:2,1:ntal,0:npert,1:nbatch+1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'tally_base(1:2,1:ntal,0:npert,1:nbatch+1)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(tally_pos, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'tally_pos(1:2,1:ntal,0:npert,1:nbatch+1)'
          stop
       end if
    end if
    if (.NOT. allocated(tally_pos)) &
         allocate(tally_pos(1:2,1:ntal,0:npert,1:nbatch+1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'tally_pos(1:2,1:ntal,0:npert,1:nbatch+1)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    if (force_alloc) then
       deallocate(tally_neg, STAT=ierr)
       if(ierr /= 0)then
          write(*,'("error. cannot deallocate array: ",a,".")') &
               'tally_neg(1:2,1:ntal,0:npert,1:nbatch+1)'
          stop
       end if
    end if
    if (.NOT. allocated(tally_neg)) &
         allocate(tally_neg(1:2,1:ntal,0:npert,1:nbatch+1), STAT=ierr)
    if(ierr /= 0)then
       write(*,'("error. cannot allocate array: ",a,".")') &
            'tally_neg(1:2,1:ntal,0:npert,1:nbatch+1)'
       write(*,'("ierror=",i3)')ierr
       stop
    end if

    return
    end subroutine

    subroutine write_head1(iu,relpert,italratio,cmdlen,version,command,bfile,nfile,pfile,cunc,ofile)
    implicit none
    integer iu,italratio,cmdlen
    real*8 relpert
    character version*8,bfile*120,nfile*120,pfile*120,cunc*120,ofile*120
    character(len=cmdlen) :: command
!
    write(iu,'("COrrelated Sampling Using Batch Statistics with MCNP TFCs, version ",a)')trim(version)
    write(iu,'("command line:   ",a)')command
    write(iu,'("base file:      ",a)')trim(bfile)
    write(iu,'("negative file:  ",a)')trim(nfile)
    write(iu,'("positive file:  ",a)')trim(pfile)
    write(iu,'("relative pert:",1pe12.4)')relpert
    if(italratio == 0)then
      write(iu,'("tally for norm: <none>")')
    else
      write(iu,'("tally for norm: ",i0)')italratio
    end if
    write(iu,'("unc. output:    ",a)')trim(cunc)
    write(iu,'("output file:    ",a)')trim(ofile)
    write(iu,'("")')
    return
    end subroutine

    subroutine write_head2(iu,iratio,lpert,lsens,ital1,itn,ipert)
    implicit none
    integer iu,iratio,ital1,itn,ipert
    logical lpert,lsens
!
    if(iratio > 0)then
      write(iu,'(/,"(tally ",i0,")/(tally ",i0,")")')ital1,itn
    else
      write(iu,'(/,"tally ",i0)')ital1
    end if
    if(lpert)then
      write(iu,'("perturbation ",i0)')ipert
    end if
    if(iratio == 1)then
      write(iu,'(10x,"Ratio     ",18x,"Uncertainty")',advance="no")
    else
      write(iu,'(10x,"Difference",18x,"Uncertainty")',advance="no")
    end if
    if(lsens)then
      write(iu,'(17x,"Rel.Sens.",19x,"Uncertainty")')
    else
      write(iu,'("")')
    end if
    if(iratio == 1)then
      write(iu,'(6x,"B",3x,"Ratio     ",4x,"Cov,batch",5x,"NoCov,batch",3x,"NoCov,nobch")',advance="no")
    else
      write(iu,'(6x,"B",3x,"Difference",4x,"Cov,batch",5x,"NoCov,batch",3x,"NoCov,nobch")',advance="no")
    end if
    if(lsens)then
      write(iu,'(3x,"Rel.Sens.",5x,"Cov,batch",5x,"NoCov,batch",3x,"NoCov,nobch")')
    else
      write(iu,'("")')
    end if
    return
    end subroutine

    subroutine calc_time(iend,ibeg,itimes)
    implicit none
    integer iend(6),ibeg(6),itimes
    integer i,idm(12) ! days in each month
    data idm/31,28,31,30,31,30,31,31,30,31,30,31/
! indices: end time             start time
!          01/02/03 04:05:06    01/02/03 04:05:06
! avoid Y3K
    if(iend(3) < ibeg(3))then
      iend(3)=iend(3)+100
    end if
! days from the beginning of the year, accounting for leap year
    iend(2)=iend(2)+sum(idm(1:iend(1)-1))
    if(mod(iend(3),4) == 0.and.iend(1).gt.2)then
      iend(2)=iend(2)+1
    end if
    ibeg(2)=ibeg(2)+sum(idm(1:ibeg(1)-1))
    if(mod(ibeg(3),4) == 0.and.ibeg(1).gt.2)then
      ibeg(2)=ibeg(2)+1
    end if
! iend(2) will be days from beginning of ibeg year
    do i=ibeg(3),iend(3)-1
      iend(2)=iend(2)+365
      if(mod(i,4) == 0)then
        iend(2)=iend(2)+1
      end if
    end do ! i
    itimes=(iend(2)-ibeg(2))*86400 &
     +(iend(4)-ibeg(4))*3600+(iend(5)-ibeg(5))*60+(iend(6)-ibeg(6))
!   write(*,'("time:",12i6,i12," s")')iend(1:6),ibeg(1:6),itimes
    return
    end subroutine
