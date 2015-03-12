PROGRAM MAIN
  !use libgomp
  IMPLICIT NONE
  integer :: I,J,K,L,IK,KMOD,N,M,A,B,C,D
  integer :: fileLu,IOS,iii,ifile
  real(8),pointer :: EVocc(:)       !do not change to allocatable
  real(8),pointer :: EVvirt(:)      !do not change to allocatable
  real(8),pointer :: UoccEOST(:,:)  !do not change to allocatable
  real(8),pointer :: UvirtT(:,:)    !do not change to allocatable
  real(8),pointer :: UoccT(:,:)  !do not change to allocatable
  real(8),pointer :: UvirtEOST(:,:)    !do not change to allocatable
  real(8),pointer :: tocc(:,:,:,:)  !do not change to allocatable
  real(8),pointer :: toccTest(:,:,:,:)  !do not change to allocatable
  real(8),pointer :: tocc2(:,:,:,:)     !do not change to allocatable
  real(8),pointer :: tocc3(:,:,:,:)     !do not change to allocatable
  real(8),pointer :: toccEOS(:,:,:,:)     !do not change to allocatable
  real(8),pointer :: tvirt(:,:,:,:)  !do not change to allocatable
  real(8),pointer :: tvirtTest(:,:,:,:)  !do not change to allocatable
  real(8),pointer :: tvirt2(:,:,:,:)     !do not change to allocatable
  real(8),pointer :: tvirt3(:,:,:,:)     !do not change to allocatable
  real(8),pointer :: tvirtEOS(:,:,:,:)     !do not change to allocatable
  real(8),pointer :: goccEOS(:,:,:,:)   !do not change to allocatable
  real(8),pointer :: gvirtEOS(:,:,:,:)   !do not change to allocatable
  real(8),pointer :: blad(:,:,:,:)   !do not change to allocatable
  real(8),pointer :: djik(:,:,:,:)   !do not change to allocatable
  real(8),pointer :: CalphaOcc(:,:,:)   !do not change to allocatable
  real(8),pointer :: CalphaVV(:,:,:)   !do not change to allocatable
  real(8),pointer :: Calpha(:,:,:)   !do not change to allocatable
  real(8),pointer :: Calpha2(:,:,:)   !do not change to allocatable
  real(8),pointer :: Calpha3(:,:,:)   !do not change to allocatable
  real(8),pointer :: Calpha4(:,:,:)   !do not change to allocatable
  integer :: nvirt,nocc,noccEOS,NBA,nvirtEOS,nocctot
  logical :: Debug,first_order
  character(len=1) :: NR(3)
  fileLu = 1345
  do ifile=1,2
     IF(ifile.EQ.1)THEN        
        OPEN(UNIT=fileLu,FILE='Calpha1.data',STATUS='OLD',&
             & FORM='UNFORMATTED',IOSTAT=IOS)
     ELSE
        OPEN(UNIT=fileLu,FILE='Calpha2.data',STATUS='OLD',&
             & FORM='UNFORMATTED',IOSTAT=IOS)
     ENDIF
     rewind fileLu
     read(fileLu) nvirt,nocc,noccEOS,nvirtEOS,NBA,nocctot
     print*,'nvirt,nocc,noccEOS,NBA',nvirt,nocc,noccEOS,nvirtEOS,NBA,nocctot

     allocate(Calpha(NBA,nvirt,nocc))
     read(fileLu) Calpha
     allocate(EVocc(nocc))
     read(fileLu) EVocc
     allocate(EVvirt(nvirt))
     read(fileLu) EVvirt

     allocate(UoccEOST(nocc,noccEOS))
     read(fileLu) UoccEOST        
     allocate(UvirtT(nvirt,nvirt))
     read(fileLu) UvirtT   

     allocate(UoccT(nocc,nocc))
     read(fileLu) UoccT        
     allocate(UvirtEOST(nvirt,nvirtEOS))
     read(fileLu) UvirtEOST   

     first_order = .TRUE.

     !=====================================================================================
     !  Major Step 5: Generate toccEOS(nvirt,noccEOS,nvirt,noccEOS)
     !=====================================================================================

     allocate(tocc(nocc,noccEOS,nvirt,nvirt))
     !Calculate and partial transform to local basis - transform 1 occupied indices (IDIAG,JLOC,ADIAG,BDIAG)
     call RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST)
     !Transform second occupied index (IDIAG,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BDIAG)
     M = noccEOS              !rows of Output Matrix
     N = noccEOS*nvirt*nvirt  !columns of Output Matrix
     K = nocc                 !summation dimension
     allocate(tocc2(noccEOS,noccEOS,nvirt,nvirt))
     call dgemm('T','N',M,N,K,1.0E0_8,UoccEOST,K,tocc,K,0.0E0_8,tocc2,M)
     deallocate(tocc)

     !Transform first Virtual index (ILOC,JLOC,ADIAG,BDIAG) => (ILOC,JLOC,ADIAG,BLOC)
     M = noccEOS*noccEOS*nvirt  !rows of Output Matrix
     N = nvirt                  !columns of Output Matrix
     K = nvirt                  !summation dimension
     allocate(tocc3(nvirt,nvirt,noccEOS,noccEOS))
     call dgemm('N','N',M,N,K,1.0E0_8,tocc2,M,UvirtT,K,0.0E0_8,tocc3,M)
     deallocate(tocc2)

     !Final virtual transformation and reorder to dimocc
     allocate(toccEOS(nvirt,noccEOS,nvirt,noccEOS))
     call RIMP2_calc_toccB(nvirt,noccEOS,tocc3,UvirtT,toccEOS)
     deallocate(tocc3)     

     !=====================================================================================
     !  Major Step 6: Generate tvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
     !=====================================================================================

     allocate(tvirt(nocc,nocc,nvirtEOS,nvirt)) !IDIAG,JDIAG,ALOC,BDIAG
     !Calculate and partial transform to local basis - transform occupied indices
     call RIMP2_calc_tvirtA(nvirt,nocc,nvirtEOS,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
     deallocate(EVocc)
     deallocate(EVvirt)

     !Transform first Virtual index (IDIAG,JDIAG,ALOC,BDIAG) => (IDIAG,JDIAG,ALOC,BLOC)
     M = nocc*nocc*nvirtEOS     !rows of Output Matrix
     N = nvirtEOS               !columns of Output Matrix
     K = nvirt                  !summation dimension
     allocate(tvirt2(nocc,nocc,nvirtEOS,nvirtEOS))
     call dgemm('N','N',M,N,K,1.0E0_8,tvirt,M,UvirtEOST,K,0.0E0_8,tvirt2,M)
     deallocate(tvirt)

     !Transform first occupied index (IDIAG,JDIAG,ALOC,BLOC) => (ILOC,JDIAG,ALOC,BLOC)
     M = nocc                    !rows of Output Matrix
     N = nocc*nvirtEOS*nvirtEOS  !columns of Output Matrix
     K = nocc                    !summation dimension
     allocate(tvirt3(nocc,nocc,nvirtEOS,nvirtEOS))
     call dgemm('T','N',M,N,K,1.0E0_8,UoccT,K,tvirt2,M,0.0E0_8,tvirt3,M)
     deallocate(tvirt2)

     !transform last occ index to local basis and reorder 
     allocate(tvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc))
     call RIMP2_calc_tvirtB(nvirtEOS,nocc,tvirt3,UoccT,tvirtEOS)
     deallocate(tvirt3)

     !=====================================================================================
     !  Major Step 7: Generate goccEOS(nvirt,noccEOS,nvirt,noccEOS)
     !=====================================================================================

     ! Transform Calpha(ALPHA,a,i) to local occupied index and local Virt
     ! Transform index delta to local occupied index 
     !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
     M = nba*nvirt  !rows of Output Matrix
     N = noccEOS          !columns of Output Matrix
     K = nocc             !summation dimension
     allocate(Calpha2(nba,nvirt,noccEOS))
     call dgemm('N','N',M,N,K,1.0E0_8,Calpha,M,UoccEOST,nocc,0.0E0_8,Calpha2,M)
     IF(.NOT.first_order)deallocate(UoccEOST)

     allocate(Calpha3(nba,nvirt,noccEOS))
     call RIMP2_TransAlpha1(nvirt,noccEOS,nba,UvirtT,Calpha2,Calpha3)
     deallocate(Calpha2)
     IF(.NOT.first_order)deallocate(UvirtT)

     allocate(goccEOS(nvirt,noccEOS,nvirt,noccEOS))
     call RIMP2_calc_gocc(nvirt,noccEOS,NBA,Calpha3,goccEOS)
     deallocate(Calpha3)

     !=====================================================================================
     !  Major Step 8: Generate gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
     !=====================================================================================

     ! Transform index delta to local occupied index 
     !(alphaAux;gamma,Jloc) = (alphaAux;gamma,J)*U(J,Jloc)     UoccEOST(iDIAG,iLOC)
     M = nba*nvirt  !rows of Output Matrix
     N = nocc             !columns of Output Matrix
     K = nocc             !summation dimension
     allocate(Calpha2(nba,nvirt,nocc))
     !OpenMP hopefully
     call dgemm('N','N',M,N,K,1.0E0_8,Calpha,M,UoccT,nocc,0.0E0_8,Calpha2,M)
     IF(.NOT.first_order)deallocate(UoccT)
     IF(.NOT.first_order)deallocate(Calpha)
     allocate(Calpha3(nba,nvirtEOS,nocc))
     call RIMP2_TransAlpha2(nocc,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha3)
     IF(.NOT.first_order)deallocate(UvirtEOST)
     deallocate(Calpha2)

     allocate(gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc))
     call RIMP2_calc_gvirt(nvirtEOS,nocc,NBA,Calpha3,gvirtEOS)
     deallocate(Calpha3)

     IF(first_order)THEN

        !=====================================================================================
        !  first_order prop: Generate djik(nvirtAOS,noccEOS,noccEOS,noccAOS=nocctot)
        !=====================================================================================
        !(alphaAux;nvirt,JnoccEOS) = (alphaAux;nvirt,J)*U(J,JnoccEOS)
        M = nba*nvirt        !rows of Output Matrix
        N = noccEOS          !columns of Output Matrix
        K = nocc             !summation dimension
        allocate(Calpha2(nba,nvirt,noccEOS))
        call dgemm('N','N',M,N,K,1.0E0_8,Calpha,M,UoccEOST,K,0.0E0_8,Calpha2,M)

        !(alphaAux,nvirtAOS,noccEOS) = (alphaAux;nvirt,noccEOS)*Uvirt(nvirt,nvirtAOS)
        allocate(Calpha3(nba,nvirt,noccEOS))
        call RIMP2_TransAlpha2(noccEOS,nvirt,nvirt,nba,UvirtT,Calpha2,Calpha3)
        deallocate(Calpha2)

        allocate(CalphaOcc(NBA,nocc,nocctot))
        read(fileLu) CalphaOcc

        !(alphaAux;nocc,noccAOS=nocctot) = (alphaAux;nocc,nocc)*UoccallT(nocctot,nocctot)
        M = nba*nocc         !rows of Output Matrix
        N = nocctot          !columns of Output Matrix
        K = nocctot          !summation dimension
        allocate(Calpha2(nba,nocctot,nocctot))
        call dgemm('N','N',M,N,K,1.0E0_8,CalphaOcc,M,UoccT,K,0.0E0_8,Calpha2,M)
        deallocate(CalphaOcc)
        !(alphaAux,noccEOS,noccAOS=nocctot) = (alphaAux;nocc,noccAOS=nocctot)*UoccEOST(nocc,noccEOS)
        allocate(Calpha4(nba,noccEOS,nocctot))
        call RIMP2_TransAlpha2(nocctot,nocc,noccEOS,nba,UoccEOST,Calpha2,Calpha4)
        deallocate(UoccEOST)
        deallocate(Calpha2)

        !  djikEOS(nvirtAOS,noccEOS,noccEOS,noccAOS)
        allocate(djik(nvirt,noccEOS,noccEOS,nocctot))
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirt,noccEOS,Calpha4,noccEOS,nocctot,djik)
        deallocate(Calpha3)
        deallocate(Calpha4)

        !=====================================================================================
        !  first_order prop: Generate blad(nvirtEOS,noccAOS,nvirtEOS,nvirtAOS)
        !=====================================================================================  

        !(alphaAux;nvirt,noccAOS) = (alphaAux;nvirt,nocc)*U(nocc,noccAOS)
        M = nba*nvirt        !rows of Output Matrix
        N = nocc             !columns of Output Matrix
        K = nocc             !summation dimension
        allocate(Calpha2(nba,nvirt,nocc))
        call dgemm('N','N',M,N,K,1.0E0_8,Calpha,M,UoccT,K,0.0E0_8,Calpha2,M)
        deallocate(UoccT)

        !(alphaAux,nvirtEOS,noccAOS) = (alphaAux;nvirt,noccAOS)*Uvirt(nvirt,nvirtEOS)
        allocate(Calpha3(nba,nvirtEOS,nocc))
        call RIMP2_TransAlpha2(nocc,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha3)
        deallocate(Calpha2)

        allocate(CalphaVV(NBA,nvirt,nvirt))
        read(fileLu) CalphaVV
        !(alphaAux;nvirt,nvirtAOS) = (alphaAux;nvirt,nvirt)*UvirtT(nvirt,nvirt)
        M = nba*nvirt        !rows of Output Matrix
        N = nvirt            !columns of Output Matrix
        K = nocc             !summation dimension
        allocate(Calpha2(nba,nvirt,nvirt))
        call dgemm('N','N',M,N,K,1.0E0_8,CalphaVV,M,UvirtT,K,0.0E0_8,Calpha2,M)
        deallocate(CalphaVV)
        deallocate(UvirtT)
        deallocate(Calpha)

        !(alphaAux,nvirtEOS,nvirtAOS) = (alphaAux;nvirt,nvirtAOS)*UvirtEOST(nvirt,nvirtEOS)
        allocate(Calpha4(nba,nvirtEOS,nvirt))
        call RIMP2_TransAlpha2(nvirt,nvirt,nvirtEOS,nba,UvirtEOST,Calpha2,Calpha4)
        deallocate(UvirtEOST)
        deallocate(Calpha2)

        !generate blad(nvirtEOS,noccAOS,nvirtEOS,nvirtAOS)
        allocate(blad(nvirtEOS,nocc,nvirtEOS,nvirt))
        call RIMP2_calc_gen4DimFO(NBA,Calpha3,nvirtEOS,nocc,Calpha4,nvirtEOS,nvirt,blad)
        deallocate(Calpha3)
        deallocate(Calpha4)
     ENDIF

     Debug = .TRUE.
     IF(Debug)THEN
        !verify that toccEOS is correct 
        allocate(toccTest(nvirt,noccEOS,nvirt,noccEOS))
        read(fileLu) toccTest
        do D=1,noccEOS
           do C=1,nvirt
              do B=1,noccEOS
                 do A=1,nvirt
                    IF(ABS(toccTest(A,B,C,D)-toccEOS(A,B,C,D)).GT.1.0E-10_8)THEN
                       print*,'Test(A,B,C,D)',toccTest(A,B,C,D)
                       print*,'toccEOS(A,B,C,D)',toccEOS(A,B,C,D)
                       print*,'Error ',A,B,C,D,'DIFF',ABS(toccTest(A,B,C,D)-toccEOS(A,B,C,D))
                       stop 'ERROR3 '
                    ENDIF
                 enddo
              enddo
           enddo
        enddo
        deallocate(toccTest)
        !verify that tvirtEOS is correct 
        allocate(toccTest(nvirtEOS,nocc,nvirtEOS,nocc))
        read(fileLu) toccTest
        do D=1,nocc
           do C=1,nvirtEOS
              do B=1,nocc
                 do A=1,nvirtEOS
                    IF(ABS(toccTest(A,B,C,D)-tvirtEOS(A,B,C,D)).GT.1.0E-10_8)THEN
                       print*,'Test(A,B,C,D)',toccTest(A,B,C,D)
                       print*,'tvirtEOS(A,B,C,D)',tvirtEOS(A,B,C,D)
                       print*,'Error ',A,B,C,D,'DIFF',ABS(toccTest(A,B,C,D)-tvirtEOS(A,B,C,D))
                       stop 'ERROR3 '
                    ENDIF
                 enddo
              enddo
           enddo
        enddo
        deallocate(toccTest)
        !verify that gvirtEOS is correct 
        allocate(toccTest(nvirtEOS,nocc,nvirtEOS,nocc))
        read(fileLu) toccTest
        do D=1,nocc
           do C=1,nvirtEOS
              do B=1,nocc
                 do A=1,nvirtEOS
                    IF(ABS(toccTest(A,B,C,D)-gvirtEOS(A,B,C,D)).GT.1.0E-10_8)THEN
                       print*,'Test(A,B,C,D)',toccTest(A,B,C,D)
                       print*,'gvirtEOS(A,B,C,D)',gvirtEOS(A,B,C,D)
                       print*,'Error ',A,B,C,D,'DIFF',ABS(toccTest(A,B,C,D)-gvirtEOS(A,B,C,D))
                       stop 'ERROR3 '
                    ENDIF
                 enddo
              enddo
           enddo
        enddo
        deallocate(toccTest)
        !verify that goccEOS is correct 
        allocate(toccTest(nvirt,noccEOS,nvirt,noccEOS))
        read(fileLu) toccTest
        do D=1,noccEOS
           do C=1,nvirt
              do B=1,noccEOS
                 do A=1,nvirt
                    IF(ABS(toccTest(A,B,C,D)-goccEOS(A,B,C,D)).GT.1.0E-10_8)THEN
                       print*,'Test(A,B,C,D)',toccTest(A,B,C,D)
                       print*,'goccEOS(A,B,C,D)',goccEOS(A,B,C,D)
                       print*,'Error ',A,B,C,D,'DIFF',ABS(toccTest(A,B,C,D)-goccEOS(A,B,C,D))
                       stop 'ERROR3 '
                    ENDIF
                 enddo
              enddo
           enddo
        enddo
        deallocate(toccTest)
        !verify that djik(nvirt,noccEOS,noccEOS,nocctot)
        allocate(toccTest(nvirt,noccEOS,noccEOS,nocctot))
        read(fileLu) toccTest
        do D=1,nocctot
           do C=1,noccEOS
              do B=1,noccEOS
                 do A=1,nvirt
                    IF(ABS(toccTest(A,B,C,D)-djik(A,B,C,D)).GT.1.0E-10_8)THEN
                       print*,'Test(A,B,C,D)',toccTest(A,B,C,D)
                       print*,'djik(A,B,C,D)',djik(A,B,C,D)
                       print*,'Error ',A,B,C,D,'DIFF',ABS(toccTest(A,B,C,D)-djik(A,B,C,D))
                       stop 'ERROR3 '
                    ENDIF
                 enddo
              enddo
           enddo
        enddo
        deallocate(toccTest)
        !verify that blad(nvirtEOS,nocc,nvirtEOS,nvirt) is correct
        allocate(toccTest(nvirtEOS,nocc,nvirtEOS,nvirt))
        read(fileLu) toccTest
        do D=1,nvirt
           do C=1,nvirtEOS
              do B=1,nocc
                 do A=1,nvirtEOS
                    IF(ABS(toccTest(A,B,C,D)-blad(A,B,C,D)).GT.1.0E-10_8)THEN
                       print*,'Test(A,B,C,D)',toccTest(A,B,C,D)
                       print*,'blad(A,B,C,D)',blad(A,B,C,D)
                       print*,'Error ',A,B,C,D,'DIFF',ABS(toccTest(A,B,C,D)-blad(A,B,C,D))
                       stop 'ERROR3 '
                    ENDIF
                 enddo
              enddo
           enddo
        enddo
        deallocate(toccTest)        
     ENDIF
     CLOSE(fileLu,STATUS='KEEP',iostat=ios)
     deallocate(toccEOS)        
     deallocate(tvirtEOS)        
     deallocate(goccEOS)        
     deallocate(gvirtEOS)        
     deallocate(djik)   
     deallocate(blad)
  enddo

CONTAINS
  !alphaCD(NBA,nvirt,nocc) is in the diagonal basis 
  subroutine RIMP2_calc_toccA(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST)
    implicit none
    integer,intent(in) :: nvirt,nocc,noccEOS,NBA
    real(8),intent(in) :: Calpha(NBA,nvirt,nocc)
    real(8),intent(in) :: EVocc(nocc),EVvirt(nvirt),UoccEOST(nocc,noccEOS)
    real(8),intent(inout) :: tocc(nocc,noccEOS,nvirt,nvirt)
    !
    integer :: BDIAG,ADIAG,IDIAG,JDIAG,ALPHAAUX,ILOC,JLOC
    real(8) :: gmocont,deltaEPS,TMP
    real(8) :: toccTMP(nocc)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
    !$OMP PRIVATE(BDIAG,ADIAG,IDIAG,JDIAG,&
    !$OMP         ALPHAAUX,ILOC,JLOC,gmocont,deltaEPS,toccTMP,TMP) &
    !$OMP SHARED(nvirt,nocc,noccEOS,NBA,Calpha,EVocc,EVvirt,tocc,UoccEOST)
    do BDIAG=1,nvirt
       do ADIAG=1,nvirt
          do IDIAG=1,nocc
             do JDIAG=1,nocc
                gmocont = 0.0E0_8  
                do ALPHAAUX=1,nba  
                   gmocont = gmocont + Calpha(ALPHAAUX,ADIAG,IDIAG)*Calpha(ALPHAAUX,BDIAG,JDIAG)
                enddo
                deltaEPS = EVocc(IDIAG)+EVocc(JDIAG)-EVvirt(BDIAG)-EVvirt(ADIAG)
                toccTMP(JDIAG)=gmocont/deltaEPS                
             enddo
             do jLOC=1,noccEOS
                TMP = 0.0E0_8
                do JDIAG=1,nocc
                   TMP = TMP + toccTMP(JDIAG)*UoccEOST(jDIAG,jLOC)
                enddo
                tocc(IDIAG,JLOC,ADIAG,BDIAG) = TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  END subroutine RIMP2_calc_toccA

  subroutine RIMP2_calc_tvirtA(nvirt,nocc,nvirtEOS,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
    implicit none
    integer,intent(in) :: nvirt,nocc,nvirtEOS,NBA
    real(8),intent(in) :: Calpha(NBA,nvirt,nocc)
    real(8),intent(in) :: EVocc(nocc),EVvirt(nvirt),UvirtEOST(nvirt,nvirtEOS)
    real(8),intent(inout) :: tvirt(nocc,nocc,nvirtEOS,nvirt)
    !
    integer :: BDIAG,ADIAG,IDIAG,JDIAG,ALPHAAUX,ALOC,BLOC
    real(8) :: gmocont,deltaEPS,TMP,tvirtTMP(nvirt)
    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
    !$OMP PRIVATE(BDIAG,ADIAG,IDIAG,JDIAG,ALPHAAUX,ALOC,BLOC,gmocont,deltaEPS,TMP,tvirtTMP) &
    !$OMP SHARED(nvirt,nocc,nvirtEOS,NBA,Calpha,EVocc,EVvirt,tvirt,UvirtEOST)
    do JDIAG=1,nocc
       do IDIAG=1,nocc
          do BDIAG=1,nvirt
             do ADIAG=1,nvirt
                gmocont = 0.0E0_8  
                do ALPHAAUX=1,nba  
                   gmocont = gmocont + Calpha(ALPHAAUX,ADIAG,IDIAG)*Calpha(ALPHAAUX,BDIAG,JDIAG)
                enddo
                deltaEPS = EVocc(IDIAG)+EVocc(JDIAG)-EVvirt(BDIAG)-EVvirt(ADIAG)
                tvirtTMP(ADIAG)=gmocont/deltaEPS                
             enddo
             do ALOC=1,nvirtEOS
                TMP = 0.0E0_8
                do ADIAG=1,nvirt
                   TMP = TMP + tvirtTMP(ADIAG)*UvirtEOST(ADIAG,ALOC)
                enddo
                tvirt(IDIAG,JDIAG,ALOC,BDIAG) = TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  END subroutine RIMP2_calc_tvirtA

  !tocc(occLOC,occLOC,virtDIAG,virtLOC)=(I,J,A,B) !Transform A
  subroutine RIMP2_calc_toccB(nvirt,noccEOS,tocc,UvirtT,toccEOS)
    implicit none
    integer,intent(in) :: nvirt,noccEOS
    real(8),intent(in) :: tocc(noccEOS,noccEOS,nvirt,nvirt),UvirtT(nvirt,nvirt)
    real(8),intent(inout) :: toccEOS(nvirt,noccEOS,nvirt,noccEOS)
    !local variables
    integer :: BLOC,JLOC,ILOC,ALOC,ADIAG
    real(8) :: TMP
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
    !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,ADIAG,TMP) &
    !$OMP SHARED(nvirt,noccEOS,tocc,UvirtT,toccEOS)
    do bLOC=1,nvirt
       do jLOC=1,noccEOS
          do aLOC=1,nvirt
             do iLOC=1,noccEOS
                TMP = 0.0E0_8
                do ADIAG=1,nvirt
                   TMP = TMP + tocc(ILOC,JLOC,ADIAG,BLOC)*UvirtT(ADIAG,aLOC)
                enddo
                toccEOS(ALOC,ILOC,BLOC,JLOC) = TMP              
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine RIMP2_calc_toccB

  subroutine RIMP2_calc_tvirtB(nvirtEOS,nocc,tvirt,UoccT,tvirtEOS)
    implicit none
    integer,intent(in) :: nvirtEOS,nocc
    real(8),intent(in) :: tvirt(nocc,nocc,nvirtEOS,nvirtEOS),UoccT(nocc,nocc)
    real(8),intent(inout) :: tvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
    !local variables
    integer :: BLOC,JLOC,ILOC,ALOC,JDIAG
    real(8) :: TMP
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
    !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,JDIAG,TMP) &
    !$OMP SHARED(nocc,nvirtEOS,tvirt,UoccT,tvirtEOS)
    do jLOC=1,nocc
       do bLOC=1,nvirtEOS
          do iLOC=1,nocc
             do aLOC=1,nvirtEOS
                TMP = 0.0E0_8
                do JDIAG=1,nocc
                   TMP = TMP + tvirt(ILOC,JDIAG,ALOC,BLOC)*UoccT(JDIAG,JLOC)
                enddo
                tvirtEOS(ALOC,ILOC,BLOC,JLOC) = TMP
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine RIMP2_calc_tvirtB

  subroutine RIMP2_TransAlpha1(nvirt,noccEOS,nba,UvirtT,AlphaCD4,AlphaCD5)
    implicit none
    integer,intent(in) :: nvirt,noccEOS,nba
    real(8),intent(in) :: UvirtT(nvirt,nvirt),AlphaCD4(nba,nvirt,noccEOS) 
    real(8),intent(inout) :: AlphaCD5(nba,nvirt,noccEOS)
    !local variables
    integer :: BLOC,JLOC,BDIAG,ALPHAAUX
    real(8) :: TMP
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
    !$OMP PRIVATE(BLOC,JLOC,BDIAG,ALPHAAUX,TMP) &
    !$OMP SHARED(nvirt,noccEOS,nba,UvirtT,AlphaCD4,AlphaCD5)
    do JLOC = 1,noccEOS
       do BLOC = 1,nvirt
          do ALPHAAUX = 1,nba
             TMP = 0.0E0_8
             do BDIAG = 1,nvirt
                TMP = TMP + UvirtT(BDIAG,BLOC)*AlphaCD4(ALPHAAUX,BDIAG,JLOC)
             enddo
             AlphaCD5(ALPHAAUX,BLOC,JLOC) = TMP
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine RIMP2_TransAlpha1

  !AlphaCD5(NBA,n3,n2) = UvirtEOST(n1,n3)*AlphaCD4(NBA,n1,n2)
  subroutine RIMP2_TransAlpha2(n2,n1,n3,nba,UvirtEOST,AlphaCD4,AlphaCD5)
    implicit none
    integer,intent(in) :: nba,n1,n2,n3
    real(8),intent(in) :: UvirtEOST(n1,n3)
    real(8),intent(in) :: AlphaCD4(NBA,n1,n2)
    real(8),intent(inout) :: AlphaCD5(NBA,n3,n2)
    !
    integer :: JLOC,BLOC,ALPHAAUX,BDIAG
    real(8) :: TMP
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
    !$OMP PRIVATE(BLOC,JLOC,BDIAG,ALPHAAUX,TMP) &
    !$OMP SHARED(n1,n2,n3,nba,UvirtEOST,AlphaCD4,AlphaCD5)
    do JLOC = 1,n2
       do BLOC = 1,n3
          do ALPHAAUX = 1,nba
             TMP = 0.0E0_8
             do BDIAG = 1,n1
                TMP = TMP + UvirtEOST(BDIAG,BLOC)*AlphaCD4(ALPHAAUX,BDIAG,JLOC)
             enddo
             AlphaCD5(ALPHAAUX,BLOC,JLOC) = TMP
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine RIMP2_TransAlpha2

  subroutine RIMP2_calc_gocc(nvirt,noccEOS,NBA,Calpha3,goccEOS)
    implicit none
    integer,intent(in) :: nvirt,noccEOS,NBA
    real(8),intent(in) :: Calpha3(NBA,nvirt,noccEOS)
    real(8),intent(inout) :: goccEOS(nvirt,noccEOS,nvirt,noccEOS)
    !local variables
    integer :: BLOC,JLOC,ILOC,ALOC,ALPHAAUX
    real(8) :: TMP
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
    !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
    !$OMP SHARED(nvirt,noccEOS,NBA,Calpha3,goccEOS)
    do jLOC=1,noccEOS
       do bLOC=1,nvirt
          do iLOC=1,noccEOS
             do aLOC=1,nvirt
                TMP = 0.0E0_8
                do ALPHAAUX = 1,nba
                   tmp = tmp + Calpha3(alphaAUX,ALOC,ILOC)*Calpha3(alphaAUX,BLOC,JLOC) 
                enddo
                goccEOS(ALOC,ILOC,BLOC,JLOC) = tmp
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine RIMP2_calc_gocc

  subroutine RIMP2_calc_gvirt(nvirtEOS,nocc,NBA,Calpha3,gvirtEOS)
    implicit none
    integer,intent(in) :: nvirtEOS,nocc,NBA
    real(8),intent(in) :: Calpha3(NBA,nvirtEOS,nocc)
    real(8),intent(inout) :: gvirtEOS(nvirtEOS,nocc,nvirtEOS,nocc)
    !local variables
    integer :: BLOC,JLOC,ILOC,ALOC,ALPHAAUX
    real(8) :: TMP
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
    !$OMP PRIVATE(BLOC,JLOC,ILOC,ALOC,ALPHAAUX,TMP) &
    !$OMP SHARED(nvirtEOS,nocc,NBA,Calpha3,gvirtEOS)
    do jLOC=1,nocc
       do bLOC=1,nvirtEOS
          do iLOC=1,nocc
             do aLOC=1,nvirtEOS
                TMP = 0.0E0_8
                do ALPHAAUX = 1,nba
                   tmp = tmp + Calpha3(alphaAUX,ALOC,ILOC)*Calpha3(alphaAUX,BLOC,JLOC) 
                enddo
                gvirtEOS(ALOC,ILOC,BLOC,JLOC) = tmp
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine RIMP2_calc_gvirt

  subroutine RIMP2_calc_gen4DimFO(NBA,Calpha3,n1,n2,Calpha4,n3,n4,djik)
    implicit none
    integer,intent(in) :: NBA,n1,n2,n3,n4
    real(8),intent(in) :: Calpha3(NBA,n1,n2),Calpha4(NBA,n3,n4)
    real(8),intent(inout) :: djik(n1,n2,n3,n4)
    !local variables
    integer :: A,B,C,D,ALPHAAUX
    real(8) :: TMP
    !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
    !$OMP PRIVATE(A,B,C,D,ALPHAAUX,TMP) &
    !$OMP SHARED(n1,n2,n3,n4,NBA,Calpha3,Calpha4,djik)
    do d=1,n4
       do c=1,n3
          do b=1,n2
             do a=1,n1
                TMP = 0.0E0_8
                do ALPHAAUX = 1,nba
                   tmp = tmp + Calpha3(alphaAUX,A,B)*Calpha4(alphaAUX,C,D)
                enddo
                djik(A,B,C,D) = tmp
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine RIMP2_calc_gen4DimFO

END PROGRAM MAIN

