C PropKa 2.00: a program for protein pKa predictions
C Copyright (c) 2007, Jan H. Jensen, Hui Li, Delphine C. Bas, and David M. Rogers
C
C * This library is free software; you can redistribute it and/or
C * modify it under the terms of the GNU Lesser General Public
C * License as published by the Free Software Foundation; either
C * version 2.1 of the License, or (at your option) any later version.
C *
C * This library is distributed in the hope that it will be useful,
C * but WITHOUT ANY WARRANTY; without even the implied warranty of
C * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C * Lesser General Public License for more details.
C
C For questions or comments contact jhjensen@kemi.ku.dk 
C
C      PROPKA: A PROTEIN PKA PREDICTOR
C      VERSION 1.0
C      BY HUI LI
C
C      04/25/2004, IOWA CITY
C
C      VERSION 2.0
C      BY DELPHINE C. BAS AND DAVID M. ROGERS 
C
C      11/05/2007, IOWA CITY, NANCY, COPENHAGEN
C
C*************************************************************************************
C REFERENCES:
C
C Very Fast Empirical Prediction and Rationalization of Protein pKa Values
C Hui Li, Andrew D. Robertson and Jan H. Jensen
C PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)
C
C Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes
C Delphine C. Bas, David M. Rogers and Jan H. Jensen
C PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)
C
C*************************************************************************************
C
C***********************************************************
C
C      THIS PROGRAM PREDICTS PROTEIN PKA VALUES 
C      ACCORDING TO THE EMPIRICAL RULES PROPOSED BY:
C
C      HUI LI, ANDREW D. ROBERTSON AND JAN H. JENSEN
C
C***********************************************************
C
       PROGRAM   PROPKA
C      *******   ******
C
       CHARACTER HEAD*2,NAMATM*5,NAMRES*3,SPACE4*4,
     $           SPACE2*2,NAMCAR*3,TYPCAR*6,AORB*1,
     $           NAMSDC*3,NAMBKB*3,NAMCOL*3,
     $           NAMPRT*3,NAMLIG*3,NBATM*4,
     $           NAMLCOL*3,NBLCOL*4,TYPCH*2,
     $           NAMt*5,NAMLN1*5,NAMLN2*5,NAMLN*5,
     $           TYPHIS*6,TYPCYS*6,TYPLYS*6,TYPTYR*6,
     $           TYPARG*6,
     $           TYPLGN3*6,TYPLGNar*6,TYPLGCAR*6,
     $           TYPLGCg*6,TYPLGCN2*6,
     $           SPACE1*1,LABEL*4,
     $           SPACE7*7, HYDRGN*1, SPAC60*60
C
       LOGICAL   CONV,TEST1,TEST2,TEST3,TEST
C
       integer   rmicar,rmiser,rmicys
C-PDB
       DIMENSION HEAD(100000),NAMATM(100000),
     $           NAMRES(100000),NUMRES(100000),NBATM(100000),
     $           X(100000),Y(100000),Z(100000),OCC(100000), 
     $           PKALIG(100000),TYPCH(100000),
C-BACKBONE O=C-N-H
     $           NUMPRT(10000),NAMPRT(10000),
     $           XPRTON(10000),YPRTON(10000),ZPRTON(10000),
     $           XNITRN(10000),YNITRN(10000),ZNITRN(10000),
     $           XCARBN(10000),YCARBN(10000),ZCARBN(10000),
     $           XOXYGN(10000),YOXYGN(10000),ZOXYGN(10000),
C - PKA SITES -
C-CARBOXYL
     $           NAMCAR(1000),TYPCAR(1000),PKACAR(1000),
     $           LCARO1(1000),LCARO2(1000),LCARRS(1000),
     $           NSDCAR(1000),NBKCAR(1000),NCLCAR(1000),
     $           PK1CAR(1000),PK2CAR(1000),NHBCAR(1000),
C-HIS
     $           XHISP1(1000),YHISP1(1000),ZHISP1(1000),
     $           XHISP2(1000),YHISP2(1000),ZHISP2(1000),
     $           LHISCG(1000),LHISND(1000),LHISCE(1000),
     $           LHISNE(1000),LHISCD(1000),LHISRS(1000),
     $           TYPHIS(1000),PKAHIS(1000),
     $           NSDHIS(1000),NBKHIS(1000),NCLHIS(1000),
     $           PK1HIS(1000),PK2HIS(1000),NHBHIS(1000),
C-CYS
     $           LCYSSG(1000),LCYSRS(1000),
     $           TYPCYS(1000),PKACYS(1000),
     $           NSDCYS(1000),NBKCYS(1000),NCLCYS(1000),
     $           PK1CYS(1000),PK2CYS(1000),
C-TYR
     $           LTYROH(1000),LTYRRS(1000),
     $           TYPTYR(1000),PKATYR(1000),
     $           NSDTYR(1000),NBKTYR(1000),NCLTYR(1000),
     $           PK1TYR(10000),PK2TYR(1000),
C-LYS
     $           LLYSNZ(1000),LLYSRS(1000),
     $           TYPLYS(1000),PKALYS(1000),
     $           NSDLYS(1000),NBKLYS(1000),NCLLYS(1000),
     $           PK1LYS(1000),PK2LYS(1000),
C-ARG
     $           XARGP1(1000),YARGP1(1000),ZARGP1(1000),
     $           XARGP2(1000),YARGP2(1000),ZARGP2(1000),
     $           XARGP3(1000),YARGP3(1000),ZARGP3(1000),
     $           XARGP4(1000),YARGP4(1000),ZARGP4(1000),
     $           XARGP5(1000),YARGP5(1000),ZARGP5(1000),
     $           LARGN1(1000),LARGN2(1000),LARGN3(1000),
     $           LARGCD(1000),LARGCZ(1000),LARGRS(1000),
     $           TYPARG(1000),PKAARG(1000),
     $           NSDARG(1000),NBKARG(1000),NCLARG(1000),
     $           PK1ARG(1000),PK2ARG(1000),
C

     $           LLIGAND(10000),
     $           NLGCAR(1000),
     $           NLGCYS(1000), 
     $           NLGTYR(1000),
     $           NLGHIS(1000),nlglys(1000),nlgarg(1000),
     $           NCLLCAR(1000),
     $           NCLLCYS(1000),
     $           NCLLTYR(1000),
     $           NCLLHIS(1000),nclllys(1000),ncllarg(1000),
C-NH groups
     $           LLIGNam(1000),
     $           XHLIG(1000),YHLIG(1000),ZHLIG(1000),
C-N aromatic
     $           LLIGNar(1000),
     $           XgLIG(1000),YgLIG(1000),ZgLIG(1000),
     $           NSDNar(1000),TYPLGNar(1000),
     $           NCLNar(1000),
     $           nbknar(1000),
     $           PK1LGNar(1000),
     $           PK2LGNar(1000),
     $           PKALGNar(1000),
C-N3 
     $           LLIGN3(1000),
     $           XN3P1(1000),YN3P1(1000),ZN3P1(1000),
     $           XN3P2(1000),YN3P2(1000),ZN3P2(1000),
     $           XN3P(1000),YN3P(1000),ZN3P(1000),
     $           NSDN3(1000),TYPLGN3(1000),
     $           NCLN3(1000),
     $           nbkn3(1000),
     $           PK1LGN3(1000),
     $           PK2LGN3(1000),
     $           PKALGN3(1000),
C-Npl
     $           LLIGNpl(1000),
     $           XLN(1000),YLN(1000),ZLN(1000),
     $           XLN1(1000),YLN1(1000),ZLN1(1000),
     $           XLN2(1000),YLN2(1000),ZLN2(1000),
     $           NAMLN1(1000),NAMLN2(1000),NAMLN(1000),
     $           XP1Np1(1000),YP1Np1(1000),ZP1Np1(1000),
     $           XP2Np1(1000),YP2Np1(1000),ZP2Np1(1000),
     $           XPNp2(1000),YPNp2(1000),ZPNp2(1000),
     $           XNg1(1000),YNg1(1000),ZNg1(1000),
     $           XNg1P(1000),YNg1P(1000),ZNg1P(1000),
     $           XNg2(1000),YNg2(1000),ZNg2(1000),
     $           XNg2P1(1000),YNg2P1(1000),ZNg2P1(1000),
     $           XNg2P2(1000),YNg2P2(1000),ZNg2P2(1000),
     $           XNg3(1000),YNg3(1000),ZNg3(1000),
     $           XNg3P1(1000),YNg3P1(1000),ZNg3P1(1000),
     $           XNg3P2(1000),YNg3P2(1000),ZNg3P2(1000),
     $           NCLCN2(1000),NCLCg(1000),NclNpl(1000),
     $           NbkCN2(1000),NbkCg(1000),NbkNpl(1000),
     $           NsdCN2(1000),NsdCg(1000),NsdNpl(1000),
     $           PK1LGC2(1000),
     $           PK1LGNp1(1000),
     $           PK1LGCg(1000),
     $           PK2LGCg(1000),
     $           PKALGCg(1000),
     $           TYPLGCg(1000),TYPLGCN2(1000),
     $           PK1LGCN2(1000),
     $           PK2LGCN2(1000),
     $           PKALGCN2(1000),
C-C3
     $           LLIGC3(1000),
C-C2
     $           LLIGC2(1000),
C-S2
     $           LLIGSo2(1000),lligs3(1000),
C-O sp2
     $           LLIGO2(1000),
C-O CARBOXYLIC ACID
     $           LLIGCAR(1000),
     $           TYPLGCAR(1000),
     $           NSDLCAR(1000),
     $           NCLLGCAR(1000),
     $           PK1LGCAR(1000),
     $           PK2LGCAR(1000),
     $           PKALGCAR(1000),
     $           XLO1(1000),YLO1(1000),ZLO1(1000),
     $           XLO2(1000),YLO2(1000),ZLO2(1000),
     $           XCac(1000),YCac(1000),ZCac(1000),
     $           NBKLCar(1000),
C-O3 groups
     $           LLIGO3(1000),
C-Cl atoms
     $           LLIGCl(1000),
C-F atoms
     $           LLIGF(1000),
C-N1 groups
     $           LLIGN1(1000),
     $           XN1G(1000),YN1G(1000),ZN1G(1000),
c dmr extra arrays
     $           lligchr(1000),lligchrval(1000),
     $           pkamod(30,1000),
c dmr
C -NONE PKA GROUPS-
C-GLN
     $           LGLNNE(1000),LGLNOE(1000),LGLNCD(1000),
     $           LGLNRS(1000),
     $           XGLNP1(1000),YGLNP1(1000),ZGLNP1(1000),
     $           XGLNP2(1000),YGLNP2(1000),ZGLNP2(1000),
C-ASN
     $           LASNND(1000),LASNOD(1000),LASNCG(1000),
     $           LASNRS(1000),
     $           XASNP1(1000),YASNP1(1000),ZASNP1(1000),
     $           XASNP2(1000),YASNP2(1000),ZASNP2(1000),
C-SER
     $           LSEROH(1000),LSERRS(1000),
C-THR
     $           LTHROH(1000),LTHRRS(1000),
C-TRP
     $           XTRPP1(1000),YTRPP1(1000),ZTRPP1(1000),
     $           LTRPCD(1000),LTRPNE(1000),LTRPCE(1000),
     $           LTRPRS(1000),
C
C-PKA DETERMINANTS
C                1: CAR (ASP/GLU)
C                2: HIS
C                3: CYS
C                4: TYR
C                5: LYS
C                6: ARG
C-PKA GROUP OF THE LIGAND
C                7: N3
C                8: Nar
C                9: Npl
C               11: Oco
     $           NMASS(20,1000),NLOCAL(20,1000),
     $           NAMSDC(20,1000,30),NUMSDC(20,1000,30),
     $           VALSDC(20,1000,30),
     $           NAMBKB(20,1000,30),NUMBKB(20,1000,30),
     $           VALBKB(20,1000,30),
     $           NAMCOL(20,1000,30),NUMCOL(20,1000,30),
     $           VALCOL(20,1000,30),
     $           TOLSDC(20,1000),TOLBKB(20,1000),TOLCOL(20,1000),
     $           TOLMAS(20,1000),TOLLOC(20,1000),
C- LIGAND PKA SITES
     $           VALLIG(20,1000,30),
     $           LABEL(20,1000,30),NAMLIG(20,1000,30),
     $           VALLCOL(20,1000,30),
     $           NAMLCOL(20,1000,30),NBLCOL(20,1000,30) 
C
C
C      **********************
C      STEP 1. READ PDB FILES
C      **********************
C
       NATOM=0
       NPRTON=1
       NCAR=0
       NLYS=0
       NTYR=0
       NSER=0
       NTHR=0
       NARG=0
       NASN=0
       NGLN=0
       NHIS=0
       NCYS=0
       NTRP=0
       NLIGAND=0
       NLIGNam=0
       NLIGC2=0
       NLIGNar=0
       NLIGO2=0
       NLIGOh=0
       NLIGN3=0
       NLIGC3=0
       NLIGNpl=0
       NLIGCAR=0
       nligchr=0
C
C      - REWRITE PDB FILE -
C
       NATOM=0
       OPEN (13, FILE='PROPKATMP2', STATUS='OLD', FORM='FORMATTED',
     $       ACCESS='SEQUENTIAL')
       DO J=1,100000
         READ (13, '(A2,A4,A7,A1,A)', END=100) 
     $        HEAD(1),NBATM(1),SPACE7,HYDRGN,SPAC60
         IF ((HEAD(1).EQ.'AT' .OR. HEAD(1).EQ.'HE'
     $                        .OR. HEAD(1).EQ.'LG') .AND.
     $        HYDRGN.NE.'H') THEN
           NATOM=NATOM+1
         END IF
       END DO
 100   CONTINUE
       CLOSE (13)
C
       IF(NATOM.LT.2)THEN
         WRITE(6,*)'PLEASE CHECK YOUR INPUT PDB FILE!'
         STOP
       END IF 
C
C      - READ INFORMATION -
C
       OPEN (10, FILE='PROPKATMP2', STATUS='OLD', FORM='FORMATTED', 
     $       ACCESS='SEQUENTIAL')
       OPEN (11, FILE='PROPKATMP3', STATUS='NEW', FORM='FORMATTED',
     $       ACCESS='SEQUENTIAL')
C
       DO I = 1, NATOM
         READ(10,'(A2,A4,I5,A5,A1,A3,a1,A1,I4,A4,F8.3,F8.3,F8.3,
     $             A2,F4.2,   F6.2)',END=200)
     $              HEAD(I),NBATM(I), NUMATM, NAMATM(I), AORB,
     $              NAMRES(I),space1,TYPCH(I), 
     $              NUMRES(I), SPACE4, X(I), Y(I), Z(I), SPACE2,
     $              OCC(I),     PKALIG(I)
C
         IF(NAMRES(I).EQ.'ASP' .OR. NAMRES(I).EQ.'GLU') THEN
           IF(NAMATM(I).EQ.'  N  ') THEN
             NCAR=NCAR+1
             LCARRS(NCAR)=NUMRES(I)
             NAMCAR(NCAR)=NAMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  OD1' .OR. NAMATM(I).EQ.'  OE1')
     $                  LCARO1(NCAR)=I
           IF(NAMATM(I).EQ.'  OD2' .OR. NAMATM(I).EQ.'  OE2')
     $                  LCARO2(NCAR)=I
         END IF
         IF(NAMATM(I).EQ.'  OXT')THEN
           NCAR=NCAR+1
           LCARRS(NCAR)=NUMRES(I)
           NAMCAR(NCAR)='C- '
           LCARO1(NCAR)=I
           DO K=1,50
             IF(NAMATM(I-K).EQ.'  O  ')THEN
               LCARO2(NCAR)=I-K
               GO TO 300
             ELSE IF(NAMATM(I+K).EQ.'  O  ')THEN
               LCARO2(NCAR)=I+K
               GO TO 300
             END IF
           END DO
         END IF
 300     CONTINUE
C
         IF((NAMRES(I).EQ.'LYS'.AND.NAMATM(I).EQ.'  NZ ') .OR.
     $                  (I.EQ.1.AND.NAMATM(I).EQ.'  N  ') )THEN
           NLYS=NLYS+1
           LLYSNZ(NLYS)=I
           LLYSRS(NLYS)=NUMRES(I)
         END IF
         IF(NAMRES(I).EQ.'TYR'.AND.NAMATM(I).EQ.'  OH ') THEN
           NTYR=NTYR+1
           LTYROH(NTYR)=I
           LTYRRS(NTYR)=NUMRES(I)
         END IF
         IF(NAMRES(I).EQ.'CYS'.AND.NAMATM(I).EQ.'  SG ') THEN
           NCYS=NCYS+1
           LCYSSG(NCYS)=I
           LCYSRS(NCYS)=NUMRES(I)
         END IF
         IF(NAMRES(I).EQ.'SER'.AND.NAMATM(I).EQ.'  OG ') THEN
           NSER=NSER+1
           LSEROH(NSER)=I
           LSERRS(NSER)=NUMRES(I)
         END IF
         IF(NAMRES(I).EQ.'THR'.AND.NAMATM(I).EQ.'  OG1') THEN
           NTHR=NTHR+1
           LTHROH(NTHR)=I
           LTHRRS(NTHR)=NUMRES(I)
         END IF
         IF(NAMRES(I).EQ.'GLN')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NGLN=NGLN+1
             LGLNRS(NGLN)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  NE2')LGLNNE(NGLN)=I
           IF(NAMATM(I).EQ.'  OE1')LGLNOE(NGLN)=I
           IF(NAMATM(I).EQ.'  CD ')LGLNCD(NGLN)=I
         END IF
         IF(NAMRES(I).EQ.'ASN')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NASN=NASN+1
             LASNRS(NASN)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  ND2')LASNND(NASN)=I
           IF(NAMATM(I).EQ.'  OD1')LASNOD(NASN)=I
           IF(NAMATM(I).EQ.'  CG ')LASNCG(NASN)=I
         END IF
         IF(NAMRES(I).EQ.'ARG')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NARG=NARG+1
             LARGRS(NARG)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  NE ')LARGN1(NARG)=I
           IF(NAMATM(I).EQ.'  NH2')LARGN2(NARG)=I
           IF(NAMATM(I).EQ.'  NH1')LARGN3(NARG)=I
           IF(NAMATM(I).EQ.'  CZ ')LARGCZ(NARG)=I
           IF(NAMATM(I).EQ.'  CD ')LARGCD(NARG)=I
         END IF
         IF(NAMRES(I).EQ.'HIS')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NHIS=NHIS+1
             LHISRS(NHIS)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  CG ')LHISCG(NHIS)=I
           IF(NAMATM(I).EQ.'  ND1')LHISND(NHIS)=I
           IF(NAMATM(I).EQ.'  CE1')LHISCE(NHIS)=I
           IF(NAMATM(I).EQ.'  NE2')LHISNE(NHIS)=I
           IF(NAMATM(I).EQ.'  CD2')LHISCD(NHIS)=I
         END IF
         IF(NAMRES(I).EQ.'TRP')THEN
           IF(NAMATM(I).EQ.'  N  ')THEN
             NTRP=NTRP+1
             LTRPRS(NTRP)=NUMRES(I)
           END IF
           IF(NAMATM(I).EQ.'  CD1')LTRPCD(NTRP)=I
           IF(NAMATM(I).EQ.'  NE1')LTRPNE(NTRP)=I
           IF(NAMATM(I).EQ.'  CE2')LTRPCE(NTRP)=I
         END IF
C
C LIGAND SITES
C
         IF(HEAD(I).EQ.'LG')THEN
           NLIGAND=NLIGAND+1
           LLIGAND(NLIGAND)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  Nam ')
     $    THEN
           NLIGNam=NLIGNam+1
           LLIGNam(NLIGNam)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  C2  ')
     $    THEN
           NLIGC2=NLIGC2+1
           LLIGC2(NLIGC2)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  Nar  ')
     $    THEN
           NLIGNar=NLIGNar+1
           LLIGNar(NLIGNar)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  O2  ')
     $    THEN
           NLIGO2=NLIGO2+1
           LLIGO2(NLIGO2)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  O3  ')
     $    THEN
           NLIGO3=NLIGO3+1
           LLIGO3(NLIGO3)=I
         END IF
         IF((HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  N3  ').OR.
     $      (HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  N4  '))
     $    THEN
           NLIGN3=NLIGN3+1
           LLIGN3(NLIGN3)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  C3  ')
     $    THEN
           NLIGC3=NLIGC3+1
           LLIGC3(NLIGC3)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  Oco ')
     $    THEN
           NLIGCAR=NLIGCAR+1
           LLIGCAR(NLIGCAR)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  Npl ')
     $    THEN
           NLIGNpl=NLIGNpl+1
           LLIGNpl(NLIGNpl)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  So2 ')
     $    THEN
           NLIGSo2=NLIGSo2+1
           LLIGSo2(NLIGSo2)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  N1  ')
     $    THEN
           NLIGN1=NLIGN1+1
           LLIGN1(NLIGN1)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  Cl  ')
     $    THEN
           NLIGCl=NLIGCl+1
           LLIGCl(NLIGCl)=I
         END IF
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  F   ')
     $    THEN
           NLIGF=NLIGF+1
           LLIGF(NLIGF)=I
         END IF
c dmr N4 to N30
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  N4  ')
     $    THEN
         NAMATM(I)='  N3  '
         END IF
c dmr Charged atoms
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  1P  ')
     $    THEN
           NLIGchr=NLIGchr+1
           LLIGchr(NLIGchr)=I
           LLIGchrval(NLIGchr)=1.0
       write(6,*)'Charged Atom ',nligchr,nbatm(i),
     $      lligchr(nligchr),lligchrval(nligchr) 
         END IF
c
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  1N  ')
     $    THEN
           NLIGchr=NLIGchr+1
           LLIGchr(NLIGchr)=I
           LLIGchrval(NLIGchr)=-1.0
       write(6,*)'Charged Atom ',nligchr,nbatm(i),
     $      lligchr(nligchr),lligchrval(nligchr) 
         END IF
c end dmr
c dcb Charged atoms (2charges)
c test Ribonuclease phosphate ligands
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  2P  ')
     $    THEN
           NLIGchr=NLIGchr+1
           LLIGchr(NLIGchr)=I
           LLIGchrval(NLIGchr)=+2.0
       write(6,*)'Charged Atom ',nligchr,nbatm(i),
     $      lligchr(nligchr),lligchrval(nligchr)
         END IF
c
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  2N  ')
     $    THEN
           NLIGchr=NLIGchr+1
           LLIGchr(NLIGchr)=I
           LLIGchrval(NLIGchr)=-2.0
       write(6,*)'Charged Atom ',nligchr,nbatm(i),
     $      lligchr(nligchr),lligchrval(nligchr)
         END IF
c end dcb
c dmr Charged atoms (3charges)
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  3P  ')
     $    THEN
           NLIGchr=NLIGchr+1
           LLIGchr(NLIGchr)=I
           LLIGchrval(NLIGchr)=+3.0
       write(6,*)'Charged Atom ',nligchr,nbatm(i),
     $      lligchr(nligchr),lligchrval(nligchr)
         END IF
c
         IF(HEAD(I).EQ.'LG'.AND.NAMATM(I).EQ.'  3N  ')
     $    THEN
           NLIGchr=NLIGchr+1
           LLIGchr(NLIGchr)=I
           LLIGchrval(NLIGchr)=-3.0
       write(6,*)'Charged Atom ',nligchr,nbatm(i),
     $      lligchr(nligchr),lligchrval(nligchr)
         END IF
c end dmr
       END DO
 200   CONTINUE
       CLOSE (10)
C
C
C
C
c dmr Check for missing side-chain atoms
       imiss=0
       do i=1,natom
c
c      if(namres(i).eq.'ASP' .or. namres(i).eq.'GLU'
c    $ .or. namres(i).eq.'HIS' .or. namres(i).eq.'TYR'
c    $ .or. namres(i).eq.'LYS' .or. namres(i).eq.'ARG')then
c      if(namatm(i).eq.'  CA ' .and. namatm(i+4).ne.'  CG ')then
c      write(6,1971)namres(i),numres(i),typch(i)
c      write(11,1971)namres(i),numres(i),typch(i)
c      imiss=imiss+1
c      end if
c      end if
c
       if(namres(i).eq.'ASP')then
       if(namatm(i).eq.'  CA ' .and. namatm(i+5).ne.'  OD1'
     $ .and. namatm(i+6).ne.'  OD2')then
       write(6,1971)namres(i),numres(i),typch(i)
       write(11,1971)namres(i),numres(i),typch(i)
       imiss=imiss+1
       end if
       end if
c
       if(namres(i).eq.'GLU')then
       if(namatm(i).eq.'  CA ' .and. namatm(i+6).ne.'  OE1'
     $ .and. namatm(i+7).ne.'  OE2')then
       write(6,1971)namres(i),numres(i),typch(i)
       write(11,1971)namres(i),numres(i),typch(i)
       imiss=imiss+1
       end if
       end if
c
       if(namres(i).eq.'HIS')then
       if(namatm(i).eq.'  CA ' .and. namatm(i+5).ne.'  ND1'
     $ .and. namatm(i+8).ne.'  NE2')then
       write(6,1971)namres(i),numres(i),typch(i)
       write(11,1971)namres(i),numres(i),typch(i)
       imiss=imiss+1
       end if
       end if
c
       if(namres(i).eq.'CYS')then
       if(namatm(i).eq.'  CA ' .and. namatm(i+4).ne.'  SG ')then
       write(6,1971)namres(i),numres(i),typch(i)
       write(11,1971)namres(i),numres(i),typch(i)
       imiss=imiss+1
       end if
       end if
c
       if(namres(i).eq.'TYR')then
       if(namatm(i).eq.'  CA ' .and. namatm(i+10).ne.'  OH ')then
       write(6,1971)namres(i),numres(i),typch(i)
       write(11,1971)namres(i),numres(i),typch(i)
       imiss=imiss+1
       end if
       end if
c
       if(namres(i).eq.'LYS')then
       if(namatm(i).eq.'  CA ' .and. namatm(i+7).ne.'  NZ ')then
       write(6,1971)namres(i),numres(i),typch(i)
       write(11,1971)namres(i),numres(i),typch(i)
       imiss=imiss+1
       end if
       end if
c
       if(namres(i).eq.'ARG')then
       if(namatm(i).eq.'  CA ' .and. namatm(i+6).ne.'  NE '
     $ .and. namatm(i+8).ne.'  NH1' .and. namatm(i+9).ne.'  NH2')then
       write(6,1971)namres(i),numres(i),typch(i)
       write(11,1971)namres(i),numres(i),typch(i)
       imiss=imiss+1
       end if
       end if
c
       end do
       if(imiss.gt.0)then
       write(6,1972)imiss
       write(11,1972)imiss
       end if
 1971  format(' Warning: Missing side-chain atoms for ',a3,i4,a2)
 1972  format(' Missing side-chain atoms for ',i4,' ionizable residues')
c
c dmr if ligand bonded to Asp or Glu Oco atom (OD/E1 or OD/E2), eg. 1BVV.pdb
         distb=1.9
       DO ICAR=1,NCAR
         XO1=X(LCARO1(ICAR))
         YO1=Y(LCARO1(ICAR))
         ZO1=Z(LCARO1(ICAR))
         XO2=X(LCARO2(ICAR))
         YO2=Y(LCARO2(ICAR))
         ZO2=Z(LCARO2(ICAR))
        do i=1,nligand
         xl=x(lligand(i))
         yl=y(lligand(i))
         zl=z(lligand(i))
              DIS1=SQRT((Xl-Xo1)**2+(Yl-Yo1)**2+(Zl-Zo1)**2)
              DIS2=SQRT((Xl-Xo2)**2+(Yl-Yo2)**2+(Zl-Zo2)**2)
 
            if(dis1.lt.distb)then

         write(6,*)namcar(icar),'',lcarrs(icar),typch(lcaro1(icar)),
     $   namatm(lcaro1(icar)),' is bound to ligand'
              rmicar=icar

      if (namcar(icar).eq.'ASP')then
              nligo3=nligo3+1
           lligo3(nligo3)=lcaro1(icar)
           NAMres(lligo3(nligo3))=namres(lligand(i))
           NbATM(lligo3(nligo3))=' OD1'
           NAMATM(lligo3(nligo3))='  O3  '

              nligo2=nligo2+1
           lligo2(nligo2)=lcaro2(icar)
           NAMres(lligo2(nligo2))=namres(lligand(i))
           NbATM(lligo2(nligo2))=' OD2'
           NAMATM(lligo2(nligo2))='  O2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)-1
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' CG '
           NAMATM(lligand(nligand))='  C2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OD1'
           NAMATM(lligand(nligand))='  O3  '

              nligand=nligand+1
           lligand(nligand)=lcaro2(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OD2'
           NAMATM(lligand(nligand))='  O2  '

      else if (namcar(icar).eq.'GLU')then
              nligo3=nligo3+1
           lligo3(nligo3)=lcaro1(icar)
           NAMres(lligo3(nligo3))=namres(lligand(i))
           NbATM(lligo3(nligo3))=' OE1'
           NAMATM(lligo3(nligo3))='  O3  '

              nligo2=nligo2+1
           lligo2(nligo2)=lcaro2(icar)
           NAMres(lligo2(nligo2))=namres(lligand(i))
           NbATM(lligo2(nligo2))=' OE2'
           NAMATM(lligo2(nligo2))='  O2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)-1
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' CD '
           NAMATM(lligand(nligand))='  C2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OE1'
           NAMATM(lligand(nligand))='  O3  '

              nligand=nligand+1
           lligand(nligand)=lcaro2(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OE2'
           NAMATM(lligand(nligand))='  O2  '
      end if

       DO ICARr=rmicar,NCAR-1
       if(icarr.ge.rmicar)then
         NAMcar(icarr)=NAMcar(icarr+1)
         lcarrs(icarr)=lcarrs(icarr+1)
         lcaro1(icarr)=lcaro1(icarr+1)
         lcaro2(icarr)=lcaro2(icarr+1)
         NAMATM(lcaro1(icarr))=NAMATM(lcaro1(icarr+1))
         NAMATM(lcaro2(icarr))=NAMATM(lcaro2(icarr+1))
       end if
       end do
       ncar=ncar-1

            end if

            if(dis2.lt.distb)then

         write(6,*)namcar(icar),'',lcarrs(icar),typch(lcaro2(icar)),
     $   namatm(lcaro2(icar)),' is bound to ligand'
              rmicar=icar

      if (namcar(icar).eq.'ASP')then
              nligo3=nligo3+1
           lligo3(nligo3)=lcaro2(icar)
           NAMres(lligo3(nligo3))=namres(lligand(i))
           NbATM(lligo3(nligo3))=' OD2'
           NAMATM(lligo3(nligo3))='  O3  '

              nligo2=nligo2+1
           lligo2(nligo2)=lcaro1(icar)
           NAMres(lligo2(nligo2))=namres(lligand(i))
           NbATM(lligo2(nligo2))=' OD1'
           NAMATM(lligo2(nligo2))='  O2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)-1
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' CG '
           NAMATM(lligand(nligand))='  C2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OD1'
           NAMATM(lligand(nligand))='  O2  '

              nligand=nligand+1
           lligand(nligand)=lcaro2(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OD2'
           NAMATM(lligand(nligand))='  O3  '

      else if (namcar(icar).eq.'GLU')then
              nligo3=nligo3+1
           lligo3(nligo3)=lcaro2(icar)
           NAMres(lligo3(nligo3))=namres(lligand(i))
           NbATM(lligo3(nligo3))=' OE2'
           NAMATM(lligo3(nligo3))='  O3  '

              nligo2=nligo2+1
           lligo2(nligo2)=lcaro1(icar)
           NAMres(lligo2(nligo2))=namres(lligand(i))
           NbATM(lligo2(nligo2))=' OE1'
           NAMATM(lligo2(nligo2))='  O2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)-1
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' CD '
           NAMATM(lligand(nligand))='  C2  '

              nligand=nligand+1
           lligand(nligand)=lcaro1(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OE1'
           NAMATM(lligand(nligand))='  O2  '

              nligand=nligand+1
           lligand(nligand)=lcaro2(icar)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OE2'
           NAMATM(lligand(nligand))='  O3  '
      end if

       DO ICARr=rmicar,NCAR-1
       if(icarr.ge.rmicar)then
         NAMcar(icarr)=NAMcar(icarr+1)
         lcarrs(icarr)=lcarrs(icarr+1)
         lcaro1(icarr)=lcaro1(icarr+1)
         lcaro2(icarr)=lcaro2(icarr+1)
         NAMATM(lcaro1(icarr))=NAMATM(lcaro1(icarr+1))
         NAMATM(lcaro2(icarr))=NAMATM(lcaro2(icarr+1))
       end if
       end do
       ncar=ncar-1
c do if both OE1 and OE2 bound to ligand?
      end if
      end do
      end do


c if Ser OG is bonded to ligand, eg. 1GGD.pdb
         distb=1.9
       DO Iser=1,Nser
         XOg=X(LSEROH(Iser))
         YOg=Y(LSEROH(Iser))
         ZOg=Z(LSEROH(Iser))
        do i=1,nligand
         xl=x(lligand(i))
         yl=y(lligand(i))
         zl=z(lligand(i))
              DIS=SQRT((Xl-Xog)**2+(Yl-Yog)**2+(Zl-Zog)**2)
             if(dis.lt.distb)then
         write(6,*)namres(lserrs(iser)),'',lserrs(iser),
     $   typch(lseroh(iser)),NAMATM(lseroh(iser)),' is bound to ligand'
              rmiser=iser

              nligo3=nligo3+1
           lligo3(nligo3)=lseroh(iser)
           NAMres(lligo3(nligo3))=namres(lligand(i))
           NbATM(lligo3(nligo3))=' OG '
           NAMATM(lligo3(nligo3))='  O3  '

              nligand=nligand+1
           lligand(nligand)=lseroh(iser)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' OG '
           NAMATM(lligand(nligand))='  O3  '

              nligand=nligand+1
           lligand(nligand)=lseroh(iser)-1
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' CB '
           NAMATM(lligand(nligand))='  C3  '

       DO Iserr=rmiser,Nser-1
       if(iserr.ge.rmiser)then
         lserrs(iserr)=lserrs(iserr+1)
         lseroh(iserr)=lseroh(iserr+1)
         NAMATM(lseroh(iserr))=NAMATM(lseroh(iserr+1))
       end if
       end do
       nser=nser-1

      end if
      end do
      end do


c if Cys SG is bonded to ligand, eg. 2GB0.pdb
         nligs3=0
         distb=2.1
       DO Icys=1,Ncys
         XSg=X(LCysSg(Icys))
         YSg=Y(LCysSg(Icys))
         ZSg=Z(LCysSg(Icys))
        do i=1,nligand
         xl=x(lligand(i))
         yl=y(lligand(i))
         zl=z(lligand(i))
              DIS=SQRT((Xl-XSg)**2+(Yl-YSg)**2+(Zl-ZSg)**2)
             if(dis.lt.distb)then
         write(6,*)namres(lcyssg(icys)),'',lcysrs(icys),
     $   typch(lcyssg(icys)),NAMATM(lcyssg(icys)),' is bound to ligand'
              rmicys=icys
              rmicysn=lcysrs(icys)

              nligs3=nligs3+1
           lligs3(nligs3)=lcyssg(icys)
           NAMres(lligs3(nligs3))=namres(lligand(i))
           NbATM(lligs3(nligs3))=' SG '
           NAMATM(lligs3(nligs3))='  S3 '

              nligand=nligand+1
           lligand(nligand)=lcyssg(icys)
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' SG '
           NAMATM(lligand(nligand))='  S3 '

       DO j=1,Natom
       if(numres(j).eq.rmicysn.and.namatm(j).eq.'  CB ' )then
              nligand=nligand+1
           lligand(nligand)=j
           NAMres(lligand(nligand))=namres(lligand(i))
           NbATM(lligand(nligand))=' CB '
           NAMATM(lligand(nligand))='  C3 '
       end if
       end do

       DO Icysr=rmicys,Ncys
       if(icys.lt.ncys)then
         lcysrs(icysr)=lcysrs(icysr+1)
         lcyssg(icysr)=lcyssg(icysr+1)
         NAMATM(lcyssg(icysr))=NAMATM(lcyssg(icysr+1))
        else
         lcysrs(icysr)=0
         lcyssg(icysr)=0
         NAMATM(lcyssg(icysr))='  X  '
       end if
       end do
       ncys=ncys-1

      end if
      end do
      end do

       write(6,*) nligand,' ligand atoms:'
      do i=1,nligand
       write(6,'(a3,i5,a5,a4,f8.3,f8.3,f8.3)')
     $ NAMres(lligand(i)),lligand(i),
     $ NAMATM(lligand(i)),NbATM(lligand(i)),
     $ X(Lligand(i)),Y(Lligand(i)),Z(Lligand(i))
      end do
c dmr
C      **********************
C      STEP 2. ASSIGN H ATOMS
C      **********************
C
C      -- DETERMINE BACKBONE PROTON POSITIONS --
C
       NUMPRT(1)=NUMRES(1)
       NAMPRT(1)=NAMRES(1)
       DO IATOM = 2, NATOM
         IF(NAMATM(IATOM).EQ.'  C  ')THEN
           XC=X(IATOM)
           YC=Y(IATOM)
           ZC=Z(IATOM)
         END IF
         IF(NAMATM(IATOM).EQ.'  O  ')THEN
           XO=X(IATOM)
           YO=Y(IATOM)
           ZO=Z(IATOM)
         END IF
         IF(NAMATM(IATOM).EQ.'  N  '.AND.
     $      NUMRES(IATOM).GT.NUMRES(IATOM-1).AND.
     $      NAMRES(IATOM).NE.'PRO')THEN
           XN=X(IATOM)
           YN=Y(IATOM)
           ZN=Z(IATOM)
           VECNRM=SQRT((XC-XO)**2+(YC-YO)**2+(ZC-ZO)**2)
           XVEC=(XC-XO)/VECNRM
           YVEC=(YC-YO)/VECNRM
           ZVEC=(ZC-ZO)/VECNRM
C          -- NOTE N-H BOND LENGTH = 1 ANGSTROM 
           XH=XN+XVEC
           YH=YN+YVEC
           ZH=ZN+ZVEC
           NPRTON=NPRTON+1
           XPRTON(NPRTON)=XH
           YPRTON(NPRTON)=YH
           ZPRTON(NPRTON)=ZH
           XNITRN(NPRTON)=XN
           YNITRN(NPRTON)=YN
           ZNITRN(NPRTON)=ZN
           XOXYGN(NPRTON)=XO
           YOXYGN(NPRTON)=YO
           ZOXYGN(NPRTON)=ZO
           XCARBN(NPRTON)=XC
           YCARBN(NPRTON)=YC
           ZCARBN(NPRTON)=ZC
           NUMPRT(NPRTON)=NUMRES(IATOM)
           NAMPRT(NPRTON)=NAMRES(IATOM)
         END IF
       END DO
C      WRITE(6,*)NPRTON,'AMIDE PROTONS ASSIGNED'
C
C      -- DETERMINE GLN PROTON POSITIONS --
C
       NGLNP=0
       DO IGLN = 1, NGLN
         XC=X(LGLNCD(IGLN))
         YC=Y(LGLNCD(IGLN))
         ZC=Z(LGLNCD(IGLN))
         XO=X(LGLNOE(IGLN))
         YO=Y(LGLNOE(IGLN))
         ZO=Z(LGLNOE(IGLN))
         XN=X(LGLNNE(IGLN))
         YN=Y(LGLNNE(IGLN))
         ZN=Z(LGLNNE(IGLN))
         VECNRM=SQRT((XC-XO)**2+(YC-YO)**2+(ZC-ZO)**2)
         XVEC=(XC-XO)/VECNRM
         YVEC=(YC-YO)/VECNRM
         ZVEC=(ZC-ZO)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN+XVEC
         YH=YN+YVEC
         ZH=ZN+ZVEC
         NGLNP=NGLNP+1
         XGLNP1(IGLN)=XH
         YGLNP1(IGLN)=YH
         ZGLNP1(IGLN)=ZH
C
         XON=(XO+XN)/2.0
         YON=(YO+YN)/2.0
         ZON=(ZO+ZN)/2.0
         VECNRM=SQRT((XC-XON)**2+(YC-YON)**2+(ZC-ZON)**2)
         XVEC=(XC-XON)/VECNRM
         YVEC=(YC-YON)/VECNRM
         ZVEC=(ZC-ZON)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN-XVEC
         YH=YN-YVEC
         ZH=ZN-ZVEC
         NGLNP=NGLNP+1
         XGLNP2(IGLN)=XH
         YGLNP2(IGLN)=YH
         ZGLNP2(IGLN)=ZH
       END DO
C      WRITE(6,*)NGLNP,'GLN PROTONS ASSIGNED'
C
C      -- DETERMINE ASN PROTON POSITIONS --
C
       NASNP=0
       DO IASN = 1, NASN
         XC=X(LASNCG(IASN))
         YC=Y(LASNCG(IASN))
         ZC=Z(LASNCG(IASN))
         XO=X(LASNOD(IASN))
         YO=Y(LASNOD(IASN))
         ZO=Z(LASNOD(IASN))
         XN=X(LASNND(IASN))
         YN=Y(LASNND(IASN))
         ZN=Z(LASNND(IASN))
         VECNRM=SQRT((XC-XO)**2+(YC-YO)**2+(ZC-ZO)**2)
         XVEC=(XC-XO)/VECNRM
         YVEC=(YC-YO)/VECNRM
         ZVEC=(ZC-ZO)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN+XVEC
         YH=YN+YVEC
         ZH=ZN+ZVEC
         NASNP=NASNP+1
         XASNP1(IASN)=XH
         YASNP1(IASN)=YH
         ZASNP1(IASN)=ZH
C
         XON=(XO+XN)/2.0
         YON=(YO+YN)/2.0
         ZON=(ZO+ZN)/2.0
         VECNRM=SQRT((XC-XON)**2+(YC-YON)**2+(ZC-ZON)**2)
         XVEC=(XC-XON)/VECNRM
         YVEC=(YC-YON)/VECNRM
         ZVEC=(ZC-ZON)/VECNRM
C        -- NOTE N-H BOND LENGTH = 1 ANGSTROM
         XH=XN-XVEC
         YH=YN-YVEC
         ZH=ZN-ZVEC
         NASNP=NASNP+1
         XASNP2(IASN)=XH
         YASNP2(IASN)=YH
         ZASNP2(IASN)=ZH
       END DO
C      WRITE(6,*)NASNP,'ASN PROTONS ASSIGNED'
C
C
C      -- DETERMINE ARG PROTON POSITIONS --
C
       NARGP=0
       DO IARG = 1, NARG
         XCD=X(LARGCD(IARG))
         YCD=Y(LARGCD(IARG))
         ZCD=Z(LARGCD(IARG))
         XCZ=X(LARGCZ(IARG))
         YCZ=Y(LARGCZ(IARG))
         ZCZ=Z(LARGCZ(IARG))
         XN1=X(LARGN1(IARG))
         YN1=Y(LARGN1(IARG))
         ZN1=Z(LARGN1(IARG))
         XN2=X(LARGN2(IARG))
         YN2=Y(LARGN2(IARG))
         ZN2=Z(LARGN2(IARG))
         XN3=X(LARGN3(IARG))
         YN3=Y(LARGN3(IARG))
         ZN3=Z(LARGN3(IARG))
C
         XO=(XCD+XCZ)/2.0
         YO=(YCD+YCZ)/2.0
         ZO=(ZCD+ZCZ)/2.0
         VECNRM=SQRT((XN1-XO)**2+(YN1-YO)**2+(ZN1-ZO)**2)
         XVEC=(XN1-XO)/VECNRM
         YVEC=(YN1-YO)/VECNRM
         ZVEC=(ZN1-ZO)/VECNRM
         NARGP=NARGP+1
         XARGP1(IARG)=XN1+XVEC
         YARGP1(IARG)=YN1+YVEC
         ZARGP1(IARG)=ZN1+ZVEC
         NARGP=NARGP+1
         XARGP2(IARG)=XN2+XVEC
         YARGP2(IARG)=YN2+YVEC
         ZARGP2(IARG)=ZN2+ZVEC
C
         XO=XN1
         YO=YN1
         ZO=ZN1
         VECNRM=SQRT((XCZ-XO)**2+(YCZ-YO)**2+(ZCZ-ZO)**2)
         XVEC=(XCZ-XO)/VECNRM
         YVEC=(YCZ-YO)/VECNRM
         ZVEC=(ZCZ-ZO)/VECNRM
         NARGP=NARGP+1
         XARGP3(IARG)=XN2+XVEC
         YARGP3(IARG)=YN2+YVEC
         ZARGP3(IARG)=ZN2+ZVEC
         NARGP=NARGP+1
         XARGP4(IARG)=XN3+XVEC
         YARGP4(IARG)=YN3+YVEC
         ZARGP4(IARG)=ZN3+ZVEC
C
         XO=XN1
         YO=YN1
         ZO=ZN1
         VECNRM=SQRT((XCD-XO)**2+(YCD-YO)**2+(ZCD-ZO)**2)
         XVEC=(XCD-XO)/VECNRM
         YVEC=(YCD-YO)/VECNRM
         ZVEC=(ZCD-ZO)/VECNRM
         NARGP=NARGP+1
         XARGP5(IARG)=XN3+XVEC
         YARGP5(IARG)=YN3+YVEC
         ZARGP5(IARG)=ZN3+ZVEC
C
       END DO
C      WRITE(6,*)NARGP,'ARG PROTONS ASSIGNED'
C
C
C      -- DETERMINE TRP PROTON POSITIONS --
C
       NTRPP=0
       DO ITRP = 1, NTRP
         XNE=X(LTRPNE(ITRP))
         YNE=Y(LTRPNE(ITRP))
         ZNE=Z(LTRPNE(ITRP))
         XCE=X(LTRPCE(ITRP))
         YCE=Y(LTRPCE(ITRP))
         ZCE=Z(LTRPCE(ITRP))
         XCD=X(LTRPCD(ITRP))
         YCD=Y(LTRPCD(ITRP))
         ZCD=Z(LTRPCD(ITRP))
C
         XO=(XCD+XCE)/2.0
         YO=(YCD+YCE)/2.0
         ZO=(ZCD+ZCE)/2.0
         VECNRM=SQRT((XNE-XO)**2+(YNE-YO)**2+(ZNE-ZO)**2)
         XVEC=(XNE-XO)/VECNRM
         YVEC=(YNE-YO)/VECNRM
         ZVEC=(ZNE-ZO)/VECNRM
         NTRPP=NTRPP+1
         XTRPP1(ITRP)=XNE+XVEC
         YTRPP1(ITRP)=YNE+YVEC
         ZTRPP1(ITRP)=ZNE+ZVEC
       END DO
C      WRITE(6,*)NTRPP,'TRP PROTONS ASSIGNED'
C
C
C      -- DETERMINE HIS PROTON POSITIONS --
C
       NHISP=0
       DO IHIS = 1, NHIS
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XO=(XCG+XCE)/2.0
         YO=(YCG+YCE)/2.0
         ZO=(ZCG+ZCE)/2.0
         VECNRM=SQRT((XND-XO)**2+(YND-YO)**2+(ZND-ZO)**2)
         XVEC=(XND-XO)/VECNRM
         YVEC=(YND-YO)/VECNRM
         ZVEC=(ZND-ZO)/VECNRM
         NHISP=NHISP+1
         XHISP1(IHIS)=XND+XVEC
         YHISP1(IHIS)=YND+YVEC
         ZHISP1(IHIS)=ZND+ZVEC
         XO=(XCD+XCE)/2.0
         YO=(YCD+YCE)/2.0
         ZO=(ZCD+ZCE)/2.0
         VECNRM=SQRT((XNE-XO)**2+(YNE-YO)**2+(ZNE-ZO)**2)
         XVEC=(XNE-XO)/VECNRM
         YVEC=(YNE-YO)/VECNRM
         ZVEC=(ZNE-ZO)/VECNRM
         NHISP=NHISP+1
         XHISP2(IHIS)=XNE+XVEC
         YHISP2(IHIS)=YNE+YVEC
         ZHISP2(IHIS)=ZNE+ZVEC
       END DO
C      WRITE(6,*)NHISP,'HIS PROTONS ASSIGNED'
C
C
C
C
C      -- DETERMINE AMIDE PROTON POSITIONS (LIGAND)--
C
       NHLIG=0
       DO INHLIG=1,NLIGNam
         XN=X(LLIGNam(INHLIG))
         YN=Y(LLIGNam(INHLIG))
         ZN=Z(LLIGNam(INHLIG))
         NB_neighbours=0
         DO I=1,NLIGAND
          IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGNam(INHLIG))) THEN
           Xb=X(LLIGAND(I))
           Yb=Y(LLIGAND(I))
           Zb=Z(LLIGAND(I))
           DISn=SQRT((Xb-XN)**2+(Yb-YN)**2+(Zb-ZN)**2)
           IF (DISn.LT.2) THEN 
              NB_neighbours=NB_neighbours+1
           ENDIF
          ENDIF
         END DO
         IF (NB_neighbours.EQ.3) THEN
          NAMATM(LLIGNam(INHLIG))='  C3 '
         ENDIF
         IF (NB_neighbours.EQ.2) THEN
          DIST=1000.0
          DO J=1,NLIGC2
            XC=X(LLIGC2(J))
            YC=Y(LLIGC2(J))
            ZC=Z(LLIGC2(J))
            DIS=SQRT((XC-XN)**2+(YC-YN)**2+(ZC-ZN)**2)
            DIST=MIN(DIST,DIS)
             IF(DIS.EQ.DIST) THEN
              XCam=XC
              YCam=YC
              ZCam=ZC
             ENDIF
          END DO
          DO K=1,NLIGO2
           XO=X(LLIGO2(K))
           YO=Y(LLIGO2(K))
           ZO=Z(LLIGO2(K))
            DIS=SQRT((XO-XCam)**2+(YO-YCam)**2+(ZO-ZCam)**2)
            DIST=MIN(DIST,DIS)
             IF(DIS.EQ.DIST) THEN
              XOcar=XO
              YOcar=YO
              ZOcar=ZO
             ENDIF
          END DO
          VECNRM=SQRT((XCam-XOcar)**2+(YCam-YOcar)**2+(ZCam-ZOcar)**2)
          XVEC=(XCam-XOcar)/VECNRM
          YVEC=(YCam-YOcar)/VECNRM
          ZVEC=(ZCam-ZOcar)/VECNRM
          NHLIG=NHLIG+1
          XHLIG(INHLIG)=XN+XVEC
          YHLIG(INHLIG)=YN+YVEC
          ZHLIG(INHLIG)=ZN+ZVEC
         ENDIF 
       END DO
C
C
C      -- DETERMINE PROTON POSITION FOR Nar ATOMS --
C
       NgLIG=0
       DO INarLIG=1,NLIGNar
         XN=X(LLIGNar(INarLIG))
         YN=Y(LLIGNar(INarLIG))
         ZN=Z(LLIGNar(INarLIG))
         DIST=1000.0
         DO J=1,NLIGAND
          IF(NAMATM(LLIGAND(J)).EQ.'  Car ')THEN
           XCar=X(LLIGAND(J))
           YCar=Y(LLIGAND(J))
           ZCar=Z(LLIGAND(J))
           DIS=SQRT((XCar-XN)**2+(YCar-YN)**2+(ZCar-ZN)**2)
           DIST=MIN(DIST,DIS)
            IF(DIS.EQ.DIST) THEN
             XC1=XCar
             YC1=YCar
             ZC1=ZCar
            ENDIF
          ENDIF
         END DO
C
C XC1,YC1,ZC1 coordinates of the closest Car atom
C
         DIST=1000.0
         DO J=1,NLIGAND
          IF(NAMATM(LLIGAND(J)).EQ.'  Car ' .AND.
     $        X(LLIGAND(J)).NE.XC1)THEN
           XCar=X(LLIGAND(J))
           YCar=Y(LLIGAND(J))
           ZCar=Z(LLIGAND(J))
           DIS=SQRT((XCar-XN)**2+(YCar-YN)**2+(ZCar-ZN)**2)
           DIST=MIN(DIST,DIS)
            IF(DIS.EQ.DIST) THEN
             XC2=XCar
             YC2=YCar
             ZC2=ZCar
            ENDIF
          ENDIF
         END DO
C
C XC2,YC2,ZC2 coordinates of the 2nd closest Car atom
C
C coordinates of the middle between C1 and C2
         Xmid=(XC1+XC2)/2
         Ymid=(YC1+YC2)/2
         Zmid=(ZC1+ZC2)/2
C Direction vector midC1C2->Nar
         VECNRM=SQRT((XN-Xmid)**2+(YN-Ymid)**2+(ZN-Zmid)**2)
         XVEC=(XN-Xmid)/VECNRM
         YVEC=(YN-Ymid)/VECNRM
         ZVEC=(ZN-Zmid)/VECNRM
         NgLIG=NgLIG+1
         XgLIG(NgLIG)=XN+XVEC
         YgLIG(NgLIG)=YN+YVEC
         ZgLIG(NgLIG)=ZN+ZVEC
       END DO
C
C
C      -- DETERMINE Npl PROTON POSITIONS --
C
       NHLIGNp1=0
       NHLIGNp2=0
       DO IN=1,NLIGNpl
         XLN(IN)=0
         YLN(IN)=0
         ZLN(IN)=0
         XLN1(IN)=0
         YLN1(IN)=0
         ZLN1(IN)=0
         XLN2(IN)=0
         YLN2(IN)=0
         ZLN2(IN)=0
         NAMLN1(IN)='00000'
         NAMLN2(IN)='00000'
         NAMLN(IN)='00000'
         XP1Np1(IN)=0
         YP1Np1(IN)=0
         ZP1Np1(IN)=0
         XP2Np1(IN)=0
         YP2Np1(IN)=0
         ZP2Np1(IN)=0
         XPNp2(IN)=0
         YPNp2(IN)=0
         ZPNp2(IN)=0
         PK1LGNp1(IN)=0
       END DO
       DO IN=1,NLIGNpl
         XN=X(LLIGNpl(IN))
         YN=Y(LLIGNpl(IN))
         ZN=Z(LLIGNpl(IN))
         NB_Npl=0
         DO I=1,NLIGAND
          IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGNpl(IN))) THEN
           Xt=X(LLIGAND(I))
           Yt=Y(LLIGAND(I))
           Zt=Z(LLIGAND(I))
           DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
           IF (DIS.LT.2) THEN
              NB_Npl=NB_Npl+1
           ENDIF
          ENDIF
         END DO
C
         IF (NB_Npl.EQ.0) THEN
           NAMATM(LLIGNpl(IN))='  Np0'
         ENDIF
C
         IF (NB_Npl.EQ.1) THEN
           NAMATM(LLIGNpl(IN))='  Np1'
           DO I=1,NLIGAND
            IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGNpl(IN))) THEN
             Xb=X(LLIGAND(I))
             Yb=Y(LLIGAND(I))
             Zb=Z(LLIGAND(I))
             NAMt=NAMATM(LLIGAND(I))
             DISn=SQRT((Xb-XN)**2+(Yb-YN)**2+(Zb-ZN)**2)
             IF (DISn.LT.2) THEN
               XLN(IN)=Xb
               YLN(IN)=Yb
               ZLN(IN)=Zb
               NAMLN(IN)=NAMt
             ENDIF
            ENDIF
           END DO
           DIST=1000.0
           DO J=1,NLIGAND
            IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGNpl(IN))
     $          .AND.X(LLIGAND(J)).NE.XLN(IN)) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              NAMt=NAMATM(LLIGAND(J))
              DIS=SQRT((Xt-XLN(IN))**2+(Yt-YLN(IN))**2+(Zt-ZLN(IN))**2)
              DIST=MIN(DIST,DIS)
              IF(DIS.EQ.DIST) THEN
                XLN1(IN)=Xt
                YLN1(IN)=Yt
                ZLN1(IN)=Zt
                NAMLN1(IN)=NAMt
              ENDIF
            ENDIF
           END DO
C
C -- X1,Y1,Z1 coordinates of the closest neighbour --
C
C
           DIST=1000.0
           DO J=1,NLIGAND
             IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGNpl(IN))
     $          .AND.X(LLIGAND(J)).NE.XLN(IN)
     $          .AND.X(LLIGAND(J)).NE.XLN1(IN)) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              NAMt=NAMATM(LLIGAND(J))
              DIS=SQRT((Xt-XLN(IN))**2+(Yt-YLN(IN))**2+(Zt-ZLN(IN))**2)
              DIST=MIN(DIST,DIS)
               IF(DIS.EQ.DIST) THEN
                XLN2(IN)=Xt
                YLN2(IN)=Yt
                ZLN2(IN)=Zt
                NAMLN2(IN)=NAMt
               ENDIF
             ENDIF
           END DO
C
C
C -- X2,Y2,Z2 coordinates of the 2nd closest neighbour --
C
C
           VECNRM1=SQRT((XLN(IN)-XLN1(IN))**2+(YLN(IN)-YLN1(IN))**2
     $                 +(ZLN(IN)-ZLN1(IN))**2)
           XVEC1=(XLN(IN)-XLN1(IN))/VECNRM1
           YVEC1=(YLN(IN)-YLN1(IN))/VECNRM1
           ZVEC1=(ZLN(IN)-ZLN1(IN))/VECNRM1
C
           VECNRM2=SQRT((XLN(IN)-XLN2(IN))**2+(YLN(IN)-YLN2(IN))**2
     $                 +(ZLN(IN)-ZLN2(IN))**2)
           XVEC2=(XLN(IN)-XLN2(IN))/VECNRM2
           YVEC2=(YLN(IN)-YLN2(IN))/VECNRM2
           ZVEC2=(ZLN(IN)-ZLN2(IN))/VECNRM2
C
           IF((NAMLN1(IN).NE.'  Npl'.OR.NAMLN1(IN).NE.'  Np1')
     $        .AND.(NAMLN2(IN).EQ.'  Npl'.OR.NAMLN2(IN).EQ.'  Np1')
     $                                                       ) THEN
             NHLIGNp1=NHLIGNp1+2
             XP1Np1(IN)=XN+XVEC1
             YP1Np1(IN)=YN+YVEC1
             ZP1Np1(IN)=ZN+ZVEC1
             XP2Np1(IN)=XN+XVEC2
             YP2Np1(IN)=YN+YVEC2
             ZP2Np1(IN)=ZN+ZVEC2
             PK1LGNp1(IN)=PKALIG(LLIGNpl(IN))
           ENDIF

           IF((NAMLN1(IN).EQ.'  Npl'.OR.NAMLN1(IN).EQ.'  Np1')
     $        .AND.(NAMLN2(IN).NE.'  Npl'.OR.NAMLN2(IN).NE.'  Np1')
     $                                                       ) THEN
             NHLIGNp1=NHLIGNp1+2
             XP1Np1(IN)=XN+XVEC2
             YP1Np1(IN)=YN+YVEC2
             ZP1Np1(IN)=ZN+ZVEC2
             XP2Np1(IN)=XN+XVEC1
             YP2Np1(IN)=YN+YVEC1
             ZP2Np1(IN)=ZN+ZVEC1
             PK1LGNp1(IN)=PKALIG(LLIGNpl(IN))
           ENDIF
C
           IF(NAMLN1(IN).NE.'  Npl'.AND.NAMLN1(IN).NE.'  Np1'
     $        .AND.NAMLN2(IN).NE.'  Npl'.AND.NAMLN2(IN).NE.'  Np1'
     $                                                       ) THEN
             NHLIGNp1=NHLIGNp1+2
             XP1Np1(IN)=XN+XVEC1
             YP1Np1(IN)=YN+YVEC1
             ZP1Np1(IN)=ZN+ZVEC1
             XP2Np1(IN)=XN+XVEC2
             YP2Np1(IN)=YN+YVEC2
             ZP2Np1(IN)=ZN+ZVEC2
           ENDIF
         ENDIF
C
         IF (NB_Npl.EQ.2) THEN
           NAMATM(LLIGNpl(IN))='  Np2'
           DIST=1000.0
           DO J=1,NLIGAND
            IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGNpl(IN))) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
              DIST=MIN(DIST,DIS)
              IF(DIS.EQ.DIST) THEN
                X1=Xt
                Y1=Yt
                Z1=Zt
              ENDIF
            ENDIF
           END DO
C
C -- X1,Y1,Z1 coordinates of the closest neighbour --
C
           DIST=1000.0
           DO J=1,NLIGAND
             IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGNpl(IN))
     $          .AND.X(LLIGAND(J)).NE.X1) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
              DIST=MIN(DIST,DIS)
               IF(DIS.EQ.DIST) THEN
                X2=Xt
                Y2=Yt
                Z2=Zt
               ENDIF
             ENDIF
           END DO
C
C
C -- X2,Y2,Z2 coordinates of the 2nd closest neighbour --
C
           XM=(X1+X2)/2.0
           YM=(Y1+Y2)/2.0
           ZM=(Z1+Z2)/2.0
C
C -- M middle between the two neighbours --
C -- M->N direction ot assign the proton --
C 
         VECNRM=SQRT((XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
         XVEC=(XN-XM)/VECNRM
         YVEC=(YN-YM)/VECNRM
         ZVEC=(ZN-ZM)/VECNRM
         NHLIGNp2=NHLIGNp2+1
         XPNp2(IN)=XN+XVEC
         YPNp2(IN)=YN+YVEC
         ZPNp2(IN)=ZN+ZVEC
         ENDIF
C
       END DO
C
C
C identification of guadininium type group (1K1O.pdb) and
C pair of Np1 nitrogen group (1K1I.pdb)
C based on the nature of the nitrogen linked (ie within a distance less
C than 2.0 angstr�ms) to the C2 carbon
C
C
       DO IC2=1,NLIGC2 
         TEST1=.FALSE.
         TEST2=.FALSE.
         TEST3=.FALSE.
         XC=X(LLIGC2(IC2))
         YC=Y(LLIGC2(IC2))
         ZC=Z(LLIGC2(IC2))
         DO IN=1,NLIGNpl
           Xt=X(LLIGNpl(IN))
           Yt=Y(LLIGNpl(IN))
           Zt=Z(LLIGNpl(IN))
           DIS=SQRT((Xt-XC)**2+(Yt-YC)**2+(Zt-ZC)**2)
           IF(NAMATM(LLIGNpl(IN)).EQ.'  Np2'
     $        .AND.DIS.LT.2.0                ) THEN
             XNg1(IC2)=Xt
             YNg1(IC2)=Yt
             ZNg1(IC2)=Zt
             XNg1P(IC2)=XPNp2(IN)
             YNg1P(IC2)=YPNp2(IN)
             ZNg1P(IC2)=ZPNp2(IN)
             TEST1=.TRUE.
           ENDIF
C
C XNg1,YNg1,ZNg1 coordinates of the closest Np2 atom
C
           IF(NAMATM(LLIGNpl(IN)).EQ.'  Np1'
     $        .AND.(Xt.NE.XNg2(IC2).AND.Xt.NE.XNg3(IC2))
     $        .AND.DIS.LT.2.0) THEN
             XNg2(IC2)=Xt
             YNg2(IC2)=Yt
             ZNg2(IC2)=Zt
             XNg2P1(IC2)=XP1Np1(IN)
             YNg2P1(IC2)=YP1Np1(IN)
             ZNg2P1(IC2)=ZP1Np1(IN)
             XNg2P2(IC2)=XP2Np1(IN)
             YNg2P2(IC2)=YP2Np1(IN)
             ZNg2P2(IC2)=ZP2Np1(IN)
             PK1LGC2(IC2)=PK1LGNp1(IN)
             TEST2=.TRUE.
C
C XNg2,YNg2,ZNg2 coordinates of the closest Np1 atom
C
             DO J=1,NLIGNpl
              IF(J.GT.IN) THEN
               Xt=X(LLIGNpl(J))
               Yt=Y(LLIGNpl(J))
               Zt=Z(LLIGNpl(J))              
               IF(NAMATM(LLIGNpl(J)).EQ.'  Np1'
     $            .AND.DIS.LT.2.0)  THEN
                XNg3(IC2)=Xt
                YNg3(IC2)=Yt
                ZNg3(IC2)=Zt
                XNg3P1(IC2)=XP1Np1(J)
                YNg3P1(IC2)=YP1Np1(J)
                ZNg3P1(IC2)=ZP1Np1(J)
                XNg3P2(IC2)=XP2Np1(J)
                YNg3P2(IC2)=YP2Np1(J)
                ZNg3P2(IC2)=ZP2Np1(J)
                TEST3=.TRUE.
               ENDIF
              ENDIF
             END DO
           ENDIF
         IF(TEST1.AND.TEST2.AND.TEST3) THEN
           NAMATM(LLIGC2(IC2))='  Cg '
           NAMLN(IN)='  Cg '
         ENDIF
         IF(.NOT.TEST1.AND.TEST2.AND.TEST3) THEN
           NAMATM(LLIGC2(IC2))='  CN2'
           NAMLN(IN)='  CN2'
         ENDIF
         END DO
C
C XNg3,YNg3,ZNg3 coordinates of the second closest Np1 atom
C
C        IF(TEST1.AND.TEST2.AND.TEST3) THEN
C          NAMATM(LLIGC2(IC2))='  Cg '
C        ENDIF
C        IF(.NOT.TEST1.AND.TEST2.AND.TEST3) THEN
C          NAMATM(LLIGC2(IC2))='  CN2'
C        ENDIF
C--------------------------------------------------------------------
C it doesn't work because the C is aromatic and not C2 in the HPV.pdb
C
C        IF(.NOT.TEST1.AND.(TEST2.NEQV.TEST3)) THEN
C          NAMATM(LLIGC2(IC2))='  CN1'
C        ENDIF
C--------------------------------------------------------------------
       END DO
C
C
C
C      -- DETERMINE N3 PROTON POSITIONS --
C
       NHLIGN3=0
       PI=3.14159265359D+00
       ALPHA=125.00*PI/180.00
       BETA=90.00*PI/180.00
C
       DO IN3=1,NLIGN3 
         XN=X(LLIGN3(IN3))
         YN=Y(LLIGN3(IN3))
         ZN=Z(LLIGN3(IN3))
         NB_neighbours=0
         DO I=1,NLIGAND
          IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGN3(IN3))) THEN
           Xb=X(LLIGAND(I))
           Yb=Y(LLIGAND(I))
           Zb=Z(LLIGAND(I))
           DISn=SQRT((Xb-XN)**2+(Yb-YN)**2+(Zb-ZN)**2)
           IF (DISn.LT.2) THEN
              NB_neighbours=NB_neighbours+1
           ENDIF
          ENDIF
         END DO
C
         IF (NB_neighbours.EQ.0) THEN
           NAMATM(LLIGN3(IN3))='  N30'
         ENDIF
         IF (NB_neighbours.EQ.1) THEN
           NAMATM(LLIGN3(IN3))='  N31'
         ENDIF
         IF (NB_neighbours.EQ.2) THEN
           NAMATM(LLIGN3(IN3))='  N32'
           DIST=1000.0
           DO J=1,NLIGAND
            IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGN3(IN3))) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
              DIST=MIN(DIST,DIS)
              IF(DIS.EQ.DIST) THEN
                X1=Xt
                Y1=Yt
                Z1=Zt
              ENDIF
            ENDIF
           END DO
C
C -- X1,Y1,Z1 coordinates of the closest neighbour --
C
           DIST=1000.0
           DO J=1,NLIGAND
             IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGN3(IN3))
     $          .AND.X(LLIGAND(J)).NE.X1) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
              DIST=MIN(DIST,DIS)
               IF(DIS.EQ.DIST) THEN
                X2=Xt
                Y2=Yt
                Z2=Zt
               ENDIF
             ENDIF
           END DO
C
C
C -- X2,Y2,Z2 coordinates of the 2nd closest neighbour --
C
C
C  --  normalized vectors V1 and V2  --
C
         VECNRM1=SQRT((X1-XN)**2+(Y1-YN)**2+(Z1-ZN)**2)
         XV1=(X1-XN)/VECNRM1
         YV1=(Y1-YN)/VECNRM1
         ZV1=(Z1-ZN)/VECNRM1
C
         VECNRM2=SQRT((X2-XN)**2+(Y2-YN)**2+(Z2-ZN)**2)
         XV2=(X2-XN)/VECNRM2
         YV2=(Y2-YN)/VECNRM2
         ZV2=(Z2-ZN)/VECNRM2
C
C   -- M middle of 1 and 2 (M->N direction of Z axe) --
C
         XM=((XN+XV1)+(XN+XV2))/2
         YM=((YN+YV1)+(YN+YV2))/2
         ZM=((ZN+ZV1)+(ZN+ZV2))/2
         VECNRMM=SQRT((XM-XN)**2+(YM-YN)**2+(ZM-ZN)**2)
         XVM=(XM-XN)/VECNRMM
         YVM=(YM-YN)/VECNRMM
         ZVM=(ZM-ZN)/VECNRMM
         XVZ=XN-XVM
         YVZ=YN-YVM
         ZVZ=ZN-ZVM
C
C
C   --  VP12 vector product of V1 and V2  -- 
C
         XVP12=YV1*ZV2-ZV1*YV2
         YVP12=ZV1*XV2-XV1*ZV2
         ZVP12=XV1*YV2-YV1*XV2
C        
C   -- V3 unit vector of VP12.THETA (direction of Y axe) --
C   -- THETA angle between V1 and V2 -- 
C 
         SCAL=XV1*XV2+YV1*YV2+ZV1*ZV2
         XV3=XVP12/SQRT(1-SCAL**2)
         YV3=YVP12/SQRT(1-SCAL**2)
         ZV3=ZVP12/SQRT(1-SCAL**2)
         XVY=XN+XV3
         YVY=YN+YV3
         ZVY=ZN+ZV3
C
C
C   --  V4 vector product of VM and V3  (direction of X axe) --
C
         XPM3=YVM*ZV3-ZVM*YV3
         YPM3=ZVM*XV3-XVM*ZV3
         ZPM3=XVM*YV3-YVM*XV3
         XVX=XN+XPM3
         YVX=YN+YPM3
         ZVX=ZN+ZPM3
C
C -- VX,VY,VZ are a new set of mutually orthogonal axes  --
C -- the coord of the protons are expressed in this new  --
C -- coordinates system using the spherical coordinates  --
C
         XVEC=XVM*COS(ALPHA)+XVPM3*SIN(ALPHA)*COS(BETA)
     $        +XV3*SIN(ALPHA)*SIN(BETA)
         YVEC=YVM*COS(ALPHA)+YVPM3*SIN(ALPHA)*COS(BETA)
     $        +YV3*SIN(ALPHA)*SIN(BETA)
         ZVEC=ZVM*COS(ALPHA)+ZVPM3*SIN(ALPHA)*COS(BETA)
     $        +ZV3*SIN(ALPHA)*SIN(BETA)
C
         NHLIGN3=NHLIGN3+1
         XN3P1(IN3)=XN+XVEC
         YN3P1(IN3)=YN+YVEC
         ZN3P1(IN3)=ZN+ZVEC
         XVECP=XVM*COS(ALPHA)+XVPM3*SIN(ALPHA)*COS(-BETA)
     $        +XV3*SIN(ALPHA)*SIN(-BETA)
         YVECP=YVM*COS(ALPHA)+YVPM3*SIN(ALPHA)*COS(-BETA)
     $        +YV3*SIN(ALPHA)*SIN(-BETA)
         ZVECP=ZVM*COS(ALPHA)+ZVPM3*SIN(ALPHA)*COS(-BETA)
     $        +ZV3*SIN(ALPHA)*SIN(-BETA)
         NHLIGN3=NHLIGN3+1 
         XN3P2(IN3)=XN+XVECP
         YN3P2(IN3)=YN+YVECP
         ZN3P2(IN3)=ZN+ZVECP
         ENDIF
C
         IF (NB_neighbours.EQ.3) THEN
           NAMATM(LLIGN3(IN3))='  N33'
           DIST=1000.0
           DO J=1,NLIGAND
            IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGN3(IN3))) THEN
            Xt=X(LLIGAND(J))
            Yt=Y(LLIGAND(J))
            Zt=Z(LLIGAND(J))
           DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
           DIST=MIN(DIST,DIS)
            IF(DIS.EQ.DIST) THEN
             X1=Xt
             Y1=Yt
             Z1=Zt
            ENDIF
            ENDIF
           END DO
C
C -- X1,Y1,Z1 coordinates of the closest neighbour --
C
           DIST=1000.0
           DO J=1,NLIGAND
             IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGN3(IN3))
     $          .AND.X(LLIGAND(J)).NE.X1) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
              DIST=MIN(DIST,DIS)
               IF(DIS.EQ.DIST) THEN
                X2=Xt
                Y2=Yt
                Z2=Zt
               ENDIF
             ENDIF
           END DO
C
C -- X2,Y2,Z2 coordinates of the 2nd closest neighbour --
C
           DIST=1000.0
           DO J=1,NLIGAND
             IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGN3(IN3))
     $          .AND.X(LLIGAND(J)).NE.X1
     $          .AND.X(LLIGAND(J)).NE.X2) THEN
              Xt=X(LLIGAND(J))
              Yt=Y(LLIGAND(J))
              Zt=Z(LLIGAND(J))
              DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
              DIST=MIN(DIST,DIS)
               IF(DIS.EQ.DIST) THEN
                X3=Xt
                Y3=Yt
                Z3=Zt
               ENDIF
             ENDIF
           END DO
C
C -- X3,Y3,Z3 coordinates of the 3rd closest neighbour --
C
C -- N1 vector N-nb1, N2 vector N-nb2, N3 vector N-nb3 --
C -- V1,V2,V3 are the corresponding normalized vectors --
C
           VECNRM1=SQRT((XN-X1)**2+(YN-Y1)**2+(ZN-Z1)**2)
           XV1=(XN-X1)/VECNRM1
           YV1=(YN-Y1)/VECNRM1
           ZV1=(ZN-Z1)/VECNRM1
C
           VECNRM2=SQRT((XN-X2)**2+(YN-Y2)**2+(ZN-Z2)**2)
           XV2=(XN-X2)/VECNRM2
           YV2=(YN-Y2)/VECNRM2
           ZV2=(ZN-Z2)/VECNRM2
C
           VECNRM3=SQRT((XN-X3)**2+(YN-Y3)**2+(ZN-Z3)**2)
           XV3=(XN-X3)/VECNRM3
           YV3=(YN-Y3)/VECNRM3
           ZV3=(ZN-Z3)/VECNRM3
C
C -- center of mass of V1,V2,V3 --
C
           XB=(3*XN+XV1+XV2+XV3)/3.0
           YB=(3*YN+YV1+YV2+YV3)/3.0
           ZB=(3*ZN+ZV1+ZV2+ZV3)/3.0
           VECNRM=SQRT((XN-XB)**2+(YN-YB)**2+(ZN-ZB)**2)
           XVEC=(XN-XB)/VECNRM
           YVEC=(YN-YB)/VECNRM
           ZVEC=(ZN-ZB)/VECNRM
           N3=N3+1
C
C -- assign a proton along the direction BN at 1.0 Angstrm --
C
           XN3P(IN3)=XN-XVEC
           YN3P(IN3)=YN-YVEC
           ZN3P(IN3)=ZN-ZVEC
         ENDIF
C
C
       END DO
C      WRITE(6,*)NHLIGN3,' N3 PROTONS ASSIGNED'
C
C
C
C      -- DETERMINE N1 GHOST POSITIONS --
C
       N1LIG=0
       DO IN1=1,NLIGN1
         XN=X(LLIGN1(IN1))
         YN=Y(LLIGN1(IN1))
         ZN=Z(LLIGN1(IN1))
         DIST=1000.0
         DO J=1,NLIGAND
          IF(NBATM(LLIGAND(J)).NE.NBATM(LLIGN1(IN1))) THEN
           Xt=X(LLIGAND(J))
           Yt=Y(LLIGAND(J))
           Zt=Z(LLIGAND(J))
           DIS=SQRT((Xt-XN)**2+(Yt-YN)**2+(Zt-ZN)**2)
           DIST=MIN(DIST,DIS)
            IF(DIS.EQ.DIST) THEN
             Xb=Xt
             Yb=Yt
             Zb=Zt
            ENDIF
          ENDIF
         END DO
C
C Xb,Yb,Zb coordinates of the atom linked to N1
C
         VECNRM=SQRT((XN-Xb)**2+(YN-Yb)**2+(ZN-Zb)**2)
         XVEC=(XN-Xb)/VECNRM
         YVEC=(YN-Yb)/VECNRM
         ZVEC=(ZN-Zb)/VECNRM
         N1LIG=N1LIG+1
C
C -- assign a proton along the direction BN at 1.0 Angstrm --
C
         XN1G(IN1)=XN+XVEC
         YN1G(IN1)=YN+YVEC
         ZN1G(IN1)=ZN+ZVEC
       END DO
C
       DO I=1,1000
         TYPCAR(I)='SUFACE'
         TYPHIS(I)='SUFACE'
         TYPCYS(I)='SUFACE'
         TYPTYR(I)='SUFACE'
         TYPLYS(I)='SUFACE'
         TYPARG(I)='SUFACE'
         TYPLGN3(I)='SUFACE'
         TYPLGNar(I)='SUFACE'
         TYPLGCAR(I)='SUFACE'
         TYPLGCg(I)='SUFACE'
         TYPLGCN2(I)='SUFACE'
       END DO
C
C
C        -- FIND DISULFIDE BONDS --
C
       DO ICYS=1,NCYS
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
         DO JCYS=1,NCYS
           IF(JCYS.NE.ICYS.AND.
     $        LCYSRS(ICYS).NE.LCYSRS(JCYS)-3 .AND.
     $        LCYSRS(ICYS).NE.LCYSRS(JCYS)+3   )THEN
             XS=X(LCYSSG(JCYS))
             YS=Y(LCYSSG(JCYS))
             ZS=Z(LCYSSG(JCYS))
             DIS=SQRT((XSG-XS)**2+(YSG-YS)**2+(ZSG-ZS)**2)
             IF(DIS.LT.2.50)THEN
               TYPCYS(ICYS)='BONDED'
             END IF
           END IF
         END DO
       END DO
C
C
C
C      INTRINSIC pKa VALUES WHICH ARE pH INDEPENDENT
C
       DO I=1,1000
         NHBCAR(I)=0
         NHBHIS(I)=0
         NSDCAR(I)=0
         NBKCAR(I)=0
         NCLCAR(I)=0
         NSDHIS(I)=0
         NBKHIS(I)=0
         NCLHIS(I)=0
         NSDCYS(I)=0
         NBKCYS(I)=0
         NCLCYS(I)=0
         NSDTYR(I)=0
         NBKTYR(I)=0
         NCLTYR(I)=0
         NSDLYS(I)=0
         NBKLYS(I)=0
         NCLLYS(I)=0
         NSDARG(I)=0
         NBKARG(I)=0
         NCLARG(I)=0
         NLGCAR(I)=0
         NCLLCAR(I)=0
         NLGHIS(I)=0
         NCLLHIS(I)=0
         NLGCYS(I)=0
         NCLLCYS(I)=0
         NLGTYR(I)=0
         NCLLTYR(I)=0
c dmr
         nlglys(i)=0
         nclllys(i)=0
         nlgarg(i)=0
         ncllarg(i)=0
c dmr
         NCLN3(I)=0
         NCLCN2(I)=0
         NCLCg(I)=0
         NCLLGCAR(I)=0 
         NCLNar(I)=0
         DO J=1,20
           NMASS(J,I)=0
           NLOCAL(J,I)=0
           TOLBKB(J,I)=0.0
           TOLSDC(J,I)=0.0
           TOLCOL(J,I)=0.0
           TOLMAS(J,I)=0.0
           TOLLOC(J,I)=0.0
           DO K=1,30
             NAMBKB(J,I,K)='000'
             NUMBKB(J,I,K)=0
             VALBKB(J,I,K)=0.0
             NAMSDC(J,I,K)='000'
             NUMSDC(J,I,K)=0
             VALSDC(J,I,K)=0.0
             NAMCOL(J,I,K)='000'
             NUMCOL(J,I,K)=0
             VALCOL(J,I,K)=0.0
             VALLCOL(J,I,K)=0.0
           END DO
         END DO
       END DO
       DO I=0,1000
         NLGCAR(I)=0
         DO J=0,10
          DO K=0,30
          VALLIG(J,I,K)=0.0
          END DO
         END DO
       END DO
C
       DO ICAR=1, NCAR
         IF(NAMCAR(ICAR).EQ.'C- ')PK1CAR(ICAR)=3.20
         IF(NAMCAR(ICAR).EQ.'ASP')PK1CAR(ICAR)=3.80
         IF(NAMCAR(ICAR).EQ.'GLU')PK1CAR(ICAR)=4.50
       END DO
       DO IHIS=1, NHIS
         PK1HIS(IHIS)=6.50
       END DO
       DO ICYS=1, NCYS
         PK1CYS(ICYS)=99.99
         IF(TYPCYS(ICYS).NE.'BONDED')PK1CYS(ICYS)=9.00
       END DO
       DO ITYR=1, NTYR
         PK1TYR(ITYR)=10.00
       END DO
       PK1LYS(1)=8.0
       DO ILYS=2, NLYS
         PK1LYS(ILYS)=10.50
       END DO
       DO IARG=1, NARG
         PK1ARG(IARG)=12.50
       END DO
       DO IN3=1, NLIGN3
         PK1LGN3(IN3)=PKALIG(LLIGN3(IN3))
         pkamod(7,IN3)=PKALIG(LLIGN3(IN3))
       END DO
       DO INar=1, NLIGNar
         PK1LGNar(INar)=PKALIG(LLIGNar(INar))
         pkamod(8,INar)=PKALIG(LLIGNar(INar))
       END DO
       DO IC2=1, NLIGC2
        IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           PK1LGCg(IC2)=PK1LGC2(IC2)
           pkamod(9,IC2)=PK1LGC2(IC2)
        ENDIF
        IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           PK1LGCN2(IC2)=PK1LGC2(IC2)
           pkamod(9,IC2)=PK1LGC2(IC2)
        ENDIF
       END DO
C
C      -- IDENTIFICATION OF THE CAR GROUPS IN THE LIGAND --
C      -- Oco identified by pair using the C2 to which they are linked--
C
       DO IC2=1,NLIGC2
         XLO1(IC2)=0
         YLO1(IC2)=0
         ZLO1(IC2)=0
         XLO2(IC2)=0
         YLO2(IC2)=0
         ZLO2(IC2)=0
         XCac(IC2)=0
         YCac(IC2)=0
         ZCac(IC2)=0
       END DO
       DO IC2=1,NLIGC2
         TEST=.FALSE.
         XC=X(LLIGC2(IC2))
         YC=Y(LLIGC2(IC2))
         ZC=Z(LLIGC2(IC2))
         DIST=1000.0
         DO ICAR=1,NLIGCAR
           Xt=X(LLIGCAR(ICAR))
           Yt=Y(LLIGCAR(ICAR))
           Zt=Z(LLIGCAR(ICAR))
           DIS=SQRT((Xt-XC)**2+(Yt-YC)**2+(Zt-ZC)**2)
           DIST=MIN(DIS,DIST)
           IF (DIS.EQ.DIST.AND.DIST.LT.2.0) THEN
             XLO1(IC2)=Xt
             YLO1(IC2)=Yt
             ZLO1(IC2)=Zt
             XCac(IC2)=XC
             YCac(IC2)=YC
             ZCac(IC2)=ZC
             PK1LGCAR(IC2)=PKALIG(LLIGCAR(ICAR))
             pkamod(11,IC2)=PKALIG(LLIGCAR(ICAR))
             TEST=.TRUE.  
           ENDIF
         END DO
C
C -- XLO1,YLO1,ZLO1 coordinates of the closest neighbour --
C
C
         IF(TEST)THEN
         DIST=1000.0
         DO JCAR=1,NLIGCAR
            IF(X(LLIGCAR(JCAR)).NE.XLO1(IC2)) THEN
              Xt=X(LLIGCAR(JCAR))
              Yt=Y(LLIGCAR(JCAR))
              Zt=Z(LLIGCAR(JCAR))
              DIS=SQRT((Xt-XC)**2+(Yt-YC)**2+(Zt-ZC)**2)
              DIST=MIN(DIST,DIS)
               IF(DIS.EQ.DIST) THEN
                XLO2(IC2)=Xt
                YLO2(IC2)=Yt
                ZLO2(IC2)=Zt
               ENDIF
            ENDIF
         END DO
         ENDIF
C
C
C -- XLO2,YLO2,ZLO2 coordinates of the 2nd closest neighbour --
C
       END DO
C
C
C      *******************
C      STEP 3. DESOLVATION
C      *******************
C
C      -- ASP/GLU --
C
       DO ICAR=1, NCAR
         XO1=X(LCARO1(ICAR))
         YO1=Y(LCARO1(ICAR))
         ZO1=Z(LCARO1(ICAR))
         XO2=X(LCARO2(ICAR))
         YO2=Y(LCARO2(ICAR))
         ZO2=Z(LCARO2(ICAR))
         XO=(XO1+XO2)/2.0
         YO=(YO1+YO2)/2.0
         ZO=(ZO1+ZO2)/2.0
C
         FMASS=0.010
         FLOCAL=0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(1,ICAR)=0
         DO IATOM=1,NATOM
           IF((HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).NE.LCARRS(ICAR)) .OR.
     $        (HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).EQ.LCARRS(ICAR) .AND.
     $        TYPCH(IATOM).NE.TYPCH(LCARO1(ICAR)))
     $        .OR. HEAD(IATOM).EQ.'LG') THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF((ABS(XO-XDS).LT.DMASS .AND.
     $         ABS(YO-YDS).LT.DMASS .AND.
     $         ABS(ZO-ZDS).LT.DMASS     )     ) THEN
             DIS=SQRT((XO-XDS)**2+(YO-YDS)**2+(ZO-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(1,ICAR)=NMASS(1,ICAR)+1 
             IF(DIS.LT.DLOCL)NLOCAL(1,ICAR)=NLOCAL(1,ICAR)+1
           END IF
           END IF
         END DO
         IF(NMASS(1,ICAR).GT.400)TYPCAR(ICAR)='BURIED'
         TOLMAS(1,ICAR)=MAX(0.0,FMASS*(NMASS(1,ICAR)-400))
         TOLLOC(1,ICAR)=FLOCAL*NLOCAL(1,ICAR)
         PK1CAR(ICAR)=PK1CAR(ICAR)+TOLMAS(1,ICAR)+TOLLOC(1,ICAR)
       END DO
C
C      -- HIS --
C
       DO IHIS=1, NHIS
         XH1=XHISP1(IHIS)
         YH1=YHISP1(IHIS)
         ZH1=ZHISP1(IHIS)
         XH2=XHISP2(IHIS)
         YH2=YHISP2(IHIS)
         ZH2=ZHISP2(IHIS)
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XCT=(XCG+XND+XCE+XNE+XCD)/5.0
         YCT=(YCG+YND+YCE+YNE+YCD)/5.0
         ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL1=4.00
         DLOCL2=6.00
         DMASS=15.50
         NMASS(2,IHIS)=0
         DO IATOM=1,NATOM
           IF((HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).NE.LHISRS(IHIS)) .OR.
     $        (HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).EQ.LHISRS(IHIS) .AND.
     $        TYPCH(IATOM).NE.TYPCH(LHISCG(IHIS)))
     $        .OR. HEAD(IATOM).EQ.'LG') THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF((ABS(XCT-XDS).LT.DMASS .AND.
     $         ABS(YCT-YDS).LT.DMASS .AND.
     $         ABS(ZCT-ZDS).LT.DMASS     )     ) THEN
             DIS=SQRT((XCT-XDS)**2+(YCT-YDS)**2+(ZCT-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(2,IHIS)=NMASS(2,IHIS)+1
             IF(DIS.LT.DLOCL1)NLOCAL(2,IHIS)=NLOCAL(2,IHIS)+1
             IF(DIS.LT.DLOCL2)NLOCAL(10,IHIS)=NLOCAL(10,IHIS)+1
           END IF
           END IF
         END DO
         IF(NMASS(2,IHIS).GT.400)TYPHIS(IHIS)='BURIED'
         TOLMAS(2,IHIS)=MIN(0.0,FMASS*(NMASS(2,IHIS)-400))
         IF(NMASS(2,IHIS).GT.400)NLOCAL(2,IHIS)=NLOCAL(10,IHIS)
         TOLLOC(2,IHIS)=FLOCAL*NLOCAL(2,IHIS)
         PK1HIS(IHIS)=PK1HIS(IHIS)+TOLMAS(2,IHIS)+TOLLOC(2,IHIS)
       END DO
C
C      -- CYS --
C
       DO ICYS=1, NCYS
       IF(TYPCYS(ICYS).NE.'BONDED')THEN
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
C
         FMASS=0.010
         FLOCAL=0.070
         DLOCL=3.50
         DMASS=15.50
         NMASS(3,ICYS)=0
         DO IATOM=1,NATOM
           IF((HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).NE.LCYSRS(ICYS)) .OR.
     $        (HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).EQ.LCYSRS(ICYS) .AND.
     $        TYPCH(IATOM).NE.TYPCH(LCYSSG(ICYS)))
     $        .OR. HEAD(IATOM).EQ.'LG') THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XSG-XDS).LT.DMASS .AND.
     $        ABS(YSG-YDS).LT.DMASS .AND.
     $        ABS(ZSG-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XSG-XDS)**2+(YSG-YDS)**2+(ZSG-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(3,ICYS)=NMASS(3,ICYS)+1
             IF(DIS.LT.DLOCL)NLOCAL(3,ICYS)=NLOCAL(3,ICYS)+1
           END IF
           END IF
         END DO
         IF(NMASS(3,ICYS).GT.400)TYPCYS(ICYS)='BURIED'
         TOLMAS(3,ICYS)=MAX(0.0,FMASS*(NMASS(3,ICYS)-400))
         TOLLOC(3,ICYS)=FLOCAL*NLOCAL(3,ICYS)
         PK1CYS(ICYS)=PK1CYS(ICYS)+TOLMAS(3,ICYS)+TOLLOC(3,ICYS)
       END IF
       END DO
C
C      -- TYR --
C
       DO ITYR=1, NTYR
         XOH=X(LTYROH(ITYR))
         YOH=Y(LTYROH(ITYR))
         ZOH=Z(LTYROH(ITYR))
C
         FMASS=0.010
         FLOCAL=0.070
         DLOCL=3.50
         DMASS=15.50
         NMASS(4,ITYR)=0
         DO IATOM=1,NATOM
           IF((HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).NE.LTYRRS(ITYR)) .OR.
     $        (HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).EQ.LTYRRS(ITYR) .AND.
     $        TYPCH(IATOM).NE.TYPCH(LTYROH(ITYR)))
     $        .OR. HEAD(IATOM).EQ.'LG') THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XOH-XDS).LT.DMASS .AND.
     $        ABS(YOH-YDS).LT.DMASS .AND.
     $        ABS(ZOH-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XOH-XDS)**2+(YOH-YDS)**2+(ZOH-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(4,ITYR)=NMASS(4,ITYR)+1
             IF(DIS.LT.DLOCL)NLOCAL(4,ITYR)=NLOCAL(4,ITYR)+1
           END IF
           END IF
         END DO
         IF(NMASS(4,ITYR).GT.400)TYPTYR(ITYR)='BURIED'
         TOLMAS(4,ITYR)=MAX(0.0,FMASS*(NMASS(4,ITYR)-400))
         TOLLOC(4,ITYR)=FLOCAL*NLOCAL(4,ITYR)
         PK1TYR(ITYR)=PK1TYR(ITYR)+TOLMAS(4,ITYR)+TOLLOC(4,ITYR)
       END DO
C
C      -- LYS --
C
       DO ILYS=1, NLYS
         XNZ=X(LLYSNZ(ILYS))
         YNZ=Y(LLYSNZ(ILYS))
         ZNZ=Z(LLYSNZ(ILYS))
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(5,ILYS)=0
         DO IATOM=1,NATOM
           IF((HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).NE.LLYSRS(ILYS)) .OR.
     $        (HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).EQ.LLYSRS(ILYS) .AND.
     $        TYPCH(IATOM).NE.TYPCH(LLYSNZ(ILYS)))
     $        .OR. HEAD(IATOM).EQ.'LG') THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XNZ-XDS).LT.DMASS .AND.
     $        ABS(YNZ-YDS).LT.DMASS .AND.
     $        ABS(ZNZ-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XNZ-XDS)**2+(YNZ-YDS)**2+(ZNZ-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(5,ILYS)=NMASS(5,ILYS)+1
             IF(DIS.LT.DLOCL)NLOCAL(5,ILYS)=NLOCAL(5,ILYS)+1
           END IF
           END IF
         END DO
         IF(NMASS(5,ILYS).GT.400)TYPLYS(ILYS)='BURIED'
         TOLMAS(5,ILYS)=MIN(0.0,FMASS*(NMASS(5,ILYS)-400))
         TOLLOC(5,ILYS)=FLOCAL*NLOCAL(5,ILYS)
         PK1LYS(ILYS)=PK1LYS(ILYS)+TOLMAS(5,ILYS)+TOLLOC(5,ILYS)
       END DO
C
C      -- ARG --
C
       DO IARG=1, NARG
         X1=X(LARGN1(IARG))
         Y1=Y(LARGN1(IARG))
         Z1=Z(LARGN1(IARG))
         X2=X(LARGN2(IARG))
         Y2=Y(LARGN2(IARG))
         Z2=Z(LARGN2(IARG))
         X3=X(LARGN3(IARG))
         Y3=Y(LARGN3(IARG))
         Z3=Z(LARGN3(IARG))
         XCZ=X(LARGCZ(IARG))
         YCZ=Y(LARGCZ(IARG))
         ZCZ=Z(LARGCZ(IARG))
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL=5.00
         DMASS=15.50
         NMASS(6,IARG)=0
         DO IATOM=1,NATOM
           IF((HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).NE.LARGRS(IARG)) .OR.
     $        (HEAD(IATOM).EQ.'AT' .AND.
     $        NUMRES(IATOM).EQ.LARGRS(IARG) .AND.
     $        TYPCH(IATOM).NE.TYPCH(LARGN1(IARG)))
     $        .OR. HEAD(IATOM).EQ.'LG') THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XCZ-XDS).LT.DMASS .AND.
     $        ABS(YCZ-YDS).LT.DMASS .AND.
     $        ABS(ZCZ-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XCZ-XDS)**2+(YCZ-YDS)**2+(ZCZ-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(6,IARG)=NMASS(6,IARG)+1
             IF(DIS.LT.DLOCL)NLOCAL(6,IARG)=NLOCAL(6,IARG)+1
           END IF
           END IF
         END DO
         IF(NMASS(6,IARG).GT.400)TYPARG(IARG)='BURIED'
         TOLMAS(6,IARG)=MIN(0.0,FMASS*(NMASS(6,IARG)-400))
         TOLLOC(6,IARG)=FLOCAL*NLOCAL(6,IARG)
         PK1ARG(IARG)=PK1ARG(IARG)+TOLMAS(6,IARG)+TOLLOC(6,IARG)
       END DO
C
C
C      -- LIGAND --
C
C       - N3 atoms -
C
       DO IN3=1, NLIGN3
C
         XN3=X(LLIGN3(IN3))
         YN3=Y(LLIGN3(IN3))
         ZN3=Z(LLIGN3(IN3))
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(7,IN3)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'AT') THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XN3-XDS).LT.DMASS .AND.
     $        ABS(YN3-YDS).LT.DMASS .AND.
     $        ABS(ZN3-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XN3-XDS)**2+(YN3-YDS)**2+(ZN3-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(7,IN3)=NMASS(7,IN3)+1
             IF(DIS.LT.DLOCL)NLOCAL(7,IN3)=NLOCAL(7,IN3)+1
           END IF
           END IF
         END DO
         IF(NMASS(7,IN3).GT.400)TYPLGN3(IN3)='BURIED'
         TOLMAS(7,IN3)=MIN(0.0,FMASS*(NMASS(7,IN3)-400))
         TOLLOC(7,IN3)=FLOCAL*NLOCAL(7,IN3)
         PK1LGN3(IN3)=PK1LGN3(IN3)+TOLMAS(7,IN3)+TOLLOC(7,IN3)
       END DO
C
C
C       - Nar atoms -
C
       DO INar=1, NLIGNar
C
         XNar=X(LLIGNar(INar))
         YNar=Y(LLIGNar(INar))
         ZNar=Z(LLIGNar(INar))
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(8,INar)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'AT' ) THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(XNar-XDS).LT.DMASS .AND.
     $        ABS(YNar-YDS).LT.DMASS .AND.
     $        ABS(ZNar-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XNar-XDS)**2+(YNar-YDS)**2+(ZNar-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(8,INar)=NMASS(8,INar)+1
             IF(DIS.LT.DLOCL)NLOCAL(8,INar)=NLOCAL(8,INar)+1
           END IF
           END IF
         END DO
         IF(NMASS(8,INar).GT.400)TYPLGNar(INar)='BURIED'
         TOLMAS(8,INar)=MIN(0.0,FMASS*(NMASS(8,INar)-400))
         TOLLOC(8,INar)=FLOCAL*NLOCAL(8,INar)
         PK1LGNar(INar)=PK1LGNar(INar)+TOLMAS(8,INar)+TOLLOC(8,INar)
       END DO
C
C
C       - Npl atoms -
C
       DO IC2=1,NLIGC2
C
        FMASS=-0.010
        FLOCAL=-0.070
        DLOCL=4.50
        DMASS=15.50
        NMASS(9,IC2)=0
C
C  -- GUADININIUM TYPE GROUP (eg. 1K1O.pdb) --
C
C
          XC=X(LLIGC2(IC2))
          YC=Y(LLIGC2(IC2))
          ZC=Z(LLIGC2(IC2))
        IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
          DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'AT' ) THEN
              XDS=X(IATOM)
              YDS=Y(IATOM)
              ZDS=Z(IATOM)
            IF(ABS(XC-XDS).LT.DMASS .AND.
     $         ABS(YC-YDS).LT.DMASS .AND.
     $         ABS(ZC-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XC-XDS)**2+(YC-YDS)**2+(ZC-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(9,IC2)=NMASS(9,IC2)+1
             IF(DIS.LT.DLOCL)NLOCAL(9,IC2)=NLOCAL(9,IC2)+1
            END IF
           END IF
          END DO
         IF(NMASS(9,IC2).GT.400)TYPLGCg(IC2)='BURIED'
         TOLMAS(9,IC2)=MIN(0.0,FMASS*(NMASS(9,IC2)-400))
         TOLLOC(9,IC2)=FLOCAL*NLOCAL(9,IC2)
         PK1LGCg(IC2)=PK1LGCg(IC2)+TOLMAS(9,IC2)+TOLLOC(9,IC2)
        ENDIF
C
C  -- Np1 pair GROUP (eg. 1K1I.pdb) --
C
        IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
          DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'AT' ) THEN
              XDS=X(IATOM)
              YDS=Y(IATOM)
              ZDS=Z(IATOM)
            IF(ABS(XC-XDS).LT.DMASS .AND.
     $         ABS(YC-YDS).LT.DMASS .AND.
     $         ABS(ZC-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((XC-XDS)**2+(YC-YDS)**2+(ZC-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(9,IC2)=NMASS(9,IC2)+1
             IF(DIS.LT.DLOCL)NLOCAL(9,IC2)=NLOCAL(9,IC2)+1
            END IF
           END IF
          END DO
         IF(NMASS(9,IC2).GT.400)TYPLGCN2(IC2)='BURIED'
         TOLMAS(9,IC2)=MIN(0.0,FMASS*(NMASS(9,IC2)-400))
         TOLLOC(9,IC2)=FLOCAL*NLOCAL(9,IC2)
         PK1LGCN2(IC2)=PK1LGCN2(IC2)+TOLMAS(9,IC2)+TOLLOC(9,IC2)
        ENDIF
       END DO
C
C
C
C
C       - CAR GROUP OF THE LIGAND (Oco atoms) -
C
C
       DO IC2=1,NLIGC2
        IF(XCac(IC2).NE.0.AND.YCac(IC2).NE.0.AND.ZCac(IC2).NE.0)THEN
         XLGO1=XLO1(IC2)
         YLGO1=YLO1(IC2)
         ZLGO1=ZLO1(IC2)
         XLGO2=XLO2(IC2)
         YLGO2=YLO2(IC2)
         ZLGO2=ZLO2(IC2)
         XO=(XLGO1+XLGO2)/2.0
         YO=(YLGO1+YLGO2)/2.0
         ZO=(ZLGO1+ZLGO2)/2.0
C
         FMASS=0.010
         FLOCAL=0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(11,IC2)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'AT' ) THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF((ABS(XO-XDS).LT.DMASS .AND.
     $         ABS(YO-YDS).LT.DMASS .AND.
     $         ABS(ZO-ZDS).LT.DMASS     )     ) THEN
             DIS=SQRT((XO-XDS)**2+(YO-YDS)**2+(ZO-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(11,IC2)=NMASS(11,IC2)+1
             IF(DIS.LT.DLOCL)NLOCAL(11,IC2)=NLOCAL(11,IC2)+1
           END IF
           END IF
         END DO
         IF(NMASS(11,IC2).GT.400)TYPLGCAR(IC2)='BURIED'
         TOLMAS(11,IC2)=MAX(0.0,FMASS*(NMASS(11,IC2)-400))
         TOLLOC(11,IC2)=FLOCAL*NLOCAL(11,IC2)
         PK1LGCAR(IC2)=PK1LGCAR(IC2)+TOLMAS(11,IC2)
     $                    +TOLLOC(11,IC2)
        ENDIF
       END DO
C
C
c dmr st. below is cp of Nar atoms
C
C       - Charged atoms -
C
       DO Ichr=1, NLIGchr
C
         Xchr=X(LLIGchr(Ichr))
         Ychr=Y(LLIGchr(Ichr))
         Zchr=Z(LLIGchr(Ichr))
C
         FMASS=-0.010
         FLOCAL=-0.070
         DLOCL=4.50
         DMASS=15.50
         NMASS(12,Ichr)=0
         DO IATOM=1,NATOM
           IF(HEAD(IATOM).EQ.'AT' ) THEN
           XDS=X(IATOM)
           YDS=Y(IATOM)
           ZDS=Z(IATOM)
           IF(ABS(Xchr-XDS).LT.DMASS .AND.
     $        ABS(Ychr-YDS).LT.DMASS .AND.
     $        ABS(Zchr-ZDS).LT.DMASS     ) THEN
             DIS=SQRT((Xchr-XDS)**2+(Ychr-YDS)**2+(Zchr-ZDS)**2)
             IF(DIS.LT.DMASS)NMASS(12,Ichr)=NMASS(12,Ichr)+1
             IF(DIS.LT.DLOCL)NLOCAL(12,Ichr)=NLOCAL(12,Ichr)+1
           END IF
           END IF
         END DO
       END DO
c dmr fin
C
C
C      **********************
C      STEP 4. ASP/GLU 
C      **********************
C
       DO ICAR=1,NCAR
C
         XO1=X(LCARO1(ICAR))
         YO1=Y(LCARO1(ICAR))
         ZO1=Z(LCARO1(ICAR))
         XO2=X(LCARO2(ICAR))
         YO2=Y(LCARO2(ICAR))
         ZO2=Z(LCARO2(ICAR))
         XO=(XO1+XO2)/2.0
         YO=(YO1+YO2)/2.0
         ZO=(ZO1+ZO2)/2.0
C
C
C        -- 2. SIDECHAIN INTERACTION --
C
C           - FIND SER -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XO1-XSER).LT.DIS2 .AND.
     $         ABS(YO1-YSER).LT.DIS2 .AND.
     $         ABS(ZO1-ZSER).LT.DIS2     ) .OR.
     $        (ABS(XO2-XSER).LT.DIS2 .AND.
     $         ABS(YO2-YSER).LT.DIS2 .AND.
     $         ABS(ZO2-ZSER).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XSER)**2+(YO1-YSER)**2+(ZO1-ZSER)**2)
             DISO2O=SQRT((XO2-XSER)**2+(YO2-YSER)**2+(ZO2-ZSER)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='SER'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FOH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XO1-XTHR).LT.DIS2 .AND.
     $         ABS(YO1-YTHR).LT.DIS2 .AND.
     $         ABS(ZO1-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XO2-XTHR).LT.DIS2 .AND.
     $         ABS(YO2-YTHR).LT.DIS2 .AND.
     $         ABS(ZO2-ZTHR).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XTHR)**2+(YO1-YTHR)**2+(ZO1-ZTHR)**2)
             DISO2O=SQRT((XO2-XTHR)**2+(YO2-YTHR)**2+(ZO2-ZTHR)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='THR'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FOH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND ASN -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN1=XASNP1(IASN)
           YASN1=YASNP1(IASN)
           ZASN1=ZASNP1(IASN)
           XASN2=XASNP2(IASN)
           YASN2=YASNP2(IASN)
           ZASN2=ZASNP2(IASN)
           IF((ABS(XO1-XASN1).LT.DIS2.AND.
     $         ABS(YO1-YASN1).LT.DIS2.AND.
     $         ABS(ZO1-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XASN1).LT.DIS2.AND.
     $         ABS(YO2-YASN1).LT.DIS2.AND.
     $         ABS(ZO2-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XASN2).LT.DIS2.AND.
     $         ABS(YO1-YASN2).LT.DIS2.AND.
     $         ABS(ZO1-ZASN2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XASN2).LT.DIS2.AND.
     $         ABS(YO2-YASN2).LT.DIS2.AND.
     $         ABS(ZO2-ZASN2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XASN1)**2+(YO1-YASN1)**2+(ZO1-ZASN1)**2)
             DISO21=SQRT((XO2-XASN1)**2+(YO2-YASN1)**2+(ZO2-ZASN1)**2)
             DISO12=SQRT((XO1-XASN2)**2+(YO1-YASN2)**2+(ZO1-ZASN2)**2)
             DISO22=SQRT((XO2-XASN2)**2+(YO2-YASN2)**2+(ZO2-ZASN2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='ASN'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN1=XGLNP1(IGLN)
           YGLN1=YGLNP1(IGLN)
           ZGLN1=ZGLNP1(IGLN)
           XGLN2=XGLNP2(IGLN)
           YGLN2=YGLNP2(IGLN)
           ZGLN2=ZGLNP2(IGLN)
           IF((ABS(XO1-XGLN1).LT.DIS2.AND.
     $         ABS(YO1-YGLN1).LT.DIS2.AND.
     $         ABS(ZO1-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XGLN1).LT.DIS2.AND.
     $         ABS(YO2-YGLN1).LT.DIS2.AND.
     $         ABS(ZO2-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XGLN2).LT.DIS2.AND.
     $         ABS(YO1-YGLN2).LT.DIS2.AND.
     $         ABS(ZO1-ZGLN2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XGLN2).LT.DIS2.AND.
     $         ABS(YO2-YGLN2).LT.DIS2.AND.
     $         ABS(ZO2-ZGLN2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XGLN1)**2+(YO1-YGLN1)**2+(ZO1-ZGLN1)**2)
             DISO21=SQRT((XO2-XGLN1)**2+(YO2-YGLN1)**2+(ZO2-ZGLN1)**2)
             DISO12=SQRT((XO1-XGLN2)**2+(YO1-YGLN2)**2+(ZO1-ZGLN2)**2)
             DISO22=SQRT((XO2-XGLN2)**2+(YO2-YGLN2)**2+(ZO2-ZGLN2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='GLN'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C           - FIND TRP -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ITRP=1,NTRP
           XTRP=XTRPP1(ITRP)
           YTRP=YTRPP1(ITRP)
           ZTRP=ZTRPP1(ITRP)
           IF((ABS(XO1-XTRP).LT.DIS2 .AND.
     $         ABS(YO1-YTRP).LT.DIS2 .AND.
     $         ABS(ZO1-ZTRP).LT.DIS2     ) .OR.
     $        (ABS(XO2-XTRP).LT.DIS2 .AND.
     $         ABS(YO2-YTRP).LT.DIS2 .AND.
     $         ABS(ZO2-ZTRP).LT.DIS2     )     ) THEN
             DISO1P=SQRT((XO1-XTRP)**2+(YO1-YTRP)**2+(ZO1-ZTRP)**2)
             DISO2P=SQRT((XO2-XTRP)**2+(YO2-YTRP)**2+(ZO2-ZTRP)**2)
             DIS=MIN(DISO1P,DISO2P)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='TRP'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LTRPRS(ITRP)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND HIS H-BONDING
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IHIS=1,NHIS
           XHIS1=XHISP1(IHIS)
           YHIS1=YHISP1(IHIS)
           ZHIS1=ZHISP1(IHIS)
           XHIS2=XHISP2(IHIS)
           YHIS2=YHISP2(IHIS)
           ZHIS2=ZHISP2(IHIS)
           IF((ABS(XO1-XHIS1).LT.DIS2.AND.
     $         ABS(YO1-YHIS1).LT.DIS2.AND.
     $         ABS(ZO1-ZHIS1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XHIS1).LT.DIS2.AND.
     $         ABS(YO2-YHIS1).LT.DIS2.AND.
     $         ABS(ZO2-ZHIS1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XHIS2).LT.DIS2.AND.
     $         ABS(YO1-YHIS2).LT.DIS2.AND.
     $         ABS(ZO1-ZHIS2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XHIS2).LT.DIS2.AND.
     $         ABS(YO2-YHIS2).LT.DIS2.AND.
     $         ABS(ZO2-ZHIS2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XHIS1)**2+(YO1-YHIS1)**2+(ZO1-ZHIS1)**2)
             DISO21=SQRT((XO2-XHIS1)**2+(YO2-YHIS1)**2+(ZO2-ZHIS1)**2)
             DISO12=SQRT((XO1-XHIS2)**2+(YO1-YHIS2)**2+(ZO1-ZHIS2)**2)
             DISO22=SQRT((XO2-XHIS2)**2+(YO2-YHIS2)**2+(ZO2-ZHIS2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='HIS'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
C
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))=NAMCAR(ICAR)
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=-FNH*MIN(1.0,VALUE)
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
               IF(NMASS(1,ICAR)+NMASS(2,IHIS).GT.900 .OR.
     $     (NMASS(1,ICAR).GT.400.AND.NMASS(2,IHIS).GT.400))THEN
                 VALSDC(1,ICAR,NSDCAR(ICAR))=-1.60
                 VALSDC(2,IHIS,NSDHIS(IHIS))=+1.60
                 NHBCAR(ICAR)=NHBCAR(ICAR) + 1
                 NHBHIS(IHIS)=NHBHIS(IHIS) + 1
                 PK1CAR(ICAR)=PK1CAR(ICAR)-6.00
                 PK1HIS(IHIS)=PK1HIS(IHIS)+6.00
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C           - FIND CYS-SH -
C
         FSH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED')THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XO1-XS).LT.DIS2.AND.
     $         ABS(YO1-YS).LT.DIS2.AND.
     $         ABS(ZO1-ZS).LT.DIS2    ) .OR.
     $        (ABS(XO2-XS).LT.DIS2.AND.
     $         ABS(YO2-YS).LT.DIS2.AND.
     $         ABS(ZO2-ZS).LT.DIS2    )     ) THEN
             DISO1S=SQRT((XO1-XS)**2+(YO1-YS)**2+(ZO1-ZS)**2)
             DISO2S=SQRT((XO2-XS)**2+(YO2-YS)**2+(ZO2-ZS)**2)
             DIS=MIN(DISO1S,DISO2S)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='CYS'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FSH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
C
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))=NAMCAR(ICAR)
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=-FSH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END IF
         END DO
C
C
C           - FIND TYR-OH -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ITYR=1,NTYR
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF((ABS(XO1-XOH).LT.DIS2.AND.
     $         ABS(YO1-YOH).LT.DIS2.AND.
     $         ABS(ZO1-ZOH).LT.DIS2    ) .OR.
     $        (ABS(XO2-XOH).LT.DIS2.AND.
     $         ABS(YO2-YOH).LT.DIS2.AND.
     $         ABS(ZO2-ZOH).LT.DIS2    )     ) THEN
             DISO1O=SQRT((XO1-XOH)**2+(YO1-YOH)**2+(ZO1-ZOH)**2)
             DISO2O=SQRT((XO2-XOH)**2+(YO2-YOH)**2+(ZO2-ZOH)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='TYR'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FOH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
C
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))=NAMCAR(ICAR)
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=-FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND LYS H-BONDING
C
         FNH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ILYS=1,NLYS
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           FNH=-0.80
           DIS2=4.00
           IF(ILYS.EQ.1)FNH=-1.20
           IF(ILYS.EQ.1)DIS2=4.50
           IF((ABS(XO1-XN).LT.DIS2.AND.
     $         ABS(YO1-YN).LT.DIS2.AND.
     $         ABS(ZO1-ZN).LT.DIS2    ) .OR.
     $        (ABS(XO2-XN).LT.DIS2.AND.
     $         ABS(YO2-YN).LT.DIS2.AND.
     $         ABS(ZO2-ZN).LT.DIS2    )     ) THEN
             DISO1N=SQRT((XO1-XN)**2+(YO1-YN)**2+(ZO1-ZN)**2)
             DISO2N=SQRT((XO2-XN)**2+(YO2-YN)**2+(ZO2-ZN)**2)
             DIS=MIN(DISO1N,DISO2N)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='LYS'
               IF(ILYS.EQ.1)NAMSDC(1,ICAR,NSDCAR(ICAR))='N+ '
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C           - FIND ARG H-BONDING
C
         FNH=-0.80
         DIS1=2.00
         DIS2=4.00
         DO IARG=1,NARG
           XH1=XARGP1(IARG)
           YH1=YARGP1(IARG)
           ZH1=ZARGP1(IARG)
           XH2=XARGP2(IARG)
           YH2=YARGP2(IARG)
           ZH2=ZARGP2(IARG)
           XH3=XARGP3(IARG)
           YH3=YARGP3(IARG)
           ZH3=ZARGP3(IARG)
           XH4=XARGP4(IARG)
           YH4=YARGP4(IARG)
           ZH4=ZARGP4(IARG)
           XH5=XARGP5(IARG)
           YH5=YARGP5(IARG)
           ZH5=ZARGP5(IARG)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH5).LT.DIS2.AND.
     $         ABS(YO1-YH5).LT.DIS2.AND.
     $         ABS(ZO1-ZH5).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH5).LT.DIS2.AND.
     $         ABS(YO2-YH5).LT.DIS2.AND.
     $         ABS(ZO2-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS15=SQRT((XO1-XH5)**2+(YO1-YH5)**2+(ZO1-ZH5)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS25=SQRT((XO2-XH5)**2+(YO2-YH5)**2+(ZO2-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             DIS=MIN(DIS,DIS25)
             IF(DIS.LT.DIS2)THEN
               NSDCAR(ICAR)=NSDCAR(ICAR)+1
               NAMSDC(1,ICAR,NSDCAR(ICAR))='ARG'
               NUMSDC(1,ICAR,NSDCAR(ICAR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
C
C     -- THERE ARE POSSIBLY 2 H-BONDS BETWEEN ARG AND CARBOXYL
               IF((DIS11.LT.2.2 .AND. DIS22.LT.2.2).OR.
     $            (DIS21.LT.2.2 .AND. DIS12.LT.2.2).OR.
     $            (DIS13.LT.2.2 .AND. DIS24.LT.2.2).OR.
     $            (DIS23.LT.2.2 .AND. DIS14.LT.2.2)    )THEN
                 VALSDC(1,ICAR,NSDCAR(ICAR))=-2.40
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C
C        -- 3. FIND BACKBONE HYDROGEN BONDING --
C
         FBKB=-1.20
         DIS1=2.0
         DIS2=3.5
         DO I=2,NPRTON
           XP=XPRTON(I)
           YP=YPRTON(I)
           ZP=ZPRTON(I)
           XN=XNITRN(I)
           YN=YNITRN(I)
           ZN=ZNITRN(I)
           IF((ABS(XO1-XP).LT.DIS2.AND.
     $         ABS(YO1-YP).LT.DIS2.AND.
     $         ABS(ZO1-ZP).LT.DIS2    ) .OR.
     $        (ABS(XO2-XP).LT.DIS2.AND.
     $         ABS(YO2-YP).LT.DIS2.AND.
     $         ABS(ZO2-ZP).LT.DIS2    )     ) THEN
             DISO1P=SQRT((XO1-XP)**2+(YO1-YP)**2+(ZO1-ZP)**2)
             DISO2P=SQRT((XO2-XP)**2+(YO2-YP)**2+(ZO2-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISO1P.LT.DISO2P) THEN
               DIS=DISO1P
               XVPO1=-(XP-XO1)/DISO1P
               YVPO1=-(YP-YO1)/DISO1P
               ZVPO1=-(ZP-ZO1)/DISO1P
               AGPO=XVNP*XVPO1 + YVNP*YVPO1 + ZVNP*ZVPO1
             ENDIF
             IF(DISO2P.LT.DISO1P) THEN
               DIS=DISO2P
               XVPO2=-(XP-XO2)/DISO2P
               YVPO2=-(YP-YO2)/DISO2P
               ZVPO2=-(ZP-ZO2)/DISO2P
               AGPO=XVNP*XVPO2 + YVNP*YVPO2 + ZVNP*ZVPO2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
               NBKCAR(ICAR)=NBKCAR(ICAR)+1
               NAMBKB(1,ICAR,NBKCAR(ICAR))=NAMPRT(I)
               NUMBKB(1,ICAR,NBKCAR(ICAR))=NUMPRT(I)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALBKB(1,ICAR,NBKCAR(ICAR))=FBKB*VALUE*AGPO
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALBKB(1,ICAR,NBKCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C        -- 4. IONIZABLE INTERACTION --
C
C
C           - FIND CYS(-)
C
         FCOUL=+2.40
         DIS1=4.00
         DIS2=7.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED')THEN
           IF(NMASS(1,ICAR)+NMASS(3,ICYS).GT.900 .OR.
     $     (NMASS(1,ICAR).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XO-XS).LT.DIS2.AND.
     $         ABS(YO-YS).LT.DIS2.AND.
     $         ABS(ZO-ZS).LT.DIS2    )     ) THEN
             DIS=SQRT((XO-XS)**2+(YO-YS)**2+(ZO-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))=NAMCAR(ICAR)
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND TYR(-)
C
C   Carboxyl-Tyr pairs is one exception as concerns the charge-charge 
C   interactions. The contribution is included no matter wether they
C   are defined as "buried" or "surface".
C   For more details see 
C   Very Fast Empirical Prediction and Rationalization of Protein pKa Values
C   Hui Li, Andrew D. Robertson and Jan H. Jensen
C   PROTEINS: Structure, Function, and Bioinformatics 61:704721 (2005)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO ITYR=1,NTYR
C          IF(NMASS(1,ICAR)+NMASS(4,ITYR).GT.900 .OR.
C    $     (NMASS(1,ICAR).GT.400.AND.NMASS(4,ITYR).GT.400))THEN
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF(ABS(XO-XOH).LT.DIS2.AND.
     $        ABS(YO-YOH).LT.DIS2.AND.
     $        ABS(ZO-ZOH).LT.DIS2    ) THEN
             DIS=SQRT((XO-XOH)**2+(YO-YOH)**2+(ZO-ZOH)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))=NAMCAR(ICAR)
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
             END IF
           END IF
C          END IF
         END DO
C
C Coulombic interaction between ASP/GLU residue and either LYS(+) or ARG(+)
C are in the non-iterative part because Lys and Arg are assumed to be  
C always positively charged when carboxyl groups titrate                    
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(1,ICAR)+NMASS(5,ILYS).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XO-XN).LT.DIS2.AND.
     $        ABS(YO-YN).LT.DIS2.AND.
     $        ABS(ZO-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XO-XN)**2+(YO-YN)**2+(ZO-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCAR(ICAR)=NCLCAR(ICAR)+1
               NAMCOL(1,ICAR,NCLCAR(ICAR))='LYS'
               IF(ILYS.EQ.1)NAMCOL(1,ICAR,NCLCAR(ICAR))='N+ '
               NUMCOL(1,ICAR,NCLCAR(ICAR))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
C
               NCLLYS(ILYS)=NCLLYS(ILYS)+1
               NAMCOL(5,ILYS,NCLLYS(ILYS))=NAMCAR(ICAR)
               NUMCOL(5,ILYS,NCLLYS(ILYS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(5,ILYS,NCLLYS(ILYS))=-FCOUL*MIN(1.0,VALUE)
               PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(1,ICAR)+NMASS(6,IARG).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF (ABS(XO-XCZ).LT.DIS2.AND.
     $         ABS(YO-YCZ).LT.DIS2.AND.
     $         ABS(ZO-ZCZ).LT.DIS2    )THEN
             DIS=SQRT((XO-XCZ)**2+(YO-YCZ)**2+(ZO-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCAR(ICAR)=NCLCAR(ICAR)+1
               NAMCOL(1,ICAR,NCLCAR(ICAR))='ARG'
               NUMCOL(1,ICAR,NCLCAR(ICAR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
C
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))=NAMCAR(ICAR)
               NUMCOL(6,IARG,NCLARG(IARG))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=-FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END DO
C
C        -- 5. LIGAND INTERACTION --
C
C           -- H-BONDING --
C
C           - FIND O3 group -
C
C        Identify the Oh group according to the number
C        of the linked atom
C        ie atom within a distance lower than 2.0 angstr�ms
C
       DO IO3=1,NLIGO3
         XO3=X(LLIGO3(IO3))
         YO3=Y(LLIGO3(IO3))
         ZO3=Z(LLIGO3(IO3))
         NB_neighbours=0
         DO I=1,NLIGAND
          IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGO3(IO3))) THEN
           Xb=X(LLIGAND(I))
           Yb=Y(LLIGAND(I))
           Zb=Z(LLIGAND(I))
           DISn=SQRT((Xb-XO3)**2+(Yb-YO3)**2+(Zb-ZO3)**2)
           IF (DISn.LT.2) THEN
              NB_neighbours=NB_neighbours+1
           ENDIF
          ENDIF
         END DO
C
C           - FIND OH group -
C
         IF (NB_neighbours.EQ.1) THEN
           FOH=-0.80
           DIS1=3.00
           DIS2=4.00
           IF((ABS(XO1-XO3).LT.DIS2 .AND.
     $         ABS(YO1-YO3).LT.DIS2 .AND.
     $         ABS(ZO1-ZO3).LT.DIS2     ) .OR.
     $        (ABS(XO2-XO3).LT.DIS2 .AND.
     $         ABS(YO2-YO3).LT.DIS2 .AND.
     $         ABS(ZO2-ZO3).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XO3)**2+(YO1-YO3)**2+(ZO1-ZO3)**2)
             DISO2O=SQRT((XO2-XO3)**2+(YO2-YO3)**2+(ZO2-ZO3)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2) THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FOH*MIN(1.0,VALUE)
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGO3(IO3))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGO3(IO3))  
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             END IF
           END IF
         ENDIF
C
C
C           - FIND O3 (non Oh) group -
C
         IF (NB_neighbours.EQ.2) THEN
           FOO=1.60
           DIS1=3.00
           DIS2=4.00
           IF((ABS(XO1-XO3).LT.DIS2 .AND.
     $         ABS(YO1-YO3).LT.DIS2 .AND.
     $         ABS(ZO1-ZO3).LT.DIS2     ) .OR.
     $        (ABS(XO2-XO3).LT.DIS2 .AND.
     $         ABS(YO2-YO3).LT.DIS2 .AND.
     $         ABS(ZO2-ZO3).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XO3)**2+(YO1-YO3)**2+(ZO1-ZO3)**2)
             DISO2O=SQRT((XO2-XO3)**2+(YO2-YO3)**2+(ZO2-ZO3)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2) THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FOO*MIN(1.0,VALUE)
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGO3(IO3))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGO3(IO3))
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             END IF
           END IF
         ENDIF
       END DO
C
C
C           - FIND Cl group -
C
           FOCl=1.60 
           DIS1=3.00
           DIS2=4.00
           DO ICl=1,NLIGCl
            XCl=X(LLIGCl(ICl))
            YCl=Y(LLIGCl(ICl))
            ZCl=Z(LLIGCl(ICl))
            IF((ABS(XO1-XCl).LT.DIS2 .AND.
     $          ABS(YO1-YCl).LT.DIS2 .AND.
     $          ABS(ZO1-ZCl).LT.DIS2     ) .OR.
     $         (ABS(XO2-XCl).LT.DIS2 .AND.
     $          ABS(YO2-YCl).LT.DIS2 .AND.
     $          ABS(ZO2-ZCl).LT.DIS2     )     ) THEN
             DISO1Cl=SQRT((XO1-XCl)**2+(YO1-YCl)**2+(ZO1-ZCl)**2)
             DISO2Cl=SQRT((XO2-XCl)**2+(YO2-YCl)**2+(ZO2-ZCl)**2)
             DIS=MIN(DISO1Cl,DISO2Cl)
             IF(DIS.LT.DIS2) THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FOCl*MIN(1.0,VALUE)
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGCl(ICl))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGCl(ICl))
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             END IF
            END IF
           END DO
C
C
C           - FIND F group -
C
           FOF=1.60 
           DIS1=3.00
           DIS2=4.00
           DO IF=1,NLIGF
            XF=X(LLIGF(IF))
            YF=Y(LLIGF(IF))
            ZF=Z(LLIGF(IF))
            IF((ABS(XO1-XF).LT.DIS2 .AND.
     $          ABS(YO1-YF).LT.DIS2 .AND.
     $          ABS(ZO1-ZF).LT.DIS2     ) .OR.
     $         (ABS(XO2-XF).LT.DIS2 .AND.
     $          ABS(YO2-YF).LT.DIS2 .AND.
     $          ABS(ZO2-ZF).LT.DIS2     )     ) THEN
             DISO1F=SQRT((XO1-XF)**2+(YO1-YF)**2+(ZO1-ZF)**2)
             DISO2F=SQRT((XO2-XF)**2+(YO2-YF)**2+(ZO2-ZF)**2)
             DIS=MIN(DISO1F,DISO2F)
             IF(DIS.LT.DIS2) THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FOF*MIN(1.0,VALUE)
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGF(IF))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGF(IF))
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             END IF
            END IF
           END DO
C
C
C           - FIND Amide group (Nam) -
C
         FNH=-1.20
         DIS1=2.0
         DIS2=3.5
         DO I=1,NLIGNam
          IF(NAMATM(LLIGNam(I)).NE.'  C3 ') THEN
           XN=X(LLIGNam(I))
           YN=Y(LLIGNam(I))
           ZN=Z(LLIGNam(I))
           XP=XHLIG(I)
           YP=YHLIG(I)
           ZP=ZHLIG(I)
           IF((ABS(XO1-XP).LT.DIS2.AND.
     $         ABS(YO1-YP).LT.DIS2.AND.
     $         ABS(ZO1-ZP).LT.DIS2    ) .OR.
     $        (ABS(XO2-XP).LT.DIS2.AND.
     $         ABS(YO2-YP).LT.DIS2.AND.
     $         ABS(ZO2-ZP).LT.DIS2    )     ) THEN
             DISO1P=SQRT((XO1-XP)**2+(YO1-YP)**2+(ZO1-ZP)**2)
             DISO2P=SQRT((XO2-XP)**2+(YO2-YP)**2+(ZO2-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISO1P.LT.DISO2P) THEN
               DIS=DISO1P
               XVPO1=-(XP-XO1)/DISO1P
               YVPO1=-(YP-YO1)/DISO1P
               ZVPO1=-(ZP-ZO1)/DISO1P
               AGPO=XVNP*XVPO1 + YVNP*YVPO1 + ZVNP*ZVPO1
             ENDIF
             IF(DISO2P.LT.DISO1P) THEN
               DIS=DISO2P
               XVPO2=-(XP-XO2)/DISO2P
               YVPO2=-(YP-YO2)/DISO2P
               ZVPO2=-(ZP-ZO2)/DISO2P
               AGPO=XVNP*XVPO2 + YVNP*YVPO2 + ZVNP*ZVPO2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*VALUE*AGPO
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGNam(I))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGNam(I))
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             END IF
           END IF
          ENDIF
         END DO
C
C
C           - FIND Aromatic Nitrogen (Nar) -
C
C
         FNg=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ILIG=1,NLIGNar
           XP=XgLIG(ILIG)
           YP=YgLIG(ILIG)
           ZP=ZgLIG(ILIG)
           IF((ABS(XO1-XP).LT.DIS2.AND.
     $         ABS(YO1-YP).LT.DIS2.AND.
     $         ABS(ZO1-ZP).LT.DIS2    ) .OR.
     $        (ABS(XO2-XP).LT.DIS2.AND.
     $         ABS(YO2-YP).LT.DIS2.AND.
     $         ABS(ZO2-ZP).LT.DIS2    )     ) THEN
             DISO1P=SQRT((XO1-XP)**2+(YO1-YP)**2+(ZO1-ZP)**2)
             DISO2P=SQRT((XO2-XP)**2+(YO2-YP)**2+(ZO2-ZP)**2)
             DIS=MIN(DISO1P,DISO2P)
             IF(DIS.LT.DIS2)THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FNg*VALUE
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGNar(ILIG))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGNar(ILIG))
C
               NSDNar(ILIG)=NSDNar(ILIG)+1
               NAMSDC(8,ILIG,NSDNar(ILIG))=NAMCAR(ICAR)
               NUMSDC(8,ILIG,NSDNar(ILIG))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(8,ILIG,NSDNar(ILIG))=-FNg*VALUE 
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
               IF(NMASS(1,ICAR)+NMASS(8,ILIG).GT.900 .OR.
     $     (NMASS(1,ICAR).GT.400.AND.NMASS(8,ILIG).GT.400))THEN
                 VALLIG(1,ICAR,NLGCAR(ICAR))=-1.60
                 VALSDC(8,ILIG,NSDNar(ILIG))=+1.60
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
               PK1LGNar(ILIG)=PK1LGNar(ILIG)
     $                       +VALSDC(8,ILIG,NSDNar(ILIG))
             END IF
           END IF
         END DO
C
C
C           - FIND C=O - 
C
         FOO=0.8
         DIS1=3.0
         DIS2=4.0
         DO I=1,NLIGO2 
           XOd=X(LLIGO2(I))
           YOd=Y(LLIGO2(I))
           ZOd=Z(LLIGO2(I))
C          DIST=1000.0
C          DO J=1,NLIGAND
C           IF(NAMATM(J).EQ.'  So2 '.OR.NAMATM(J).EQ.'  C2  ') THEN 
C            Xb=X(LLIGAND(J))
C            Yb=Y(LLIGAND(J))
C            Zb=Z(LLIGAND(J))
C            DIS=SQRT((Xb-XOd)**2+(Yb-YOd)**2+(Zb-ZOd)**2)
C            DIST=MIN(DIST,DIS)
C             IF(DIS.EQ.DIST) THEN
C              Xn=Xb
C              Yn=Yb
C              Zn=Zb
C             ENDIF
C           ENDIF
C          END DO
           IF((ABS(XO1-XOd).LT.DIS2.AND.
     $         ABS(YO1-YOd).LT.DIS2.AND.
     $         ABS(ZO1-ZOd).LT.DIS2    ) .OR.
     $        (ABS(XO2-XOd).LT.DIS2.AND.
     $         ABS(YO2-YOd).LT.DIS2.AND.
     $         ABS(ZO2-ZOd).LT.DIS2    )     ) THEN
             DISO1O=SQRT((XO1-XOd)**2+(YO1-YOd)**2+(ZO1-ZOd)**2)
             DISO2O=SQRT((XO2-XOd)**2+(YO2-YOd)**2+(ZO2-ZOd)**2)
             DIS=MIN(DISO1O,DISO2O)
C            VECNRM=SQRT((XOd-Xn)**2+(YOd-Yn)**2+(ZOd-Zn)**2)
C            XVCO=(XOd-Xn)/VECNRM
C            YVCO=(YOd-Yn)/VECNRM
C            ZVCO=(ZOd-Zn)/VECNRM
C            IF(DISO1O.LT.DISO2O) THEN
C              DIS=DISO1O
C              XVOO1=-(XOd-XO1)/DISO1O
C              YVOO1=-(YOd-YO1)/DISO1O
C              ZVOO1=-(ZOd-ZO1)/DISO1O
C              AGOO=XVCO*XVOO1 + YVCO*YVOO1 + ZVCO*ZVOO1
C            ENDIF
C            IF(DISO2O.LT.DISO1O) THEN
C              DIS=DISO2O
C              XVOO2=-(XOd-XO2)/DISO2O
C              YVOO2=-(YOd-YO2)/DISO2O
C              ZVOO2=-(ZOd-ZO2)/DISO2O
C              AGOO=XVCO*XVOO2 + YVCO*YVOO2 + ZVCO*ZVOO2
C            ENDIF
C            IF(DIS.LT.DIS2 .AND. AGOO.GT.0.001)THEN
             IF(DIS.LT.DIS2)THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FOO*VALUE
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGO2(I))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGO2(I))
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             END IF
           END IF
         END DO
C
C
C
C           - FIND Npl -        
C
         FNH=-0.80
         DIS1=2.00
         DIS2=4.00
         DO IC2=1,NLIGC2
C
C  -- GUADININIUM TYPE GROUP (eg. 1K1O.pdb) --
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg1P(IC2)
           YH1=YNg1P(IC2)
           ZH1=ZNg1P(IC2)
           XH2=XNg2P2(IC2)  
           YH2=YNg2P2(IC2)
           ZH2=ZNg2P2(IC2)
           XH3=XNg2P1(IC2)
           YH3=YNg2P1(IC2)
           ZH3=ZNg2P1(IC2)
           XH4=XNg3P2(IC2)
           YH4=YNg3P2(IC2)
           ZH4=ZNg3P2(IC2)
           XH5=XNg3P1(IC2)
           YH5=YNg3P1(IC2)
           ZH5=ZNg3P1(IC2)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH5).LT.DIS2.AND.
     $         ABS(YO1-YH5).LT.DIS2.AND.
     $         ABS(ZO1-ZH5).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH5).LT.DIS2.AND.
     $         ABS(YO2-YH5).LT.DIS2.AND.
     $         ABS(ZO2-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS15=SQRT((XO1-XH5)**2+(YO1-YH5)**2+(ZO1-ZH5)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS25=SQRT((XO2-XH5)**2+(YO2-YH5)**2+(ZO2-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             DIS=MIN(DIS,DIS25)
             IF(DIS.LT.DIS2)THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*MIN(1.0,VALUE)
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGC2(IC2))
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGC2(IC2))

                 NSDcg(Ic2)=NSDcg(Ic2)+1
                 NAMSDC(9,Ic2,NSDcg(Ic2))=NAMCAR(ICAR)
                 NUMSDC(9,Ic2,NSDcg(Ic2))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(9,Ic2,NSDcg(Ic2))=-FNH*MIN(1.0,VALUE)
C
C
C     -- THERE ARE POSSIBLY 2 H-BONDS BETWEEN GUADININIUM GROUP AND CARBOXYL
               IF((DIS11.LT.2.2 .AND. DIS22.LT.2.2).OR.
     $            (DIS21.LT.2.2 .AND. DIS12.LT.2.2).OR.
     $            (DIS13.LT.2.2 .AND. DIS24.LT.2.2).OR.
     $            (DIS23.LT.2.2 .AND. DIS14.LT.2.2)    )THEN
               VALLIG(1,ICAR,NLGCAR(ICAR))=-2.40
               VALSDC(9,IC2,NSDCg(IC2))=2.40
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALSDC(9,IC2,NsdCg(IC2))
             END IF
           END IF
          ENDIF
C
C  -- Np1 pair GROUP (eg. 1K1I.pdb) --
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg2P2(IC2)  
           YH1=YNg2P2(IC2)
           ZH1=ZNg2P2(IC2)
           XH2=XNg2P1(IC2)
           YH2=YNg2P1(IC2)
           ZH2=ZNg2P1(IC2)
           XH3=XNg3P2(IC2)
           YH3=YNg3P2(IC2)
           ZH3=ZNg3P2(IC2)
           XH4=XNg3P1(IC2)
           YH4=YNg3P1(IC2)
           ZH4=ZNg3P1(IC2)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             IF(DIS.LT.DIS2)THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*MIN(1.0,VALUE)
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGC2(IC2))
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGC2(IC2))

                 NSDcn2(Ic2)=NSDcn2(Ic2)+1
                 NAMSDC(9,Ic2,NSDcn2(Ic2))=NAMCAR(ICAR)
                 NUMSDC(9,Ic2,NSDcn2(Ic2))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(9,Ic2,NSDcn2(Ic2))=-FNH*MIN(1.0,VALUE)
C
C
C     -- THERE ARE POSSIBLY 2 H-BONDS BETWEEN C2N GROUP AND CARBOXYL
               IF((DIS12.LT.2.2 .AND. DIS23.LT.2.2).OR.
     $            (DIS22.LT.2.2 .AND. DIS13.LT.2.2)    )THEN
               VALLIG(1,ICAR,NLGCAR(ICAR))=-2.40
               VALSDC(9,IC2,NSDCn2(IC2))=2.40
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NsdCn2(IC2))
             END IF
           END IF
          ENDIF
         END DO
C
C
C  -- NH2 GROUP (eg. 1HPV.pdb) --
C
C     This groups is assumed to be non ionizable
C
         FNH=-0.80
         DIS1=2.00
         DIS2=4.00
         DO IN=1,NLIGNpl
          IF(NAMATM(LLIGNpl(IN)).EQ.'  Np1'.AND.
     $       NAMLN(IN).NE.'  Cg '.AND.NAMLN(IN).NE.'  CN2') THEN
           XN=X(LLIGNpl(IN))
           YN=Y(LLIGNpl(IN))
           ZN=Z(LLIGNpl(IN))
           XH1=XP1Np1(IN)
           YH1=YP1Np1(IN)
           ZH1=ZP1Np1(IN)
           XH2=XP2Np1(IN)
           YH2=YP2Np1(IN)
           ZH2=ZP2Np1(IN)
            IF(
     $         (ABS(XO1-XH1).LT.DIS2.AND.
     $          ABS(YO1-YH1).LT.DIS2.AND.
     $          ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $         (ABS(XO2-XH1).LT.DIS2.AND.
     $          ABS(YO2-YH1).LT.DIS2.AND.
     $          ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $         (ABS(XO1-XH2).LT.DIS2.AND.
     $          ABS(YO1-YH2).LT.DIS2.AND.
     $          ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $         (ABS(XO2-XH2).LT.DIS2.AND.
     $          ABS(YO2-YH2).LT.DIS2.AND.
     $          ABS(ZO2-ZH2).LT.DIS2    )     ) THEN
              DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
              DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
              DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
              DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
              DIS=MIN(DIS11,DIS12)
              DIS=MIN(DIS,DIS21)
              DIS=MIN(DIS,DIS22)
             IF(DIS.LT.DIS2)THEN  
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*MIN(1.0,VALUE)
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGNpl(IN))
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGNpl(IN))
c dmr inserted line below for PK1CAR
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             ENDIF
            ENDIF
          ENDIF
         END DO
C
C           - FIND N3 atoms -
C
         FNH=-0.80
         DIS1=2.0
         DIS2=3.5
         DO IN3=1,NLIGN3 
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N30')then
           FNH=-0.80
           DIS1=3.00
           DIS2=4.00
             XN=X(LLIGN3(IN3))
             YN=Y(LLIGN3(IN3))
             ZN=Z(LLIGN3(IN3))
             IF((ABS(XO1-XN).LT.DIS2.AND.
     $           ABS(YO1-YN).LT.DIS2.AND.
     $           ABS(ZO1-ZN).LT.DIS2    ) .OR.
     $          (ABS(XO2-XN).LT.DIS2.AND.
     $           ABS(YO2-YN).LT.DIS2.AND.
     $           ABS(ZO2-ZN).LT.DIS2    )     ) THEN
               DISO1N=SQRT((XO1-XN)**2+(YO1-YN)**2+(ZO1-ZN)**2)
               DISO2N=SQRT((XO2-XN)**2+(YO2-YN)**2+(ZO2-ZN)**2)
               DIS=MIN(DISO1N,DISO2N)
               IF(DIS.LT.DIS2)THEN
                 NLGCAR(ICAR)=NLGCAR(ICAR)+1
                 NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGN3(IN3))
                 LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*MIN(1.0,VALUE)
                 PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
C
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))=NAMCAR(ICAR)
                 NUMSDC(7,IN3,NSDN3(IN3))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
             END IF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N31')then
           FNH=-0.80
           DIS1=3.00
           DIS2=4.00
             XN=X(LLIGN3(IN3))
             YN=Y(LLIGN3(IN3))
             ZN=Z(LLIGN3(IN3))
             IF((ABS(XO1-XN).LT.DIS2.AND.
     $           ABS(YO1-YN).LT.DIS2.AND.
     $           ABS(ZO1-ZN).LT.DIS2    ) .OR.
     $          (ABS(XO2-XN).LT.DIS2.AND.
     $           ABS(YO2-YN).LT.DIS2.AND.
     $           ABS(ZO2-ZN).LT.DIS2    )     ) THEN
               DISO1N=SQRT((XO1-XN)**2+(YO1-YN)**2+(ZO1-ZN)**2)
               DISO2N=SQRT((XO2-XN)**2+(YO2-YN)**2+(ZO2-ZN)**2)
               DIS=MIN(DISO1N,DISO2N)
               IF(DIS.LT.DIS2)THEN
                 NLGCAR(ICAR)=NLGCAR(ICAR)+1
                 NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGN3(IN3))
                 LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*MIN(1.0,VALUE)
                 PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
C
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))=NAMCAR(ICAR)
                 NUMSDC(7,IN3,NSDN3(IN3))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
             END IF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N32') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XH1=XN3P1(IN3)
           YH1=YN3P1(IN3)
           ZH1=ZN3P1(IN3)
           XH2=XN3P2(IN3)
           YH2=YN3P2(IN3)
           ZH2=ZN3P2(IN3)
           IF((ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             IF (DIS11.LT.DIS12) THEN
               XP=XH1
               YP=YH1
               ZP=ZH1
               DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
               IF (DIS21.LT.DIS11) THEN
                 XP=XH1
                 YP=YH1
                 ZP=ZH1
                 DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
                 IF (DIS22.LT.DIS21) THEN
                   XP=XH2
                   YP=YH2
                   ZP=ZH2
                 ENDIF
               ENDIF
               DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
               IF (DIS22.LT.DIS11) THEN
                 XP=XH2
                 YP=YH2
                 ZP=ZH2
               ENDIF
             ELSE
               XP=XH2
               YP=YH2
               ZP=ZH2
               DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
               IF (DIS21.LT.DIS12) THEN
                 XP=XH1
                 YP=YH1
                 ZP=ZH1
                 DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
                 IF (DIS22.LT.DIS21) THEN
                   XP=XH2
                   YP=YH2
                   ZP=ZH2
                 ENDIF
               ENDIF
               DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
               IF (DIS22.LT.DIS12) THEN
                 XP=XH2
                 YP=YH2
                 ZP=ZH2
              ENDIF
             ENDIF
             DISO1P=SQRT((XO1-XP)**2+(YO1-YP)**2+(ZO1-ZP)**2)
             DISO2P=SQRT((XO2-XP)**2+(YO2-YP)**2+(ZO2-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISO1P.LT.DISO2P) THEN
               DIS=DISO1P
               XVPO1=-(XP-XO1)/DISO1P
               YVPO1=-(YP-YO1)/DISO1P
               ZVPO1=-(ZP-ZO1)/DISO1P
               AGPO=XVNP*XVPO1 + YVNP*YVPO1 + ZVNP*ZVPO1
             ENDIF
             IF(DISO2P.LT.DISO1P) THEN
               DIS=DISO2P  
               XVPO2=-(XP-XO2)/DISO2P
               YVPO2=-(YP-YO2)/DISO2P
               ZVPO2=-(ZP-ZO2)/DISO2P
               AGPO=XVNP*XVPO2 + YVNP*YVPO2 + ZVNP*ZVPO2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*VALUE*AGPO
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGN3(IN3))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGN3(IN3))
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
C
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))=NAMCAR(ICAR)
               NUMSDC(7,IN3,NSDN3(IN3))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)*AGPO
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             ENDIF
           ENDIF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N33') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XP=XN3P(IN3)
           YP=YN3P(IN3)
           ZP=ZN3P(IN3)
           IF((ABS(XO1-XP).LT.DIS2.AND.
     $         ABS(YO1-YP).LT.DIS2.AND.
     $         ABS(ZO1-ZP).LT.DIS2    ) .OR.
     $        (ABS(XO2-XP).LT.DIS2.AND.
     $         ABS(YO2-YP).LT.DIS2.AND.
     $         ABS(ZO2-ZP).LT.DIS2    )     ) THEN
             DISO1P=SQRT((XO1-XP)**2+(YO1-YP)**2+(ZO1-ZP)**2)
             DISO2P=SQRT((XO2-XP)**2+(YO2-YP)**2+(ZO2-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISO1P.LT.DISO2P) THEN
               DIS=DISO1P
               XVPO1=-(XP-XO1)/DISO1P
               YVPO1=-(YP-YO1)/DISO1P
               ZVPO1=-(ZP-ZO1)/DISO1P
               AGPO=XVNP*XVPO1 + YVNP*YVPO1 + ZVNP*ZVPO1
             ENDIF
             IF(DISO2P.LT.DISO1P) THEN
               DIS=DISO2P
               XVPO2=-(XP-XO2)/DISO2P
               YVPO2=-(YP-YO2)/DISO2P
               ZVPO2=-(ZP-ZO2)/DISO2P
               AGPO=XVNP*XVPO2 + YVNP*YVPO2 + ZVNP*ZVPO2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
                NLGCAR(ICAR)=NLGCAR(ICAR)+1
                VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                VALUE=MIN(1.0,VALUE)
                VALLIG(1,ICAR,NLGCAR(ICAR))=FNH*VALUE*AGPO
                LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGN3(IN3))
                NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGN3(IN3))
                PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
C
                NSDN3(IN3)=NSDN3(IN3)+1
                NAMSDC(7,IN3,NSDN3(IN3))=NAMCAR(ICAR)
                NUMSDC(7,IN3,NSDN3(IN3))=LCARRS(ICAR)
                VALSDC(7,IN3,NSDN3(IN3))=-FNH*VALUE*AGPO
                PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
              ENDIF
            ENDIF
           ENDIF
C
         END DO
C
C
C
C
C           - FIND N1 atoms -
C
C           The ghost atom is used to take into
C           account the orientation of the lone pair
C
         FNG=-1.20
         DIS1=2.0
         DIS2=3.5
         DO IN1=1,NLIGN1
           XG=XN1G(IN1)
           YG=YN1G(IN1)
           ZG=ZN1G(IN1)
           XN=X(LLIGN1(IN1))
           YN=Y(LLIGN1(IN1))
           ZN=Z(LLIGN1(IN1))
           IF((ABS(XO1-XG).LT.DIS2.AND.
     $         ABS(YO1-YG).LT.DIS2.AND.
     $         ABS(ZO1-ZG).LT.DIS2    ) .OR.
     $        (ABS(XO2-XG).LT.DIS2.AND.
     $         ABS(YO2-YG).LT.DIS2.AND.
     $         ABS(ZO2-ZG).LT.DIS2    )     ) THEN
             DISO1G=SQRT((XO1-XG)**2+(YO1-YG)**2+(ZO1-ZG)**2)
             DISO2G=SQRT((XO2-XG)**2+(YO2-YG)**2+(ZO2-ZG)**2)
             VECNRM=SQRT((XG-XN)**2+(YG-YN)**2+(ZG-ZN)**2)
             XVNG=(XG-XN)/VECNRM
             YVNG=(YG-YN)/VECNRM
             ZVNG=(ZG-ZN)/VECNRM
             IF(DISO1G.LT.DISO2G) THEN
               DIS=DISO1G
               XVGO1=-(XG-XO1)/DISO1G
               YVGO1=-(YG-YO1)/DISO1G
               ZVGO1=-(ZG-ZO1)/DISO1G
               AGGO=XVNG*XVGO1 + YVNG*YVGO1 + ZVNG*ZVGO1
             ENDIF
             IF(DISO2G.LT.DISO1G) THEN
               DIS=DISO2G
               XVGO2=-(XG-XO2)/DISO2G
               YVGO2=-(YG-YO2)/DISO2G
               ZVGO2=-(ZG-ZO2)/DISO2G
               AGGO=XVNG*XVGO2 + YVNG*YVGO2 + ZVNG*ZVGO2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGGO.GT.0.001)THEN
               NLGCAR(ICAR)=NLGCAR(ICAR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(1,ICAR,NLGCAR(ICAR))=FNG*VALUE*AGGO
               LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGN1(IN1))
               NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGN1(IN1))
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
             ENDIF
           ENDIF
         END DO
C
C
C           -- IONIZABLE INTERACTION --
C
C           N3+, Npl+ are considered in the iterative part (STEP 8)
C
C
C           - FIND Charged Atom -
C
      call charatm(1,icar,xo,yo,zo,ncllcar,namlcol,namres,
     $ nblcol,nbatm,vallcol,pk1car,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
      END DO
C
C
C
C      **********************
C      STEP 5. HIS
C      **********************
C
C
       DO IHIS=1, NHIS
         XH1=XHISP1(IHIS)
         YH1=YHISP1(IHIS)
         ZH1=ZHISP1(IHIS)
         XH2=XHISP2(IHIS)
         YH2=YHISP2(IHIS)
         ZH2=ZHISP2(IHIS)
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XCT=(XCG+XND+XCE+XNE+XCD)/5.0
         YCT=(YCG+YND+YCE+YNE+YCD)/5.0
         ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
C
C
C           - FIND CYS-HIS H-BONDING
C
C        IT SEEMS THAT CYS-S...H...N-HIS SYSTEM
C        TENDS TO BE CYS-S(-)...(+)H-N-HIS
C        PK(CYS) < PK(HIS)
C
         FSN=1.60
         DIS1=3.00
         DIS2=4.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED')THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XH1-XS).LT.DIS2.AND.
     $         ABS(YH1-YS).LT.DIS2.AND.
     $         ABS(ZH1-ZS).LT.DIS2    ) .OR.
     $        (ABS(XH2-XS).LT.DIS2.AND.
     $         ABS(YH2-YS).LT.DIS2.AND.
     $         ABS(ZH2-ZS).LT.DIS2    )     ) THEN
             DISH1S=SQRT((XH1-XS)**2+(YH1-YS)**2+(ZH1-ZS)**2)
             DISH2S=SQRT((XH2-XS)**2+(YH2-YS)**2+(ZH2-ZS)**2)
             DIS=MIN(DISH1S,DISH2S)
             IF(DIS.LT.DIS2)THEN
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))='CYS'
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=FSN*MIN(1.0,VALUE)
C
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='HIS'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=-FSN*MIN(1.0,VALUE)
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
               IF(NMASS(2,IHIS)+NMASS(3,ICYS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
                 VALSDC(2,IHIS,NSDHIS(IHIS))=+3.60
                 VALSDC(3,ICYS,NSDCYS(ICYS))=-3.60
               END IF
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END IF
         END DO
C
C
C
C           - FIND TYR-OH -
C             HIS-NH...(-)O-TYR DECREASES TYR'S PK
C
         FOH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ITYR=1,NTYR
           XO=X(LTYROH(ITYR))
           YO=Y(LTYROH(ITYR))
           ZO=Z(LTYROH(ITYR))
           IF((ABS(XH1-XO).LT.DIS2.AND.
     $         ABS(YH1-YO).LT.DIS2.AND.
     $         ABS(ZH1-ZO).LT.DIS2    ) .OR.
     $        (ABS(XH2-XO).LT.DIS2.AND.
     $         ABS(YH2-YO).LT.DIS2.AND.
     $         ABS(ZH2-ZO).LT.DIS2    )     ) THEN
             DISH1O=SQRT((XH1-XO)**2+(YH1-YO)**2+(ZH1-ZO)**2)
             DISH2O=SQRT((XH2-XO)**2+(YH2-YO)**2+(ZH2-ZO)**2)
             DIS=MIN(DISH1O,DISH2O)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='HIS'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND ASN -
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF((ABS(XH1-XASN).LT.DIS2.AND.
     $         ABS(YH1-YASN).LT.DIS2.AND.
     $         ABS(ZH1-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XASN).LT.DIS2.AND.
     $         ABS(YH2-YASN).LT.DIS2.AND.
     $         ABS(ZH2-ZASN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XASN)**2+(YH1-YASN)**2+(ZH1-ZASN)**2)
             DISH2=SQRT((XH2-XASN)**2+(YH2-YASN)**2+(ZH2-ZASN)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))='ASN'
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=FOH*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF((ABS(XH1-XGLN).LT.DIS2.AND.
     $         ABS(YH1-YGLN).LT.DIS2.AND.
     $         ABS(ZH1-ZGLN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XGLN).LT.DIS2.AND.
     $         ABS(YH2-YGLN).LT.DIS2.AND.
     $         ABS(ZH2-ZGLN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XGLN)**2+(YH1-YGLN)**2+(ZH1-ZGLN)**2)
             DISH2=SQRT((XH2-XGLN)**2+(YH2-YGLN)**2+(ZH2-ZGLN)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDHIS(IHIS)=NSDHIS(IHIS)+1
               NAMSDC(2,IHIS,NSDHIS(IHIS))='GLN'
               NUMSDC(2,IHIS,NSDHIS(IHIS))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(2,IHIS,NSDHIS(IHIS))=FOH*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALSDC(2,IHIS,NSDHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C        -- 3. FIND BACKBONE INTERACTION --
C
         FBKB=1.20
         DIS1=2.00
         DIS2=3.50
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(XH1-XBKO).LT.DIS2.AND.
     $         ABS(YH1-YBKO).LT.DIS2.AND.
     $         ABS(ZH1-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH2-XBKO).LT.DIS2.AND.
     $         ABS(YH2-YBKO).LT.DIS2.AND.
     $         ABS(ZH2-ZBKO).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XBKO)**2+(YH1-YBKO)**2+(ZH1-ZBKO)**2)
             DISH2=SQRT((XH2-XBKO)**2+(YH2-YBKO)**2+(ZH2-ZBKO)**2)
             VECNM=SQRT((XBKC-XBKO)**2+(YBKC-YBKO)**2+(ZBKC-ZBKO)**2)
             XVCO=(XBKO-XBKC)/VECNM
             YVCO=(YBKO-YBKC)/VECNM
             ZVCO=(ZBKO-ZBKC)/VECNM
             IF(DISH1.LT.DISH2) THEN
               DIS=DISH1
               XVOH1=(XH1-XBKO)/DISH1
               YVOH1=(YH1-YBKO)/DISH1
               ZVOH1=(ZH1-ZBKO)/DISH1
               AGOH=XVCO*XVOH1+YVCO*YVOH1+ZVCO*ZVOH1
             ENDIF
             IF(DISH2.LT.DISH1) THEN
               DIS=DISH2
               XVOH2=(XH2-XBKO)/DISH2
               YVOH2=(YH2-YBKO)/DISH2
               ZVOH2=(ZH2-ZBKO)/DISH2
               AGOH=XVCO*XVOH2+YVCO*YVOH2+ZVCO*ZVOH2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NBKHIS(IHIS)=NBKHIS(IHIS)+1
               NAMBKB(2,IHIS,NBKHIS(IHIS))=NAMPRT(I-1)
               NUMBKB(2,IHIS,NBKHIS(IHIS))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(2,IHIS,NBKHIS(IHIS))=FBKB*AGOH*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALBKB(2,IHIS,NBKHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C        -- 4. IONIZABLE INTERACTION   --
C
C Coulombic interaction between HIS residue and either LYS(+) or ARG(+) 
C are in the non-iterative part because Lys and Arg are assumed to be  
C always positively charged when His groups titrate                   
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(2,IHIS)+NMASS(5,ILYS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF((ABS(XCT-XN).LT.DIS2.AND.
     $         ABS(YCT-YN).LT.DIS2.AND.
     $         ABS(ZCT-ZN).LT.DIS2    )     ) THEN
             DIS=SQRT((XCT-XN)**2+(YCT-YN)**2+(ZCT-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='LYS'
               IF(ILYS.EQ.1)NAMCOL(2,IHIS,NCLHIS(IHIS))='N+ '
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(2,IHIS)+NMASS(6,IARG).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XCT-XCZ).LT.DIS2.AND.
     $        ABS(YCT-YCZ).LT.DIS2.AND.
     $        ABS(ZCT-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XCT-XCZ)**2+(YCT-YCZ)**2+(ZCT-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='ARG'
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END DO
C
C
C
C           -- H-BONDING --
C
C           - FIND O3 group -
C
C        Identify the Oh group according to the number
C        of the linked atom , ie atoms within a distance lower than
C        2.0 angstr�ms)
C
C        NB: AT THIS POINT WE IGNORE THE OH GROUPS FROM THE LIGAND
C
       DO IO3=1,NLIGO3
         XO3=X(LLIGO3(IO3))
         YO3=Y(LLIGO3(IO3))
         ZO3=Z(LLIGO3(IO3))
         NB_neighbours=0
         DO I=1,NLIGAND
          IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGO3(IO3))) THEN
           Xb=X(LLIGAND(I))
           Yb=Y(LLIGAND(I))
           Zb=Z(LLIGAND(I))
           DISn=SQRT((Xb-XO3)**2+(Yb-YO3)**2+(Zb-ZO3)**2)
           IF (DISn.LT.2) THEN
              NB_neighbours=NB_neighbours+1
           ENDIF
          ENDIF
         END DO
C
C
C           - FIND O3 (non Oh) group -
C
         IF (NB_neighbours.EQ.2) THEN
           FNO=1.60
           DIS1=3.00
           DIS2=4.00
           IF((ABS(XH1-XO3).LT.DIS2 .AND.
     $         ABS(YH1-YO3).LT.DIS2 .AND.
     $         ABS(ZH1-ZO3).LT.DIS2     ) .OR.
     $        (ABS(XH2-XO3).LT.DIS2 .AND.
     $         ABS(YH2-YO3).LT.DIS2 .AND.
     $         ABS(ZH2-ZO3).LT.DIS2     )     ) THEN
             DISH1O=SQRT((XH1-XO3)**2+(YH1-YO3)**2+(ZH1-ZO3)**2)
             DISH2O=SQRT((XH2-XO3)**2+(YH2-YO3)**2+(ZH2-ZO3)**2)
             IF(DISH1O.LT.DISH2O) THEN
               DIS=DISH1O
               VECNM1=SQRT((XND-XH1)**2+(YND-YH1)**2+(ZND-ZH1)**2)
               VECNM2=SQRT((XNE-XH2)**2+(YNE-YH2)**2+(ZNE-ZH2)**2)
               XVNH1=(XH1-XND)/VECNM1
               YVNH1=(YH1-YND)/VECNM1
               ZVNH1=(ZH1-ZND)/VECNM1
               XVH1O=(XO3-XH1)/DISH1O
               YVH1O=(YO3-YH1)/DISH1O
               ZVH1O=(ZO3-ZH1)/DISH1O
               AGOH=XVNH1*XVH1O+YVNH1*YVH1O+ZVNH1*ZVH1O
             ENDIF
             IF(DISH2O.LT.DISH1O) THEN
               DIS=DISH2O
               VECNM2=SQRT((XNE-XH2)**2+(YNE-YH2)**2+(ZNE-ZH2)**2)
               XVNH2=(XH2-XNE)/VECNM2
               YVNH2=(YH2-YNE)/VECNM2
               ZVNH2=(ZH2-ZNE)/VECNM2
               XVH2O=(XO3-XH2)/DISH2O
               YVH2O=(YO3-YH2)/DISH2O
               ZVH2O=(ZO3-ZH2)/DISH2O
               AGOH2=XVNH2*XVH2O+YVNH2*YVH2O+ZVNH2*ZVH2O
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNO*AGOH*MIN(1.0,VALUE)
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGO3(IO3))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGO3(IO3))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             END IF
           END IF
        ENDIF
       END DO
C
C
C           - FIND Cl group -
C
           FNCl=1.60 
           DIS1=3.00
           DIS2=4.00
           DO ICl=1,NLIGCl
            XCl=X(LLIGCl(ICl))
            YCl=Y(LLIGCl(ICl))
            ZCl=Z(LLIGCl(ICl))
           IF((ABS(XH1-XCl).LT.DIS2 .AND.
     $         ABS(YH1-YCl).LT.DIS2 .AND.
     $         ABS(ZH1-ZCl).LT.DIS2     ) .OR.
     $        (ABS(XH2-XCl).LT.DIS2 .AND.
     $         ABS(YH2-YCl).LT.DIS2 .AND.
     $         ABS(ZH2-ZCl).LT.DIS2     )     ) THEN
             DISH1Cl=SQRT((XH1-XCl)**2+(YH1-YCl)**2+(ZH1-ZCl)**2)
             DISH2Cl=SQRT((XH2-XCl)**2+(YH2-YCl)**2+(ZH2-ZCl)**2)
             IF(DISH1Cl.LT.DISH2Cl) THEN
               DIS=DISH1Cl
               VECNM1=SQRT((XND-XH1)**2+(YND-YH1)**2+(ZND-ZH1)**2)
               XVNH1=(XH1-XND)/VECNM1
               YVNH1=(YH1-YND)/VECNM1
               ZVNH1=(ZH1-ZND)/VECNM1
               XVH1Cl=(XCl-XH1)/DISH1Cl
               YVH1Cl=(YCl-YH1)/DISH1Cl
               ZVH1Cl=(ZCl-ZH1)/DISH1Cl
               AGClH=XVNH1*XVH1Cl+YVNH1*YVH1Cl+ZVNH1*ZVH1Cl
             ENDIF
             IF(DISH2Cl.LT.DISH1Cl) THEN
               DIS=DISH2Cl
               VECNM2=SQRT((XNE-XH2)**2+(YNE-YH2)**2+(ZNE-ZH2)**2)
               XVNH2=(XH2-XNE)/VECNM2
               YVNH2=(YH2-YNE)/VECNM2
               ZVNH2=(ZH2-ZNE)/VECNM2
               XVH2Cl=(XCl-XH2)/DISH2Cl
               YVH2Cl=(YCl-YH2)/DISH2Cl
               ZVH2Cl=(ZCl-ZH2)/DISH2Cl
               AGClH=XVNH2*XVH2Cl+YVNH2*YVH2Cl+ZVNH2*ZVH2Cl
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGClH.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNCl*AGClH*MIN(1.0,VALUE)
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGCl(ICl))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGCl(ICl))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             END IF
           END IF
           END DO
C
C
C           - FIND F group -
C
           FNF=1.60 
           DIS1=3.00
           DIS2=4.00
           DO IF=1,NLIGF
            XF=X(LLIGF(IF))
            YF=Y(LLIGF(IF))
            ZF=Z(LLIGF(IF))
           IF((ABS(XH1-XF).LT.DIS2 .AND.
     $         ABS(YH1-YF).LT.DIS2 .AND.
     $         ABS(ZH1-ZF).LT.DIS2     ) .OR.
     $        (ABS(XH2-XF).LT.DIS2 .AND.
     $         ABS(YH2-YF).LT.DIS2 .AND.
     $         ABS(ZH2-ZF).LT.DIS2     )     ) THEN
             DISH1F=SQRT((XH1-XF)**2+(YH1-YF)**2+(ZH1-ZF)**2)
             DISH2F=SQRT((XH2-XF)**2+(YH2-YF)**2+(ZH2-ZF)**2)
             IF(DISH1F.LT.DISH2F) THEN
               DIS=DISH1F
               VECNM1=SQRT((XND-XH1)**2+(YND-YH1)**2+(ZND-ZH1)**2)
               XVNH1=(XH1-XND)/VECNM1
               YVNH1=(YH1-YND)/VECNM1
               ZVNH1=(ZH1-ZND)/VECNM1
               XVH1F=(XF-XH1)/DISH1F
               YVH1F=(YF-YH1)/DISH1F
               ZVH1F=(ZF-ZH1)/DISH1F
               AGFH=XVNH1*XVH1F+YVNH1*YVH1F+ZVNH1*ZVH1F
             ENDIF
             IF(DISH2F.LT.DISH1F) THEN
               DIS=DISH2F
               VECNM2=SQRT((XNE-XH2)**2+(YNE-YH2)**2+(ZNE-ZH2)**2)
               XVNH2=(XH2-XNE)/VECNM2
               YVNH2=(YH2-YNE)/VECNM2
               ZVNH2=(ZH2-ZNE)/VECNM2
               XVH2F=(XF-XH2)/DISH2F
               YVH2F=(YF-YH2)/DISH2F
               ZVH2F=(ZF-ZH2)/DISH2F
               AGFH=XVNH2*XVH2F+YVNH2*YVH2F+ZVNH2*ZVH2F
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGFH.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNF*AGFH*MIN(1.0,VALUE)
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGF(IF))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGF(IF))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             END IF
           END IF
           END DO
C
C
C           - FIND Amide group (Nam) -
C
         FNH=-1.20
         DIS1=2.0
         DIS2=3.5
         DO I=1,NLIGNam
          IF(NAMATM(LLIGNam(I)).NE.'  C3 ') THEN
           XN=X(LLIGNam(I))
           YN=Y(LLIGNam(I))
           ZN=Z(LLIGNam(I))
           XP=XHLIG(I)
           YP=YHLIG(I)
           ZP=ZHLIG(I)
           IF((ABS(XNE-XP).LT.DIS2.AND.
     $         ABS(YNE-YP).LT.DIS2.AND.
     $         ABS(ZNE-ZP).LT.DIS2    ) .OR.
     $        (ABS(XND-XP).LT.DIS2.AND.
     $         ABS(YND-YP).LT.DIS2.AND.
     $         ABS(ZND-ZP).LT.DIS2    )     ) THEN
             DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
             DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISNEP.LT.DISNDP) THEN
               DIS=DISNEP
               XVPNE=-(XP-XNE)/DISNEP
               YVPNE=-(YP-YNE)/DISNEP
               ZVPNE=-(ZP-ZNE)/DISNEP
               AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
             ENDIF
             IF(DISNDP.LT.DISNEP) THEN
               DIS=DISNDP
               XVPND=-(XP-XND)/DISNDP
               YVPND=-(YP-YND)/DISNDP
               ZVPND=-(ZP-ZND)/DISNDP
               AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPN.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*VALUE*AGPN
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGNam(I))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGNam(I))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             END IF
           END IF
         ENDIF
         END DO
C
C *************************************************
C after debate in during my stay in Copenhagen
C interaction hereafter moved to hte iterative part
C *************************************************
C
C           - FIND Aromatic Nitrogen (Nar) -
C
C
C        FNg=-0.80
C        DIS1=2.00
C        DIS2=3.00
C        DO ILIG=1,NLIGNar
C          XN=X(LLIGNar(ILIG))
C          YN=Y(LLIGNar(ILIG))
C          ZN=Z(LLIGNar(ILIG))
C          XP=XgLIG(ILIG)
C          YP=YgLIG(ILIG)
C          ZP=ZgLIG(ILIG)
C          IF((ABS(XNE-XP).LT.DIS2.AND.
C    $         ABS(YNE-YP).LT.DIS2.AND.
C    $         ABS(ZNE-ZP).LT.DIS2    ) .OR.
C    $        (ABS(XND-XP).LT.DIS2.AND.
C    $         ABS(YND-YP).LT.DIS2.AND.
C    $         ABS(ZND-ZP).LT.DIS2    )     ) THEN
C            DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
C            DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
C            VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
C            XVNP=(XP-XN)/VECNRM
C            YVNP=(YP-YN)/VECNRM
C            ZVNP=(ZP-ZN)/VECNRM
C            IF(DISNEP.LT.DISNDP) THEN
C              DIS=DISNEP
C              XVPNE=-(XP-XNE)/DISNEP
C              YVPNE=-(YP-YNE)/DISNEP
C              ZVPNE=-(ZP-ZNE)/DISNEP
C              AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
C            ENDIF
C            IF(DISNDP.LT.DISNEP) THEN
C              DIS=DISNDP    
C              XVPND=-(XP-XND)/DISNDP
C              YVPND=-(YP-YND)/DISNDP
C              ZVPND=-(ZP-ZND)/DISNDP
C              AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
C            ENDIF
C            IF(DIS.LT.DIS2 .AND. AGPN.GT.0.001)THEN
C              NLGHIS(IHIS)=NLGHIS(IHIS)+1
C              VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
C              VALLIG(2,IHIS,NLGHIS(IHIS))=FNg*MIN(1.0,VALUE)*AGPN
C              LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGNar(ILIG))
C              NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGNar(ILIG))
C
C              NSDNar(ILIG)=NSDNar(ILIG)+1
C              NAMSDC(8,ILIG,NSDNar(ILIG))='HIS'
C              NUMSDC(8,ILIG,NSDNar(ILIG))=LHISRS(IHIS)
C              VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
C              VALSDC(8,ILIG,NSDNar(ILIG))=-FNg*MIN(1.0,VALUE)*AGPN
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
C              IF(NMASS(1,ICAR)+NMASS(8,ILIG).GT.900 .OR.
C    $     (NMASS(1,ICAR).GT.400.AND.NMASS(8,ILIG).GT.400))THEN
C                VALLIG(1,ICAR,NLGCAR(ICAR))=-1.60
C                VALSDC(8,ILIG,NSDNar(ILIG))=+1.60
C                NHBCAR(ICAR)=NHBCAR(ICAR) + 1
C                NHBHIS(IHIS)=NHBHIS(ILIG) + 1
C                PK1CAR(ICAR)=PK1CAR(ICAR)-6.00
C                PK1LGNar(ILIG)=PK1LGNar(ILIG)+6.00
C              END IF
C              PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
C              PK1LGNar(ILIG)=PK1LGNar(ILIG)
C    $                       +VALSDC(8,ILIG,NSDNar(ILIG))
C            END IF
C          END IF
C        END DO
C
C
C           - FIND C=O - 
C
C*****************************************************
C july 2007 we decided to include the angle dependence 
C for this interaction
C the angle considered is the one formed by O from the
C carbonyl group and the NH atoms from the His resiude
C*****************************************************
C
         FON=0.8
         DIS1=3.0
         DIS2=4.0
         DO I=1,NLIGO2 
           XOd=X(LLIGO2(I))
           YOd=Y(LLIGO2(I))
           ZOd=Z(LLIGO2(I))
           IF((ABS(XH1-XOd).LT.DIS2.AND.
     $         ABS(YH1-YOd).LT.DIS2.AND.
     $         ABS(ZH1-ZOd).LT.DIS2    ) .OR.
     $        (ABS(XH2-XOd).LT.DIS2.AND.
     $         ABS(YH2-YOd).LT.DIS2.AND.
     $         ABS(ZH2-ZOd).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XOd)**2+(YH1-YOd)**2+(ZH1-ZOd)**2)
             DISH2=SQRT((XH2-XOd)**2+(YH2-YOd)**2+(ZH2-ZOd)**2)
             DIS=MIN(DISH1,DISH2)
C
C NB : We need to identify the closest NH group from the His
C residue to evaluate the angle NHO (either NDH1 or NEH2)
C
            IF(DISH1.LT.DISH2)THEN
             VECNM=SQRT((XH1-XND)**2+(YH1-YND)**2+(ZH1-ZND)**2)
             XNDH1=(XH1-XND)/VECNM
             YNDH1=(YH1-YND)/VECNM
             ZNDH1=(ZH1-ZND)/VECNM
C
             XVH1O=-(XH1-XOd)/DISH1
             YVH1O=-(YH1-YOd)/DISH1
             ZVH1O=-(ZH1-ZOd)/DISH1
             AGOH=XNDH1*XVH1O+YNDH1*YVH1O+ZNDH1*ZVH1O
            ELSE 
             VECNM=SQRT((XH2-XNE)**2+(YH2-YNE)**2+(ZH2-ZNE)**2)
             XNEH2=(XH2-XNE)/VECNM
             YNEH2=(YH2-YNE)/VECNM
             ZNEH2=(ZH2-ZNE)/VECNM
C
             XVH2O=-(XH2-XOd)/DISH2
             YVH2O=-(YH2-YOd)/DISH2
             ZVH2O=-(ZH2-ZOd)/DISH2
             AGOH=XNEH2*XVH2O+YNEH2*YVH2O+ZNEH2*ZVH2O
            ENDIF
            IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FON*MIN(1.0,VALUE)*AGOH
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGO2(I))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGO2(I))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             END IF
           END IF
         END DO
C
C
C
C           - FIND Npl -        
C
         FNH=-0.80
         DIS1=2.00
         DIS2=4.00
         DO IC2=1,NLIGC2
C
C  -- GUADININIUM TYPE GROUP (eg. 1K1O.pdb) --
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XCgN1=XNg1(IC2)
           YCgN1=YNg1(IC2)
           ZCgN1=ZNg1(IC2)
           XCgN2=XNg2(IC2)
           YCgN2=YNg2(IC2)
           ZCgN2=ZNg2(IC2)
           XCgN3=XNg3(IC2)
           YCgN3=YNg3(IC2)
           ZCgN3=ZNg3(IC2)
           XCgH1=XNg1P(IC2)
           YCgH1=YNg1P(IC2)
           ZCgH1=ZNg1P(IC2)
           XCgH2=XNg2P2(IC2)  
           YCgH2=YNg2P2(IC2)
           ZCgH2=ZNg2P2(IC2)
           XCgH3=XNg2P1(IC2)
           YCgH3=YNg2P1(IC2)
           ZCgH3=ZNg2P1(IC2)
           XCgH4=XNg3P2(IC2)
           YCgH4=YNg3P2(IC2)
           ZCgH4=ZNg3P2(IC2)
           XCgH5=XNg3P1(IC2)
           YCgH5=YNg3P1(IC2)
           ZCgH5=ZNg3P1(IC2)
           IF(
     $        (ABS(XNE-XCgH1).LT.DIS2.AND.
     $         ABS(YNE-YCgH1).LT.DIS2.AND.
     $         ABS(ZNE-ZCgH1).LT.DIS2    ) .OR.
     $        (ABS(XND-XCgH1).LT.DIS2.AND.
     $         ABS(YND-YCgH1).LT.DIS2.AND.
     $         ABS(ZND-ZCgH1).LT.DIS2    ) .OR.
     $        (ABS(XNE-XCgH2).LT.DIS2.AND.
     $         ABS(YNE-YCgH2).LT.DIS2.AND.
     $         ABS(ZNE-ZCgH2).LT.DIS2    ) .OR.
     $        (ABS(XND-XCgH2).LT.DIS2.AND.
     $         ABS(YND-YCgH2).LT.DIS2.AND.
     $         ABS(ZND-ZCgH2).LT.DIS2    ) .OR.
     $        (ABS(XNE-XCgH3).LT.DIS2.AND.
     $         ABS(YNE-YCgH3).LT.DIS2.AND.
     $         ABS(ZNE-ZCgH3).LT.DIS2    ) .OR.
     $        (ABS(XND-XCgH3).LT.DIS2.AND.
     $         ABS(YND-YCgH3).LT.DIS2.AND.
     $         ABS(ZND-ZCgH3).LT.DIS2    ) .OR.
     $        (ABS(XNE-XCgH4).LT.DIS2.AND.
     $         ABS(YNE-YCgH4).LT.DIS2.AND.
     $         ABS(ZNE-ZCgH4).LT.DIS2    ) .OR.
     $        (ABS(XND-XCgH4).LT.DIS2.AND.
     $         ABS(YND-YCgH4).LT.DIS2.AND.
     $         ABS(ZND-ZCgH4).LT.DIS2    ) .OR.
     $        (ABS(XNE-XCgH5).LT.DIS2.AND.
     $         ABS(YNE-YCgH5).LT.DIS2.AND.
     $         ABS(ZNE-ZCgH5).LT.DIS2    ) .OR.
     $        (ABS(XND-XCgH5).LT.DIS2.AND.
     $         ABS(YND-YCgH5).LT.DIS2.AND.
     $         ABS(ZND-ZCgH5).LT.DIS2    )     ) THEN
             DISNE1=SQRT((XNE-XCgH1)**2+(YNE-YCgH1)**2+(ZNE-ZCgH1)**2)
             DISNE2=SQRT((XNE-XCgH2)**2+(YNE-YCgH2)**2+(ZNE-ZCgH2)**2)
             IF (DISNE1.LT.DISNE2) THEN
               XP1=XCgH1
               YP1=YCgH1
               ZP1=ZCgH1
               XN1=XCgN1
               YN1=YCgN1
               ZN1=ZCgN1
               DISNE3=SQRT((XNE-XCgH3)**2+(YNE-YCgH3)**2+(ZNE-ZCgH3)**2)
               IF (DISNE3.LT.DISNE1) THEN
                 XP1=XCgH3
                 YP1=YCgH3
                 ZP1=ZCgH3
                 XN1=XCgN2
                 YN1=YCgN2
                 ZN1=ZCgN2
                 DISNE4=SQRT((XNE-XCgH4)**2+(YNE-YCgH4)**2
     $                      +(ZNE-ZCgH4)**2)
                 IF (DISNE4.LT.DISNE3) THEN
                   XP1=XCgH4
                   YP1=YCgH4
                   ZP1=ZCgH4
                   XN1=XCgN3
                   YN1=YCgN3
                   ZN1=ZCgN3
                   DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2
     $                        +(ZNE-ZCgH5)**2)
                   IF (DISNE5.LT.DISNE4) THEN
                     XP1=XCgH5
                     YP1=YCgH5
                     ZP1=ZCgH5
                     XN1=XCgN3
                     YN1=YCgN3
                     ZN1=ZCgN3
                   ENDIF
                 ENDIF
                 DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2
     $                      +(ZNE-ZCgH5)**2)
                 IF (DISNE5.LT.DISNE4) THEN
                   XP1=XCgH5
                   YP1=YCgH5
                   ZP1=ZCgH5
                   XN1=XCgN3
                   YN1=YCgN3
                   ZN1=ZCgN3
                 ENDIF
               ENDIF
               DISNE4=SQRT((XNE-XCgH4)**2+(YNE-YCgH4)**2
     $                    +(ZNE-ZCgH4)**2)
               IF (DISNE4.LT.DISNE3) THEN
                 XP1=XCgH4
                 YP1=YCgH4
                 ZP1=ZCgH4
                 XN1=XCgN3
                 YN1=YCgN3
                 ZN1=ZCgN3
                 DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2
     $                      +(ZNE-ZCgH5)**2)
                 IF (DISNE5.LT.DISNE4) THEN
                   XP1=XCgH5
                   YP1=YCgH5
                   ZP1=ZCgH5
                   XN1=XCgN3
                   YN1=YCgN3
                   ZN1=ZCgN3
                 ENDIF
               ENDIF
               DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2
     $                    +(ZNE-ZCgH5)**2)
               IF (DISNE5.LT.DISNE4) THEN
                 XP1=XCgH5
                 YP1=YCgH5
                 ZP1=ZCgH5
                 XN1=XCgN3
                 YN1=YCgN3
                 ZN1=ZCgN3
               ENDIF
             ELSE
               XP1=XCgH2
               YP1=YCgH2
               ZP1=ZCgH2
               XN1=XCgN2
               YN1=YCgN2
               ZN1=ZCgN2
               DISNE3=SQRT((XNE-XCgH3)**2+(YNE-YCgH3)**2
     $                    +(ZNE-ZCgH3)**2)
               IF (DISNE3.LT.DISNE2) THEN
                 XP1=XCgH3
                 YP1=YCgH3
                 ZP1=ZCgH3
                 XN1=XCgN2
                 YN1=YCgN2
                 ZN1=ZCgN2
                 DISNE4=SQRT((XNE-XCgH4)**2+(YNE-YCgH4)**2
     $                      +(ZNE-ZCgH4)**2)
                 IF (DISNE4.LT.DISNE3) THEN
                   XP1=XCgH4
                   YP1=YCgH4
                   ZP1=ZCgH4
                   XN1=XCgN3
                   YN1=YCgN3
                   ZN1=ZCgN3
                   DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2
     $                        +(ZNE-ZCgH5)**2)
                   IF (DISNE5.LT.DISNE4) THEN
                     XP1=XCgH5
                     YP1=YCgH5
                     ZP1=ZCgH5
                     XN1=XCgN3
                     YN1=YCgN3
                     ZN1=ZCgN3
                   ENDIF
                 ENDIF
                 DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2
     $                      +(ZNE-ZCgH5)**2)
                 IF (DISNE5.LT.DISNE4) THEN
                   XP1=XCgH5
                   YP1=YCgH5
                   ZP1=ZCgH5
                   XN1=XCgN3
                   YN1=YCgN3
                   ZN1=ZCgN3
                 ENDIF
               ENDIF
               DISNE4=SQRT((XNE-XCgH4)**2+(YNE-YCgH4)**2
     $                    +(ZNE-ZCgH4)**2)
               IF (DISNE4.LT.DISNE3) THEN
                 XP1=XCgH4
                 YP1=YCgH4
                 ZP1=ZCgH4
                 XN1=XCgN3
                 YN1=YCgN3
                 ZN1=ZCgN3
                 DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2
     $                      +(ZNE-ZCgH5)**2)
                 IF (DISNE5.LT.DISNE4) THEN
                   XP1=XCgH5
                   YP1=YCgH5
                   ZP1=ZCgH5
                   XN1=XCgN3
                   YN1=YCgN3
                   ZN1=ZCgN3
                 ENDIF
               ENDIF
               DISNE5=SQRT((XNE-XCgH5)**2+(YNE-YCgH5)**2+(ZNE-ZCgH5)**2)
               IF (DISNE5.LT.DISNE4) THEN
                 XP1=XCgH5
                 YP1=YCgH5
                 ZP1=ZCgH5
                 XN1=XCgN3
                 YN1=YCgN3
                 ZN1=ZCgN3
               ENDIF
             ENDIF
C
             IF (DISND1.LT.DISND2) THEN
               XP2=XCgH1
               YP2=YCgH1
               ZP2=ZCgH1
               XN2=XCgN1
               YN2=YCgN1
               ZN2=ZCgN1
               DISND3=SQRT((XND-XCgH3)**2+(YND-YCgH3)**2+(ZND-ZCgH3)**2)
               IF (DISND3.LT.DISND1) THEN
                 XP2=XCgH3
                 YP2=YCgH3
                 ZP2=ZCgH3
                 XN2=XCgN2
                 YN2=YCgN2
                 ZN2=ZCgN2
                 DISND4=SQRT((XND-XCgH4)**2+(YND-YCgH4)**2
     $                      +(ZND-ZCgH4)**2)
                 IF (DISND4.LT.DISND3) THEN
                   XP2=XCgH4
                   YP2=YCgH4
                   ZP2=ZCgH4
                   XN2=XCgN3
                   YN2=YCgN3
                   ZN2=ZCgN3
                   DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2
     $                        +(ZND-ZCgH5)**2)
                   IF (DISND5.LT.DISND4) THEN
                     XP2=XCgH5
                     YP2=YCgH5
                     ZP2=ZCgH5
                     XN2=XCgN3
                     YN2=YCgN3
                     ZN2=ZCgN3
                   ENDIF
                 ENDIF
                 DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2
     $                      +(ZND-ZCgH5)**2)
                 IF (DISND5.LT.DISND4) THEN
                   XP2=XCgH5
                   YP2=YCgH5
                   ZP2=ZCgH5
                   XN2=XCgN3
                   YN2=YCgN3
                   ZN2=ZCgN3
                 ENDIF
               ENDIF
               DISND4=SQRT((XND-XCgH4)**2+(YND-YCgH4)**2+(ZND-ZCgH4)**2)
               IF (DISND4.LT.DISND3) THEN
                 XP2=XCgH4
                 YP2=YCgH4
                 ZP2=ZCgH4
                 XN2=XCgN3
                 YN2=YCgN3
                 ZN2=ZCgN3
                 DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2
     $                      +(ZND-ZCgH5)**2)
                 IF (DISND5.LT.DISND4) THEN
                   XP2=XCgH5
                   YP2=YCgH5
                   ZP2=ZCgH5
                   XN2=XCgN3
                   YN2=YCgN3
                   ZN2=ZCgN3
                 ENDIF
               ENDIF
               DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2+(ZND-ZCgH5)**2)
               IF (DISND5.LT.DISND4) THEN
                 XP2=XCgH5
                 YP2=YCgH5
                 ZP2=ZCgH5
                 XN2=XCgN3
                 YN2=YCgN3
                 ZN2=ZCgN3
               ENDIF
             ELSE
               XP2=XCgH2
               YP2=YCgH2
               ZP2=ZCgH2
               XN2=XCgN2
               YN2=YCgN2
               ZN2=ZCgN2
               DISND3=SQRT((XND-XCgH3)**2+(YND-YCgH3)**2+(ZND-ZCgH3)**2)
               IF (DISND3.LT.DISND2) THEN
                 XP2=XCgH3
                YP2=YCgH3
                 ZP2=ZCgH3
                 XN2=XCgN2
                 YN2=YCgN2
                 ZN2=ZCgN2
                 DISND4=SQRT((XND-XCgH4)**2+(YND-YCgH4)**2
     $                      +(ZND-ZCgH4)**2)
                 IF (DISND4.LT.DISND3) THEN
                   XP2=XCgH4
                   YP2=YCgH4
                   ZP2=ZCgH4
                   XN2=XCgN3
                   YN2=YCgN3
                   ZN2=ZCgN3
                   DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2
     $                        +(ZND-ZCgH5)**2)
                   IF (DISND5.LT.DISND4) THEN
                     XP2=XCgH5
                     YP2=YCgH5
                     ZP2=ZCgH5
                     XN2=XCgN3
                     YN2=YCgN3
                     ZN2=ZCgN3
                   ENDIF
                 ENDIF
                 DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2
     $                      +(ZND-ZCgH5)**2)
                 IF (DISND5.LT.DISND4) THEN
                   XP2=XCgH5
                   YP2=YCgH5
                   ZP2=ZCgH5
                   XN2=XCgN3
                   YN2=YCgN3
                   ZN2=ZCgN3
                 ENDIF
               ENDIF
               DISND4=SQRT((XND-XCgH4)**2+(YND-YCgH4)**2+(ZND-ZCgH4)**2)
               IF (DISND4.LT.DISND3) THEN
                 XP2=XCgH4
                 YP2=YCgH4
                 ZP2=ZCgH4
                 XN2=XCgN3
                 YN2=YCgN3
                 ZN2=ZCgN3
                 DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2
     $                      +(ZND-ZCgH5)**2)
                 IF (DISND5.LT.DISND4) THEN
                   XP2=XCgH5
                   YP2=YCgH5
                   ZP2=ZCgH5
                   XN2=XCgN3
                   YN2=YCgN3
                   ZN2=ZCgN3
                 ENDIF
               ENDIF
               DISND5=SQRT((XND-XCgH5)**2+(YND-YCgH5)**2+(ZND-ZCgH5)**2)
               IF (DISND5.LT.DISND4) THEN
                 XP2=XCgH5
                 YP2=YCgH5
                 ZP2=ZCgH5
                 XN2=XCgN3
                 YN2=YCgN3
                 ZN2=ZCgN3
               ENDIF
             ENDIF
             DISNEP1=SQRT((XNE-XP1)**2+(YNE-YP1)**2+(ZNE-ZP1)**2)
             DISNEP2=SQRT((XNE-XP2)**2+(YNE-YP2)**2+(ZNE-ZP2)**2)
             IF (DISNEP1.LT.DISNEP2) THEN
               XP=XP1
               YP=YP1
               ZP=ZP1
               XN=XN1
               YN=YN1
               ZN=ZN1
               DISNDP1=SQRT((XND-XP1)**2+(YND-YP1)**2+(ZND-ZP1)**2)
               IF (DISNDP1.LT.DISNEP1) THEN
                 XP=XP1
                 YP=YP1
                 ZP=ZP1
                 XN=XN1
                 YN=YN1
                 ZN=ZN1
                 DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
                 IF (DISNDP2.LT.DISNDP1) THEN
                   XP=XP2
                   YP=YP2
                   ZP=ZP2
                   XN=XN2
                   YN=YN2
                   ZN=ZN2
                 ENDIF
               ENDIF
               DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
               IF (DISNDP2.LT.DISNEP1) THEN
                 XP=XP2
                 YP=YP2
                 ZP=ZP2
                 XN=XN2
                 YN=YN2
                 ZN=ZN2
               ENDIF
             ELSE
               XP=XP2
               YP=YP2
               ZP=ZP2
               XN=XN2
               YN=YN2
               ZN=ZN2
               DISNDP1=SQRT((XND-XP1)**2+(YND-YP1)**2+(ZND-ZP1)**2)
               IF (DISNDP1.LT.DISNEP2) THEN
                 XP=XP1
                 YP=YP1
                 ZP=ZP1
                 XN=XN1
                 YN=YN1
                 ZN=ZN1
                 DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
                 IF (DISNDP2.LT.DISNDP1) THEN
                   XP=XP2
                   YP=YP2
                   ZP=ZP2
                   XN=XN2
                   YN=YN2
                   ZN=ZN2
                 ENDIF
               ENDIF
               DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
               IF (DISNDP2.LT.DISNEP2) THEN
                 XP=XP2
                 YP=YP2
                 ZP=ZP2
                 XN=XN2
                 YN=YN2
                 ZN=ZN2
              ENDIF
             ENDIF
             DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
             DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
             VECNRMP=SQRT((XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2)
             XVNP=-(XN-XP)/VECNRMNP
             YVNP=-(YN-YP)/VECNRMNP
             ZVNP=-(ZN-ZP)/VECNRMNP
             IF(DISNEP.LT.DISNDP)THEN
               DIS=DISNEP
               XVPNE=-(XP-XNE)/DISNEP
               YVPNE=-(YP-YNE)/DISNEP
               ZVPNE=-(ZP-ZNE)/DISNEP
               AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
             ENDIF
             IF(DISNDP.LT.DISNEP)THEN
               DIS=DISNDP
               XVPND=-(XP-XND)/DISNDP
               YVPND=-(YP-YND)/DISNDP
               ZVPND=-(ZP-ZND)/DISNDP
               AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
             ENDIF
             IF(DIS.LT.DIS2.AND.AGPN.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*MIN(1.0,VALUE)*AGPN
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGC2(IC2))
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGC2(IC2))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
c dmr is below valid?
               NCLCg(IC2)=NCLCg(IC2)+1
               NAMsdc(9,IC2,NCLCg(IC2))='HIS'
               NUMsdc(9,IC2,NCLCg(IC2))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALsdc(9,IC2,NCLCg(IC2))=-Fnh*MIN(1.0,VALUE)*agpn
               PK1LGCg(IC2)=PK1LGCg(IC2)
     $                        +VALsdc(9,IC2,NCLCg(IC2))
             END IF
           END IF
          ENDIF
C
C  -- Np1 pair GROUP (eg. 1K1I.pdb) --
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XCN1=XNg2(IC2)
           YCN1=YNg2(IC2)
           ZCN1=ZNg2(IC2)
           XCN2=XNg3(IC2)
           YCN2=YNg3(IC2)
           ZCN2=ZNg3(IC2)
           XCH1=XNg2P2(IC2)  
           YCH1=YNg2P2(IC2)
           ZCH1=ZNg2P2(IC2)
           XCH2=XNg2P1(IC2)
           YCH2=YNg2P1(IC2)
           ZCH2=ZNg2P1(IC2)
           XCH3=XNg3P2(IC2)
           YCH3=YNg3P2(IC2)
           ZCH3=ZNg3P2(IC2)
           XCH4=XNg3P1(IC2)
           YCH4=YNg3P1(IC2)
           ZCH4=ZNg3P1(IC2)
           IF(
     $        (ABS(XNE-XCH1).LT.DIS2.AND.
     $         ABS(YNE-YCH1).LT.DIS2.AND.
     $         ABS(ZNE-ZCH1).LT.DIS2    ) .OR.
     $        (ABS(XND-XCH1).LT.DIS2.AND.
     $         ABS(YND-YCH1).LT.DIS2.AND.
     $         ABS(ZND-ZCH1).LT.DIS2    ) .OR.
     $        (ABS(XNE-XCH2).LT.DIS2.AND.
     $         ABS(YNE-YCH2).LT.DIS2.AND.
     $         ABS(ZNE-ZCH2).LT.DIS2    ) .OR.
     $        (ABS(XND-XCH2).LT.DIS2.AND.
     $         ABS(YND-YCH2).LT.DIS2.AND.
     $         ABS(ZND-ZCH2).LT.DIS2    ) .OR.
     $        (ABS(XNE-XCH3).LT.DIS2.AND.
     $         ABS(YNE-YCH3).LT.DIS2.AND.
     $         ABS(ZNE-ZCH3).LT.DIS2    ) .OR.
     $        (ABS(XND-XCH3).LT.DIS2.AND.
     $         ABS(YND-YCH3).LT.DIS2.AND.
     $         ABS(ZND-ZCH3).LT.DIS2    ) .OR.
     $        (ABS(XNE-XCH4).LT.DIS2.AND.
     $         ABS(YNE-YCH4).LT.DIS2.AND.
     $         ABS(ZNE-ZCH4).LT.DIS2    ) .OR.
     $        (ABS(XND-XCH4).LT.DIS2.AND.
     $         ABS(YND-YCH4).LT.DIS2.AND.
     $         ABS(ZND-ZCH4).LT.DIS2    )     ) THEN
             DISNE1=SQRT((XNE-XCH1)**2+(YNE-YCH1)**2+(ZNE-ZCH1)**2)
             DISNE2=SQRT((XNE-XCH2)**2+(YNE-YCH2)**2+(ZNE-ZCH2)**2)
             IF (DISNE1.LT.DISNE2) THEN
               XP1=XCH1
               YP1=YCH1
               ZP1=ZCH1
               XN1=XCN1
               YN1=YCN1
               ZN1=ZCN1
               DISNE3=SQRT((XNE-XCH3)**2+(YNE-YCH3)**2+(ZNE-ZCH3)**2)
               IF (DISNE3.LT.DISNE1) THEN
                 XP1=XCH3
                 YP1=YCH3
                 ZP1=ZCH3
                 XN1=XCN2
                 YN1=YCN2
                 ZN1=ZCN2
                 DISNE4=SQRT((XNE-XCH4)**2+(YNE-YCH4)**2+(ZNE-ZCH4)**2)
                 IF (DISNE4.LT.DISNE3) THEN
                   XP1=XCH4
                   YP1=YCH4
                   ZP1=ZCH4
                   XN1=XCN2
                   YN1=YCN2
                   ZN1=ZCN2
                 ENDIF
               ENDIF
               DISNE4=SQRT((XNE-XCH4)**2+(YNE-YCH4)**2+(ZNE-ZCH4)**2)
               IF (DISNE4.LT.DISNE3) THEN
                 XP1=XCH4
                 YP1=YCH4
                 ZP1=ZCH4
                 XN1=XCN2
                 YN1=YCN2
                 ZN1=ZCN2
               ENDIF
             ELSE
               XP1=XCH2
               YP1=YCH2
               ZP1=ZCH2
               XN1=XCN1
               YN1=YCN1
               ZN1=ZCN1
               DISNE3=SQRT((XNE-XCH3)**2+(YNE-YCH3)**2+(ZNE-ZCH3)**2)
               IF (DISNE3.LT.DISNE2) THEN
                 XP1=XCH3
                 YP1=YCH3
                 ZP1=ZCH3
                 XN1=XCN2
                 YN1=YCN2
                 ZN1=ZCN2
                 DISNE4=SQRT((XNE-XCH4)**2+(YNE-YCH4)**2+(ZNE-ZCH4)**2)
                 IF (DISNE4.LT.DISNE3) THEN
                   XP1=XCH4
                   YP1=YCH4
                   ZP1=ZCH4
                   XN1=XCN2
                   YN1=YCN2
                   ZN1=ZCN2
                 ENDIF
               ENDIF
               DISNE4=SQRT((XNE-XCH4)**2+(YNE-YCH4)**2+(ZNE-ZCH4)**2)
               IF (DISNE4.LT.DISNE3) THEN
                 XP1=XCH4
                 YP1=YCH4
                 ZP1=ZCH4
                 XN1=XCN2
                 YN1=YCN2
                 ZN1=ZCN2
               ENDIF
             ENDIF
             IF (DISND1.LT.DISND2) THEN
               XP2=XCH1
               YP2=YCH1
               ZP2=ZCH1
               XN2=XCN1
               YN2=YCN1
               ZN2=ZCN1
               DISND3=SQRT((XND-XCH3)**2+(YND-YCH3)**2+(ZND-ZCH3)**2)
               IF (DISND3.LT.DISND1) THEN
                 XP2=XCH3
                 YP2=YCH3
                 ZP2=ZCH3
                 XN2=XCN2
                 YN2=YCN2
                 ZN2=ZCN2
                 DISND4=SQRT((XND-XCH4)**2+(YND-YCH4)**2+(ZND-ZCH4)**2)
                 IF (DISND4.LT.DISND3) THEN
                   XP2=XCH4
                   YP2=YCH4
                   ZP2=ZCH4
                   XN2=XCN2
                   YN2=YCN2
                   ZN2=ZCN2
                 ENDIF
               ENDIF
               DISND4=SQRT((XND-XCH4)**2+(YND-YCH4)**2+(ZND-ZCH4)**2)
               IF (DISND4.LT.DISND3) THEN
                 XP2=XCH4
                 YP2=YCH4
                 ZP2=ZCH4
                 XN2=XCN2
                 YN2=YCN2
                 ZN2=ZCN2
               ENDIF
             ELSE
               XP2=XCH2
               YP2=YCH2
               ZP2=ZCH2
               XN2=XCN1
               YN2=YCN1
               ZN2=ZCN1
               DISND3=SQRT((XND-XCH3)**2+(YND-YCH3)**2+(ZND-ZCH3)**2)
               IF (DISND3.LT.DISND2) THEN
                 XP2=XCH3
                 YP2=YCH3
                 ZP2=ZCH3
                 XN2=XCN2
                 YN2=YCN2
                 ZN2=ZCN2
                 DISND4=SQRT((XND-XCH4)**2+(YND-YCH4)**2+(ZND-ZCH4)**2)
                 IF (DISND4.LT.DISND3) THEN
                   XP2=XCH4
                   YP2=YCH4
                   ZP2=ZCH4
                   XN2=XCN2
                   YN2=YCN2
                   ZN2=ZCN2
                 ENDIF
               ENDIF
               DISND4=SQRT((XND-XCH4)**2+(YND-YCH4)**2+(ZND-ZCH4)**2)
               IF (DISND4.LT.DISND3) THEN
                 XP2=XCH4
                 YP2=YCH4
                 ZP2=ZCH4
                 XN2=XCN2
                 YN2=YCN2
                 ZN2=ZCN2
               ENDIF
             ENDIF
             DISNEP1=SQRT((XNE-XP1)**2+(YNE-YP1)**2+(ZNE-ZP1)**2)
             DISNEP2=SQRT((XNE-XP2)**2+(YNE-YP2)**2+(ZNE-ZP2)**2)
             IF (DISNEP1.LT.DISNEP2) THEN
               XP=XP1
               YP=YP1
               ZP=ZP1
               XN=XN1
               YN=YN1
               ZN=ZN1
               DISNDP1=SQRT((XND-XP1)**2+(YND-YP1)**2+(ZND-ZP1)**2)
               IF (DISNDP1.LT.DISNEP1) THEN
                 XP=XP1
                 YP=YP1
                 ZP=ZP1
                 XN=XN1
                 YN=YN1
                 ZN=ZN1
                 DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
                 IF (DISNDP2.LT.DISNDP1) THEN
                   XP=XP2
                   YP=YP2
                   ZP=ZP2
                   XN=XN2
                   YN=YN2
                   ZN=ZN2
                 ENDIF
               ENDIF
               DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
               IF (DISNDP2.LT.DISNEP1) THEN
                 XP=XP2
                 YP=YP2
                 ZP=ZP2
                 XN=XN2
                 YN=YN2
                 ZN=ZN2
               ENDIF
             ELSE
               XP=XP2
               YP=YP2
               ZP=ZP2
               XN=XN2
               YN=YN2
               ZN=ZN2
               DISNDP1=SQRT((XND-XP1)**2+(YND-YP1)**2+(ZND-ZP1)**2)
               IF (DISNDP1.LT.DISNEP2) THEN
                 XP=XP1
                 YP=YP1
                 ZP=ZP1
                 XN=XN1
                 YN=YN1
                 ZN=ZN1
                 DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
                 IF (DISNDP2.LT.DISNDP1) THEN
                   XP=XP2
                   YP=YP2
                   ZP=ZP2
                   XN=XN2
                   YN=YN2
                   ZN=ZN2
                 ENDIF
               ENDIF
               DISNDP2=SQRT((XND-XP2)**2+(YND-YP2)**2+(ZND-ZP2)**2)
               IF (DISNDP2.LT.DISNEP2) THEN
                 XP=XP2
                 YP=YP2
                 ZP=ZP2
                 XN=XN2
                 YN=YN2
                 ZN=ZN2
              ENDIF
             ENDIF
             DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
             DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
             VECNRMP=SQRT((XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2)
             XVNP=-(XN-XP)/VECNRMNP
             YVNP=-(YN-YP)/VECNRMNP
             ZVNP=-(ZN-ZP)/VECNRMNP
             IF(DISNEP.LT.DISNDP)THEN
               DIS=DISNEP
               XVPNE=-(XP-XNE)/DISNEP
               YVPNE=-(YP-YNE)/DISNEP
               ZVPNE=-(ZP-ZNE)/DISNEP
               AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
             ENDIF
             IF(DISNDP.LT.DISNEP)THEN
               DIS=DISNDP
               XVPND=-(XP-XND)/DISNDP
               YVPND=-(YP-YND)/DISNDP
               ZVPND=-(ZP-ZND)/DISNDP
               AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
             ENDIF
             IF(DIS.LT.DIS2.AND.AGPN.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*MIN(1.0,VALUE)*AGPN
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGC2(IC2))
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGC2(IC2))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
c dmr is below valid?
               NCLCn2(IC2)=NCLCn2(IC2)+1
               NAMsdc(9,IC2,NCLCn2(IC2))='HIS'
               NUMsdc(9,IC2,NCLCn2(IC2))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALsdc(9,IC2,NCLCn2(IC2))=-Fnh*MIN(1.0,VALUE)*agpn
               PK1LGCn2(IC2)=PK1LGCn2(IC2)
     $                        +VALsdc(9,IC2,NCLCn2(IC2))
             END IF
           END IF
          ENDIF
         END DO
C
C
C  -- NH2 GROUP (eg. 1HPV.pdb) --
C
C     This group is assumed to be non ionizable
C
         FNH=-0.80
         DIS1=2.00
         DIS2=4.00
         DO IN=1,NLIGNpl
          IF(NAMATM(LLIGNpl(IN)).EQ.'  Np1'.AND.
     $       NAMLN(IN).NE.'  Cg '.AND.NAMLN(IN).NE.'  CN2') THEN
           XN=X(LLIGNpl(IN))
           YN=Y(LLIGNpl(IN))
           ZN=Z(LLIGNpl(IN))
           XNH1=XP1Np1(IN)
           YNH1=YP1Np1(IN)
           ZNH1=ZP1Np1(IN)
           XNH2=XP2Np1(IN)
           YNH2=YP2Np1(IN)
           ZNH2=ZP2Np1(IN)
            IF(
     $         (ABS(XNE-XNH1).LT.DIS2.AND.
     $          ABS(YNE-YNH1).LT.DIS2.AND.
     $          ABS(ZNE-ZNH1).LT.DIS2    ) .OR.
     $         (ABS(XND-XNH1).LT.DIS2.AND.
     $          ABS(YND-YNH1).LT.DIS2.AND.
     $          ABS(ZND-ZNH1).LT.DIS2    ) .OR.
     $         (ABS(XNE-XNH2).LT.DIS2.AND.
     $          ABS(YNE-YNH2).LT.DIS2.AND.
     $          ABS(ZNE-ZNH2).LT.DIS2    ) .OR.
     $         (ABS(XND-XNH2).LT.DIS2.AND.
     $          ABS(YND-YNH2).LT.DIS2.AND.
     $          ABS(ZND-ZNH2).LT.DIS2    )     ) THEN
             DISNE1=SQRT((XNE-XNH1)**2+(YNE-YNH1)**2+(ZNE-ZNH1)**2)
             DISNE2=SQRT((XNE-XNH2)**2+(YNE-YNH2)**2+(ZNE-ZNH2)**2)
             IF (DISNE1.LT.DISNE2) THEN
               XP=XNH1
               YP=YNH1
               ZP=ZNH1
               DISND1=SQRT((XND-XNH1)**2+(YND-YNH1)**2+(ZND-ZNH1)**2)
               IF (DISND1.LT.DISNE1) THEN
                 XP=XNH1
                 YP=YNH1
                 ZP=ZNH1
                 DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
                 IF (DISND2.LT.DISND1) THEN
                   XP=XNH2
                   YP=YNH2
                   ZP=ZNH2
                 ENDIF
               ENDIF
               DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
               IF (DISND2.LT.DISNE1) THEN
                 XP=XNH2
                 YP=YNH2
                 ZP=ZNH2
               ENDIF
             ELSE
               XP=XNH2
               YP=YNH2
               ZP=ZNH2
               DISND1=SQRT((XND-XNH1)**2+(YND-YNH1)**2+(ZND-ZNH1)**2)
               IF (DISND1.LT.DISNE2) THEN
                 XP=XNH1
                 YP=YNH1
                 ZP=ZNH1
                 DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
                 IF (DISND2.LT.DISND1) THEN
                   XP=XNH2
                   YP=YNH2
                   ZP=ZNH2
                 ENDIF
               ENDIF
               DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
               IF (DISND2.LT.DISNE2) THEN
                 XP=XNH2
                 YP=YNH2
                 ZP=ZNH2
              ENDIF
             ENDIF
             DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
             DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
             VECNRMP=SQRT((XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2)
             XVNP=-(XN-XP)/VECNRMNP
             YVNP=-(YN-YP)/VECNRMNP
             ZVNP=-(ZN-ZP)/VECNRMNP
             IF(DISNEP.LT.DISNDP)THEN
               DIS=DISNEP
               XVPNE=-(XP-XNE)/DISNEP
               YVPNE=-(YP-YNE)/DISNEP
               ZVPNE=-(ZP-ZNE)/DISNEP
               AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
             ENDIF
             IF(DISNDP.LT.DISNEP)THEN
               DIS=DISNDP
               XVPND=-(XP-XND)/DISNDP
               YVPND=-(YP-YND)/DISNDP
               ZVPND=-(ZP-ZND)/DISNDP
               AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
             ENDIF
             IF(DIS.LT.DIS2.AND.AGPN.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*MIN(1.0,VALUE)*AGPN
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGNpl(IN))
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGNpl(IN))
c dmr added line below
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             ENDIF
            ENDIF
          ENDIF
         END DO
C
C
C
C           - FIND N3 atoms -
C
         FNH=-0.80
         DIS1=2.0
         DIS2=3.5
         DO IN3=1,NLIGN3 
C
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N30')THEN
           FNH=-0.80
           DIS1=3.00
           DIS2=4.00
             XN=X(LLIGN3(IN3))
             YN=Y(LLIGN3(IN3))
             ZN=Z(LLIGN3(IN3))
             IF((ABS(XNE-XN).LT.DIS2.AND.
     $           ABS(YNE-YN).LT.DIS2.AND.
     $           ABS(ZNE-ZN).LT.DIS2    ) .OR.
     $          (ABS(XND-XN).LT.DIS2.AND.
     $           ABS(YND-YN).LT.DIS2.AND.
     $           ABS(ZND-ZN).LT.DIS2    )     ) THEN
               DISNEN=SQRT((XNE-XN)**2+(YNE-YN)**2+(ZNE-ZN)**2)
               DISNDN=SQRT((XND-XN)**2+(YND-YN)**2+(ZND-ZN)**2)
               DIS=MIN(DISNEN,DISNDN)
               IF(DIS.LT.DIS2)THEN
                 NLGHIS(IHIS)=NLGHIS(IHIS)+1
                 NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGN3(IN3))
                 LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*MIN(1.0,VALUE)
                 PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(1,IHIS,NLGHIS(IHIS))
 
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))='HIS'
                 NUMSDC(7,IN3,NSDN3(IN3))=LHISRS(IHIS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
             END IF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N31')THEN
           FNH=-0.80
           DIS1=3.00
           DIS2=4.00
             XN=X(LLIGN3(IN3))
             YN=Y(LLIGN3(IN3))
             ZN=Z(LLIGN3(IN3))
             IF((ABS(XNE-XN).LT.DIS2.AND.
     $           ABS(YNE-YN).LT.DIS2.AND.
     $           ABS(ZNE-ZN).LT.DIS2    ) .OR.
     $          (ABS(XND-XN).LT.DIS2.AND.
     $           ABS(YND-YN).LT.DIS2.AND.
     $           ABS(ZND-ZN).LT.DIS2    )     ) THEN
               DISNEN=SQRT((XNE-XN)**2+(YNE-YN)**2+(ZNE-ZN)**2)
               DISNDN=SQRT((XND-XN)**2+(YND-YN)**2+(ZND-ZN)**2)
               DIS=MIN(DISNEN,DISNDN)
               IF(DIS.LT.DIS2)THEN
                 NLGHIS(IHIS)=NLGHIS(IHIS)+1
                 NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGN3(IN3))
                 LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*MIN(1.0,VALUE)
                 PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(1,IHIS,NLGHIS(IHIS))
 
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))='HIS'
                 NUMSDC(7,IN3,NSDN3(IN3))=LHISRS(IHIS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
             END IF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N32') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XNH1=XN3P1(IN3)
           YNH1=YN3P1(IN3)
           ZNH1=ZN3P1(IN3)
           XNH2=XN3P2(IN3)
           YNH2=YN3P2(IN3)
           ZNH2=ZN3P2(IN3)
           IF((ABS(XNE-XNH1).LT.DIS2.AND.
     $         ABS(YNE-YNH1).LT.DIS2.AND.
     $         ABS(ZNE-ZNH1).LT.DIS2    ) .OR.
     $        (ABS(XND-XNH1).LT.DIS2.AND.
     $         ABS(YND-YNH1).LT.DIS2.AND.
     $         ABS(ZND-ZNH1).LT.DIS2    ) .OR.
     $        (ABS(XNE-XNH2).LT.DIS2.AND.
     $         ABS(YNE-YNH2).LT.DIS2.AND.
     $         ABS(ZNE-ZNH2).LT.DIS2    ) .OR.
     $        (ABS(XND-XNH2).LT.DIS2.AND.
     $         ABS(YND-YNH2).LT.DIS2.AND.
     $         ABS(ZND-ZNH2).LT.DIS2    )     ) THEN
             DISNE1=SQRT((XNE-XNH1)**2+(YNE-YNH1)**2+(ZNE-ZNH1)**2)
             DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
             IF (DISNE1.LT.DISND2) THEN
               XP=XNH1
               YP=YNH1
               ZP=ZNH1
               DISND1=SQRT((XND-XNH1)**2+(YND-YNH1)**2+(ZND-ZNH1)**2)
               IF (DISND1.LT.DISNE1) THEN
                 XP=XNH1
                 YP=YNH1
                 ZP=ZNH1
                 DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
                 IF (DISND2.LT.DISND1) THEN
                   XP=XNH2
                   YP=YNH2
                   ZP=ZNH2
                 ENDIF
               ENDIF
               DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
               IF (DISND2.LT.DISNE1) THEN
                 XP=XNH2
                 YP=YNH2
                 ZP=ZNH2
               ENDIF
             ELSE
               XP=XNH2
               YP=YNH2
               ZP=ZNH2
               DISND1=SQRT((XND-XNH1)**2+(YND-YNH1)**2+(ZND-ZNH1)**2)
               IF (DISND1.LT.DISNE2) THEN
                 XP=XNH1
                 YP=YNH1
                 ZP=ZNH1
                 DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
                 IF (DISND2.LT.DISND1) THEN
                   XP=XNH2
                   YP=YNH2
                   ZP=ZNH2
                 ENDIF
               ENDIF
               DISND2=SQRT((XND-XNH2)**2+(YND-YNH2)**2+(ZND-ZNH2)**2)
               IF (DISND2.LT.DISNE2) THEN
                 XP=XNH2
                 YP=YNH2
                 ZP=ZNH2
              ENDIF
             ENDIF
             DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
             DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISNEP.LT.DISNDP) THEN
               DIS=DISNEP
               XVPNE=-(XP-XNE)/DISNEP
               YVPNE=-(YP-YNE)/DISNEP
               ZVPNE=-(ZP-ZNE)/DISNEP
               AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
             ENDIF
             IF(DISNDP.LT.DISNEP) THEN
               DIS=DISNDP
               XVPND=-(XP-XND)/DISNDP
               YVPND=-(YP-YND)/DISNDP
               ZVPND=-(ZP-ZND)/DISNDP
               AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPN.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*VALUE*AGPN
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGN3(IN3))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGN3(IN3))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
C
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='HIS'
               NUMSDC(7,IN3,NSDN3(IN3))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)*AGPN
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             ENDIF
           ENDIF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N33') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XP=XN3P(IN3)
           YP=YN3P(IN3)
           ZP=ZN3P(IN3)
           IF((ABS(XNE-XP).LT.DIS2.AND.
     $         ABS(YNE-YP).LT.DIS2.AND.
     $         ABS(ZNE-ZP).LT.DIS2    ) .OR.
     $        (ABS(XND-XP).LT.DIS2.AND.
     $         ABS(YND-YP).LT.DIS2.AND.
     $         ABS(ZND-ZP).LT.DIS2    )     ) THEN
             DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
             DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISNEP.LT.DISNDP) THEN
               DIS=DISNEP
               XVPNE=-(XP-XNE)/DISNEP
               YVPNE=-(YP-YNE)/DISNEP
               ZVPNE=-(ZP-ZNE)/DISNEP
               AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
             ENDIF
             IF(DISNDP.LT.DISNEP) THEN
               DIS=DISNDP
               XVPND=-(XP-XND)/DISNDP
               YVPND=-(YP-YND)/DISNDP
               ZVPND=-(ZP-ZND)/DISNDP
               AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
             ENDIF 
             IF(DIS.LT.DIS2 .AND. AGPN.GT.0.001)THEN
                NLGHIS(IHIS)=NLGHIS(IHIS)+1
                VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                VALUE=MIN(1.0,VALUE)
                VALLIG(2,IHIS,NLGHIS(IHIS))=FNH*VALUE*AGPN
                LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGN3(IN3))
                NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGN3(IN3))
                PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
C
                NSDN3(IN3)=NSDN3(IN3)+1
                NAMSDC(7,IN3,NSDN3(IN3))='HIS'
                NUMSDC(7,IN3,NSDN3(IN3))=LHISRS(IHIS)
                VALSDC(7,IN3,NSDN3(IN3))=-FNH*VALUE*AGPN
                PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
              ENDIF
            ENDIF
           ENDIF
C
         END DO
C
C
C
C
C           - FIND N1 atoms -
C
C           The ghost atom is used to take into
C           account the orientation of the lone pair
C
         FNG=-1.20
         DIS1=2.0
         DIS2=3.5
         DO IN1=1,NLIGN1
           XG=XN1G(IN1)
           YG=YN1G(IN1)
           ZG=ZN1G(IN1)
           XN=X(LLIGN1(IN1))
           YN=Y(LLIGN1(IN1))
           ZN=Z(LLIGN1(IN1))
           IF((ABS(XNE-XG).LT.DIS2.AND.
     $         ABS(YNE-YG).LT.DIS2.AND.
     $         ABS(ZNE-ZG).LT.DIS2    ) .OR.
     $        (ABS(XND-XG).LT.DIS2.AND.
     $         ABS(YND-YG).LT.DIS2.AND.
     $         ABS(ZND-ZG).LT.DIS2    )     ) THEN
             DISNEG=SQRT((XNE-XG)**2+(YNE-YG)**2+(ZNE-ZG)**2)
             DISNDG=SQRT((XND-XG)**2+(YND-YG)**2+(ZND-ZG)**2)
             VECNRM=SQRT((XG-XN)**2+(YG-YN)**2+(ZG-ZN)**2)
             XVNG=(XG-XN)/VECNRM
             YVNG=(YG-YN)/VECNRM
             ZVNG=(ZG-ZN)/VECNRM
             IF(DISNEG.LT.DISNDG) THEN
               DIS=DISNEG
               XVGNE=-(XG-XNE)/DISNEG
               YVGNE=-(YG-YNE)/DISNEG
               ZVGNE=-(ZG-ZNE)/DISNEG
               AGGN=XVNG*XVGNE + YVNG*YVGNE + ZVNG*ZVGNE
             ENDIF
             IF(DISNDG.LT.DISNEG) THEN
               DIS=DISNDG
               XVGND=-(XG-XND)/DISNDG
               YVGND=-(YG-YND)/DISNDG
               ZVGND=-(ZG-ZND)/DISNDG
               AGGN=XVNG*XVGND + YVNG*YVGND + ZVNG*ZVGND
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGGN.GT.0.001)THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNG*MIN(1.0,VALUE)*AGGN
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGN1(IN1))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGN1(IN1))
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             ENDIF
           ENDIF
         END DO
C
C
C
C
C           - FIND Oco H-BONDING
C
         FNO=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ILIG=1,NLIGC2
          IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                     .AND.ZCac(ILIG).NE.0)THEN
           XO1=XLO1(ILIG)
           YO1=YLO1(ILIG)
           ZO1=ZLO1(ILIG)
           XO2=XLO2(ILIG)
           YO2=YLO2(ILIG)
           ZO2=ZLO2(ILIG)
           IF((ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DISO12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DISO21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DISO22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
              NSDLCAR(ILIG)=NSDLCAR(ILIG)+1
              NAMSDC(11,ILIG,NSDLCAR(ILIG))='HIS'
              NUMSDC(11,ILIG,NSDLCAR(ILIG))=LHISRS(IHIS)
              VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
              VALSDC(11,ILIG,NSDLCAR(ILIG))=FNO*MIN(1.0,VALUE)
C
              NLGHIS(IHIS)=NLGHIS(IHIS)+1
              LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGC2(ILIG))
              NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGC2(ILIG))
              VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
              VALLIG(2,IHIS,NLGHIS(IHIS))=-FNO*MIN(1.0,VALUE) 
C
C             STRONG H-BONDING IF BURIED
C             (MOSTLY COULOMBIC INTERACTION)
C
              IF(NMASS(11,ILIG)+NMASS(2,IHIS).GT.900 .OR.
     $     (NMASS(11,ILIG).GT.400.AND.NMASS(2,IHIS).GT.400))THEN
                 VALSDC(11,ILIG,NSDlCAR(ILIG))=-1.60
                 VALlig(2,IHIS,NlgHIS(IHIS))=+1.60
              END IF
              PK1LGCAR(ILIG)=PK1LGCAR(ILIG)
     $                      +VALSDC(11,ILIG,NSDLCAR(ILIG))
              PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             END IF
           END IF
          ENDIF
         END DO
C
C
C           -- IONIZABLE INTERACTION --
C
C           N3+, Npl+ in the iterative part STEP 8
C
C
C           - FIND Charged Atom -
C
      call charatm(2,ihis,xct,yct,zct,ncllhis,namlcol,namres,
     $ nblcol,nbatm,vallcol,pk1his,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
       END DO
C
C
C
C      **********************
C      STEP 6. CYS
C      **********************
C
C
       DO ICYS=1, NCYS
       IF(TYPCYS(ICYS).NE.'BONDED')THEN
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
C
C
C
C        -- 2. SIDECHAIN INTERACTION --
C
C
C           - FIND SER -
C
         FOH=-1.60
         DIS1=3.50
         DIS2=4.50
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XSG-XSER).LT.DIS2 .AND.
     $         ABS(YSG-YSER).LT.DIS2 .AND.
     $         ABS(ZSG-ZSER).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XSER)**2+(YSG-YSER)**2+(ZSG-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='SER'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FOH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=-1.60
         DIS1=3.50
         DIS2=4.50
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XSG-XTHR).LT.DIS2 .AND.
     $         ABS(YSG-YTHR).LT.DIS2 .AND.
     $         ABS(ZSG-ZTHR).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XTHR)**2+(YSG-YTHR)**2+(ZSG-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='THR'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FOH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FNH=-1.60
         DIS1=2.50
         DIS2=3.50
         DO IGLN=1,NGLN
           XGLN1=XGLNP1(IGLN)
           YGLN1=YGLNP1(IGLN)
           ZGLN1=ZGLNP1(IGLN)
           XGLN2=XGLNP2(IGLN)
           YGLN2=YGLNP2(IGLN)
           ZGLN2=ZGLNP2(IGLN)
           IF((ABS(XSG-XGLN1).LT.DIS2.AND.
     $         ABS(YSG-YGLN1).LT.DIS2.AND.
     $         ABS(ZSG-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XSG-XGLN2).LT.DIS2.AND.
     $         ABS(YSG-YGLN2).LT.DIS2.AND.
     $         ABS(ZSG-ZGLN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XSG-XGLN1)**2+(YSG-YGLN1)**2+(ZSG-ZGLN1)**2)
             DIS12=SQRT((XSG-XGLN2)**2+(YSG-YGLN2)**2+(ZSG-ZGLN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='GLN'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C           - FIND ASN -
C
         FNH=-1.60
         DIS1=2.50
         DIS2=3.50
         DO IASN=1,NASN
           XASN1=XASNP1(IASN)
           YASN1=YASNP1(IASN)
           ZASN1=ZASNP1(IASN)
           XASN2=XASNP2(IASN)
           YASN2=YASNP2(IASN)
           ZASN2=ZASNP2(IASN)
           IF((ABS(XSG-XASN1).LT.DIS2.AND.
     $         ABS(YSG-YASN1).LT.DIS2.AND.
     $         ABS(ZSG-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XSG-XASN2).LT.DIS2.AND.
     $         ABS(YSG-YASN2).LT.DIS2.AND.
     $         ABS(ZSG-ZASN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XSG-XASN1)**2+(YSG-YASN1)**2+(ZSG-ZASN1)**2)
             DIS12=SQRT((XSG-XASN2)**2+(YSG-YASN2)**2+(ZSG-ZASN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='ASN'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C           - FIND TRP -
C
         FNH=-1.60
         DIS1=2.50
         DIS2=3.50
         DO ITRP=1,NTRP
           XTRP1=XTRPP1(ITRP)
           YTRP1=YTRPP1(ITRP)
           ZTRP1=ZTRPP1(ITRP)
           IF((ABS(XSG-XTRP1).LT.DIS2.AND.
     $         ABS(YSG-YTRP1).LT.DIS2.AND.
     $         ABS(ZSG-ZTRP1).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XTRP1)**2+(YSG-YTRP1)**2+(ZSG-ZTRP1)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='TRP'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LTRPRS(ITRP)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C        -- 3. FIND BACKBONE INTERACTION --
C
         FBKB=-2.40
         DIS1=3.5
         DIS2=4.5
         DO I=2,NPRTON
           XP=XPRTON(I)
           YP=YPRTON(I)
           ZP=ZPRTON(I)
           XN=XNITRN(I)
           YN=YNITRN(I)
           ZN=ZNITRN(I)
           IF((ABS(XSG-XP).LT.DIS2.AND.
     $         ABS(YSG-YP).LT.DIS2.AND.
     $         ABS(ZSG-ZP).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XP)**2+(YSG-YP)**2+(ZSG-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPS=-(XP-XSG)/DIS
             YVPS=-(YP-YSG)/DIS
             ZVPS=-(ZP-ZSG)/DIS
             AGPS=XVNP*XVPS + YVNP*YVPS + ZVNP*ZVPS
             IF(DIS.LT.DIS2 .AND. AGPS.GT.0.001)THEN
               NBKCYS(ICYS)=NBKCYS(ICYS)+1
               NAMBKB(3,ICYS,NBKCYS(ICYS))=NAMPRT(I)
               NUMBKB(3,ICYS,NBKCYS(ICYS))=NUMPRT(I)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALBKB(3,ICYS,NBKCYS(ICYS))=FBKB*VALUE*AGPS
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALBKB(3,ICYS,NBKCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C        -- COULOMBIC INTERACTIONS --
C
C
C
C           - FIND CYS PAIRS -
C
         FSS=-1.6
         DIS1=3.00
         DIS2=5.00
         DO JCYS=1,NCYS
           IF(LCYSRS(ICYS).EQ.LCYSRS(JCYS)-3)THEN
             XS=X(LCYSSG(JCYS))
             YS=Y(LCYSSG(JCYS))
             ZS=Z(LCYSSG(JCYS))
             DIS=SQRT((XSG-XS)**2+(YSG-YS)**2+(ZSG-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='CYS'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LCYSRS(JCYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FSS*MIN(1.0,VALUE)
C
               NSDCYS(JCYS)=NSDCYS(JCYS)+1
               NAMSDC(3,JCYS,NSDCYS(JCYS))='CYS'
               NUMSDC(3,JCYS,NSDCYS(JCYS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,JCYS,NSDCYS(JCYS))=-FSS*MIN(1.0,VALUE)
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
               IF(NMASS(3,ICYS)+NMASS(3,JCYS).GT.900 .OR.
     $    (NMASS(3,ICYS).GT.400.AND.NMASS(3,JCYS).GT.400))THEN
                 VALSDC(3,ICYS,NSDCYS(ICYS))=-3.60
                 VALSDC(3,JCYS,NSDCYS(JCYS))=+3.60
               END IF
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
               PK1CYS(JCYS)=PK1CYS(JCYS)+VALSDC(3,JCYS,NSDCYS(JCYS))
             END IF
           END IF
         END DO
C
C
C           - FIND TYR-OH -
C
         FOH=-0.80
         DIS1=3.50
         DIS2=4.50
         DO ITYR=1,NTYR
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF(ABS(XSG-XOH).LT.DIS2.AND.
     $        ABS(YSG-YOH).LT.DIS2.AND.
     $        ABS(ZSG-ZOH).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XOH)**2+(YSG-YOH)**2+(ZSG-ZOH)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='TYR'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FOH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
C
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='CYS'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=-FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C           - FIND TYR(-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO ITYR=1,NTYR
           IF(NMASS(3,ICYS)+NMASS(4,ITYR).GT.900 .OR.
     $     (NMASS(3,ICYS).GT.400.AND.NMASS(4,ITYR).GT.400))THEN
           XOH=X(LTYROH(ITYR))
           YOH=Y(LTYROH(ITYR))
           ZOH=Z(LTYROH(ITYR))
           IF(ABS(XSG-XOH).LT.DIS2.AND.
     $        ABS(YSG-YOH).LT.DIS2.AND.
     $        ABS(ZSG-ZOH).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XOH)**2+(YSG-YOH)**2+(ZSG-ZOH)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))='CYS'
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND LYS H-BONDING
C
         FSN=-1.60
         DIS1=3.00
         DIS2=4.00
         DO ILYS=1,NLYS
           FSN=-1.60
           DIS2=4.00
           IF(ILYS.EQ.1)FSN=-2.40
           IF(ILYS.EQ.1)DIS2=4.50
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XSG-XN).LT.DIS2.AND.
     $        ABS(YSG-YN).LT.DIS2.AND.
     $        ABS(ZSG-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='LYS'
               IF(ILYS.EQ.1)NAMSDC(3,ICYS,NSDCYS(ICYS))='N+ '
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FSN*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C
C Coulombic interaction between CYS residue and LYS(+)
C is in the non-iterative part because Lys is assumed to be
C always positively charged when Cys group titrates
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(3,ICYS)+NMASS(5,ILYS).GT.900 .OR.
     $     (NMASS(3,ICYS).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XSG-XN).LT.DIS2.AND.
     $        ABS(YSG-YN).LT.DIS2.AND.
     $        ABS(ZSG-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='LYS'
               IF(ILYS.EQ.1)NAMCOL(3,ICYS,NCLCYS(ICYS))='N+ '
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LLYSRS(ILYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
C
               NCLLYS(ILYS)=NCLLYS(ILYS)+1
               NAMCOL(5,ILYS,NCLLYS(ILYS))='CYS'
               NUMCOL(5,ILYS,NCLLYS(ILYS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(5,ILYS,NCLLYS(ILYS))=-FCOUL*MIN(1.0,VALUE)
               PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND ARG(+) H-BONDING
C
         FSN=-1.60
         DIS1=2.50
         DIS2=4.00
         DO IARG=1,NARG
           XH1=XARGP1(IARG)
           YH1=YARGP1(IARG)
           ZH1=ZARGP1(IARG)
           XH2=XARGP2(IARG)
           YH2=YARGP2(IARG)
           ZH2=ZARGP2(IARG)
           XH3=XARGP3(IARG)
           YH3=YARGP3(IARG)
           ZH3=ZARGP3(IARG)
           XH4=XARGP4(IARG)
           YH4=YARGP4(IARG)
           ZH4=ZARGP4(IARG)
           XH5=XARGP5(IARG)
           YH5=YARGP5(IARG)
           ZH5=ZARGP5(IARG)
           IF(
     $        (ABS(XSG-XH1).LT.DIS2.AND.
     $         ABS(YSG-YH1).LT.DIS2.AND.
     $         ABS(ZSG-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH2).LT.DIS2.AND.
     $         ABS(YSG-YH2).LT.DIS2.AND.
     $         ABS(ZSG-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH3).LT.DIS2.AND.
     $         ABS(YSG-YH3).LT.DIS2.AND.
     $         ABS(ZSG-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH4).LT.DIS2.AND.
     $         ABS(YSG-YH4).LT.DIS2.AND.
     $         ABS(ZSG-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH5).LT.DIS2.AND.
     $         ABS(YSG-YH5).LT.DIS2.AND.
     $         ABS(ZSG-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XSG-XH1)**2+(YSG-YH1)**2+(ZSG-ZH1)**2)
             DIS12=SQRT((XSG-XH2)**2+(YSG-YH2)**2+(ZSG-ZH2)**2)
             DIS13=SQRT((XSG-XH3)**2+(YSG-YH3)**2+(ZSG-ZH3)**2)
             DIS14=SQRT((XSG-XH4)**2+(YSG-YH4)**2+(ZSG-ZH4)**2)
             DIS15=SQRT((XSG-XH5)**2+(YSG-YH5)**2+(ZSG-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             IF(DIS.LT.DIS2)THEN
               NSDCYS(ICYS)=NSDCYS(ICYS)+1
               NAMSDC(3,ICYS,NSDCYS(ICYS))='ARG'
               NUMSDC(3,ICYS,NSDCYS(ICYS))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(3,ICYS,NSDCYS(ICYS))=FSN*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NSDCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C Coulombic interaction between CYS residue and ARG(+)
C is in the non-iterative part because Arg is assumed to be
C always positively charged when Cys group titrates
C
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(3,ICYS)+NMASS(6,IARG).GT.900 .OR.
     $     (NMASS(3,ICYS).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XSG-XCZ).LT.DIS2.AND.
     $        ABS(YSG-YCZ).LT.DIS2.AND.
     $        ABS(ZSG-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XSG-XCZ)**2+(YSG-YCZ)**2+(ZSG-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='ARG'
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
C
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))='CYS'
               NUMCOL(6,IARG,NCLARG(IARG))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=-FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END DO
C
C
       END IF
C
C
C
C        -- 5. LIGAND INTERACTION --
C
C           -- H-BONDING --
C
C
C           - FIND O3 group -
C
C           We identify the Oh group according 
C           to the number of the linked atoms
C           (ie atoms within a distance lower than 
C           2.0 angstr�ms)
C
       DO IO3=1,NLIGO3
         XO3=X(LLIGO3(IO3))
         YO3=Y(LLIGO3(IO3))
         ZO3=Z(LLIGO3(IO3))
         DO I=1,NLIGAND
          IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGO3(IO3))) THEN
           Xb=X(LLIGAND(I))
           Yb=Y(LLIGAND(I))
           Zb=Z(LLIGAND(I))
           DISn=SQRT((Xb-XO3)**2+(Yb-YO3)**2+(Zb-ZO3)**2)
           IF (DISn.LT.2) THEN
              NB_neighbours=NB_neighbours+1
           ENDIF
          ENDIF
         END DO
C
C           - FIND OH group -
C
         IF (NB_neighbours.EQ.1) THEN
C
           FOH=-1.60
           DIS1=3.50
           DIS2=4.50
           IF((ABS(XSG-XO3).LT.DIS2 .AND.
     $         ABS(YSG-YO3).LT.DIS2 .AND.
     $         ABS(ZSG-ZO3).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XO3)**2+(YSG-YO3)**2+(ZSG-ZO3)**2)
             IF(DIS.LT.DIS2) THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FOH*MIN(1.0,VALUE)
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGO3(IO3))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGO3(IO3))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             END IF
           END IF
         ENDIF
C
C
C           - FIND O3 (non Oh) group -
C
         IF (NB_neighbours.EQ.2) THEN
           FOO=3.20
           DIS1=3.50
           DIS2=4.50
           IF((ABS(XSG-XO3).LT.DIS2 .AND.
     $         ABS(YSG-YO3).LT.DIS2 .AND.
     $         ABS(ZSG-ZO3).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XO3)**2+(YSG-YO3)**2+(ZSG-ZO3)**2)
             IF(DIS.LT.DIS2) THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FOO*MIN(1.0,VALUE)
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGO3(IO3))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGO3(IO3))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             END IF
           END IF
         ENDIF
C
       END DO
C
C
C           - FIND Cl group -
C
         FOCl=3.20
         DIS1=3.50
         DIS2=4.50
C
         DO ICl=1,NLIGCl
           XCl=X(LLIGCl(ICl))
           YCl=Y(LLIGCl(ICl))
           ZCl=Z(LLIGCl(ICl))
           IF((ABS(XSG-XCl).LT.DIS2 .AND.
     $         ABS(YSG-YCl).LT.DIS2 .AND.
     $         ABS(ZSG-ZCl).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XCl)**2+(YSG-YCl)**2+(ZSG-ZCl)**2)
             IF(DIS.LT.DIS2) THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FOCl*MIN(1.0,VALUE)
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGCl(ICl))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGCl(ICl))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C           - FIND F group -
C
         FOF=3.20
         DIS1=3.50
         DIS2=4.50
         DO IF=1,NLIGF
           XF=X(LLIGF(IF))
           YF=Y(LLIGF(IF))
           ZF=Z(LLIGF(IF))
           IF((ABS(XSG-XF).LT.DIS2 .AND.
     $         ABS(YSG-YF).LT.DIS2 .AND.
     $         ABS(ZSG-ZF).LT.DIS2     )     ) THEN
             DIS=SQRT((XSG-XF)**2+(YSG-YF)**2+(ZSG-ZF)**2)
             IF(DIS.LT.DIS2) THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FOF*MIN(1.0,VALUE)
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGF(IF))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGF(IF))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C           - FIND Amide group (Nam) -
C
         FNH=-2.40
         DIS1=3.5
         DIS2=4.5
         DO I=1,NLIGNam
          IF(NAMATM(LLIGNam(I)).NE.'  C3 ') THEN
           XN=X(LLIGNam(I))
           YN=Y(LLIGNam(I))
           ZN=Z(LLIGNam(I))
           XP=XHLIG(I)
           YP=YHLIG(I)
           ZP=ZHLIG(I)
           IF((ABS(XSG-XP).LT.DIS2.AND.
     $         ABS(YSG-YP).LT.DIS2.AND.
     $         ABS(ZSG-ZP).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XP)**2+(YSG-YP)**2+(ZSG-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPS=-(XP-XSG)/DIS
             YVPS=-(YP-YSG)/DIS
             ZVPS=-(ZP-ZSG)/DIS
             AGPS=XVNP*XVPS + YVNP*YVPS + ZVNP*ZVPS
             IF(DIS.LT.DIS2 .AND. AGPS.GT.0.001)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FNH*VALUE*AGPO
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGNam(I))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGNam(I))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             END IF
           END IF
          ENDIF
         END DO
C
C
C           - FIND Aromatic Nitrogen (Nar) -
C
C
         FNg=-1.60
         DIS1=3.00
         DIS2=4.00
         DO ILIG=1,NLIGNar
           XP=XgLIG(ILIG)
           YP=YgLIG(ILIG)
           ZP=ZgLIG(ILIG)
           IF((ABS(XSG-XP).LT.DIS2.AND.
     $         ABS(YSG-YP).LT.DIS2.AND.
     $         ABS(ZSG-ZP).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XP)**2+(YSG-YP)**2+(ZSG-ZP)**2)
             IF(DIS.LT.DIS2)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FNg*VALUE
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGNar(ILIG))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGNar(ILIG))
C
               NSDNar(ILIG)=NSDNar(ILIG)+1
               NAMSDC(8,ILIG,NSDNar(ILIG))='CYS'
               NUMSDC(8,ILIG,NSDNar(ILIG))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(8,ILIG,NSDNar(ILIG))=-FNg*VALUE 
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
               IF(NMASS(3,ICAR)+NMASS(8,ILIG).GT.900 .OR.
     $     (NMASS(3,ICYS).GT.400.AND.NMASS(8,ILIG).GT.400))THEN
                 VALLIG(3,ICYS,NLGCYS(ICYS))=-3.60
                 VALSDC(8,ILIG,NSDNar(ILIG))=3.60
               END IF
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
               PK1LGNar(ILIG)=PK1LGNar(ILIG)
     $                       +VALSDC(8,ILIG,NSDNar(ILIG))
             END IF
           END IF
         END DO
C
C
C
C           - FIND C=O - 
C
C
C dcb
C angle dependence removed july 2007
C dcb
C
         FOO=1.6
         DIS1=3.5
         DIS2=4.5
         DO I=1,NLIGO2 
           XO=X(LLIGO2(I))
           YO=Y(LLIGO2(I))
           ZO=Z(LLIGO2(I))
           IF((ABS(XSG-XO).LT.DIS2.AND.
     $         ABS(YSG-YO).LT.DIS2.AND.
     $         ABS(ZSG-ZO).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XO)**2+(YSG-YO)**2+(ZSG-ZO)**2)
             IF(DIS.LT.DIS2)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FOO*VALUE
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGo2(i))
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGO2(I))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             END IF
           END IF
         END DO
C
C
C 
C           - FIND Npl -        
C
         FNH=-1.60
         DIS1=2.50
         DIS2=4.00
         DO IC2=1,NLIGC2
C
C  -- GUADININIUM TYPE GROUP (eg. 1K1O.pdb) --
C
C
C dcb july 2007
C There is no angle dependence when the H-bonding interaction
C between Cys and Arg group is considered so for the moment
C we don't inlcude the angle dependence for the interaction
C between Cys residues and guadininium group from the ligand
C dcb
C
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg1P(IC2)
           YH1=YNg1P(IC2)
           ZH1=ZNg1P(IC2)
           XH2=XNg2P2(IC2)  
           YH2=YNg2P2(IC2)
           ZH2=ZNg2P2(IC2)
           XH3=XNg2P1(IC2)
           YH3=YNg2P1(IC2)
           ZH3=ZNg2P1(IC2)
           XH4=XNg3P2(IC2)
           YH4=YNg3P2(IC2)
           ZH4=ZNg3P2(IC2)
           XH5=XNg3P1(IC2)
           YH5=YNg3P1(IC2)
           ZH5=ZNg3P1(IC2)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH5).LT.DIS2.AND.
     $         ABS(YO1-YH5).LT.DIS2.AND.
     $         ABS(ZO1-ZH5).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH5).LT.DIS2.AND.
     $         ABS(YO2-YH5).LT.DIS2.AND.
     $         ABS(ZO2-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS15=SQRT((XO1-XH5)**2+(YO1-YH5)**2+(ZO1-ZH5)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS25=SQRT((XO2-XH5)**2+(YO2-YH5)**2+(ZO2-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             DIS=MIN(DIS,DIS25)
             IF(DIS.LT.DIS2)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FNH*MIN(1.0,VALUE)
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGNpl(IN))
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGNpl(IN))
                 NSDcg(Ic2)=NSDcg(Ic2)+1
                 NAMSDC(9,Ic2,NSDcg(Ic2))='CYS'
                 NUMSDC(9,Ic2,NSDcg(Ic2))=LCYSRS(ICYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(9,Ic2,NSDcg(Ic2))=-FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALSDC(9,IC2,NSDCg(IC2))
             END IF
           END IF
          ENDIF
C
C  -- Np1 pair GROUP (eg. 1K1I.pdb) --
C
C dcb july 2007
C There is no angle dependence when the H-bonding interaction
C between Cys and Arg group is considered so for the moment
C we don't include the angle dependence for the interaction
C between Cys residues and Np1 pair group from the ligand
C dcb
C
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg2P2(IC2)  
           YH1=YNg2P2(IC2)
           ZH1=ZNg2P2(IC2)
           XH2=XNg2P1(IC2)
           YH2=YNg2P1(IC2)
           ZH2=ZNg2P1(IC2)
           XH3=XNg3P2(IC2)
           YH3=YNg3P2(IC2)
           ZH3=ZNg3P2(IC2)
           XH4=XNg3P1(IC2)
           YH4=YNg3P1(IC2)
           ZH4=ZNg3P1(IC2)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             IF(DIS.LT.DIS2)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FNH*MIN(1.0,VALUE)
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGNpl(IN))
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGNpl(IN))
                 NSDcn2(Ic2)=NSDcn2(Ic2)+1
                 NAMSDC(9,Ic2,NSDcn2(Ic2))='CYS'
                 NUMSDC(9,Ic2,NSDcn2(Ic2))=LCYSRS(ICYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(9,Ic2,NSDcn2(Ic2))=-FNH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NSDCN2(IC2))
             END IF
           END IF
          ENDIF
         END DO
C
C
C  -- NH2 GROUP (eg. 1HPV.pdb) --
C
C     This groups is assumed to be non ionizable
C
         DO IN=1,NLIGNpl
           XN=X(LLIGNpl(IN))
           YN=Y(LLIGNpl(IN))
           ZN=Z(LLIGNpl(IN))
          DO IC2=1,NLIGC2
           XC2=X(LLIGC2(IC2))
           YC2=Y(LLIGC2(IC2))
           ZC2=Z(LLIGC2(IC2))
           DIS=SQRT((XC2-XN)**2+(YC2-YN)**2+(ZC2-ZN)**2)
           IF (DIS.LT.2.AND. NAMATM(LLIGC2(IC2)).EQ.'  CN1') THEN
             XH1=XP1Np1(IN)
             YH1=YP1Np1(IN)
             ZH1=ZP1Np1(IN)
             XH2=XP2Np1(IN)
             YH2=YP2Np1(IN)
             ZH2=ZP2Np1(IN)
            IF(
     $         (ABS(XO1-XH1).LT.DIS2.AND.
     $          ABS(YO1-YH1).LT.DIS2.AND.
     $          ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $         (ABS(XO2-XH1).LT.DIS2.AND.
     $          ABS(YO2-YH1).LT.DIS2.AND.
     $          ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $         (ABS(XO1-XH2).LT.DIS2.AND.
     $          ABS(YO1-YH2).LT.DIS2.AND.
     $          ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $         (ABS(XO2-XH2).LT.DIS2.AND.
     $          ABS(YO2-YH2).LT.DIS2.AND.
     $          ABS(ZO2-ZH2).LT.DIS2    )     ) THEN
              DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
              DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
              DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
              DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
              DIS=MIN(DIS11,DIS12)
              DIS=MIN(DIS,DIS21)
              DIS=MIN(DIS,DIS22)
             IF(DIS.LT.DIS2)THEN  
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FNH*MIN(1.0,VALUE)
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGNpl(IN))
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGNpl(IN))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             ENDIF
            ENDIF
           ENDIF
          END DO
         END DO
C
C
C
C           - FIND N3 atoms -
C
         FNH=-1.60
         DIS1=2.0
         DIS2=3.5
         DO IN3=1,NLIGN3 
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N30')then
           FNH=-0.80
           DIS1=3.00
           DIS2=4.00
             XN=X(LLIGN3(IN3))
             YN=Y(LLIGN3(IN3))
             ZN=Z(LLIGN3(IN3))
             IF((ABS(XSG-XN).LT.DIS2.AND.
     $           ABS(YSG-YN).LT.DIS2.AND.
     $           ABS(ZSG-ZN).LT.DIS2    )     ) THEN
               DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
               IF(DIS.LT.DIS2)THEN
                 NLGCYS(ICYS)=NLGCYS(ICYS)+1
                 NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGN3(IN3))
                 LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(3,ICYS,NLGCYS(ICYS))=FNH*MIN(1.0,VALUE)
                 PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
C
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))='CYS'
                 NUMSDC(7,IN3,NSDN3(IN3))=LCYSRS(ICYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
             END IF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N31')then
           FNH=-0.80
           DIS1=3.00
           DIS2=4.00
             XN=X(LLIGN3(IN3))
             YN=Y(LLIGN3(IN3))
             ZN=Z(LLIGN3(IN3))
             IF((ABS(XSG-XN).LT.DIS2.AND.
     $           ABS(YSG-YN).LT.DIS2.AND.
     $           ABS(ZSG-ZN).LT.DIS2    )     ) THEN
               DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
               IF(DIS.LT.DIS2)THEN
                 NLGCYS(ICYS)=NLGCYS(ICYS)+1
                 NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGN3(IN3))
                 LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(3,ICYS,NLGCYS(ICYS))=FNH*MIN(1.0,VALUE)
                 PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
C
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))='CYS'
                 NUMSDC(7,IN3,NSDN3(IN3))=LCYSRS(ICYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
             END IF
          ENDIF
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N32') THEN
C
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XH1=XN3P1(IN3)
           YH1=YN3P1(IN3)
           ZH1=ZN3P1(IN3)
           XH2=XN3P2(IN3)
           YH2=YN3P2(IN3)
           ZH2=ZN3P2(IN3)
           IF((ABS(XSG-XH1).LT.DIS2.AND.
     $         ABS(YSG-YH1).LT.DIS2.AND.
     $         ABS(ZSG-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XSG-XH2).LT.DIS2.AND.
     $         ABS(YSG-YH2).LT.DIS2.AND.
     $         ABS(ZSG-ZH2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XSG-XH1)**2+(YSG-YH1)**2+(ZSG-ZH1)**2)
             DIS12=SQRT((XSG-XH2)**2+(YSG-YH2)**2+(ZSG-ZH2)**2)
             IF (DIS11.LT.DIS12) THEN
               XP=XH1
               YP=YH1
               ZP=ZH1
               DIS21=SQRT((XSG-XH1)**2+(YSG-YH1)**2+(ZSG-ZH1)**2)
             ELSE
               XP=XH2
               YP=YH2
               ZP=ZH2
             ENDIF
             DIS=SQRT((XSG-XP)**2+(YSG-YP)**2+(ZSG-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPS=-(XP-XSG)/DIS
             YVPS=-(YP-YSG)/DIS
             ZVPS=-(ZP-ZSG)/DIS
             AGPS=XVNP*XVPS + YVNP*YVPS + ZVNP*ZVPS
             IF(DIS.LT.DIS2 .AND. AGPS.GT.0.001)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FNH*VALUE*AGPS
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGN3(IN3))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGN3(IN3))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
C
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='CYS'       
               NUMSDC(7,IN3,NSDN3(IN3))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)*AGPS
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             ENDIF
           ENDIF
          ENDIF
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N33') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XP=XN3P(IN3)
           YP=YN3P(IN3)
           ZP=ZN3P(IN3)
           IF((ABS(XO1-XP).LT.DIS2.AND.
     $         ABS(YO1-YP).LT.DIS2.AND.
     $         ABS(ZO1-ZP).LT.DIS2    ) .OR.
     $        (ABS(XO2-XP).LT.DIS2.AND.
     $         ABS(YO2-YP).LT.DIS2.AND.
     $         ABS(ZO2-ZP).LT.DIS2    )     ) THEN
             DISO1P=SQRT((XO1-XP)**2+(YO1-YP)**2+(ZO1-ZP)**2)
             DISO2P=SQRT((XO2-XP)**2+(YO2-YP)**2+(ZO2-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISO1P.LT.DISO2P) THEN
               DIS=DISO1P
               XVPO1=-(XP-XO1)/DISO1P
               YVPO1=-(YP-YO1)/DISO1P
               ZVPO1=-(ZP-ZO1)/DISO1P
               AGPO=XVNP*XVPO1 + YVNP*YVPO1 + ZVNP*ZVPO1
             ENDIF
             IF(DISO2P.LT.DISO1P) THEN
               DIS=DISO2P
               XVPO2=-(XP-XO2)/DISO2P
               YVPO2=-(YP-YO2)/DISO2P
               ZVPO2=-(ZP-ZO2)/DISO2P
               AGPO=XVNP*XVPO2 + YVNP*YVPO2 + ZVNP*ZVPO2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
                NLGCYS(ICYS)=NLGCYS(ICYS)+1
                VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                VALUE=MIN(1.0,VALUE)
                VALLIG(1,ICYS,NLGCYS(ICYS))=FNH*VALUE*AGPO
                LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGN3(IN3))
                NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGN3(IN3))
                PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
C
                NSDN3(IN3)=NSDN3(IN3)+1
                NAMSDC(7,IN3,NSDN3(IN3))='CYS'
                NUMSDC(7,IN3,NSDN3(IN3))=LCYSRS(ICYS)
                VALSDC(7,IN3,NSDN3(IN3))=-FNH*VALUE*AGPO
                PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
              ENDIF
            ENDIF
           ENDIF
C
 
         END DO
C
C
C
C
C           - FIND N1 atoms -
C
C           The ghost atom is used to take into
C           account the orientation of the lone pair
C
         FSG=-1.20
         DIS1=2.0
         DIS2=3.5
         DO IN1=1,NLIGN1
           XG=XN1G(IN1)
           YG=YN1G(IN1)
           ZG=ZN1G(IN1)
           XN=X(LLIGN1(IN1))
           YN=Y(LLIGN1(IN1))
           ZN=Z(LLIGN1(IN1))
           IF((ABS(XSG-XG).LT.DIS2.AND.
     $         ABS(YSG-YG).LT.DIS2.AND.
     $         ABS(ZSG-ZG).LT.DIS2    )     ) THEN
             DISSGG=SQRT((XSG-XG)**2+(YSG-YG)**2+(ZSG-ZG)**2)
             VECNRM=SQRT((XG-XN)**2+(YG-YN)**2+(ZG-ZN)**2)
             XVNG=(XG-XN)/VECNRM
             YVNG=(YG-YN)/VECNRM
             ZVNG=(ZG-ZN)/VECNRM
             XVGSG=-(XG-XSG)/DISSGG
             YVGSG=-(YG-YSG)/DISSGG
             ZVGSG=-(ZG-ZSG)/DISSGG
             AGGSG=XVNG*XVGSG + YVNG*YVGSG + ZVNG*ZVGSG
             IF(DIS.LT.DIS2 .AND. AGGSG.GT.0.001)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=FSG*MIN(1.0,VALUE)*AGGSG
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGN1(IN1))
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGN1(IN1))
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLIG(3,ICYS,NLGCYS(ICYS))
             ENDIF
           ENDIF
         END DO
C
C
C
C
C           - FIND Oco group -
C
         FSH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ILIG=1,NLIGC2
          IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                     .AND.ZCac(ILIG).NE.0)THEN
           XO1=XLO1(ILIG)
           YO1=YLO1(ILIG)
           ZO1=ZLO1(ILIG)
           XO2=XLO2(ILIG)
           YO2=YLO2(ILIG)
           ZO2=ZLO2(ILIG)
           IF((ABS(XO1-XSG).LT.DIS2.AND.
     $         ABS(YO1-YSG).LT.DIS2.AND.
     $         ABS(ZO1-ZSG).LT.DIS2    ) .OR.
     $        (ABS(XO2-XSG).LT.DIS2.AND.
     $         ABS(YO2-YSG).LT.DIS2.AND.
     $         ABS(ZO2-ZSG).LT.DIS2    )     ) THEN
             DISO1S=SQRT((XO1-XSG)**2+(YO1-YSG)**2+(ZO1-ZSG)**2)
             DISO2S=SQRT((XO2-XSG)**2+(YO2-YSG)**2+(ZO2-ZSG)**2)
             DIS=MIN(DISO1S,DISO2S)
             IF(DIS.LT.DIS2)THEN
               NLGCYS(ICYS)=NLGCYS(ICYS)+1
               NAMLIG(3,ICYS,NLGCYS(ICYS))=NAMRES(LLIGC2(ILIG))
               LABEL(3,ICYS,NLGCYS(ICYS))=NBATM(LLIGC2(ILIG))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(3,ICYS,NLGCYS(ICYS))=-FSH*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALSDC(3,ICYS,NLGCYS(ICYS))
C
               NSDLCAR(ILIG)=NSDLCAR(ILIG)+1
               NAMSDC(11,ILIG,NSDLCAR(ILIG))='CYS'
               NUMSDC(11,ILIG,NSDLCAR(ILIG))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,ILIG,NSDLCAR(ILIG))=FSH*MIN(1.0,VALUE)
               PK1LGCAR(ILIG)=PK1LGCAR(ILIG)
     $                         +VALSDC(11,ILIG,NSDLCAR(ILIG))
C
             END IF
           END IF
          ENDIF
         END DO
C
C
C
C           -- IONIZABLE INTERACTION --
C
C           N3+, Npl+, Oco- are in the iterative part STEP 8
C
C
C
C           - FIND Charged Atom -
C
      call charatm(3,icys,xsg,ysg,zsg,ncllcys,namlcol,namres,
     $ nblcol,nbatm,vallcol,pk1cys,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
C
       END DO
C
C
C
C      **********************
C      STEP 7. TYR
C      **********************
C
C
       DO ITYR=1, NTYR
         XOH=X(LTYROH(ITYR))
         YOH=Y(LTYROH(ITYR))
         ZOH=Z(LTYROH(ITYR))
C
C
C
C        -- 2. SIDECHAIN INTERACTION --
C
C
C           - FIND SER -
C
         FOH=-0.80
         DIS1=3.50
         DIS2=4.50
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XOH-XSER).LT.DIS2 .AND.
     $         ABS(YOH-YSER).LT.DIS2 .AND.
     $         ABS(ZOH-ZSER).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XSER)**2+(YOH-YSER)**2+(ZOH-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='SER'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND THR -
C
         FOH=-0.80
         DIS1=3.50
         DIS2=4.50
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XOH-XTHR).LT.DIS2 .AND.
     $         ABS(YOH-YTHR).LT.DIS2 .AND.
     $         ABS(ZOH-ZTHR).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XTHR)**2+(YOH-YTHR)**2+(ZOH-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='THR'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND ASN -
C
         FNH=-0.80
         DIS1=2.50
         DIS2=3.50
         DO IASN=1,NASN
           XASN1=XASNP1(IASN)
           YASN1=YASNP1(IASN)
           ZASN1=ZASNP1(IASN)
           XASN2=XASNP2(IASN)
           YASN2=YASNP2(IASN)
           ZASN2=ZASNP2(IASN)
           IF((ABS(XOH-XASN1).LT.DIS2.AND.
     $         ABS(YOH-YASN1).LT.DIS2.AND.
     $         ABS(ZOH-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XOH-XASN2).LT.DIS2.AND.
     $         ABS(YOH-YASN2).LT.DIS2.AND.
     $         ABS(ZOH-ZASN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XOH-XASN1)**2+(YOH-YASN1)**2+(ZOH-ZASN1)**2)
             DIS12=SQRT((XOH-XASN2)**2+(YOH-YASN2)**2+(ZOH-ZASN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='ASN'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FNH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FNH=-0.80
         DIS1=2.50
         DIS2=3.50
         DO IGLN=1,NGLN
           XGLN1=XGLNP1(IGLN)
           YGLN1=YGLNP1(IGLN)
           ZGLN1=ZGLNP1(IGLN)
           XGLN2=XGLNP2(IGLN)
           YGLN2=YGLNP2(IGLN)
           ZGLN2=ZGLNP2(IGLN)
           IF((ABS(XOH-XGLN1).LT.DIS2.AND.
     $         ABS(YOH-YGLN1).LT.DIS2.AND.
     $         ABS(ZOH-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XOH-XGLN2).LT.DIS2.AND.
     $         ABS(YOH-YGLN2).LT.DIS2.AND.
     $         ABS(ZOH-ZGLN2).LT.DIS2    )     ) THEN
             DIS11=SQRT((XOH-XGLN1)**2+(YOH-YGLN1)**2+(ZOH-ZGLN1)**2)
             DIS12=SQRT((XOH-XGLN2)**2+(YOH-YGLN2)**2+(ZOH-ZGLN2)**2)
             DIS=MIN(DIS11,DIS12)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='GLN'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FNH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C
C           - FIND TRP -
C
         FNH=-0.80
         DIS1=2.50
         DIS2=3.50
         DO ITRP=1,NTRP
           XTRP1=XTRPP1(ITRP)
           YTRP1=YTRPP1(ITRP)
           ZTRP1=ZTRPP1(ITRP)
           IF((ABS(XOH-XTRP1).LT.DIS2.AND.
     $         ABS(YOH-YTRP1).LT.DIS2.AND.
     $         ABS(ZOH-ZTRP1).LT.DIS2    )     ) THEN
             DIS=SQRT((XOH-XTRP1)**2+(YOH-YTRP1)**2+(ZOH-ZTRP1)**2)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='TRP'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LTRPRS(ITRP)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FNH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C
C        -- 3. FIND BACKBONE INTERACTION --
C
         FBKB=-1.20
         DIS1=3.5
         DIS2=4.5
         DO I=2,NPRTON
           XP=XPRTON(I)
           YP=YPRTON(I)
           ZP=ZPRTON(I)
           XN=XNITRN(I)
           YN=YNITRN(I)
           ZN=ZNITRN(I)
           IF((ABS(XOH-XP).LT.DIS2.AND.
     $         ABS(YOH-YP).LT.DIS2.AND.
     $         ABS(ZOH-ZP).LT.DIS2    )     ) THEN
             DIS=SQRT((XOH-XP)**2+(YOH-YP)**2+(ZOH-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPS=-(XP-XOH)/DIS
             YVPS=-(YP-YOH)/DIS
             ZVPS=-(ZP-ZOH)/DIS
             AGPS=XVNP*XVPS + YVNP*YVPS + ZVNP*ZVPS
             IF(DIS.LT.DIS2 .AND. AGPS.GT.0.001)THEN
               NBKTYR(ITYR)=NBKTYR(ITYR)+1
               NAMBKB(4,ITYR,NBKTYR(ITYR))=NAMPRT(I)
               NUMBKB(4,ITYR,NBKTYR(ITYR))=NUMPRT(I)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALBKB(4,ITYR,NBKTYR(ITYR))=FBKB*VALUE*AGPS
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALBKB(4,ITYR,NBKTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C        -- COULOMBIC INTERACTIONS --
C
C
C           - FIND ARG H-BONDING
C
         FOH=-0.80
         DIS1=2.50
         DIS2=4.00
         DO IARG=1,NARG
           XH1=XARGP1(IARG)
           YH1=YARGP1(IARG)
           ZH1=ZARGP1(IARG)
           XH2=XARGP2(IARG)
           YH2=YARGP2(IARG)
           ZH2=ZARGP2(IARG)
           XH3=XARGP3(IARG)
           YH3=YARGP3(IARG)
           ZH3=ZARGP3(IARG)
           XH4=XARGP4(IARG)
           YH4=YARGP4(IARG)
           ZH4=ZARGP4(IARG)
           XH5=XARGP5(IARG)
           YH5=YARGP5(IARG)
           ZH5=ZARGP5(IARG)
           IF(
     $        (ABS(XOH-XH1).LT.DIS2.AND.
     $         ABS(YOH-YH1).LT.DIS2.AND.
     $         ABS(ZOH-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH2).LT.DIS2.AND.
     $         ABS(YOH-YH2).LT.DIS2.AND.
     $         ABS(ZOH-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH3).LT.DIS2.AND.
     $         ABS(YOH-YH3).LT.DIS2.AND.
     $         ABS(ZOH-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH4).LT.DIS2.AND.
     $         ABS(YOH-YH4).LT.DIS2.AND.
     $         ABS(ZOH-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XOH-XH5).LT.DIS2.AND.
     $         ABS(YOH-YH5).LT.DIS2.AND.
     $         ABS(ZOH-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XOH-XH1)**2+(YOH-YH1)**2+(ZOH-ZH1)**2)
             DIS12=SQRT((XOH-XH2)**2+(YOH-YH2)**2+(ZOH-ZH2)**2)
             DIS13=SQRT((XOH-XH3)**2+(YOH-YH3)**2+(ZOH-ZH3)**2)
             DIS14=SQRT((XOH-XH4)**2+(YOH-YH4)**2+(ZOH-ZH4)**2)
             DIS15=SQRT((XOH-XH5)**2+(YOH-YH5)**2+(ZOH-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             IF(DIS.LT.DIS2)THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='ARG'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C Coulombic interaction between TYR residue and ARG(+)
C is in the non-iterative part because Arg is assumed to be
C always positively charged when Tyr group titrates
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
           IF(NMASS(4,ITYR)+NMASS(6,IARG).GT.900 .OR.
     $     (NMASS(4,ITYR).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XOH-XCZ).LT.DIS2.AND.
     $        ABS(YOH-YCZ).LT.DIS2.AND.
     $        ABS(ZOH-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XOH-XCZ)**2+(YOH-YCZ)**2+(ZOH-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))='ARG'
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LARGRS(IARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
C
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))='TYR'
               NUMCOL(6,IARG,NCLARG(IARG))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=-FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END DO
C
C
C
C
C
C        -- 5. LIGAND INTERACTION --
C
C           -- H-BONDING --
C
C
C           - FIND O3 group -
C
C identify the Oh group according to the number
C of the linked atom
C
       DO IO3=1,NLIGO3
         XO3=X(LLIGO3(IO3))
         YO3=Y(LLIGO3(IO3))
         ZO3=Z(LLIGO3(IO3))
        NB_neighbours=0
         DO I=1,NLIGAND
          IF(NBATM(LLIGAND(I)).NE.NBATM(LLIGO3(IO3))) THEN
           Xb=X(LLIGAND(I))
           Yb=Y(LLIGAND(I))
           Zb=Z(LLIGAND(I))
           DISn=SQRT((Xb-XO3)**2+(Yb-YO3)**2+(Zb-ZO3)**2)
           IF (DISn.LT.2.0) THEN
              NB_neighbours=NB_neighbours+1
           ENDIF
          ENDIF
         END DO
C
C           - FIND OH group -
C
         IF (NB_neighbours.EQ.1) THEN
 
           FOH=-0.80
           DIS1=3.50
           DIS2=4.50
           IF((ABS(XOH-XO3).LT.DIS2 .AND.
     $         ABS(YOH-YO3).LT.DIS2 .AND.
     $         ABS(ZOH-ZO3).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XO3)**2+(YOH-YO3)**2+(ZOH-ZO3)**2)
             IF(DIS.LT.DIS2) THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FOH*MIN(1.0,VALUE)
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGO3(IO3))
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGO3(IO3))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             END IF
           END IF
         ENDIF
C
C
C           - FIND O3 (non Oh) group -
C
         IF (NB_neighbours.EQ.2) THEN
C
           FOO=1.60
           DIS1=3.50
           DIS2=4.50
           IF((ABS(XOH-XO3).LT.DIS2 .AND.
     $         ABS(YOH-YO3).LT.DIS2 .AND.
     $         ABS(ZOH-ZO3).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XO3)**2+(YOH-YO3)**2+(ZOH-ZO3)**2)
             IF(DIS.LT.DIS2) THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FOO*MIN(1.0,VALUE)
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGO3(IO3))
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGO3(IO3))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             END IF
           END IF
         ENDIF
       END DO
C
C
C
C
C           - FIND Cl group -
C
         FOCl=1.60
         DIS1=3.50
         DIS2=4.50
C
         DO ICl=1,NLIGCl
           XCl=X(LLIGCl(ICl))
           YCl=Y(LLIGCl(ICl))
           ZCl=Z(LLIGCl(ICl))
           IF((ABS(XOH-XCl).LT.DIS2 .AND.
     $         ABS(YOH-YCl).LT.DIS2 .AND.
     $         ABS(ZOH-ZCl).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XCl)**2+(YOH-YCl)**2+(ZOH-ZCl)**2)
             IF(DIS.LT.DIS2) THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FOCl*MIN(1.0,VALUE)
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGCl(ICl))
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGCl(ICl))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C
C           - FIND F group -
C
         FOF=1.60
         DIS1=3.00
         DIS2=4.00
C
         DO IF=1,NLIGF
           XF=X(LLIGF(IF))
           YF=Y(LLIGF(IF))
           ZF=Z(LLIGF(IF))
           IF((ABS(XOH-XF).LT.DIS2 .AND.
     $         ABS(YOH-YF).LT.DIS2 .AND.
     $         ABS(ZOH-ZF).LT.DIS2     )     ) THEN
             DIS=SQRT((XOH-XF)**2+(YOH-YF)**2+(ZOH-ZF)**2)
             IF(DIS.LT.DIS2) THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FOF*MIN(1.0,VALUE)
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGF(IF))
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGF(IF))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C           - FIND Amide group (Nam) -
C
         FNH=-1.20
         DIS1=3.5
         DIS2=4.5
         DO I=1,NLIGNam
          IF(NAMATM(LLIGNam(I)).NE.'  C3 ') THEN
           XN=X(LLIGNam(I))
           YN=Y(LLIGNam(I))
           ZN=Z(LLIGNam(I))
           XP=XHLIG(I)
           YP=YHLIG(I)
           ZP=ZHLIG(I)
           IF((ABS(XOH-XP).LT.DIS2.AND.
     $         ABS(YOH-YP).LT.DIS2.AND.
     $         ABS(ZOH-ZP).LT.DIS2    )     ) THEN
             DIS=SQRT((XOH-XP)**2+(YOH-YP)**2+(ZOH-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             XVPO=-(XP-XOH)/DIS
             YVPO=-(YP-YOH)/DIS
             ZVPO=-(ZP-ZOH)/DIS
             AGPO=XVNP*XVPO + YVNP*YVPO + ZVNP*ZVPO
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FNH*VALUE*AGPO
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGNam(I))
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGNam(I))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             END IF
           END IF
         ENDIF
         END DO
C
C
C           - FIND Aromatic Nitrogen (Nar) -
C
c dmr start
         FNh=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ILIG=1,NLIGNar
           Xp=XgLIG(ILIG)
           Yp=YgLIG(ILIG)
           Zp=ZgLIG(ILIG)
            if(ABS(XOh-Xp).LT.DIS2.AND.
     $         ABS(YOh-Yp).LT.DIS2.AND.
     $         ABS(ZOh-Zp).LT.DIS2          ) THEN
             DIS=SQRT((XOh-Xp)**2+(YOh-Yp)**2+(ZOh-Zp)**2)
             IF(DIS.LT.DIS2)THEN
               NLGtyr(Ityr)=NLGtyr(Ityr)+1
               VALLIG(4,Ityr,NLGtyr(Ityr))=FNh*MIN(1.0,VALUE)
               LABEL(4,Ityr,NLGtyr(Ityr))=NBATM(LLIGNar(ILIG))
               NAMLIG(4,Ityr,NLGtyr(Ityr))=NAMRES(LLIGNar(ILIG))
               PK1tyr(Ityr)=PK1tyr(Ityr)+VALLIG(4,Ityr,NLGtyr(Ityr))
 
               NSDNar(ILIG)=NSDNar(ILIG)+1
               NAMSDC(8,ILIG,NSDNar(ILIG))='TYR'
               NUMSDC(8,ILIG,NSDNar(ILIG))=LtyrRS(Ityr)
               VALSDC(8,ILIG,NSDNar(ILIG))=-FNh*MIN(1.0,VALUE)
               PK1LGNar(ILIG)=PK1LGNar(ILIG)
     $                       +VALSDC(8,ILIG,NSDNar(ILIG))
 
c              STRONG H-BONDING IF BURIED
c              (MOSTLY COULOMBIC INTERACTION)
C
C              IF(NMASS(4,Ityr)+NMASS(8,ILIG).GT.900 .OR.
C    $     (NMASS(4,Ityr).GT.400.AND.NMASS(8,ILIG).GT.400))THEN
C                VALLIG(4,Ityr,NLGtyr(Ityr))=-1.60
C                VALSDC(8,ILIG,NSDNar(ILIG))=+1.60
C              END IF
c              PK1tyr(Ityr)=PK1tyr(Ityr)+VALLIG(4,Ityr,NLGtyr(Ityr))
c              PK1LGNar(ILIG)=PK1LGNar(ILIG)
c    $                       +VALSDC(8,ILIG,NSDNar(ILIG))
             END IF
           END IF
         END DO
c dmr end
C
C
C           - FIND C=O - 
C dcb
C angle dependence removed july 2007
C dcb
C
         FOO=0.8
         DIS1=3.5
         DIS2=4.5
         DO I=1,NLIGO2 
           XO=X(LLIGO2(I))
           YO=Y(LLIGO2(I))
           ZO=Z(LLIGO2(I))
           IF((ABS(XOH-XO).LT.DIS2.AND.
     $         ABS(YOH-YO).LT.DIS2.AND.
     $         ABS(ZOH-ZO).LT.DIS2    )     ) THEN
             DIS=SQRT((XOH-XO)**2+(YOH-YO)**2+(ZOH-ZO)**2)
             IF(DIS.LT.DIS2)THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FOO*VALUE
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGO2(I))
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGO2(I))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             END IF
           END IF
         END DO
C
C
C 
C
C           - FIND Npl -        
C
         FNH=-0.80
         DIS1=2.50
         DIS2=4.00
         DO IC2=1,NLIGC2
C
C  -- GUADININIUM TYPE GROUP (eg. 1K1O.pdb) --
C
C dcb july 2007
C There is no angle dependence when the H-bonding interaction
C between Cys and Arg group is considered so for the moment
C we don't include the angle dependence for the interaction
C between Cys residues and Np1 pair group from the ligand
C dcb
C
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg1P(IC2)
           YH1=YNg1P(IC2)
           ZH1=ZNg1P(IC2)
           XH2=XNg2P2(IC2)  
           YH2=YNg2P2(IC2)
           ZH2=ZNg2P2(IC2)
           XH3=XNg2P1(IC2)
           YH3=YNg2P1(IC2)
           ZH3=ZNg2P1(IC2)
           XH4=XNg3P2(IC2)
           YH4=YNg3P2(IC2)
           ZH4=ZNg3P2(IC2)
           XH5=XNg3P1(IC2)
           YH5=YNg3P1(IC2)
           ZH5=ZNg3P1(IC2)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH5).LT.DIS2.AND.
     $         ABS(YO1-YH5).LT.DIS2.AND.
     $         ABS(ZO1-ZH5).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH5).LT.DIS2.AND.
     $         ABS(YO2-YH5).LT.DIS2.AND.
     $         ABS(ZO2-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS15=SQRT((XO1-XH5)**2+(YO1-YH5)**2+(ZO1-ZH5)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS25=SQRT((XO2-XH5)**2+(YO2-YH5)**2+(ZO2-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             DIS=MIN(DIS,DIS25)
             IF(DIS.LT.DIS2)THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FNH*MIN(1.0,VALUE)
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGNpl(IN))
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGNpl(IN))
                 NSDcg(Ic2)=NSDcg(Ic2)+1
                 NAMSDC(9,Ic2,NSDcg(Ic2))='TYR'
                 NUMSDC(9,Ic2,NSDcg(Ic2))=LTYRRS(ITYR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(9,Ic2,NSDcg(Ic2))=-FNH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALSDC(9,IC2,NSDCg(IC2))
             END IF
           END IF
          ENDIF
C
C  -- Np1 pair GROUP (eg. 1K1I.pdb) --
C
C dcb july 2007
C There is no angle dependence when the H-bonding interaction
C between Cys and Arg group is considered so for the moment
C we don't include the angle dependence for the interaction
C between Cys residues and Np1 pair group from the ligand
C dcb
C
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg2P2(IC2)  
           YH1=YNg2P2(IC2)
           ZH1=ZNg2P2(IC2)
           XH2=XNg2P1(IC2)
           YH2=YNg2P1(IC2)
           ZH2=ZNg2P1(IC2)
           XH3=XNg3P2(IC2)
           YH3=YNg3P2(IC2)
           ZH3=ZNg3P2(IC2)
           XH4=XNg3P1(IC2)
           YH4=YNg3P1(IC2)
           ZH4=ZNg3P1(IC2)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             IF(DIS.LT.DIS2)THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FNH*MIN(1.0,VALUE)
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGNpl(IN))
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGNpl(IN))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))

                 NSDcn2(Ic2)=NSDcn2(Ic2)+1
                 NAMSDC(9,Ic2,NSDcn2(Ic2))='TYR'
                 NUMSDC(9,Ic2,NSDcn2(Ic2))=LTYRRS(ITYR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(9,Ic2,NSDcn2(Ic2))=-FNH*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NSDCN2(IC2))
C
C
C     -- THERE ARE POSSIBLY 2 H-BONDS BETWEEN C2N GROUP AND CARBOXYL
C              IF((DIS12.LT.2.2 .AND. DIS23.LT.2.2).OR.
C    $            (DIS22.LT.2.2 .AND. DIS13.LT.2.2)    )THEN
C                VALLIG(4,ITYR,NLGTYR(ITYR))=-2.40
C                VALSDC(9,IC2,NSDCn2(IC2))=2.40
C              END IF
c              PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
c              PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NSDCN2(IC2))
             END IF
           END IF
          ENDIF
         END DO
C
C
C  -- NH2 GROUP (eg. 1HPV.pdb) --
C
C     This group is assuemd to be non ionizable
C
         DO IN=1,NLIGNpl
           XN=X(LLIGNpl(IN))
           YN=Y(LLIGNpl(IN))
           ZN=Z(LLIGNpl(IN))
          DO IC2=1,NLIGC2
           XC2=X(LLIGC2(IC2))
           YC2=Y(LLIGC2(IC2))
           ZC2=Z(LLIGC2(IC2))
           DIS=SQRT((XC2-XN)**2+(YC2-YN)**2+(ZC2-ZN)**2)
           IF (DIS.LT.2.AND. NAMATM(LLIGC2(IC2)).EQ.'  CN1') THEN
             XH1=XP1Np1(IN)
             YH1=YP1Np1(IN)
             ZH1=ZP1Np1(IN)
             XH2=XP2Np1(IN)
             YH2=YP2Np1(IN)
             ZH2=ZP2Np1(IN)
            IF(
     $         (ABS(XO1-XH1).LT.DIS2.AND.
     $          ABS(YO1-YH1).LT.DIS2.AND.
     $          ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $         (ABS(XO2-XH1).LT.DIS2.AND.
     $          ABS(YO2-YH1).LT.DIS2.AND.
     $          ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $         (ABS(XO1-XH2).LT.DIS2.AND.
     $          ABS(YO1-YH2).LT.DIS2.AND.
     $          ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $         (ABS(XO2-XH2).LT.DIS2.AND.
     $          ABS(YO2-YH2).LT.DIS2.AND.
     $          ABS(ZO2-ZH2).LT.DIS2    )     ) THEN
              DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
              DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
              DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
              DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
              DIS=MIN(DIS11,DIS12)
              DIS=MIN(DIS,DIS21)
              DIS=MIN(DIS,DIS22)
             IF(DIS.LT.DIS2)THEN  
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FNH*MIN(1.0,VALUE)
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGNpl(IN))
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGNpl(IN))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             ENDIF
            ENDIF
           ENDIF
          END DO
         END DO
C
C
C
C           - FIND N3 atoms -
C
C           in the iterative part as for the Tyr-Lys interaction
C
c        END DO
C
C
C
C
C           - FIND N1 atoms -
C
C           The ghost atom is used to take into
C           account the orientation of the lone pair
C
         FOG=-1.20
         DIS1=2.0
         DIS2=3.5
         DO IN1=1,NLIGN1
           XG=XN1G(IN1)
           YG=YN1G(IN1)
           ZG=ZN1G(IN1)
           XN=X(LLIGN1(IN1))
           YN=Y(LLIGN1(IN1))
           ZN=Z(LLIGN1(IN1))
           IF((ABS(XOH-XG).LT.DIS2.AND.
     $         ABS(YOH-YG).LT.DIS2.AND.
     $         ABS(ZOH-ZG).LT.DIS2    )     ) THEN
             DISOHG=SQRT((XOH-XG)**2+(YOH-YG)**2+(ZOH-ZG)**2)
             VECNRM=SQRT((XG-XN)**2+(YG-YN)**2+(ZG-ZN)**2)
             XVNG=(XG-XN)/VECNRM
             YVNG=(YG-YN)/VECNRM
             ZVNG=(ZG-ZN)/VECNRM
             XVGOH=-(XG-XOH)/DISOHG
             YVGOH=-(YG-YOH)/DISOHG
             ZVGOH=-(ZG-ZOH)/DISOHG
             AGGOH=XVNG*XVGOH + YVNG*YVGOH + ZVNG*ZVGOH
             IF(DIS.LT.DIS2 .AND. AGGOH.GT.0.001)THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=FOG*MIN(1.0,VALUE)*AGGOH
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGN1(IN1))
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGN1(IN1))
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
             ENDIF
           ENDIF
         END DO
C
C
C           - FIND Oco group -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ILIG=1,NLIGC2
          IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                     .AND.ZCac(ILIG).NE.0)THEN
           XO1=XLO1(ILIG)
           YO1=YLO1(ILIG)
           ZO1=ZLO1(ILIG)
           XO2=XLO2(ILIG)
           YO2=YLO2(ILIG)
           ZO2=ZLO2(ILIG)
           IF((ABS(XO1-XOH).LT.DIS2.AND.
     $         ABS(YO1-YOH).LT.DIS2.AND.
     $         ABS(ZO1-ZOH).LT.DIS2    ) .OR.
     $        (ABS(XO2-XOH).LT.DIS2.AND.
     $         ABS(YO2-YOH).LT.DIS2.AND.
     $         ABS(ZO2-ZOH).LT.DIS2    )     ) THEN
             DISO1O=SQRT((XO1-XOH)**2+(YO1-YOH)**2+(ZO1-ZOH)**2)
             DISO2O=SQRT((XO2-XOH)**2+(YO2-YOH)**2+(ZO2-ZOH)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NLGTYR(ITYR)=NLGTYR(ITYR)+1
               NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGC2(ILIG))
               LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGC2(ILIG))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(4,ITYR,NLGTYR(ITYR))=-FOH*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NLGTYR(ITYR))
C
               NSDLCAR(ILIG)=NSDLCAR(ILIG)+1
               NAMSDC(11,ILIG,NSDLCAR(ILIG))='TYR'
               NUMSDC(11,ILIG,NSDLCAR(ILIG))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,ILIG,NSDLCAR(ILIG))=FOH*MIN(1.0,VALUE)
               PK1LGCAR(ILIG)=PK1LGCAR(ILIG)
     $                         +VALSDC(11,ILIG,NSDLCAR(ILIG))
C
             END IF
           END IF
          ENDIF
         END DO
C
C
C
C           -- IONIZABLE INTERACTION --
C           N3+ Npl+, Oco- in the iterative part STEP 8
C
C
C
C           - FIND Charged Atom -
C
      call charatm(4,ityr,xoh,yoh,zoh,nclltyr,namlcol,namres,
     $ nblcol,nbatm,vallcol,pk1tyr,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
       END DO
C
C
C
C
C
c dmr H-bonds: ligand w/ backbone; ligand w/ non-ionizable grps.
c      **********************
c      STEP 7A. -CN2 and guadinium grp (Cg).
c      **********************
c
c 1HPV -NH2 , Np1
c         
c  -CN2
C  -- Np1 pair GROUP (eg. 1K1I.pdb) --
C
         DO ic2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg2P2(IC2)  
           YH1=YNg2P2(IC2)
           ZH1=ZNg2P2(IC2)
           XH2=XNg2P1(IC2)
           YH2=YNg2P1(IC2)
           ZH2=ZNg2P1(IC2)
           XH3=XNg3P2(IC2)
           YH3=YNg3P2(IC2)
           ZH3=ZNg3P2(IC2)
           XH4=XNg3P1(IC2)
           YH4=YNg3P1(IC2)
           ZH4=ZNg3P1(IC2)
c Lys params
         FBKB=0.80
         DIS1=2.00
         DIS2=4.00
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(XH1-XBKO).LT.DIS2.AND.
     $         ABS(YH1-YBKO).LT.DIS2.AND.
     $         ABS(ZH1-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH2-XBKO).LT.DIS2.AND.
     $         ABS(YH2-YBKO).LT.DIS2.AND.
     $         ABS(ZH2-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH3-XBKO).LT.DIS2.AND.
     $         ABS(YH3-YBKO).LT.DIS2.AND.
     $         ABS(ZH3-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH4-XBKO).LT.DIS2.AND.
     $         ABS(YH4-YBKO).LT.DIS2.AND.
     $         ABS(ZH4-ZBKO).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XBKO)**2+(YH1-YBKO)**2+(ZH1-ZBKO)**2)
             DISH2=SQRT((XH2-XBKO)**2+(YH2-YBKO)**2+(ZH2-ZBKO)**2)
             DISH3=SQRT((XH3-XBKO)**2+(YH3-YBKO)**2+(ZH3-ZBKO)**2)
             DISH4=SQRT((XH4-XBKO)**2+(YH4-YBKO)**2+(ZH4-ZBKO)**2)
             VECNM=SQRT((XBKC-XBKO)**2+(YBKC-YBKO)**2+(ZBKC-ZBKO)**2)
             XVCO=(XBKO-XBKC)/VECNM
             YVCO=(YBKO-YBKC)/VECNM
             ZVCO=(ZBKO-ZBKC)/VECNM
             dis=min(dish1,dish2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             if(dis.eq.dish1)then
               XVOH1=(XH1-XBKO)/DISH1
               YVOH1=(YH1-YBKO)/DISH1
               ZVOH1=(ZH1-ZBKO)/DISH1
               AGOH=XVCO*XVOH1+YVCO*YVOH1+ZVCO*ZVOH1
             ENDIF
             if(dis.eq.dish2)then
               XVOH2=(XH2-XBKO)/DISH2
               YVOH2=(YH2-YBKO)/DISH2
               ZVOH2=(ZH2-ZBKO)/DISH2
               AGOH=XVCO*XVOH2+YVCO*YVOH2+ZVCO*ZVOH2
             ENDIF
             if(dis.eq.dish3)then
               XVOH3=(XH3-XBKO)/DISH3
               YVOH3=(YH3-YBKO)/DISH3
               ZVOH3=(ZH3-ZBKO)/DISH3
               AGOH=XVCO*XVOH3+YVCO*YVOH3+ZVCO*ZVOH3
             ENDIF
             if(dis.eq.dish4)then
               XVOH4=(XH4-XBKO)/DISH4
               YVOH4=(YH4-YBKO)/DISH4
               ZVOH4=(ZH4-ZBKO)/DISH4
               AGOH=XVCO*XVOH4+YVCO*YVOH4+ZVCO*ZVOH4
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NBKCN2(IC2)=NBKCN2(IC2)+1
               NAMBKB(9,IC2,NBKCN2(IC2))=NAMPRT(I-1)
               NUMBKB(9,IC2,NBKCN2(IC2))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(9,IC2,NBKCN2(IC2))=FBKB*AGOH*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALBKB(9,IC2,NBKCN2(IC2))
             END IF
           END IF
         END DO
c
         end if
         END DO

c non-ionizable, from His step
         DO ic2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg2P2(IC2)  
           YH1=YNg2P2(IC2)
           ZH1=ZNg2P2(IC2)
           XH2=XNg2P1(IC2)
           YH2=YNg2P1(IC2)
           ZH2=ZNg2P1(IC2)
           XH3=XNg3P2(IC2)
           YH3=YNg3P2(IC2)
           ZH3=ZNg3P2(IC2)
           XH4=XNg3P1(IC2)
           YH4=YNg3P1(IC2)
           ZH4=ZNg3P1(IC2)
C
C           - FIND ASN -
C
C dcb july 2007
C At first we don't include the angle dependence 
C dcb end
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF((ABS(XH1-XASN).LT.DIS2.AND.
     $         ABS(YH1-YASN).LT.DIS2.AND.
     $         ABS(ZH1-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XASN).LT.DIS2.AND.
     $         ABS(YH2-YASN).LT.DIS2.AND.
     $         ABS(ZH2-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH3-XASN).LT.DIS2.AND.
     $         ABS(YH3-YASN).LT.DIS2.AND.
     $         ABS(ZH3-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH4-XASN).LT.DIS2.AND.
     $         ABS(YH4-YASN).LT.DIS2.AND.
     $         ABS(ZH4-ZASN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XASN)**2+(YH1-YASN)**2+(ZH1-ZASN)**2)
             DISH2=SQRT((XH2-XASN)**2+(YH2-YASN)**2+(ZH2-ZASN)**2)
             DISH3=SQRT((XH3-XASN)**2+(YH3-YASN)**2+(ZH3-ZASN)**2)
             DISH4=SQRT((XH4-XASN)**2+(YH4-YASN)**2+(ZH4-ZASN)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             IF(DIS.LT.DIS2)THEN
               NSDCN2(IC2)=NSDCN2(IC2)+1
               NAMSDC(9,IC2,NSDCN2(IC2))='ASN'
               NUMSDC(9,IC2,NSDCN2(IC2))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCN2(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NSDCN2(IC2))
             END IF
           END IF
         END DO
c
c           - FIND GLN -
C
C dcb july 2007
C At first we don't include the angle dependence
C dcb end
C
c
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF((ABS(XH1-XGLN).LT.DIS2.AND.
     $         ABS(YH1-YGLN).LT.DIS2.AND.
     $         ABS(ZH1-ZGLN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XGLN).LT.DIS2.AND.
     $         ABS(YH2-YGLN).LT.DIS2.AND.
     $         ABS(ZH2-ZGLN).LT.DIS2    ) .OR.           
     $        (ABS(XH3-XGLN).LT.DIS2.AND.
     $         ABS(YH3-YGLN).LT.DIS2.AND.
     $         ABS(ZH3-ZGLN).LT.DIS2    ) .OR.           
     $        (ABS(XH4-XGLN).LT.DIS2.AND.
     $         ABS(YH4-YGLN).LT.DIS2.AND.
     $         ABS(ZH4-ZGLN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XGLN)**2+(YH1-YGLN)**2+(ZH1-ZGLN)**2)
             DISH2=SQRT((XH2-XGLN)**2+(YH2-YGLN)**2+(ZH2-ZGLN)**2)
             DISH3=SQRT((XH3-XGLN)**2+(YH3-YGLN)**2+(ZH3-ZGLN)**2)
             DISH4=SQRT((XH4-XGLN)**2+(YH4-YGLN)**2+(ZH4-ZGLN)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             IF(DIS.LT.DIS2)THEN
               NSDCN2(IC2)=NSDCN2(IC2)+1
               NAMSDC(9,IC2,NSDCN2(IC2))='GLN'
               NUMSDC(9,IC2,NSDCN2(IC2))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCN2(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NSDCN2(IC2))
             END IF
           END IF
         END DO
C
C           - FIND SER -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XH1-XSER).LT.DIS2 .AND.
     $         ABS(YH1-YSER).LT.DIS2 .AND.
     $         ABS(ZH1-ZSER).LT.DIS2     ) .OR.
     $        (ABS(XH2-XSER).LT.DIS2 .AND.
     $         ABS(YH2-YSER).LT.DIS2 .AND.
     $         ABS(ZH2-ZSER).LT.DIS2     ) .OR.           
     $        (ABS(XH3-XSER).LT.DIS2 .AND.
     $         ABS(YH3-YSER).LT.DIS2 .AND.
     $         ABS(ZH3-ZSER).LT.DIS2     ) .OR.           
     $        (ABS(XH4-XSER).LT.DIS2 .AND.
     $         ABS(YH4-YSER).LT.DIS2 .AND.
     $         ABS(ZH4-ZSER).LT.DIS2     )     ) THEN
             DISH1=SQRT((XH1-XSER)**2+(YH1-YSER)**2+(ZH1-ZSER)**2)
             DISH2=SQRT((XH2-XSER)**2+(YH2-YSER)**2+(ZH2-ZSER)**2)
             DISH3=SQRT((XH3-XSER)**2+(YH3-YSER)**2+(ZH3-ZSER)**2)
             DISH4=SQRT((XH4-XSER)**2+(YH4-YSER)**2+(ZH4-ZSER)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             IF(DIS.LT.DIS2)THEN
               NSDCN2(IC2)=NSDCN2(IC2)+1
               NAMSDC(9,IC2,NSDCN2(IC2))='SER'
               NUMSDC(9,IC2,NSDCN2(IC2))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCN2(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NSDCN2(IC2))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XH1-XTHR).LT.DIS2 .AND.
     $         ABS(YH1-YTHR).LT.DIS2 .AND.
     $         ABS(ZH1-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XH2-XTHR).LT.DIS2 .AND.
     $         ABS(YH2-YTHR).LT.DIS2 .AND.
     $         ABS(ZH2-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XH3-XTHR).LT.DIS2 .AND.
     $         ABS(YH3-YTHR).LT.DIS2 .AND.
     $         ABS(ZH3-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XH4-XTHR).LT.DIS2 .AND.
     $         ABS(YH4-YTHR).LT.DIS2 .AND.
     $         ABS(ZH4-ZTHR).LT.DIS2     )     ) THEN
             DISH1=SQRT((XH1-XTHR)**2+(YH1-YTHR)**2+(ZH1-ZTHR)**2)
             DISH2=SQRT((XH2-XTHR)**2+(YH2-YTHR)**2+(ZH2-ZTHR)**2)
             DISH3=SQRT((XH3-XTHR)**2+(YH3-YTHR)**2+(ZH3-ZTHR)**2)
             DISH4=SQRT((XH4-XTHR)**2+(YH4-YTHR)**2+(ZH4-ZTHR)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             IF(DIS.LT.DIS2)THEN
               NSDCN2(IC2)=NSDCN2(IC2)+1
               NAMSDC(9,IC2,NSDCN2(IC2))='THR'
               NUMSDC(9,IC2,NSDCN2(IC2))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCN2(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)+VALSDC(9,IC2,NSDCN2(IC2))
             END IF
           END IF
         END DO
c end non-ionizable
C
C           - FIND Charged Atom -
C
      call charatml(9,ic2,lligc2,xc,yc,zc,nclcn2,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgcn2,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO

c  Cg
C  -- GUADININIUM TYPE GROUP (eg. 1K1O.pdb) --
C
         DO IC2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg1P(IC2)
           YH1=YNg1P(IC2)
           ZH1=ZNg1P(IC2)
           XH2=XNg2P2(IC2)  
           YH2=YNg2P2(IC2)
           ZH2=ZNg2P2(IC2)
           XH3=XNg2P1(IC2)
           YH3=YNg2P1(IC2)
           ZH3=ZNg2P1(IC2)
           XH4=XNg3P2(IC2)
           YH4=YNg3P2(IC2)
           ZH4=ZNg3P2(IC2)
           XH5=XNg3P1(IC2)
           YH5=YNg3P1(IC2)
           ZH5=ZNg3P1(IC2)
c Arg params
         FBKB=0.80
         DIS1=2.00
         DIS2=4.00
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(XH1-XBKO).LT.DIS2.AND.
     $         ABS(YH1-YBKO).LT.DIS2.AND.
     $         ABS(ZH1-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH2-XBKO).LT.DIS2.AND.
     $         ABS(YH2-YBKO).LT.DIS2.AND.
     $         ABS(ZH2-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH3-XBKO).LT.DIS2.AND.
     $         ABS(YH3-YBKO).LT.DIS2.AND.
     $         ABS(ZH3-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH4-XBKO).LT.DIS2.AND.
     $         ABS(YH4-YBKO).LT.DIS2.AND.
     $         ABS(ZH4-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XH5-XBKO).LT.DIS2.AND.
     $         ABS(YH5-YBKO).LT.DIS2.AND.
     $         ABS(ZH5-ZBKO).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XBKO)**2+(YH1-YBKO)**2+(ZH1-ZBKO)**2)
             DISH2=SQRT((XH2-XBKO)**2+(YH2-YBKO)**2+(ZH2-ZBKO)**2)
             DISH3=SQRT((XH3-XBKO)**2+(YH3-YBKO)**2+(ZH3-ZBKO)**2)
             DISH4=SQRT((XH4-XBKO)**2+(YH4-YBKO)**2+(ZH4-ZBKO)**2)
             DISH5=SQRT((XH5-XBKO)**2+(YH5-YBKO)**2+(ZH5-ZBKO)**2)
             VECNM=SQRT((XBKC-XBKO)**2+(YBKC-YBKO)**2+(ZBKC-ZBKO)**2)
             XVCO=(XBKO-XBKC)/VECNM
             YVCO=(YBKO-YBKC)/VECNM
             ZVCO=(ZBKO-ZBKC)/VECNM
             dis=min(dish1,dish2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             dis=min(dis,dish5)
             if(dis.eq.dish1)then
               XVOH1=(XH1-XBKO)/DISH1
               YVOH1=(YH1-YBKO)/DISH1
               ZVOH1=(ZH1-ZBKO)/DISH1
               AGOH=XVCO*XVOH1+YVCO*YVOH1+ZVCO*ZVOH1
             ENDIF
             if(dis.eq.dish2)then
               XVOH2=(XH2-XBKO)/DISH2
               YVOH2=(YH2-YBKO)/DISH2
               ZVOH2=(ZH2-ZBKO)/DISH2
               AGOH=XVCO*XVOH2+YVCO*YVOH2+ZVCO*ZVOH2
             ENDIF
             if(dis.eq.dish3)then
               XVOH3=(XH3-XBKO)/DISH3
               YVOH3=(YH3-YBKO)/DISH3
               ZVOH3=(ZH3-ZBKO)/DISH3
               AGOH=XVCO*XVOH3+YVCO*YVOH3+ZVCO*ZVOH3
             ENDIF
             if(dis.eq.dish4)then
               XVOH4=(XH4-XBKO)/DISH4
               YVOH4=(YH4-YBKO)/DISH4
               ZVOH4=(ZH4-ZBKO)/DISH4
               AGOH=XVCO*XVOH4+YVCO*YVOH4+ZVCO*ZVOH4
             ENDIF
             if(dis.eq.dish5)then
               XVOH5=(XH5-XBKO)/DISH5
               YVOH5=(YH5-YBKO)/DISH5
               ZVOH5=(ZH5-ZBKO)/DISH5
               AGOH=XVCO*XVOH5+YVCO*YVOH5+ZVCO*ZVOH5
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NBKCg(IC2)=NBKCg(IC2)+1
               NAMBKB(9,IC2,NBKCg(IC2))=NAMPRT(I-1)
               NUMBKB(9,IC2,NBKCg(IC2))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(9,IC2,NBKCg(IC2))=FBKB*AGOH*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALBKB(9,IC2,NBKCg(IC2))
             END IF
           END IF
         END DO
c
         end if
         END DO

c non-ionizable, from CN2 (above)
         DO IC2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           XH1=XNg1P(IC2)
           YH1=YNg1P(IC2)
           ZH1=ZNg1P(IC2)
           XH2=XNg2P2(IC2)  
           YH2=YNg2P2(IC2)
           ZH2=ZNg2P2(IC2)
           XH3=XNg2P1(IC2)
           YH3=YNg2P1(IC2)
           ZH3=ZNg2P1(IC2)
           XH4=XNg3P2(IC2)
           YH4=YNg3P2(IC2)
           ZH4=ZNg3P2(IC2)
           XH5=XNg3P1(IC2)
           YH5=YNg3P1(IC2)
           ZH5=ZNg3P1(IC2)
C
C           - FIND ASN -
C
C dcb july 2007
C At first we don't include the angle dependence
C dcb end
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF((ABS(XH1-XASN).LT.DIS2.AND.
     $         ABS(YH1-YASN).LT.DIS2.AND.
     $         ABS(ZH1-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XASN).LT.DIS2.AND.
     $         ABS(YH2-YASN).LT.DIS2.AND.
     $         ABS(ZH2-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH3-XASN).LT.DIS2.AND.
     $         ABS(YH3-YASN).LT.DIS2.AND.
     $         ABS(ZH3-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH4-XASN).LT.DIS2.AND.
     $         ABS(YH4-YASN).LT.DIS2.AND.
     $         ABS(ZH4-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XH5-XASN).LT.DIS2.AND.
     $         ABS(YH5-YASN).LT.DIS2.AND.
     $         ABS(ZH5-ZASN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XASN)**2+(YH1-YASN)**2+(ZH1-ZASN)**2)
             DISH2=SQRT((XH2-XASN)**2+(YH2-YASN)**2+(ZH2-ZASN)**2)
             DISH3=SQRT((XH3-XASN)**2+(YH3-YASN)**2+(ZH3-ZASN)**2)
             DISH4=SQRT((XH4-XASN)**2+(YH4-YASN)**2+(ZH4-ZASN)**2)
             DISH5=SQRT((XH5-XASN)**2+(YH5-YASN)**2+(ZH5-ZASN)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             dis=min(dis,dish5)
             IF(DIS.LT.DIS2)THEN
               NSDCg(IC2)=NSDCg(IC2)+1
               NAMSDC(9,IC2,NSDCg(IC2))='ASN'
               NUMSDC(9,IC2,NSDCg(IC2))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCg(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALSDC(9,IC2,NSDCg(IC2))
             END IF
           END IF
         END DO
c
c           - FIND GLN -
C
C dcb july 2007
C At first we don't include the angle dependence
C dcb end
c
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF((ABS(XH1-XGLN).LT.DIS2.AND.
     $         ABS(YH1-YGLN).LT.DIS2.AND.
     $         ABS(ZH1-ZGLN).LT.DIS2    ) .OR.
     $        (ABS(XH2-XGLN).LT.DIS2.AND.
     $         ABS(YH2-YGLN).LT.DIS2.AND.
     $         ABS(ZH2-ZGLN).LT.DIS2    ) .OR.           
     $        (ABS(XH3-XGLN).LT.DIS2.AND.
     $         ABS(YH3-YGLN).LT.DIS2.AND.
     $         ABS(ZH3-ZGLN).LT.DIS2    ) .OR.           
     $        (ABS(XH4-XGLN).LT.DIS2.AND.
     $         ABS(YH4-YGLN).LT.DIS2.AND.
     $         ABS(ZH4-ZGLN).LT.DIS2    ) .OR.           
     $        (ABS(XH5-XGLN).LT.DIS2.AND.
     $         ABS(YH5-YGLN).LT.DIS2.AND.
     $         ABS(ZH5-ZGLN).LT.DIS2    )     ) THEN
             DISH1=SQRT((XH1-XGLN)**2+(YH1-YGLN)**2+(ZH1-ZGLN)**2)
             DISH2=SQRT((XH2-XGLN)**2+(YH2-YGLN)**2+(ZH2-ZGLN)**2)
             DISH3=SQRT((XH3-XGLN)**2+(YH3-YGLN)**2+(ZH3-ZGLN)**2)
             DISH4=SQRT((XH4-XGLN)**2+(YH4-YGLN)**2+(ZH4-ZGLN)**2)
             DISH5=SQRT((XH5-XGLN)**2+(YH5-YGLN)**2+(ZH5-ZGLN)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             dis=min(dis,dish5)
             IF(DIS.LT.DIS2)THEN
               NSDCg(IC2)=NSDCg(IC2)+1
               NAMSDC(9,IC2,NSDCg(IC2))='GLN'
               NUMSDC(9,IC2,NSDCg(IC2))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCg(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALSDC(9,IC2,NSDCg(IC2))
             END IF
           END IF
         END DO
C
C           - FIND SER -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XH1-XSER).LT.DIS2 .AND.
     $         ABS(YH1-YSER).LT.DIS2 .AND.
     $         ABS(ZH1-ZSER).LT.DIS2     ) .OR.
     $        (ABS(XH2-XSER).LT.DIS2 .AND.
     $         ABS(YH2-YSER).LT.DIS2 .AND.
     $         ABS(ZH2-ZSER).LT.DIS2     ) .OR.           
     $        (ABS(XH3-XSER).LT.DIS2 .AND.
     $         ABS(YH3-YSER).LT.DIS2 .AND.
     $         ABS(ZH3-ZSER).LT.DIS2     ) .OR.           
     $        (ABS(XH4-XSER).LT.DIS2 .AND.
     $         ABS(YH4-YSER).LT.DIS2 .AND.
     $         ABS(ZH4-ZSER).LT.DIS2     ) .OR.           
     $        (ABS(XH5-XSER).LT.DIS2 .AND.
     $         ABS(YH5-YSER).LT.DIS2 .AND.
     $         ABS(ZH5-ZSER).LT.DIS2     )     ) THEN
             DISH1=SQRT((XH1-XSER)**2+(YH1-YSER)**2+(ZH1-ZSER)**2)
             DISH2=SQRT((XH2-XSER)**2+(YH2-YSER)**2+(ZH2-ZSER)**2)
             DISH3=SQRT((XH3-XSER)**2+(YH3-YSER)**2+(ZH3-ZSER)**2)
             DISH4=SQRT((XH4-XSER)**2+(YH4-YSER)**2+(ZH4-ZSER)**2)
             DISH5=SQRT((XH5-XSER)**2+(YH5-YSER)**2+(ZH5-ZSER)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             dis=min(dis,dish5)
             IF(DIS.LT.DIS2)THEN
               NSDCg(IC2)=NSDCg(IC2)+1
               NAMSDC(9,IC2,NSDCg(IC2))='SER'
               NUMSDC(9,IC2,NSDCg(IC2))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCg(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALSDC(9,IC2,NSDCg(IC2))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XH1-XTHR).LT.DIS2 .AND.
     $         ABS(YH1-YTHR).LT.DIS2 .AND.
     $         ABS(ZH1-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XH2-XTHR).LT.DIS2 .AND.
     $         ABS(YH2-YTHR).LT.DIS2 .AND.
     $         ABS(ZH2-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XH3-XTHR).LT.DIS2 .AND.
     $         ABS(YH3-YTHR).LT.DIS2 .AND.
     $         ABS(ZH3-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XH4-XTHR).LT.DIS2 .AND.
     $         ABS(YH4-YTHR).LT.DIS2 .AND.
     $         ABS(ZH4-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XH5-XTHR).LT.DIS2 .AND.
     $         ABS(YH5-YTHR).LT.DIS2 .AND.
     $         ABS(ZH5-ZTHR).LT.DIS2     )     ) THEN
             DISH1=SQRT((XH1-XTHR)**2+(YH1-YTHR)**2+(ZH1-ZTHR)**2)
             DISH2=SQRT((XH2-XTHR)**2+(YH2-YTHR)**2+(ZH2-ZTHR)**2)
             DISH3=SQRT((XH3-XTHR)**2+(YH3-YTHR)**2+(ZH3-ZTHR)**2)
             DISH4=SQRT((XH4-XTHR)**2+(YH4-YTHR)**2+(ZH4-ZTHR)**2)
             DISH5=SQRT((XH5-XTHR)**2+(YH5-YTHR)**2+(ZH5-ZTHR)**2)
             DIS=MIN(DISH1,DISH2)
             dis=min(dis,dish3)
             dis=min(dis,dish4)
             dis=min(dis,dish5)
             IF(DIS.LT.DIS2)THEN
               NSDCg(IC2)=NSDCg(IC2)+1
               NAMSDC(9,IC2,NSDCg(IC2))='THR'
               NUMSDC(9,IC2,NSDCg(IC2))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(9,IC2,NSDCg(IC2))=Foh*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)+VALSDC(9,IC2,NSDCg(IC2))
             END IF
           END IF
         END DO
c end non-ionizable
C
C           - FIND Charged Atom -
C
      call charatml(9,ic2,lligc2,xc,yc,zc,nclcg,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgcg,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO

c      **********************
c      STEP 7B. N3 atoms:
c      **********************
c N30
         FBKB=0.80
         DIS1=2.0
         DIS2=4.0
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N30') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(Xn-XBKO).LT.DIS2.AND.
     $         ABS(Yn-YBKO).LT.DIS2.AND.
     $         ABS(Zn-ZBKO).LT.DIS2    )) then
             DIS=SQRT((Xn-XBKO)**2+(Yn-YBKO)**2+(Zn-ZBKO)**2)
             IF(DIS.LT.DIS2)THEN
               NBKN3(IN3)=NBKN3(IN3)+1
               NAMBKB(7,IN3,NBKN3(IN3))=NAMPRT(I-1)
               NUMBKB(7,IN3,NBKN3(IN3))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(7,IN3,NBKN3(IN3))=FBKB*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALBKB(7,IN3,NBKN3(IN3))
             END IF
           END IF
         END DO
         end if
         END DO

c non-ionizable
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N30') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
C
C           - FIND ASN -
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF(ABS(Xn-XASN).LT.DIS2.AND.
     $        ABS(Yn-YASN).LT.DIS2.AND.
     $        ABS(Zn-ZASN).LT.DIS2    ) then
             DIS=SQRT((Xn-XASN)**2+(Yn-YASN)**2+(Zn-ZASN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='ASN'
               NUMSDC(7,IN3,NSDN3(IN3))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=FOH*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c
c           - FIND GLN -
c
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF(ABS(Xn-XGLN).LT.DIS2.AND.
     $        ABS(Yn-YGLN).LT.DIS2.AND.
     $        ABS(Zn-ZGLN).LT.DIS2    ) then
             DIS=SQRT((Xn-XGLN)**2+(Yn-YGLN)**2+(Zn-ZGLN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='GLN'
               NUMSDC(7,IN3,NSDN3(IN3))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND SER -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF(ABS(Xn-XSER).LT.DIS2 .AND.
     $        ABS(Yn-YSER).LT.DIS2 .AND.
     $        ABS(Zn-ZSER).LT.DIS2     ) then
             DIS=SQRT((Xn-XSER)**2+(Yn-YSER)**2+(Zn-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='SER'
               NUMSDC(7,IN3,NSDN3(IN3))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF(ABS(Xn-XTHR).LT.DIS2 .AND.
     $        ABS(Yn-YTHR).LT.DIS2 .AND.
     $        ABS(Zn-ZTHR).LT.DIS2     ) then
             DIS=SQRT((Xn-XTHR)**2+(Yn-YTHR)**2+(Zn-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='THR'
               NUMSDC(7,IN3,NSDN3(IN3))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c end non-ionizable
C
C           - FIND Charged Atom -
C
      call charatml(7,in3,llign3,xn,yn,zn,ncln3,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgn3,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO
         
c N31
         FBKB=0.80
         DIS1=2.0
         DIS2=4.0
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N31') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(Xn-XBKO).LT.DIS2.AND.
     $         ABS(Yn-YBKO).LT.DIS2.AND.
     $         ABS(Zn-ZBKO).LT.DIS2    )) then
             DIS=SQRT((Xn-XBKO)**2+(Yn-YBKO)**2+(Zn-ZBKO)**2)
             IF(DIS.LT.DIS2)THEN
               NBKN3(IN3)=NBKN3(IN3)+1
               NAMBKB(7,IN3,NBKN3(IN3))=NAMPRT(I-1)
               NUMBKB(7,IN3,NBKN3(IN3))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(7,IN3,NBKN3(IN3))=FBKB*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALBKB(7,IN3,NBKN3(IN3))
             END IF
           END IF
         END DO
         end if
         END DO
         
c non-ionizable
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N31') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
C
C           - FIND ASN -
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF(ABS(Xn-XASN).LT.DIS2.AND.
     $        ABS(Yn-YASN).LT.DIS2.AND.
     $        ABS(Zn-ZASN).LT.DIS2    ) then
             DIS=SQRT((Xn-XASN)**2+(Yn-YASN)**2+(Zn-ZASN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='ASN'
               NUMSDC(7,IN3,NSDN3(IN3))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=FOH*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c
c           - FIND GLN -
c
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF(ABS(Xn-XGLN).LT.DIS2.AND.
     $        ABS(Yn-YGLN).LT.DIS2.AND.
     $        ABS(Zn-ZGLN).LT.DIS2    ) then
             DIS=SQRT((Xn-XGLN)**2+(Yn-YGLN)**2+(Zn-ZGLN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='GLN'
               NUMSDC(7,IN3,NSDN3(IN3))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND SER -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF(ABS(Xn-XSER).LT.DIS2 .AND.
     $        ABS(Yn-YSER).LT.DIS2 .AND.
     $        ABS(Zn-ZSER).LT.DIS2     ) then
             DIS=SQRT((Xn-XSER)**2+(Yn-YSER)**2+(Zn-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='SER'
               NUMSDC(7,IN3,NSDN3(IN3))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF(ABS(Xn-XTHR).LT.DIS2 .AND.
     $        ABS(Yn-YTHR).LT.DIS2 .AND.
     $        ABS(Zn-ZTHR).LT.DIS2     ) then
             DIS=SQRT((Xn-XTHR)**2+(Yn-YTHR)**2+(Zn-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='THR'
               NUMSDC(7,IN3,NSDN3(IN3))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c end non-ionizable
C
C           - FIND Charged Atom -
C
      call charatml(7,in3,llign3,xn,yn,zn,ncln3,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgn3,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO
         
c N32
         FBKB=0.80
         DIS1=2.0
         DIS2=3.5
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N32') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XNH1=XN3P1(IN3)
           YNH1=YN3P1(IN3)
           ZNH1=ZN3P1(IN3)
           XNH2=XN3P2(IN3)
           YNH2=YN3P2(IN3)
           ZNH2=ZN3P2(IN3)
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(XnH1-XBKO).LT.DIS2.AND.
     $         ABS(YnH1-YBKO).LT.DIS2.AND.
     $         ABS(ZnH1-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XnH2-XBKO).LT.DIS2.AND.
     $         ABS(YnH2-YBKO).LT.DIS2.AND.
     $         ABS(ZnH2-ZBKO).LT.DIS2    )     ) THEN
             DISH1=SQRT((XnH1-XBKO)**2+(YnH1-YBKO)**2+(ZnH1-ZBKO)**2)
             DISH2=SQRT((XnH2-XBKO)**2+(YnH2-YBKO)**2+(ZnH2-ZBKO)**2)
             VECNM=SQRT((XBKC-XBKO)**2+(YBKC-YBKO)**2+(ZBKC-ZBKO)**2)
             XVCO=(XBKO-XBKC)/VECNM
             YVCO=(YBKO-YBKC)/VECNM
             ZVCO=(ZBKO-ZBKC)/VECNM
             IF(DISH1.LT.DISH2) THEN
               DIS=DISH1
               XVOH1=(XnH1-XBKO)/DISH1
               YVOH1=(YnH1-YBKO)/DISH1
               ZVOH1=(ZnH1-ZBKO)/DISH1
               AGOH=XVCO*XVOH1+YVCO*YVOH1+ZVCO*ZVOH1
             ENDIF
             IF(DISH2.LT.DISH1) THEN
               DIS=DISH2
               XVOH2=(XnH2-XBKO)/DISH2
               YVOH2=(YnH2-YBKO)/DISH2
               ZVOH2=(ZnH2-ZBKO)/DISH2
               AGOH=XVCO*XVOH2+YVCO*YVOH2+ZVCO*ZVOH2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NBKN3(IN3)=NBKN3(IN3)+1
               NAMBKB(7,IN3,NBKN3(IN3))=NAMPRT(I-1)
               NUMBKB(7,IN3,NBKN3(IN3))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(7,IN3,NBKN3(IN3))=FBKB*AGOH*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALBKB(7,IN3,NBKN3(IN3))
             END IF
           END IF
         END DO
         end if
         END DO
         
c non-ionizable
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N32') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XNH1=XN3P1(IN3)
           YNH1=YN3P1(IN3)
           ZNH1=ZN3P1(IN3)
           XNH2=XN3P2(IN3)
           YNH2=YN3P2(IN3)
           ZNH2=ZN3P2(IN3)
C
C           - FIND ASN -
C
C dcb july 2007
C At first we don't include the angle dependence
C dcb end
C
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF((ABS(XnH1-XASN).LT.DIS2.AND.
     $         ABS(YnH1-YASN).LT.DIS2.AND.
     $         ABS(ZnH1-ZASN).LT.DIS2    ) .OR.
     $        (ABS(XnH2-XASN).LT.DIS2.AND.
     $         ABS(YnH2-YASN).LT.DIS2.AND.
     $         ABS(ZnH2-ZASN).LT.DIS2    )) then
             DISH1=SQRT((XnH1-XASN)**2+(YnH1-YASN)**2+(ZnH1-ZASN)**2)
             DISH2=SQRT((XnH2-XASN)**2+(YnH2-YASN)**2+(ZnH2-ZASN)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='ASN'
               NUMSDC(7,IN3,NSDN3(IN3))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c
c           - FIND GLN -
C
C dcb july 2007
C At first we don't include the angle dependence
C dcb end
C
c
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF((ABS(XnH1-XGLN).LT.DIS2.AND.
     $         ABS(YnH1-YGLN).LT.DIS2.AND.
     $         ABS(ZnH1-ZGLN).LT.DIS2    ) .OR.
     $        (ABS(XnH2-XGLN).LT.DIS2.AND.
     $         ABS(YnH2-YGLN).LT.DIS2.AND.
     $         ABS(ZnH2-ZGLN).LT.DIS2    )) then          
             DISH1=SQRT((XnH1-XGLN)**2+(YnH1-YGLN)**2+(ZnH1-ZGLN)**2)
             DISH2=SQRT((XnH2-XGLN)**2+(YnH2-YGLN)**2+(ZnH2-ZGLN)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='GLN'
               NUMSDC(7,IN3,NSDN3(IN3))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND SER -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XnH1-XSER).LT.DIS2 .AND.
     $         ABS(YnH1-YSER).LT.DIS2 .AND.
     $         ABS(ZnH1-ZSER).LT.DIS2     ) .OR.
     $        (ABS(XnH2-XSER).LT.DIS2 .AND.
     $         ABS(YnH2-YSER).LT.DIS2 .AND.
     $         ABS(ZnH2-ZSER).LT.DIS2     )) then           
             DISH1=SQRT((XnH1-XSER)**2+(YnH1-YSER)**2+(ZnH1-ZSER)**2)
             DISH2=SQRT((XnH2-XSER)**2+(YnH2-YSER)**2+(ZnH2-ZSER)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='SER'
               NUMSDC(7,IN3,NSDN3(IN3))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XnH1-XTHR).LT.DIS2 .AND.
     $         ABS(YnH1-YTHR).LT.DIS2 .AND.
     $         ABS(ZnH1-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XnH2-XTHR).LT.DIS2 .AND.
     $         ABS(YnH2-YTHR).LT.DIS2 .AND.
     $         ABS(ZnH2-ZTHR).LT.DIS2     )) then
             DISH1=SQRT((XnH1-XTHR)**2+(YnH1-YTHR)**2+(ZnH1-ZTHR)**2)
             DISH2=SQRT((XnH2-XTHR)**2+(YnH2-YTHR)**2+(ZnH2-ZTHR)**2)
             DIS=MIN(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='THR'
               NUMSDC(7,IN3,NSDN3(IN3))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c end non-ionizable
C
C           - FIND Charged Atom -
C
      call charatml(7,in3,llign3,xn,yn,zn,ncln3,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgn3,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO

c N33
         FBKB=0.80
         DIS1=2.0
         DIS2=3.5
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N33') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XNH1=XN3P(IN3)
           YNH1=YN3P(IN3)
           ZNH1=ZN3P(IN3)
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(XnH1-XBKO).LT.DIS2.AND.
     $         ABS(YnH1-YBKO).LT.DIS2.AND.
     $         ABS(ZnH1-ZBKO).LT.DIS2    )) then
             DISH1=SQRT((XnH1-XBKO)**2+(YnH1-YBKO)**2+(ZnH1-ZBKO)**2)
             VECNM=SQRT((XBKC-XBKO)**2+(YBKC-YBKO)**2+(ZBKC-ZBKO)**2)
             XVCO=(XBKO-XBKC)/VECNM
             YVCO=(YBKO-YBKC)/VECNM
             ZVCO=(ZBKO-ZBKC)/VECNM
               DIS=DISH1
               XVOH1=(XnH1-XBKO)/DISH1
               YVOH1=(YnH1-YBKO)/DISH1
               ZVOH1=(ZnH1-ZBKO)/DISH1
               AGOH=XVCO*XVOH1+YVCO*YVOH1+ZVCO*ZVOH1
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NBKN3(IN3)=NBKN3(IN3)+1
               NAMBKB(7,IN3,NBKN3(IN3))=NAMPRT(I-1)
               NUMBKB(7,IN3,NBKN3(IN3))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(7,IN3,NBKN3(IN3))=FBKB*AGOH*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALBKB(7,IN3,NBKN3(IN3))
             END IF
           END IF
         END DO
         end if
         END DO

c non-ionizable
         DO IN3=1,NLIGN3 
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N33') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XNH1=XN3P(IN3)
           YNH1=YN3P(IN3)
           ZNH1=ZN3P(IN3)
C
C           - FIND ASN -
C
C dcb july 2007
C At first we don't include the angle dependence
C dcb end
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF(ABS(XnH1-XASN).LT.DIS2.AND.
     $        ABS(YnH1-YASN).LT.DIS2.AND.
     $        ABS(ZnH1-ZASN).LT.DIS2    ) then
             DIS=SQRT((XnH1-XASN)**2+(YnH1-YASN)**2+(ZnH1-ZASN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='ASN'
               NUMSDC(7,IN3,NSDN3(IN3))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c
c           - FIND GLN -
C
C dcb july 2007
C At first we don't include the angle dependence
C dcb end
C
c
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF(ABS(XnH1-XGLN).LT.DIS2.AND.
     $        ABS(YnH1-YGLN).LT.DIS2.AND.
     $        ABS(ZnH1-ZGLN).LT.DIS2    ) then
             DIS=SQRT((XnH1-XGLN)**2+(YnH1-YGLN)**2+(ZnH1-ZGLN)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='GLN'
               NUMSDC(7,IN3,NSDN3(IN3))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND SER -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF(ABS(XnH1-XSER).LT.DIS2 .AND.
     $        ABS(YnH1-YSER).LT.DIS2 .AND.
     $        ABS(ZnH1-ZSER).LT.DIS2     ) then
             DIS=SQRT((XnH1-XSER)**2+(YnH1-YSER)**2+(ZnH1-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='SER'
               NUMSDC(7,IN3,NSDN3(IN3))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF(ABS(XnH1-XTHR).LT.DIS2 .AND.
     $        ABS(YnH1-YTHR).LT.DIS2 .AND.
     $        ABS(ZnH1-ZTHR).LT.DIS2     ) then
             DIS=SQRT((XnH1-XTHR)**2+(YnH1-YTHR)**2+(ZnH1-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDN3(IN3)=NSDN3(IN3)+1
               NAMSDC(7,IN3,NSDN3(IN3))='THR'
               NUMSDC(7,IN3,NSDN3(IN3))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(7,IN3,NSDN3(IN3))=Foh*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
             END IF
           END IF
         END DO
c end non-ionizable
C
C           - FIND Charged Atom -
C
      call charatml(7,in3,llign3,xn,yn,zn,ncln3,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgn3,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO

c      **********************
c      STEP 7C. Oco atoms:
c      **********************
       DO IC2=1,NLIGC2
        IF(XCac(IC2).NE.0.AND.YCac(IC2).NE.0.AND.ZCac(IC2).NE.0)THEN
         XLGO1=XLO1(IC2)
         YLGO1=YLO1(IC2)
         ZLGO1=ZLO1(IC2)
         XLGO2=XLO2(IC2)
         YLGO2=YLO2(IC2)
         ZLGO2=ZLO2(IC2)
         XO=(XLGO1+XLGO2)/2.0
         YO=(YLGO1+YLGO2)/2.0
         ZO=(ZLGO1+ZLGO2)/2.0
 
c w/ O=C-
         FBKB=1.20
         DIS1=3.0
         DIS2=4.0
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(XlgO1-XBKO).LT.DIS2.AND.
     $         ABS(YlgO1-YBKO).LT.DIS2.AND.
     $         ABS(ZlgO1-ZBKO).LT.DIS2    ) .OR.
     $        (ABS(XlgO2-XBKO).LT.DIS2.AND.
     $         ABS(YlgO2-YBKO).LT.DIS2.AND.
     $         ABS(ZlgO2-ZBKO).LT.DIS2    )     ) THEN
             DISH1=SQRT((XlgO1-XBKO)**2+(YlgO1-YBKO)**2+(ZlgO1-ZBKO)**2)
             DISH2=SQRT((XlgO2-XBKO)**2+(YlgO2-YBKO)**2+(ZlgO2-ZBKO)**2)
               DIS=min(DISH1,DISH2)
             IF(DIS.LT.DIS2)THEN
               NBKLCar(Ic2)=NBKLCar(Ic2)+1
               NAMBKB(11,Ic2,NBKLCar(Ic2))=NAMPRT(I-1)
               NUMBKB(11,Ic2,NBKLCar(Ic2))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(11,Ic2,NBKLCar(Ic2))=FBKB*MIN(1.0,VALUE)
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALBKB(11,Ic2,NBKLCar(Ic2))
             END IF
           END IF
         END DO
C
c below w/ H-N- from Asp/Glu step 4
C
         FBKB=-1.20
         DIS1=2.0
         DIS2=3.5
         DO I=2,NPRTON
           XP=XPRTON(I)
           YP=YPRTON(I)
           ZP=ZPRTON(I)
           XN=XNITRN(I)
           YN=YNITRN(I)
           ZN=ZNITRN(I)
           IF((ABS(XlgO1-XP).LT.DIS2.AND.
     $         ABS(YlgO1-YP).LT.DIS2.AND.
     $         ABS(ZlgO1-ZP).LT.DIS2    ) .OR.
     $        (ABS(XlgO2-XP).LT.DIS2.AND.
     $         ABS(YlgO2-YP).LT.DIS2.AND.
     $         ABS(ZlgO2-ZP).LT.DIS2    )     ) THEN
             DISO1P=SQRT((XlgO1-XP)**2+(YlgO1-YP)**2+(ZlgO1-ZP)**2)
             DISO2P=SQRT((XlgO2-XP)**2+(YlgO2-YP)**2+(ZlgO2-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISO1P.LT.DISO2P) THEN
               DIS=DISO1P
               XVPO1=-(XP-XlgO1)/DISO1P
               YVPO1=-(YP-YlgO1)/DISO1P
               ZVPO1=-(ZP-ZlgO1)/DISO1P
               AGPO=XVNP*XVPO1 + YVNP*YVPO1 + ZVNP*ZVPO1
             ENDIF
             IF(DISO2P.LT.DISO1P) THEN
               DIS=DISO2P
               XVPO2=-(XP-XlgO2)/DISO2P
               YVPO2=-(YP-YlgO2)/DISO2P
               ZVPO2=-(ZP-ZlgO2)/DISO2P
               AGPO=XVNP*XVPO2 + YVNP*YVPO2 + ZVNP*ZVPO2
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPO.GT.0.001)THEN
               NBKLCar(Ic2)=NBKLCar(Ic2)+1
               NAMBKB(11,Ic2,NBKLCar(Ic2))=NAMPRT(I)
               NUMBKB(11,Ic2,NBKLCar(Ic2))=NUMPRT(I)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALUE=MIN(1.0,VALUE)
               VALBKB(11,Ic2,NBKLCar(Ic2))=FBKB*VALUE*AGPO
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALBKB(11,Ic2,NBKLCar(Ic2))
             END IF
           END IF
         END DO
C
         end if
         END DO
         
c non-ionizable
       DO IC2=1,NLIGC2
        IF(XCac(IC2).NE.0.AND.YCac(IC2).NE.0.AND.ZCac(IC2).NE.0)THEN
         XO1=XLO1(IC2)
         YO1=YLO1(IC2)
         ZO1=ZLO1(IC2)
         XO2=XLO2(IC2)
         YO2=YLO2(IC2)
         ZO2=ZLO2(IC2)
c below from step 4 Asp/Glu 
C
C           - FIND SER -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF((ABS(XO1-XSER).LT.DIS2 .AND.
     $         ABS(YO1-YSER).LT.DIS2 .AND.
     $         ABS(ZO1-ZSER).LT.DIS2     ) .OR.
     $        (ABS(XO2-XSER).LT.DIS2 .AND.
     $         ABS(YO2-YSER).LT.DIS2 .AND.
     $         ABS(ZO2-ZSER).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XSER)**2+(YO1-YSER)**2+(ZO1-ZSER)**2)
             DISO2O=SQRT((XO2-XSER)**2+(YO2-YSER)**2+(ZO2-ZSER)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDLCar(Ic2)=NSDLCar(Ic2)+1
               NAMSDC(11,Ic2,NSDLCar(Ic2))='SER'
               NUMSDC(11,Ic2,NSDLCar(Ic2))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,Ic2,NSDLCar(Ic2))=Foh*MIN(1.0,VALUE)
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALSDC(11,Ic2,NSDLCar(Ic2))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF((ABS(XO1-XTHR).LT.DIS2 .AND.
     $         ABS(YO1-YTHR).LT.DIS2 .AND.
     $         ABS(ZO1-ZTHR).LT.DIS2     ) .OR.
     $        (ABS(XO2-XTHR).LT.DIS2 .AND.
     $         ABS(YO2-YTHR).LT.DIS2 .AND.
     $         ABS(ZO2-ZTHR).LT.DIS2     )     ) THEN
             DISO1O=SQRT((XO1-XTHR)**2+(YO1-YTHR)**2+(ZO1-ZTHR)**2)
             DISO2O=SQRT((XO2-XTHR)**2+(YO2-YTHR)**2+(ZO2-ZTHR)**2)
             DIS=MIN(DISO1O,DISO2O)
             IF(DIS.LT.DIS2)THEN
               NSDLCar(Ic2)=NSDLCar(Ic2)+1
               NAMSDC(11,Ic2,NSDLCar(Ic2))='THR'
               NUMSDC(11,Ic2,NSDLCar(Ic2))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,Ic2,NSDLCar(Ic2))=Foh*MIN(1.0,VALUE)
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALSDC(11,Ic2,NSDLCar(Ic2))
             END IF
           END IF
         END DO
C
C
C           - FIND ASN -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN1=XASNP1(IASN)
           YASN1=YASNP1(IASN)
           ZASN1=ZASNP1(IASN)
           XASN2=XASNP2(IASN)
           YASN2=YASNP2(IASN)
           ZASN2=ZASNP2(IASN)
           IF((ABS(XO1-XASN1).LT.DIS2.AND.
     $         ABS(YO1-YASN1).LT.DIS2.AND.
     $         ABS(ZO1-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XASN1).LT.DIS2.AND.
     $         ABS(YO2-YASN1).LT.DIS2.AND.
     $         ABS(ZO2-ZASN1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XASN2).LT.DIS2.AND.
     $         ABS(YO1-YASN2).LT.DIS2.AND.
     $         ABS(ZO1-ZASN2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XASN2).LT.DIS2.AND.
     $         ABS(YO2-YASN2).LT.DIS2.AND.
     $         ABS(ZO2-ZASN2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XASN1)**2+(YO1-YASN1)**2+(ZO1-ZASN1)**2)
             DISO21=SQRT((XO2-XASN1)**2+(YO2-YASN1)**2+(ZO2-ZASN1)**2)
             DISO12=SQRT((XO1-XASN2)**2+(YO1-YASN2)**2+(ZO1-ZASN2)**2)
             DISO22=SQRT((XO2-XASN2)**2+(YO2-YASN2)**2+(ZO2-ZASN2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDLCar(Ic2)=NSDLCar(Ic2)+1
               NAMSDC(11,Ic2,NSDLCar(Ic2))='ASN'
               NUMSDC(11,Ic2,NSDLCar(Ic2))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,Ic2,NSDLCar(Ic2))=Fnh*MIN(1.0,VALUE)
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALSDC(11,Ic2,NSDLCar(Ic2))
             END IF
           END IF
         END DO
C
C
C           - FIND GLN -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN1=XGLNP1(IGLN)
           YGLN1=YGLNP1(IGLN)
           ZGLN1=ZGLNP1(IGLN)
           XGLN2=XGLNP2(IGLN)
           YGLN2=YGLNP2(IGLN)
           ZGLN2=ZGLNP2(IGLN)
           IF((ABS(XO1-XGLN1).LT.DIS2.AND.
     $         ABS(YO1-YGLN1).LT.DIS2.AND.
     $         ABS(ZO1-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XGLN1).LT.DIS2.AND.
     $         ABS(YO2-YGLN1).LT.DIS2.AND.
     $         ABS(ZO2-ZGLN1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XGLN2).LT.DIS2.AND.
     $         ABS(YO1-YGLN2).LT.DIS2.AND.
     $         ABS(ZO1-ZGLN2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XGLN2).LT.DIS2.AND.
     $         ABS(YO2-YGLN2).LT.DIS2.AND.
     $         ABS(ZO2-ZGLN2).LT.DIS2    )     ) THEN
             DISO11=SQRT((XO1-XGLN1)**2+(YO1-YGLN1)**2+(ZO1-ZGLN1)**2)
             DISO21=SQRT((XO2-XGLN1)**2+(YO2-YGLN1)**2+(ZO2-ZGLN1)**2)
             DISO12=SQRT((XO1-XGLN2)**2+(YO1-YGLN2)**2+(ZO1-ZGLN2)**2)
             DISO22=SQRT((XO2-XGLN2)**2+(YO2-YGLN2)**2+(ZO2-ZGLN2)**2)
             DIS=MIN(DISO11,DISO21)
             DIS=MIN(DIS,DISO12)
             DIS=MIN(DIS,DISO22)
             IF(DIS.LT.DIS2)THEN
               NSDLCar(Ic2)=NSDLCar(Ic2)+1
               NAMSDC(11,Ic2,NSDLCar(Ic2))='GLN'
               NUMSDC(11,Ic2,NSDLCar(Ic2))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,Ic2,NSDLCar(Ic2))=Fnh*MIN(1.0,VALUE)
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALSDC(11,Ic2,NSDLCar(Ic2))
             END IF
           END IF
         END DO
C
C           - FIND TRP -
C
         FNH=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ITRP=1,NTRP
           XTRP=XTRPP1(ITRP)
           YTRP=YTRPP1(ITRP)
           ZTRP=ZTRPP1(ITRP)
           IF((ABS(XO1-XTRP).LT.DIS2 .AND.
     $         ABS(YO1-YTRP).LT.DIS2 .AND.
     $         ABS(ZO1-ZTRP).LT.DIS2     ) .OR.
     $        (ABS(XO2-XTRP).LT.DIS2 .AND.
     $         ABS(YO2-YTRP).LT.DIS2 .AND.
     $         ABS(ZO2-ZTRP).LT.DIS2     )     ) THEN
             DISO1P=SQRT((XO1-XTRP)**2+(YO1-YTRP)**2+(ZO1-ZTRP)**2)
             DISO2P=SQRT((XO2-XTRP)**2+(YO2-YTRP)**2+(ZO2-ZTRP)**2)
             DIS=MIN(DISO1P,DISO2P)
             IF(DIS.LT.DIS2)THEN
               NSDLCar(Ic2)=NSDLCar(Ic2)+1
               NAMSDC(11,Ic2,NSDLCar(Ic2))='TRP'
               NUMSDC(11,Ic2,NSDLCar(Ic2))=LTRPRS(ITRP)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,Ic2,NSDLCar(Ic2))=Fnh*MIN(1.0,VALUE)
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALSDC(11,Ic2,NSDLCar(Ic2))
             END IF
           END IF
         END DO
C
c end non-ionizable
c w/ Lys and Arg from Step 4
C
C           - FIND LYS H-BONDING
C
         FNH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ILYS=1,NLYS
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           DIS2=4.00
           IF(ILYS.EQ.1)FNH=-1.20
           IF(ILYS.EQ.1)DIS2=4.50
           IF((ABS(XO1-XN).LT.DIS2.AND.
     $         ABS(YO1-YN).LT.DIS2.AND.
     $         ABS(ZO1-ZN).LT.DIS2    ) .OR.
     $        (ABS(XO2-XN).LT.DIS2.AND.
     $         ABS(YO2-YN).LT.DIS2.AND.
     $         ABS(ZO2-ZN).LT.DIS2    )     ) THEN
             DISO1N=SQRT((XO1-XN)**2+(YO1-YN)**2+(ZO1-ZN)**2)
             DISO2N=SQRT((XO2-XN)**2+(YO2-YN)**2+(ZO2-ZN)**2)
             DIS=MIN(DISO1N,DISO2N)
             IF(DIS.LT.DIS2)THEN
c              NSDCAR(ICAR)=NSDCAR(ICAR)+1
c              NAMSDC(1,ICAR,NSDCAR(ICAR))='LYS'
c              IF(ILYS.EQ.1)NAMSDC(1,ICAR,NSDCAR(ICAR))='N+ '
c              NUMSDC(1,ICAR,NSDCAR(ICAR))=LLYSRS(ILYS)
c              VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
c              VALSDC(1,ICAR,NSDCAR(ICAR))=FNH*MIN(1.0,VALUE)
c              PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
               NSDLCar(Ic2)=NSDLCar(Ic2)+1
               NAMSDC(11,Ic2,NSDLCar(Ic2))='LYS'
               IF(ILYS.EQ.1)NAMSDC(11,IC2,NSDlCAR(IC2))='N+ '
               NUMSDC(11,Ic2,NSDLCar(Ic2))=LlysRS(Ilys)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,Ic2,NSDLCar(Ic2))=Fnh*MIN(1.0,VALUE)
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALSDC(11,Ic2,NSDLCar(Ic2))
             END IF
           END IF
         END DO
C
C           - FIND ARG H-BONDING
C
         FNH=-0.80
         DIS1=2.00
         DIS2=4.00
         DO IARG=1,NARG
           XH1=XARGP1(IARG)
           YH1=YARGP1(IARG)
           ZH1=ZARGP1(IARG)
           XH2=XARGP2(IARG)
           YH2=YARGP2(IARG)
           ZH2=ZARGP2(IARG)
           XH3=XARGP3(IARG)
           YH3=YARGP3(IARG)
           ZH3=ZARGP3(IARG)
           XH4=XARGP4(IARG)
           YH4=YARGP4(IARG)
           ZH4=ZARGP4(IARG)
           XH5=XARGP5(IARG)
           YH5=YARGP5(IARG)
           ZH5=ZARGP5(IARG)
           IF(
     $        (ABS(XO1-XH1).LT.DIS2.AND.
     $         ABS(YO1-YH1).LT.DIS2.AND.
     $         ABS(ZO1-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH1).LT.DIS2.AND.
     $         ABS(YO2-YH1).LT.DIS2.AND.
     $         ABS(ZO2-ZH1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH2).LT.DIS2.AND.
     $         ABS(YO1-YH2).LT.DIS2.AND.
     $         ABS(ZO1-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH2).LT.DIS2.AND.
     $         ABS(YO2-YH2).LT.DIS2.AND.
     $         ABS(ZO2-ZH2).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH3).LT.DIS2.AND.
     $         ABS(YO1-YH3).LT.DIS2.AND.
     $         ABS(ZO1-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH3).LT.DIS2.AND.
     $         ABS(YO2-YH3).LT.DIS2.AND.
     $         ABS(ZO2-ZH3).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH4).LT.DIS2.AND.
     $         ABS(YO1-YH4).LT.DIS2.AND.
     $         ABS(ZO1-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH4).LT.DIS2.AND.
     $         ABS(YO2-YH4).LT.DIS2.AND.
     $         ABS(ZO2-ZH4).LT.DIS2    ) .OR.
     $        (ABS(XO1-XH5).LT.DIS2.AND.
     $         ABS(YO1-YH5).LT.DIS2.AND.
     $         ABS(ZO1-ZH5).LT.DIS2    ) .OR.
     $        (ABS(XO2-XH5).LT.DIS2.AND.
     $         ABS(YO2-YH5).LT.DIS2.AND.
     $         ABS(ZO2-ZH5).LT.DIS2    )     ) THEN
             DIS11=SQRT((XO1-XH1)**2+(YO1-YH1)**2+(ZO1-ZH1)**2)
             DIS12=SQRT((XO1-XH2)**2+(YO1-YH2)**2+(ZO1-ZH2)**2)
             DIS13=SQRT((XO1-XH3)**2+(YO1-YH3)**2+(ZO1-ZH3)**2)
             DIS14=SQRT((XO1-XH4)**2+(YO1-YH4)**2+(ZO1-ZH4)**2)
             DIS15=SQRT((XO1-XH5)**2+(YO1-YH5)**2+(ZO1-ZH5)**2)
             DIS21=SQRT((XO2-XH1)**2+(YO2-YH1)**2+(ZO2-ZH1)**2)
             DIS22=SQRT((XO2-XH2)**2+(YO2-YH2)**2+(ZO2-ZH2)**2)
             DIS23=SQRT((XO2-XH3)**2+(YO2-YH3)**2+(ZO2-ZH3)**2)
             DIS24=SQRT((XO2-XH4)**2+(YO2-YH4)**2+(ZO2-ZH4)**2)
             DIS25=SQRT((XO2-XH5)**2+(YO2-YH5)**2+(ZO2-ZH5)**2)
             DIS=MIN(DIS11,DIS12)
             DIS=MIN(DIS,DIS13)
             DIS=MIN(DIS,DIS14)
             DIS=MIN(DIS,DIS15)
             DIS=MIN(DIS,DIS21)
             DIS=MIN(DIS,DIS22)
             DIS=MIN(DIS,DIS23)
             DIS=MIN(DIS,DIS24)
             DIS=MIN(DIS,DIS25)
             IF(DIS.LT.DIS2)THEN
               NSDLCar(Ic2)=NSDLCar(Ic2)+1
               NAMSDC(11,Ic2,NSDLCar(Ic2))='ARG'
               NUMSDC(11,Ic2,NSDLCar(Ic2))=LargRS(Iarg)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(11,Ic2,NSDLCar(Ic2))=Fnh*MIN(1.0,VALUE)
C
C     -- THERE ARE POSSIBLY 2 H-BONDS BETWEEN ARG AND CARBOXYL
               IF((DIS11.LT.2.2 .AND. DIS22.LT.2.2).OR.
     $            (DIS21.LT.2.2 .AND. DIS12.LT.2.2).OR.
     $            (DIS13.LT.2.2 .AND. DIS24.LT.2.2).OR.
     $            (DIS23.LT.2.2 .AND. DIS14.LT.2.2)    )THEN
               VALSDC(11,Ic2,NSDLCar(Ic2))=-2.40
               END IF
               PK1LGcar(Ic2)=PK1LGcar(Ic2)+VALSDC(11,Ic2,NSDLCar(Ic2))
             END IF
           END IF
         END DO
C
C           - FIND Charged Atom -
C
      call charatml(11,ic2,lligc2,xo,yo,zo,ncllgcar,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgcar,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO
C
c      **********************
c      STEP 7D. Nar atoms:
c      **********************
c Nar+ (basic w/ H) i.e. protonated basic. Non-basic not titratable.
c        FBKB=0.80
c FBKB * 2 to account for N: O=C-N-H interaction.
         FBKB=1.60
         DIS1=2.0
         DIS2=3.5
         DO INar=1,NLIGNar
          IF (NAMATM(LLIGNar(INar)).EQ.'  Nar  ') THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           XP=XgLIG(INar)
           YP=YgLIG(INar)
           ZP=ZgLIG(INar)
         DO I=2,NPRTON
           XBKO=XOXYGN(I)
           YBKO=YOXYGN(I)
           ZBKO=ZOXYGN(I)
           XBKC=XCARBN(I)
           YBKC=YCARBN(I)
           ZBKC=ZCARBN(I)
           IF((ABS(Xp-XBKO).LT.DIS2.AND.
     $         ABS(Yp-YBKO).LT.DIS2.AND.
     $         ABS(Zp-ZBKO).LT.DIS2    )) then
             DISH1=SQRT((Xp-XBKO)**2+(Yp-YBKO)**2+(Zp-ZBKO)**2)
             VECNM=SQRT((XBKC-XBKO)**2+(YBKC-YBKO)**2+(ZBKC-ZBKO)**2)
             XVCO=(XBKO-XBKC)/VECNM
             YVCO=(YBKO-YBKC)/VECNM
             ZVCO=(ZBKO-ZBKC)/VECNM
               DIS=DISH1
               XVOH1=(Xp-XBKO)/DISH1
               YVOH1=(Yp-YBKO)/DISH1
               ZVOH1=(Zp-ZBKO)/DISH1
               AGOH=XVCO*XVOH1+YVCO*YVOH1+ZVCO*ZVOH1
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NBKNar(INar)=NBKNar(INar)+1
               NAMBKB(8,INar,NBKNar(INar))=NAMPRT(I-1)
               NUMBKB(8,INar,NBKNar(INar))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(8,INar,NBKNar(INar))=FBKB*AGOH*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALBKB(8,INar,NBKNar(INar))
             END IF
           END IF
         END DO
         end if
         END DO

c N-H w/ H-N-C=O bb
         FBKB=-0.80
         DIS1=1.0
         DIS2=2.0
         DO INar=1,NLIGNar
          IF (NAMATM(LLIGNar(INar)).EQ.'  Nar  ') THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           XP=XgLIG(INar)
           YP=YgLIG(INar)
           ZP=ZgLIG(INar)
         DO I=2,NPRTON
           XBKN=XNITRN(I)
           YBKN=YNITRN(I)
           ZBKN=ZNITRN(I)
           XBKH=XPRTON(I)
           YBKH=YPRTON(I)
           ZBKH=ZPRTON(I)
           IF((ABS(Xp-XBKH).LT.DIS2.AND.
     $         ABS(Yp-YBKH).LT.DIS2.AND.
     $         ABS(Zp-ZBKH).LT.DIS2    )) then
             DIS=SQRT((Xp-XBKH)**2+(Yp-YBKH)**2+(Zp-ZBKH)**2)
             IF(DIS.LT.DIS2)THEN
               NBKNar(INar)=NBKNar(INar)+1
               NAMBKB(8,INar,NBKNar(INar))=NAMPRT(I-1)
               NUMBKB(8,INar,NBKNar(INar))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALBKB(8,INar,NBKNar(INar))=FBKB*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALBKB(8,INar,NBKNar(INar))
             END IF
           END IF
         END DO
         end if
         END DO

c N: w/ H-N-C=O bb
         FBKB=-0.80
         DIS1=2.0
         DIS2=3.5
         DO INar=1,NLIGNar
          IF (NAMATM(LLIGNar(INar)).EQ.'  Nar  ') THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           XP=XgLIG(INar)
           YP=YgLIG(INar)
           ZP=ZgLIG(INar)
         DO I=2,NPRTON
           XBKN=XNITRN(I)
           YBKN=YNITRN(I)
           ZBKN=ZNITRN(I)
           XBKH=XPRTON(I)
           YBKH=YPRTON(I)
           ZBKH=ZPRTON(I)
           IF((ABS(XN-XBKH).LT.DIS2.AND.
     $         ABS(YN-YBKH).LT.DIS2.AND.
     $         ABS(ZN-ZBKH).LT.DIS2    )) then
             DISH1=SQRT((XN-XBKH)**2+(YN-YBKH)**2+(ZN-ZBKH)**2)
c            VECNM=SQRT((XBKn-XBKh)**2+(YBKn-YBKh)**2+(ZBKn-ZBKh)**2)
c            XVnh=(XBKh-XBKn)/VECNM
c            YVnh=(YBKh-YBKn)/VECNM
c            ZVnh=(ZBKh-ZBKn)/VECNM
               DIS=DISH1
c              XVhH1=(Xn-XBKh)/DISH1
c              YVhH1=(Yn-YBKh)/DISH1
c              ZVhH1=(Zn-ZBKh)/DISH1
c              AGhH=XVnh*XVhH1+YVnh*YVhH1+ZVnh*ZVhH1
c            IF(DIS.LT.DIS2 .AND. AGhH.GT.0.001)THEN
             IF(DIS.LT.DIS2)THEN
               NBKNar(INar)=NBKNar(INar)+1
               NAMBKB(8,INar,NBKNar(INar))=NAMPRT(I-1)
               NUMBKB(8,INar,NBKNar(INar))=NUMPRT(I-1)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
c              VALBKB(8,INar,NBKNar(INar))=FBKB*AGhH*MIN(1.0,VALUE)
               VALBKB(8,INar,NBKNar(INar))=FBKB*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALBKB(8,INar,NBKNar(INar))
             END IF
           END IF
         END DO
         end if
         END DO

c non-ionizable
         DO INar=1,NLIGNar
          IF (NAMATM(LLIGNar(INar)).EQ.'  Nar  ') THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           XP=XgLIG(INar)
           YP=YgLIG(INar)
           ZP=ZgLIG(INar)
C
C           - FIND ASN -
C
C dcb july 2007
C The angle contribution is now included
C dcb july 2007
C
C
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IASN=1,NASN
           XASN=X(LASNOD(IASN))
           YASN=Y(LASNOD(IASN))
           ZASN=Z(LASNOD(IASN))
           IF(ABS(Xp-XASN).LT.DIS2.AND.
     $        ABS(Yp-YASN).LT.DIS2.AND.
     $        ABS(Zp-ZASN).LT.DIS2    ) then
             DIS=SQRT((Xp-XASN)**2+(Yp-YASN)**2+(Zp-ZASN)**2)
             VECNM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XNP=(XP-XN)/VECNM
             YNP=(YP-YN)/VECNM
             ZNP=(ZP-ZN)/VECNM
C
             XVPOD=-(XASN-XP)/DIS
             YVPOD=-(YASN-YP)/DIS
             ZVPOD=-(ZASN-ZP)/DIS
             AGOH=XNP*XVPOD+YNP*YVPOD+ZNP*ZVPOD
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NSDNar(INar)=NSDNar(INar)+1
               NAMSDC(8,INar,NSDNar(INar))='ASN'
               NUMSDC(8,INar,NSDNar(INar))=LASNRS(IASN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(8,INar,NSDNar(INar))=Foh*MIN(1.0,VALUE)*AGOH
               PK1LGNar(INar)=PK1LGNar(INar)+VALSDC(8,INar,NSDNar(INar))
             END IF
           END IF
         END DO
c
c           - FIND GLN -
C
C dcb july 2007
C The angle contribution is now included
C dcb july 2007
C
c
         FOH=0.80
         DIS1=2.00
         DIS2=3.00
         DO IGLN=1,NGLN
           XGLN=X(LGLNOE(IGLN))
           YGLN=Y(LGLNOE(IGLN))
           ZGLN=Z(LGLNOE(IGLN))
           IF(ABS(Xp-XGLN).LT.DIS2.AND.
     $        ABS(Yp-YGLN).LT.DIS2.AND.
     $        ABS(Zp-ZGLN).LT.DIS2    ) then
             DIS=SQRT((Xp-XGLN)**2+(Yp-YGLN)**2+(Zp-ZGLN)**2)
              VECNM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
              XNP=(XP-XN)/VECNM
              YNP=(YP-YN)/VECNM
              ZNP=(ZP-ZN)/VECNM
C
              XVPOE=-(XGLN-XP)/DIS
              YVPOE=-(YGLN-YP)/DIS
              ZVPOE=-(ZGLN-ZP)/DIS
             AGOH=XNP*XVPOE+YNP*YVPOE+ZNP*ZVPOE
             IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               NSDNar(INar)=NSDNar(INar)+1
               NAMSDC(8,INar,NSDNar(INar))='GLN'
               NUMSDC(8,INar,NSDNar(INar))=LGLNRS(IGLN)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(8,INar,NSDNar(INar))=Foh*MIN(1.0,VALUE)*AGOH
               PK1LGNar(INar)=PK1LGNar(INar)+VALSDC(8,INar,NSDNar(INar))
             END IF
           END IF
         END DO
C
C           - FIND SER -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ISER=1,NSER
           XSER=X(LSEROH(ISER))
           YSER=Y(LSEROH(ISER))
           ZSER=Z(LSEROH(ISER))
           IF(ABS(Xp-XSER).LT.DIS2 .AND.
     $        ABS(Yp-YSER).LT.DIS2 .AND.
     $        ABS(Zp-ZSER).LT.DIS2     ) then
             DIS=SQRT((Xp-XSER)**2+(Yp-YSER)**2+(Zp-ZSER)**2)
             IF(DIS.LT.DIS2)THEN
               NSDNar(INar)=NSDNar(INar)+1
               NAMSDC(8,INar,NSDNar(INar))='SER'
               NUMSDC(8,INar,NSDNar(INar))=LSERRS(ISER)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(8,INar,NSDNar(INar))=Foh*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALSDC(8,INar,NSDNar(INar))
             END IF
           END IF
         END DO
C
C           - FIND THR -
C
         FOH=0.80
         DIS1=3.00
         DIS2=4.00
         DO ITHR=1,NTHR
           XTHR=X(LTHROH(ITHR))
           YTHR=Y(LTHROH(ITHR))
           ZTHR=Z(LTHROH(ITHR))
           IF(ABS(Xp-XTHR).LT.DIS2 .AND.
     $        ABS(Yp-YTHR).LT.DIS2 .AND.
     $        ABS(Zp-ZTHR).LT.DIS2     ) then
             DIS=SQRT((Xp-XTHR)**2+(Yp-YTHR)**2+(Zp-ZTHR)**2)
             IF(DIS.LT.DIS2)THEN
               NSDNar(INar)=NSDNar(INar)+1
               NAMSDC(8,INar,NSDNar(INar))='THR'
               NUMSDC(8,INar,NSDNar(INar))=LTHRRS(ITHR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(8,INar,NSDNar(INar))=Foh*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALSDC(8,INar,NSDNar(INar))
             END IF
           END IF
         END DO
c end non-ionizable
C
C           - FIND Charged Atom -
C
      call charatml(8,inar,llignar,xn,yn,zn,nclnar,namcol,namres,
     $ nblcol,nbatm,valcol,pk1lgnar,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
         end if
         END DO
c dmr end of ligand-backbone; ligand w/ non-ionizable 
c dmr str
c      **********************
c      STEP 7E. LYS 
c      **********************
       DO ILYS=1, NLYS
         XNZ=X(LLYSNZ(ILYS))
         YNZ=Y(LLYSNZ(ILYS))
         ZNZ=Z(LLYSNZ(ILYS))
C
C           - FIND Charged Atom -
C
      call charatm(5,ilys,xnz,ynz,znz,nclllys,namlcol,namres,
     $ nblcol,nbatm,vallcol,pk1lys,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
       END DO
c
c      **********************
c      STEP 7F. ARG 
c      **********************
       DO IARG=1, NARG
         X1=X(LARGN1(IARG))
         Y1=Y(LARGN1(IARG))
         Z1=Z(LARGN1(IARG))
         X2=X(LARGN2(IARG))
         Y2=Y(LARGN2(IARG))
         Z2=Z(LARGN2(IARG))
         X3=X(LARGN3(IARG))
         Y3=Y(LARGN3(IARG))
         Z3=Z(LARGN3(IARG))
         XCZ=X(LARGCZ(IARG))
         YCZ=Y(LARGCZ(IARG))
         ZCZ=Z(LARGCZ(IARG))
C
C           - FIND Charged Atom -
C
      call charatm(6,iarg,xcz,ycz,zcz,ncllarg,namlcol,namres,
     $ nblcol,nbatm,vallcol,pk1arg,nmass,nligchr,lligchr,
     $ lligchrval,x,y,z)
C
       END DO
c dmr fin
C
C
       DO I=1, 1000
         PK2CAR(I)=PK1CAR(I)
         PK2HIS(I)=PK1HIS(I)
         PK2CYS(I)=PK1CYS(I)
         PK2TYR(I)=PK1TYR(I)
         PK2LYS(I)=PK1LYS(I)
         PK2ARG(I)=PK1ARG(I)
         PK2LGN3(I)=PK1LGN3(I)
         PK2LGNar(I)=PK1LGNar(I)
         PK2LGCN2(I)=PK1LGCN2(I)
         PK2LGCg(I)=PK1LGCg(I)
         PK2LGCAR(I)=PK1LGCAR(I)
         PKACAR(I)=PK1CAR(I)
         PKAHIS(I)=PK1HIS(I)
         PKACYS(I)=PK1CYS(I)
         PKATYR(I)=PK1TYR(I)
         PKALYS(I)=PK1LYS(I)
         PKAARG(I)=PK1ARG(I)
         PKALGN3(I)=PK1LGN3(I)
         PKALGNar(I)=PK1LGNar(I)
         PKALGCg(I)=PK1LGCg(I)
         PKALGCN2(I)=PK1LGCN2(I)
         PKALGCAR(I)=PK1LGCAR(I)
       END DO
       DO I=501,1000
         NSDCAR(I)=NSDCAR(I-500)
         NBKCAR(I)=NBKCAR(I-500)
         NCLCAR(I)=NCLCAR(I-500)
         NSDHIS(I)=NSDHIS(I-500)
         NBKHIS(I)=NBKHIS(I-500)
         NCLHIS(I)=NCLHIS(I-500)
         NSDCYS(I)=NSDCYS(I-500)
         NBKCYS(I)=NBKCYS(I-500)
         NCLCYS(I)=NCLCYS(I-500)
         NSDTYR(I)=NSDTYR(I-500)
         NBKTYR(I)=NBKTYR(I-500)
         NCLTYR(I)=NCLTYR(I-500)
         NSDLYS(I)=NSDLYS(I-500)
         NBKLYS(I)=NBKLYS(I-500)
         NCLLYS(I)=NCLLYS(I-500)
         NSDARG(I)=NSDARG(I-500)
         NBKARG(I)=NBKARG(I-500)
         NCLARG(I)=NCLARG(I-500)
         NlgCAR(I)=NlgCAR(I-500)
         NSDLCAR(I)=NSDLCAR(I-500)
         NCLLCAR(I)=NCLLCAR(I-500)
         NCLLCYS(I)=NCLLCYS(I-500)
         NCLLTYR(I)=NCLLTYR(I-500)
         NCLLHIS(I)=NCLLHIS(I-500)
c
         NCLLlys(I)=NCLLlys(I-500)
         NCLLarg(I)=NCLLarg(I-500)
c
         NCLLGCAR(I)=NCLLGCAR(I-500)
         NSDN3(I)=NSDN3(I-500)
         NCLN3(I)=NCLN3(I-500)
         NSDNar(I)=NSDNar(I-500)
         NCLNar(I)=NCLNar(I-500)
         NCLCN2(I)=NCLCN2(I-500)
         NCLCg(I)=NCLCg(I-500)
c
         NclNpl(I)=NclNpl(I-500)
         NsdCN2(I)=NsdCN2(I-500)
         NsdCg(I)=NsdCg(I-500)
         NsdNpl(I)=NsdNpl(I-500)
         NbkCN2(I)=NbkCN2(I-500)
         NbkCg(I)=NbkCg(I-500)
         NbkNpl(I)=NbkNpl(I-500)
       END DO

C
C
C      **********************
C      STEP 8. PKA ITERATION
C      **********************
C         - USUALLY 1 TO 3 ITERATIONS ARE REQUIRED
C
       DO 500 ITER=1,10
C
       DO I=1, 1000
         PK1CAR(I)=PK2CAR(I)
         PK1HIS(I)=PK2HIS(I)
         PK1CYS(I)=PK2CYS(I)
         PK1TYR(I)=PK2TYR(I)
         PK1LYS(I)=PK2LYS(I)
         PK1ARG(I)=PK2ARG(I)
         PK1LGN3(I)=PK2LGN3(I)
         PK1LGNar(I)=PK2LGNar(I)
         PK1LGCN2(I)=PK2LGCN2(I)
         PK1LGCg(I)=PK2LGCg(I)
         PK1LGCAR(I)=PK2LGCAR(I)
       END DO
       DO I=1,500
         NSDCAR(I)=NSDCAR(I+500)
         NBKCAR(I)=NBKCAR(I+500)
         NCLCAR(I)=NCLCAR(I+500)
         NSDHIS(I)=NSDHIS(I+500)
         NBKHIS(I)=NBKHIS(I+500)
         NCLHIS(I)=NCLHIS(I+500)
         NSDCYS(I)=NSDCYS(I+500)
         NBKCYS(I)=NBKCYS(I+500)
         NCLCYS(I)=NCLCYS(I+500)
         NSDTYR(I)=NSDTYR(I+500)
         NBKTYR(I)=NBKTYR(I+500)
         NCLTYR(I)=NCLTYR(I+500)
         NSDLYS(I)=NSDLYS(I+500)
         NBKLYS(I)=NBKLYS(I+500)
         NCLLYS(I)=NCLLYS(I+500)
         NSDARG(I)=NSDARG(I+500)
         NBKARG(I)=NBKARG(I+500)
         NCLARG(I)=NCLARG(I+500)
         NlgCAR(I)=NlgCAR(I+500)
         NSDLCAR(I)=NSDLCAR(I+500)
         NCLLCAR(I)=NCLLCAR(I+500)
         NCLLCYS(I)=NCLLCYS(I+500)
         NCLLTYR(I)=NCLLTYR(I+500)
         NCLLHIS(I)=NCLLHIS(I+500)
c
         NCLLlys(I)=NCLLlys(I+500)
         NCLLarg(I)=NCLLarg(I+500)
c
         NCLLGCAR(I)=NCLLGCAR(I+500)
         NSDN3(I)=NSDN3(I+500)
         NCLN3(I)=NCLN3(I+500)
         NSDNar(I)=NSDNar(I+500)
         NCLNar(I)=NCLNar(I+500)
         NCLCN2(I)=NCLCN2(I+500)
         NCLCg(I)=NCLCg(I+500)
c
         NclNpl(I)=NclNpl(I+500)
         NsdCN2(I)=NsdCN2(I+500)
         NsdCg(I)=NsdCg(I+500)
         NsdNpl(I)=NsdNpl(I+500)
         NbkCN2(I)=NbkCN2(I+500)
         NbkCg(I)=NbkCg(I+500)
         NbkNpl(I)=NbkNpl(I+500)
       END DO
C
C      **********************
C      ASP/GLU 
C      **********************
C
C
       DO ICAR=1, NCAR
         K=NSDCAR(ICAR)
         DO J=K+1,30
           NAMSDC(1,ICAR,J)='000'
           NUMSDC(1,ICAR,J)=0
           VALSDC(1,ICAR,J)=0.0
         END DO
C
         K=NBKCAR(ICAR)
         DO J=K+1,30
           NAMBKB(1,ICAR,J)='000'
           NUMBKB(1,ICAR,J)=0
           VALBKB(1,ICAR,J)=0.0
         END DO
C
         K=NCLCAR(ICAR)
         DO J=K+1,30
           NAMCOL(1,ICAR,J)='000'
           NUMCOL(1,ICAR,J)=0
           VALCOL(1,ICAR,J)=0.0
         END DO
C
         K=NLGCAR(ICAR)
         DO J=K+1,30
           NAMLIG(1,ICAR,J)='000'
           LABEL(1,ICAR,J)='0000'
           VALLIG(1,ICAR,J)=0.0
         END DO
C
         K=NCLLCAR(ICAR)
         DO J=K+1,30
           NAMLCOL(1,ICAR,J)='000'
           NBLCOL(1,ICAR,J)='0000'
           VALLCOL(1,ICAR,J)=0.0
         END DO
C
         XO1=X(LCARO1(ICAR))
         YO1=Y(LCARO1(ICAR))
         ZO1=Z(LCARO1(ICAR))
         XO2=X(LCARO2(ICAR))
         YO2=Y(LCARO2(ICAR))
         ZO2=Z(LCARO2(ICAR))
         XO=(XO1+XO2)/2.0
         YO=(YO1+YO2)/2.0
         ZO=(ZO1+ZO2)/2.0
C
C           - FIND HIS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IHIS=1,NHIS
         IF(PKACAR(ICAR).LT.PKAHIS(IHIS))THEN
           IF(NMASS(1,ICAR)+NMASS(2,IHIS).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(2,IHIS).GT.400))THEN
           XCG=X(LHISCG(IHIS))
           YCG=Y(LHISCG(IHIS))
           ZCG=Z(LHISCG(IHIS))
           XND=X(LHISND(IHIS))
           YND=Y(LHISND(IHIS))
           ZND=Z(LHISND(IHIS))
           XCE=X(LHISCE(IHIS))
           YCE=Y(LHISCE(IHIS))
           ZCE=Z(LHISCE(IHIS))
           XNE=X(LHISNE(IHIS))
           YNE=Y(LHISNE(IHIS))
           ZNE=Z(LHISNE(IHIS))
           XCD=X(LHISCD(IHIS))
           YCD=Y(LHISCD(IHIS))
           ZCD=Z(LHISCD(IHIS))
           XCT=(XCG+XND+XCE+XNE+XCD)/5.0
           YCT=(YCG+YND+YCE+YNE+YCD)/5.0
           ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
           IF(ABS(XO-XCT).LT.DIS2.AND.
     $        ABS(YO-YCT).LT.DIS2.AND.
     $        ABS(ZO-ZCT).LT.DIS2    )  THEN
             DIS=SQRT((XO-XCT)**2+(YO-YCT)**2+(ZO-ZCT)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCAR(ICAR)=NCLCAR(ICAR)+1
               NAMCOL(1,ICAR,NCLCAR(ICAR))='HIS'
               NUMCOL(1,ICAR,NCLCAR(ICAR))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
C
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))=NAMCAR(ICAR)
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=-FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND ASP/GLU
C
         FOO=-0.80
         DIS1=2.50
         DIS2=3.50
         IF(ICAR.LT.NCAR)THEN
         DO JCAR=ICAR+1,NCAR
             XJO1=X(LCARO1(JCAR))
             YJO1=Y(LCARO1(JCAR))
             ZJO1=Z(LCARO1(JCAR))
             XJO2=X(LCARO2(JCAR))
             YJO2=Y(LCARO2(JCAR))
             ZJO2=Z(LCARO2(JCAR))
             IF((ABS(XO1-XJO1).LT.DIS2.AND.
     $         ABS(YO1-YJO1).LT.DIS2.AND.
     $         ABS(ZO1-ZJO1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XJO1).LT.DIS2.AND.
     $         ABS(YO2-YJO1).LT.DIS2.AND.
     $         ABS(ZO2-ZJO1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XJO2).LT.DIS2.AND.
     $         ABS(YO1-YJO2).LT.DIS2.AND.
     $         ABS(ZO1-ZJO2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XJO2).LT.DIS2.AND.
     $         ABS(YO2-YJO2).LT.DIS2.AND.
     $         ABS(ZO2-ZJO2).LT.DIS2    )     ) THEN
               DIS11=SQRT((XO1-XJO1)**2
     $                   +(YO1-YJO1)**2
     $                   +(ZO1-ZJO1)**2)
               DIS12=SQRT((XO1-XJO2)**2
     $                   +(YO1-YJO2)**2
     $                   +(ZO1-ZJO2)**2)
               DIS21=SQRT((XO2-XJO1)**2
     $                   +(YO2-YJO1)**2
     $                   +(ZO2-ZJO1)**2)
               DIS22=SQRT((XO2-XJO2)**2
     $                   +(YO2-YJO2)**2
     $                   +(ZO2-ZJO2)**2)
               DIS=MIN(DIS11,DIS12)
               DIS=MIN(DIS,DIS21)
               DIS=MIN(DIS,DIS22)
               IF(DIS.LT.DIS2)THEN
               IF(PKACAR(ICAR).LT.PKACAR(JCAR))THEN
                 NSDCAR(ICAR)=NSDCAR(ICAR)+1
                 NAMSDC(1,ICAR,NSDCAR(ICAR))=NAMCAR(JCAR)
                 NUMSDC(1,ICAR,NSDCAR(ICAR))=LCARRS(JCAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,ICAR,NSDCAR(ICAR))=FOO*MIN(1.0,VALUE)
C
                 NSDCAR(JCAR)=NSDCAR(JCAR)+1
                 NAMSDC(1,JCAR,NSDCAR(JCAR))=NAMCAR(ICAR)
                 NUMSDC(1,JCAR,NSDCAR(JCAR))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,JCAR,NSDCAR(JCAR))=-FOO*MIN(1.0,VALUE)
C
C                STRONG H-BONDING IF BURIED
C
                 IF(NMASS(1,ICAR)+NMASS(1,JCAR).GT.900 .OR.
     $           (NMASS(1,ICAR).GT.400.AND.NMASS(1,JCAR).GT.400))THEN
                   VALSDC(1,ICAR,NSDCAR(ICAR))=-1.60
                   VALSDC(1,JCAR,NSDCAR(JCAR))=+1.60
                 END IF
               ELSE
                 NSDCAR(ICAR)=NSDCAR(ICAR)+1
                 NAMSDC(1,ICAR,NSDCAR(ICAR))=NAMCAR(JCAR)
                 NUMSDC(1,ICAR,NSDCAR(ICAR))=LCARRS(JCAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,ICAR,NSDCAR(ICAR))=-FOO*MIN(1.0,VALUE)
C
                 NSDCAR(JCAR)=NSDCAR(JCAR)+1
                 NAMSDC(1,JCAR,NSDCAR(JCAR))=NAMCAR(ICAR)
                 NUMSDC(1,JCAR,NSDCAR(JCAR))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(1,JCAR,NSDCAR(JCAR))=FOO*MIN(1.0,VALUE)
C
C                STRONG H-BONDING IF BURIED
C
                 IF(NMASS(1,ICAR)+NMASS(1,JCAR).GT.900 .OR.
     $           (NMASS(1,ICAR).GT.400.AND.NMASS(1,JCAR).GT.400))THEN
                   VALSDC(1,ICAR,NSDCAR(ICAR))=+1.60
                   VALSDC(1,JCAR,NSDCAR(JCAR))=-1.60
                 END IF
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALSDC(1,ICAR,NSDCAR(ICAR))
               PK1CAR(JCAR)=PK1CAR(JCAR)+VALSDC(1,JCAR,NSDCAR(JCAR))
C
               END IF
             END IF
         END DO
         END IF
C
C
C           - FIND ASP/GLU (-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO JCAR=1,NCAR
           IF(JCAR.NE.ICAR .AND. PKACAR(ICAR).GT.PKACAR(JCAR))THEN
             IF(NMASS(1,ICAR)+NMASS(1,JCAR).GT.900 .OR.
     $         (NMASS(1,ICAR).GT.400.AND.NMASS(1,JCAR).GT.400))THEN
             XJO1=X(LCARO1(JCAR))
             YJO1=Y(LCARO1(JCAR))
             ZJO1=Z(LCARO1(JCAR))
             XJO2=X(LCARO2(JCAR))
             YJO2=Y(LCARO2(JCAR))
             ZJO2=Z(LCARO2(JCAR))
             XJO=(XJO1+XJO2)/2.0
             YJO=(YJO1+YJO2)/2.0
             ZJO=(ZJO1+ZJO2)/2.0
             IF(ABS(XO-XJO).LT.DIS2.AND.
     $          ABS(YO-YJO).LT.DIS2.AND.
     $          ABS(ZO-ZJO).LT.DIS2    ) THEN
               DIS=SQRT((XO-XJO)**2
     $                 +(YO-YJO)**2
     $                 +(ZO-ZJO)**2)
               IF(DIS.LT.DIS2)THEN
                 NCLCAR(ICAR)=NCLCAR(ICAR)+1
                 NAMCOL(1,ICAR,NCLCAR(ICAR))=NAMCAR(JCAR)
                 NUMCOL(1,ICAR,NCLCAR(ICAR))=LCARRS(JCAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(1,ICAR,NCLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
                 PK1CAR(ICAR)=PK1CAR(ICAR)+VALCOL(1,ICAR,NCLCAR(ICAR))
               END IF
             END IF
             END IF
           END IF
         END DO
C
C
C
C
C           - FIND CAR GROUP OF THE LIGAND (Oco)
C
         FOO=-0.80
         DIS1=2.50
         DIS2=3.50
         DO ILGCAR=1,NLIGC2
          IF(XCac(ILGCAR).NE.0.AND.YCac(ILGCAR).NE.0
     $                        .AND.ZCac(ILGCAR).NE.0)THEN
             XLGO1=XLO1(ILGCAR)
             YLGO1=YLO1(ILGCAR)
             ZLGO1=ZLO1(ILGCAR)
             XLGO2=XLO2(ILGCAR)
             YLGO2=YLO2(ILGCAR)
             ZLGO2=ZLO2(ILGCAR)
             IF((ABS(XO1-XLGO1).LT.DIS2.AND.
     $         ABS(YO1-YLGO1).LT.DIS2.AND.
     $         ABS(ZO1-ZLGO1).LT.DIS2    ) .OR.
     $        (ABS(XO2-XLGO1).LT.DIS2.AND.
     $         ABS(YO2-YLGO1).LT.DIS2.AND.
     $         ABS(ZO2-ZLGO1).LT.DIS2    ) .OR.
     $        (ABS(XO1-XLGO2).LT.DIS2.AND.
     $         ABS(YO1-YLGO2).LT.DIS2.AND.
     $         ABS(ZO1-ZLGO2).LT.DIS2    ) .OR.
     $        (ABS(XO2-XLGO2).LT.DIS2.AND.
     $         ABS(YO2-YLGO2).LT.DIS2.AND.
     $         ABS(ZO2-ZLGO2).LT.DIS2    )     ) THEN
               DIS11=SQRT((XO1-XLGO1)**2
     $                   +(YO1-YLGO1)**2
     $                   +(ZO1-ZLGO1)**2)
               DIS12=SQRT((XO1-XLGO2)**2
     $                   +(YO1-YLGO2)**2
     $                   +(ZO1-ZLGO2)**2)
               DIS21=SQRT((XO2-XLGO1)**2
     $                   +(YO2-YLGO1)**2
     $                   +(ZO2-ZLGO1)**2)
               DIS22=SQRT((XO2-XLGO2)**2
     $                   +(YO2-YLGO2)**2
     $                   +(ZO2-ZLGO2)**2)
               DIS=MIN(DIS11,DIS12)
               DIS=MIN(DIS,DIS21)
               DIS=MIN(DIS,DIS22)
               IF(DIS.LT.DIS2)THEN
C
               IF(PKACAR(ICAR).LT.PKALGCAR(ILGCAR))THEN
                 NLGCAR(ICAR)=NLGCAR(ICAR)+1
                 LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGC2(ILGCAR))
                 NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGC2(ILGCAR))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(1,ICAR,NLGCAR(ICAR))=FOO*MIN(1.0,VALUE)
C
                 NSDLCAR(ILGCAR)=NSDLCAR(ILGCAR)+1
                 NAMSDC(11,ILGCAR,NSDLCAR(ILGCAR))=NAMCAR(ICAR)
                 NUMSDC(11,ILGCAR,NSDLCAR(ILGCAR))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(11,ILGCAR,NSDLCAR(ILGCAR))=-FOO*MIN(1.0,VALUE)
C
C                STRONG H-BONDING IF BURIED
C
                 IF(NMASS(1,ICAR)+NMASS(11,ILGCAR).GT.900 .OR.
     $           (NMASS(1,ICAR).GT.400.AND.NMASS(11,ILGCAR).GT.400))THEN
                   VALLIG(1,ICAR,NLGCAR(ICAR))=-1.60
                   VALSDC(11,ILGCAR,NSDLCAR(ILGCAR))=+1.60
                 END IF
               ELSE
                 NLGCAR(ICAR)=NLGCAR(ICAR)+1
                 LABEL(1,ICAR,NLGCAR(ICAR))=NBATM(LLIGC2(ILGCAR))
                 NAMLIG(1,ICAR,NLGCAR(ICAR))=NAMRES(LLIGC2(ILGCAR))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(1,ICAR,NLGCAR(ICAR))=-FOO*MIN(1.0,VALUE)
C
                 NSDLCAR(ILGCAR)=NSDLCAR(ILGCAR)+1
                 NAMSDC(11,ILGCAR,NSDLCAR(ILGCAR))=NAMCAR(ICAR)
                 NUMSDC(11,ILGCAR,NSDLCAR(ILGCAR))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(11,ILGCAR,NSDLCAR(ILGCAR))=FOO*MIN(1.0,VALUE)
C
C                STRONG H-BONDING IF BURIED
C
                 IF(NMASS(1,ICAR)+NMASS(11,ILGCAR).GT.900 .OR.
     $           (NMASS(1,ICAR).GT.400.AND.NMASS(11,ILGCAR).GT.400))THEN
                   VALLIG(1,ICAR,NLGCAR(ICAR))=+1.60
                   VALSDC(11,ILGCAR,NSDLCAR(ILGCAR))=-1.60
                 END IF
               END IF
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLIG(1,ICAR,NLGCAR(ICAR))
               PK1LGCAR(ILGCAR)=PK1LGCAR(ILGCAR)
     $                         +VALSDC(11,ILGCAR,NSDLCAR(ILGCAR))
C
C
C
               END IF
             END IF
          ENDIF
         END DO
C
C
C
C           - FIND Oco (-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO ILGCAR=1,NLIGC2
          IF(XCac(ILGCAR).NE.0.AND.YCac(ILGCAR).NE.0
     $                        .AND.ZCac(ILGCAR).NE.0)THEN
             IF(NMASS(1,ICAR)+NMASS(11,ILGCAR).GT.900 .OR.
     $         (NMASS(1,ICAR).GT.400.AND.NMASS(11,ILGCAR).GT.400))THEN
             XLGO1=XLO1(ILGCAR)
             YLGO1=YLO1(ILGCAR)
             ZLGO1=ZLO1(ILGCAR)
             XLGO2=XLO2(ILGCAR)
             YLGO2=YLO2(ILGCAR)
             ZLGO2=ZLO2(ILGCAR)
             XLGO=(XLGO1+XLGO2)/2.0
             YLGO=(YLGO1+YLGO2)/2.0
             ZLGO=(ZLGO1+ZLGO2)/2.0
             IF(ABS(XO-XLGO).LT.DIS2.AND.
     $          ABS(YO-YLGO).LT.DIS2.AND.
     $          ABS(ZO-ZLGO).LT.DIS2    ) THEN
               DIS=SQRT((XO-XLGO)**2
     $                 +(YO-YLGO)**2
     $                 +(ZO-ZLGO)**2)
               IF(DIS.LT.DIS2)THEN
           IF(PKACAR(ICAR).GT.PKALGCAR(ILGCAR))THEN
                 NCLLCAR(ICAR)=NCLLCAR(ICAR)+1
                 NAMLCOL(1,ICAR,NCLLCAR(ICAR))=NAMRES(LLIGC2(ILGCAR))
                 NBLCOL(1,ICAR,NCLLCAR(ICAR))=NBATM(LLIGC2(ILGCAR))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLCOL(1,ICAR,NCLLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
                 PK1CAR(ICAR)=PK1CAR(ICAR)+VALLCOL(1,ICAR,NCLLCAR(ICAR))
C
           else
                 NCLLGCAR(ILGCAR)=NCLLGCAR(ILGCAR)+1
                 NAMCOL(11,ILGCAR,NCLLGCAR(ILGCAR))=NAMCAR(ICAR)
                 NUMCOL(11,ILGCAR,NCLLGCAR(ILGCAR))=LCARRS(ICAR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(11,ILGCAR,NCLLGCAR(ILGCAR))=
     $                                          FCOUL*MIN(1.0,VALUE)
                 PK1LGCAR(ILGCAR)=PK1LGCAR(ILGCAR)
     $                           +VALCOL(11,ILGCAR,NCLLGCAR(ILGCAR))
               END IF
             END IF
             END IF
           END IF
          ENDIF
         END DO
C
C
C           - FIND Nar(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO INar=1,NLIGNar
         IF(PKACAR(ICAR).LT.PKALGNar(INar))THEN
           IF(NMASS(1,ICAR)+NMASS(8,INar).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(8,INar).GT.400))THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           IF(ABS(XO-XN).LT.DIS2.AND.
     $        ABS(YO-YN).LT.DIS2.AND.
     $        ABS(ZO-ZN).LT.DIS2    )  THEN
             DIS=SQRT((XO-XN)**2+(YO-YN)**2+(ZO-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLCAR(ICAR)=NCLLCAR(ICAR)+1
               NAMLCOL(1,ICAR,NCLLCAR(ICAR))=NAMRES(LLIGNar(INar))
               NBLCOL(1,ICAR,NCLLCAR(ICAR))=NBATM(LLIGNar(INar))  
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(1,ICAR,NCLLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)
     $                     +VALLCOL(1,ICAR,NCLLCAR(ICAR))
C
               NCLNar(INar)=NCLNar(INar)+1
               NAMCOL(8,INar,NCLNar(INar))=NAMCAR(ICAR)
               NUMCOL(8,INar,NCLNar(INar))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(8,INar,NCLNar(INar))=-FCOUL*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALCOL(8,INar,NCLNar(INar))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND N3+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
C
         DO IN3=1,NLIGN3
         IF(PKACAR(ICAR).LT.PKALGN3(IN3))THEN
           IF(NMASS(1,ICAR)+NMASS(7,IN3).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(7,IN3).GT.400))THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           IF(ABS(XO-XN).LT.DIS2.AND.
     $        ABS(YO-YN).LT.DIS2.AND.
     $        ABS(ZO-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XO-XN)**2+(YO-YN)**2+(ZO-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLCAR(ICAR)=NCLLCAR(ICAR)+1
               NAMLCOL(1,ICAR,NCLLCAR(ICAR))=NAMRES(LLIGN3(IN3))
               NBLCOL(1,ICAR,NCLLCAR(ICAR))=NBATM(LLIGN3(IN3))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(1,ICAR,NCLLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLCOL(1,ICAR,NCLLCAR(ICAR))
C
               NCLN3(IN3)=NCLN3(IN3)+1
               NAMCOL(7,IN3,NCLN3(IN3))=NAMCAR(ICAR)
               NUMCOL(7,IN3,NCLN3(IN3))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(7,IN3,NCLN3(IN3))=-FCOUL*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALCOL(7,IN3,NCLN3(IN3))
             END IF
           END IF
           END IF
          END IF
         END DO
C
C
C           - FIND Npl+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IC2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
          IF(PKACAR(ICAR).LT.PKALGCg(Ic2))THEN
           IF(NMASS(1,ICAR)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(9,IC2).GT.400))THEN
C
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
            IF (ABS(XO-XC).LT.DIS2.AND.
     $          ABS(YO-YC).LT.DIS2.AND.
     $          ABS(ZO-ZC).LT.DIS2    )THEN
             DIS=SQRT((XO-XC)**2+(YO-YC)**2+(ZO-ZC)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLCAR(ICAR)=NCLLCAR(ICAR)+1
               NAMLCOL(1,ICAR,NCLLCAR(ICAR))=NAMRES(LLIGC2(IC2))
               NBLCOL(1,ICAR,NCLLCAR(ICAR))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(1,ICAR,NCLLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLCOL(1,ICAR,NCLLCAR(ICAR))
C
               NCLCg(IC2)=NCLCg(IC2)+1
               NAMCOL(9,IC2,NCLCg(IC2))=NAMCAR(ICAR)
               NUMCOL(9,IC2,NCLCg(IC2))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCg(IC2))=-FCOUL*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)
     $                        +VALCOL(9,IC2,NCLCg(IC2))
             END IF
            END IF
           END IF
          END IF
          END IF
C
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2' ) THEN
          IF(PKACAR(ICAR).LT.PKALGCg(Ic2))THEN
           IF(NMASS(1,ICAR)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(1,ICAR).GT.400.AND.NMASS(9,IC2).GT.400))THEN
C
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
            IF (ABS(XO-XC).LT.DIS2.AND.
     $          ABS(YO-YC).LT.DIS2.AND.
     $          ABS(ZO-ZC).LT.DIS2    )THEN
             DIS=SQRT((XO-XC)**2+(YO-YC)**2+(ZO-ZC)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLCAR(ICAR)=NCLLCAR(ICAR)+1
               NAMLCOL(1,ICAR,NCLLCAR(ICAR))=NAMRES(LLIGC2(IC2))
               NBLCOL(1,ICAR,NCLLCAR(ICAR))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(1,ICAR,NCLLCAR(ICAR))=FCOUL*MIN(1.0,VALUE)
               PK1CAR(ICAR)=PK1CAR(ICAR)+VALLCOL(1,ICAR,NCLLCAR(ICAR))
C
               NCLCN2(IC2)=NCLCN2(IC2)+1
               NAMCOL(9,IC2,NCLCN2(IC2))=NAMCAR(ICAR)
               NUMCOL(9,IC2,NCLCN2(IC2))=LCARRS(ICAR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCN2(IC2))=-FCOUL*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)
     $                        +VALCOL(9,IC2,NCLCN2(IC2))
             END IF
            END IF
           END IF
          END IF
          END IF
         END DO
C
C
C
       END DO
C
C
C
C      **********************
C      HIS
C      **********************
C
C
       DO IHIS=1, NHIS
         K=NSDHIS(IHIS)
         DO J=K+1,30
           NAMSDC(2,IHIS,J)='000'
           NUMSDC(2,IHIS,J)=0
           VALSDC(2,IHIS,J)=0.0
         END DO
C
         K=NBKHIS(IHIS)
         DO J=K+1,30
           NAMBKB(2,IHIS,J)='000'
           NUMBKB(2,IHIS,J)=0
           VALBKB(2,IHIS,J)=0.0
         END DO
C
         K=NCLHIS(IHIS)
         DO J=K+1,30
           NAMCOL(2,IHIS,J)='000'
           NUMCOL(2,IHIS,J)=0
           VALCOL(2,IHIS,J)=0.0
         END DO
C
         K=NLGHIS(IHIS)
         DO J=K+1,30
           NAMLIG(2,IHIS,J)='000'
           LABEL(2,IHIS,J)='0000'
           VALLIG(2,IHIS,J)=0.0
         END DO
C
         K=NCLLHIS(IHIS)
         DO J=K+1,30
           NAMLCOL(2,IHIS,J)='000'
           NBLCOL(2,IHIS,J)='0000'
           VALLCOL(2,IHIS,J)=0.0
         END DO
C
         XH1=XHISP1(IHIS)
         YH1=YHISP1(IHIS)
         ZH1=ZHISP1(IHIS)
         XH2=XHISP2(IHIS)
         YH2=YHISP2(IHIS)
         ZH2=ZHISP2(IHIS)
         XCG=X(LHISCG(IHIS))
         YCG=Y(LHISCG(IHIS))
         ZCG=Z(LHISCG(IHIS))
         XND=X(LHISND(IHIS))
         YND=Y(LHISND(IHIS))
         ZND=Z(LHISND(IHIS))
         XCE=X(LHISCE(IHIS))
         YCE=Y(LHISCE(IHIS))
         ZCE=Z(LHISCE(IHIS))
         XNE=X(LHISNE(IHIS))
         YNE=Y(LHISNE(IHIS))
         ZNE=Z(LHISNE(IHIS))
         XCD=X(LHISCD(IHIS))
         YCD=Y(LHISCD(IHIS))
         ZCD=Z(LHISCD(IHIS))
         XCT=(XCG+XND+XCE+XNE+XCD)/5.0
         YCT=(YCG+YND+YCE+YNE+YCD)/5.0
         ZCT=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
C
C
C
C           - FIND HIS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO JHIS=1,NHIS
         IF(IHIS.NE.JHIS.AND.PKAHIS(IHIS).LT.PKAHIS(JHIS))THEN
           IF(NMASS(2,IHIS)+NMASS(2,JHIS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(2,JHIS).GT.400))THEN
           XCG=X(LHISCG(JHIS))
           YCG=Y(LHISCG(JHIS))
           ZCG=Z(LHISCG(JHIS))
           XND=X(LHISND(JHIS))
           YND=Y(LHISND(JHIS))
           ZND=Z(LHISND(JHIS))
           XCE=X(LHISCE(JHIS))
           YCE=Y(LHISCE(JHIS))
           ZCE=Z(LHISCE(JHIS))
           XNE=X(LHISNE(JHIS))
           YNE=Y(LHISNE(JHIS))
           ZNE=Z(LHISNE(JHIS))
           XCD=X(LHISCD(JHIS))
           YCD=Y(LHISCD(JHIS))
           ZCD=Z(LHISCD(JHIS))
           XCTJ=(XCG+XND+XCE+XNE+XCD)/5.0
           YCTJ=(YCG+YND+YCE+YNE+YCD)/5.0
           ZCTJ=(ZCG+ZND+ZCE+ZNE+ZCD)/5.0
           IF(ABS(XCT-XCTJ).LT.DIS2.AND.
     $        ABS(YCT-YCTJ).LT.DIS2.AND.
     $        ABS(ZCT-ZCTJ).LT.DIS2    ) THEN
             DIS=SQRT((XCT-XCTJ)**2+(YCT-YCTJ)**2+(ZCT-ZCTJ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='HIS'
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LHISRS(JHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND CYS(-) COULOMBICS
C
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO ICYS=1,NCYS
         IF(TYPCYS(ICYS).NE.'BONDED'.AND.
     $                   PKAHIS(IHIS).GT.PKACYS(ICYS))THEN
           IF(NMASS(2,IHIS)+NMASS(3,ICYS).GT.900 .OR.
     $     (NMASS(2,IHIS).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
           XS=X(LCYSSG(ICYS))
           YS=Y(LCYSSG(ICYS))
           ZS=Z(LCYSSG(ICYS))
           IF((ABS(XCT-XS).LT.DIS2.AND.
     $         ABS(YCT-YS).LT.DIS2.AND.
     $         ABS(ZCT-ZS).LT.DIS2    )     ) THEN
             DIS=SQRT((XCT-XS)**2+(YCT-YS)**2+(ZCT-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NCLHIS(IHIS)=NCLHIS(IHIS)+1
               NAMCOL(2,IHIS,NCLHIS(IHIS))='CYS'
               NUMCOL(2,IHIS,NCLHIS(IHIS))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(2,IHIS,NCLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALCOL(2,IHIS,NCLHIS(IHIS))
C
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='HIS'
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=-FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C      END DO
C
C
C dcb
C moved from STEP 5
C
C           - FIND Aromatic Nitrogen (Nar) -
C
C
         FNg=-0.80
         DIS1=2.00
         DIS2=3.00
         DO ILIG=1,NLIGNar
           XN=X(LLIGNar(ILIG))
           YN=Y(LLIGNar(ILIG))
           ZN=Z(LLIGNar(ILIG))
           XP=XgLIG(ILIG)
           YP=YgLIG(ILIG)
           ZP=ZgLIG(ILIG)
           IF((ABS(XNE-XP).LT.DIS2.AND.
     $         ABS(YNE-YP).LT.DIS2.AND.
     $         ABS(ZNE-ZP).LT.DIS2    ) .OR.
     $        (ABS(XND-XP).LT.DIS2.AND.
     $         ABS(YND-YP).LT.DIS2.AND.
     $         ABS(ZND-ZP).LT.DIS2    )     ) THEN
             DISNEP=SQRT((XNE-XP)**2+(YNE-YP)**2+(ZNE-ZP)**2)
             DISNDP=SQRT((XND-XP)**2+(YND-YP)**2+(ZND-ZP)**2)
             VECNRM=SQRT((XP-XN)**2+(YP-YN)**2+(ZP-ZN)**2)
             XVNP=(XP-XN)/VECNRM
             YVNP=(YP-YN)/VECNRM
             ZVNP=(ZP-ZN)/VECNRM
             IF(DISNEP.LT.DISNDP) THEN
               DIS=DISNEP
               XVPNE=-(XP-XNE)/DISNEP
               YVPNE=-(YP-YNE)/DISNEP
               ZVPNE=-(ZP-ZNE)/DISNEP
               AGPN=XVNP*XVPNE + YVNP*YVPNE + ZVNP*ZVPNE
             ENDIF
             IF(DISNDP.LT.DISNEP) THEN
               DIS=DISNDP    
               XVPND=-(XP-XND)/DISNDP
               YVPND=-(YP-YND)/DISNDP
               ZVPND=-(ZP-ZND)/DISNDP
               AGPN=XVNP*XVPND + YVNP*YVPND + ZVNP*ZVPND
             ENDIF
             IF(DIS.LT.DIS2 .AND. AGPN.GT.0.001)THEN
C
             IF(PK1HIS(IHIS).LT.PK1LGNar(ILIG))THEN
               NLGHIS(IHIS)=NLGHIS(IHIS)+1
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLIG(2,IHIS,NLGHIS(IHIS))=FNg*MIN(1.0,VALUE)*AGPN
               LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGNar(ILIG))
               NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGNar(ILIG))
C
               NSDNar(ILIG)=NSDNar(ILIG)+1
               NAMSDC(8,ILIG,NSDNar(ILIG))='HIS'
               NUMSDC(8,ILIG,NSDNar(ILIG))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(8,ILIG,NSDNar(ILIG))=-FNg*MIN(1.0,VALUE)*AGPN
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
               IF(NMASS(1,ICAR)+NMASS(8,ILIG).GT.900 .OR.
     $     (NMASS(1,ICAR).GT.400.AND.NMASS(8,ILIG).GT.400))THEN
                 VALLIG(1,ICAR,NLGCAR(ICAR))=-1.60
                 VALSDC(8,ILIG,NSDNar(ILIG))=+1.60
               END IF
             ELSE
                NLGHIS(IHIS)=NLGHIS(IHIS)+1
                VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                VALLIG(2,IHIS,NLGHIS(IHIS))=-FNg*MIN(1.0,VALUE)*AGPN
                LABEL(2,IHIS,NLGHIS(IHIS))=NBATM(LLIGNar(ILIG))
                NAMLIG(2,IHIS,NLGHIS(IHIS))=NAMRES(LLIGNar(ILIG))
C
                NSDNar(ILIG)=NSDNar(ILIG)+1
                NAMSDC(8,ILIG,NSDNar(ILIG))='HIS'
                NUMSDC(8,ILIG,NSDNar(ILIG))=LHISRS(IHIS)
                VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                VALSDC(8,ILIG,NSDNar(ILIG))=FNg*MIN(1.0,VALUE)*AGPN
C
C              STRONG H-BONDING IF BURIED
C              (MOSTLY COULOMBIC INTERACTION)
C
              IF(NMASS(1,ICAR)+NMASS(8,ILIG).GT.900 .OR.
     $     (NMASS(1,ICAR).GT.400.AND.NMASS(8,ILIG).GT.400))THEN
                  VALLIG(1,ICAR,NLGCAR(ICAR))=+1.60
                  VALSDC(8,ILIG,NSDNar(ILIG))=-1.60
              END IF 
             END IF 
             PK1HIS(IHIS)=PK1HIS(IHIS)+VALLIG(2,IHIS,NLGHIS(IHIS))
             PK1LGNar(ILIG)=PK1LGNar(ILIG)
     $                     +VALSDC(8,ILIG,NSDNar(ILIG))
             END IF
           END IF
         END DO
C
C
C
C           - FIND Nar(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO INar=1,NLIGNar
           IF(NMASS(2,IHIS)+NMASS(8,INar).GT.900 .OR.
     $        (NMASS(2,IHIS).GT.400.AND.NMASS(8,INar).GT.400))THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           IF(ABS(XCT-XN).LT.DIS2.AND.
     $        ABS(YCT-YN).LT.DIS2.AND.
     $        ABS(ZCT-ZN).LT.DIS2    )  THEN
             DIS=SQRT((XCT-XN)**2+(YCT-YN)**2+(ZCT-ZN)**2)
             IF(DIS.LT.DIS2)THEN
         IF(PKAHIS(IHIS).LT.PKALGNar(INar))THEN
               NCLLHIS(IHIS)=NCLLHIS(IHIS)+1
               NAMLCOL(2,IHIS,NCLLHIS(IHIS))=NAMRES(LLIGNar(INar))
               NBLCOL(2,IHIS,NCLLHIS(IHIS))=NBATM(LLIGNar(INar))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(2,IHIS,NCLLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)
     $                     +VALLCOL(2,IHIS,NCLLHIS(IHIS))
C
         ELSE
               NCLNar(INar)=NCLNar(INar)+1
               NAMCOL(8,INar,NCLNar(INar))='HIS'
               NUMCOL(8,INar,NCLNar(INar))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(8,INar,NCLNar(INar))=FCOUL*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALCOL(8,INar,NCLNar(INar))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C
C           - FIND Oco(-)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILIG=1,NLIGC2
         IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                     .AND.ZCac(ILIG).NE.0)THEN
         IF(PKALGCAR(ILIG).LT.PKAHIS(IHIS))THEN
           IF(NMASS(11,ILIG)+NMASS(2,IHIS).GT.900 .OR.
     $        (NMASS(11,ILIG).GT.400.AND.NMASS(2,IHIS).GT.400))THEN
           XO1=XLO1(ILIG)
           YO1=YLO1(ILIG)
           ZO1=ZLO1(ILIG)
           XO2=XLO2(ILIG)
           YO2=YLO2(ILIG)
           ZO2=ZLO2(ILIG)
           XO=(XO1+XO2)/2.0
           YO=(YO1+YO2)/2.0
           ZO=(ZO1+ZO2)/2.0
           IF(ABS(XO-XCT).LT.DIS2.AND.
     $        ABS(YO-YCT).LT.DIS2.AND.
     $        ABS(ZO-ZCT).LT.DIS2    )  THEN
             DIS=SQRT((XO-XCT)**2+(YO-YCT)**2+(ZO-ZCT)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLHIS(IHIS)=NCLLHIS(IHIS)+1
               NAMLCOL(2,IHIS,NCLLHIS(IHIS))=NAMRES(LLIGC2(ILIG))
               NBLCOL(2,IHIS,NCLLHIS(IHIS))=NBATM(LLIGC2(ILIG))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(2,IHIS,NCLLHIS(IHIS))=-FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLCOL(2,IHIS,NCLLHIS(IHIS))
C
               NCLLGCAR(ILIG)=NCLLGCAR(ILIG)+1
               NAMCOL(11,ILIG,NCLLGCAR(ILIG))='HIS'
               NUMCOL(11,ILIG,NCLLGCAR(ILIG))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(11,ILIG,NCLLGCAR(ILIG))=FCOUL*MIN(1.0,VALUE)
               PK1LGCAR(ILIG)=PK1LGCAR(ILIG)
     $                       +VALCOL(11,ILIG,NCLLGCAR(ILIG))
C
             END IF
           END IF
           END IF
         END IF
         END IF
         END DO
C
C
C           - FIND N3+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
C
         DO IN3=1,NLIGN3
C
           IF(NMASS(2,Ihis)+NMASS(7,IN3).GT.900 .OR.
     $        (NMASS(2,Ihis).GT.400.AND.NMASS(7,IN3).GT.400))THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           IF(ABS(Xct-XN).LT.DIS2.AND.
     $        ABS(Yct-YN).LT.DIS2.AND.
     $        ABS(Zct-ZN).LT.DIS2    ) THEN
             DIS=SQRT((Xct-XN)**2+(Yct-YN)**2+(Zct-ZN)**2)
             IF(DIS.LT.DIS2)THEN
              IF(PKAhis(Ihis).LT.PKALGN3(IN3))THEN
               NCLLhis(Ihis)=NCLLhis(Ihis)+1
               NAMLCOL(2,Ihis,NCLLhis(Ihis))=NAMRES(LLIGN3(IN3))
               NBLCOL(2,Ihis,NCLLhis(Ihis))=NBATM(LLIGN3(IN3))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(2,Ihis,NCLLhis(Ihis))=FCOUL*MIN(1.0,VALUE)
               PK1his(Ihis)=PK1his(Ihis)+VALLCOL(2,Ihis,NCLLhis(Ihis))
C
              ELSE
               NCLN3(IN3)=NCLN3(IN3)+1
               NAMCOL(7,IN3,NCLN3(IN3))='HIS'
               NUMCOL(7,IN3,NCLN3(IN3))=LhisRS(Ihis)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(7,IN3,NCLN3(IN3))=FCOUL*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALCOL(7,IN3,NCLN3(IN3))
              END IF
             END IF
           END IF
           END IF
         END DO

C
C
C           - FIND Npl+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IC2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           IF(NMASS(2,IHIS)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(2,IHIS).GT.400.AND.NMASS(9,IC2).GT.400))THEN
 
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
            IF (ABS(XCT-XC).LT.DIS2.AND.
     $          ABS(YCT-YC).LT.DIS2.AND.
     $          ABS(ZCT-ZC).LT.DIS2    )THEN
             DIS=SQRT((XCT-XC)**2+(YCT-YC)**2+(ZCT-ZC)**2)
             IF(DIS.LT.DIS2)THEN
              IF(PKAhis(Ihis).LT.PKALGCg(Ic2))THEN
               NCLLHIS(IHIS)=NCLLHIS(IHIS)+1
               NAMLCOL(2,IHIS,NCLLHIS(IHIS))=NAMRES(LLIGC2(IC2))
               NBLCOL(2,IHIS,NCLLHIS(IHIS))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(2,IHIS,NCLLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLCOL(2,IHIS,NCLLHIS(IHIS))
c
              ELSE 
               NCLCg(IC2)=NCLCg(IC2)+1
               NAMCOL(9,IC2,NCLCg(IC2))='HIS'
               NUMCOL(9,IC2,NCLCg(IC2))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCg(IC2))=FCOUL*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)
     $                        +VALCOL(9,IC2,NCLCg(IC2))
              END IF
             END IF
            END IF
           END IF
          END IF
c
c
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2' ) THEN
           IF(NMASS(2,IHIS)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(2,IHIS).GT.400.AND.NMASS(9,IC2).GT.400))THEN
 
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
            IF (ABS(XCT-XC).LT.DIS2.AND.
     $          ABS(YCT-YC).LT.DIS2.AND.
     $          ABS(ZCT-ZC).LT.DIS2    )THEN
             DIS=SQRT((XCT-XC)**2+(YCT-YC)**2+(ZCT-ZC)**2)
             IF(DIS.LT.DIS2)THEN
              IF(PKAhis(Ihis).LT.PKALGCn2(Ic2))THEN
               NCLLHIS(IHIS)=NCLLHIS(IHIS)+1
               NAMLCOL(2,IHIS,NCLLHIS(IHIS))=NAMRES(LLIGC2(IC2))
               NBLCOL(2,IHIS,NCLLHIS(IHIS))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(2,IHIS,NCLLHIS(IHIS))=FCOUL*MIN(1.0,VALUE)
               PK1HIS(IHIS)=PK1HIS(IHIS)+VALLCOL(2,IHIS,NCLLHIS(IHIS))
c
              ELSE
               NCLCN2(IC2)=NCLCN2(IC2)+1
               NAMCOL(9,IC2,NCLCN2(IC2))='HIS'
               NUMCOL(9,IC2,NCLCN2(IC2))=LHISRS(IHIS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCN2(IC2))=FCOUL*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)
     $                        +VALCOL(9,IC2,NCLCN2(IC2))
              END IF
             END IF
            END IF
           END IF
          END IF
         END DO
c dmr finish
C
       END DO
C
C
C
C
C      **********************
C      CYS
C      **********************
C
C
       DO ICYS=1, NCYS
         K=NSDCYS(ICYS)
         DO J=K+1,30
           NAMSDC(3,ICYS,J)='000'
           NUMSDC(3,ICYS,J)=0
           VALSDC(3,ICYS,J)=0.0
         END DO
C
         K=NBKCYS(ICYS)
         DO J=K+1,30
           NAMBKB(3,ICYS,J)='000'
           NUMBKB(3,ICYS,J)=0
           VALBKB(3,ICYS,J)=0.0
         END DO
C
         K=NCLCYS(ICYS)
         DO J=K+1,30
           NAMCOL(3,ICYS,J)='000'
           NUMCOL(3,ICYS,J)=0
           VALCOL(3,ICYS,J)=0.0
         END DO
C
         K=NLGCYS(ICYS)
         DO J=K+1,30
           NAMLIG(3,ICYS,J)='000'
           LABEL(3,ICYS,J)='0000'
           VALLIG(3,ICYS,J)=0.0
         END DO
C
         K=NCLLCYS(ICYS)
         DO J=K+1,30
           NAMLCOL(3,ICYS,J)='000'
           NBLCOL(3,ICYS,J)='0000'
           VALLCOL(3,ICYS,J)=0.0
         END DO
C
       IF(TYPCYS(ICYS).NE.'BONDED')THEN
         XSG=X(LCYSSG(ICYS))
         YSG=Y(LCYSSG(ICYS))
         ZSG=Z(LCYSSG(ICYS))
C
C
C           - FIND CYS(-)
C
         FCOUL=+2.40
         DIS1=4.00
         DIS2=7.00
         DO JCYS=1,NCYS
         IF(TYPCYS(JCYS).NE.'BONDED'.AND.
     $            PKACYS(ICYS).GT.PKACYS(JCYS))THEN
           IF(NMASS(3,ICYS)+NMASS(3,JCYS).GT.900 .OR.
     $    (NMASS(3,ICYS).GT.400.AND.NMASS(3,JCYS).GT.400))THEN
             XS=X(LCYSSG(JCYS))
             YS=Y(LCYSSG(JCYS))
             ZS=Z(LCYSSG(JCYS))
             DIS=SQRT((XSG-XS)**2+(YSG-YS)**2+(ZSG-ZS)**2)
             IF(DIS.LT.DIS2)THEN
               NCLCYS(ICYS)=NCLCYS(ICYS)+1
               NAMCOL(3,ICYS,NCLCYS(ICYS))='CYS'
               NUMCOL(3,ICYS,NCLCYS(ICYS))=LCYSRS(JCYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(3,ICYS,NCLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALCOL(3,ICYS,NCLCYS(ICYS))
             END IF
           END IF
         END IF
         END DO
C
C
C           - FIND Nar(+) 
C
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO INar=1,NLIGNar
         IF( PKALGNar(INar).GT.PKACYS(ICYS))THEN
           IF(NMASS(8,INar)+NMASS(3,ICYS).GT.900 .OR.
     $     (NMASS(8,INar).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           IF((ABS(XSG-XN).LT.DIS2.AND.
     $         ABS(YSG-YN).LT.DIS2.AND.
     $         ABS(ZSG-ZN).LT.DIS2    )     ) THEN
             DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLCYS(ICYS)=NCLLCYS(ICYS)+1
               NAMLCOL(3,ICYS,NCLLCYS(ICYS))=NAMRES(LLIGNar(INar))
               NBLCOL(3,ICYS,NCLLCYS(ICYS))=NBATM(LLIGNar(INar))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(3,ICYS,NCLLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLCOL(3,ICYS,NCLLCYS(ICYS))
C
               NCLNar(INar)=NCLNar(INar)+1
               NAMCOL(8,INar,NCLNar(INar))='CYS'
               NUMCOL(8,INar,NCLNar(INar))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(8,INar,NCLNar(INar))=-FCOUL*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALCOL(8,INar,NCLNar(INar))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C
Cdmr st
C           - FIND Oco(-)
C
         FCOUL=+2.40
         DIS1=4.00
         DIS2=7.00
         DO ILIG=1,NLIGC2
          IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                     .AND.ZCac(ILIG).NE.0)THEN
          IF(NMASS(11,ILIG)+NMASS(3,ICYS).GT.900 .OR.
     $      (NMASS(11,ILIG).GT.400.AND.NMASS(3,ICYS).GT.400))THEN
            XO1=XLO1(ILIG)
            YO1=YLO1(ILIG)
            ZO1=ZLO1(ILIG)
            XO2=XLO2(ILIG)
            YO2=YLO2(ILIG)
            ZO2=ZLO2(ILIG)
            XO=(XO1+XO2)/2.0
            YO=(YO1+YO2)/2.0
            ZO=(ZO1+ZO2)/2.0
           IF((ABS(XO-XSG).LT.DIS2.AND.
     $         ABS(YO-YSG).LT.DIS2.AND.
     $         ABS(ZO-ZSG).LT.DIS2    )     ) THEN
             DIS=SQRT((XO-XSG)**2+(YO-YSG)**2+(ZO-ZSG)**2)
             IF(DIS.LT.DIS2)THEN
                   ilgcar=ilig
           IF(PKAcys(Icys).GT.PKALGCAR(ILGCAR))THEN
               NCLLCYS(ICYS)=NCLLCYS(ICYS)+1
               NAMLCOL(3,ICYS,NCLLCYS(ICYS))=NAMRES(LLIGC2(ILIG))
               NBLCOL(3,ICYS,NCLLCYS(ICYS))=NBATM(LLIGC2(ILIG))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(3,ICYS,NCLLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLCOL(3,ICYS,NCLLCYS(ICYS))
           else
                 NCLLGCAR(ILGCAR)=NCLLGCAR(ILGCAR)+1
                 NAMCOL(11,ILGCAR,NCLLGCAR(ILGCAR))='CYS'
                 NUMCOL(11,ILGCAR,NCLLGCAR(ILGCAR))=LcysRS(Icys)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(11,ILGCAR,NCLLGCAR(ILGCAR))=
     $                                          FCOUL*MIN(1.0,VALUE)
                 PK1LGCAR(ILGCAR)=PK1LGCAR(ILGCAR)
     $                           +VALCOL(11,ILGCAR,NCLLGCAR(ILGCAR))
           end if
             END IF
            END IF
          END IF
          ENDIF
         END DO
C
C
C           - FIND N3+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
C
         DO IN3=1,NLIGN3
         IF( PKALGN3(IN3).GT.PKACYS(ICYS))THEN
           IF(NMASS(3,ICYS)+NMASS(7,IN3).GT.900 .OR.
     $        (NMASS(3,ICYS).GT.400.AND.NMASS(7,IN3).GT.400))THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           IF(ABS(XSG-XN).LT.DIS2.AND.
     $        ABS(YSG-YN).LT.DIS2.AND.
     $        ABS(ZSG-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XSG-XN)**2+(YSG-YN)**2+(ZSG-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLCYS(ICYS)=NCLLCYS(ICYS)+1
               NAMLCOL(3,ICYS,NCLLCYS(ICYS))=NAMRES(LLIGN3(IN3))
               NBLCOL(3,ICYS,NCLLCYS(ICYS))=NBATM(LLIGN3(IN3))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(3,ICYS,NCLLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLCOL(3,ICYS,NCLLCYS(ICYS))
C
               NCLN3(IN3)=NCLN3(IN3)+1
               NAMCOL(7,IN3,NCLN3(IN3))='CYS'
               NUMCOL(7,IN3,NCLN3(IN3))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(7,IN3,NCLN3(IN3))=-FCOUL*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALCOL(7,IN3,NCLN3(IN3))
             END IF
           END IF
           END IF
          END IF
         END DO
C
C
C           - FIND Npl+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IC2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
          IF( PKALGCg(Ic2).GT.PKACYS(ICYS))THEN
           IF(NMASS(3,ICYS)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(3,ICYS).GT.400.AND.NMASS(9,IC2).GT.400))THEN
C
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
            IF (ABS(XSG-XC).LT.DIS2.AND.
     $          ABS(YSG-YC).LT.DIS2.AND.
     $          ABS(ZSG-ZC).LT.DIS2    )THEN
              DIS=SQRT((XSG-XC)**2+(YSG-YC)**2+(ZSG-ZC)**2)
              IF(DIS.LT.DIS2)THEN
               NCLLCYS(ICYS)=NCLLCYS(ICYS)+1
               NAMLCOL(3,ICYS,NCLLCYS(ICYS))=NAMRES(LLIGC2(IC2))
               NBLCOL(3,ICYS,NCLLCYS(ICYS))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(3,ICYS,NCLLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLCOL(3,ICYS,NCLLCYS(ICYS))
C
               NCLCg(IC2)=NCLCg(IC2)+1
               NAMCOL(9,IC2,NCLCg(IC2))='CYS'
               NUMCOL(9,IC2,NCLCg(IC2))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCg(IC2))=-FCOUL*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)
     $                        +VALCOL(9,IC2,NCLCg(IC2))
              END IF
            END IF
           END IF
          END IF
          END IF
C
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2' ) THEN
          IF( PKALGCg(Ic2).GT.PKACYS(ICYS))THEN
           IF(NMASS(3,ICYS)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(3,ICYS).GT.400.AND.NMASS(9,IC2).GT.400))THEN
C
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
            IF (ABS(XSG-XC).LT.DIS2.AND.
     $          ABS(YSG-YC).LT.DIS2.AND.
     $          ABS(ZSG-ZC).LT.DIS2    )THEN
              DIS=SQRT((XSG-XC)**2+(YSG-YC)**2+(ZSG-ZC)**2)
              IF(DIS.LT.DIS2)THEN
               NCLLCYS(ICYS)=NCLLCYS(ICYS)+1
               NAMLCOL(3,ICYS,NCLLCYS(ICYS))=NAMRES(LLIGC2(IC2))
               NBLCOL(3,ICYS,NCLLCYS(ICYS))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(3,ICYS,NCLLCYS(ICYS))=FCOUL*MIN(1.0,VALUE)
               PK1CYS(ICYS)=PK1CYS(ICYS)+VALLCOL(3,ICYS,NCLLCYS(ICYS))
C
               NCLCN2(IC2)=NCLCN2(IC2)+1
               NAMCOL(9,IC2,NCLCN2(IC2))='CYS'
               NUMCOL(9,IC2,NCLCN2(IC2))=LCYSRS(ICYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCN2(IC2))=-FCOUL*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)
     $                        +VALCOL(9,IC2,NCLCN2(IC2))
              END IF
            END IF
           END IF
          END IF
          END IF
         END DO
C
C
Cdmr fin
C
C
C
       END IF
       END DO
C
C
C
C      **********************
C      TYR
C      **********************
C
C
       DO ITYR=1, NTYR
         K=NSDTYR(ITYR)
         DO J=K+1,30
           NAMSDC(4,ITYR,J)='000'
           NUMSDC(4,ITYR,J)=0
           VALSDC(4,ITYR,J)=0.0
         END DO
C
         K=NBKTYR(ITYR)
         DO J=K+1,30
           NAMBKB(4,ITYR,J)='000'
           NUMBKB(4,ITYR,J)=0
           VALBKB(4,ITYR,J)=0.0
         END DO
C
         K=NCLTYR(ITYR)
         DO J=K+1,30
           NAMCOL(4,ITYR,J)='000'
           NUMCOL(4,ITYR,J)=0
           VALCOL(4,ITYR,J)=0.0
         END DO
C
         K=NLGTYR(ITYR)
         DO J=K+1,30
           NAMLIG(4,ITYR,J)='000'
           LABEL(4,ITYR,J)='0000'
           VALLIG(4,ITYR,J)=0.0
         END DO
C
         K=NCLLTYR(ITYR)
         DO J=K+1,30
           NAMLCOL(4,ITYR,J)='000'
           NBLCOL(4,ITYR,J)='0000'
           VALLCOL(4,ITYR,J)=0.0
         END DO
C
         XOH=X(LTYROH(ITYR))
         YOH=Y(LTYROH(ITYR))
         ZOH=Z(LTYROH(ITYR))
C
C           - FIND TYR-OH
C
         FOO=0.80
         DIS1=3.50
         DIS2=4.50
         IF(ITYR.LT.NTYR)THEN
         DO JTYR=ITYR+1,NTYR
           XO=X(LTYROH(JTYR))
           YO=Y(LTYROH(JTYR))
           ZO=Z(LTYROH(JTYR))
           IF(ABS(XOH-XO).LT.DIS2.AND.
     $        ABS(YOH-YO).LT.DIS2.AND.
     $        ABS(ZOH-ZO).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XO)**2+(YOH-YO)**2+(ZOH-ZO)**2)
             IF(DIS.LT.DIS2)THEN
             IF(PKATYR(ITYR).GT.PKATYR(JTYR))THEN
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='TYR'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LTYRRS(JTYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=FOO*MIN(1.0,VALUE)
C
               NSDTYR(JTYR)=NSDTYR(JTYR)+1
               NAMSDC(4,JTYR,NSDTYR(JTYR))='TYR'
               NUMSDC(4,JTYR,NSDTYR(JTYR))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,JTYR,NSDTYR(JTYR))=-FOO*MIN(1.0,VALUE)
             ELSE
               NSDTYR(ITYR)=NSDTYR(ITYR)+1
               NAMSDC(4,ITYR,NSDTYR(ITYR))='TYR'
               NUMSDC(4,ITYR,NSDTYR(ITYR))=LTYRRS(JTYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,ITYR,NSDTYR(ITYR))=-FOO*MIN(1.0,VALUE)
C
               NSDTYR(JTYR)=NSDTYR(JTYR)+1
               NAMSDC(4,JTYR,NSDTYR(JTYR))='TYR'
               NUMSDC(4,JTYR,NSDTYR(JTYR))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALSDC(4,JTYR,NSDTYR(JTYR))=FOO*MIN(1.0,VALUE)
             END IF
             PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
             PK1TYR(JTYR)=PK1TYR(JTYR)+VALSDC(4,JTYR,NSDTYR(JTYR))
             END IF
           END IF
         END DO
         END IF
C
C           - FIND TYR(-)
C
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO JTYR=1,NTYR
         IF(ITYR.NE.JTYR .AND. PKATYR(ITYR).GT.PKATYR(JTYR))THEN
           IF(NMASS(4,ITYR)+NMASS(4,JTYR).GT.900 .OR.
     $       (NMASS(4,ITYR).GT.400.AND.NMASS(4,JTYR).GT.400))THEN
           XOJ=X(LTYROH(JTYR))
           YOJ=Y(LTYROH(JTYR))
           ZOJ=Z(LTYROH(JTYR))
           IF(ABS(XOH-XOJ).LT.DIS2.AND.
     $        ABS(YOH-YOJ).LT.DIS2.AND.
     $        ABS(ZOH-ZOJ).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XOJ)**2+(YOH-YOJ)**2+(ZOH-ZOJ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLTYR(ITYR)=NCLTYR(ITYR)+1
               NAMCOL(4,ITYR,NCLTYR(ITYR))='TYR'
               NUMCOL(4,ITYR,NCLTYR(ITYR))=LTYRRS(JTYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND LYS H-BONDING
C
         FOH=-0.80
         DIS1=3.00
         DIS2=4.00
         DO ILYS=1,NLYS
           FOH=-0.80
           DIS2=4.00
           IF(ILYS.EQ.1)FOH=-1.20
           IF(ILYS.EQ.1)DIS2=4.50
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XOH-XN).LT.DIS2.AND.
     $        ABS(YOH-YN).LT.DIS2.AND.
     $        ABS(ZOH-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XN)**2+(YOH-YN)**2+(ZOH-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               IF(PKATYR(ITYR).LT.PKALYS(ILYS))THEN
                 NSDTYR(ITYR)=NSDTYR(ITYR)+1
                 NAMSDC(4,ITYR,NSDTYR(ITYR))='LYS'
                 IF(ILYS.EQ.1)NAMSDC(4,ITYR,NSDTYR(ITYR))='N+ '
                 NUMSDC(4,ITYR,NSDTYR(ITYR))=LLYSRS(ILYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(4,ITYR,NSDTYR(ITYR))=FOH*MIN(1.0,VALUE)
                 PK1TYR(ITYR)=PK1TYR(ITYR)+VALSDC(4,ITYR,NSDTYR(ITYR))
               END IF
             END IF
           END IF
         END DO
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILYS=1,NLYS
           IF(NMASS(4,ITYR)+NMASS(5,ILYS).GT.900 .OR.
     $        (NMASS(4,ITYR).GT.400.AND.NMASS(5,ILYS).GT.400))THEN
           XN=X(LLYSNZ(ILYS))
           YN=Y(LLYSNZ(ILYS))
           ZN=Z(LLYSNZ(ILYS))
           IF(ABS(XOH-XN).LT.DIS2.AND.
     $        ABS(YOH-YN).LT.DIS2.AND.
     $        ABS(ZOH-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XN)**2+(YOH-YN)**2+(ZOH-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               IF(PKATYR(ITYR).LT.PKALYS(ILYS))THEN
                 NCLTYR(ITYR)=NCLTYR(ITYR)+1
                 NAMCOL(4,ITYR,NCLTYR(ITYR))='LYS'
                 IF(ILYS.EQ.1)NAMCOL(4,ITYR,NCLTYR(ITYR))='N+ '
                 NUMCOL(4,ITYR,NCLTYR(ITYR))=LLYSRS(ILYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(4,ITYR,NCLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
                 PK1TYR(ITYR)=PK1TYR(ITYR)+VALCOL(4,ITYR,NCLTYR(ITYR))
C
                 NCLLYS(ILYS)=NCLLYS(ILYS)+1
                 NAMCOL(5,ILYS,NCLLYS(ILYS))='TYR'
                 NUMCOL(5,ILYS,NCLLYS(ILYS))=LTYRRS(ITYR)
                 VALCOL(5,ILYS,NCLLYS(ILYS))=-FCOUL*MIN(1.0,VALUE)
                 PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
               END IF
             END IF
           END IF
           END IF
         END DO
C
C
C
C
C           - FIND N3 atoms -
C
         FNH=-0.80
         DIS1=2.0
         DIS2=3.5
         DO IN3=1,NLIGN3
C
C
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N31'.OR.
     $        NAMATM(LLIGN3(IN3)).EQ.'  N30') THEN
         FNH=-0.80
         DIS1=3.0
         DIS2=4.0
             XN=X(LLIGN3(IN3))
             YN=Y(LLIGN3(IN3))
             ZN=Z(LLIGN3(IN3))
             IF((ABS(XOH-XN).LT.DIS2.AND.
     $           ABS(YOH-YN).LT.DIS2.AND.
     $           ABS(ZOH-ZN).LT.DIS2    )     ) THEN
               DIS=SQRT((XOH-XN)**2+(YOH-YN)**2+(ZOH-ZN)**2)
               IF(DIS.LT.DIS2)THEN
               IF(PKATYR(ITYR).LT.PKALGN3(IN3))THEN
                 NLGTYR(ITYR)=NLGTYR(ITYR)+1
                 NAMLIG(4,ITYR,NLGTYR(ITYR))=NAMRES(LLIGN3(IN3))
                 LABEL(4,ITYR,NLGTYR(ITYR))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(4,ITYR,NLGTYR(ITYR))=FNH*MIN(1.0,VALUE)
                 PK1TYR(ITYR)=PK1TYR(ITYR)+VALLIG(4,ITYR,NLGTYR(ITYR))
C
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))='TYR'
                 NUMSDC(7,IN3,NSDN3(IN3))=LTYRRS(ITYR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
               END IF
             END IF
          ENDIF
c
c dmr start
C
C dcb
C angle contribution added in july 2007
C dcb
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N32') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XNH1=XN3P1(IN3)
           YNH1=YN3P1(IN3)
           ZNH1=ZN3P1(IN3)
           XNH2=XN3P2(IN3)
           YNH2=YN3P2(IN3)
           ZNH2=ZN3P2(IN3)
           IF((ABS(XnH1-Xoh).LT.DIS2.AND.
     $         ABS(YnH1-Yoh).LT.DIS2.AND.
     $         ABS(ZnH1-Zoh).LT.DIS2    ) .OR.
     $        (ABS(XnH2-Xoh).LT.DIS2.AND.
     $         ABS(YnH2-Yoh).LT.DIS2.AND.
     $         ABS(ZnH2-Zoh).LT.DIS2    )     ) THEN
             DISH1=SQRT((XnH1-Xoh)**2+(YnH1-Yoh)**2+(ZnH1-Zoh)**2)
             DISH2=SQRT((XnH2-Xoh)**2+(YnH2-Yoh)**2+(ZnH2-Zoh)**2)
             IF(DISH1.LT.DISH2) THEN
               DIS=DISH1
               VECNM=SQRT((XNH1-XN)**2+(YNH1-YN)**2+(ZNH1-ZN)**2)
               XNP=(XNH1-XN)/VECNM
               YNP=(YNH1-YN)/VECNM
               ZNP=(ZNH1-ZN)/VECNM
C
               XVPO=-(XOH-XNH1)/DIS
               YVPO=-(YOH-YNH1)/DIS
               ZVPO=-(ZOH-ZNH1)/DIS
               AGOH=XNP*XVPO+YNP*YVPO+ZNP*ZVPO
             ENDIF
             IF(DISH2.LT.DISH1) THEN
               DIS=DISH2
               VECNM=SQRT((XNH2-XN)**2+(YNH2-YN)**2+(ZNH2-ZN)**2)
               XNP=(XNH2-XN)/VECNM
               YNP=(YNH2-YN)/VECNM
               ZNP=(ZNH2-ZN)/VECNM
C
               XVPO=-(XOH-XNH2)/DIS
               YVPO=-(YOH-YNH2)/DIS
               ZVPO=-(ZOH-ZNH2)/DIS
               AGOH=XNP*XVPO+YNP*YVPO+ZNP*ZVPO
             ENDIF
               IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               IF(PKAtyr(Ityr).LT.PKALGN3(IN3))THEN
                 NLGtyr(Ityr)=NLGtyr(Ityr)+1
                 NAMLIG(4,Ityr,NLGtyr(Ityr))=NAMRES(LLIGN3(IN3))
                 LABEL(4,Ityr,NLGtyr(Ityr))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(4,Ityr,NLGtyr(Ityr))=FNH*MIN(1.0,VALUE)*AGOH
                 PK1tyr(Ityr)=PK1tyr(Ityr)+VALLIG(4,Ityr,NLGtyr(Ityr))
 
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))='TYR'
                 NUMSDC(7,IN3,NSDN3(IN3))=LTYRRS(ITYR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)*AGOH
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
               END IF
           END IF
          ENDIF

C
C dcb
C angle contribution added in july 2007
C dcb
          IF (NAMATM(LLIGN3(IN3)).EQ.'  N33') THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           XNH1=XN3P(IN3)
           YNH1=YN3P(IN3)
           ZNH1=ZN3P(IN3)
           IF( ABS(XnH1-Xoh).LT.DIS2.AND.
     $         ABS(YnH1-Yoh).LT.DIS2.AND.
     $         ABS(ZnH1-Zoh).LT.DIS2    ) then
             DIS=SQRT((XnH1-Xoh)**2+(YnH1-Yoh)**2+(ZnH1-Zoh)**2)
             VECNM=SQRT((XNH1-XN)**2+(YNH1-YN)**2+(ZNH1-ZN)**2)
             XNP=(XNH1-XN)/VECNM
             YNP=(YNH1-YN)/VECNM
             ZNP=(ZNH1-ZN)/VECNM
C
             XVPO=-(XOH-XNH1)/DIS
             YVPO=-(YOH-YNH1)/DIS
             ZVPO=-(ZOH-ZNH1)/DIS
             AGOH=XNP*XVPO+YNP*YVPO+ZNP*ZVPO
               IF(DIS.LT.DIS2 .AND. AGOH.GT.0.001)THEN
               IF(PKAtyr(Ityr).LT.PKALGN3(IN3))THEN
                 NLGtyr(Ityr)=NLGtyr(Ityr)+1
                 NAMLIG(2,Ityr,NLGtyr(Ityr))=NAMRES(LLIGN3(IN3))
                 LABEL(2,Ityr,NLGtyr(Ityr))=NBATM(LLIGN3(IN3))
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALLIG(2,Ityr,NLGtyr(Ityr))=FNH*MIN(1.0,VALUE)*AGOH
                 PK1tyr(Ityr)=PK1tyr(Ityr)+VALLIG(2,Ityr,NLGtyr(Ityr))
 
                 NSDN3(IN3)=NSDN3(IN3)+1
                 NAMSDC(7,IN3,NSDN3(IN3))='TYR'
                 NUMSDC(7,IN3,NSDN3(IN3))=LTYRRS(ITYR)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALSDC(7,IN3,NSDN3(IN3))=-FNH*MIN(1.0,VALUE)*AGOH
                 PK1LGN3(IN3)=PK1LGN3(IN3)+VALSDC(7,IN3,NSDN3(IN3))
               END IF
               END IF
             END IF
          ENDIF
         END DO
C
C           - FIND N3+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
C
         DO IN3=1,NLIGN3
          IF(PKATYR(ITYR).LT.PKALGN3(IN3))THEN
           IF(NMASS(4,ITYR)+NMASS(7,IN3).GT.900 .OR.
     $        (NMASS(4,ITYR).GT.400.AND.NMASS(7,IN3).GT.400))THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           IF(ABS(XOH-XN).LT.DIS2.AND.
     $        ABS(YOH-YN).LT.DIS2.AND.
     $        ABS(ZOH-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XOH-XN)**2+(YOH-YN)**2+(ZOH-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLTYR(ITYR)=NCLLTYR(ITYR)+1
               NAMLCOL(4,ITYR,NCLLTYR(ITYR))=NAMRES(LLIGN3(IN3))
               NBLCOL(4,ITYR,NCLLTYR(ITYR))=NBATM(LLIGN3(IN3))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(4,ITYR,NCLLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLCOL(4,ITYR,NCLLTYR(ITYR))
C
               NCLN3(IN3)=NCLN3(IN3)+1
               NAMCOL(7,IN3,NCLN3(IN3))='TYR'
               NUMCOL(7,IN3,NCLN3(IN3))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(7,IN3,NCLN3(IN3))=-FCOUL*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALCOL(7,IN3,NCLN3(IN3))
              END IF
             END IF
           END IF
           END IF
         END DO
C
C
C           - FIND Oco(-)
C
         FCOUL=+2.40
         DIS1=4.00
         DIS2=7.00
         DO ILIG=1,NLIGC2
          IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                      .AND.ZCac(ILIG).NE.0)THEN
          IF(NMASS(11,ILIG)+NMASS(4,ITYR).GT.900 .OR.
     $      (NMASS(11,ILIG).GT.400.AND.NMASS(4,ITYR).GT.400))THEN
            XO1=XLO1(ILIG)
            YO1=YLO1(ILIG)
            ZO1=ZLO1(ILIG)
            XO2=XLO2(ILIG)
            YO2=YLO2(ILIG)
            ZO2=ZLO2(ILIG)
            XO=(XO1+XO2)/2.0
            YO=(YO1+YO2)/2.0
            ZO=(ZO1+ZO2)/2.0
           IF((ABS(XO-XOH).LT.DIS2.AND.
     $         ABS(YO-YOH).LT.DIS2.AND.
     $         ABS(ZO-ZOH).LT.DIS2    )     ) THEN
             DIS=SQRT((XO-XOH)**2+(YO-YOH)**2+(ZO-ZOH)**2)
             IF(DIS.LT.DIS2)THEN
                   ilgcar=ilig
           IF(PKAtyr(Ityr).GT.PKALGCAR(ILGCAR))THEN
               NCLLTYR(ITYR)=NCLLTYR(ITYR)+1
               NAMLCOL(4,ITYR,NCLLTYR(ITYR))=NAMRES(LLIGC2(ILIG))
               NBLCOL(4,ITYR,NCLLTYR(ITYR))=NBATM(LLIGC2(ILIG))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(4,ITYR,NCLLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLCOL(4,ITYR,NCLLTYR(ITYR))
           else
                 NCLLGCAR(ILGCAR)=NCLLGCAR(ILGCAR)+1
                 NAMCOL(11,ILGCAR,NCLLGCAR(ILGCAR))='TYR'
                 NUMCOL(11,ILGCAR,NCLLGCAR(ILGCAR))=LtyrRS(Ityr)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(11,ILGCAR,NCLLGCAR(ILGCAR))=
     $                                          FCOUL*MIN(1.0,VALUE)
                 PK1LGCAR(ILGCAR)=PK1LGCAR(ILGCAR)
     $                           +VALCOL(11,ILGCAR,NCLLGCAR(ILGCAR))
           end if
             END IF
            END IF
          END IF
          ENDIF
         END DO
C
C
C           - FIND Nar(+) 
C
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO INar=1,NLIGNar
         IF( PKALGNar(INar).GT.PKAtyr(Ityr))THEN
           IF(NMASS(8,INar)+NMASS(4,Ityr).GT.900 .OR.
     $     (NMASS(8,INar).GT.400.AND.NMASS(4,Ityr).GT.400))THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           IF((ABS(Xoh-XN).LT.DIS2.AND.
     $         ABS(Yoh-YN).LT.DIS2.AND.
     $         ABS(Zoh-ZN).LT.DIS2    )     ) THEN
             DIS=SQRT((Xoh-XN)**2+(Yoh-YN)**2+(Zoh-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLtyr(Ityr)=NCLLtyr(Ityr)+1
               NAMLCOL(3,Ityr,NCLLtyr(Ityr))=NAMRES(LLIGNar(INar))
               NBLCOL(3,Ityr,NCLLtyr(Ityr))=NBATM(LLIGNar(INar))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(3,Ityr,NCLLtyr(Ityr))=FCOUL*MIN(1.0,VALUE)
               PK1tyr(Ityr)=PK1tyr(Ityr)+VALLCOL(3,Ityr,NCLLtyr(Ityr))
C
               NCLNar(INar)=NCLNar(INar)+1
               NAMCOL(8,INar,NCLNar(INar))='TYR'
               NUMCOL(8,INar,NCLNar(INar))=LtyrRS(Ityr)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(8,INar,NCLNar(INar))=-FCOUL*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALCOL(8,INar,NCLNar(INar))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C
C           - FIND Npl+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IC2=1,NLIGC2
         IF( PKALGCg(Ic2).GT.PKAtyr(Ityr))THEN
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           IF(NMASS(4,ITYR)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(4,ITYR).GT.400.AND.NMASS(9,IC2).GT.400))THEN
C
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
           IF (ABS(XOH-XC).LT.DIS2.AND.
     $         ABS(YOH-YC).LT.DIS2.AND.
     $         ABS(ZOH-ZC).LT.DIS2    )THEN
             DIS=SQRT((XOH-XC)**2+(YOH-YC)**2+(ZOH-ZC)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLTYR(ITYR)=NCLLTYR(ITYR)+1
               NAMLCOL(4,ITYR,NCLLTYR(ITYR))=NAMRES(LLIGC2(IC2))
               NBLCOL(4,ITYR,NCLLTYR(ITYR))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(4,ITYR,NCLLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLCOL(4,ITYR,NCLLTYR(ITYR))
C
               NCLCg(IC2)=NCLCg(IC2)+1
               NAMCOL(9,IC2,NCLCg(IC2))='TYR'
               NUMCOL(9,IC2,NCLCg(IC2))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCg(IC2))=-FCOUL*MIN(1.0,VALUE)
               PK1LGCg(IC2)=PK1LGCg(IC2)
     $                        +VALCOL(9,IC2,NCLCg(IC2))
             END IF
           END IF
           END IF
          END IF
          END IF
C
C
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2' ) THEN
          IF( PKALGCn2(Ic2).GT.PKAtyr(Ityr))THEN
           IF(NMASS(4,ITYR)+NMASS(9,IC2).GT.900 .OR.
     $        (NMASS(4,ITYR).GT.400.AND.NMASS(9,IC2).GT.400))THEN
C
            XC=X(LLIGC2(IC2))
            YC=Y(LLIGC2(IC2))
            ZC=Z(LLIGC2(IC2))
           IF (ABS(XOH-XC).LT.DIS2.AND.
     $         ABS(YOH-YC).LT.DIS2.AND.
     $         ABS(ZOH-ZC).LT.DIS2    )THEN
             DIS=SQRT((XOH-XC)**2+(YOH-YC)**2+(ZOH-ZC)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLTYR(ITYR)=NCLLTYR(ITYR)+1
               NAMLCOL(4,ITYR,NCLLTYR(ITYR))=NAMRES(LLIGC2(IC2))
               NBLCOL(4,ITYR,NCLLTYR(ITYR))=NBATM(LLIGC2(IC2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(4,ITYR,NCLLTYR(ITYR))=FCOUL*MIN(1.0,VALUE)
               PK1TYR(ITYR)=PK1TYR(ITYR)+VALLCOL(4,ITYR,NCLLTYR(ITYR))
C
               NCLCN2(IC2)=NCLCN2(IC2)+1
               NAMCOL(9,IC2,NCLCN2(IC2))='TYR'
               NUMCOL(9,IC2,NCLCN2(IC2))=LTYRRS(ITYR)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,IC2,NCLCN2(IC2))=-FCOUL*MIN(1.0,VALUE)
               PK1LGCN2(IC2)=PK1LGCN2(IC2)
     $                        +VALCOL(9,IC2,NCLCN2(IC2))
             END IF
           END IF
           END IF
          END IF
          END IF
         END DO
C
Cdmr fin
C
C
       END DO
C
C
C
C
C
C
C
C      **********************
C      LYS
C      **********************
C
C
       DO ILYS=1, NLYS
         K=NSDLYS(ILYS)
         DO J=K+1,30
           NAMSDC(5,ILYS,J)='000'
           NUMSDC(5,ILYS,J)=0
           VALSDC(5,ILYS,J)=0.0
         END DO
C
         K=NBKLYS(ILYS)
         DO J=K+1,30
           NAMBKB(5,ILYS,J)='000'
           NUMBKB(5,ILYS,J)=0
           VALBKB(5,ILYS,J)=0.0
         END DO
C
         K=NCLLYS(ILYS)
         DO J=K+1,30
           NAMCOL(5,ILYS,J)='000'
           NUMCOL(5,ILYS,J)=0
           VALCOL(5,ILYS,J)=0.0
         END DO
C
c dmr str
         K=NCLLlys(Ilys)
         DO J=K+1,30
           NAMLCOL(5,Ilys,J)='000'
           NBLCOL(5,Ilys,J)='0000'
           VALLCOL(5,Ilys,J)=0.0
         END DO
c dmr fin
C
         XNZ=X(LLYSNZ(ILYS))
         YNZ=Y(LLYSNZ(ILYS))
         ZNZ=Z(LLYSNZ(ILYS))
C
C
C           - FIND LYS(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO JLYS=1,NLYS
         IF(ILYS.NE.JLYS .AND. PKALYS(ILYS).LT.PKALYS(JLYS))THEN
           IF(NMASS(5,ILYS)+NMASS(5,JLYS).GT.900 .OR.
     $     (NMASS(5,ILYS).GT.400.AND.NMASS(5,JLYS).GT.400))THEN
           XN=X(LLYSNZ(JLYS))
           YN=Y(LLYSNZ(JLYS))
           ZN=Z(LLYSNZ(JLYS))
           IF(ABS(XNZ-XN).LT.DIS2.AND.
     $        ABS(YNZ-YN).LT.DIS2.AND.
     $        ABS(ZNZ-ZN).LT.DIS2    ) THEN
             DIS=SQRT((XNZ-XN)**2+(YNZ-YN)**2+(ZNZ-ZN)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLYS(ILYS)=NCLLYS(ILYS)+1
               NAMCOL(5,ILYS,NCLLYS(ILYS))='LYS'
               IF(JLYS.EQ.1)NAMCOL(5,ILYS,NCLLYS(ILYS))='N+ '
               NUMCOL(5,ILYS,NCLLYS(ILYS))=LLYSRS(JLYS)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(5,ILYS,NCLLYS(ILYS))=FCOUL*MIN(1.0,VALUE)
               PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
             END IF
           END IF
           END IF
         END IF
         END DO
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO IARG=1,NARG
         IF(NMASS(5,ILYS)+NMASS(6,IARG).GT.900 .OR.
     $   (NMASS(5,ILYS).GT.400.AND.NMASS(6,IARG).GT.400))THEN
           XCZ=X(LARGCZ(IARG))
           YCZ=Y(LARGCZ(IARG))
           ZCZ=Z(LARGCZ(IARG))
           IF(ABS(XNZ-XCZ).LT.DIS2.AND.
     $        ABS(YNZ-YCZ).LT.DIS2.AND.
     $        ABS(ZNZ-ZCZ).LT.DIS2    )  THEN
             DIS=SQRT((XNZ-XCZ)**2+(YNZ-YCZ)**2+(ZNZ-ZCZ)**2)
             IF(DIS.LT.DIS2)THEN
               IF(PKALYS(ILYS).LT.PKAARG(IARG))THEN
                 NCLLYS(ILYS)=NCLLYS(ILYS)+1
                 NAMCOL(5,ILYS,NCLLYS(ILYS))='ARG'
                 NUMCOL(5,ILYS,NCLLYS(ILYS))=LARGRS(IARG)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(5,ILYS,NCLLYS(ILYS))=FCOUL*MIN(1.0,VALUE)
                 PK1LYS(ILYS)=PK1LYS(ILYS)+VALCOL(5,ILYS,NCLLYS(ILYS))
               ELSE
                 NCLARG(IARG)=NCLARG(IARG)+1
                 NAMCOL(6,IARG,NCLARG(IARG))='LYS'
                 IF(ILYS.EQ.1)NAMCOL(6,IARG,NCLARG(IARG))='N+ '
                 NUMCOL(6,IARG,NCLARG(IARG))=LLYSRS(ILYS)
                 VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
                 VALCOL(6,IARG,NCLARG(IARG))=FCOUL*MIN(1.0,VALUE)
                 PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
               END IF
             END IF
           END IF
         END IF
         END DO
c dmr start
C
C           - FIND N3+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
C
         DO IN3=1,NLIGN3
C
           IF(NMASS(5,Ilys)+NMASS(7,IN3).GT.900 .OR.
     $        (NMASS(5,Ilys).GT.400.AND.NMASS(7,IN3).GT.400))THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           IF(ABS(Xnz-XN).LT.DIS2.AND.
     $        ABS(Ynz-YN).LT.DIS2.AND.
     $        ABS(Znz-ZN).LT.DIS2    ) THEN
             DIS=SQRT((Xnz-XN)**2+(Ynz-YN)**2+(Znz-ZN)**2)
             IF(DIS.LT.DIS2)THEN
              IF(PKAlys(Ilys).LT.PKALGN3(IN3))THEN
               NCLllys(Ilys)=NCLllys(Ilys)+1
               NAMLCOL(5,Ilys,NCLllys(Ilys))=NAMRES(LLIGN3(IN3))
               NBLCOL(5,Ilys,NCLllys(Ilys))=NBATM(LLIGN3(IN3))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(5,Ilys,NCLllys(Ilys))=FCOUL*MIN(1.0,VALUE)
               PK1lys(Ilys)=PK1lys(Ilys)+VALLCOL(5,Ilys,NCLllys(Ilys))
C
              ELSE
               NCLN3(IN3)=NCLN3(IN3)+1
               NAMCOL(7,IN3,NCLN3(IN3))='LYS'
               NUMCOL(7,IN3,NCLN3(IN3))=LlysRS(Ilys)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(7,IN3,NCLN3(IN3))=FCOUL*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALCOL(7,IN3,NCLN3(IN3))
             END IF
             END IF
           END IF
           END IF
         END DO
c
c
C           - FIND Npl+ -
C  -- Np1 pair GROUP -CN2 (eg. 1K1I.pdb) --
C
         DO ic2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           IF(NMASS(5,Ilys)+NMASS(9,Ic2).GT.900 .OR.
     $        (NMASS(5,Ilys).GT.400.AND.NMASS(9,Ic2).GT.400))THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           IF(ABS(Xnz-Xc).LT.DIS2.AND.
     $        ABS(Ynz-Yc).LT.DIS2.AND.
     $        ABS(Znz-Zc).LT.DIS2    ) THEN
             DIS=SQRT((Xnz-Xc)**2+(Ynz-Yc)**2+(Znz-Zc)**2)
               IF(DIS.LT.DIS2)THEN
              IF(PKAlys(Ilys).LT.PKALGcn2(Ic2))THEN
               NCLllys(Ilys)=NCLllys(Ilys)+1
               NAMLCOL(5,Ilys,NCLllys(Ilys))=NAMRES(LLIGNpl(Ic2))
               NBLCOL(5,Ilys,NCLllys(Ilys))=NBATM(LLIGNpl(Ic2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(5,Ilys,NCLllys(Ilys))=FCOUL*MIN(1.0,VALUE)
               PK1lys(Ilys)=PK1lys(Ilys)+VALLCOL(5,Ilys,NCLllys(Ilys))
C
              ELSE
               NCLcn2(Ic2)=NCLcn2(Ic2)+1
               NAMCOL(9,Ic2,NCLcn2(Ic2))='LYS'
               NUMCOL(9,Ic2,NCLcn2(Ic2))=LlysRS(Ilys)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,Ic2,NCLcn2(Ic2))=FCOUL*MIN(1.0,VALUE)
               PK1LGcn2(Ic2)=PK1LGcn2(Ic2)+VALCOL(9,Ic2,NCLcn2(Ic2))
              END IF
               END IF
             END IF
          ENDIF
          ENDIF
         END DO
c
C  -- GUADININIUM TYPE GROUP Cg (eg. 1K1O.pdb) --
C
         DO IC2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           IF(NMASS(5,Ilys)+NMASS(9,Ic2).GT.900 .OR.
     $        (NMASS(5,Ilys).GT.400.AND.NMASS(9,Ic2).GT.400))THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           IF(ABS(Xnz-Xc).LT.DIS2.AND.
     $        ABS(Ynz-Yc).LT.DIS2.AND.
     $        ABS(Znz-Zc).LT.DIS2    ) THEN
             DIS=SQRT((Xnz-Xc)**2+(Ynz-Yc)**2+(Znz-Zc)**2)
               IF(DIS.LT.DIS2)THEN
              IF(PKAlys(Ilys).LT.PKALGcg(Ic2))THEN
               NCLllys(Ilys)=NCLllys(Ilys)+1
               NAMLCOL(5,Ilys,NCLllys(Ilys))=NAMRES(LLIGNpl(Ic2))
               NBLCOL(5,Ilys,NCLllys(Ilys))=NBATM(LLIGNpl(Ic2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(5,Ilys,NCLllys(Ilys))=FCOUL*MIN(1.0,VALUE)
               PK1lys(Ilys)=PK1lys(Ilys)+VALLCOL(5,Ilys,NCLllys(Ilys))
C
              ELSE
               NCLcg(Ic2)=NCLcg(Ic2)+1
               NAMCOL(9,Ic2,NCLcg(Ic2))='LYS'
               NUMCOL(9,Ic2,NCLcg(Ic2))=LlysRS(Ilys)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,Ic2,NCLcg(Ic2))=FCOUL*MIN(1.0,VALUE)
               PK1LGcg(Ic2)=PK1LGcg(Ic2)+VALCOL(9,Ic2,NCLcg(Ic2))
              END IF
               END IF
             END IF
          ENDIF
          ENDIF
         END DO
C
C
C           - FIND Oco(-)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILIG=1,NLIGC2
         IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                     .AND.ZCac(ILIG).NE.0)THEN
         IF(PKALGCAR(ILIG).LT.PKAlys(Ilys))THEN
           IF(NMASS(11,ILIG)+NMASS(5,Ilys).GT.900 .OR.
     $        (NMASS(11,ILIG).GT.400.AND.NMASS(5,Ilys).GT.400))THEN
           XO1=XLO1(ILIG)
           YO1=YLO1(ILIG)
           ZO1=ZLO1(ILIG)
           XO2=XLO2(ILIG)
           YO2=YLO2(ILIG)
           ZO2=ZLO2(ILIG)
           XO=(XO1+XO2)/2.0
           YO=(YO1+YO2)/2.0
           ZO=(ZO1+ZO2)/2.0
           IF(ABS(XO-Xnz).LT.DIS2.AND.
     $        ABS(YO-Ynz).LT.DIS2.AND.
     $        ABS(ZO-Znz).LT.DIS2    )  THEN
             DIS=SQRT((XO-Xnz)**2+(YO-Ynz)**2+(ZO-Znz)**2)
             IF(DIS.LT.DIS2)THEN
               NCLllys(Ilys)=NCLllys(Ilys)+1
               NAMLCOL(5,Ilys,NCLllys(Ilys))=NAMRES(LLIGC2(ILIG))
               NBLCOL(5,Ilys,NCLllys(Ilys))=NBATM(LLIGC2(ILIG))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(5,Ilys,NCLllys(Ilys))=-FCOUL*MIN(1.0,VALUE)
               PK1lys(Ilys)=PK1lys(Ilys)+VALLCOL(5,Ilys,NCLllys(Ilys))
C
               NCLLGCAR(ILIG)=NCLLGCAR(ILIG)+1
               NAMCOL(11,ILIG,NCLLGCAR(ILIG))='LYS'
               NUMCOL(11,ILIG,NCLLGCAR(ILIG))=LlysRS(Ilys)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(11,ILIG,NCLLGCAR(ILIG))=FCOUL*MIN(1.0,VALUE)
               PK1LGCAR(ILIG)=PK1LGCAR(ILIG)
     $                       +VALCOL(11,ILIG,NCLLGCAR(ILIG))
C
             END IF
           END IF
           END IF
         END IF
         END IF
         END DO
C
C
C           - FIND Nar(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO INar=1,NLIGNar
           IF(NMASS(5,Ilys)+NMASS(8,INar).GT.900 .OR.
     $        (NMASS(5,Ilys).GT.400.AND.NMASS(8,INar).GT.400))THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           IF(ABS(Xnz-XN).LT.DIS2.AND.
     $        ABS(Ynz-YN).LT.DIS2.AND.
     $        ABS(Znz-ZN).LT.DIS2    )  THEN
             DIS=SQRT((Xnz-XN)**2+(Ynz-YN)**2+(Znz-ZN)**2)
             IF(DIS.LT.DIS2)THEN
         IF(PKAlys(Ilys).LT.PKALGNar(INar))THEN
               NCLllys(Ilys)=NCLllys(Ilys)+1
               NAMLCOL(5,Ilys,NCLLlys(Ilys))=NAMRES(LLIGNar(INar))
               NBLCOL(5,Ilys,NCLLlys(Ilys))=NBATM(LLIGNar(INar))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(5,Ilys,NCLLlys(Ilys))=FCOUL*MIN(1.0,VALUE)
               PK1lys(Ilys)=PK1lys(Ilys)
     $                     +VALLCOL(5,Ilys,NCLLlys(Ilys))
C
         ELSE
               NCLNar(INar)=NCLNar(INar)+1
               NAMCOL(8,INar,NCLNar(INar))='LYS'
               NUMCOL(8,INar,NCLNar(INar))=LlysRS(Ilys)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(8,INar,NCLNar(INar))=FCOUL*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALCOL(8,INar,NCLNar(INar))
             END IF
           END IF
           END IF
         END IF
         END DO
c dmr finish
C
       END DO
C
C
C      **********************
C      ARG
C      **********************
C
C
       DO IARG=1, NARG
         K=NSDARG(IARG)
         DO J=K+1,30
           NAMSDC(6,IARG,J)='000'
           NUMSDC(6,IARG,J)=0
           VALSDC(6,IARG,J)=0.0
         END DO
C
         K=NBKARG(IARG)
         DO J=K+1,30
           NAMBKB(6,IARG,J)='000'
           NUMBKB(6,IARG,J)=0
           VALBKB(6,IARG,J)=0.0
         END DO
C
         K=NCLARG(IARG)
         DO J=K+1,30
           NAMCOL(6,IARG,J)='000'
           NUMCOL(6,IARG,J)=0
           VALCOL(6,IARG,J)=0.0
         END DO
C
c dmr str
         K=NCLLarg(Iarg)
         DO J=K+1,30
           NAMLCOL(6,Iarg,J)='000'
           NBLCOL(6,Iarg,J)='0000'
           VALLCOL(6,Iarg,J)=0.0
         END DO
c dmr fin
C
         X1=X(LARGN1(IARG))
         Y1=Y(LARGN1(IARG))
         Z1=Z(LARGN1(IARG))
         X2=X(LARGN2(IARG))
         Y2=Y(LARGN2(IARG))
         Z2=Z(LARGN2(IARG))
         X3=X(LARGN3(IARG))
         Y3=Y(LARGN3(IARG))
         Z3=Z(LARGN3(IARG))
         XCZ=X(LARGCZ(IARG))
         YCZ=Y(LARGCZ(IARG))
         ZCZ=Z(LARGCZ(IARG))
C
C           - FIND ARG(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO JARG=1,NARG
         IF(IARG.NE.JARG.AND.PKAARG(IARG).LT.PKAARG(JARG))THEN
           IF(NMASS(6,IARG)+NMASS(6,JARG).GT.900 .OR.
     $     (NMASS(6,IARG).GT.400.AND.NMASS(6,JARG).GT.400))THEN
           XCZJ=X(LARGCZ(JARG))
           YCZJ=Y(LARGCZ(JARG))
           ZCZJ=Z(LARGCZ(JARG))
           IF(ABS(XCZ-XCZJ).LT.DIS2.AND.
     $        ABS(YCZ-YCZJ).LT.DIS2.AND.
     $        ABS(ZCZ-ZCZJ).LT.DIS2    )  THEN
             DIS=SQRT((XCZ-XCZJ)**2+(YCZ-YCZJ)**2+(ZCZ-ZCZJ)**2)
             IF(DIS.LT.DIS2)THEN
               NCLARG(IARG)=NCLARG(IARG)+1
               NAMCOL(6,IARG,NCLARG(IARG))='ARG'
               NUMCOL(6,IARG,NCLARG(IARG))=LARGRS(JARG)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(6,IARG,NCLARG(IARG))=FCOUL*MIN(1.0,VALUE)
               PK1ARG(IARG)=PK1ARG(IARG)+VALCOL(6,IARG,NCLARG(IARG))
             END IF
           END IF
           END IF
         END IF
         END DO
c dmr start
C
C           - FIND N3+ -
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
C
         DO IN3=1,NLIGN3
C
           IF(NMASS(6,Iarg)+NMASS(7,IN3).GT.900 .OR.
     $        (NMASS(6,Iarg).GT.400.AND.NMASS(7,IN3).GT.400))THEN
           XN=X(LLIGN3(IN3))
           YN=Y(LLIGN3(IN3))
           ZN=Z(LLIGN3(IN3))
           IF(ABS(Xcz-XN).LT.DIS2.AND.
     $        ABS(Ycz-YN).LT.DIS2.AND.
     $        ABS(Zcz-ZN).LT.DIS2    ) THEN
             DIS=SQRT((Xcz-XN)**2+(Ycz-YN)**2+(Zcz-ZN)**2)
             IF(DIS.LT.DIS2)THEN
              IF(PKAarg(Iarg).LT.PKALGN3(IN3))THEN
               NCLlarg(Iarg)=NCLlarg(Iarg)+1
               NAMLCOL(6,Iarg,NCLlarg(Iarg))=NAMRES(LLIGN3(IN3))
               NBLCOL(6,Iarg,NCLlarg(Iarg))=NBATM(LLIGN3(IN3))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(6,Iarg,NCLlarg(Iarg))=FCOUL*MIN(1.0,VALUE)
               PK1arg(Iarg)=PK1arg(Iarg)+VALLCOL(6,Iarg,NCLlarg(Iarg))
C
              ELSE
               NCLN3(IN3)=NCLN3(IN3)+1
               NAMCOL(7,IN3,NCLN3(IN3))='ARG'
               NUMCOL(7,IN3,NCLN3(IN3))=LargRS(Iarg)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(7,IN3,NCLN3(IN3))=FCOUL*MIN(1.0,VALUE)
               PK1LGN3(IN3)=PK1LGN3(IN3)+VALCOL(7,IN3,NCLN3(IN3))
              END IF
             END IF
           END IF
           END IF
         END DO
c
c
C           - FIND Npl+ -
C  -- Np1 pair GROUP -CN2 (eg. 1K1I.pdb) --
C
         DO ic2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
           IF(NMASS(6,Iarg)+NMASS(9,Ic2).GT.900 .OR.
     $        (NMASS(6,Iarg).GT.400.AND.NMASS(9,Ic2).GT.400))THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           IF(ABS(Xcz-Xc).LT.DIS2.AND.
     $        ABS(Ycz-Yc).LT.DIS2.AND.
     $        ABS(Zcz-Zc).LT.DIS2    ) THEN
             DIS=SQRT((Xcz-Xc)**2+(Ycz-Yc)**2+(Zcz-Zc)**2)
               IF(DIS.LT.DIS2)THEN
              IF(PKAarg(Iarg).LT.PKALGcn2(Ic2))THEN
               NCLlarg(Iarg)=NCLlarg(Iarg)+1
               NAMLCOL(6,Iarg,NCLlarg(Iarg))=NAMRES(LLIGNpl(Ic2))
               NBLCOL(6,Iarg,NCLlarg(Iarg))=NBATM(LLIGNpl(Ic2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(6,Iarg,NCLlarg(Iarg))=FCOUL*MIN(1.0,VALUE)
               PK1arg(Iarg)=PK1arg(Iarg)+VALLCOL(6,Iarg,NCLlarg(Iarg))
C
              ELSE
               NCLcn2(Ic2)=NCLcn2(Ic2)+1
               NAMCOL(9,Ic2,NCLcn2(Ic2))='ARG'
               NUMCOL(9,Ic2,NCLcn2(Ic2))=LargRS(Iarg)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,Ic2,NCLcn2(Ic2))=FCOUL*MIN(1.0,VALUE)
               PK1LGcn2(Ic2)=PK1LGcn2(Ic2)+VALCOL(9,Ic2,NCLcn2(Ic2))
              END IF
               END IF
             END IF
          ENDIF
          ENDIF
         END DO
c
C  -- GUADININIUM TYPE GROUP Cg (eg. 1K1O.pdb) --
C
         DO IC2=1,NLIGC2
          IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
           IF(NMASS(6,Iarg)+NMASS(9,Ic2).GT.900 .OR.
     $        (NMASS(6,Iarg).GT.400.AND.NMASS(9,Ic2).GT.400))THEN
           XC=X(LLIGC2(IC2))
           YC=Y(LLIGC2(IC2))
           ZC=Z(LLIGC2(IC2))
           IF(ABS(Xcz-Xc).LT.DIS2.AND.
     $        ABS(Ycz-Yc).LT.DIS2.AND.
     $        ABS(Zcz-Zc).LT.DIS2    ) THEN
             DIS=SQRT((Xcz-Xc)**2+(Ycz-Yc)**2+(Zcz-Zc)**2)
               IF(DIS.LT.DIS2)THEN
              IF(PKAarg(Iarg).LT.PKALGcg(Ic2))THEN
               NCLlarg(Iarg)=NCLlarg(Iarg)+1
               NAMLCOL(6,Iarg,NCLlarg(Iarg))=NAMRES(LLIGNpl(Ic2))
               NBLCOL(6,Iarg,NCLlarg(Iarg))=NBATM(LLIGNpl(Ic2))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(6,Iarg,NCLlarg(Iarg))=FCOUL*MIN(1.0,VALUE)
               PK1arg(Iarg)=PK1arg(Iarg)+VALLCOL(6,Iarg,NCLlarg(Iarg))
C
              ELSE
               NCLcg(Ic2)=NCLcg(Ic2)+1
               NAMCOL(9,Ic2,NCLcg(Ic2))='ARG'
               NUMCOL(9,Ic2,NCLcg(Ic2))=LargRS(Iarg)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(9,Ic2,NCLcg(Ic2))=FCOUL*MIN(1.0,VALUE)
               PK1LGcg(Ic2)=PK1LGcg(Ic2)+VALCOL(9,Ic2,NCLcg(Ic2))
              END IF
               END IF
             END IF
          ENDIF
          ENDIF
         END DO
C
C
C           - FIND Oco(-)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO ILIG=1,NLIGC2
         IF(XCac(ILIG).NE.0.AND.YCac(ILIG).NE.0
     $                     .AND.ZCac(ILIG).NE.0)THEN
         IF(PKALGCAR(ILIG).LT.PKAarg(Iarg))THEN
           IF(NMASS(11,ILIG)+NMASS(6,Iarg).GT.900 .OR.
     $        (NMASS(11,ILIG).GT.400.AND.NMASS(6,Iarg).GT.400))THEN
           XO1=XLO1(ILIG)
           YO1=YLO1(ILIG)
           ZO1=ZLO1(ILIG)
           XO2=XLO2(ILIG)
           YO2=YLO2(ILIG)
           ZO2=ZLO2(ILIG)
           XO=(XO1+XO2)/2.0
           YO=(YO1+YO2)/2.0
           ZO=(ZO1+ZO2)/2.0
           IF(ABS(XO-Xcz).LT.DIS2.AND.
     $        ABS(YO-Ycz).LT.DIS2.AND.
     $        ABS(ZO-Zcz).LT.DIS2    )  THEN
             DIS=SQRT((XO-Xcz)**2+(YO-Ycz)**2+(ZO-Zcz)**2)
             IF(DIS.LT.DIS2)THEN
               NCLlarg(Iarg)=NCLlarg(Iarg)+1
               NAMLCOL(6,Iarg,NCLlarg(Iarg))=NAMRES(LLIGC2(ILIG))
               NBLCOL(6,Iarg,NCLlarg(Iarg))=NBATM(LLIGC2(ILIG))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(6,Iarg,NCLlarg(Iarg))=-FCOUL*MIN(1.0,VALUE)
               PK1arg(Iarg)=PK1arg(Iarg)+VALLCOL(6,Iarg,NCLlarg(Iarg))
C
               NCLLGCAR(ILIG)=NCLLGCAR(ILIG)+1
               NAMCOL(11,ILIG,NCLLGCAR(ILIG))='ARG'
               NUMCOL(11,ILIG,NCLLGCAR(ILIG))=LargRS(Iarg)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(11,ILIG,NCLLGCAR(ILIG))=FCOUL*MIN(1.0,VALUE)
               PK1LGCAR(ILIG)=PK1LGCAR(ILIG)
     $                       +VALCOL(11,ILIG,NCLLGCAR(ILIG))
C
             END IF
           END IF
           END IF
         END IF
         END IF
         END DO
C
C
C           - FIND Nar(+)
C
         FCOUL=-2.40
         DIS1=4.00
         DIS2=7.00
         DO INar=1,NLIGNar
           IF(NMASS(6,Iarg)+NMASS(8,INar).GT.900 .OR.
     $        (NMASS(6,Iarg).GT.400.AND.NMASS(8,INar).GT.400))THEN
           XN=X(LLIGNar(INar))
           YN=Y(LLIGNar(INar))
           ZN=Z(LLIGNar(INar))
           IF(ABS(Xcz-XN).LT.DIS2.AND.
     $        ABS(Ycz-YN).LT.DIS2.AND.
     $        ABS(Zcz-ZN).LT.DIS2    )  THEN
             DIS=SQRT((Xcz-XN)**2+(Ycz-YN)**2+(Zcz-ZN)**2)
             IF(DIS.LT.DIS2)THEN
         IF(PKAarg(Iarg).LT.PKALGNar(INar))THEN
               NCLlarg(Iarg)=NCLlarg(Iarg)+1
               NAMLCOL(6,Iarg,NCLlarg(Iarg))=NAMRES(LLIGNar(INar))
               NBLCOL(6,Iarg,NCLlarg(Iarg))=NBATM(LLIGNar(INar))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(6,Iarg,NCLlarg(Iarg))=FCOUL*MIN(1.0,VALUE)
               PK1arg(Iarg)=PK1arg(Iarg)
     $                     +VALLCOL(6,Iarg,NCLlarg(Iarg))
C
         ELSE
               NCLNar(INar)=NCLNar(INar)+1
               NAMCOL(8,INar,NCLNar(INar))='ARG'
               NUMCOL(8,INar,NCLNar(INar))=LargRS(Iarg)
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALCOL(8,INar,NCLNar(INar))=FCOUL*MIN(1.0,VALUE)
               PK1LGNar(INar)=PK1LGNar(INar)+VALCOL(8,INar,NCLNar(INar))
             END IF
           END IF
           END IF
         END IF
         END DO
c dmr finish
C
       END DO
C
C
C
C
C
       CONV=.TRUE.
       DO I=1,1000
C        IF(PK1CAR(I).NE.PKACAR(I))
C    $ WRITE(11,*)NAMCAR(I),LCARRS(I),PK1CAR(I),PKACAR(I)
C        IF(PK1HIS(I).NE.PKAHIS(I))
C    $ WRITE(11,*)'HIS',LHISRS(I),PK1HIS(I),PKAHIS(I)
C        IF(PK1CYS(I).NE.PKACYS(I))
C    $ WRITE(11,*)'CYS',LCYSRS(I),PK1CYS(I),PKACYS(I)
C        IF(PK1TYR(I).NE.PKATYR(I))
C    $ WRITE(11,*)'TYR',LTYRRS(I),PK1TYR(I),PKATYR(I)
C        IF(PK1LYS(I).NE.PKALYS(I))
C    $ WRITE(11,*)'LYS',LLYSRS(I),PK1LYS(I),PKALYS(I)
C        IF(PK1ARG(I).NE.PKAARG(I))
C    $ WRITE(11,*)'ARG',LARGRS(I),PK1ARG(I),PKAARG(I)
         IF(PK1CAR(I).NE.PKACAR(I))CONV=.FALSE.
         IF(PK1HIS(I).NE.PKAHIS(I))CONV=.FALSE.
         IF(PK1CYS(I).NE.PKACYS(I))CONV=.FALSE.
         IF(PK1TYR(I).NE.PKATYR(I))CONV=.FALSE.
         IF(PK1LYS(I).NE.PKALYS(I))CONV=.FALSE.
         IF(PK1ARG(I).NE.PKAARG(I))CONV=.FALSE.
         IF(PK1LGN3(I).NE.PKALGN3(I))CONV=.FALSE. 
         IF(PK1LGNar(I).NE.PKALGNar(I))CONV=.FALSE. 
         IF(PK1LGCg(I).NE.PKALGCg(I))CONV=.FALSE.
         IF(PK1LGCN2(I).NE.PKALGCN2(I))CONV=.FALSE.
         IF(PK1LGCAR(I).NE.PKALGCAR(I))CONV=.FALSE.
         PKACAR(I)=PK1CAR(I)
         PKAHIS(I)=PK1HIS(I)
         PKACYS(I)=PK1CYS(I)
         PKATYR(I)=PK1TYR(I)
         PKALYS(I)=PK1LYS(I)
         PKAARG(I)=PK1ARG(I)
         PKALGN3(I)=PK1LGN3(I)
         PKALGNar(I)=PK1LGNar(I)
         PKALGCN2(I)=PK1LGCN2(I)
         PKALGCg(I)=PK1LGCg(I)
         PKALGCAR(I)=PK1LGCAR(I)
       END DO
C
       WRITE(6,*) 'ITER=',ITER
       IF(CONV)GOTO 510
 500   CONTINUE
 510   CONTINUE
C
C
C      **********************
C      STEP 9. PRINT RESULTS
C      **********************
C
       DO ICAR=1, NCAR
          PKACAR(ICAR)=PKACAR(ICAR)+6.00*NHBCAR(ICAR)
       END DO
       DO IHIS=1, NHIS
          PKAHIS(IHIS)=PKAHIS(IHIS)-6.00*NHBHIS(IHIS)
       END DO
C
C
       DO I=1,2
       IF(I.EQ.1)K=6
       IF(I.EQ.2)K=11
C
       WRITE(K,*)' '
       WRITE(K,'(95(1H-))')
       WRITE(K,'(2(1H-),30X,A,30X,2(1H-))')
     $          'PROPKA: A PROTEIN PKA PREDICTOR'
       WRITE(K,'(2(1H-),28X,A,28X,2(1H-))')
     $          'VERSION 1.0,  04/25/2004, IOWA CITY'
       WRITE(K,'(2(1H-),41X,A,41X,2(1H-))')
     $          'BY HUI LI'
C dcb
c      WRITE(K,'(2(1H-),28X,A,28X,2(1H-))')
       WRITE(K,'(2(1H-),28X,A,17X,2(1H-))')
     $          'VERSION 2.0,  11/05/2007, IOWA CITY/COPENHAGEN'
c      WRITE(K,'(2(1H-),41X,A,41X,2(1H-))')
       WRITE(K,'(2(1H-),28X,A,25X,2(1H-))')
     $          'BY DELPHINE C. BAS AND DAVID M. ROGERS'
C dcb end
       WRITE(K,'(2(1H-),25X,41X,25X,2(1H-))')
       WRITE(K,'(2(1H-),22X,A,22X,2(1H-))')
     $          'PROPKA PREDICTS PROTEIN PKA VALUES ACCORDING TO'
       WRITE(K,'(2(1H-),30X,A,30X,2(1H-))')
     $          'THE EMPIRICAL RULES PROPOSED BY'
       WRITE(K,'(2(1H-),23X,A,23X,2(1H-))')
     $          'HUI LI, ANDREW D. ROBERTSON AND JAN H. JENSEN'
       WRITE(K,'(95(1H-))')
       WRITE(K,*)' '
C
C dcb add references in the output file
       WRITE(K,*)' '
       WRITE(K,'(95(1H-))')
       WRITE(K,*)'References: '
       WRITE(K,*)''
       WRITE(K,'(1X,A51,1X,A21)')
     $ 'Very Fast Empirical Prediction and Rationalization',
     $ 'of Protein pKa Values'
       WRITE(K,*)' Hui Li, Andrew D. Robertson and Jan H. Jensen'
       WRITE(K,'(A36,1X,A14,1X,A17)')
     $ 'PROTEINS: Structure, Function, and',
     $ 'Bioinformatics',
     $ '61:704-721 (2005)'
       WRITE(K,*)''
       WRITE(K,'(2X,A47,1X,A35)')
     $ 'Very Fast Prediction and Rationalization of pKa',
     $ 'Values for Protein-Ligand Complexes'
       WRITE(K,*)' Delphine C. Bas, David M. Rogers and Jan H. Jensen'
       WRITE(K,'(A36,1X,A14,1X,A17)')
     $ 'PROTEINS: Structure, Function, and',
     $ 'Bioinformatics',
     $ '73:765-783 (2008)'
       WRITE(K,'(95(1H-))')
       WRITE(K,*)' '
C dcb add warning for identification of groups in the ouput
       WRITE(K,*)' '
       WRITE(K,'(95(1H-))')
       WRITE(K,*)'WARNING ! '
       WRITE(K,*)''
       WRITE(K,'(1X,A57)')
     $ 'CARBOXYL, GUADININIUM, AND AMIDINO GROUPS FROM THE LIGAND'
       WRITE(K,'(1X,A59)')  
     $ 'ARE IDENTIFIED IN THE OUTPUT USING THE CENTRAL CARBON ATOM.'
       WRITE(K,'(95(1H-))')
       WRITE(K,*)' '
C
       WRITE(K,'(7(1H-),3X,5(1H-),3X,6(1H-),3X,9(1H-),
     $        11(1H-),3(3X,13(1H-)))')
       WRITE(K,'(15X,12X,A,3(3X,A13))')
     $   'DESOLVATION  EFFECTS',
     $   '    SIDECHAIN',
     $   '     BACKBONE',
     $   '    COULOMBIC'
       WRITE(K,'(A15,6(3X,A))')
     $   'RESIDUE     pKa ',
     $   'LOCATE',
     $   '  MASSIVE',
     $   '   LOCAL',
     $   'HYDROGEN BOND',
     $   'HYDROGEN BOND',
     $   '  INTERACTION'
       WRITE(K,'(7(1H-),3X,5(1H-),3X,6(1H-),3X,9(1H-),
     $        3X,8(1H-),3(3X,13(1H-)))')
       WRITE(K,*)' '
C
C
       DO ICAR=1,NCAR
         WRITE(K,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   NAMCAR(ICAR), LCARRS(ICAR), PKACAR(ICAR),
     $   TYPCAR(ICAR),TOLMAS(1,ICAR),NMASS(1,ICAR),
     $   TOLLOC(1,ICAR),NLOCAL(1,ICAR),
     $   VALSDC(1,ICAR,1),NAMSDC(1,ICAR,1),NUMSDC(1,ICAR,1),
     $   VALBKB(1,ICAR,1),NAMBKB(1,ICAR,1),NUMBKB(1,ICAR,1),
     $   VALCOL(1,ICAR,1),NAMCOL(1,ICAR,1),NUMCOL(1,ICAR,1)
         DO J=2,30
           IF(J.LE.NSDCAR(ICAR) .OR.
     $        J.LE.NBKCAR(ICAR) .OR.
     $        J.LE.NCLCAR(ICAR)  )THEN
             WRITE(K,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       NAMCAR(ICAR), LCARRS(ICAR),
     $       VALSDC(1,ICAR,J),NAMSDC(1,ICAR,J),NUMSDC(1,ICAR,J),
     $       VALBKB(1,ICAR,J),NAMBKB(1,ICAR,J),NUMBKB(1,ICAR,J),
     $       VALCOL(1,ICAR,J),NAMCOL(1,ICAR,J),NUMCOL(1,ICAR,J)
           END IF
         END DO
         DO IT=1,NLGCAR(ICAR)
         WRITE(K,'(A3,I4,40X,3X,F5.2,1X,A3,A4)')
     $   NAMCAR(ICAR), LCARRS(ICAR),
     $   VALLIG(1,ICAR,IT),NAMLIG(1,ICAR,IT),
     $   LABEL(1,ICAR,IT)
         END DO
         DO IC=1,NCLLCAR(ICAR)
          WRITE(K,'(A3,I4,72X,3X,F5.2,1X,A3,A4)')
     $    NAMCAR(ICAR), LCARRS(ICAR),
     $    VALLCOL(1,ICAR,IC),NAMLCOL(1,ICAR,IC),
     $    NBLCOL(1,ICAR,IC)
         END DO
         WRITE(K,*)' '
       END DO
C
C
       DO IHIS=1,NHIS
         WRITE(K,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'HIS', LHISRS(IHIS), PKAHIS(IHIS),
     $   TYPHIS(IHIS),TOLMAS(2,IHIS),NMASS(2,IHIS),
     $   TOLLOC(2,IHIS),NLOCAL(2,IHIS),
     $   VALSDC(2,IHIS,1),NAMSDC(2,IHIS,1),NUMSDC(2,IHIS,1),
     $   VALBKB(2,IHIS,1),NAMBKB(2,IHIS,1),NUMBKB(2,IHIS,1),
     $   VALCOL(2,IHIS,1),NAMCOL(2,IHIS,1),NUMCOL(2,IHIS,1)
         DO J=2,30
           IF(J.LE.NSDHIS(IHIS) .OR.
     $        J.LE.NBKHIS(IHIS) .OR.
     $        J.LE.NCLHIS(IHIS)    )THEN
             WRITE(K,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'HIS', LHISRS(IHIS),
     $       VALSDC(2,IHIS,J),NAMSDC(2,IHIS,J),NUMSDC(2,IHIS,J),
     $       VALBKB(2,IHIS,J),NAMBKB(2,IHIS,J),NUMBKB(2,IHIS,J),
     $       VALCOL(2,IHIS,J),NAMCOL(2,IHIS,J),NUMCOL(2,IHIS,J)
           END IF
         END DO
         DO IT=1,NLGHIS(IHIS)
         WRITE(K,'(A3,I4,40X,3X,F5.2,1X,A3,A4)')
     $   'HIS', LHISRS(IHIS),
     $   VALLIG(2,IHIS,IT),NAMLIG(2,IHIS,IT),
     $   LABEL(2,IHIS,IT)
         END DO
         DO IC=1,NCLLHIS(IHIS)
          WRITE(K,'(A3,I4,72X,3X,F5.2,1X,A3,A4)')
     $    'HIS', LHISRS(IHIS),
     $    VALLCOL(2,IHIS,IC),NAMLCOL(2,IHIS,IC),
     $    NBLCOL(2,IHIS,IC)
         END DO
         WRITE(K,*)' '
       END DO
C
C
       DO ICYS=1,NCYS
         WRITE(K,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'CYS', LCYSRS(ICYS), PKACYS(ICYS),
     $   TYPCYS(ICYS),TOLMAS(3,ICYS),NMASS(3,ICYS),
     $   TOLLOC(3,ICYS),NLOCAL(3,ICYS),
     $   VALSDC(3,ICYS,1),NAMSDC(3,ICYS,1),NUMSDC(3,ICYS,1),
     $   VALBKB(3,ICYS,1),NAMBKB(3,ICYS,1),NUMBKB(3,ICYS,1),
     $   VALCOL(3,ICYS,1),NAMCOL(3,ICYS,1),NUMCOL(3,ICYS,1)
         DO J=2,30
           IF(J.LE.NSDCYS(ICYS) .OR.
     $        J.LE.NBKCYS(ICYS) .OR.
     $        J.LE.NCLCYS(ICYS)    )THEN
             WRITE(K,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'CYS', LCYSRS(ICYS),
     $       VALSDC(3,ICYS,J),NAMSDC(3,ICYS,J),NUMSDC(3,ICYS,J),
     $       VALBKB(3,ICYS,J),NAMBKB(3,ICYS,J),NUMBKB(3,ICYS,J),
     $       VALCOL(3,ICYS,J),NAMCOL(3,ICYS,J),NUMCOL(3,ICYS,J)
           END IF
         END DO
         DO IT=1,NLGCYS(ICYS)
         WRITE(K,'(A3,I4,40X,3X,F5.2,1X,A3,A4)')
     $   'CYS', LCYSRS(ICYS),
     $   VALLIG(3,ICYS,IT),NAMLIG(3,ICYS,IT),
     $   LABEL(3,ICYS,IT)
         END DO
         DO IC=1,NCLLCYS(ICYS)
          WRITE(K,'(A3,I4,72X,3X,F5.2,1X,A3,A4)')
     $    'CYS', LCYSRS(ICYS),
     $    VALLCOL(3,ICYS,IC),NAMLCOL(3,ICYS,IC),
     $    NBLCOL(3,ICYS,IC)
         END DO
         WRITE(K,*)' '
       END DO
C
C
       DO ITYR=1,NTYR
         WRITE(K,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'TYR', LTYRRS(ITYR), PKATYR(ITYR),
     $   TYPTYR(ITYR),TOLMAS(4,ITYR),NMASS(4,ITYR),
     $   TOLLOC(4,ITYR),NLOCAL(4,ITYR),
     $   VALSDC(4,ITYR,1),NAMSDC(4,ITYR,1),NUMSDC(4,ITYR,1),
     $   VALBKB(4,ITYR,1),NAMBKB(4,ITYR,1),NUMBKB(4,ITYR,1),
     $   VALCOL(4,ITYR,1),NAMCOL(4,ITYR,1),NUMCOL(4,ITYR,1)
         DO J=2,30
           IF(J.LE.NSDTYR(ITYR) .OR.
     $        J.LE.NBKTYR(ITYR) .OR.
     $        J.LE.NCLTYR(ITYR)    )THEN
             WRITE(K,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'TYR', LTYRRS(ITYR),
     $       VALSDC(4,ITYR,J),NAMSDC(4,ITYR,J),NUMSDC(4,ITYR,J),
     $       VALBKB(4,ITYR,J),NAMBKB(4,ITYR,J),NUMBKB(4,ITYR,J),
     $       VALCOL(4,ITYR,J),NAMCOL(4,ITYR,J),NUMCOL(4,ITYR,J)
           END IF
         END DO
         DO IT=1,NLGTYR(ITYR)
         WRITE(K,'(A3,I4,40X,3X,F5.2,1X,A3,A4)')
     $   'TYR', LTYRRS(ITYR),
     $   VALLIG(4,ITYR,IT),NAMLIG(4,ITYR,IT),
     $   LABEL(4,ITYR,IT)
         END DO
         DO IC=1,NCLLTYR(ITYR)
          WRITE(K,'(A3,I4,72X,3X,F5.2,1X,A3,A4)')
     $    'TYR', LTYRRS(ITYR),
     $    VALLCOL(4,ITYR,IC),NAMLCOL(4,ITYR,IC),
     $    NBLCOL(4,ITYR,IC)
         END DO
         WRITE(K,*)' '
       END DO
C
C
       DO ILYS=1,NLYS
       IF(ILYS.EQ.1)THEN
         WRITE(K,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'N+ ', LLYSRS(ILYS), PKALYS(ILYS),
     $   TYPLYS(ILYS),TOLMAS(5,ILYS),NMASS(5,ILYS),
     $   TOLLOC(5,ILYS),NLOCAL(5,ILYS),
     $   VALSDC(5,ILYS,1),NAMSDC(5,ILYS,1),NUMSDC(5,ILYS,1),
     $   VALBKB(5,ILYS,1),NAMBKB(5,ILYS,1),NUMBKB(5,ILYS,1),
     $   VALCOL(5,ILYS,1),NAMCOL(5,ILYS,1),NUMCOL(5,ILYS,1)
         DO J=2,30
           IF(J.LE.NSDLYS(ILYS) .OR.
     $        J.LE.NBKLYS(ILYS) .OR.
     $        J.LE.NCLLYS(ILYS)    )THEN
             WRITE(K,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'N+ ', LLYSRS(ILYS),
     $       VALSDC(5,ILYS,J),NAMSDC(5,ILYS,J),NUMSDC(5,ILYS,J),
     $       VALBKB(5,ILYS,J),NAMBKB(5,ILYS,J),NUMBKB(5,ILYS,J),
     $       VALCOL(5,ILYS,J),NAMCOL(5,ILYS,J),NUMCOL(5,ILYS,J)
           END IF
         END DO
         WRITE(K,*)' '
       ELSE
         WRITE(K,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'LYS', LLYSRS(ILYS), PKALYS(ILYS),
     $   TYPLYS(ILYS),TOLMAS(5,ILYS),NMASS(5,ILYS),
     $   TOLLOC(5,ILYS),NLOCAL(5,ILYS),
     $   VALSDC(5,ILYS,1),NAMSDC(5,ILYS,1),NUMSDC(5,ILYS,1),
     $   VALBKB(5,ILYS,1),NAMBKB(5,ILYS,1),NUMBKB(5,ILYS,1),
     $   VALCOL(5,ILYS,1),NAMCOL(5,ILYS,1),NUMCOL(5,ILYS,1)
         DO J=2,30
           IF(J.LE.NSDLYS(ILYS) .OR.
     $        J.LE.NBKLYS(ILYS) .OR.
     $        J.LE.NCLLYS(ILYS)    )THEN
             WRITE(K,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'LYS', LLYSRS(ILYS),
     $       VALSDC(5,ILYS,J),NAMSDC(5,ILYS,J),NUMSDC(5,ILYS,J),
     $       VALBKB(5,ILYS,J),NAMBKB(5,ILYS,J),NUMBKB(5,ILYS,J),
     $       VALCOL(5,ILYS,J),NAMCOL(5,ILYS,J),NUMCOL(5,ILYS,J)
           END IF
         END DO
c dmr str; cp from Tyr above
         DO IC=1,NCLLlys(Ilys)
          WRITE(K,'(A3,I4,72X,3X,F5.2,1X,A3,A4)')
     $    'LYS', LlysRS(Ilys),
     $    VALLCOL(5,Ilys,IC),NAMLCOL(5,Ilys,IC),
     $    NBLCOL(5,Ilys,IC)
         END DO
c dmr fin
         WRITE(K,*)' '
       END IF
       END DO
C
C
       DO IARG=1,NARG
         WRITE(K,'(A3,I4,F8.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   'ARG', LARGRS(IARG), PKAARG(IARG),
     $   TYPARG(IARG),TOLMAS(6,IARG),NMASS(6,IARG),
     $   TOLLOC(6,IARG),NLOCAL(6,IARG),
     $   VALSDC(6,IARG,1),NAMSDC(6,IARG,1),NUMSDC(6,IARG,1),
     $   VALBKB(6,IARG,1),NAMBKB(6,IARG,1),NUMBKB(6,IARG,1),
     $   VALCOL(6,IARG,1),NAMCOL(6,IARG,1),NUMCOL(6,IARG,1)
         DO J=2,30
           IF(J.LE.NSDARG(IARG) .OR.
     $        J.LE.NBKARG(IARG) .OR.
     $        J.LE.NCLARG(IARG)    )THEN
             WRITE(K,'(A3,I4,40X,3(3X,F5.2,1X,A3,I4))')
     $       'ARG', LARGRS(IARG),
     $       VALSDC(6,IARG,J),NAMSDC(6,IARG,J),NUMSDC(6,IARG,J),
     $       VALBKB(6,IARG,J),NAMBKB(6,IARG,J),NUMBKB(6,IARG,J),
     $       VALCOL(6,IARG,J),NAMCOL(6,IARG,J),NUMCOL(6,IARG,J)
           END IF
         END DO
c dmr str
         DO IC=1,NCLLarg(Iarg)
          WRITE(K,'(A3,I4,72X,3X,F5.2,1X,A3,A4)')
     $    'ARG', LargRS(Iarg),
     $    VALLCOL(6,Iarg,IC),NAMLCOL(6,Iarg,IC),
     $    NBLCOL(6,Iarg,IC)
         END DO
c dmr fin
         WRITE(K,*)' '
       END DO
C
C
c dmr start
       DO IN3=1,NLIGN3
         WRITE(K,'(A3,1X,A4,F7.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   NAMRES(LLIGN3(IN3)),NBATM(LLIGN3(IN3)), PKALGN3(IN3),
     $   TYPLGN3(IN3),TOLMAS(7,IN3),NMASS(7,IN3),
     $   TOLLOC(7,IN3),NLOCAL(7,IN3), 
     $   VALSDC(7,IN3,1),
     $   NAMSDC(7,IN3,1),NUMSDC(7,IN3,1),
     $   VALbkb(7,IN3,1),NAMbkb(7,IN3,1),
     $   NUMbkb(7,IN3,1),
     $   VALCOL(7,IN3,1),NAMCOL(7,IN3,1),
     $   NUMCOL(7,IN3,1)
         DO J=2,30
           IF(J.LE.NSDN3(IN3) .OR.
     $        J.LE.NBKN3(IN3) .OR.
     $        J.LE.NCLN3(IN3)     )THEN
             WRITE(K,'(A3,1X,A4,39X,3(3X,F5.2,1X,A3,I4))')
     $       NAMRES(LLIGN3(IN3)),NBATM(LLIGN3(IN3)), 
     $       VALSDC(7,IN3,J),NAMSDC(7,IN3,J),NUMSDC(7,IN3,J),
     $       VALbkb(7,IN3,j),NAMbkb(7,IN3,j),
     $       NUMbkb(7,IN3,j),
     $       VALCOL(7,IN3,J),NAMCOL(7,IN3,J),
     $       NUMCOL(7,IN3,J) 
           END IF
         END DO
         WRITE(K,*)' '
       END DO
c dmr finish
C
C
c dmr start
       DO INar=1,NLIGNar
         WRITE(K,'(A3,1X,A4,F7.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   NAMRES(LLIGNar(INar)),NBATM(LLIGNar(INar)), PKALGNar(INar),
     $   TYPLGNar(INar),TOLMAS(8,INar),NMASS(8,INar),
     $   TOLLOC(8,INar),NLOCAL(8,INar),
     $   VALSDC(8,INar,1),
     $   NAMSDC(8,INar,1),NUMSDC(8,INar,1),
     $   VALbkb(8,INar,1),NAMbkb(8,INar,1),
     $   NUMbkb(8,INar,1),
     $   VALCOL(8,INar,1),NAMCOL(8,INar,1),
     $   NUMCOL(8,INar,1)
         DO J=2,30
           IF(J.LE.NSDNar(INar) .OR.
     $        J.LE.NBKNar(INar) .OR.
     $        J.LE.NCLNar(INar)     ) THEN
             WRITE(K,'(A3,1X,A4,39X,3(3X,F5.2,1X,A3,I4))')
     $       NAMRES(LLIGNar(INar)),NBATM(LLIGNar(INar)),
     $       VALSDC(8,INar,J),NAMSDC(8,INar,J),NUMSDC(8,INar,J),
     $       VALbkb(8,INar,j),NAMbkb(8,INar,j),
     $       NUMbkb(8,INar,j),
     $       VALCOL(8,INar,J),NAMCOL(8,INar,J),
     $       NUMCOL(8,INar,J)
           END IF
         END DO
         WRITE(K,*)' '
       END DO
c dmr finish
C
C
c dmr start
       DO IC2=1,NLIGC2
        IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
         WRITE(K,'(A3,1X,A4,F7.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   NAMRES(LLIGC2(IC2)),NBATM(LLIGC2(IC2)), PKALGCg(IC2),
     $   TYPLGCg(IC2),TOLMAS(9,IC2),NMASS(9,IC2),
     $   TOLLOC(9,IC2),NLOCAL(9,IC2),
     $   VALsdc(9,IC2,1),
     $   NAMsdc(9,IC2,1),NUMsdc(9,IC2,1),
     $   VALbkb(9,IC2,1),NAMbkb(9,IC2,1),
     $   NUMbkb(9,IC2,1),
     $   VALCOL(9,IC2,1),NAMCOL(9,IC2,1),
     $   NUMCOL(9,IC2,1)
         DO J=2,30
           IF(J.LE.NSDCg(IC2) .OR.
     $        J.LE.NBKCg(IC2) .OR.
     $        J.LE.NCLCg(IC2)     ) THEN
             WRITE(K,'(A3,1X,A4,39X,3(3X,F5.2,1X,A3,I4))')
     $       NAMRES(LLIGC2(IC2)),NBATM(LLIGC2(IC2)),
     $       VALsdc(9,IC2,J),
     $       NAMsdc(9,IC2,J),NUMsdc(9,IC2,J),
     $       VALbkb(9,IC2,J),NAMbkb(9,IC2,J),
     $       NUMbkb(9,IC2,J),
     $       VALCOL(9,IC2,J),NAMCOL(9,IC2,J),
     $       NUMCOL(9,IC2,J)
           END IF
         END DO
         WRITE(K,*)' '
        END IF
C
        IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
         WRITE(K,'(A3,1X,A4,F7.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   NAMRES(LLIGC2(IC2)),NBATM(LLIGC2(IC2)), PKALGCN2(IC2),
     $   TYPLGCN2(IC2),TOLMAS(9,IC2),NMASS(9,IC2),
     $   TOLLOC(9,IC2),NLOCAL(9,IC2),
     $   VALsdc(9,IC2,1),
     $   NAMsdc(9,IC2,1),NUMsdc(9,IC2,1),
     $   VALbkb(9,IC2,1),NAMbkb(9,IC2,1),
     $   NUMbkb(9,IC2,1),
     $   VALCOL(9,IC2,1),NAMCOL(9,IC2,1),
     $   NUMCOL(9,IC2,1)
         DO J=2,30
           IF(J.LE.NSDCN2(IC2) .OR.
     $        J.LE.NBKCN2(IC2) .OR.
     $        J.LE.NCLCN2(IC2)     ) THEN
             WRITE(K,'(A3,1X,A4,39X,3(3X,F5.2,1X,A3,I4))')
     $       NAMRES(LLIGC2(IC2)),NBATM(LLIGC2(IC2)),
     $       VALsdc(9,IC2,J),
     $       NAMsdc(9,IC2,J),NUMsdc(9,IC2,J),
     $       VALbkb(9,IC2,J),NAMbkb(9,IC2,J),
     $       NUMbkb(9,IC2,J),
     $       VALCOL(9,IC2,J),NAMCOL(9,IC2,J),
     $       NUMCOL(9,IC2,J)
           END IF
         END DO
         WRITE(K,*)' '
        END IF
       END DO
c dmr finish
C
C
C
C
c dmr below does what? Oco ligand atoms
c
       DO ILCAR=1,NLIGC2
        IF(XCac(ILCAR).NE.0.AND.YCac(ILCAR).NE.0
     $                     .AND.ZCac(ILCAR).NE.0)THEN
         WRITE(K,'(A3,1X,A4,F7.2,3X,A6,2X,F5.2,I5,2X,F5.2,I4,
     $           3(3X,F5.2,1X,A3,I4))')
     $   NAMRES(LLIGC2(ILCAR)),NBATM(LLIGC2(ILCAR)), PKALGCAR(ILCAR),
     $   TYPLGCAR(ILCAR),TOLMAS(11,ILCAR),NMASS(11,ILCAR),
     $   TOLLOC(11,ILCAR),NLOCAL(11,ILCAR),
     $   VALSDC(11,ILCAR,1),
     $   NAMSDC(11,ILCAR,1),NUMSDC(11,ILCAR,1),
     $   VALbkb(11,ILCar,1),NAMbkb(11,ILCar,1),
     $   NUMbkb(11,ILCar,1),
     $   VALCOL(11,ILCAR,1),NAMCOL(11,ILCAR,1),
     $   NUMCOL(11,ILCAR,1)
         DO J=2,30
           IF(J.LE.NSDLCAR(ILCAR) .OR.
     $        J.LE.NBKLCar(ILCAR) .OR.
     $        J.LE.NCLLGCAR(ILCAR)     ) THEN
             WRITE(K,'(A3,1X,A4,39X,3(3X,F5.2,1X,A3,I4))')
     $       NAMRES(LLIGC2(ILCAR)),NBATM(LLIGC2(ILCAR)),
     $       VALSDC(11,ILCAR,J),NAMSDC(11,ILCAR,J),
     $       NUMSDC(11,ILCAR,J),
     $       VALbkb(11,ILCar,J),NAMbkb(11,ILCar,J),
     $       NUMbkb(11,ILCar,J),
     $       VALCOL(11,ILCAR,J),NAMCOL(11,ILCAR,J),
     $       NUMCOL(11,ILCAR,J)
           END IF
         END DO
         WRITE(K,*)' '
        ENDIF
       END DO
c dmr
C
C
C
C
       WRITE(K,'(95(1H-))')
       WRITE(K,'(A)')'SUMMARY OF THIS PREDICTION'
       WRITE(K,'(3x,a8,a7,3x,a7,3x,a16)')'RESIDUE','pKa','pKmodel'
     $   ,'ligand atom-type'
       DO ICAR=1,NCAR
         IF(NAMCAR(ICAR).EQ.'ASP')THEN
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $     NAMCAR(ICAR), LCARRS(ICAR), typch(lcaro1(icar)), PKACAR(ICAR)
     $    ,3.80
         END IF
       END DO
       DO ICAR=1,NCAR
         IF(NAMCAR(ICAR).EQ.'GLU')THEN
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $     NAMCAR(ICAR), LCARRS(ICAR), typch(lcaro1(icar)), PKACAR(ICAR)
     $    ,4.50
         END IF
       END DO
       DO ICAR=1,NCAR
         IF(NAMCAR(ICAR).EQ.'C- ')THEN
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $     NAMCAR(ICAR), LCARRS(ICAR), typch(lcaro1(icar)), PKACAR(ICAR)
     $    ,3.20
         END IF
       END DO
       DO IHIS=1,NHIS
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $   'HIS', LHISRS(IHIS), typch(lhiscg(ihis)), PKAHIS(IHIS),6.50
       END DO
       DO ICYS=1,NCYS
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $   'CYS', LCYSRS(ICYS), typch(lcyssg(icys)), PKACYS(ICYS),9.00
       END DO
       DO ITYR=1,NTYR
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $   'TYR', LTYRRS(ITYR), typch(ltyroh(ityr)), PKATYR(ITYR),10.00
       END DO
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $   'N+ ', LLYSRS(1), typch(llysnz(1)), PKALYS(1),8.00
       DO ILYS=2,NLYS
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $   'LYS', LLYSRS(ILYS), typch(llysnz(ilys)), PKALYS(ILYS),10.50
       END DO
       DO IARG=1,NARG
           WRITE(K,'(3X,A3,I4,a1,F7.2,3x,f7.2)')
     $   'ARG', LARGRS(IARG), typch(largn1(iarg)), PKAARG(IARG),12.50
       END DO

       DO IN3=1,NLIGN3
         WRITE(K,'(3X,A3,1a,A4,F6.2,3x,f7.2,3x,a3)')
     $   NAMRES(LLIGN3(IN3)),typch(llign3(in3)),NBATM(LLIGN3(IN3))
     $  ,PKALGN3(IN3),pkamod(7,in3),'N3'
       END DO
       DO INar=1,NLIGNar
         WRITE(K,'(3X,A3,1a,A4,F6.2,3x,f7.2,3x,a3)')
     $   NAMRES(LLIGNar(INar)),typch(llignar(inar)),NBATM(LLIGNar(INar))
     $  ,PKALGNar(INar),pkamod(8,inar),'Nar'
       END DO
       DO IC2=1,NLIGC2 
        IF(NAMATM(LLIGC2(IC2)).EQ.'  CN2') THEN
         WRITE(K,'(3X,A3,1a,A4,F6.2,3x,f7.2,3x,a3)')
     $   NAMRES(LLIGC2(IC2)),typch(lligc2(ic2)),NBATM(LLIGC2(IC2))
     $  ,PKALGCN2(IC2),pkamod(9,ic2),'Npl'
        END IF 
        IF(NAMATM(LLIGC2(IC2)).EQ.'  Cg ') THEN
         WRITE(K,'(3X,A3,1a,A4,F6.2,3x,f7.2,3x,a3)')
     $   NAMRES(LLIGC2(IC2)),typch(lligc2(ic2)),NBATM(LLIGC2(IC2))
     $  ,PKALGCg(IC2),pkamod(9,ic2),'Npl'
        END IF
       END DO
       DO ICAR=1,NLIGC2
        IF(XCac(ICAR).NE.0.AND.YCac(ICAR).NE.0
     $                    .AND.ZCac(ICAR).NE.0)THEN
         WRITE(K,'(3X,A3,1a,A4,F6.2,3x,f7.2,3x,a3)')
     $   NAMRES(LLIGC2(ICAR)),typch(lligc2(icar)),NBATM(LLIGC2(ICAR))
     $  ,PKALGCAR(ICAR),pkamod(11,icar),'Oco'
        END IF
       END DO
       WRITE(K,'(95(1H-))')
C
       END DO
C
       CLOSE(11)
C
C
       END

      subroutine charatm(grpid,igrp,xgrp,ygrp,zgrp,
     $ NCLLgrp,NAMLCOL,NAMRES,NBLCOL,NBATM,VALLCOL,PK1grp,
     $ nmass,nligchr,lligchr,lligchrval,x,y,z)
C
C           - FIND Charged Atom -
C
      parameter(np=1000,npl=100000,idchr=12)
      dimension ncllgrp(np),namlcol(20,np,30),namres(npl),
     $          nblcol(20,np,30),nbatm(npl),vallcol(20,np,30),
     $          pk1grp(np),nmass(20,np),lligchr(np),lligchrval(np),
     $          x(npl),y(npl),z(npl)
      character namres*3,nbatm*4,namlcol*3,nblcol*4
      integer grpid
c
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO Ica=1,NLIGchr
           IF(NMASS(grpid,igrp).GT.400.and.nmass(idchr,ica).gt.400
     $     .or.(NMASS(grpid,igrp)+nmass(idchr,ica).gt.900))THEN
           Xa=X(LLIGchr(Ica))
           Ya=Y(LLIGchr(Ica))
           Za=Z(LLIGchr(Ica))
           IF(ABS(Xgrp-Xa).LT.DIS2.AND.
     $        ABS(Ygrp-Ya).LT.DIS2.AND.
     $        ABS(Zgrp-Za).LT.DIS2    ) THEN
             DIS=SQRT((Xgrp-Xa)**2+(Ygrp-Ya)**2+(Zgrp-Za)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLgrp(Igrp)=NCLLgrp(Igrp)+1
               NAMLCOL(grpid,Igrp,NCLLgrp(Igrp))=NAMRES(LLIGchr(Ica))
               NBLCOL(grpid,Igrp,NCLLgrp(Igrp))=NBATM(LLIGchr(Ica))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(grpid,Igrp,NCLLgrp(Igrp))=lligchrval(ica)*
     $         (-FCOUL)*MIN(1.0,VALUE)
             PK1grp(Igrp)=PK1grp(Igrp)+VALLCOL(grpid,Igrp,NCLLgrp(Igrp))
             END IF
           END IF
           END IF
         END DO
         return
         end

      subroutine charatml(grpid,igrp,jgrp,xgrp,ygrp,zgrp,
     $ NCLLgrp,NAMLCOL,NAMRES,NBLCOL,NBATM,VALLCOL,PK1grp,
     $ nmass,nligchr,lligchr,lligchrval,x,y,z)
C
C           - FIND Charged Atom -
C
      parameter(np=1000,npl=100000,idchr=12)
      dimension ncllgrp(np),namlcol(20,np,30),namres(npl),
     $          nblcol(20,np,30),nbatm(npl),vallcol(20,np,30),
     $          pk1grp(np),nmass(20,np),lligchr(np),lligchrval(np),
     $          x(npl),y(npl),z(npl),jgrp(np)
      character namres*3,nbatm*4,namlcol*3,nblcol*4
      integer grpid
c
         FCOUL=2.40
         DIS1=4.00
         DIS2=7.00
         DO Ica=1,NLIGchr
         if(namres(jgrp(igrp)).ne.namres(lligchr(ica)))then
           IF(NMASS(grpid,igrp).GT.400.and.nmass(idchr,ica).gt.400
     $     .or.(NMASS(grpid,igrp)+nmass(idchr,ica).gt.900))THEN
           Xa=X(LLIGchr(Ica))
           Ya=Y(LLIGchr(Ica))
           Za=Z(LLIGchr(Ica))
           IF(ABS(Xgrp-Xa).LT.DIS2.AND.
     $        ABS(Ygrp-Ya).LT.DIS2.AND.
     $        ABS(Zgrp-Za).LT.DIS2    ) THEN
             DIS=SQRT((Xgrp-Xa)**2+(Ygrp-Ya)**2+(Zgrp-Za)**2)
             IF(DIS.LT.DIS2)THEN
               NCLLgrp(Igrp)=NCLLgrp(Igrp)+1
               NAMLCOL(grpid,Igrp,NCLLgrp(Igrp))=NAMRES(LLIGchr(Ica))
               NBLCOL(grpid,Igrp,NCLLgrp(Igrp))=NBATM(LLIGchr(Ica))
               VALUE=1.0-(DIS-DIS1)/(DIS2-DIS1)
               VALLCOL(grpid,Igrp,NCLLgrp(Igrp))=lligchrval(ica)*
     $         (-FCOUL)*MIN(1.0,VALUE)
             PK1grp(Igrp)=PK1grp(Igrp)+VALLCOL(grpid,Igrp,NCLLgrp(Igrp))
             END IF
           END IF
           END IF
         end if
         END DO
         return
         end
