      subroutine tmatrix(number_of_angles)
c       #########################################################################
c       Date: 2/3/2003:
c       PCK mods started this date to implement the Maxwell-Garnett mixing formula
c       when calculating the complex dielectric function for an ice+water combination.
c
c
c       Date: 1/28/2003:
c       Mods started on this date to make the program simulate a single layer
c       (one refractive index) particle.
c       New refractive index subroutine added to simulate ice / air mixture.
c       These mods were done with advice of Larry Carey, Walt Petersen, and Bringi.
c       This source code version re-named to tmatrix1_pck.f as of 1/28/2003.
c
c       Original PCK documentation from NOvember 2002.
c       Parameters of primary interest (axis ratio, temperature, radar wavelength,
c       etc.) now read in from an external initialization file.
c       NOTE: In subroutine rddata, the refractive index calls have been set
c       so that the inner portion of the scattering particle is ice, and
c       the outer layer is water.
c       #########################################################################

c       Electromagnetic scattering from a lossy dielectric imbedded within
c       another lossy dielectric.
c       exp(-jwt) time convention is used.
c       Perform the numerical integration and fill the a,b,c,d matrices
c       for the outer surface and x,y matrices for the inner surface.

c       Added implicit double precision statements 11/98 wap

c       Unit 8 outputs just the size information & RCS (not used)
c       Unit 7 outputs final T-matrices
c       Unit 6 outputs intermediate parameters for checking convergence
c       after processing each mode for each size

        implicit none

c
c       Program arguments:
c

        integer number_of_angles

c
c       function declarations
c
        double precision cdmpyr, cdmpyi

c
c       program "global" "structs" if you will . . .
c

        double precision dtr, rtd, cpi
        common dtr, rtd, cpi

        double precision dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1                      aovrb, sigma, ib
        common /bdycom/ dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1          aovrb, sigma, ib

        integer nm, kmv
        double precision cmi(10), cmv, cm2, twm, prodm
        common /cmvcom/ nm, kmv, cmi, cmv, cm2, twm, prodm

        double precision qem, b89
        common /variam/ qem, b89

        double precision dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2
        common /bdycm2/ dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2

        integer nrank, nranki
        double precision a(40,40,2), b(40,40,2), c(40,40,2), d(40,40,2),
     1                      x(40,40,2), y(40,40,2), cmxnrm(40)
        common /mtxcom/ nrank, nranki, a, b, c, d, x, y, cmxnrm

        double precision anginc, acans(181,2,2), uang(181), rtsfct,
     1                      dltang
        integer nuang
        common /uvccom/ anginc, acans, uang, rtsfct, dltang, nuang

        integer irow, irow1, icol, icol1
        double precision crow, crowm, crow1, ccol, ccolm, ccol1, crij,
     1                      crssij
        common /rowcol/ irow, irow1, crow, crowm, crow1, icol, icol1,
     1                  ccol, ccolm, ccol1, crij, crssij

        double precision height
        integer ml, ijk
        common /brin/ ml, ijk, height

        double precision epps(4)
        integer nsect
        common /endpnt/ epps, nsect

        double precision summ1, summ2
        common /scatt/ summ1, summ2

        double precision ptopdia, prhoice, pzhgt, ptempt, palamd, paovrb
        common /pck_values/ ptopdia, prhoice, pzhgt, ptempt, palamd,
     1                      paovrb

        double precision temp(48, 8, 4)
        double precision clrmtx(19200), clrtot(724), result(48, 8)
        double precision buff(5)

c
c       local vars
c

        double precision bdyfct, ckrow, em, fctki, fprod, quanm, rhoice
        double precision scale, temati, tematr, topdia, zhgt

        integer i, ibfct, icom1, icom2, iefct, ifct, iflag, im, ith, j
        integer js, k, kk, l, lfct, ll, nn, npts, nsect1

        equivalence (a(1, 1, 1), clrmtx(1)), (acans(1, 1, 1), clrtot(1))

        character * 36 runid
c
c
c       Get the values from the PCK input file.
c       These will be transferred into the "pck_values" named common block.
c       At the appropriate points in the program, the pck_values
c       (which start with "p")
c       will be copied into the variables as they were originally named.
        call pck_inp
c
c       #############################
c       PCK mods:
c       Open a file to hold some run-identifying items.
        open(22, file='pck_tmatrix2b.out', form='formatted')

c       Make a header line entry to designate the start of a new model run
        write(22, '(a56)')
     1  'XT1   KM MSL   INITIAL DIA (CM) LABEL DENSITY (G/CM**3)'

c       Copy in the run-identification variables.
        zhgt = pzhgt
        topdia = ptopdia
        rhoice = prhoice

c       Put the the PCK items into the run identification line in the output file.
        write(runid, 8600) zhgt, topdia, rhoice
 8600   format('RI',3x,f4.1,5x,f4.1,15x,f3.2)
        write(22, '(a36)') runid

c       #############################
c       Open Units 7 and 8 for output
c       Open (unit=7, file='tmatrix2_7')   -smg08
        Open (unit=8, file='spongyice.out')

c       Read the number of input sizes
c       A T matrix is written to unit 7 for each size
c       read(5,*) ml
        read(20, *) ml
c       write(7,*)ml    -smg08
c       ml=1

c       Set data counter to 0
        ijk = 0
        cpi = 4.0 * atan(1.0)
        dtr = cpi / 180.0
        rtd = 180.0 / cpi

c       Read input data and print headings for output files
   20   call rddata

c        print *,'H A'

        iflag = 0
c       Set the number of scattering angles (nuang = 181 maximum and nuang-1
c       should be divisible into 180 by a whole number)
c         nuang= 2
        nuang = number_of_angles
        dltang = DBLE(180 / (nuang - 1))
        uang(1) = anginc

c        print *,'H B'

        do 30 i = 2, nuang
          uang(i) = uang(i - 1) + dltang
   30   continue

c        print *,'H C'

c       Clear the accumulating answer register (used in addprc).
        do 40 j = 1, 724
          clrtot(j) = 0.0
   40   continue

c        print *,'H D'

        summ1 = 0.0
        summ2 = 0.0
        rtsfct = 8.0 / conk

c       Set multiplier b89 dependent on ib value (symmetry indicator).
        b89 = 1.0
        if (ib.eq.8) b89 = 2.0
        bdyfct = 1.0

c        print *,'H E'

c       Set up a loop for each m value.
        do 900 im = 1, nm

c          print *,'H E A'

c         set m dependent variables.
          cmv = cmi(im)
          kmv = cmv
          cm2 = cmv**2
          prodm = 1.0
          if(kmv.gt.0) go to 44

c          print *,'H E B'

          em = 1.0
          go to 60

c          print *,'H E C'

   44     em = 2.0
          quanm = cmv
          do 52 ifct = 1, kmv
            quanm = quanm + 1.0
            prodm = quanm * prodm / 2.0
   52     continue

c          print *,'H E D'

   60     qem = -2.0 / em
          twm = 2.0 * cmv

c         Initialize all matrix areas to zero
          do 80 i = 1, 19200
            clrmtx(i) = 0.0
   80     continue

c          print *,'H E E'

c         Set up a loop for each row of the matrices.
          crow = 0.0
          crowm = cmv
          do 600 irow = 1, nrank
            irow1 = irow + nrank
            crow = crow + 1.0
            crowm = crowm + 1.0
            crow1 = crow + 1.0

c           Set up a loop for each column of the matrices.
            ccol = 0.0
            ccolm = cmv
            do 400 icol = 1, nrank
              icol1 = icol + nrank
              ccol = ccol + 1.0
              ccolm = ccolm + 1.0
              ccol1 = ccol + 1.0
c
c             Calculate matrices a,b,c,d associated with the outer surface,
c             following notation by Peterson and Strom, q1(out,re) is stored
c             in a, q1(re,re) is stored in b,q1(re,out) is stored in c and
c             q1(out,out) is stored in d. This completes the matrices required
c             for the outer surface. Note: T-matrix is not calculated for the
c             outer surface. All matrices are transposed in the following code.
c
c             Calculate matrices x,y associated with the embedded or inner
c             surface. q2(out,re) is stored in x and q2(re,re) is stored in y.
c
c             Perform integration using a sequence of 1,3,7,15,31,63,127 and 255
c             point extended gauss-type quadrature formulae.
c             result(48,k) contains the values of the integrals .
c             There are 48 integrations to
c             be performed for each looping through irow and icol. These correspond
c             to 4 sub-matrix elements for each of the 6 matrices (a,b,c,d,x,y)
c             and assiciated real and imaginary parts.
c
              nsect1 = nsect - 1
              do 301 j=1, nsect1
                js = j
                ith = 0
                call quad(epps(j), epps(j + 1), k, result, npts, ith,
     1                      js)
                do 401 i = 1, 48
                  temp(i, k, j) = result(i, k)
  401           continue
  301         continue

              do 501 j = 1, nsect1
                a(icol, irow1, 1) = temp(1, k, j) + a(icol, irow1, 1)
                a(icol, irow1, 2) = temp(2, k, j) + a(icol, irow1, 2)
                b(icol, irow1, 1) = temp(3, k, j) + b(icol, irow1, 1)
                b(icol, irow1, 2) = temp(4, k, j) + b(icol, irow1, 2)
                c(icol, irow1, 1) = temp(5, k, j) + c(icol, irow1, 1)
                c(icol, irow1, 2) = temp(6, k, j) + c(icol, irow1, 2)
                d(icol, irow1, 1) = temp(7, k, j) + d(icol, irow1, 1)
                d(icol, irow1, 2) = temp(8, k, j) + d(icol, irow1, 2)
                x(icol, irow1, 1) = temp(9, k, j) + x(icol, irow1, 1)
                x(icol, irow1, 2) = temp(10, k, j) + x(icol, irow1, 2)
                y(icol, irow1, 1) = temp(11, k, j) + y(icol, irow1, 1)
                y(icol, irow1, 2) = temp(12, k, j) + y(icol, irow1, 2)
                a(icol1, irow, 1) = temp(13, k, j) + a(icol1, irow, 1)
                a(icol1, irow, 2) = temp(14, k, j) + a(icol1, irow, 2)
                b(icol1, irow, 1) = temp(15, k, j) + b(icol1, irow, 1)
                b(icol1, irow, 2) = temp(16, k, j) + b(icol1, irow, 2)
                c(icol1, irow, 1) = temp(17, k, j) + c(icol1, irow, 1)
                c(icol1, irow, 2) = temp(18, k, j) + c(icol1, irow, 2)
                d(icol1, irow, 1) = temp(19, k, j) + d(icol1, irow, 1)
                d(icol1, irow, 2) = temp(20, k, j) + d(icol1, irow, 2)
                x(icol1, irow, 1) = temp(21, k, j) + x(icol1, irow, 1)
                x(icol1, irow, 2) = temp(22, k, j) + x(icol1, irow, 2)
                y(icol1, irow, 1) = temp(23, k, j) + y(icol1, irow, 1)
                y(icol1, irow, 2) = temp(24, k, j) + y(icol1, irow, 2)
                a(icol1, irow1, 1) = temp(25, k, j) + a(icol1, irow1, 1)
                a(icol1, irow1, 2) = temp(26, k, j) + a(icol1, irow1, 2)
                b(icol1, irow1, 1) = temp(27, k, j) + b(icol1, irow1, 1)
                b(icol1, irow1, 2) = temp(28, k, j) + b(icol1, irow1, 2)
                c(icol1, irow1, 1) = temp(29, k, j) + c(icol1, irow1, 1)
                c(icol1, irow1, 2) = temp(30, k, j) + c(icol1, irow1, 2)
                d(icol1, irow1, 1) = temp(31, k, j) + d(icol1, irow1, 1)
                d(icol1, irow1, 2) = temp(32, k, j) + d(icol1, irow1, 2)
                x(icol1, irow1, 1) = temp(33, k, j) + x(icol1, irow1, 1)
                x(icol1, irow1, 2) = temp(34, k, j) + x(icol1, irow1, 2)
                y(icol1, irow1, 1) = temp(35, k, j) + y(icol1, irow1, 1)
                y(icol1, irow1, 2) = temp(36, k, j) + y(icol1, irow1, 2)
                a(icol, irow, 1) = temp(37, k, j) + a(icol, irow, 1)
                a(icol, irow, 2) = temp(38, k, j) + a(icol, irow, 2)
                b(icol, irow, 1) = temp(39, k, j) + b(icol, irow, 1)
                b(icol, irow, 2) = temp(40, k, j) + b(icol, irow, 2)
                c(icol, irow, 1) = temp(41, k, j) + c(icol, irow, 1)
                c(icol, irow, 2) = temp(42, k, j) + c(icol, irow, 2)
                d(icol, irow, 1) = temp(43, k, j) + d(icol, irow, 1)
                d(icol, irow, 2) = temp(44, k, j) + d(icol, irow, 2)
                x(icol, irow, 1) = temp(45, k, j) + x(icol, irow, 1)
                x(icol, irow, 2) = temp(46, k, j) + x(icol, irow, 2)
                y(icol, irow, 1) = temp(47, k, j) + y(icol, irow, 1)
                y(icol, irow, 2) = temp(48, k, j) + y(icol, irow, 2)
  501         continue

c             if (iflag.eq.0) write(6,101) npts
c 101         format(1x,'Number of Gauss points used= ',i5)
              iflag = iflag + 1
  400       continue

c           Calculate the normalization factor (used in addprc)
            ckrow = irow
            if (kmv.gt.0) go to 426

            fctki = 1.0
            goto 440

  426       if (irow.ge.kmv) go to 430

            cmxnrm(irow) = 1.0
            goto 600

  430       ibfct = irow - kmv + 1
            iefct = irow + kmv
            fprod = ibfct
            fctki = 1.0

            do 432 lfct = ibfct, iefct
              fctki = fctki * fprod
              fprod = fprod + 1.0
  432       continue

  440       cmxnrm(irow) = 4. * ckrow * (ckrow + 1.) * fctki / (em * (2.
     1      * ckrow + 1.))
  600     continue

c          print *,'H E F'

c         write(*,*) 'after 600'
          nn = 2 * nrank

c         Process computed matrices

c         Calculate t(2) and store in y. x is left free. also print y.
          call prcssm(x, y, nrank, nranki)

c          print *,'H E G'

c         write(*,*) 'after prcssm'
c         call printm(y,nn,40)
c         normalise t - matrix  of the inner body before coupling
c         with the outer body matrices. this change made on nov 21/84.
          do 1008 l = 1, nrank
            ll=l+nrank

            do 1008 k = 1, nrank
              kk = k + nrank
              scale = cmxnrm(l) / cmxnrm(k)

              do 1008 j = 1, 2
                y(l, k, j) = y(l, k, j) * scale
                y(l, kk, j) = y(l, kk, j) * scale
                y(ll, k, j) = y(ll, k, j) * scale
                y(ll, kk, j) = y(ll, kk, j) * scale
 1008     continue

c          print *,'H E H'

c         Calculate b-t(2)*c and store in x.
          do 1001 l = 1, nn
            do 1002 k = 1 , nn
              tematr = 0.0
              temati = 0.0

              do 1003 j = 1, nn
                tematr = cdmpyr(y(l, j, 1), y(l, j, 2), c(j, k, 1), c(j,
     1                          k, 2)) + tematr

                temati = cdmpyi(y(l, j, 1), y(l, j, 2), c(j, k, 1), c(j,
     1                          k, 2)) + temati
 1003         continue

            x(l, k, 1) = b(l, k, 1) - tematr
            x(l, k, 2) = b(l, k, 2) - temati

 1002       continue
 1001     continue

c          print *,'H E I'

c         Calculate a-t(2)*d and store in b.
          do 1004 l = 1, nn
            do 1005 k = 1, nn
              tematr = 0.0
              temati = 0.0

              do 1006 j = 1, nn
                tematr = cdmpyr(y(l, j, 1), y(l, j, 2), d(j, k, 1),
     1                          d(j, k, 2)) + tematr
                temati = cdmpyi(y(l, j, 1), y(l, j, 2), d(j, k, 1),
     1                          d(j, k, 2)) + temati
 1006         continue

              b(l, k, 1) = a(l, k, 1) - tematr
              b(l, k, 2) = a(l, k, 2) - temati

 1005       continue
 1004     continue

c          print *,'H E J'

c         Calculate t(1,2)=b(inverse)*x and store in x.
c         write(*,*) 'prcssm'
          call prcssm(b, x, nrank, nranki)

c          print *,'H E K'

c         Print T matrix for this mode
c         write(7,11) cmi(im)  -smg08
c 11       format(2x,f10.5)
c          do 7  icom1 = 1, 2
c            do 7  icom2 = 1, nn
c          write(7,9) (x(icom2,icom3,icom1),icom3=1,nn)  -smg08
c 9        format(2x,8e15.7) -smg08
c 7        continue

c         call printm(x,nn,40)
          call addprc

c          print *,'H E L'

  900   continue

c        print *,'H F'

        if (ijk.le.ml) go to 20
c       closing particle model file:
        close(20)
      end

c--------------------------------------------------------------
      subroutine rddata
c       Subroutine to read input data for the scattering program

        implicit none

        double precision dtr, rtd, cpi
        common dtr, rtd, cpi

        integer kgauss
        common /gauss/ kgauss

        double precision a(40, 40, 2), b(40, 40, 2), c(40, 40, 2)
        double precision d(40, 40, 2), x(40, 40, 2), y(40, 40, 2)
        double precision cmxnrm(40)
        integer nrank, nranki
        common /mtxcom/ nrank, nranki, a, b, c, d, x, y, cmxnrm

        double precision cmi(10), cmv, cm2, twm, prodm
        integer nm, kmv
        common /cmvcom/ nm, kmv, cmi, cmv, cm2, twm, prodm

        double precision dcnr, dcni, ckprr, ckpri, ckr, dckr, conk
        double precision aovrb, sigma
        integer ib
        common /bdycom/ dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1                  aovrb,sigma,ib

        double precision anginc, acans(181, 2, 2), uang(181), rtsfct
        double precision dtlang
        integer nuang
        common /uvccom/ anginc, acans, uang, rtsfct, dtlang, nuang

        double precision dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2
        common /bdycm2/ dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2

        double precision height
        integer ml, ijk
        common /brin/ ml, ijk, height

        double precision epps(4)
        integer nsect
        common /endpnt/ epps, nsect

        double precision dpart
        common /vivek2/ dpart

        double precision ptopdia, prhoice, pzhgt, ptempt, palamd, paovrb
        common /pck_values/ ptopdia, prhoice, pzhgt, ptempt, palamd,
     1                      paovrb

        double precision pvw
        common /sponge/ pvw

        double precision epdeg(4), alamd, dcore, epsinnerimag, vi, vw
        double precision epsinnerreal, epsouterimag, epsouterreal, tempt

        integer i, iwet

c       Check to see if all input data sets are read
        if (ijk.eq.ml) go to 200

        print *,'on particle ', ijk + 1,' of ', ml

c       Read necessary input data
c       nm    = number of m values (modes)
c       nrank = n value (matrix order)
c       nsect = number of sections in the body
c       ib    = symmetry code  ib = 8 for mirror symmetry about
c       theta= 90 degrees; ib = 9 for general shaped body
c       anginc= angle of incidence of the incident wave.

c       read(5,*) nm,nrank,nsect,ib,anginc

        nrank = 17
        nsect = 2
        ib = 8
        nm = 7
        anginc = 90.
        nranki = nrank + 1

c        print *,'A'

c       Read in parameters for the outer body
c       conk  = ka of body (size parameter)
c       aovrb = a/b ratio of outer spheroid (<1 for oblate, >1 for prolate)
c       sigma = 4/k (used with scattered field subroutine addprc)
c       dcnr  = real part of dielectric constant of outer body
c       dcni  = imaginary part of dielectric constant of outer body
c       dcnr and dcni are both positive since exp(-jwt)
c       time convention is assumed.

c       Read in parameters for the imbedded body
c       conk2 = ka of embedded body (size parameter)
c       aovrb2= a/b ratio of inner spheroid
c       dcnr2 = real part of dielectric constant of embedded body
c       dcni2 = imaginary part of dielectric constant of embedded body
c       dcnr2 and dcni2 are both positive since exp(-jwt)
c       time convention is assumed.
c       dpart= outer diameter
c       dcore=inner diameter
c       alamd= radar wavelength
c       tempt= temperature of both inner and outer

c       read(5,*) dpart,dcore,aovrb,aovrb2,tempt,alamd
        read(20,*) dpart, dcore, aovrb, aovrb2, epsouterreal,
     1               epsouterimag, epsinnerreal, epsinnerimag

c        print *,'B'

c       For freezing drops first assumption will be that dpart=dcore+1mm; Bringi et
c       al. 1997;  something to test in sensitivity studies. Drop diameters less
c       than 2 mm will be assumed to be uniformly frozen.  This seems reasonable
c       given the r**2 dependence of freezing time and that 1.5-3 minutes seems about
c       right for 4mm drops if we interpret the results of Murray and List and
C       Blanchard (1957).  Hence, small drops are going to be almost all ice
c       in a reasonably short (say 10-30 s) amount of time.
c       Note- for freezing of raindrops where rapid freezing is assumed (outer
c       ice shell, inner water core)- an oblateness relationship based on teh
c       shape of a raindrop of similar diameter is used (WAP 1/99)

c       Prupacher and Beard (1975) (for raindrops)

c       aovrb=1.03-.62*dpart

c       Green (1975) relationship (for raindrops) (WAP; 1/99)

C       aovrb=1.6993*dpart**5-5.3728*dpart**4+6.4728*dpart**3-
C    1       +3.3406*dpart**2+0.033*dpart+1.0001
c
c       ###################
c       PCK MOD
c       Axis ratio, temperature, and radar wavelength transferred in
c       by the following few lines:
c       aovrb=paovrb
c       aovrb2 = aovrb
c       Walt's original comment lines:
c       Use temperature of -10 C for rapid freezing (-3 C for slow freezing and
C       sloshing water case. (0 for Pat's case of 0.8 aovrb, 15-.45 water fractions)
        tempt = ptempt
        alamd = palamd
c       #####################
        conk = cpi * dpart / (alamd * (aovrb**(1. / 3.)))
        conk2 = cpi * dcore / (alamd * (aovrb2**(1. / 3.)))
c       write (*,*) 'conks', conk, conk2
        sigma = 2.0 * alamd / (cpi * 100.0)

c        print *,'C'

c       figos computes conk, conk2, aovrb, and aovrb2 for a given
c       dpart, dcore, and thickness of water coating for the
c       special case of an ice sphere inside of a oblate raindrop
c       See J Atmos Sci, V47, 1990, pp. 549-564 for model geometry
c       call figos(thick,alamd,dpart,dcore,aovrb,conk,conk2)

c       NOTE on 2/3/2003:
c       Following variables used by snwep subroutine are no longer used.
c       Spongy ice dielectric function calculation now none by sbr die_funct,
c       with the water volume fraction passed from the input file through pvw.
c       Dielectric const for dry or wet snow
        iwet = 0
        vi = 1.0
        vw = 0.0
c       #############################################################################
c       PCK Mod on 2/3/2003:
c
c        print *,'D'

c       Larry's subroutine called once here, should return all necessary
c       complex dielectric function values (water, ice, spongy mixture).
c       call die_funct(watr,wati,xicer,xicei,sponjr,sponji)
c       Transfer the various complex dielectric functions to the output file.
c       write(22,8860)
c8860   format('XT',4x,'COMPLX DIELCT FUNCTIONS
c     > VIA LARRY CAREY + MAXWELL-GARNETT')
c       write(22,8861)watr,wati
c8861   format('DW',11x,'WATER      (RE,IM)=',2(2x,f6.3))
c       write(22,8862)xicer,xicei
c8862   format('DI',11x,'ICE        (RE,IM)=',2(2x,f6.3))
c       write(22,8863)sponjr,sponji
c8863   format('DS',11x,'SPONGY ICE (RE,IM)=',2(2x,f6.3))
c
c       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       THIS IS WHERE THE DIELECTRIC PROPERTIES OF THE TWO LAYERS ARE ASSIGNED.
c       The following are for Tracy's spongy hail calculations:
c       Inner (suffix 2) layer is spongy ice:
        dcnr2 = epsinnerreal
        dcni2 = epsinnerimag

c       Outer layer is water:
        dcnr = epsouterreal
        dcni = epsouterimag
c       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c        print *,'E'

c       read(25,*)vi,dcnr2,dcni2
c       ##############################################################################
c
c       Unit 6 outputs intermediate parameters after each mode and
c       for each size; used for checking convergence
c       write(6,"(i15,  5x,'Number of sizes')")ml
c       write(6,"(i15,  5x,'Symmetry code')")ib
c       write(6,"(i15,  5x,'Number of sections')")nsect
c       write(6,"(i15,  5x,'Nrank')")nrank
c       write(6,"(i15,  5x,'Nmode')")nm
c       write(6,"(e15.7,5x,'Outer size parameter')")conk
c       write(6,"(e15.7,5x,'Inner size parameter')")conk2
c       write(6,"(e15.7,5x,'Wavelength in cm')")alamd
c       write(6,"(e15.7,5x,'Temperature in C')")tempt
c       write(6,"(e15.7,5x,'Equivalent particle diameter in cm')")dpart
c       write(6,"(e15.7,5x,'Equivalent inner diameter in cm')")dcore
c       write(6,"(e15.7,5x,'Outer A/B ratio')")aovrb
c       write(6,"(e15.7,5x,'Inner A/B ratio')")aovrb2
c       write(6,"(e15.7,5x,'Inner Re(dielectric const)')")dcnr2
c       write(6,"(e15.7,5x,'Inner Im(dielectric const)')")dcni2
c       write(6,"(e15.7,5x,'Outer Re(dielectric const)')")dcnr
c       write(6,"(e15.7,5x,'Outer Im(dielectric const)')")dcni
c       write(6,"(e15.7,5x,'DENSITY FRACTION OF ICE >NOT USED PCK 2/3/2003')") vi
c
c       #######################
c       Mod on 2/3/2003:
c       Original output list modified to recognize externally-input volume fraction of
c       water (pvw).
c       Put selected items into the PCK output file.
        write(22, 8879)
 8879   format('XT', 4x, 'TEMP C', 2x, 'WAVELENGTH CM', 2x,
     >  'VOLUME FRACTION OF WATER IN SPONGY ICE')

        write(22, 8880) tempt, alamd, pvw
 8880   format('EV',2(2x,f7.3),9x,f4.2)

        write(22, 8881)
 8881   format('XT', 4x, 'AXIS', 5x, 'OUT DIA', 2x, 'IN DIA', 2x,
     >  'OUT DILECT RE IM',3x,'IN DILECT RE IM')

        write(22, 8882) aovrb, dpart, dcore, dcnr, dcni, dcnr2, dcni2
 8882   format('PA', 4(2x, f7.3), 4x, f7.3, 1x, f7.3, 2x, f7.3)
c       #######################

c        print *,'F'

c       Unit 7 eventually outputs the T-matrix for each size
c       write(7,"(i15,  5x,'NRANK')")nrank   -smg08
c       write(7,"(i15,  5x,'NMODE')")nm      -smg08
c       write(7,"(e15.7,5x,'OUTER SIZE PARAMETER')")conk   -smg08
c       write(7,"(e15.7,5x,'WAVELENGTH IN CM')")alamd    -smg08
c       write(7,"(e15.7,5x,'EQUIVALENT OUTER DIAMETER IN CM')")dpart  -smg08
c       write(7,"(e15.7,5x,'AXIS RATIO')")aovrb   -smg08
c       write(7,"(e15.7,5x,'SIGMA')")sigma   -smg08
c       write(7,"(e15.7,5x,'OUTER RE(DIELECTRIC CONSTANT)')") dcnr   -smg08
c       write(7,"(e15.7,5x,'OUTER IM(DIELECTRIC CONSTANT)')") dcni    -smg08

c       cmi = m values (only m= 1 is required for anginc= 0,but
c       m= 0,1,2,... is required for general incidence angle)
        if(nm.eq.1) then
          cmi(1) = 1.0
        else
          do 7  i = 1, nm
 7          cmi(i) = DBLE(i - 1)
        endif

c        print *,'G'

c       Preferable to start with kgauss=3 or more
        kgauss = 5
        call calenp

        do 140 i = 1, nsect
          epdeg(i) = rtd * epps(i)
  140   continue

c        print *,'H'

c       write(6,148) (epdeg(i),i=1,nsect)
c148    format(1x,'End points= ',8e12.4,/(1x,23x,8e12.4))
c       Jump up data set counter
        ijk = ijk + 1
        return

c        print *,'I'

  200   write(6, *) 'END OF INPUT DATA'
        print *, ' '
        print *, '@@@ Input values as read from file tmatrix2_pck.inp1:'

        write(6, 991)ptopdia, prhoice
  991   format(4x, 'Initial dia (cm)=    ', f6.2,
     >  '    Ice density (gm/cm**3)=', f6.3)

        write(6, 992) pzhgt, ptempt
  992   format(4x, 'Height (km MSL)=     ', f6.2,
     >  '    Scatterer temp (C)=   ', f5.1)

        write(6, 993) palamd, paovrb
  993   format(4x, 'Radar wavelength (cm)= ', f4.1,
     >  '    Scatterer axis ratio=   ', f5.3)
c       print *,'ptopdia, prhoice= ',ptopdia,prhoice
c       print *,'pzhgt, ptempt= ',pzhgt,ptempt
c       print *,'palamd,paovrb= ',palamd,paovrb
        print *, ' '
        print *, 'NOTE: Rename output file "tmatrix2_7"'
        print *, '      as required for nmueller program input.'
        print *, '      For water-coated hail, the required nmueller'
        print *, '      input file name is:'
        print *, '      "hail_wet"'
        print *, ' '
        write(6, *) 'PROGRAM COMPLETE: tmatrix2b_pck
     >  (ice dielectric=die_funct)'
c       close (7)
       close (8)
       close (22)
c        print *,'J'
c        stop
       ijk = ijk + 1
      end
c-----------------------------------------------------------------
      subroutine gener(xxx, t, ith, js)

      IMPLICIT DOUBLE PRECISION (A-H)
      IMPLICIT DOUBLE PRECISION (O-Z)
      IMPLICIT INTEGER (I-N)

        double precision dtr, rtd, cpi
        common dtr, rtd, cpi

        double precision theta, sinth, costh
        common/thtcom/theta,sinth,costh

        double precision a(40, 40, 2), b(40, 40, 2), c(40, 40, 2)
        double precision d(40, 40, 2), x(40, 40, 2), y(40, 40, 2)
        double precision cmxnrm(40)
        integer nrank,nranki
        common /mtxcom/ nrank, nranki,a, b, c, d, x, y, cmxnrm

        double precision dcnr, dcni, ckprr, ckpri, ckr, dckr, conk
        double precision aovrb, sigma
        integer ib
        common /bdycom/ dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1                  aovrb, sigma, ib

        double precision pnmllg(41), bsslsp(41, 31, 3)
        double precision cneumn(41, 31, 3)
        double precision bslkpr(41, 31, 3), bslkpi(41, 31, 3)
        double precision cneumr(41, 31, 3), cneumi(41, 31, 3)
        common /fnccom/ pnmllg, bsslsp, cneumn, bslkpr, bslkpi,
     1                  cneumr, cneumi

        double precision bkpr2(41, 31, 3), bkpi2(41, 31, 3)
        double precision cnpr2(41, 31, 3), cnpi2(41, 31, 3)
        double precision bkprr2(41, 31, 3), bkpii2(41,31,3)
        common /bringi/ bkpr2, bkpi2, cnpr2, cnpi2, bkprr2, bkpii2

        double precision dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2
        common /bdycm2/ dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2

        double precision crow, crowm, crow1, ccolm, ccol1, crij, crssij
        integer irow, irow1, icol,icol1
        common /rowcol/ irow, irow1, crow, crowm, crow1, icol, icol1,
     1                  ccol, ccolm, ccol1, crij, crssij

        double precision cmi(10), cmv, cm2, twm, prodm
        integer nm, kmv
        common /cmvcom/ nm, kmv, cmi, cmv, cm2, twm, prodm

        double precision qem, b89
        common /variam/ qem, b89

        double precision rhankl(40, 2), rbskpr(40, 2), rbessl(40)
        double precision hanklr(41), hankli(41), rhskpr(40, 2)
        double precision rbess2(40, 2), hankr2(41), hanki2(41)
        double precision rbskp2(40, 2), rhank2(40, 2)
        double precision t(4, 6, 2)

c        print *,'START gener'

        xxxr = xxx * rtd

        do 1101 i = 1, 4
          do 1101 j = 1, 6
            do 1101 k = 1, 2
              t(i, j, k) = 0.0
 1101   continue

        sqreal = crootr(dcnr, dcni)
        sqimag = crooti(dcnr, dcni)
        qsreal = cddvdr(1.D0, 0.D0, sqreal, sqimag)
        qsimag = cddvdi(1.D0, 0.D0, sqreal, sqimag)

c       Calculate sqrt(eps-2)
c
        sqrs2 = crootr(dcnr2, dcni2)
        sqis2 = crooti(dcnr2, dcni2)

c       Calculate eps-2/eps-1

        rat12r = cddvdr(dcnr2, dcni2, dcnr, dcni)
        rat12i = cddvdi(dcnr2, dcni2, dcnr, dcni)

c       calculate sqrt(eps-2/eps-1)
        sqrea2 = crootr(rat12r, rat12i)
        sqima2 = crooti(rat12r, rat12i)

c       calculate sqrt(eps-1/eps-2)
        qsrea2 = cddvdr(1.D0, 0.D0, sqrea2, sqima2)
        qsima2 = cddvdi(1.D0, 0.D0, sqrea2, sqima2)
        theta = xxx
        costh = dcos(theta)
        sinth = dsin(theta)
        srmsin = sinth
        tempp = cpi - theta

        if (dabs(tempp).lt.1.0d-8) sinth = 0.0

c       Generate the Legendre polynomials.
        call genlgp

c       Evaluate kr and its derivative as a function of theta.
  348   call genkr

c       Generate all necessary Bessel and Neumann function and their ratios.
        if((irow.eq.1).and.(icol.eq.1)) call genbsl(ith, js)

c       Calculate k*r for outer surface.
        ckprr = sqreal * ckr
        ckpri = sqimag * ckr
        iswt = 1

        if((irow.eq.1).and.(icol.eq.1)) call genbkr(ckprr,
     1                                              ckpri,
     2                                              iswt,
     3                                              ith,
     4                                              js)

        do 347 k = 1, nranki
          hanklr(k) = bslkpr(k, ith, js) - cneumi(k, ith, js)
          hankli(k) = bslkpi(k, ith, js) + cneumr(k, ith, js)
  347   continue

        do 350 k = 1, nrank
          rbessl(k) = bsslsp(k, ith, js) / bsslsp(k + 1, ith, js)
          rhankl(k, 1)=cddvdr(bsslsp(k, ith, js),
     1                          cneumn(k, ith, js),
     2                          bsslsp(k + 1, ith, js),
     3                          cneumn(k + 1, ith, js))

          rhankl(k,2) = cddvdi(bsslsp(k,ith,js),cneumn(k,ith,js),
     1    bsslsp(k+1,ith,js),cneumn(k+1,ith,js))
          rbskpr(k,1) = cddvdr(bslkpr(k,ith,js),bslkpi(k,ith,js),
     1    bslkpr(k+1,ith,js),bslkpi(k+1,ith,js))
          rbskpr(k,2) = cddvdi(bslkpr(k,ith,js),bslkpi(k,ith,js),
     1    bslkpr(k+1,ith,js),bslkpi(k+1,ith,js))
          bkr = cdmpyr(sqreal,sqimag,rbskpr(k,1),rbskpr(k,2))
          bki = cdmpyi(sqreal,sqimag,rbskpr(k,1),rbskpr(k,2))
          rbskpr(k,1) = bkr
          rbskpr(k,2) = bki
          tempnr=bslkpr(k,ith,js)-cneumi(k,ith,js)
          tempni=bslkpi(k,ith,js)+cneumr(k,ith,js)
          tempdr=bslkpr(k+1,ith,js)-cneumi(k+1,ith,js)
          tempdi=bslkpi(k+1,ith,js)+cneumr(k+1,ith,js)
          rhskpr(k,1)=cddvdr(tempnr,tempni,tempdr,tempdi)
          rhskpr(k,2)=cddvdi(tempnr,tempni,tempdr,tempdi)
          hkr=cdmpyr(sqreal,sqimag,rhskpr(k,1),rhskpr(k,2))
          hki=cdmpyi(sqreal,sqimag,rhskpr(k,1),rhskpr(k,2))
          rhskpr(k,1)=hkr
          rhskpr(k,2)=hki
  350   continue

      call genkr2

c      print *,'gener 1'

c     calculate k*r2 for inner surface.
      ckprr2=ckr2*sqreal
      ckpri2=ckr2*sqimag
c     calculate derivative of k*r2
      dkpr2=dckr2*sqreal
      dkpi2=dckr2*sqimag
      iswt=2
      if((irow.eq.1).and.(icol.eq.1))call genbkr(ckprr2,ckpri2,iswt,
     1ith,js)

c        print *,'gener 2'

c     calculate k*r2 for inner surface.
      ckhr2=ckr2*sqrs2
      ckhi2=ckr2*sqis2
      iswt=3
      if((irow.eq.1).and.(icol.eq.1)) call genbkr(ckhr2,ckhi2,iswt,
     1ith,js)

c        print *,'gener 3'

      do 349 k=1,nranki
      hankr2(k)=bkpr2(k,ith,js)-cnpi2(k,ith,js)
      hanki2(k)=bkpi2(k,ith,js)+cnpr2(k,ith,js)
  349   continue

c        print *,'gener 4'

      do 351 k=1,nrank
      rbess2(k,1)=cddvdr(bkpr2(k,ith,js),bkpi2(k,ith,js),
     1bkpr2(k+1,ith,js),bkpi2(k+1,ith,js))
      rbess2(k,2)=cddvdi(bkpr2(k,ith,js),bkpi2(k,ith,js),
     1bkpr2(k+1,ith,js),bkpi2(k+1,ith,js))
      tempnr=bkpr2(k,ith,js)-cnpi2(k,ith,js)
      tempni=bkpi2(k,ith,js)+cnpr2(k,ith,js)
      tempdr=bkpr2(k+1,ith,js)-cnpi2(k+1,ith,js)
      tempdi=bkpi2(k+1,ith,js)+cnpr2(k+1,ith,js)
      rhank2(k,1)=cddvdr(tempnr,tempni,tempdr,tempdi)
      rhank2(k,2)=cddvdi(tempnr,tempni,tempdr,tempdi)
  351   continue

c        print *,'gener 5'

      do 352 k=1,nrank
      rbskp2(k,1)=cddvdr(bkprr2(k,ith,js),bkpii2(k,ith,js),
     1bkprr2(k+1,ith,js),bkpii2(k+1,ith,js))
      rbskp2(k,2)=cddvdi(bkprr2(k,ith,js),bkpii2(k,ith,js),
     1bkprr2(k+1,ith,js),bkpii2(k+1,ith,js))
      bkr2=cdmpyr(sqrea2,sqima2,rbskp2(k,1),rbskp2(k,2))
      bki2=cdmpyi(sqrea2,sqima2,rbskp2(k,1),rbskp2(k,2))
      rbskp2(k,1)=bkr2
      rbskp2(k,2)=bki2
  352   continue

c        print *,'gener 6'

      br = rbessl(irow)
      hr = rhankl(irow,1)
      hi = rhankl(irow,2)
      hr2=rhank2(irow,1)
      hi2=rhank2(irow,2)
      br2=rbess2(irow,1)
      bi2=rbess2(irow,2)

c      print *,'gener 7'

c     calculate frequently used variable combinations for use in a,b,c,
c     d matrices.
      crij = crow+ccol
      crssij = crow*ccol
      cmcrco = cm2-qem*crssij*costh**2
      pnr0c0 = pnmllg(irow)*pnmllg(icol)
      pnr0c1 = pnmllg(irow)*pnmllg(icol+1)
      pnr1c0 = pnmllg(irow+1)*pnmllg(icol)
      pnr1c1 = pnmllg(irow+1)*pnmllg(icol+1)
      b1a = crow*costh*pnr1c1-crowm*pnr0c1
      b1b = ccol*costh*pnr1c1-ccolm*pnr1c0
      bkr = rbskpr(icol,1)
      bki = rbskpr(icol,2)
      hkr=rhskpr(icol,1)
      hki=rhskpr(icol,2)
      hbkmlr=cdmpyr(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),
     1bslkpr(icol+1,ith,js),bslkpi(icol+1,ith,js))
      hbkmli=cdmpyi(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),
     1bslkpr(icol+1,ith,js),bslkpi(icol+1,ith,js))
      bbkmlr = bsslsp(irow+1,ith,js)*bslkpr(icol+1,ith,js)
      bbkmli = bsslsp(irow+1,ith,js)*bslkpi(icol+1,ith,js)
      hhkmlr=cdmpyr(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),
     1hanklr(icol+1),hankli(icol+1))
      hhkmli=cdmpyi(bsslsp(irow+1,ith,js),cneumn(irow+1,ith,js),
     1hanklr(icol+1),hankli(icol+1))
      bhkmlr=bsslsp(irow+1,ith,js)*hanklr(icol+1)
      bhkmli=bsslsp(irow+1,ith,js)*hankli(icol+1)
      hepsr = cdmpyr(qsreal,qsimag,hbkmlr,hbkmli)
      hepsi = cdmpyi(qsreal,qsimag,hbkmlr,hbkmli)
      bepsr = cdmpyr(qsreal,qsimag,bbkmlr,bbkmli)
      bepsi = cdmpyi(qsreal,qsimag,bbkmlr,bbkmli)
      hhepsr=cdmpyr(qsreal,qsimag,hhkmlr,hhkmli)
      hhepsi=cdmpyi(qsreal,qsimag,hhkmlr,hhkmli)
      bbepsr=cdmpyr(qsreal,qsimag,bhkmlr,bhkmli)
      bbepsi=cdmpyi(qsreal,qsimag,bhkmlr,bhkmli)

c      print *,'gener 8'

c     calculate frequently used variable combinations for use in
c     x,y matrices.
      bkr2=rbskp2(icol,1)
      bki2=rbskp2(icol,2)
      hb2mlr=cdmpyr(hankr2(irow+1),hanki2(irow+1),
     1bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
      hb2mli=cdmpyi(hankr2(irow+1),hanki2(irow+1),
     1bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
      bb2mlr=cdmpyr(bkpr2(irow+1,ith,js),bkpi2(irow+1,ith,js),
     1bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
      bb2mli=cdmpyi(bkpr2(irow+1,ith,js),bkpi2(irow+1,ith,js),
     1bkprr2(icol+1,ith,js),bkpii2(icol+1,ith,js))
 902   format(2x,4(2x,i3),4(4x,e10.3))
       ftheta=theta*rtd

c       print *,'gener 9'

c      write(6,902) js,ith,irow,icol,ftheta,bkpr2(irow+1,ith,js),
c    1 bkpi2(irow+1,ith,js),bsslsp(irow+1,ith,js)
      heps2r=cdmpyr(qsrea2,qsima2,hb2mlr,hb2mli)
      heps2i=cdmpyi(qsrea2,qsima2,hb2mlr,hb2mli)
      beps2r=cdmpyr(qsrea2,qsima2,bb2mlr,bb2mli)
      beps2i=cdmpyi(qsrea2,qsima2,bb2mlr,bb2mli)
      if(ib.eq.9) go to 380

c      print *,'gener 9'

c     if ib = 8 (mirror symmetry body), i=l=0 if irow and icol are both
c     odd or both even, j=k=0 if irow and icol are odd,even or even,odd.
      if((irow+icol).eq.((irow+icol)/2)*2) go to 392

c      print *,'gener 10'

c     test for m=0 (if m=0 the i and l submatrices are zero).
  380 if(kmv.eq.0) go to 390

c        print *,'gener 11'

c     fill out elements for equivalent i-submatrix position.
      b1 = b1a+b1b
      htbkr = cdmpyr(hr,hi,bkr,bki)
      htbki = cdmpyi(hr,hi,bkr,bki)
      btbkr = br*bkr
      btbki = br*bki
      tempp=(crow*crow1*bkr+ccol*ccol1*hr-crssij*(crij+2.0  )/ckr
     1)*dckr*sinth
      sumar=tempp*pnr1c1
      sumr=(ckr*(1.0  +htbkr)-ccol*hr-crow*bkr+crssij/ckr)*b1*ckr+sumar
      sumai = pnr1c1*(crow*crow1*bki+ccol*ccol1*hi)*dckr*sinth
      sumi = (ckr*htbki-ccol*hi-crow*bki)*b1*ckr+sumai
      t(1,1,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,hbkmlr,hbkmli)
      t(1,1,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,hbkmlr,hbkmli)
      sumbr = pnr1c1*(crow*crow1*bkr+ccol*ccol1*br-crssij*(crij+2.0d0)/c
     1kr)*dckr*sinth
      sumr=(ckr*(1.0  +btbkr)-ccol*br-crow*bkr+crssij/ckr)*b1*ckr+sumbr
      sumbi = pnr1c1*crow*crow1*bki*dckr*sinth
      sumi = (ckr*btbki-crow*bki)*b1*ckr+sumbi
      t(1,2,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,bbkmlr,bbkmli)
      t(1,2,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,bbkmlr,bbkmli)
      bthkr=br*hkr
      bthki=br*hki
      hthkr=cdmpyr(hr,hi,hkr,hki)
      hthki=cdmpyi(hr,hi,hkr,hki)
      sumcr=pnr1c1*(crow*crow1*hkr+ccol*ccol1*br-crssij*(crij+2.0d0)/ckr
     1)*dckr*sinth
      sumr=(ckr*(1.0  +bthkr)-ccol*br-crow*hkr+crssij/ckr)*b1*ckr+sumcr
      sumci=pnr1c1*(crow*crow1*hki)*dckr*sinth
      sumi=(ckr*bthki-crow*hki)*b1*ckr+sumci
      t(1,3,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,bhkmlr,bhkmli)
      t(1,3,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,bhkmlr,bhkmli)
      sumdr=pnr1c1*(crow*crow1*hkr+ccol*ccol1*hr-crssij*(crij+2.0d0)/ckr
     1)*dckr*sinth
      sumr=(ckr*(1.0  +hthkr)-ccol*hr-crow*hkr+crssij/ckr)*b1*ckr+sumdr
      sumdi=pnr1c1*(crow*crow1*hki+ccol*ccol1*hi)*dckr*sinth
      sumi=(ckr*hthki-ccol*hi-crow*hki)*b1*ckr+sumdi
      t(1,4,1)=b89*cmv*srmsin*cdmpyr(sumr,sumi,hhkmlr,hhkmli)
      t(1,4,2)=b89*cmv*srmsin*cdmpyi(sumr,sumi,hhkmlr,hhkmli)

c      print *,'gener 12'

c     place q2(out,re)in x-matrix and q2(re,re) in y-matrix.
      htbkr2=cdmpyr(hr2,hi2,bkr2,bki2)
      htbki2=cdmpyi(hr2,hi2,bkr2,bki2)
      btbkr2=cdmpyr(br2,bi2,bkr2,bki2)
      btbki2=cdmpyi(br2,bi2,bkr2,bki2)
      tempr=cddvdr(1.D0  ,0.D0  ,ckprr2,ckpri2)
      tempi=cddvdi(1.D0  ,0.D0  ,ckprr2,ckpri2)
      temr=1.0  +htbkr2
      temrz=cdmpyr(temr,htbki2,ckprr2,ckpri2)
      temiz=cdmpyi(temr,htbki2,ckprr2,ckpri2)
      sumxr=pnr1c1*(crow*crow1*bkr2+ccol*ccol1*hr2-crssij*(crij+2.0d0)*
     1tempr)*sinth
      sumxi=pnr1c1*(crow*crow1*bki2+ccol*ccol1*hi2-crssij*(crij+2.0d0)*
     1tempi)*sinth
      sumgr=cdmpyr(sumxr,sumxi,dkpr2,dkpi2)
      sumgi=cdmpyi(sumxr,sumxi,dkpr2,dkpi2)
      sumr=(temrz-ccol*hr2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*hi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumgr+sumr2
      six=sumgi+sumi2
      t(1,5,1)=b89*cmv*srmsin*cdmpyr(srx,six,hb2mlr,hb2mli)
      t(1,5,2)=b89*cmv*srmsin*cdmpyi(srx,six,hb2mlr,hb2mli)
      sumyr=pnr1c1*(crow*crow1*bkr2+ccol*ccol1*br2-crssij*(crij+2.0d0)*
     1tempr)*sinth
      sumyi=pnr1c1*(crow*crow1*bki2+ccol*ccol1*bi2-crssij*(crij+2.0d0)*
     1tempi)*sinth
      sumhr=cdmpyr(sumyr,sumyi,dkpr2,dkpi2)
      sumhi=cdmpyi(sumyr,sumyi,dkpr2,dkpi2)
      temr=1.0+btbkr2
      temrz=cdmpyr(temr,btbki2,ckprr2,ckpri2)
      temiz=cdmpyi(temr,btbki2,ckprr2,ckpri2)
      sumr=(temrz-ccol*br2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*bi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumhr+sumr2
      six=sumhi+sumi2
      t(1,6,1)=b89*cmv*srmsin*cdmpyr(srx,six,bb2mlr,bb2mli)
      t(1,6,2)=b89*cmv*srmsin*cdmpyi(srx,six,bb2mlr,bb2mli)

c      print *,'gener 13'

c     fill out elements for equivalent l-submatrix position.
      sumr=(ckr*(dcnr+htbkr)-ccol*hr-crow*bkr+crssij/ckr)*b1*ckr+sumar
      sumi=(ckr*(dcni+htbki)-ccol*hi-crow*bki)*b1*ckr+sumai
      t(2,1,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,hepsr,hepsi)
      t(2,1,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,hepsr,hepsi)
      sumr = (ckr*(dcnr+btbkr)-ccol*br-crow*bkr+crssij/ckr)*b1*ckr+sumbr
      sumi = (ckr*(dcni+btbki)-crow*bki)*b1*ckr+sumbi
      t(2,2,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,bepsr,bepsi)
      t(2,2,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,bepsr,bepsi)
      sumr=(ckr*(dcnr+bthkr)-ccol*br-crow*hkr+crssij/ckr)*b1*ckr+sumcr
      sumi=(ckr*(dcni+bthki)-crow*hki)*b1*ckr+sumci
      t(2,3,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,bbepsr,bbepsi)
      t(2,3,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,bbepsr,bbepsi)
      sumr=(ckr*(dcnr+hthkr)-ccol*hr-crow*hkr+crssij/ckr)*b1*ckr+sumdr
      sumi=(ckr*(dcni+hthki)-ccol*hi-crow*hki)*b1*ckr+sumdi
      t(2,4,1)=-b89*cmv*srmsin*cdmpyr(sumr,sumi,hhepsr,hhepsi)
      t(2,4,2)=-b89*cmv*srmsin*cdmpyi(sumr,sumi,hhepsr,hhepsi)
      temr=rat12r+htbkr2
      temi=rat12i+htbki2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz-ccol*hr2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*hi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumr2+sumgr
      six=sumi2+sumgi
      t(2,5,1)=-b89*cmv*srmsin*cdmpyr(srx,six,heps2r,heps2i)
      t(2,5,2)=-b89*cmv*srmsin*cdmpyi(srx,six,heps2r,heps2i)
      temr=rat12r+btbkr2
      temi=rat12i+btbki2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz-ccol*br2-crow*bkr2+crssij*tempr)*b1
      sumi=(temiz-ccol*bi2-crow*bki2+crssij*tempi)*b1
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumr2+sumhr
      six=sumi2+sumhi
      t(2,6,1)=-b89*cmv*srmsin*cdmpyr(srx,six,beps2r,beps2i)
      t(2,6,2)=-b89*cmv*srmsin*cdmpyi(srx,six,beps2r,beps2i)

c      print *,'gener 14'

390   if (ib.eq.8) go to 400
c     fill out elements for eqiivalent j-submatrix position.
  392 a12=cmcrco*pnr1c1+qem*(crow*ccolm*costh*pnr1c0+ccol*crowm*costh*pn
     1r0c1-crowm*ccolm*pnr0c0)
      b1a = ccol*ccol1*b1a
      b1b = crow*crow1*b1b
      b1 = (b1a-b1b)*sinth
      dd=-qem*dckr
      cr = cdmpyr(dcnr,dcni,hr,hi)
      ci = cdmpyi(dcnr,dcni,hr,hi)
      sumr=(ckr*(bkr-cr)+dcnr*crow-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(bki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,1,1)=b89*srmsin*cdmpyr(sumr,sumi,hepsr,hepsi)
      t(3,1,2)=b89*srmsin*cdmpyi(sumr,sumi,hepsr,hepsi)
      cr = br*dcnr
      ci = br*dcni
      sumr=(ckr*(bkr-cr)+dcnr*crow-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(bki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,2,1)=b89*srmsin*cdmpyr(sumr,sumi,bepsr,bepsi)
      t(3,2,2)=b89*srmsin*cdmpyi(sumr,sumi,bepsr,bepsi)
      sumr=(ckr*(hkr-cr)+crow*dcnr-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*dd
      sumi=(ckr*(hki-ci)+dcni*crow)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,3,1)=b89*srmsin*cdmpyr(sumr,sumi,bbepsr,bbepsi)
      t(3,3,2)=b89*srmsin*cdmpyi(sumr,sumi,bbepsr,bbepsi)
      crh=cdmpyr(dcnr,dcni,hr,hi)
      cih=cdmpyi(dcnr,dcni,hr,hi)
      sumr=(ckr*(hkr-crh)+crow*dcnr-ccol)*a12*ckr+(b1a-dcnr*b1b)*sinth*d
     1d
      sumi=(ckr*(hki-cih)+crow*dcni)*a12*ckr-(dcni*b1b)*sinth*dd
      t(3,4,1)=b89*srmsin*cdmpyr(sumr,sumi,hhepsr,hhepsi)
      t(3,4,2)=b89*srmsin*cdmpyi(sumr,sumi,hhepsr,hhepsi)
      d2r=-qem*dkpr2
      d2i=-qem*dkpi2
      cr2=cdmpyr(rat12r,rat12i,hr2,hi2)
      ci2=cdmpyi(rat12r,rat12i,hr2,hi2)
      temr=bkr2-cr2
      temi=bki2-ci2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+crow*rat12r-ccol)*a12
      sumi=(temiz+crow*rat12i)*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      sumxr=(b1a-rat12r*b1b)*sinth
      sumxi=-rat12i*b1b*sinth
      sumxr2=cdmpyr(sumxr,sumxi,d2r,d2i)
      sumxi2=cdmpyi(sumxr,sumxi,d2r,d2i)
      srx=sumr2+sumxr2
      six=sumi2+sumxi2
      t(3,5,1)=b89*srmsin*cdmpyr(srx,six,heps2r,heps2i)
      t(3,5,2)=b89*srmsin*cdmpyi(srx,six,heps2r,heps2i)
      cr2=cdmpyr(rat12r,rat12i,br2,bi2)
      ci2=cdmpyi(rat12r,rat12i,br2,bi2)
      temr=bkr2-cr2
      temi=bki2-ci2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+rat12r*crow-ccol)*a12
      sumi=(temiz+rat12i*crow)*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      sumxr=(b1a-rat12r*b1b)*sinth
      sumxi=-rat12i*b1b*sinth
      sumxr2=cdmpyr(sumxr,sumxi,d2r,d2i)
      sumxi2=cdmpyi(sumxr,sumxi,d2r,d2i)
      srx=sumr2+sumxr2
      six=sumi2+sumxi2
      t(3,6,1)=b89*srmsin*cdmpyr(srx,six,beps2r,beps2i)
      t(3,6,2)=b89*srmsin*cdmpyi(srx,six,beps2r,beps2i)

c      print *,'gener 15'

c     fill out elements for equivalent k-submatrix.
      sumr = (ckr*(bkr-hr)+crow-ccol)*a12*ckr+b1*dd
      sumi = (ckr*(bki-hi))*a12*ckr
      t(4,1,1)=b89*srmsin*cdmpyr(sumr,sumi,hbkmlr,hbkmli)
      t(4,1,2)=b89*srmsin*cdmpyi(sumr,sumi,hbkmlr,hbkmli)
      sumr = (ckr*(bkr-br)+crow-ccol)*a12*ckr+b1*dd
      sumi = bki*a12*ckr**2
      t(4,2,1)=b89*srmsin*cdmpyr(sumr,sumi,bbkmlr,bbkmli)
      t(4,2,2)=b89*srmsin*cdmpyi(sumr,sumi,bbkmlr,bbkmli)
      sumr=(ckr*(hkr-br)+crow-ccol)*a12*ckr+b1*dd
      sumi=(ckr*hki)*a12*ckr
      t(4,3,1)=b89*srmsin*cdmpyr(sumr,sumi,bhkmlr,bhkmli)
      t(4,3,2)=b89*srmsin*cdmpyi(sumr,sumi,bhkmlr,bhkmli)
      sumr=(ckr*(hkr-hr)+crow-ccol)*a12*ckr+b1*dd
      sumi=(ckr*(hki-hi))*a12*ckr
      t(4,4,1)=b89*srmsin*cdmpyr(sumr,sumi,hhkmlr,hhkmli)
      t(4,4,2)=b89*srmsin*cdmpyi(sumr,sumi,hhkmlr,hhkmli)
      temr=bkr2-hr2
      temi=bki2-hi2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+crow-ccol)*a12
      sumi=temiz*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      sumxr=b1*d2r
      sumxi=b1*d2i
      srx=sumr2+sumxr
      six=sumi2+sumxi
      t(4,5,1)=b89*srmsin*cdmpyr(srx,six,hb2mlr,hb2mli)
      t(4,5,2)=b89*srmsin*cdmpyi(srx,six,hb2mlr,hb2mli)
      temr=bkr2-br2
      temi=bki2-bi2
      temrz=cdmpyr(temr,temi,ckprr2,ckpri2)
      temiz=cdmpyi(temr,temi,ckprr2,ckpri2)
      sumr=(temrz+crow-ccol)*a12
      sumi=temiz*a12
      sumr2=cdmpyr(sumr,sumi,ckprr2,ckpri2)
      sumi2=cdmpyi(sumr,sumi,ckprr2,ckpri2)
      srx=sumr2+sumxr
      six=sumi2+sumxi
      t(4,6,1)=b89*srmsin*cdmpyr(srx,six,bb2mlr,bb2mli)
      t(4,6,2)=b89*srmsin*cdmpyi(srx,six,bb2mlr,bb2mli)

c      print *,'END gener'

400   return
      end
c-----------------------------------------------------------------
      subroutine quad(a,b,k,result,npts,ith,js)
c     Numerical integration using Gauss-Legendre method.

      IMPLICIT DOUBLE PRECISION (A-H)
      IMPLICIT DOUBLE PRECISION (O-Z)
      IMPLICIT INTEGER (I-N)

      double precision p,d1,d2,d3,d4,d5,d6,d7,d8
      common/gauss/kgauss
      dimension p(381),funct(48,127),fzero(4,6,2),acum(48),test(48),
     1t(4,6,2),result(48,8)
      dimension d1(54),d2(54),d3(54),d4(54),d5(54),d6(54),d7(54),d8(3)
      equivalence (p(1),d1(1)),(p(55),d2(1)),(p(109),d3(1)),(p(163),d4
     1(1)),(p(217),d5(1)),(p(271),d6(1)),(p(325),d7(1)),(p(379),d8(1))
      data d1/
     1 7.74596669241483d-01, 5.55555555555557d-01, 8.88888888888889d-01,
     2 2.68488089868333d-01, 9.60491268708019d-01, 1.04656226026467d-01,
     3 4.34243749346802d-01, 4.01397414775962d-01, 4.50916538658474d-01,
     4 1.34415255243784d-01, 5.16032829970798d-02, 2.00628529376989d-01,
     5 9.93831963212756d-01, 1.70017196299402d-02, 8.88459232872258d-01,
     6 9.29271953151245d-02, 6.21102946737228d-01, 1.71511909136392d-01,
     7 2.23386686428967d-01, 2.19156858401588d-01, 2.25510499798206d-01,
     8 6.72077542959908d-02, 2.58075980961766d-02, 1.00314278611795d-01,
     9 8.43456573932111d-03, 4.64628932617579d-02, 8.57559200499902d-02,
     + 1.09578421055925d-01, 9.99098124967666d-01, 2.54478079156187d-03,
     1 9.81531149553739d-01, 1.64460498543878d-02, 9.29654857429739d-01,
     2 3.59571033071293d-02, 8.36725938168868d-01, 5.69795094941234d-02,
     3 7.02496206491528d-01, 7.68796204990037d-02, 5.31319743644374d-01,
     4 9.36271099812647d-02, 3.31135393257977d-01, 1.05669893580235d-01,
     5 1.12488943133187d-01, 1.11956873020953d-01, 1.12755256720769d-01,
     6 3.36038771482077d-02, 1.29038001003512d-02, 5.01571393058995d-02,
     7 4.21763044155885d-03, 2.32314466399103d-02, 4.28779600250078d-02,
     8 5.47892105279628d-02, 1.26515655623007d-03, 8.22300795723591d-03/
      data d2/
     1 1.79785515681282d-02, 2.84897547458336d-02, 3.84398102494556d-02,
     2 4.68135549906281d-02, 5.28349467901166d-02, 5.59784365104763d-02,
     3 9.99872888120358d-01, 3.63221481845531d-04, 9.97206259372224d-01,
     4 2.57904979468569d-03, 9.88684757547428d-01, 6.11550682211726d-03,
     5 9.72182874748583d-01, 1.04982469096213d-02, 9.46342858373402d-01,
     6 1.54067504665595d-02, 9.10371156957005d-01, 2.05942339159128d-02,
     7 8.63907938193691d-01, 2.58696793272147d-02, 8.06940531950218d-01,
     8 3.10735511116880d-02, 7.39756044352696d-01, 3.60644327807826d-02,
     9 6.62909660024781d-01, 4.07155101169443d-02, 5.77195710052045d-01,
     + 4.49145316536321d-02, 4.83618026945841d-01, 4.85643304066732d-02,
     1 3.83359324198731d-01, 5.15832539520484d-02, 2.77749822021825d-01,
     2 5.39054993352661d-02, 1.68235251552208d-01, 5.54814043565595d-02,
     3 5.63443130465928d-02, 5.62776998312542d-02, 5.63776283603847d-02,
     4 1.68019385741038d-02, 6.45190005017574d-03, 2.50785696529497d-02,
     5 2.10881524572663d-03, 1.16157233199551d-02, 2.14389800125039d-02,
     6 2.73946052639814d-02, 6.32607319362634d-04, 4.11150397865470d-03,
     7 8.98927578406411d-03, 1.42448773729168d-02, 1.92199051247278d-02,
     8 2.34067774953141d-02, 2.64174733950583d-02, 2.79892182552381d-02/
      data d3/
     1 1.80739564445388d-04, 1.28952408261042d-03, 3.05775341017553d-03,
     2 5.24912345480885d-03, 7.70337523327974d-03, 1.02971169579564d-02,
     3 1.29348396636074d-02, 1.55367755558440d-02, 1.80322163903913d-02,
     4 2.03577550584721d-02, 2.24572658268161d-02, 2.42821652033366d-02,
     5 2.57916269760242d-02, 2.69527496676331d-02, 2.77407021782797d-02,
     6 2.81388499156271d-02, 9.99982430354891d-01, 5.05360952078625d-05,
     7 9.99598799671912d-01, 3.77746646326985d-04, 9.98316635318407d-01,
     8 9.38369848542380d-04, 9.95724104698407d-01, 1.68114286542147d-03,
     9 9.91495721178104d-01, 2.56876494379402d-03, 9.85371499598521d-01,
     + 3.57289278351730d-03, 9.77141514639705d-01, 4.67105037211432d-03,
     1 9.66637851558417d-01, 5.84344987583563d-03, 9.53730006425761d-01,
     2 7.07248999543356d-03, 9.38320397779592d-01, 8.34283875396818d-03,
     3 9.20340025470011d-01, 9.64117772970252d-03, 8.99744899776941d-01,
     4 1.09557333878379d-02, 8.76513414484705d-01, 1.22758305600827d-02,
     5 8.50644494768350d-01, 1.35915710097655d-02, 8.22156254364980d-01,
     6 1.48936416648152d-02, 7.91084933799848d-01, 1.61732187295777d-02,
     7 7.57483966380512d-01, 1.74219301594641d-02, 7.21423085370098d-01,
     8 1.86318482561388d-02, 6.82987431091078d-01, 1.97954950480975d-02/
      data d4/
     1 6.42276642509760d-01, 2.09058514458120d-02, 5.99403930242243d-01,
     2 2.19563663053178d-02, 5.54495132631931d-01, 2.29409642293877d-02,
     3 5.07687757533716d-01, 2.38540521060385d-02, 4.59130011989833d-01,
     4 2.46905247444876d-02, 4.08979821229888d-01, 2.54457699654648d-02,
     5 3.57403837831532d-01, 2.61156733767061d-02, 3.04576441556714d-01,
     6 2.66966229274503d-02, 2.50678730303482d-01, 2.71855132296248d-02,
     7 1.95897502711100d-01, 2.75797495664819d-02, 1.40424233152560d-01,
     8 2.78772514766137d-02, 8.44540400837110d-02, 2.80764557938172d-02,
     9 2.81846489497457d-02, 2.81763190330167d-02, 2.81888141801924d-02,
     + 8.40096928705192d-03, 3.22595002508787d-03, 1.25392848264749d-02,
     1 1.05440762286332d-03, 5.80786165997757d-03, 1.07194900062519d-02,
     2 1.36973026319907d-02, 3.16303660822264d-04, 2.05575198932735d-03,
     3 4.49463789203206d-03, 7.12243868645840d-03, 9.60995256236391d-03,
     4 1.17033887476570d-02, 1.32087366975291d-02, 1.39946091276191d-02,
     5 9.03727346587510d-05, 6.44762041305726d-04, 1.52887670508776d-03,
     6 2.62456172740443d-03, 3.85168761663987d-03, 5.14855847897819d-03,
     7 6.46741983180368d-03, 7.76838777792199d-03, 9.01610819519566d-03,
     8 1.01788775292361d-02, 1.12286329134080d-02, 1.21410826016683d-02/
      data d5/
     1 1.28958134880121d-02, 1.34763748338165d-02, 1.38703510891399d-02,
     2 1.40694249578135d-02, 2.51578703842806d-05, 1.88873264506505d-04,
     3 4.69184924247851d-04, 8.40571432710723d-04, 1.28438247189701d-03,
     4 1.78644639175865d-03, 2.33552518605716d-03, 2.92172493791781d-03,
     5 3.53624499771678d-03, 4.17141937698409d-03, 4.82058886485126d-03,
     6 5.47786669391895d-03, 6.13791528004137d-03, 6.79578550488277d-03,
     7 7.44682083240758d-03, 8.08660936478883d-03, 8.71096507973207d-03,
     8 9.31592412806942d-03, 9.89774752404876d-03, 1.04529257229060d-02,
     9 1.09781831526589d-02, 1.14704821146939d-02, 1.19270260530193d-02,
     + 1.23452623722438d-02, 1.27228849827324d-02, 1.30578366883530d-02,
     1 1.33483114637252d-02, 1.35927566148124d-02, 1.37898747832410d-02,
     2 1.39386257383068d-02, 1.40382278969086d-02, 1.40881595165083d-02,
     3 9.99997596379750d-01, 6.93793643241083d-06, 9.99943996207055d-01,
     4 5.32752936697805d-05, 9.99760490924434d-01, 1.35754910949228d-04,
     5 9.99380338025023d-01, 2.49212400482998d-04, 9.98745614468096d-01,
     6 3.89745284473282d-04, 9.97805354495956d-01, 5.54295314930373d-04,
     7 9.96514145914890d-01, 7.40282804244503d-04, 9.94831502800622d-01,
     8 9.45361516858527d-04, 9.92721344282788d-01, 1.16748411742996d-03/
      data d6/
     1 9.90151370400771d-01, 1.40490799565515d-03, 9.87092527954033d-01,
     2 1.65611272815445d-03, 9.83518657578632d-01, 1.91971297101387d-03,
     3 9.79406281670862d-01, 2.19440692536384d-03, 9.74734459752401d-01,
     4 2.47895822665757d-03, 9.69484659502459d-01, 2.77219576459345d-03,
     5 9.63640621569812d-01, 3.07301843470258d-03, 9.57188216109859d-01,
     6 3.38039799108691d-03, 9.50115297521293d-01, 3.69337791702565d-03,
     7 9.42411565191083d-01, 4.01106872407503d-03, 9.34068436157727d-01,
     8 4.33264096809299d-03, 9.25078932907077d-01, 4.65731729975685d-03,
     9 9.15437587155765d-01, 4.98436456476553d-03, 9.05140358813263d-01,
     + 5.31308660518706d-03, 8.94184568335557d-01, 5.64281810138445d-03,
     1 8.82568840247341d-01, 5.97291956550816d-03, 8.70293055548114d-01,
     2 6.30277344908575d-03, 8.57358310886234d-01, 6.63178124290190d-03,
     3 8.43766882672707d-01, 6.95936140939044d-03, 8.29522194637402d-01,
     4 7.28494798055382d-03, 8.14628787655138d-01, 7.60798966571904d-03,
     5 7.99092290960843d-01, 7.92794933429486d-03, 7.82919394118284d-01,
     6 8.24430376303287d-03, 7.66117819303759d-01, 8.55654356130769d-03,
     7 7.48696293616938d-01, 8.86417320948252d-03, 7.30664521242183d-01,
     8 9.16671116356077d-03, 7.12033155362253d-01, 9.46368999383007d-03/
      data d7/
     1 6.92813769779114d-01, 9.75465653631741d-03, 6.73018830230419d-01,
     2 1.00391720440569d-02, 6.52661665410019d-01, 1.03168123309476d-02,
     3 6.31756437711193d-01, 1.05871679048852d-02, 6.10318113715188d-01,
     4 1.08498440893373d-02, 5.88362434447664d-01, 1.11044611340069d-02,
     5 5.65905885423653d-01, 1.13506543159806d-02, 5.42965666498311d-01,
     6 1.15880740330440d-02, 5.19559661537457d-01, 1.18163858908302d-02,
     7 4.95706407918762d-01, 1.20352707852796d-02, 4.71425065871658d-01,
     8 1.22444249816120d-02, 4.46735387662029d-01, 1.24435601907140d-02,
     9 4.21657686626164d-01, 1.26324036435421d-02, 3.96212806057616d-01,
     + 1.28106981638774d-02, 3.70422087950079d-01, 1.29782022395374d-02,
     1 3.44307341599437d-01, 1.31346900919602d-02, 3.17890812068477d-01,
     2 1.32799517439305d-02, 2.91195148518247d-01, 1.34137930851101d-02,
     3 2.64243372410927d-01, 1.35360359349562d-02, 2.37058845589829d-01,
     4 1.36465181025713d-02, 2.09665238243181d-01, 1.37450934430019d-02,
     5 1.82086496759252d-01, 1.38316319095064d-02, 1.54346811481378d-01,
     6 1.39060196013255d-02, 1.26470584372302d-01, 1.39681588065169d-02,
     7 9.84823965981194d-02, 1.40179680394566d-02, 7.04069760428552d-02,
     8 1.40553820726499d-02, 4.22691647653637d-02, 1.40803519625536d-02/
      data d8/
     1 1.40938864107825d-02, 1.40928450691604d-02, 1.40944070900962d-02/

c        print *,'START quad'

      if(a.eq.b)go to 107
      sum=(b+a)/2.0
      diff=(b-a)/2.0
c     one point formula.
c     set up variable combinations for use in evaluation of integrands
       ith=ith+1
      call gener(sum,t,ith,js)
      do 1010 ii=1,4
      do 1020 jj=1,6
      do 1030 kk=1,2
c     jj=1,6 corresponds to matrices a thro d and x thro y. ii=1,4
c     corresponds to filling out eq. i,l,j,k positoons or sub-matrices
c     kk=1,2 corresponds to real and imaginary parts.
      fzero(ii,jj,kk)=t(ii,jj,kk)
      l=(ii-1)*12+(jj-1)*2+kk
      result(l,1)=2.0  *t(ii,jj,kk)*diff
1030  continue
1020  continue
1010  continue
      i=0
      iold=0
      inew=1
      k=2
      do 1040 n=1,48
      acum(n)=0.0
1040   continue
      go to 103
101   continue
      if(k.eq.kgauss) go to 105
      k=k+1
      do 1050 n=1,48
      acum(n)=0.0
1050  continue
c     contribution from function values already computed.
      do 102 j=1,iold
      i=i+1
      do 1060 n=1,48
      acum(n)=acum(n)+p(i)*funct(n,j)
1060  continue
102   continue
c     contribution from new values.
103   continue
      iold=iold+inew
      do 104 j=inew,iold
      i=i+1
      x=p(i)*diff
      temp1=sum+x
        ith=ith+1
      call gener(temp1,t,ith,js)
      do 1070 ii=1,4
      do 1080 jj=1,6
      do 1090 kk=1,2
      l=(ii-1)*12+(jj-1)*2+kk
c     l goes from 1 to 48.
      funct(l,j)=t(ii,jj,kk)
1090  continue
1080  continue
1070  continue
      temp2=sum-x
        ith=ith+1
      call gener(temp2,t,ith,js)
      do 1071 ii=1,4
      do 1081 jj=1,6
      do 1091 kk=1,2
      l=(ii-1)*12+(jj-1)*2+kk
      funct(l,j)=funct(l,j)+t(ii,jj,kk)
1091  continue
1081  continue
1071  continue
      i=i+1
      do 1100 n=1,48
      acum(n)=acum(n)+p(i)*funct(n,j)
1100  continue
104   continue
      inew=iold+1
      i=i+1
      do 1200 ii=1,4
      do 1300 jj=1,6
      do 1400 kk=1,2
      n=(ii-1)*12+(jj-1)*2+kk
      result(n,k)=(acum(n)+p(i)*fzero(ii,jj,kk))*diff
1400  continue
1300  continue
1200  continue
      go to 101
105   continue
c     normal termination.
106   continue
      npts=inew+iold
      return
c     trivial case
 107  continue
      k=2
      do 1600 m=1,48
      do 1700 n=1,2
      result(m,n)=0.0
1700  continue
1600  continue
      npts=0

c      print *,'END quad'

      return
      end
c-----------------------------------------------------------------
      subroutine genlgp

c     A routine to generate Legendre polynomials.
c     The index on the function is incremented by one.

      IMPLICIT DOUBLE PRECISION (A-H)
      IMPLICIT DOUBLE PRECISION (O-Z)
      IMPLICIT INTEGER (I-N)

      common dtr,rtd,cpi
      common /mtxcom/ nrank,nranki,a(40,40,2),b(40,40,2),c(40,40,2),
     1d(40,40,2),x(40,40,2),y(40,40,2),cmxnrm(40)
      common/fnccom/pnmllg(41),bsslsp(41,31,3),cneumn(41,31,3),
     1bslkpr(41,31,3),bslkpi(41,31,3),cneumr(41,31,3),cneumi(41,31,3)
      common /cmvcom/ nm,kmv,cmi(10),cmv,cm2,twm,prodm
      common /thtcom/ theta,sinth,costh

c      print *,'START genlgp'

      dtwm=twm+1.0
c     this is special case when theta equals cpiand m=0.
c     when theta equals cpi all integrands are 0 and any values can be
c     put in pnmllg(41).here we have put them equal to 0.
      if((sinth.eq.0.0  ).and.(kmv.eq.0)) go to 6
c     at this point theta lies strictly between 0 and cpi.
      if(theta)16,4,16
    4 if(kmv-1)6,12,6
    6 do 8 ilg = 1,nranki
      pnmllg(ilg)=0.0
    8 continue
      go to 88
12    pnmllg(1)=0.0
      pnmllg(2)=1.0
      pla=1.0
      go to 48
   16 if(kmv)20,20,40
c     the special case when m = 0.
20    pla=1.0/sinth
      plb = costh*pla
      pnmllg(1) = pla
      pnmllg(2) = plb
      ibeg = 3
      go to 60
c     general case for m not equal to 0.
   40 do 44 ilg = 1,kmv
      pnmllg(ilg)=0.0
   44 continue
      if((sinth.eq.0.0  ).and.(kmv.eq.1)) go to 1001
      pla = prodm*sinth**(kmv-1)
      go to 1002
1001  pla=0.0
1002  continue
      pnmllg(kmv+1) = pla
   48 plb = dtwm*costh*pla
      pnmllg(kmv+2) = plb
      ibeg = kmv+3
c     do recursion formula for all remaining legendre polynomials.
   60 cnmul = ibeg+ibeg-3
      cnm=2.0
      cnmm = dtwm
      do 80 ilgr = ibeg,nranki
      plc = (cnmul*costh*plb-cnmm*pla)/cnm
      pnmllg(ilgr) = plc
      pla = plb
      plb = plc
      cnmul=cnmul+2.0
      cnm=cnm+1.0
      cnmm=cnmm+1.0
   80 continue

c        print *,'END genlgp'

   88 return
      end
c-----------------------------------------------------------------
      double precision function cdabx(a, b)

        implicit none

        double precision a, b, e, f, g

        if(a) 4, 22, 4

    4   if(b) 8, 30, 8
    8   e = dmax1(a, b)
        f = dmin1(a, b)
        g = f / e
        cdabx = dabs(e) * dsqrt(1.0 + g * g)
        return

   22   if(b) 28, 26, 28
   26   cdabx = 0.0

        return

   28   cdabx = dabs(b)
        return

   30   cdabx = dabs(a)

        return
      end
c-----------------------------------------------------------------
      double precision function cdmpyr(a, b, c, d)

        implicit none

        double precision a, b, c, d

        cdmpyr = a * c - b * d

        return
      end
c-----------------------------------------------------------------
      double precision function cdmpyi(a, b, c, d)

        implicit none

        double precision a, b, c, d

        cdmpyi = b * c + a * d

        return
      end
c-----------------------------------------------------------------
      double precision function cddvdr(a, b, c, d)

        implicit none

        double precision a, b, c, d, e, f

        e = c * c + d * d
        f = a * c + b * d
        cddvdr = f / e

        return
      end
c-----------------------------------------------------------------
      double precision function cddvdi(a, b, c, d)

        implicit none

        double precision a, b, c, d, e, f

        e = c * c + d * d
        f = b * c - a * d
        cddvdi = f / e

        return
      end
c-----------------------------------------------------------------
      subroutine genbsl(ith,js)

c     Generate Bessel and Neumann functions for real arguments.
c     The index on the function is incremented by one.

      IMPLICIT DOUBLE PRECISION (A-H)
      IMPLICIT DOUBLE PRECISION (O-Z)
      IMPLICIT INTEGER (I-N)

      common dtr,rtd,cpi
      common /mtxcom/ nrank,nranki,a(40,40,2),b(40,40,2),c(40,40,2),
     1d(40,40,2),x(40,40,2),y(40,40,2),cmxnrm(40)
      common /thtcom/ theta,sinth,costh
      common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
      common/fnccom/pnmllg(41),bsslsp(41,31,3),cneumn(41,31,3),
     1bslkpr(41,31,3),bslkpi(41,31,3),cneumr(41,31,3),cneumi(41,31,3)
c     set up a loop to get 2 successive bessel functions
      nval=nrank-1
      pckr=ckr
      do 40 i=1,4
      call bessel(nval,pckr,answr,ierror)
      if(ierror) 20,20,32
20    ansa=answr
      nval=nval+1
      call bessel(nval,pckr,answr,ierror)
      if(ierror) 24,24,28
24    ansb=answr
      go to 60
28    nval=nval-1
32    nval=nval+nrank
40    continue
c     program unable to generate bessel functions
      write(6,1001)
1001  format(///,5x,'unable to generate bessel functions',//)
c     set up for proper recursion of the bessel functons
60    if(nval-nrank)100,100,64
64    iend=nval-nrank
      conn=2*(nval-1)+1.0
      do 72 ip=1,iend
      ansc=conn*ansa/pckr-ansb
      conn=conn-2.0
      ansb=ansa
      ansa=ansc
72     continue
c     program is ready to recurse downward into bessel function
100   bsslsp(nranki,ith,js)=ansb
      bsslsp(nranki-1,ith,js)=ansa
      conn=    (DBLE(nrank+nrank-1))
      ie=nranki-2
      je=ie
      do 120 jb=1,je
      ansc=conn*ansa/pckr-ansb
      bsslsp(ie,ith,js)=ansc
      ansb=ansa
      ansa=ansc
      ie=ie-1
      conn=conn-2.0
120   continue
c     generate neumann functions
      cskrx= dcos(pckr)/pckr
      snkrx= dsin(pckr)/pckr
      ckr2=pckr**2
      cmuln=3.0
      snsa=-cskrx
      snsb=-cskrx/pckr-snkrx
      cneumn(1,ith,js)=snsa
      cneumn(2,ith,js)=snsb
      do 280 i=3,nranki
      snsc=cmuln*snsb/pckr-snsa
      cneumn(i,ith,js)=snsc
      snsa=snsb
      snsb=snsc
      cmuln=cmuln+2.0
280   continue
c     perform wronskian test on orders 0 and 1 and orders nrank-1 and nrank.
      quanbt=dabs(ckr2*(bsslsp(2,ith,js)*cneumn(1,ith,js)-
     1bsslsp(1,ith,js)*cneumn(2,ith,js))-1.0)
      quannt=dabs(ckr2*(bsslsp(nranki,ith,js)*cneumn(nrank,ith,js)-
     1bsslsp(nrank,ith,js)*cneumn(nranki,ith,js))-1.0  )
      if(quanbt-1.0e-10)360,352,352
  352 thtprt = rtd*theta
      write(6,356) thtprt,pckrr,quanbt,quannt
356   format(/,10x,'theta=',f9.4,'kr=',f10.4,'bessel test=',e12.5,
     1'neumann test=',e12.5)
      go to 362
  360 if(quannt-1.0e-10)362,352,352
  362 return
      end
c-----------------------------------------------------------------
      subroutine bessel(norder, argmnt, answr, ierror)

        implicit none

        double precision argmnt, answr, acr, apr, ci, cn, cni, fact
        double precision prod, sum, topr, x
        integer norder, ierror, ifct, n, i

        ierror = 0
        n = norder
        x = argmnt
        cn = n
        sum = 1.0
        apr = 1.0
        topr = -0.5 * x * x
        ci = 1.0
        cni = DBLE(2 * n) + 3.0
        cni = DBLE(2 * n) + 3.0

        do 60 i = 1, 100
          acr = topr * apr / (ci * cni)
          sum = sum + acr
          if(dabs(acr / sum) - 1.0e-20) 100, 100, 40

   40     apr = acr
          ci = ci + 1.0
          cni = cni + 2.0
   60   continue

        ierror = 1
        go to 200

c       the series has converged
  100   prod = (DBLE(2 * n)) + 1.0
        fact = 1.0d0

        if(n) 160, 160, 120
  120   do 140 ifct = 1, n
          fact = fact * x / prod
          prod = prod - 2.0
  140   continue

  160   answr = fact * sum
  200   return
      end
c-----------------------------------------------------------------
      subroutine genbkr(xxr,xxi,iswt,ith,js)

c     Generate Bessel functions for complex arguments.
c     The index on the function is incremented by one.

      IMPLICIT DOUBLE PRECISION (A-H)
      IMPLICIT DOUBLE PRECISION (O-Z)
      IMPLICIT INTEGER (I-N)

      common dtr,rtd,cpi
      common /mtxcom/ nrank,nranki,a(40,40,2),b(40,40,2),c(40,40,2),
     1d(40,40,2),x(40,40,2),y(40,40,2),cmxnrm(40)
      common/fnccom/pnmllg(41),bsslsp(41,31,3),cneumn(41,31,3),
     1bslkpr(41,31,3),bslkpi(41,31,3),cneumr(41,31,3),cneumi(41,31,3)
      common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
      common /thtcom/ theta,sinth,costh
      common/bringi/bkpr2(41,31,3),bkpi2(41,31,3),cnpr2(41,31,3),
     1cnpi2(41,31,3),bkprr2(41,31,3),bkpii2(41,31,3)
      dimension bsr(41),bsi(41),cnr(41),cni(41)
      dimension rjr(301),rji(301)
      dimension pr(20),pi(20),ipntr(20)

c      print *,'START genbkr'

      alarge=1.0e30
      pckrr=xxr
      pckri=xxi
      rm = cdabx(pckrr,pckri)
      nval=50
      if((rm.gt.25.0  ).and.(rm.le.150.0  )) nval=2*rm
      if(rm.gt.150.0  ) nval=300
c
c        print *,'genbkr 1'

c     generate bessel functions.
c
      rjr(nval+1)=0.0
      rji(nval+1)=0.0
      rji(nval)=0.0
      rjr(nval)=1.0
      if(rm.gt.2.0  ) rjr(nval) = 1.0e-10
      if(rm.gt.10.0  ) rjr(nval) = 1.0e-20
      if(rm.gt.25.0  ) rjr(nval) = 1.0e-30
      ie=nval+2
      eposx=exp(pckri)
      enegx=exp(-pckri)
      ex1 = (enegx+eposx)/2.0
      ex2 = (enegx-eposx)/2.0
      csr =  dsin(pckrr)*ex1
      csi = - dcos(pckrr)*ex2
      ar = cddvdr(csr,csi,pckrr,pckri)
      ai = cddvdi(csr,csi,pckrr,pckri)
      k=0

c      print *,'genbkr 2'

      do 10 i=2,nval
      ij=ie-i
      f=DBLE(2*ij-1)
      cfr=cddvdr(f,0.D0,pckrr,pckri)
      cfi=cddvdi(f,0.D0,pckrr,pckri)
      rjr(ij-1) = cdmpyr(rjr(ij),rji(ij),cfr,cfi)-rjr(ij+1)
      rji(ij-1) = cdmpyi(rjr(ij),rji(ij),cfr,cfi)-rji(ij+1)
      temr= dabs(rjr(ij-1))
      temi= dabs(rji(ij-1))
      if (temr.le.alarge.and.temi.le.alarge)  go to 10
c
c******* bump up the pointer k and store the current value
c******* of rj at p(k), and its position at ipntr(k)
c******* (taking care of real and im. parts)
c
      ii=ij+1
      if (ii.le.nranki) go to 8010
c
c******* (ij+1).gt.nranki  so normalize without storing
c******* rj and its position
c
      tr=cddvdr(rjr(ij),rji(ij),rjr(ij-1),rji(ij-1))
      ti=cddvdi(rjr(ij),rji(ij),rjr(ij-1),rji(ij-1))
      rjr(ij)=tr
      rji(ij)=ti
      rjr(ij-1)=1.0
      rji(ij-1)=0.0
      go to 10
 8010 continue
c
c******* (ij+1).le.nranki  so normalize and store rj and its position
c
      k=k+1
      if (k.le.20) go to 7000
      write (6,7003)
 7003 format (' pr,pi,ipntr arrays are too small')
      stop
 7000 continue
      pr(k)=rjr(ij-1)
      pi(k)=rji(ij-1)
      ipntr(k)=ij+1
      rjr(ij-1)=1.0
      rji(ij-1)=0.0
      tr=cddvdr(rjr(ij),rji(ij),pr(k),pi(k))
      ti=cddvdi(rjr(ij),rji(ij),pr(k),pi(k))
      rjr(ij)=tr
      rji(ij)=ti
      ii=ipntr(k)-2
   10 continue

c        print *,'genbkr 3'

c
c******* backsolution
c
      pxr = cddvdr(ar,ai,rjr(1),rji(1))
      pxi = cddvdi(ar,ai,rjr(1),rji(1))
      do 20 i=1,nranki
      if (ipntr(k).eq.i) go to 5001
      bsr(i)=cdmpyr(rjr(i),rji(i),pxr,pxi)
      bsi(i)=cdmpyi(rjr(i),rji(i),pxr,pxi)
      go to 20
 5001 continue
c
c******* make correction for px  ...  divide   by p(k) (real and im)
c
      tr=cddvdr(pxr,pxi,pr(k),pi(k))
      ti=cddvdi(pxr,pxi,pr(k),pi(k))
      pxr=tr
      pxi=ti
      k=k-1
      bsr(i)=cdmpyr(rjr(i),rji(i),pxr,pxi)
      bsi(i)=cdmpyi(rjr(i),rji(i),pxr,pxi)
   20 continue

c        print *,'genbkr 4'

c     generate neumann functions for test.
c
      ccr=dcos(pckrr)*ex1
      cci=dsin(pckrr)*ex2
      cnr(1)=-cddvdr(ccr,cci,pckrr,pckri)
      cni(1)=-cddvdi(ccr,cci,pckrr,pckri)
      cnr(2)=cddvdr(cnr(1),cni(1),pckrr,pckri)-ar
      cni(2)=cddvdi(cnr(1),cni(1),pckrr,pckri)-ai

c      print *,'genbkr 5'

      do 30 i=3,nranki
      f=DBLE(2*i-3)
      cfr=cddvdr(f,0.D0,pckrr,pckri)
      cfi=cddvdi(f,0.D0,pckrr,pckri)
      cnr(i)=cdmpyr(cnr(i-1),cni(i-1),cfr,cfi)-cnr(i-2)
      cni(i)=cdmpyi(cnr(i-1),cni(i-1),cfr,cfi)-cni(i-2)
   30 continue

c        print *,'genbkr 6'

      if(iswt.ne.1) go to 101

c      print *,'genbkr 7'

      do 201 i=1,nranki
      bslkpr(i,ith,js)=bsr(i)
      bslkpi(i,ith,js)=bsi(i)
      cneumr(i,ith,js)=cnr(i)
      cneumi(i,ith,js)=cni(i)
  201    continue

c        print *,'genbkr 8'

      go to  8001

c        print *,'genbkr 9'

  101   if(iswt.ne.2) go to 102

c        print *,'genbkr 10'

      do 202 i=1,nranki
      bkpr2(i,ith,js)=bsr(i)
      bkpi2(i,ith,js)=bsi(i)
      cnpr2(i,ith,js)=cnr(i)
      cnpi2(i,ith,js)=cni(i)
  202   continue

c        print *,'genbkr 11'

      go to  8001

c        print *,'genbkr 12'

  102   do 203 i=1,nranki
c        print *,'102 - i =', i, ' ith = ', ith, ' js = ', js

      bkprr2(i,ith,js)=bsr(i)
      bkpii2(i,ith,js)=bsi(i)
  203   continue

c        print *,'genbkr 13'


c     Perform two tests on Bessel and Neumann functions.  First test is
c     most accurate for large arguments and the second is most accurate
c     for smaller arguments.  If either test is passed, functions are good.

c     For large arguments abs(bessel) should equal abs(neumann).
8001   continue
      quabt=cdabx(bsr(1),bsi(1) )/cdabx(cnr(1),cni(1))-1.0
      quant=cdabx(bsr(nranki),bsi(nranki))/cdabx(cnr(nranki),cni(nranki)
     1)-1.0
      if(( dabs(quabt).gt.1.0d-8).or.( dabs(quant).gt.1.0d-8)) go to 32

c      print *,'END genbkr'

      return

c      print *,'genbkr 14'

c     Perform Wronskian test if large argument test fails.
   32 pckr2r = pckrr**2-pckri**2
      pckr2i = 2.0  *pckrr*pckri

c      print *,'genbkr 15'

c     bessel test
      itest = 1
      t1r=cdmpyr(bsr(2),bsi(2),cnr(1),cni(1))
      t1i=cdmpyi(bsr(2),bsi(2),cnr(1),cni(1))
      t2r=cdmpyr(bsr(1),bsi(1),cnr(2),cni(2))
      t2i=cdmpyi(bsr(1),bsi(1),cnr(2),cni(2))
   35 t3r = cddvdr(t1r,t1i,t2r,t2i)-1.0
      t3i = cddvdi(t1r,t1i,t2r,t2i)
      t4r = cdmpyr(t2r,t2i,t3r,t3i)
      t4i = cdmpyi(t2r,t2i,t3r,t3i)
      t5r = cdmpyr(pckr2r,pckr2i,t4r,t4i)-1.0
      t5i = cdmpyi(pckr2r,pckr2i,t4r,t4i)
      if(itest.eq.2) go to 40

c      print *,'genbkr 16'

      quanbt = cdabx(t5r,t5i)

c      print *,'genbkr 17'

c     neumann test
      t1r=cdmpyr(bsr(nranki),bsi(nranki),cnr(nrank),cni(nrank))
      t1i=cdmpyi(bsr(nranki),bsi(nranki),cnr(nrank),cni(nrank))
      t2r=cdmpyr(bsr(nrank),bsi(nrank),cnr(nranki),cni(nranki))
      t2i=cdmpyi(bsr(nrank),bsi(nrank),cnr(nranki),cni(nranki))
      itest = 2
      go to 35

c      print *,'genbkr 18'

   40 quannt = cdabx(t5r,t5i)
      if((quanbt.gt.1.0e-08).or.(quannt.gt.1.0e-08)) go to 45

c      print *,'END genbkr'

      return

c      print *,'genbkr 19'

   45 thtprt = rtd*theta
      write(6,50) thtprt,pckrr,pckri,quabt,quant,quanbt,quannt
   50 format(//,2x,8h theta =,f9.4,6h, kr =,2f10.4,7h, abt =,f12.5,7h, a
     1nt =,e12.5,7h, nbt =,e12.5,7h, nnt =,e12.5)

c        print *,'END genbkr'

      return
      end
c-----------------------------------------------------------------
      double precision function crootr(a, b)

        implicit none

        double precision a, b, dmag, angle

        dmag = (a * a + b * b)**0.25
        angle = 0.5 * atan2(b, a)
        crootr = dmag * dcos(angle)
        return
      end
c-----------------------------------------------------------------
      double precision function crooti(a, b)

        implicit none

        double precision a, b, dmag, angle

        dmag = (a * a + b * b)**0.25
        angle = 0.5 * atan2(b, a)
        crooti = dmag * dsin(angle)
        return
      end
c-----------------------------------------------------------------
      subroutine prcssm(a,b,nr,nri)

c     A routine to solve the equation t = (a-inverse)*b  ( all matrices
c     are transposed) using Gauss-Jordan elimination.

      IMPLICIT DOUBLE PRECISION (A-H)
      IMPLICIT DOUBLE PRECISION (O-Z)
      IMPLICIT INTEGER (I-N)

      dimension a(40,40,2),b(40,40,2)
      dimension aijmax(2),arat(2)
      equivalence (l,fl),(k,fk)

c      write(*,*) 'In prcssm'
      n = 2*nr
c     start reduction of the a matrix.
      do 80 i = 1,n
c     search for the maximum element in the ith row of the a-matrix.
      aijmax(1) = a(i,1,1)
      aijmax(2) = a(i,1,2)
      jmax = 1
      do 10 j = 2,n
      if(cdabx(a(i,j,1),a(i,j,2)).le.cdabx(aijmax(1),aijmax(2))) goto 10
      aijmax(1) = a(i,j,1)
      aijmax(2) = a(i,j,2)
      jmax = j
   10 continue
c     if aijmax is zero ( as it will be for any row (or column) where the
c     index m is .gt. the index n, i.e., the legendre functions force those
c     matrix elements to zero),then the matrix is singular so solve the
c     reduced matrix (order = 2*(nrank-m)).
      if(cdabx(aijmax(1),aijmax(2)).gt.0.0  ) go to 20
      jmax = i
      go to 75
c     normalize the ith row by aijmax (jmax element of the ith row).
   20 do 30 j = 1,n
      t1 = a(i,j,1)
      t2 = a(i,j,2)
      a(i,j,1) = cddvdr(t1,t2,aijmax(1),aijmax(2))
      a(i,j,2) = cddvdi(t1,t2,aijmax(1),aijmax(2))
c     normalize the ith row of b.
      t1 = b(i,j,1)
      t2 = b(i,j,2)
      b(i,j,1) = cddvdr(t1,t2,aijmax(1),aijmax(2))
      b(i,j,2) = cddvdi(t1,t2,aijmax(1),aijmax(2))
   30 continue
c     use row transformations to get zeros above and below the jmax
c     element of the ith row of a.  apply same row transformations
c     to the b matrix.
      do 70 k = 1,n
      if(k.eq.i) go to 70
      arat(1) = -a(k,jmax,1)
      arat(2) = -a(k,jmax,2)
      do 50 j = 1,n
      if(cdabx(a(i,j,1),a(i,j,2)).le.0.0  ) go to 50
      a(k,j,1) = cdmpyr(arat(1),arat(2),a(i,j,1),a(i,j,2))+a(k,j,1)
      a(k,j,2) = cdmpyi(arat(1),arat(2),a(i,j,1),a(i,j,2))+a(k,j,2)
   50 continue
      a(k,jmax,1)=0.0
      a(k,jmax,2)=0.0
      do 60 j=1,n
      if(cdabx(b(i,j,1),b(i,j,2)).le.0.0  ) go to 60
      b(k,j,1) = cdmpyr(arat(1),arat(2),b(i,j,1),b(i,j,2))+b(k,j,1)
      b(k,j,2) = cdmpyi(arat(1),arat(2),b(i,j,1),b(i,j,2))+b(k,j,2)
   60 continue
   70 continue
c     store row counter (i) in top element of jmax column.  thus,
c     the top row of a will contain the location of the pivot
c     (unity) element of each column (after reduction).
   75 l = i
c     store the integer i in the top row of a.
      a(1,jmax,1) = fl
   80 continue
c      write(*,*) 'after line 80 prcssm'
c     the reduction of a is complete.  perform row interchanges
c     as indicated in the first row of a.
      do 120 i = 1,n
      k=i
c     put the integer value in a into fk.
c     write(*,*) n, '=after 120'
      if (i .eq. 1)then
         do 44 jj=1,34
c           write(*,*)a(1,jj,1)
 44         continue
            endif
   90 fk = a(1,k,1)
c      write(*,*) k,i,fk,a(1,k,1)
      if(k-i) 90,120,100
c     if k(1,i) is less than i, then that row has already been
c     involved in an interchange, and we use k(1,k) until we get
c     a value of k greater than i (corresponding to a row stored
c     below the ith row).
  100 do 110 j=1,n
      arat(1) = b(i,j,1)
      arat(2) = b(i,j,2)
      b(i,j,1) = b(k,j,1)
      b(i,j,2) = b(k,j,2)
      b(k,j,1) = arat(1)
      b(k,j,2) = arat(2)
  110 continue
  120 continue
c      write (*,*) 'back to main'
      return
      end
c-----------------------------------------------------------------
      subroutine addprc

c     A routine to obtain the scattered field coefficients and calculate
c     the differential scattering cross section in the azimuthal plane.

      IMPLICIT DOUBLE PRECISION (A-H)
      IMPLICIT DOUBLE PRECISION (O-Z)
      IMPLICIT INTEGER (I-N)

      complex ci,cim
      common dtr,rtd,cpi
      common/mtxcom/nrank,nranki,a(40,40,2),b(40,40,2),c(40,40,2),
     1d(40,40,2),tmat(40,40,2),y(40,40,2),cmxnrm(40)
      common/bdycm2/dcnr2,dcni2,ckr2,dckr2,conk2,aovrb2
      common/fnccom/pnmllg(41),bsslsp(41,31,3),cneumn(41,31,3),
     1bslkpr(41,31,3),bslkpi(41,31,3),cneumr(41,31,3),cneumi(41,31,3)
      common /bdycom/ dcnr,dcni,ckprr,ckpri,ckr,dckr,conk,aovrb,sigma,ib
      common /cmvcom/ nm,kmv,cmi(10),cmv,cm2,twm,prodm
      common /thtcom/ theta,sinth,costh
      common /uvccom/anginc,acans(181,2,2),uang(181),rtsfct,dltang,nuang
      common/scatt/summ1,summ2
      common/vivek2/twoa
c########################################################################
c Mod on 11/4/2002:
c Next two lines were the original code; compile failure due to more than
c one array declaration for arrays "a" and "b":
c     dimension zxold(181),zyold(181),ab1(40,2),ab2(40,2),fg1(40,2),fg2(
c    140,2),fgans(181,2,2),a(40,40,2), b(40,40,2)
c Next two lines started as duplicates of the originals above.
c For PCK test, the "a" and "b" declarations were removed:
      dimension zxold(181),zyold(181),ab1(40,2),ab2(40,2),fg1(40,2),fg2(
     140,2),fgans(181,2,2)
c########################################################################
      logical test
      data test/.true./
      ci = (0.0,1.0)
c     Generate the legendre functions for the incident angle.
      if(anginc) 15,5,15
5     costh=1.0
10    sinth=0.0
      theta=0.0
      go to 30
   15 if (anginc-180.0) 25,20,25
20    costh=-1.0
      go to 10
   25 theta = dtr*anginc
      sinth=dsin(theta)
      costh=dcos(theta)
   30 call genlgp

c     Generate the incident field coefficients -- ab1 = theta polarization
c     and ab2 = phi polarization.

      cn=0.0
      do 35 n=1,nrank
      np = n+nrank
      cn=cn+1.0
      n1 = n+1
      ci1r=real(ci**n)
      ci1i=aimag(ci**n)
      ci2r=real(ci**n1)
      ci2i=aimag(ci**n1)
      p1 = cn*costh*pnmllg(n1)-(cn+cmv)*pnmllg(n)
      p2 = cmv*pnmllg(n1)
      ab1(n,1) = -ci1r*p2
      ab1(n,2) = -ci1i*p2
      ab1(np,1) = ci2r*p1
      ab1(np,2) = ci2i*p1
      ab2(n,1) = ci1r*p1
      ab2(n,2) = ci1i*p1
      ab2(np,1) = -ci2r*p2
      ab2(np,2) = -ci2i*p2
35    continue

c     The scattered field coefficients = the row vector of incident field
c     coefficients times the t-transposed matrix.

      nr2 = 2*nrank
      do 45 j = 1,nr2
      s1r=0.0
      s1i=0.0
      s2r=0.0
      s2i=0.0
      do 40 i = 1,nr2
      s1r = s1r+cdmpyr(ab1(i,1),ab1(i,2),tmat(i,j,1),tmat(i,j,2))
      s1i = s1i+cdmpyi(ab1(i,1),ab1(i,2),tmat(i,j,1),tmat(i,j,2))
      s2r = s2r+cdmpyr(ab2(i,1),ab2(i,2),tmat(i,j,1),tmat(i,j,2))
      s2i = s2i+cdmpyi(ab2(i,1),ab2(i,2),tmat(i,j,1),tmat(i,j,2))
   40 continue
      fg1(j,1) = s1r
      fg1(j,2) = s1i
      fg2(j,1) = s2r
      fg2(j,2) = s2i
   45 continue

c     Calculate scattering cossections normalized for parallel and
c     perpendicular polz
      sum1=0.0
      sum2=0.0
      do 1001 i=1,nrank
      ii=i+nrank
      temp1=fg1(i,1)**2+fg1(i,2)**2+fg1(ii,1)**2+fg1(ii,2)**2
      temp2=fg2(i,1)**2+fg2(i,2)**2+fg2(ii,1)**2+fg2(ii,2)**2
      sum1=sum1+temp1/cmxnrm(i)
      sum2=sum2+temp2/cmxnrm(i)
1001  continue
c     normalize scattering crossections.
      sum1=(rtsfct*2.0/conk)*sum1
      sum2=(rtsfct*2.0/conk)*sum2
c         normalize w.r.t. eq. spherical dia.
          cnorm=aovrb**(-2./3.)
          sum1=sum1*cnorm
          sum2=sum2*cnorm
c     accumulate results for each m value
      summ1=sum1+summ1
      summ2=sum2+summ2

c     Evaluate the scattered field at each scattering angle.

      do 170 iu = 1,nuang
c     Generate the Legendre multipliers.
      if(uang(iu)) 95,85,95
85    costh=1.0
90    sinth=0.0
      theta=0.0
      go to 110
   95 if(uang(iu)-180.0  ) 105,100,105
100   costh=-1.0
      go to 90
  105 theta = dtr*uang(iu)
      sinth=dsin(theta)
      costh=dcos(theta)
  110 call genlgp
      fgans(iu,1,1)=0.0
      fgans(iu,1,2)=0.0
      fgans(iu,2,1)=0.0
      fgans(iu,2,2)=0.0
      cn=0.0
      do 160 n = 1,nrank
      np = n+nrank
      n1 = n+1
      cn=cn+1.0
      p1 = cn*costh*pnmllg(n1)-(cn+cmv)*pnmllg(n)
      p2 = cmv*pnmllg(n1)
      cim = (-ci)**n1
      cir=real(cim)
      cii=aimag(cim)
      f1r = fg1(n,1)*p2
      f1i = fg1(n,2)*p2
      g1r = -fg1(np,2)*p1
      g1i = fg1(np,1)*p1
      fgans(iu,1,1) = fgans(iu,1,1)+cdmpyr(cir,cii,f1r+g1r,f1i+g1i)/cmxn
     1rm(n)
      fgans(iu,1,2) = fgans(iu,1,2)+cdmpyi(cir,cii,f1r+g1r,f1i+g1i)/cmxn
     1rm(n)
      f2r = fg2(n,1)*p1
      f2i = fg2(n,2)*p1
      g2r = -fg2(np,2)*p2
      g2i = fg2(np,1)*p2
      fgans(iu,2,1) = fgans(iu,2,1)-cdmpyr(cir,cii,f2r+g2r,f2i+g2i)/cmxn
     1rm(n)
      fgans(iu,2,2) = fgans(iu,2,2)-cdmpyi(cir,cii,f2r+g2r,f2i+g2i)/cmxn
     1rm(n)
  160 continue

c     The normalized diff.scat.cross sect. is given by ((8/ka)*fgans)**2
c     scale fgans to calculate diff. scat. cross sect. (rtsfct = 8/ka)

      fgans(iu,1,1) = rtsfct*fgans(iu,1,1)
      fgans(iu,1,2) = rtsfct*fgans(iu,1,2)
      fgans(iu,2,1) = rtsfct*fgans(iu,2,1)
      fgans(iu,2,2) = rtsfct*fgans(iu,2,2)
  170 continue

c     Accumulate the results for each m value.

c      write(6,175) kmv,anginc
c 175  format(30x,'ACCUMULATED SUMS FOR m= ',i3,
c     1 '  ANGLE OF INCIDENCE= ',f6.2,' DEGREES')
c      write(6,*)

      do 172 iup = 1,nuang
        acans(iup,1,1) = acans(iup,1,1)+fgans(iup,1,1)
        acans(iup,1,2) = acans(iup,1,2)+fgans(iup,1,2)
        acans(iup,2,1) = acans(iup,2,1)+fgans(iup,2,1)
        acans(iup,2,2) = acans(iup,2,2)+fgans(iup,2,2)
c        write(6,202) uang(iup),acans(iup,1,1),acans(iup,1,2)
c        write(6,203) uang(iup),acans(iup,2,1),acans(iup,2,2)
 172  continue
c 202  format(1x,'Parallel polz ',f7.3,2x,2e15.7)
c 203  format(1x,'Perpend  polz ',f7.3,2x,2e15.7)

c     Calculate the extinction crossections.
      extpp=acans(1,1,2)*rtsfct/4.0
      extper=acans(1,2,2)*rtsfct/4.0

c     Normalize wrt equivalent spherical diameter
      extpp=extpp*cnorm
      extper=extper*cnorm

c     Calculate forward and backward amplitude in far zone.
c     sigma equals 4.0/k
      forrp=sigma*acans(1,1,1)/rtsfct
      forip=sigma*acans(1,1,2)/rtsfct
      forpe=sigma*acans(1,2,1)/rtsfct
      foripe=sigma*acans(1,2,2)/rtsfct
      borrp=sigma*acans(nuang,1,1)/rtsfct
      borip=sigma*acans(nuang,1,2)/rtsfct
      borrpe=sigma*acans(nuang,2,1)/rtsfct
      boripe=sigma*acans(nuang,2,2)/rtsfct

c     Calculate normalized radar crossections for both polarizations
      xhor=acans(nuang,1,1)**2+acans(nuang,1,2)**2
      yver=acans(nuang,2,1)**2+acans(nuang,2,2)**2

c     Normalize wrt equivalent spherical diameter
      xhor= xhor*cnorm
      yver= yver*cnorm

c     Print the scattering results

c      write(6,2001)
c      write(6,2002) summ1,extpp,xhor,forrp,forip,borrp,borip
c      write(6,2003) summ2,extper,yver,forpe,foripe,borrpe,boripe
c      write(6,*)
c2001  format(18x,'SCAT',11x,'EXT',12x,'RADAR',10x,'RE(FORW)',7x,
c     1'IM(FORW)',7x,'RE(BACK)',7x,'IM(BACK)')
c2002  format(1x,'Parallel polz ',7(3x,e12.5))
c2003  format(1x,'Perpend  polz ',7(3x,e12.5))

c     Unit 8 writes the forward and backward amplitudes for
c     parallel/perp polz

      if(kmv.eq.6) then
        print *,'writing'
        write(8,2020) borrp,borip,borrpe,boripe,forrp,forip,forpe,foripe
c        write(8,2020) forpe,foripe,borrpe,boripe   -SMG 09 changed to write just backwards amps
2020    format(8(e15.7))
      endif

      return
      end

c-----------------------------------------------------------------
      subroutine genkr

        implicit none

        double precision dtr, rtd, cpi
        common dtr, rtd, cpi

c       calculate ckr and dckr as a function of theta for a oblate spheroid
        double precision dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1                      aovrb, sigma
        integer ib
        common /bdycom/ dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1                  aovrb, sigma, ib

        double precision theta, sinth, costh
        common /thtcom/ theta, sinth, costh

        double precision bovra, qb

        bovra = 1.0 / aovrb
        qb = 1.000 / dsqrt((bovra * costh)**2 + sinth**2)
        ckr = conk * qb
        dckr = conk * costh * sinth * (bovra**2 - 1.000) * qb**3
c       thet1=theta*rtd

        return
      end
c-----------------------------------------------------------------
      subroutine genkr2

        implicit none

        double precision dtr, rtd, cpi
        common dtr, rtd, cpi

c       calculate ckr2 and dckr2 for sph-cone-oblate
        double precision dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2
        common /bdycm2/ dcnr2, dcni2, ckr2, dckr2, conk2, aovrb2

        double precision theta, sinth, costh
        common /thtcom/ theta, sinth, costh

        double precision bovra2, qb2

        bovra2 = 1.0 / aovrb2

c       inner particle is an oblate spheroid
        qb2 = 1. / dsqrt((bovra2 * costh)**2 + sinth**2)
        ckr2 = conk2 * qb2
        dckr2 = conk2 * costh * sinth * (bovra2**2 - 1.0) * qb2**3

c       inner particle is a sphere
c       ckr2=conk2
c       dckr2=0.0
c       thet1=theta*rtd
c 10     format(2x,3(e15.7,4x))

        return
      end
c----------------------------------------------------------------
      subroutine calenp

        implicit none

        double precision dtr, rtd, cpi
        common dtr, rtd, cpi

        double precision epps(4)
        integer nsect
        common /endpnt/ epps, nsect

        double precision dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1                      aovrb, sigma
        integer ib
        common /bdycom/ dcnr, dcni, ckprr, ckpri, ckr, dckr, conk,
     1                  aovrb, sigma, ib

        epps(1) = 0.
        epps(2) = cpi / 2.

        return
      end

      subroutine pck_inp
c       Date: 12/11/2002
c       Reads program initialization values from an input file.
c       "Namelist" input procedure is used.
c       See p221 of DAB's Sun FORTRAN Reference Manual.
c
        implicit none
        namelist /ingest/ ptopdia, prhoice, pzhgt, ptempt, palamd,
     1                      paovrb, pvw
        common /pck_values/ ptopdia, prhoice, pzhgt, ptempt, palamd,
     1                      paovrb
        common /sponge/ pvw
        DOUBLE PRECISION ptopdia, prhoice, pzhgt, ptempt, palamd, paovrb
        DOUBLE PRECISION pvw
c
        open(60,file='spongyice.inp')
c
        read(unit=60, nml=ingest)
c
        close(60)
        return
c
      end
