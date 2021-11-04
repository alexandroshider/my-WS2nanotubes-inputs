!  by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

!  structure for 1T-MoS2   

      program main
      implicit none
      integer, parameter:: i4  = 4
      integer, parameter:: dp  = kind(1.0d0)

      integer(i4), parameter:: natomcell = 12_i4
      integer(i4), parameter:: natype = 3_i4
      integer(i4), parameter:: maxneigh = 6_i4
      real(dp), parameter:: bigspace = 100.0_dp
      real(dp), parameter:: pi = 3.141593_dp
      character(len=1):: orientation
      integer(i4):: nlayer
      integer(i4):: ntot
      integer(i4):: atnumcell(natomcell)
      integer(i4):: atypecell(natomcell)
      integer(i4), allocatable:: atnum(:)
      integer(i4), allocatable:: atype(:)
      real(dp):: bond
      real(dp):: latconst
      real(dp):: thicknessiner
      real(dp):: thickness
      real(dp):: lenxyz(3)
      real(dp):: mdbox(3)
      real(dp), allocatable:: xalat(:)
      real(dp), allocatable:: yalat(:)
      real(dp), allocatable:: zalat(:)
      real(dp), allocatable:: Trasla(:)
      real(dp), allocatable:: Trasla4(:)
      real(dp), allocatable:: coord(:,:,:)

      integer(i4):: jxsf
      integer(i4):: jlammps
      integer(i4):: jovito
      integer(i4):: i
      integer(i4):: j
      integer(i4):: k
      integer(i4):: l
      integer(i4):: ii
      integer(i4):: il
      integer(i4):: nx
      integer(i4):: ny
      integer(i4):: nz
      integer(i4):: N
      integer(i4):: quit
      integer(i4):: num
      logical:: ltube
      real(dp):: box_x
      real(dp):: box_y
      real(dp):: a
      real(dp):: b
      real(dp):: c
      real(dp):: ci
      real(dp):: d
      real(dp):: h
      real(dp):: x0
      real(dp):: y0
      real(dp):: z0
      real(dp):: zav
      real(dp):: radius
      real(dp):: radius0
      real(dp):: radius1
      real(dp):: radius2
      real(dp):: diameter
      real(dp):: theta
      real(dp):: xmin
      real(dp):: ymin
      real(dp):: ly
      real(dp):: diam
      real(dp):: diametro
      real(dp):: YM
      real(dp):: ZM
      real(dp):: Yme
      real(dp):: Zme
      real(dp):: ang
      real(dp):: distancia
      real(dp):: FG
      real(dp):: idvacan
      character(len=1)::ori
      integer(i4)::forma
      OPEN( 10,  FILE='datos.dat' )
      open(20, file='Fuera.dat')
      open(30,file='Fueralammps.dat')

   ! SelecciÃ³n de nano tubo o nano capa
!10   write(*,*) 'Nanotubo o Nanocapa, (1) o (2)'
!   read(*,*) forma
!   If (forma.eq.1) then
!   ltube = .true.
!   Else if (forma.eq.2) then
!   ltube = .false.
!   Else
!       write(*,*) 'No se encontro ningun tipo'
!       go to 10
!       write(*,*)
!   End if

   ! tube(tubo) or sheet
      ltube = .true.
      idvacan=0
   ! Orientation, z=zigzag, a=armchair
      write(*,*) 'Ingrese la orientacion (a) o (z)'
      read(*,*) ori
      orientation = ori
   ! Numero de monocapas
      nlayer = 1

!   lenxyz(1) = Distancia ligada a la longitud
!   lenxyz(2) = Distancia ligada al diametro
      write(*,*) 'Ingrese el diametro'
      read(*,*) diametro
      write(*,*)'Ingrese el n£mero de vectores de Burguers que desea'
      write(*,*)'recorrer la estructura (0 para no recorrerla)'
      read(*,*) FG
      diam=(diametro)*pi
      lenxyz(1) = 103
      lenxyz(2) = diam


   ! for tube, this version only works for single-layer
      if(ltube .and. nlayer.ne.1)then
        write(*,*)'MUST set nlayer = 1 for tube, input nlayer = ',nlayer
      stop
      end if

      latconst = 3.1998
      bond = 2.4193

   ! derived parameters
   ! inplane 'bond', mimic graphene                                                                    
      b = latconst/sqrt(3.0_dp)
   ! heighth of X atom
      h = sqrt(bond*bond-b*b)
      thicknessiner = 2.0*h
   ! thickness, arbitrary value, should be revised forfew-layer or bulk system
      thickness = 6.0

      a = latconst
      c = thickness
      ci = thicknessiner
      allocate(coord(3,natomcell,nlayer))
   ! armchair and zigzag are different by the followingcoordinate transformation.
   !            0  1  0
   !           -1  0  0
   !            0  0  1
   ! i.e., x --> y, y--> -x, and z --> z
   ! armchair
      if(orientation.eq.'z')then
      do i = 1, nlayer, 1
         x0 = 0.5*b
         y0 = 0.0
         z0 = real(i-1) * thickness
         coord(1,1,i) = 0.0_dp + x0
         coord(2,1,i) = 0.0_dp + y0
         coord(3,1,i) = 0.0_dp + z0
         coord(1,2,i) = -0.5_dp*b + x0
         coord(2,2,i) = 0.5_dp*sqrt(3.0)*b + y0
         coord(3,2,i) = h + z0
         coord(1,3,i) = 0.5*b + x0
         coord(2,3,i) = 0.5*sqrt(3.0)*b + y0
         coord(3,3,i) = -h + z0
         coord(1,4,i) = 0.0_dp + x0
         coord(2,4,i) = a + y0
         coord(3,4,i) = 0.0_dp + z0
         coord(1,5,i) = -0.5_dp*b + x0
         coord(2,5,i) = a + 0.5_dp*sqrt(3.0)*b + y0
         coord(3,5,i) = h + z0
         coord(1,6,i) = 0.5*b + x0
         coord(2,6,i) = a + 0.5*sqrt(3.0)*b + y0
         coord(3,6,i) = -h + z0
         coord(1,7,i) = 1.5*b + x0
         coord(2,7,i) = a + 0.5*sqrt(3.0)*b + y0
         coord(3,7,i) = 0.0 + z0
         coord(1,8,i) = b + x0
         coord(2,8,i) = a + y0
         coord(3,8,i) = h + z0
         coord(1,9,i) = 2.0*b + x0
         coord(2,9,i) = a + y0
         coord(3,9,i) = -h + z0
         coord(1,10,i) = 1.5*b + x0
         coord(2,10,i) = 0.5*sqrt(3.0)*b + y0
         coord(3,10,i) = 0.0 + z0
         coord(1,11,i) = b + x0
         coord(2,11,i) = 0.0 + y0
         coord(3,11,i) = h + z0
         coord(1,12,i) = 2.0*b + x0
         coord(2,12,i) = 0.0 + y0
         coord(3,12,i) = -h + z0
      end do !i
      box_x = 3.0*b
      box_y = 2.0*a
   ! zigzag
      else if(orientation.eq.'a')then
      do i = 1, nlayer, 1
         x0 = 0.5*b
         y0 = 0.0
         z0 = real(i-1) * thickness
         coord(1,1,i) = 0.0_dp + y0
         coord(2,1,i) = 0.0_dp + x0
         coord(3,1,i) = 0.0_dp + z0
         coord(1,2,i) = -0.5_dp*sqrt(3.0)*b - y0
         coord(2,2,i) = -0.5_dp*b + x0
         coord(3,2,i) = h + z0
         coord(1,3,i) = -0.5*sqrt(3.0)*b - y0
         coord(2,3,i) = 0.5*b + x0
         coord(3,3,i) = -h + z0
         coord(1,4,i) = -a - y0
         coord(2,4,i) = 0.0_dp + x0
         coord(3,4,i) = 0.0_dp + z0
         coord(1,5,i) = -a - 0.5_dp*sqrt(3.0)*b - y0
         coord(2,5,i) = -0.5_dp*b + x0
         coord(3,5,i) = h + z0
         coord(1,6,i) = -a - 0.5*sqrt(3.0)*b - y0
         coord(2,6,i) = 0.5*b + x0
         coord(3,6,i) = -h + z0
         coord(1,7,i) = -a - 0.5*sqrt(3.0)*b - y0
         coord(2,7,i) = 1.5*b + x0
         coord(3,7,i) = 0.0 + z0
         coord(1,8,i) = -a - y0
         coord(2,8,i) = b + x0
         coord(3,8,i) = h + z0
         coord(1,9,i) = -a - y0
         coord(2,9,i) = 2.0*b + x0
         coord(3,9,i) = -h + z0
         coord(1,10,i) = -0.5*sqrt(3.0)*b - y0
         coord(2,10,i) = 1.5*b + x0
         coord(3,10,i) = 0.0 + z0
         coord(1,11,i) = -0.0 - y0
         coord(2,11,i) = b + x0
         coord(3,11,i) = h + z0
         coord(1,12,i) = -0.0 - y0
         coord(2,12,i) = 2.0*b + x0
         coord(3,12,i) = -h + z0
      end do !i
      box_x = sqrt(3.0_dp)*b * 2.0
      box_y = 3.0_dp*b
      else
      write(*,*)"unknown orientation type",orientation
      stop
      end if

   ! set atom number
      atnumcell(1) =  42
      atnumcell(2) =  16
      atnumcell(3) =  16
      atnumcell(4) =  42
      atnumcell(5) =  16
      atnumcell(6) =  16
      atnumcell(7) =  42
      atnumcell(8) =  16
      atnumcell(9) =  16
      atnumcell(10) =  42
      atnumcell(11) =  16
      atnumcell(12) =  16

   ! set atom type
      atypecell(1) = 1
      atypecell(2) = 2
      atypecell(3) = 3
      atypecell(4) = 1
      atypecell(5) = 2
      atypecell(6) = 3
      atypecell(7) = 1
      atypecell(8) = 2
      atypecell(9) = 3
      atypecell(10) = 1
      atypecell(11) = 2
      atypecell(12) = 3

   ! reset size
      nx = lenxyz(1)/box_x
      lenxyz(1) = real(nx) * box_x
      ny = lenxyz(2)/box_y
      lenxyz(2) = real(ny) * box_y
      nz = nlayer
      lenxyz(3) = real(nz)*thickness

   ! total atom number
      N = nlayer * natomcell * nx * ny
      ntot = N

   ! allocate ntot-related arrays
      allocate(xalat(ntot), yalat(ntot), zalat(ntot),Trasla(ntot))
      allocate(atnum(ntot), atype(ntot),Trasla4(ntot))

   ! generate structure
      do i = 1, nx
      do j = 1, ny
         x0 = real(i-1) * box_x
         y0 = real(j-1) * box_y
         z0 = 0.0_dp
         k = (i-1)*ny*natomcell*nlayer + (j-1)*natomcell*nlayer
         do l = 1, nlayer
            do ii = 1, natomcell
               k = k + 1
               xalat(k) = x0 + coord(1,ii,l)
               yalat(k) = y0 + coord(2,ii,l)
               zalat(k) = z0 + coord(3,ii,l)
               atnum(k) = atnumcell(ii)
               atype(k) = atypecell(ii)
            end do !ii
         end do !l
      end do !j
      end do !i
      if(k.ne.N)then
       write(*,*)"inconsistent for atom number",k,N
      stop
      end if

   ! roll up plane into tube, only works for single-layer
   ! the middle layer is purely bent, the inner layer is compressed, 
   ! the outer layer is tensiled
      if(ltube)then
      zav = 0.0
      do i = 1, ntot
         zav = zav + zalat(i)
      end do !i
      zav = zav/real(ntot)

      diameter = lenxyz(2)/pi
      radius0 = 0.5_dp*diameter
      radius1 = 0.5_dp*diameter - 0.5*thicknessiner
      radius2 = 0.5_dp*diameter + 0.5*thicknessiner
      ly = pi * diameter
      do i = 1, ntot, 1
         ! inner layer
         if(zalat(i).lt.(zav-0.1))then
            radius = radius1
         ! outer layer
         else if(zalat(i).gt.(zav+0.1))then
            radius = radius2
         ! middle layer
         else
            radius = radius0
         end if
         theta = 2.0_dp*pi*yalat(i)/ly
         xalat(i) = xalat(i)
         yalat(i) = radius * cos(theta)
         zalat(i) = radius * sin(theta)
      end do !i
      end if

      mdbox(:) = lenxyz(:)
      if(ltube)then
      mdbox(2) = mdbox(2) + bigspace
      mdbox(3) = mdbox(2)
      else
      mdbox(3) = mdbox(3) + bigspace
      end if

   ! shift to center of simulation box
      if(ltube)then
      do i = 1, ntot
         yalat(i) = yalat(i) + 0.5_dp * mdbox(2)
         zalat(i) = zalat(i) + 0.5_dp * mdbox(3)
      end do !i
      else
      do i = 1, ntot
         zalat(i) = zalat(i) + 0.5_dp * mdbox(3)
      end do !i
  ! small shift in y, avoid possible boundary effects for plane layer
      ymin = 1000.0
      do i = 1, ntot
      if(yalat(i).lt.ymin) ymin = yalat(i)
      end do !i
      do i = 1, ntot
       yalat(i) = yalat(i) - ymin + 0.1
      end do !i
      end if

   ! small shift in x, avoid possible boundary effects
      xmin = 1000.0
      do i = 1, ntot
      if(xalat(i).lt.xmin) xmin = xalat(i)
      end do !i
      do i = 1, ntot
         xalat(i) = xalat(i) - xmin + 0.1
      end do !i

   ! output structure information
      write(*,'(''# ntot = '', i12)')ntot
      if(ltube)then
      write(*,'(''# system = MoS2 tube, length, diameter = '',2f16.4)')
     &lenxyz(1),diameter
      write(*,*)
      write(*,'(''# system = MoS2 tube, rad0, rad1, rad2 = '',2f16.4)')
     &radius0,radius1,radius2
      else
      write(*,'(''# system = MoS2 sheet, lenxyz(1:3) = '',3f16.4 )')
     &lenxyz(1:3)
      end if
      write(*,'(''# mdbox(1:3) = '',3f16.4 )')mdbox(1:3)
      write(*,'(''# nx, ny, nz = '',3i12 )')nx, ny, nz
   !Indices n,m   
      If (orientation.eq.'a')then
      write(*,*) '# Indices n, m','    ', '(',ny,',',ny,')'
      Else
      write(*,*) '# Indices n, m','    ', '(',ny*2,',',0,')'
      End if

   !-------------------------------------------------------------------------------------------------
      ! Informacion de salida de a.out
      write(10,'(''# ntot = '', i12)')ntot
      if(ltube)then
      write(10,'(''# system = MoS2 tube, length, diameter = '',2f16.4)')
     &lenxyz(1),diameter
      write(10,*)
      write(10,'(''# system = MoS2 tube, rad0, rad1, rad2 = '',2f16.4)')
     &radius0,radius1,radius2
      else
      write(10,'(''# system = MoS2 sheet, lenxyz(1:3) = '',3f16.4 )')
     &lenxyz(1:3)
      end if
      write(10,'(''# mdbox(1:3) = '',3f16.4 )')mdbox(1:3)
      write(10,'(''# nx, ny, nz = '',3i12 )')nx, ny, nz
      !Indices n,m
      If (orientation.eq.'a')then
      write(10,*) '# Indices n, m','    ', '(',ny,',',ny,')'
      Else
      write(10,*) '# Indices n, m','    ', '(',2*ny,',',0,')'
      End if
   !-------------------------------------------------------------------------------------------------

   ! out .xsf file
      jxsf = 11
      open(jxsf,file='xyz.xsf',status='unknown',form='formatted')
      write(jxsf,*)"# structure data, .xsf formate, viewed byXCRYSDEN"
      write(jxsf,*)"# ntot =",ntot
      write(jxsf,*)"ATOMS"
      do i = 1, ntot
      write(jxsf,'(i12,3f20.8)')atnum(i), xalat(i), yalat(i), zalat(i)                                    
      end do !i
      close(jxsf)

   ! out lammps.dat file
      jlammps = 12
      open(jlammps,file='lammps.dat',status='unknown',form='formatted')
      write(jlammps,*)"# Lammps sturcture data file"
      write(jlammps,*)"# ntot =",ntot
      write(jlammps,*)ntot,"   atoms"
      write(jlammps,*)natype," atom types"
      write(jlammps,*)
      write(jlammps,*)0.0_dp,mdbox(1),"   xlo xhi"
      write(jlammps,*)0.0_dp,mdbox(2),"   ylo yhi"
      write(jlammps,*)0.0_dp,mdbox(3),"   zlo zhi"
      write(jlammps,*)
      write(jlammps,*)"Atoms"
      write(jlammps,*)
      write(*,*)
      do i = 1, ntot
      write(jlammps,'(2i12,3f12.4)')i,atype(i),
     &xalat(i),yalat(i),zalat(i)
      end do !i
      close(jlammps)


   !Archivo para h‚lice
      open(32, file='Helice.dat')
      open(33, file='Traslacion.dat')
            YM=yalat(1)
            ZM=zalat(1)
            Yme=yalat(1)
            Zme=zalat(1)
      do i=1,ntot

      if  (yalat(i).gt.YM) then
        YM=yalat(i)
        end if

        if (zalat(i).gt.ZM) then
           ZM=zalat(i)
           end if

           if (yalat(i).lt.Yme) then
           Yme=yalat(i)
           end if

           if(zalat(i).lt.Zme) then
           Zme=zalat(i)
           end if

           end do
        write(*,*) 'La coordenada y mayor es:' ,YM
        write(*,*) 'La coordenada y menor es:' ,Yme
        write(*,*) 'La coordenada z mayor es:' ,ZM
        write(*,*) 'La coordenada z menor es:' ,Zme
        write(*,*) 'El centro est  en',(YM+Yme)/2,(ZM+Zme)/2,'->(y,z)'
        write(32,*) 0,(YM+Yme)/2,(ZM+Zme)/2
         distancia=3.119*FG

         do i=1, ntot

      ang=ATAN((((ZM+Zme)/2)-zalat(i))/(((YM+Yme)/2)-yalat(i)))

      !Cuadrantes
      if( (-(YM+Yme)/2+yalat(i)).ge.0 .and. (-(ZM+Zme)/2+zalat(i)).ge.0)         !Primer cuadrante
     &then
       Trasla(i)=xalat(i)-distancia*(ang/(2*pi))
      end if

      if((-(YM+Yme)/2+yalat(i)).le.0 .and. (-(ZM+Zme)/2+zalat(i)).ge.0)         !Segundo cuadrante
     &then
       Trasla(i)=xalat(i)-distancia*((ang+pi)/(2*pi))
      end if

      if((-(YM+Yme)/2+yalat(i)).le.0 .and. (-(ZM+Zme)/2+zalat(i)).le.0)       !Tercer cuadrante
     &then
       Trasla(i)=xalat(i)-distancia*((ang+pi)/(2*pi))
      end if

      if((-(YM+Yme)/2+yalat(i)).ge.0 .and. (-(ZM+Zme)/2+zalat(i)).le.0)       !Cuarto cuadrante
     &then
       Trasla(i)=xalat(i)-distancia*((ang+2*pi)/(2*pi))
      end if

      !Condiciones para que no salga de la caja
      if (Trasla(i).gt.mdbox(1)) then
      Trasla4(i)=Trasla(i)-mdbox(1)
      Trasla(i)=Trasla4(i)
      end if

      if (Trasla(i).lt.0.0) then
      if ((-(YM+Yme)/2+yalat(i)).le.0 .and. (-(ZM+Zme)/2+zalat(i)).le.0
     &.or.((-(YM+Yme)/2+yalat(i)).ge.0 .and. (-(ZM+Zme)/2+zalat(i)).le.0
     &))then
      !Trasla4(i)=Trasla(i)+mdbox(1)
      !Trasla(i)=Trasla4(i)
      end if
      end if

      if (Trasla(i).lt.0.0) then
      if(((-(YM+Yme)/2+yalat(i)).ge.0 .and. (-(ZM+Zme)/2+zalat(i)).ge.0)
     &.or.((-(YM+Yme)/2+yalat(i)).le.0 .and. (-(ZM+Zme)/2+zalat(i)).ge.0
     &))then
      !Trasla4(i)=Trasla(i)+mdbox(1)
      !Trasla(i)=Trasla4(i)
      end if
      end if

      write(33,*)Trasla(i), yalat(i), zalat(i)
      write(32,*)xalat(i), yalat(i),zalat(i)

      end do
      close(32)
      close(33)

      jovito = 55
      open(jovito,file='xyzhelic.ovito',status='unknown',
     &form='formatted')
      write(jovito,'(''ITEM: TIMESTEP'')')
      write(jovito,'(i12)')1
      write(jovito,'(''ITEM: NUMBER OF ATOMS'')')
      write(jovito,'(i12)')ntot
      write(jovito,'(''ITEM: BOX BOUNDS pp pp pp'')')
      write(jovito,'(2f20.8,''   xlo xhi'')')0.0_dp,mdbox(1)
      write(jovito,'(2f20.8,''   ylo yhi'')')0.0_dp,mdbox(2)
      write(jovito,'(2f20.8,''   zlo zhi'')')0.0_dp,mdbox(3)
      write(jovito,'(''ITEM: ATOMS id type x y z c_csym'')')
      
      quit=1
      num=0
      do i = 1, ntot
       num=num+1
       if (quit.eq.idvacan) then
       !'(2i12,3f20.8)'
       write(20,*)i,atype(i),Trasla(i),yalat(i)
     &,zalat(i),0.0
       num=num-1

       else

       write(jovito,'(2i12,4f20.8)')num,atype(i),Trasla(i),yalat(i)
     &,zalat(i),0.0

       end if
       quit=quit+1
       end do !i
      close(jovito)
   
   ! out ovito file
      jovito = 13
      open(jovito,file='xyz.ovito',status='unknown',form='formatted')
      write(jovito,'(''ITEM: TIMESTEP'')')
      write(jovito,'(i12)')1
      write(jovito,'(''ITEM: NUMBER OF ATOMS'')')
      write(jovito,'(i12)')ntot
      write(jovito,'(''ITEM: BOX BOUNDS pp pp pp'')')
      write(jovito,'(2f20.8,''   xlo xhi'')')0.0_dp,mdbox(1)
      write(jovito,'(2f20.8,''   ylo yhi'')')0.0_dp,mdbox(2)
      write(jovito,'(2f20.8,''   zlo zhi'')')0.0_dp,mdbox(3)
      write(jovito,'(''ITEM: ATOMS id type x y z c_csym'')')
      do i = 1, ntot
      write(jovito,'(2i12,4f20.8)')i,atype(i),xalat(i),yalat(i),
     &zalat(i),0.0
      end do !i
      close(jovito)
      
      !Archivo xyz
      open(56,file='crystal_helic.xyz',status='unknown',
     &form='formatted')
      write(56,*)ntot
      write(56,*) ''
      do i=1,ntot
      if(atype(i).eq.1)then
      write(56,*)'Mo',Trasla(i),yalat(i),zalat(i)
      else
      write(56,*)'S',Trasla(i),yalat(i), zalat(i)
      end if
      end do
      close(56)
      
      
      jlammps = 21
      open(jlammps,file='lammpshelic.dat',status='unknown',
     &form='formatted')
      write(jlammps,'(''# Lammps sturcture data file'')')
      write(jlammps,'(''# ntot = '',i12)')ntot
      write(jlammps,'(i12,''   atoms'')')ntot
      write(jlammps,'(i12,'' atom types'')')natype
      write(jlammps,*)
      write(jlammps,'(2f20.8,''   xlo xhi'')')0.0_dp,mdbox(1)
      write(jlammps,'(2f20.8,''   ylo yhi'')')0.0_dp,mdbox(2)
      write(jlammps,'(2f20.8,''   zlo zhi'')')0.0_dp,mdbox(3)
      write(jlammps,*)
      write(jlammps,'(''Atoms'')')
      write(jlammps,*)
      quit=1
      num=0
      do i = 1, ntot
       num=num+1
       if (quit.eq.idvacan) then
      ! '(2i12,3f20.8)'
       write(30,*)i,atype(i),Trasla(i),yalat(i)
     &,zalat(i),0.0
       num=num-1
       
       else
       
       write(jlammps,'(2i12,4f20.8)')num,atype(i),Trasla(i),yalat(i)
     &,zalat(i)
     
       end if
       quit=quit+1
       end do !i
       
      close(jlammps)
      
      

      write(*,*)"Job done, Sir!"
       pause
      end
