  program CORTOPDB
!**************************************************************************************************
!a. Author - Amitava Roy                                                                          *
!b. Contact information - amitava.roy@nih.gov, amitroy74@gmail.com                                *
!c. Date - August 3rd, 2017                                                                       *
!d. Version - 1.00                                                                                *
!e. Description - The F90 code reads a CHARMM formatted coordinate file and write a PDB formatted *
!   coordinate file. The program can read extended formatted CHARMM coordinate file and can       *
!   preserve the residue ID column in the CHARMM coordinate file as the residue numbers in the    *
!   PDB file. HSD, HSN and HSE residues are renamed as HIS in the PDB file.                       *
!f. Detailed usage instructions -                                                                 *
!    cor_to_pdb  xxx.pdb xxx.cor [ext] [res]                                                      *
!    The input can be in any order. Files are recognized by the extensions. Consequently the PDB  *
!    file has to have an extension .pdb, the coordinate files has to have extension .cor.         *
!    xxx.pdb – Output PDB file name. Required. Maximum length of file name is 255 character.      *
!    xxx.cor – Input CHARMM formated coordinate file. Required.                                   *
!              Maximum length of file name is 255 character.                                      *
!    ext - Flag to indicate the input coordinate file is in extended format. Optional.            *
!    res – Flag to indicate that the residue numbers will be taken from RESID column of the       *
!          coordinate files. Otherwise the residue numbers will start from 1. Optional.           *
!g. Example of usage -                                                                            *
!    Compilation: ifort -mcmodel=medium -shared-intel -o cor_to_pdb cor_to_pdb.f90                *
!                 pgf90 -o cor_to_pdb cor_to_pdb.f90                                              *
!    Change the values below in the code as needed. Higher values will require higher memory.     *
!    mxcr (default=150000)- Maximum number of of residues.                                        *
!    mxca (default=50)- Maximum number of atoms per residue.                                      *
!   Usage: cor_to_pdb 5iq9.cor 5iq9.pdb ext res                                                   *
!h. Keywords - PDB, CHARMM, FORTRAN, COORDINATE                                                   *
!**************************************************************************************************
! INPUTS                                                               

  implicit none
  integer, parameter   :: mxcr=150000,mxca=50
  character(len=6), parameter :: duma="ATOM  "
  integer              :: i,j,k,l,m,maxa,mxar(mxcr),maxr,nremark
  real*4               :: x(mxcr,mxca),y(mxcr,mxca),z(mxcr,mxca),tfact(mxcr,mxca)
  real*4               :: lx,ly,lz,ltfact
  character(len=255)   :: arg,pdb,cor,larg,remarks(50)
  character(len=4)     :: rnam(mxcr,mxca),anam(mxcr,mxca),cid(mxcr,mxca),presid(mxcr,mxca)
  character(len=4)     :: lrnam,lanam,lcid,lresid
  logical              :: lpdb,lcor,lext,lres,cter,nace,c1,c2,n1,n2,n3

! Charmm coordinate format
100 format(i5,i5,1x,a4,1x,a4,3f10.5,1x,a4,1x,a4,f10.5)
101 format(i10,i10,2x,a4,2x,a8,4x,3f20.10,2x,a4,6x,a4,f20.10)
! PDB file format
110 format(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3(f8.3),2(f6.2))  


   pdb=" "
   cor=" "
   lpdb=.false.
   lcor=.false.
   lext=.false.
   lres=.false.
   mxar=0

! Parse the input
   i=0
   do
    call get_command_argument(i,arg)
    if (len_trim(arg).eq.0) exit
    j=len_trim(arg)
    larg=" "
    larg(1:j)=trim(arg)
    if((larg(j-3:j).eq.'.pdb').or.(larg(j-3:j).eq.'.PDB')) then
     pdb(1:j)=larg
     lpdb=.true.
     open(unit=21,file=pdb(1:j),status='new',err=902)
    endif

    if((larg(j-3:j).eq.'.cor').or.(larg(j-3:j).eq.'.COR')) then
     cor(1:j)=larg
     lcor=.true.
     open(unit=11,file=cor(1:j),status='old',err=901)
    endif 

    if((larg.eq.'ext').or.(larg.eq.'EXT')) lext=.true.
    if((larg.eq.'res').or.(larg.eq.'RES')) lres=.true.

    i=i+1
   enddo

   if(.not. lcor) then
    print '(a25)', "No input coordinate file!"
    goto 999
   endif

   if(.not. lpdb) then
    print '(a19)', "No output PDB file!"
   endif

!Read charmm coordinate file
   if(lcor) then
    nremark=0
    do
     read(11,'(a)') arg
     if (arg(1:1).ne."*") exit
     nremark=nremark+1
     j=len_trim(arg)
     remarks(nremark)=arg(2:j)
    enddo
    if(lext) then
     read(arg(1:10),'(i10)') maxa
    else
     read(arg(1:5),'(i5)') maxa
    endif

    j=0
    m=0
    maxr=0
    do i=1,maxa
     if(lext) then
      read(11,101) l,k,lrnam,lanam,lx,ly,lz,lcid,lresid,ltfact
     else
      read(11,100) l,k,lrnam,lanam,lx,ly,lz,lcid,lresid,ltfact
     endif
     if(k.ne.j) then
      if(k.gt.1) mxar(maxr)=m
      m=0
      maxr=maxr+1
      if(maxr.gt.mxcr) then
       print '(a50)', "Number of residue is larger than allocated memory."
       print '(a28)', "Increase MXCR and recompile."
       goto 999
      endif
     endif
     m=m+1
     if(m.gt.mxca) then
      print '(a60)', "Number of atom in a residue is larger than allocated memory."
      print '(a28)', "Increase MXCA and recompile."
      goto 999
     endif

     rnam(k,m)=lrnam
     anam(k,m)=lanam
     x(k,m)=lx
     y(k,m)=ly
     z(k,m)=lz
     cid(k,m)=lcid
     presid(k,m)=lresid
     tfact(k,m)=ltfact
     j=k
    enddo
    mxar(maxr)=m
   endif

! Write pdb file
   if(lpdb) then
    do i=1,nremark
     write(21,'(a6,1x,a)') "REMARKS",trim(remarks(i))
    enddo

    do i=1,maxr
     cter=.false.
     c1=.false.
     c2=.false.
     do j=1,mxar(i)
      if(anam(i,j)(1:3).eq.'OT1') c1=.true.
      if(anam(i,j)(1:3).eq.'OT2') c2=.true.
     enddo
     if(c1.and.c2) cter=.true.
     do j=1,mxar(i)
      if(i.gt.1) then
       k=sum(mxar(1:i-1))+j
      else
       k=j
      endif
      lcid=cid(i,j)(1:1)
      lanam="    "
      if(anam(i,j)(4:4).eq." ") then
       lanam(2:4)=anam(i,j)(1:3)
      else
       lanam(1:1)=anam(i,j)(4:4)
       lanam(2:4)=anam(i,j)(1:3)
      endif
      lrnam=rnam(i,j)(1:3)
      if(lrnam.eq."HSD") lrnam="HIS"
      if(lrnam.eq."HSE") lrnam="HIS"
      if(lrnam.eq."HSN") lrnam="HIS"
      if((lrnam.eq."ILE").and.(lanam.eq." CD ")) lanam=" CD1"
      if((cter).and.(lanam.eq." OT1")) lanam=" OXT"
      if((cter).and.(lanam.eq." OT2")) lanam=" O  "
      if(lres) then
       read(presid(i,j),*) l
      else
       l=i
      endif
      write(21,110) duma,k,lanam," ",lrnam,lcid,l," ",x(i,j),y(i,j),z(i,j),1.00,tfact(i,j) 
     enddo
    enddo
    write(21,'(a3)') "END"

   endif

!   Error messages
    goto 999
901 print *, 'error in openning file, file may not exist!'
    print *, cor(1:j)
    goto 999
902 print *, 'error in writing file, file may already exist!'
    print *, pdb(1:j)
    goto 999
999 continue
  end program
