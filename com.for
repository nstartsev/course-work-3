*
      subroutine com(epsy)
      parameter (Jmax=64)
      parameter (Kmax=64)
      implicit real*8 (a-h,o-z)
        real*8
     >  q(Jmax,Kmax), rwork(10000),dpr(Jmax,Kmax)
      common
     >/dimy/yn(0:1024),yn1(0:1024),ym(0:1025),ym1(0:1025),Jm
     >/dimz/zn(0:1024),zn1(0:1024),zm(0:1025),zm1(0:1025),Km
     >/pry/apym(256),bpym(256),cpym(256)
     >/prz/apzm(256),bpzm(256),cpzm(256)
     >/pyp/apyp(256),bpyp(256),cpyp(256)
     >/pzp/apzp(256),bpzp(256),cpzp(256)
     >/pyb/apyb(256),bpyb(256),cpyb(256)
     >/pzb/apzb(256),bpzb(256),cpzb(256)
     >/pw/pwork(10000)
     >/alre/al, Re
     >/basefl/ub(0:257,0:257),Oxb(0:257,0:257),vb(0:257,0:257)
     >,Oyb(0:257,0:257)
     >,wb(0:257,0:257),Ozb(0:257,0:257)
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
*


        hy=2.d0/Jm
        do j=0,Jm
          yn(j)=hy*j-1
          yn1(j)=hy
        end do
        do j=0,Jm+1
          ym(j)=(j-0.5d0)*hy-1
          ym1(j)=hy
        end do

        hz=2.d0/Km
        do k=0,Km
          zn(k)=hz*k-1
          zn1(k)=hz
        end do
        do k=0,Km+1
          zm(k)=(k-0.5d0)*hz-1
          zm1(k)=hz
        end do

* race coefficients

        do j=1,Jm
         apym(j) = 1.d0/(ym1(j)*yn1(j-1))
         cpym(j) = 1.d0/(ym1(j)*yn1(j))
         bpym(j) = -apym(j)- cpym(j)
        end do
        do k=1,Km
         apzm(k) = 1.d0/(zm1(k)*zn1(k-1))
         cpzm(k) = 1.d0/(zm1(k)*zn1(k))
         bpzm(k) = -apzm(k)- cpzm(k)
        end do
*

* pressure coefficients
      do j = 1, Jm
        apyp(j) = apym(j)
        bpyp(j) = bpym(j)-al**2
        cpyp(j) = cpym(j)
      end do

      do k = 1 , Km
        apzp(k) =  apzm(k)
        bpzp(k) =  bpzm(k)
        cpzp(k) =  cpzm(k)
      end do
        bpyp(1)=   apyp(1)+bpyp(1)
        apyp(1) = 0.
        bpyp(Jm) = bpyp(Jm)+cpyp(Jm)
        cpyp(Jm) = 0.

        bpzp(1)=  apzp(1)+bpzp(1)
        apzp(1) = 0
        bpzp(Km) = bpzp(Km)+cpzp(Km)
        cpzp(Km) = 0


      do j = 1, Jm
        apyb(j) = apym(j)
        bpyb(j) = bpym(j)
        cpyb(j) = cpym(j)
      end do

      do k = 1 , Km
        apzb(k) =  apzm(k)
        bpzb(k) =  bpzm(k)
        cpzb(k) =  cpzm(k)
      end do
        bpyb(1)=  -apyb(1)+bpyb(1)
        apyb(1) = 0.
        bpyb(Jm) = bpyb(Jm)-cpyb(Jm)
        cpyb(Jm) = 0.

        bpzb(1)=  -apzb(1)+bpzb(1)
        apzb(1) = 0
        bpzb(Km) = bpzb(Km)-cpzb(Km)
        cpzb(Km) = 0
*

* Base flow
        do j = 1, Jm
            do k =1 , Km
                q(j,k) = -1
            end do
        end do
      call blktri(0,1,Jm,apyb,bpyb,cpyb,1,Km,apzb,bpzb,cpzb,Jmax,
     >  q,ierror,rwork)
      call blktri(1,1,Jm,apyb,bpyb,cpyb,1,Km,apzb,bpzb,cpzb,Jmax,
     >  q,ierror,rwork)
      umax = 0
      do j = 1 , Jm
        do k = 1, Km
            umax = max(abs(q(j,k)), umax)
        end do
      end do
      do j = 1 , Jm
        do k = 1 , Km
            q(j,k) = q(j,k)/umax
        end do
      end do
        do j=1,Jm
            do k=1,Km
*                ub(j,k)=(1.-ym(j)**2)*(1.-zm(k)**2)
                 ub(j,k) = q(j,k)
            end do
        end do
        do j=0,Jm+1
            do k=0,Km+1
                vb(j,k)=0
                wb(j,k)=0
            end do
        end do


        do k=0,Km+1
            ub(0,k)=-ub(1,k)
            ub(Jm+1,k)=-ub(Jm,k)
        end do
        do j=0,Jm+1
            ub(j,0)=-ub(j,1)
            ub(j,Km+1)=-ub(j,Km)
        end do

        do j=0,Jm+1
            vb(j,0)=-vb(j,1)
            vb(j,Km+1)=-vb(j,Km)
        end do

        do k=0,Km+1
            wb(0,k)=-wb(1,k)
            wb(Jm+1,k)=-wb(Jm,k)
        end do



        do j=0,Jm
            do k=0,Km
                oxb(j,k)=-(vb(j,k+1)-vb(j,k))/zn1(k)
     >          +(wb(j+1,k)-wb(j,k))/yn1(j)
            end do
        end do
        do j=0,Jm+1
            do k=0,Km
                oyb(j,k)=(ub(j,k+1)-ub(j,k))/zn1(k)
     >          -cial*wb(j,k)
            end do
        end do
        do j=0,Jm
            do k=0,Km+1
                ozb(j,k)=-(ub(j+1,k)-ub(j,k))/yn1(j)
     >          +cial*vb(j,k)
            end do
        end do
*


* fill pwork
        call blktri(0,1,Jm,apyp,bpyp,cpyp,1,Km,apzp,bpzp,cpzp,Jmax,
     >               dpr,ierror,pwork)
*


      return
      end
