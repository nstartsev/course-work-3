*
      subroutine rp(t,u,v,w,ut,vt,wt,ox,oy,oz,p,Jmax, Kmax)
      implicit real*8 (a-h,o-z)
	     complex*16
     > u(0:Jmax,0:Kmax),v(0:Jmax,0:Kmax),w(0:Jmax,0:Kmax)
     >,ut(0:Jmax,0:Kmax),vt(0:Jmax,0:Kmax),wt(0:Jmax,0:Kmax)
     >,ox(0:Jmax,0:Kmax),oy(0:Jmax,0:Kmax),oz(0:Jmax,0:Kmax)
     >,p(Jmax,Kmax),dp(Jmax,Kmax)
     >,ci,alci,cial
     >,tmp1 , tmp2 ,tmp3 ,tmp4, tmp5, tmp6
     >,tmp11, tmp12, tmp21,tmp22, tmp31, tmp32
     >,tmp41, tmp42, tmp51, tmp52, tmp61,tmp62
         real*8
     >  dpr(Jmax,Kmax), dpi(Jmax,Kmax)
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
*velocity

      do k=0,Km
        w(0,k)=-w(1,k)
        w(Jm+1,k)=-w(Jm,k)
      end do
      do j=0,Jm
        v(j,0)=-v(j,1)
        v(j,Km+1)=-v(j,Km)
      end do
      do j = 1 , Jm
        do k = 1, Km
            u(j,k) = -(v(j,k)-v(j-1,k))/ym1(j)/cial
     >               -(w(j,k)-w(j,k-1))/zm1(k)/cial
        end do
      end do

      do k=0,Km+1
        u(0,k)=-u(1,k)
        u(Jm+1,k)=-u(Jm,k)
      end do
      do j=0,Jm+1
        u(j,0)=-u(j,1)
        u(j,Km+1)=-u(j,Km)
      end do


*
* Vorticities

      do j=0,Jm
        do k=0,Km
            ox(j,k)=(w(j+1,k)-w(j,k))/yn1(j)
     >            - (v(j,k+1)-v(j,k))/zn1(k)

            oy(j,k)= -cial*w(j,k)
     >             +(u(j,k+1)-u(j,k))/zn1(k)

            oz(j,k)=-(u(j+1,k)-u(j,k))/yn1(j)
     >               +cial*v(j,k)
        end do
      end do
****


      do j = 1 , Jm
        do k = 1, Km
            tmp1 = 0.5d0*(v(j-1,k)*ozb(j-1,k)+v(j,k)*ozb(j,k))
            tmp2 = 0.5d0*(vb(j-1,k)*oz(j-1,k)+vb(j,k)*oz(j,k))
            tmp3 = - 0.5d0*(wb(j,k-1)*oy(j,k-1)+wb(j,k)*oy(j,k))
            tmp4 = - 0.5d0*(w(j,k-1)*oyb(j,k-1)+w(j,k)*oyb(j,k))
            tmp5 = - (oz(j,k)-oz(j-1,k))/(ym1(j)*Re)
            tmp6 = (oy(j,k)-oy(j,k-1))/(zm1(k)*Re)
            ut(j,k) = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6

        end do
      end do

      do j=0, Jm-1
        do k=1, Km
            tmp11 = 0.5d0*(wb(j+1,k)+wb(j,k))*ox(j,k)
            tmp12 = 0.5d0*(wb(j+1,k-1)+wb(j,k-1))*ox(j,k-1)

            tmp1 = 0.5d0*(tmp11+tmp12)

            tmp21 = 0.5d0*(w(j+1,k)+w(j,k))*oxb(j,k)
            tmp22 = 0.5d0*(w(j+1,k-1)+w(j,k-1))*oxb(j,k-1)

            tmp2 = 0.5d0*(tmp21+tmp22)

            tmp3 = - 0.5d0*(ub(j+1,k)+ub(j,k))*oz(j,k)
            tmp4 = - 0.5d0*(u(j+1,k)+u(j,k))*ozb(j,k)
            tmp5 = - (ox(j,k)-ox(j,k-1))/(zm1(j)*Re)
            tmp6 =  cial*oz(j,k)/Re
            vt(j,k) = tmp1+tmp2+tmp3+tmp4 +tmp5+tmp6
        end do
      end do

*      do k = 1, Km
*        write(*,*) vt(0,k)
*      end do
      do j=1,Jm
        do k = 0,Km-1
            tmp1 =   0.5d0*(ub(j,k+1)+ub(j,k))*oy(j,k)
            tmp2 =   0.5d0*(u(j,k+1)+u(j,k))*oyb(j,k)
            tmp31 =  0.5d0*(vb(j,k+1)+vb(j,k))*ox(j,k)
            tmp32 =  0.5d0*(vb(j-1,k+1)+vb(j-1,k))*ox(j-1,k)

            tmp3 = - 0.5d0*(tmp31+tmp32)

            tmp41 = 0.5d0*(v(j,k+1)+v(j,k))*oxb(j,k)
            tmp42 = 0.5d0*(v(j-1,k+1)+v(j-1,k))*oxb(j-1,k)


            tmp4 = -0.5d0*(tmp41 + tmp42)

            tmp5 = - cial*oy(j,k)/Re
            tmp6 = (ox(j,k)-ox(j-1,k))/(ym1(j)*Re)
            wt(j,k) = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6
        end do
      end do


*      do k = 0 , Km+1
*        vt(0,k)  = 0
*        vt(Jm,k) = 0
*      end do

*      do j = 0, Jm+1
*        wt(j,0) = 0
*        wt(j,Km) = 0
*      end do
*
* Pressure


       do j=1,Jm
            do k = 1,Km
                dp(j,k) = cial*ut(j,k) + (vt(j,k)-vt(j-1,k))/(ym1(j))
     >          + (wt(j,k)-wt(j,k-1))/(zm1(k))
            end do
       end do
        do k = 1 , Km
            dp(1,k)=dp(1,k)+apym(1)*yn1(0)*vt(0,k)
            dp(Jm,k) = dp(Jm,k)-cpym(Jm)*vt(Jm,k)*yn1(Jm)
        end do
        do j = 1, Jm
            dp(j,1) = dp(j,1)+apzm(1)*zn1(0)*wt(j,0)
            dp(j,Km) = dp(j,Km)-cpzm(Km)*wt(j,Km)*zn1(Km)
        end do


        do j = 1 , Jm
            do k = 1 , Km
                dpr(j,k) = dreal(dp(j,k))
                dpi(j,k) = dimag(dp(j,k))
            end do
        end do

      call blktri(1,1,Jm,apyp,bpyp,cpyp,1,Km,apzp,bpzp,cpzp,Jmax,
     >  dpr,ierror,pwork)
      call blktri(1,1,Jm,apyp,bpyp,cpyp,1,Km,apzp,bpzp,cpzp,Jmax,
     >  dpi,ierror,pwork)




        do j = 1 , Jm
            do k = 1 , Km
                p(j,k) = dpr(j,k)+ci*dpi(j,k)
            end do
        end do
        do j = 1, Jm-1
            do k = 1 , Km
                vt(j,k) = vt(j,k) - (p(j+1,k)-p(j,k))/yn1(j)
            end do
        end do
        do j = 1, Jm
            do k = 1 , Km-1
                wt(j,k) = wt(j,k) - (p(j,k+1)-p(j,k))/zn1(k)
            end do
        end do
      return
      end
