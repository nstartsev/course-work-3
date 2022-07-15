       implicit real*8 (a-h,o-z)

      parameter (Jmax=64)
      parameter (Kmax=64)

      integer LWMAX
      parameter ( LWMAX = 10000 )
      parameter  (KJmax = 5000)
      real*8
     >rwork(LWMAX)
      complex*16
     > u(0:Jmax,0:Kmax),v(0:Jmax,0:Kmax),w(0:Jmax,0:Kmax)
     >,u1(0:Jmax,0:Kmax),v1(0:Jmax,0:Kmax),w1(0:Jmax,0:Kmax)
     >,p(Jmax,Kmax)
     >,ox(0:Jmax,0:Kmax), oy(0:Jmax,0:Kmax),oz(0:Jmax,0:Kmax)
     >,Amat(Jmax*Kmax,Jmax*Kmax), evecl(Jmax*Kmax,Jmax*Kmax)
     >,evec(Jmax*Kmax,Jmax*Kmax),eval(Jmax*Kmax),work(LWMAX)
     >,Amat1(Jmax*Kmax,Jmax*Kmax), tmp1
     >, ppert(0:Jmax, 0:Kmax)
     >, ut(0:Jmax,0:Kmax),vt(0:Jmax,0:Kmax),wt(0:Jmax,0:Kmax)
     >,ci,alci,cial
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
	     integer info, lwork
     	INTRINSIC INT, MIN
         ci=(0.d0,1.d0)
	     alci=al*ci
	     cial=alci
       open(5,file='oz.car')
        read(5,*)Jm
        read(5,*)Km
        read(5,*) al
        read(5,*) Re
        read(5,*) eps
        close(5)
      write(*,*) 'Jm=',Jm,'Km=', Km,'al=', al ,'Re=', Re
       call Com(eps)





300   format('  ZONE  I=',i4,'  J=',i4,'  F=POINT')
      open(4, File='streamfun.dat')
      write(4,300)Jm,Km
        do j = 1, Jm
            do k = 1, Km
                write(4,*) ym(j), zm(k), ub(j,k)
            end do
        end do
      close(4)


* Main loop
      t=0.

      Nm=(Jm-1)*Km+Jm*(Km-1)
*        do n=1,Nm

        n = 0

        do j = 1, Jm-1
            do k = 1 , Km
                n = n+1
                do j1 = 0 , Jm+1
                    do k1 = 0 , Km+1
                        v(j1,k1)=0
                        w(j1,k1)=0
                    end do
                end do
                v(j,k)=1
                call rp(t,u,v,w,u1,v1,w1,ox,oy,oz,p,Jmax, Kmax)

                m=0
                do j1 = 1, Jm-1
                    do k1 = 1, Km
                        m = m+1
                        amat(m,n) = v1(j1,k1)
                    end do
                end do

                do j1 = 1, Jm
                    do k1 = 1 , Km-1
                        m=m+1
                        amat(m,n) = w1(j1,k1)
                    end do
                end do
            end do
        end do

        do j = 1, Jm
            do k = 1 , Km-1
                n=n+1
                do j1 = 0 , Jm+1
                    do k1 = 0 , Km+1
                        v(j1,k1)=0
                        w(j1,k1)=0
                    end do
                end do
                w(j,k)=1
                call rp(t,u,v,w,u1,v1,w1,ox,oy,oz,p,Jmax, Kmax)
                m=0
                do j1 = 1, Jm-1
                    do k1 = 1, Km
                        m = m+1
                        amat(m,n) = v1(j1,k1)
                    end do
                end do
                do j1 = 1, Jm
                    do k1 = 1 , Km-1
                        m = m+1
                        amat(m,n) = w1(j1,k1)
                    end do
                end do
            end do
        end do



        do j = 1 , Nm
            do k = 1 , Nm
                amat1(j,k)= amat(j,k)
            end do
        end do

* Eigenvectors and eigenvalues
          write(*,*) 'ok'
          info = 2
          lwork=-1
          CALL ZGEEV( 'N', 'V', Nm, amat, Jmax*Kmax, eval
     >, evecl, Jmax*Kmax, evec, Jmax*Kmax, work, lwork, rwork, info )
          lwork = MIN( LWMAX, INT( work( 1 ) ) )
          CALL ZGEEV( 'N', 'V', Nm, amat, Jmax*Kmax, eval
     >, evecl, Jmax*Kmax, evec, Jmax*Kmax, work, lwork, rwork, info )
*

      do k=1,Nm-1
         do n=k+1,Nm
             if(dreal(eval(n)).gt.dreal(eval(k))) then
              eval(Nm+1)=eval(k)
              eval(k)=eval(n)
              eval(n)=eval(Nm+1)
              do j=1,Nm
                evec(j,Nm+1)=evec(j,k)
                evec(j,k)=evec(j,n)
                evec(j,n)=evec(j,Nm+1)
              end do
             end if
         end do
      end do
*
      tmp3 = 0
      do l = 1, Nm
          do j=1, Nm
              tmp1 = (0.d0, 0.d0)
              tmp2 = 0
              do  k = 1, Nm
                  tmp1 = tmp1 + amat1(j,k)*evec(k,l)
              end do
              tmp2 = max(abs(eval(l)*evec(j,l) - tmp1), tmp2)
          end do
              tmp3 = max(tmp2, tmp3)
      end do
      write(*,*) tmp3



*выписывание главной моды и нахождение давления
      n = 0
      do j = 1, Jm-1
        do k = 1, Km
            n = n+1
            v(j,k) = evec(n,1)

        end do
      end do

      do j = 1, Jm
        do k = 1 , Km-1
            n=n+1
            w(j,k) = evec(n,1)

        end do
       end do

      ci=(0.d0,1.d0)
      alci=al*ci
      cial=alci

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



      do j = 1 , Jm
        do k = 1, Km
            tmp1 = 0.5d0*(v(j-1,k)*ozb(j-1,k)+v(j,k)*ozb(j,k))
            tmp2 = 0.5d0*(vb(j-1,k)*oz(j-1,k)+vb(j,k)*oz(j,k))
            tmp3 = - 0.5d0*(wb(j,k-1)*oy(j,k-1)+wb(j,k)*oy(j,k))
            tmp4 = - 0.5d0*(w(j,k-1)*oyb(j,k-1)+w(j,k)*oyb(j,k-1))
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




      do j = 1, Jm
        do k = 1, Km
            ppert(j,k) = (ut(j,k)-eval(1)*u(j,k))/cial
        end do
      end do


      diff1 = 0.

      do j = 1, Jm-1
        do k = 1, Km
            tmp1 = vt(j,k)-v(j,k)*eval(1)
     >      -(ppert(j+1,k)-ppert(j,k))/yn1(j)
            write(*,*) abs(tmp1)
            diff1 = max(diff1, abs(tmp1))
        end do
      end do
      write(*,*) diff1

*


* выписывание массива собственных векторов
        open(6, FILE ='c.dat')
        do n=1,Nm
*            write(6,*) -dreal(eval(n)/(al*(0.,1.))),
*     >               -dimag(eval(n)/(al*(0.,1.)))
        end do
        do n=1, Nm
             write(6,*) dreal(eval(n)), dimag(eval(n))
        end do
        close(6)
*

      end
