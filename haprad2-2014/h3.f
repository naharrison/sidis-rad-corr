      double precision function h3(x,q2,z)
      implicit none
      double precision x,q2,z
      double precision q0,lambda
      double precision a,a1,a2,b1,b2,bb
      data q0,lambda/1.d0,0.25d0/
      data a/-0.36544D-03/,a1/-2.1855d0/,a2/3.4176d0/,
     &b1/-1.7567d0/,b2/1.1272d0/,bb/8.9985d0/
      if(q2.gt.q0) then
      h3=a*x**a1*(1.d0-x)**a2*z**b1*(1.d0-z)**b2*
     &(log(q2/(lambda**2))/log(q0/(lambda**2)))**bb
      else
      h3=a*x**a1*(1.d0-x)**a2*z**b1*(1.d0-z)**b2
      endif
      return
      end
