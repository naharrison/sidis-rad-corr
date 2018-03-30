      double precision function h4(x,q2,z)
      implicit none
      double precision x,q2,z
      double precision q0,lambda
      double precision a,a1,a2,b1,b2,bb
      data q0,lambda/1.d0,0.25d0/
      data a/0.10908D-02/,a1/-0.35265D-06/,a2/0.30276D-07/,
     &b1/-0.66787d0/,b2/3.5868d0/,bb/6.8777d0/
      if(q2.gt.q0) then
      h4=a*x**a1*(1.d0-x)**a2*z**b1*(1.d0-z)**b2*
     &(log(q2/(lambda**2))/log(q0/(lambda**2)))**(bb/x)
      else
      h4=a*x**a1*(1.d0-x)**a2*z**b1*(1.d0-z)**b2
      endif
      return
      end
