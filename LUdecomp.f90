subroutine LUDecomp()
    use globalConstants

    ! nanoHUB.org
    ! Numerical Solution of Poisson's Equation
    ! Dragica Vasileska
    ! LU Decomposition Slides

    !do i=2, nmax
    !    if(a(i).eq.-c(i-1)) then
    !    a(i)=a(i)+0.01
    !    endif
    !enddo
    ! a, b, c, f are defined coming into LU Decomp
    alpha(1) = a(1)
    g(1) = f(1)

    ! From 2 to nmax by step of 1
    ! Evaluate expressions for beta, alpha, and g
    ! Obtain beta, alpha, g from 1 to nmax

    do i = 2, nmax
       if (b(i).eq.0) then
          beta(i)=0
       else
          beta(i)  = b(i) / (alpha(i-1))
       endif
       alpha(i) = a(i) - (beta(i)*c(i-1))
       g(i)     = f(i) - (beta(i)*g(i-1))

       if (isnan(beta(i))) then
          print *, "Error: NaN problem with beta(", i, ")"
          print *, "alpha(i-1)=", alpha(i-1)
          print *, "b(i)=", b(i)
       endif

       if (alpha(i).eq.0) then
          print *, "Error: alpha is 0"
          print *, "a(",i,")=", a(i)
          print *, "beta(",i,")=", beta(i)
          print *, "c(",i-1,")=", c(i-1)
          print *, "b(",i,")=", b(i)
          print *, "alpha(",i-1,")=", alpha(i-1)
          print *, "psi(",i+1,")=", psi(i+1)
          print *, "psi(",i,")=", psi(i)
          print *, "psi(",i-1,")=", psi(i-1)
          print *, "AA=", (psi(i+1)-psi(i))
          print *, "CC=", (psi(i-1)-psi(i))
          print *, "ZZ=", (psi(i)-psi(i+1))
          print *, "YY=", (psi(i)-psi(i-1))
       endif

       if (a(i).eq.0) then
          print *, "Error: a is 0"
       endif
    enddo

    ! Evaluate x at nmax
    ! Use alpha and g at nmax from previous do loop
    x(nmax)=g(nmax)/alpha(nmax)

    ! From n-max-1 to 1 by step of -1
    ! Evaluate expression for x
    do i = nmax-1, 1, -1
        x(i)=(g(i)-(c(i)*x(i+1)))/alpha(i)
     enddo

end subroutine LUDecomp