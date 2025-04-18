      subroutine roptim(n,S,L,thr,maxIt,msg,warm,X,W,Wd,WXj,info)
!
!     .. Scalar Arguments ..
      implicit double precision(a-h, o-z)
      integer n, warm, msg, maxIt
!     ..
!     .. Array Arguments ..
      double precision S(n,n), L(n,n), X(n,n), W(n,n), Wd(n), WXj(n)
!     ..
!
!  Purpose
!  =======
!
!  This subroutine computes the L1 regularized covariance matrix estimate
!  using the algorithm described in the paper:
!    J. Friedman, T. Hastie, R. Tibshirani:
!    Sparse inverse covariance estimation with the graphical lasso
!    Biostatistics, 9(3):432-441, July 2008.
!  This implementation is documented in the technical report:
!    M. A. Sustik B. Calderhead:
!    GLASSOFAST: An efficient GLASSO implementation
!    Technical Report TR-12-29, University of Texas at Austin
!
!  Arguments
!  =========
!
!  n      (input) integer
!         The dimension of the input matrix s
!
!  S      (input) double precision array, dimension n x n
!         The empirical covariance matrix
!
!  L      (input) double precision array, dimension n x n
!         Regularization matrix (symmetric); U NAS off-diagonal ze stala wartoscia
!
!  thr    (input) double precision
!         Convergence threshold
!
!  maxIt  (input/output) integer
!         Maximum number of whole matrix sweeps on input, actual number
!         of sweeps on output
!
!  msg    (input) integer
!         Controls amount of messages printed
!
!  warm   (input) integer
!         flag indicating cold versus warm start, see also X, W
!
!  X      (input/output) double precision array, dimension n x n
!         Inverse covariance matrix estimate
!
!  W      (input/output) double precision array, dimension n x n
!         Covariance matrix estimate
!
!  Wd     (input) double precision array, dimension n
!         Internally used workspace
!
!  WXj    (input) double precision array, dimension n
!         Internally used workspace
!
!  info   (output) integer
!         Indicates errors
!
!  Wd, WDXj are not allocatable arrays, because gfortran/gdb has trouble
!  displaying those.

integer iter
double precision EPS
parameter (EPS = 1.1e-16)
info = 0
shr = sum(abs(S))
do i = 1,n
   shr = shr - abs(S(i, i))
enddo
if (shr .eq. 0.0) then
!  S is diagonal.
   W = 0.0
   do i = 1,n
      W(i,i) = 1.0
   enddo
   X = 0.0
   do i = 1,n
      X(i,i) = 1.0
   enddo
   return
endif
shr = thr*shr/(n-1)
thrLasso = shr/n
if (thrLasso .lt. 2*EPS) then
   thrLasso = 2*EPS
end if
if (warm .eq. 0) then
   X = 0.0 ! X to kwadracik
   W = 0.0 ! TODO(Sprawdzic, czy nie lepiej z innego)
   do i = 1,n
      W(i,i) = 1.0
   enddo
else
   do i = 1,n
     X(1:n,i) = -X(1:n,i)/X(i,i) ! To wazne -- X zmienia interpretacje; teraz to bedzie kwadracik
     X(i,i) = 0
  end do
end if
!do i = 1,n ! My tego nie mamy
!   Wd(i) = S(i,i) + L(i,i)
!   W(i,i) = Wd(i)
!end do
do iter = 1,maxIt
! if (msg .ne. 0) write(6,*) "iteration =", iter
!   print iterations to the console (ok with CRAN)
if (msg .ne. 0)  call intpr('iter:',-1,iter,1)
   dw = 0.0
   do j = 1,n
      WXj(1:n) = 0.0
!     We exploit sparsity of X when computing column j of W*X*D:
      do i = 1,n
         if (X(i,j) .ne. 0.0) then
            WXj = WXj + W(:,i)*X(i,j)
         endif
      enddo
      do
         dlx = 0.0
         do i = 1,n
            if (i .ne. j) then !ten kawałek to jest zwykly lasso
               a = S(i,j) - WXj(i) + W(i,i)*X(i,j)
               b = abs(a) - L(i,j)
               if (b .gt. 0.0) then
                  c = sign(b, a)/W(i,i)
               else
                  c = 0.0
               endif
               delta = c - X(i,j)
               if (delta .ne. 0.0) then
                  X(i,j) = c
                  WXj(1:n) = WXj(1:n) + W(:,i)*delta
                  dlx = max(dlx, abs(delta))
               endif
            endif
         enddo
         if (dlx .lt. thrLasso) then
            exit
         endif
      enddo
      WXj(j) = W(j,j)
      dw = max(dw, sum(abs(WXj(1:n) - W(:,j))))
      W(:,j) = WXj(1:n)
      W(j,:) = WXj(1:n)
      W(j,j) = 1 + sum(X(:,j)*W(:,j))
   enddo
!   write(6,*) "  dw =", dw
   if (dw .le. shr) then
      exit
   endif
enddo
do i = 1,n
   X(1:n,i) = -X(1:n,i)
   X(i,i) = 1
enddo
do i = 1,n-1
   X(i+1:n,i) = (X(i+1:n,i) + X(i,i+1:n))/2;
   X(i,i+1:n) = X(i+1:n,i)
enddo
maxIt = iter
return
end subroutine roptim
