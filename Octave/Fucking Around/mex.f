      #include "fintrf.h"
C======================================================================
#if 0
C     
C     timestwo.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     timestwo.f
C
C     Computational function that takes a scalar and doubles it.
      
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2018 The MathWorks, Inc.
C     
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:

#if MX_HAS_INTERLEAVED_COMPLEX
      mwPointer mxGetDoubles
#else
      mwPointer mxGetPr
#endif

      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN

C     Pointers to input/output mxArrays:
      mwPointer x_ptr, y_ptr

C     Array information:
      mwPointer mrows, ncols
      mwSize size

C     Arguments for computational routine:
      real*8  x_input, y_output

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput',
     +                           'One input required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nOutput',
     +                           'Too many output arguments.')
      endif

C     Validate inputs
C     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:NonNumeric',
     +                           'Input must be a number.')
      endif

C     Get the size of the input array.
      mrows = mxGetM(prhs(1))
      ncols = mxGetN(prhs(1))
      size = mrows*ncols

C     Create Fortran array from the input argument.

#if MX_HAS_INTERLEAVED_COMPLEX
      x_ptr = mxGetDoubles(prhs(1))
#else
      x_ptr = mxGetPr(prhs(1))
#endif
      call mxCopyPtrToReal8(x_ptr,x_input,size)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(mrows,ncols,0)
#if MX_HAS_INTERLEAVED_COMPLEX
      y_ptr = mxGetDoubles(plhs(1))
#else
      y_ptr = mxGetPr(plhs(1))
#endif

C     Call the computational subroutine.
      call timestwo(y_output, x_input)

C     Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(y_output,y_ptr,size)     

      return
      end

C-----------------------------------------------------------------------
      subroutine RK4(t,y,h,alpha,y0,v0,tf,n)
      implicit none
      real*8::t,y,h,f1,f2,K,alpha,y0,v0,tf
      dimension::t(0:1000),y(1:2,0:1000),K(1:2,1:4)
      integer::n,i,m
      n=tf/h
      do i=0,n
        t(i)= h * i
      end do
      y(1,0) = y0
      y(2,0)= v0
c     Note: \ddot{y} = -sin(y) - alpha*\dot{y}, i.e. \dot{y} = v and \dot{v} = -sin(y)- alpha *v
c    i.e, our diff equation will be: f1 = y2,f2 = -sin(y1) - alpha y2
      do i=0,n-1
        
      K(1,1)=f1(t(i),y(1,i),y(2,i),alpha)
      K(2,1)=f2(t(i),y(1,i),y(2,i),alpha)
      K(1,2) = f1(t(i)+0.5*h,y(1,i)+0.5*h*K(1,1),y(2,i)+0.5*h*K(2,1)
     c ,alpha)
      K(2,2) = f2(t(i)+0.5*h,y(1,i)+0.5*h*K(1,1),y(2,i)+0.5*h*K(2,1)
     c,alpha)
      K(1,3) = f1(t(i)+0.5*h,y(1,i)+0.5*h*K(1,2),y(2,i)+0.5*h*K(2,2)
      
     c,alpha)
      K(2,3) = f2(t(i)+0.5*h,y(1,i)+0.5*h*K(1,2),y(2,i)+0.5*h*K(2,2)
     c,alpha)
      K(1,4) = f1(t(i)+0.5*h,y(1,i)+h*K(1,3),y(2,i)+h*K(2,3),alpha)
      K(2,4) = f2(t(i)+0.5*h,y(1,i)+h*K(1,3),y(2,i)+h*K(2,3),alpha)
      y(1,i+1) = y(1,i) +(h/6)*(K(1,1)+2*K(1,2)+2*K(1,3)+K(1,4)) 
      y(2,i+1) = y(2,i) +(h/6)*(K(2,1)+2*K(2,2)+2*K(2,3)+K(2,4)) 
      end do
      end subroutine
      
      function f1(t,y,v,alpha)
      implicit none
      real*8::t,y,v,alpha,f1
      f1 = v
      end function
      
      function f2(t,y,v,alpha)
      implicit none
      real*8::t,y,v,f2,alpha
      f2 = -sin(y)- alpha*v
      end function
