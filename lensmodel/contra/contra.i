%module contra
%{
// Headers from contra.c
#define nmax 1001
int pymain(int MACin, int BARin, double cin, double fbin, double rbin, double innerSlopein, double nuin, double Ain, double win, int nrad, double rfout[nmax], double rhofout[nmax]);
%}

// The arg list for pymain in contra has arrays meant to be filled as return values
// These typemaps eliminate the need to send these arrays as inputs from a python call
//   instead filling a blank array to send to the C code.
%typemap(in,numinputs=0) double rfout[ANY] (double temp[$1_dim0]) {
  int i;
  for (i = 0; i < $1_dim0; i++) {
    temp[i] = 0;
  }
  $1 = temp;
}
%typemap(in,numinputs=0) double rhofout[ANY] (double temp[$1_dim0]) {
  int i;
  for (i = 0; i < $1_dim0; i++) {
    temp[i] = 0;
  }
  $1 = temp;
}

// The arg list for pymain in contra has arrays meant to be filled as return values
// These typemaps add those arrays to the list of outputs returned by the C code
//   so they are received as lists by the python call.
%typemap(argout) (double rfout[nmax], double rhofout[nmax]) {
    PyObject *o = PyList_New($1_dim0+$2_dim0);
    int i;
    for(i=0; i<$1_dim0; i++)
    {
        PyList_SetItem(o, i, PyFloat_FromDouble($1[i]));
    }
    for(i=0; i<$2_dim0; i++)
    {
        PyList_SetItem(o, i+$1_dim0, PyFloat_FromDouble($2[i]));
    }
    $result = o;
}
int pymain(int MACin, int BARin, double cin, double fbin, double rbin, double innerSlopein, double nuin, double Ain, double win, int nrad, double rfout[nmax], double rhofout[nmax]);

%pythoncode{
  def cleanOutput(out,nrad):
     import numpy as np
     arr=np.array(out)
     return arr.reshape((2,len(arr)/2))[:,:nrad]

  def contra(MAC, BAR, conc, fb, rb, isl, nu, A, w, nrad):
     out=pymain(MAC,BAR,conc,fb,rb,isl,nu,A,w,nrad)
     return cleanOutput(out,nrad)
}
