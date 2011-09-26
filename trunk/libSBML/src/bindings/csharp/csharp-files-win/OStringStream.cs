/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

namespace libsbml {

using System;
using System.Runtime.InteropServices;

public class OStringStream : OStream {
	private HandleRef swigCPtr;
	
	internal OStringStream(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.OStringStreamUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.OStringStreamUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(OStringStream obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (OStringStream obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~OStringStream() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_OStringStream(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public OStringStream() : this(libsbmlPINVOKE.new_OStringStream(), true) {
  }

  public string str() {
    string ret = libsbmlPINVOKE.OStringStream_str__SWIG_0(swigCPtr);
    return ret;
  }

  public void str(string s) {
    libsbmlPINVOKE.OStringStream_str__SWIG_1(swigCPtr, s);
  }

}

}
