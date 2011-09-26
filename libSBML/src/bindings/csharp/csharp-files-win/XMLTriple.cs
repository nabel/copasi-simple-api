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

public class XMLTriple : IDisposable {
	private HandleRef swigCPtr;
	protected bool swigCMemOwn;
	
	internal XMLTriple(IntPtr cPtr, bool cMemoryOwn)
	{
		swigCMemOwn = cMemoryOwn;
		swigCPtr    = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(XMLTriple obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (XMLTriple obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~XMLTriple() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_XMLTriple(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
    }
  }

  public XMLTriple() : this(libsbmlPINVOKE.new_XMLTriple__SWIG_0(), true) {
  }

  public XMLTriple(string name, string uri, string prefix) : this(libsbmlPINVOKE.new_XMLTriple__SWIG_1(name, uri, prefix), true) {
  }

  public XMLTriple(string triplet, char sepchar) : this(libsbmlPINVOKE.new_XMLTriple__SWIG_2(triplet, sepchar), true) {
  }

  public XMLTriple(string triplet) : this(libsbmlPINVOKE.new_XMLTriple__SWIG_3(triplet), true) {
  }

  public XMLTriple(XMLTriple orig) : this(libsbmlPINVOKE.new_XMLTriple__SWIG_4(XMLTriple.getCPtr(orig)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public XMLTriple clone() {
    IntPtr cPtr = libsbmlPINVOKE.XMLTriple_clone(swigCPtr);
    XMLTriple ret = (cPtr == IntPtr.Zero) ? null : new XMLTriple(cPtr, true);
    return ret;
  }

  public string getName() {
    string ret = libsbmlPINVOKE.XMLTriple_getName(swigCPtr);
    return ret;
  }

  public string getPrefix() {
    string ret = libsbmlPINVOKE.XMLTriple_getPrefix(swigCPtr);
    return ret;
  }

  public string getURI() {
    string ret = libsbmlPINVOKE.XMLTriple_getURI(swigCPtr);
    return ret;
  }

  public string getPrefixedName() {
    string ret = libsbmlPINVOKE.XMLTriple_getPrefixedName(swigCPtr);
    return ret;
  }

  public bool isEmpty() {
    bool ret = libsbmlPINVOKE.XMLTriple_isEmpty(swigCPtr);
    return ret;
  }

}

}
