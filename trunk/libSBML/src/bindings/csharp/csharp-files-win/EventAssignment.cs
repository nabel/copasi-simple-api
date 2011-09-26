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

public class EventAssignment : SBase {
	private HandleRef swigCPtr;
	
	internal EventAssignment(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.EventAssignmentUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.EventAssignmentUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(EventAssignment obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (EventAssignment obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~EventAssignment() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_EventAssignment(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public EventAssignment(long level, long version) : this(libsbmlPINVOKE.new_EventAssignment__SWIG_0(level, version), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public EventAssignment(SBMLNamespaces sbmlns) : this(libsbmlPINVOKE.new_EventAssignment__SWIG_1(SBMLNamespaces.getCPtr(sbmlns)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public EventAssignment(EventAssignment orig) : this(libsbmlPINVOKE.new_EventAssignment__SWIG_2(EventAssignment.getCPtr(orig)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public new EventAssignment clone() {
    IntPtr cPtr = libsbmlPINVOKE.EventAssignment_clone(swigCPtr);
    EventAssignment ret = (cPtr == IntPtr.Zero) ? null : new EventAssignment(cPtr, true);
    return ret;
  }

  public string getVariable() {
    string ret = libsbmlPINVOKE.EventAssignment_getVariable(swigCPtr);
    return ret;
  }

  public ASTNode getMath() {
    IntPtr cPtr = libsbmlPINVOKE.EventAssignment_getMath(swigCPtr);
    ASTNode ret = (cPtr == IntPtr.Zero) ? null : new ASTNode(cPtr, false);
    return ret;
  }

  public bool isSetVariable() {
    bool ret = libsbmlPINVOKE.EventAssignment_isSetVariable(swigCPtr);
    return ret;
  }

  public bool isSetMath() {
    bool ret = libsbmlPINVOKE.EventAssignment_isSetMath(swigCPtr);
    return ret;
  }

  public int setVariable(string sid) {
    int ret = libsbmlPINVOKE.EventAssignment_setVariable(swigCPtr, sid);
    return ret;
  }

  public int setMath(ASTNode math) {
    int ret = libsbmlPINVOKE.EventAssignment_setMath(swigCPtr, ASTNode.getCPtr(math));
    return ret;
  }

  public UnitDefinition getDerivedUnitDefinition() {
    IntPtr cPtr = libsbmlPINVOKE.EventAssignment_getDerivedUnitDefinition__SWIG_0(swigCPtr);
    UnitDefinition ret = (cPtr == IntPtr.Zero) ? null : new UnitDefinition(cPtr, false);
    return ret;
  }

  public bool containsUndeclaredUnits() {
    bool ret = libsbmlPINVOKE.EventAssignment_containsUndeclaredUnits__SWIG_0(swigCPtr);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.EventAssignment_getTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.EventAssignment_getElementName(swigCPtr);
    return ret;
  }

  public override bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.EventAssignment_hasRequiredAttributes(swigCPtr);
    return ret;
  }

  public override bool hasRequiredElements() {
    bool ret = libsbmlPINVOKE.EventAssignment_hasRequiredElements(swigCPtr);
    return ret;
  }

  public new string getId() {
    string ret = libsbmlPINVOKE.EventAssignment_getId(swigCPtr);
    return ret;
  }

}

}
