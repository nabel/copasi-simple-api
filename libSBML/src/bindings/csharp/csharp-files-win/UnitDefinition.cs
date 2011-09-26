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

public class UnitDefinition : SBase {
	private HandleRef swigCPtr;
	
	internal UnitDefinition(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.UnitDefinitionUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.UnitDefinitionUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(UnitDefinition obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (UnitDefinition obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~UnitDefinition() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_UnitDefinition(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public UnitDefinition(long level, long version) : this(libsbmlPINVOKE.new_UnitDefinition__SWIG_0(level, version), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public UnitDefinition(SBMLNamespaces sbmlns) : this(libsbmlPINVOKE.new_UnitDefinition__SWIG_1(SBMLNamespaces.getCPtr(sbmlns)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public UnitDefinition(UnitDefinition orig) : this(libsbmlPINVOKE.new_UnitDefinition__SWIG_2(UnitDefinition.getCPtr(orig)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public new UnitDefinition clone() {
    IntPtr cPtr = libsbmlPINVOKE.UnitDefinition_clone(swigCPtr);
    UnitDefinition ret = (cPtr == IntPtr.Zero) ? null : new UnitDefinition(cPtr, true);
    return ret;
  }

  public new string getId() {
    string ret = libsbmlPINVOKE.UnitDefinition_getId(swigCPtr);
    return ret;
  }

  public new string getName() {
    string ret = libsbmlPINVOKE.UnitDefinition_getName(swigCPtr);
    return ret;
  }

  public new bool isSetId() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isSetId(swigCPtr);
    return ret;
  }

  public new bool isSetName() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isSetName(swigCPtr);
    return ret;
  }

  public new int setId(string sid) {
    int ret = libsbmlPINVOKE.UnitDefinition_setId(swigCPtr, sid);
    return ret;
  }

  public new int setName(string name) {
    int ret = libsbmlPINVOKE.UnitDefinition_setName(swigCPtr, name);
    return ret;
  }

  public int unsetName() {
    int ret = libsbmlPINVOKE.UnitDefinition_unsetName(swigCPtr);
    return ret;
  }

  public bool isVariantOfArea() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfArea(swigCPtr);
    return ret;
  }

  public bool isVariantOfLength() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfLength(swigCPtr);
    return ret;
  }

  public bool isVariantOfSubstance() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfSubstance(swigCPtr);
    return ret;
  }

  public bool isVariantOfTime() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfTime(swigCPtr);
    return ret;
  }

  public bool isVariantOfVolume() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfVolume(swigCPtr);
    return ret;
  }

  public bool isVariantOfDimensionless() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfDimensionless(swigCPtr);
    return ret;
  }

  public bool isVariantOfMass() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfMass(swigCPtr);
    return ret;
  }

  public bool isVariantOfSubstancePerTime() {
    bool ret = libsbmlPINVOKE.UnitDefinition_isVariantOfSubstancePerTime(swigCPtr);
    return ret;
  }

  public int addUnit(Unit u) {
    int ret = libsbmlPINVOKE.UnitDefinition_addUnit(swigCPtr, Unit.getCPtr(u));
    return ret;
  }

  public Unit createUnit() {
    IntPtr cPtr = libsbmlPINVOKE.UnitDefinition_createUnit(swigCPtr);
    Unit ret = (cPtr == IntPtr.Zero) ? null : new Unit(cPtr, false);
    return ret;
  }

  public ListOfUnits getListOfUnits() {
    IntPtr cPtr = libsbmlPINVOKE.UnitDefinition_getListOfUnits__SWIG_0(swigCPtr);
    ListOfUnits ret = (cPtr == IntPtr.Zero) ? null : new ListOfUnits(cPtr, false);
    return ret;
  }

  public Unit getUnit(long n) {
    IntPtr cPtr = libsbmlPINVOKE.UnitDefinition_getUnit__SWIG_0(swigCPtr, n);
    Unit ret = (cPtr == IntPtr.Zero) ? null : new Unit(cPtr, false);
    return ret;
  }

  public long getNumUnits() { return (long)libsbmlPINVOKE.UnitDefinition_getNumUnits(swigCPtr); }

  public Unit removeUnit(long n) {
    IntPtr cPtr = libsbmlPINVOKE.UnitDefinition_removeUnit(swigCPtr, n);
    Unit ret = (cPtr == IntPtr.Zero) ? null : new Unit(cPtr, true);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.UnitDefinition_getTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.UnitDefinition_getElementName(swigCPtr);
    return ret;
  }

  public static void simplify(UnitDefinition ud) {
    libsbmlPINVOKE.UnitDefinition_simplify(UnitDefinition.getCPtr(ud));
  }

  public static void reorder(UnitDefinition ud) {
    libsbmlPINVOKE.UnitDefinition_reorder(UnitDefinition.getCPtr(ud));
  }

  public static UnitDefinition convertToSI(UnitDefinition arg0) {
    IntPtr cPtr = libsbmlPINVOKE.UnitDefinition_convertToSI(UnitDefinition.getCPtr(arg0));
    UnitDefinition ret = (cPtr == IntPtr.Zero) ? null : new UnitDefinition(cPtr, true);
    return ret;
  }

  public static bool areIdentical(UnitDefinition ud1, UnitDefinition ud2) {
    bool ret = libsbmlPINVOKE.UnitDefinition_areIdentical(UnitDefinition.getCPtr(ud1), UnitDefinition.getCPtr(ud2));
    return ret;
  }

  public static bool areEquivalent(UnitDefinition ud1, UnitDefinition ud2) {
    bool ret = libsbmlPINVOKE.UnitDefinition_areEquivalent(UnitDefinition.getCPtr(ud1), UnitDefinition.getCPtr(ud2));
    return ret;
  }

  public static UnitDefinition combine(UnitDefinition ud1, UnitDefinition ud2) {
    IntPtr cPtr = libsbmlPINVOKE.UnitDefinition_combine(UnitDefinition.getCPtr(ud1), UnitDefinition.getCPtr(ud2));
    UnitDefinition ret = (cPtr == IntPtr.Zero) ? null : new UnitDefinition(cPtr, true);
    return ret;
  }

  public static string printUnits(UnitDefinition ud, bool compact) {
    string ret = libsbmlPINVOKE.UnitDefinition_printUnits__SWIG_0(UnitDefinition.getCPtr(ud), compact);
    return ret;
  }

  public static string printUnits(UnitDefinition ud) {
    string ret = libsbmlPINVOKE.UnitDefinition_printUnits__SWIG_1(UnitDefinition.getCPtr(ud));
    return ret;
  }

  public override bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.UnitDefinition_hasRequiredAttributes(swigCPtr);
    return ret;
  }

  public override bool hasRequiredElements() {
    bool ret = libsbmlPINVOKE.UnitDefinition_hasRequiredElements(swigCPtr);
    return ret;
  }

}

}
