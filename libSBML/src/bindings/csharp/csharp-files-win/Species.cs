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

public class Species : SBase {
	private HandleRef swigCPtr;
	
	internal Species(IntPtr cPtr, bool cMemoryOwn) : base(libsbmlPINVOKE.SpeciesUpcast(cPtr), cMemoryOwn)
	{
		//super(libsbmlPINVOKE.SpeciesUpcast(cPtr), cMemoryOwn);
		swigCPtr = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(Species obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (Species obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~Species() {
    Dispose();
  }

  public override void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_Species(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
      base.Dispose();
    }
  }

  public Species(long level, long version) : this(libsbmlPINVOKE.new_Species__SWIG_0(level, version), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public Species(SBMLNamespaces sbmlns) : this(libsbmlPINVOKE.new_Species__SWIG_1(SBMLNamespaces.getCPtr(sbmlns)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public Species(Species orig) : this(libsbmlPINVOKE.new_Species__SWIG_2(Species.getCPtr(orig)), true) {
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
  }

  public new Species clone() {
    IntPtr cPtr = libsbmlPINVOKE.Species_clone(swigCPtr);
    Species ret = (cPtr == IntPtr.Zero) ? null : new Species(cPtr, true);
    return ret;
  }

  public void initDefaults() {
    libsbmlPINVOKE.Species_initDefaults(swigCPtr);
  }

  public new string getId() {
    string ret = libsbmlPINVOKE.Species_getId(swigCPtr);
    return ret;
  }

  public new string getName() {
    string ret = libsbmlPINVOKE.Species_getName(swigCPtr);
    return ret;
  }

  public string getSpeciesType() {
    string ret = libsbmlPINVOKE.Species_getSpeciesType(swigCPtr);
    return ret;
  }

  public string getCompartment() {
    string ret = libsbmlPINVOKE.Species_getCompartment(swigCPtr);
    return ret;
  }

  public double getInitialAmount() {
    double ret = libsbmlPINVOKE.Species_getInitialAmount(swigCPtr);
    return ret;
  }

  public double getInitialConcentration() {
    double ret = libsbmlPINVOKE.Species_getInitialConcentration(swigCPtr);
    return ret;
  }

  public string getSubstanceUnits() {
    string ret = libsbmlPINVOKE.Species_getSubstanceUnits(swigCPtr);
    return ret;
  }

  public string getSpatialSizeUnits() {
    string ret = libsbmlPINVOKE.Species_getSpatialSizeUnits(swigCPtr);
    return ret;
  }

  public string getUnits() {
    string ret = libsbmlPINVOKE.Species_getUnits(swigCPtr);
    return ret;
  }

  public bool getHasOnlySubstanceUnits() {
    bool ret = libsbmlPINVOKE.Species_getHasOnlySubstanceUnits(swigCPtr);
    return ret;
  }

  public bool getBoundaryCondition() {
    bool ret = libsbmlPINVOKE.Species_getBoundaryCondition(swigCPtr);
    return ret;
  }

  public int getCharge() {
    int ret = libsbmlPINVOKE.Species_getCharge(swigCPtr);
    return ret;
  }

  public bool getConstant() {
    bool ret = libsbmlPINVOKE.Species_getConstant(swigCPtr);
    return ret;
  }

  public string getConversionFactor() {
    string ret = libsbmlPINVOKE.Species_getConversionFactor(swigCPtr);
    return ret;
  }

  public new bool isSetId() {
    bool ret = libsbmlPINVOKE.Species_isSetId(swigCPtr);
    return ret;
  }

  public new bool isSetName() {
    bool ret = libsbmlPINVOKE.Species_isSetName(swigCPtr);
    return ret;
  }

  public bool isSetSpeciesType() {
    bool ret = libsbmlPINVOKE.Species_isSetSpeciesType(swigCPtr);
    return ret;
  }

  public bool isSetCompartment() {
    bool ret = libsbmlPINVOKE.Species_isSetCompartment(swigCPtr);
    return ret;
  }

  public bool isSetInitialAmount() {
    bool ret = libsbmlPINVOKE.Species_isSetInitialAmount(swigCPtr);
    return ret;
  }

  public bool isSetInitialConcentration() {
    bool ret = libsbmlPINVOKE.Species_isSetInitialConcentration(swigCPtr);
    return ret;
  }

  public bool isSetSubstanceUnits() {
    bool ret = libsbmlPINVOKE.Species_isSetSubstanceUnits(swigCPtr);
    return ret;
  }

  public bool isSetSpatialSizeUnits() {
    bool ret = libsbmlPINVOKE.Species_isSetSpatialSizeUnits(swigCPtr);
    return ret;
  }

  public bool isSetUnits() {
    bool ret = libsbmlPINVOKE.Species_isSetUnits(swigCPtr);
    return ret;
  }

  public bool isSetCharge() {
    bool ret = libsbmlPINVOKE.Species_isSetCharge(swigCPtr);
    return ret;
  }

  public bool isSetConversionFactor() {
    bool ret = libsbmlPINVOKE.Species_isSetConversionFactor(swigCPtr);
    return ret;
  }

  public bool isSetBoundaryCondition() {
    bool ret = libsbmlPINVOKE.Species_isSetBoundaryCondition(swigCPtr);
    return ret;
  }

  public bool isSetHasOnlySubstanceUnits() {
    bool ret = libsbmlPINVOKE.Species_isSetHasOnlySubstanceUnits(swigCPtr);
    return ret;
  }

  public bool isSetConstant() {
    bool ret = libsbmlPINVOKE.Species_isSetConstant(swigCPtr);
    return ret;
  }

  public new int setId(string sid) {
    int ret = libsbmlPINVOKE.Species_setId(swigCPtr, sid);
    return ret;
  }

  public new int setName(string name) {
    int ret = libsbmlPINVOKE.Species_setName(swigCPtr, name);
    return ret;
  }

  public int setSpeciesType(string sid) {
    int ret = libsbmlPINVOKE.Species_setSpeciesType(swigCPtr, sid);
    return ret;
  }

  public int setCompartment(string sid) {
    int ret = libsbmlPINVOKE.Species_setCompartment(swigCPtr, sid);
    return ret;
  }

  public int setInitialAmount(double value) {
    int ret = libsbmlPINVOKE.Species_setInitialAmount(swigCPtr, value);
    return ret;
  }

  public int setInitialConcentration(double value) {
    int ret = libsbmlPINVOKE.Species_setInitialConcentration(swigCPtr, value);
    return ret;
  }

  public int setSubstanceUnits(string sid) {
    int ret = libsbmlPINVOKE.Species_setSubstanceUnits(swigCPtr, sid);
    return ret;
  }

  public int setSpatialSizeUnits(string sid) {
    int ret = libsbmlPINVOKE.Species_setSpatialSizeUnits(swigCPtr, sid);
    return ret;
  }

  public int setUnits(string sname) {
    int ret = libsbmlPINVOKE.Species_setUnits(swigCPtr, sname);
    return ret;
  }

  public int setHasOnlySubstanceUnits(bool value) {
    int ret = libsbmlPINVOKE.Species_setHasOnlySubstanceUnits(swigCPtr, value);
    return ret;
  }

  public int setBoundaryCondition(bool value) {
    int ret = libsbmlPINVOKE.Species_setBoundaryCondition(swigCPtr, value);
    return ret;
  }

  public int setCharge(int value) {
    int ret = libsbmlPINVOKE.Species_setCharge(swigCPtr, value);
    return ret;
  }

  public int setConstant(bool value) {
    int ret = libsbmlPINVOKE.Species_setConstant(swigCPtr, value);
    return ret;
  }

  public int setConversionFactor(string sid) {
    int ret = libsbmlPINVOKE.Species_setConversionFactor(swigCPtr, sid);
    return ret;
  }

  public int unsetName() {
    int ret = libsbmlPINVOKE.Species_unsetName(swigCPtr);
    return ret;
  }

  public int unsetSpeciesType() {
    int ret = libsbmlPINVOKE.Species_unsetSpeciesType(swigCPtr);
    return ret;
  }

  public int unsetInitialAmount() {
    int ret = libsbmlPINVOKE.Species_unsetInitialAmount(swigCPtr);
    return ret;
  }

  public int unsetInitialConcentration() {
    int ret = libsbmlPINVOKE.Species_unsetInitialConcentration(swigCPtr);
    return ret;
  }

  public int unsetSubstanceUnits() {
    int ret = libsbmlPINVOKE.Species_unsetSubstanceUnits(swigCPtr);
    return ret;
  }

  public int unsetSpatialSizeUnits() {
    int ret = libsbmlPINVOKE.Species_unsetSpatialSizeUnits(swigCPtr);
    return ret;
  }

  public int unsetUnits() {
    int ret = libsbmlPINVOKE.Species_unsetUnits(swigCPtr);
    return ret;
  }

  public int unsetCharge() {
    int ret = libsbmlPINVOKE.Species_unsetCharge(swigCPtr);
    return ret;
  }

  public int unsetConversionFactor() {
    int ret = libsbmlPINVOKE.Species_unsetConversionFactor(swigCPtr);
    return ret;
  }

  public UnitDefinition getDerivedUnitDefinition() {
    IntPtr cPtr = libsbmlPINVOKE.Species_getDerivedUnitDefinition__SWIG_0(swigCPtr);
    UnitDefinition ret = (cPtr == IntPtr.Zero) ? null : new UnitDefinition(cPtr, false);
    return ret;
  }

  public override int getTypeCode() {
    int ret = libsbmlPINVOKE.Species_getTypeCode(swigCPtr);
    return ret;
  }

  public override string getElementName() {
    string ret = libsbmlPINVOKE.Species_getElementName(swigCPtr);
    return ret;
  }

  public override bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.Species_hasRequiredAttributes(swigCPtr);
    return ret;
  }

}

}
