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

public class SBase : IDisposable {
	private HandleRef swigCPtr;
	protected bool swigCMemOwn;
	
	internal SBase(IntPtr cPtr, bool cMemoryOwn)
	{
		swigCMemOwn = cMemoryOwn;
		swigCPtr    = new HandleRef(this, cPtr);
	}
	
	internal static HandleRef getCPtr(SBase obj)
	{
		return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
	}
	
	internal static HandleRef getCPtrAndDisown (SBase obj)
	{
		HandleRef ptr = new HandleRef(null, IntPtr.Zero);
		
		if (obj != null)
		{
			ptr             = obj.swigCPtr;
			obj.swigCMemOwn = false;
		}
		
		return ptr;
	}

  ~SBase() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libsbmlPINVOKE.delete_SBase(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
    }
  }

  public static bool operator==(SBase lhs, SBase rhs)
  {
    if((Object)lhs == (Object)rhs)
    {
      return true;
    }

    if( ((Object)lhs == null) || ((Object)rhs == null) )
    {
      return false;
    }

    return (getCPtr(lhs).Handle.ToString() == getCPtr(rhs).Handle.ToString());
  }

  public static bool operator!=(SBase lhs, SBase rhs)
  {
    return !(lhs == rhs);
  }

  public override bool Equals(Object sb)
  {
    if ( ! (sb is SBase) )
    {
      return false;
    }

    return this == (SBase)sb;
  }

  public override int GetHashCode()
  {
    return swigCPtr.Handle.ToInt32();
  }

  public virtual SBase clone() {
	return libsbml.DowncastSBase(libsbmlPINVOKE.SBase_clone(swigCPtr), true);
}

  public string getMetaId() {
    string ret = libsbmlPINVOKE.SBase_getMetaId(swigCPtr);
    return ret;
  }

  public string getId() {
    string ret = libsbmlPINVOKE.SBase_getId(swigCPtr);
    return ret;
  }

  public string getName() {
    string ret = libsbmlPINVOKE.SBase_getName(swigCPtr);
    return ret;
  }

  public XMLNode getNotes() {
    IntPtr cPtr = libsbmlPINVOKE.SBase_getNotes(swigCPtr);
    XMLNode ret = (cPtr == IntPtr.Zero) ? null : new XMLNode(cPtr, false);
    return ret;
  }

  public string getNotesString() {
    string ret = libsbmlPINVOKE.SBase_getNotesString(swigCPtr);
    return ret;
  }

  public XMLNode getAnnotation() {
    IntPtr cPtr = libsbmlPINVOKE.SBase_getAnnotation(swigCPtr);
    XMLNode ret = (cPtr == IntPtr.Zero) ? null : new XMLNode(cPtr, false);
    return ret;
  }

  public string getAnnotationString() {
    string ret = libsbmlPINVOKE.SBase_getAnnotationString(swigCPtr);
    return ret;
  }

  public virtual XMLNamespaces getNamespaces() {
    IntPtr cPtr = libsbmlPINVOKE.SBase_getNamespaces(swigCPtr);
    XMLNamespaces ret = (cPtr == IntPtr.Zero) ? null : new XMLNamespaces(cPtr, false);
    return ret;
  }

  public SBMLDocument getSBMLDocument() {
    IntPtr cPtr = libsbmlPINVOKE.SBase_getSBMLDocument__SWIG_0(swigCPtr);
    SBMLDocument ret = (cPtr == IntPtr.Zero) ? null : new SBMLDocument(cPtr, false);
    return ret;
  }

  public SBase getParentSBMLObject() {
	return libsbml.DowncastSBase(libsbmlPINVOKE.SBase_getParentSBMLObject(swigCPtr), false);
}

  public SBase getAncestorOfType(int type) {
	return libsbml.DowncastSBase(libsbmlPINVOKE.SBase_getAncestorOfType(swigCPtr, type), false);
}

  public int getSBOTerm() {
    int ret = libsbmlPINVOKE.SBase_getSBOTerm(swigCPtr);
    return ret;
  }

  public string getSBOTermID() {
    string ret = libsbmlPINVOKE.SBase_getSBOTermID(swigCPtr);
    return ret;
  }

  public long getLine() { return (long)libsbmlPINVOKE.SBase_getLine(swigCPtr); }

  public long getColumn() { return (long)libsbmlPINVOKE.SBase_getColumn(swigCPtr); }

  public ModelHistory getModelHistory() {
    IntPtr cPtr = libsbmlPINVOKE.SBase_getModelHistory__SWIG_0(swigCPtr);
    ModelHistory ret = (cPtr == IntPtr.Zero) ? null : new ModelHistory(cPtr, false);
    return ret;
  }

  public bool isSetMetaId() {
    bool ret = libsbmlPINVOKE.SBase_isSetMetaId(swigCPtr);
    return ret;
  }

  public bool isSetId() {
    bool ret = libsbmlPINVOKE.SBase_isSetId(swigCPtr);
    return ret;
  }

  public bool isSetName() {
    bool ret = libsbmlPINVOKE.SBase_isSetName(swigCPtr);
    return ret;
  }

  public bool isSetNotes() {
    bool ret = libsbmlPINVOKE.SBase_isSetNotes(swigCPtr);
    return ret;
  }

  public bool isSetAnnotation() {
    bool ret = libsbmlPINVOKE.SBase_isSetAnnotation(swigCPtr);
    return ret;
  }

  public bool isSetSBOTerm() {
    bool ret = libsbmlPINVOKE.SBase_isSetSBOTerm(swigCPtr);
    return ret;
  }

  public int setMetaId(string metaid) {
    int ret = libsbmlPINVOKE.SBase_setMetaId(swigCPtr, metaid);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public bool isSetModelHistory() {
    bool ret = libsbmlPINVOKE.SBase_isSetModelHistory(swigCPtr);
    return ret;
  }

  public int setId(string sid) {
    int ret = libsbmlPINVOKE.SBase_setId(swigCPtr, sid);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int setName(string name) {
    int ret = libsbmlPINVOKE.SBase_setName(swigCPtr, name);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public virtual int setAnnotation(XMLNode annotation) {
    int ret = libsbmlPINVOKE.SBase_setAnnotation__SWIG_0(swigCPtr, XMLNode.getCPtr(annotation));
    return ret;
  }

  public virtual int setAnnotation(string annotation) {
    int ret = libsbmlPINVOKE.SBase_setAnnotation__SWIG_1(swigCPtr, annotation);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public virtual int appendAnnotation(XMLNode annotation) {
    int ret = libsbmlPINVOKE.SBase_appendAnnotation__SWIG_0(swigCPtr, XMLNode.getCPtr(annotation));
    return ret;
  }

  public virtual int appendAnnotation(string annotation) {
    int ret = libsbmlPINVOKE.SBase_appendAnnotation__SWIG_1(swigCPtr, annotation);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int setNotes(XMLNode notes) {
    int ret = libsbmlPINVOKE.SBase_setNotes__SWIG_0(swigCPtr, XMLNode.getCPtr(notes));
    return ret;
  }

  public int setNotes(string notes) {
    int ret = libsbmlPINVOKE.SBase_setNotes__SWIG_1(swigCPtr, notes);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int appendNotes(XMLNode notes) {
    int ret = libsbmlPINVOKE.SBase_appendNotes__SWIG_0(swigCPtr, XMLNode.getCPtr(notes));
    return ret;
  }

  public int appendNotes(string notes) {
    int ret = libsbmlPINVOKE.SBase_appendNotes__SWIG_1(swigCPtr, notes);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int setModelHistory(ModelHistory history) {
    int ret = libsbmlPINVOKE.SBase_setModelHistory(swigCPtr, ModelHistory.getCPtr(history));
    return ret;
  }

  public virtual int setSBOTerm(int value) {
    int ret = libsbmlPINVOKE.SBase_setSBOTerm__SWIG_0(swigCPtr, value);
    return ret;
  }

  public virtual int setSBOTerm(string sboid) {
    int ret = libsbmlPINVOKE.SBase_setSBOTerm__SWIG_1(swigCPtr, sboid);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int setNamespaces(XMLNamespaces xmlns) {
    int ret = libsbmlPINVOKE.SBase_setNamespaces(swigCPtr, XMLNamespaces.getCPtr(xmlns));
    return ret;
  }

  public int unsetMetaId() {
    int ret = libsbmlPINVOKE.SBase_unsetMetaId(swigCPtr);
    return ret;
  }

  public int unsetNotes() {
    int ret = libsbmlPINVOKE.SBase_unsetNotes(swigCPtr);
    return ret;
  }

  public int unsetAnnotation() {
    int ret = libsbmlPINVOKE.SBase_unsetAnnotation(swigCPtr);
    return ret;
  }

  public int unsetSBOTerm() {
    int ret = libsbmlPINVOKE.SBase_unsetSBOTerm(swigCPtr);
    return ret;
  }

  public int addCVTerm(CVTerm term, bool newBag) {
    int ret = libsbmlPINVOKE.SBase_addCVTerm__SWIG_0(swigCPtr, CVTerm.getCPtr(term), newBag);
    return ret;
  }

  public int addCVTerm(CVTerm term) {
    int ret = libsbmlPINVOKE.SBase_addCVTerm__SWIG_1(swigCPtr, CVTerm.getCPtr(term));
    return ret;
  }

  public  CVTermList  getCVTerms() { 
  IntPtr cPtr = libsbmlPINVOKE.SBase_getCVTerms__SWIG_0(swigCPtr);
  return (cPtr == IntPtr.Zero) ? null : new CVTermList(cPtr, true);
}

  public long getNumCVTerms() { return (long)libsbmlPINVOKE.SBase_getNumCVTerms(swigCPtr); }

  public CVTerm getCVTerm(long n) {
    IntPtr cPtr = libsbmlPINVOKE.SBase_getCVTerm(swigCPtr, n);
    CVTerm ret = (cPtr == IntPtr.Zero) ? null : new CVTerm(cPtr, false);
    return ret;
  }

  public int unsetCVTerms() {
    int ret = libsbmlPINVOKE.SBase_unsetCVTerms(swigCPtr);
    return ret;
  }

  public int unsetModelHistory() {
    int ret = libsbmlPINVOKE.SBase_unsetModelHistory(swigCPtr);
    return ret;
  }

  public int getResourceBiologicalQualifier(string resource) {
    int ret = libsbmlPINVOKE.SBase_getResourceBiologicalQualifier(swigCPtr, resource);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public int getResourceModelQualifier(string resource) {
    int ret = libsbmlPINVOKE.SBase_getResourceModelQualifier(swigCPtr, resource);
    if (libsbmlPINVOKE.SWIGPendingException.Pending) throw libsbmlPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public Model getModel() {
    IntPtr cPtr = libsbmlPINVOKE.SBase_getModel(swigCPtr);
    Model ret = (cPtr == IntPtr.Zero) ? null : new Model(cPtr, false);
    return ret;
  }

  public long getLevel() { return (long)libsbmlPINVOKE.SBase_getLevel(swigCPtr); }

  public long getVersion() { return (long)libsbmlPINVOKE.SBase_getVersion(swigCPtr); }

  public virtual int getTypeCode() {
    int ret = libsbmlPINVOKE.SBase_getTypeCode(swigCPtr);
    return ret;
  }

  public bool hasValidLevelVersionNamespaceCombination() {
    bool ret = libsbmlPINVOKE.SBase_hasValidLevelVersionNamespaceCombination(swigCPtr);
    return ret;
  }

  public virtual string getElementName() {
    string ret = libsbmlPINVOKE.SBase_getElementName(swigCPtr);
    return ret;
  }

  public string toSBML() {
    string ret = libsbmlPINVOKE.SBase_toSBML(swigCPtr);
    return ret;
  }

  public virtual bool hasRequiredAttributes() {
    bool ret = libsbmlPINVOKE.SBase_hasRequiredAttributes(swigCPtr);
    return ret;
  }

  public virtual bool hasRequiredElements() {
    bool ret = libsbmlPINVOKE.SBase_hasRequiredElements(swigCPtr);
    return ret;
  }

}

}
