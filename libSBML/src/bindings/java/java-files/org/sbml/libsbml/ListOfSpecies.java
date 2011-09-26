/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

/** 
 * LibSBML implementation of SBML Level&nbsp;2's ListOfSpecies construct.
 * <p>
 * The various ListOf___ classes in SBML are merely containers used for
 * organizing the main components of an SBML model.  All are derived from
 * the abstract class {@link SBase}, and inherit the various attributes and
 * subelements of {@link SBase}, such as 'metaid' as and 'annotation'.  The
 * ListOf___ classes do not add any attributes of their own.
 * <p>
 * The relationship between the lists and the rest of an SBML model is
 * illustrated by the following (for SBML Level&nbsp;2 Version&nbsp;4):
 * <p>
 * <center><img src='listof-illustration.jpg'></center><br>
 * 
 * <p>
 * Readers may wonder about the motivations for using the ListOf___
 * containers.  A simpler approach in XML might be to place the components
 * all directly at the top level of the model definition.  We chose instead
 * to group them within XML elements named after {@link ListOf}<em>Classname</em>,
 * in part because we believe this helps organize the components and makes
 * visual reading of models in XML easier.  More importantly, the fact that
 * the container classes are derived from {@link SBase} means that software tools
 * can add information about the lists themselves into each list
 * container's 'annotation'.
 * <p>
 * @see ListOfFunctionDefinitions
 * @see ListOfUnitDefinitions
 * @see ListOfCompartmentTypes
 * @see ListOfSpeciesTypes
 * @see ListOfCompartments
 * @see ListOfSpecies
 * @see ListOfParameters
 * @see ListOfInitialAssignments
 * @see ListOfRules
 * @see ListOfConstraints
 * @see ListOfReactions
 * @see ListOfEvents
 */

public class ListOfSpecies extends ListOf {
   private long swigCPtr;

   protected ListOfSpecies(long cPtr, boolean cMemoryOwn)
   {
     super(libsbmlJNI.SWIGListOfSpeciesUpcast(cPtr), cMemoryOwn);
     swigCPtr = cPtr;
   }

   protected static long getCPtr(ListOfSpecies obj)
   {
     return (obj == null) ? 0 : obj.swigCPtr;
   }

   protected static long getCPtrAndDisown (ListOfSpecies obj)
   {
     long ptr = 0;

     if (obj != null)
     {
       ptr             = obj.swigCPtr;
       obj.swigCMemOwn = false;
     }

     return ptr;
   }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        libsbmlJNI.delete_ListOfSpecies(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  
  /**
   * Creates and returns a deep copy of this ListOfSpeciess instance.
   * <p>
   * @return a (deep) copy of this ListOfSpeciess.
   */
 public ListOfSpecies cloneObject() {
    long cPtr = libsbmlJNI.ListOfSpecies_cloneObject(swigCPtr, this);
    return (cPtr == 0) ? null : new ListOfSpecies(cPtr, true);
  }

  
  /**
   * Returns the libSBML type code for this SBML object.
   * <p>
   * LibSBML attaches an
   * identifying code to every kind of SBML object.  These are known as
   * <em>SBML type codes</em>.  In other languages, the set of type codes
   * is stored in an enumeration; in the Java language interface for
   * libSBML, the type codes are defined as static integer constants in
   * interface class {@link libsbmlConstants}.  The names of the type codes
   * all begin with the characters <code>SBML_</code>. 
   * <p>
   * @return the SBML type code for this object, or @link SBMLTypeCode_t#SBML_UNKNOWN SBML_UNKNOWN@endlink (default).
   * <p>
   * @see #getElementName()
   */
 public int getTypeCode() {
    return libsbmlJNI.ListOfSpecies_getTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the libSBML type code for the objects contained in this {@link ListOf}
   * (i.e., {@link Species} objects, if the list is non-empty).
   * <p>
   * LibSBML attaches an
   * identifying code to every kind of SBML object.  These are known as
   * <em>SBML type codes</em>.  In other languages, the set of type codes
   * is stored in an enumeration; in the Java language interface for
   * libSBML, the type codes are defined as static integer constants in
   * interface class {@link libsbmlConstants}.  The names of the type codes
   * all begin with the characters <code>SBML_</code>. 
   * <p>
   * @return the SBML type code for the objects contained in this {@link ListOf}
   * instance, or @link SBMLTypeCode_t#SBML_UNKNOWN SBML_UNKNOWN@endlink (default).
   * <p>
   * @see #getElementName()
   */
 public int getItemTypeCode() {
    return libsbmlJNI.ListOfSpecies_getItemTypeCode(swigCPtr, this);
  }

  
  /**
   * Returns the XML element name of this object.
   * <p>
   * For ListOfSpeciess, the XML element name is <code>'listOfSpeciess'</code>.
   * <p>
   * @return the name of this element, i.e., <code>'listOfSpeciess'</code>.
   */
 public String getElementName() {
    return libsbmlJNI.ListOfSpecies_getElementName(swigCPtr, this);
  }

  
  /**
   * Get a {@link Species} from the {@link ListOfSpecies}.
   * <p>
   * @param n the index number of the {@link Species} to get.
   * <p>
   * @return the nth {@link Species} in this {@link ListOfSpecies}.
   * <p>
   * @see #size()
   */
 public Species get(long n) {
    long cPtr = libsbmlJNI.ListOfSpecies_get__SWIG_0(swigCPtr, this, n);
    return (cPtr == 0) ? null : new Species(cPtr, false);
  }

  
  /**
   * Get a {@link Species} from the {@link ListOfSpecies}
   * based on its identifier.
   * <p>
   * @param sid a string representing the identifier 
   * of the {@link Species} to get.
   * <p>
   * @return {@link Species} in this {@link ListOfSpecies}
   * with the given id or <code>NULL</code> if no such
   * {@link Species} exists.
   * <p>
   * @see #get(long n)
   * @see #size()
   */
 public Species get(String sid) {
    long cPtr = libsbmlJNI.ListOfSpecies_get__SWIG_2(swigCPtr, this, sid);
    return (cPtr == 0) ? null : new Species(cPtr, false);
  }

  
  /**
   * Removes the nth item from this ListOfSpeciess items and returns a pointer to
   * it.
   * <p>
   * The caller owns the returned item and is responsible for deleting it.
   * <p>
   * @param n the index of the item to remove
   * <p>
   * @see #size()
   */
 public Species remove(long n) {
    long cPtr = libsbmlJNI.ListOfSpecies_remove__SWIG_0(swigCPtr, this, n);
    return (cPtr == 0) ? null : new Species(cPtr, true);
  }

  
  /**
   * Removes item in this ListOfSpeciess items with the given identifier.
   * <p>
   * The caller owns the returned item and is responsible for deleting it.
   * If none of the items in this list have the identifier <code>sid</code>, then 
   * <code>NULL</code> is returned.
   * <p>
   * @param sid the identifier of the item to remove
   * <p>
   * @return the item removed.  As mentioned above, the caller owns the
   * returned item.
   */
 public Species remove(String sid) {
    long cPtr = libsbmlJNI.ListOfSpecies_remove__SWIG_1(swigCPtr, this, sid);
    return (cPtr == 0) ? null : new Species(cPtr, true);
  }

  public ListOfSpecies() {
    this(libsbmlJNI.new_ListOfSpecies(), true);
  }

}
