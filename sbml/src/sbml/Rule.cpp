/**
 * @file    Rule.cpp
 * @brief   Implementations of Rule, ListOfRules, AlgebraicRule, AssignmentRule
 *          and RateRule.
 * @author  Ben Bornstein
 *
 * $Id: Rule.cpp 10156 2009-09-01 12:09:35Z sarahkeating $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/sbml/Rule.cpp $
 *
 * <!--------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright (C) 2009-2011 jointly by the following organizations: 
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. EMBL European Bioinformatics Institute (EBML-EBI), Hinxton, UK
 *  
 * Copyright (C) 2006-2008 by the California Institute of Technology,
 *     Pasadena, CA, USA 
 *  
 * Copyright (C) 2002-2005 jointly by the following organizations: 
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. Japan Science and Technology Agency, Japan
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 * ---------------------------------------------------------------------- -->*/

#include <sbml/xml/XMLNode.h>
#include <sbml/xml/XMLAttributes.h>
#include <sbml/xml/XMLInputStream.h>
#include <sbml/xml/XMLOutputStream.h>
#include <sbml/xml/XMLNamespaces.h>

#include <sbml/math/FormulaFormatter.h>
#include <sbml/math/FormulaParser.h>
#include <sbml/math/MathML.h>
#include <sbml/math/ASTNode.h>

#include <sbml/SBO.h>
#include <sbml/SBMLTypeCodes.h>
#include <sbml/SBMLVisitor.h>
#include <sbml/SBMLError.h>
#include <sbml/SBMLDocument.h>
#include <sbml/Model.h>
#include <sbml/Rule.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */
/*
 * Only subclasses may create Rules.
 */
Rule::Rule (int type, unsigned int level, unsigned int version)
 :
   SBase   ( level, version)
 , mFormula( ""  )
 , mMath   (  NULL       )
 , mType   ( type     )
 , mL1Type ( SBML_UNKNOWN )
 , mInternalId ( "" )
{
}

Rule::Rule (int type, SBMLNamespaces * sbmlns) :
   SBase   ( sbmlns )
 , mFormula( ""       )
 , mMath   (  NULL       )
 , mType   ( type     )
 , mL1Type ( SBML_UNKNOWN )
 , mInternalId ( "" )
{
}
/** @endcond */


/*
 * Destroys this Rule.
 */
Rule::~Rule ()
{
  delete mMath;
}


/*
 * Copy constructor. Creates a copy of this Rule.
 */
Rule::Rule (const Rule& orig) :
   SBase   ( orig          )
 , mMath   ( NULL            )
{
  if (&orig == NULL)
  {
    throw SBMLConstructorException("Null argument to copy constructor");
  }
  else
  {
    mVariable   = orig.mVariable;
    mFormula    = orig.mFormula ;
    mUnits      = orig.mUnits   ;
    mType       = orig.mType    ;
    mL1Type     = orig.mL1Type  ;
    mInternalId = orig.mInternalId;

    if (orig.mMath != NULL) 
    {
      mMath = orig.mMath->deepCopy();
      mMath->setParentSBMLObject(this);
    }
  }
}


/*
 * Assignment operator.
 */
Rule& Rule::operator=(const Rule& rhs)
{
  if (&rhs == NULL)
  {
    throw SBMLConstructorException("Null argument to assignment operator");
  }
  else if(&rhs!=this)
  {
    this->SBase::operator =(rhs);
    mVariable = rhs.mVariable;
    mFormula = rhs.mFormula ;
    mUnits   = rhs.mUnits   ;
    mType    = rhs.mType    ;
    mL1Type  = rhs.mL1Type  ;
    mInternalId = rhs.mInternalId;

    delete mMath;
    if (rhs.mMath != NULL) 
    {
      mMath = rhs.mMath->deepCopy();
      mMath->setParentSBMLObject(this);
    }
    else
    {
      mMath = NULL;
    }
   
  }

  return *this;
}


/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Model's next Rule
 * (if available).
 */
bool
Rule::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * @return a (deep) copy of this Rule.
 */
Rule*
Rule::clone () const
{
  return new Rule(*this);
}


/*
 * @return the formula for this Rule.
 */
const string&
Rule::getFormula () const
{
  if (mFormula.empty() == true && mMath != NULL)
  {
    char* s  = SBML_formulaToString(mMath);
    mFormula = s;

    free(s);
  }

  return mFormula;
}


/*
 * @return the math for this Rule.
 */
const ASTNode*
Rule::getMath () const
{
  if (mMath == NULL && mFormula.empty() == false)
  {
    mMath = SBML_parseFormula( mFormula.c_str() );
  }

  return mMath;
}


/*
 * @return the variable for this Rule.
 */
const string&
Rule::getVariable () const
{
  return mVariable;
}


/**
 * @return the string of variable attribute of this object.
 * 
 * @note this function is an alias for getVariable()
 * 
 * @see getVariable()
 */
const std::string& 
Rule::getId() const
{
  return getVariable();
}


/*
 * @return the units for this Rule (L1 ParameterRules only).
 */
const string&
Rule::getUnits () const
{
  return mUnits;
}


/*
 * @return true if the formula (or equivalently the math) for this Rule is
 * set, false otherwise.
 */
bool
Rule::isSetFormula () const
{
  return (mFormula.empty() == false) || (mMath != NULL);
}


/*
 * @return true if the math (or equivalently the formula) for this Rule is
 * set, false otherwise.
 */
bool
Rule::isSetMath () const
{
  /* if the formula has been set but it is not a correct formula
   * it cannot be correctly transferred to an ASTNode so in fact
   * getMath will return @c NULL
   *
   * this function needs to test for this
   */
  bool formula = isSetFormula();
  
  if (formula)
  {
    const ASTNode *temp = getMath();
    if (temp == NULL)
      formula = false;
  }
    
  return formula;
}


/*
 * @return true if the variable of this Rule is set, false
 * otherwise.
 */
bool
Rule::isSetVariable () const
{
  return (mVariable.empty() == false);
}


/*
 * @return true if the units for this Rule is set, false otherwise
 * (L1 ParameterRules only).
 */
bool
Rule::isSetUnits () const
{
  return (mUnits.empty() == false);
}


/*
 * Sets the formula of this Rule to a copy of string.
 */
int
Rule::setFormula (const std::string& formula)
{
  if (&(formula) == NULL)
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else 
  {
    ASTNode * math = SBML_parseFormula(formula.c_str());
    if (formula == "")
    {
      mFormula.erase();
      delete mMath;
      mMath = NULL;
      return LIBSBML_OPERATION_SUCCESS;
    }
    else if (math == NULL || !(math->isWellFormedASTNode()))
    {
      return LIBSBML_INVALID_OBJECT;
    }
    else
    {
      mFormula = formula;

      if (mMath != NULL)
      {
        delete mMath;
        mMath = NULL;
      }
      return LIBSBML_OPERATION_SUCCESS;
    }
  }
}


/*
 * Sets the math of this Rule to a copy of the given ASTNode.
 */
int
Rule::setMath (const ASTNode* math)
{
  if (mMath == math) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else if (math == NULL)
  {
    delete mMath;
    mMath = NULL;
    mFormula.erase();
    return LIBSBML_OPERATION_SUCCESS;
  }
  else if (!(math->isWellFormedASTNode()))
  {
    return LIBSBML_INVALID_OBJECT;
  }
  else
  {
    delete mMath;
    mMath = (math != NULL) ? math->deepCopy() : NULL;
    if (mMath != NULL) mMath->setParentSBMLObject(this);
    mFormula.erase();
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the variable of this Rule to a copy of sid.
 */
int
Rule::setVariable (const std::string& sid)
{
  if (&(sid) == NULL)
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else if (isAlgebraic())
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else if (!(SyntaxChecker::isValidInternalSId(sid)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mVariable = sid;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Sets the units for this Rule to a copy of sname (L1 ParameterRules
 * only).
 */
int
Rule::setUnits (const std::string& sname)
{
  /* only in L1 ParameterRule */
  if (&(sname) == NULL)
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else if (getLevel() > 1)
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else if ( !isParameter())
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else if (!(SyntaxChecker::isValidInternalUnitSId(sname)))
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  else
  {
    mUnits = sname;
    return LIBSBML_OPERATION_SUCCESS;
  }
}


/*
 * Unsets the units for this Rule (L1 ParameterRules only).
 */
int
Rule::unsetUnits ()
{
  /* only in L1 Parameter rule */
  if (getLevel() > 1)
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }
  else if ( !isParameter())
  {
    return LIBSBML_UNEXPECTED_ATTRIBUTE;
  }

  mUnits.erase();

  if (mUnits.empty()) 
  {
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


/*
  * Calculates and returns a UnitDefinition that expresses the units
  * returned by the math expression of this Rule.
  */
UnitDefinition * 
Rule::getDerivedUnitDefinition()
{
  if (!isSetMath())
    return NULL;
  /* if we have the whole model but it is not in a document
   * it is still possible to determine the units
   */
  Model * m = static_cast <Model *> (getAncestorOfType(SBML_MODEL));

  if (m != NULL)
  {
    if (!m->isPopulatedListFormulaUnitsData())
    {
      m->populateListFormulaUnitsData();
    }
    
    if (isAlgebraic())
    {
      if (m->getFormulaUnitsData(getInternalId(), getTypeCode()) != NULL)
      {
        return m->getFormulaUnitsData(getInternalId(), getTypeCode())
                                               ->getUnitDefinition();
      }
      else
      {
        return NULL;
      }
    }
    else
    {
      if (m->getFormulaUnitsData(getVariable(), getTypeCode()) != NULL)
      {
        return m->getFormulaUnitsData(getVariable(), getTypeCode())
                                              ->getUnitDefinition();
      }
      else
      {
        return NULL;
      }
    }
  }
  else
  {
    return NULL;
  }

}


/*
  * Constructs and returns a UnitDefinition that expresses the units of this 
  * Compartment.
  */
const UnitDefinition *
Rule::getDerivedUnitDefinition() const
{
  return const_cast <Rule *> (this)->getDerivedUnitDefinition();
}


/*
 * Predicate returning @c true if 
 * the math expression of this Rule contains
 * parameters/numbers with undeclared units that cannot be ignored.
 */
bool 
Rule::containsUndeclaredUnits()
{
  if (!isSetMath())
    return false;
  /* if we have the whole model but it is not in a document
   * it is still possible to determine the units
   */
  Model * m = static_cast <Model *> (getAncestorOfType(SBML_MODEL));

  if (m != NULL)
  {
    if (!m->isPopulatedListFormulaUnitsData())
    {
      m->populateListFormulaUnitsData();
    }
    
    if (isAlgebraic())
    {
      if (m->getFormulaUnitsData(getInternalId(), getTypeCode()) != NULL)
      {
        return m->getFormulaUnitsData(getInternalId(), getTypeCode())
          ->getContainsUndeclaredUnits();
      }
      else
      {
        return false;
      }
    }
    else
    {
      if (m->getFormulaUnitsData(getVariable(), getTypeCode()) != NULL)
      {
        return m->getFormulaUnitsData(getVariable(), getTypeCode())
          ->getContainsUndeclaredUnits();
      }
      else
      {
        return false;
      }
    }
  }
  else
  {
    return false;
  }
}



bool 
Rule::containsUndeclaredUnits() const
{
  return const_cast<Rule *> (this)->containsUndeclaredUnits();
}


/*
 * @return the type of this Rule, either RULE_TYPE_RATE or
 * RULE_TYPE_SCALAR.
 */
RuleType_t
Rule::getType () const
{
  if (mType == SBML_ASSIGNMENT_RULE) return RULE_TYPE_SCALAR;
  if (mType == SBML_RATE_RULE)       return RULE_TYPE_RATE;
  return RULE_TYPE_INVALID;
}


/*
 * @return true if this Rule is an AlgebraicRule, false otherwise.
 */
bool
Rule::isAlgebraic () const
{
  return (mType == SBML_ALGEBRAIC_RULE);
}


/*
 * @return true if this Rule is an AssignmentRule, false otherwise.
 */
bool
Rule::isAssignment () const
{
  return (mType == SBML_ASSIGNMENT_RULE);
}


/*
 * @return true if this Rule is a CompartmentVolumeRule, false otherwise.
 */
bool
Rule::isCompartmentVolume () const
{
  if (mL1Type == SBML_COMPARTMENT_VOLUME_RULE)
  {
    return true;
  }
  else
  {
    const Model* model = getModel();
    return (model == NULL) ? false : model->getCompartment( getVariable() ) != NULL;
  }
}


/*
 * @return true if this Rule is a ParameterRule, false otherwise.
 */
bool
Rule::isParameter () const
{
  if (mL1Type == SBML_PARAMETER_RULE)
  {
    return true;
  }
  else
  {
    const Model* model = getModel();
    return (model == NULL) ? false : model->getParameter( getVariable() ) != NULL;
  }
}


/*
 * @return true if this Rule is a RateRule (L2) or has type="rate" (L1),
 * false otherwise.
 */
bool
Rule::isRate () const
{
  return (mType == SBML_RATE_RULE);
}


/*
 * @return true if this Rule is an AssignmentRule (L2) has type="scalar"
 * (L1), false otherwise.
 */
bool
Rule::isScalar () const
{
  return (mType == SBML_ASSIGNMENT_RULE);
}


/*
 * @return true if this Rule is a SpeciesConcentrationRule, false
 * otherwise.
 */
bool
Rule::isSpeciesConcentration () const
{
  if (mL1Type == SBML_SPECIES_CONCENTRATION_RULE)
  {
    return true;
  }
  else
  {
    const Model* model = getModel();
    return (model == NULL) ? false : model->getSpecies( getVariable() ) != NULL;
  }
}


/*
 * @return the typecode (int) of this SBML object or SBML_UNKNOWN
 * (default).
 *
 * @see getElementName()
 */
int
Rule::getTypeCode () const
{
  return mType;
}


/*
 * @return the SBML Level 1 typecode for this Rule or SBML_UNKNOWN
 * (default).
 */
int
Rule::getL1TypeCode () const
{
  return mL1Type;
}


/*
 * @return the name of this element eg "algebraicRule".
 
 */
const string&
Rule::getElementName () const
{
  static const string algebraic   = "algebraicRule";
  static const string specie      = "specieConcentrationRule";
  static const string species     = "speciesConcentrationRule";
  static const string compartment = "compartmentVolumeRule";
  static const string parameter   = "parameterRule";
  static const string assignment  = "assignmentRule";
  static const string rate        = "rateRule";
  static const string unknown     = "unknownRule";

  if ( isAlgebraic() )
  {
    return algebraic;
  }
  else if (getLevel() == 1)
  {
    if ( isSpeciesConcentration() )
    {
      return (getVersion() == 2) ? species : specie;
    }
    else if ( isCompartmentVolume() )
    {
      return compartment;
    }
    else if ( isParameter() )
    {
      return parameter;
    }
  }
  else
  {
    if ( isAssignment() )
    {
      return assignment;
    }
    else if ( isRate() )
    {
      return rate;
    }
  }

  return unknown;
}


bool 
Rule::hasRequiredElements() const
{
  bool allPresent = true;

  /* required attributes for rule: math */

  if (!isSetMath())
    allPresent = false;

  return allPresent;
}


/** @cond doxygen-libsbml-internal */
bool 
Rule::hasRequiredAttributes() const
{
  bool allPresent = true;

  /* required attributes for rules:(formula in L1) */

  if (getLevel() == 1 && !isSetFormula())
    allPresent = false;

  return allPresent;
}


/*
 * Subclasses should override this method to write out their contained
 * SBML objects as XML elements.  Be sure to call your parents
 * implementation of this method as well.
 */
void
Rule::writeElements (XMLOutputStream& stream) const
{
  SBase::writeElements(stream);
  if ( getLevel() > 1 && isSetMath() ) writeMathML(getMath(), stream);

  //
  // (EXTENSION)
  //
  SBase::writeExtensionElements(stream);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read (and store) XHTML,
 * MathML, etc. directly from the XMLInputStream.
 *
 * @return true if the subclass read from the stream, false otherwise.
 */
bool
Rule::readOtherXML (XMLInputStream& stream)
{
  bool          read = false;
  const string& name = stream.peek().getName();

  if (name == "math")
  {
    // if this is level 1 there shouldnt be any math!!!
    if (getLevel() == 1) 
    {
      logError(NotSchemaConformant, getLevel(), getVersion(),
	       "SBML Level 1 does not support MathML.");
      delete mMath;
      return false;
    }

    if (mMath != NULL)
    {
      if (getLevel() < 3) 
      {
        logError(NotSchemaConformant, getLevel(), getVersion(),
	        "Only one <math> element is permitted inside a "
	        "particular containing element.");
      }
      else
      {
        logError(OneMathElementPerRule, getLevel(), getVersion());
      }
    }
    delete mMath;

    /* check for MathML namespace 
     * this may be explicitly declared here
     * or implicitly declared on the whole document
     */
    const XMLToken elem = stream.peek();
    const std::string prefix = checkMathMLNamespace(elem);

    mMath = readMathML(stream, prefix);
    if (mMath != NULL) mMath->setParentSBMLObject(this);
    read  = true;
  }

  /* ------------------------------
   *
   *   (EXTENSION)
   *
   * ------------------------------ */
  if ( SBase::readOtherXML(stream) )
    read = true;

  return read;
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/**
 * Subclasses should override this method to get the list of
 * expected attributes.
 * This function is invoked from corresponding readAttributes()
 * function.
 */
void
Rule::addExpectedAttributes(ExpectedAttributes& attributes)
{
  SBase::addExpectedAttributes(attributes);

  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();

  if (level == 1)
  {
    attributes.add("formula");
    const string s = (level == 1 && version == 1) ? "specie" : "species";
    attributes.add(s);
    attributes.add("compartment");
    attributes.add("name");
    attributes.add("units");
    attributes.add("type");
  }
  else
  {
    attributes.add("variable");
    if (level == 2 && version == 2)
    {
      attributes.add("sboTerm");
    }
  }
}


/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Rule::readAttributes (const XMLAttributes& attributes,
                      const ExpectedAttributes& expectedAttributes)
{
  const unsigned int level   = getLevel  ();

  SBase::readAttributes(attributes, expectedAttributes);

  switch (level)
  {
  case 1:
    readL1Attributes(attributes);
    break;
  case 2:
    readL2Attributes(attributes);
    break;
  case 3:
  default:
    readL3Attributes(attributes);
    break;
  }

}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Rule::readL1Attributes (const XMLAttributes& attributes)
{
  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();

  //
  // formula: string  { use="required" }  (L1v1, L1v2)
  //
  attributes.readInto("formula", mFormula, getErrorLog(), true);

  //
  // type { use="optional" default="scalar" }  (L1v1, L1v2)
  //
  // This attribute is handled by ListOfRules::createObject();
  //

  if ( isSpeciesConcentration() )
  {
    //
    // specie : SName   { use="required" }  (L1v1)
    // species: SName   { use="required" }  (L1v2)
    //
    const string s = (level == 1 && version == 1) ? "specie" : "species";
    bool assigned = attributes.readInto(s, mVariable, getErrorLog(), true);
    if (assigned && mVariable.size() == 0)
    {
      logEmptyString(s, level, version, "<rule>");
    }
    if (!SyntaxChecker::isValidInternalSId(mVariable)) 
      logError(InvalidIdSyntax);
  }
  else if ( isCompartmentVolume() )
  {
    //
    // compartment: SName  { use="required" }  (L1v1, L1v2)
    //
    bool assigned = attributes.readInto("compartment", mVariable, getErrorLog(), true);
    if (assigned && mVariable.size() == 0)
    {
      logEmptyString("compartment", level, version, "<rule>");
    }
    if (!SyntaxChecker::isValidInternalSId(mVariable)) 
      logError(InvalidIdSyntax);
  }
  else if ( isParameter() )
  {
    //
    // name: SName  { use="required" } (L1v1, L1v2)
    //
    bool assigned = attributes.readInto("name", mVariable, getErrorLog(), true);
    if (assigned && mVariable.size() == 0)
    {
      logEmptyString("name", level, version, "<rule>");
    }
    if (!SyntaxChecker::isValidInternalSId(mVariable)) 
      logError(InvalidIdSyntax);

    //
    // units  { use="optional" }  (L1v1, L1v2);
    //
    attributes.readInto("units", mUnits);

  }
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Rule::readL2Attributes (const XMLAttributes& attributes)
{
  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();

  if (isAssignment() || isRate())
  {
    //
    // variable: SId  { use="required" }  (L2v1 ->)
    //
    bool assigned = attributes.readInto("variable", mVariable, getErrorLog(), true);
    if (assigned && mVariable.size() == 0)
    {
      logEmptyString("variable", level, version, "<rule>");
    }
    if (!SyntaxChecker::isValidInternalSId(mVariable)) 
      logError(InvalidIdSyntax);
  }

  //
  // sboTerm: SBOTerm { use="optional" }  (L2v2 ->)
  //
  if (version == 2) 
    mSBOTerm = SBO::readTerm(attributes, this->getErrorLog(), level, version);
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to read values from the given
 * XMLAttributes set into their specific fields.  Be sure to call your
 * parents implementation of this method as well.
 */
void
Rule::readL3Attributes (const XMLAttributes& attributes)
{
  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();

  if (isAssignment() || isRate())
  {
    //
    // variable: SId  { use="required" }  (L2v1 ->)
    //
    bool assigned = attributes.readInto("variable", mVariable, getErrorLog());
    if (!assigned)
    {
      if (isAssignment())
        getErrorLog()->logError(AllowedAttributesOnAssignRule, level, 
                                version);
      else
        getErrorLog()->logError(AllowedAttributesOnRateRule, level, 
                                version);

    }
    if (assigned && mVariable.size() == 0)
    {
      logEmptyString("variable", level, version, "<rule>");
    }
    if (!SyntaxChecker::isValidInternalSId(mVariable)) 
      logError(InvalidIdSyntax);
  }
}
/** @endcond */


/** @cond doxygen-libsbml-internal */
/*
 * Subclasses should override this method to write their XML attributes
 * to the XMLOutputStream.  Be sure to call your parents implementation
 * of this method as well.
 */
void
Rule::writeAttributes (XMLOutputStream& stream) const
{
  SBase::writeAttributes(stream);

  const unsigned int level   = getLevel  ();
  const unsigned int version = getVersion();


  if (level == 1)
  {
    //
    // formula: string  { use="required" }  (L1v1, L1v2)
    //
    stream.writeAttribute("formula", getFormula());

    //
    // type { use="optional" default="scalar" }  (L1v1, L1v2)
    //
    if (getType() == RULE_TYPE_RATE)
    {
      const string rate = "rate";
      stream.writeAttribute("type", rate);
    }

    //
    // specie : SName   { use="required" }  (L1v1)
    // species: SName   { use="required" }  (L1v2)
    //
    if ( isSpeciesConcentration() )
    {
      const string species = (version == 1) ? "specie" : "species";
      stream.writeAttribute(species, mVariable);
    }

    //
    // compartment: SName  { use="required" }  (L1v1, L1v2)
    //
    else if ( isCompartmentVolume() )
    {
      stream.writeAttribute("compartment", mVariable);
    }

    else if ( isParameter() )
    {
      //
      // name: SName  { use="required" } (L1v1, L1v2)
      //
      stream.writeAttribute("name", mVariable);

      //
      // units  { use="optional" }  (L1v1, L1v2);
      //
      stream.writeAttribute("units", mUnits);
    }
  }
  else if (level > 1)
  {
    //
    // sboTerm: SBOTerm { use="optional" }  (L2v2->)
    //
    // sboTerm for L2V3 or later is written in SBase::writeAttributes()
    //
    if ( (level == 2) && (version == 2) )
    {
      SBO::writeTerm(stream, mSBOTerm);
    }

    //
    // variable: SId  { use="required" }  (L2v1-> )
    //
    if(!isAlgebraic())
      stream.writeAttribute("variable", mVariable);
  }

  //
  // (EXTENSION)
  //
  SBase::writeExtensionAttributes(stream);
}
/** @endcond */


/*
 * Sets the SBML Level 1 typecode for this Rule.
 */
int
Rule::setL1TypeCode (int type)
{
  if (    (type == SBML_PARAMETER_RULE) 
       || (type == SBML_COMPARTMENT_VOLUME_RULE) 
       || (type == SBML_SPECIES_CONCENTRATION_RULE) 
     )
  {
    mSBMLNamespaces->setLevel(1);
    mL1Type = type;
  }
  else
  {
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
  return LIBSBML_OPERATION_SUCCESS;
}



AlgebraicRule::AlgebraicRule (unsigned int level, unsigned int version) :
  Rule(SBML_ALGEBRAIC_RULE, level, version)
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();

  mInternalIdOnly = false;
}


AlgebraicRule::AlgebraicRule (SBMLNamespaces * sbmlns) :
  Rule(SBML_ALGEBRAIC_RULE, sbmlns)
{
  if (!hasValidLevelVersionNamespaceCombination())
  {
    std::string err(getElementName());
    XMLNamespaces* xmlns = sbmlns->getNamespaces();
    if (xmlns)
    {
      std::ostringstream oss;
      XMLOutputStream xos(oss);
      xos << *xmlns;
      err.append(oss.str());
    }
    throw SBMLConstructorException(err);
  }

  mInternalIdOnly = false;

  loadPlugins(sbmlns);
}


/*
 * Destroys this AlgebraicRule.
 */
AlgebraicRule::~AlgebraicRule ()
{
}

/*
 * @return a (deep) copy of this Rule.
 */
AlgebraicRule*
AlgebraicRule::clone () const
{
  return new AlgebraicRule(*this);
}


/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Model's next Rule
 * (if available).
 */
bool
AlgebraicRule::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}

bool 
AlgebraicRule::hasRequiredAttributes() const
{
  bool allPresent = Rule::hasRequiredAttributes();

  return allPresent;
}



/** @cond doxygen-libsbml-internal */

/*
 * sets the mInternalIdOnly flag
 */
void 
AlgebraicRule::setInternalIdOnly()
{
  mInternalIdOnly = true;
}

/*
 * gets the mInternalIdOnly flag
 */
bool 
AlgebraicRule::getInternalIdOnly() const
{
  return mInternalIdOnly;
}

  /** @endcond */



AssignmentRule::AssignmentRule (unsigned int level, unsigned int version) :
  Rule(SBML_ASSIGNMENT_RULE, level, version)
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}

AssignmentRule::AssignmentRule (SBMLNamespaces *sbmlns) :
  Rule(SBML_ASSIGNMENT_RULE, sbmlns)
{
  if (!hasValidLevelVersionNamespaceCombination())
  {
    std::string err(getElementName());
    XMLNamespaces* xmlns = sbmlns->getNamespaces();
    if (xmlns)
    {
      std::ostringstream oss;
      XMLOutputStream xos(oss);
      xos << *xmlns;
      err.append(oss.str());
    }
    throw SBMLConstructorException(err);
  }

  loadPlugins(sbmlns);
}


/*
 * Destroys this AssignmentRule.
 */
AssignmentRule::~AssignmentRule ()
{
}


/*
 * @return a (deep) copy of this Rule.
 */
AssignmentRule*
AssignmentRule::clone () const
{
  return new AssignmentRule(*this);
}


/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Model's next Rule
 * (if available).
 */
bool
AssignmentRule::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


bool 
AssignmentRule::hasRequiredAttributes() const
{
  bool allPresent = Rule::hasRequiredAttributes();

  /* required attributes for assignment: variable (comp/species/name in L1) */

  if (!isSetVariable())
    allPresent = false;

  return allPresent;
}


RateRule::RateRule (unsigned int level, unsigned int version) :
  Rule(SBML_RATE_RULE, level, version)
{
  if (!hasValidLevelVersionNamespaceCombination())
    throw SBMLConstructorException();
}

RateRule::RateRule (SBMLNamespaces *sbmlns) :
  Rule(SBML_RATE_RULE, sbmlns)
{
  if (!hasValidLevelVersionNamespaceCombination())
  {
    std::string err(getElementName());
    XMLNamespaces* xmlns = sbmlns->getNamespaces();
    if (xmlns)
    {
      std::ostringstream oss;
      XMLOutputStream xos(oss);
      xos << *xmlns;
      err.append(oss.str());
    }
    throw SBMLConstructorException(err);
  }

  loadPlugins(sbmlns);
}


/*
 * Destroys this RateRule.
 */
RateRule::~RateRule ()
{
}

/*
 * @return a (deep) copy of this Rule.
 */
RateRule*
RateRule::clone () const
{
  return new RateRule(*this);
}


bool 
RateRule::hasRequiredAttributes() const
{
  bool allPresent = Rule::hasRequiredAttributes();

  /* required attributes for rateRule: variable (comp/species/name in L1) */

  if (!isSetVariable())
    allPresent = false;

  return allPresent;
}


/*
 * Accepts the given SBMLVisitor.
 *
 * @return the result of calling <code>v.visit()</code>, which indicates
 * whether or not the Visitor would like to visit the Model's next Rule
 * (if available).
 */
bool
RateRule::accept (SBMLVisitor& v) const
{
  return v.visit(*this);
}


/*
 * Creates a new ListOfRules items.
 */
ListOfRules::ListOfRules (unsigned int level, unsigned int version)
  : ListOf(level,version)
{
}


/*
 * Creates a new ListOfRules items.
 */
ListOfRules::ListOfRules (SBMLNamespaces* sbmlns)
  : ListOf(sbmlns)
{
}


/*
 * @return a (deep) copy of this ListOfRules.
 */
ListOfRules*
ListOfRules::clone () const
{
  return new ListOfRules(*this);
}


/*
 * @return the typecode (int) of this SBML object or SBML_UNKNOWN
 * SBML_UNKNOWN (default).
 */
int
ListOfRules::getItemTypeCode () const
{
  return SBML_RULE;
}


/*
 * @return the name of this element ie "listOfRules".
 */
const string&
ListOfRules::getElementName () const
{
  static const string name = "listOfRules";
  return name;
}


/* return nth item in list */
Rule *
ListOfRules::get(unsigned int n)
{
  return static_cast<Rule*>(ListOf::get(n));
}


/* return nth item in list */
const Rule *
ListOfRules::get(unsigned int n) const
{
  return static_cast<const Rule*>(ListOf::get(n));
}

/**
 * Used by ListOf::get() to lookup an SBase based by its id.
 */
struct IdEqRule : public unary_function<SBase*, bool>
{
  const string& id;

  IdEqRule (const string& id) : id(id) { }
  bool operator() (SBase* sb) 
       { return static_cast <Rule *> (sb)->getVariable() == id; }
};


/* return item by id */
Rule*
ListOfRules::get (const std::string& sid)
{
  return const_cast<Rule*>( 
    static_cast<const ListOfRules&>(*this).get(sid) );
}


/* return item by id */
const Rule*
ListOfRules::get (const std::string& sid) const
{
  vector<SBase*>::const_iterator result;

  if (&(sid) == NULL)
  {
    return NULL;
  }
  else
  {
    result = find_if( mItems.begin(), mItems.end(), IdEqRule(sid) );
    return (result == mItems.end()) ? NULL : static_cast <Rule*> (*result);
  }
}


/* Removes the nth item from this list */
Rule*
ListOfRules::remove (unsigned int n)
{
   return static_cast<Rule*>(ListOf::remove(n));
}


/* Removes item in this list by id */
Rule*
ListOfRules::remove (const std::string& sid)
{
  SBase* item = NULL;
  vector<SBase*>::iterator result;

  if (&(sid) != NULL)
  {
    result = find_if( mItems.begin(), mItems.end(), IdEqRule(sid) );

    if (result != mItems.end())
    {
      item = *result;
      mItems.erase(result);
    }
  }

  return static_cast <Rule*> (item);
}


/** @cond doxygen-libsbml-internal */
/*
 * @return the ordinal position of the element with respect to its siblings
 * or -1 (default) to indicate the position is not significant.
 */
int
ListOfRules::getElementPosition () const
{
  return 9;
}
/** @endcond */

/** @cond doxygen-libsbml-internal */

bool 
ListOfRules::isValidTypeForList(SBase * item)
{
  int tc = item->getTypeCode();
  return ((tc == SBML_ALGEBRAIC_RULE )
    ||    (tc == SBML_ASSIGNMENT_RULE )
    ||    (tc == SBML_RATE_RULE )
    ||    (tc == SBML_SPECIES_CONCENTRATION_RULE )
    ||    (tc == SBML_COMPARTMENT_VOLUME_RULE )
    ||    (tc == SBML_PARAMETER_RULE));
}

/** @endcond */

/** @cond doxygen-libsbml-internal */
/*
 * @return the SBML object corresponding to next XMLToken in the
 * XMLInputStream or @c NULL if the token was not recognized.
 */
SBase*
ListOfRules::createObject (XMLInputStream& stream)
{
  const unsigned int level  = getLevel();
  const string&      name   = stream.peek().getName();
  Rule*              object = NULL;


  if (name == "algebraicRule")
  {
    try
    {
      object = new AlgebraicRule(getSBMLNamespaces());
    }
    catch (SBMLConstructorException*)
    {
      object = new AlgebraicRule(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
    catch ( ... )
    {
      object = new AlgebraicRule(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
    }
  }
  else if (level == 1)
  {
    string type = "scalar";
    stream.peek().getAttributes().readInto("type", type);

    if (type == "scalar")
    {
      try
      {
        object = new AssignmentRule(getSBMLNamespaces());
      }
      catch (SBMLConstructorException*)
      {
        object = new AssignmentRule(1, 2);
      }
      catch ( ... )
      {
        object = new AssignmentRule(1, 2);
      }
    }
    else if (type == "rate")
    {
      try
      {
        object = new RateRule(getSBMLNamespaces());
      }
      catch (SBMLConstructorException*)
      {
        object = new RateRule(1, 2);
      }
      catch ( ... )
      {
        object = new RateRule(1, 2);
      }
    }

    if (object)
    {
      if ( name == "speciesConcentrationRule" ||
           name == "specieConcentrationRule" )
      {
        object->setL1TypeCode(SBML_SPECIES_CONCENTRATION_RULE);
      }
      else if (name == "compartmentVolumeRule")
      {
        object->setL1TypeCode(SBML_COMPARTMENT_VOLUME_RULE);
      }
      else if (name == "parameterRule")
      {
        object->setL1TypeCode(SBML_PARAMETER_RULE);
      }
      else
      {
        delete object;
        object = NULL;
      }
    }
  }
  else
  {
    if (name == "assignmentRule")
    {
      try
      {
        object = new AssignmentRule(getSBMLNamespaces());
      }
      catch (SBMLConstructorException*)
      {
        object = new AssignmentRule(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
      }
      catch ( ... )
      {
        object = new AssignmentRule(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
      }
    }
    else if (name == "rateRule")
    {
      try
      {
        object = new RateRule(getSBMLNamespaces());
      }
      catch (SBMLConstructorException*)
      {
        object = new RateRule(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
      }
      catch ( ... )
      {
        object = new RateRule(SBMLDocument::getDefaultLevel(),
        SBMLDocument::getDefaultVersion());
      }
    }
  }

  if (object) mItems.push_back(object);

  return object;
}
/** @endcond */


/** @cond doxygen-c-only */


/**
 * Creates a new AlgebraicRule (Rule_t) structure using the given SBML 
 * @p level and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * AlgebraicRule
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * AlgebraicRule
 *
 * @return a pointer to the newly created Rule_t structure.
 *
 * @note Once a AlgebraicRule has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the AlgebraicRule.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
Rule_t *
Rule_createAlgebraic (unsigned int level, unsigned int version)
{
  try
  {
    AlgebraicRule* obj = new AlgebraicRule(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new AlgebraicRule (Rule_t) structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this AlgebraicRule
 *
 * @return a pointer to the newly created Rule_t structure.
 *
 * @note Once a AlgebraicRule has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the AlgebraicRule.  Despite this, the ability to supply the values at creation 
 * time is an important aid to creating valid SBML.  Knowledge of the intended 
 * SBML Level and Version determine whether it is valid to assign a particular 
 * value to an attribute, or whether it is valid to add an object to an existing
 * SBMLDocument.
 */
LIBSBML_EXTERN
Rule_t *
Rule_createAlgebraicWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    AlgebraicRule* obj = new AlgebraicRule(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new AssignmentRule (Rule_t) structure using the given SBML
 * @p level and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * AssignmentRule
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * AssignmentRule
 *
 * @return a pointer to the newly created Rule_t structure.
 *
 * @note Once a AssignmentRule has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the AssignmentRule.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
Rule_t *
Rule_createAssignment (unsigned int level, unsigned int version)
{
  try
  {
    AssignmentRule* obj = new AssignmentRule(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new AssignmentRule (Rule_t) structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this AssignmentRule
 *
 * @return a pointer to the newly created Rule_t structure.
 *
 * @note Once a AssignmentRule has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the AssignmentRule.  Despite this, the ability to supply the values at creation
 * time is an important aid to creating valid SBML.  Knowledge of the intended
 * SBML Level and Version determine whether it is valid to assign a particular
 * value to an attribute, or whether it is valid to add an object to an existing
 * SBMLDocument.
 */
LIBSBML_EXTERN
Rule_t *
Rule_createAssignmentWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    AssignmentRule* obj = new AssignmentRule(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new RateRule (Rule_t) structure using the given SBML
 * @p level and @p version values.
 *
 * @param level an unsigned int, the SBML Level to assign to this
 * RateRule
 *
 * @param version an unsigned int, the SBML Version to assign to this
 * RateRule
 *
 * @return a pointer to the newly created Rule_t structure.
 *
 * @note Once a RateRule has been added to an SBMLDocument, the @p
 * level and @p version for the document @em override those used to create
 * the RateRule.  Despite this, the ability to supply the values at
 * creation time is an important aid to creating valid SBML.  Knowledge of
 * the intended SBML Level and Version  determine whether it is valid to
 * assign a particular value to an attribute, or whether it is valid to add
 * an object to an existing SBMLDocument.
 */
LIBSBML_EXTERN
Rule_t *
Rule_createRate (unsigned int level, unsigned int version)
{
  try
  {
    RateRule* obj = new RateRule(level,version);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Creates a new RateRule (Rule_t) structure using the given
 * SBMLNamespaces_t structure.
 *
 * @param sbmlns SBMLNamespaces, a pointer to an SBMLNamespaces structure
 * to assign to this RateRule
 *
 * @return a pointer to the newly created Rule_t structure.
 *
 * @note Once a RateRule has been added to an SBMLDocument, the
 * @p sbmlns namespaces for the document @em override those used to create
 * the RateRule.  Despite this, the ability to supply the values at creation
 * time is an important aid to creating valid SBML.  Knowledge of the intended
 * SBML Level and Version determine whether it is valid to assign a particular
 * value to an attribute, or whether it is valid to add an object to an existing
 * SBMLDocument.
 */
LIBSBML_EXTERN
Rule_t *
Rule_createRateWithNS (SBMLNamespaces_t* sbmlns)
{
  try
  {
    RateRule* obj = new RateRule(sbmlns);
    return obj;
  }
  catch (SBMLConstructorException)
  {
    return NULL;
  }
}


/**
 * Destroys this Rule.
 */
LIBSBML_EXTERN
void
Rule_free (Rule_t *r)
{
  if (r != NULL)
  delete r;
}


/**
 * @return a (deep) copy of this Rule.
 */
LIBSBML_EXTERN
Rule_t *
Rule_clone (const Rule_t *r)
{
  return  (r != NULL) ? static_cast<Rule*>( r->clone() ) : NULL;
}


/**
 * Returns a list of XMLNamespaces_t associated with this Rule_t
 * structure.
 *
 * @param r the Rule_t structure
 * 
 * @return pointer to the XMLNamespaces_t structure associated with 
 * this SBML object
 */
LIBSBML_EXTERN
const XMLNamespaces_t *
Rule_getNamespaces(Rule_t *r)
{
  return (r != NULL) ? r->getNamespaces() : NULL;
}

/**
 * @note SBML Level 1 uses a text-string format for mathematical formulas.
 * SBML Level 2 uses MathML, an XML format for representing mathematical
 * expressions.  LibSBML provides an Abstract Syntax Tree API for working
 * with mathematical expressions; this API is more powerful than working
 * with formulas directly in text form, and ASTs can be translated into
 * either MathML or the text-string syntax.  The libSBML methods that
 * accept text-string formulas directly (such as this one) are
 * provided for SBML Level 1 compatibility, but developers are encouraged
 * to use the AST mechanisms.  
 *
 * @return the formula for this Rule.
 */
LIBSBML_EXTERN
const char *
Rule_getFormula (const Rule_t *r)
{
  return (r != NULL && r->isSetFormula()) ? r->getFormula().c_str() : NULL;
}


/**
 * @return the math for this Rule.
 */
LIBSBML_EXTERN
const ASTNode_t *
Rule_getMath (const Rule_t *r)
{
  return (r != NULL) ? r->getMath() : NULL;
}


/**
 * @return the type of this Rule, either RULE_TYPE_RATE or
 * RULE_TYPE_SCALAR.
 */
LIBSBML_EXTERN
RuleType_t
Rule_getType (const Rule_t *r)
{
  return (r != NULL) ? r->getType() : RULE_TYPE_INVALID ;
}


/**
 * @return the variable for this Rule.
 */
LIBSBML_EXTERN
const char *
Rule_getVariable (const Rule_t *r)
{
  return (r != NULL && r->isSetVariable()) ? r->getVariable().c_str() : NULL;
}


/**
 * @return the units for this Rule (L1 ParameterRules only).
 */
LIBSBML_EXTERN
const char *
Rule_getUnits (const Rule_t *r)
{
  return (r != NULL && r->isSetUnits()) ? r->getUnits().c_str() : NULL;
}


/**
 * @note SBML Level 1 uses a text-string format for mathematical formulas.
 * SBML Level 2 uses MathML, an XML format for representing mathematical
 * expressions.  LibSBML provides an Abstract Syntax Tree API for working
 * with mathematical expressions; this API is more powerful than working
 * with formulas directly in text form, and ASTs can be translated into
 * either MathML or the text-string syntax.  The libSBML methods that
 * accept text-string formulas directly (such as this one) are
 * provided for SBML Level 1 compatibility, but developers are encouraged
 * to use the AST mechanisms.
 *
 * @return true (non-zero) if the formula (or equivalently the math) for
 * this Rule is set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Rule_isSetFormula (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isSetFormula() ) : 0;
}


/**
 * @return true (non-zero) if the math (or equivalently the formula) for
 * this Rule is set, false (0) otherwise.
 */
LIBSBML_EXTERN
int
Rule_isSetMath (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isSetMath() ) : 0;
}


/**
 * @return true (non-zero) if the variable of this Rule is set, false
 * (0) otherwise.
 */
LIBSBML_EXTERN
int
Rule_isSetVariable (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isSetVariable() ) : 0;
}


/**
 * @return true (non-zero) if the units for this Rule is set, false
 * (0) otherwise (L1 ParameterRules only).
 */
LIBSBML_EXTERN
int
Rule_isSetUnits (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isSetUnits() ) : 0;
}


/**
 * Sets the formula of this Rule to a copy of string.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_OBJECT
 *
 * @note SBML Level 1 uses a text-string format for mathematical formulas.
 * SBML Level 2 uses MathML, an XML format for representing mathematical
 * expressions.  LibSBML provides an Abstract Syntax Tree API for working
 * with mathematical expressions; this API is more powerful than working
 * with formulas directly in text form, and ASTs can be translated into
 * either MathML or the text-string syntax.  The libSBML methods that
 * accept text-string formulas directly (such as this one) are
 * provided for SBML Level 1 compatibility, but developers are encouraged
 * to use the AST mechanisms.
 */
LIBSBML_EXTERN
int
Rule_setFormula (Rule_t *r, const char *formula)
{
  if (r != NULL)
    return (formula == NULL) ? r->setMath(NULL) : r->setFormula(formula);
  else
    return LIBSBML_INVALID_OBJECT;
}


/**
 * Sets the math of this Rule to a copy of the given ASTNode.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_OBJECT
 */
LIBSBML_EXTERN
int
Rule_setMath (Rule_t *r, const ASTNode_t *math)
{
  if (r != NULL)
    return r->setMath(math);
  else
    return LIBSBML_INVALID_OBJECT;
}


/**
 * Sets the variable of this RateRule to a copy of sid.
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "variable" attribute.
 */
LIBSBML_EXTERN
int
Rule_setVariable (Rule_t *r, const char *sid)
{
  if (r != NULL)
    return (sid == NULL) ? r->setVariable("") : r->setVariable(sid);
  else
    return LIBSBML_INVALID_OBJECT;
}


/**
 * Sets the units for this Rule to a copy of sname (L1 ParameterRules
 * only).
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 * @li LIBSBML_UNEXPECTED_ATTRIBUTE
 *
 * @note Using this function with an id of NULL is equivalent to
 * unsetting the "units" attribute.
 */
LIBSBML_EXTERN
int
Rule_setUnits (Rule_t *r, const char *sname)
{
  if (r != NULL)
    return (sname == NULL) ? r->unsetUnits() : r->setUnits(sname);
  else
    return LIBSBML_INVALID_OBJECT;
}


/**
 * Unsets the units for this Rule (L1 ParameterRules only).
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_OPERATION_FAILED
 */
LIBSBML_EXTERN
int
Rule_unsetUnits (Rule_t *r)
{
  if (r != NULL)
    return r->unsetUnits();
  else
    return LIBSBML_INVALID_OBJECT;
}


/**
 * @return true (non-zero) if this Rule is an AlgebraicRule, false (0)
 * otherwise.
 */
LIBSBML_EXTERN
int
Rule_isAlgebraic (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isAlgebraic() ) : 0;
}


/**
 * @return true (non-zero) if this Rule is an AssignmentRule, false (0)
 * otherwise.
 */
LIBSBML_EXTERN
int
Rule_isAssignment (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isAssignment() ) : 0;
}


/**
 * This method attempts to lookup the Rule's variable in the Model's list
 * of Compartments.
 *
 * @return true (non-zero) if this Rule is a CompartmentVolumeRule, false
 * (0) otherwise.
 */
LIBSBML_EXTERN
int
Rule_isCompartmentVolume (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isCompartmentVolume() ) : 0;
}


/**
 * This method attempts to lookup the Rule's variable in the Model's list
 * of Parameters.
 *
 * @return true (non-zero) if this Rule is a ParameterRule, false (0)
 * otherwise.
 */
LIBSBML_EXTERN
int
Rule_isParameter (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isParameter() ) : 0;
}


/**
 * @return true (non-zero) if this Rule is a RateRule (L2) or has
 * type="rate" (L1), false (0) otherwise.
 */
LIBSBML_EXTERN
int
Rule_isRate (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isRate() ) : 0;
}


/**
 * @return true (non-zero) if this Rule is an AssignmentRule (L2) has
 * type="scalar" (L1), false (0) otherwise.
 */
LIBSBML_EXTERN
int
Rule_isScalar (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isScalar() ) : 0;
}


/**
 * This method attempts to lookup the Rule's variable in the Model's list
 * of Species.
 *
 * @return true (non-zero) if this Rule is a SpeciesConcentrationRule, false
 * (0) otherwise.
 */
LIBSBML_EXTERN
int
Rule_isSpeciesConcentration (const Rule_t *r)
{
  return (r != NULL) ? static_cast<int>( r->isSpeciesConcentration() ) : 0;
}


/**
 * @return the typecode (int) of SBML objects contained in this ListOf or
 * (default).
 */
LIBSBML_EXTERN
int
Rule_getTypeCode (const Rule_t *r)
{
  return (r != NULL) ? r->getTypeCode() : SBML_UNKNOWN;
}


/**
 * @return the SBML Level 1 typecode for this Rule or SBML_UNKNOWN
 * (default).
 */
LIBSBML_EXTERN
int
Rule_getL1TypeCode (const Rule_t *r)
{
  return (r != NULL) ? r->getL1TypeCode() : SBML_UNKNOWN;
}

/**
 * Sets the SBML Level&nbsp;1 typecode for this Rule.
 *
 * @param r the Rule_t structure
 * @param type the SBML Level&nbsp;1 typecode for this Rule
 * (@c SBML_COMPARTMENT_VOLUME_RULE, @c SBML_PARAMETER_RULE,
 * or @c SBML_SPECIES_CONCENTRATION_RULE).
 *
 * @return integer value indicating success/failure of the
 * function.  @if clike The value is drawn from the
 * enumeration #OperationReturnValues_t. @endif The possible values
 * returned by this function are:
 * @li LIBSBML_OPERATION_SUCCESS
 * @li LIBSBML_INVALID_ATTRIBUTE_VALUE
 *
 */
LIBSBML_EXTERN
int
Rule_setL1TypeCode (Rule_t *r, int L1Type)
{
  return (r != NULL) ? r->setL1TypeCode(L1Type) : LIBSBML_INVALID_OBJECT;
}

/**
  * Calculates and returns a UnitDefinition_t that expresses the units
  * returned by the math expression of this Rule_t.
  *
  * @return a UnitDefinition_t that expresses the units of the math 
  * expression of this Rule_t.
  *
  * Note that the functionality that facilitates unit analysis depends 
  * on the model as a whole.  Thus, in cases where the object has not 
  * been added to a model or the model itself is incomplete,
  * unit analysis is not possible and this method will return @c NULL.
  *
  * @note The units are calculated by applying the mathematics 
  * from the expression to the units of the <ci> elements used 
  * within the expression. Where there are parameters/numbers
  * with undeclared units the UnitDefinition_t returned by this
  * function may not accurately represent the units of the expression.
  * 
  * @see Rule_containsUndeclaredUnits()
  */
LIBSBML_EXTERN
UnitDefinition_t * 
Rule_getDerivedUnitDefinition(Rule_t *r)
{
  return (r != NULL) ? r->getDerivedUnitDefinition() : NULL;
}


/**
  * Predicate returning @c true or @c false depending on whether 
  * the math expression of this Rule_t contains
  * parameters/numbers with undeclared units.
  * 
  * @return @c true if the math expression of this Rule_t
  * includes parameters/numbers 
  * with undeclared units, @c false otherwise.
  *
  * @note a return value of @c true indicates that the UnitDefinition_t
  * returned by the getDerivedUnitDefinition function may not 
  * accurately represent the units of the expression.
  *
  * @see Rule_getDerivedUnitDefinition()
  */
LIBSBML_EXTERN
int 
Rule_containsUndeclaredUnits(Rule_t *r)
{
  return (r != NULL) ? static_cast<int>(r->containsUndeclaredUnits()) : 0;
}


/**
 * @return item in this ListOfRule with the given id or @c NULL if no such
 * item exists.
 */
LIBSBML_EXTERN
Rule_t *
ListOfRules_getById (ListOf_t *lo, const char *sid)
{
  if (lo != NULL)
      return (sid != NULL) ? 
      static_cast <ListOfRules *> (lo)->get(sid) : NULL;
  else
    return NULL;
}


/**
 * Removes item in this ListOf items with the given id or @c NULL if no such
 * item exists.  The caller owns the returned item and is responsible for
 * deleting it.
 */
LIBSBML_EXTERN
Rule_t *
ListOfRules_removeById (ListOf_t *lo, const char *sid)
{
  if (lo != NULL)
    return (sid != NULL) ? 
      static_cast <ListOfRules *> (lo)->remove(sid) : NULL;
  else
    return NULL;
}

/** @endcond */
LIBSBML_CPP_NAMESPACE_END
