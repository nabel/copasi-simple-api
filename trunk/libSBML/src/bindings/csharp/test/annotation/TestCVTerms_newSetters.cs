///  @file    TestCVTerms_newSetters.cs
///  @brief   CVTerms unit tests
///  @author  Frank Bergmann (Csharp conversion)
///  @author  Akiya Jouraku (Csharp conversion)
///  @author  Sarah Keating 
/// 
///  $Id: TestCVTerms_newSetters.cs 11545 2010-07-23 02:19:10Z mhucka $
///  $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/bindings/csharp/test/annotation/TestCVTerms_newSetters.cs $
/// 
///  ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
/// 
///  DO NOT EDIT THIS FILE.
/// 
///  This file was generated automatically by converting the file located at
///  src/annotation/test/TestCVTerms_newSetters.c
///  using the conversion program dev/utilities/translateTests/translateTests.pl.
///  Any changes made here will be lost the next time the file is regenerated.
/// 
///  -----------------------------------------------------------------------------
///  This file is part of libSBML.  Please visit http://sbml.org for more
///  information about SBML, and the latest version of libSBML.
/// 
///  Copyright 2005-2010 California Institute of Technology.
///  Copyright 2002-2005 California Institute of Technology and
///                      Japan Science and Technology Corporation.
///  
///  This library is free software; you can redistribute it and/or modify it
///  under the terms of the GNU Lesser General Public License as published by
///  the Free Software Foundation.  A copy of the license agreement is provided
///  in the file named "LICENSE.txt" included with this software distribution
///  and also available online as http://sbml.org/software/libsbml/license.html
///  -----------------------------------------------------------------------------


namespace LibSBMLCSTest {

  using libsbml;

  using System;

  using System.IO;

  public class TestCVTerms_newSetters {
    public class AssertionError : System.Exception 
    {
      public AssertionError() : base()
      {
        
      }
    }


    static void assertTrue(bool condition)
    {
      if (condition == true)
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        return;
      }
      else if ( (a == null) || (b == null) )
      {
        throw new AssertionError();
      }
      else if (a.Equals(b))
      {
        return;
      }
  
      throw new AssertionError();
    }

    static void assertNotEquals(object a, object b)
    {
      if ( (a == null) && (b == null) )
      {
        throw new AssertionError();
      }
      else if ( (a == null) || (b == null) )
      {
        return;
      }
      else if (a.Equals(b))
      {
        throw new AssertionError();
      }
    }

    static void assertEquals(bool a, bool b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(bool a, bool b)
    {
      if ( a != b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertEquals(int a, int b)
    {
      if ( a == b )
      {
        return;
      }
      throw new AssertionError();
    }

    static void assertNotEquals(int a, int b)
    {
      if ( a != b )
      {
        return;
      }
      throw new AssertionError();
    }


    public void test_CVTerm_addResource()
    {
      CVTerm term = new  CVTerm(libsbml.MODEL_QUALIFIER);
      string resource =  "GO6666";
      XMLAttributes xa;
      assertTrue( term != null );
      assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
      int i = term.addResource( "");
      assertTrue( i == libsbml.LIBSBML_OPERATION_FAILED );
      xa = term.getResources();
      assertTrue( xa.getLength() == 0 );
      i = term.addResource(resource);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      xa = term.getResources();
      assertTrue( xa.getLength() == 1 );
      assertTrue((  "rdf:resource" == xa.getName(0) ));
      assertTrue((  "GO6666" == xa.getValue(0) ));
      term = null;
    }

    public void test_CVTerm_removeResource()
    {
      CVTerm term = new  CVTerm(libsbml.MODEL_QUALIFIER);
      string resource =  "GO6666";
      XMLAttributes xa;
      assertTrue( term != null );
      assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
      term.addResource(resource);
      xa = term.getResources();
      assertTrue( xa.getLength() == 1 );
      int i = term.removeResource( "CCC");
      assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
      xa = term.getResources();
      assertTrue( xa.getLength() == 1 );
      i = term.removeResource(resource);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      xa = term.getResources();
      assertTrue( xa.getLength() == 0 );
      term = null;
    }

    public void test_CVTerm_setBiolQualifierType()
    {
      CVTerm term = new  CVTerm(libsbml.BIOLOGICAL_QUALIFIER);
      assertTrue( term != null );
      assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
      int i = term.setBiologicalQualifierType(libsbml.BQB_IS);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_IS );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
      i = term.setQualifierType(libsbml.MODEL_QUALIFIER);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
      i = term.setBiologicalQualifierType(libsbml.BQB_IS);
      assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
      assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
      term = null;
    }

    public void test_CVTerm_setModelQualifierType()
    {
      CVTerm term = new  CVTerm(libsbml.MODEL_QUALIFIER);
      assertTrue( term != null );
      assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
      int i = term.setModelQualifierType(libsbml.BQM_IS);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( term.getQualifierType() == libsbml.MODEL_QUALIFIER );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_IS );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
      i = term.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER);
      assertTrue( i == libsbml.LIBSBML_OPERATION_SUCCESS );
      assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
      i = term.setModelQualifierType(libsbml.BQM_IS);
      assertTrue( i == libsbml.LIBSBML_INVALID_ATTRIBUTE_VALUE );
      assertTrue( term.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER );
      assertTrue( term.getBiologicalQualifierType() == libsbml.BQB_UNKNOWN );
      assertTrue( term.getModelQualifierType() == libsbml.BQM_UNKNOWN );
      term = null;
    }

  }
}
