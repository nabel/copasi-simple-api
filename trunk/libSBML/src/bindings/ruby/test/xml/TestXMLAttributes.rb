# @file    TestXMLAttributes.rb
# @brief   TestXMLAttributes unit tests
#
# @author  Akiya Jouraku (Ruby conversion)
# @author  Ben Bornstein 
#
# $Id$
# $HeadURL$
#
# ====== WARNING ===== WARNING ===== WARNING ===== WARNING ===== WARNING ======
#
# DO NOT EDIT THIS FILE.
#
# This file was generated automatically by converting the file located at
# src/xml/test/TestXMLAttributes.cpp
# using the conversion program dev/utilities/translateTests/translateTests.pl.
# Any changes made here will be lost the next time the file is regenerated.
#
# -----------------------------------------------------------------------------
# This file is part of libSBML.  Please visit http://sbml.org for more
# information about SBML, and the latest version of libSBML.
#
# Copyright 2005-2010 California Institute of Technology.
# Copyright 2002-2005 California Institute of Technology and
#                     Japan Science and Technology Corporation.
# 
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation.  A copy of the license agreement is provided
# in the file named "LICENSE.txt" included with this software distribution
# and also available online as http://sbml.org/software/libsbml/license.html
# -----------------------------------------------------------------------------
require 'test/unit'
require 'libSBML'

class TestXMLAttributes < Test::Unit::TestCase

  def util_NaN
    z = 0.0
    return 0.0/z
  end

  def util_PosInf
    z = 0.0
    return 1.0/z
  end

  def util_NegInf
    z = 0.0
    return -1.0/z
  end

  def equals(*x)
    case x.size
    when 2
      e, s = x
      return e == s
    when 1
      e, = x
      return e == @@oss.str()
    end
  end

  def test_XMLAttributes_add_get
    attrs = LibSBML::XMLAttributes.new()
    assert( attrs.getLength() == 0 )
    assert_equal true, attrs.isEmpty()
    attrs.add("xmlns", "http://foo.org/")
    assert( attrs.getLength() == 1 )
    assert( attrs.isEmpty() == false )
    attrs.add("foo", "bar")
    assert( attrs.getLength() == 2 )
    assert( attrs.isEmpty() == false )
    assert( attrs.getIndex("xmlns") == 0 )
    assert( attrs.getIndex("foo"  ) == 1 )
    assert( attrs.getIndex("bar"  ) == -1 )
    assert( attrs.getValue("xmlns") ==  "http://foo.org/"  )
    assert( attrs.getValue("foo"  ) ==  "bar"              )
    assert( attrs.getValue("bar"  ) ==  ""                 )
    assert( attrs.getName(0) ==  "xmlns"  )
    assert( attrs.getName(1) ==  "foo"    )
    assert( attrs.getName(2) ==  ""       )
  end

  def test_XMLAttributes_assignment
    att1 = LibSBML::XMLAttributes.new()
    att1.add("xmlns", "http://foo.org/")
    assert( att1.getLength() == 1 )
    assert( att1.isEmpty() == false )
    assert( att1.getIndex("xmlns") == 0 )
    assert( att1.getName(0) ==   "xmlns"  )
    assert( att1.getValue("xmlns") ==  "http://foo.org/"  )
    att2 = LibSBML::XMLAttributes.new()
    att2 = att1
    assert( att2.getLength() == 1 )
    assert( att2.isEmpty() == false )
    assert( att2.getIndex("xmlns") == 0 )
    assert( att2.getName(0) ==   "xmlns"  )
    assert( att2.getValue("xmlns") ==  "http://foo.org/"  )
    att2 = nil
    att1 = nil
  end

  def test_XMLAttributes_clone
    att1 = LibSBML::XMLAttributes.new()
    att1.add("xmlns", "http://foo.org/")
    assert( att1.getLength() == 1 )
    assert( att1.isEmpty() == false )
    assert( att1.getIndex("xmlns") == 0 )
    assert( att1.getName(0) ==   "xmlns"  )
    assert( att1.getValue("xmlns") ==  "http://foo.org/"  )
    att2 = att1.clone()
    assert( att2.getLength() == 1 )
    assert( att2.isEmpty() == false )
    assert( att2.getIndex("xmlns") == 0 )
    assert( att2.getName(0) ==   "xmlns"  )
    assert( att2.getValue("xmlns") ==  "http://foo.org/"  )
    att2 = nil
    att1 = nil
  end

  def test_XMLAttributes_copy
    att1 = LibSBML::XMLAttributes.new()
    att1.add("xmlns", "http://foo.org/")
    assert( att1.getLength() == 1 )
    assert( att1.isEmpty() == false )
    assert( att1.getIndex("xmlns") == 0 )
    assert( att1.getName(0) ==   "xmlns"  )
    assert( att1.getValue("xmlns") ==  "http://foo.org/"  )
    att2 = LibSBML::XMLAttributes.new(att1)
    assert( att2.getLength() == 1 )
    assert( att2.isEmpty() == false )
    assert( att2.getIndex("xmlns") == 0 )
    assert( att2.getName(0) ==   "xmlns"  )
    assert( att2.getValue("xmlns") ==  "http://foo.org/"  )
    att2 = nil
    att1 = nil
  end

end