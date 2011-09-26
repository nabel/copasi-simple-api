/**
 * @file    XMLHandler.cpp
 * @brief   XMLHandler interface
 * @author  Ben Bornstein
 *
 * $Id: XMLHandler.cpp 11633 2010-08-03 03:53:27Z mhucka $
 * $HeadURL: https://sbml.svn.sourceforge.net/svnroot/sbml/trunk/libsbml/src/xml/XMLHandler.cpp $
 *
 *<!---------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright 2005-2010 California Institute of Technology.
 * Copyright 2002-2005 California Institute of Technology and
 *                     Japan Science and Technology Corporation.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution and
 * also available online as http://sbml.org/software/libsbml/license.html
 *----------------------------------------------------------------------- -->*/

#include <sbml/xml/XMLHandler.h>
#include <sbml/xml/XMLToken.h>

/** @cond doxygen-ignored */

using namespace std;

/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN

/** @cond doxygen-libsbml-internal */

/*
 * Creates a new XMLHandler.
 */
XMLHandler::XMLHandler ()
{
}


/*
 * Destroys this XMLHandler.
 */
XMLHandler::~XMLHandler ()
{
}


/*
 * Receive notification of the beginning of the document.
 *
 * By default, do nothing. Application writers may override this method
 * in a subclass to take specific actions at the start of the document.
 */
void
XMLHandler::startDocument ()
{
}


/*
 * Receive notification of the XML declaration,
 * i.e.  <?xml version="1.0" encoding="UTF-8"?>
 *
 * By default, do nothing. Application writers may override this method in
 * a subclass to take specific actions at the declaration.
 */
void
XMLHandler::XML (const string& version, const string& encoding)
{
}


/*
 * Receive notification of the start of an element.
 *
 * By default, do nothing. Application writers may override this method
 * in a subclass to take specific actions at the start of each element.
 */
void
XMLHandler::startElement (const XMLToken& element)
{
}


/*
 * Receive notification of the end of the document.
 *
 * By default, do nothing. Application writers may override this method
 * in a subclass to take specific actions at the end of the document.
 */
void
XMLHandler::endDocument ()
{
}


/*
 * Receive notification of the end of an element.
 *
 * By default, do nothing. Application writers may override this method
 * in a subclass to take specific actions at the end of each element.
 */
void
XMLHandler::endElement (const XMLToken& element)
{
}


/*
 * Receive notification of character data inside an element.
 *
 * By default, do nothing. Application writers may override this method
 * to take specific actions for each chunk of character data.
 */
void
XMLHandler::characters (const XMLToken& data)
{
}


/** @endcond */

LIBSBML_CPP_NAMESPACE_END
