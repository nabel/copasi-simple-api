# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.40
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_copasi', [dirname(__file__)])
        except ImportError:
            import _copasi
            return _copasi
        if fp is not None:
            try:
                _mod = imp.load_module('_copasi', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _copasi = swig_import_helper()
    del swig_import_helper
else:
    import _copasi
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class c_strings(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, c_strings, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, c_strings, name)
    __repr__ = _swig_repr
    __swig_setmethods__["length"] = _copasi.c_strings_length_set
    __swig_getmethods__["length"] = _copasi.c_strings_length_get
    if _newclass:length = _swig_property(_copasi.c_strings_length_get, _copasi.c_strings_length_set)
    __swig_setmethods__["strings"] = _copasi.c_strings_strings_set
    __swig_getmethods__["strings"] = _copasi.c_strings_strings_get
    if _newclass:strings = _swig_property(_copasi.c_strings_strings_get, _copasi.c_strings_strings_set)
    def __init__(self): 
        this = _copasi.new_c_strings()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _copasi.delete_c_strings
    __del__ = lambda self : None;
c_strings_swigregister = _copasi.c_strings_swigregister
c_strings_swigregister(c_strings)

class c_matrix(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, c_matrix, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, c_matrix, name)
    __repr__ = _swig_repr
    __swig_setmethods__["rows"] = _copasi.c_matrix_rows_set
    __swig_getmethods__["rows"] = _copasi.c_matrix_rows_get
    if _newclass:rows = _swig_property(_copasi.c_matrix_rows_get, _copasi.c_matrix_rows_set)
    __swig_setmethods__["cols"] = _copasi.c_matrix_cols_set
    __swig_getmethods__["cols"] = _copasi.c_matrix_cols_get
    if _newclass:cols = _swig_property(_copasi.c_matrix_cols_get, _copasi.c_matrix_cols_set)
    __swig_setmethods__["values"] = _copasi.c_matrix_values_set
    __swig_getmethods__["values"] = _copasi.c_matrix_values_get
    if _newclass:values = _swig_property(_copasi.c_matrix_values_get, _copasi.c_matrix_values_set)
    __swig_setmethods__["rownames"] = _copasi.c_matrix_rownames_set
    __swig_getmethods__["rownames"] = _copasi.c_matrix_rownames_get
    if _newclass:rownames = _swig_property(_copasi.c_matrix_rownames_get, _copasi.c_matrix_rownames_set)
    __swig_setmethods__["colnames"] = _copasi.c_matrix_colnames_set
    __swig_getmethods__["colnames"] = _copasi.c_matrix_colnames_get
    if _newclass:colnames = _swig_property(_copasi.c_matrix_colnames_get, _copasi.c_matrix_colnames_set)
    def __init__(self): 
        this = _copasi.new_c_matrix()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _copasi.delete_c_matrix
    __del__ = lambda self : None;
c_matrix_swigregister = _copasi.c_matrix_swigregister
c_matrix_swigregister(c_matrix)


def c_createMatrix(*args):
  return _copasi.c_createMatrix(*args)
c_createMatrix = _copasi.c_createMatrix

def c_createStringsArray(*args):
  return _copasi.c_createStringsArray(*args)
c_createStringsArray = _copasi.c_createStringsArray

def c_getMatrixValue(*args):
  return _copasi.c_getMatrixValue(*args)
c_getMatrixValue = _copasi.c_getMatrixValue

def c_setMatrixValue(*args):
  return _copasi.c_setMatrixValue(*args)
c_setMatrixValue = _copasi.c_setMatrixValue

def c_getRowName(*args):
  return _copasi.c_getRowName(*args)
c_getRowName = _copasi.c_getRowName

def c_setRowName(*args):
  return _copasi.c_setRowName(*args)
c_setRowName = _copasi.c_setRowName

def c_getColumnName(*args):
  return _copasi.c_getColumnName(*args)
c_getColumnName = _copasi.c_getColumnName

def c_setColumnName(*args):
  return _copasi.c_setColumnName(*args)
c_setColumnName = _copasi.c_setColumnName

def c_getString(*args):
  return _copasi.c_getString(*args)
c_getString = _copasi.c_getString

def c_setString(*args):
  return _copasi.c_setString(*args)
c_setString = _copasi.c_setString

def c_getStringIndex(*args):
  return _copasi.c_getStringIndex(*args)
c_getStringIndex = _copasi.c_getStringIndex

def c_getRowIndex(*args):
  return _copasi.c_getRowIndex(*args)
c_getRowIndex = _copasi.c_getRowIndex

def c_getColumnIndex(*args):
  return _copasi.c_getColumnIndex(*args)
c_getColumnIndex = _copasi.c_getColumnIndex

def c_deleteMatrix(*args):
  return _copasi.c_deleteMatrix(*args)
c_deleteMatrix = _copasi.c_deleteMatrix

def c_deleteStringsArray(*args):
  return _copasi.c_deleteStringsArray(*args)
c_deleteStringsArray = _copasi.c_deleteStringsArray

def c_appendColumns(*args):
  return _copasi.c_appendColumns(*args)
c_appendColumns = _copasi.c_appendColumns

def c_appendRows(*args):
  return _copasi.c_appendRows(*args)
c_appendRows = _copasi.c_appendRows

def c_printMatrixToFile(*args):
  return _copasi.c_printMatrixToFile(*args)
c_printMatrixToFile = _copasi.c_printMatrixToFile

def c_printOutMatrix(*args):
  return _copasi.c_printOutMatrix(*args)
c_printOutMatrix = _copasi.c_printOutMatrix
class copasi_model(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, copasi_model, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, copasi_model, name)
    __repr__ = _swig_repr
    __swig_setmethods__["CopasiModelPtr"] = _copasi.copasi_model_CopasiModelPtr_set
    __swig_getmethods__["CopasiModelPtr"] = _copasi.copasi_model_CopasiModelPtr_get
    if _newclass:CopasiModelPtr = _swig_property(_copasi.copasi_model_CopasiModelPtr_get, _copasi.copasi_model_CopasiModelPtr_set)
    __swig_setmethods__["CopasiDataModelPtr"] = _copasi.copasi_model_CopasiDataModelPtr_set
    __swig_getmethods__["CopasiDataModelPtr"] = _copasi.copasi_model_CopasiDataModelPtr_get
    if _newclass:CopasiDataModelPtr = _swig_property(_copasi.copasi_model_CopasiDataModelPtr_get, _copasi.copasi_model_CopasiDataModelPtr_set)
    __swig_setmethods__["qHash"] = _copasi.copasi_model_qHash_set
    __swig_getmethods__["qHash"] = _copasi.copasi_model_qHash_get
    if _newclass:qHash = _swig_property(_copasi.copasi_model_qHash_get, _copasi.copasi_model_qHash_set)
    __swig_setmethods__["errorMessage"] = _copasi.copasi_model_errorMessage_set
    __swig_getmethods__["errorMessage"] = _copasi.copasi_model_errorMessage_get
    if _newclass:errorMessage = _swig_property(_copasi.copasi_model_errorMessage_get, _copasi.copasi_model_errorMessage_set)
    __swig_setmethods__["warningMessage"] = _copasi.copasi_model_warningMessage_set
    __swig_getmethods__["warningMessage"] = _copasi.copasi_model_warningMessage_get
    if _newclass:warningMessage = _swig_property(_copasi.copasi_model_warningMessage_get, _copasi.copasi_model_warningMessage_set)
    def __init__(self): 
        this = _copasi.new_copasi_model()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _copasi.delete_copasi_model
    __del__ = lambda self : None;
copasi_model_swigregister = _copasi.copasi_model_swigregister
copasi_model_swigregister(copasi_model)

class copasi_reaction(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, copasi_reaction, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, copasi_reaction, name)
    __repr__ = _swig_repr
    __swig_setmethods__["CopasiReactionPtr"] = _copasi.copasi_reaction_CopasiReactionPtr_set
    __swig_getmethods__["CopasiReactionPtr"] = _copasi.copasi_reaction_CopasiReactionPtr_get
    if _newclass:CopasiReactionPtr = _swig_property(_copasi.copasi_reaction_CopasiReactionPtr_get, _copasi.copasi_reaction_CopasiReactionPtr_set)
    __swig_setmethods__["CopasiModelPtr"] = _copasi.copasi_reaction_CopasiModelPtr_set
    __swig_getmethods__["CopasiModelPtr"] = _copasi.copasi_reaction_CopasiModelPtr_get
    if _newclass:CopasiModelPtr = _swig_property(_copasi.copasi_reaction_CopasiModelPtr_get, _copasi.copasi_reaction_CopasiModelPtr_set)
    __swig_setmethods__["qHash"] = _copasi.copasi_reaction_qHash_set
    __swig_getmethods__["qHash"] = _copasi.copasi_reaction_qHash_get
    if _newclass:qHash = _swig_property(_copasi.copasi_reaction_qHash_get, _copasi.copasi_reaction_qHash_set)
    def __init__(self): 
        this = _copasi.new_copasi_reaction()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _copasi.delete_copasi_reaction
    __del__ = lambda self : None;
copasi_reaction_swigregister = _copasi.copasi_reaction_swigregister
copasi_reaction_swigregister(copasi_reaction)

class copasi_compartment(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, copasi_compartment, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, copasi_compartment, name)
    __repr__ = _swig_repr
    __swig_setmethods__["CopasiCompartmentPtr"] = _copasi.copasi_compartment_CopasiCompartmentPtr_set
    __swig_getmethods__["CopasiCompartmentPtr"] = _copasi.copasi_compartment_CopasiCompartmentPtr_get
    if _newclass:CopasiCompartmentPtr = _swig_property(_copasi.copasi_compartment_CopasiCompartmentPtr_get, _copasi.copasi_compartment_CopasiCompartmentPtr_set)
    __swig_setmethods__["CopasiModelPtr"] = _copasi.copasi_compartment_CopasiModelPtr_set
    __swig_getmethods__["CopasiModelPtr"] = _copasi.copasi_compartment_CopasiModelPtr_get
    if _newclass:CopasiModelPtr = _swig_property(_copasi.copasi_compartment_CopasiModelPtr_get, _copasi.copasi_compartment_CopasiModelPtr_set)
    __swig_setmethods__["qHash"] = _copasi.copasi_compartment_qHash_set
    __swig_getmethods__["qHash"] = _copasi.copasi_compartment_qHash_get
    if _newclass:qHash = _swig_property(_copasi.copasi_compartment_qHash_get, _copasi.copasi_compartment_qHash_set)
    def __init__(self): 
        this = _copasi.new_copasi_compartment()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _copasi.delete_copasi_compartment
    __del__ = lambda self : None;
copasi_compartment_swigregister = _copasi.copasi_compartment_swigregister
copasi_compartment_swigregister(copasi_compartment)


def copasi_init():
  return _copasi.copasi_init()
copasi_init = _copasi.copasi_init

def copasi_end():
  return _copasi.copasi_end()
copasi_end = _copasi.copasi_end

def cRemoveModel(*args):
  return _copasi.cRemoveModel(*args)
cRemoveModel = _copasi.cRemoveModel

def cSetSBMLLevelAndVersion(*args):
  return _copasi.cSetSBMLLevelAndVersion(*args)
cSetSBMLLevelAndVersion = _copasi.cSetSBMLLevelAndVersion

def cReadAntimonyFile(*args):
  return _copasi.cReadAntimonyFile(*args)
cReadAntimonyFile = _copasi.cReadAntimonyFile

def cReadAntimonyString(*args):
  return _copasi.cReadAntimonyString(*args)
cReadAntimonyString = _copasi.cReadAntimonyString

def cReadSBMLFile(*args):
  return _copasi.cReadSBMLFile(*args)
cReadSBMLFile = _copasi.cReadSBMLFile

def cReadSBMLString(*args):
  return _copasi.cReadSBMLString(*args)
cReadSBMLString = _copasi.cReadSBMLString

def cWriteSBMLFile(*args):
  return _copasi.cWriteSBMLFile(*args)
cWriteSBMLFile = _copasi.cWriteSBMLFile

def cWriteAntimonyFile(*args):
  return _copasi.cWriteAntimonyFile(*args)
cWriteAntimonyFile = _copasi.cWriteAntimonyFile

def cCreateModel(*args):
  return _copasi.cCreateModel(*args)
cCreateModel = _copasi.cCreateModel

def cCompileModel(*args):
  return _copasi.cCompileModel(*args)
cCompileModel = _copasi.cCompileModel

def cCreateCompartment(*args):
  return _copasi.cCreateCompartment(*args)
cCreateCompartment = _copasi.cCreateCompartment

def cSetVolume(*args):
  return _copasi.cSetVolume(*args)
cSetVolume = _copasi.cSetVolume

def cGetCompartments(*args):
  return _copasi.cGetCompartments(*args)
cGetCompartments = _copasi.cGetCompartments

def cSetCompartmentVolumes(*args):
  return _copasi.cSetCompartmentVolumes(*args)
cSetCompartmentVolumes = _copasi.cSetCompartmentVolumes

def cGetNumberOfCompartments(*args):
  return _copasi.cGetNumberOfCompartments(*args)
cGetNumberOfCompartments = _copasi.cGetNumberOfCompartments

def cSetAssignmentRule(*args):
  return _copasi.cSetAssignmentRule(*args)
cSetAssignmentRule = _copasi.cSetAssignmentRule

def cCreateVariable(*args):
  return _copasi.cCreateVariable(*args)
cCreateVariable = _copasi.cCreateVariable

def cCreateEvent(*args):
  return _copasi.cCreateEvent(*args)
cCreateEvent = _copasi.cCreateEvent

def cCreateReaction(*args):
  return _copasi.cCreateReaction(*args)
cCreateReaction = _copasi.cCreateReaction

def cAddReactant(*args):
  return _copasi.cAddReactant(*args)
cAddReactant = _copasi.cAddReactant

def cAddProduct(*args):
  return _copasi.cAddProduct(*args)
cAddProduct = _copasi.cAddProduct

def cSetReactionRate(*args):
  return _copasi.cSetReactionRate(*args)
cSetReactionRate = _copasi.cSetReactionRate

def cGetReactionRate(*args):
  return _copasi.cGetReactionRate(*args)
cGetReactionRate = _copasi.cGetReactionRate

def cGetReactionRatesEx(*args):
  return _copasi.cGetReactionRatesEx(*args)
cGetReactionRatesEx = _copasi.cGetReactionRatesEx

def cGetFlux(*args):
  return _copasi.cGetFlux(*args)
cGetFlux = _copasi.cGetFlux

def cGetParticleFlux(*args):
  return _copasi.cGetParticleFlux(*args)
cGetParticleFlux = _copasi.cGetParticleFlux

def cGetReactionRates(*args):
  return _copasi.cGetReactionRates(*args)
cGetReactionRates = _copasi.cGetReactionRates

def cGetNumberOfReactions(*args):
  return _copasi.cGetNumberOfReactions(*args)
cGetNumberOfReactions = _copasi.cGetNumberOfReactions

def cGetNumberOfSpecies(*args):
  return _copasi.cGetNumberOfSpecies(*args)
cGetNumberOfSpecies = _copasi.cGetNumberOfSpecies

def cGetNumberOfFloatingSpecies(*args):
  return _copasi.cGetNumberOfFloatingSpecies(*args)
cGetNumberOfFloatingSpecies = _copasi.cGetNumberOfFloatingSpecies

def cGetNumberOfBoundarySpecies(*args):
  return _copasi.cGetNumberOfBoundarySpecies(*args)
cGetNumberOfBoundarySpecies = _copasi.cGetNumberOfBoundarySpecies

def cCreateSpecies(*args):
  return _copasi.cCreateSpecies(*args)
cCreateSpecies = _copasi.cCreateSpecies

def cSetSpeciesType(*args):
  return _copasi.cSetSpeciesType(*args)
cSetSpeciesType = _copasi.cSetSpeciesType

def cSetSpeciesConcentration(*args):
  return _copasi.cSetSpeciesConcentration(*args)
cSetSpeciesConcentration = _copasi.cSetSpeciesConcentration

def cSetInitialConcentration(*args):
  return _copasi.cSetInitialConcentration(*args)
cSetInitialConcentration = _copasi.cSetInitialConcentration

def cSetSpeciesAmount(*args):
  return _copasi.cSetSpeciesAmount(*args)
cSetSpeciesAmount = _copasi.cSetSpeciesAmount

def cSetSpeciesConcentrations(*args):
  return _copasi.cSetSpeciesConcentrations(*args)
cSetSpeciesConcentrations = _copasi.cSetSpeciesConcentrations

def cGetFloatingSpeciesIntitialConcentrations(*args):
  return _copasi.cGetFloatingSpeciesIntitialConcentrations(*args)
cGetFloatingSpeciesIntitialConcentrations = _copasi.cGetFloatingSpeciesIntitialConcentrations

def cSetFloatingSpeciesIntitialConcentrations(*args):
  return _copasi.cSetFloatingSpeciesIntitialConcentrations(*args)
cSetFloatingSpeciesIntitialConcentrations = _copasi.cSetFloatingSpeciesIntitialConcentrations

def cGetAllSpecies(*args):
  return _copasi.cGetAllSpecies(*args)
cGetAllSpecies = _copasi.cGetAllSpecies

def cGetFloatingSpeciesConcentrations(*args):
  return _copasi.cGetFloatingSpeciesConcentrations(*args)
cGetFloatingSpeciesConcentrations = _copasi.cGetFloatingSpeciesConcentrations

def cGetBoundarySpecies(*args):
  return _copasi.cGetBoundarySpecies(*args)
cGetBoundarySpecies = _copasi.cGetBoundarySpecies

def cGetAmounts(*args):
  return _copasi.cGetAmounts(*args)
cGetAmounts = _copasi.cGetAmounts

def cGetConcentration(*args):
  return _copasi.cGetConcentration(*args)
cGetConcentration = _copasi.cGetConcentration

def cGetAmount(*args):
  return _copasi.cGetAmount(*args)
cGetAmount = _copasi.cGetAmount

def cSetValue(*args):
  return _copasi.cSetValue(*args)
cSetValue = _copasi.cSetValue

def cSetGlobalParameter(*args):
  return _copasi.cSetGlobalParameter(*args)
cSetGlobalParameter = _copasi.cSetGlobalParameter

def cGetNumberOfGlobalParameters(*args):
  return _copasi.cGetNumberOfGlobalParameters(*args)
cGetNumberOfGlobalParameters = _copasi.cGetNumberOfGlobalParameters

def cGetGlobalParameters(*args):
  return _copasi.cGetGlobalParameters(*args)
cGetGlobalParameters = _copasi.cGetGlobalParameters

def cSetGlobalParameterValues(*args):
  return _copasi.cSetGlobalParameterValues(*args)
cSetGlobalParameterValues = _copasi.cSetGlobalParameterValues

def cSetValues(*args):
  return _copasi.cSetValues(*args)
cSetValues = _copasi.cSetValues

def cGetRatesOfChange(*args):
  return _copasi.cGetRatesOfChange(*args)
cGetRatesOfChange = _copasi.cGetRatesOfChange

def cGetRateOfChange(*args):
  return _copasi.cGetRateOfChange(*args)
cGetRateOfChange = _copasi.cGetRateOfChange

def cGetRatesOfChangeEx(*args):
  return _copasi.cGetRatesOfChangeEx(*args)
cGetRatesOfChangeEx = _copasi.cGetRatesOfChangeEx

def cSimulateDeterministic(*args):
  return _copasi.cSimulateDeterministic(*args)
cSimulateDeterministic = _copasi.cSimulateDeterministic

def cOneStep(*args):
  return _copasi.cOneStep(*args)
cOneStep = _copasi.cOneStep

def cSimulateStochastic(*args):
  return _copasi.cSimulateStochastic(*args)
cSimulateStochastic = _copasi.cSimulateStochastic

def cSimulateHybrid(*args):
  return _copasi.cSimulateHybrid(*args)
cSimulateHybrid = _copasi.cSimulateHybrid

def cSimulateTauLeap(*args):
  return _copasi.cSimulateTauLeap(*args)
cSimulateTauLeap = _copasi.cSimulateTauLeap

def cResetState(*args):
  return _copasi.cResetState(*args)
cResetState = _copasi.cResetState

def cGetReactionRatesFromTimeCourse(*args):
  return _copasi.cGetReactionRatesFromTimeCourse(*args)
cGetReactionRatesFromTimeCourse = _copasi.cGetReactionRatesFromTimeCourse

def cGetDerivativesFromTimeCourse(*args):
  return _copasi.cGetDerivativesFromTimeCourse(*args)
cGetDerivativesFromTimeCourse = _copasi.cGetDerivativesFromTimeCourse

def cGetCCFromTimeCourse(*args):
  return _copasi.cGetCCFromTimeCourse(*args)
cGetCCFromTimeCourse = _copasi.cGetCCFromTimeCourse

def cGetElasticitiesFromTimeCourse(*args):
  return _copasi.cGetElasticitiesFromTimeCourse(*args)
cGetElasticitiesFromTimeCourse = _copasi.cGetElasticitiesFromTimeCourse

def cFilterTimeCourseResults(*args):
  return _copasi.cFilterTimeCourseResults(*args)
cFilterTimeCourseResults = _copasi.cFilterTimeCourseResults

def cGetSteadyState(*args):
  return _copasi.cGetSteadyState(*args)
cGetSteadyState = _copasi.cGetSteadyState

def cGetSteadyStateUsingSimulation(*args):
  return _copasi.cGetSteadyStateUsingSimulation(*args)
cGetSteadyStateUsingSimulation = _copasi.cGetSteadyStateUsingSimulation

def cGetJacobian(*args):
  return _copasi.cGetJacobian(*args)
cGetJacobian = _copasi.cGetJacobian

def cGetEigenvalues(*args):
  return _copasi.cGetEigenvalues(*args)
cGetEigenvalues = _copasi.cGetEigenvalues

def cGetUnscaledFluxControlCoeffs(*args):
  return _copasi.cGetUnscaledFluxControlCoeffs(*args)
cGetUnscaledFluxControlCoeffs = _copasi.cGetUnscaledFluxControlCoeffs

def cGetScaledFluxControlCoeffs(*args):
  return _copasi.cGetScaledFluxControlCoeffs(*args)
cGetScaledFluxControlCoeffs = _copasi.cGetScaledFluxControlCoeffs

def cGetUnscaledConcentrationControlCoeffs(*args):
  return _copasi.cGetUnscaledConcentrationControlCoeffs(*args)
cGetUnscaledConcentrationControlCoeffs = _copasi.cGetUnscaledConcentrationControlCoeffs

def cGetScaledConcentrationConcentrationCoeffs(*args):
  return _copasi.cGetScaledConcentrationConcentrationCoeffs(*args)
cGetScaledConcentrationConcentrationCoeffs = _copasi.cGetScaledConcentrationConcentrationCoeffs

def cGetUnscaledElasticities(*args):
  return _copasi.cGetUnscaledElasticities(*args)
cGetUnscaledElasticities = _copasi.cGetUnscaledElasticities

def cGetScaledElasticities(*args):
  return _copasi.cGetScaledElasticities(*args)
cGetScaledElasticities = _copasi.cGetScaledElasticities

def cGetFullStoichiometryMatrix(*args):
  return _copasi.cGetFullStoichiometryMatrix(*args)
cGetFullStoichiometryMatrix = _copasi.cGetFullStoichiometryMatrix

def cGetReducedStoichiometryMatrix(*args):
  return _copasi.cGetReducedStoichiometryMatrix(*args)
cGetReducedStoichiometryMatrix = _copasi.cGetReducedStoichiometryMatrix

def cGetElementaryFluxModes(*args):
  return _copasi.cGetElementaryFluxModes(*args)
cGetElementaryFluxModes = _copasi.cGetElementaryFluxModes

def cGetGammaMatrix(*args):
  return _copasi.cGetGammaMatrix(*args)
cGetGammaMatrix = _copasi.cGetGammaMatrix

def cGetKMatrix(*args):
  return _copasi.cGetKMatrix(*args)
cGetKMatrix = _copasi.cGetKMatrix

def cGetK0Matrix(*args):
  return _copasi.cGetK0Matrix(*args)
cGetK0Matrix = _copasi.cGetK0Matrix

def cGetLinkMatrix(*args):
  return _copasi.cGetLinkMatrix(*args)
cGetLinkMatrix = _copasi.cGetLinkMatrix

def cGetL0Matrix(*args):
  return _copasi.cGetL0Matrix(*args)
cGetL0Matrix = _copasi.cGetL0Matrix

def cOptimize(*args):
  return _copasi.cOptimize(*args)
cOptimize = _copasi.cOptimize

def cSetOptimizerIterations(*args):
  return _copasi.cSetOptimizerIterations(*args)
cSetOptimizerIterations = _copasi.cSetOptimizerIterations

def cSetOptimizerSize(*args):
  return _copasi.cSetOptimizerSize(*args)
cSetOptimizerSize = _copasi.cSetOptimizerSize

def cSetOptimizerMutationRate(*args):
  return _copasi.cSetOptimizerMutationRate(*args)
cSetOptimizerMutationRate = _copasi.cSetOptimizerMutationRate

def cSetOptimizerCrossoverRate(*args):
  return _copasi.cSetOptimizerCrossoverRate(*args)
cSetOptimizerCrossoverRate = _copasi.cSetOptimizerCrossoverRate

def cDisableAssignmentRuleReordering():
  return _copasi.cDisableAssignmentRuleReordering()
cDisableAssignmentRuleReordering = _copasi.cDisableAssignmentRuleReordering

def cEnableAssignmentRuleReordering():
  return _copasi.cEnableAssignmentRuleReordering()
cEnableAssignmentRuleReordering = _copasi.cEnableAssignmentRuleReordering

def toItems(array):
    n = len(array);
    A = c_createItemsArray(n);
    for i in range(0, n):
        c_setItem(A, i, array[i]);

    return A;

def fromItems(array):
    n = array.length;
    A = range(0,n);
    for i in range(0, n):
        A[i] = c_getItem(array,i);

    #c_deleteItemsArray(array);
    return A;

def toStrings(array):
    n = len(array);
    A = c_createStringsArray(n);
    for i in range(0, n):
        c_setString(A, i, array[i]);

    return A;

def fromStrings(array):
    n = array.length;
    A = range(0,n);
    for i in range(0, n):
        A[i] = c_getString(array,i);

    #c_deleteStringsArray(array);
    return A;

def fromMatrix(matrix, row_wise = False):
    n = matrix.rows;
    m = matrix.cols;
    cols = fromStrings(matrix.colnames);
    rows = fromStrings(matrix.rownames);
    if row_wise:
        A = range(0,n);
        for i in range(0, n):
            A[i] = range(0,m);
            for j in range(0,m):
                A[i][j] = c_getMatrixValue(matrix,i,j);
    else:
        A = range(0,m);
        for i in range(0, m):
            A[i] = range(0,n);
            for j in range(0,n):
                A[i][j] = c_getMatrixValue(matrix,j,i);

    #c_deleteMatrix(matrix);
    return [rows, cols, A];

def toMatrix(lists, row_wise = False , rows = [], cols = []):
    n = len(lists);
    m = len(lists[0]);
    A = c_createMatrix(0,0);
    if row_wise:
        A = c_createMatrix(n,m);
    else:
        A = c_createMatrix(m,n);
    for i in range(0, n):
        for j in range(0,m):
            if row_wise:
                c_setMatrixValue(A,i,j,lists[i][j]);
            else:
                c_setMatrixValue(A,j,i,lists[i][j]);
    n = len(rows);
    m = len(cols);

    for i in range(0,n):
        c_setRowName(A,i,rows[i]);

    for i in range(0,m):
        c_setColumnName(A,i,cols[i]);

    return A;

def fromTable(table, row_wise = False):
    n = table.rows;
    m = table.cols;
    cols = fromStrings(table.colnames);
    rows = fromStrings(table.rownames);
    if row_wise:
        A = range(0,n);
        for i in range(0, n):
            A[i] = range(0,m);
            for j in range(0,m):
                A[i][j] = c_getTableValue(table,i,j);
    else:
        A = range(0,m);
        for i in range(0, m):
            A[i] = range(0,n);
            for j in range(0,n):
                A[i][j] = c_getTableValue(table,j,i);

    return [rows, cols, A];

def toTable(lists, row_wise = False , rows = [], cols = []):
    n = len(lists);
    m = len(lists[0]);
    
    A = c_createTable(0,0);
    if row_wise:
        A = c_createTable(n,m);
    else:
        A = c_createTable(m,n);

    for i in range(0, n):
        for j in range(0,m):
            if row_wise:
                c_setTableValue(A,i,j,lists[i][j]);
            else:
                c_setTableValue(A,j,i,lists[i][j]);
    n = len(rows);
    m = len(cols);

    for i in range(0,n):
        c_setString(A.rownames,i,rows[i]);

    for i in range(0,m):
        c_setString(A.colnames,i,cols[i]);

    return A;

