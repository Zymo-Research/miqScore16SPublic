import os
from pytest import mark

def setEnvironmentForTest():
    os.environ["INTEGER"] = "42"
    os.environ["FLOAT"] = "4.2"
    os.environ["STRING"] = "text"
    os.environ["NEGBOOL"] = "False"
    os.environ["POSBOOL"] = "True"

@mark.build
@mark.parameters
def test_parameterParse():
    from .. import parameters
    parameterSet = parameters.EnvParameters()
    parameterSet.addParameter("integer", int, default=0)
    parameterSet.addParameter("float", float, default=0.0)
    parameterSet.addParameter("string", str, default="override")
    parameterSet.addParameter("negbool", bool, default=True)
    parameterSet.addParameter("posbool", bool, default=False)
    assert parameterSet.integer == 42
    assert parameterSet.float == 4.2
    assert parameterSet.string == "text"
    assert parameterSet.negbool == False
    assert parameterSet.posbool == True
