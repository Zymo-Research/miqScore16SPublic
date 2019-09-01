from . import masterTable
from . import pipelineParameters

__all__ = ["masterTable",
           "pipelineParameters"]




def crossValidationPassed(sampleData: masterTable.MasterTable, parameters: pipelineParameters.PipelineParameters):
    for line in sampleData:
        if not line.seqType in parameters:
            return False
    return True