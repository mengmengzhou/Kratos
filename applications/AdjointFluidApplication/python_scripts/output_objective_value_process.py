import KratosMultiphysics
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if settings["process_name"].GetString() == "OutputObjectiveValueProcess2D":
        return OutputObjectiveValueProcess2D(Model, settings["Parameters"])
    elif settings["process_name"].GetString() == "OutputObjectiveValueProcess3D":
        return OutputObjectiveValueProcess3D(Model, settings["Parameters"])
    else:
        raise Exception("process_name %s not found." % settings["process_name"].GetString())

class OutputObjectiveValueProcess2D(AdjointFluidApplication.OutputObjectiveValueProcess2D):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        AdjointFluidApplication.OutputObjectiveValueProcess2D.__init__(self, self.model_part, settings["objective_settings"])

class OutputObjectiveValueProcess3D(AdjointFluidApplication.OutputObjectiveValueProcess3D):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        AdjointFluidApplication.OutputObjectiveValueProcess3D.__init__(self, self.model_part, settings["objective_settings"])
