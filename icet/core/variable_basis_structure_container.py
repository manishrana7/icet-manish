from icet import StructureContainer, VariableBasisClusterSpace

class VariableBasisStructureContainer(StructureContainer):
    cluster_space_type = VariableBasisClusterSpace

    def __init__(self, cluster_space: VariableBasisClusterSpace):
        super().__init__(cluster_space=cluster_space)