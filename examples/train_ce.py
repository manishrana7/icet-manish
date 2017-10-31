import clusterspace

prim = 1

clusterspace = clusterspace(prim, cutoffs, components)
sc = StructureContainer(clusterspace, supercells)
optimizer = Optimizer('LASSO',
                      sc.get_fit_data(property=..., structures=...),
                      validation_data=sc.get_fit_data(property=..., structures=...))

ce_quatuplets = cluster_expansion(optimizer, fix_ce=[ce_pairs, ce_triplets])
ce_quatuplets.train(training_fraction=0.5, validation_fraction=0.9)

model = clusters.get_model(fj.parameters)


def fix_ce(self, cluster_expansion):
    """ 
    Changes self so that the eci's of the input cluster_expansion is kept fixed
    during training.
    """

    for i, 

