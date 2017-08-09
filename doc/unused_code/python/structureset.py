class StructureSet(object):

    def __init__(self, clusters):

        self.clusters = clusters
        self.structures = []
        self.ids = []
        self.cluster_vectors = []
        self.fit_properties = []


    def load_structures_from_database(self, database,
    								  fit_property_field=None,
    								  id_type='running',
    								  cluster_vector_field=None,
    								  map_to_prototype=True,
    								  save_cluster_vector_in_database=False):
	    from ase.db import connect
         
        try:
            db = connect(database)
        except:
             raise Exception()
        
        for entry in db.select():
        	self.structures.append(entry.toatoms())

        	# Structure id
            if id_type == 'running':
        		if self.ids == []:
        			structure_id = 0
        		else:
        			structure_id = max(self.ids) + 1
        	else:
        		try:
        			structure_id = entry.get(id_type)
        		except:
        			raise Exception()
        	self.ids.append(structure_id)


        	# Cluster vector
            if cluster_vector_field != None:
                try:
                    cv = entry.get(cluster_vector_field)
                except:
                    raise Exception()
                assert len(cv) == len(self.clusters)
            else:
                cv = clusters.calculate_cluster_vector(entry.to_atoms())
            self.cvs.append(cv)
            if save_cluster_vector_in_databse:
               assert cluster_vector_field == None
            	    db.update(entry.id, cluster_vector=cv)

             # Fit property
             if fit_property_field != None:
             	try:
             		fit_property = entry.get(fit_property_field)
             	except:
             		raise Exception()
             else:
             	fit_property = None
             self.fit_properties.append(fit_property)



    def add_structure(self, structure, fit_property=None, structure_id=None,
                      cluster_vector=[], map_to_prototype=True):
        self.structures.append(structure)

        # Structure id
        if structure_id == None:
        	if self.ids == []:
        		structure_id = 0
        	else:
        		structure_id = max(self.ids) + 1
        self.ids.append(structure_id)

        # Cluster vector
        if cluster_vector != []:
            assert len(cluster_vector) == len(self.clusters)
        elif calculate_cluster_vector:
            cluster_vector = clusters.calculate_cluster_vector(structure)
        self.cluster_vectors.append(cluster_vector)
        
        self.fit_properties.append(fit_property)


    def save(self, filename, save_clusters=True):
        return 0

    def fit(self, subset=[], exclude=[], method='Split-Bregman', **kwargs):
    	# check fit proporty existence etc
        return ClusterExpansion(self, subset=subset, exclude=exclude,
                                method=method, **kwargs)

    def get_loo_error(self, method='Split-Bregman', **kwargs):
        return 0

    def get_lpo_error(self, p=20, repetitions=10, 
                      method='Split-Bregman', **kwargs):
        return 0

    def get_number_of_structures(self):
        return len(self.structures)


    def copy(self):
        return copy.copy(self)

    def copy_with_filter(self, cluster_filter):
        new_cluster = self.clusters.copy_with_filter(cluster_filter)
        new_structureset = self.copy()
        new_structureset.clusters = new_cluster
        # modify cluster vectors etc
        return new_structureset

    def __len__(self):
        return len(self.structures)

    def __str__(self):
        return('Structure set based on {}'.format(str(self.clusters)))
