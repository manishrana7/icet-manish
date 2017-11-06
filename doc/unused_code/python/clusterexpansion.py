class ClusterExpansion(object):

    def __init__(self, structureset, subset=None, exclude=None,
                 method='Split-Bregman', **kwargs):
        self.method = 'Split-Bregman'

    def get_ecis(self):
        return self.ecis

    def get_errors(self):
        return 0

    def get_true_and_fitted_values(self):
        return 0

    def get_training_rms_error(self):
        return 0

    def get_rms_error(self, structureset):
        return 0
