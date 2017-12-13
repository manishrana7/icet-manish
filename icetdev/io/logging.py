import logging
import sys


# This is the root logger of icet
icetdev_logger = logging.getLogger('icetdev')

# Will process all levels of INFO or higher
icetdev_logger.setLevel(logging.INFO)

# If you know what you are doing you may set this to True
icetdev_logger.propagate = False

# The HiPhive logger will collect events from childs and the default behaviour
# is to print it directly to stdout
ch = logging.StreamHandler(sys.stdout)
icetdev_logger.addHandler(ch)


# continuous = True

def set_config(filename=None, level=None, continuous=None):

    # If a filename is provided a logfile will be created
    if filename is not None:
        fh = logging.FileHandler(filename)
        icetdev_logger.addHandler(fh)

    if level is not None:
        icetdev_logger.setLevel(level)

#    if continuous is not None:
#        sys.modules[__name__].continuous = continuous
