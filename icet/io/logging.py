import logging
import sys

# This is the root logger of icet
logger = logging.getLogger('icet')

# Set the format which is simpler for the console
FORMAT = '%(name)s: %(levelname)s  %(message)s'
formatter = logging.Formatter(FORMAT)

# Will process all levels of WARNING or higher
logger.setLevel(logging.WARNING)

# If you know what you are doing you may set this to True
logger.propagate = False

# The logger will collect events from children and the default
# behaviour is to print it directly to stdout
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(formatter)
logger.addHandler(handler)


def set_config(filename=None, level=None):

    # If a filename is provided a logfile will be created
    if filename is not None:
        fh = logging.FileHandler(filename)
        logger.addHandler(fh)

    if level is not None:
        logger.setLevel(level)
