import sys
import numpy
def get_data():
  for line in sys.stdin:
    yield float(line)


print numpy.median(list(get_data()))
