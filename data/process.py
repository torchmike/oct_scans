import sys
import numpy
def get_data():
  data = []
  fail_count = 0
  for line in sys.stdin:
    if float(line) != -1:
      data.append(float(line))
    else:
      fail_count +=1
  return [data,fail_count]


data, fail_count = get_data()
median = numpy.median(data)
name = sys.argv[1]
mean = numpy.mean(data)
stdev = numpy.std(data)
print "{name}\t{median}\t{mean}\t{stdev}\t{fail_count}\t{count}".format(name=name, mean=mean, stdev=stdev, median=median, fail_count=fail_count, count = len(data))
