from pybedtools import BedTool
import sys
bam = BedTool(sys.argv[1])
bf = bam.filter(lambda x: (x.end - x.start >= 200000))
bf.saveas(sys.argv[1]+'_filtered.bed')
