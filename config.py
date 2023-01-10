host = ''

port = ''

user = ''

passwd = ''

db = ''

db_cosmic = 'COSMIC'

db_chip = ''

charset=''

# Fraction of samples of the cohort to set the threshold used to discard mutations present
# in too many samples as presumable germinal variant or artifact
TOTAL_SAMPLES_FRACTION = 0.1
# VAF threshold to emit a mutacion
VAF_THRESHOLD = 0.01
# MAF threshold to consider a mutation as a populaton variant (germinal)
MAX_MAF_THRESHOLD = 0.01
#  Threshold of samples with mutation a vaf compatible with germinal variants
MAX_NUM_SAMPLES_GERMINAL_VAF = 10
# Max occurrences consider common in COSMIC (used in heuristic filter)
COMMON_COSMIC_THRESHOLD = 3
# VAF threshold for the heuristic filter
HEURISTIC_VAF_THRESHOLD = 0.02
# Fraction of samples of the cohort for the heuristic filter
HEURISTIC_FRACTION_SAMPLES = 0.02
# Max occurrences in COSMIC CMC for the heuristic filter
HEURISTIC_COSMIC_OCURRENCES_THRESHOLD = 3
