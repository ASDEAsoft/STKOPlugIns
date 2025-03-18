
# beam_column_elements elasticBeamColumn
# Geometric transformation command
geomTransf Linear 1 0.0 1.0 0.0
element elasticBeamColumn 1 1 2 558.0 493634.48 214623.687 14784.12 4994.272 5481.24 1
# Geometric transformation command
geomTransf Linear 2 0.0 1.0 0.0
element elasticBeamColumn 2 2 3 558.0 4936.3448 2146.23687 14784.12 4994.272 5481.24 2
# Geometric transformation command
geomTransf Linear 3 0.0 1.0 0.0
element elasticBeamColumn 3 3 4 558.0 493634.48 214623.687 14784.12 4994.272 5481.24 3

# zero_length_elements zeroLength
element zeroLength 4 5 7 -mat 7 7 7 7 8 9 -dir 1 2 3 4 5 6 -orient 0.0 0.0 1.0 1.0 0.0 0.0

# zero_length_elements zeroLengthSection
element zeroLengthSection 5 7 1 3000001 -orient 0.0 0.0 1.0 -1.0 0.0 0.0

# zero_length_elements zeroLength
element zeroLength 6 6 8 -mat 7 7 7 7 10 11 -dir 1 2 3 4 5 6 -orient 0.0 0.0 1.0 1.0 0.0 0.0

# zero_length_elements zeroLengthSection
element zeroLengthSection 7 8 4 3000002 -orient 0.0 0.0 -1.0 -1.0 0.0 0.0
