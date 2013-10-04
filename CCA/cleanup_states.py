with open('states.txt') as ofile:
	with open('states_clean.txt', 'w') as wfile:
		for line in ofile:
			states = line.strip().split(' ')

			for state in states:
				wfile.write(state)

			wfile.write('\n')


