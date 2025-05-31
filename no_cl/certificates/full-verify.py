from subprocess import run

# Proposition 3.1
run('python3 validator.py prop_3_1_c5k4.pickle', shell=True, check=True)  # 0m
run('python3 verifier.py prop_3_1_c5k4.pickle', shell=True, check=True)   # 19m
run('python3 validator.py prop_3_1_c7.pickle', shell=True, check=True)    # 0m
run('python3 verifier.py prop_3_1_c7.pickle', shell=True, check=True)     # 10m

# Proposition 3.2
run('python3 validator.py prop_3_2_c5k4.pickle', shell=True, check=True)  # 25m
run('python3 verifier.py prop_3_2_c5k4.pickle', shell=True, check=True)   # 19m
run('python3 validator.py prop_3_2_c7.pickle', shell=True, check=True)    # 10m
run('python3 verifier.py prop_3_2_c7.pickle', shell=True, check=True)     # 10m

# Proposition 3.3
run('python3 validator.py prop_3_3_c5k4.pickle', shell=True, check=True)  # 28m
run('python3 verifier.py prop_3_3_c5k4.pickle', shell=True, check=True)   # 57h
run('python3 validator.py prop_3_3_c7.pickle', shell=True, check=True)    # 16m
run('python3 verifier.py prop_3_3_c7.pickle', shell=True, check=True)     # 26h
