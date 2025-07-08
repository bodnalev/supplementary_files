from subprocess import run

# Proposition 3.1
run('python3 validator.py proposition_3_1.pickle', shell=True, check=True) # 0m
run('python3 verifier.py proposition_3_1.pickle', shell=True, check=True)  # 90h

# Proposition 3.2
run('python3 validator.py proposition_3_2.pickle', shell=True, check=True) # 0m
run('python3 verifier.py proposition_3_2.pickle', shell=True, check=True)  # 90h

# Proposition 3.3
run('python3 validator.py proposition_3_3.pickle', shell=True, check=True) # 0m
run('python3 verifier.py proposition_3_3.pickle', shell=True, check=True)  # 3h

# Proposition 3.4
run('python3 validator.py proposition_3_4.pickle', shell=True, check=True) # 0m
run('python3 verifier.py proposition_3_4.pickle', shell=True, check=True)  # 0m
