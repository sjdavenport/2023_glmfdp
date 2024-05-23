import sys
sys.path.insert(0, '/home/mind/sdavenpo/2023_glmfdp/Simulations/JER_control/EJER_runs')
from EJER_runs import run_nsubj_dim_review
run_nsubj_dim_review(41, 100, simtype = 1, niters = 5000, df=5)
run_nsubj_dim_review(41, 100, simtype = 0, niters = 5000, df=5)