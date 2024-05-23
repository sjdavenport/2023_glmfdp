import sys
sys.path.insert(0, '/home/mind/sdavenpo/2023_glmfdp/Simulations/JER_control/')
from EJER_runs import run_nsubj_dim_review
run_nsubj_dim_review(26, 100, simtype = 1, niters = 5000, df=3)
run_nsubj_dim_review(26, 100, simtype = -1, niters = 5000, df=3)