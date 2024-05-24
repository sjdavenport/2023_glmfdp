import sys
sys.path.insert(0, '/home/mind/sdavenpo/2023_glmfdp/Simulations/JER_control/')
from EJER_runs import run_nsubj_dim_review
run_nsubj_dim_review(107, 100, simtype = 1, niters = 1000, df=5)
run_nsubj_dim_review(107, 100, simtype = -1, niters = 1000, df=5)