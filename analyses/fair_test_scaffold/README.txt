Fair test of the effect of large scaffold groups on training (same test under conditions A, B, and C).
Prerequisites: copy features.py, modeling.py, and modeling_tox21.py from the pipeline folder into work/; datasets in
work/data_ic50, work/data_layers, and work/data_tox21; results from the old regime in work/out/old_<group>_<endpoint>.json
(from the resultados_brutos package of the previous reassessment).
Order: (1) cd work; python fair_test.py <ic50|layers|tox21> <endpoint> <random_seeds...>  (fair_tasks.txt lists the calls;
the sweep used 4 random seeds on ChEMBL and 3 on Tox21: the argument for the fifth random seed was cut from the queue);
(2) python fair_test_analysis.py  -> out/fair_long.csv, out/fair_summary_paired.csv, out/fair_consistency_random_seeds.csv.
Environment: edcqsar-ml (rdkit, scikit-learn, xgboost per the pipeline's requirements.txt).
