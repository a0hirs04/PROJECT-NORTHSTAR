# RC2 True Biological Pass Baseline

This directory freezes the first clean RC2 true biological pass produced from the `drug_kill_multiplier = 0.01` seed-42 run.

## Winning Run

- Output directory: `/home/a0hirs04/PROJECT-NORTHSTAR/build/rc2_full_seed42/replicate_01_seed42/output`
- Frozen run config: [`winning_config.xml`](/home/a0hirs04/PROJECT-NORTHSTAR/artifacts/rc2_true_biological_pass_km0p01_seed42/winning_config.xml)
- Evaluator report: [`evaluate_rc2_true_biological_pass_km0p01_seed42.txt`](/home/a0hirs04/PROJECT-NORTHSTAR/artifacts/rc2_true_biological_pass_km0p01_seed42/evaluate_rc2_true_biological_pass_km0p01_seed42.txt)

## Key Results

- Final classification: `TRUE BIOLOGICAL PASS`
- Hard score: `5/5`
- Day 14 live tumor: `723`
- Day 28 live tumor: `301`
- Day 42 live tumor: `431`
- Reduction at day 28: `58.37%`
- `peri_ecm` day 14: `0.9825`
- `peri_ecm` day 28: `0.9995`
- `ecm_at_survivors` day 14: `0.5594`
- `ecm_at_survivors` day 28: `0.6577`

## Interpretation

Treatment produced a partial response without eradication, the ECM barrier remained intact through treatment, survivors occupied higher-ECM sanctuary zones, ABCB1 resistance emerged, and the tumor regrew after withdrawal. This is the first clean RC2 true biological pass for this branch using the post-fix ECM pipeline and `drug_kill_multiplier = 0.01`.
