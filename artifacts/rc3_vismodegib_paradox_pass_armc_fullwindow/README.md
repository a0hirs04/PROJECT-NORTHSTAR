# RC3 Vismodegib Paradox Pass

Verdict: `RC3 = PASS`  
`C3.5 = PASS, but provisional`

## What Failed Originally

The first RC3 run reproduced the intended SHH-paradox first half:

- `C3.1` passed: `B > A`
- `C3.2` passed: `radius(B) > radius(A)`
- `C3.3` passed: `periECM(B) < periECM(A)`

But the second half failed:

- `C3.4` failed: `C` was not better than `A`
- `C3.5` failed against the available RC2 comparator
- `C3.6` failed: `C < A < B` was broken

## Schedule Audit Finding

Arm C was inheriting the RC2 drug window unchanged:

- `drug_start_time = 20160`
- `drug_end_time = 40320`
- `max_time = 80640`

That left Arm C off drug during weeks 6-8 after the barrier had already been thinned. The mechanism audit showed Arm C had better drug access and stronger kill than the frozen RC2 gold baseline during active dosing, so the dominant failure mode was schedule mismatch with post-withdrawal rebound.

## Rescue

Arm C was rerun with a schedule-only change:

- `drug_start_time = 20160`
- `drug_end_time = 80640`

All other biology and model code were left unchanged.

Old vs new Arm C day-56 median live tumor: `1111 -> 127`

## Final RC3 Result

- `C3.1`: PASS
- `C3.2`: PASS
- `C3.3`: PASS
- `C3.4`: PASS
- `C3.5`: PASS, provisional
- `C3.6`: PASS

Median endpoint tumor counts:

- `A = 842`
- `B = 2772`
- `C = 127`

Final rank order:

- `C < A < B`

## Reproducibility

This tag includes the tracked launcher and evaluator used to produce and assess this pass:

- [launch_rc3.py](/home/a0hirs04/PROJECT-NORTHSTAR/launch_rc3.py)
- [evaluate_rc3.py](/home/a0hirs04/PROJECT-NORTHSTAR/evaluate_rc3.py)

The passing Arm C full-window configs are saved under [`arm_c_fullwindow/`](./arm_c_fullwindow), and the mixed-arm evaluation report is saved as [`rc3_armc_fullwindow_summary.txt`](./rc3_armc_fullwindow_summary.txt).

## Provisional Note

`C3.5` is still provisional because the clean RC2 gold comparator currently covers only seed `42` in:

`/home/a0hirs04/PROJECT-NORTHSTAR/build/rc2_full_seed42`

Next cleanup task:

Create a clean 5-seed RC2 gold benchmark set for seeds `42-46` with unchanged gold-master biology so `C3.5` can be finalized.
