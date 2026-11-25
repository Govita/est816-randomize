# EST-816 randomization list

This program builds the allocation sequence for the EST-816 trial: 60 participants, 40 EST-816 and 20 placebo. It uses permuted blocks. A block of six, with four active and two placebo assignments, is one allowed block. The program also uses blocks of nine and mixes the two sizes, so the sequence cannot be predicted from a single known block length.

The realized block-size schedule, the random seed and the treatment assigned to each kit are written only to restricted files. The command used at enrollment prints the next kit number and does not print the treatment or any later kit.

No stratification is the default. The primary analysis adjusts for baseline biological age, and sex and chronological age are sensitivity covariates. With 20 placebo participants, splitting the list by sex and age would leave very small groups. Stratification also cannot keep an exact 40/20 total unless each stratum size is fixed before enrollment. Use strata only when those sizes are known in advance.

## Requirements

Python 3.11 or newer. The program uses the Python standard library only. The generation record stores the Python version and a hash of this program. Reproduce a sequence with that same Python version.

## Create the list

Run this once, from a computer that is not used for recruitment or outcome assessment. Choose a seed and store it only with the restricted files. Do not send it by email.

```bash
python3 est816_randomize.py generate \
  --output restricted \
  --prepared-by "Name of generator" \
  --seed 123456789
```

The seed in this example is not the trial seed. Choose a different one and keep it with the restricted files. Omit `--seed` to let the program draw one. It is printed once and stored in the restricted record. The command checks that the file contains 40 EST-816 and 20 placebo allocations before it writes anything.

If sex-specific lists are required, fix the sizes first. Each size must add up to 60, and each must itself split into the allowed blocks:

```bash
python3 est816_randomize.py generate \
  --output restricted \
  --prepared-by "Name of generator" \
  --seed 123456789 \
  --stratum female=36 \
  --stratum male=24
```

That example gives 24 and 12 for women and 16 and 8 for men.

## Files

| File | Who may see it | Contents |
| --- | --- | --- |
| `RESTRICTED_generation_record.json` | Generator and independent checker | Seed, Python version, program hash, block schedule, approval |
| `RESTRICTED_randomization_master.csv` | Allocation-list holders | Randomization number, kit, treatment and block |
| `RESTRICTED_kit_packing_list.csv` | People who pack the kits | Kit number and treatment, with no block column |
| `site_dispensing_instructions.txt` | Recruiting staff | How to ask for the next kit; no sequence and no treatments |
| `enrollment_log.csv` | Allocation-list holder | Created as kits are released; no treatment column |
| `unblinding_log.csv` | Allocation-list holder | Created only for an emergency disclosure |

Pack the kits from the packing list, then keep that list away from the recruiting staff. Kits with the same number must look the same, including labelling, as required by the protocol.

## Check and approve

A second person, not the person who generated the list, repeats the checks. Approval is recorded only when the totals, the file hash and a fresh regeneration from the seed all agree.

```bash
python3 est816_randomize.py verify \
  --output restricted \
  --approve-by "Name of independent checker"
```

## Release kits

After consent and eligibility are confirmed, the list holder runs one command and tells the recruiter only the printed kit number.

```bash
python3 est816_randomize.py assign-next \
  --output restricted \
  --participant-id GZ-0001
```

For a stratified list, add `--stratum female` or `--stratum male`. The holder can print the number released, but not the upcoming kits:

```bash
python3 est816_randomize.py status --output restricted
```

## Emergency unblinding

This reveals one kit that has already been assigned. It does not reveal a future kit. Follow the protocol's emergency-unblinding rules before using it.

```bash
python3 est816_randomize.py unblind \
  --output restricted \
  --randomization-number R001 \
  --requested-by "Study physician" \
  --reason "Urgent clinical care"
```

## Tests

```bash
python3 -m unittest discover -s tests
```