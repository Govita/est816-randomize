import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import est816_randomize as randomize


PROGRAM = ROOT / "est816_randomize.py"


def run_program(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(PROGRAM), *args],
        check=False,
        capture_output=True,
        text=True,
    )


class SequenceTests(unittest.TestCase):
    def test_default_trial_totals_and_blocks(self) -> None:
        rows, schedule = randomize.build_sequence(60, 40, 20, seed=20260215, strata=[("none", 60)])
        self.assertEqual([row["treatment"] for row in rows].count(randomize.ACTIVE), 40)
        self.assertEqual([row["treatment"] for row in rows].count(randomize.PLACEBO), 20)
        sizes = schedule[0]["block_sizes"]
        self.assertEqual(sum(sizes), 60)
        self.assertTrue(all(size in (6, 9) for size in sizes))
        self.assertEqual(len({row["randomization_number"] for row in rows}), 60)
        self.assertEqual(rows[0]["randomization_number"], "R001")
        self.assertEqual(rows[0]["kit_number"], "EST816-001")

    def test_same_seed_is_reproduced_and_other_seeds_differ(self) -> None:
        first, _ = randomize.build_sequence(60, 40, 20, 99, [("none", 60)])
        second, _ = randomize.build_sequence(60, 40, 20, 99, [("none", 60)])
        self.assertEqual(first, second)
        others = {
            tuple(row["treatment"] for row in randomize.build_sequence(60, 40, 20, seed, [("none", 60)])[0])
            for seed in range(20)
        }
        self.assertGreater(len(others), 1)

    def test_only_the_last_slot_in_a_block_is_forced(self) -> None:
        rows, _ = randomize.build_sequence(60, 40, 20, 7, [("none", 60)])
        blocks: dict[str, list[str]] = {}
        for row in rows:
            blocks.setdefault(row["block_id"], []).append(row["treatment"])
        for block in blocks.values():
            for index in range(len(block) - 1):
                self.assertEqual(set(block[index:]), {randomize.ACTIVE, randomize.PLACEBO})

    def test_many_seeds_pass_the_independent_check(self) -> None:
        for seed in range(30):
            rows, _ = randomize.build_sequence(60, 40, 20, seed, [("none", 60)])
            randomize.check_rows(rows, 60, 40, 20, [("none", 60)])

    def test_stratum_order_does_not_change_the_sequence(self) -> None:
        forward, _ = randomize.build_sequence(60, 40, 20, 5, [("female", 36), ("male", 24)])
        reverse, _ = randomize.build_sequence(60, 40, 20, 5, [("male", 24), ("female", 36)])
        self.assertEqual(forward, reverse)
        female = [row for row in forward if row["stratum"] == "female"]
        male = [row for row in forward if row["stratum"] == "male"]
        self.assertEqual([row["treatment"] for row in female].count(randomize.ACTIVE), 24)
        self.assertEqual([row["treatment"] for row in male].count(randomize.ACTIVE), 16)

    def test_invalid_strata_are_rejected(self) -> None:
        with self.assertRaises(randomize.RandomizationError):
            randomize.parse_strata(["female=36", "male=20"], 60)
        with self.assertRaises(randomize.RandomizationError):
            randomize.parse_strata(["female=3", "male=57"], 60)

    def test_tampered_treatment_fails_verification(self) -> None:
        rows, _ = randomize.build_sequence(60, 40, 20, 3, [("none", 60)])
        rows[0]["treatment"] = randomize.PLACEBO
        with self.assertRaises(randomize.RandomizationError):
            randomize.check_rows(rows, 60, 40, 20, [("none", 60)])


class CommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name) / "restricted"

    def tearDown(self) -> None:
        self.temp.cleanup()

    def generate(self, *extra: str) -> subprocess.CompletedProcess[str]:
        return run_program(
            "generate",
            "--output",
            str(self.directory),
            "--prepared-by",
            "Test Generator",
            "--seed",
            "12345",
            *extra,
        )

    def test_generate_verify_and_release_one_kit_at_a_time(self) -> None:
        generated = self.generate()
        self.assertEqual(generated.returncode, 0, generated.stderr)
        self.assertIn("pending independent verification", generated.stdout)
        self.assertNotIn("block_sizes", generated.stdout)

        packing = (self.directory / randomize.PACKING_NAME).read_text(encoding="utf-8")
        self.assertNotIn("block_id", packing)
        instructions = (self.directory / randomize.INSTRUCTIONS_NAME).read_text(encoding="utf-8")
        self.assertNotIn(randomize.ACTIVE, instructions)
        self.assertNotIn(randomize.PLACEBO, instructions)

        early = run_program(
            "assign-next",
            "--output",
            str(self.directory),
            "--participant-id",
            "P001",
        )
        self.assertEqual(early.returncode, 1)
        self.assertIn("not approved", early.stderr)

        verified = run_program(
            "verify",
            "--output",
            str(self.directory),
            "--approve-by",
            "Independent Reviewer",
        )
        self.assertEqual(verified.returncode, 0, verified.stderr)
        record = json.loads((self.directory / randomize.RECORD_NAME).read_text(encoding="utf-8"))
        self.assertEqual(record["seed"], 12345)
        self.assertEqual(record["python_version"], sys.version.split()[0])
        self.assertEqual(record["code_sha256"], randomize.file_sha256(PROGRAM))
        self.assertEqual(record["approval"]["approved_by"], "Independent Reviewer")
        self.assertEqual(len(record["block_schedule_restricted"][0]["block_sizes"]), record["block_count"])

        first = run_program(
            "assign-next",
            "--output",
            str(self.directory),
            "--participant-id",
            "P001",
        )
        second = run_program(
            "assign-next",
            "--output",
            str(self.directory),
            "--participant-id",
            "P002",
        )
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("Dispense kit: EST816-001", first.stdout)
        self.assertIn("Dispense kit: EST816-002", second.stdout)
        self.assertNotIn(randomize.ACTIVE, first.stdout)
        self.assertNotIn(randomize.PLACEBO, first.stdout)
        duplicate = run_program(
            "assign-next",
            "--output",
            str(self.directory),
            "--participant-id",
            "P001",
        )
        self.assertEqual(duplicate.returncode, 1)

        status = run_program("status", "--output", str(self.directory))
        self.assertIn("Assigned: 2 / 60", status.stdout)
        self.assertNotIn("EST816-003", status.stdout)

        hidden = run_program(
            "unblind",
            "--output",
            str(self.directory),
            "--randomization-number",
            "R003",
            "--requested-by",
            "Study physician",
            "--reason",
            "future kit",
        )
        self.assertEqual(hidden.returncode, 1)
        opened = run_program(
            "unblind",
            "--output",
            str(self.directory),
            "--randomization-number",
            "R001",
            "--requested-by",
            "Study physician",
            "--reason",
            "urgent clinical care",
        )
        self.assertEqual(opened.returncode, 0, opened.stderr)
        self.assertIn("Treatment:", opened.stdout)
        with (self.directory / randomize.ENROLLMENT_NAME).open(encoding="utf-8") as handle:
            log = list(csv.DictReader(handle))
        self.assertEqual([row["participant_id"] for row in log], ["P001", "P002"])
        self.assertNotIn("treatment", log[0])

    def test_existing_output_is_not_replaced(self) -> None:
        self.assertEqual(self.generate().returncode, 0)
        second = self.generate()
        self.assertEqual(second.returncode, 1)
        self.assertIn("Refusing to overwrite", second.stderr)

    def test_stratified_release_follows_each_stratum(self) -> None:
        generated = self.generate("--stratum", "male=24", "--stratum", "female=36")
        self.assertEqual(generated.returncode, 0, generated.stderr)
        approved = run_program(
            "verify",
            "--output",
            str(self.directory),
            "--approve-by",
            "Independent Reviewer",
        )
        self.assertEqual(approved.returncode, 0, approved.stderr)
        female = run_program(
            "assign-next",
            "--output",
            str(self.directory),
            "--participant-id",
            "F001",
            "--stratum",
            "female",
        )
        male = run_program(
            "assign-next",
            "--output",
            str(self.directory),
            "--participant-id",
            "M001",
            "--stratum",
            "male",
        )
        self.assertIn("EST816-female-001", female.stdout)
        self.assertIn("EST816-male-001", male.stdout)
        missing = run_program(
            "assign-next",
            "--output",
            str(self.directory),
            "--participant-id",
            "X001",
        )
        self.assertEqual(missing.returncode, 1)


if __name__ == "__main__":
    unittest.main()
