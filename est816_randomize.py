#!/usr/bin/env python3
"""Permuted-block randomization for EST-816.

Builds a reproducible 2:1 allocation list, checks it, and releases one kit
number at a time. The realized block-size schedule is written only to the
restricted generation record. Recruiting staff must not receive that record,
the master list, or the packing list.

The allowed block sizes live in this program because the retained code has to
specify the algorithm. They are not printed in site files or in the usual
on-screen summary.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import secrets
import stat
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from random import Random

ACTIVE = "EST-816"
PLACEBO = "Placebo"
# Sizes that keep a 2:1 ratio inside every block. A block of 6 is the
# illustration used in the protocol (4 active and 2 placebo). A block of 9
# is 6 active and 3 placebo. The generated schedule mixes these sizes.
BLOCK_SIZES = (6, 9)

MASTER_NAME = "RESTRICTED_randomization_master.csv"
PACKING_NAME = "RESTRICTED_kit_packing_list.csv"
RECORD_NAME = "RESTRICTED_generation_record.json"
INSTRUCTIONS_NAME = "site_dispensing_instructions.txt"
ENROLLMENT_NAME = "enrollment_log.csv"
UNBLIND_NAME = "unblinding_log.csv"

MASTER_FIELDS = [
    "stratum",
    "randomization_number",
    "kit_number",
    "block_id",
    "treatment",
]
PACKING_FIELDS = ["kit_number", "treatment"]
ENROLLMENT_FIELDS = [
    "assigned_at_utc",
    "participant_id",
    "stratum",
    "randomization_number",
    "kit_number",
]
UNBLIND_FIELDS = [
    "unblinded_at_utc",
    "randomization_number",
    "kit_number",
    "treatment",
    "requested_by",
    "reason",
]


class RandomizationError(Exception):
    """Invalid settings or a sequence that failed its checks."""


def program_path() -> Path:
    return Path(__file__).resolve()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def restrict(path: Path) -> None:
    """Limit a file to the account that generated it."""
    os.chmod(path, stat.S_IRUSR | stat.S_IWUSR)


def fisher_yates(items: list[str], rng: Random) -> list[str]:
    """Shuffle with randrange so the result depends only on the recorded seed."""
    shuffled = list(items)
    for index in range(len(shuffled) - 1, 0, -1):
        swap = rng.randrange(index + 1)
        shuffled[index], shuffled[swap] = shuffled[swap], shuffled[index]
    return shuffled


def compositions(total: int, sizes: tuple[int, ...] = BLOCK_SIZES) -> list[tuple[int, ...]]:
    """Every ordered way to sum exactly to total using the allowed block sizes."""
    ordered = tuple(sorted(sizes))
    found: list[tuple[int, ...]] = []

    def walk(remaining: int, prefix: list[int]) -> None:
        if remaining == 0:
            found.append(tuple(prefix))
            return
        for size in ordered:
            if size <= remaining:
                prefix.append(size)
                walk(remaining - size, prefix)
                prefix.pop()

    walk(total, [])
    if not found:
        raise RandomizationError(
            f"{total} participants cannot be split into blocks of {ordered[0]} and {ordered[-1]}."
        )
    return found


def one_block(size: int, rng: Random) -> list[str]:
    """Build one 2:1 block whose last assignment is the only forced slot.

    One active and one placebo assignment are reserved for the final two
    positions, then those two positions are shuffled. Until the final slot,
    both treatments remain possible.
    """
    if size % 3 != 0 or size < 6:
        raise RandomizationError(f"Block size {size} does not keep a 2:1 ratio.")
    placebo = size // 3
    active = size - placebo
    head = fisher_yates([ACTIVE] * (active - 1) + [PLACEBO] * (placebo - 1), rng)
    tail = fisher_yates([ACTIVE, PLACEBO], rng)
    return head + tail


def build_stratum(name: str, total: int, rng: Random, block_offset: int) -> tuple[list[dict[str, str]], list[int]]:
    schedule = list(rng.choice_index(compositions(total)))
    rows: list[dict[str, str]] = []
    position = 1
    for offset, size in enumerate(schedule):
        block_id = str(block_offset + offset + 1)
        for treatment in one_block(size, rng):
            if name == "none":
                number = f"R{position:03d}"
                kit = f"EST816-{position:03d}"
            else:
                number = f"R-{name}-{position:03d}"
                kit = f"EST816-{name}-{position:03d}"
            rows.append(
                {
                    "stratum": name,
                    "randomization_number": number,
                    "kit_number": kit,
                    "block_id": block_id,
                    "treatment": treatment,
                }
            )
            position += 1
    return rows, schedule


class SeededRandom(Random):
    """Random source with an explicit, testable way to pick a schedule."""

    def choice_index(self, options: list[tuple[int, ...]]) -> tuple[int, ...]:
        return options[self.randrange(len(options))]


def parse_strata(values: list[str] | None, total: int) -> list[tuple[str, int]]:
    if not values:
        return [("none", total)]
    parsed: list[tuple[str, int]] = []
    seen: set[str] = set()
    for value in values:
        name, separator, raw_size = value.partition("=")
        name = name.strip()
        if (
            separator != "="
            or not name
            or not raw_size.strip().isdigit()
            or not name.replace("-", "").isalnum()
        ):
            raise RandomizationError(
                f"Stratum {value!r} must look like female=36. Use letters, numbers and a size."
            )
        if name in seen or name == "none":
            raise RandomizationError(f"Stratum name {name!r} is repeated or reserved.")
        seen.add(name)
        parsed.append((name, int(raw_size)))
    sizes = [size for _, size in parsed]
    if sum(sizes) != total:
        raise RandomizationError(
            f"Stratum sizes add up to {sum(sizes)}, but the trial total is {total}."
        )
    for name, size in parsed:
        if size < min(BLOCK_SIZES) or compositions(size) == []:
            raise RandomizationError(
                f"Stratum {name} has {size} participants. Each stratum must be at least "
                f"{min(BLOCK_SIZES)} and divisible into the allowed block sizes."
            )
    return parsed


def build_sequence(
    total: int,
    active: int,
    placebo: int,
    seed: int,
    strata: list[tuple[str, int]],
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    if total != active + placebo:
        raise RandomizationError("Active and placebo counts must add up to the number of participants.")
    if active != 2 * placebo:
        raise RandomizationError("This program builds a 2:1 sequence. Active must be twice placebo.")
    if seed < 0 or seed.bit_length() > 63:
        raise RandomizationError("The seed must be a whole number from 0 through 2^63 - 1.")

    rng = SeededRandom(seed)
    rows: list[dict[str, str]] = []
    schedule: list[dict[str, object]] = []
    block_offset = 0
    # Alphabetical stratum order makes the seed independent of argument order.
    for name, size in sorted(strata, key=lambda item: item[0]):
        stratum_rows, stratum_schedule = build_stratum(name, size, rng, block_offset)
        expected_active = size * active // total
        if sum(row["treatment"] == ACTIVE for row in stratum_rows) != expected_active:
            raise RandomizationError(f"Stratum {name} did not keep the 2:1 ratio.")
        rows.extend(stratum_rows)
        schedule.append({"stratum": name, "participants": size, "block_sizes": stratum_schedule})
        block_offset += len(stratum_schedule)
    check_rows(rows, total, active, placebo, strata)
    return rows, schedule


def check_rows(
    rows: list[dict[str, str]],
    total: int,
    active: int,
    placebo: int,
    strata: list[tuple[str, int]],
) -> None:
    """Independent checks of a written or newly built sequence."""
    if len(rows) != total:
        raise RandomizationError(f"The sequence has {len(rows)} rows, not {total}.")
    treatments = [row["treatment"] for row in rows]
    active_count = treatments.count(ACTIVE)
    placebo_count = treatments.count(PLACEBO)
    if active_count != active or placebo_count != placebo:
        raise RandomizationError(
            f"Totals are {active_count} {ACTIVE} and {placebo_count} {PLACEBO}, "
            f"not {active} and {placebo}."
        )
    numbers = [row["randomization_number"] for row in rows]
    kits = [row["kit_number"] for row in rows]
    if len(set(numbers)) != len(numbers) or len(set(kits)) != len(kits):
        raise RandomizationError("A randomization number or kit number is repeated.")
    if any(row["treatment"] not in (ACTIVE, PLACEBO) for row in rows):
        raise RandomizationError("A row has a treatment other than EST-816 or Placebo.")

    expected = dict(strata)
    seen_sizes: dict[str, int] = {}
    by_block: dict[tuple[str, str], list[str]] = {}
    previous_block: dict[str, str | None] = {name: None for name in expected}
    for row in rows:
        stratum = row["stratum"]
        if stratum not in expected:
            raise RandomizationError(f"Unexpected stratum {stratum!r}.")
        seen_sizes[stratum] = seen_sizes.get(stratum, 0) + 1
        block_key = (stratum, row["block_id"])
        if previous_block[stratum] not in (None, row["block_id"]) and block_key in by_block:
            raise RandomizationError(f"Block {row['block_id']} in {stratum} is not contiguous.")
        previous_block[stratum] = row["block_id"]
        by_block.setdefault(block_key, []).append(row["treatment"])

    for name, size in expected.items():
        if seen_sizes.get(name, 0) != size:
            raise RandomizationError(f"Stratum {name} has {seen_sizes.get(name, 0)} rows, not {size}.")
    for (stratum, block_id), block in by_block.items():
        if len(block) not in BLOCK_SIZES:
            raise RandomizationError(
                f"Block {block_id} in {stratum} has {len(block)} assignments, which is not an allowed size."
            )
        if block.count(PLACEBO) * 2 != block.count(ACTIVE):
            raise RandomizationError(f"Block {block_id} in {stratum} is not a 2:1 block.")
        for index in range(len(block) - 1):
            if set(block[index:]) != {ACTIVE, PLACEBO}:
                raise RandomizationError(
                    f"Block {block_id} in {stratum} forces an assignment before its final slot."
                )


def canonical_master(rows: list[dict[str, str]]) -> str:
    lines = [",".join(MASTER_FIELDS)]
    for row in rows:
        lines.append(",".join(row[field] for field in MASTER_FIELDS))
    return "\n".join(lines) + "\n"


def site_instructions(stratified: bool) -> str:
    stratum_step = (
        "Tell the list holder the participant's stratum when you ask for the next kit.\n"
        if stratified
        else ""
    )
    return (
        "Kit dispensing instructions\n"
        "\n"
        "Kits are identical in appearance. Take them only in the order released by the\n"
        "allocation-list holder, after written consent and confirmation of eligibility.\n"
        "\n"
        "1. Confirm consent and eligibility.\n"
        f"{stratum_step}"
        "2. Ask the list holder for the next kit. You receive one randomization number\n"
        "   and one kit number.\n"
        "3. Dispense that sealed kit. Record the participant ID, randomization number\n"
        "   and kit number on the study forms.\n"
        "4. Do not open a later kit, do not choose a kit out of order, and do not ask\n"
        "   what later participants will receive.\n"
        "\n"
        "The treatment list, random seed and packing list are not site documents.\n"
        "Emergency unblinding follows the protocol and is logged separately.\n"
    )


def write_outputs(
    directory: Path,
    rows: list[dict[str, str]],
    schedule: list[dict[str, object]],
    seed: int,
    total: int,
    active: int,
    placebo: int,
    strata: list[tuple[str, int]],
    prepared_by: str,
) -> dict[str, object]:
    directory.mkdir(parents=True, exist_ok=True)
    master_text = canonical_master(rows)
    master_path = directory / MASTER_NAME
    master_path.write_text(master_text, encoding="utf-8")
    packing_rows = [
        {"kit_number": row["kit_number"], "treatment": row["treatment"]}
        for row in sorted(rows, key=lambda row: row["kit_number"])
    ]
    write_csv(directory / PACKING_NAME, PACKING_FIELDS, packing_rows)
    instructions = site_instructions(any(name != "none" for name, _ in strata))
    (directory / INSTRUCTIONS_NAME).write_text(instructions, encoding="utf-8")

    record: dict[str, object] = {
        "study": "EST-816",
        "algorithm": "permuted-block randomization, 2:1, variable block size",
        "rng": "Python random.Random with an integer seed; Fisher-Yates uses randrange",
        "python_version": sys.version.split()[0],
        "python_implementation": sys.implementation.name,
        "code_file": program_path().name,
        "code_sha256": file_sha256(program_path()),
        "generated_at_utc": utc_now(),
        "prepared_by": prepared_by,
        "seed": seed,
        "participants": total,
        "est816": active,
        "placebo": placebo,
        "stratification": "none" if strata == [("none", total)] else "fixed stratum sizes",
        "strata": [{"name": name, "participants": size} for name, size in strata],
        "block_count": sum(len(item["block_sizes"]) for item in schedule),
        "block_schedule_restricted": schedule,
        "master_sha256": text_sha256(master_text),
        "verification": {
            "checked": True,
            "est816": active,
            "placebo": placebo,
            "participants": total,
        },
        "approval": {
            "status": "pending independent verification",
            "approved_by": None,
            "approved_at_utc": None,
        },
    }
    record_path = directory / RECORD_NAME
    record_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    for name in (MASTER_NAME, PACKING_NAME, RECORD_NAME):
        restrict(directory / name)
    return record


def load_record(directory: Path) -> dict[str, object]:
    path = directory / RECORD_NAME
    if not path.exists():
        raise RandomizationError(f"No generation record in {directory}.")
    return json.loads(path.read_text(encoding="utf-8"))


def strata_from_record(record: dict[str, object]) -> list[tuple[str, int]]:
    return [(item["name"], int(item["participants"])) for item in record["strata"]]  # type: ignore[index]


def verify_directory(directory: Path) -> dict[str, object]:
    record = load_record(directory)
    rows = read_csv(directory / MASTER_NAME)
    strata = strata_from_record(record)
    check_rows(
        rows,
        int(record["participants"]),
        int(record["est816"]),
        int(record["placebo"]),
        strata,
    )
    master_hash = text_sha256((directory / MASTER_NAME).read_text(encoding="utf-8"))
    if master_hash != record["master_sha256"]:
        raise RandomizationError("The master file hash does not match the generation record.")
    packing = read_csv(directory / PACKING_NAME)
    master_pairs = sorted((row["kit_number"], row["treatment"]) for row in rows)
    packing_pairs = sorted((row["kit_number"], row["treatment"]) for row in packing)
    if master_pairs != packing_pairs:
        raise RandomizationError("The packing list does not match the master list.")
    if "block_id" in (directory / PACKING_NAME).read_text(encoding="utf-8"):
        raise RandomizationError("The packing list contains block information.")
    instructions = (directory / INSTRUCTIONS_NAME).read_text(encoding="utf-8")
    if ACTIVE in instructions or PLACEBO in instructions or "block" in instructions.lower():
        raise RandomizationError("The site instructions reveal treatment or blocking.")
    regenerated, regenerated_schedule = build_sequence(
        int(record["participants"]),
        int(record["est816"]),
        int(record["placebo"]),
        int(record["seed"]),
        strata,
    )
    if canonical_master(regenerated) != canonical_master(rows):
        raise RandomizationError("Regenerating the seed does not reproduce the approved master list.")
    if regenerated_schedule != record["block_schedule_restricted"]:
        raise RandomizationError("Regenerating the seed does not reproduce the block schedule.")
    return record


def approve(directory: Path, approved_by: str) -> None:
    record = verify_directory(directory)
    approval = record["approval"]
    if not isinstance(approval, dict):
        raise RandomizationError("The generation record has no approval section.")
    if approval["status"] == "approved" and approval["approved_by"] != approved_by:
        raise RandomizationError(
            f"Already approved by {approval['approved_by']}. Refusing to replace that name."
        )
    approval["status"] = "approved"
    approval["approved_by"] = approved_by
    approval["approved_at_utc"] = utc_now()
    path = directory / RECORD_NAME
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    restrict(path)


def read_enrollment(directory: Path) -> list[dict[str, str]]:
    path = directory / ENROLLMENT_NAME
    if not path.exists():
        return []
    return read_csv(path)


@contextmanager
def assignment_lock(directory: Path):
    """One kit can be released at a time, even if the command is started twice."""
    directory.mkdir(parents=True, exist_ok=True)
    handle = (directory / ".assignment.lock").open("a", encoding="utf-8")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def require_approved(directory: Path) -> dict[str, object]:
    record = verify_directory(directory)
    approval = record["approval"]
    if not isinstance(approval, dict) or approval["status"] != "approved":
        raise RandomizationError(
            "The sequence is not approved. A second person must run verify --approve-by first."
        )
    return record


def assign_next(directory: Path, participant_id: str, stratum: str | None) -> dict[str, str]:
    with assignment_lock(directory):
        return _assign_next(directory, participant_id, stratum)


def _assign_next(directory: Path, participant_id: str, stratum: str | None) -> dict[str, str]:
    if not participant_id.strip() or "," in participant_id or "\n" in participant_id:
        raise RandomizationError("The participant ID must be a single non-empty value without commas.")
    record = require_approved(directory)
    rows = read_csv(directory / MASTER_NAME)
    log = read_enrollment(directory)
    if any(entry["participant_id"] == participant_id for entry in log):
        raise RandomizationError(f"Participant {participant_id} already has a kit.")
    strata = [item[0] for item in strata_from_record(record)]
    if strata == ["none"]:
        if stratum not in (None, "none"):
            raise RandomizationError("This sequence has no strata. Omit --stratum.")
        selected = "none"
    else:
        if stratum is None or stratum not in strata:
            raise RandomizationError(
                "This sequence is stratified. Pass --stratum with one of: " + ", ".join(strata)
            )
        selected = stratum
    used = {entry["randomization_number"] for entry in log}
    chosen = next(
        (
            row
            for row in rows
            if row["stratum"] == selected and row["randomization_number"] not in used
        ),
        None,
    )
    if chosen is None:
        raise RandomizationError(f"No kits remain for stratum {selected}.")
    entry = {
        "assigned_at_utc": utc_now(),
        "participant_id": participant_id,
        "stratum": selected,
        "randomization_number": chosen["randomization_number"],
        "kit_number": chosen["kit_number"],
    }
    log.append(entry)
    write_csv(directory / ENROLLMENT_NAME, ENROLLMENT_FIELDS, log)
    restrict(directory / ENROLLMENT_NAME)
    return entry


def unblind(directory: Path, randomization_number: str, requested_by: str, reason: str) -> dict[str, str]:
    with assignment_lock(directory):
        return _unblind(directory, randomization_number, requested_by, reason)


def _unblind(directory: Path, randomization_number: str, requested_by: str, reason: str) -> dict[str, str]:
    if not reason.strip() or not requested_by.strip():
        raise RandomizationError("Emergency unblinding needs the requester's name and a reason.")
    require_approved(directory)
    log = read_enrollment(directory)
    if not any(entry["randomization_number"] == randomization_number for entry in log):
        raise RandomizationError(
            "That randomization number has not been assigned. Refusing to reveal a future kit."
        )
    rows = read_csv(directory / MASTER_NAME)
    matched = next(row for row in rows if row["randomization_number"] == randomization_number)
    entry = {
        "unblinded_at_utc": utc_now(),
        "randomization_number": matched["randomization_number"],
        "kit_number": matched["kit_number"],
        "treatment": matched["treatment"],
        "requested_by": requested_by,
        "reason": reason,
    }
    path = directory / UNBLIND_NAME
    existing = read_csv(path) if path.exists() else []
    existing.append(entry)
    write_csv(path, UNBLIND_FIELDS, existing)
    restrict(path)
    return entry


def allocation_status(directory: Path) -> dict[str, object]:
    record = verify_directory(directory)
    log = read_enrollment(directory)
    approval = record["approval"]
    return {
        "approval": approval["status"] if isinstance(approval, dict) else "unknown",
        "assigned": len(log),
        "participants": int(record["participants"]),
        "stratification": record["stratification"],
    }


def command_generate(args: argparse.Namespace) -> int:
    directory = Path(args.output)
    if (directory / RECORD_NAME).exists() and not args.force:
        raise RandomizationError(
            f"{directory} already contains a sequence. Refusing to overwrite it without --force."
        )
    seed = secrets.randbits(63) if args.seed is None else args.seed
    strata = parse_strata(args.stratum, args.participants)
    rows, schedule = build_sequence(args.participants, args.active, args.placebo, seed, strata)
    record = write_outputs(
        directory,
        rows,
        schedule,
        seed,
        args.participants,
        args.active,
        args.placebo,
        strata,
        args.prepared_by,
    )
    print(f"Restricted files written to: {directory.resolve()}")
    print(f"Stratification: {record['stratification']}")
    print(f"Allocations checked: {args.active} {ACTIVE}, {args.placebo} {PLACEBO}")
    print(f"Blocks: {record['block_count']} (the size schedule is only in {RECORD_NAME})")
    print(f"Seed: {seed}")
    print(f"Python: {record['python_version']}")
    print(f"Code SHA-256: {record['code_sha256']}")
    print(f"Master SHA-256: {record['master_sha256']}")
    print("Status: pending independent verification")
    print("Do not send the RESTRICTED files, the seed, or this screen output to recruiting staff.")
    return 0


def command_verify(args: argparse.Namespace) -> int:
    directory = Path(args.output)
    record = verify_directory(directory)
    print(f"Participants: {record['participants']}")
    print(f"{ACTIVE}: {record['est816']}")
    print(f"{PLACEBO}: {record['placebo']}")
    print(f"Stratification: {record['stratification']}")
    print("Duplicate numbers: none")
    print("Every block keeps a 2:1 ratio: yes")
    print("Master hash matches the record: yes")
    print("Seed reproduces the master list: yes")
    if args.approve_by:
        approve(directory, args.approve_by)
        print(f"Approved by: {args.approve_by}")
    else:
        approval = record["approval"]
        status = approval["status"] if isinstance(approval, dict) else "unknown"
        print(f"Approval status: {status}")
    return 0


def command_assign(args: argparse.Namespace) -> int:
    entry = assign_next(Path(args.output), args.participant_id, args.stratum)
    print(f"Participant: {entry['participant_id']}")
    print(f"Randomization number: {entry['randomization_number']}")
    print(f"Dispense kit: {entry['kit_number']}")
    return 0


def command_status(args: argparse.Namespace) -> int:
    status = allocation_status(Path(args.output))
    print(f"Approval: {status['approval']}")
    print(f"Assigned: {status['assigned']} / {status['participants']}")
    print(f"Stratification: {status['stratification']}")
    return 0


def command_unblind(args: argparse.Namespace) -> int:
    entry = unblind(Path(args.output), args.randomization_number, args.requested_by, args.reason)
    print(f"Randomization number: {entry['randomization_number']}")
    print(f"Kit: {entry['kit_number']}")
    print(f"Treatment: {entry['treatment']}")
    print("This disclosure was added to the unblinding log.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate and release the EST-816 randomization list. "
            "No stratification is the default and is appropriate for this trial: "
            "the primary analysis adjusts for baseline biological age, and 20 placebo "
            "participants would make several strata very small."
        )
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    generate = subcommands.add_parser("generate", help="create a restricted randomization list")
    generate.add_argument("--output", required=True, help="directory for the restricted files")
    generate.add_argument("--prepared-by", required=True, help="name of the person generating the list")
    generate.add_argument("--seed", type=int, help="integer seed; a random seed is chosen if omitted")
    generate.add_argument("--participants", type=int, default=60)
    generate.add_argument("--active", type=int, default=40)
    generate.add_argument("--placebo", type=int, default=20)
    generate.add_argument(
        "--stratum",
        action="append",
        help="optional fixed stratum, for example female=36; repeat once per stratum",
    )
    generate.add_argument("--force", action="store_true", help="replace an existing sequence")
    generate.set_defaults(func=command_generate)

    verify = subcommands.add_parser("verify", help="check totals, hash and seed reproduction")
    verify.add_argument("--output", required=True)
    verify.add_argument("--approve-by", help="record a second person's approval if the checks pass")
    verify.set_defaults(func=command_verify)

    assign = subcommands.add_parser("assign-next", help="release the next kit number only")
    assign.add_argument("--output", required=True)
    assign.add_argument("--participant-id", required=True)
    assign.add_argument("--stratum", help="required only when the sequence was stratified")
    assign.set_defaults(func=command_assign)

    status = subcommands.add_parser("status", help="show how many kits have been released")
    status.add_argument("--output", required=True)
    status.set_defaults(func=command_status)

    unblind = subcommands.add_parser("unblind", help="reveal one assigned kit for an emergency")
    unblind.add_argument("--output", required=True)
    unblind.add_argument("--randomization-number", required=True)
    unblind.add_argument("--requested-by", required=True)
    unblind.add_argument("--reason", required=True)
    unblind.set_defaults(func=command_unblind)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except RandomizationError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
