# Running a major refactor across several sessions

Read this only when setting up work of that shape; the standing rules
that apply to ordinary work are in `CLAUDE.md`.


For work too large for one session, split by **file ownership**, not by task: several
worker sessions edit disjoint file sets concurrently, one coordinator session holds the
decisions. 

**Roles.** Workers own files and do the engineering. The coordinator owns no files. It
dispatches, arbitrates ownership, rules on cross-cutting questions, submits cluster
jobs, and records decisions. A worker never edits another's file, never runs a long job
(it writes the qsub), and never commits.

**Model choice.** The coordinator needs the strongest model available: its job is
*judgment at decision points*, and the expensive mistakes are there. But its token burn
is dominated by routing and log-reading, not by judgment. So keep its I/O disciplined
(below) rather than downgrading it. Workers also need a strong model -- every real bug in
the rebuild (bracket-H hydrogen counts, an mmCIF quoting bug, a shared-area credit that
contradicted its own reference state, an error handler that did not cover the write) was
found by a worker verifying on real data in past projects, not by the coordinator reading 
a report.

**The coordinator must not hold the project in its context -- it holds it on disk.**
Three files, all under `~/docking/scratch/claude_plan_notes/`:
- `README_sessions.md` -- one paragraph: who owns what, read first in every session.
- `decision_records.md` -- one record per settled decision: problem, decision, why,
  evidence, rejected alternatives, cost accepted, **when to revisit**. Write these
  *before* dispatch. Decisions made mid-flight get relitigated; pre-recorded ones do not.
- `inbox.md` -- append-only, one section per session. A worker messages another by
  appending under its heading. This removes the user from the relay loop entirely.
Plus an ownership table (function-level where two sessions share a file). With those on
disk a fresh coordinator is current in minutes, so **restart the coordinator when its
context gets heavy** instead of nursing it.

**Starting the workers.** The user opens the worker sessions manually and tells the
coordinator their names (`ListAgents` lists peers; `SendMessage` reaches them). Do not
have the coordinator spawn them as subagents. First prompt to each worker: read the 
three files above and your plan file, then reply with the files you own, your first 
three tasks, and any conflict with the decision records -- and wait. Coordinator should 
check that reply before saying go; it catches bad ownership assumptions and genuinely 
caught issues for past projects.

**Keeping the bug rate down.** Most defects found during the 2026-09 rebuild were
pre-existing ones the refactor surfaced (a latent `makedirs` race, a stale contact cutoff,
a permutation group that ignored the key) or new-code slips the author caught in its own
verification pass. Those are cheap and expected. The **expensive** ones all came from the
coordinator, and they are the ones to design against:
- **An underspecified spec.** "Credit each shared point 1/k to its k occluders" did not say
  whether the ligand counts in k. A worker implemented it faithfully, it went through a
  calibration sweep and a recall test, and only then was it found to contradict the
  leave-one-out reference it sat next to. **A spec for anything numeric must state the
  invariant it has to satisfy, not just the formula** -- the closure identity that caught it
  (`sum over partners == area occluded by >=1 partner and no ligand atom`) belonged in the
  spec, where it would have failed on day one.
- **A constant fixed off a thin stratum.** A threshold was set from a bound resting on 11
  pairs; at 5x the sample it moved by 5x and two full calibration runs were wasted.
  **Before quoting any bound, state n for the smallest cell it depends on**, and treat a
  cell under ~30 as directional only.
- **Re-running a job against code that just changed.** Check the script's mtime against the
  submit time before trusting a result, and never resubmit a failed job until the fix that
  caused the failure has actually landed.
- **Reasoning from a variable's NAME instead of its definition.** Three coordinator specs
  failed this way in one night, each caught by the implementer in one arithmetic check:
  banding a size-dependent quantity and reading composition off it (it measured residue
  volume); dividing by `cg_buriable` believing it was the fragment's surface (it is the
  area partners took, so the ratio's dominant axis was partner count); and asserting "0.25
  is about one area quantum" from oxygen alone (the quantum scales with radius, so
  phosphorus and the heavy halogens do not clear it). **Before specifying a comparison, ask
  what the quantity is mechanically a function of, and write its formula into the spec. If
  you cannot write the formula, you do not know it well enough to specify a statistic on
  it.** For any "X is approximately Y" claim, evaluate Y at the extremes of its domain, not
  at the typical case.
- **The cheapest fix for all of the above is a spec echo.** Before building, the implementer
  restates the spec as a formula plus the invariant it must satisfy, and sends it back. It
  is one message; it would have caught all three errors above before a single calibration
  run. Cheaper still: prefer specs the implementer can falsify on day one (state the
  closure identity, the control that must come out flat, the case that must fail).

**Verify continuously; the user is usually not watching.** A coordinator running unattended
must re-check its own state rather than trust its last report: confirm a job is running the
file version you think it is, re-read a file's mtime before acting on a stale belief about
it, and treat "a worker said it is done" as a claim to spot-check on disk, not a fact. Every
worker claim in the rebuild that was checked on disk held up -- but two stopped without
doing the work they had acknowledged, and only a mtime check caught it. Prefer one cheap
verifying command per exchange over a batch of trust.

**Rules that earned their place.**
- Hand-offs get *pinned* (exact names, dtypes, shapes) at dispatch, not negotiated
  later. An unpinned hand-off makes one session wait on another for no reason.
- Each worker runs only its own test subset; the full suite is a coordinator gate on a
  settled tree, because three sessions mid-edit make it red for uninteresting reasons.
- Workers report verdicts and measurements, never raw output. The coordinator greps for
  the verdict line and never reads a raw log.
- The user makes every commit. Commit per session so a clobber is attributable, and
  commit a baseline *before* any worker starts.
- When a worker pushes back on a ruling with evidence, it is usually right: it has the
  file open. Ask for the missing measurement rather than overruling.

