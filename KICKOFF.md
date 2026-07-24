# Kickoff — read this, then start

You are picking up mzLibR with no memory of the session that planned it. Everything you need is on
disk.

## 1. Read, in this order

1. **`E:\CodeReview\mzLibRust\STATUS.md`** — what was just built and what it cost. Start here.
2. **`E:\CodeReview\mzLibR\PLAN.md`** — your brief. Long on purpose; it is meant to be the only thing
   you need. §6 (every trap, with numbers) and §8 (the documentation doctrine) are the parts two
   previous builds actually paid for.

## 2. Check what moved before you write anything

Two pull requests were open when the plan was written, and `PLAN.md` describes their subjects as
*current* behaviour. If they merged, some of it is now wrong:

```
gh pr view 1114 --repo smith-chem-wisc/mzLib --json state,title
gh pr view 9 --repo smith-chem-wisc/pyMzLib --json state,title
```

- **mzLib #1114** removes `y` ions from ETD/ECD. If merged, §6.3's first row is history, not a
  live trap.
- **pyMzLib #9** fixes `Peptide.intensity()` returning `None`. If merged, that back-port is done.

Also worth a glance: mzLib issues #1110, #1111, #1112, #1113 and pyMzLib #8.

## 3. Stage a bridge

```
E:\CodeReview\pyMzLib\code\pyMzLib\pkg\python\src\pymzlib\_dotnet\win-x64\mzlib-bridge.exe
```

Set `MZLIB_BRIDGE` to it, or copy it where the package will look. `E:\CodeReview\mzLibRust\scripts\stage-bridge.ps1`
is the reference for doing this properly, including the version probe.

## 4. Build it

Follow `PLAN.md` §13 milestones: M0 transport → M1 PRIDE → M2 peptidoform → M3 FlashLFQ.

**Work locally. Do not create the public repo until M0–M3 are green** — the first public commit
should be a working package, which is how mzLibRust was done.

`E:\CodeReview\mzLibRust` is the closest model for everything: module layout, the pure-function test
seam, doc voice, licence files, findings and test-parity docs.

## 5. Two standing rules

- **An mzLib bug becomes an mzLib issue.** Verify it first — one finding last session was
  investigated and correctly dropped as intended behaviour.
- **A documentation lesson is back-ported to every binding.** The traps live in mzLib's behaviour,
  not in one language.

## 6. Before you claim it works

Offline tests must pass with **no network and no bridge staged**. Live tests must **skip, not fail**,
when a service is down. See `PLAN.md` §9.
