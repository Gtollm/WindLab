"""
Shape optimizer for drag coefficient extremization.

Uses the CMA-ES algorithm (Covariance Matrix Adaptation - Evolution Strategy)
to search for a 3-D body of revolution that minimizes or maximizes Cd at a
given Reynolds number, subject to:
  - fixed length along the flow direction (L_FLOW)
  - fixed volume (V_TARGET)
  - shape is far from domain boundaries (enforced by eval-stl padding)

Usage:
  # minimize Cd (most streamlined shape)
  python optimize.py --re 50 --mode min

  # maximize Cd (bluffest shape)
  python optimize.py --re 50 --mode max

  # override volume / length
  python optimize.py --re 100 --mode min --volume 2e-5 --length 0.08

Requirements:
  pip install cma numpy scipy
  cargo build --release   (windlab binary must be up-to-date)
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import cma
import numpy as np


# ── progress display ──────────────────────────────────────────────────────────

class ProgressTracker:
    """Prints a compact progress line that overwrites itself each update."""

    def __init__(self, max_iter: int, pop_size: int, mode: str):
        self.max_iter  = max_iter
        self.pop_size  = pop_size
        self.mode      = mode
        self.total_est = max_iter * pop_size   # rough upper bound on evals
        self.iter      = 0
        self.eval      = 0
        self.best_cd   = float("inf") if mode == "min" else float("-inf")
        self.t_start   = time.time()
        self.t_iter    = time.time()   # time of last iter start
        self.secs_per_eval: float | None = None

    def record_eval(self, cd: float) -> bool:
        """Call after each windlab run. Returns True if new best."""
        self.eval += 1
        is_best = (
            (self.mode == "min" and cd < self.best_cd) or
            (self.mode == "max" and np.isfinite(cd) and cd > self.best_cd)
        )
        if is_best:
            self.best_cd = cd
        return is_best

    def record_iter(self, sigma: float) -> None:
        """Call after each CMA-ES iteration (after es.tell)."""
        now = time.time()
        elapsed_iter = now - self.t_iter
        self.secs_per_eval = elapsed_iter / max(self.pop_size, 1)
        self.t_iter = now
        self.iter  += 1
        self._print_iter(sigma)

    def _fmt_time(self, seconds: float) -> str:
        seconds = int(seconds)
        if seconds < 60:
            return f"{seconds}s"
        if seconds < 3600:
            return f"{seconds//60}m{seconds%60:02d}s"
        return f"{seconds//3600}h{(seconds%3600)//60:02d}m"

    def _print_iter(self, sigma: float) -> None:
        elapsed = time.time() - self.t_start
        pct     = self.iter / self.max_iter * 100.0

        # ETA from average time-per-eval
        if self.secs_per_eval is not None:
            remaining_evals = (self.max_iter - self.iter) * self.pop_size
            eta = self._fmt_time(remaining_evals * self.secs_per_eval)
        else:
            eta = "?"

        cd_str = (
            f"{self.best_cd:.4f}" if np.isfinite(self.best_cd) else "  n/a"
        )
        bar_w = 20
        filled = int(bar_w * self.iter / self.max_iter)
        bar = "█" * filled + "░" * (bar_w - filled)

        line = (
            f"\r  [{bar}] iter {self.iter:3d}/{self.max_iter}"
            f"  eval {self.eval:4d}"
            f"  best Cd={cd_str}"
            f"  σ={sigma:.4f}"
            f"  {self._fmt_time(elapsed)} / ETA {eta}"
            f"   "          # padding to overwrite previous longer line
        )
        print(line, end="", flush=True)

# ── project root relative to this script ─────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
_exe = "windlab.exe" if sys.platform == "win32" else "windlab"

def _pick_binary() -> Path:
    """Return the release binary. Debug is forbidden (15-20x slower)."""
    import subprocess as _sp
    candidate = ROOT / "target" / "release" / _exe
    if candidate.exists():
        try:
            out = _sp.run([str(candidate), "--help"],
                          capture_output=True, text=True, timeout=5).stdout
            if "eval-stl" in out:
                return candidate
        except Exception:
            pass
    raise RuntimeError(
        f"Release binary not found or missing 'eval-stl'. Run: cargo build --release"
    )

WINDLAB = _pick_binary()

# ── import shape generator ────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
import shape_gen as sg


# ── evaluator ────────────────────────────────────────────────────────────────

def evaluate(
    params: np.ndarray,
    *,
    re: float,
    tau: float,
    cpd: int,
    steps: int,
    V_target: float,
    L: float,
    stl_dir: str,
) -> float:
    """
    Generate STL → call windlab eval-stl → parse Cd.
    Returns +inf on any failure so CMA-ES treats it as a very bad candidate.
    """
    fd, stl_path = tempfile.mkstemp(suffix=".stl", dir=stl_dir)
    os.close(fd)
    try:
        sg.save_stl(params, stl_path, V_target=V_target, L=L)

        result = subprocess.run(
            [
                str(WINDLAB), "eval-stl",
                "--stl",   stl_path,
                "--re",    str(re),
                "--tau",   str(tau),
                "--cpd",   str(cpd),
                "--steps", str(steps),
                "--no-progress",
            ],
            capture_output=True,
            text=True,
            timeout=900,
        )

        for line in result.stdout.splitlines():
            line = line.strip()
            if line.startswith("Cd="):
                cd = float(line.split("=", 1)[1])
                if np.isfinite(cd) and cd > 0.0:
                    return cd
                break

        print(f"  [WARN] no valid Cd — stdout: {result.stdout!r}  "
              f"stderr: {result.stderr[-300:]!r}", file=sys.stderr)
        return float("inf")

    except subprocess.TimeoutExpired:
        print("  [WARN] windlab timed out", file=sys.stderr)
        return float("inf")
    except Exception as exc:
        print(f"  [ERR] {exc}", file=sys.stderr)
        return float("inf")
    finally:
        try:
            os.unlink(stl_path)
        except OSError:
            pass


# ── optimizer ─────────────────────────────────────────────────────────────────

def run(
    mode: str,          # "min" or "max"
    re: float,
    V_target: float,
    L: float,
    cpd: int,
    steps: int,
    n_ctrl: int,
    max_iter: int,
    seed: int,
) -> tuple[np.ndarray, float]:

    tau = 0.5 + 4.5 / re   # keeps u_inlet = 0.05·D/D = 0.05 at any Re

    sign = -1.0 if mode == "max" else 1.0   # CMA-ES always minimizes

    x0 = sg.default_params(n_ctrl=n_ctrl, V_target=V_target, L=L)
    lb, ub = sg.param_bounds(n_ctrl=n_ctrl, V_target=V_target, L=L)
    sigma0 = float((ub[0] - lb[0]) * 0.25)

    stl_dir = tempfile.mkdtemp(prefix="windlab_opt_")

    opts = cma.CMAOptions()
    opts["seed"]     = seed
    opts["maxiter"]  = max_iter
    opts["bounds"]   = [lb.tolist(), ub.tolist()]
    opts["verbose"]  = -9          # suppress cma internal output
    opts["tolx"]     = 1e-5
    opts["tolfun"]   = 1e-4

    es = cma.CMAEvolutionStrategy(x0.tolist(), sigma0, opts)

    pop_size    = es.popsize
    best_params = x0.copy()
    tracker     = ProgressTracker(max_iter=max_iter, pop_size=pop_size, mode=mode)

    print(f"\n{'─'*60}")
    print(f"  {'MINIMIZE' if mode == 'min' else 'MAXIMIZE'} Cd")
    print(f"  Re={re}  tau={tau:.3f}  V={V_target:.2e} m³  L={L:.3f} m")
    print(f"  params={len(x0)}  pop={pop_size}  cpd={cpd}  steps/eval={steps}")
    print(f"{'─'*60}")

    while not es.stop():
        solutions = es.ask()
        fitnesses = []

        for e_idx, params in enumerate(solutions):
            print(f"\r  iter {tracker.iter+1}/{max_iter}  eval {tracker.eval+1}"
                  f"/{tracker.eval + len(solutions) - e_idx}  running...   ",
                  end="", flush=True)
            cd = evaluate(
                np.asarray(params),
                re=re, tau=tau, cpd=cpd, steps=steps,
                V_target=V_target, L=L, stl_dir=stl_dir,
            )
            fitnesses.append(sign * cd)

            is_best = tracker.record_eval(cd)
            if is_best:
                best_params = np.asarray(params).copy()

        es.tell(solutions, fitnesses)
        tracker.record_iter(es.sigma)   # prints progress line

    # newline after the last \r progress line
    print()

    best_cd = tracker.best_cd
    print(f"\n{'─'*60}")
    print(f"  RESULT ({mode}): Cd = {best_cd:.4f}")
    print(f"  Total evaluations: {tracker.eval}  iterations: {tracker.iter}")
    print(f"{'─'*60}\n")

    # Clean up temp dir
    try:
        import shutil
        shutil.rmtree(stl_dir, ignore_errors=True)
    except Exception:
        pass

    return best_params, best_cd


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="WindLab shape optimizer")
    parser.add_argument("--mode",    choices=["min", "max"], default="min",
                        help="Minimize or maximize Cd (default: min)")
    parser.add_argument("--re",      type=float, default=50.0,
                        help="Reynolds number (default: 50)")
    parser.add_argument("--volume",  type=float, default=sg.V_TARGET,
                        help=f"Target volume in m³ (default: {sg.V_TARGET:.2e})")
    parser.add_argument("--length",  type=float, default=sg.L_FLOW,
                        help=f"Fixed shape length in m (default: {sg.L_FLOW})")
    parser.add_argument("--cpd",     type=int,   default=30,
                        help="Cells per cross-stream diameter (default: 30)")
    parser.add_argument("--steps",   type=int,   default=3000,
                        help="LBM steps per evaluation (default: 3000)")
    parser.add_argument("--n-ctrl",  type=int,   default=sg.N_CTRL,
                        help=f"Interior control points (default: {sg.N_CTRL})")
    parser.add_argument("--max-iter",type=int,   default=60,
                        help="CMA-ES max iterations (default: 60)")
    parser.add_argument("--seed",    type=int,   default=42)
    parser.add_argument("--out-stl", type=str,   default=None,
                        help="Save best shape to this STL path")
    args = parser.parse_args()

    if not WINDLAB.exists():
        print(f"ERROR: windlab binary not found at {WINDLAB}", file=sys.stderr)
        print("Run:  cargo build --release", file=sys.stderr)
        sys.exit(1)

    best_params, best_cd = run(
        mode=args.mode,
        re=args.re,
        V_target=args.volume,
        L=args.length,
        cpd=args.cpd,
        steps=args.steps,
        n_ctrl=args.n_ctrl,
        max_iter=args.max_iter,
        seed=args.seed,
    )

    # Save best shape
    out_path = args.out_stl or f"best_{args.mode}_cd_Re{args.re:.0f}.stl"
    sg.save_stl(best_params, out_path, V_target=args.volume, L=args.length)
    print(f"Best shape saved → {out_path}")
    print(f"Final Cd = {best_cd:.4f}")


if __name__ == "__main__":
    main()
