#!/usr/bin/env python3
"""Plant known geodesy, LOD and rendering bugs one at a time and see what notices.

Each mutant is applied to a scratch copy of the repo (the working tree is
never touched) and built in its own build dir. Then the unit tests and the
engine checks (terra --check) run against it.

  python3 tools/mutants.py                 all mutants
  python3 tools/mutants.py e2-sign     one mutant
  python3 tools/mutants.py --list

Extra configure flags can be passed in TERRA_CMAKE_ARGS.
Exit code 0 if every mutant was caught, 1 otherwise.
"""
import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MUTANTS = {
    "e2-sign": (
        "src/ellipsoid.cpp",
        "(n * (1.0 - e2_) + g.h) * s",
        "(n * (1.0 + e2_) + g.h) * s",
        "ECEF z uses 1 + e^2 instead of 1 - e^2",
    ),
    "geodetic-one-iteration": (
        "src/ellipsoid.cpp",
        "for (int k = 0; k < 30; ++k) {",
        "for (int k = 0; k < 1; ++k) {",
        "ECEF to geodetic stops after one iteration",
    ),
    "vincenty-c-term": (
        "src/ellipsoid.cpp",
        "const double C = f / 16 * cos2_alpha",
        "const double C = f / 8 * cos2_alpha",
        "Vincenty's C coefficient doubled",
    ),
    "child-quadrant-mixup": (
        "src/tiles.cpp",
        "2 * k.y + static_cast<std::uint32_t>(quadrant >> 1)",
        "2 * k.y + static_cast<std::uint32_t>(quadrant & 1)",
        "child row taken from the column bit",
    ),
    "no-balance": (
        "src/select.cpp",
        "for (bool changed = opt.balance; changed;) {",
        "for (bool changed = false; changed;) {",
        "quadtree never balanced",
    ),
    "horizon-plane-ignored": (
        "src/select.cpp",
        "  if (dot(sc, c_hat) + rs >= 1.0 / d) return false;  // reaches the camera side of the horizon plane\n",
        "",
        "horizon test only checks the cone, so spheres in front of the globe get culled",
    ),
    "bounding-sphere-halved": (
        "src/select.cpp",
        "  s.radius *= 1.02;",
        "  s.radius *= 0.5;",
        "tile bounding spheres too small, so visible tiles get culled",
    ),
    "shallow-skirts": (
        "src/mesh.cpp",
        "return 0.5 * tile_span(k.level) * e.a() / grid + 50.0;",
        "return 0.01 * tile_span(k.level) * e.a() / grid;",
        "skirts too short to cover T-junction cracks",
    ),
    "fill-rule-inclusive": (
        "src/render.cpp",
        "if ((w0 == 0 && !tl0) || (w1 == 0 && !tl1) || (w2 == 0 && !tl2)) continue;",
        "if (false && ((w0 == 0 && !tl0) || (w1 == 0 && !tl1) || (w2 == 0 && !tl2))) continue;",
        "pixels on shared edges drawn by both triangles",
    ),
    "no-near-clip": (
        "src/render.cpp",
        "const bool ain = a.p.z >= near, bin = b.p.z >= near;",
        "const bool ain = true, bin = true;",
        "triangles behind the camera are projected anyway",
    ),
    "no-depth-test": (
        "src/render.cpp",
        "      if (iz <= img.inv_depth[idx]) continue;\n",
        "",
        "later triangles always overwrite nearer ones",
    ),
    "focal-from-width": (
        "src/camera.cpp",
        "return 0.5 * height / std::tan(0.5 * fov_y);",
        "return 0.5 * width / std::tan(0.5 * fov_y);",
        "focal length computed from the width with a vertical field of view",
    ),
}



def run(cmd, cwd, timeout=900):
    try:
        return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        return subprocess.CompletedProcess(cmd, -999, e.stdout or "", "timeout")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    args = ap.parse_args()
    if args.list:
        for n, m in MUTANTS.items():
            print(f"{n:26} {m[3]}")
        return 0
    names = args.names or list(MUTANTS)
    for n in names:
        if n not in MUTANTS:
            print(f"unknown mutant {n}", file=sys.stderr)
            return 2

    extra = shlex.split(os.environ.get("TERRA_CMAKE_ARGS", ""))
    ignore = shutil.ignore_patterns("build*", ".git", "__pycache__")
    survivors = []
    print(f"{'mutant':26} {'unit tests':>11} {'checks':>10}")
    with tempfile.TemporaryDirectory(prefix="terra-mutants-") as tmp:
        for name in names:
            path, old, new, what = MUTANTS[name]
            src = os.path.join(tmp, name)
            shutil.copytree(ROOT, src, ignore=ignore)
            target = os.path.join(src, path)
            text = open(target).read()
            if text.count(old) != 1:
                print(f"{name}: anchor not found exactly once in {path}; update tools/mutants.py")
                return 2
            open(target, "w").write(text.replace(old, new))
            build = os.path.join(src, "build")
            c = run(["cmake", "-S", src, "-B", build, "-DCMAKE_BUILD_TYPE=Release", *extra], src)
            b = run(["cmake", "--build", build, "-j", str(args.jobs)], src)
            if c.returncode or b.returncode:
                print(f"{name}: does not compile\n{c.stderr[-1500:]}{b.stdout[-2500:]}{b.stderr[-1500:]}")
                return 2
            t = run([os.path.join(build, "terra_tests"), "--gtest_brief=1"], src)
            failed_tests = sorted(set(re.findall(r"\[  FAILED  \] (\S+) \(", t.stdout)))
            unit = "pass" if t.returncode == 0 else (f"{len(failed_tests)} fail" if failed_tests else "crash")
            k = run([os.path.join(build, "terra"), "--check"], src)
            failed_checks = [l.split()[0] for l in k.stdout.splitlines() if re.match(r"^\S+\s+FAIL\b", l)]
            checks = "pass" if k.returncode == 0 else (f"{len(failed_checks)} fail" if failed_checks else "crash")
            caught = t.returncode != 0 or k.returncode != 0
            print(f"{name:26} {unit:>11} {checks:>10}   {'' if caught else 'SURVIVED  '}{what}")
            if failed_checks:
                print(f"{'':26} checks: {', '.join(failed_checks)}")
            if failed_tests:
                print(f"{'':26} tests: {', '.join(failed_tests[:4])}{' ...' if len(failed_tests) > 4 else ''}")
            if not caught:
                survivors.append(name)
    if survivors:
        print(f"\n{len(survivors)} mutant(s) survived: {', '.join(survivors)}")
        return 1
    print("\nall mutants caught")
    return 0


if __name__ == "__main__":
    sys.exit(main())
