# === FIX: pin JAX to a version ColabDesign supports ===
# ColabDesign calls jax.lib.xla_bridge, removed in modern JAX (this image got 0.11.1).
# Try candidate versions until BOTH xla_bridge exists AND a binder model preps.
import subprocess, time
t0 = time.time()

PROBE = (
    "import jax, jax.lib; "
    "assert hasattr(jax.lib,'xla_bridge'), 'no xla_bridge'; "
    "from colabdesign import mk_afdesign_model; "
    "m = mk_afdesign_model(protocol='binder', data_dir='/kaggle/working'); "
    "m.prep_inputs(pdb_filename='/kaggle/working/targets/EGFR_A_domainI.pdb', "
    "chain='A', binder_len=60, hotspot='20,21,22,25,28,50,51'); "
    "print('PROBE OK', jax.__version__, jax.devices())"
)

def probe():
    r = subprocess.run(["python", "-c", PROBE], capture_output=True, text=True,
                       cwd="/kaggle/working")
    return r.returncode == 0, (r.stdout + r.stderr)[-500:]

ok, msg = probe()
print("current:", "OK" if ok else "fails ->", msg.strip()[-200:], flush=True)

if not ok:
    for ver in ["0.4.38", "0.4.35", "0.5.3", "0.4.30"]:
        print("\n--- trying jax %s ---" % ver, flush=True)
        r = subprocess.run(
            "pip install -q 'jax[cuda12]==%s' 'jaxlib==%s'" % (ver, ver),
            shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            print("  install failed:", (r.stdout + r.stderr)[-300:], flush=True)
            continue
        ok, msg = probe()
        print("  probe:", "OK" if ok else "fail", msg.strip()[-250:], flush=True)
        if ok:
            print("\nPINNED JAX %s" % ver, flush=True)
            break

print("\n" + ("GENERATOR READY" if ok else "GENERATOR STILL BROKEN"))
print("jax pin step done in %.1f min" % ((time.time() - t0) / 60))
