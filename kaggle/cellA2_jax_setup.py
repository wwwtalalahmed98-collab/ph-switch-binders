# === CHALLENGE 1 : JAX GENERATOR SETUP (replaces RFdiffusion, which needs dgl) ===
# dgl has no build that imports on this image's Python 3.13, so RFdiffusion cannot run
# here. ColabDesign/AfDesign is JAX-based, needs no dgl, and is the engine BindCraft
# itself is built on.
import os, subprocess, time, urllib.request
t0 = time.time()

def sh(cmd, check=False):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print("FAILED:", cmd)
        print((r.stdout + r.stderr)[-1200:])
        if check:
            raise SystemExit("setup failed")
    return r

import sys
print("python", sys.version.split()[0], flush=True)

# ---- JAX with CUDA, then ColabDesign ----
# Pin JAX: ColabDesign needs jax.lib.xla_bridge, removed in modern JAX. 0.4.38 is the
# newest version verified working with ColabDesign on this image's Python 3.13.
PROBE_V = "python -c \"import jax, jax.lib; assert hasattr(jax.lib,'xla_bridge'); print(jax.__version__)\""
r = subprocess.run(PROBE_V, shell=True, capture_output=True, text=True)
if r.returncode != 0 or not r.stdout.strip().startswith("0.4.38"):
    print("installing jax 0.4.38 (ColabDesign-compatible) ...", flush=True)
    sh("pip install -q 'jax[cuda12]==0.4.38' 'jaxlib==0.4.38'")
r = subprocess.run("python -c \"import jax; print('JAX', jax.__version__, jax.devices())\"",
                   shell=True, capture_output=True, text=True)
print(r.stdout.strip() or (r.stderr[-600:]), flush=True)

sh("pip install -q git+https://github.com/sokrypton/ColabDesign.git")

# ---- AlphaFold2 parameters (~3.5 GB) ----
os.makedirs("/kaggle/working/params", exist_ok=True)
if not os.path.exists("/kaggle/working/params/params_model_1_ptm.npz"):
    print("downloading AF2 params (~3.5 GB) ...", flush=True)
    sh("wget -q -O /kaggle/working/af2.tar "
       "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar")
    sh("tar -xf /kaggle/working/af2.tar -C /kaggle/working/params")
    sh("rm -f /kaggle/working/af2.tar")
print("AF2 params:",
      len([f for f in os.listdir('/kaggle/working/params') if f.endswith('.npz')]),
      "files", flush=True)

# ---- verify the generator actually imports and can see the target ----
chk = subprocess.run(
    "cd /kaggle/working && python -c \""
    "from colabdesign import mk_afdesign_model; "
    "m = mk_afdesign_model(protocol='binder', data_dir='/kaggle/working'); "
    "print('COLABDESIGN IMPORT OK')\"",
    shell=True, capture_output=True, text=True)
print(chk.stdout[-600:])
if chk.returncode != 0:
    print("COLABDESIGN IMPORT FAILED:")
    print(chk.stderr[-1800:])
else:
    print("GENERATOR READY")

print("JAX SETUP DONE in %.1f min" % ((time.time() - t0) / 60))
