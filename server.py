import os
import uuid
import sys
from flask import Flask, request, jsonify, send_from_directory
import shutil
import subprocess
from werkzeug.utils import secure_filename

app = Flask(__name__)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_DIR = os.path.join(BASE_DIR, "workspace")
STATIC_DIR = os.path.join(BASE_DIR, "static")
os.makedirs(WORKSPACE_DIR, exist_ok=True)

@app.route("/")
def index():
    return send_from_directory(STATIC_DIR, "index.html")

@app.route("/static/<path:filename>")
def static_files(filename):
    return send_from_directory(STATIC_DIR, filename)

@app.route("/api/status")
def status():
    return jsonify({"ok": True})

@app.route("/api/env")
def env():
    return jsonify({"ok": True, "python": sys.version, "executable": sys.executable, "venv": os.path.exists(os.path.join(BASE_DIR, ".venv"))})

@app.route("/api/upload", methods=["POST"])
def upload():
    files = {}
    for key in ["top", "trj", "mapping"]:
        f = request.files.get(key)
        if f:
            files[key] = f
    if not files:
        return jsonify({"ok": False, "error": "no files"}), 400
    rid = str(uuid.uuid4())
    out_dir = os.path.join(WORKSPACE_DIR, rid)
    os.makedirs(out_dir, exist_ok=True)
    saved = {}
    for k, f in files.items():
        name = secure_filename(f.filename)
        path = os.path.join(out_dir, name)
        f.save(path)
        saved[k] = path
    return jsonify({"ok": True, "rid": rid, "saved": saved})

@app.route("/api/votca/check")
def votca_check():
    cmds = ["csg_map", "csg_boltzmann", "csg_stat", "csg_fmatch", "csg_call"]
    avail = {}
    for c in cmds:
        avail[c] = {"system_path": shutil.which(c)}
    conda_avail = {}
    if shutil.which("conda"):
        for c in cmds:
            try:
                res = subprocess.run(["conda", "run", "-n", "cg", "bash", "-lc", f"which {c}"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                conda_avail[c] = res.stdout.strip() if res.returncode == 0 else None
            except Exception:
                conda_avail[c] = None
    return jsonify({"ok": True, "available": avail, "conda_cg": conda_avail})

@app.route("/api/votca/map", methods=["POST"])
def votca_map():
    rid = request.form.get("rid")
    if not rid:
        return jsonify({"ok": False, "error": "missing rid"}), 400
    out_dir = os.path.join(WORKSPACE_DIR, rid)
    if not os.path.isdir(out_dir):
        return jsonify({"ok": False, "error": "rid not found"}), 404
    mapping = None
    top = None
    trj = None
    for name in os.listdir(out_dir):
        p = os.path.join(out_dir, name)
        if name.endswith(".xml"):
            mapping = p
        elif name.endswith(".tpr") or name.endswith(".top") or name.endswith(".gro"):
            top = p
        elif name.endswith(".xtc") or name.endswith(".trr"):
            trj = p
    if not mapping or not top or not trj:
        return jsonify({"ok": False, "error": "missing inputs"}), 400
    use_conda = shutil.which("conda") is not None
    try:
        if use_conda:
            cmd = ["conda", "run", "-n", "cg", "bash", "-lc", f"csg_map --top '{top}' --cg '{mapping}' --trj '{trj}'"]
        else:
            if shutil.which("csg_map") is None:
                return jsonify({"ok": False, "error": "csg_map not available"}), 500
            cmd = ["csg_map", "--top", top, "--cg", mapping, "--trj", trj]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return jsonify({"ok": True, "returncode": res.returncode, "stdout": res.stdout, "stderr": res.stderr})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500

@app.route("/api/workflow/ethylene", methods=["POST"])
def workflow_ethylene():
    tpr = os.path.join(BASE_DIR, "example/ethylene/eth.tpr")
    trr = os.path.join(BASE_DIR, "example/ethylene/eth.trr")
    xtc = os.path.join(BASE_DIR, "example/ethylene/eth.xtc")
    gro = os.path.join(BASE_DIR, "example/ethylene/eth.gro")
    top = os.path.join(BASE_DIR, "example/ethylene/eth.top")
    mapping = os.path.join(BASE_DIR, "example/ethylene/mapping_two_beads.xml")
    logs = []
    if shutil.which("conda"):
        cmd1 = ["conda", "run", "-n", "cg", "bash", "-lc", f"gmx grompp -f '{BASE_DIR}/example/ethylene/min.mdp' -c '{gro}' -p '{top}' -o '{tpr}'"]
        res1 = subprocess.run(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "grompp", "rc": res1.returncode, "out": res1.stdout, "err": res1.stderr})
        cmd2 = ["conda", "run", "-n", "cg", "bash", "-lc", f"gmx mdrun -s '{tpr}' -o '{trr}' -x '{xtc}' -nsteps 10"]
        res2 = subprocess.run(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "mdrun", "rc": res2.returncode, "out": res2.stdout, "err": res2.stderr})
    else:
        return jsonify({"ok": False, "error": "conda cg env required"}), 500
    cg_out = os.path.join(BASE_DIR, "example/ethylene/cg_conf.gro")
    if shutil.which("docker"):
        cmd3 = ["bash", "-lc", f"docker run --rm -v '{BASE_DIR}':'/work' votca/votca:latest csg_map --top /work/example/ethylene/eth.tpr --trj /work/example/ethylene/eth.trr --cg /work/example/ethylene/mapping_two_beads.xml --out /work/example/ethylene/cg_conf.gro"]
        res3 = subprocess.run(cmd3, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_map_docker", "rc": res3.returncode, "out": res3.stdout, "err": res3.stderr})
    else:
        cmd3 = ["conda", "run", "-n", "cg", "bash", "-lc", f"csg_map --top '{tpr}' --trj '{trr}' --cg '{mapping}' --out '{cg_out}'"]
        res3 = subprocess.run(cmd3, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_map_local", "rc": res3.returncode, "out": res3.stdout, "err": res3.stderr})
    exists = os.path.exists(cg_out)
    return jsonify({"ok": exists, "cg": "example/ethylene/cg_conf.gro" if exists else None, "logs": logs})

@app.route("/api/file")
def read_file():
    p = request.args.get("path")
    if not p:
        return jsonify({"ok": False, "error": "missing path"}), 400
    abspath = os.path.abspath(os.path.join(BASE_DIR, p)) if not os.path.isabs(p) else p
    if not abspath.startswith(BASE_DIR):
        return jsonify({"ok": False, "error": "forbidden"}), 403
    if not os.path.isfile(abspath):
        return jsonify({"ok": False, "error": "not found"}), 404
    with open(abspath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    return jsonify({"ok": True, "path": abspath.replace(BASE_DIR+"/", ""), "content": content})

@app.route("/api/workflow/ff", methods=["POST"])
def workflow_ff():
    tpr = os.path.join(BASE_DIR, "example/ethylene/eth.tpr")
    trr = os.path.join(BASE_DIR, "example/ethylene/eth.trr")
    mapping = os.path.join(BASE_DIR, "example/ethylene/mapping_two_beads.xml")
    out_dir = os.path.join(BASE_DIR, "example/ethylene")
    logs = []
    if not (os.path.isfile(tpr) and os.path.isfile(trr)):
        return jsonify({"ok": False, "error": "missing tpr/trr, please run ethylene workflow first"}), 400
    if shutil.which("docker"):
        cmd = ["bash", "-lc", f"docker run --rm -v '{BASE_DIR}':'/work' votca/votca:latest csg_boltzmann --top /work/example/ethylene/eth.tpr --trj /work/example/ethylene/eth.trr --cg /work/example/ethylene/mapping_two_beads.xml"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_boltzmann", "rc": res.returncode, "out": res.stdout, "err": res.stderr})
    else:
        cmd = ["conda", "run", "-n", "cg", "bash", "-lc", f"csg_boltzmann --top '{tpr}' --trj '{trr}' --cg '{mapping}'"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_boltzmann", "rc": res.returncode, "out": res.stdout, "err": res.stderr})
    produced = []
    for name in os.listdir(out_dir):
        if name.endswith('.pot') or name.endswith('.force') or name.endswith('.dist'):
            produced.append(os.path.join('example/ethylene', name))
    return jsonify({"ok": len(produced) > 0, "files": produced, "logs": logs})

@app.route("/api/ff/list")
def ff_list():
    out_dir = os.path.join(BASE_DIR, "example/ethylene")
    files = []
    for name in os.listdir(out_dir):
        if name.endswith('.pot') or name.endswith('.force') or name.endswith('.dist'):
            files.append(os.path.join('example/ethylene', name))
    return jsonify({"ok": True, "files": sorted(files)})

@app.route("/api/ff/parse", methods=["POST"])
def ff_parse():
    path = request.form.get("path") or (request.json.get("path") if request.is_json else None)
    if not path:
        return jsonify({"ok": False, "error": "missing path"}), 400
    abspath = os.path.abspath(os.path.join(BASE_DIR, path)) if not os.path.isabs(path) else path
    if not abspath.startswith(BASE_DIR):
        return jsonify({"ok": False, "error": "forbidden"}), 403
    if not os.path.isfile(abspath):
        return jsonify({"ok": False, "error": "not found"}), 404
    xs = []
    ys = []
    with open(abspath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(";") or s.startswith("#"):
                continue
            parts = s.split()
            nums = []
            for p in parts:
                try:
                    nums.append(float(p))
                except Exception:
                    nums = []
                    break
            if len(nums) >= 2:
                xs.append(nums[0])
                ys.append(nums[1])
    return jsonify({"ok": True, "x": xs, "y": ys, "path": path})

@app.route("/api/rdf/run", methods=["POST"])
def rdf_run():
    tpr = os.path.join(BASE_DIR, "example/ethylene/eth.tpr")
    trr = os.path.join(BASE_DIR, "example/ethylene/eth.trr")
    mapping = os.path.join(BASE_DIR, "example/ethylene/mapping_two_beads.xml")
    options = os.path.join(BASE_DIR, "example/ethylene/stat_options.xml")
    out_dir = os.path.join(BASE_DIR, "example/ethylene")
    if not (os.path.isfile(tpr) and os.path.isfile(trr)):
        return jsonify({"ok": False, "error": "missing tpr/trr"}), 400
    cmd = ["bash", "-lc", f"docker run --rm -v '{BASE_DIR}':'/work' -w /work/example/ethylene votca/votca:latest csg_stat --top /work/example/ethylene/eth.tpr --trj /work/example/ethylene/eth.trr --cg /work/example/ethylene/mapping_two_beads.xml --options /work/example/ethylene/stat_options.xml --ext dist.tgt"]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    ab_path = os.path.join(out_dir, "AB.dist.tgt")
    xs, ys = [], []
    if os.path.isfile(ab_path):
        with open(ab_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        xs.append(float(parts[0]))
                        ys.append(float(parts[1]))
                    except Exception:
                        pass
    return jsonify({"ok": True, "stdout": res.stdout, "stderr": res.stderr, "x": xs, "y": ys})

@app.route("/api/gmx/table", methods=["POST"])
def gmx_table():
    out_dir = os.path.join(BASE_DIR, "example/ethylene")
    produced = []
    for name in os.listdir(out_dir):
        if name.endswith('.pot'):
            base = os.path.splitext(name)[0]
            in_path = os.path.join(out_dir, name)
            out_name = f"table_{base}.xvg"
            out_path = os.path.join(out_dir, out_name)
            try:
                res = subprocess.run([sys.executable, os.path.join(BASE_DIR, 'scripts/pot_to_table.py'), in_path, out_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if os.path.isfile(out_path):
                    produced.append(os.path.join('example/ethylene', out_name))
            except Exception:
                pass
    return jsonify({"ok": len(produced) > 0, "tables": produced})

@app.route("/api/cg/top", methods=["POST"])
def cg_top():
    tpr = os.path.join(BASE_DIR, "example/ethylene/eth.tpr")
    mapping = os.path.join(BASE_DIR, "example/ethylene/mapping_two_beads.xml")
    out_top = os.path.join(BASE_DIR, "example/ethylene/cg_conf.top")
    if not os.path.isfile(tpr):
        return jsonify({"ok": False, "error": "missing tpr"}), 400
    cmd = [
        "bash",
        "-lc",
        f"docker run --rm -v '{BASE_DIR}':'/work' -w /work/example/ethylene votca/votca:latest csg_gmxtopol --top /work/example/ethylene/eth.tpr --out /work/example/ethylene/cg_conf.top --cg /work/example/ethylene/mapping_two_beads.xml",
    ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    ok = os.path.isfile(out_top)
    return jsonify({"ok": ok, "stdout": res.stdout, "stderr": res.stderr, "top": "example/ethylene/cg_conf.top" if ok else None})

@app.route("/api/convert/pdb", methods=["POST"])
def convert_pdb():
    path = request.form.get("path") or (request.json.get("path") if request.is_json else None)
    if not path:
        return jsonify({"ok": False, "error": "missing path"}), 400
    in_abs = os.path.abspath(os.path.join(BASE_DIR, path)) if not os.path.isabs(path) else path
    if not in_abs.startswith(BASE_DIR):
        return jsonify({"ok": False, "error": "forbidden"}), 403
    if not os.path.isfile(in_abs):
        return jsonify({"ok": False, "error": "not found"}), 404
    out_abs = os.path.splitext(in_abs)[0] + ".pdb"
    if shutil.which("conda"):
        cmd = ["conda", "run", "-n", "cg", "bash", "-lc", f"gmx editconf -f '{in_abs}' -o '{out_abs}'"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        ok = res.returncode == 0 and os.path.isfile(out_abs)
        return jsonify({"ok": ok, "stdout": res.stdout, "stderr": res.stderr, "out": out_abs.replace(BASE_DIR+"/", "")})
    return jsonify({"ok": False, "error": "conda cg env required"}), 500

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5050)