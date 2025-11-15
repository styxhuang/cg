import os
import uuid
import sys
from flask import Flask, request, jsonify, send_from_directory
import logging
import logging.handlers
import shutil
import subprocess
from werkzeug.utils import secure_filename

app = Flask(__name__)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("backend")
log_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "backend.log")
file_handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=1048576, backupCount=3, encoding="utf-8")
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
logger.addHandler(file_handler)
logger.propagate = False

@app.before_request
def _log_in():
    try:
        logger.info(f"[backend] {request.method} {request.path} ct={request.content_type}")
    except Exception:
        pass

@app.after_request
def _log_out(resp):
    try:
        logger.info(f"[backend] {resp.status_code} {request.path}")
    except Exception:
        pass
    return resp
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)
WORKSPACE_DIR = os.path.join(BASE_DIR, "workspace")
STATIC_DIR = os.path.join(ROOT_DIR, "static")
os.makedirs(WORKSPACE_DIR, exist_ok=True)

@app.route("/")
def index():
    return send_from_directory(STATIC_DIR, "index.html")

@app.route("/static/<path:filename>")
def static_files(filename):
    return send_from_directory(STATIC_DIR, filename)

@app.route("/status")
def status():
    return jsonify({"ok": True})

@app.route("/env")
def env():
    return jsonify({"ok": True, "python": sys.version, "executable": sys.executable, "venv": os.path.exists(os.path.join(BASE_DIR, ".venv"))})

@app.route("/upload", methods=["POST"])
def upload():
    files = {}
    try:
        tp = request.form.get("top_path")
        trp = request.form.get("trj_path")
        mp = request.form.get("mapping_path")
        logger.info(f"[api/upload] incoming top_path={tp} trj_path={trp} mapping_path={mp}")
    except Exception:
        logger.info("[api/upload] incoming error")
    try:
        logger.info(f"[api/upload] form_keys={list(request.form.keys())} file_keys={list(request.files.keys())}")
    except Exception:
        pass
    for key in ["top", "trj", "mapping"]:
        f = request.files.get(key)
        if f:
            files[key] = f
    if not files:
        try:
            logger.warning("[api/upload] no files in request")
        except Exception:
            pass
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
        try:
            logger.info(f"[api/upload] saved {k} name={name} path={path}")
        except Exception:
            pass
        try:
            logger.info(f"[api/upload] saved {k} name={name} path={path}")
        except Exception:
            pass
    try:
        logger.info(f"[api/upload] rid={rid} saved={saved}")
    except Exception:
        pass
    return jsonify({"ok": True, "rid": rid, "saved": saved})


@app.errorhandler(404)
def _handle_404(e):
    try:
        logger.warning(f"[404] {request.method} {request.path}")
    except Exception:
        pass
    return jsonify({"ok": False, "error": "not found", "path": request.path}), 404

@app.errorhandler(405)
def _handle_405(e):
    try:
        logger.warning(f"[405] {request.method} {request.path}")
    except Exception:
        pass
    return jsonify({"ok": False, "error": "method not allowed", "path": request.path}), 405

def _find_workspace_inputs(rid: str):
    out_dir = os.path.join(WORKSPACE_DIR, rid)
    if not os.path.isdir(out_dir):
        return None
    mapping = None
    top = None
    trj = None
    targets = []
    stat_options = None
    for name in os.listdir(out_dir):
        p = os.path.join(out_dir, name)
        if name.endswith(".xml"):
            if name == "stat_options.xml":
                stat_options = p
            else:
                mapping = p
        elif name.endswith(".tpr") or name.endswith(".top") or name.endswith(".gro"):
            top = p
        elif name.endswith(".xtc") or name.endswith(".trr"):
            trj = p
        elif name.endswith(".dist.tgt"):
            targets.append(p)
    return {"dir": out_dir, "mapping": mapping, "top": top, "trj": trj, "targets": sorted(targets), "stat_options": stat_options}

@app.route("/votca/check")
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

@app.route("/votca/map", methods=["POST"])
def votca_map():
    rid = request.form.get("rid")
    if not rid:
        return jsonify({"ok": False, "error": "missing rid"}), 400
    info = _find_workspace_inputs(rid)
    if not info:
        return jsonify({"ok": False, "error": "rid not found"}), 404
    mapping = info["mapping"]
    top = info["top"]
    trj = info["trj"]
    if not mapping or not top or not trj:
        return jsonify({"ok": False, "error": "missing inputs (top/trj/mapping)"}), 400
    use_conda = shutil.which("conda") is not None
    out_gro = os.path.join(info["dir"], "cg_conf.gro")
    try:
        if use_conda:
            cmd = ["conda", "run", "-n", "cg", "bash", "-lc", f"csg_map --top '{top}' --cg '{mapping}' --trj '{trj}' --out '{out_gro}'"]
        else:
            if shutil.which("csg_map") is None:
                return jsonify({"ok": False, "error": "csg_map not available"}), 500
            cmd = ["csg_map", "--top", top, "--cg", mapping, "--trj", trj, "--out", out_gro]
        res = subprocess.run(cmd, cwd=info["dir"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        rel = os.path.join(os.path.relpath(info["dir"], ROOT_DIR), "cg_conf.gro")
        return jsonify({"ok": os.path.isfile(out_gro), "cg_conf": rel, "returncode": res.returncode, "stdout": res.stdout, "stderr": res.stderr, "rid": rid})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500

@app.route("/workflow/ethylene", methods=["POST"])
def workflow_ethylene():
    tpr = os.path.join(ROOT_DIR, "examples/ethylene/eth.tpr")
    trr = os.path.join(ROOT_DIR, "examples/ethylene/eth.trr")
    xtc = os.path.join(ROOT_DIR, "examples/ethylene/eth.xtc")
    gro = os.path.join(ROOT_DIR, "examples/ethylene/eth.gro")
    top = os.path.join(ROOT_DIR, "examples/ethylene/eth.top")
    mapping = os.path.join(ROOT_DIR, "examples/ethylene/mapping_two_beads.xml")
    logs = []
    if shutil.which("conda"):
        cmd1 = ["conda", "run", "-n", "cg", "bash", "-lc", f"gmx grompp -f '{ROOT_DIR}/examples/ethylene/min.mdp' -c '{gro}' -p '{top}' -o '{tpr}' -maxwarn 100"]
        res1 = subprocess.run(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "grompp", "rc": res1.returncode, "out": res1.stdout, "err": res1.stderr})
        cmd2 = ["conda", "run", "-n", "cg", "bash", "-lc", f"gmx mdrun -s '{tpr}' -o '{trr}' -x '{xtc}' -nsteps 10"]
        res2 = subprocess.run(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "mdrun", "rc": res2.returncode, "out": res2.stdout, "err": res2.stderr})
    else:
        return jsonify({"ok": False, "error": "conda cg env required"}), 500
    cg_out = os.path.join(ROOT_DIR, "examples/ethylene/cg_conf.gro")
    if shutil.which("docker"):
        cmd3 = ["bash", "-lc", f"docker run --rm -v '{ROOT_DIR}':'/work' votca/votca:latest csg_map --top /work/examples/ethylene/eth.tpr --trj /work/examples/ethylene/eth.trr --cg /work/examples/ethylene/mapping_two_beads.xml --out /work/examples/ethylene/cg_conf.gro"]
        res3 = subprocess.run(cmd3, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_map_docker", "rc": res3.returncode, "out": res3.stdout, "err": res3.stderr})
    else:
        cmd3 = ["conda", "run", "-n", "cg", "bash", "-lc", f"csg_map --top '{tpr}' --trj '{trr}' --cg '{mapping}' --out '{cg_out}'"]
        res3 = subprocess.run(cmd3, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_map_local", "rc": res3.returncode, "out": res3.stdout, "err": res3.stderr})
    exists = os.path.exists(cg_out)
    rel = "examples/ethylene/cg_conf.gro" if exists else None
    return jsonify({"ok": exists, "cg": rel, "logs": logs})

@app.route("/file")
def read_file():
    p = request.args.get("path")
    if not p:
        return jsonify({"ok": False, "error": "missing path"}), 400
    abspath = os.path.abspath(os.path.join(ROOT_DIR, p)) if not os.path.isabs(p) else p
    if not abspath.startswith(ROOT_DIR):
        return jsonify({"ok": False, "error": "forbidden"}), 403
    if not os.path.isfile(abspath):
        return jsonify({"ok": False, "error": "not found"}), 404
    with open(abspath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    return jsonify({"ok": True, "path": abspath.replace(BASE_DIR+"/", ""), "content": content})

@app.route("/workflow/ff", methods=["POST"])
def workflow_ff():
    tpr_candidates = [
        os.path.join(ROOT_DIR, "examples/ethylene/eth.tpr"),
        os.path.join(ROOT_DIR, "examples/ethylene/min.tpr"),
        os.path.join(ROOT_DIR, "examples/ethylene/min-1.tpr"),
        os.path.join(ROOT_DIR, "external/cg/example/ethylene/eth.tpr"),
        os.path.join(ROOT_DIR, "external/cg/example/ethylene/topol.tpr"),
    ]
    trr_candidates = [
        os.path.join(ROOT_DIR, "examples/ethylene/eth.trr"),
        os.path.join(ROOT_DIR, "examples/ethylene/min.trr"),
        os.path.join(ROOT_DIR, "examples/ethylene/min-1.trr"),
        os.path.join(ROOT_DIR, "external/cg/example/ethylene/eth.trr"),
    ]
    mapping_candidates = [
        os.path.join(ROOT_DIR, "examples/ethylene/mapping_two_beads.xml"),
        os.path.join(ROOT_DIR, "external/cg/example/ethylene/mapping_two_beads.xml"),
    ]
    def pick(cands):
        for p in cands:
            if os.path.isfile(p):
                return p
        return None
    tpr = pick(tpr_candidates)
    trr = pick(trr_candidates)
    mapping = pick(mapping_candidates)
    out_dir = os.path.join(ROOT_DIR, "examples/ethylene")
    logs = []
    if not (tpr and trr and os.path.isfile(tpr) and os.path.isfile(trr)):
        return jsonify({"ok": False, "error": "missing tpr/trr, please run ethylene workflow first"}), 400
    if shutil.which("docker"):
        cmd = ["bash", "-lc", f"docker run --rm -v '{ROOT_DIR}':'/work' votca/votca:latest csg_boltzmann --top /work/{os.path.relpath(tpr, ROOT_DIR)} --trj /work/{os.path.relpath(trr, ROOT_DIR)} --cg /work/{os.path.relpath(mapping, ROOT_DIR)}"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_boltzmann", "rc": res.returncode, "out": res.stdout, "err": res.stderr})
    else:
        cmd = ["conda", "run", "-n", "cg", "bash", "-lc", f"csg_boltzmann --top '{tpr}' --trj '{trr}' --cg '{mapping}'"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_boltzmann", "rc": res.returncode, "out": res.stdout, "err": res.stderr})
    produced = []
    for name in os.listdir(out_dir):
        if name.endswith('.pot') or name.endswith('.force') or name.endswith('.dist'):
            produced.append(os.path.join('examples/ethylene', name))
    return jsonify({"ok": len(produced) > 0, "files": produced, "logs": logs})

@app.route("/ff/list")
def ff_list():
    rid = request.args.get("rid")
    if not rid:
        return jsonify({"ok": False, "error": "missing rid"}), 400
    info = _find_workspace_inputs(rid)
    if not info:
        return jsonify({"ok": False, "error": "rid not found"}), 404
    out_dir = info["dir"]
    files = []
    for name in os.listdir(out_dir):
        if name.endswith('.pot') or name.endswith('.force') or name.endswith('.dist') or name.endswith('.dist.tgt'):
            files.append(os.path.join(os.path.relpath(out_dir, ROOT_DIR), name))
    return jsonify({"ok": True, "files": sorted(files), "rid": rid})

def _pick(*paths):
    for p in paths:
        if p and os.path.isfile(p):
            return p
    return None

@app.route("/inverse/config", methods=["POST"])
def inverse_config():
    rid = request.form.get("rid")
    if not rid:
        return jsonify({"ok": False, "error": "missing rid"}), 400
    info = _find_workspace_inputs(rid)
    if not info:
        return jsonify({"ok": False, "error": "rid not found"}), 404
    mapping = info["mapping"]
    tpr = info["top"]
    trr = info["trr"]
    targets = info["targets"]
    if not (tpr and trr and mapping and targets):
        return jsonify({"ok": False, "error": "missing inputs (top/trj/mapping/targets)"}), 400
    out_dir = info["dir"]
    opt_path = os.path.join(out_dir, "inverse_options.xml")
    lines = ["<inverse>", "  <general>"]
    lines += [f"    <cg>{mapping}</cg>", f"    <topology>{tpr}</topology>", f"    <trajectory>{trr}</trajectory>"]
    lines += ["  </general>", "  <potentials>", "    <nonbonded>"]
    for tgt in targets:
        base = os.path.basename(tgt).split(".")[0]
        lines += ["      <rdf>", f"        <name>{base}</name>", f"        <target>{tgt}</target>", "      </rdf>"]
    lines += ["    </nonbonded>", "  </potentials>", "  <tasks>", "    <type>ibi</type>", "  </tasks>", "</inverse>"]
    with open(opt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return jsonify({"ok": True, "options": opt_path, "targets": targets, "rid": rid})

@app.route("/inverse/run", methods=["POST"])
def inverse_run():
    rid = request.form.get("rid")
    if not rid:
        return jsonify({"ok": False, "error": "missing rid"}), 400
    info = _find_workspace_inputs(rid)
    if not info:
        return jsonify({"ok": False, "error": "rid not found"}), 404
    opt = os.path.join(info["dir"], "inverse_options.xml")
    if not os.path.isfile(opt):
        return jsonify({"ok": False, "error": "options not found, run inverse/config first"}), 400
    work = info["dir"]
    logs = []
    if shutil.which("docker"):
        cmd = ["bash", "-lc", f"docker run --rm -v '{ROOT_DIR}':'/work' -w /work/{os.path.relpath(work, ROOT_DIR)} votca/votca:latest csg_inverse --options inverse_options.xml"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_inverse_docker", "rc": res.returncode, "out": res.stdout, "err": res.stderr})
    else:
        cmd = ["conda", "run", "-n", "cg", "bash", "-lc", "csg_inverse --options inverse_options.xml"]
        res = subprocess.run(cmd, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logs.append({"step": "csg_inverse_local", "rc": res.returncode, "out": res.stdout, "err": res.stderr})
    produced = []
    for name in os.listdir(work):
        if name.endswith('.pot') or name.endswith('.force') or name.endswith('.dist'):
            produced.append(os.path.join(os.path.relpath(work, ROOT_DIR), name))
    return jsonify({"ok": len(produced) > 0, "files": produced, "logs": logs, "rid": rid})

@app.route("/ff/parse", methods=["POST"])
def ff_parse():
    path = request.form.get("path") or (request.json.get("path") if request.is_json else None)
    if not path:
        return jsonify({"ok": False, "error": "missing path"}), 400
    abspath = os.path.abspath(os.path.join(ROOT_DIR, path)) if not os.path.isabs(path) else path
    if not abspath.startswith(ROOT_DIR):
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

@app.route("/rdf/run", methods=["POST"])
def rdf_run():
    rid = request.form.get("rid")
    if not rid:
        return jsonify({"ok": False, "error": "missing rid"}), 400
    info = _find_workspace_inputs(rid)
    if not info:
        return jsonify({"ok": False, "error": "rid not found"}), 404
    tpr = info["top"]
    trr = info["trj"]
    mapping = info["mapping"]
    options = info["stat_options"]
    if not (tpr and trr and mapping and options):
        return jsonify({"ok": False, "error": "missing inputs (top/trj/mapping/stat_options.xml)"}), 400
    work = info["dir"]
    if shutil.which("docker"):
        cmd = ["bash", "-lc", f"docker run --rm -v '{ROOT_DIR}':'/work' -w /work/{os.path.relpath(work, ROOT_DIR)} votca/votca:latest csg_stat --top '{os.path.relpath(tpr, ROOT_DIR)}' --trj '{os.path.relpath(trr, ROOT_DIR)}' --cg '{os.path.relpath(mapping, ROOT_DIR)}' --options '{os.path.relpath(options, ROOT_DIR)}' --ext dist.tgt"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        cmd = ["conda", "run", "-n", "cg", "bash", "-lc", f"csg_stat --top '{tpr}' --trj '{trr}' --cg '{mapping}' --options '{options}' --ext dist.tgt"]
        res = subprocess.run(cmd, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    produced = []
    for name in os.listdir(work):
        if name.endswith(".dist.tgt"):
            produced.append(os.path.join(os.path.relpath(work, ROOT_DIR), name))
    return jsonify({"ok": len(produced) > 0, "files": produced, "stdout": res.stdout, "stderr": res.stderr, "rid": rid})

@app.route("/ff/table", methods=["POST"])
def ff_table():
    rid = request.form.get("rid")
    if not rid:
        return jsonify({"ok": False, "error": "missing rid"}), 400
    info = _find_workspace_inputs(rid)
    if not info:
        return jsonify({"ok": False, "error": "rid not found"}), 404
    out_dir = info["dir"]
    produced = []
    script = os.path.join(ROOT_DIR, 'external/cg/scripts/pot_to_table.py')
    for name in os.listdir(out_dir):
        if name.endswith('.pot'):
            base = os.path.splitext(name)[0]
            in_path = os.path.join(out_dir, name)
            out_name = f"table_{base}.xvg"
            out_path = os.path.join(out_dir, out_name)
            try:
                res = subprocess.run([sys.executable, script, in_path, out_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if os.path.isfile(out_path):
                    produced.append(os.path.join(os.path.relpath(out_dir, ROOT_DIR), out_name))
            except Exception as e:
                pass
    return jsonify({"ok": len(produced) > 0, "tables": produced, "rid": rid})

@app.route("/cg/top", methods=["POST"])
def cg_top():
    tpr = os.path.join(ROOT_DIR, "examples/ethylene/eth.tpr")
    mapping = os.path.join(ROOT_DIR, "examples/ethylene/mapping_two_beads.xml")
    out_top = os.path.join(ROOT_DIR, "examples/ethylene/cg_conf.top")
    if not os.path.isfile(tpr):
        return jsonify({"ok": False, "error": "missing tpr"}), 400
        cmd = [
        "bash",
        "-lc",
        f"docker run --rm -v '{ROOT_DIR}':'/work' -w /work/examples/ethylene votca/votca:latest csg_gmxtopol --top /work/examples/ethylene/eth.tpr --out /work/examples/ethylene/cg_conf.top --cg /work/examples/ethylene/mapping_two_beads.xml",
        ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    ok = os.path.isfile(out_top)
    return jsonify({"ok": ok, "stdout": res.stdout, "stderr": res.stderr, "top": "example/ethylene/cg_conf.top" if ok else None})

@app.route("/convert/pdb", methods=["POST"])
def convert_pdb():
    path = request.form.get("path") or (request.json.get("path") if request.is_json else None)
    if not path:
        return jsonify({"ok": False, "error": "missing path"}), 400
    in_abs = os.path.abspath(os.path.join(ROOT_DIR, path)) if not os.path.isabs(path) else path
    if not in_abs.startswith(ROOT_DIR):
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

@app.route("/routes")
def api_routes():
    try:
        routes = []
        for r in app.url_map.iter_rules():
            routes.append({"rule": str(r), "methods": list(r.methods)})
        return jsonify({"ok": True, "routes": routes})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500

if __name__ == "__main__":
    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", "5050"))
    app.run(host=host, port=port)
