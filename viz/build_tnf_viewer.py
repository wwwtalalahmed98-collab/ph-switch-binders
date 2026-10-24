# Build a self-contained viewer for the 20 submitted TNF-alpha designs.
#
# Why this exists: about 250 complexes have been predicted in this project and not one had
# ever been looked at. Every selection decision was made from numbers in a table. A table
# cannot show you that a binder grips with a floppy tail, that an "interface" is one
# contact pretending to be ten, or that a histidine we called buried sits in a gap.
#
# Output goes to incoming/ because the structures embed the design sequences, and those
# must not enter the public repository. This script is the committed part; its output is
# not. Open the HTML directly in a browser - it needs no server.
import csv, json, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "incoming", "tnf_structures")
OUT = os.path.join(ROOT, "incoming", "tnf_viewer.html")

met = {r["name"]: r for r in csv.DictReader(
    open(os.path.join(ROOT, "submission", "challenge2_metrics_v4.csv")))}
seqs = {r["name"]: r["sequence"] for r in csv.DictReader(
    open(os.path.join(ROOT, "submission", "challenge2_proposed_v4.csv")))}
cmap = json.load(open(os.path.join(ROOT, "tnf", "results", "construct_map.json")))["order"]

# The target's two chains are concatenated into one chain renumbered from 1, so a residue
# quoted in UniProt numbering (the way every one of our result files quotes it) has to be
# mapped back through the construct map or it will point at the wrong residue entirely.
uni2idx = {(d["chain"], d["uniprot"]): i + 1 for i, d in enumerate(cmap)}

T31 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G",
           HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S",
           THR="T", TRP="W", TYR="Y", VAL="V")


def load_first_model(path):
    """Keep only the first model. The control structures carry two, which doubles the
    file and makes every residue appear twice to a viewer."""
    out, started = [], False
    for ln in open(path):
        if ln.startswith("MODEL"):
            if started:
                break
            started = True
            continue
        if ln.startswith(("ATOM", "TER")):
            out.append(ln.rstrip("\n"))
        elif ln.startswith("ENDMDL"):
            break
    return "\n".join(out)


def chains_of(text):
    chains = {}
    for ln in text.split("\n"):
        if not ln.startswith("ATOM"):
            continue
        ch, num, rn = ln[21], int(ln[22:26]), ln[17:20].strip()
        chains.setdefault(ch, {}).setdefault(num, rn)
    return chains


def binder_chain_and_seq(chains, want_len):
    for c, v in chains.items():
        if len(v) == want_len:
            return c, "".join(T31.get(v[k], "X") for k in sorted(v))
    return None, None


designs = []
for name in sorted(seqs):
    path = os.path.join(SRC, name + ".pdb")
    if not os.path.exists(path):
        print("skipping %s: no structure" % name)
        continue
    m = met[name]
    want = seqs[name]
    pdb = load_first_model(path)
    chains = chains_of(pdb)
    chain, got = binder_chain_and_seq(chains, len(want))
    if chain is None:
        print("skipping %s: no chain of length %d" % (name, len(want)))
        continue

    # Say plainly where the structure is not the molecule that was submitted. Five of the
    # eight controls were predicted before the repair substitutions were applied, so the
    # structure is the parent, not the final sequence. Showing it silently would make the
    # switch-versus-control comparison look cleaner than it is.
    diffs = [(i + 1, got[i], want[i]) for i in range(len(want)) if got[i] != want[i]]

    d = dict(name=name, role=m["role"], pdb=pdb, chain=chain, length=int(m["length"]),
             paired_with=m["paired_with"], internal=m["internal_id"],
             human_iptm=m["human_i_ptm"], mouse_iptm=m["mouse_i_ptm"],
             human_plddt=m["human_plddt"],
             mismatch=["%s%d%s" % (a, i, b) for i, a, b in diffs])

    if m["role"] == "switch":
        sub = m["his_substitutions"]
        d["sub"] = sub
        d["his_resi"] = int(sub[1:-1])
        d["burial"] = m["burial_on_binding"]
        d["charge"] = m["charge_lys_minus_gln"]
        d["anchor"] = m["his_partner"]
        mt = re.match(r"^([A-Z])(\d+)([A-Z])$", m["his_partner"] or "")
        if mt:
            # The target chain keeps the gapped numbering of the original crop - 122
            # residues spread over 1..322 - so the Nth residue of the construct map is
            # NOT residue number N in the file. Converting the ordinal through the
            # structure's own residue numbers is the difference between pointing at
            # Arg108 and pointing at an unrelated lysine eighty residues away.
            tnums = sorted(chains["A"])
            idx = uni2idx.get((mt.group(3), int(mt.group(2))))
            resi = tnums[idx - 1] if idx and idx <= len(tnums) else None
            d["anchor_resi"] = resi
            # Assert against the STRUCTURE, not against the map. Checking the map
            # against itself is what let the wrong residue through the first time.
            d["anchor_ok"] = bool(resi) and T31.get(chains["A"][resi]) == mt.group(1)
    designs.append(d)

print("packed %d designs (%d switches, %d controls)"
      % (len(designs), sum(1 for d in designs if d["role"] == "switch"),
         sum(1 for d in designs if d["role"] == "control")))
bad = [d["name"] for d in designs if d.get("anchor_resi") and not d.get("anchor_ok")]
if bad:
    raise SystemExit("anchor mapping disagrees with the construct map for %s" % bad)
flagged = [d["name"] for d in designs if d["mismatch"]]
print("structures that are NOT the submitted sequence: %s"
      % (", ".join(flagged) if flagged else "none"))

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TNF binders</title>
<style>
:root{--bg:#faf9f7;--panel:#fff;--ink:#23231f;--dim:#6b6a64;--line:#e2e0d9;
--accent:#1d9e75;--warn:#993c1d;--sw:#534ab7;--ct:#888780}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#17171a;
--panel:#1f1f23;--ink:#e8e7e2;--dim:#9c9a92;--line:#33333a;--accent:#5dcaa5;
--warn:#f0997b;--sw:#afa9ec;--ct:#b4b2a9}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.6 ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
header{padding:18px 20px 12px;border-bottom:1px solid var(--line)}
h1{margin:0 0 4px;font-size:19px;font-weight:500}
.sub{color:var(--dim);font-size:13px;max-width:70ch}
.wrap{display:flex;gap:0;align-items:stretch;min-height:560px}
.side{width:196px;flex:none;border-right:1px solid var(--line);overflow-y:auto;
max-height:78vh;padding:8px}
.grp{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--dim);
padding:10px 8px 4px}
.item{padding:6px 8px;border-radius:6px;cursor:pointer;font-size:13px;
display:flex;justify-content:space-between;gap:6px;align-items:baseline}
.item:hover{background:var(--line)}
.item.on{background:var(--sw);color:#fff}
.item.on .meta{color:#fff;opacity:.85}
.meta{font-size:11px;color:var(--dim);font-variant-numeric:tabular-nums}
.main{flex:1;display:flex;flex-direction:column;min-width:0}
.stage{display:flex;gap:10px;padding:10px}
.vwrap{position:relative;flex:1;min-width:0}
.vlab{position:absolute;top:8px;left:10px;z-index:5;font-size:12px;color:var(--dim);
background:var(--panel);padding:2px 7px;border-radius:5px;border:1px solid var(--line)}
.v{width:100%;height:460px;position:relative;border:1px solid var(--line);
border-radius:10px;overflow:hidden;background:var(--panel)}
.bar{display:flex;flex-wrap:wrap;gap:14px;padding:4px 20px 12px;align-items:center;
font-size:13px;color:var(--dim)}
label.ck{display:flex;gap:5px;align-items:center;cursor:pointer;user-select:none}
.facts{padding:0 20px 20px;display:flex;flex-wrap:wrap;gap:8px}
.f{border:1px solid var(--line);border-radius:9px;padding:7px 11px;background:var(--panel)}
.f b{display:block;font-size:11px;font-weight:400;color:var(--dim);
letter-spacing:.04em;text-transform:uppercase}
.f span{font-size:15px;font-variant-numeric:tabular-nums}
.note{margin:0 20px 18px;padding:9px 12px;border:1px solid var(--line);
border-left:3px solid var(--warn);border-radius:0;background:var(--panel);
font-size:13px;color:var(--dim)}
.key{display:flex;gap:12px;font-size:12px;color:var(--dim);padding:0 20px 14px;
flex-wrap:wrap}
.sw{display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:4px;
vertical-align:-1px}
select,button{font:inherit;font-size:13px;padding:3px 7px;border-radius:6px;
border:1px solid var(--line);background:var(--panel);color:var(--ink);cursor:pointer}
</style></head><body>
<header>
<h1>TNF-&alpha; designs, as predicted</h1>
<div class="sub">The twenty structures submitted for Challenge 2. Confidence colouring is
the model's own per-residue pLDDT, not a measurement. Nothing here has been tested in a
laboratory.</div>
</header>
<div class="wrap">
  <div class="side" id="side"></div>
  <div class="main">
    <div class="stage">
      <div class="vwrap"><div class="vlab" id="labA"></div><div class="v" id="vA"></div></div>
      <div class="vwrap" id="wrapB" style="display:none">
        <div class="vlab" id="labB"></div><div class="v" id="vB"></div></div>
    </div>
    <div class="bar">
      <label class="ck"><input type="checkbox" id="twin"> show its control twin</label>
      <label class="ck"><input type="checkbox" id="his" checked> histidine</label>
      <label class="ck"><input type="checkbox" id="anch" checked> anchor it must bury against</label>
      <label class="ck"><input type="checkbox" id="surf"> target surface</label>
      <label class="ck"><input type="checkbox" id="spin"> spin</label>
      <span>colour <select id="mode">
        <option value="plddt">by confidence</option>
        <option value="chain">target vs binder</option>
        <option value="spectrum">binder N&rarr;C</option>
      </select></span>
    </div>
    <div class="key" id="key"></div>
    <div class="facts" id="facts"></div>
    <div class="note" id="note" style="display:none"></div>
  </div>
</div>
<script>window.__D = __DATA__;</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<script>
if (typeof $3Dmol === "undefined") {
  document.write('<scr'+'ipt src="https://3dmol.org/build/3Dmol-min.js"></scr'+'ipt>');
}
</script>
<script>
var D = window.__D, byName = {}, vA = null, vB = null;
// open on a switch, not whichever name sorts first - the switches are the
// point of the set and the controls only make sense beside one
var cur = D.filter(function(x){ return x.role === "switch"; })[0] || D[0];
D.forEach(function(d){ byName[d.name] = d; });

function el(id){ return document.getElementById(id); }

function buildSide(){
  var h = '', roles = [["switch","acid switches"],["control","controls"]];
  roles.forEach(function(r){
    h += '<div class="grp">' + r[1] + '</div>';
    D.filter(function(d){ return d.role === r[0]; }).forEach(function(d){
      var m = d.role === "switch" ? ("burial " + d.burial) : (d.length + " aa");
      h += '<div class="item" data-n="' + d.name + '"><span>' + d.name +
           '</span><span class="meta">' + m + '</span></div>';
    });
  });
  el("side").innerHTML = h;
  el("side").querySelectorAll(".item").forEach(function(n){
    n.onclick = function(){ cur = byName[n.dataset.n]; render(); };
  });
}

function styleOne(v, d){
  var mode = el("mode").value;
  v.setStyle({}, {});
  v.setStyle({chain: "A"}, {cartoon: {color: "#9b9a93", opacity: 0.95}});
  if (mode === "plddt"){
    // red is low confidence, blue is high, matching the convention AlphaFold's own
    // viewers use. Getting this the wrong way round would make the worst-modelled
    // loops look like the best part of the design.
    v.setStyle({chain: d.chain}, {cartoon: {colorscheme:
      {prop: "b", gradient: "rwb", min: 45, max: 95}}});
  } else if (mode === "chain"){
    v.setStyle({chain: d.chain}, {cartoon: {color: "#534ab7"}});
  } else {
    v.setStyle({chain: d.chain}, {cartoon: {color: "spectrum"}});
  }
  if (el("surf").checked){
    v.addSurface($3Dmol.SurfaceType.VDW,
      {opacity: 0.55, color: "#c9c7bf"}, {chain: "A"});
  }
  if (el("his").checked && d.his_resi){
    v.addStyle({chain: d.chain, resi: d.his_resi},
      {stick: {radius: 0.22, colorscheme: "greenCarbon"}});
    v.addResLabels({chain: d.chain, resi: d.his_resi},
      {font: "sans-serif", fontSize: 11, fontColor: "white",
       backgroundColor: "#0f6e56", backgroundOpacity: 0.85, showBackground: true});
  }
  if (el("anch").checked && d.anchor_resi){
    v.addStyle({chain: "A", resi: d.anchor_resi},
      {stick: {radius: 0.2, colorscheme: "yellowCarbon"}});
    v.addResLabels({chain: "A", resi: d.anchor_resi},
      {font: "sans-serif", fontSize: 11, fontColor: "black",
       backgroundColor: "#FAC775", backgroundOpacity: 0.85, showBackground: true});
  }
  v.setBackgroundColor("white");
}

function mount(id, d){
  var v = $3Dmol.createViewer(el(id), {backgroundColor: "white"});
  v.addModel(d.pdb, "pdb");
  styleOne(v, d);
  v.zoomTo({chain: d.chain});
  v.zoom(0.75);
  v.render();
  v.spin(el("spin").checked ? "y" : false);
  return v;
}

function facts(d){
  function f(lab, val){
    return '<div class="f"><b>' + lab + '</b><span>' + val + '</span></div>';
  }
  var h = f("length", d.length + " aa") +
          f("confidence, human", d.human_iptm) +
          f("confidence, mouse", d.mouse_iptm);
  if (d.role === "switch"){
    h += f("substitution", d.sub) +
         f("buried on binding", d.burial) +
         f("charge effect", d.charge) +
         f("anchor", d.anchor || "-");
  } else {
    h += f("role", "control, no histidine") + f("twin", d.paired_with || "-");
  }
  el("facts").innerHTML = h;

  var n = el("note");
  if (d.mismatch && d.mismatch.length){
    n.style.display = "";
    n.innerHTML = "<b>This structure is not the molecule we submitted.</b> It was " +
      "predicted before the repair substitutions were applied, so it differs at " +
      d.mismatch.length + " position" + (d.mismatch.length > 1 ? "s" : "") + ": " +
      d.mismatch.join(", ") + ". The shape is a close stand-in, not the final sequence.";
  } else { n.style.display = "none"; }

  var k = el("key");
  if (el("mode").value === "plddt"){
    k.innerHTML = '<span><span class="sw" style="background:#e24b4a"></span>' +
      'low confidence (pLDDT 45)</span>' +
      '<span><span class="sw" style="background:#378add"></span>high (95)</span>' +
      '<span><span class="sw" style="background:#9b9a93"></span>TNF target</span>';
  } else {
    k.innerHTML = '<span><span class="sw" style="background:#9b9a93"></span>' +
      'TNF target</span><span><span class="sw" style="background:#534ab7"></span>' +
      'designed binder</span>';
  }
}

function render(){
  el("side").querySelectorAll(".item").forEach(function(n){
    n.classList.toggle("on", n.dataset.n === cur.name);
  });
  el("vA").innerHTML = ""; vA = mount("vA", cur);
  el("labA").textContent = cur.name + (cur.role === "switch" ? "  (switch)" : "  (control)");

  var twin = el("twin").checked && cur.paired_with && byName[cur.paired_with];
  el("wrapB").style.display = twin ? "" : "none";
  if (twin){
    var t = byName[cur.paired_with];
    el("vB").innerHTML = ""; vB = mount("vB", t);
    el("labB").textContent = t.name + (t.role === "switch" ? "  (switch)" : "  (control)");
  } else { vB = null; }
  facts(cur);
}

["his","anch","surf","twin","mode"].forEach(function(id){
  el(id).onchange = render;
});
el("spin").onchange = function(){
  [vA, vB].forEach(function(v){ if (v) v.spin(el("spin").checked ? "y" : false); });
};

buildSide();
render();
</script></body></html>"""

html = HTML.replace("__DATA__", json.dumps(designs, separators=(",", ":")))
with open(OUT, "w", encoding="utf-8") as fh:
    fh.write(html)
print("wrote %s  (%.2f MB)" % (OUT, os.path.getsize(OUT) / 1e6))
