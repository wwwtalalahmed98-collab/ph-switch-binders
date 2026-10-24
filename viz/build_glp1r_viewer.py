# Build a self-contained viewer for GLP1R complexes.
#
# This one has to show something the TNF viewer did not: DEPTH. The whole challenge turns
# on how far a design reaches into the helical bundle, and on 23 Oct we established that
# the depth NUMBER cannot be trusted to rank candidates - the real hormone reached 34.1 A,
# the blocker-shaped truncation 33.3 A, and a scrambled nonsense sequence 37.0 A. So the
# page draws the depth instead of reporting it: a ruler down the receptor axis with the
# hormone's own reach marked, and the N-terminus of each design marked, because polarity
# (which end points down) is the one directional signal that did separate sense from
# nonsense.
#
# It also draws where our crop stops. glp1r_thread_001 scored the highest i_ptm this
# project has produced by threading into space that only exists because we cut the bottom
# off the receptor; with the ruler on screen that is obvious at a glance instead of
# needing a separate analysis.
#
# Runs unchanged on the design structures when they land - it globs whatever complexes are
# present. Output goes to incoming/ because structures embed sequences.
import glob, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "incoming", "glp1r_structures")
OUT = os.path.join(ROOT, "incoming", "glp1r_viewer.html")

T31 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G",
           HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S",
           THR="T", TRP="W", TYR="Y", VAL="V")

M = json.load(open(os.path.join(SRC, "glp1r_map.json")))
ORDER = M["order"]
AXIS = M["axis"]
ECD_C = M["ecd_centroid"]
NATIVE_DEEP = max(M["native_glp1"]["depth"].values())

ctrl = {}
cpath = os.path.join(SRC, "glp1r_control.json")
if os.path.exists(cpath):
    for r in json.load(open(cpath)):
        if "seq" in r:
            ctrl[r["seq"]] = r

LABELS = {
    "GLP-1_7-36": ("The real hormone", "activates the receptor"),
    "GLP-1_9-36": ("Blocker shape", "exendin(9-39) is this molecule, a drug that blocks"),
    "scrambled": ("Scrambled", "same amino acids, shuffled: meaningless"),
    "poly-alanine,": ("Poly-alanine", "a floor to measure against"),
}


def load(path):
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
    ch = {}
    for ln in text.split("\n"):
        if ln.startswith("ATOM"):
            ch.setdefault(ln[21], {}).setdefault(int(ln[22:26]), ln[17:20].strip())
    return ch


items = []
for path in sorted(glob.glob(os.path.join(SRC, "*.pdb"))):
    base = os.path.basename(path)[:-4]
    if base.startswith("GLP1R_") and "crop" in base:
        continue                                   # the bare target, not a complex
    pdb = load(path)
    ch = chains_of(pdb)
    # The receptor is the chain whose residue count matches the map. Never "the bigger
    # chain": that assumption is what put the TNF viewer's anchor marker eighty residues
    # off the residue it was meant to point at.
    rid = next((c for c, v in ch.items() if len(v) == len(ORDER)), None)
    if rid is None:
        print("skipping %s: no chain of %d residues (have %s)"
              % (base, len(ORDER), {c: len(v) for c, v in ch.items()}))
        continue
    bid = [c for c in ch if c != rid]
    if len(bid) != 1:
        print("skipping %s: expected one binder chain, got %s" % (base, bid))
        continue
    bid = bid[0]
    bnums = sorted(ch[bid])
    seq = "".join(T31.get(ch[bid][k], "X") for k in bnums)

    # ordinal -> real residue number, taken from the structure, for the receptor regions
    rnums = sorted(ch[rid])
    ecd = [rnums[i] for i, d in enumerate(ORDER) if d["region"] == "ecd"]
    tm = [rnums[i] for i, d in enumerate(ORDER) if d["region"] == "tm"]
    # sanity: the map's amino acids must match the structure, or the regions are wrong
    bad = [i for i, d in enumerate(ORDER) if T31.get(ch[rid][rnums[i]]) != d["aa"]]
    if bad:
        raise SystemExit("%s: receptor sequence disagrees with the map at %d positions "
                         "(first at ordinal %d)" % (base, len(bad), bad[0] + 1))

    title, note = LABELS.get(base, (base, ""))
    c = ctrl.get(seq, {})
    items.append(dict(key=base, title=title, note=note, pdb=pdb, rid=rid, bid=bid,
                      seq=seq, length=len(seq), ecd=ecd, tm=tm,
                      nterm=bnums[0], cterm=bnums[-1],
                      i_ptm=c.get("i_ptm"), plddt=c.get("plddt"),
                      deepest=c.get("deepest"), nterm_depth=c.get("nterm_depth"),
                      tm_contacts=c.get("tm_contacts"),
                      ecd_contacts=c.get("ecd_contacts")))

if not items:
    raise SystemExit("no complexes found in %s" % SRC)
print("packed %d complexes: %s" % (len(items), ", ".join(i["key"] for i in items)))

DATA = dict(items=items, axis=AXIS, ecd_c=ECD_C, native_deep=round(NATIVE_DEEP, 1),
            crop_bottom=46.8)

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GLP-1R complexes</title>
<style>
:root{--bg:#faf9f7;--panel:#fff;--ink:#23231f;--dim:#6b6a64;--line:#e2e0d9;
--accent:#1d9e75;--warn:#993c1d}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#17171a;
--panel:#1f1f23;--ink:#e8e7e2;--dim:#9c9a92;--line:#33333a;--accent:#5dcaa5;
--warn:#f0997b}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.6 ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
header{padding:18px 20px 12px;border-bottom:1px solid var(--line)}
h1{margin:0 0 4px;font-size:19px;font-weight:500}
.sub{color:var(--dim);font-size:13px;max-width:76ch}
.wrap{display:flex;min-height:560px}
.side{width:230px;flex:none;border-right:1px solid var(--line);overflow-y:auto;
max-height:80vh;padding:8px}
.item{padding:7px 9px;border-radius:7px;cursor:pointer;font-size:13px;margin-bottom:2px}
.item:hover{background:var(--line)}
.item.on{background:var(--accent);color:#fff}
.item.on .m{color:#fff;opacity:.9}
.m{display:block;font-size:11px;color:var(--dim)}
.main{flex:1;display:flex;flex-direction:column;min-width:0}
.stage{display:flex;gap:10px;padding:10px}
.vwrap{position:relative;flex:1;min-width:0}
.vlab{position:absolute;top:8px;left:10px;z-index:5;font-size:12px;color:var(--dim);
background:var(--panel);padding:2px 7px;border-radius:5px;border:1px solid var(--line)}
.v{width:100%;height:470px;border:1px solid var(--line);border-radius:10px;
overflow:hidden;background:var(--panel);position:relative}
.bar{display:flex;flex-wrap:wrap;gap:14px;padding:4px 20px 10px;align-items:center;
font-size:13px;color:var(--dim)}
label.ck{display:flex;gap:5px;align-items:center;cursor:pointer;user-select:none}
.key{display:flex;gap:12px;font-size:12px;color:var(--dim);padding:0 20px 12px;
flex-wrap:wrap}
.sw{display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:4px;
vertical-align:-1px}
.facts{padding:0 20px 18px;display:flex;flex-wrap:wrap;gap:8px}
.f{border:1px solid var(--line);border-radius:9px;padding:7px 11px;background:var(--panel)}
.f b{display:block;font-size:11px;font-weight:400;color:var(--dim);
letter-spacing:.04em;text-transform:uppercase}
.f span{font-size:15px;font-variant-numeric:tabular-nums}
.note{margin:0 20px 18px;padding:9px 12px;border:1px solid var(--line);
border-left:3px solid var(--warn);background:var(--panel);font-size:13px;color:var(--dim)}
select{font:inherit;font-size:13px;padding:3px 7px;border-radius:6px;
border:1px solid var(--line);background:var(--panel);color:var(--ink);cursor:pointer}
</style></head><body>
<header>
<h1>GLP-1R, and how deep things reach</h1>
<div class="sub">The receptor is cut off at the dashed marker because we cropped it to
make the design job affordable. Anything reaching past that is in space our target
invented. Depth numbers alone were shown on 23 October not to separate an activating
peptide from a blocking one &mdash; which is why this draws them instead.</div>
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
      <span>compare with <select id="cmp"><option value="">nothing</option></select></span>
      <label class="ck"><input type="checkbox" id="ruler" checked> depth ruler</label>
      <label class="ck"><input type="checkbox" id="ends" checked> mark N and C ends</label>
      <label class="ck"><input type="checkbox" id="spin"> spin</label>
      <span>binder <select id="mode">
        <option value="solid">one colour</option>
        <option value="plddt">by confidence</option>
        <option value="spectrum">N&rarr;C</option>
      </select></span>
    </div>
    <div class="key">
      <span><span class="sw" style="background:#b4b2a9"></span>extracellular domain</span>
      <span><span class="sw" style="background:#5dcaa5"></span>seven-helix bundle</span>
      <span><span class="sw" style="background:#534ab7"></span>peptide</span>
      <span><span class="sw" style="background:#1d9e75"></span>N-terminus, the end that must go deep</span>
      <span><span class="sw" style="background:#BA7517"></span>C-terminus</span>
    </div>
    <div class="facts" id="facts"></div>
    <div class="note" id="note" style="display:none"></div>
  </div>
</div>
<script>window.__G = __DATA__;</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<script>
if (typeof $3Dmol === "undefined") {
  document.write('<scr'+'ipt src="https://3dmol.org/build/3Dmol-min.js"></scr'+'ipt>');
}
</script>
<script>
var G = window.__G, D = G.items, byKey = {}, cur = D[0], vA = null, vB = null;
D.forEach(function(d){ byKey[d.key] = d; });
function el(id){ return document.getElementById(id); }
function pt(depth){
  return {x: G.ecd_c[0] + G.axis[0]*depth,
          y: G.ecd_c[1] + G.axis[1]*depth,
          z: G.ecd_c[2] + G.axis[2]*depth};
}

function build(){
  el("side").innerHTML = D.map(function(d){
    var sub = (d.deepest != null) ? ("reaches " + d.deepest + " Å") : (d.length + " aa");
    return '<div class="item" data-k="' + d.key + '"><span>' + d.title +
           '</span><span class="m">' + sub + '</span></div>';
  }).join("");
  el("side").querySelectorAll(".item").forEach(function(n){
    n.onclick = function(){ cur = byKey[n.dataset.k]; render(); };
  });
  el("cmp").innerHTML = '<option value="">nothing</option>' +
    D.map(function(d){ return '<option value="' + d.key + '">' + d.title + '</option>'; })
     .join("");
}

function styleOne(v, d){
  var mode = el("mode").value;
  v.setStyle({}, {});
  v.setStyle({chain: d.rid, resi: d.ecd}, {cartoon: {color: "#b4b2a9", opacity: 0.9}});
  v.setStyle({chain: d.rid, resi: d.tm}, {cartoon: {color: "#5dcaa5", opacity: 0.9}});
  if (mode === "plddt"){
    v.setStyle({chain: d.bid}, {cartoon: {colorscheme:
      {prop: "b", gradient: "rwb", min: 45, max: 95}}, stick: {radius: 0.12}});
  } else if (mode === "spectrum"){
    v.setStyle({chain: d.bid}, {cartoon: {color: "spectrum"}, stick: {radius: 0.12}});
  } else {
    v.setStyle({chain: d.bid}, {cartoon: {color: "#534ab7"}, stick: {radius: 0.12}});
  }
  if (el("ends").checked){
    v.addStyle({chain: d.bid, resi: d.nterm},
               {sphere: {radius: 1.9, color: "#1d9e75"}});
    v.addStyle({chain: d.bid, resi: d.cterm},
               {sphere: {radius: 1.9, color: "#BA7517"}});
  }
  if (el("ruler").checked){
    v.addCylinder({start: pt(-6), end: pt(G.crop_bottom + 6), radius: 0.16,
                   color: "#888780", dashed: true, fromCap: 1, toCap: 1});
    // where the natural hormone stops reaching
    v.addSphere({center: pt(G.native_deep), radius: 1.3, color: "#378add"});
    v.addLabel("hormone stops  " + G.native_deep + " Å",
      {position: pt(G.native_deep), fontSize: 11, fontColor: "white",
       backgroundColor: "#185FA5", backgroundOpacity: 0.85, inFront: true});
    // where OUR CROP stops - past this the target is empty because we cut it off
    v.addSphere({center: pt(G.crop_bottom), radius: 1.3, color: "#D85A30"});
    v.addLabel("our crop ends  " + G.crop_bottom + " Å",
      {position: pt(G.crop_bottom), fontSize: 11, fontColor: "white",
       backgroundColor: "#993C1D", backgroundOpacity: 0.85, inFront: true});
  }
  v.setBackgroundColor("white");
}

function mount(id, d){
  var v = $3Dmol.createViewer(el(id), {backgroundColor: "white"});
  v.addModel(d.pdb, "pdb");
  styleOne(v, d);
  v.zoomTo({chain: d.rid});
  v.zoom(0.9);
  v.render();
  v.spin(el("spin").checked ? "y" : false);
  return v;
}

function facts(d){
  function f(lab, val){
    return '<div class="f"><b>' + lab + '</b><span>' + val + '</span></div>';
  }
  var h = f("length", d.length + " aa");
  if (d.i_ptm != null) h += f("confidence", d.i_ptm);
  if (d.deepest != null){
    h += f("deepest point", d.deepest + " Å");
    h += f("N-terminus depth", d.nterm_depth + " Å");
    h += f("threaded?", (d.nterm_depth >= d.deepest - 5) ? "yes, N-end is deepest"
                                                         : "NO, it went in backwards");
  }
  if (d.tm_contacts != null)
    h += f("contacts", d.tm_contacts + " bundle / " + d.ecd_contacts + " surface");
  el("facts").innerHTML = h;

  var n = el("note"), msgs = [];
  if (d.note) msgs.push(d.note);
  if (d.deepest != null && d.nterm_depth < d.deepest - 5)
    msgs.push("<b>This one is in backwards.</b> Its deepest atom sits at " + d.deepest +
      " Å but its N-terminus is at " + d.nterm_depth + " Å. Depth alone would " +
      "have scored it as reaching deeper than the real hormone.");
  if (d.deepest != null && d.deepest > G.crop_bottom - 2)
    msgs.push("<b>This reaches the bottom edge of our crop.</b> Past " + G.crop_bottom +
      " Å the target is empty only because we cut the receptor off there.");
  n.style.display = msgs.length ? "" : "none";
  n.innerHTML = msgs.join("<br><br>");
}

function render(){
  el("side").querySelectorAll(".item").forEach(function(n){
    n.classList.toggle("on", n.dataset.k === cur.key);
  });
  el("vA").innerHTML = ""; vA = mount("vA", cur);
  el("labA").textContent = cur.title;
  var k = el("cmp").value, other = k && byKey[k];
  el("wrapB").style.display = other ? "" : "none";
  if (other){
    el("vB").innerHTML = ""; vB = mount("vB", other);
    el("labB").textContent = other.title;
  } else { vB = null; }
  facts(cur);
}

["ruler","ends","mode","cmp"].forEach(function(id){ el(id).onchange = render; });
el("spin").onchange = function(){
  [vA, vB].forEach(function(v){ if (v) v.spin(el("spin").checked ? "y" : false); });
};
build();
render();
</script></body></html>"""

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write(HTML.replace("__DATA__", json.dumps(DATA, separators=(",", ":"))))
print("wrote %s  (%.2f MB)" % (OUT, os.path.getsize(OUT) / 1e6))
