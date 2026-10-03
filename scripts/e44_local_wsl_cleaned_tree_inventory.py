#!/usr/bin/env python3
from pathlib import Path
import json, os, time

ROOTS = [
    Path.home() / "stress-energy-closure" / "eboss_workspace",
    Path.home() / "stress-energy-closure",
]

patterns = (
    "cleaned_halo_info",
    "cleaning",
    "merger",
    "tree",
    "association",
    "mainprog",
)

seen=set()
hits=[]

for root in ROOTS:
    if not root.exists():
        continue
    for dirpath, dirnames, filenames in os.walk(root):
        p=Path(dirpath)
        # avoid obvious large irrelevant/generated dirs
        dirnames[:] = [d for d in dirnames if d not in {".git",".venv","venv","__pycache__"}]
        # cap depth relative to root
        try:
            depth=len(p.relative_to(root).parts)
        except ValueError:
            depth=0
        if depth>10:
            dirnames[:] = []
            continue

        lowdir=str(p).lower()
        if any(x in lowdir for x in patterns):
            key=("dir",str(p))
            if key not in seen:
                seen.add(key)
                hits.append({"type":"dir","path":str(p)})

        for fn in filenames:
            low=fn.lower()
            if any(x in low for x in patterns):
                fp=p/fn
                try:
                    st=fp.stat()
                    rec={"type":"file","path":str(fp),"bytes":st.st_size}
                except OSError as e:
                    rec={"type":"file","path":str(fp),"stat_error":str(e)}
                key=("file",str(fp))
                if key not in seen:
                    seen.add(key)
                    hits.append(rec)

# prioritize exact cleaned ASDFs and files with mainprog/tree naming
def score(r):
    s=r["path"].lower()
    return (
        0 if "cleaned_halo_info" in s else
        1 if "mainprog" in s else
        2 if "merger" in s or "tree" in s else
        3
    )
hits=sorted(hits,key=lambda r:(score(r),r["path"]))

out={
    "stage":"E44_LOCAL_FILESYSTEM_INVENTORY_ONLY",
    "status":"PASS",
    "roots":[str(r) for r in ROOTS],
    "hit_count":len(hits),
    "hits":hits[:500],
    "arrays_opened":False,
    "asdf_decompression_performed":False,
    "next_action":"If a cleaned_halo_info ASDF exists, stop and upload this inventory. Do NOT open it yet."
}

q=Path("source_data/e44_local_cleaned_tree_inventory.json")
q.parent.mkdir(exist_ok=True)
q.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")

print("E44_WSL_INVENTORY_PASS")
print("HIT_COUNT",len(hits))
for r in hits[:40]:
    if r["type"]=="file":
        print("FILE",r["bytes"],r["path"])
    else:
        print("DIR",r["path"])
print("ARRAYS_OPENED",False)
print("REPORT",q)
