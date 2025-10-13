# ----- 106: scan aggregator (parse-time setup) -------------------------------
from pathlib import Path
import json

CONFIG_DIR_106 = "config"
DUMP_DIR = globals().get("DUMP_DIR", "dump")
BIN_DIR  = globals().get("BIN_DIR",  "build/bin/scripts")
LOG_DIR  = globals().get("LOG_DIR",  "log")

def _read_json(p):
    with open(p, "r") as fh:
        return json.load(fh)

# Map output_file_name -> config path (ensure uniqueness)
OUTFILE_TO_CFG = {}
for p in sorted(Path(CONFIG_DIR_106).glob("106_*.json")):
    cfg = _read_json(p)
    out_name = cfg.get("output_file_name")
    if not out_name:
        # fallback: use stem + .root if not provided
        out_name = f"{p.stem}.root"
    if out_name in OUTFILE_TO_CFG and str(OUTFILE_TO_CFG[out_name]) != str(p):
        raise ValueError(f"Duplicate output_file_name '{out_name}' in {p} and {OUTFILE_TO_CFG[out_name]}")
    OUTFILE_TO_CFG[out_name] = str(p)

def _inputs_for_outfile(outfile_name: str):
    """Get upstream inputs for a given final outfile name, using the mapped config."""
    cfg_path = OUTFILE_TO_CFG[outfile_name]
    cfg = _read_json(cfg_path)
    ins = []
    for n in cfg.get("run_numbers", []):
        rn = f"{int(n):04d}"
        ins.append(f"{DUMP_DIR}/102_EventMatch/beamtests/Run{rn}.root")
    return ins

# All final outputs discovered now (parse-time)
ALL_106_OUTPUTS = [f"{DUMP_DIR}/106_ADC_Compare/{name}" for name in OUTFILE_TO_CFG.keys()]

# ----- 106: main rule --------------------------------------------------------
rule ADC_Compare_106:
    """
    Build one scan result from its output filename.
    The corresponding config is looked up by OUTFILE_TO_CFG[outfile].
    """
    input:
        exe  = f"{BIN_DIR}/106_ADC_Compare",
        cfg  = lambda wc: OUTFILE_TO_CFG[f"{wc.outfile}.root"],
        runs = lambda wc: _inputs_for_outfile(f"{wc.outfile}.root")
    output:
        out  = f"{DUMP_DIR}/106_ADC_Compare/{{outfile}}.root"
    log:
        f"{LOG_DIR}/106_ADC_Compare/{{outfile}}.log"
    shell:
        r"""
        mkdir -p "{DUMP_DIR}/106_ADC_Compare" "{LOG_DIR}/106_ADC_Compare"
        "{input.exe}" -c "{input.cfg}" -i {input.runs} -o "{output.out}" > "{log}" 2>&1
        """

# ----- 106: convenience target ----------------------------------------------
rule run_all_106:
    """
    Run all configs starting with 106_*.json in config/.
    """
    input: ALL_106_OUTPUTS