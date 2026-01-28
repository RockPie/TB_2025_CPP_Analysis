# ----- 106: scan aggregator (parse-time setup) -------------------------------
from pathlib import Path
import json

CONFIG_DIR_106 = "config"
CONFIG_DIR_305 = "config"
CONFIG_DIR_312 = "config"
CONFIG_DIR_317 = "config"
DUMP_DIR = globals().get("DUMP_DIR", "dump")
BIN_DIR  = globals().get("BIN_DIR",  "build/bin/scripts")
LOG_DIR  = globals().get("LOG_DIR",  "log")

def _read_json(p):
    with open(p, "r") as fh:
        return json.load(fh)

# Map output_file_name -> config path (ensure uniqueness)
OUTFILE_TO_CFG_106 = {}
for p in sorted(Path(CONFIG_DIR_106).glob("106_*.json")):
    cfg = _read_json(p)
    out_name = cfg.get("output_file_name")
    if not out_name:
        # fallback: use stem + .root if not provided
        out_name = f"{p.stem}.root"
    if out_name in OUTFILE_TO_CFG_106 and str(OUTFILE_TO_CFG_106[out_name]) != str(p):
        raise ValueError(f"Duplicate output_file_name '{out_name}' in {p} and {OUTFILE_TO_CFG_106[out_name]}")
    OUTFILE_TO_CFG_106[out_name] = str(p)

def _inputs_for_outfile(outfile_name: str):
    """Get upstream inputs for a given final outfile name, using the mapped config."""
    cfg_path = OUTFILE_TO_CFG_106[outfile_name]
    cfg = _read_json(cfg_path)
    ins = []
    for n in cfg.get("run_numbers", []):
        rn = f"{int(n):04d}"
        # ins.append(f"{DUMP_DIR}/102_EventMatch/beamtests/Run{rn}.root")
        ins.append(f"{DUMP_DIR}/302_EventMatchX/beamtests/Run{rn}_em.root")
    return ins

OUTFILE_TO_CFG_305 = {}
for p in sorted(Path(CONFIG_DIR_305).glob("305_*.json")):
    cfg = _read_json(p)
    out_name = cfg.get("output_file_name")
    if not out_name:
        # fallback: use stem + .root if not provided
        out_name = f"{p.stem}.root"
    if out_name in OUTFILE_TO_CFG_305 and str(OUTFILE_TO_CFG_305[out_name]) != str(p):
        raise ValueError(f"Duplicate output_file_name '{out_name}' in {p} and {OUTFILE_TO_CFG_305[out_name]}")
    OUTFILE_TO_CFG_305[out_name] = str(p)

def _inputs_for_outfile_305(outfile_name: str):
    """Get upstream inputs for a given final outfile name, using the mapped config."""
    cfg_path = OUTFILE_TO_CFG_305[outfile_name]
    cfg = _read_json(cfg_path)
    ins = []
    for n in cfg.get("run_numbers", []):
        rn = f"{int(n):04d}"
        ins.append(f"{DUMP_DIR}/304_RawADC/beamtests/Run{rn}.root")
    return ins

OUTFILE_TO_CFG_312 = {}
for p in sorted(Path(CONFIG_DIR_312).glob("312_*.json")):
    cfg = _read_json(p)
    out_name = cfg.get("output_file_name")
    if not out_name:
        # fallback: use stem + .root if not provided
        out_name = f"{p.stem}.root"
    if out_name in OUTFILE_TO_CFG_312 and str(OUTFILE_TO_CFG_312[out_name]) != str(p):
        raise ValueError(f"Duplicate output_file_name '{out_name}' in {p} and {OUTFILE_TO_CFG_312[out_name]}")
    OUTFILE_TO_CFG_312[out_name] = str(p)

def _inputs_for_outfile_312(outfile_name: str):
    """Get upstream inputs for a given final outfile name, using the mapped config."""
    cfg_path = OUTFILE_TO_CFG_312[outfile_name]
    cfg = _read_json(cfg_path)
    ins = []
    for n in cfg.get("run_numbers", []):
        rn = f"{int(n):04d}"
        ins.append(f"{DUMP_DIR}/311_RawToT/beamtests/Run{rn}.root")
    return ins

OUTFILE_TO_CFG_317 = {}
for p in sorted(Path(CONFIG_DIR_317).glob("317_*.json")):
    cfg = _read_json(p)
    out_name = cfg.get("output_file_name")
    if not out_name:
        # fallback: use stem + .root if not provided
        out_name = f"{p.stem}.root"
    if out_name in OUTFILE_TO_CFG_317 and str(OUTFILE_TO_CFG_317[out_name]) != str(p):
        raise ValueError(f"Duplicate output_file_name '{out_name}' in {p} and {OUTFILE_TO_CFG_317[out_name]}")
    OUTFILE_TO_CFG_317[out_name] = str(p)

def _inputs_for_outfile_317(outfile_name: str):
    """Get upstream inputs for a given final outfile name, using the mapped config."""
    cfg_path = OUTFILE_TO_CFG_317[outfile_name]
    cfg = _read_json(cfg_path)
    ins = []
    for n in cfg.get("run_numbers", []):
        rn = f"{int(n):04d}"
        ins.append(f"{DUMP_DIR}/316_ADC_ToT/beamtests/Run{rn}.root")
    return ins

# All final outputs discovered now (parse-time)
ALL_106_OUTPUTS = [f"{DUMP_DIR}/106_ADC_Compare/{name}" for name in OUTFILE_TO_CFG_106.keys()]
ALL_305_OUTPUTS = [f"{DUMP_DIR}/305_ADC_Fit_Compare/{name}" for name in OUTFILE_TO_CFG_305.keys()]
ALL_312_OUTPUTS = [f"{DUMP_DIR}/312_ToT_Fit_Compare/{name}" for name in OUTFILE_TO_CFG_312.keys()]
ALL_317_OUTPUTS = [f"{DUMP_DIR}/317_ADC_ToT_Fit_Compare/{name}" for name in OUTFILE_TO_CFG_317.keys()]

# ----- 106: main rule --------------------------------------------------------
rule ADC_Compare_106:
    """
    Build one scan result from its output filename.
    The corresponding config is looked up by OUTFILE_TO_CFG_106[outfile].
    """
    input:
        exe  = f"{BIN_DIR}/106_ADC_Compare",
        cfg  = lambda wc: OUTFILE_TO_CFG_106[f"{wc.outfile}.root"],
        runs = lambda wc: _inputs_for_outfile(f"{wc.outfile}.root")
    output:
        out  = f"{DUMP_DIR}/106_ADC_Compare/{{outfile}}.root"
    log:
        f"{LOG_DIR}/106_ADC_Compare/{{outfile}}.log"
    shell:
        r"""
        mkdir -p "{DUMP_DIR}/106_ADC_Compare" "{LOG_DIR}/106_ADC_Compare"
        "{input.exe}" -f "{input.cfg}" -o "{output.out}" > "{log}" 2>&1
        """

# ----- 305: main rule --------------------------------------------------------
rule ADC_Fit_Compare_305:
    """
    Build one fit comparison result from its output filename.
    The corresponding config is looked up by OUTFILE_TO_CFG_106[outfile].
    """
    input:
        exe  = f"{BIN_DIR}/305_ADC_Fit_Compare",
        cfg  = lambda wc: OUTFILE_TO_CFG_305[f"{wc.outfile}.root"],
        runs = lambda wc: _inputs_for_outfile_305(f"{wc.outfile}.root")
    output:
        out  = f"{DUMP_DIR}/305_ADC_Fit_Compare/{{outfile}}.root"
    log:
        f"{LOG_DIR}/305_ADC_Fit_Compare/{{outfile}}.log"
    shell:
        r"""
        mkdir -p "{DUMP_DIR}/305_ADC_Fit_Compare" "{LOG_DIR}/305_ADC_Fit_Compare"
        "{input.exe}" -f "{input.cfg}" -o "{output.out}" > "{log}" 2>&1
        """

# ----- 312: main rule --------------------------------------------------------
rule ToT_Fit_Compare_312:
    """
    Build one ToT fit comparison result from its output filename.
    The corresponding config is looked up by OUTFILE_TO_CFG_312[outfile].
    """
    input:
        exe  = f"{BIN_DIR}/312_ToT_Fit_Compare",
        cfg  = lambda wc: OUTFILE_TO_CFG_312[f"{wc.outfile}.root"],
        runs = lambda wc: _inputs_for_outfile_312(f"{wc.outfile}.root")
    output:
        out  = f"{DUMP_DIR}/312_ToT_Fit_Compare/{{outfile}}.root"
    log:
        f"{LOG_DIR}/312_ToT_Fit_Compare/{{outfile}}.log"
    shell:
        r"""
        mkdir -p "{DUMP_DIR}/312_ToT_Fit_Compare" "{LOG_DIR}/312_ToT_Fit_Compare"
        "{input.exe}" -f "{input.cfg}" -o "{output.out}" > "{log}" 2>&1
        """

# ----- 317: main rule --------------------------------------------------------
rule ADC_ToT_Fit_Compare_317:
    """
    Build one ADC-ToT fit comparison result from its output filename.
    The corresponding config is looked up by OUTFILE_TO_CFG_317[outfile].
    """
    input:
        exe  = f"{BIN_DIR}/317_ADC_ToT_Fit_Compare",
        cfg  = lambda wc: OUTFILE_TO_CFG_317[f"{wc.outfile}.root"],
        runs = lambda wc: _inputs_for_outfile_317(f"{wc.outfile}.root")
    output:
        out  = f"{DUMP_DIR}/317_ADC_ToT_Fit_Compare/{{outfile}}.root"
    log:
        f"{LOG_DIR}/317_ADC_ToT_Fit_Compare/{{outfile}}.log"
    shell:
        r"""
        mkdir -p "{DUMP_DIR}/317_ADC_ToT_Fit_Compare" "{LOG_DIR}/317_ADC_ToT_Fit_Compare"
        "{input.exe}" -f "{input.cfg}" -o "{output.out}" > "{log}" 2>&1
        """

# ----- 106: convenience target ----------------------------------------------
rule run_all_106:
    """
    Run all configs starting with 106_*.json in config/.
    """
    input: ALL_106_OUTPUTS

# ----- 305: convenience target ----------------------------------------------
rule run_all_305:
    """
    Run all configs starting with 305_*.json in config/.
    """
    input: ALL_305_OUTPUTS

# ----- 312: convenience target ----------------------------------------------
rule run_all_312:
    """
    Run all configs starting with 312_*.json in config/.
    """
    input: ALL_312_OUTPUTS